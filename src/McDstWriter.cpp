#include "McDstWriter.h"

#include "McEvent.h"
#include "McParticle.h"
#include "McRun.h"
#include <TClonesArray.h>
#include <TFile.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <limits>
#include <system_error>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

struct McDstWriter::Impl {
    std::string destination, temporary, definitions, error;
    std::unique_ptr<TFile> file;
    std::unique_ptr<TTree> tree;
    std::unique_ptr<TClonesArray> events, particles;
    TClonesArray* event_array = nullptr;
    TClonesArray* particle_array = nullptr;
    int ensemble = 0, initial = 0, spectators = 0, final_count = 0;
    int participants = 0, binary = 0, nn_elastic = 0, nn_inelastic = 0;
    std::vector<int> oscar_ids;
    bool finished = false;

    ~Impl() {
        tree.reset();
        if (file)
            file->Close();
        if (!temporary.empty()) {
            std::error_code ignored;
            fs::remove(temporary, ignored);
        }
    }

    bool Fail(const std::string& message) {
        error = message;
        return false;
    }
};

McDstWriter::McDstWriter(const std::string& filename, const std::string& definitions)
    : impl_(std::make_unique<Impl>()) {
    auto& w = *impl_;
    w.destination = filename;
    w.definitions = definitions;
    if (filename.empty() || fs::is_directory(filename)) {
        w.Fail("output must be a file path");
        return;
    }
    const std::string pattern = filename + ".tmp.XXXXXX";
    std::vector<char> temporary_name(pattern.begin(), pattern.end());
    temporary_name.push_back('\0');
    const int descriptor = mkstemp(temporary_name.data());
    if (descriptor < 0) {
        w.Fail("cannot create output beside " + filename + ": " + std::strerror(errno));
        return;
    }
    w.temporary = temporary_name.data();
    if (close(descriptor) != 0) {
        w.Fail("cannot close temporary output: " + std::string(std::strerror(errno)));
        return;
    }
    w.file.reset(TFile::Open(w.temporary.c_str(), "RECREATE", "Converted McDst"));
    if (!w.file || w.file->IsZombie()) {
        w.Fail("cannot create output beside " + filename);
        return;
    }
    w.file->SetCompressionLevel(1);
    w.tree = std::make_unique<TTree>("McDst", "Converted particle data");
    w.tree->SetDirectory(w.file.get());
    w.events = std::make_unique<TClonesArray>("McEvent", 1);
    w.particles = std::make_unique<TClonesArray>("McParticle", 1000);
    w.event_array = w.events.get();
    w.particle_array = w.particles.get();
    auto* event_branch = w.tree->Branch("Event", &w.event_array, 65536, 99);
    auto* particle_branch = w.tree->Branch("Particle", &w.particle_array, 65536, 99);
    if (!event_branch || !particle_branch) {
        w.Fail("cannot create McDst branches");
        return;
    }
    event_branch->SetAutoDelete(kFALSE);
    particle_branch->SetAutoDelete(kFALSE);
    w.tree->Branch("OscarEnsemble", &w.ensemble);
    w.tree->Branch("InitialNucleons", &w.initial);
    w.tree->Branch("SpectatorNucleons", &w.spectators);
    w.tree->Branch("FinalParticles", &w.final_count);
    w.tree->Branch("TransportParticipants", &w.participants);
    w.tree->Branch("BinaryCollisions", &w.binary);
    w.tree->Branch("NNElasticCollisions", &w.nn_elastic);
    w.tree->Branch("NNInelasticCollisions", &w.nn_inelastic);
    w.tree->Branch("OscarParticleId", &w.oscar_ids);
}

McDstWriter::~McDstWriter() = default;

const std::string& McDstWriter::Error() const {
    return impl_->error;
}

bool McDstWriter::WriteEvent(const OscarEvent& source) {
    auto& w = *impl_;
    if (!w.error.empty())
        return false;
    if (w.finished)
        return w.Fail("writer already finished");
    if (source.particles.size() > std::numeric_limits<unsigned short>::max())
        return w.Fail("too many particles for McParticle indices");
    // Delete invokes destructors before placement new, including McEvent's TString.
    w.events->Delete();
    w.particles->Delete();
    w.oscar_ids.clear();
    int child[2] = {-1, -1};
    for (std::size_t i = 0; i < source.particles.size(); ++i) {
        const auto& p = source.particles[i];
        new ((*w.particles)[static_cast<int>(i)])
            McParticle(static_cast<int>(i), p.pdg, p.spectator ? 1 : 0, -1, 0, -1, -1, child, p.px,
                       p.py, p.pz, p.p0, p.x, p.y, p.z, p.t);
        w.oscar_ids.push_back(p.id);
    }
    auto* event = new ((*w.events)[0]) McEvent();
    event->setEventNr(source.number);
    event->setB(source.impact_parameter);
    event->setPhi(0.);
    event->setNes(1);
    event->setStepNr(1);
    event->setStepT(source.final_time);
    // McEvent uses Short_t; exact counters are also stored in Int_t branches.
    const int limit = std::numeric_limits<short>::max();
    event->setNpart(source.participants <= limit ? source.participants : -1);
    event->setNcoll(source.binary_collisions <= limit ? source.binary_collisions : -1);
    event->setComment(w.definitions.c_str());
    w.ensemble = source.ensemble;
    w.initial = source.initial_nucleons;
    w.spectators = source.spectators;
    w.final_count = static_cast<int>(source.particles.size());
    w.participants = source.participants;
    w.binary = source.binary_collisions;
    w.nn_elastic = source.nn_elastic_collisions;
    w.nn_inelastic = source.nn_inelastic_collisions;
    if (w.tree->Fill() < 0 || w.file->TestBit(TFile::kWriteError))
        return w.Fail("failed to write event");
    return true;
}

bool McDstWriter::Finish(long long input_events, long long skipped_events) {
    auto& w = *impl_;
    if (!w.error.empty())
        return false;
    if (w.finished)
        return w.Fail("writer already finished");
    w.file->cd();
    const auto entries = w.tree->GetEntries();
    const int run_count =
        entries <= std::numeric_limits<int>::max() ? static_cast<int>(entries) : 0;
    // OSCAR alone does not describe the beam setup or generated cross section.
    McRun run("SMASH", "OSCAR conversion; beam setup and cross section unspecified", 0, 0, 0., 0, 0,
              0., 0., 0., 0, 0., 0., 0., run_count);
    TNamed definitions("ConverterDefinitions", w.definitions.c_str());
    TParameter<Long64_t> input_count("InputEvents", input_events);
    TParameter<Long64_t> skipped_count("SkippedEmptyEvents", skipped_events);
    TParameter<Long64_t> output_count("OutputEvents", entries);
    if (run.Write() <= 0 || definitions.Write() <= 0 || input_count.Write() <= 0 ||
        skipped_count.Write() <= 0 || output_count.Write() <= 0 || w.file->Write() <= 0 ||
        w.file->TestBit(TFile::kWriteError))
        return w.Fail("failed to finalize ROOT output");
    w.tree.reset();
    w.file->Close();
    if (w.file->TestBit(TFile::kWriteError))
        return w.Fail("failed to close ROOT output");
    w.file.reset();
    std::error_code error;
    fs::rename(w.temporary, w.destination, error);
    if (error)
        return w.Fail("cannot replace output: " + error.message());
    w.temporary.clear();
    w.finished = true;
    return true;
}
