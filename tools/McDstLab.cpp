#include "McDstLab.h"
#include "McEvent.h"
#include "McParticle.h"
#include "McRun.h"

#include <TBranchElement.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THashList.h>
#include <TKey.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TParticlePDG.h>
#include <TROOT.h>
#include <TTree.h>

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string_view>
#include <thread>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace mcdst_lab {
namespace {
void Require(bool condition, std::string_view message) {
    if (!condition)
        throw std::runtime_error(std::string(message));
}

// Destruction order keeps ROOT objects closed before removing a failed output.
struct TemporaryFile {
    std::string name;
    explicit TemporaryFile(const std::string& destination) : name(destination + ".tmp.XXXXXX") {
        const int fd = mkstemp(name.data());
        Require(fd >= 0, "cannot create temporary output: " + std::string(std::strerror(errno)));
        if (close(fd) != 0) {
            std::filesystem::remove(name);
            throw std::runtime_error("cannot close temporary output");
        }
    }
    ~TemporaryFile() {
        std::error_code ignored;
        if (!name.empty())
            std::filesystem::remove(name, ignored);
    }
};

void Write(TObject& object, const char* name = nullptr) {
    Require(object.Write(name, TObject::kOverwrite) > 0,
            "cannot write ROOT object: " + std::string(object.GetName()));
}

template <class T> void Parameter(const char* name, T value) {
    TParameter<T> parameter(name, value);
    Write(parameter);
}

void CheckFloat(double value) {
    Require(std::isfinite(value) && std::abs(value) <= std::numeric_limits<float>::max(),
            "non-finite value or value outside McDst Float_t range");
}

struct Histogram {
    TH2D values;
    Long64_t undefined = 0, outside = 0, axial = 0;
    Histogram(const char* name, const char* title, const Options& options)
        : values(name, title, 240, -options.eta_max, options.eta_max, 150, 0., options.pt_max) {
        values.SetDirectory(nullptr);
    }
    void Fill(double longitudinal, double pt) {
        double eta;
        if (pt == 0.) {
            if (longitudinal == 0.) {
                ++undefined;
                return;
            }
            ++axial;
            // Represent infinite eta in the proper flow bin, without INF in ROOT statistics.
            eta = std::copysign(values.GetXaxis()->GetXmax() + 1., longitudinal);
        } else {
            eta = std::asinh(longitudinal / pt);
        }
        const int x = values.GetXaxis()->FindFixBin(eta);
        const int y = values.GetYaxis()->FindFixBin(pt);
        if (x == 0 || x > values.GetNbinsX() || y == 0 || y > values.GetNbinsY())
            ++outside;
        values.Fill(eta, pt);
    }
    void Save() {
        Write(values);
        const std::string name = values.GetName();
        Parameter((name + "_undefined_eta").c_str(), undefined);
        Parameter((name + "_outside_plot").c_str(), outside);
        Parameter((name + "_axial").c_str(), axial);
    }
};

// -1: unknown PDG, 0: not a charged hadron, 1: charged hadron/nuclear fragment.
const std::unordered_map<int, int>& KnownCharges() {
    // Immutable after initialization: no shared lazy PDG lookups in worker threads.
    static const auto charges = [] {
        auto* database = TDatabasePDG::Instance();
        database->GetParticle(211);
        std::unordered_map<int, int> result;
        TIter particles(database->ParticleList());
        while (auto* particle = dynamic_cast<TParticlePDG*>(particles())) {
            const std::string type = particle->ParticleClass();
            const bool hadron = type.find("Meson") != std::string::npos ||
                                type.find("Baryon") != std::string::npos || type == "Ion" ||
                                type == "Nucleus";
            result.emplace(particle->PdgCode(), hadron && particle->Charge() != 0. ? 1 : 0);
        }
        return result;
    }();
    return charges;
}

void Initialize() {
    static std::once_flag initialized;
    std::call_once(initialized, [] {
        ROOT::EnableThreadSafety();
        // Never split a large output into untracked files. Set once, before workers start.
        TTree::SetMaxTreeSize(std::numeric_limits<Long64_t>::max());
        McEvent::Class();
        McParticle::Class();
        McRun::Class();
        KnownCharges();
    });
}

int ChargedHadron(int pdg) {
    const long long code = std::abs(static_cast<long long>(pdg));
    if (code >= 1000000000LL && code < 1100000000LL) {
        const int z = (code / 10000) % 1000, a = (code / 10) % 1000;
        if (a > 0 && z <= a)
            return z > 0 ? 1 : 0;
    }
    const auto& known = KnownCharges();
    const auto found = known.find(pdg);
    return found == known.end() ? -1 : found->second;
}

void CheckArray(TTree& tree, const char* name, const char* type) {
    auto* branch = tree.GetBranch(name);
    auto* element = dynamic_cast<TBranchElement*>(branch);
    // Unsplit TClonesArray branches may not expose the element class until read.
    Require(branch && std::string(branch->GetClassName()) == "TClonesArray" &&
                (!element || !element->GetClonesName() || !*element->GetClonesName() ||
                 std::string(element->GetClonesName()) == type),
            std::string("expected TClonesArray<") + type + "> branch " + name);
}

std::string OutputFilename(const std::filesystem::path& input, bool hist_only) {
    std::string stem = input.stem().string();
    constexpr std::string_view suffix = ".mcdst";
    if (stem.size() >= suffix.size()) {
        std::string ending = stem.substr(stem.size() - suffix.size());
        std::transform(ending.begin(), ending.end(), ending.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        if (ending == suffix)
            stem.resize(stem.size() - suffix.size());
    }
    return stem + (hist_only ? ".qa.root" : ".McDst.root");
}

} // namespace

Frame::Frame(const Options& options) {
    Require(std::isfinite(options.sqrtsnn) && std::isfinite(options.nucleon_mass) &&
                options.nucleon_mass > 0. && options.sqrtsnn > 2. * options.nucleon_mass,
            "sqrtsnn must exceed twice the positive nucleon mass (GeV)");
    gamma = options.sqrtsnn / (2. * options.nucleon_mass);
    beta = std::sqrt(1. - 1. / (gamma * gamma));
    Require(std::isfinite(gamma) && beta < 1., "beam energy is too large for this boost");
}

double Frame::Energy(double energy, double pz) const {
    return gamma * (energy + beta * pz);
}

double Frame::Longitudinal(double energy, double pz) const {
    return -gamma * (pz + beta * energy);
}

Result Process(const std::string& input_name, const std::string& output_name,
               const Options& options) {
    Initialize();
    const Frame frame(options);
    Require(std::isfinite(options.eta_max) && options.eta_max > 0. && options.eta_max <= 100. &&
                std::isfinite(options.pt_max) && options.pt_max > 0.,
            "plot limits must be positive and finite, eta-max <= 100");
    Require(!output_name.empty() && !std::filesystem::is_directory(output_name),
            "output must be a file path");
    Require(std::filesystem::weakly_canonical(input_name) !=
                    std::filesystem::weakly_canonical(output_name) &&
                (!std::filesystem::exists(output_name) ||
                 !std::filesystem::equivalent(input_name, output_name)),
            "input and output must be different files");
    std::unique_ptr<TFile> input(TFile::Open(input_name.c_str(), "READ"));
    Require(input && !input->IsZombie() && !input->TestBit(TFile::kRecovered),
            "cannot read a complete input ROOT file");
    Require(!input->Get("LabFrame"), "input was already processed by mcdst_lab");
    auto* source = dynamic_cast<TTree*>(input->Get("McDst"));
    Require(source != nullptr, "missing McDst tree");
    CheckArray(*source, "Particle", "McParticle");
    CheckArray(*source, "Event", "McEvent");

    TClonesArray particles("McParticle"), events("McEvent");
    auto* particle_ptr = &particles;
    auto* event_ptr = &events;
    Require(source->SetBranchAddress("Particle", &particle_ptr) >= 0 &&
                source->SetBranchAddress("Event", &event_ptr) >= 0,
            "cannot bind McDst branches");
    source->GetBranch("Particle")->SetAutoDelete(kFALSE);
    source->GetBranch("Event")->SetAutoDelete(kFALSE);
    // Reset before the arrays leave scope, including on failure.
    struct ResetAddresses {
        TTree* tree;
        ~ResetAddresses() {
            tree->ResetBranchAddresses();
        }
    } reset{source};
    if (options.hist_only) {
        source->SetBranchStatus("*", false);
        source->SetBranchStatus("Particle*", true);
        source->SetBranchStatus("Event*", true);
    }
    source->SetCacheSize(8 * 1024 * 1024);

    TemporaryFile temporary(output_name);
    std::unique_ptr<TFile> output(TFile::Open(temporary.name.c_str(), "RECREATE"));
    Require(output && !output->IsZombie(), "cannot create output ROOT file");
    output->SetCompressionLevel(1);
    std::unique_ptr<TTree> destination;
    if (!options.hist_only) {
        destination.reset(source->CloneTree(0));
        Require(destination != nullptr, "cannot clone McDst structure");
        destination->SetDirectory(output.get());
        destination->SetAutoFlush(64);
        destination->SetMaxVirtualSize(4 * 1024 * 1024);
    }

    TTree multiplicity("LabMultiplicity", "Counts without acceptance or spectator cuts");
    multiplicity.SetDirectory(output.get());
    multiplicity.SetAutoFlush(1024);
    Long64_t entry = 0, total_particles = 0, total_unknown = 0, masked_positions = 0;
    int event_number = 0, n_all = 0, n_ch = 0, n_unknown = 0;
    multiplicity.Branch("InputEntry", &entry);
    multiplicity.Branch("EventNr", &event_number);
    multiplicity.Branch("NAll", &n_all);
    multiplicity.Branch("NCh", &n_ch);
    multiplicity.Branch("NUnknownPDG", &n_unknown);
    Histogram before("hEtaPtBefore", "CM;#eta (axis opposite Projectile);p_{T} [GeV/c]", options);
    Histogram after("hEtaPtAfter",
                    "Target rest frame;#eta (axis opposite Projectile);p_{T} [GeV/c]", options);
    TH1D all_hist("hNAll", "All final particles;N_{all};Events", 65536, -.5, 65535.5);
    TH1D charged_hist("hNCh", "Charged final hadrons;N_{ch};Events", 65536, -.5, 65535.5);
    all_hist.SetDirectory(nullptr);
    charged_hist.SetDirectory(nullptr);
    std::unordered_map<int, int> charge_cache;
    const Long64_t entries = source->GetEntries();
    for (entry = 0; entry < entries; ++entry) {
        Require(source->GetEntry(entry) > 0, "failed to read entry " + std::to_string(entry));
        Require(particles.GetClass() == McParticle::Class() &&
                    events.GetClass() == McEvent::Class(),
                "unexpected class in McDst arrays");
        Require(events.GetEntriesFast() == 1, "expected one McEvent per entry");
        auto* event = static_cast<McEvent*>(events.At(0));
        Require(event != nullptr, "null McEvent");
        event_number = event->eventNr();
        n_all = particles.GetEntriesFast();
        n_ch = n_unknown = 0;
        for (int i = 0; i < n_all; ++i) {
            auto* p = static_cast<McParticle*>(particles.At(i));
            Require(p != nullptr, "null McParticle at entry " + std::to_string(entry));
            const double e = p->e(), px = p->px(), py = p->py(), pz = p->pz();
            for (double value : {e, px, py, pz})
                CheckFloat(value);
            Require(e > 0., "nonpositive particle energy");
            const double transformed_e = frame.Energy(e, pz);
            const double transformed_pz = frame.Longitudinal(e, pz);
            CheckFloat(transformed_e);
            CheckFloat(transformed_pz);
            Require(transformed_e > 0., "nonpositive transformed energy");
            const double pt = std::hypot(px, py);
            before.Fill(-pz, pt);
            after.Fill(transformed_pz, pt);
            auto found = charge_cache.find(p->pdg());
            if (found == charge_cache.end())
                found = charge_cache.emplace(p->pdg(), ChargedHadron(p->pdg())).first;
            n_ch += found->second == 1;
            n_unknown += found->second == -1;

            // A partially masked four-position cannot be boosted. Preserve it verbatim.
            const double x = p->x(), y = p->y(), z = p->z(), t = p->t();
            const bool masked = x == -999. || y == -999. || z == -999. || t == -999.;
            masked_positions += masked;
            if (!options.hist_only) {
                p->setMomentum(-px, py, transformed_pz, transformed_e);
                if (!masked) {
                    const double new_t = frame.Energy(t, z), new_z = frame.Longitudinal(t, z);
                    for (double value : {x, y, new_z, new_t})
                        CheckFloat(value);
                    p->setPosition(-x, y, new_z, new_t);
                }
            }
        }
        if (destination) {
            // Rotate the reaction-plane direction with the transverse axes.
            if (event->phi() != -999.) {
                Require(std::isfinite(event->phi()), "non-finite event phi");
                const double pi = std::acos(-1.);
                double phi = std::fmod(pi - event->phi(), 2. * pi);
                if (phi < 0.)
                    phi += 2. * pi;
                event->setPhi(phi);
            }
            Require(destination->Fill() >= 0, "cannot write transformed event");
        }
        // Multiplicity is invariant: a single set of counts applies to both frames.
        total_particles += n_all;
        total_unknown += n_unknown;
        all_hist.Fill(n_all);
        charged_hist.Fill(n_ch);
        Require(multiplicity.Fill() >= 0 && !output->TestBit(TFile::kWriteError),
                "cannot write event counts");
    }
    output->cd();
    if (destination)
        Write(*destination);
    Write(multiplicity);
    before.Save();
    after.Save();
    Write(all_hist);
    Write(charged_hist);
    TNamed description(
        "LabFrame",
        "v1: symmetric SMASH CM -> Target rest; x'=-x,y'=y,z'=-z after boost; "
        "all saved particles, no cuts; NCh=charged hadrons including nuclear fragments; "
        "NUnknownPDG>0 means NCh is incomplete; masked positions unchanged; "
        "Event.stepT remains generator CM time; extra branches retain source meanings; "
        "source file-level metadata is in SourceMetadata");
    Write(description);
    Parameter("LabSqrtsNN", options.sqrtsnn);
    Parameter("LabNucleonMass", options.nucleon_mass);
    Parameter("LabBeta", frame.beta);
    Parameter("LabGamma", frame.gamma);
    Parameter("LabEvents", entries);
    Parameter("LabParticles", total_particles);
    Parameter("LabUnknownPDG", total_unknown);
    Parameter("LabMaskedPositions", masked_positions);
    Parameter("LabHistOnly", static_cast<int>(options.hist_only));
    TNamed source_name("LabSource", std::filesystem::absolute(input_name).c_str());
    Write(source_name);

    // Copy small run metadata, but never unrelated input trees or old histograms.
    // Those retain their CM interpretation and cannot be relabeled as lab observables.
    auto* source_metadata = output->mkdir("SourceMetadata");
    Require(source_metadata != nullptr, "cannot create source metadata directory");
    TIter keys(input->GetListOfKeys());
    std::unordered_set<std::string> copied;
    while (auto* key = dynamic_cast<TKey*>(keys())) {
        const std::string class_name = key->GetClassName();
        if (class_name != "McRun" && class_name != "TNamed" &&
            class_name.rfind("TParameter<", 0) != 0)
            continue;
        if (!copied.insert(key->GetName()).second)
            continue;
        std::unique_ptr<TObject> object(key->ReadObj());
        Require(object != nullptr, "cannot read source metadata");
        source_metadata->cd();
        Write(*object, key->GetName());
        if (!options.hist_only) {
            if (auto* run = dynamic_cast<McRun*>(object.get())) {
                // McRun stores momentum magnitudes per nucleon; axis sign is in LabFrame.
                const double beam_p = 2. * frame.gamma *
                                      std::sqrt(std::pow(options.sqrtsnn / 2., 2) -
                                                std::pow(options.nucleon_mass, 2));
                CheckFloat(beam_p);
                run->setPProj(beam_p);
                run->setPTarg(0.);
                output->cd();
                Write(*run, key->GetName());
            }
        }
    }
    output->cd();
    // Detach owned trees before closing the file; avoid ROOT double ownership.
    multiplicity.SetDirectory(nullptr);
    destination.reset();
    output->Close();
    Require(!output->TestBit(TFile::kWriteError), "failed to finalize ROOT output");
    output.reset();
    std::filesystem::rename(temporary.name, output_name);
    temporary.name.clear();
    return {entries, total_particles, total_unknown};
}

BatchResult ProcessDirectory(
    const std::string& input_name, 
    const std::string& output_name,
    const Options& options, 
    unsigned threads,
    std::size_t max_files) {
    const Frame frame(options);
    Require(std::filesystem::is_directory(input_name), "input must be a directory");
    const auto input = std::filesystem::canonical(input_name);
    const auto output = std::filesystem::weakly_canonical(output_name);
    Require(!output_name.empty() && input != output, "input and output directories must differ");
    if (std::filesystem::exists(output))
        Require(std::filesystem::is_directory(output), "output must be a directory");
    std::vector<std::filesystem::path> files;
    for (const auto& file : std::filesystem::directory_iterator(input)) {
        if (file.is_regular_file() && file.path().extension() == ".root")
            files.push_back(file.path());
    }
    Require(!files.empty(), "no .root files in input directory (scan is not recursive)");
    std::sort(files.begin(), files.end());

    if (files.size() > max_files) files.resize(max_files);

    std::vector<std::filesystem::path> destinations;
    std::unordered_set<std::string> output_names;
    for (const auto& file : files) {
        const auto filename = OutputFilename(file, options.hist_only);
        Require(output_names.insert(filename).second,
                "multiple input files map to the same output filename: " + filename);
        destinations.push_back(output / filename);
    }
    std::filesystem::create_directories(output);
    Initialize();
    const unsigned available =
        threads ? threads : std::max(1u, std::thread::hardware_concurrency());
    BatchResult result;
    result.workers = static_cast<unsigned>(std::min<std::size_t>(available, files.size()));
    std::cout << "Files: " << files.size() << "; workers: " << result.workers
              << "; beta=" << frame.beta << "; gamma=" << frame.gamma << '\n';
    std::atomic<std::size_t> next{0};
    std::mutex reporting;
    auto work = [&] {
        while (true) {
            const std::size_t index = next.fetch_add(1);
            if (index >= files.size())
                return;
            try {
                const auto report =
                    Process(files[index].string(), destinations[index].string(), options);
                std::lock_guard<std::mutex> lock(reporting);
                ++result.succeeded;
                std::cout << files[index].filename().string() << " -> "
                          << destinations[index].filename().string() << ": " << report.events
                          << " events, " << report.particles << " particles";
                if (report.unknown_pdg)
                    std::cout << "; unknown PDG records=" << report.unknown_pdg
                              << " (NCh incomplete; see NUnknownPDG)";
                std::cout << '\n';
            } catch (const std::exception& error) {
                std::lock_guard<std::mutex> lock(reporting);
                ++result.failed;
                std::cerr << files[index].filename().string() << ": " << error.what() << '\n';
            }
        }
    };
    std::vector<std::thread> workers;
    // If the OS refuses another thread, join all already started workers before unwinding.
    try {
        for (unsigned i = 1; i < result.workers; ++i)
            workers.emplace_back(work);
    } catch (...) {
        for (auto& worker : workers)
            worker.join();
        throw;
    }
    work();
    for (auto& worker : workers)
        worker.join();
    std::cout << "Completed: " << result.succeeded << "; failed: " << result.failed << '\n';
    return result;
}
} // namespace mcdst_lab
