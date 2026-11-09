#include "McDstLab.h"
#include "McEvent.h"
#include "McParticle.h"
#include "McRun.h"

#include <TClonesArray.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TNamed.h>
#include <TParameter.h>
#include <TTree.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

namespace fs = std::filesystem;
namespace {
void Check(bool condition, const std::string& message) {
    if (!condition)
        throw std::runtime_error(message);
}
void Near(double actual, double expected, double tolerance, const std::string& message) {
    Check(std::abs(actual - expected) <= tolerance, message);
}
struct Temp {
    fs::path path =
        fs::temp_directory_path() / ("mcdst-lab-" + std::to_string(std::random_device{}()));
    Temp() {
        fs::create_directories(path);
    }
    ~Temp() {
        std::error_code error;
        fs::remove_all(path, error);
    }
};
std::string Bytes(const fs::path& path) {
    std::ifstream file(path, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(file), {});
}
template <class T> T* Get(TFile& file, const char* name) {
    auto* object = dynamic_cast<T*>(file.Get(name));
    Check(object != nullptr, std::string("missing object: ") + name);
    return object;
}

// Mixed species, an empty event, a 394-particle event, then one particle.
// Extra scalar/vector branches exercise the CloneTree address propagation.
void Fixture(const fs::path& filename, int repetitions = 1, int split = 99, bool bad = false) {
    TFile file(filename.c_str(), "RECREATE");
    TClonesArray particles("McParticle"), events("McEvent");
    auto* particle_pointer = &particles;
    auto* event_pointer = &events;
    TTree tree("McDst", "test CM input");
    tree.Branch("Particle", &particle_pointer, 32000, split);
    tree.Branch("Event", &event_pointer, 32000, split);
    int tag = 0;
    std::vector<int> ids;
    tree.Branch("UserTag", &tag);
    tree.Branch("OscarParticleId", &ids);
    const double beam = std::sqrt(1.5 * 1.5 - .938 * .938);
    auto add = [&](int pdg, double mass, double px, double py, double pz, int mask = 0) {
        const int index = particles.GetEntriesFast();
        int child[2] = {5, 7};
        double energy = std::sqrt(mass * mass + px * px + py * py + pz * pz);
        if (bad && index == 3)
            energy = std::numeric_limits<double>::quiet_NaN();
        new (particles[index]) McParticle(index, pdg, index % 2, 3, 2, 4, 0, child, px, py, pz,
                                          energy, mask ? -999. : 1., mask ? -999. : 2.,
                                          mask ? -999. : 3., mask == 1 ? -999. : 10.);
        ids.push_back(1000 + index);
    };
    for (int repeat = 0; repeat < repetitions; ++repeat) {
        for (int type = 0; type < 4; ++type) {
            events.Delete();
            particles.Delete();
            ids.clear();
            auto* event = new (events[0]) McEvent();
            tag = repeat * 4 + type;
            event->setEventNr(tag + 100);
            event->setB(7.5);
            event->setStepT(200.);
            event->setPhi(.3);
            event->setNpart(12);
            event->setNcoll(19);
            event->setComment("source event comment");
            if (type == 0) {
                add(2212, .938, 0., 0., -beam, 1);
                add(2212, .938, 0., 0., beam, 2);
                add(-211, .138, .3, .4, .8);
                add(2112, .938, 0., 0., 0.);
                add(11, .000511, .2, 0., .1);
                add(22, 0., 0., 0., 1.);
                add(1000010020, 1.876, .1, 0., -.2);
                add(987654321, .5, .2, 0., -.3);
            } else if (type == 2) {
                for (int i = 0; i < 394; ++i)
                    add(2212, .938, .05 * (i % 5), .03 * (i % 7), i % 2 ? beam : -beam);
            } else if (type == 3) {
                add(211, .138, .2, .3, -.1);
            }
            tree.Fill();
        }
    }
    tree.Write();
    McRun run("SMASH", "CM fixture", 197, 79, beam, 197, 79, beam, 0., 16., 0, 0., 0., 0.,
              4 * repetitions);
    run.Write();
    TNamed metadata("Definitions", "test source metadata");
    metadata.Write();
}

void CheckHistograms(TFile& file, int repetitions = 1) {
    const int particles = 403 * repetitions;
    auto* before = Get<TH2>(file, "hEtaPtBefore");
    auto* after = Get<TH2>(file, "hEtaPtAfter");
    Check(std::string(before->GetXaxis()->GetTitle()).find("eta") != std::string::npos &&
              std::string(before->GetYaxis()->GetTitle()).find("p_{T}") != std::string::npos,
          "eta and pt axes were swapped");
    for (const char* name : {"hEtaPtBefore", "hEtaPtAfter"}) {
        auto* hist = Get<TH2>(file, name);
        const auto undefined =
            Get<TParameter<Long64_t>>(file, (std::string(name) + "_undefined_eta").c_str())
                ->GetVal();
        Near(hist->Integral(0, hist->GetNbinsX() + 1, 0, hist->GetNbinsY() + 1) + undefined,
             particles, 0., "histogram lost or duplicated particles");
    }
    Check(Get<TParameter<Long64_t>>(file, "hEtaPtBefore_undefined_eta")->GetVal() == repetitions,
          "zero momentum particle incorrectly has finite eta");
    Check(Get<TParameter<Long64_t>>(file, "LabParticles")->GetVal() == particles,
          "wrong total particle count");
    Check(Get<TParameter<Long64_t>>(file, "LabMaskedPositions")->GetVal() == 2 * repetitions,
          "masked positions not accounted for");
    const auto* hist = Get<TH1>(file, "hNAll");
    Check(hist->GetBinContent(hist->FindFixBin(394.)) == repetitions,
          "394-particle event was filtered again");
    auto* counts = Get<TTree>(file, "LabMultiplicity");
    Check(counts->GetEntries() == 4 * repetitions, "empty event lost from multiplicity");
    int all = 0, charged = 0, unknown = 0, event = 0;
    Long64_t input_entry = 0;
    counts->SetBranchAddress("NAll", &all);
    counts->SetBranchAddress("NCh", &charged);
    counts->SetBranchAddress("NUnknownPDG", &unknown);
    counts->SetBranchAddress("EventNr", &event);
    counts->SetBranchAddress("InputEntry", &input_entry);
    const int expected_all[] = {8, 0, 394, 1}, expected_ch[] = {4, 0, 394, 1};
    for (int i = 0; i < 4 * repetitions; ++i) {
        counts->GetEntry(i);
        Check(all == expected_all[i % 4] && charged == expected_ch[i % 4] &&
                  unknown == (i % 4 == 0) && event == i + 100 && input_entry == i,
              "incorrect multiplicity or charged-hadron classification");
    }
    counts->ResetBranchAddresses();
    // The isolated pion must be in bins for pseudorapidity, not rapidity, with reversed axis.
    const double pt = .5, e = std::sqrt(.138 * .138 + pt * pt + .8 * .8);
    const mcdst_lab::Frame frame({});
    Check(before->GetBinContent(before->FindFixBin(std::asinh(-.8 / pt), pt)) >= repetitions,
          "CM eta does not use the opposite beam axis");
    Check(after->GetBinContent(after->FindFixBin(std::asinh(frame.Longitudinal(e, .8) / pt), pt)) >=
              repetitions,
          "lab eta bin is wrong");
}

void CheckTransformed(const fs::path& input_path, const fs::path& output_path) {
    TFile input(input_path.c_str(), "READ"), output(output_path.c_str(), "READ");
    CheckHistograms(output);
    auto* source = Get<TTree>(input, "McDst");
    auto* target = Get<TTree>(output, "McDst");
    Check(source->GetEntries() == target->GetEntries() &&
              source->GetListOfBranches()->GetEntries() ==
                  target->GetListOfBranches()->GetEntries(),
          "tree shape changed");
    TClonesArray in_particles("McParticle"), out_particles("McParticle"), in_events("McEvent"),
        out_events("McEvent");
    auto* ip = &in_particles;
    auto* op = &out_particles;
    auto* ie = &in_events;
    auto* oe = &out_events;
    source->SetBranchAddress("Particle", &ip);
    target->SetBranchAddress("Particle", &op);
    source->SetBranchAddress("Event", &ie);
    target->SetBranchAddress("Event", &oe);
    for (auto* tree : {source, target}) {
        tree->GetBranch("Particle")->SetAutoDelete(kFALSE);
        tree->GetBranch("Event")->SetAutoDelete(kFALSE);
    }
    int source_tag = -1, target_tag = -2;
    std::vector<int> source_ids, target_ids;
    auto* si = &source_ids;
    auto* ti = &target_ids;
    source->SetBranchAddress("UserTag", &source_tag);
    target->SetBranchAddress("UserTag", &target_tag);
    source->SetBranchAddress("OscarParticleId", &si);
    target->SetBranchAddress("OscarParticleId", &ti);
    const mcdst_lab::Frame frame({});
    for (int entry = 0; entry < 4; ++entry) {
        source->GetEntry(entry);
        target->GetEntry(entry);
        Check(source_tag == target_tag && source_ids == target_ids, "extra branches changed");
        const auto* a_event = static_cast<McEvent*>(in_events.At(0));
        const auto* b_event = static_cast<McEvent*>(out_events.At(0));
        Near(b_event->phi(), std::acos(-1.) - a_event->phi(), 1e-6,
             "reaction plane rotation wrong");
        Check(a_event->eventNr() == b_event->eventNr() && a_event->b() == b_event->b() &&
                  b_event->stepT() == 200. && b_event->npart() == 12 && b_event->ncoll() == 19,
              "event metadata changed");
        Check(in_particles.GetEntriesFast() == out_particles.GetEntriesFast(),
              "particle count changed");
        for (int i = 0; i < in_particles.GetEntriesFast(); ++i) {
            const auto* a = static_cast<McParticle*>(in_particles.At(i));
            const auto* b = static_cast<McParticle*>(out_particles.At(i));
            Check(a->pdg() == b->pdg() && a->index() == b->index() && a->status() == b->status() &&
                      a->parent() == b->parent() && a->mate() == b->mate() &&
                      a->firstChild() == b->firstChild() && a->lastChild() == b->lastChild(),
                  "identity, spectator marker or ancestry changed");
            Near(a->momentum().M2(), b->momentum().M2(), 5e-6, "mass invariant failed");
            Near(a->pt(), b->pt(), 1e-7, "pt invariant failed");
            Near(b->px(), -a->px(), 1e-7, "transverse orientation wrong");
            Near(b->e(), frame.Energy(a->e(), a->pz()), 5e-7, "energy boost wrong");
            Near(b->pz(), frame.Longitudinal(a->e(), a->pz()), 5e-7, "longitudinal boost wrong");
            if (a->x() == -999.) {
                Check(a->x() == b->x() && a->y() == b->y() && a->z() == b->z() && a->t() == b->t(),
                      "sentinel coordinates were transformed");
            } else {
                Near(b->position().M2(), a->position().M2(), 5e-5, "spacetime interval changed");
                Near(b->z(), frame.Longitudinal(a->t(), a->z()), 2e-6, "position boost wrong");
            }
            if (entry == 0 && i == 0) {
                Near(b->pz(), 0., 2e-7, "Target is not stopped");
                Near(b->e(), .938, 2e-7, "Target rest energy wrong");
            }
            if (entry == 0 && i == 1)
                Near(b->pz(), -3.743720562065003, 5e-7, "Projectile momentum/direction wrong");
        }
    }
    Near(Get<McRun>(output, "run")->pTarg(), 0., 0., "run target momentum wrong");
    Near(Get<McRun>(output, "run")->pProj(), 3.743720562065003, 5e-7,
         "run projectile momentum wrong");
    Check(Get<McRun>(output, "SourceMetadata/run")->pTarg() > 0., "source run metadata lost");
    Get<TNamed>(output, "SourceMetadata/Definitions");
    source->ResetBranchAddresses();
    target->ResetBranchAddresses();
}

template <class Action> void Reject(Action action, const std::string& expected) {
    try {
        action();
    } catch (const std::exception& error) {
        Check(std::string(error.what()).find(expected) != std::string::npos,
              "wrong rejection: " + std::string(error.what()));
        return;
    }
    throw std::runtime_error("expected rejection: " + expected);
}
} // namespace

int main(int argc, char** argv) {
    try {
        if (argc >= 3 && std::string(argv[1]) == "--fixture") {
            Fixture(argv[2], argc > 3 ? std::stoi(argv[3]) : 1, argc > 4 ? std::stoi(argv[4]) : 99);
            return 0;
        }
        Temp temp;
        const auto input = temp.path / "input.root", output = temp.path / "output.root";
        Fixture(input);
        const auto input_bytes = Bytes(input);
        mcdst_lab::Process(input.string(), output.string());
        Check(Bytes(input) == input_bytes, "input modified");
        CheckTransformed(input, output);

        mcdst_lab::Options hist_only;
        hist_only.hist_only = true;
        const auto hist_path = temp.path / "hist.root";
        mcdst_lab::Process(input.string(), hist_path.string(), hist_only);
        {
            TFile hist(hist_path.c_str(), "READ"), full(output.c_str(), "READ");
            Check(!hist.Get("McDst"), "hist-only wrote a particle tree");
            CheckHistograms(hist);
            for (const char* name : {"hEtaPtBefore", "hEtaPtAfter", "hNAll", "hNCh"}) {
                auto* a = Get<TH1>(hist, name);
                auto* b = Get<TH1>(full, name);
                for (int bin = 0; bin < a->GetNcells(); ++bin)
                    Check(a->GetBinContent(bin) == b->GetBinContent(bin),
                          "hist-only differs from full output");
            }
        }
        const auto old_output = Bytes(output);
        const auto bad_input = temp.path / "bad.root";
        Fixture(bad_input, 1, 99, true);
        Reject([&] { mcdst_lab::Process(bad_input.string(), output.string()); }, "non-finite");
        Check(Bytes(output) == old_output, "failed processing replaced the existing output");
        Reject([&] { mcdst_lab::Process(input.string(), input.string()); }, "different files");
        Reject([&] { mcdst_lab::Process(output.string(), hist_path.string()); },
               "already processed");
        auto invalid = hist_only;
        invalid.sqrtsnn = 1.;
        Reject([&] { mcdst_lab::Process(input.string(), output.string(), invalid); }, "sqrtsnn");
        invalid.sqrtsnn = 3.;
        invalid.eta_max = 0.;
        Reject([&] { mcdst_lab::Process(input.string(), output.string(), invalid); },
               "plot limits");
        Fixture(temp.path / "unsplit.root", 1, 0);
        mcdst_lab::Process((temp.path / "unsplit.root").string(), output.string());
        CheckTransformed(temp.path / "unsplit.root", output);
        Fixture(temp.path / "empty.root", 0);
        mcdst_lab::Process((temp.path / "empty.root").string(), output.string());
        {
            TFile empty(output.c_str(), "READ");
            Check(Get<TTree>(empty, "McDst")->GetEntries() == 0, "empty tree processing failed");
        }

        fs::create_directory(temp.path / "batch");
        fs::copy_file(input, temp.path / "batch/a.root");
        fs::copy_file(input, temp.path / "batch/b.root");
        fs::copy_file(input, temp.path / "batch/ignore.dat");
        fs::create_directory(temp.path / "batch/nested");
        fs::copy_file(input, temp.path / "batch/nested/ignore.root");
        auto report = mcdst_lab::ProcessDirectory((temp.path / "batch").string(),
                                                  (temp.path / "batch-out").string(), {}, 8);
        Check(report.workers == 2 && report.succeeded == 2 && report.failed == 0,
              "parallel directory scheduling failed");
        CheckTransformed(input, temp.path / "batch-out/a.root");
        CheckTransformed(input, temp.path / "batch-out/b.root");
        Check(!fs::exists(temp.path / "batch-out/ignore.dat"), "unexpected non-ROOT output");
        fs::copy_file(bad_input, temp.path / "batch/bad.root");
        report = mcdst_lab::ProcessDirectory((temp.path / "batch").string(),
                                             (temp.path / "batch-hist").string(), hist_only, 0);
        Check(report.succeeded == 2 && report.failed == 1, "batch failure handling failed");
        Reject(
            [&] {
                mcdst_lab::ProcessDirectory((temp.path / "batch").string(),
                                            (temp.path / "batch").string());
            },
            "differ");
        for (const auto& file : fs::recursive_directory_iterator(temp.path))
            Check(file.path().filename().string().find(".tmp.") == std::string::npos,
                  "staging file leaked");
        std::cout << "Lab kinematics, ROOT round trips, multiplicity, histogram and parallel I/O "
                     "checks passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "TEST FAILURE: " << error.what() << '\n';
        return 1;
    }
}
