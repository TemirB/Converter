#include "McDstLabInternal.h"
#include "McEvent.h"
#include "McParticle.h"
#include "McRun.h"

#include <TBranchElement.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TH2D.h>
#include <THashList.h>
#include <TKey.h>
#include <TNamed.h>
#include <TParticlePDG.h>
#include <TROOT.h>

#include <filesystem>
#include <mutex>
#include <unordered_map>
#include <unordered_set>

namespace mcdst_lab::detail {
namespace {

struct Charge {
    bool known = false;
    bool charged = false;
    bool hadron = false;
};

using ChargeTable = std::unordered_map<int, Charge>;

const ChargeTable& KnownCharges() {
    // Initialize once before workers start; ROOT PDG lookups never run concurrently.
    static const auto charges = [] {
        auto* database = TDatabasePDG::Instance();
        database->GetParticle(211);
        ChargeTable table;
        TIter particles(database->ParticleList());
        while (auto* particle = dynamic_cast<TParticlePDG*>(particles())) {
            const std::string type = particle->ParticleClass();
            const bool hadron = type.find("Meson") != std::string::npos ||
                                type.find("Baryon") != std::string::npos || type == "Ion" ||
                                type == "Nucleus";
            table.emplace(particle->PdgCode(), Charge{true, particle->Charge() != 0.0, hadron});
        }
        return table;
    }();
    return charges;
}

Charge ParticleCharge(int pdg) {
    const long long code = std::abs(static_cast<long long>(pdg));
    if (code >= 1000000000LL && code < 1100000000LL) {
        const int charge = (code / 10000) % 1000;
        const int mass_number = (code / 10) % 1000;
        if (mass_number > 0 && charge <= mass_number) {
            return {true, charge > 0, true};
        }
    }
    const auto& table = KnownCharges();
    const auto found = table.find(pdg);
    return found == table.end() ? Charge{} : found->second;
}

void CheckArray(TTree& tree, const char* name, const char* type) {
    auto* branch = tree.GetBranch(name);
    auto* element = dynamic_cast<TBranchElement*>(branch);
    // Unsplit arrays may expose their element class only after the first read.
    const bool valid_element = !element || !element->GetClonesName() ||
                               !*element->GetClonesName() ||
                               std::string(element->GetClonesName()) == type;
    Require(branch && std::string(branch->GetClassName()) == "TClonesArray" && valid_element,
            std::string("expected TClonesArray<") + type + "> branch " + name);
}

std::unique_ptr<TFile> OpenInput(const std::string& path) {
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "READ"));
    Require(file && !file->IsZombie() && !file->TestBit(TFile::kRecovered),
            "cannot read a complete input ROOT file");
    Require(!file->Get("LabFrame"), "input was already processed by mcdst_lab");
    return file;
}

TTree& ReadMcDstTree(TFile& file) {
    auto* tree = dynamic_cast<TTree*>(file.Get("McDst"));
    Require(tree != nullptr, "missing McDst tree");
    CheckArray(*tree, "Particle", "McParticle");
    CheckArray(*tree, "Event", "McEvent");
    return *tree;
}

struct InputBuffers {
    TClonesArray particles{"McParticle"};
    TClonesArray events{"McEvent"};
    TClonesArray* particle_pointer = &particles;
    TClonesArray* event_pointer = &events;
    ResetTreeAddresses reset;

    InputBuffers(TTree& tree, bool hist_only) : reset(tree) {
        Require(tree.SetBranchAddress("Particle", &particle_pointer) >= 0 &&
                    tree.SetBranchAddress("Event", &event_pointer) >= 0,
                "cannot bind McDst branches");
        tree.GetBranch("Particle")->SetAutoDelete(kFALSE);
        tree.GetBranch("Event")->SetAutoDelete(kFALSE);
        if (hist_only) {
            tree.SetBranchStatus("*", false);
            tree.SetBranchStatus("Particle*", true);
            tree.SetBranchStatus("Event*", true);
        }
        tree.SetCacheSize(8 * 1024 * 1024);
    }

    McEvent& ReadEvent(TTree& tree, Long64_t entry) {
        Require(tree.GetEntry(entry) > 0, "failed to read entry " + std::to_string(entry));
        Require(particles.GetClass() == McParticle::Class() &&
                    events.GetClass() == McEvent::Class(),
                "unexpected class in McDst arrays");
        Require(events.GetEntriesFast() == 1, "expected one McEvent per entry");
        auto* event = static_cast<McEvent*>(events.At(0));
        Require(event != nullptr, "null McEvent");
        return *event;
    }
};

std::unique_ptr<TFile> OpenOutput(const std::string& path) {
    std::unique_ptr<TFile> file(TFile::Open(path.c_str(), "RECREATE"));
    Require(file && !file->IsZombie(), "cannot create output ROOT file");
    file->SetCompressionLevel(1);
    return file;
}

std::unique_ptr<TTree> CloneOutputTree(TTree& source, TFile& output, bool hist_only) {
    if (hist_only) {
        return nullptr;
    }
    std::unique_ptr<TTree> tree(source.CloneTree(0));
    Require(tree != nullptr, "cannot clone McDst structure");
    tree->SetDirectory(&output);
    tree->SetAutoFlush(64);
    tree->SetMaxVirtualSize(4 * 1024 * 1024);
    return tree;
}

struct EventCounts {
    Long64_t input_entry = 0;
    int event_number = 0;
    int all = 0;
    int charged = 0;
    int unknown_pdg = 0;
    int fxt = 0;

    void Bind(TTree& tree) {
        tree.Branch("InputEntry", &input_entry);
        tree.Branch("EventNr", &event_number);
        tree.Branch("NAll", &all);
        tree.Branch("NCh", &charged);
        tree.Branch("NUnknownPDG", &unknown_pdg);
        tree.Branch("NFxt", &fxt);
    }
};

struct EtaPtHistogram {
    TH2D values;
    Long64_t undefined = 0;
    Long64_t outside = 0;
    Long64_t axial = 0;

    EtaPtHistogram(const char* name, const char* title, const Options& options)
        : values(name, title, 240, -options.eta_max, options.eta_max, 150, 0.0, options.pt_max) {
        values.SetDirectory(nullptr);
    }

    void Fill(double longitudinal, double pt) {
        double eta;
        if (pt == 0.0) {
            if (longitudinal == 0.0) {
                ++undefined;
                return;
            }
            ++axial;
            // Put infinite eta in a flow bin without INF in ROOT statistics.
            eta = std::copysign(values.GetXaxis()->GetXmax() + 1.0, longitudinal);
        } else {
            eta = std::asinh(longitudinal / pt);
        }

        const int x = values.GetXaxis()->FindFixBin(eta);
        const int y = values.GetYaxis()->FindFixBin(pt);
        if (x == 0 || x > values.GetNbinsX() || y == 0 || y > values.GetNbinsY()) {
            ++outside;
        }
        values.Fill(eta, pt);
    }

    void Save() {
        Write(values);
        const std::string name = values.GetName();
        WriteParameter((name + "_undefined_eta").c_str(), undefined);
        WriteParameter((name + "_outside_plot").c_str(), outside);
        WriteParameter((name + "_axial").c_str(), axial);
    }
};

struct QaHistograms {
    EtaPtHistogram before;
    EtaPtHistogram after;
    std::unique_ptr<TH1D> all;
    std::unique_ptr<TH1D> charged;
    std::unique_ptr<TH1D> fxt;

    explicit QaHistograms(const Options& options)
        : before("hEtaPtBefore", "CM;#eta (axis opposite Projectile);p_{T} [GeV/c]", options),
          after("hEtaPtAfter", "Target rest frame;#eta (axis opposite Projectile);p_{T} [GeV/c]",
                options),
          all(MakeMultiplicityHistogram("hNAll", "All final particles;N_{all};Events")),
          charged(MakeMultiplicityHistogram("hNCh", "Charged final hadrons;N_{ch};Events")),
          fxt(MakeMultiplicityHistogram(
              "hNFxt",
              "Charged particles, -2 < #eta_{lab} < 0, p_{T} > 0.15 GeV/c;N_{FXT};Events")) {}

    void FillCounts(const EventCounts& counts) {
        all->Fill(counts.all);
        charged->Fill(counts.charged);
        fxt->Fill(counts.fxt);
    }

    void Save() {
        before.Save();
        after.Save();
        Write(*all);
        Write(*charged);
        Write(*fxt);
    }
};

struct Momentum {
    double energy;
    double px;
    double py;
    double pz;
    double pt;
    double lab_energy;
    double lab_pz;
};

Momentum TransformMomentum(const McParticle& particle, const Frame& frame) {
    const double energy = particle.e();
    const double px = particle.px();
    const double py = particle.py();
    const double pz = particle.pz();
    for (double component : {energy, px, py, pz}) {
        CheckFloat(component);
    }
    Require(energy > 0.0, "nonpositive particle energy");

    const double lab_energy = frame.Energy(energy, pz);
    const double lab_pz = frame.Longitudinal(energy, pz);
    CheckFloat(lab_energy);
    CheckFloat(lab_pz);
    Require(lab_energy > 0.0, "nonpositive transformed energy");
    return {energy, px, py, pz, std::hypot(px, py), lab_energy, lab_pz};
}

bool HasMaskedPosition(const McParticle& particle) {
    return particle.x() == -999.0 || particle.y() == -999.0 || particle.z() == -999.0 ||
           particle.t() == -999.0;
}

void TransformPosition(McParticle& particle, const Frame& frame) {
    const double x = particle.x();
    const double y = particle.y();
    const double time = frame.Energy(particle.t(), particle.z());
    const double z = frame.Longitudinal(particle.t(), particle.z());
    for (double component : {x, y, z, time}) {
        CheckFloat(component);
    }
    particle.setPosition(-x, y, z, time);
}

void RotateReactionPlane(McEvent& event) {
    if (event.phi() == -999.0) {
        return;
    }
    Require(std::isfinite(event.phi()), "non-finite event phi");
    const double pi = std::acos(-1.0);
    double phi = std::fmod(pi - event.phi(), 2.0 * pi);
    if (phi < 0.0) {
        phi += 2.0 * pi;
    }
    event.setPhi(phi);
}

bool PassesFxtCuts(const Momentum& momentum) {
    if (momentum.pt <= kFxtPtMin) {
        return false;
    }
    const double eta = std::asinh(momentum.lab_pz / momentum.pt);
    return eta > kFxtEtaMin && eta < kFxtEtaMax;
}

class FileAnalysis {
  public:
    EventCounts event_counts;
    Result totals;
    Long64_t masked_positions = 0;
    FxtCounts fxt_counts;

    FileAnalysis(const Frame& frame, const Options& options)
        : frame_(frame), hist_only_(options.hist_only), histograms_(options) {}

    void ProcessEvent(McEvent& event, TClonesArray& particles, Long64_t entry) {
        event_counts = {entry, event.eventNr(), particles.GetEntriesFast(), 0, 0, 0};
        for (int index = 0; index < event_counts.all; ++index) {
            auto* particle = static_cast<McParticle*>(particles.At(index));
            Require(particle != nullptr, "null McParticle at entry " + std::to_string(entry));
            ProcessParticle(*particle);
        }
        if (!hist_only_) {
            RotateReactionPlane(event);
        }

        ++totals.events;
        totals.particles += event_counts.all;
        totals.unknown_pdg += event_counts.unknown_pdg;
        ++fxt_counts[event_counts.fxt];
        histograms_.FillCounts(event_counts);
    }

    void SaveHistograms() {
        histograms_.Save();
    }

  private:
    void CountParticle(int pdg, const Momentum& momentum) {
        auto found = charges_.find(pdg);
        if (found == charges_.end()) {
            found = charges_.emplace(pdg, ParticleCharge(pdg)).first;
        }
        const Charge charge = found->second;
        event_counts.charged += charge.charged && charge.hadron;
        event_counts.unknown_pdg += !charge.known;
        event_counts.fxt += charge.charged && PassesFxtCuts(momentum);
    }

    void ProcessParticle(McParticle& particle) {
        const auto momentum = TransformMomentum(particle, frame_);
        histograms_.before.Fill(-momentum.pz, momentum.pt);
        histograms_.after.Fill(momentum.lab_pz, momentum.pt);
        CountParticle(particle.pdg(), momentum);

        const bool masked = HasMaskedPosition(particle);
        masked_positions += masked;
        if (hist_only_) {
            return;
        }
        particle.setMomentum(-momentum.px, momentum.py, momentum.lab_pz, momentum.lab_energy);
        if (!masked) {
            TransformPosition(particle, frame_);
        }
    }

    const Frame& frame_;
    bool hist_only_;
    QaHistograms histograms_;
    ChargeTable charges_;
};

void WriteFxtMetadata() {
    TNamed definition(
        "FxtDefinition",
        "NFxt=all known charged particles including leptons and nuclear fragments; "
        "Target rest frame, +z opposite Projectile; strict -2<eta<0 and pT>0.15 GeV/c; "
        "no spectator cut; unknown PDG excluded; plotted eta/pT limits do not affect NFxt");
    Write(definition);
    WriteParameter("FxtEtaMin", kFxtEtaMin);
    WriteParameter("FxtEtaMax", kFxtEtaMax);
    WriteParameter("FxtPtMin", kFxtPtMin);
}

void WriteFrameMetadata(const std::string& input, const Options& options, const Frame& frame) {
    TNamed description(
        "LabFrame",
        "v2: symmetric SMASH CM -> Target rest; x'=-x,y'=y,z'=-z after boost; "
        "all saved particles, no McDst cuts; NFxt uses lab acceptance, see FxtDefinition; "
        "NCh=charged hadrons including nuclear fragments; "
        "NUnknownPDG>0 means NCh is incomplete; masked positions unchanged; "
        "Event.stepT remains generator CM time; extra branches retain source meanings; "
        "source file-level metadata is in SourceMetadata");
    TNamed source("LabSource", std::filesystem::absolute(input).c_str());
    Write(description);
    WriteParameter("LabSqrtsNN", options.sqrtsnn);
    WriteParameter("LabNucleonMass", options.nucleon_mass);
    WriteParameter("LabBeta", frame.beta);
    WriteParameter("LabGamma", frame.gamma);
    Write(source);
}

void WriteCounters(const FileAnalysis& analysis, const Options& options) {
    WriteParameter("LabEvents", static_cast<Long64_t>(analysis.totals.events));
    WriteParameter("LabParticles", static_cast<Long64_t>(analysis.totals.particles));
    WriteParameter("LabUnknownPDG", static_cast<Long64_t>(analysis.totals.unknown_pdg));
    WriteParameter("LabMaskedPositions", analysis.masked_positions);
    WriteParameter("LabHistOnly", static_cast<int>(options.hist_only));
}

bool IsSourceMetadata(const std::string& class_name) {
    return class_name == "McRun" || class_name == "TNamed" ||
           class_name.rfind("TParameter<", 0) == 0;
}

void WriteTransformedRun(TObject& object, const char* name, TFile& output, const Options& options,
                         const Frame& frame) {
    auto* run = dynamic_cast<McRun*>(&object);
    if (options.hist_only || !run) {
        return;
    }
    // McRun stores magnitudes per nucleon; the axis sign is described by LabFrame.
    const double beam_momentum =
        2.0 * frame.gamma *
        std::sqrt(std::pow(options.sqrtsnn / 2.0, 2) - std::pow(options.nucleon_mass, 2));
    CheckFloat(beam_momentum);
    run->setPProj(beam_momentum);
    run->setPTarg(0.0);
    output.cd();
    Write(*run, name);
}

void CopySourceMetadata(TFile& input, TFile& output, const Options& options, const Frame& frame) {
    auto* metadata = output.mkdir("SourceMetadata");
    Require(metadata != nullptr, "cannot create source metadata directory");
    TIter keys(input.GetListOfKeys());
    std::unordered_set<std::string> copied;
    while (auto* key = dynamic_cast<TKey*>(keys())) {
        if (!IsSourceMetadata(key->GetClassName()) || !copied.insert(key->GetName()).second) {
            continue;
        }
        std::unique_ptr<TObject> object(key->ReadObj());
        Require(object != nullptr, "cannot read source metadata");
        metadata->cd();
        Write(*object, key->GetName());
        WriteTransformedRun(*object, key->GetName(), output, options, frame);
    }
    output.cd();
}

void CheckPlotLimits(const Options& options) {
    Require(std::isfinite(options.eta_max) && options.eta_max > 0.0 && options.eta_max <= 100.0 &&
                std::isfinite(options.pt_max) && options.pt_max > 0.0,
            "plot limits must be positive and finite, eta-max <= 100");
}

} // namespace

void InitializeRoot() {
    static std::once_flag initialized;
    std::call_once(initialized, [] {
        ROOT::EnableThreadSafety();
        TTree::SetMaxTreeSize(std::numeric_limits<Long64_t>::max());
        McEvent::Class();
        McParticle::Class();
        McRun::Class();
        KnownCharges();
    });
}

StagedResult WriteTransformedFile(const std::string& input_path, const std::string& output_path,
                                  const Options& options) {
    InitializeRoot();
    const Frame frame(options);
    CheckPlotLimits(options);
    auto input = OpenInput(input_path);
    auto& source = ReadMcDstTree(*input);
    InputBuffers buffers(source, options.hist_only);
    auto output = OpenOutput(output_path);
    auto destination = CloneOutputTree(source, *output, options.hist_only);

    FileAnalysis analysis(frame, options);
    TTree multiplicity("LabMultiplicity",
                       "Full multiplicities and FXT charged-particle acceptance");
    multiplicity.SetDirectory(output.get());
    multiplicity.SetAutoFlush(1024);
    analysis.event_counts.Bind(multiplicity);

    for (Long64_t entry = 0; entry < source.GetEntries(); ++entry) {
        auto& event = buffers.ReadEvent(source, entry);
        analysis.ProcessEvent(event, buffers.particles, entry);
        if (destination) {
            Require(destination->Fill() >= 0, "cannot write transformed event");
        }
        Require(multiplicity.Fill() >= 0 && !output->TestBit(TFile::kWriteError),
                "cannot write event counts");
    }

    output->cd();
    if (destination) {
        Write(*destination);
    }
    Write(multiplicity);
    analysis.SaveHistograms();
    WriteFxtMetadata();
    WriteFrameMetadata(input_path, options, frame);
    WriteCounters(analysis, options);
    CopySourceMetadata(*input, *output, options, frame);

    // Detach stack-owned trees before closing the file to avoid double ownership.
    multiplicity.SetDirectory(nullptr);
    destination.reset();
    CloseOutput(*output, "failed to finalize ROOT output");
    return {analysis.totals, std::move(analysis.fxt_counts)};
}

} // namespace mcdst_lab::detail
