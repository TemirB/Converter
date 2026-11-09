#include "ConverterValidation.h"
#include "OscarConverter.h"
#include "OscarConverterV1.h"

#include "McEvent.h"
#include "McParticle.h"
#include "TClonesArray.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TTree.h"

#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <locale>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using EventKey = std::pair<int, int>; // Event number and occurrence, starting at zero.
using ParticlePair = std::pair<const McParticle*, const McParticle*>;
using Getter = Double_t (McParticle::*)() const;
constexpr std::array<Getter, 6> kGetters = {&McParticle::px, &McParticle::py, &McParticle::pz,
                                            &McParticle::x,  &McParticle::y,  &McParticle::z};
constexpr std::array<const char*, 8> kComponentNames = {
    "px", "py", "pz", "x", "y", "z", "multiplicity", "impact_parameter"};

class ValidationHistograms {
  public:
    explicit ValidationHistograms(const std::string& filename)
        : file_(TFile::Open(filename.c_str(), "RECREATE", "Current converter vs v1")) {
        eventStatus_.SetDirectory(nullptr);
        particleStatus_.SetDirectory(nullptr);
        if (!file_ || file_->IsZombie())
            throw std::runtime_error("cannot create histogram file: " + filename);
        file_->SetCompressionLevel(1);
        for (std::size_t i = 0; i < histograms_.size(); ++i) {
            const std::string name = std::string("diff_") + kComponentNames[i];
            const char* unit = i < 3 ? " [GeV/c]" : i == 6 ? "" : " [fm]";
            const std::string title = std::string(kComponentNames[i]) + " (current - v1);Delta " +
                                      kComponentNames[i] + unit +
                                      (i < 6 ? ";Matched particles" : ";Matched events");
            const double halfRange = i == 6 ? 100.5 : 1.005e-6;
            histograms_[i] =
                std::make_unique<TH1D>(name.c_str(), title.c_str(), 201, -halfRange, halfRange);
            histograms_[i]->SetDirectory(nullptr);
            histograms_[i]->SetCanExtend(TH1::kXaxis);
        }
        for (auto* histogram : {&eventStatus_, &particleStatus_}) {
            int bin = 1;
            for (const char* label : {"equal", "different", "current_only", "v1_only"})
                histogram->GetXaxis()->SetBinLabel(bin++, label);
        }
    }

    void FillParticle(const ParticlePair& pair, bool equal) {
        particleStatus_.Fill(StatusBin(pair.first, pair.second, equal));
        if (pair.first && pair.second)
            for (std::size_t i = 0; i < kGetters.size(); ++i)
                histograms_[i]->Fill((pair.first->*kGetters[i])() - (pair.second->*kGetters[i])());
    }

    void FillEvent(const McEvent* current, const McEvent* v1, int currentMultiplicity,
                   int v1Multiplicity, bool equal) {
        eventStatus_.Fill(StatusBin(current, v1, equal));
        if (current && v1) {
            histograms_[6]->Fill(currentMultiplicity - v1Multiplicity);
            histograms_[7]->Fill(current->b() - v1->b());
        }
    }

    void Write(const std::string& currentFilename, const std::string& v1Filename) {
        file_->cd();
        for (const auto& histogram : histograms_)
            WriteObject(*histogram);
        WriteObject(eventStatus_);
        WriteObject(particleStatus_);
        TNamed convention("comparison",
                          "Differences: current - v1. Particles matched by OSCAR ID; "
                          "events by event number and occurrence. Particle histograms "
                          "contain matched IDs; multiplicity and impact parameter "
                          "are filled once per matched event. Missing entries are "
                          "counted in event_status and particle_status.");
        TNamed currentInput("current_input", currentFilename.c_str());
        TNamed v1Input("v1_input", v1Filename.c_str());
        WriteObject(convention);
        WriteObject(currentInput);
        WriteObject(v1Input);
        file_->Write();
        file_->Close();
        if (file_->TestBit(TFile::kWriteError))
            throw std::runtime_error("failed to finalize histogram file");
    }

  private:
    static int StatusBin(bool current, bool v1, bool equal) {
        return !current ? 4 : !v1 ? 3 : equal ? 1 : 2;
    }

    static void WriteObject(TObject& object) {
        if (object.Write() <= 0)
            throw std::runtime_error(std::string("failed to write ROOT object: ") +
                                     object.GetName());
    }

    std::unique_ptr<TFile> file_;
    std::array<std::unique_ptr<TH1D>, 8> histograms_;
    TH1D eventStatus_{"event_status", "Event comparison;Status;Events", 4, .5, 4.5};
    TH1D particleStatus_{"particle_status", "Particle comparison;Status;Particles", 4, .5, 4.5};
};

class McDstInput {
  public:
    explicit McDstInput(const std::string& filename)
        : file_(TFile::Open(filename.c_str(), "READ")) {
        if (!file_ || file_->IsZombie())
            throw std::runtime_error("cannot open McDst file: " + filename);
        tree_ = dynamic_cast<TTree*>(file_->Get("McDst"));
        if (!tree_ || tree_->SetBranchAddress("Event", &eventPointer_) < 0 ||
            tree_->SetBranchAddress("Particle", &particlePointer_) < 0 ||
            tree_->SetBranchAddress("SourceId", &idPointer_) < 0)
            throw std::runtime_error("missing McDst Event/Particle/SourceId branches: " + filename);
        tree_->GetBranch("Event")->SetAutoDelete(kFALSE);
        tree_->GetBranch("Particle")->SetAutoDelete(kFALSE);
        tree_->GetBranch("SourceId")->SetAutoDelete(kFALSE);
    }

    ~McDstInput() {
        if (tree_)
            tree_->ResetBranchAddresses();
    }

    std::map<EventKey, Long64_t> Index() {
        std::map<EventKey, Long64_t> result;
        std::map<int, int> occurrences;
        for (Long64_t entry = 0; entry < tree_->GetEntries(); ++entry) {
            if (tree_->GetBranch("Event")->GetEntry(entry) <= 0 || events_.GetEntriesFast() != 1)
                throw std::runtime_error("cannot read McDst event metadata");
            const int number = static_cast<McEvent*>(events_.At(0))->eventNr();
            result.emplace(EventKey{number, occurrences[number]++}, entry);
        }
        return result;
    }

    void Read(Long64_t entry) {
        if (tree_->GetEntry(entry) <= 0 || events_.GetEntriesFast() != 1 ||
            ids_.size() != static_cast<std::size_t>(particles_.GetEntriesFast()))
            throw std::runtime_error("invalid McDst event or SourceId count");
    }

    const McEvent* Event() const {
        return static_cast<McEvent*>(events_.At(0));
    }
    int Multiplicity() const {
        return particles_.GetEntriesFast();
    }

    void AddParticles(std::map<int, ParticlePair>& result, bool current) const {
        for (int i = 0; i < particles_.GetEntriesFast(); ++i) {
            auto& pair = result[ids_[i]];
            auto& slot = current ? pair.first : pair.second;
            if (slot)
                throw std::runtime_error("duplicate SourceId in McDst event");
            slot = static_cast<McParticle*>(particles_.At(i));
        }
    }

  private:
    std::unique_ptr<TFile> file_;
    TClonesArray events_{"McEvent", 1}, particles_{"McParticle", 1000};
    TClonesArray* eventPointer_ = &events_;
    TClonesArray* particlePointer_ = &particles_;
    std::vector<int> ids_;
    std::vector<int>* idPointer_ = &ids_;
    TTree* tree_ = nullptr;
};

bool EqualParticles(const ParticlePair& pair) {
    if (!pair.first || !pair.second || pair.first->pdg() != pair.second->pdg())
        return false;
    for (auto getter : kGetters)
        if ((pair.first->*getter)() != (pair.second->*getter)())
            return false;
    return true;
}

void WriteValue(std::ostream& output, std::optional<double> value) {
    output << ',';
    if (value)
        output << *value;
}

void WritePair(std::ostream& output, std::optional<double> current, std::optional<double> v1) {
    WriteValue(output, current);
    WriteValue(output, v1);
    WriteValue(output, current && v1 ? std::optional<double>(*current - *v1) : std::nullopt);
}

void WriteRow(std::ostream& output, const EventKey& key, std::optional<int> id, const char* status,
              const ParticlePair& pair, const McDstInput* current, const McDstInput* v1) {
    output << (id ? "particle" : "event") << ',' << key.first << ',' << key.second << ',';
    if (id)
        output << *id;
    output << ',' << status;
    WriteValue(output, pair.first ? std::optional<double>(pair.first->index()) : std::nullopt);
    WriteValue(output, pair.second ? std::optional<double>(pair.second->index()) : std::nullopt);
    WriteValue(output, pair.first ? std::optional<double>(pair.first->pdg()) : std::nullopt);
    WriteValue(output, pair.second ? std::optional<double>(pair.second->pdg()) : std::nullopt);
    for (auto getter : kGetters)
        WritePair(output,
                  pair.first ? std::optional<double>((pair.first->*getter)()) : std::nullopt,
                  pair.second ? std::optional<double>((pair.second->*getter)()) : std::nullopt);
    WritePair(output, current ? std::optional<double>(current->Multiplicity()) : std::nullopt,
              v1 ? std::optional<double>(v1->Multiplicity()) : std::nullopt);
    WritePair(output, current ? std::optional<double>(current->Event()->b()) : std::nullopt,
              v1 ? std::optional<double>(v1->Event()->b()) : std::nullopt);
    output << '\n';
}

void CheckDistinct(const std::filesystem::path& source, const std::filesystem::path& target) {
    if (std::filesystem::weakly_canonical(source) == std::filesystem::weakly_canonical(target) ||
        (std::filesystem::exists(source) && std::filesystem::exists(target) &&
         std::filesystem::equivalent(source, target)))
        throw std::runtime_error("input and output must be different files: " + target.string());
}

} // namespace

ValidationResult CompareMcDstVersions(const std::string& currentFilename,
                                      const std::string& v1Filename,
                                      const std::string& reportFilename) {
    TDirectory::TContext context;
    try {
        CheckDistinct(currentFilename, reportFilename);
        CheckDistinct(v1Filename, reportFilename);
        auto histogramFilename = std::filesystem::path(reportFilename);
        histogramFilename.replace_extension(".root");
        CheckDistinct(currentFilename, histogramFilename);
        CheckDistinct(v1Filename, histogramFilename);
        CheckDistinct(reportFilename, histogramFilename);
        McDstInput current(currentFilename), v1(v1Filename);
        const auto currentIndex = current.Index(), v1Index = v1.Index();
        std::map<EventKey, std::pair<Long64_t, Long64_t>> index;
        for (const auto& entry : currentIndex)
            index.emplace(entry.first, std::make_pair(entry.second, -1));
        for (const auto& entry : v1Index) {
            const auto found = index.find(entry.first);
            if (found == index.end())
                index.emplace(entry.first, std::make_pair(-1, entry.second));
            else
                found->second.second = entry.second;
        }

        ValidationHistograms histograms(histogramFilename.string());
        std::ofstream report(reportFilename);
        if (!report)
            throw std::runtime_error("cannot create diff report: " + reportFilename);
        report.imbue(std::locale::classic());
        report << std::setprecision(std::numeric_limits<double>::max_digits10);
        report << "row_type,event_nr,event_occurrence,particle_id,status,"
                  "current_index,v1_index,current_pdg,v1_pdg";
        for (const char* name : kComponentNames)
            report << ",current_" << name << ",v1_" << name << ",diff_" << name;
        report << '\n';

        std::size_t differingEvents = 0, matchedParticles = 0, differingParticles = 0;
        std::size_t currentOnlyParticles = 0, v1OnlyParticles = 0;
        for (const auto& entry : index) {
            const auto* currentEvent = entry.second.first >= 0 ? &current : nullptr;
            const auto* v1Event = entry.second.second >= 0 ? &v1 : nullptr;
            if (currentEvent)
                current.Read(entry.second.first);
            if (v1Event)
                v1.Read(entry.second.second);
            std::map<int, ParticlePair> particles;
            if (currentEvent)
                current.AddParticles(particles, true);
            if (v1Event)
                v1.AddParticles(particles, false);
            bool equal = currentEvent && v1Event && current.Multiplicity() == v1.Multiplicity() &&
                         current.Event()->b() == v1.Event()->b();
            for (const auto& particle : particles) {
                const auto& pair = particle.second;
                const bool same = EqualParticles(pair);
                histograms.FillParticle(pair, same);
                equal = equal && same;
                if (pair.first && pair.second) {
                    ++matchedParticles;
                    if (!same)
                        ++differingParticles;
                } else if (pair.first)
                    ++currentOnlyParticles;
                else
                    ++v1OnlyParticles;
            }
            if (!equal)
                ++differingEvents;
            histograms.FillEvent(
                currentEvent ? current.Event() : nullptr, v1Event ? v1.Event() : nullptr,
                currentEvent ? current.Multiplicity() : 0, v1Event ? v1.Multiplicity() : 0, equal);
            const char* status = !currentEvent ? "v1_only"
                                 : !v1Event    ? "current_only"
                                 : equal       ? "equal"
                                               : "different";
            WriteRow(report, entry.first, std::nullopt, status, {}, currentEvent, v1Event);
            for (const auto& particle : particles) {
                const auto& pair = particle.second;
                const char* particleStatus = !pair.first            ? "v1_only"
                                             : !pair.second         ? "current_only"
                                             : EqualParticles(pair) ? "equal"
                                                                    : "different";
                WriteRow(report, entry.first, particle.first, particleStatus, pair, currentEvent,
                         v1Event);
            }
        }
        report.close();
        if (!report)
            throw std::runtime_error("failed to write diff report: " + reportFilename);
        histograms.Write(currentFilename, v1Filename);
        std::cout << "Events: current=" << currentIndex.size() << " v1=" << v1Index.size()
                  << " differing=" << differingEvents << '\n'
                  << "Particles: matched=" << matchedParticles
                  << " differing=" << differingParticles << " current_only=" << currentOnlyParticles
                  << " v1_only=" << v1OnlyParticles << '\n'
                  << "Diff report (current - v1): " << reportFilename << '\n'
                  << "ROOT histograms (current - v1): " << histogramFilename.string() << '\n';
        return differingEvents == 0 ? ValidationResult::Equal : ValidationResult::Different;
    } catch (const std::exception& error) {
        std::cerr << "Comparison failed: " << error.what() << '\n';
        return ValidationResult::Error;
    }
}

ValidationResult ValidateOscarVersions(const std::string& inputFilename,
                                       const std::string& outputDirectory) {
    TDirectory::TContext context;
    namespace fs = std::filesystem;
    try {
        if (!fs::is_regular_file(inputFilename))
            throw std::runtime_error("input must be an OSCAR file: " + inputFilename);
        const fs::path directory =
            outputDirectory.empty() ? inputFilename + ".validation" : outputDirectory;
        const auto current = directory / "current.McDst.root";
        const auto v1 = directory / "v1.McDst.root";
        const auto report = directory / "diff.csv";
        const auto histograms = directory / "diff.root";
        for (const auto& output : {current, v1, report, histograms})
            CheckDistinct(inputFilename, output);
        fs::create_directories(directory);
        ConversionOptions options;
        options.record_source_ids = true;
        std::cout << "Converting with current version...\n";
        if (!OscarConverter().Convert(inputFilename, current.string(), options))
            return ValidationResult::Error;
        std::cout << "Converting with v1...\n";
        if (!OscarConverterV1().Convert(inputFilename, v1.string(), true))
            return ValidationResult::Error;
        return CompareMcDstVersions(current.string(), v1.string(), report.string());
    } catch (const std::exception& error) {
        std::cerr << "Validation failed: " << error.what() << '\n';
        return ValidationResult::Error;
    }
}
