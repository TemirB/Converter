#include "McDstLabInternal.h"

#include <TBranch.h>
#include <TNamed.h>

namespace mcdst_lab::detail {
namespace {

int FindCentralityClass(double percentile) {
    for (std::size_t index = 0; index < kCentralityClasses.size(); ++index) {
        if (percentile < kCentralityClasses[index].upper_percentile) {
            return static_cast<int>(index);
        }
    }
    return -1;
}

struct CentralityRow {
    double low = 0.0;
    double high = 0.0;
    double percentile = 0.0;
    int class_index = -1;

    void Assign(const CentralityBin& bin) {
        low = bin.percentile_low;
        high = bin.percentile_high;
        percentile = bin.percentile;
        class_index = bin.class_index;
    }

    std::array<TBranch*, 4> AddEventBranches(TTree& tree) {
        return {
            tree.Branch("FxtCentralityLow", &low),
            tree.Branch("FxtCentralityHigh", &high),
            tree.Branch("FxtCentralityPercentile", &percentile),
            tree.Branch("FxtCentralityClass", &class_index),
        };
    }
};

class ClassHistograms {
  public:
    ClassHistograms() : classes_("hFxtCentrality", "FXT centrality;Class [%];Events", 5, 0., 5.) {
        classes_.SetDirectory(nullptr);
        for (std::size_t index = 0; index < kCentralityClasses.size(); ++index) {
            const auto& definition = kCentralityClasses[index];
            const std::string title =
                std::string("FXT centrality ") + definition.label + "% ;N_{FXT};Events";
            slices_[index] = MakeMultiplicityHistogram(definition.histogram_name, title.c_str());
            classes_.GetXaxis()->SetBinLabel(index + 1, definition.label);
        }
        classes_.GetXaxis()->SetBinLabel(5, "80-100");
    }

    void Fill(int multiplicity, int class_index) {
        if (class_index < 0) {
            classes_.Fill(4.5);
            return;
        }
        classes_.Fill(class_index + 0.5);
        slices_[class_index]->Fill(multiplicity);
    }

    void Save() {
        Write(classes_);
        for (auto& slice : slices_) {
            Write(*slice);
        }
    }

  private:
    TH1D classes_;
    std::array<std::unique_ptr<TH1D>, kCentralityClasses.size()> slices_;
};

void WriteEventCentrality(TTree& tree, const CentralityCalibration& calibration) {
    int multiplicity = 0;
    CentralityRow row;
    ResetTreeAddresses reset(tree);
    tree.SetBranchStatus("*", false);
    tree.SetBranchStatus("NFxt", true);
    Require(tree.SetBranchAddress("NFxt", &multiplicity) >= 0, "cannot read NFxt");

    const auto branches = row.AddEventBranches(tree);
    for (auto* branch : branches) {
        Require(branch != nullptr, "cannot create FXT centrality branch");
    }

    ClassHistograms histograms;
    auto* multiplicity_branch = tree.GetBranch("NFxt");
    for (Long64_t entry = 0; entry < tree.GetEntries(); ++entry) {
        Require(multiplicity_branch->GetEntry(entry) > 0, "cannot read staged FXT multiplicity");
        const auto found = calibration.bins.find(multiplicity);
        Require(found != calibration.bins.end(), "NFxt absent from centrality calibration");
        row.Assign(found->second);
        for (auto* branch : branches) {
            Require(branch->Fill() >= 0, "cannot write FXT centrality");
        }
        histograms.Fill(multiplicity, row.class_index);
    }

    tree.SetBranchStatus("*", true);
    Write(tree);
    histograms.Save();
}

void AddFrequency(TH1D& histogram, int multiplicity, Long64_t events) {
    const int bin = histogram.FindFixBin(multiplicity);
    const double frequency = histogram.GetBinContent(bin) + events;
    histogram.SetBinContent(bin, frequency);
    histogram.SetBinError(bin, std::sqrt(frequency));
}

void WriteCalibrationMap(const CentralityCalibration& calibration) {
    auto reference =
        MakeMultiplicityHistogram("hNFxtReference", "FXT calibration sample;N_{FXT};Events");
    TTree mapping("FxtCentralityMap", "Shared right-tail integral over all calibration events");
    int multiplicity = 0;
    Long64_t events = 0;
    CentralityRow row;
    mapping.Branch("NFxt", &multiplicity);
    mapping.Branch("Events", &events);
    mapping.Branch("PercentileLow", &row.low);
    mapping.Branch("PercentileHigh", &row.high);
    mapping.Branch("Percentile", &row.percentile);
    mapping.Branch("Class", &row.class_index);

    // Use exact frequencies, including multiplicities outside the plotted range.
    for (const auto& [value, bin] : calibration.bins) {
        multiplicity = value;
        events = bin.events;
        row.Assign(bin);
        Require(mapping.Fill() >= 0, "cannot write FXT centrality map");
        AddFrequency(*reference, multiplicity, events);
    }

    reference->ResetStats();
    reference->SetEntries(calibration.events);
    Write(*reference);
    Write(mapping);
    mapping.SetDirectory(nullptr);
}

void WriteCalibrationMetadata(const CentralityCalibration& calibration,
                              const CalibrationScope& scope) {
    TNamed definition(
        "FxtCentralityDefinition",
        "right-tail event integral including NFxt=0; percentile=midpoint of tied-bin interval; "
        "equal NFxt stays together; classes [0,10),[10,30),[30,50),[50,80); "
        "class=-1 for [80,100); calibration includes successfully staged selected files");
    TNamed description("FxtCentralityScope", scope.description);
    Write(definition);
    Write(description);
    WriteParameter("FxtCentralityEvents", calibration.events);
    WriteParameter("FxtCentralityFiles", static_cast<Long64_t>(scope.files));
    WriteParameter("FxtCentralityFailedFiles", static_cast<Long64_t>(scope.failed_files));
}

} // namespace

CentralityCalibration CalibrateCentrality(const FxtCounts& counts) {
    CentralityCalibration calibration;
    for (const auto& [multiplicity, events] : counts) {
        calibration.events += events;
    }

    Long64_t events_above = 0;
    for (auto entry = counts.rbegin(); entry != counts.rend(); ++entry) {
        const auto [multiplicity, events] = *entry;
        const double low = 100.0 * static_cast<double>(events_above) / calibration.events;
        const double high = 100.0 * static_cast<double>(events_above + events) / calibration.events;
        const double midpoint = (low + high) / 2.0;
        const CentralityBin bin{events, low, high, midpoint, FindCentralityClass(midpoint)};
        calibration.bins.emplace(multiplicity, bin);
        events_above += events;
    }
    return calibration;
}

void AddCentrality(const std::string& path, const CentralityCalibration& calibration,
                   const CalibrationScope& scope) {
    TFile file(path.c_str(), "UPDATE");
    Require(!file.IsZombie(), "cannot reopen staged output for FXT centrality");
    auto* multiplicity = dynamic_cast<TTree*>(file.Get("LabMultiplicity"));
    Require(multiplicity != nullptr, "missing LabMultiplicity for FXT centrality");

    // Only the small event-count tree is read again; particle baskets stay untouched.
    WriteEventCentrality(*multiplicity, calibration);
    WriteCalibrationMap(calibration);
    WriteCalibrationMetadata(calibration, scope);
    CloseOutput(file, "failed to finalize FXT centrality");
}

} // namespace mcdst_lab::detail
