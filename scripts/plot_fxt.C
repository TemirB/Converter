#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TParameter.h>
#include <TStyle.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

namespace fxt_plot {

constexpr std::array<const char*, 4> kLabels = {"0-10%", "10-30%", "30-50%", "50-80%"};
constexpr std::array<int, 4> kColors = {kRed - 7, kOrange - 2, kGreen + 1, kAzure - 4};
using Slices = std::array<std::unique_ptr<TH1D>, kLabels.size()>;

struct Calibration {
    TH1D& reference;
    TTree& mapping;
    Long64_t events;
};

struct PlotData {
    Slices slices;
    std::array<Long64_t, kLabels.size()> totals{};
    int maximum_multiplicity = 0;
};

Calibration ReadCalibration(TFile& file) {
    auto* reference = dynamic_cast<TH1D*>(file.Get("hNFxtReference"));
    auto* mapping = dynamic_cast<TTree*>(file.Get("FxtCentralityMap"));
    auto* events = dynamic_cast<TParameter<Long64_t>*>(file.Get("FxtCentralityEvents"));
    if (file.IsZombie() || !reference || !mapping || !events) {
        throw std::runtime_error("missing FXT calibration; rerun mcdst_lab with FXT support");
    }
    if (events->GetVal() == 0) {
        throw std::runtime_error("cannot plot FXT centrality for an empty sample");
    }
    return {*reference, *mapping, events->GetVal()};
}

Slices CreateSlices(const TH1D& reference) {
    Slices slices;
    for (std::size_t index = 0; index < slices.size(); ++index) {
        const std::string name = "fxt_slice_" + std::to_string(index);
        slices[index].reset(static_cast<TH1D*>(reference.Clone(name.c_str())));
        slices[index]->SetDirectory(nullptr);
        slices[index]->Reset();
        slices[index]->SetFillColor(kColors[index]);
        slices[index]->SetLineColor(kColors[index]);
    }
    return slices;
}

void AddFrequency(TH1D& histogram, int multiplicity, Long64_t events) {
    const int bin = histogram.FindFixBin(multiplicity);
    const double frequency = histogram.GetBinContent(bin) + events;
    histogram.SetBinContent(bin, frequency);
    histogram.SetBinError(bin, std::sqrt(frequency));
}

PlotData ReadPlotData(const Calibration& calibration) {
    PlotData data;
    data.slices = CreateSlices(calibration.reference);
    int multiplicity = 0;
    int class_index = -1;
    Long64_t events = 0;
    auto& mapping = calibration.mapping;
    mapping.SetBranchAddress("NFxt", &multiplicity);
    mapping.SetBranchAddress("Events", &events);
    mapping.SetBranchAddress("Class", &class_index);

    for (Long64_t entry = 0; entry < mapping.GetEntries(); ++entry) {
        mapping.GetEntry(entry);
        data.maximum_multiplicity = std::max(data.maximum_multiplicity, multiplicity);
        if (class_index < 0) {
            continue;
        }
        if (class_index >= static_cast<int>(data.slices.size())) {
            throw std::runtime_error("unexpected FXT centrality class");
        }
        AddFrequency(*data.slices[class_index], multiplicity, events);
        data.totals[class_index] += events;
    }
    mapping.ResetBranchAddresses();
    return data;
}

void ConfigureAxes(TH1D& reference, int maximum_multiplicity) {
    reference.SetTitle("FXT multiplicity: charged, -2 < #eta_{lab} < 0, p_{T} > 0.15 GeV/c");
    reference.GetXaxis()->SetTitle("N_{FXT}");
    reference.GetYaxis()->SetTitle("Events / multiplicity bin");
    reference.GetXaxis()->SetRangeUser(-0.5, std::min(65535.5, maximum_multiplicity * 1.08 + 1.0));
    reference.SetMinimum(0.5);
    reference.SetMaximum(std::max(1.0, reference.GetMaximum()) * 4.0);
    reference.SetLineColor(kBlack);
    reference.SetLineWidth(2);
}

void DrawClasses(const Calibration& calibration, const PlotData& data, TLegend& legend) {
    legend.SetBorderSize(0);
    const std::string header = "Calibration: " + std::to_string(calibration.events) + " events";
    legend.SetHeader(header.c_str());
    legend.AddEntry(&calibration.reference, "Full distribution (including 80-100%)", "l");
    for (std::size_t index = 0; index < data.slices.size(); ++index) {
        data.slices[index]->Draw("HIST SAME");
        const std::string label =
            std::string(kLabels[index]) + ": " + std::to_string(data.totals[index]) + " events";
        legend.AddEntry(data.slices[index].get(), label.c_str(), "f");
    }
    calibration.reference.Draw("HIST SAME");
    legend.Draw();
}

} // namespace fxt_plot

// Uses the whole calibration sample, even when invoked on just one output file.
// root -l -b -q 'scripts/plot_fxt.C("lab_files/sample.McDst.root","fxt_multi")'
void plot_fxt(const char* input, const char* prefix = "fxt_multi") {
    TFile file(input, "READ");
    const auto calibration = fxt_plot::ReadCalibration(file);
    const auto data = fxt_plot::ReadPlotData(calibration);

    gStyle->SetOptStat(0);
    TCanvas canvas("fxt_centrality", "FXT multiplicity and centrality", 1100, 750);
    canvas.SetLogy();
    canvas.SetLeftMargin(0.12);
    canvas.SetBottomMargin(0.13);
    fxt_plot::ConfigureAxes(calibration.reference, data.maximum_multiplicity);
    calibration.reference.Draw("HIST");
    TLegend legend(0.54, 0.59, 0.89, 0.88);
    fxt_plot::DrawClasses(calibration, data, legend);

    canvas.RedrawAxis();
    canvas.SaveAs((std::string(prefix) + ".png").c_str());
    canvas.SaveAs((std::string(prefix) + ".pdf").c_str());
}
