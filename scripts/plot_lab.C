#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TStyle.h>

#include <algorithm>
#include <stdexcept>
#include <string>

// root -l -b -q 'scripts/plot_lab.C("lab.root","eta_pt")'
void plot_lab(const char* input, const char* prefix = "eta_pt") {
    TFile file(input, "READ");
    auto* before = dynamic_cast<TH2*>(file.Get("hEtaPtBefore"));
    auto* after = dynamic_cast<TH2*>(file.Get("hEtaPtAfter"));
    if (file.IsZombie() || !before || !after)
        throw std::runtime_error("missing eta-pT histograms");
    gStyle->SetOptStat(0);
    TCanvas canvas("eta_pt_comparison", "CM and Target rest frame", 1400, 600);
    canvas.Divide(2, 1);
    const double maximum = std::max({1., before->GetMaximum(), after->GetMaximum()});
    int pad = 0;
    for (auto* histogram : {before, after}) {
        canvas.cd(++pad);
        gPad->SetRightMargin(.20);
        gPad->SetBottomMargin(.13);
        gPad->SetLogz();
        histogram->SetMinimum(.5);
        histogram->SetMaximum(maximum);
        histogram->GetZaxis()->SetTitle("Particles / bin");
        histogram->GetZaxis()->SetLabelSize(.032);
        histogram->GetZaxis()->SetTitleSize(.032);
        histogram->GetZaxis()->SetTitleOffset(1.5);
        histogram->Draw("COLZ");
    }
    canvas.SaveAs((std::string(prefix) + ".png").c_str());
    canvas.SaveAs((std::string(prefix) + ".pdf").c_str());
}
