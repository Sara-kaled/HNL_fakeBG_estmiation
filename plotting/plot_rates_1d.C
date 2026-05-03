// =============================================================================
// FILE:    plotting/plot_rates_1d.C
// PURPOSE: Self-contained 1D fake-rate and real-efficiency plotting.
//          Replaces the old wrapper that called legacy plotSubtractedFake...
//          and plot_real.C — those are now NOT required.
//
//          Produces clean ATLAS Internal-style plots:
//            - 1D fake rate ε_f(pT) and ε_f(|η|), per CR region
//            - 1D real efficiency ε_r(pT) and ε_r(|η|), per CR region
//          With proper √s, luminosity, and ATLAS-style colors and markers.
//          Clopper-Pearson errors via TGraphAsymmErrors.
//
// INPUTS:  ../efficiency/outputs/MC_histograms.root      (prompt MC)
//          ../efficiency/outputs/Data_histograms.root    (data)
//          ../efficiency/outputs/real_eff_histograms.root (real e_r MC)
// OUTPUTS: outputs/02_1d_rates/fake/<region>_<obs>.pdf
//          outputs/02_1d_rates/real/<region>_<obs>.pdf
// =============================================================================
#include <TFile.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include <iostream>
#include <vector>
#include <string>

#include "../utils/AtlasStyle.h"

// ---------------------------------------------------------------------------
// Compute ε from (pass, total) with Clopper-Pearson errors.
// Both inputs may be data or (data − MC_prompt) bin-by-bin.
// ---------------------------------------------------------------------------
static TGraphAsymmErrors* makeEffGraph(const TH1D* hPass, const TH1D* hTotal)
{
    auto* g = new TGraphAsymmErrors();
    g->Divide(hPass, hTotal, "cl=0.683 b(1,1) mode");
    return g;
}

// ---------------------------------------------------------------------------
// Subtract MC histogram bin-by-bin, clamping negative bins to 0 to avoid
// pathological efficiencies. Returns a new histogram caller owns.
// ---------------------------------------------------------------------------
static TH1D* subtract(const TH1D* data, const TH1D* mc, const char* name)
{
    auto* h = (TH1D*)data->Clone(name);
    h->SetDirectory(nullptr);
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double v = data->GetBinContent(i) - mc->GetBinContent(i);
        if (v < 0) v = 0.;
        double e = std::sqrt(data->GetBinError(i)*data->GetBinError(i)
                           + mc  ->GetBinError(i)*mc  ->GetBinError(i));
        h->SetBinContent(i, v);
        h->SetBinError  (i, e);
    }
    return h;
}

// ---------------------------------------------------------------------------
// Draw one fake-rate graph + one real-efficiency graph in ATLAS Internal style
// for a given (region, observable). yLabel = "Fake efficiency" or "Real efficiency".
// ---------------------------------------------------------------------------
static void drawOne(TGraphAsymmErrors* g,
                    const char* xLabel, const char* yLabel,
                    const char* regionLabel, const char* outPath,
                    double yMax = 0.6)
{
    SetATLASStyle();
    TCanvas c("c", "", 700, 600);

    g->SetMarkerStyle(20);
    g->SetMarkerSize (1.2);
    g->SetMarkerColor(kATLAS_data);
    g->SetLineColor  (kATLAS_data);
    g->SetLineWidth  (2);
    g->Draw("AP");
    g->GetXaxis()->SetTitle(xLabel);
    g->GetYaxis()->SetTitle(yLabel);
    g->GetYaxis()->SetRangeUser(0., yMax);

    DrawAtlasLabel(0.16, 0.88);

    TLatex lat;
    lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.038);
    lat.DrawLatex(0.16, 0.78, regionLabel);

    c.SaveAs(outPath);
}

// ---------------------------------------------------------------------------
void plot_rates_1d()
{
    const std::string fMC_path   = "../efficiency/outputs/MC_histograms.root";
    const std::string fData_path = "../efficiency/outputs/Data_histograms.root";
    const std::string fReal_path = "../efficiency/outputs/real_eff_histograms.root";

    TFile* fMC   = TFile::Open(fMC_path.c_str(),   "READ");
    TFile* fData = TFile::Open(fData_path.c_str(), "READ");
    TFile* fReal = TFile::Open(fReal_path.c_str(), "READ");
    if (!fMC || !fData || !fReal ||
        fMC->IsZombie() || fData->IsZombie() || fReal->IsZombie()) {
        std::cerr << "ERROR: cannot open one of the input ROOT files.\n";
        return;
    }

    gSystem->mkdir("outputs/02_1d_rates/fake", true);
    gSystem->mkdir("outputs/02_1d_rates/real", true);

    const std::vector<std::string> regions = {"CR", "CR_MET_low", "CR_MET_high"};
    // pt and eta drive the AMM map; d0sig and MET are diagnostic / motivate
    // the prompt-subtraction and MET systematics. d0 is the bare impact
    // parameter (vs significance) — useful for cross-check with d0sig.
    const std::vector<std::string> obs = {"pt", "eta", "d0sig", "d0", "met"};

    auto regionLatex = [](const std::string& r) {
        if (r == "CR")          return "Fake CR (SS #mue)";
        if (r == "CR_MET_low")  return "Fake CR, E_{T}^{miss} < 50 GeV";
        if (r == "CR_MET_high") return "Fake CR, E_{T}^{miss} > 50 GeV";
        return r.c_str();
    };
    auto axisLatex = [](const std::string& o) {
        if (o == "pt")    return "Electron p_{T} [GeV]";
        if (o == "eta")   return "Electron |#eta|";
        if (o == "d0sig") return "Electron |d_{0}/#sigma|";
        if (o == "d0")    return "Electron |d_{0}|";
        if (o == "met")   return "E_{T}^{miss} [GeV]";
        return o.c_str();
    };

    // ---------- 1D Fake rate ----------
    std::cout << "\n=== 1D fake rate ε_f = (data−MC_prompt)_T / (data−MC_prompt)_L ===\n";
    for (const auto& reg : regions) {
        for (const auto& o : obs) {
            const std::string nT = "h_tight_" + o + "_" + reg;
            const std::string nL = "h_loose_" + o + "_" + reg;

            auto* dT = (TH1D*)fData->Get(nT.c_str());
            auto* dL = (TH1D*)fData->Get(nL.c_str());
            auto* mT = (TH1D*)fMC  ->Get(nT.c_str());
            auto* mL = (TH1D*)fMC  ->Get(nL.c_str());
            if (!dT || !dL || !mT || !mL) {
                std::cerr << "  missing hist for " << reg << " / " << o << "\n"; continue;
            }
            auto* sT = subtract(dT, mT, "sub_T");
            auto* sL = subtract(dL, mL, "sub_L");
            auto* g  = makeEffGraph(sT, sL);

            const std::string outPath = "outputs/02_1d_rates/fake/" + reg + "_" + o + ".pdf";
            drawOne(g, axisLatex(o), "Fake efficiency #varepsilon_{f}",
                    regionLatex(reg), outPath.c_str(), 0.6);
            std::cout << "  wrote " << outPath << "\n";
        }
    }

    // ---------- 1D Real efficiency (MC truth-matched) ----------
    std::cout << "\n=== 1D real efficiency ε_r = MC_T / MC_L (truth-prompt) ===\n";
    for (const auto& reg : regions) {
        for (const auto& o : obs) {
            const std::string nT = "h_tight_" + o + "_" + reg;
            const std::string nL = "h_loose_" + o + "_" + reg;

            auto* T = (TH1D*)fReal->Get(nT.c_str());
            auto* L = (TH1D*)fReal->Get(nL.c_str());
            if (!T || !L) {
                std::cerr << "  missing real hist for " << reg << " / " << o << "\n"; continue;
            }
            auto* g = makeEffGraph(T, L);
            const std::string outPath = "outputs/02_1d_rates/real/" + reg + "_" + o + ".pdf";
            drawOne(g, axisLatex(o), "Real efficiency #varepsilon_{r}",
                    regionLatex(reg), outPath.c_str(), 1.0);
            std::cout << "  wrote " << outPath << "\n";
        }
    }

    fMC->Close(); fData->Close(); fReal->Close();
    std::cout << "\nplot_rates_1d done.\n";
}
