// =============================================================================
// FILE:    plotting/plot_real.C
// PURPOSE: Real-electron efficiency plots from the MC tight/loose histograms
//          produced by efficiency/real_eff.C.
//
//          Reconstructed after the file was accidentally overwritten by the
//          contents of efficiency/real_eff.C — the wrapper macros
//          plot_eff_maps.C and plot_rates_1d.C #include this file and call
//          plot_real(...) with the signature below.
//
// INPUTS:  MC_histograms.root produced by real_eff.C, containing
//             h2_tight_eta_pt_<region>   (TH2D, X = pT [GeV], Y = |η|)
//             h2_loose_eta_pt_<region>
//             h_tight_pt_<region>, h_loose_pt_<region>     (1D)
//             h_tight_eta_<region>, h_loose_eta_<region>   (1D)
//
// HELPER NAMES: makeEffGraphReal / divide2DReal — distinct from the
//          identically-shaped helpers in plotSubtractedFakeElectronEfficiencies.c,
//          so that wrappers including both files do not get a duplicate-
//          definition error from CINT/Cling.
// =============================================================================
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include "TSystem.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

// --------------------------------------------------------------------
// 1D efficiency graph with Clopper–Pearson errors (file-local helper).
// Renamed from `makeEffGraph` because plotSubtractedFakeElectronEfficiencies.c
// already provides one — including both files in the same wrapper would
// otherwise hit a redefinition error under ACLiC.
// --------------------------------------------------------------------
static TGraphAsymmErrors* makeEffGraphReal(TH1D* hNum, TH1D* hDen,
                                           const char* name,
                                           Color_t color = kBlack)
{
    if (!hNum || !hDen) {
        cout << "Error: null histogram in makeEffGraphReal\n";
        return nullptr;
    }
    auto* g = new TGraphAsymmErrors(hNum, hDen, "n");   // CP intervals
    g->SetName(name);
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.0);
    g->SetLineWidth(2);
    return g;
}

// --------------------------------------------------------------------
// 2D bin-by-bin efficiency with symmetric Clopper–Pearson errors.
// --------------------------------------------------------------------
static TH2D* divide2DReal(TH2D* hNum, TH2D* hDen, const char* name)
{
    if (!hNum || !hDen) {
        cout << "Error: null 2D histogram in divide2DReal\n";
        return nullptr;
    }
    auto* hEff = (TH2D*)hDen->Clone(name);
    hEff->Reset();
    for (int ix = 1; ix <= hEff->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hEff->GetNbinsY(); ++iy) {
            double num = hNum->GetBinContent(ix, iy);
            double den = hDen->GetBinContent(ix, iy);
            if (den > 0 && num >= 0 && num <= den) {
                double eff   = num / den;
                double up    = TEfficiency::ClopperPearson(den, num, 0.683, true)  - eff;
                double down  = eff - TEfficiency::ClopperPearson(den, num, 0.683, false);
                double err   = 0.5 * (up + down);
                hEff->SetBinContent(ix, iy, eff);
                hEff->SetBinError  (ix, iy, err);
            } else {
                hEff->SetBinContent(ix, iy, 0.);
                hEff->SetBinError  (ix, iy, 0.);
            }
        }
    }
    return hEff;
}

// ---------------------------------------------------------------------------
// Main entry point — signature matches what plot_eff_maps.C / plot_rates_1d.C
// expect.  Default args allow a quick interactive `plot_real()` invocation.
// ---------------------------------------------------------------------------
void plot_real(const char* mc_filename = "MC_histograms.root",
               const char* region      = "CR",
               bool        plot2D      = true,
               bool        logx        = false,
               float       yMin        = 0.0f,
               float       yMax        = 1.0f,
               const char* outDir      = ".",
               const char* fmt         = "pdf")
{
    gStyle->SetOptStat(0);
    gSystem->mkdir(outDir, kTRUE);

    TFile* fIn = TFile::Open(mc_filename, "READ");
    if (!fIn || fIn->IsZombie()) {
        cerr << "ERROR: cannot open " << mc_filename << "\n";
        return;
    }

    auto get = [&](const char* prefix, const char* axis) {
        TString nm = TString::Format("%s_%s_%s", prefix, axis, region);
        return (TH1D*)fIn->Get(nm);
    };
    auto get2 = [&](const char* prefix) {
        TString nm = TString::Format("%s_eta_pt_%s", prefix, region);
        return (TH2D*)fIn->Get(nm);
    };

    TH1D* hMCTightPt   = get ("h_tight",  "pt");
    TH1D* hMCLoosePt   = get ("h_loose",  "pt");
    TH1D* hMCTightEta  = get ("h_tight",  "eta");
    TH1D* hMCLooseEta  = get ("h_loose",  "eta");
    TH2D* hMCTight2D   = get2("h2_tight");
    TH2D* hMCLoose2D   = get2("h2_loose");

    if (!hMCTightPt || !hMCLoosePt) {
        cerr << "ERROR: 1D pT histograms missing for region " << region << "\n";
        fIn->Close();
        return;
    }

    // ---- 1D pT efficiency -------------------------------------------------
    TGraphAsymmErrors* gEffPt  = makeEffGraphReal(hMCTightPt,  hMCLoosePt,  "gEffPt");
    TGraphAsymmErrors* gEffEta = makeEffGraphReal(hMCTightEta, hMCLooseEta, "gEffEta");

    {
        TCanvas c("c_pt", "real eff vs pT", 800, 700);
        if (logx) c.SetLogx();
        gEffPt->SetTitle((TString::Format(
            "Real-electron efficiency vs p_{T} (%s);Electron p_{T} [GeV];#varepsilon_{r}",
            region)));
        gEffPt->GetYaxis()->SetRangeUser(yMin, yMax);
        gEffPt->Draw("AP");
        TLatex t; t.SetNDC(); t.SetTextFont(72); t.SetTextSize(0.04);
        t.DrawLatex(0.16, 0.93, "ATLAS");
        t.SetTextFont(42); t.DrawLatex(0.26, 0.93, "Internal");
        c.SaveAs(TString::Format("%s/eff_real_pt_%s.%s", outDir, region, fmt));
    }

    // ---- 1D |η| efficiency ------------------------------------------------
    if (gEffEta) {
        TCanvas c("c_eta", "real eff vs |eta|", 800, 700);
        gEffEta->SetTitle((TString::Format(
            "Real-electron efficiency vs |#eta| (%s);|Electron #eta|;#varepsilon_{r}",
            region)));
        gEffEta->GetYaxis()->SetRangeUser(yMin, yMax);
        gEffEta->Draw("AP");
        TLatex t; t.SetNDC(); t.SetTextFont(72); t.SetTextSize(0.04);
        t.DrawLatex(0.16, 0.93, "ATLAS");
        t.SetTextFont(42); t.DrawLatex(0.26, 0.93, "Internal");
        c.SaveAs(TString::Format("%s/eff_real_eta_%s.%s", outDir, region, fmt));
    }

    // ---- 2D pT × |η| map --------------------------------------------------
    if (plot2D && hMCTight2D && hMCLoose2D) {
        TH2D* hEff2D = divide2DReal(hMCTight2D, hMCLoose2D, "hEff2D");
        TCanvas c("c2d", "real eff 2D", 820, 720);
        c.SetRightMargin(0.18);
        if (logx) c.SetLogx();
        hEff2D->SetTitle((TString::Format(
            "Real-electron efficiency (%s);Electron p_{T} [GeV];|Electron #eta|;#varepsilon_{r}",
            region)));
        hEff2D->GetZaxis()->SetRangeUser(yMin, yMax);
        hEff2D->Draw("colz text");
        TLatex t; t.SetNDC(); t.SetTextFont(72); t.SetTextSize(0.04);
        t.DrawLatex(0.16, 0.94, "ATLAS");
        t.SetTextFont(42); t.DrawLatex(0.26, 0.94, "Internal");
        c.SaveAs(TString::Format("%s/eff_real_2d_%s.%s", outDir, region, fmt));
    }

    fIn->Close();
}
