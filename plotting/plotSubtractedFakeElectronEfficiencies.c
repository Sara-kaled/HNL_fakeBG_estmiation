#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TEfficiency.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

// Helper function to subtract histograms (bin-by-bin, with error propagation, floor to 0)
TH1D* subtractHist(TH1D* hData, TH1D* hMC, const char* name) {
    if (!hData || !hMC) return nullptr;
    TH1D* hSub = (TH1D*)hData->Clone(name);
    for (int i = 1; i <= hSub->GetNbinsX(); ++i) {
        double val = max(0., hData->GetBinContent(i) - hMC->GetBinContent(i));
        double err = sqrt(pow(hData->GetBinError(i), 2) + pow(hMC->GetBinError(i), 2));
        hSub->SetBinContent(i, val);
        hSub->SetBinError(i, err);
    }
    return hSub;
}

TH2D* subtractHist2D(TH2D* hData, TH2D* hMC, const char* name) {
    if (!hData || !hMC) return nullptr;
    TH2D* hSub = (TH2D*)hData->Clone(name);
    for (int ix = 1; ix <= hSub->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hSub->GetNbinsY(); ++iy) {
            double val = max(0., hData->GetBinContent(ix, iy) - hMC->GetBinContent(ix, iy));
            double err = sqrt(pow(hData->GetBinError(ix, iy), 2) + pow(hMC->GetBinError(ix, iy), 2));
            hSub->SetBinContent(ix, iy, val);
            hSub->SetBinError(ix, iy, err);
        }
    }
    return hSub;
}

// Helper to create efficiency graph with Clopper-Pearson errors
TGraphAsymmErrors* makeEffGraph(TH1D* hNum, TH1D* hDen, const char* name, Color_t color = kBlack) {
    if (!hNum || !hDen) return nullptr;
    TGraphAsymmErrors* g = new TGraphAsymmErrors(hNum, hDen, "n");  // Clopper-Pearson
    g->SetName(name);
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(1.0);
    g->SetLineWidth(2);
    return g;
}

// Helper to divide 2D histograms for efficiency (with errors)
TH2D* divide2D(TH2D* hNum, TH2D* hDen, const char* name) {
    if (!hNum || !hDen) return nullptr;
    TH2D* hEff = (TH2D*)hDen->Clone(name);
    hEff->Reset();
    for (int ix = 1; ix <= hEff->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hEff->GetNbinsY(); ++iy) {
            double num = hNum->GetBinContent(ix, iy);
            double den = hDen->GetBinContent(ix, iy);
            double eff = (den > 0) ? num / den : 0;
            // Clopper-Pearson error (approx symmetric for simplicity; use TEfficiency for exact)
            double err = (den > 0) ? sqrt(eff * (1 - eff) / den) : 0;
            hEff->SetBinContent(ix, iy, eff);
            hEff->SetBinError(ix, iy, err);
        }
    }
    return hEff;
}

// Main function to plot 1D and 2D fake electron efficiencies with subtraction
void plotSubtractedFakeElectronEfficiencies(
    const char* mc_filename   = "../MC_eff/MC_histograms.root",
    const char* data_filename = "../data_eff/Data_histograms.root",
    const char* region        = "CR",  // e.g., "CR", "CR_MET_high", "SR", "incl"
    bool        plot2D        = true,           // Set to true for 2D plots
    bool        logx          = true,           // Log x-axis for pT
    float       yMin          = 0.0,            // Efficiency y-range min
    float       yMax          = 1.0,            // Efficiency y-range max
    const char* outDir        = ".",            // Output directory
    const char* figType       = "pdf"           // Output format
) {
    TFile* fMC   = TFile::Open(mc_filename, "READ");
    TFile* fData = TFile::Open(data_filename, "READ");
    if (!fMC || !fData) {
        cout << "Error opening files!" << endl;
        return;
    }

    string reg(region);

    // Histogram names
    string loose_pt  = "h_loose_pt_"  + reg;
    string tight_pt  = "h_tight_pt_"  + reg;
    string loose_eta = "h_loose_eta_" + reg;
    string tight_eta = "h_tight_eta_" + reg;
    string loose_2d  = "h2_loose_eta_pt_" + reg;
    string tight_2d  = "h2_tight_eta_pt_" + reg;

    // ── Load histograms ──────────────────────────────────────────────────────
    TH1D* hDataLoosePt  = (TH1D*)fData->Get(loose_pt.c_str());
    TH1D* hDataTightPt  = (TH1D*)fData->Get(tight_pt.c_str());
    TH1D* hMCLoosePt    = (TH1D*)fMC->Get(loose_pt.c_str());
    TH1D* hMCTightPt    = (TH1D*)fMC->Get(tight_pt.c_str());

    TH1D* hDataLooseEta = (TH1D*)fData->Get(loose_eta.c_str());
    TH1D* hDataTightEta = (TH1D*)fData->Get(tight_eta.c_str());
    TH1D* hMCLooseEta   = (TH1D*)fMC->Get(loose_eta.c_str());
    TH1D* hMCTightEta   = (TH1D*)fMC->Get(tight_eta.c_str());

    TH2D* hDataLoose2D  = nullptr;
    TH2D* hDataTight2D  = nullptr;
    TH2D* hMCLoose2D    = nullptr;
    TH2D* hMCTight2D    = nullptr;
    if (plot2D) {
        hDataLoose2D = (TH2D*)fData->Get(loose_2d.c_str());
        hDataTight2D = (TH2D*)fData->Get(tight_2d.c_str());
        hMCLoose2D   = (TH2D*)fMC->Get(loose_2d.c_str());
        hMCTight2D   = (TH2D*)fMC->Get(tight_2d.c_str());
    }

    // ── Subtracted numerators and denominators ───────────────────────────────
    TH1D* hNumPt  = subtractHist(hDataTightPt,  hMCTightPt,  "hNumPt");
    TH1D* hDenPt  = subtractHist(hDataLoosePt,  hMCLoosePt,  "hDenPt");
    TH1D* hNumEta = subtractHist(hDataTightEta, hMCTightEta, "hNumEta");
    TH1D* hDenEta = subtractHist(hDataLooseEta, hMCLooseEta, "hDenEta");

    // ── Efficiencies ─────────────────────────────────────────────────────────
    TGraphAsymmErrors* gEffPt  = makeEffGraph(hNumPt,  hDenPt,  "gEffPt");
    TGraphAsymmErrors* gEffEta = makeEffGraph(hNumEta, hDenEta, "gEffEta");

    TH2D* hEff2D = nullptr;
    if (plot2D && hDataLoose2D && hDataTight2D && hMCLoose2D && hMCTight2D) {
        TH2D* hNum2D = subtractHist2D(hDataTight2D, hMCTight2D, "hNum2D");
        TH2D* hDen2D = subtractHist2D(hDataLoose2D, hMCLoose2D, "hDen2D");
        hEff2D = divide2D(hNum2D, hDen2D, "hEff2D");
    }

    // ── Style setup ──────────────────────────────────────────────────────────
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPalette(104);  // Nice color palette for 2D

    // ── Plot 1D pT ───────────────────────────────────────────────────────────
    string cnamePt = "fake_el_eff_pt_" + reg;
    TCanvas* cPt = new TCanvas(cnamePt.c_str(), cnamePt.c_str(), 770, 560);
   // cPt->SetLogx(logx);
    TH1F* framePt = new TH1F("framePt", "", 100, 1, 80);  // Adjust bins as needed
    framePt->GetYaxis()->SetRangeUser(yMin, yMax);
    framePt->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
    framePt->GetYaxis()->SetTitle("Fake efficiency");
    framePt->Draw("AXIS");
    if (gEffPt) gEffPt->Draw("P SAME");

    TLegend* leg = new TLegend(0.58, 0.80, 0.88, 0.89);
    leg->AddEntry(gEffPt, "Data - MC (fake CR)", "P");
    leg->Draw();

    TLatex lat; lat.SetNDC(); lat.SetTextSize(0.03);
    lat.DrawLatex(0.18, 0.92, ("Fake electron efficiency - " + reg).c_str());

    cPt->SaveAs((string(outDir) + "/" + cnamePt + "." + figType).c_str());

    // ── Plot 1D eta ──────────────────────────────────────────────────────────
    string cnameEta = "fake_el_eff_eta_" + reg;
    TCanvas* cEta = new TCanvas(cnameEta.c_str(), cnameEta.c_str(), 770, 560);
    TH1F* frameEta = new TH1F("frameEta", "", 100, 0, 3.0);  // Adjust for |eta|
    frameEta->GetYaxis()->SetRangeUser(yMin, yMax);
    frameEta->GetXaxis()->SetTitle("Electron |#eta|");
    frameEta->GetYaxis()->SetTitle("Fake efficiency");
    frameEta->Draw("AXIS");
    if (gEffEta) gEffEta->Draw("P SAME");

    leg->Draw();
    lat.DrawLatex(0.18, 0.92, ("Fake electron efficiency - " + reg).c_str());

    cEta->SaveAs((string(outDir) + "/" + cnameEta + "." + figType).c_str());

    // ── Plot 2D ──────────────────────────────────────────────────────────────
    if (plot2D && hEff2D) {
        string cname2D = "fake_el_eff_2D_" + reg;
        TCanvas* c2D = new TCanvas(cname2D.c_str(), cname2D.c_str(), 770, 560);
        c2D->SetRightMargin(0.15);
        c2D->SetLogx(logx);

//        hEff2D->SetMinimum(0.001); 
        hEff2D->GetZaxis()->SetRangeUser(yMin, yMax);
        hEff2D->GetXaxis()->SetTitleOffset(1.25);   // 25 % farther from the axis

        hEff2D->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
        hEff2D->GetYaxis()->SetTitle("Electron |#eta|");
        hEff2D->GetZaxis()->SetTitle("Fake efficiency");
        hEff2D->Draw("COLZ TEXT45 E");

        gStyle->SetPaintTextFormat(".3f");


        lat.DrawLatex(0.18, 0.92, ("ATLAS #sqrt{s}=13 TeV, 140 fb^{-1} Fake electron efficiencies"));

        c2D->SaveAs((string(outDir) + "/" + cname2D + "." + figType).c_str());
    }

    cout << "Plots saved in " << outDir << " for region: " << region << endl;

    // Cleanup (optional)
    delete fMC;
    delete fData;
}

// Example usage:
// int main() {
//     plotSubtractedFakeElectronEfficiencies();
//     return 0;
// }

