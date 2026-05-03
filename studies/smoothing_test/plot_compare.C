// =============================================================================
// FILE:    studies/smoothing_test/plot_compare.C
// PURPOSE: Plot the per-event AMM weight distribution from data16, overlaying
//          raw-ε_f and smoothed-ε_f predictions. Also produces:
//            - Weight vs leading electron pT (profile of mean weight)
//            - Cumulative integral comparison (predicted yield vs pT cut)
//
// INPUT:   outputs/data16_amm_compare.root
// OUTPUT:  outputs/compare_smoothing.pdf
// =============================================================================
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <iostream>

void plot_compare()
{
    TFile* fIn = TFile::Open("outputs/data16_amm_compare.root", "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "ERROR: cannot open outputs/data16_amm_compare.root\n";
        return;
    }
    TTree* t = (TTree*)fIn->Get("cmp");
    if (!t) { std::cerr << "ERROR: tree 'cmp' missing\n"; return; }

    gStyle->SetOptStat(0);

    TCanvas c("c","smoothing comparison", 1400, 900);
    c.Divide(2, 2);

    // ----- Pad 1: per-event weight distribution -----
    c.cd(1);
    gPad->SetLogy();
    TH1D* hRaw    = new TH1D("hRaw",    ";AMM event weight;events", 100, -0.4, 0.6);
    TH1D* hSmooth = new TH1D("hSmooth", ";AMM event weight;events", 100, -0.4, 0.6);
    t->Draw("w_amm_raw    >> hRaw",    "", "goff");
    t->Draw("w_amm_smooth >> hSmooth", "", "goff");
    hRaw   ->SetLineColor(kBlack); hRaw   ->SetLineWidth(2);
    hSmooth->SetLineColor(kRed);   hSmooth->SetLineWidth(2);
    hRaw->Draw("HIST");
    hSmooth->Draw("HIST SAME");
    {
        TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
        leg->AddEntry(hRaw,    Form("raw     (#sum w = %.1f)", hRaw->Integral(0,-1)),    "l");
        leg->AddEntry(hSmooth, Form("smoothed (#sum w = %.1f)", hSmooth->Integral(0,-1)), "l");
        leg->SetBorderSize(0);
        leg->Draw();
    }

    // ----- Pad 2: profile of mean weight vs pT -----
    c.cd(2);
    TProfile* pRaw = new TProfile("pRaw", ";leading e p_{T} [GeV];#LTw_{AMM}#GT",
                                  20, 15., 200.);
    TProfile* pSm  = new TProfile("pSm",  ";leading e p_{T} [GeV];#LTw_{AMM}#GT",
                                  20, 15., 200.);
    t->Draw("w_amm_raw    : l2_pt >> pRaw", "", "prof goff");
    t->Draw("w_amm_smooth : l2_pt >> pSm",  "", "prof goff");
    pRaw->SetLineColor(kBlack); pRaw->SetMarkerColor(kBlack); pRaw->SetMarkerStyle(20);
    pSm ->SetLineColor(kRed);   pSm ->SetMarkerColor(kRed);   pSm ->SetMarkerStyle(20);
    pRaw->Draw("E1");
    pSm ->Draw("E1 SAME");
    {
        TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
        leg->AddEntry(pRaw, "raw",      "lp");
        leg->AddEntry(pSm,  "smoothed", "lp");
        leg->SetBorderSize(0);
        leg->Draw();
    }

    // ----- Pad 3: weight vs pT, restricted to loose-only events (positive weights) -----
    c.cd(3);
    TH1D* hPosRaw = new TH1D("hPosRaw", "loose-only events;leading e p_{T} [GeV];sum w_{L}",
                             20, 15., 200.);
    TH1D* hPosSm  = new TH1D("hPosSm",  "loose-only events;leading e p_{T} [GeV];sum w_{L}",
                             20, 15., 200.);
    t->Draw("l2_pt >> hPosRaw", "w_amm_raw    * (l2_loose && !l2_tight)", "goff");
    t->Draw("l2_pt >> hPosSm",  "w_amm_smooth * (l2_loose && !l2_tight)", "goff");
    hPosRaw->SetLineColor(kBlack); hPosRaw->SetLineWidth(2);
    hPosSm ->SetLineColor(kRed);   hPosSm ->SetLineWidth(2);
    hPosRaw->Draw("HIST");
    hPosSm ->Draw("HIST SAME");
    {
        TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
        leg->AddEntry(hPosRaw, "raw",      "l");
        leg->AddEntry(hPosSm,  "smoothed", "l");
        leg->SetBorderSize(0);
        leg->Draw();
    }

    // ----- Pad 4: cumulative-integral vs pT-cut -----
    c.cd(4);
    TH1D* hCumRaw = (TH1D*)hPosRaw->Clone("hCumRaw");
    TH1D* hCumSm  = (TH1D*)hPosSm ->Clone("hCumSm");
    // build right-cumulative: integral from pT bin onwards
    {
        const int n = hCumRaw->GetNbinsX();
        for (int i = 1; i <= n; ++i) {
            hCumRaw->SetBinContent(i, hPosRaw->Integral(i, n));
            hCumSm ->SetBinContent(i, hPosSm ->Integral(i, n));
        }
    }
    hCumRaw->SetTitle(";leading e p_{T} cut [GeV];loose-only fake yield above cut");
    hCumRaw->Draw("HIST");
    hCumSm ->Draw("HIST SAME");
    {
        TLegend* leg = new TLegend(0.55, 0.75, 0.88, 0.88);
        leg->AddEntry(hCumRaw, "raw",      "l");
        leg->AddEntry(hCumSm,  "smoothed", "l");
        leg->SetBorderSize(0);
        leg->Draw();
    }

    c.SaveAs("outputs/compare_smoothing.pdf");
    std::cout << "\nIntegrated AMM yield (data16):\n";
    std::cout << "  raw     : " << hRaw   ->GetMean() * hRaw   ->GetEntries() << "\n";
    std::cout << "  smoothed: " << hSmooth->GetMean() * hSmooth->GetEntries() << "\n";
    std::cout << "\n✓ Wrote outputs/compare_smoothing.pdf\n";

    fIn->Close();
}
