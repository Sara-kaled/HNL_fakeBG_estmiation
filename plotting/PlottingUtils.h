// PlottingUtils.h — shared colors, labels, canvas helpers for all study macros
#pragma once
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TLine.h>
#include <TFile.h>
#include <TSystem.h>
#include <TH1D.h>
#include <TH2F.h>
#include <THStack.h>
#include <TGraphErrors.h>
#include <map>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/stat.h>

// -------------------------------------------------------------------------
// Sample metadata
// -------------------------------------------------------------------------
struct SampleInfo { int color; const char* title; };

inline std::map<std::string, SampleInfo> sampleMeta() {
    return {
        {"ttbar",     {kGray+2,      "#it{t#bar{t}}"}},
        {"ztautau",   {kAzure-4,     "Z#rightarrow#tau#tau"}},
        {"diboson",   {kGreen-6,     "Diboson"}},
        {"singletop", {kOrange+1,    "Single top"}},
        {"wjets",     {kYellow-3,    "W+jets"}},
        {"higgs",     {kViolet-4,    "Higgs"}},
        {"fake",      {kRed-7,       "Fake/Non-Prompt"}},
    };
}

// -------------------------------------------------------------------------
// Minimal style
// -------------------------------------------------------------------------
inline void setStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetPadLeftMargin(0.13);
    gStyle->SetLabelSize(0.04,"xyz");
    gStyle->SetTitleSize(0.045,"xyz");
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetHistLineWidth(2);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
}

// -------------------------------------------------------------------------
// ATLAS label
// -------------------------------------------------------------------------
inline void drawATLASLabel(double x=0.18, double y=0.88, bool internal=true) {
    TLatex l; l.SetNDC();
    l.SetTextFont(72); l.SetTextSize(0.045);
    l.DrawLatex(x, y, "ATLAS");
    l.SetTextFont(42);
    if (internal) l.DrawLatex(x+0.12, y, "Internal");
    l.SetTextSize(0.035);
    l.DrawLatex(x, y-0.06, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
}

// -------------------------------------------------------------------------
// Make output directory
// -------------------------------------------------------------------------
inline void makeDir(const char* path) {
    mkdir(path, 0755);
}

// -------------------------------------------------------------------------
// Build stacked canvas with ratio panel
// Returns: {upper pad, lower pad} — caller must switch between them
// -------------------------------------------------------------------------
inline std::pair<TPad*,TPad*> makeRatioCanvas(TCanvas* c) {
    TPad* upper = new TPad("upper","upper",0,0.28,1,1);
    upper->SetBottomMargin(0.02);
    upper->SetTopMargin(0.06);
    upper->SetRightMargin(0.05);
    upper->SetLeftMargin(0.13);
    upper->Draw(); upper->cd();
    c->cd();
    TPad* lower = new TPad("lower","lower",0,0,1,0.28);
    lower->SetTopMargin(0.03);
    lower->SetBottomMargin(0.38);
    lower->SetRightMargin(0.05);
    lower->SetLeftMargin(0.13);
    lower->Draw();
    return {upper, lower};
}

// -------------------------------------------------------------------------
// Draw stack + data + uncertainty band on current pad
// sortedNames: process names from largest to smallest (largest drawn first = bottom)
// hMC: map of process → histogram (already weighted/scaled)
// hData: observed data (nullptr if blinded)
// hFake: fake histogram (nullptr if not available)
// hTotalUnc: total MC uncertainty band (stat only for postfit)
// -------------------------------------------------------------------------
inline TH1D* drawStack(
    const std::vector<std::string>& sortedNames,
    std::map<std::string, TH1D*>& hMC,
    TH1D* hFake,
    TH1D* hData,
    TH1D* hTotalUnc,
    TLegend* leg,
    const char* xTitle,
    const char* yTitle = "Events")
{
    auto meta = sampleMeta();
    THStack* stack = new THStack("stack","");

    // Add in reverse order so largest ends up at bottom
    for (int i = (int)sortedNames.size()-1; i >= 0; --i) {
        const std::string& n = sortedNames[i];
        if (hMC.count(n) && hMC[n]) {
            hMC[n]->SetFillColor(meta[n].color);
            hMC[n]->SetLineColor(kBlack);
            hMC[n]->SetLineWidth(1);
            stack->Add(hMC[n]);
        }
    }
    // Fake always on top of MC stack
    if (hFake) {
        hFake->SetFillColor(meta["fake"].color);
        hFake->SetLineColor(kBlack);
        hFake->SetLineWidth(1);
        stack->Add(hFake);
    }

    stack->Draw("HIST");
    stack->GetXaxis()->SetTitle(xTitle);
    stack->GetYaxis()->SetTitle(yTitle);
    stack->GetXaxis()->SetLabelSize(0.);  // hidden if ratio pad below

    // Uncertainty band
    if (hTotalUnc) {
        hTotalUnc->SetFillStyle(3354);
        hTotalUnc->SetFillColor(kBlack);
        hTotalUnc->SetMarkerSize(0);
        hTotalUnc->SetLineWidth(0);
        hTotalUnc->Draw("E2 SAME");
        leg->AddEntry(hTotalUnc, "Stat. unc.", "f");
    }

    // Data
    if (hData) {
        hData->SetMarkerStyle(20);
        hData->SetMarkerSize(0.9);
        hData->SetLineColor(kBlack);
        hData->Draw("E SAME");
        leg->AddEntry(hData, "Data", "pe");
    }

    // Legend entries (top to bottom = smallest to largest)
    if (hFake) leg->AddEntry(hFake, meta["fake"].title, "f");
    for (const auto& n : sortedNames) {
        if (hMC.count(n) && hMC[n] && hMC[n]->Integral() > 0)
            leg->AddEntry(hMC[n], meta[n].title, "f");
    }
    leg->Draw();

    return (TH1D*)stack->GetHistogram();
}

// -------------------------------------------------------------------------
// Draw ratio pad
// -------------------------------------------------------------------------
inline void drawRatio(TPad* lower, TH1D* hData, TH1D* hTotalMC, const char* xTitle) {
    lower->cd();
    TH1D* hRatio = (TH1D*)hData->Clone("ratio");
    hRatio->Divide(hTotalMC);
    hRatio->SetMarkerStyle(20);
    hRatio->SetMarkerSize(0.8);
    hRatio->SetLineColor(kBlack);
    hRatio->GetYaxis()->SetTitle("Data / Pred.");
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.13);
    hRatio->GetYaxis()->SetLabelSize(0.12);
    hRatio->GetYaxis()->SetTitleOffset(0.4);
    hRatio->GetXaxis()->SetTitle(xTitle);
    hRatio->GetXaxis()->SetTitleSize(0.14);
    hRatio->GetXaxis()->SetLabelSize(0.12);
    hRatio->GetYaxis()->SetRangeUser(0.5, 1.5);
    hRatio->Draw("E");

    TLine* line = new TLine(hRatio->GetXaxis()->GetXmin(),1,hRatio->GetXaxis()->GetXmax(),1);
    line->SetLineStyle(2); line->SetLineColor(kRed);
    line->Draw();
}
