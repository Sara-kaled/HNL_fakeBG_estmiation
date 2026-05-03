// =============================================================================
// FILE:    utils/AtlasStyle.h
// PURPOSE: Single source of ATLAS plotting style + labels for the framework.
//          Replaces the duplicated SetATLASStyle / DrawAtlasLabel definitions
//          scattered across plotting/* and studies/* macros.
//
// USAGE:   #include "../utils/AtlasStyle.h"
//          SetATLASStyle();           // call at the top of each macro
//          DrawAtlasLabel(0.16, 0.88, "Internal");   // ATLAS Internal banner
//
// COLORS:  ATLAS-recommended palette for stack plots and overlays.
//          Use kATLAS<Process> below instead of raw kRed/kBlue/etc.
// =============================================================================
#ifndef ATLAS_STYLE_H_
#define ATLAS_STYLE_H_

#include <TStyle.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TColor.h>
#include <TPad.h>
#include <TH1.h>

// Full Run-2 13 TeV luminosity in fb^-1
constexpr const char* kATLASEnergyLumi = "#sqrt{s} = 13 TeV, 140 fb^{-1}";
constexpr const char* kATLASStatus     = "Internal";

// ATLAS-style colors for common processes / sources.
// Tuned for visibility on white backgrounds and consistency with ATLAS public
// plots (ttbar = orange/yellow, Z = green, single-top = light blue, diboson =
// magenta, Higgs = red, Wjets = teal, QCD/multijet = grey).
inline int kATLAS_ttbar    = TColor::GetColor("#F8C471");   // light orange
inline int kATLAS_singletop= TColor::GetColor("#85C1E9");   // light blue
inline int kATLAS_Wjets    = TColor::GetColor("#76D7C4");   // teal
inline int kATLAS_Zee      = TColor::GetColor("#82E0AA");   // green
inline int kATLAS_Zmm      = TColor::GetColor("#58D68D");   // brighter green
inline int kATLAS_Ztautau  = TColor::GetColor("#52BE80");   // dark green
inline int kATLAS_diboson  = TColor::GetColor("#BB8FCE");   // light purple
inline int kATLAS_Higgs    = TColor::GetColor("#F1948A");   // salmon
inline int kATLAS_fake     = TColor::GetColor("#5DADE2");   // bright blue
inline int kATLAS_data     = kBlack;

// Truth fake-source colors (B / C / Conv / LF / Tau / Other)
inline int kATLAS_src_B    = TColor::GetColor("#E74C3C");   // red
inline int kATLAS_src_C    = TColor::GetColor("#E67E22");   // orange
inline int kATLAS_src_Conv = TColor::GetColor("#3498DB");   // blue
inline int kATLAS_src_LF   = TColor::GetColor("#27AE60");   // green
inline int kATLAS_src_Tau  = TColor::GetColor("#9B59B6");   // purple
inline int kATLAS_src_Other= TColor::GetColor("#95A5A6");   // grey

inline void SetATLASStyle() {
    static bool done = false;
    if (done) return;
    done = true;

    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetCanvasColor(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameLineWidth(1);

    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    gStyle->SetPadLeftMargin  (0.14);
    gStyle->SetPadRightMargin (0.05);
    gStyle->SetPadTopMargin   (0.07);
    gStyle->SetPadBottomMargin(0.13);

    gStyle->SetTextFont(42);
    gStyle->SetTextSize(0.045);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetLabelSize(0.045, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetTitleSize(0.05, "xyz");
    gStyle->SetTitleOffset(1.1, "x");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.3, "z");

    gStyle->SetEndErrorSize(2);
    gStyle->SetMarkerStyle(20);
    gStyle->SetMarkerSize(1.0);

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillColor(0);
    gStyle->SetLegendFont(42);
    gStyle->SetLegendTextSize(0.04);

    gStyle->SetPalette(57);   // viridis (good for 2D maps)

    gROOT->ForceStyle();
}

inline void DrawAtlasLabel(double x = 0.16, double y = 0.88,
                           const char* status = kATLASStatus,
                           const char* energyLumi = kATLASEnergyLumi)
{
    TLatex l;
    l.SetNDC();
    l.SetTextFont(72);
    l.SetTextSize(0.05);
    l.DrawLatex(x, y, "ATLAS");
    l.SetTextFont(42);
    l.DrawLatex(x + 0.11, y, status);
    l.SetTextSize(0.04);
    l.DrawLatex(x, y - 0.05, energyLumi);
}

// Convenience: histogram styling shortcut
inline void StyleHist(TH1* h, int color, int markerStyle = 20,
                      double markerSize = 1.0, int lineWidth = 2)
{
    if (!h) return;
    h->SetLineColor  (color);
    h->SetMarkerColor(color);
    h->SetMarkerStyle(markerStyle);
    h->SetMarkerSize (markerSize);
    h->SetLineWidth  (lineWidth);
}

#endif // ATLAS_STYLE_H_
