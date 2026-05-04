// =============================================================================
// FILE:    efficiency/plot_composition_cr_sr.C
// PURPOSE: Produce two stacked-fraction plots (one for CR, one for SR) of
//          the truth-source composition of loose fakes vs pT, using the
//          frac_<src>_CR / frac_<src>_SR histograms already written by
//          build_composition_systematic.C into composition_systematic.root.
//
// OUTPUT:  outputs/composition_cr.pdf
//          outputs/composition_sr.pdf
// =============================================================================
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>

#include "../utils/AtlasStyle.h"

static const char* kSrcName[] = {"B","C","Conv","LF","Tau","Unknown","Other"};
static const int   kNSrc      = 7;

static int colorOf(int s) {
    switch (s) {
        case 0: return kATLAS_src_B;
        case 1: return kATLAS_src_C;
        case 2: return kATLAS_src_Conv;
        case 3: return kATLAS_src_LF;
        case 4: return kATLAS_src_Tau;
        default: return kATLAS_src_Other;
    }
}

static void drawOne(TFile* f, const char* tag, const char* title,
                    const char* outPath)
{
    SetATLASStyle();
    TCanvas c("c","",700,550);

    THStack* hs = new THStack(Form("hs_%s",tag),
                              Form(";p_{T}^{e} [GeV];fraction"));

    // We display fractions; one TH1D per source built from frac_<src>_<tag>.
    TH1D* hOne[kNSrc] = {nullptr};
    TLegend* leg = new TLegend(0.75, 0.55, 0.95, 0.92);
    leg->SetBorderSize(0); leg->SetFillStyle(0);

    for (int s = 0; s < kNSrc; ++s) {
        TString name = TString::Format("frac_%s_%s", kSrcName[s], tag);
        auto* h2 = (TH2D*)f->Get(name);
        if (!h2) continue;
        // Project onto pT (X) integrating |η| (Y); average over η bins
        TH1D* h = h2->ProjectionX(Form("p_%s_%s", kSrcName[s], tag));
        h->Scale(1.0 / h2->GetNbinsY());     // mean over η bins
        h->SetFillColor(colorOf(s));
        h->SetLineColor(colorOf(s));
        h->SetMarkerColor(colorOf(s));
        hOne[s] = h;
        hs->Add(h);
        leg->AddEntry(h, kSrcName[s], "f");
    }

    hs->Draw("HIST");
    hs->SetMaximum(1.2);
    hs->GetYaxis()->SetTitleOffset(1.2);
    leg->Draw();

    DrawAtlasLabel(0.18, 0.88);
    TLatex lat; lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.04);
    lat.DrawLatex(0.18, 0.80, title);

    c.SaveAs(outPath);
}

void plot_composition_cr_sr()
{
    gSystem->mkdir("outputs", true);

    TFile* f = TFile::Open("outputs/composition_systematic.root", "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "ERROR: outputs/composition_systematic.root not found.\n"
                     "       Run efficiency/build_composition_systematic.C first.\n";
        return;
    }

    drawOne(f, "CR", "SS #mue Fake CR (truth fakes)",     "outputs/composition_cr.pdf");
    drawOne(f, "SR", "OS VBF SR (truth fakes)",            "outputs/composition_sr.pdf");

    f->Close();
    std::cout << "Wrote outputs/composition_cr.pdf and outputs/composition_sr.pdf\n";
}
