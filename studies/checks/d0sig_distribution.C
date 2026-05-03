// =============================================================================
// FILE:    studies/d0sig_distribution.C
// PURPOSE: Compare |d0sig| distributions of loose electrons in MC, split by
//          truth source, per pT bin. Heavy-flavor real e's (b/c → e) have
//          large |d0sig| due to parent-hadron flight; light-flavor jet→e
//          fakes have small |d0sig|. The shape of the high-pT loose
//          population's |d0sig| tells us if "fakes" there are mostly real
//          heavy-flavor electrons.
//
// OUTPUT:  outputs/d0sig_distribution.root
//          outputs/d0sig_distribution.pdf  (one panel per pT bin)
//          stdout: per-pT mean |d0sig| split by truth source
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include "ROOT/RVec.hxx"
#include <iostream>
#include <cstdio>

#include "../../utils/SelectionUtils.h"

using ROOT::VecOps::RVec;

enum Cat { kPrompt=0, kHF=1, kLF=2, kConv=3, kOther=4, kNCat=5 };
static const char* kCatName[kNCat] = {"prompt","HF (b/c)","LF","Conv","Other"};
static const int   kCatColor[kNCat] = {kBlack, kRed+1, kGreen+2, kBlue+1, kGray+1};

static int classify(int iff, bool isPrompt) {
    if (isPrompt) return kPrompt;
    if (iff == 8 || iff == 9) return kHF;
    if (iff == 10) return kLF;
    if (iff == 5)  return kConv;
    return kOther;
}

void d0sig_distribution()
{
    gStyle->SetOptStat(0);
    gSystem->mkdir("outputs", true);

    TChain mc("tree");
    const std::string base = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    for (auto p : {"ttbar","ztautau","wjets","singletop","diboson","higgs"})
        mc.Add((base + p + "/*.root").c_str());

    EventData ev;
    SetupBranches(&mc, ev);
    RVecI* el_IFFClass = nullptr;
    mc.SetBranchAddress("el_IFFClass", &el_IFFClass);

    // h[pT-bin][cat]
    std::vector<std::vector<TH1D*>> h(kNPtBins, std::vector<TH1D*>(kNCat, nullptr));
    for (int ix = 0; ix < kNPtBins; ++ix) {
        for (int s = 0; s < kNCat; ++s) {
            h[ix][s] = new TH1D(Form("h_d0sig_pt%d_%s", ix, kCatName[s]),
                                Form(";|d_{0}/#sigma| (loose e, p_{T} #in [%.0f,%.0f]);electrons",
                                     kPtBins[ix], kPtBins[ix+1]),
                                40, 0., 8.);
            h[ix][s]->Sumw2();
            h[ix][s]->SetLineColor(kCatColor[s]);
            h[ix][s]->SetLineWidth(2);
        }
    }

    const Long64_t n = mc.GetEntries();
    std::cout << "MC entries: " << n << "\n";
    for (Long64_t i = 0; i < n; ++i) {
        if (i % 1000000 == 0) std::cout << "  " << i << " / " << n << "\n";
        mc.GetEntry(i);
        ComputeDerived(ev);
        if (!passBasePresel_FakeMeasurement(ev)) continue;
        if (!isFakeCR(ev)) continue;
        if (!ev.no_bjets) continue;
        if (!ev.el_is_loose) continue;
        if (!ev.el_d0sig || ev.el_d0sig->empty()) continue;
        if (!el_IFFClass || el_IFFClass->empty()) continue;

        const double d0sig = std::fabs((*ev.el_d0sig)[0]);
        const double pt    = ev.l2_pt * 1e-3;
        int ix = -1;
        for (int b = 0; b < kNPtBins; ++b)
            if (pt >= kPtBins[b] && pt < kPtBins[b+1]) { ix = b; break; }
        if (ix < 0) continue;

        const int cat = classify((*el_IFFClass)[0], isMCPromptElectron(ev));
        h[ix][cat]->Fill(d0sig, mcWeightLoose(ev));
    }

    // ----- Per pT mean |d0sig| table -----
    printf("\n========================================================================\n");
    printf("  Mean |d0sig| of loose electrons split by truth category\n");
    printf("========================================================================\n");
    printf("  pT [GeV]      ");
    for (int s = 0; s < kNCat; ++s) printf("%9s   ", kCatName[s]);
    printf("\n  -------------------------------------------------------------------\n");
    for (int ix = 0; ix < kNPtBins; ++ix) {
        printf("  [%4.0f,%4.0f]  ", kPtBins[ix], kPtBins[ix+1]);
        for (int s = 0; s < kNCat; ++s)
            printf("  %.2f      ", h[ix][s]->GetMean());
        printf("\n");
    }
    printf("\n  Interpretation:\n");
    printf("    HF (b/c) typically has <|d0sig|> ~ 2-4 (displaced from PV).\n");
    printf("    LF mis-ID and prompt have <|d0sig|> < 1 (compatible with PV).\n");
    printf("    If your high-pT bin's mean shifts toward HF, those 'fakes' are real e's.\n");

    // ----- Plot: one panel per pT bin -----
    TCanvas c("c","|d0sig| split", 1500, 1000);
    c.Divide(3, (kNPtBins+2)/3);
    for (int ix = 0; ix < kNPtBins; ++ix) {
        c.cd(ix+1); gPad->SetLogy();
        // Stack just by drawing on top
        double maxY = 0.;
        for (int s = 0; s < kNCat; ++s) maxY = std::max(maxY, h[ix][s]->GetMaximum());
        for (int s = 0; s < kNCat; ++s) {
            h[ix][s]->GetYaxis()->SetRangeUser(0.5, 5.0 * maxY);
            h[ix][s]->Draw(s == 0 ? "HIST" : "HIST SAME");
        }
        if (ix == 0) {
            auto* leg = new TLegend(0.55, 0.55, 0.88, 0.88);
            for (int s = 0; s < kNCat; ++s) leg->AddEntry(h[ix][s], kCatName[s], "l");
            leg->SetBorderSize(0); leg->Draw();
        }
    }
    c.SaveAs("outputs/d0sig_distribution.pdf");

    TFile* fo = TFile::Open("outputs/d0sig_distribution.root", "RECREATE");
    for (int ix = 0; ix < kNPtBins; ++ix)
        for (int s = 0; s < kNCat; ++s) h[ix][s]->Write();
    fo->Close();
    printf("\nWrote outputs/d0sig_distribution.{root,pdf}\n");
}
