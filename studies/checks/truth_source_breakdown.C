// =============================================================================
// FILE:    studies/truth_source_breakdown.C
// PURPOSE: Per-pT truth-source composition of the loose fake-electron
//          population in the SS μe fake CR. Uses el_IFFClass to classify
//          each non-prompt electron into:
//             B (8), C (9), Conv (5), LF (10), Tau (7), Other
//          Tells us at high pT whether "fakes" are mostly heavy-flavor real
//          electrons (in which case the matrix method is fundamentally
//          struggling) or genuine jet→e mis-IDs.
//
// OUTPUT:  outputs/truth_source_breakdown.root
//          outputs/truth_source_breakdown.pdf
//          stdout: per-pT-bin fraction table
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include "ROOT/RVec.hxx"
#include <iostream>
#include <cstdio>

#include "../../utils/SelectionUtils.h"

using ROOT::VecOps::RVec;

enum FakeSrc { kB=0, kC=1, kConv=2, kLF=3, kTau=4, kOther=5, kNSrc=6 };
static const char* kSrcName[kNSrc] = {"B","C","Conv","LF","Tau","Other"};
static const int   kSrcColor[kNSrc] = {kRed+1, kOrange+1, kBlue+1, kGreen+2, kMagenta+1, kGray+1};

static int classifyIFF(int iff) {
    switch(iff) {
        case 8:  return kB;
        case 9:  return kC;
        case 5:  return kConv;
        case 10: return kLF;
        case 7:  return kTau;
        default: return kOther;
    }
}

void truth_source_breakdown()
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

    // One histogram per truth source, x = pT
    TH1D* h[kNSrc];
    for (int s = 0; s < kNSrc; ++s) {
        h[s] = new TH1D(Form("h_%s", kSrcName[s]),
                        ";p_{T}^{e} [GeV];loose fakes (weighted)",
                        kNPtBins, kPtBins);
        h[s]->Sumw2();
        h[s]->SetFillColor(kSrcColor[s]);
        h[s]->SetLineColor(kSrcColor[s]);
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
        if (isMCPromptElectron(ev)) continue;   // restrict to truth fakes only
        if (!el_IFFClass || el_IFFClass->empty()) continue;

        const int src = classifyIFF((*el_IFFClass)[0]);
        h[src]->Fill(ev.l2_pt * 1e-3, mcWeightLoose(ev));
    }

    // ----- Per-pT fraction table -----
    printf("\n========================================================================\n");
    printf("  Truth-source composition of loose fakes per pT bin\n");
    printf("========================================================================\n");
    printf("  pT bin       ");
    for (int s = 0; s < kNSrc; ++s) printf("  %-6s", kSrcName[s]);
    printf("    total\n");
    printf("  -------------------------------------------------------------------\n");
    for (int ix = 1; ix <= kNPtBins; ++ix) {
        double tot = 0.;
        for (int s = 0; s < kNSrc; ++s) tot += h[s]->GetBinContent(ix);
        printf("  [%4.0f,%4.0f]  ", kPtBins[ix-1], kPtBins[ix]);
        for (int s = 0; s < kNSrc; ++s) {
            double f = (tot > 0) ? h[s]->GetBinContent(ix)/tot : 0.;
            printf("  %5.1f%%", 100.*f);
        }
        printf("    %.1f\n", tot);
    }
    printf("\n  Interpretation:\n");
    printf("    B/C dominant at high pT  → heavy-flavor real e's; matrix method assumes\n");
    printf("                                they're fakes but they pass tight at near-100%%.\n");
    printf("                                CAP pT or use a fake-factor with looser ID.\n");
    printf("    LF dominant at low pT    → genuine jet→e; matrix method clean.\n");
    printf("    Conv dominant somewhere  → photon-to-e^{±}; consider explicit handling.\n");

    // ----- Stacked histogram -----
    TCanvas c("c","Truth-source", 900, 600);
    THStack* hs = new THStack("hs", ";p_{T}^{e} [GeV];fraction of loose fakes");
    // Normalize each pT bin to unity to show fractions
    TH1D* hSum = (TH1D*)h[0]->Clone("hSum"); hSum->Reset();
    for (int s = 0; s < kNSrc; ++s) hSum->Add(h[s]);
    for (int s = 0; s < kNSrc; ++s) {
        TH1D* hN = (TH1D*)h[s]->Clone(Form("hN_%s", kSrcName[s]));
        hN->Divide(hSum);
        hs->Add(hN);
    }
    hs->Draw("HIST");
    hs->GetYaxis()->SetRangeUser(0., 1.0);
    auto* leg = new TLegend(0.78, 0.55, 0.95, 0.88);
    for (int s = 0; s < kNSrc; ++s) leg->AddEntry(h[s], kSrcName[s], "f");
    leg->SetBorderSize(0); leg->Draw();
    c.SaveAs("outputs/truth_source_breakdown.pdf");

    TFile* fo = TFile::Open("outputs/truth_source_breakdown.root", "RECREATE");
    for (int s = 0; s < kNSrc; ++s) h[s]->Write();
    fo->Close();
    printf("\nWrote outputs/truth_source_breakdown.{root,pdf}\n");
}
