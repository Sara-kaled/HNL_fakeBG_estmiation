// =============================================================================
// FILE:    studies/data_mc_loose_ratio.C
// PURPOSE: For each (pT, |η|) bin in the SS μe fake CR, compare the data
//          loose-electron yield to the MC prompt-only loose yield.
//          The ratio data_L / MC_L tells us how much "extra" the data has
//          relative to MC's prompt model.
//
//          Interpretation:
//            ratio ~ 1   : MC tracks data → prompt subtraction is essentially
//                          all of data, and the "fake residual" is small.
//                          Over-statement of ε_f then comes from a real
//                          composition difference, not unmodeled fakes.
//            ratio >> 1  : data has much more loose than MC → unmodeled
//                          fake population (likely QCD multijet) inflating
//                          the "data−MC" denominator. Apply a MET cut or
//                          tighten the CR to reduce this contamination.
//            ratio < 1   : MC over-predicts loose population (overestimated
//                          xs, wrong scale factors, etc.) — investigate.
//
// OUTPUT:  outputs/data_mc_loose_ratio.root
//          outputs/data_mc_loose_ratio.pdf
//          stdout: per-bin yield + ratio table
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include <cstdio>

#include "../../utils/SelectionUtils.h"

void data_mc_loose_ratio()
{
    gStyle->SetOptStat(0);
    gSystem->mkdir("outputs", true);

    TChain dataCh("tree");
    dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data15.root");
    dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data16.root");
    dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_d/data17.root");
    dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/data18.root");

    TChain mcCh("tree");
    const std::string base = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    for (auto p : {"ttbar","ztautau","wjets","singletop","diboson","higgs"})
        mcCh.Add((base + p + "/*.root").c_str());

    auto book = [](const char* nm) {
        return new TH2D(nm, ";p_{T}^{e} [GeV];|#eta^{e}|",
                        kNPtBins, kPtBins, kNEtaBins, kEtaBins);
    };
    TH2D *hData = book("hData_loose"), *hMC = book("hMC_loose");
    hData->Sumw2(); hMC->Sumw2();

    auto loop = [&](TTree* t, bool isMC, const char* tag) {
        std::cout << "\n=== Loop " << tag << " (entries " << t->GetEntries() << ") ===\n";
        EventData ev;
        SetupBranches(t, ev);
        const Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            if (i % 1000000 == 0) std::cout << "  " << i << " / " << n << "\n";
            t->GetEntry(i);
            ComputeDerived(ev);
            if (!passBasePresel_FakeMeasurement(ev)) continue;
            if (!isFakeCR(ev)) continue;
            if (!ev.no_bjets) continue;
            if (!ev.el_is_loose) continue;
            if (isMC && !isMCPromptElectron(ev)) continue;
            const double w = isMC ? mcWeightLoose(ev) : 1.0;
            (isMC ? hMC : hData)->Fill(ev.l2_pt * 1e-3, ev.el_aeta, w);
        }
    };
    loop(&mcCh,   true,  "MC");
    loop(&dataCh, false, "DATA");

    // ----- Per-bin table -----
    printf("\n========================================================================\n");
    printf("  data_loose / MC_prompt_loose per (pT, |η|) bin\n");
    printf("========================================================================\n");
    printf("  pT [GeV]      |η|        N_data        N_MC          ratio\n");
    printf("  ----------------------------------------------------------------\n");
    for (int ix = 1; ix <= kNPtBins; ++ix) {
        for (int iy = 1; iy <= kNEtaBins; ++iy) {
            double d = hData->GetBinContent(ix, iy);
            double m = hMC  ->GetBinContent(ix, iy);
            double r = (m > 0) ? d / m : 0.;
            printf("  [%4.0f,%4.0f]  [%.2f,%.2f]  %10.1f  %10.1f   %6.2f\n",
                   kPtBins[ix-1], kPtBins[ix],
                   kEtaBins[iy-1], kEtaBins[iy], d, m, r);
        }
    }

    // ----- Average per pT (collapse over η) -----
    printf("\n  Averages per pT bin:\n");
    printf("  pT [GeV]      <data>      <MC_prompt>   <ratio>\n");
    printf("  ----------------------------------------------------------------\n");
    TH1D* hRatPt = new TH1D("hRatPt", ";p_{T} [GeV];data_{loose}/MC_{loose}",
                            kNPtBins, kPtBins);
    for (int ix = 1; ix <= kNPtBins; ++ix) {
        double d = 0., m = 0.;
        for (int iy = 1; iy <= kNEtaBins; ++iy) {
            d += hData->GetBinContent(ix, iy);
            m += hMC  ->GetBinContent(ix, iy);
        }
        double r = (m > 0) ? d / m : 0.;
        printf("  [%4.0f,%4.0f]  %10.1f  %10.1f   %6.2f\n",
               kPtBins[ix-1], kPtBins[ix], d, m, r);
        hRatPt->SetBinContent(ix, r);
    }

    // ----- Plot -----
    TCanvas c("c","data/MC ratio", 900, 600);
    hRatPt->SetMarkerStyle(20); hRatPt->SetMarkerSize(1.4);
    hRatPt->SetLineWidth(2);
    hRatPt->Draw("PL");
    hRatPt->GetYaxis()->SetRangeUser(0., 1.5 * hRatPt->GetMaximum());
    c.SaveAs("outputs/data_mc_loose_ratio.pdf");

    TFile* fo = TFile::Open("outputs/data_mc_loose_ratio.root", "RECREATE");
    hData->Write(); hMC->Write(); hRatPt->Write();
    fo->Close();
    printf("\nWrote outputs/data_mc_loose_ratio.{root,pdf}\n");
}
