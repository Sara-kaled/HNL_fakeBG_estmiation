// =============================================================================
// FILE:    studies/multijet_test.C
// PURPOSE: Test whether multijet (QCD) contamination is inflating the
//          measured ε_f. Idea:
//            - "Standard CR"   : SS μe with TIGHT muon  (current measurement)
//            - "Anti-iso CR"   : SS μe with LOOSE-ONLY muon (anti-isolated)
//                                → suppresses prompt μ from W/Z, enriches QCD
//                                  bb̄ → μe events.
//          Measure ε_f(pT, |η|) in BOTH regions on (data − MC_prompt).
//          If both ε_f match → no significant QCD contamination, your
//          standard CR ε_f is clean.
//          If anti-iso CR ε_f is much higher → QCD multijet really IS in your
//          standard CR; the data-MC residual contains an unmodeled multijet
//          fake population that inflates ε_f.
//
// INPUTS:  Data ntuples (data15-18)
//          MC ntuples (ttbar, wjets, ztautau, singletop, diboson, higgs)
// OUTPUT:  studies/outputs/multijet_test.root
//          studies/outputs/multijet_test.pdf  (per-η ε_f overlay)
//          stdout: per-bin ε_f comparison table
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include "../../utils/SelectionUtils.h"

// ---------------------------------------------------------------------------
// Mirror of passBasePresel WITHOUT the mu_is_tight requirement, so we can
// classify events into tight-mu vs loose-only-mu downstream. Keeps every
// other cut identical to the standard preselection.
// ---------------------------------------------------------------------------
static inline bool passPreselNoMuTight(const EventData& ev)
{
    return passTopology(ev)
        && passTrigger(ev)
        && passTriggerMatching(ev)
        && passMuonIP(ev)
        && passElectronIP(ev)
        && passCrackVeto(ev)
        && passLeptonPt(ev);
}

// ---------------------------------------------------------------------------
// Helper: fill a (pT, |η|) histogram for tight / loose electron at the bin
// containing ev.l2_pt (GeV) × ev.el_aeta. Both 2D histograms have the same
// kPtBins × kEtaBins binning so they divide cleanly.
// ---------------------------------------------------------------------------
static inline void fillEle(TH2D* hT, TH2D* hL, const EventData& ev, double w)
{
    const double pt_gev = ev.l2_pt * 1e-3;
    if (ev.el_is_tight) hT->Fill(pt_gev, ev.el_aeta, w);
    if (ev.el_is_loose) hL->Fill(pt_gev, ev.el_aeta, w);
}

// ---------------------------------------------------------------------------
// Compute ε_f = (data_T − mc_T) / (data_L − mc_L) per (pT, |η|) bin.
// Returns a TH2D with content + symmetric error.
// ---------------------------------------------------------------------------
static TH2D* computeEpsF(TH2D* dT, TH2D* dL, TH2D* mT, TH2D* mL,
                         const char* name)
{
    auto* h = (TH2D*)dT->Clone(name);
    h->Reset();
    for (int ix = 1; ix <= dT->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= dT->GetNbinsY(); ++iy) {
            double nT = dT->GetBinContent(ix, iy) - mT->GetBinContent(ix, iy);
            double nL = dL->GetBinContent(ix, iy) - mL->GetBinContent(ix, iy);
            double eT = std::sqrt(dT->GetBinError(ix,iy) * dT->GetBinError(ix,iy)
                                + mT->GetBinError(ix,iy) * mT->GetBinError(ix,iy));
            double eL = std::sqrt(dL->GetBinError(ix,iy) * dL->GetBinError(ix,iy)
                                + mL->GetBinError(ix,iy) * mL->GetBinError(ix,iy));
            if (nL <= 0 || nT < 0) {
                h->SetBinContent(ix, iy, 0.);
                h->SetBinError  (ix, iy, 0.);
                continue;
            }
            double r = nT / nL;
            double rel2 = (nT > 0 ? eT*eT/(nT*nT) : 0.) + eL*eL/(nL*nL);
            h->SetBinContent(ix, iy, r);
            h->SetBinError  (ix, iy, r * std::sqrt(rel2));
        }
    }
    return h;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void multijet_test()
{
    gStyle->SetOptStat(0);
    gSystem->mkdir("outputs", true);

    // ----- Build chains -----
    TChain dataCh("tree");
    dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data15.root");
    dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data16.root");
    dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_d/data17.root");
    dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/data18.root");

    TChain mcCh("tree");
    const std::string base = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    for (auto p : {"ttbar","ztautau","wjets","singletop","diboson","higgs"})
        mcCh.Add((base + p + "/*.root").c_str());

    // ----- 2D histograms: tight & loose, data & MC, standard & anti-iso -----
    auto book = [](const char* name) {
        return new TH2D(name, ";p_{T}^{e} [GeV];|#eta^{e}|",
                        kNPtBins, kPtBins, kNEtaBins, kEtaBins);
    };
    TH2D* dT_std = book("dT_std");
    TH2D* dL_std = book("dL_std");
    TH2D* mT_std = book("mT_std");
    TH2D* mL_std = book("mL_std");
    TH2D* dT_aiso = book("dT_aiso");
    TH2D* dL_aiso = book("dL_aiso");
    TH2D* mT_aiso = book("mT_aiso");
    TH2D* mL_aiso = book("mL_aiso");
    for (auto* h : {dT_std,dL_std,mT_std,mL_std,dT_aiso,dL_aiso,mT_aiso,mL_aiso})
        h->Sumw2();

    // ----- Loop helper: works on either chain -----
    auto loop = [&](TTree* t, bool isMC, const char* tag) {
        std::cout << "\n=== Loop " << tag << " (entries " << t->GetEntries() << ") ===\n";
        EventData ev;
        SetupBranches(t, ev);
        const Long64_t n = t->GetEntries();
        Long64_t nStd = 0, nAiso = 0;
        for (Long64_t i = 0; i < n; ++i) {
            if (i % 1000000 == 0) std::cout << "  " << i << " / " << n << "\n";
            t->GetEntry(i);
            ComputeDerived(ev);

            if (!passPreselNoMuTight(ev)) continue;

            // SS fake-CR-like cuts (no VBF / b-veto here — measurement region)
            if (!isFakeCR(ev)) continue;
            if (!(ev.jet1_pt > 0. && ev.jet2_pt > 0.)) continue;
            if (!ev.no_bjets) continue;

            // Truth-prompt cut on MC for subtraction (same as fake-rate measurement)
            if (isMC && !isMCPromptElectron(ev)) continue;

            const double w = isMC ? mcWeightLoose(ev) : 1.0;

            const bool muTight  = ev.mu_is_tight;
            const bool muLoose  = ev.mu_is_loose;
            const bool muAiso   = muLoose && !muTight;   // loose but failing tight = anti-iso

            if (muTight) {
                fillEle(isMC ? mT_std : dT_std,
                        isMC ? mL_std : dL_std, ev, w);
                ++nStd;
            } else if (muAiso) {
                fillEle(isMC ? mT_aiso : dT_aiso,
                        isMC ? mL_aiso : dL_aiso, ev, w);
                ++nAiso;
            }
        }
        std::cout << "  filled  std-mu: " << nStd
                  << "  anti-iso-mu: " << nAiso << "\n";
    };

    loop(&mcCh,   true,  "MC");
    loop(&dataCh, false, "DATA");

    // ----- Compute ε_f maps -----
    TH2D* eps_std  = computeEpsF(dT_std,  dL_std,  mT_std,  mL_std,  "eps_std");
    TH2D* eps_aiso = computeEpsF(dT_aiso, dL_aiso, mT_aiso, mL_aiso, "eps_aiso");

    // ----- Per-bin comparison table -----
    printf("\n========================================================================\n");
    printf("  Multijet-enrichment test: ε_f in standard CR vs anti-iso-μ CR\n");
    printf("========================================================================\n");
    printf("  pT [GeV]      |η|         ε_f (std)        ε_f (anti-iso)   ratio\n");
    printf("  ----------------------------------------------------------------\n");
    double maxRatio = 0., sumRatio = 0.; int nRatio = 0;
    for (int ix = 1; ix <= eps_std->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= eps_std->GetNbinsY(); ++iy) {
            double e1 = eps_std ->GetBinContent(ix, iy);
            double s1 = eps_std ->GetBinError  (ix, iy);
            double e2 = eps_aiso->GetBinContent(ix, iy);
            double s2 = eps_aiso->GetBinError  (ix, iy);
            double r  = (e1 > 0. && e2 > 0.) ? e2 / e1 : 0.;
            double ptLo = eps_std->GetXaxis()->GetBinLowEdge(ix);
            double ptHi = eps_std->GetXaxis()->GetBinUpEdge(ix);
            double etLo = eps_std->GetYaxis()->GetBinLowEdge(iy);
            double etHi = eps_std->GetYaxis()->GetBinUpEdge(iy);
            printf("  [%4.0f,%4.0f]  [%.2f,%.2f]   %.4f±%.4f   %.4f±%.4f   %.2f\n",
                   ptLo, ptHi, etLo, etHi, e1, s1, e2, s2, r);
            if (r > 0) { sumRatio += r; ++nRatio; if (r > maxRatio) maxRatio = r; }
        }
    }
    printf("  ----------------------------------------------------------------\n");
    if (nRatio > 0) {
        printf("  Average ratio (anti-iso / std): %.2f   max bin: %.2f\n",
               sumRatio/nRatio, maxRatio);
    }
    printf("\n  Interpretation:\n");
    printf("    ratio ~ 1: no detectable multijet contamination in std CR.\n");
    printf("    ratio > 1.3: anti-iso CR ε_f much larger → QCD multijet enriches\n");
    printf("                 the loose-only population there. Likely some of that\n");
    printf("                 contamination is also in the std CR, inflating its ε_f.\n");
    printf("    ratio < 0.7: std CR has a population (real prompt residual?) that the\n");
    printf("                 anti-iso CR removes — investigate.\n");

    // ----- Plot per-η slice: ε_f vs pT, std vs anti-iso -----
    TCanvas c("c","Multijet test", 1500, 900);
    c.Divide(3, 2);
    for (int iy = 1; iy <= kNEtaBins; ++iy) {
        c.cd(iy);
        TH1D* h1 = eps_std ->ProjectionX(Form("std_etabin%d", iy),  iy, iy);
        TH1D* h2 = eps_aiso->ProjectionX(Form("aiso_etabin%d", iy), iy, iy);
        // ProjectionX SUMS the y-bin, but for ε_f map (already a ratio per bin) we
        // want the actual per-bin values — use SetBinContent to copy directly
        for (int ix = 1; ix <= eps_std->GetNbinsX(); ++ix) {
            h1->SetBinContent(ix, eps_std ->GetBinContent(ix, iy));
            h1->SetBinError  (ix, eps_std ->GetBinError  (ix, iy));
            h2->SetBinContent(ix, eps_aiso->GetBinContent(ix, iy));
            h2->SetBinError  (ix, eps_aiso->GetBinError  (ix, iy));
        }
        h1->SetTitle(Form("|#eta| in [%.2f, %.2f];p_{T} [GeV];#varepsilon_{f}",
                          eps_std->GetYaxis()->GetBinLowEdge(iy),
                          eps_std->GetYaxis()->GetBinUpEdge(iy)));
        h1->SetMarkerStyle(20); h1->SetMarkerColor(kBlack); h1->SetLineColor(kBlack);
        h2->SetMarkerStyle(21); h2->SetMarkerColor(kRed);   h2->SetLineColor(kRed);
        h1->GetYaxis()->SetRangeUser(0., std::max(0.5, h1->GetMaximum()*1.5));
        h1->Draw("E1");
        h2->Draw("E1 SAME");
        TLegend* l = new TLegend(0.55, 0.75, 0.88, 0.88);
        l->AddEntry(h1, "std CR (tight #mu)",     "lp");
        l->AddEntry(h2, "anti-iso CR (loose #mu)", "lp");
        l->SetBorderSize(0);
        l->Draw();
    }
    c.SaveAs("outputs/multijet_test.pdf");

    // ----- Save -----
    TFile* fo = TFile::Open("outputs/multijet_test.root", "RECREATE");
    eps_std ->Write(); eps_aiso->Write();
    dT_std ->Write(); dL_std ->Write(); mT_std ->Write(); mL_std ->Write();
    dT_aiso->Write(); dL_aiso->Write(); mT_aiso->Write(); mL_aiso->Write();
    fo->Close();

    printf("\nWrote outputs/multijet_test.root\n");
    printf("Wrote outputs/multijet_test.pdf  (6-panel: ε_f vs pT per |η| slice)\n");
}
