// =============================================================================
// FILE:    studies/eps_f_vs_met_cut.C
// PURPOSE: Scan ε_f as a function of MET-cut value to find the threshold that
//          best removes multijet contamination without crashing statistics.
//
//          For each MET threshold in {0, 10, 20, ..., 100} GeV, computes
//          ε_f(pT, |η|) using only events with MET > threshold. Plots:
//            - ε_f vs MET-cut, one curve per pT bin (averaged over |η|)
//            - Integrated ε_f over the [50,200] high-pT bin (where multijet
//              is most contaminating) vs MET-cut → use this to pick the cut
//            - Surviving event count vs MET-cut → trade-off check
//
// The "best" MET cut is the smallest value at which ε_f stabilises (plateaus).
// Below that the multijet leak distorts ε_f; above that you only lose stat.
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

#include "../../utils/SelectionUtils.h"

// MET thresholds to scan (in GeV)
static const std::vector<double> kMetCuts = {0, 10, 20, 30, 40, 50, 60, 70, 80, 100};

static inline void fillEle(TH2D* hT, TH2D* hL, const EventData& ev, double w)
{
    const double pt_gev = ev.l2_pt * 1e-3;
    if (ev.el_is_tight) hT->Fill(pt_gev, ev.el_aeta, w);
    if (ev.el_is_loose) hL->Fill(pt_gev, ev.el_aeta, w);
}

static double divEpsF(double dT, double dL, double mT, double mL)
{
    double nT = dT - mT;
    double nL = dL - mL;
    if (nL <= 0 || nT < 0) return 0.;
    return std::min(1.0, nT / nL);
}

static double divEpsErr(double dT, double dL, double mT, double mL,
                        double edT, double edL, double emT, double emL)
{
    double nT = dT - mT;
    double nL = dL - mL;
    if (nL <= 0 || nT <= 0) return 0.;
    double eT = std::sqrt(edT*edT + emT*emT);
    double eL = std::sqrt(edL*edL + emL*emL);
    double r  = nT / nL;
    return r * std::sqrt(eT*eT/(nT*nT) + eL*eL/(nL*nL));
}

// ---------------------------------------------------------------------------
void eps_f_vs_met_cut()
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

    // For each MET cut, one set of (data_T, data_L, mc_T, mc_L) histograms.
    const int K = (int)kMetCuts.size();
    auto book = [](const char* nm) {
        return new TH2D(nm, ";p_{T}^{e} [GeV];|#eta^{e}|",
                        kNPtBins, kPtBins, kNEtaBins, kEtaBins);
    };
    std::vector<TH2D*> dT(K), dL(K), mT(K), mL(K);
    for (int k = 0; k < K; ++k) {
        dT[k] = book(Form("dT_%d", k));
        dL[k] = book(Form("dL_%d", k));
        mT[k] = book(Form("mT_%d", k));
        mL[k] = book(Form("mL_%d", k));
        for (auto* h : {dT[k], dL[k], mT[k], mL[k]}) h->Sumw2();
    }
    // Surviving event count per MET cut
    std::vector<long> nSurvData(K, 0), nSurvMC(K, 0);

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
            if (!(ev.jet1_pt > 0. && ev.jet2_pt > 0.)) continue;
            if (!ev.no_bjets) continue;

            if (isMC && !isMCPromptElectron(ev)) continue;
            const double w = isMC ? mcWeightLoose(ev) : 1.0;

            // Fill EVERY histogram whose MET threshold is satisfied
            for (int k = 0; k < K; ++k) {
                if (ev.MET < kMetCuts[k] * 1000.) continue;   // MET in MeV
                fillEle(isMC ? mT[k] : dT[k], isMC ? mL[k] : dL[k], ev, w);
                if (isMC) ++nSurvMC[k]; else ++nSurvData[k];
            }
        }
    };

    loop(&mcCh,   true,  "MC");
    loop(&dataCh, false, "DATA");

    // ----- Per pT bin: ε_f vs MET cut, averaged over |η| -----
    // Use simple unweighted average over η (5 bins)
    std::vector<std::vector<double>> epsByPt(kNPtBins, std::vector<double>(K, 0.));
    std::vector<std::vector<double>> errByPt(kNPtBins, std::vector<double>(K, 0.));
    for (int ix = 1; ix <= kNPtBins; ++ix) {
        for (int k = 0; k < K; ++k) {
            double sumE = 0., sumE2 = 0.; int n = 0;
            for (int iy = 1; iy <= kNEtaBins; ++iy) {
                double e = divEpsF(dT[k]->GetBinContent(ix,iy),
                                   dL[k]->GetBinContent(ix,iy),
                                   mT[k]->GetBinContent(ix,iy),
                                   mL[k]->GetBinContent(ix,iy));
                double s = divEpsErr(dT[k]->GetBinContent(ix,iy),
                                     dL[k]->GetBinContent(ix,iy),
                                     mT[k]->GetBinContent(ix,iy),
                                     mL[k]->GetBinContent(ix,iy),
                                     dT[k]->GetBinError  (ix,iy),
                                     dL[k]->GetBinError  (ix,iy),
                                     mT[k]->GetBinError  (ix,iy),
                                     mL[k]->GetBinError  (ix,iy));
                if (e > 0) { sumE += e; sumE2 += s*s; ++n; }
            }
            if (n > 0) {
                epsByPt[ix-1][k] = sumE / n;
                errByPt[ix-1][k] = std::sqrt(sumE2) / n;
            }
        }
    }

    // ----- Print scan table -----
    printf("\n========================================================================\n");
    printf("  ε_f vs MET cut  (averaged over |η|)\n");
    printf("========================================================================\n");
    printf("  MET cut [GeV]   ");
    for (int ix = 1; ix <= kNPtBins; ++ix) {
        printf("  pT[%3.0f-%3.0f]   ",
               kPtBins[ix-1], kPtBins[ix]);
    }
    printf("  N_data   N_MC\n");
    printf("  ----------------------------------------------------------------\n");
    for (int k = 0; k < K; ++k) {
        printf("       > %3.0f       ", kMetCuts[k]);
        for (int ix = 1; ix <= kNPtBins; ++ix) {
            printf(" %.4f±%.4f ", epsByPt[ix-1][k], errByPt[ix-1][k]);
        }
        printf("  %7ld  %7ld\n", nSurvData[k], nSurvMC[k]);
    }

    // ----- Plot -----
    TCanvas c("c","ε_f vs MET cut", 1200, 800);
    auto* leg = new TLegend(0.62, 0.55, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    int colors[] = {kBlue, kRed, kGreen+2, kMagenta, kOrange+1, kCyan+1, kBlack};
    std::vector<TGraphErrors*> g(kNPtBins);

    // First pass: build all graphs and find the global y-max so the y-range
    // accommodates EVERY pT-bin curve, not just the first one drawn.
    double yMax = 0.;
    for (int ix = 0; ix < kNPtBins; ++ix) {
        g[ix] = new TGraphErrors(K);
        for (int k = 0; k < K; ++k) {
            g[ix]->SetPoint     (k, kMetCuts[k], epsByPt[ix][k]);
            g[ix]->SetPointError(k, 0,           errByPt[ix][k]);
            yMax = std::max(yMax, epsByPt[ix][k] + errByPt[ix][k]);
        }
        g[ix]->SetTitle(";MET threshold [GeV];#LT#varepsilon_{f}#GT (averaged over |#eta|)");
        g[ix]->SetMarkerStyle(20);
        g[ix]->SetMarkerColor(colors[ix % 7]);
        g[ix]->SetLineColor  (colors[ix % 7]);
        g[ix]->SetLineWidth(2);
    }

    // Second pass: draw with consistent y-range covering all curves.
    for (int ix = 0; ix < kNPtBins; ++ix) {
        if (ix == 0) {
            g[ix]->Draw("APL");
            g[ix]->GetYaxis()->SetRangeUser(0., 1.25 * yMax);
        } else {
            g[ix]->Draw("PL SAME");
        }
        leg->AddEntry(g[ix], Form("p_{T} #in [%.0f, %.0f] GeV",
                                  kPtBins[ix], kPtBins[ix+1]), "lp");
    }
    leg->Draw();
    c.SaveAs("outputs/eps_f_vs_met_cut.pdf");

    // ----- Recommendation: where does the highest-pT bin plateau? -----
    int hiBin = kNPtBins - 1;   // last (highest-pT) bin
    printf("\n========================================================================\n");
    printf("  Plateau search: pT bin [%.0f, %.0f] GeV\n",
           kPtBins[hiBin], kPtBins[hiBin+1]);
    printf("========================================================================\n");
    int plateauK = -1;
    for (int k = 1; k < K; ++k) {
        double e0 = epsByPt[hiBin][k-1], e1 = epsByPt[hiBin][k];
        double s0 = errByPt[hiBin][k-1], s1 = errByPt[hiBin][k];
        double pull = (s0 + s1 > 0) ? std::fabs(e1 - e0) / std::sqrt(s0*s0 + s1*s1) : 0.;
        printf("  step %3.0f → %3.0f: ε_f %.3f → %.3f  (Δ pull = %.2fσ)\n",
               kMetCuts[k-1], kMetCuts[k], e0, e1, pull);
        if (plateauK < 0 && pull < 1.0 && k >= 2) plateauK = k - 1;
    }
    if (plateauK >= 0) {
        printf("\n  ε_f stabilises starting at MET > %.0f GeV "
               "(further cuts only cost statistics).\n", kMetCuts[plateauK]);
        printf("  RECOMMENDED CUT: MET > %.0f GeV\n", kMetCuts[plateauK]);
    } else {
        printf("\n  ε_f never plateaus — keeps moving with MET cut.\n");
        printf("  Consider widening the scan to higher MET, or pick a cut where the\n");
        printf("  trade-off (statistics vs ε_f stability) is acceptable.\n");
    }

    TFile* fo = TFile::Open("outputs/eps_f_vs_met_cut.root", "RECREATE");
    for (int k = 0; k < K; ++k) {
        dT[k]->Write(); dL[k]->Write(); mT[k]->Write(); mL[k]->Write();
    }
    for (int ix = 0; ix < kNPtBins; ++ix)
        g[ix]->Write(Form("g_eps_f_pt_bin_%d", ix));
    fo->Close();
    printf("\nWrote outputs/eps_f_vs_met_cut.root\nWrote outputs/eps_f_vs_met_cut.pdf\n");
}
