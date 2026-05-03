// =============================================================================
// FILE:    studies/met_systematic_investigation.C
// PURPOSE: Test whether the fake-rate's MET-dependence is statistically
//          significant. Splits the SS fake CR at MET = 40 GeV and at 50 GeV
//          (in addition to the nominal 30 GeV split), measures the per-pT-bin
//          fake rate (tight/loose, after prompt-MC subtraction), and reports
//          the per-bin compatibility |Δ|/σ_stat.
//
// DECISION RULE: if max |Δ|/σ_stat < 1 across pT bins, the MET-dependence is
//                statistically compatible with zero — the MET systematic can
//                be dropped or absorbed into the (larger) composition or
//                prompt-subtraction uncertainty.
//
// INPUTS:  /eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/  (prompt MC)
//          AMM/data_with_electron_fake_weights.root or raw data ntuples
// OUTPUTS: studies/outputs/met_systematic_investigation.root
//          stdout: per-bin table
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "../utils/SelectionUtils.h"

// ---------------------------------------------------------------------------
// One per-pT-bin counter: numerator (tight) and denominator (loose), with sumw2.
// ---------------------------------------------------------------------------
struct Counter {
    std::vector<double> num, num2, den, den2;
    explicit Counter(int n) : num(n,0.), num2(n,0.), den(n,0.), den2(n,0.) {}
    void addNum(int b, double w) { if(b>=0&&b<(int)num.size()){ num[b]+=w; num2[b]+=w*w; } }
    void addDen(int b, double w) { if(b>=0&&b<(int)den.size()){ den[b]+=w; den2[b]+=w*w; } }
};

// ---------------------------------------------------------------------------
// pT bin lookup using kPtBins from SelectionUtils.h (GeV)
// ---------------------------------------------------------------------------
static int ptBin(double pt_gev) {
    for (int i = 0; i < kNPtBins; ++i)
        if (pt_gev >= kPtBins[i] && pt_gev < kPtBins[i+1]) return i;
    return -1;
}

// ---------------------------------------------------------------------------
// Compute fake rate ε = (data_T − mc_T) / (data_L − mc_L) per bin.
// Statistical uncertainty: standard ratio propagation with sumw2.
// ---------------------------------------------------------------------------
struct Rate { double val=0., err=0.; };
static Rate divide(double nT, double n2T, double nL, double n2L) {
    Rate r;
    if (nL <= 0.) return r;
    r.val = nT / nL;
    double rel2 = (n2T > 0. ? n2T/(nT*nT) : 0.) + (n2L > 0. ? n2L/(nL*nL) : 0.);
    r.err = std::fabs(r.val) * std::sqrt(rel2);
    return r;
}

static void fillRateHist(TH1D* h, const Counter& dat, const Counter& mc) {
    for (int b = 0; b < kNPtBins; ++b) {
        double nT  = dat.num[b]  - mc.num[b];
        double n2T = dat.num2[b] + mc.num2[b];   // add MC stat in quadrature
        double nL  = dat.den[b]  - mc.den[b];
        double n2L = dat.den2[b] + mc.den2[b];
        Rate r = divide(nT, n2T, nL, n2L);
        h->SetBinContent(b+1, r.val);
        h->SetBinError  (b+1, r.err);
    }
}

// ---------------------------------------------------------------------------
// Helper: classify MET (in GeV) into (low,high) at threshold T_GeV
// returns 0 for low, 1 for high, -1 if MET unavailable (treat as low).
// ---------------------------------------------------------------------------
static int metSplit(double met_gev, double T) { return met_gev < T ? 0 : 1; }

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void met_systematic_investigation()
{
    // -- 1. Set up data + MC chains
    const std::string base = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    TChain mcCh("tree");
    for (auto p : {"ttbar","ztautau","wjets","singletop","diboson","higgs"})
        mcCh.Add((base + p + "/*.root").c_str());

    TChain dataCh("tree");
    dataCh.Add("../AMM/data_with_electron_fake_weights.root");   // run from studies/
    if (dataCh.GetEntries() == 0) {
        std::cerr << "WARNING: AMM data ntuple not found, falling back to raw data tree\n";
        dataCh.Reset();
        // Same data list used in efficiency/data_eff.C
        dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data15.root");
        dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data16.root");
        dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_d/data17.root");
        dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/data18.root");
    }

    // -- 2. Counters: nominal + 2 split points × {low,high} = 5 buckets
    Counter d_nom(kNPtBins), m_nom(kNPtBins);
    Counter d_lo40(kNPtBins), m_lo40(kNPtBins);
    Counter d_hi40(kNPtBins), m_hi40(kNPtBins);
    Counter d_lo50(kNPtBins), m_lo50(kNPtBins);
    Counter d_hi50(kNPtBins), m_hi50(kNPtBins);

    // -- 3. Loop helper (works on any chain that uses EventData)
    auto loop = [&](TTree* t, bool isMC) {
        EventData ev;
        SetupBranches(t, ev);
        Long64_t n = t->GetEntries();
        std::cout << "  entries: " << n << (isMC ? "  (MC)\n" : "  (data)\n");

        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            // No data15 dispatch any more — passTriggerMatching ORs the
            // global flag with the per-chain flags, so missing branches
            // don't matter.
            ev.isData   = !isMC;
            ev.isData15 = false;
            ComputeDerived(ev);

            if (!passBasePresel(ev))         continue;
            if (!isFakeCR(ev))               continue;          // SS fake CR
            if (!ev.no_bjets)                continue;          // CR is b-veto
            if (!passVBFJets(ev))            continue;
            if (isMC && !isMCPromptEvent(ev)) continue;          // prompt-only MC

            const double pt_gev  = ev.l2_pt * 1e-3;
            const double met_gev = ev.MET   * 1e-3;
            const int    pb      = ptBin(pt_gev);
            if (pb < 0) continue;

            const double w = isMC ? mcWeightLoose(ev) : 1.0;

            const bool tight = ev.el_is_tight;
            const bool loose = ev.el_is_loose;
            if (!loose) continue;

            auto fill = [&](Counter& dC, Counter& mC) {
                Counter& C = isMC ? mC : dC;
                if (tight) C.addNum(pb, w);
                C.addDen(pb, w);
            };
            fill(d_nom, m_nom);
            if (metSplit(met_gev, 40.) == 0) fill(d_lo40, m_lo40); else fill(d_hi40, m_hi40);
            if (metSplit(met_gev, 50.) == 0) fill(d_lo50, m_lo50); else fill(d_hi50, m_hi50);
        }
    };

    std::cout << "\n=== MC loop ===\n";   loop(&mcCh,   true);
    std::cout << "\n=== Data loop ===\n"; loop(&dataCh, false);

    // -- 4. Build TH1Ds (bin axis = pT)
    auto mkHist = [](const char* name) {
        TH1D* h = new TH1D(name, ";p_{T}^{e} [GeV]; fake rate", kNPtBins, kPtBins);
        h->Sumw2();
        return h;
    };
    TH1D* h_nom    = mkHist("h_rate_nominal");
    TH1D* h_lo40   = mkHist("h_rate_met_low_40");
    TH1D* h_hi40   = mkHist("h_rate_met_high_40");
    TH1D* h_lo50   = mkHist("h_rate_met_low_50");
    TH1D* h_hi50   = mkHist("h_rate_met_high_50");

    fillRateHist(h_nom , d_nom , m_nom );
    fillRateHist(h_lo40, d_lo40, m_lo40);
    fillRateHist(h_hi40, d_hi40, m_hi40);
    fillRateHist(h_lo50, d_lo50, m_lo50);
    fillRateHist(h_hi50, d_hi50, m_hi50);

    // -- 5. Per-bin compatibility table
    std::cout << "\n=========================================================================================================\n";
    std::cout << "  pT bin       nominal           MET<40            MET>=40           |Δ|/σ40   MET<50         MET>=50         |Δ|/σ50\n";
    std::cout << "-----------------------------------------------------------------------------------------------------------\n";

    double max_pull40 = 0., max_pull50 = 0.;
    for (int b = 0; b < kNPtBins; ++b) {
        const int B = b + 1;
        double v0 = h_nom->GetBinContent(B), e0 = h_nom->GetBinError(B);
        double v_lo40 = h_lo40->GetBinContent(B), e_lo40 = h_lo40->GetBinError(B);
        double v_hi40 = h_hi40->GetBinContent(B), e_hi40 = h_hi40->GetBinError(B);
        double v_lo50 = h_lo50->GetBinContent(B), e_lo50 = h_lo50->GetBinError(B);
        double v_hi50 = h_hi50->GetBinContent(B), e_hi50 = h_hi50->GetBinError(B);

        double d40 = v_hi40 - v_lo40;
        double s40 = std::sqrt(e_lo40*e_lo40 + e_hi40*e_hi40);
        double p40 = (s40 > 0.) ? std::fabs(d40)/s40 : 0.;
        if (p40 > max_pull40) max_pull40 = p40;

        double d50 = v_hi50 - v_lo50;
        double s50 = std::sqrt(e_lo50*e_lo50 + e_hi50*e_hi50);
        double p50 = (s50 > 0.) ? std::fabs(d50)/s50 : 0.;
        if (p50 > max_pull50) max_pull50 = p50;

        printf("  [%3.0f-%3.0f]  %.3f±%.3f   %.3f±%.3f   %.3f±%.3f   %5.2fσ    %.3f±%.3f   %.3f±%.3f   %5.2fσ\n",
               kPtBins[b], kPtBins[b+1],
               v0, e0,
               v_lo40, e_lo40, v_hi40, e_hi40, p40,
               v_lo50, e_lo50, v_hi50, e_hi50, p50);
    }
    std::cout << "-----------------------------------------------------------------------------------------------------------\n";
    printf("  Max |Δ|/σ over pT bins: 40-GeV split = %.2f,  50-GeV split = %.2f\n", max_pull40, max_pull50);
    if (std::max(max_pull40, max_pull50) < 1.0) {
        std::cout << "\n  ==> MET-dependence is statistically compatible with zero across pT bins.\n";
        std::cout << "      The MET systematic can be dropped or replaced by the larger composition /\n";
        std::cout << "      prompt-subtraction uncertainty.\n";
    } else {
        std::cout << "\n  ==> MET-dependence is statistically resolved in at least one bin.\n";
        std::cout << "      Keep the MET systematic; consider symmetrizing using the larger of the two splits.\n";
    }

    // -- 6. Save
    gSystem->mkdir("outputs", true);
    TFile* fout = TFile::Open("outputs/met_systematic_investigation.root", "RECREATE");
    h_nom->Write(); h_lo40->Write(); h_hi40->Write(); h_lo50->Write(); h_hi50->Write();
    fout->Close();
    std::cout << "\nWrote studies/outputs/met_systematic_investigation.root\n";
}
