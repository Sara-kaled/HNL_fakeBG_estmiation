// =============================================================================
// FILE:       studies/mc_closure_fake.C
// PURPOSE:    Single-pass MC closure for the AMM fake estimate, plus five
//             diagnostics to localize any closure failure:
//
//   D1. Effective ε_f from THIS MC chain (data-driven recipe applied to MC):
//       counts fake-tight vs fake-loose in the SS CR and prints the per-pT-bin
//       ε_f. Compare to the stored fake_rate_nominal map — if they differ,
//       the stored map is biased relative to what the truth says it should be.
//
//   D2. mcWeightLoose / mcWeightTight ratio histogram. The closure compares
//       N_truth (weighted by mcWeightTight) against N_MM (weighted by
//       mcWeightLoose × amm_weight). If the two SF combinations differ
//       systematically, the closure is offset by their average ratio.
//
//   D3. MM positive / negative breakdown per region. The MM total is the sum
//       of negative contributions from tight events (w_T) and positive
//       contributions from loose-only events (w_L). If one side dominates
//       unphysically, you've got a measurement bias.
//
//   D4. Truth fake counts split by tight vs loose-only. ε_f is supposed to
//       describe the rate at which loose-only fakes get promoted to tight.
//       Printing the two counts side-by-side tells you the "ground truth" ε_f.
//
//   D5. Stored vs effective ε_f comparison table per pT bin.
//
// INPUTS:     MC ntuples (TChain)
//             ../efficiency/outputs/fake_eff_with_systematics.root
//             ../efficiency/outputs/real_eff_histograms.root
// OUTPUTS:    outputs/mc_closure_fake_SR_shapes.root
// RUNS AFTER: build_fake_systematics.C, real_eff.C
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TSystem.h>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <string>

using std::cout;
using std::endl;

#include "../utils/SelectionUtils.h"

enum Region { kPre=0, kCR=1, kVR=2, kSR=3, kNRegion=4 };
static const char* regionNames[kNRegion] = { "Preselection", "CR", "VR", "SR" };

// ---------------------------------------------------------------------------
static double amm_weight(bool isTight, double real_eff, double fake_eff)
{
    real_eff = std::max(1e-5, std::min(1.0 - 1e-5, real_eff));
    fake_eff = std::max(1e-5, std::min(real_eff - 1e-5, fake_eff));
    const double denom = real_eff - fake_eff;
    if (std::abs(denom) < 1e-9) return 0.0;
    return isTight ? (fake_eff / denom) * (real_eff - 1.0)
                   : (fake_eff / denom) *  real_eff;
}

static double get_eff(TH2F* h, double pt_gev, double eta)
{
    if (!h) return 0.0;
    int bx = std::max(1, std::min(h->GetXaxis()->FindBin(pt_gev), h->GetNbinsX()));
    int by = std::max(1, std::min(h->GetYaxis()->FindBin(std::abs(eta)), h->GetNbinsY()));
    return h->GetBinContent(bx, by);
}

// ---------------------------------------------------------------------------
void mc_closure_fake()
{
    // -----------------------------------------------------------------------
    // 0. Open efficiency maps
    // -----------------------------------------------------------------------
    TFile* fFake = TFile::Open("../efficiency/outputs/fake_eff_with_systematics.root", "READ");
    TFile* fReal = TFile::Open("../efficiency/outputs/real_eff_histograms.root", "READ");
    if (!fFake || fFake->IsZombie() || !fReal || fReal->IsZombie()) {
        std::cerr << "ERROR: cannot open efficiency input files.\n"; return;
    }
    TH2F* hFakeNom   = (TH2F*)fFake->Get("fake_rate_nominal");
    TH2F* hRealTight = (TH2F*)fReal->Get("h2_tight_eta_pt_CR");
    TH2F* hRealLoose = (TH2F*)fReal->Get("h2_loose_eta_pt_CR");
    if (!hFakeNom || !hRealTight || !hRealLoose) {
        std::cerr << "ERROR: efficiency histograms missing\n"; return;
    }
    TH2F* hRealEff = (TH2F*)hRealTight->Clone("real_eff_el");
    hRealEff->Divide(hRealLoose);

    // -----------------------------------------------------------------------
    // 1. MC chain
    // -----------------------------------------------------------------------
    const std::string base = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    TChain mc("tree");
    for (auto p : {"ttbar","ztautau","wjets","singletop","diboson","higgs"})
        mc.Add((base + p + "/*.root").c_str());

    // -----------------------------------------------------------------------
    // 2. Yield accumulators + diagnostic accumulators
    // -----------------------------------------------------------------------
    double Ntruth[kNRegion]={}, Ntruth2[kNRegion]={};
    double Nmm   [kNRegion]={}, Nmm2   [kNRegion]={};

    // D3 — MM split by sign of weight (tight vs loose-only contributions)
    double NmmPos[kNRegion]={}, NmmNeg[kNRegion]={};
    long   nMmPos[kNRegion]={}, nMmNeg[kNRegion]={};

    // D4 — Truth fakes split by tight / loose-only in CR
    double Tfake_tight_CR  = 0., Tfake_loose_CR  = 0.;
    double Tprompt_tight_CR= 0., Tprompt_loose_CR= 0.;

    // D1 — Effective ε_f per pT bin from THIS MC, mirroring the data recipe
    //      (loose-only denominator includes ALL fakes; tight numerator =
    //       fake & tight). No prompt-MC subtraction here — pure truth.
    const int nPtDiag = hFakeNom->GetNbinsX();
    std::vector<double> num_diag(nPtDiag, 0.), den_diag(nPtDiag, 0.);

    // D2 — weight ratio collector
    long   nWeightSamples = 0;
    double sumW_tight = 0., sumW_loose = 0.;
    TH1D* hWRatio = new TH1D("hWRatio",
        "mcWeightLoose / mcWeightTight; ratio; events", 200, 0., 2.);

    TH1D* hTruth_SR = new TH1D("hTruth_SR","Truth fake; m_{jj} [GeV]; Events", 10, 0., 2000.);
    TH1D* hMM_SR    = new TH1D("hMM_SR",   "MM fake;    m_{jj} [GeV]; Events", 10, 0., 2000.);
    hTruth_SR->Sumw2(); hMM_SR->Sumw2();

    // -----------------------------------------------------------------------
    // 3. Single-pass loop
    // -----------------------------------------------------------------------
    cout << "\n=== Single-pass closure loop ===\n";
    EventData ev;
    SetupBranches(&mc, ev);
    const Long64_t n = mc.GetEntries();
    cout << "  Entries: " << n << "\n";

    for (Long64_t i = 0; i < n; ++i) {
        if (i % 1000000 == 0) cout << "  " << i << " / " << n << "\n";
        mc.GetEntry(i);
        ComputeDerived(ev);

        if (!passBasePresel(ev)) continue;

        const double w_tight = mcWeightTight(ev);
        const double w_loose = mcWeightLoose(ev);
        const bool   is_fake = isMCFakeElectron(ev);

        // D2: collect weight ratio for ALL events (any region)
        if (w_tight > 0. && std::isfinite(w_loose / w_tight)) {
            ++nWeightSamples;
            sumW_tight += w_tight;
            sumW_loose += w_loose;
            hWRatio->Fill(w_loose / w_tight);
        }

        // ---- Truth side ----
        // The AMM predicts the TIGHT fake yield in each region.
        // Therefore truth must count tight fakes only (is_fake AND el_is_tight),
        // otherwise we're comparing apples to oranges.
        if (is_fake && ev.el_is_tight) {
            Ntruth[kPre]  += w_tight;
            Ntruth2[kPre] += w_tight*w_tight;
        }
        if (isFakeCR(ev) && passVBFJets(ev) && ev.no_bjets) {
            // D4: split CR truth fakes by tight / loose-only
            if (ev.el_is_tight) {
                if (is_fake) Tfake_tight_CR  += w_tight;
                else         Tprompt_tight_CR+= w_tight;
            } else if (ev.el_is_loose) {
                if (is_fake) Tfake_loose_CR  += w_loose;
                else         Tprompt_loose_CR+= w_loose;
            }
            // D1: effective ε_f numerator/denominator per pT bin
            if (ev.el_is_loose && is_fake) {
                int b = std::max(1, std::min(hFakeNom->GetXaxis()->FindBin(ev.l2_pt*1e-3),
                                              hFakeNom->GetNbinsX())) - 1;
                den_diag[b] += w_loose;
                if (ev.el_is_tight) num_diag[b] += w_tight;
            }
            if (is_fake && ev.el_is_tight) {
                Ntruth[kCR]  += w_tight;
                Ntruth2[kCR] += w_tight*w_tight;
            }
        }
        if (ev.el_is_tight) {
            if (isSR(ev) && is_fake) {
                Ntruth[kSR]  += w_tight;
                Ntruth2[kSR] += w_tight*w_tight;
                hTruth_SR->Fill(ev.mjj * 1e-3, w_tight);
            }
            if (isVR(ev) && is_fake) {
                Ntruth[kVR]  += w_tight;
                Ntruth2[kVR] += w_tight*w_tight;
            }
        }

        // ---- MM side ----
        if (!ev.el_is_loose) continue;
        const double pt_gev = ev.l2_pt * 1e-3;
        const double eta    = ev.l2_eta;
        const double er     = get_eff(hRealEff, pt_gev, eta);
        const double ef     = get_eff(hFakeNom, pt_gev, eta);
        const double w_amm  = amm_weight(ev.el_is_tight, er, ef);
        const double wf     = w_amm * w_loose;

        auto bookMM = [&](int r) {
            Nmm [r]  += wf;
            Nmm2[r]  += wf*wf;
            if (wf >= 0) { NmmPos[r] += wf; ++nMmPos[r]; }
            else         { NmmNeg[r] += wf; ++nMmNeg[r]; }
        };
        bookMM(kPre);
        if (isFakeCR(ev) && passVBFJets(ev) && ev.no_bjets) bookMM(kCR);
        if (isSR(ev)) { bookMM(kSR); hMM_SR->Fill(ev.mjj * 1e-3, wf); }
        if (isVR(ev))   bookMM(kVR);
    }

    // -----------------------------------------------------------------------
    // 4. Closure summary
    // -----------------------------------------------------------------------
    cout << "\n========================================================\n";
    cout << "  MC closure: truth fake vs MM fake\n";
    cout << "========================================================\n";
    for (int r = 0; r < kNRegion; ++r) {
        double Nt = Ntruth[r], eNt = std::sqrt(Ntruth2[r]);
        double Nm = Nmm[r],    eNm = std::sqrt(Nmm2[r]);
        double diff = Nm - Nt;
        double rel  = (Nt != 0.) ? diff/Nt : 0.;
        cout << "\nRegion: " << regionNames[r] << "\n";
        cout << "  N_fake_truth = " << Nt << " ± " << eNt << "\n";
        cout << "  N_fake_MM    = " << Nm << " ± " << eNm << "\n";
        cout << "  Diff (MM-truth) = " << diff;
        if (Nt != 0.) cout << "  (" << 100.*rel << "% of truth)";
        cout << "\n";
    }

    // -----------------------------------------------------------------------
    // D2 — weight-SF ratio diagnostic
    // -----------------------------------------------------------------------
    cout << "\n========================================================\n";
    cout << "  D2 — mcWeightLoose / mcWeightTight\n";
    cout << "========================================================\n";
    cout << "  events sampled:           " << nWeightSamples << "\n";
    cout << "  sum mcWeightTight:        " << sumW_tight << "\n";
    cout << "  sum mcWeightLoose:        " << sumW_loose << "\n";
    if (sumW_tight > 0.) {
        cout << "  global ratio loose/tight: " << sumW_loose/sumW_tight << "\n";
        cout << "  (if << 1 → MM is suppressed by missing/different lepton-SF;"
                "  loose SF should be used on BOTH sides for a fair closure)\n";
    }
    cout << "  ratio mean:               " << hWRatio->GetMean()
         << "  RMS: " << hWRatio->GetRMS() << "\n";

    // -----------------------------------------------------------------------
    // D3 — MM positive vs negative breakdown
    // -----------------------------------------------------------------------
    cout << "\n========================================================\n";
    cout << "  D3 — MM positive (loose-only) vs negative (tight) yield per region\n";
    cout << "========================================================\n";
    cout << std::left << std::setw(15) << "region"
         << std::right << std::setw(14) << "ΣwMM(+)"
         << std::setw(10) << "n(+)"
         << std::setw(14) << "ΣwMM(-)"
         << std::setw(10) << "n(-)"
         << std::setw(14) << "net" << "\n";
    for (int r = 0; r < kNRegion; ++r) {
        cout << std::left << std::setw(15) << regionNames[r]
             << std::right << std::setw(14) << NmmPos[r]
             << std::setw(10) << nMmPos[r]
             << std::setw(14) << NmmNeg[r]
             << std::setw(10) << nMmNeg[r]
             << std::setw(14) << (NmmPos[r] + NmmNeg[r]) << "\n";
    }

    // -----------------------------------------------------------------------
    // D4 — CR truth fake split (ground truth ε_f from this MC, integrated)
    // -----------------------------------------------------------------------
    cout << "\n========================================================\n";
    cout << "  D4 — Truth tight vs loose-only in SS CR\n";
    cout << "========================================================\n";
    cout << "  truth fakes:    tight = " << Tfake_tight_CR
         << ",  loose-only = " << (Tfake_loose_CR - Tfake_tight_CR) << "\n";
    cout << "  truth prompts:  tight = " << Tprompt_tight_CR
         << ",  loose-only = " << (Tprompt_loose_CR - Tprompt_tight_CR) << "\n";
    if (Tfake_loose_CR > 0.)
        cout << "  effective ε_f (truth, integrated) = "
             << Tfake_tight_CR / Tfake_loose_CR << "\n";

    // -----------------------------------------------------------------------
    // D1 + D5 — per-pT-bin effective ε_f vs stored map
    // -----------------------------------------------------------------------
    cout << "\n========================================================\n";
    cout << "  D1/D5 — Per-pT effective ε_f from THIS MC vs stored map\n";
    cout << "========================================================\n";
    cout << std::left << std::setw(20) << "pT [GeV]"
         << std::setw(15) << "ε_f (truth)"
         << std::setw(15) << "ε_f (stored)"
         << std::setw(12) << "ratio" << "\n";
    for (int b = 0; b < nPtDiag; ++b) {
        double lo = hFakeNom->GetXaxis()->GetBinLowEdge(b+1);
        double hi = hFakeNom->GetXaxis()->GetBinUpEdge(b+1);
        double effTruth  = (den_diag[b] > 0.) ? num_diag[b] / den_diag[b] : 0.;
        // Stored ε_f averaged over |η| (just bin 1 for the diagnostic)
        double effStored = hFakeNom->GetBinContent(b+1, 1);
        cout << std::left << "  ["
             << std::setw(5) << lo << "," << std::setw(8) << hi << "]  "
             << std::right << std::setw(12) << effTruth
             << std::setw(15) << effStored
             << std::setw(12) << (effStored>0 ? effTruth/effStored : 0.) << "\n";
    }
    cout << "\n  Interpretation (ratio = truth ε_f / stored ε_f):\n"
         << "    ratio ~ 1   : ε_f matches what truth says → closure should hold.\n"
         << "    ratio < 1   : stored ε_f is BIGGER than truth → MM over-predicts.\n"
         << "                  (data has more tight fakes per loose than MC predicts;\n"
         << "                   could be data-MC fake-mix difference or under-subtracted prompt MC).\n"
         << "    ratio > 1   : stored ε_f is SMALLER than truth → MM under-predicts.\n";

    // -----------------------------------------------------------------------
    // Save shape histograms + diagnostics
    // -----------------------------------------------------------------------
    gSystem->mkdir("outputs", true);
    TFile* fout = new TFile("outputs/mc_closure_fake_SR_shapes.root","RECREATE");
    hTruth_SR->Write();
    hMM_SR->Write();
    hWRatio->Write();
    fout->Close();
    cout << "\nWrote outputs/mc_closure_fake_SR_shapes.root\n";

    fFake->Close();
    fReal->Close();
}
