// =============================================================================
// FILE:    AMM/AMM_compare_yaml.C
// PURPOSE: Independent AMM weight calculation, applying the SAME selection as
//          the YAML EventSelection 'only_emu' block:
//             JET_N baselineJvtLoose 20000 >= 2
//             EL_N 4500 == 1
//             MU_N 15000 == 1
//
//          Reads the flat ntuple that already contains:
//             - el_pt_NOSYS, el_eta, el_select_loose_NOSYS, el_select_tight_NOSYS
//             - mu_pt_NOSYS, mu_select_loose_NOSYS, mu_select_tight_NOSYS
//             - jet_pt_NOSYS, jet_select_baselineJvtLoose_NOSYS
//             - fake_weight_nominal  (AMM weight from AMM_ElectronOnly.C)
//
//          Re-computes the per-event AMM weight from fake_eff_input.root
//          and writes a small TTree with one row per passing event:
//             eventNumber, runNumber
//             w_amm_macro       (recomputed here)
//             w_amm_framework   (= fake_weight_nominal from the input)
//             w_diff            = macro − framework
//          Plus 1D histograms of both weights and their per-event difference.
//
// USAGE:
//   root -b -q -l 'AMM_compare_yaml.C("input.root","fake_eff_input.root","out.root")'
//
// OUTPUT:
//   <out.root>:
//     - tree comparison
//     - h_w_macro, h_w_framework, h_w_diff
//     - h_w_macro_vs_framework  (TH2D)
// =============================================================================
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TBranch.h>
#include "ROOT/RVec.hxx"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

using ROOT::VecOps::RVec;

// ---------------------------------------------------------------------------
// AMM per-electron weight (matches AMM_ElectronOnly.C clamps + formulas)
// ---------------------------------------------------------------------------
static double amm_weight(bool isTight, double er, double ef)
{
    er = std::max(1e-5, std::min(1.0 - 1e-5, er));
    ef = std::max(1e-5, std::min(er - 1e-5, ef));
    const double denom = er - ef;
    if (std::abs(denom) < 1e-9) return 0.0;
    return isTight ? (ef / denom) * (er - 1.0)   // negative
                   : (ef / denom) *  er;          // positive
}

// ---------------------------------------------------------------------------
// 2D efficiency lookup with last-bin clamping
// ---------------------------------------------------------------------------
static double get_eff(TH2D* h, double pt_gev, double eta)
{
    if (!h) return 0.0;
    int bx = std::max(1, std::min(h->GetXaxis()->FindBin(pt_gev), h->GetNbinsX()));
    int by = std::max(1, std::min(h->GetYaxis()->FindBin(std::abs(eta)), h->GetNbinsY()));
    return h->GetBinContent(bx, by);
}

// ---------------------------------------------------------------------------
// Count entries above pT threshold in an RVec (in MeV)
// ---------------------------------------------------------------------------
template <class V>
static int count_above(const V* v, double thresh_mev) {
    if (!v) return 0;
    int n = 0;
    for (size_t i = 0; i < v->size(); ++i) if ((*v)[i] > thresh_mev) ++n;
    return n;
}

// ---------------------------------------------------------------------------
// Count jets passing baselineJvtLoose AND pT threshold
// ---------------------------------------------------------------------------
static int count_good_jets(const RVec<float>* pt, const RVec<char>* jvt,
                           double pt_thresh_mev) {
    if (!pt || !jvt) return 0;
    int n = 0;
    const size_t N = std::min(pt->size(), jvt->size());
    for (size_t i = 0; i < N; ++i)
        if ((*pt)[i] > pt_thresh_mev && (*jvt)[i]) ++n;
    return n;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void AMM_compare_yaml(const char* inFile  = "data_with_electron_fake_weights.root",
                      const char* effFile = "fake_eff_input.root",
                      const char* outFile = "amm_comparison.root")
{
    // ----- 1. Open efficiency file -----
    TFile* fe = TFile::Open(effFile, "READ");
    if (!fe || fe->IsZombie()) {
        std::cerr << "ERROR: cannot open " << effFile << "\n"; return;
    }
    TH2D* hRealE = (TH2D*)fe->Get("RealEfficiency2D_el_pt_eta");
    TH2D* hFakeE = (TH2D*)fe->Get("FakeEfficiency2D_el_pt_eta");
    if (!hRealE || !hFakeE) {
        std::cerr << "ERROR: efficiency histograms not found in " << effFile << "\n";
        std::cerr << "       Expected: RealEfficiency2D_el_pt_eta, FakeEfficiency2D_el_pt_eta\n";
        return;
    }

    // ----- 2. Open input ntuple -----
    TFile* fi = TFile::Open(inFile, "READ");
    if (!fi || fi->IsZombie()) {
        std::cerr << "ERROR: cannot open " << inFile << "\n"; return;
    }
    TTree* t = (TTree*)fi->Get("tree");
    if (!t) {
        std::cerr << "ERROR: 'tree' not found in " << inFile << "\n"; return;
    }

    // ----- 3. Set up branches -----
    UInt_t   runNumber   = 0;
    ULong64_t eventNumber = 0;
    Double_t fw_framework = 0.;

    RVec<float>* el_pt           = nullptr;
    RVec<float>* el_eta          = nullptr;
    RVec<char>*  el_loose        = nullptr;
    RVec<char>*  el_tight        = nullptr;

    RVec<float>* mu_pt           = nullptr;
    RVec<char>*  mu_loose        = nullptr;
    RVec<char>*  mu_tight        = nullptr;

    RVec<float>* jet_pt          = nullptr;
    RVec<char>*  jet_jvt_loose   = nullptr;

    t->SetBranchAddress("runNumber",                          &runNumber);
    t->SetBranchAddress("eventNumber",                        &eventNumber);
    t->SetBranchAddress("fake_weight_nominal",                &fw_framework);

    t->SetBranchAddress("el_pt_NOSYS",                        &el_pt);
    t->SetBranchAddress("el_eta",                             &el_eta);
    t->SetBranchAddress("el_select_loose_NOSYS",              &el_loose);
    t->SetBranchAddress("el_select_tight_NOSYS",              &el_tight);

    t->SetBranchAddress("mu_pt_NOSYS",                        &mu_pt);
    t->SetBranchAddress("mu_select_loose_NOSYS",              &mu_loose);
    t->SetBranchAddress("mu_select_tight_NOSYS",              &mu_tight);

    t->SetBranchAddress("jet_pt_NOSYS",                       &jet_pt);
    t->SetBranchAddress("jet_select_baselineJvtLoose_NOSYS",  &jet_jvt_loose);

    // ----- 4. Output tree + histograms -----
    TFile* fo = TFile::Open(outFile, "RECREATE");
    TTree* o  = new TTree("comparison", "AMM weight comparison");
    Double_t w_macro = 0., w_diff = 0.;
    o->Branch("runNumber",       &runNumber,       "runNumber/i");
    o->Branch("eventNumber",     &eventNumber,     "eventNumber/l");
    o->Branch("w_amm_macro",     &w_macro,         "w_amm_macro/D");
    o->Branch("w_amm_framework", &fw_framework,    "w_amm_framework/D");
    o->Branch("w_diff",          &w_diff,          "w_diff/D");

    TH1D* hM   = new TH1D("h_w_macro",    ";AMM weight (macro);events", 200, -0.5, 0.5);
    TH1D* hF   = new TH1D("h_w_framework",";AMM weight (framework);events", 200, -0.5, 0.5);
    TH1D* hD   = new TH1D("h_w_diff",     ";macro − framework;events", 200, -0.05, 0.05);
    TH2D* hMF  = new TH2D("h_w_macro_vs_framework",
                          ";framework weight;macro weight",
                          100, -0.5, 0.5, 100, -0.5, 0.5);

    // ----- 5. Event loop with YAML 'only_emu' selection -----
    const Long64_t n = t->GetEntries();
    std::cout << "Total events: " << n << "\n";

    Long64_t nPass = 0;
    double sumMacro = 0., sumFW = 0., sumAbsDiff = 0.;

    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        if (i % 100000 == 0) std::cout << "  entry " << i << " / " << n << "\n";

        // YAML 'only_emu' cuts:
        //   JET_N baselineJvtLoose 20000 >= 2
        //   EL_N  4500  == 1
        //   MU_N  15000 == 1
        const int nJ = count_good_jets(jet_pt, jet_jvt_loose, 20000.);
        if (nJ < 2) continue;

        const int nE = count_above(el_pt, 4500.);
        if (nE != 1) continue;

        const int nM = count_above(mu_pt, 15000.);
        if (nM != 1) continue;

        // AMM weight: sum over loose electrons (only loose=true counted, tight = w_T,
        // loose-only = w_L)
        w_macro = 0.;
        if (el_pt && el_loose && el_tight && el_eta) {
            const size_t N = el_pt->size();
            for (size_t j = 0; j < N; ++j) {
                const bool isLoose = (j < el_loose->size()) && (*el_loose)[j];
                if (!isLoose) continue;
                const bool isTight = (j < el_tight->size()) && (*el_tight)[j];
                const double pt_gev = (*el_pt)[j] * 1e-3;
                const double eta    = (*el_eta)[j];
                const double er = get_eff(hRealE, pt_gev, eta);
                const double ef = get_eff(hFakeE, pt_gev, eta);
                w_macro += amm_weight(isTight, er, ef);
            }
        }

        w_diff = w_macro - fw_framework;

        ++nPass;
        sumMacro   += w_macro;
        sumFW      += fw_framework;
        sumAbsDiff += std::abs(w_diff);

        hM ->Fill(w_macro);
        hF ->Fill(fw_framework);
        hD ->Fill(w_diff);
        hMF->Fill(fw_framework, w_macro);

        o->Fill();
    }

    // ----- 6. Summary -----
    printf("\n=========================================================\n");
    printf("  Events passing 'only_emu' selection: %lld\n", nPass);
    printf("  Sum of weights (macro):     %.4f   (predicted fake yield, this macro)\n", sumMacro);
    printf("  Sum of weights (framework): %.4f   (from fake_weight_nominal)\n", sumFW);
    printf("  Δ yield = macro − framework: %+.4f  (%.2f%%)\n",
           sumMacro - sumFW,
           100. * (sumMacro - sumFW) / (std::abs(sumFW) > 0 ? std::abs(sumFW) : 1.));
    printf("  Mean |Δ| per event:         %.6f\n", nPass ? sumAbsDiff/double(nPass) : 0.);
    printf("=========================================================\n");

    // ----- 7. Save -----
    fo->cd();
    o ->Write();
    hM ->Write();
    hF ->Write();
    hD ->Write();
    hMF->Write();
    fo->Close();
    fi->Close();
    fe->Close();
    printf("\n✓ Wrote %s\n", outFile);
    printf("    tree:   comparison  (per-event row)\n");
    printf("    hists:  h_w_macro, h_w_framework, h_w_diff, h_w_macro_vs_framework\n");
}
