// =============================================================================
// FILE:    studies/smoothing_test/amm_compare_smoothed.C
// PURPOSE: Apply the AMM weight to data16 events with TWO ε_f maps in one
//          pass: the raw binned map and the smoothed map. Persist both
//          per-event weights so the next macro can plot them side-by-side.
//
// INPUT:   /eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data16.root
//          ../../efficiency/outputs/fake_eff_with_systematics.root  (raw)
//          outputs/fake_eff_smoothed.root                            (smooth)
//          ../../efficiency/outputs/real_eff_histograms.root         (ε_r)
// OUTPUT:  outputs/data16_amm_compare.root
//             tree 'cmp' with branches:
//               eventNumber, l2_pt, l2_eta, l2_loose, l2_tight,
//               w_amm_raw, w_amm_smooth
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TTree.h>
#include "ROOT/RVec.hxx"
#include <cmath>
#include <iostream>

using ROOT::VecOps::RVec;

// ---------------------------------------------------------------------------
static double amm_weight(bool isTight, double er, double ef)
{
    er = std::max(1e-5, std::min(1.0 - 1e-5, er));
    ef = std::max(1e-5, std::min(er - 1e-5, ef));
    const double denom = er - ef;
    if (std::abs(denom) < 1e-9) return 0.0;
    return isTight ? (ef / denom) * (er - 1.0)
                   : (ef / denom) *  er;
}

static double get_eff(TH2F* h, double pt_gev, double eta)
{
    if (!h) return 0.0;
    int bx = std::max(1, std::min(h->GetXaxis()->FindBin(pt_gev), h->GetNbinsX()));
    int by = std::max(1, std::min(h->GetYaxis()->FindBin(std::abs(eta)), h->GetNbinsY()));
    return h->GetBinContent(bx, by);
}

// ---------------------------------------------------------------------------
void amm_compare_smoothed()
{
    TFile* fRaw = TFile::Open("../../efficiency/outputs/fake_eff_with_systematics.root", "READ");
    TFile* fSm  = TFile::Open("outputs/fake_eff_smoothed.root", "READ");
    TFile* fRe  = TFile::Open("../../efficiency/outputs/real_eff_histograms.root", "READ");
    if (!fRaw || !fSm || !fRe || fRaw->IsZombie() || fSm->IsZombie() || fRe->IsZombie()) {
        std::cerr << "ERROR: cannot open one of the efficiency files\n"; return;
    }
    TH2F* hFakeRaw = (TH2F*)fRaw->Get("fake_rate_nominal");
    TH2F* hFakeSm  = (TH2F*)fSm ->Get("fake_rate_nominal");
    TH2F* hRT      = (TH2F*)fRe ->Get("h2_tight_eta_pt_CR");
    TH2F* hRL      = (TH2F*)fRe ->Get("h2_loose_eta_pt_CR");
    if (!hFakeRaw || !hFakeSm || !hRT || !hRL) {
        std::cerr << "ERROR: required histograms missing\n"; return;
    }
    TH2F* hRealEff = (TH2F*)hRT->Clone("real_eff_el");
    hRealEff->Divide(hRL);

    // ----- data16 ntuple -----
    TChain tin("tree");
    tin.Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data16.root");
    if (tin.GetEntries() == 0) {
        std::cerr << "ERROR: data16 chain empty\n"; return;
    }

    // ----- Read only what we need -----
    RVec<float>* el_pt    = nullptr;
    RVec<float>* el_eta   = nullptr;
    RVec<char>*  el_loose = nullptr;
    RVec<char>*  el_tight = nullptr;
    ULong64_t    eventNumber = 0;
    tin.SetBranchAddress("el_pt_NOSYS",           &el_pt);
    tin.SetBranchAddress("el_eta",                &el_eta);
    tin.SetBranchAddress("el_select_loose_NOSYS", &el_loose);
    tin.SetBranchAddress("el_select_tight_NOSYS", &el_tight);
    tin.SetBranchAddress("eventNumber",           &eventNumber);

    // ----- Output -----
    gSystem->mkdir("outputs", true);
    TFile* fOut = TFile::Open("outputs/data16_amm_compare.root", "RECREATE");
    fOut->SetCompressionLevel(1);
    TTree* tout = new TTree("cmp", "data16 AMM comparison: raw vs smoothed ε_f");
    Double_t w_raw=0., w_smooth=0., l2_pt_out=0., l2_eta_out=0.;
    Bool_t   l2_loose_out=false, l2_tight_out=false;
    ULong64_t evt_out = 0;
    tout->Branch("eventNumber", &evt_out,        "eventNumber/l");
    tout->Branch("l2_pt",       &l2_pt_out,      "l2_pt/D");
    tout->Branch("l2_eta",      &l2_eta_out,     "l2_eta/D");
    tout->Branch("l2_loose",    &l2_loose_out,   "l2_loose/O");
    tout->Branch("l2_tight",    &l2_tight_out,   "l2_tight/O");
    tout->Branch("w_amm_raw",   &w_raw,          "w_amm_raw/D");
    tout->Branch("w_amm_smooth",&w_smooth,       "w_amm_smooth/D");

    const Long64_t n = tin.GetEntries();
    std::cout << "data16 entries: " << n << "\n";

    long nPass = 0;
    double sumRaw = 0., sumSm = 0.;

    for (Long64_t i = 0; i < n; ++i) {
        if (i % 1000000 == 0) std::cout << "  " << i << " / " << n << "\n";
        tin.GetEntry(i);

        if (!el_pt || el_pt->empty()) continue;

        // Sum AMM over loose electrons (mirrors AMM_ElectronOnly.C)
        w_raw = 0.; w_smooth = 0.;
        for (size_t j = 0; j < el_pt->size(); ++j) {
            bool isLoose = (j < el_loose->size()) && (*el_loose)[j];
            if (!isLoose) continue;
            bool isTight = (j < el_tight->size()) && (*el_tight)[j];
            double pt_gev = (*el_pt)[j] * 1e-3;
            double eta    = (*el_eta)[j];
            double er     = get_eff(hRealEff, pt_gev, eta);
            double efRaw  = get_eff(hFakeRaw, pt_gev, eta);
            double efSm   = get_eff(hFakeSm,  pt_gev, eta);
            w_raw    += amm_weight(isTight, er, efRaw);
            w_smooth += amm_weight(isTight, er, efSm);
        }

        // Store leading-electron kinematics (or 0 if none)
        l2_pt_out  = (el_pt    && el_pt->size())    ? (double)(*el_pt)[0]    * 1e-3 : 0.;
        l2_eta_out = (el_eta   && el_eta->size())   ? (double)(*el_eta)[0]          : 0.;
        l2_loose_out = (el_loose && el_loose->size()) && (*el_loose)[0];
        l2_tight_out = (el_tight && el_tight->size()) && (*el_tight)[0];
        evt_out      = eventNumber;

        if (w_raw == 0. && w_smooth == 0.) continue;  // skip events with no loose el
        ++nPass;
        sumRaw += w_raw;
        sumSm  += w_smooth;
        tout->Fill();
    }

    fOut->cd();
    tout->Write();
    fOut->Close();

    std::cout << "\n=== Summary ===\n";
    std::cout << "Events with non-zero weight: " << nPass << "\n";
    std::cout << "Sum w_raw    = " << sumRaw << "\n";
    std::cout << "Sum w_smooth = " << sumSm  << "\n";
    std::cout << "Δ yield (smooth - raw) = " << (sumSm - sumRaw)
              << "  (" << (sumRaw != 0. ? 100.*(sumSm-sumRaw)/std::abs(sumRaw) : 0.)
              << "% of raw)\n";
    std::cout << "\n✓ Wrote outputs/data16_amm_compare.root\n";
}
