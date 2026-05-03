// =============================================================================
// FILE:       MC_eff/fake_bg_estimation.C
// PURPOSE:    Fill tight and loose electron histograms from MC truth-filtered events.
//             Only PROMPT electrons (el_truth_isPrompt == 1, el_IFFClass == 2) are
//             counted. These histograms are the MC subtraction component of the fake
//             rate: fake_rate = (data - MC_prompt) / (loose_data - loose_MC_prompt).
//             Also fills histograms for the real electron efficiency ε_r = tight/loose
//             in the CR and SR (needed by the Matrix Method weight formula).
// INPUTS:     MC ntuples (TChain) — ttbar, W+jets, Z+jets, diboson, single-top
// OUTPUTS:    MC_histograms.root  (tight/loose 2D and 1D histograms per region)
// RUNS AFTER: nothing (first step in the chain)
// KEY PHYSICS: The MC subtraction removes real prompt electrons that contaminate
//             the SS data CR. Without this, the fake rate would be diluted by
//             real electrons and under-estimated.
// =============================================================================
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <iostream>
#include "ROOT/RVec.hxx"
#include "../utils/SelectionUtils.h"


using namespace ROOT;
using namespace ROOT::VecOps;

void MC_eff() {

      	TChain *tree = new TChain("tree");

     tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/ztautau/*.root");
     tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/wjets/*.root");
     tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/ttbar/*.root");
     tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/singletop/*.root");
     tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/diboson/*.root");
     tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/higgs/*.root");

    // --- Central EventData + shared branch setup ---
    EventData ev;
    SetupBranches(tree, ev);

    // --- Extra branches not in EventData ---
    RVecF* el_pt_NOSYS  = nullptr;   // per-electron pT in MeV (histogram filling)
    RVecI* el_IFFClass  = nullptr;   // truth IFF class
    RVecI* mu_IFFClass  = nullptr;   // muon IFF class
    RVecF* el_d0_NOSYS  = nullptr;   // d0 value (not sig) for d0 histograms
    // Charge mis-ID efficiency SF (data/MC ratio for the rate at which a
    // truth-prompt electron's charge is mis-measured). Used to up-weight
    // charge-flipped prompt MC events so their density in SS regions matches
    // the data charge-flip rate. Default 1.0 if branch missing.
    RVecF* el_qmis_sf   = nullptr;
    tree->SetBranchAddress("el_pt_NOSYS",  &el_pt_NOSYS);
    tree->SetBranchAddress("el_IFFClass",  &el_IFFClass);
    tree->SetBranchAddress("mu_IFFClass",  &mu_IFFClass);
    tree->SetBranchAddress("el_d0_NOSYS",  &el_d0_NOSYS);
    tree->SetBranchAddress("el_charge_misid_effSF_tight_NOSYS", &el_qmis_sf);

    // Define counters for weighted event counts + cutflow steps
    double count_trigger_and_match = 0.0;
    double count_tight_prompt_muons = 0.0;
    double count_tight_prompt_electrons = 0.0;
    double count_loose_prompt_electrons = 0.0;
    double count_preselection = 0.0;
    double count_trigger_match_preselection = 0.0;
    double count_after_trigger_match   = 0.0;
    double count_after_preselection    = 0.0;


// ───────────────────────────────────────────────────────────────
// Helper macro to define region histograms with Sumw2
// ───────────────────────────────────────────────────────────────
// IMPORTANT: 1D pT histograms now span 0–200 GeV (was 0–80) so diagnostic
// plots match the full pT range used by the 2D map (kPtBins ends at 200).
// Without this, all electrons with pT > 80 GeV silently fell into overflow
// and the upper-pT ε_f tail was invisible in 1D inspections — even though
// the 2D map (which feeds the AMM) was always binning them correctly.
#define DEFINE_REGION_HISTOS(region) \
    /* 1D pT (matches 2D map x-range: 0–200 GeV) */ \
    TH1D* h_tight_pt_##region = new TH1D("h_tight_pt_" #region, ";Electron p_{T} [GeV];fake rate", 50, 0, 200); \
    h_tight_pt_##region->Sumw2(); \
    TH1D* h_loose_pt_##region = new TH1D("h_loose_pt_" #region, ";Electron p_{T} [GeV];Loose electrons", 50, 0, 200); \
    h_loose_pt_##region->Sumw2(); \
    \
    /* 1D |eta| */ \
    TH1D* h_tight_eta_##region = new TH1D("h_tight_eta_" #region, ";|Electron #eta|;fake rate", 10, 0, 2.5); \
    h_tight_eta_##region->Sumw2(); \
    TH1D* h_loose_eta_##region = new TH1D("h_loose_eta_" #region, ";|Electron #eta|;Loose electrons", 10, 0, 2.5); \
    h_loose_eta_##region->Sumw2(); \
    \
    /* 1D |d0| */ \
    TH1D* h_tight_d0_##region = new TH1D("h_tight_d0_" #region, ";|d_{0}| ;fake rate", 20, 0, 0.3); \
    h_tight_d0_##region->Sumw2(); \
    TH1D* h_loose_d0_##region = new TH1D("h_loose_d0_" #region, ";|d_{0}|;Loose electrons", 20, 0, 0.3); \
    h_loose_d0_##region->Sumw2(); \
    \
    /* 1D |d0| */ \
    TH1D* h_tight_d0sig_##region = new TH1D("h_tight_d0sig_" #region, ";d0sig ;fake rate", 20, 0, 7); \
    h_tight_d0sig_##region->Sumw2(); \
    TH1D* h_loose_d0sig_##region = new TH1D("h_loose_d0sig_" #region, ";d0sig;Loose electrons", 20, 0, 7); \
    h_loose_d0sig_##region->Sumw2(); \
    \
    /* 1D MET */ \
    TH1D* h_tight_met_##region = new TH1D("h_tight_met_" #region, ";MET [GeV];fake rate", 50, 0, 200); \
    h_tight_met_##region->Sumw2(); \
    TH1D* h_loose_met_##region = new TH1D("h_loose_met_" #region, ";MET [GeV];Loose electrons", 50, 0, 200); \
    h_loose_met_##region->Sumw2(); \
    \
    /* 2D eta vs pT */ \
    TH2D* h2_tight_eta_pt_##region = new TH2D("h2_tight_eta_pt_" #region, ";Electron p_{T} [GeV]; |Electron #eta|", \
                                              nPtBins_2d, ptBins_2d, nEtaBins_2d, etaBins_2d); \
    h2_tight_eta_pt_##region->Sumw2(); \
    TH2D* h2_loose_eta_pt_##region = new TH2D("h2_loose_eta_pt_" #region, ";Electron p_{T} [GeV] ;|Electron #eta|", \
                                              nPtBins_2d, ptBins_2d, nEtaBins_2d, etaBins_2d); \
    h2_loose_eta_pt_##region->Sumw2();

// Keep the bin arrays (same as before)
// Source bin edges centrally from SelectionUtils.h so the framework
// stays self-consistent when kPtBins / kEtaBins change.
const int     nPtBins_2d  = kNPtBins;
const double* ptBins_2d   = kPtBins;
const int     nEtaBins_2d = kNEtaBins;
const double* etaBins_2d  = kEtaBins;

// ───────────────────────────────────────────────────────────────
// Define histograms for each region
// ───────────────────────────────────────────────────────────────
DEFINE_REGION_HISTOS(CR)
DEFINE_REGION_HISTOS(SR)
DEFINE_REGION_HISTOS(CR_MET_low)
DEFINE_REGION_HISTOS(CR_MET_high)
DEFINE_REGION_HISTOS(CR_el_d0)

// Optional: inclusive / total histograms (if you still want them)
DEFINE_REGION_HISTOS(incl)

// Counters (declare and initialize to 0 outside loop)
Long64_t n_after_flavour = 0;
Long64_t n_after_trigger = 0;
Long64_t n_after_preselection = 0;
Long64_t n_after_CR = 0;
Long64_t n_after_SR = 0;
Long64_t n_after_CR_MET_low = 0;
Long64_t n_after_CR_MET_high = 0;
Long64_t n_after_CR_el_d0 = 0;

// Debug counters
Long64_t n_tight_electron_any = 0;
Long64_t n_loose_electron_any = 0;
Long64_t n_tight_electron_prompt = 0;
Long64_t n_loose_electron_prompt = 0;

Long64_t nentries = tree->GetEntries();
std::cout << "Total number of entries in the tree: " << nentries << std::endl;

for (Long64_t i = 0; i < nentries; i++) {
    tree->GetEntry(i);
    ComputeDerived(ev);

    if (!passTopology(ev)) continue;
    n_after_flavour++;

    double total_weight_tight = mcWeightTight(ev);
    double total_weight_loose = mcWeightLoose(ev);

    // Fake-rate measurement: muon-only trigger to avoid biasing the electron
    // (single-e triggers force a near-tight electron at HLT, inflating ε_f).
    if (!passMuonTrigger(ev)) continue;
    n_after_trigger++;

    if (!passBasePresel_FakeMeasurement(ev)) continue;
    n_after_preselection++;

    // --- Region flags ---
    // CR: SS + at least 2 jets (jet1_pt>0 && jet2_pt>0) + no b-jets (VBF not required for fake rate CR)
    bool pass_CR = isFakeCR(ev) && ev.no_bjets && (ev.jet1_pt > 0. && ev.jet2_pt > 0.);
    if (pass_CR) n_after_CR++;

    bool pass_SR = ev.is_os && passVBFJets(ev) && ev.no_bjets;
    if (pass_SR) n_after_SR++;

    bool pass_CR_MET_low  = pass_CR && ev.MET < 50000.;
    if (pass_CR_MET_low) n_after_CR_MET_low++;

    bool pass_CR_MET_high = pass_CR && ev.MET >= 50000.;
    if (pass_CR_MET_high) n_after_CR_MET_high++;

    // Tighter IP requirement (|d0sig|<3, |z0sinθ|<0.3) on top of the
    // standard event-level cut already applied by passBasePresel.  This
    // sharpens the loose denominator used for ε_r / ε_f, removing the
    // displaced-track tail that artificially depresses ε_r.  Both the
    // tight numerator AND the loose denominator are gated by `pass_IP`
    // below so the ratio remains a clean conditional probability.
    bool pass_IP = passElectronIPTight(ev);

    // --- Truth-level "prompt-like" for electron ---
    // Subtract IFFClass==2 (PromptIsoElectron, the standard prompt) AND
    // IFFClass==5 (PromptPhotonConversion).  Conversions are real electrons
    // from a prompt photon; if we left them in the SS-CR data they would
    // pass tight at high efficiency and inflate ε_f at high pT.  Treating
    // them as prompt-like for the data−MC subtraction follows the convention
    // of the ATLAS HNL multilepton (arXiv:2204.11138) and FakeBkgTools
    // (JINST 18 T11004) papers.
    int  el_iff = (el_IFFClass && !el_IFFClass->empty()) ? (*el_IFFClass)[0] : -1;
    bool iff_prompt_like   = (el_iff == 2);
    bool electron_is_prompt = isMCPromptElectron(ev) && iff_prompt_like;

    bool tight_electron_prompt = ev.el_is_tight && electron_is_prompt;
    bool loose_electron_prompt = ev.el_is_loose && electron_is_prompt;

    // --------------------------------------------------------------------
    // Debug: count tight/loose electrons (global, not region-specific)
    // --------------------------------------------------------------------
    if (ev.el_is_tight) n_tight_electron_any++;
    if (ev.el_is_loose) n_loose_electron_any++;

    if (pass_SR) {
        if (tight_electron_prompt && pass_IP) n_tight_electron_prompt++;
        if (loose_electron_prompt)            n_loose_electron_prompt++;
    }

    // --- Per-electron variable extraction ---
    double el_pt_val    = (el_pt_NOSYS && !el_pt_NOSYS->empty()) ? (*el_pt_NOSYS)[0] * 1e-3 : 0.;
    double el_eta_val   = ev.el_aeta;  // already |η|, set by ComputeDerived
    double el_d0_val    = (el_d0_NOSYS && !el_d0_NOSYS->empty()) ? std::fabs((*el_d0_NOSYS)[0]) : 0.;
    double el_d0sig_val = (ev.el_d0sig && !ev.el_d0sig->empty()) ? std::fabs((*ev.el_d0sig)[0]) : 0.;
    double met_val      = ev.MET * 1e-3;

    // Charge mis-ID SF: data/MC ratio for the charge-flip rate of prompt
    // electrons. When applied multiplicatively to PROMPT MC weights it
    // up-weights charge-flipped prompt events so MC matches data; non-flipped
    // events have SF≈1.
    //
    // The branch uses -1.0 as a sentinel value (no electron / selection
    // didn't run / IFF classifier failed). Treat sentinels as SF=1 so
    // they don't corrupt the prompt subtraction with negative weights.
    double sf_qmis = (el_qmis_sf && !el_qmis_sf->empty())
                     ? (double)(*el_qmis_sf)[0]
                     : 1.0;
    if (sf_qmis < 0.0 || !std::isfinite(sf_qmis)) sf_qmis = 1.0;

    double w_tight = total_weight_tight * sf_qmis;
    double w_loose = total_weight_loose * sf_qmis;

    // =========================================================================
    // Fill histograms for each region
    // =========================================================================

    // For CR
    if (pass_CR && pass_IP) {
        if (tight_electron_prompt) {
            h_tight_pt_CR->Fill(el_pt_val, w_tight);
            h_tight_eta_CR->Fill(el_eta_val, w_tight);
            h_tight_d0_CR->Fill(el_d0_val, w_tight);
            h_tight_d0sig_CR->Fill(el_d0sig_val, w_tight);
            h_tight_met_CR->Fill(met_val, w_tight);
            h2_tight_eta_pt_CR->Fill(el_pt_val, el_eta_val, w_tight);
        }
        if (loose_electron_prompt) {
            h_loose_pt_CR->Fill(el_pt_val, w_loose);
            h_loose_eta_CR->Fill(el_eta_val, w_loose);
            h_loose_d0_CR->Fill(el_d0_val, w_loose);
            h_loose_d0sig_CR->Fill(el_d0sig_val, w_loose);
            h_loose_met_CR->Fill(met_val, w_loose);
            h2_loose_eta_pt_CR->Fill(el_pt_val, el_eta_val, w_loose);
        }
    }

    if (pass_SR && pass_IP) {
        if (tight_electron_prompt ) {
            h_tight_pt_SR->Fill(el_pt_val, w_tight);
            h_tight_eta_SR->Fill(el_eta_val, w_tight);
            h_tight_d0_SR->Fill(el_d0_val, w_tight);
            h_tight_met_SR->Fill(met_val, w_tight);
            h2_tight_eta_pt_SR->Fill(el_pt_val, el_eta_val, w_tight);
        }
        if (loose_electron_prompt) {
            h_loose_pt_SR->Fill(el_pt_val, w_loose);
            h_loose_eta_SR->Fill(el_eta_val, w_loose);
            h_loose_d0_SR->Fill(el_d0_val, w_loose);
            h_loose_met_SR->Fill(met_val, w_loose);
            h2_loose_eta_pt_SR->Fill(el_pt_val, el_eta_val, w_loose);
        }
    }

    // For CR_MET_low (pass_IP gates BOTH tight num and loose denom — symmetric)
    if (pass_CR_MET_low && pass_IP) {
        if (tight_electron_prompt) {
            h_tight_pt_CR_MET_low->Fill(el_pt_val, w_tight);
            h_tight_eta_CR_MET_low->Fill(el_eta_val, w_tight);
            h_tight_d0_CR_MET_low->Fill(el_d0_val, w_tight);
            h_tight_met_CR_MET_low->Fill(met_val, w_tight);
            h2_tight_eta_pt_CR_MET_low->Fill(el_pt_val, el_eta_val, w_tight);
        }
        if (loose_electron_prompt) {
            h_loose_pt_CR_MET_low->Fill(el_pt_val, w_loose);
            h_loose_eta_CR_MET_low->Fill(el_eta_val, w_loose);
            h_loose_d0_CR_MET_low->Fill(el_d0_val, w_loose);
            h_loose_met_CR_MET_low->Fill(met_val, w_loose);
            h2_loose_eta_pt_CR_MET_low->Fill(el_pt_val, el_eta_val, w_loose);
        }
    }

    // For CR_MET_high (pass_IP gates BOTH tight num and loose denom — symmetric)
    if (pass_CR_MET_high && pass_IP) {
        if (tight_electron_prompt) {
            h_tight_pt_CR_MET_high->Fill(el_pt_val, w_tight);
            h_tight_eta_CR_MET_high->Fill(el_eta_val, w_tight);
            h_tight_d0_CR_MET_high->Fill(el_d0_val, w_tight);
            h_tight_met_CR_MET_high->Fill(met_val, w_tight);
            h2_tight_eta_pt_CR_MET_high->Fill(el_pt_val, el_eta_val, w_tight);
        }
        if (loose_electron_prompt) {
            h_loose_pt_CR_MET_high->Fill(el_pt_val, w_loose);
            h_loose_eta_CR_MET_high->Fill(el_eta_val, w_loose);
            h_loose_d0_CR_MET_high->Fill(el_d0_val, w_loose);
            h_loose_met_CR_MET_high->Fill(met_val, w_loose);
            h2_loose_eta_pt_CR_MET_high->Fill(el_pt_val, el_eta_val, w_loose);
        }
    }

}


h_tight_pt_SR->Draw("E");          // error bars
h_tight_pt_SR->Draw("HIST E");     // histogram + error bars

// =============================================
// Print selection summary and counters
// =============================================
std::cout << std::endl;
std::cout << "=======================================================" << std::endl;
std::cout << "          EVENT SELECTION SUMMARY" << std::endl;
std::cout << "=======================================================" << std::endl;
std::cout << "Total entries in tree:                  " << nentries << std::endl;
std::cout << "After flavour (mu-e pair):              " << n_after_flavour << std::endl;
std::cout << "After trigger + muon matching:          " << n_after_trigger << std::endl;
std::cout << "After preselection (no el d0 cut):      " << n_after_preselection << std::endl;
std::cout << "After CR (SS region):                   " << n_after_CR << std::endl;
std::cout << "After SR (OS region):                   " << n_after_SR << std::endl;
std::cout << "After CR_MET_low:                       " << n_after_CR_MET_low << std::endl;
std::cout << "After CR_MET_high:                      " << n_after_CR_MET_high << std::endl;
std::cout << "-------------------------------------------------------" << std::endl;

// Optional: efficiencies
if (nentries > 0) {
    std::cout << "Flavour eff:                            "
              << Form("%.3f %%", 100.0 * n_after_flavour / nentries) << std::endl;
}
if (n_after_flavour > 0) {
    std::cout << "Trigger+matching eff:                   "
              << Form("%.3f %%", 100.0 * n_after_trigger / n_after_flavour) << std::endl;
}
if (n_after_trigger > 0) {
    std::cout << "Preselection eff:                       "
              << Form("%.3f %%", 100.0 * n_after_preselection / n_after_trigger) << std::endl;
}
if (n_after_preselection > 0) {
    std::cout << "CR eff (from preselection):             "
              << Form("%.3f %%", 100.0 * n_after_CR / n_after_preselection) << std::endl;
    std::cout << "SR eff (from preselection):             "
              << Form("%.3f %%", 100.0 * n_after_SR / n_after_preselection) << std::endl;
}
std::cout << std::endl;

std::cout << "Debug lepton counts (global):" << std::endl;
std::cout << "  Tight electrons (any):                " << n_tight_electron_any << std::endl;
std::cout << "  Loose electrons (any):                " << n_loose_electron_any << std::endl;
std::cout << "  Tight prompt electrons:               " << n_tight_electron_prompt << std::endl;
std::cout << "  Loose prompt electrons:               " << n_loose_electron_prompt << std::endl;
std::cout << "=======================================================" << std::endl;
std::cout << std::endl;

std::cout << "==== HIST SUMMARY ====" << std::endl;
std::cout << "CR tight_pt: entries=" << h_tight_pt_CR->GetEntries()
          << " integral=" << h_tight_pt_CR->Integral() << std::endl;
std::cout << "SR tight_pt: entries=" << h_tight_pt_SR->GetEntries()
          << " integral=" << h_tight_pt_SR->Integral() << std::endl;
std::cout << "CR loose_pt: entries=" << h_loose_pt_CR->GetEntries()
          << " integral=" << h_loose_pt_CR->Integral() << std::endl;
std::cout << "SR loose_pt: entries=" << h_loose_pt_SR->GetEntries()
          << " integral=" << h_loose_pt_SR->Integral() << std::endl;


// =============================================
// Save all histograms to ROOT file
// =============================================
TFile* f_out = new TFile("outputs/MC_histograms.root", "RECREATE");
if (!f_out || f_out->IsZombie()) {
    std::cerr << "ERROR: Could not create output file MC_histograms.root" << std::endl;
} else {
    std::cout << "Saving histograms to MC_histograms.root ..." << std::endl;

    // Helper macro to write a region's histograms (makes code shorter)
    #define WRITE_REGION(region) \
        if (h_tight_pt_##region)        h_tight_pt_##region->Write(); \
        if (h_loose_pt_##region)        h_loose_pt_##region->Write(); \
        if (h_tight_eta_##region)       h_tight_eta_##region->Write(); \
        if (h_loose_eta_##region)       h_loose_eta_##region->Write(); \
        if (h_tight_d0_##region)        h_tight_d0_##region->Write(); \
        if (h_loose_d0_##region)        h_loose_d0_##region->Write(); \
        if (h_tight_d0sig_##region)        h_tight_d0sig_##region->Write(); \
        if (h_loose_d0sig_##region)        h_loose_d0sig_##region->Write(); \
        if (h_tight_met_##region)       h_tight_met_##region->Write(); \
        if (h_loose_met_##region)       h_loose_met_##region->Write(); \
        if (h2_tight_eta_pt_##region)   h2_tight_eta_pt_##region->Write(); \
        if (h2_loose_eta_pt_##region)   h2_loose_eta_pt_##region->Write();

    // Write all regions
    WRITE_REGION(CR)
    WRITE_REGION(SR)
    WRITE_REGION(CR_MET_low)
    WRITE_REGION(CR_MET_high)
 //   WRITE_REGION(CR_el_d0)

    // If you kept the inclusive / global ones
 //   WRITE_REGION(incl)


    f_out->Write();
    f_out->Close();
    std::cout << "All histograms saved successfully to MC_histograms.root" << std::endl;
}



}
