// =============================================================================
// FILE:       data_eff/fake_bg_estimation.C
// PURPOSE:    Fill tight and loose electron histograms from data in the same-sign CR.
//             No truth filtering is applied — all electrons passing selection are
//             counted (both real and fake). After subtracting the MC prompt histograms
//             (from MC_eff/fake_bg_estimation.C), the residual is the fake population.
//             MET-split sub-regions (CR_MET_low, CR_MET_high) are filled to evaluate
//             the MET-dependence systematic on the fake rate.
// INPUTS:     Data ntuples (TChain) — periods A through I of Run 2 data
// OUTPUTS:    Data_histograms.root  (tight/loose 2D and 1D histograms per region)
// RUNS AFTER: nothing (runs in parallel with MC_eff/fake_bg_estimation.C)
// KEY PHYSICS: The inclusive SS data sample is enriched in fake electrons because
//             SS μe events cannot come from prompt OS sources at leading order.
//             The MET split separates QCD multijet (low MET) from W+jets (high MET)
//             fake contributions to assess systematic uncertainty.
// =============================================================================
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <iostream>
#include "ROOT/RVec.hxx"
#include "../utils/SelectionUtils.h"
#include <string>

using namespace ROOT;
using namespace ROOT::VecOps;

void data_eff() {
    TChain *tree = new TChain("tree");
     tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data15.root");
     tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data16.root");
     tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_d/data17.root");
     tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/data18.root");

    // --- Central EventData + shared branch setup ---
    EventData ev;
    SetupBranches(tree, ev);

    // --- Extra branches not in EventData ---
    RVecF* el_pt_NOSYS = nullptr;   // per-electron pT in MeV (histogram filling)
    RVecF* el_d0_NOSYS = nullptr;   // d0 value (not sig) for d0 histograms
    tree->SetBranchAddress("el_pt_NOSYS",  &el_pt_NOSYS);
    tree->SetBranchAddress("el_d0_NOSYS",  &el_d0_NOSYS);

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
#define DEFINE_REGION_HISTOS(region) \
    /* 1D pT (matches MC_eff.C and 2D map x-range: 0–200 GeV) */ \
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
    TH1D* h_tight_d0_##region = new TH1D("h_tight_d0_" #region, ";|d_{0}| [mm];fake rate", 20, 0, 0.3); \
    h_tight_d0_##region->Sumw2(); \
    TH1D* h_loose_d0_##region = new TH1D("h_loose_d0_" #region, ";|d_{0}| [mm];Loose electrons", 20, 0, 0.3); \
    h_loose_d0_##region->Sumw2(); \
    \
    /* 1D d0sig */ \
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
    /* 1D Mjj */ \
    TH1D* h_tight_mjj_##region = new TH1D("h_tight_mjj_" #region, ";Mjj [GeV];fake rate", 60, 0, 2000); \
    h_tight_mjj_##region->Sumw2(); \
    TH1D* h_loose_mjj_##region = new TH1D("h_loose_mjj_" #region, ";Mjj [GeV];Loose electrons", 60, 0, 2000); \
    h_loose_mjj_##region->Sumw2(); \
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
DEFINE_REGION_HISTOS(CR_el_d0ig)
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
Long64_t n_after_VR_mjj_inv = 0;
// Debug counters
Long64_t n_tight_electron_any = 0;
Long64_t n_loose_electron_any = 0;
Long64_t n_tight_electron_prompt = 0;
Long64_t n_loose_electron_prompt = 0;

Long64_t nentries = tree->GetEntries();
std::cout << "Total number of entries in the tree: " << nentries << std::endl;

for (Long64_t i = 0; i < nentries; i++) {
    tree->GetEntry(i);
    std::string currentFile = tree->GetCurrentFile()->GetName();

    ev.isData   = true;
    // ev.isData15 is no longer needed for trigger-matching dispatch
    // (passTriggerMatching now ORs globalTriggerMatch with per-chain flags).

    ComputeDerived(ev);

    // Weight = 1.0 for data
    double w_tight = 1.0;
    double w_loose = 1.0;

    if (!passTopology(ev)) continue;
    n_after_flavour++;

    // Fake-rate measurement: muon-only trigger (electron is the unbiased probe).
    if (!passMuonTrigger(ev)) continue;
    n_after_trigger++;

    if (!passBasePresel_FakeMeasurement(ev)) continue;
    n_after_preselection++;

    // --- Region flags ---
    bool pass_CR = isFakeCR(ev) && ev.no_bjets && (ev.jet1_pt > 0. && ev.jet2_pt > 0.);
    if (pass_CR) n_after_CR++;

    bool pass_SR = ev.is_os && passVBFJets(ev) && ev.no_bjets;
    if (pass_SR) n_after_SR++;

    bool pass_CR_MET_low  = pass_CR && ev.MET < 50000.;
    if (pass_CR_MET_low) n_after_CR_MET_low++;

    bool pass_CR_MET_high = pass_CR && ev.MET > 50000.;
    if (pass_CR_MET_high) n_after_CR_MET_high++;

    // --- Explicit IP flag for tight histogram filling ---
    bool pass_IP = passElectronIP(ev);

    // --- No truth filter for data — just use el ID ---
    bool tight_electron_prompt = ev.el_is_tight;
    bool loose_electron_prompt = ev.el_is_loose;

    // --------------------------------------------------------------------
    // Debug: count tight/loose electrons (global, not region-specific)
    // --------------------------------------------------------------------
    if (pass_CR) {
        if (ev.el_is_tight && pass_IP) n_tight_electron_any++;
        if (el_pt_NOSYS && !el_pt_NOSYS->empty() && (*el_pt_NOSYS)[0]) n_loose_electron_any++;
    }
    if (pass_SR) {
        if (tight_electron_prompt && pass_IP) n_tight_electron_prompt++;
        if (loose_electron_prompt && pass_IP) n_loose_electron_prompt++;
    }

    // --- Per-electron variable extraction ---
    double el_pt_val    = (el_pt_NOSYS && !el_pt_NOSYS->empty()) ? (*el_pt_NOSYS)[0] * 1e-3 : 0.;
    double el_eta_val   = ev.el_aeta;  // already |η|, set by ComputeDerived
    double el_d0_val    = (el_d0_NOSYS && !el_d0_NOSYS->empty()) ? std::fabs((*el_d0_NOSYS)[0]) : 0.;
    double el_d0sig_val = (ev.el_d0sig && !ev.el_d0sig->empty()) ? std::fabs((*ev.el_d0sig)[0]) : 0.;
    double met_val      = ev.MET * 1e-3;

    // =========================================================================
    // Fill histograms for each region if the region flag is true
    // =========================================================================

    // For CR
    if (pass_CR && pass_IP) {
        if (tight_electron_prompt && pass_IP) {
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
        if (tight_electron_prompt && pass_IP) {
            h_tight_pt_SR->Fill(el_pt_val, w_tight);
            h_tight_eta_SR->Fill(el_eta_val, w_tight);
            h_tight_d0_SR->Fill(el_d0_val, w_tight);
            h_tight_d0sig_SR->Fill(el_d0sig_val, w_tight);
            h_tight_met_SR->Fill(met_val, w_tight);
            h2_tight_eta_pt_SR->Fill(el_pt_val, el_eta_val, w_tight);
        }
        if (loose_electron_prompt) {
            h_loose_pt_SR->Fill(el_pt_val, w_loose);
            h_loose_eta_SR->Fill(el_eta_val, w_loose);
            h_loose_d0_SR->Fill(el_d0_val, w_loose);
            h_loose_d0sig_SR->Fill(el_d0sig_val, w_loose);
            h_loose_met_SR->Fill(met_val, w_loose);
            h2_loose_eta_pt_SR->Fill(el_pt_val, el_eta_val, w_loose);
        }
    }

    // For CR_MET_low
    if (pass_CR_MET_low && pass_IP) {
        if (tight_electron_prompt && pass_IP) {
            h_tight_pt_CR_MET_low->Fill(el_pt_val, w_tight);
            h_tight_eta_CR_MET_low->Fill(el_eta_val, w_tight);
            h_tight_d0_CR_MET_low->Fill(el_d0_val, w_tight);
            h_tight_d0sig_CR_MET_low->Fill(el_d0sig_val, w_tight);
            h_tight_met_CR_MET_low->Fill(met_val, w_tight);
            h2_tight_eta_pt_CR_MET_low->Fill(el_pt_val, el_eta_val, w_tight);
        }
        if (loose_electron_prompt) {
            h_loose_pt_CR_MET_low->Fill(el_pt_val, w_loose);
            h_loose_eta_CR_MET_low->Fill(el_eta_val, w_loose);
            h_loose_d0_CR_MET_low->Fill(el_d0_val, w_loose);
            h_loose_d0sig_CR_MET_low->Fill(el_d0sig_val, w_loose);
            h_loose_met_CR_MET_low->Fill(met_val, w_loose);
            h2_loose_eta_pt_CR_MET_low->Fill(el_pt_val, el_eta_val, w_loose);
        }
    }

    // For CR_MET_high
    if (pass_CR_MET_high && pass_IP) {
        if (tight_electron_prompt ) {
            h_tight_pt_CR_MET_high->Fill(el_pt_val, w_tight);
            h_tight_eta_CR_MET_high->Fill(el_eta_val, w_tight);
            h_tight_d0_CR_MET_high->Fill(el_d0_val, w_tight);
            h_tight_d0sig_CR_MET_high->Fill(el_d0sig_val, w_tight);
            h_tight_met_CR_MET_high->Fill(met_val, w_tight);
            h2_tight_eta_pt_CR_MET_high->Fill(el_pt_val, el_eta_val, w_tight);
        }
        if (loose_electron_prompt) {
            h_loose_pt_CR_MET_high->Fill(el_pt_val, w_loose);
            h_loose_eta_CR_MET_high->Fill(el_eta_val, w_loose);
            h_loose_d0_CR_MET_high->Fill(el_d0_val, w_loose);
            h_loose_d0sig_CR_MET_high->Fill(el_d0sig_val, w_loose);
            h_loose_met_CR_MET_high->Fill(met_val, w_loose);
            h2_loose_eta_pt_CR_MET_high->Fill(el_pt_val, el_eta_val, w_loose);
        }
    }
    }


h_tight_pt_CR->Draw("E");          // error bars
h_tight_pt_CR->Draw("HIST E");     // histogram + error bars

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
}
std::cout << std::endl;

std::cout << "Debug lepton counts (global):" << std::endl;
std::cout << "  Tight electrons (any):                " << n_tight_electron_any << std::endl;
std::cout << "  Loose electrons (any):                " << n_loose_electron_any << std::endl;
std::cout << "  Tight prompt electrons:               " << n_tight_electron_prompt << std::endl;
std::cout << "  Loose prompt electrons:               " << n_loose_electron_prompt << std::endl;
std::cout << "=======================================================" << std::endl;
std::cout << std::endl;

// =============================================
// Save all histograms to ROOT file
// =============================================
TFile* f_out = new TFile("outputs/Data_histograms.root", "RECREATE");
if (!f_out || f_out->IsZombie()) {
    std::cerr << "ERROR: Could not create output file data_histograms.root" << std::endl;
} else {
    std::cout << "Saving histograms to data_histograms.root ..." << std::endl;

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
    WRITE_REGION(CR_el_d0)
//    WRITE_REGION(CR_el_d0sig)
  //  WRITE_REGION(VR_mjj_inv)
    // If you kept the inclusive / global ones
    WRITE_REGION(incl)


    f_out->Write();
    f_out->Close();
    std::cout << "All histograms saved successfully to Data_histograms.root" << std::endl;
}



}
