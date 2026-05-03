// =============================================================================
// FILE:       MC_eff/real/real_eff.C
// PURPOSE:    Measure the real electron efficiency epsilon_r = N(tight)/N(loose)
//             from MC, using prompt-truth-matched electrons in the OS prompt CR.
// INPUTS:     MC ntuples (TChain) — ttbar, W+jets, Z+jets, diboson, single-top, higgs
// OUTPUTS:    MC_histograms.root  (tight/loose 1D and 2D histograms per region)
// PATTERN:    Centralized EventData + SelectionUtils.h, mirrors fake_bg_estimation.C
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

void real_eff() {

    TChain* tree = new TChain("tree");
    tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/ztautau/*.root");
    tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/wjets/*.root");
    tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/ttbar/*.root");
    tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/singletop/*.root");
    tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/diboson/*.root");
    tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/higgs/*.root");

    // --- Central EventData + shared branch setup ---
    EventData ev;
    SetupBranches(tree, ev);

    // --- Extra branches not in EventData (needed for histogram filling) ---
    RVecF* el_pt_NOSYS = nullptr;   // per-electron pT in MeV
    RVecF* el_d0_NOSYS = nullptr;   // per-electron d0 (not significance) for d0 histograms
    tree->SetBranchAddress("el_pt_NOSYS", &el_pt_NOSYS);
    tree->SetBranchAddress("el_d0_NOSYS", &el_d0_NOSYS);

    // ---------------------------------------------------------------
    // Histogram definitions — same names/binning as before so plot_real.C is unchanged
    // ---------------------------------------------------------------
    const int    nPtBins_2d  = kNPtBins;
    const double* ptBins_2d  = kPtBins;
    const int    nEtaBins_2d = kNEtaBins;
    const double* etaBins_2d = kEtaBins;

#define DEFINE_REGION_HISTOS(region) \
    TH1D* h_tight_pt_##region   = new TH1D("h_tight_pt_"  #region, ";Electron p_{T} [GeV];real eff",     20, 25, 80); h_tight_pt_##region->Sumw2(); \
    TH1D* h_loose_pt_##region   = new TH1D("h_loose_pt_"  #region, ";Electron p_{T} [GeV];Loose electrons", 20, 25, 80); h_loose_pt_##region->Sumw2(); \
    TH1D* h_tight_eta_##region  = new TH1D("h_tight_eta_" #region, ";|Electron #eta|;real eff",          10, 0, 2.5); h_tight_eta_##region->Sumw2(); \
    TH1D* h_loose_eta_##region  = new TH1D("h_loose_eta_" #region, ";|Electron #eta|;Loose electrons",   10, 0, 2.5); h_loose_eta_##region->Sumw2(); \
    TH1D* h_tight_d0_##region   = new TH1D("h_tight_d0_"  #region, ";|d_{0}|;real eff",                  20, 0, 0.3); h_tight_d0_##region->Sumw2(); \
    TH1D* h_loose_d0_##region   = new TH1D("h_loose_d0_"  #region, ";|d_{0}|;Loose electrons",           20, 0, 0.3); h_loose_d0_##region->Sumw2(); \
    TH1D* h_tight_d0sig_##region= new TH1D("h_tight_d0sig_" #region, ";d0sig;real eff",                  20, 0, 7);   h_tight_d0sig_##region->Sumw2(); \
    TH1D* h_loose_d0sig_##region= new TH1D("h_loose_d0sig_" #region, ";d0sig;Loose electrons",           20, 0, 7);   h_loose_d0sig_##region->Sumw2(); \
    TH1D* h_tight_met_##region  = new TH1D("h_tight_met_" #region, ";MET [GeV];real eff",                20, 0, 200); h_tight_met_##region->Sumw2(); \
    TH1D* h_loose_met_##region  = new TH1D("h_loose_met_" #region, ";MET [GeV];Loose electrons",         20, 0, 200); h_loose_met_##region->Sumw2(); \
    TH2D* h2_tight_eta_pt_##region = new TH2D("h2_tight_eta_pt_" #region, ";Electron p_{T} [GeV];|Electron #eta|", \
                                              nPtBins_2d, ptBins_2d, nEtaBins_2d, etaBins_2d); h2_tight_eta_pt_##region->Sumw2(); \
    TH2D* h2_loose_eta_pt_##region = new TH2D("h2_loose_eta_pt_" #region, ";Electron p_{T} [GeV];|Electron #eta|", \
                                              nPtBins_2d, ptBins_2d, nEtaBins_2d, etaBins_2d); h2_loose_eta_pt_##region->Sumw2();

    DEFINE_REGION_HISTOS(CR)
    DEFINE_REGION_HISTOS(SR)
    DEFINE_REGION_HISTOS(CR_MET_low)
    DEFINE_REGION_HISTOS(CR_MET_high)

    // ---------------------------------------------------------------
    // Counters
    // ---------------------------------------------------------------
    Long64_t n_after_flavour       = 0;
    Long64_t n_after_trigger       = 0;
    Long64_t n_after_preselection  = 0;
    Long64_t n_after_CR            = 0;
    Long64_t n_after_SR            = 0;
    Long64_t n_after_CR_MET_low    = 0;
    Long64_t n_after_CR_MET_high   = 0;

    Long64_t n_tight_electron_any     = 0;
    Long64_t n_loose_electron_any     = 0;
    Long64_t n_tight_electron_prompt  = 0;
    Long64_t n_loose_electron_prompt  = 0;

    Long64_t nentries = tree->GetEntries();
    std::cout << "Total number of entries in the tree: " << nentries << std::endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        ComputeDerived(ev);

        if (!passTopology(ev)) continue;
        n_after_flavour++;

        double w_tight = mcWeightTight(ev);
        double w_loose = mcWeightLoose(ev);

        // Real-efficiency measurement: muon-only trigger to avoid biasing
        // the electron (the probe). Same convention as ε_f measurement.
        if (!passMuonTriggerOnly(ev)) continue;
        n_after_trigger++;

        if (!passBasePresel_FakeMeasurement(ev)) continue;
        n_after_preselection++;

        // --- Region flags (real-efficiency CR is OS = isPromptCR) ---
        bool pass_CR          = isPromptCR(ev) && ev.no_bjets && (ev.jet1_pt > 0. && ev.jet2_pt > 0.);
        bool pass_SR          = isPromptCR(ev) && passVBFJets(ev) && ev.no_bjets;
        bool pass_CR_MET_low  = pass_CR && ev.MET <  50000.;
        bool pass_CR_MET_high = pass_CR && ev.MET >= 50000.;

        if (pass_CR)          n_after_CR++;
        if (pass_SR)          n_after_SR++;
        if (pass_CR_MET_low)  n_after_CR_MET_low++;
        if (pass_CR_MET_high) n_after_CR_MET_high++;

        // --- Prompt-electron + ID flags ---
        // Tighter IP requirement (|d0sig|<3, |z0sinθ|<0.3) on top of the
        // standard event-level cut.  Used to sharpen the loose denominator
        // in the prompt-electron sample so ε_r is unbiased by the
        // displaced-track tail.
        bool pass_IP                = passElectronIPTight(ev);
        bool electron_is_prompt     = isMCPromptEvent(ev);  // both legs prompt
        bool tight_electron_prompt  = ev.el_is_tight && electron_is_prompt;
        bool loose_electron_prompt  = ev.el_is_loose && electron_is_prompt;

        if (ev.el_is_tight) n_tight_electron_any++;
        if (ev.el_is_loose) n_loose_electron_any++;
        if (pass_SR) {
            if (tight_electron_prompt && pass_IP) n_tight_electron_prompt++;
            if (loose_electron_prompt)            n_loose_electron_prompt++;
        }

        // --- Per-electron variables ---
        double el_pt_val    = (el_pt_NOSYS && !el_pt_NOSYS->empty()) ? (*el_pt_NOSYS)[0] * 1e-3 : 0.;
        double el_eta_val   = ev.el_aeta;
        double el_d0_val    = (el_d0_NOSYS && !el_d0_NOSYS->empty()) ? std::fabs((*el_d0_NOSYS)[0]) : 0.;
        double el_d0sig_val = (ev.el_d0sig && !ev.el_d0sig->empty()) ? std::fabs((*ev.el_d0sig)[0]) : 0.;
        double met_val      = ev.MET * 1e-3;

#define FILL_REGION(region) \
        if (tight_electron_prompt && pass_IP) { \
            h_tight_pt_##region   ->Fill(el_pt_val,    w_tight); \
            h_tight_eta_##region  ->Fill(el_eta_val,   w_tight); \
            h_tight_d0_##region   ->Fill(el_d0_val,    w_tight); \
            h_tight_d0sig_##region->Fill(el_d0sig_val, w_tight); \
            h_tight_met_##region  ->Fill(met_val,      w_tight); \
            h2_tight_eta_pt_##region->Fill(el_pt_val, el_eta_val, w_tight); \
        } \
        if (loose_electron_prompt && pass_IP) { \
            h_loose_pt_##region   ->Fill(el_pt_val,    w_loose); \
            h_loose_eta_##region  ->Fill(el_eta_val,   w_loose); \
            h_loose_d0_##region   ->Fill(el_d0_val,    w_loose); \
            h_loose_d0sig_##region->Fill(el_d0sig_val, w_loose); \
            h_loose_met_##region  ->Fill(met_val,      w_loose); \
            h2_loose_eta_pt_##region->Fill(el_pt_val, el_eta_val, w_loose); \
        }

        if (pass_CR)          { FILL_REGION(CR) }
        if (pass_SR)          { FILL_REGION(SR) }
        if (pass_CR_MET_low)  { FILL_REGION(CR_MET_low) }
        if (pass_CR_MET_high) { FILL_REGION(CR_MET_high) }
    }

    // =============================================
    // Print selection summary
    // =============================================
    std::cout << "\n=======================================================\n"
              << "          REAL-EFFICIENCY EVENT SELECTION SUMMARY\n"
              << "=======================================================\n"
              << "Total entries in tree:                  " << nentries << "\n"
              << "After flavour (mu-e pair):              " << n_after_flavour << "\n"
              << "After trigger + matching:               " << n_after_trigger << "\n"
              << "After preselection:                     " << n_after_preselection << "\n"
              << "After CR (OS):                          " << n_after_CR << "\n"
              << "After SR (OS+VBF):                      " << n_after_SR << "\n"
              << "After CR_MET_low:                       " << n_after_CR_MET_low << "\n"
              << "After CR_MET_high:                      " << n_after_CR_MET_high << "\n"
              << "-------------------------------------------------------\n";

    std::cout << "Debug lepton counts (global):\n"
              << "  Tight electrons (any):                " << n_tight_electron_any    << "\n"
              << "  Loose electrons (any):                " << n_loose_electron_any    << "\n"
              << "  Tight prompt electrons (SR+IP):       " << n_tight_electron_prompt << "\n"
              << "  Loose prompt electrons (SR):          " << n_loose_electron_prompt << "\n"
              << "=======================================================\n\n";

    std::cout << "==== HIST SUMMARY ====\n"
              << "CR tight_pt: entries=" << h_tight_pt_CR->GetEntries() << " integral=" << h_tight_pt_CR->Integral() << "\n"
              << "SR tight_pt: entries=" << h_tight_pt_SR->GetEntries() << " integral=" << h_tight_pt_SR->Integral() << "\n"
              << "CR loose_pt: entries=" << h_loose_pt_CR->GetEntries() << " integral=" << h_loose_pt_CR->Integral() << "\n"
              << "SR loose_pt: entries=" << h_loose_pt_SR->GetEntries() << " integral=" << h_loose_pt_SR->Integral() << "\n";

    // =============================================
    // Save histograms
    // =============================================
    TFile* f_out = new TFile("outputs/real_eff_histograms.root", "RECREATE");
    if (!f_out || f_out->IsZombie()) {
        std::cerr << "ERROR: Could not create output file MC_histograms.root" << std::endl;
        return;
    }
    std::cout << "Saving histograms to MC_histograms.root ..." << std::endl;

#define WRITE_REGION(region) \
    h_tight_pt_##region->Write(); \
    h_loose_pt_##region->Write(); \
    h_tight_eta_##region->Write(); \
    h_loose_eta_##region->Write(); \
    h_tight_d0_##region->Write(); \
    h_loose_d0_##region->Write(); \
    h_tight_d0sig_##region->Write(); \
    h_loose_d0sig_##region->Write(); \
    h_tight_met_##region->Write(); \
    h_loose_met_##region->Write(); \
    h2_tight_eta_pt_##region->Write(); \
    h2_loose_eta_pt_##region->Write();

    WRITE_REGION(CR)
    WRITE_REGION(SR)
    WRITE_REGION(CR_MET_low)
    WRITE_REGION(CR_MET_high)

    f_out->Write();
    f_out->Close();
    std::cout << "All histograms saved successfully to MC_histograms.root" << std::endl;
}
