// =============================================================================
// FILE:       convert_fake_ntuple.C
// PURPOSE:    Convert the AMM-weighted data ntuple to the stat-framework
//             format.  Applies base preselection (`passBasePresel` =
//             topology + trigger fired + trigger matching + muon tight +
//             muon/electron IP cuts + crack veto + loose-pT thresholds)
//             and writes a new TTree with branches renamed and units
//             converted (MeV → GeV) to match the TRExFitter convention.
//             Output tree name = TMVATree to match the config NtupleName.
//
// TRIGGER MATCHING — DATA15 NOTE:
//   data15 ntuples lack the per-chain matching branches and use a single
//   global flag (globalTriggerMatch_NOSYS).  AMM_ElectronOnly.C writes a
//   per-event `isData15/O` branch into AMM/data_with_electron_fake_weights.root,
//   computed from the input filename.  SetupBranches() picks it up via
//   SafeSet, and passTriggerMatching() automatically dispatches to
//   `globalTriggerMatch` for data15 events (and to the per-chain matching
//   for everything else).  No special handling is needed in this macro.
//
// STAT-FRAMEWORK REGION SELECTIONS (applied by TRExFitter, not here):
//   SR      : Nbjet==0 && channel==2 && 0<lep_invMass<150 && mjj>400
//             && lep_1stPt>45 && lep_2ndPt>15 && jet_1stPt>35 && jet_2ndPt>20
//   VR      : same as SR but mjj<400
//   CR_tt   : same as SR but Nbjet>=1 (no b-veto)
//   CR_ztt  : Nbjet==0 && 35<lep_1stPt<45 && ... (Z→ττ enriched)
//
// INPUTS:     AMM/data_with_electron_fake_weights.root
// OUTPUTS:    fake_for_stat.root (tree: TMVATree)
// RUNS AFTER: AMM/AMM_ElectronOnly.C
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <cmath>
#include <iostream>

#include "utils/SelectionUtils.h"

void convert_fake_ntuple()
{
    // -------------------------------------------------------------------------
    // Input
    // -------------------------------------------------------------------------
    TChain in("tree");
    in.Add("AMM/data_with_electron_fake_weights.root");
    Long64_t nEntries = in.GetEntries();
    if (nEntries == 0) { std::cerr << "ERROR: empty input tree\n"; return; }
    std::cout << "Input entries: " << nEntries << "\n";

    // -------------------------------------------------------------------------
    // Branch setup via SelectionUtils EventData
    // -------------------------------------------------------------------------
    EventData ev;
    SetupBranches(&in, ev);

    // -------------------------------------------------------------------------
    // Output branches — stat-framework naming, units in GeV
    // -------------------------------------------------------------------------
    TFile* fout = TFile::Open("fake_for_stat.root","RECREATE");
    TTree* out  = new TTree("TMVATree","Fake/Non-Prompt background for stat framework");

    Float_t lep_invMass  = 0.f;   // mll [GeV]
    Float_t mjj_out      = 0.f;   // mjj [GeV]
    Float_t lep_1stPt    = 0.f;   // leading lepton pT [GeV]
    Float_t lep_2ndPt    = 0.f;   // subleading lepton pT [GeV]
    Float_t jet_1stPt    = 0.f;   // leading jet pT [GeV]
    Float_t jet_2ndPt    = 0.f;   // subleading jet pT [GeV]
    Int_t   Nbjet_out    = 0;     // number of b-tagged jets
    Int_t   channel_out  = 0;     // dilepton channel (always 2 = μe after presel)
    Int_t   lep_charge   = 0;     // l1_charge * l2_charge (<0 = OS, >0 = SS)

    Float_t FinalWeight           = 0.f;
    Float_t FinalWeight_MET_up    = 0.f;
    Float_t FinalWeight_MET_down  = 0.f;
    Float_t FinalWeight_comp_up   = 0.f;
    Float_t FinalWeight_comp_down = 0.f;
    Float_t FinalWeight_real_up   = 0.f;   // ±20 % prompt-MC normalisation
    Float_t FinalWeight_real_down = 0.f;
    Float_t FinalWeight_stat_up   = 0.f;
    Float_t FinalWeight_stat_down = 0.f;

    out->Branch("lep_invMass",           &lep_invMass);
    out->Branch("mjj",                   &mjj_out);
    out->Branch("lep_1stPt",             &lep_1stPt);
    out->Branch("lep_2ndPt",             &lep_2ndPt);
    out->Branch("jet_1stPt",             &jet_1stPt);
    out->Branch("jet_2ndPt",             &jet_2ndPt);
    out->Branch("Nbjet",                 &Nbjet_out);
    out->Branch("channel",               &channel_out);
    out->Branch("lep_charge",            &lep_charge);
    out->Branch("FinalWeight",           &FinalWeight);
    out->Branch("FinalWeight_MET_up",    &FinalWeight_MET_up);
    out->Branch("FinalWeight_MET_down",  &FinalWeight_MET_down);
    out->Branch("FinalWeight_comp_up",   &FinalWeight_comp_up);
    out->Branch("FinalWeight_comp_down", &FinalWeight_comp_down);
    out->Branch("FinalWeight_real_up",   &FinalWeight_real_up);
    out->Branch("FinalWeight_real_down", &FinalWeight_real_down);
    out->Branch("FinalWeight_stat_up",   &FinalWeight_stat_up);
    out->Branch("FinalWeight_stat_down", &FinalWeight_stat_down);

    // -------------------------------------------------------------------------
    // Event loop — apply base preselection, convert units, write
    // -------------------------------------------------------------------------
    Long64_t nPassed       = 0;
    Long64_t n_input_d15   = 0;   // bookkeeping: how many input events were data15
    Long64_t n_passed_d15  = 0;   // and how many passed preselection

    for (Long64_t i = 0; i < nEntries; ++i) {
        in.GetEntry(i);
        ComputeDerived(ev);

        if (ev.isData15) ++n_input_d15;

        // Base preselection: topology + trigger fired + matching + muon tight
        //                  + IP cuts + crack veto + loose pT (mu>20, e>=15 GeV)
        // For data15 events, passTriggerMatching() automatically uses
        // ev.globalTriggerMatch instead of the per-chain matching, since
        // ev.isData15 was loaded by SetupBranches() from the AMM-cloned tree.
        // For data16/17/18 and MC, the per-chain matching is used.
        // NOTE: region cuts (Nbjet, mjj, lep_1stPt>45, mll) are applied by
        //       TRExFitter using the branch values written below — do NOT pre-cut here.
        if (!passBasePresel(ev)) continue;

        // Δη_jj > 3 is a hard baseline cut (implicit VBF requirement)
        if (std::fabs(ev.deta_jj) <= 3.0) continue;

        if (ev.isData15) ++n_passed_d15;

        // Fill output branches (MeV → GeV)
        lep_invMass    = static_cast<float>(ev.mll    * 1e-3);
        mjj_out        = static_cast<float>(ev.mjj    * 1e-3);
        lep_1stPt      = ev.l1_pt * 1e-3f;
        lep_2ndPt      = ev.l2_pt * 1e-3f;
        jet_1stPt      = static_cast<float>(ev.jet1_pt * 1e-3);
        jet_2ndPt      = static_cast<float>(ev.jet2_pt * 1e-3);
        Nbjet_out      = ev.Nbjet;
        channel_out    = 2;   // always μe after topology cut in passBasePresel
        lep_charge     = ev.l1_charge * ev.l2_charge;

        FinalWeight           = static_cast<float>(safe(ev.fw_nom));
        FinalWeight_MET_up    = static_cast<float>(safe(ev.fw_met_up));
        FinalWeight_MET_down  = static_cast<float>(safe(ev.fw_met_dn));
        FinalWeight_comp_up   = static_cast<float>(safe(ev.fw_comp_up));
        FinalWeight_comp_down = static_cast<float>(safe(ev.fw_comp_dn));
        FinalWeight_real_up   = static_cast<float>(safe(ev.fw_real_up));
        FinalWeight_real_down = static_cast<float>(safe(ev.fw_real_dn));
        FinalWeight_stat_up   = static_cast<float>(safe(ev.fw_stat_up));
        FinalWeight_stat_down = static_cast<float>(safe(ev.fw_stat_dn));

        out->Fill();
        ++nPassed;
    }

    fout->Write();
    fout->Close();

    std::cout << "Passed preselection : " << nPassed << " / " << nEntries << "\n";
    std::cout << "Input  data15 events : " << n_input_d15
              << "    (using globalTriggerMatch)\n";
    std::cout << "Passed data15 events : " << n_passed_d15
              << "    (data15 contribution to fake estimate)\n";
    if (n_input_d15 == 0) {
        std::cout << "  WARNING: no events tagged as data15 — either the AMM "
                     "ntuple wasn't re-run after the data15-flag changes, "
                     "or data15.root wasn't in the AMM input list.\n";
    }
    std::cout << "Output: fake_for_stat.root (tree: TMVATree)\n\n";
    std::cout << "Branch mapping:\n";
    std::cout << "  lep_invMass  = mll [GeV]    | mjj         = mjj [GeV]\n";
    std::cout << "  lep_1stPt    = l1 pT [GeV]  | lep_2ndPt   = l2 pT [GeV]\n";
    std::cout << "  jet_1stPt    = j1 pT [GeV]  | jet_2ndPt   = j2 pT [GeV]\n";
    std::cout << "  Nbjet        = b-jet count  | channel     = 2 (μe)\n";
    std::cout << "  FinalWeight  = AMM weight (nominal)\n";
    std::cout << "  FinalWeight_{MET,comp,real,stat}_{up,down} = systematics\n";
    std::cout << "  Preselection: topology + trigger fired + trigger matching "
                 "(per-chain or globalTriggerMatch for data15)\n";
    std::cout << "                + muon tight + IP cuts + crack veto + |Δη_jj|>3\n";
    std::cout << "  DO NOT cut on el tight/loose — FinalWeight encodes that via AMM.\n";
}
