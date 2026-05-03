// =============================================================================
// utils/SelectionUtils.h
// Central selection library for the HNL μe VBF fake background analysis.
//
// USAGE in any macro:
//   #include "../utils/SelectionUtils.h"   (adjust path as needed)
//
//   EventData ev;
//   SetupBranches(&chain, ev);      // once, before event loop
//   for (Long64_t i = 0; i < n; ++i) {
//       chain.GetEntry(i);
//       ComputeDerived(ev);
//       if (!passBasePresel(ev)) continue;
//       if (isSR(ev)) { /* fill SR */ }
//   }
//
// All pT / MET / mll / mjj values in the ntuple are in MeV.
// Region cut comparisons convert to GeV internally where the stat-framework
// config specifies GeV numbers.
// =============================================================================
#pragma once
#include "ROOT/RVec.hxx"
#include <TTree.h>
#include <cmath>

using RVecC = ROOT::VecOps::RVec<char>;
using RVecF = ROOT::VecOps::RVec<float>;
using RVecI = ROOT::VecOps::RVec<int>;

// ---------------------------------------------------------------------------
// pT / eta bin edges — keep in sync with fake_bg_estimation.C / real_eff.C
// ---------------------------------------------------------------------------
// NOTE: top pT bin merged ([50,80] ⊕ [80,200] → [50,200]) to absorb the
// low-stat / unphysically-high-rate behaviour we saw in the original
// [80,200] GeV bin (~30 % central value, ~30 % stat).  Merging widens the
// last bin's denominator, stabilises the rate, and removes the need for
// the flat 100 % top-bin treatment that build_fake_systematics.C used to
// inject.  Keep kPtBins in sync with any plotting-side bin definitions.
constexpr int    kNPtBins  = 5;
constexpr int    kNEtaBins = 5;
constexpr double kPtBins [kNPtBins  + 1] = {5, 15, 25, 35, 50, 200};     // GeV
constexpr double kEtaBins[kNEtaBins + 1] = {0.0, 0.7, 1.4, 1.6, 2.0, 2.5};

// ---------------------------------------------------------------------------
// safe() — discard events with non-finite generator weights
// ---------------------------------------------------------------------------
inline double safe(double x) { return std::isfinite(x) ? x : 0.0; }

// =============================================================================
// EventData — holds every branch pointer for one event.
// Declare one instance per tree, call SetupBranches once, then GetEntry + ComputeDerived.
// =============================================================================
struct EventData {
    // ---- Flavor / charge ----
    Int_t   l1_flavor = 0, l2_flavor = 0;
    Int_t   l1_charge = 0, l2_charge = 0;

    // ---- Kinematics (MeV, as stored in ntuple) ----
    Float_t  l1_pt = 0.f, l2_pt = 0.f, l2_eta = 0.f;
    Double_t mll = 0., mjj = 0., deta_jj = 0.;
    Double_t jet1_pt = 0., jet2_pt = 0., MET = 0.;

    // ---- MC event weights (absent in data — default to 1 / 1 / 1) ----
    Double_t total_weight       = 1.;
    Float_t  globalTriggerEffSF = 1.f;
    Float_t  weight_lepSF_loose = 1.f;
    Float_t  weight_lepSF_tight = 1.f;

    // ---- Lepton / jet ID (vector branches — nullable) ----
    RVecC* mu_loose      = nullptr;  // mu_select_loose_NOSYS  (for anti-iso CR studies)
    RVecC* mu_tight      = nullptr;  // mu_select_tight_NOSYS
    RVecC* el_loose      = nullptr;  // el_select_loose_NOSYS
    RVecC* el_tight      = nullptr;  // el_select_tight_NOSYS
    RVecC* jet_btag      = nullptr;  // jet_GN2v01_FixedCutBEff_85_select

    // ---- Impact parameters ----
    RVecF* mu_d0sig      = nullptr;
    RVecF* el_d0sig      = nullptr;
    RVecF* mu_z0sintheta = nullptr;
    RVecF* el_z0sintheta = nullptr;
    RVecF* el_eta_vec    = nullptr;

    // ---- MC truth prompt flags (absent in data) ----
    RVecI* mu_truth_isPrompt = nullptr;
    RVecI* el_truth_isPrompt = nullptr;

    // ---- Per-chain trigger matching (vector branches — absent in data15
    //      and possibly in some MC; SafeSet leaves them at nullptr).
    //      passTriggerMatching() ORs these with globalTriggerMatch_NOSYS,
    //      so a missing branch never makes the matching fail spuriously. ----
    RVecC* el_match_e26  = nullptr;
    RVecC* el_match_e60  = nullptr;
    RVecC* el_match_e140 = nullptr;
    RVecC* mu_match_mu26 = nullptr;
    RVecC* mu_match_mu50 = nullptr;

    // ---- Trigger fired (every chain across 2015 / 2016+ / 2022 — any
    //      branch that isn't in a given ntuple stays at default `false`
    //      thanks to SafeSet, so the OR below works for all years) ----
    Bool_t trig_e24            = false;   // 2015
    Bool_t trig_e26_nod0       = false;   // 2016+
    Bool_t trig_e26_22         = false;   // 2022
    Bool_t trig_e60            = false;   // 2015
    Bool_t trig_e60_nod0       = false;   // 2016+
    Bool_t trig_e60_22         = false;   // 2022
    Bool_t trig_e120           = false;   // 2015
    Bool_t trig_e140_nod0      = false;   // 2016+
    Bool_t trig_e140_22        = false;   // 2022
    Bool_t trig_mu20           = false;   // 2015
    Bool_t trig_mu24           = false;   // 2022
    Bool_t trig_mu26           = false;   // 2016+
    Bool_t trig_mu40           = false;   // 2015
    Bool_t trig_mu50           = false;   // 2016+
    Bool_t trig_mu50_22        = false;   // 2022

    // Aliases kept for any code that still references the old names —
    // bound below to the same branches as the matching new fields.
    Bool_t trig_e26  = false;             // alias of trig_e26_nod0
    Bool_t trig_e140 = false;             // alias of trig_e140_nod0

    // ---- Fake weights (attached by AMM_ElectronOnly.C — absent in raw MC/data) ----
    Double_t fw_nom     = 0., fw_met_up  = 0., fw_met_dn  = 0.;
    Double_t fw_comp_up = 0., fw_comp_dn = 0.;
    Double_t fw_real_up = 0., fw_real_dn = 0.;
    Double_t fw_stat_up = 0., fw_stat_dn = 0.;

    // ---- Derived quantities (filled by ComputeDerived, call after GetEntry) ----
    int   Nbjet      = 0;
    bool  no_bjets   = true;
    bool  is_os      = false;   // opposite-sign
    bool  mu_is_loose   = false;
    bool  mu_is_tight   = false;
    bool  el_is_loose   = false;
    bool  el_is_tight   = false;
    float el_aeta       = 999.f;  // |η| of leading electron
    
    // ---- Dataset flags (set from macro) ----
    bool isData   = false;
    bool isData15 = false;

    // ---- Global trigger matching (data15 only) ----
    Bool_t globalTriggerMatch = false;
    
    
    
};

// ---------------------------------------------------------------------------
// SafeSet — only call SetBranchAddress if the branch exists in the tree.
// Required for optional branches (MC weights, truth, trigger matching,
// fake weights, BDT score) so the same SetupBranches works for all tree types.
// ---------------------------------------------------------------------------
inline void SafeSet(TTree* t, const char* name, void* addr) {
    if (t->GetBranch(name)) t->SetBranchAddress(name, addr);
}

// ---------------------------------------------------------------------------
// SetupBranches — connect all EventData members to a TTree/TChain.
// Call once per tree before the event loop; re-call after switching trees.
// ---------------------------------------------------------------------------
inline void SetupBranches(TTree* t, EventData& ev) {
    // Core kinematics (always present)
    // All bindings via SafeSet — missing branches simply leave the field at
    // its EventData default value, instead of segfaulting in GetEntry.
    SafeSet(t, "l1_flavor",  &ev.l1_flavor);
    SafeSet(t, "l2_flavor",  &ev.l2_flavor);
    SafeSet(t, "l1_charge",  &ev.l1_charge);
    SafeSet(t, "l2_charge",  &ev.l2_charge);
    SafeSet(t, "l1_pt",      &ev.l1_pt);
    SafeSet(t, "l2_pt",      &ev.l2_pt);
    SafeSet(t, "l2_eta",     &ev.l2_eta);
    SafeSet(t, "mll",        &ev.mll);
    SafeSet(t, "mjj",        &ev.mjj);
    SafeSet(t, "deta_jj",    &ev.deta_jj);
    SafeSet(t, "jet1_pt",    &ev.jet1_pt);
    SafeSet(t, "jet2_pt",    &ev.jet2_pt);
    SafeSet(t, "MET", &ev.MET);

    // MC event weights (absent in data)
    SafeSet(t, "total_weight",                &ev.total_weight);
    SafeSet(t, "globalTriggerEffSF_NOSYS",    &ev.globalTriggerEffSF);
    SafeSet(t, "weight_leptonSF_loose_NOSYS", &ev.weight_lepSF_loose);
    SafeSet(t, "weight_leptonSF_tight_NOSYS", &ev.weight_lepSF_tight);
    SafeSet(t, "globalTriggerMatch_NOSYS", &ev.globalTriggerMatch);

    // isData15 — present only on AMM-cloned data ntuples
    // (AMM_ElectronOnly.C writes it from the input filename).
    // For raw MC / raw data the branch is missing and the field stays at
    // its default (false); macros looping raw data must still set it
    // manually from the filename.
    SafeSet(t, "isData15", &ev.isData15);

    // Lepton / jet ID
    SafeSet(t, "mu_select_loose_NOSYS",              &ev.mu_loose);
    SafeSet(t, "mu_select_tight_NOSYS",              &ev.mu_tight);
    SafeSet(t, "el_select_loose_NOSYS",              &ev.el_loose);
    SafeSet(t, "el_select_tight_NOSYS",              &ev.el_tight);
    SafeSet(t, "jet_GN2v01_FixedCutBEff_85_select",  &ev.jet_btag);

    // Impact parameters
    SafeSet(t, "mu_d0sig_NOSYS",      &ev.mu_d0sig);
    SafeSet(t, "el_d0sig_NOSYS",      &ev.el_d0sig);
    SafeSet(t, "mu_z0sintheta_NOSYS", &ev.mu_z0sintheta);
    SafeSet(t, "el_z0sintheta_NOSYS", &ev.el_z0sintheta);
    SafeSet(t, "el_eta",              &ev.el_eta_vec);

    // MC truth (absent in data)
    SafeSet(t, "mu_truth_isPrompt", &ev.mu_truth_isPrompt);
    SafeSet(t, "el_truth_isPrompt", &ev.el_truth_isPrompt);

    // Triggers fired
    // ----- Trigger fired -----------------------------------------------
    // Use SafeSet so a missing chain (e.g. data15 lacks 2016 names; data16+
    // lacks the 2015 e24 chain) just leaves the Bool_t at default `false`.
    // 2015 chains
    SafeSet(t, "trigPassed_HLT_e24_lhmedium_L1EM20VH",        &ev.trig_e24);
    SafeSet(t, "trigPassed_HLT_e60_lhmedium",                 &ev.trig_e60);
    SafeSet(t, "trigPassed_HLT_e120_lhloose",                 &ev.trig_e120);
    SafeSet(t, "trigPassed_HLT_mu20_iloose_L1MU15",           &ev.trig_mu20);
    SafeSet(t, "trigPassed_HLT_mu40",                         &ev.trig_mu40);
    // 2016+ chains (also bound to the legacy aliases trig_e26/trig_e140)
    SafeSet(t, "trigPassed_HLT_e26_lhtight_nod0_ivarloose",   &ev.trig_e26_nod0);
    SafeSet(t, "trigPassed_HLT_e26_lhtight_nod0_ivarloose",   &ev.trig_e26);
    SafeSet(t, "trigPassed_HLT_e60_lhmedium_nod0",            &ev.trig_e60_nod0);
    SafeSet(t, "trigPassed_HLT_e140_lhloose_nod0",            &ev.trig_e140_nod0);
    SafeSet(t, "trigPassed_HLT_e140_lhloose_nod0",            &ev.trig_e140);
    SafeSet(t, "trigPassed_HLT_mu26_ivarmedium",              &ev.trig_mu26);
    SafeSet(t, "trigPassed_HLT_mu50",                         &ev.trig_mu50);
    // 2022 chains
    SafeSet(t, "trigPassed_HLT_e26_lhtight_ivarloose_L1EM22VHI", &ev.trig_e26_22);
    SafeSet(t, "trigPassed_HLT_e60_lhmedium_L1EM22VHI",          &ev.trig_e60_22);
    SafeSet(t, "trigPassed_HLT_e140_lhloose_L1EM22VHI",          &ev.trig_e140_22);
    SafeSet(t, "trigPassed_HLT_mu24_ivarmedium_L1MU14FCH",       &ev.trig_mu24);
    SafeSet(t, "trigPassed_HLT_mu50_L1MU14FCH",                  &ev.trig_mu50_22);

    // Trigger matching (absent in data ntuples)
    SafeSet(t, "el_trigMatched_HLT_e26_lhtight_nod0_ivarloose", &ev.el_match_e26);
    SafeSet(t, "el_trigMatched_HLT_e60_lhmedium_nod0",           &ev.el_match_e60);
    SafeSet(t, "el_trigMatched_HLT_e140_lhloose_nod0",           &ev.el_match_e140);
    SafeSet(t, "mu_trigMatched_HLT_mu26_ivarmedium",              &ev.mu_match_mu26);
    SafeSet(t, "mu_trigMatched_HLT_mu50",                         &ev.mu_match_mu50);

    // Fake weights (only in AMM-weighted ntuples)
    SafeSet(t, "fake_weight_nominal",   &ev.fw_nom);
    SafeSet(t, "fake_weight_MET_up",    &ev.fw_met_up);
    SafeSet(t, "fake_weight_MET_down",  &ev.fw_met_dn);
    SafeSet(t, "fake_weight_comp_up",   &ev.fw_comp_up);
    SafeSet(t, "fake_weight_comp_down", &ev.fw_comp_dn);
    SafeSet(t, "fake_weight_real_up",   &ev.fw_real_up);
    SafeSet(t, "fake_weight_real_down", &ev.fw_real_dn);
    SafeSet(t, "fake_weight_statUp",    &ev.fw_stat_up);
    SafeSet(t, "fake_weight_statDown",  &ev.fw_stat_dn);

}

// ResetTrigMatchPointers — call before switching to a new tree so that
// dangling pointers from the old tree don't cause crashes in SafeSet-guarded branches.
inline void ResetTrigMatchPointers(EventData& ev) {
    ev.el_match_e26  = nullptr;
    ev.el_match_e60  = nullptr;
    ev.el_match_e140 = nullptr;
    ev.mu_match_mu26 = nullptr;
    ev.mu_match_mu50 = nullptr;
}

// ---------------------------------------------------------------------------
// ComputeDerived — fill derived fields after each GetEntry().
// Must be called before any region function.
// ---------------------------------------------------------------------------
inline void ComputeDerived(EventData& ev) {
    // b-jet count
    ev.Nbjet    = 0;
    ev.no_bjets = true;
    if (ev.jet_btag && !ev.jet_btag->empty()) {
        for (char b : *ev.jet_btag) {
            if (b) { ++ev.Nbjet; ev.no_bjets = false; }
        }
    }
    // Charge
    ev.is_os = (ev.l1_charge * ev.l2_charge < 0);
    // Lepton ID
    ev.mu_is_loose = ev.mu_loose && !ev.mu_loose->empty() && (*ev.mu_loose)[0];
    ev.mu_is_tight = ev.mu_tight && !ev.mu_tight->empty() && (*ev.mu_tight)[0];
    ev.el_is_loose = ev.el_loose && !ev.el_loose->empty() && (*ev.el_loose)[0];
    ev.el_is_tight = ev.el_tight && !ev.el_tight->empty() && (*ev.el_tight)[0];
    // Electron |η|
    ev.el_aeta = (ev.el_eta_vec && !ev.el_eta_vec->empty()) ?
                 std::fabs((*ev.el_eta_vec)[0]) : 999.f;
}

// =============================================================================
// Individual selection functions — each takes a fully-populated EventData.
// All pT values in MeV.
// =============================================================================

// ---- Topology ----
inline bool passTopology(const EventData& ev) {
    return ev.l1_flavor == 1 && ev.l2_flavor == 0;  // l1=muon, l2=electron
}

// ---- Trigger ----
inline bool passMuonTrigger(const EventData& ev) {
    return ev.trig_mu20    // 2015
        || ev.trig_mu40    // 2015
        || ev.trig_mu26    // 2016+
        || ev.trig_mu50    // 2016+
        || ev.trig_mu24    // 2022
        || ev.trig_mu50_22;// 2022
}
inline bool passElectronTrigger(const EventData& ev) {
    return ev.trig_e24                                            // 2015
        || ev.trig_e60   || ev.trig_e120                         // 2015
        || ev.trig_e26   || ev.trig_e26_nod0  || ev.trig_e26_22  // 2016+ / 2022
        || ev.trig_e60_nod0 || ev.trig_e60_22
        || ev.trig_e140  || ev.trig_e140_nod0 || ev.trig_e140_22;
}
// ---- Trigger logic ----
// Two distinct policies, picked by the calling region:
//
//  passTrigger:
//      Used by SR / VR / CR_tt / CR_ztautau — analysis triggers a muon OR
//      an electron. Maximum acceptance.
//
//  passMuonTriggerOnly:
//      Used by the fake CR (ε_f measurement) and prompt CR (ε_r measurement).
//      The electron is the unbiased probe and must NOT have fired the trigger,
//      otherwise the loose population is HLT-biased toward tight-passing.
//      References: HNL multilepton arXiv:2204.11138, SUSY 2-l SS arXiv:1909.08457.
inline bool passTrigger(const EventData& ev) {
    return passMuonTrigger(ev) || passElectronTrigger(ev);
}
inline bool passMuonTriggerOnly(const EventData& ev) {
    return passMuonTrigger(ev);
}

// ---- Trigger matching ----
// Same asymmetry as the trigger pass: full matching for SR/VR/CR_tt/CR_ztautau,
// muon-only matching for fake CR / prompt CR.
//
// All routes work without a year-specific dispatch — branches that aren't
// in the ntuple stay nullptr and are silently skipped:
//   1. globalTriggerMatch_NOSYS  (universal; ONLY route in data15)
//   2. mu_trigMatched_HLT_*       (per-muon flags; data16+ / MC)
//   3. el_trigMatched_HLT_*       (per-electron flags; data16+ / MC)

inline bool passTriggerMatching(const EventData& ev) {
    if (ev.globalTriggerMatch) return true;
    auto matched = [](const RVecC* p) {
        return p && !p->empty() && (*p)[0];
    };
    // muon side
    if (ev.trig_mu26 && matched(ev.mu_match_mu26)) return true;
    if (ev.trig_mu50 && matched(ev.mu_match_mu50)) return true;
    // electron side (only relevant for SR/VR/CR_tt/CR_ztautau)
    if (ev.trig_e26  && matched(ev.el_match_e26))  return true;
    if (ev.trig_e60  && matched(ev.el_match_e60))  return true;
    if (ev.trig_e140 && matched(ev.el_match_e140)) return true;
    return false;
}

// Muon-only matching, for the fake-rate measurement only.
inline bool passMuonTriggerMatchingOnly(const EventData& ev) {
    if (ev.globalTriggerMatch) return true;
    auto matched = [](const RVecC* p) {
        return p && !p->empty() && (*p)[0];
    };
    if (ev.trig_mu26 && matched(ev.mu_match_mu26)) return true;
    if (ev.trig_mu50 && matched(ev.mu_match_mu50)) return true;
    return false;
}

// ---- Impact parameters ----
inline bool passMuonIP(const EventData& ev) {
    return ev.mu_d0sig      && !ev.mu_d0sig->empty()      && std::fabs((*ev.mu_d0sig)[0])      < 3.0f
        && ev.mu_z0sintheta && !ev.mu_z0sintheta->empty() && std::fabs((*ev.mu_z0sintheta)[0]) < 0.5f;
}
inline bool passElectronIP(const EventData& ev) {
    return ev.el_d0sig      && !ev.el_d0sig->empty()      && std::fabs((*ev.el_d0sig)[0])      < 5.0f
        && ev.el_z0sintheta && !ev.el_z0sintheta->empty() && std::fabs((*ev.el_z0sintheta)[0]) < 0.5f;
}
// Tighter electron-IP cut (|d0sig|<3, |z0sinθ|<0.3) used in MC_eff.C and
// real_eff.C only.  Purpose: enrich the loose denominator with prompt-like
// electrons before computing ε_r (and ε_f).  Removes the long-lived /
// displaced tail that dilutes the prompt sample, raising ε_r toward the
// ~0.85–0.95 regime quoted by ATLAS.  Do NOT use this on the AMM
// application path — the rate maps and the application must use compatible
// IP requirements; the standard `passElectronIP` is the analysis baseline.
inline bool passElectronIPTight(const EventData& ev) {
    return ev.el_d0sig      && !ev.el_d0sig->empty()      && std::fabs((*ev.el_d0sig)[0])      < 3.0f
        && ev.el_z0sintheta && !ev.el_z0sintheta->empty() && std::fabs((*ev.el_z0sintheta)[0]) < 0.3f;
}

// ---- Electron crack veto (|η| in [1.37, 1.52]) ----
inline bool passCrackVeto(const EventData& ev) {
    return !(ev.el_aeta > 1.37f && ev.el_aeta < 1.52f);
}

// ---- Lepton pT ----
// Loose (fake CR / prompt CR): mu > 26 GeV (trigger-level), e >= 5 GeV.
// The lowered e>=5 cut lets the [5,15] GeV bin in kPtBins populate. The
// SR/VR/CR_tt/CR_ztautau add `passSignalLeptonPt` on top, which restores
// the e>15 GeV requirement, so signal selections are unaffected.
inline bool passLeptonPt(const EventData& ev) {
    return ev.l1_pt > 26000.f && ev.l2_pt >= 5000.f;
}
// Signal (ttbar CR / Ztt CR / VR / SR): mu > 45 GeV, e > 15 GeV
inline bool passSignalLeptonPt(const EventData& ev) {
    return ev.l1_pt > 45000.f && ev.l2_pt > 15000.f;
}

// ---- VBF jet requirements: j1 > 35 GeV, j2 > 20 GeV, |Δη| > 3 ----
inline bool passVBFJets(const EventData& ev) {
    return ev.jet1_pt > 35000. && ev.jet2_pt > 20000. && std::fabs(ev.deta_jj) > 3.0;
}

// ---- MC truth ----
inline bool isMCPromptMuon(const EventData& ev) {
    return ev.mu_truth_isPrompt && !ev.mu_truth_isPrompt->empty() && (*ev.mu_truth_isPrompt)[0] == 1;
}
inline bool isMCPromptElectron(const EventData& ev) {
    return ev.el_truth_isPrompt && !ev.el_truth_isPrompt->empty() && (*ev.el_truth_isPrompt)[0] == 1;
}
inline bool isMCFakeElectron(const EventData& ev) { return !isMCPromptElectron(ev); }
inline bool isMCPromptEvent(const EventData& ev)  { return isMCPromptMuon(ev) && isMCPromptElectron(ev); }

// ---- MC event weight helpers ----
inline double mcWeightTight(const EventData& ev) {
    double w = safe(ev.total_weight) * safe(ev.globalTriggerEffSF) * safe(ev.weight_lepSF_tight);
    return std::isfinite(w) ? w : 0.0;
}
inline double mcWeightLoose(const EventData& ev) {
    double w = safe(ev.total_weight) * safe(ev.globalTriggerEffSF) * safe(ev.weight_lepSF_loose);
    return std::isfinite(w) ? w : 0.0;
}

// =============================================================================
// Base preselection — common to all downstream selections.
// Does NOT include electron ID (tight/loose) — callers add that explicitly.
// =============================================================================
inline bool passBasePresel(const EventData& ev) {
    return passTopology(ev)
        && passTrigger(ev)
        && passTriggerMatching(ev)
        && ev.mu_is_tight
        && passMuonIP(ev)
        && passElectronIP(ev)
        && passCrackVeto(ev)
        && passLeptonPt(ev);   // loose: mu>20 GeV, e>=15 GeV
}

// passBasePresel_FakeMeasurement — used in MC_eff.C / data_eff.C / real_eff.C
// and fake-CR studies. Identical to passBasePresel but with MUON-ONLY trigger
// and matching, so the electron is unbiased ("probe").
//
// SR / VR / CR_tt / CR_ztautau use passBasePresel directly (mu OR e trigger).
inline bool passBasePresel_FakeMeasurement(const EventData& ev) {
    return passTopology(ev)
        && passMuonTriggerOnly(ev)
        && passMuonTriggerMatchingOnly(ev)
        && ev.mu_is_tight
        && passMuonIP(ev)
        && passElectronIP(ev)
        && passCrackVeto(ev)
        && passLeptonPt(ev);
}

// =============================================================================
// passBasePreselNoElectronIP — same as passBasePresel but WITHOUT the electron
// d0sig / z0sinθ requirement.  Used by MC_eff.C and real_eff.C when measuring
// ε_f and ε_r so that the loose denominator and the tight numerator are taken
// from the natural loose population (not artificially trimmed by IP cuts).
//
// Rationale: applying the electron-IP cut on the loose denominator removes
// part of the population the AMM is supposed to extrapolate, biasing ε_r low.
// By dropping it here we recover the "raw" ε_r that ATLAS HNL/SUSY analyses
// quote.  Muon IP is kept (the muon side is required tight everywhere).
// =============================================================================
inline bool passBasePreselNoElectronIP(const EventData& ev) {
    return passTopology(ev)
        && passTrigger(ev)
        && passTriggerMatching(ev)
        && ev.mu_is_tight
        && passMuonIP(ev)
        && passCrackVeto(ev)
        && passLeptonPt(ev);
}

// =============================================================================
// Region functions — each encodes the EXACT stat-framework selection.
// All mll / mjj comparisons are in GeV (divided from MeV inside each function).
// =============================================================================

// ---- Fake CR ----
// Same-sign (SS) μe.  pT: passLeptonPt (mu>20 GeV, e>=15 GeV).
// Used for fake-rate measurement (ε_f) and as AMM application region.
inline bool isFakeCR(const EventData& ev) {
    return !ev.is_os;   // SS; pT / ID enforced by passBasePresel
}

// ---- Prompt CR ----
// OS .  Same loose pT as fake CR (mu>20 GeV, e>=15 GeV).
// Used for MC prompt subtraction inside the SS CR fake-rate measurement.
inline bool isPromptCR(const EventData& ev) {
    return ev.is_os;
}

// ---- ttbar CR (HNL_CR_tt) ----
// OS + Nbjet>=1 + mu>45 + e>15 + j1>35 + j2>20 + mjj>400 + 0<mll<150  [GeV]
inline bool isTtbarCR(const EventData& ev) {
    if (!ev.is_os || ev.no_bjets) return false;
    if (!passSignalLeptonPt(ev)) return false;
    if (!passVBFJets(ev)) return false;
    double mjj_gev = ev.mjj * 1e-3, mll_gev = ev.mll * 1e-3;
    return mjj_gev > 400. && mll_gev > 0. && mll_gev < 150.;
}

// ---- Z→ττ CR (HNL_CR_ztautau) ----
// OS + Nbjet==0 + 35<mu<45 + e>15 + j1>35 + j2>20 + mjj>400 + 0<mll<150  [GeV]
inline bool isZtautauCR(const EventData& ev) {
    if (!ev.is_os || !ev.no_bjets) return false;
    float l1_gev = ev.l1_pt * 1e-3f;
    if (l1_gev <= 35.f || l1_gev >= 45.f) return false;   // 35 < l1_pt < 45 GeV
    if (ev.l2_pt <= 15000.f) return false;
    if (!passVBFJets(ev)) return false;
    double mjj_gev = ev.mjj * 1e-3, mll_gev = ev.mll * 1e-3;
    return mjj_gev > 400. && mll_gev > 0. && mll_gev < 150.;
}

// ---- Validation Region (HNL_VR) ----
// OS + Nbjet==0 + mu>45 + e>15 + j1>35 + j2>20 + mjj<400 + 0<mll<150  [GeV]
inline bool isVR(const EventData& ev) {
    if (!ev.is_os || !ev.no_bjets) return false;
    if (!passSignalLeptonPt(ev)) return false;
    if (ev.jet1_pt <= 35000. || ev.jet2_pt <= 20000. || std::fabs(ev.deta_jj) <= 3.0) return false;
    double mjj_gev = ev.mjj * 1e-3, mll_gev = ev.mll * 1e-3;
    return mjj_gev > 0. && mjj_gev < 400. && mll_gev > 0. && mll_gev < 150.;
}

// ---- Signal Region (HNL_SR) ----
// OS + Nbjet==0 + mu>45 + e>15 + j1>35 + j2>20 + mjj>400 + 0<mll<150  [GeV]
inline bool isSR(const EventData& ev) {
    if (!ev.is_os || !ev.no_bjets) return false;
    if (!passSignalLeptonPt(ev)) return false;
    if (!passVBFJets(ev)) return false;
    double mjj_gev = ev.mjj * 1e-3, mll_gev = ev.mll * 1e-3;
    return mjj_gev > 400. && mll_gev > 0. && mll_gev < 150.;
}
