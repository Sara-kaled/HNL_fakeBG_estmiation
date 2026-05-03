// =============================================================================
// FILE:    studies/fake_cutflow.C
// PURPOSE: Fake-background cutflow — applies the SAME analysis cuts as the
//          SR/VR (full mu OR e trigger), but reads the AMM-weighted ntuple
//          (data_with_electron_fake_weights.root) and uses fake_weight_nominal
//          as the per-event weight. The final yield is the predicted fake
//          contribution at each cut step.
//
// NOTE: The fake-rate measurement itself (in MC_eff/data_eff/real_eff) uses
//       muon-only triggering. Here we APPLY the AMM weight in the analysis
//       region — so the trigger asymmetry is "muon-only at measurement,
//       full mu||e at application". That mismatch is the standard ATLAS
//       practice (HNL multilepton arXiv:2204.11138, SUSY 2-l SS arXiv:1909.08457).
// =============================================================================
#include <TChain.h>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include "../utils/SelectionUtils.h"

void fake_cutflow()
{
    TChain* tree = new TChain("tree");
    tree->Add("../AMM/data_with_electron_fake_weights.root");
    if (tree->GetEntries() == 0) {
        std::cerr << "ERROR: AMM ntuple not found — run AMM/AMM_ElectronOnly.C first.\n";
        return;
    }

    EventData ev;
    SetupBranches(tree, ev);

    // The AMM weight branch isn't part of EventData; read it directly.
    Double_t fake_weight_nominal = 0.;
    tree->SetBranchAddress("fake_weight_nominal", &fake_weight_nominal);

    const std::vector<std::string> cuts = {
        "Initial",
        "mu-e pair",
        "Trigger + match (mu||e)",
        "Tight muon",
        "Muon IP",
        "Electron IP",
        "Crack veto",
        "Signal lepton pT (mu>45, e>15)",
        "Opposite-sign",
        "b-veto",
        "VBF jets (j1>35, j2>20, |Deta|>3)",
        "0 < mll < 150 GeV",
        "VR  (mjj < 400 GeV)",
        "SR  (mjj > 400 GeV)"
    };
    std::map<std::string, double> Y;
    for (const auto& c : cuts) Y[c] = 0.;

    Long64_t n = tree->GetEntries();
    std::cout << "AMM ntuple entries: " << n << "\n";
    for (Long64_t i = 0; i < n; ++i) {
        if (i % 1000000 == 0) std::cout << "  " << i << " / " << n << "\n";
        tree->GetEntry(i);
        ComputeDerived(ev);

        const double w = fake_weight_nominal;
        Y["Initial"] += w;

        if (!passTopology(ev))                                            continue;  Y["mu-e pair"] += w;
        if (!passTrigger(ev) || !passTriggerMatching(ev))                 continue;  Y["Trigger + match (mu||e)"] += w;
        if (!ev.mu_is_tight)                                              continue;  Y["Tight muon"] += w;
        if (!passMuonIP(ev))                                              continue;  Y["Muon IP"] += w;
        if (!passElectronIP(ev))                                          continue;  Y["Electron IP"] += w;
        if (!passCrackVeto(ev))                                           continue;  Y["Crack veto"] += w;
        if (!passSignalLeptonPt(ev))                                      continue;  Y["Signal lepton pT (mu>45, e>15)"] += w;
        if (!ev.is_os)                                                    continue;  Y["Opposite-sign"] += w;
        if (!ev.no_bjets)                                                 continue;  Y["b-veto"] += w;
        if (!passVBFJets(ev))                                             continue;  Y["VBF jets (j1>35, j2>20, |Deta|>3)"] += w;

        const double mll_gev = ev.mll * 1e-3;
        if (!(mll_gev > 0. && mll_gev < 150.))                            continue;  Y["0 < mll < 150 GeV"] += w;

        const double mjj_gev = ev.mjj * 1e-3;
        if (mjj_gev < 400.) Y["VR  (mjj < 400 GeV)"] += w;
        else                Y["SR  (mjj > 400 GeV)"] += w;
    }

    std::cout << "\n-------------------- Fake cutflow (HNL μe VBF) --------------------\n";
    std::cout << std::left  << std::setw(40) << "Cut"
              << std::right << std::setw(18) << "Σ fake_weight_nominal" << "\n";
    std::cout << "-----------------------------------------------------------------\n";
    for (const auto& c : cuts) {
        std::cout << std::left  << std::setw(40) << c
                  << std::right << std::setw(18) << std::fixed << std::setprecision(2) << Y[c] << "\n";
    }
    std::cout << "-----------------------------------------------------------------\n";
}
