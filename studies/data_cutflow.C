// =============================================================================
// FILE:    studies/data_cutflow.C
// PURPOSE: Data cutflow for the HNL μe VBF SR / VR. Uses
//          utils/SelectionUtils.h. SR/VR analysis triggers: mu OR e.
// =============================================================================
#include <TChain.h>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include "../utils/SelectionUtils.h"

void data_cutflow()
{
    TChain* tree = new TChain("tree");
    tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data15.root");
    tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data16.root");
    tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_d/data17.root");
    tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/data18.root");

    EventData ev;
    SetupBranches(tree, ev);

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
        "SR  (mjj > 400 GeV)  [BLINDED]"
    };
    std::map<std::string, double> Y;
    for (const auto& c : cuts) Y[c] = 0.;

    Long64_t n = tree->GetEntries();
    std::cout << "Data entries: " << n << "\n";
    for (Long64_t i = 0; i < n; ++i) {
        if (i % 1000000 == 0) std::cout << "  " << i << " / " << n << "\n";
        tree->GetEntry(i);
        ComputeDerived(ev);

        const double w = 1.0;
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
        else                Y["SR  (mjj > 400 GeV)  [BLINDED]"] += w;
    }

    std::cout << "\n-------------------- Data cutflow (HNL μe VBF) --------------------\n";
    std::cout << std::left  << std::setw(40) << "Cut"
              << std::right << std::setw(18) << "Events" << "\n";
    std::cout << "------------------------------------------------------------------\n";
    for (const auto& c : cuts) {
        std::cout << std::left  << std::setw(40) << c
                  << std::right << std::setw(18) << std::fixed << std::setprecision(0) << Y[c] << "\n";
    }
    std::cout << "------------------------------------------------------------------\n";
}
