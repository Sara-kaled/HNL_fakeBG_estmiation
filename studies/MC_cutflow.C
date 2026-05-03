// =============================================================================
// FILE:    studies/MC_cutflow.C
// PURPOSE: MC cutflow for the HNL μe VBF SR / VR. Uses utils/SelectionUtils.h
//          throughout — no duplicated branch handling. Truth-prompt cut is
//          applied at the end so the cutflow shows pure prompt MC. Final
//          split: SR (m_jj > 400 GeV) vs VR (m_jj < 400 GeV).
//
// SR/VR analysis triggers: mu OR e (passTrigger + passTriggerMatching).
// =============================================================================
#include <TChain.h>
#include <TH1D.h>
#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>

#include "../utils/SelectionUtils.h"

void MC_cutflow()
{
    TChain* tree = new TChain("tree");
    const std::string base = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    for (auto p : {"ttbar","wjets","ztautau","singletop","diboson","higgs"})
        tree->Add((base + p + "/*.root").c_str());

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
        "Truth-prompt mu-e",
        "VR  (mjj < 400 GeV)",
        "SR  (mjj > 400 GeV)"
    };
    std::map<std::string, double> Y;
    for (const auto& c : cuts) Y[c] = 0.;

    Long64_t n = tree->GetEntries();
    std::cout << "MC entries: " << n << "\n";
    for (Long64_t i = 0; i < n; ++i) {
        if (i % 1000000 == 0) std::cout << "  " << i << " / " << n << "\n";
        tree->GetEntry(i);
        ComputeDerived(ev);

        const double w = mcWeightTight(ev);
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

        if (!isMCPromptEvent(ev))                                         continue;  Y["Truth-prompt mu-e"] += w;

        const double mjj_gev = ev.mjj * 1e-3;
        if (mjj_gev < 400.) Y["VR  (mjj < 400 GeV)"] += w;
        else                Y["SR  (mjj > 400 GeV)"] += w;
    }

    std::cout << "\n--------------------- MC cutflow (HNL μe VBF) ---------------------\n";
    std::cout << std::left  << std::setw(40) << "Cut"
              << std::right << std::setw(18) << "Yield"
              << std::setw(12) << "Pass eff %\n";
    std::cout << "------------------------------------------------------------------\n";
    double prev = -1.;
    for (const auto& c : cuts) {
        double y = Y[c];
        double eff = (prev > 0.) ? 100. * y / prev : 100.;
        std::cout << std::left  << std::setw(40) << c
                  << std::right << std::setw(18) << std::fixed << std::setprecision(2) << y
                  << std::setw(12) << eff << "\n";
        if (c == "Initial" || c == "Truth-prompt mu-e" ||
            c == "VR  (mjj < 400 GeV)" || c == "SR  (mjj > 400 GeV)") prev = y;
        else prev = y;
    }
    std::cout << "------------------------------------------------------------------\n";
}
