// =============================================================================
// FILE:       yields/fake_yields_correct.C
// PURPOSE:    Final yield table: prompt MC backgrounds + data-driven fake prediction
//             vs observed data in Preselection, CR (SS), VR, and SR.
//             MC loop:   sums prompt-only events (isMCPromptEvent) with tight SF.
//             Data loop: sums fake_weight_nominal (and systematics) from AMM-weighted
//             data ntuple. Region definitions from SelectionUtils.h match the
//             stat-framework config exactly.
// INPUTS:     MC ntuples, AMM/data_with_electron_fake_weights.root
// OUTPUTS:    Yield table printed to stdout
// RUNS AFTER: AMM/AMM_ElectronOnly.C
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

using std::cout;
using std::endl;

#include "../utils/SelectionUtils.h"

// ---------------------------------------------------------------------------
struct Yield { double sumw=0., sumw2=0.; };

struct FakeInfo {
    double nom=0., sumw2_nom=0.;
    double met_up=0., met_dn=0.;
    double comp_up=0., comp_dn=0.;
    double stat_up=0., stat_dn=0.;
};

enum Region { kPre=0, kCR=1, kVR=2, kSR=3, kNRegion=4 };
static const char* regionNames[kNRegion] = { "Preselection","CR","VR","SR" };

// ---------------------------------------------------------------------------
inline void addFake(FakeInfo& F, const EventData& ev) {
    F.nom      += ev.fw_nom;
    F.sumw2_nom+= ev.fw_nom * ev.fw_nom;
    F.met_up   += ev.fw_met_up;
    F.met_dn   += ev.fw_met_dn;
    F.comp_up  += ev.fw_comp_up;
    F.comp_dn  += ev.fw_comp_dn;
    F.stat_up  += ev.fw_stat_up;
    F.stat_dn  += ev.fw_stat_dn;
}

// ---------------------------------------------------------------------------
void fake_yields_correct()
{
    const std::string base = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";

    // -----------------------------------------------------------------------
    // 1. Define MC processes
    // -----------------------------------------------------------------------
    std::map<std::string, std::string> mcPaths = {
        {"ttbar",     base+"ttbar/*.root"},
        {"ztautau",   base+"ztautau/*.root"},
        {"wjets",     base+"wjets/*.root"},
        {"singletop", base+"singletop/*.root"},
        {"diboson",   base+"diboson/*.root"},
        {"higgs",     base+"higgs/*.root"},
    };
    std::vector<std::string> mcNames = {
        "ttbar","ztautau","wjets","singletop","diboson","higgs"
    };

    std::map<std::string, Yield[kNRegion]> mcY;
    FakeInfo fakeReg[kNRegion];
    double   dataObs[kNRegion] = {};

    // -----------------------------------------------------------------------
    // 2. MC loop (prompt-only)
    // -----------------------------------------------------------------------
    cout << "\n=== MC (prompt only) ===\n";

    for (const auto& pname : mcNames) {
        cout << "  " << pname << "\n";
        TChain ch("tree");
        ch.Add(mcPaths[pname].c_str());

        EventData ev;
        SetupBranches(&ch, ev);
        Long64_t n = ch.GetEntries();

        for (Long64_t i = 0; i < n; ++i) {
            ch.GetEntry(i);
            ComputeDerived(ev);

            if (!passBasePresel(ev))     continue;
            if (!isMCPromptEvent(ev))    continue;   // prompt-only (no double-counting with fake)
            if (!ev.mu_is_tight)         continue;
            if (!ev.el_is_tight)         continue;

            double wt = mcWeightTight(ev);

            // kPre
            mcY[pname][kPre].sumw  += wt;
            mcY[pname][kPre].sumw2 += wt*wt;

            if (!passVBFJets(ev)) continue;

            // kCR — SS, use loose SF (AMM application CR)
            if (isFakeCR(ev) && ev.no_bjets) {
                double wl = mcWeightLoose(ev);
                mcY[pname][kCR].sumw  += wl;
                mcY[pname][kCR].sumw2 += wl*wl;
            }

            // kSR / kVR — signal kinematics (tight)
            if (isSR(ev)) {
                mcY[pname][kSR].sumw  += wt;
                mcY[pname][kSR].sumw2 += wt*wt;
            }
            if (isVR(ev)) {
                mcY[pname][kVR].sumw  += wt;
                mcY[pname][kVR].sumw2 += wt*wt;
            }
        }
    }

    // -----------------------------------------------------------------------
    // 3. Data loop (fake prediction + raw observed)
    // -----------------------------------------------------------------------
    cout << "\n=== Data (fake weights) ===\n";

    TChain dch("tree");
    dch.Add("../AMM/data_with_electron_fake_weights.root");

    EventData ev;
    SetupBranches(&dch, ev);
    Long64_t nData = dch.GetEntries();

    for (Long64_t i = 0; i < nData; ++i) {
        dch.GetEntry(i);
        ComputeDerived(ev);

        if (!passBasePresel(ev)) continue;

        // kPre — before any VBF / b-jet cut
        dataObs[kPre] += 1.;
        addFake(fakeReg[kPre], ev);

        if (!ev.no_bjets) continue;   // b-jet veto for CR/VR/SR

        // kCR — SS + VBF + no b-jets (AMM sums over loose, giving fake prediction)
        if (isFakeCR(ev) && passVBFJets(ev)) {
            dataObs[kCR] += 1.;
            addFake(fakeReg[kCR], ev);
        }

        // kSR / kVR — OS + signal kinematics (AMM: sum over loose → fake prediction)
        if (ev.el_is_loose) {
            if (isSR(ev)) {
                dataObs[kSR] += 1.;
                addFake(fakeReg[kSR], ev);
            }
            if (isVR(ev)) {
                dataObs[kVR] += 1.;
                addFake(fakeReg[kVR], ev);
            }
        }
    }

    // -----------------------------------------------------------------------
    // 4. Print yield table
    // -----------------------------------------------------------------------
    cout << "\n=== Yields per region ===\n";

    for (int r = 0; r < kNRegion; ++r) {
        cout << "\nRegion: " << regionNames[r] << "\n";
        printf("+----------------+------------+\n");
        printf("| %-14s | %10s |\n","Process","Yield");
        printf("+----------------+------------+\n");

        double Nmc=0., NmcErr2=0.;
        for (const auto& pname : mcNames) {
            printf("| %-14s | %10.2f |\n", pname.c_str(), mcY[pname][r].sumw);
            Nmc     += mcY[pname][r].sumw;
            NmcErr2 += mcY[pname][r].sumw2;
        }

        const FakeInfo& F = fakeReg[r];
        printf("| %-14s | %10.2f |\n","Fake", F.nom);
        printf("+----------------+------------+\n");

        double errMC   = NmcErr2  >0. ? std::sqrt(NmcErr2)  : 0.;
        double errFStat= F.sumw2_nom>0.? std::sqrt(F.sumw2_nom) : 0.;

        double dMET  = std::max(std::fabs(F.met_up -F.nom), std::fabs(F.met_dn -F.nom));
        double dComp = std::max(std::fabs(F.comp_up-F.nom), std::fabs(F.comp_dn-F.nom));
        double dEStat= std::max(std::fabs(F.stat_up-F.nom), std::fabs(F.stat_dn-F.nom));
        double errFSyst= std::sqrt(dMET*dMET + dComp*dComp + dEStat*dEStat);

        double Nbkg  = Nmc + F.nom;
        double errTot= std::sqrt(errMC*errMC + errFStat*errFStat + errFSyst*errFSyst);

        printf("| %-14s | %10.2f |\n","TotalMC",  Nmc);
        printf("| %-14s | %10.2f |\n","TotalFake", F.nom);
        printf("| %-14s | %10.2f |\n","TotalBkg",  Nbkg);
        printf("| %-14s | %10.0f |\n","DataObs",   dataObs[r]);
        printf("+----------------+------------+\n");

        double fn = std::fabs(F.nom);
        printf("Total uncertainty: ±%.2f (%.1f%%)\n", errTot, Nbkg>0.?100.*errTot/Nbkg:0.);
        printf("  MC stat   : %.2f (%.1f%% of MC)\n",    errMC,   Nmc>0.?100.*errMC/Nmc:0.);
        printf("  Fake stat : %.2f (%.1f%% of |Fake|)\n",errFStat, fn>0.?100.*errFStat/fn:0.);
        printf("  Fake syst : %.2f (%.1f%% of |Fake|)\n",errFSyst, fn>0.?100.*errFSyst/fn:0.);
    }

    cout << "\nDone.\n";
}
