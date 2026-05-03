// =============================================================================
// FILE:    AMM/AMM_test_DAOD.C
// PURPOSE: Standalone, simple AMM weight calculator that reads a DAOD_PHYS
//          file directly and computes per-event fake weights — for comparison
//          with the FakeBkgCalculatorAlg output (weight_fake_electronAMM)
//          produced by the main framework.
//
// SELECTIONS APPLIED:
//   - exactly 1 tight muon (Tight quality), pT > 26 GeV, |η| < 2.5
//   - ≥ 1 electron passing LooseLH, pT > 15 GeV, |η| < 2.47
//   - ≥ 2 jets, pT > 20 GeV
//
// AMM FORMULA (per loose electron, summed):
//   w_T = ε_f (ε_r − 1) / (ε_r − ε_f)     if electron is tight
//   w_L = ε_f  ε_r      / (ε_r − ε_f)     if electron is loose-only
//
// USAGE (with AnalysisBase set up):
//   root -b -q -l 'AMM_test_DAOD.C("DAOD_PHYS.39629931._000161.pool.root.1")'
//
// OUTPUT:
//   AMM_test_output.root            (TH1D h_amm_weight)
//   stdout: per-event prints for first 20 passing events + summary
// =============================================================================
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODElectron/ElectronContainer.h"
#include "xAODMuon/MuonContainer.h"
#include "xAODJet/JetContainer.h"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <string>

// ---------------------------------------------------------------------------
// AMM weight for one electron given (tight/loose, ε_r, ε_f)
// ---------------------------------------------------------------------------
static double amm_weight(bool isTight, double er, double ef)
{
    // same clamps as AMM_ElectronOnly.C
    er = std::max(1e-5, std::min(1.0 - 1e-5, er));
    ef = std::max(1e-5, std::min(er - 1e-5, ef));
    const double denom = er - ef;
    if (std::abs(denom) < 1e-9) return 0.0;
    return isTight ? (ef / denom) * (er - 1.0)   // negative
                   : (ef / denom) *  er;          // positive
}

// ---------------------------------------------------------------------------
// 2D efficiency lookup with last-bin clamping (matches AMM_ElectronOnly.C)
// ---------------------------------------------------------------------------
static double get_eff(TH2D* h, double pt_gev, double eta)
{
    if (!h) return 0.0;
    const double abs_eta = std::abs(eta);
    int bx = std::max(1, std::min(h->GetXaxis()->FindBin(pt_gev), h->GetNbinsX()));
    int by = std::max(1, std::min(h->GetYaxis()->FindBin(abs_eta), h->GetNbinsY()));
    return h->GetBinContent(bx, by);
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void AMM_test_DAOD(const char* daodPath = "DAOD_PHYS.39629931._000161.pool.root.1",
                   const char* effPath  = "fake_eff_input.root")
{
    // ----- 1. Set up xAOD reading -----
    xAOD::Init().ignore();

    // ----- 2. Open efficiency file -----
    TFile* fe = TFile::Open(effPath, "READ");
    if (!fe || fe->IsZombie()) {
        std::cerr << "ERROR: cannot open " << effPath << "\n";
        return;
    }
    TH2D* hRealE = (TH2D*)fe->Get("RealEfficiency2D_el_pt_eta");
    TH2D* hFakeE = (TH2D*)fe->Get("FakeEfficiency2D_el_pt_eta");
    if (!hRealE || !hFakeE) {
        std::cerr << "ERROR: efficiency histograms not found in " << effPath << "\n";
        std::cerr << "       Expected: RealEfficiency2D_el_pt_eta, FakeEfficiency2D_el_pt_eta\n";
        return;
    }
    std::cout << "ε_r map: " << hRealE->GetNbinsX() << " × " << hRealE->GetNbinsY() << " bins\n";
    std::cout << "ε_f map: " << hFakeE->GetNbinsX() << " × " << hFakeE->GetNbinsY() << " bins\n";

    // ----- 3. Open DAOD -----
    TFile* fd = TFile::Open(daodPath, "READ");
    if (!fd || fd->IsZombie()) {
        std::cerr << "ERROR: cannot open " << daodPath << "\n";
        return;
    }
    xAOD::TEvent event(xAOD::TEvent::kClassAccess);
    if (event.readFrom(fd).isFailure()) {
        std::cerr << "ERROR: TEvent::readFrom failed\n";
        return;
    }

    const Long64_t n = event.getEntries();
    std::cout << "Total events in DAOD: " << n << "\n";

    // ----- 4. Event loop -----
    Long64_t nPass = 0, nPrinted = 0;
    double sumW = 0.0, sumWPos = 0.0, sumWNeg = 0.0;
    TH1D* hW = new TH1D("h_amm_weight", ";AMM event weight;events", 200, -0.5, 0.5);

    for (Long64_t i = 0; i < n; ++i) {
        event.getEntry(i);

        const xAOD::ElectronContainer* electrons = nullptr;
        const xAOD::MuonContainer*     muons     = nullptr;
        const xAOD::JetContainer*      jets      = nullptr;

        if (event.retrieve(electrons, "Electrons").isFailure())             continue;
        if (event.retrieve(muons,     "Muons").isFailure())                 continue;
        if (event.retrieve(jets,      "AntiKt4EMPFlowJets").isFailure())    continue;

        // -- Muons: exactly 1 Tight, pT > 26 GeV, |η| < 2.5
        int nGoodMu = 0;
        for (const xAOD::Muon* mu : *muons) {
            if (mu->pt() < 26000.) continue;
            if (std::abs(mu->eta()) > 2.5) continue;
            if (mu->quality() != xAOD::Muon::Tight) continue;
            ++nGoodMu;
        }
        if (nGoodMu != 1) continue;

        // -- Electrons: pT > 15, |η| < 2.47, must pass LooseLH;
        //    the AMM tight/loose flag = TightLH (which is what your ε_f
        //    map was measured at).
        std::vector<std::pair<bool, std::pair<double,double>>> looseE;
        for (const xAOD::Electron* el : *electrons) {
            if (el->pt() < 15000.) continue;
            if (std::abs(el->eta()) > 2.47) continue;

            char looseLH = 0, tightLH = 0;
            if (el->isAvailable<char>("DFCommonElectronsLHLoose"))
                looseLH = el->auxdataConst<char>("DFCommonElectronsLHLoose");
            if (el->isAvailable<char>("DFCommonElectronsLHTight"))
                tightLH = el->auxdataConst<char>("DFCommonElectronsLHTight");

            if (!looseLH) continue;
            looseE.push_back({(bool)tightLH, {el->pt(), el->eta()}});
        }
        if (looseE.empty()) continue;

        // -- Jets: ≥ 2 with pT > 20 GeV
        int nGoodJ = 0;
        for (const xAOD::Jet* j : *jets) {
            if (j->pt() > 20000.) ++nGoodJ;
        }
        if (nGoodJ < 2) continue;

        // -- Compute per-event AMM weight (sum over loose electrons)
        double w = 0.0;
        for (auto& el : looseE) {
            const bool   tight   = el.first;
            const double pt_gev  = el.second.first  * 1e-3;
            const double eta     = el.second.second;
            const double er      = get_eff(hRealE, pt_gev, eta);
            const double ef      = get_eff(hFakeE, pt_gev, eta);
            w += amm_weight(tight, er, ef);
        }

        ++nPass;
        sumW += w;
        if (w > 0) sumWPos += w; else sumWNeg += w;
        hW->Fill(w);

        if (nPrinted < 20) {
            std::cout << "Event " << i
                      << "  nLooseE=" << looseE.size()
                      << "  nTightE=" << std::count_if(looseE.begin(), looseE.end(),
                                                       [](auto& p){ return p.first; })
                      << "  weight=" << w << "\n";
            ++nPrinted;
        }
    }

    // ----- 5. Summary -----
    std::cout << "\n=== Summary ===\n";
    std::cout << "Events processed:           " << n        << "\n";
    std::cout << "Events passing selection:   " << nPass    << "\n";
    std::cout << "Sum of AMM weights:         " << sumW     << "  (predicted fake yield)\n";
    std::cout << "  positive contributions:   " << sumWPos  << "  (loose-only events)\n";
    std::cout << "  negative contributions:   " << sumWNeg  << "  (tight events)\n";
    std::cout << "Mean weight per event:      "
              << (nPass ? sumW/double(nPass) : 0.) << "\n";

    // ----- 6. Persist histogram -----
    TFile* fout = TFile::Open("AMM_test_output.root", "RECREATE");
    hW->Write();
    fout->Close();
    std::cout << "\nWrote AMM_test_output.root  (h_amm_weight)\n";
}
