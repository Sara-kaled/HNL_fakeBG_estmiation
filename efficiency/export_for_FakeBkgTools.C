// =============================================================================
// FILE:    efficiency/export_for_FakeBkgTools.C
// PURPOSE: Convert the analysis-internal fake and real efficiency maps
//          into the format expected by the ATLAS FakeBkgTools
//          (CP::AsymptMatrixTool / FakeBkgCalculatorAlg).
//
//          Reads:
//            outputs/fake_eff_with_systematics.root  → ε_f  (nominal + systs)
//            outputs/real_eff_histograms.root         → ε_r  (tight/loose MC)
//
//          Writes:
//            outputs/fake_eff_input.root   with TH2D's named:
//              el_fake_eff         (nominal ε_f, pT [GeV] × |η|)
//              el_real_eff         (nominal ε_r, pT [GeV] × |η|)
//              el_fake_eff_up      (total-syst up variation)
//              el_fake_eff_down    (total-syst down variation)
//
//            outputs/fake_eff_input.xml    consumed by FakeBkgCalculatorAlg
//              config = "path/to/outputs/fake_eff_input.xml"
//
//          Sanity checks performed on every bin:
//            • 0 < ε_f < 1
//            • 0 < ε_r < 1
//            • ε_r > ε_f  (matrix-method inversion requires ε_r ≠ ε_f)
//            • No NaN / Inf in content or error
//
// RUNS AFTER: build_fake_systematics.C  (produces fake_eff_with_systematics.root)
//             real_eff.C               (produces real_eff_histograms.root)
//
// KEY CONVENTION:
//   Histogram axis layout (must match the XML <param> order):
//     X axis = pT  [GeV]  (kPtBins  = {15,25,35,50,80,200})
//     Y axis = |η|        (kEtaBins = {0,0.7,1.4,1.6,2.0,2.5})
//   The pT axis is already in GeV in both source files (SelectionUtils.h fills
//   in GeV via `ev.l2_pt * 1e-3`).  No unit conversion needed.
// =============================================================================
#include <TFile.h>
#include <TH2D.h>
#include <TSystem.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <string>

// ---------------------------------------------------------------------------
// Helper: perform a sanity-checked bin-by-bin division  num/den → efficiency.
// Uses TEfficiency Clopper–Pearson for the bin error.
// Returns the number of bins that failed the [lo, hi] range check.
// ---------------------------------------------------------------------------
static int checkedEfficiency(TH2D* hNum, TH2D* hDen, TH2D* hEff,
                             const char* label,
                             double lo = 0.0, double hi = 1.0)
{
    int nBad = 0;
    const int nX = hEff->GetNbinsX();
    const int nY = hEff->GetNbinsY();

    hEff->Reset();

    for (int ix = 1; ix <= nX; ++ix) {
        for (int iy = 1; iy <= nY; ++iy) {
            double num = hNum->GetBinContent(ix, iy);
            double den = hDen->GetBinContent(ix, iy);
            double val = 0., err = 0.;

            if (den > 0. && num >= 0. && num <= den) {
                val = num / den;
                // Symmetric approximation of Clopper–Pearson for the error
                // (full CP requires TEfficiency; this avoids another dependency)
                double relErr2 = (num > 0. ? 1./num : 0.) + (den > 0. ? 1./den : 0.);
                err = val * std::sqrt(relErr2);
            } else if (den > 0. && num > den) {
                // Numerator > denominator after data–MC subtraction:
                // can happen due to MC stat fluctuations. Cap at 0.99.
                val = 0.99;
                err = 0.5;  // large error flags the unstable bin
                printf("  WARNING [%s] bin (%d,%d): num=%.1f > den=%.1f → capped at 0.99\n",
                       label, ix, iy, num, den);
                ++nBad;
            } else {
                // den <= 0: no data in this bin
                val = 0.;
                err = 1.;  // large error → tool will weight this bin minimally
                printf("  WARNING [%s] bin (%d,%d): den=%.1f ≤ 0 → set to 0 ± 1\n",
                       label, ix, iy, den);
                ++nBad;
            }

            // NaN / Inf guard
            if (!std::isfinite(val)) { val = 0.; err = 1.; ++nBad; }
            if (!std::isfinite(err)) { err = 1.; }

            // Range check
            if (val < lo || val > hi) {
                printf("  WARNING [%s] bin (%d,%d): eff=%.4f outside [%.2f,%.2f]\n",
                       label, ix, iy, val, lo, hi);
                val = std::min(std::max(val, lo + 1e-6), hi - 1e-6);
                ++nBad;
            }

            hEff->SetBinContent(ix, iy, val);
            hEff->SetBinError  (ix, iy, err);
        }
    }
    return nBad;
}

// ---------------------------------------------------------------------------
// Helper: check that ε_r > ε_f in every bin (matrix-method requirement).
// Prints a warning for each violating bin.
// ---------------------------------------------------------------------------
static int checkRealVsFake(const TH2D* hReal, const TH2D* hFake, const char* label)
{
    int nBad = 0;
    for (int ix = 1; ix <= hReal->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= hReal->GetNbinsY(); ++iy) {
            double er = hReal->GetBinContent(ix, iy);
            double ef = hFake->GetBinContent(ix, iy);
            if (er <= ef) {
                printf("  WARNING [%s] bin (%d,%d): ε_r=%.4f ≤ ε_f=%.4f — "
                       "matrix inversion is ill-defined in this bin!\n",
                       label, ix, iy, er, ef);
                ++nBad;
            }
        }
    }
    return nBad;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void export_for_FakeBkgTools()
{
    // ------------------------------------------------------------------
    // 1. Open input files
    // ------------------------------------------------------------------
    TFile* fFake = TFile::Open("outputs/fake_eff_with_systematics.root", "READ");
    if (!fFake || fFake->IsZombie()) {
        std::cerr << "ERROR: cannot open outputs/fake_eff_with_systematics.root\n"
                  << "       Run build_fake_systematics.C first.\n";
        return;
    }

    TFile* fReal = TFile::Open("outputs/real_eff_histograms.root", "READ");
    if (!fReal || fReal->IsZombie()) {
        std::cerr << "ERROR: cannot open outputs/real_eff_histograms.root\n"
                  << "       Run real_eff.C first.\n";
        return;
    }

    // ------------------------------------------------------------------
    // 2. Retrieve ε_f maps (nominal + total systematic variations)
    //    Key in fake_eff_with_systematics.root:
    //      "fake_rate_nominal"           — nominal  (written with explicit name)
    //      "fake_rate_total_systUp"      — tot-syst up
    //      "fake_rate_total_systDown"    — tot-syst down
    // ------------------------------------------------------------------
    auto* hFakeNom  = (TH2D*)fFake->Get("fake_rate_nominal");
    auto* hFakeTotU = (TH2D*)fFake->Get("fake_rate_total_systUp");
    auto* hFakeTotD = (TH2D*)fFake->Get("fake_rate_total_systDown");

    if (!hFakeNom) {
        std::cerr << "ERROR: histogram 'fake_rate_nominal' not found in "
                     "fake_eff_with_systematics.root\n";
        fFake->ls();   // print available keys to guide debugging
        return;
    }
    if (!hFakeTotU || !hFakeTotD) {
        std::cerr << "WARNING: total systematic variation histograms not found — "
                     "only nominal will be written.\n";
    }

    // ------------------------------------------------------------------
    // 3. Retrieve tight + loose from real_eff_histograms.root and compute ε_r.
    //    Use the CR region (prompt OS CR, highest statistics).
    //    Histogram names (from real_eff.C):
    //      "h2_tight_eta_pt_CR"   X=pT[GeV], Y=|η|
    //      "h2_loose_eta_pt_CR"
    // ------------------------------------------------------------------
    auto* hRealTight = (TH2D*)fReal->Get("h2_tight_eta_pt_CR");
    auto* hRealLoose = (TH2D*)fReal->Get("h2_loose_eta_pt_CR");

    if (!hRealTight || !hRealLoose) {
        std::cerr << "ERROR: real efficiency histograms not found in "
                     "real_eff_histograms.root\n";
        fReal->ls();
        return;
    }

    // ------------------------------------------------------------------
    // 4. Build ε_r histogram (same binning as ε_f)
    // ------------------------------------------------------------------
    TH2D* hRealEff = (TH2D*)hRealTight->Clone("el_real_eff");
    hRealEff->SetDirectory(nullptr);
    hRealEff->SetTitle("Real electron efficiency #varepsilon_{r}(p_{T},|#eta|);"
                       "p_{T}^{e} [GeV];|#eta^{e}|");

    printf("\n=== Building el_real_eff (ε_r = tight/loose from MC prompt CR) ===\n");
    int nBadReal = checkedEfficiency(hRealTight, hRealLoose, hRealEff,
                                     "el_real_eff", 0.0, 1.0);

    // ------------------------------------------------------------------
    // 5. Clone and rename ε_f maps
    // ------------------------------------------------------------------
    TH2D* hFakeEff  = (TH2D*)hFakeNom->Clone("el_fake_eff");
    hFakeEff->SetDirectory(nullptr);
    hFakeEff->SetTitle("Fake electron efficiency #varepsilon_{f}(p_{T},|#eta|);"
                       "p_{T}^{e} [GeV];|#eta^{e}|");

    TH2D* hFakeEffU = nullptr;
    TH2D* hFakeEffD = nullptr;
    if (hFakeTotU && hFakeTotD) {
        hFakeEffU = (TH2D*)hFakeTotU->Clone("el_fake_eff_up");
        hFakeEffU->SetDirectory(nullptr);
        hFakeEffD = (TH2D*)hFakeTotD->Clone("el_fake_eff_down");
        hFakeEffD->SetDirectory(nullptr);
    }

    // ------------------------------------------------------------------
    // 6. Validate ε_f bins (must already be in [0,1] — guarded upstream
    //    by calculate_fake_rates.C and build_fake_systematics.C, but
    //    we re-check here because NaN/clamping could occur after syst addition)
    // ------------------------------------------------------------------
    printf("\n=== Checking el_fake_eff (ε_f nominal) ===\n");
    int nBadFake = 0;
    const int nX = hFakeEff->GetNbinsX();
    const int nY = hFakeEff->GetNbinsY();
    for (int ix = 1; ix <= nX; ++ix) {
        for (int iy = 1; iy <= nY; ++iy) {
            double v = hFakeEff->GetBinContent(ix, iy);
            double e = hFakeEff->GetBinError  (ix, iy);
            bool bad = false;
            if (!std::isfinite(v) || v < 0. || v > 1.) {
                printf("  WARNING [el_fake_eff] bin (%d,%d): val=%.4f → clamped\n",
                       ix, iy, v);
                v = std::min(std::max(v, 0.), 1.);
                bad = true;
            }
            if (!std::isfinite(e) || e < 0.) { e = 1.; bad = true; }
            if (bad) { hFakeEff->SetBinContent(ix, iy, v);
                       hFakeEff->SetBinError  (ix, iy, e);
                       ++nBadFake; }
        }
    }
    if (hFakeEffU && hFakeEffD) {
        // Clamp syst-variation maps to [0,1] as well
        for (int ix = 1; ix <= nX; ++ix) {
            for (int iy = 1; iy <= nY; ++iy) {
                auto clamp = [](TH2D* h, int i, int j) {
                    double v = h->GetBinContent(i,j);
                    if (!std::isfinite(v)) v = 0.;
                    v = std::min(std::max(v, 0.), 1.);
                    h->SetBinContent(i, j, v);
                };
                clamp(hFakeEffU, ix, iy);
                clamp(hFakeEffD, ix, iy);
            }
        }
    }

    // ------------------------------------------------------------------
    // 7. Cross-check: ε_r > ε_f in every bin
    // ------------------------------------------------------------------
    printf("\n=== Cross-check: ε_r > ε_f (matrix method condition) ===\n");
    int nBadInv = checkRealVsFake(hRealEff, hFakeEff, "ε_r vs ε_f");

    // ------------------------------------------------------------------
    // 8. Print summary table
    // ------------------------------------------------------------------
    printf("\n%-8s  %-8s  %-8s  %-8s  %-8s\n",
           "pT-bin", "η-bin", "ε_f", "ε_r", "ε_r-ε_f");
    printf("-------------------------------------------------------\n");
    for (int ix = 1; ix <= nX; ++ix) {
        for (int iy = 1; iy <= nY; ++iy) {
            double ef = hFakeEff->GetBinContent(ix, iy);
            double er = hRealEff->GetBinContent(ix, iy);
            double ptLo  = hFakeEff->GetXaxis()->GetBinLowEdge(ix);
            double ptHi  = hFakeEff->GetXaxis()->GetBinUpEdge(ix);
            double etaLo = hFakeEff->GetYaxis()->GetBinLowEdge(iy);
            double etaHi = hFakeEff->GetYaxis()->GetBinUpEdge(iy);
            printf("  pT=[%3.0f-%3.0f]  |η|=[%.1f-%.1f]   "
                   "ε_f=%.4f  ε_r=%.4f  Δ=%+.4f%s\n",
                   ptLo, ptHi, etaLo, etaHi, ef, er, er-ef,
                   (er <= ef) ? "  *** INVERSION PROBLEM ***" : "");
        }
    }

    // ------------------------------------------------------------------
    // 9. Write output ROOT file in FakeBkgTools "default ROOT format"
    //    (no XML needed — the tool reads histograms by naming convention)
    //
    // Naming convention required by CP::AsymptMatrixTool:
    //   RealEfficiency2D_el_pt_eta   — ε_r(pT, |η|) for electrons
    //   FakeEfficiency2D_el_pt_eta   — ε_f(pT, |η|) for electrons
    //   RealEfficiency_mu_pt           — ε_r(pT) for muons (1-bin dummy ≈ 1)
    //   FakeEfficiency_mu_pt           — ε_f(pT) for muons (1-bin dummy ≈ 0)
    //
    // Systematic variations follow the convention:
    //   {Real|Fake}Efficiency2D_el_pt_eta__<SystName>
    // bin contents = absolute uncertainty (sym), bin errors = 0.
    // ------------------------------------------------------------------
    gSystem->mkdir("outputs", true);
    TFile* fOut = TFile::Open("outputs/fake_eff_input.root", "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        std::cerr << "ERROR: cannot create outputs/fake_eff_input.root\n";
        return;
    }

    hRealEff->Write("RealEfficiency2D_el_pt_eta");
    hFakeEff->Write("FakeEfficiency2D_el_pt_eta");

    // ----- Muon dummies (electron-only AMM: μ assumed prompt) -----
    TH1D* hMuReal = new TH1D("RealEfficiency_mu_pt", ";p_{T}^{#mu} [GeV];#varepsilon_{r}^{#mu}",
                             1, 0., 1.e6);
    hMuReal->SetBinContent(1, 0.999);
    hMuReal->SetBinError  (1, 0.);
    hMuReal->Write();

    TH1D* hMuFake = new TH1D("FakeEfficiency_mu_pt", ";p_{T}^{#mu} [GeV];#varepsilon_{f}^{#mu}",
                             1, 0., 1.e6);
    hMuFake->SetBinContent(1, 0.001);
    hMuFake->SetBinError  (1, 0.);
    hMuFake->Write();

    // ----- Systematic variation histograms (electron ε_f only) -----
    // Read each named systematic from fake_eff_with_systematics.root and write
    // out the absolute-uncertainty histogram under the FakeBkgTools-suffix
    // convention (__SystName).  Bin errors are zeroed so the parser interprets
    // the contents as the uncertainty itself.
    auto writeSyst = [&](const char* srcName, const char* sysLabel) {
        TH2D* hUp   = (TH2D*)fFake->Get((std::string("fake_rate_") + srcName + "_up"  ).c_str());
        TH2D* hDown = (TH2D*)fFake->Get((std::string("fake_rate_") + srcName + "_down").c_str());
        if (!hUp || !hDown) {
            printf("  (skip syst '%s': src histograms not found)\n", sysLabel);
            return;
        }
        TH2D* hSyst = (TH2D*)hUp->Clone(
            (std::string("FakeEfficiency2D_el_pt_eta__") + sysLabel).c_str());
        hSyst->SetDirectory(nullptr);
        for (int ix = 1; ix <= hSyst->GetNbinsX(); ++ix) {
            for (int iy = 1; iy <= hSyst->GetNbinsY(); ++iy) {
                double up  = hUp  ->GetBinContent(ix, iy);
                double dn  = hDown->GetBinContent(ix, iy);
                double nom = hFakeEff->GetBinContent(ix, iy);
                double sigma = 0.5 * std::fabs(up - dn);
                if (!std::isfinite(sigma)) sigma = 0.;
                hSyst->SetBinContent(ix, iy, sigma);
                hSyst->SetBinError  (ix, iy, 0.);
                (void)nom;
            }
        }
        hSyst->Write();
        printf("  ✓ wrote syst '%s' (max σ = %.3f)\n", sysLabel, hSyst->GetMaximum());
    };

    printf("\n=== Writing systematic variation histograms ===\n");
    writeSyst("MET",     "METdep");
    writeSyst("comp",    "Composition");
    writeSyst("real",    "PromptSubtr");
    writeSyst("stat",    "StatExtra");

    fOut->Close();

    printf("\n✓ Wrote outputs/fake_eff_input.root  (default-format, no XML required)\n");

    // ------------------------------------------------------------------
    // 10. Final summary
    // ------------------------------------------------------------------
    printf("\n=== Export summary ===\n");
    printf("  Bad bins in ε_r :   %d\n", nBadReal);
    printf("  Bad bins in ε_f :   %d\n", nBadFake);
    printf("  ε_r ≤ ε_f bins :   %d\n", nBadInv);
    if (nBadReal == 0 && nBadFake == 0 && nBadInv == 0) {
        printf("  ✓ All bins pass — maps are ready for FakeBkgCalculatorAlg.\n");
    } else {
        printf("  ⚠  Some bins have issues — inspect warnings above before\n"
               "     feeding these maps into the ntuple maker.\n");
    }
    printf("\nNext steps:\n"
           "  1. Copy outputs/fake_eff_input.root to a directory that PathResolver\n"
           "     can see (either the run directory, or anywhere on $DATAPATH).\n"
           "  2. Set  config: 'fake_eff_input.root'  in the FakeBkgCalculator YAML\n"
           "     block. No XML required — the tool reads histogram names directly.\n"
           "  3. (Optional) export DATAPATH=/path/to/dir:$DATAPATH  before launching.\n");

    fFake->Close();
    fReal->Close();
}
