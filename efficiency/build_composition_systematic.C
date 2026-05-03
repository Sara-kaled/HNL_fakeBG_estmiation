// =============================================================================
// FILE:    efficiency/build_composition_systematic.C
// PURPOSE: Build the COMPOSITION systematic on the fake electron efficiency
//          using MC truth-based fake-source reweighting.
//
// Categories are now detailed:
//   B       = b-hadron decays              IFFClass 8
//   C       = c-hadron decays              IFFClass 9
//   Conv    = photon conversions           IFFClass 5
//   LF      = light-flavour mis-ID         IFFClass 10
//   Tau     = tau decays                   IFFClass 7
//   Unknown = unknown                      IFFClass 0
//   Other   = known-unknown / unexpected   IFFClass 1, etc.
//
// Only prompt isolated electrons are skipped:
//   IFFClass 2 = IsoElectron
// =============================================================================

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TSystem.h>
#include <TString.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

#include "ROOT/RVec.hxx"
#include "../utils/SelectionUtils.h"

using namespace ROOT;
using namespace ROOT::VecOps;

// ---------------------------------------------------------------------------
// IFFClass -> detailed fake-source category.
// NOTE: prefix every enumerator with `kFS_` to avoid clashes with ROOT's own
// global enums.  In particular, ROOT pulls `kUnknown=1` from
// TStructNode.h into the global scope when CINT/Cling parses macros, so any
// bare `kUnknown` here would redefine that ROOT enumerator and fail to
// compile under ACLiC.
// ---------------------------------------------------------------------------
enum FakeSrc {
    kFS_B       = 0,
    kFS_C       = 1,
    kFS_Conv    = 2,
    kFS_LF      = 3,
    kFS_Tau     = 4,
    kFS_Unknown = 5,
    kFS_Other   = 6,
    kNSrc       = 7
};

static const char* kSrcName[kNSrc] = {
    "B",
    "C",
    "Conv",
    "LF",
    "Tau",
    "Unknown",
    "Other"
};

static int classifyFake(int iff)
{
    switch (iff) {
        case 8:  return kFS_B;        // BHadronDecay
        case 9:  return kFS_C;        // CHadronDecay
        case 5:  return kFS_Conv;     // PromptPhotonConversion
        case 10: return kFS_LF;       // LightFlavorDecay
        case 7:  return kFS_Tau;      // TauDecay
        case 0:  return kFS_Unknown;  // Unknown
        default: return kFS_Other;    // 1 = KnownUnknown, plus unexpected codes
    }
}

void build_composition_systematic()
{
    gSystem->mkdir("outputs", true);

    // ------------------------------------------------------------------
    // 1. MC TChain
    // ------------------------------------------------------------------
    TChain* tree = new TChain("tree");

    const std::string base =
        "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";

    for (auto p : {"ttbar", "wjets", "ztautau", "singletop", "diboson", "higgs"}) {
        tree->Add((base + p + "/*.root").c_str());
    }

    EventData ev;
    SetupBranches(tree, ev);

    RVecF* el_pt_NOSYS = nullptr;
    RVecI* el_IFFClass = nullptr;

    tree->SetBranchAddress("el_pt_NOSYS", &el_pt_NOSYS);
    tree->SetBranchAddress("el_IFFClass", &el_IFFClass);

    // ------------------------------------------------------------------
    // 2. Per-source 2D histograms in CR and SR
    // ------------------------------------------------------------------
    TH2D* h_tight_CR[kNSrc] = {nullptr};
    TH2D* h_loose_CR[kNSrc] = {nullptr};
    TH2D* h_tight_SR[kNSrc] = {nullptr};
    TH2D* h_loose_SR[kNSrc] = {nullptr};

    for (int s = 0; s < kNSrc; ++s) {
        auto mk = [&](const char* tag) {
            TString nm = TString::Format("h2_%s_%s", tag, kSrcName[s]);

            auto* h = new TH2D(
                nm,
                ";p_{T}^{e} [GeV];|#eta^{e}|",
                kNPtBins,
                kPtBins,
                kNEtaBins,
                kEtaBins
            );

            h->Sumw2();
            h->SetDirectory(nullptr);
            return h;
        };

        h_tight_CR[s] = mk("tight_CR");
        h_loose_CR[s] = mk("loose_CR");
        h_tight_SR[s] = mk("tight_SR");
        h_loose_SR[s] = mk("loose_SR");
    }

    // ------------------------------------------------------------------
    // 3. Event loop
    // ------------------------------------------------------------------
    Long64_t n = tree->GetEntries();
    std::cout << "Total MC entries: " << n << std::endl;

    Long64_t n_fakes_CR = 0;
    Long64_t n_fakes_SR = 0;

    Long64_t n_src[kNSrc] = {0};

    for (Long64_t i = 0; i < n; ++i) {
        if (i % 200000 == 0) {
            std::cout << "  entry " << i << " / " << n << "\n";
        }

        tree->GetEntry(i);
        ComputeDerived(ev);

        if (!passBasePresel(ev)) continue;
        if (!ev.el_is_loose) continue;
        if (!ev.no_bjets) continue;

        if (!el_IFFClass || el_IFFClass->empty()) continue;

        const int iff = (*el_IFFClass)[0];

        // Skip only prompt isolated electrons.
        // Photon conversions, IFFClass 5, are now kept as Conv.
        if (iff == 2) continue;

        const int src = classifyFake(iff);
        ++n_src[src];

        const bool pass_CR =
            isFakeCR(ev) &&
            (ev.jet1_pt > 0.) &&
            (ev.jet2_pt > 0.);

        const bool pass_SR =
            ev.is_os &&
            passVBFJets(ev);

        if (!pass_CR && !pass_SR) continue;

        const double pt_gev = ev.l2_pt * 1e-3;
        const double aeta   = ev.el_aeta;
        const double w      = mcWeightLoose(ev);

        if (!std::isfinite(w) || w == 0.0) continue;

        const bool tight = ev.el_is_tight;

        if (pass_CR) {
            ++n_fakes_CR;
            h_loose_CR[src]->Fill(pt_gev, aeta, w);
            if (tight) {
                h_tight_CR[src]->Fill(pt_gev, aeta, w);
            }
        }

        if (pass_SR) {
            ++n_fakes_SR;
            h_loose_SR[src]->Fill(pt_gev, aeta, w);
            if (tight) {
                h_tight_SR[src]->Fill(pt_gev, aeta, w);
            }
        }
    }

    std::cout << "MC fake-electron events: CR=" << n_fakes_CR
              << "   SR=" << n_fakes_SR << "\n";

    std::cout << "\nDetailed fake-source counts before region split:\n";
    for (int s = 0; s < kNSrc; ++s) {
        std::cout << "  " << kSrcName[s] << " = " << n_src[s] << "\n";
    }

    // ------------------------------------------------------------------
    // 4. Build per-source efficiencies and source fractions
    // ------------------------------------------------------------------
    auto cloneEmpty = [&](const TH2D* tmpl, const char* nm) {
        auto* h = (TH2D*)tmpl->Clone(nm);
        h->Reset();
        h->SetDirectory(nullptr);
        return h;
    };

    TH2D* eps_s[kNSrc] = {nullptr};
    TH2D* frac_s_CR[kNSrc] = {nullptr};
    TH2D* frac_s_SR[kNSrc] = {nullptr};

    for (int s = 0; s < kNSrc; ++s) {
        eps_s[s] = cloneEmpty(
            h_loose_CR[s],
            (std::string("eps_") + kSrcName[s]).c_str()
        );

        frac_s_CR[s] = cloneEmpty(
            h_loose_CR[s],
            (std::string("frac_") + kSrcName[s] + "_CR").c_str()
        );

        frac_s_SR[s] = cloneEmpty(
            h_loose_SR[s],
            (std::string("frac_") + kSrcName[s] + "_SR").c_str()
        );
    }

    const int nX = h_loose_CR[0]->GetNbinsX();
    const int nY = h_loose_CR[0]->GetNbinsY();

    const double kMinStat = 5.0;

    for (int ix = 1; ix <= nX; ++ix) {
        for (int iy = 1; iy <= nY; ++iy) {

            for (int s = 0; s < kNSrc; ++s) {
                const double t = h_tight_CR[s]->GetBinContent(ix, iy);
                const double l = h_loose_CR[s]->GetBinContent(ix, iy);

                const double e =
                    (l > 0.)
                    ? std::min(std::max(t / l, 0.), 1.)
                    : 0.;

                eps_s[s]->SetBinContent(ix, iy, e);

                eps_s[s]->SetBinError(
                    ix,
                    iy,
                    (l > 0.) ? std::sqrt(std::abs(t) * (1.0 - e) / l) : 0.
                );
            }

            double tot_loose_CR = 0.;
            double tot_loose_SR = 0.;

            for (int s = 0; s < kNSrc; ++s) {
                tot_loose_CR += h_loose_CR[s]->GetBinContent(ix, iy);
                tot_loose_SR += h_loose_SR[s]->GetBinContent(ix, iy);
            }

            for (int s = 0; s < kNSrc; ++s) {
                const double fCR =
                    (tot_loose_CR > 0.)
                    ? h_loose_CR[s]->GetBinContent(ix, iy) / tot_loose_CR
                    : 0.;

                const double fSR =
                    (tot_loose_SR > 0.)
                    ? h_loose_SR[s]->GetBinContent(ix, iy) / tot_loose_SR
                    : 0.;

                frac_s_CR[s]->SetBinContent(ix, iy, fCR);
                frac_s_SR[s]->SetBinContent(ix, iy, fSR);
            }
        }
    }

    // ------------------------------------------------------------------
    // 5. Reweighted predictions and fractional composition systematic
    // ------------------------------------------------------------------
    TH2D* eps_pred_CR = cloneEmpty(h_loose_CR[0], "eps_pred_CR");
    TH2D* eps_pred_SR = cloneEmpty(h_loose_CR[0], "eps_pred_SR");
    TH2D* frac_comp   = cloneEmpty(h_loose_CR[0], "frac_comp");

    int nLowStat = 0;
    int nValid = 0;

    for (int ix = 1; ix <= nX; ++ix) {
        for (int iy = 1; iy <= nY; ++iy) {

            double tot_CR = 0.;
            double tot_SR = 0.;

            for (int s = 0; s < kNSrc; ++s) {
                tot_CR += h_loose_CR[s]->GetBinContent(ix, iy);
                tot_SR += h_loose_SR[s]->GetBinContent(ix, iy);
            }

            double pCR = 0.;
            double pSR = 0.;

            for (int s = 0; s < kNSrc; ++s) {
                pCR += frac_s_CR[s]->GetBinContent(ix, iy)
                     * eps_s[s]->GetBinContent(ix, iy);

                pSR += frac_s_SR[s]->GetBinContent(ix, iy)
                     * eps_s[s]->GetBinContent(ix, iy);
            }

            eps_pred_CR->SetBinContent(ix, iy, pCR);
            eps_pred_SR->SetBinContent(ix, iy, pSR);

            double fc = 0.;

            if (tot_CR < kMinStat || tot_SR < kMinStat || (pCR + pSR) <= 0.) {
                fc = 0.;
                ++nLowStat;
            } else {
                // Symmetric percent-difference: bounded in [-2, +2] when one of
                // pSR or pCR vanishes, no asymmetric artifact from putting CR
                // in the denominator. Avoids the previous (pSR-pCR)/pCR formula
                // which blew up for compositions wildly different between CR
                // and SR and silently saturated against a hard ±1 clamp.
                fc = 2.0 * (pSR - pCR) / (pSR + pCR);

                if (!std::isfinite(fc)) {
                    fc = 0.;
                }

                fc = std::min(std::max(fc, -2.0), 2.0);
                ++nValid;
            }

            frac_comp->SetBinContent(ix, iy, fc);
        }
    }

    std::cout << "\n=== Composition systematic summary ===\n"
              << "  Bins with valid composition variation: " << nValid << "\n"
              << "  Bins below stat threshold, set to 0:   " << nLowStat << "\n";

    // ------------------------------------------------------------------
    // 6. Per-bin printout
    // ------------------------------------------------------------------
    printf("\n  pT-bin       |eta|-bin       ");

    for (int s = 0; s < kNSrc; ++s) {
        printf("f_%s^CR  ", kSrcName[s]);
    }

    printf("| ");

    for (int s = 0; s < kNSrc; ++s) {
        printf("f_%s^SR  ", kSrcName[s]);
    }

    printf("| frac_comp\n");

    for (int ix = 1; ix <= nX; ++ix) {
        for (int iy = 1; iy <= nY; ++iy) {
            const double ptLo = frac_comp->GetXaxis()->GetBinLowEdge(ix);
            const double ptHi = frac_comp->GetXaxis()->GetBinUpEdge(ix);
            const double etLo = frac_comp->GetYaxis()->GetBinLowEdge(iy);
            const double etHi = frac_comp->GetYaxis()->GetBinUpEdge(iy);

            printf("  [%3.0f,%3.0f]  [%.1f,%.1f]   ",
                   ptLo, ptHi, etLo, etHi);

            for (int s = 0; s < kNSrc; ++s) {
                printf("%.3f    ", frac_s_CR[s]->GetBinContent(ix, iy));
            }

            printf("| ");

            for (int s = 0; s < kNSrc; ++s) {
                printf("%.3f    ", frac_s_SR[s]->GetBinContent(ix, iy));
            }

            printf("| %+.3f\n", frac_comp->GetBinContent(ix, iy));
        }
    }

    // ------------------------------------------------------------------
    // 7. Save
    // ------------------------------------------------------------------
    TFile* fOut = TFile::Open(
        "outputs/composition_systematic.root",
        "RECREATE"
    );

    if (!fOut || fOut->IsZombie()) {
        std::cerr << "ERROR: cannot create outputs/composition_systematic.root\n";
        return;
    }

    frac_comp->Write("frac_comp");
    eps_pred_CR->Write("eps_pred_CR");
    eps_pred_SR->Write("eps_pred_SR");

    for (int s = 0; s < kNSrc; ++s) {
        eps_s[s]->Write();
        frac_s_CR[s]->Write();
        frac_s_SR[s]->Write();

        h_tight_CR[s]->Write();
        h_loose_CR[s]->Write();
        h_tight_SR[s]->Write();
        h_loose_SR[s]->Write();
    }

    fOut->Close();

    std::cout << "\nWrote outputs/composition_systematic.root\n"
              << "Detailed composition categories saved:\n";

    for (int s = 0; s < kNSrc; ++s) {
        std::cout << "  eps_" << kSrcName[s] << "\n"
                  << "  frac_" << kSrcName[s] << "_CR\n"
                  << "  frac_" << kSrcName[s] << "_SR\n";
    }
}

