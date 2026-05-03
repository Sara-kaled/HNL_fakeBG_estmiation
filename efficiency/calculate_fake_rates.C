// =============================================================================
// FILE:       calculate_fake_rates.C
// PURPOSE:    Compute the data-driven fake electron rate map f(pT, |eta|).
//             f = P(tight | fake electron) = (data-MC)_tight / (data-MC)_loose
//             per pT x |eta| bin in the same-sign CR.
// INPUTS:     ../data_eff/Data_histograms.root  (loose/tight data histograms)
//             ../MC_eff/MC_histograms.root       (prompt-only MC histograms)
// OUTPUTS:    fake_rates.root  (2D and 1D fake rate histograms + MC prompt efficiency)
// RUNS AFTER: MC_eff/fake_bg_estimation.C, data_eff/fake_bg_estimation.C
// KEY PHYSICS: The numerator (data_tight - mc_tight) counts fake electrons passing the
//             tight WP in data after subtracting real contamination. The denominator
//             (data_loose - mc_loose) counts all fake electrons passing loose. Their
//             ratio is the fake efficiency ε_f used by the Matrix Method.
// =============================================================================
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <cstdio>
#include <iostream>
#include <vector>

void calculate_fake_rates() {
    // Open input files
    TFile* f_data = new TFile("outputs/Data_histograms.root", "READ");
    TFile* f_mc = new TFile("outputs/MC_histograms.root", "READ");

    if (!f_data || f_data->IsZombie() || !f_mc || f_mc->IsZombie()) {
        std::cerr << "Error opening input files!" << std::endl;
        return;
    }

    // Output file
    TFile* f_out = new TFile("outputs/fake_rates.root", "RECREATE");

    // List of regions
    std::vector<TString> regions = {"CR", "CR_MET_low", "CR_MET_high"};

    // List of 1D variables
    std::vector<TString> vars_1d = {"pt", "eta", "d0", "met"};

    // =====================================
    // Compute 2D fake rates for all regions
    // =====================================
    for (const auto& region : regions) {
        // Get histograms
        TH2D* h2_tight_data = (TH2D*)f_data->Get("h2_tight_eta_pt_" + region);
        TH2D* h2_loose_data = (TH2D*)f_data->Get("h2_loose_eta_pt_" + region);
        TH2D* h2_tight_mc = (TH2D*)f_mc->Get("h2_tight_eta_pt_" + region);
        TH2D* h2_loose_mc = (TH2D*)f_mc->Get("h2_loose_eta_pt_" + region);

        if (!h2_tight_data || !h2_loose_data || !h2_tight_mc || !h2_loose_mc) {
            std::cerr << "Missing 2D histograms for region: " << region << std::endl;
            continue;
        }

        // Compute numerator: data_tight - mc_tight
        TH2D* num = (TH2D*)h2_tight_data->Clone("num_" + region);
        num->Add(h2_tight_mc, -1.0);

        // Compute denominator: data_loose - mc_loose
        TH2D* den = (TH2D*)h2_loose_data->Clone("den_" + region);
        den->Add(h2_loose_mc, -1.0);

        // Compute fake rate: (data_tight - mc_tight) / (data_loose - mc_loose)
        // This gives f = P(tight | fake electron) in each pT x |eta| bin.
        TH2D* fake_rate_2d = (TH2D*)num->Clone("fake_rate_2d_" + region);
        fake_rate_2d->Divide(den);
        fake_rate_2d->SetTitle("2D Fake Rate (pt, eta) - " + region);

        // Guard against zero/negative denominator bins (data-MC subtraction can go negative).
        // Set problematic bins to 0 with a large error so they do not silently bias the result.
        int nWarn2D = 0;
        for (int ix = 1; ix <= den->GetNbinsX(); ++ix) {
            for (int iy = 1; iy <= den->GetNbinsY(); ++iy) {
                double d = den->GetBinContent(ix, iy);
                double r = fake_rate_2d->GetBinContent(ix, iy);
                if (d <= 0.0 || r < 0.0 || r > 1.0) {
                    fake_rate_2d->SetBinContent(ix, iy, 0.0);
                    fake_rate_2d->SetBinError(ix, iy, 1.0);
                    ++nWarn2D;
                }
            }
        }
        if (nWarn2D > 0)
            std::cerr << "WARNING: " << nWarn2D << " bad 2D bins in region " << region
                      << " (den<=0, rate<0, or rate>1) — set to 0 ± 1\n";

        // Save
        fake_rate_2d->Write();

        // -----------------------------------------------------------------
        // Per-bin ε_f printout for inspection (region-by-region).
        // After applying the charge-misID SF in MC_eff.C, compare these
        // numbers with your previous fake_rates.root output to see how the
        // additional charge-flip subtraction shifts the measured ε_f.
        // Expected direction: ε_f decreases (more prompt is subtracted).
        // -----------------------------------------------------------------
        printf("\n  ε_f per bin (region: %s) [after charge-misID SF on MC]\n",
               region.Data());
        printf("  ----------------------------------------------------------------\n");
        printf("  pT [GeV]            |eta|             ε_f          ± stat\n");
        printf("  ----------------------------------------------------------------\n");
        const int nXp = fake_rate_2d->GetNbinsX();
        const int nYp = fake_rate_2d->GetNbinsY();
        for (int ix = 1; ix <= nXp; ++ix) {
            const double ptLo = fake_rate_2d->GetXaxis()->GetBinLowEdge(ix);
            const double ptHi = fake_rate_2d->GetXaxis()->GetBinUpEdge(ix);
            for (int iy = 1; iy <= nYp; ++iy) {
                const double etLo = fake_rate_2d->GetYaxis()->GetBinLowEdge(iy);
                const double etHi = fake_rate_2d->GetYaxis()->GetBinUpEdge(iy);
                printf("  [%6.1f,%6.1f]  [%4.2f,%4.2f]   %8.4f   ± %.4f\n",
                       ptLo, ptHi, etLo, etHi,
                       fake_rate_2d->GetBinContent(ix, iy),
                       fake_rate_2d->GetBinError  (ix, iy));
            }
        }

        // Clean up temps
        delete num;
        delete den;
    }

    // =====================================
    // Compute 1D fake rates for each variable in each region
    // =====================================
    for (const auto& region : regions) {
        for (const auto& var : vars_1d) {
            // Get histograms
            TH1D* h_tight_data = (TH1D*)f_data->Get("h_tight_" + var + "_" + region);
            TH1D* h_loose_data = (TH1D*)f_data->Get("h_loose_" + var + "_" + region);
            TH1D* h_tight_mc = (TH1D*)f_mc->Get("h_tight_" + var + "_" + region);
            TH1D* h_loose_mc = (TH1D*)f_mc->Get("h_loose_" + var + "_" + region);

            if (!h_tight_data || !h_loose_data || !h_tight_mc || !h_loose_mc) {
                std::cerr << "Missing 1D histograms for var: " << var << " in region: " << region << std::endl;
                continue;
            }

            // Compute numerator: data_tight - mc_tight
            TH1D* num = (TH1D*)h_tight_data->Clone("num_1d_" + var + "_" + region);
            num->Add(h_tight_mc, -1.0);

            // Compute denominator: data_loose - mc_loose
            TH1D* den = (TH1D*)h_loose_data->Clone("den_1d_" + var + "_" + region);
            den->Add(h_loose_mc, -1.0);

            // Compute fake rate: num / den
            TH1D* fake_rate_1d = (TH1D*)num->Clone("fake_rate_1d_" + var + "_" + region);
            fake_rate_1d->Divide(den);
            fake_rate_1d->SetTitle("1D Fake Rate (" + var + ") - " + region);

            // Save
            fake_rate_1d->Write();

            // Clean up temps
            delete num;
            delete den;
        }
    }

    // =====================================
    // For MC only: 2D efficiency in CR and SR (tight / loose for prompt)
    // =====================================
    std::vector<TString> mc_regions = {"CR", "SR"};
    for (const auto& region : mc_regions) {
        // Get histograms (MC only, prompt)
        TH2D* h2_tight_mc = (TH2D*)f_mc->Get("h2_tight_eta_pt_" + region);
        TH2D* h2_loose_mc = (TH2D*)f_mc->Get("h2_loose_eta_pt_" + region);

        if (!h2_tight_mc || !h2_loose_mc) {
            std::cerr << "Missing MC 2D histograms for region: " << region << std::endl;
            continue;
        }

        // Compute eff: tight / loose
        TH2D* eff_2d_mc = (TH2D*)h2_tight_mc->Clone("eff_2d_mc_" + region);
        eff_2d_mc->Divide(h2_loose_mc);
        eff_2d_mc->SetTitle("2D Prompt Efficiency (pt, eta) MC - " + region);

        // Save
        eff_2d_mc->Write();
    }

    // Close files
    f_out->Write();
    f_out->Close();
    f_data->Close();
    f_mc->Close();

    std::cout << "All fake rates and efficiencies saved to fake_rates.root" << std::endl;
}
