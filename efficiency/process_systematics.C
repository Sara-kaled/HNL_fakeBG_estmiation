#include "TFile.h"
#include "TH2D.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TGraphAsymmErrors.h"
#include "TROOT.h"
#include <iostream>
#include <string>
#include <vector>
#include <memory>

// Helper function to make the plots
void make_syst_plot(TH2D* h_nominal, TH2D* h_syst_relative, const std::string& title, const std::string& output_filename);

// Main function
void process_systematics() {
    // --- Setup: Basic styling and prevent ROOT from auto-deleting histograms ---
    gROOT->SetStyle("gStyle");
    gStyle->SetOptStat(0);
    TH1::SetDefaultSumw2();
    TH1::AddDirectory(kFALSE); // CRITICAL: Prevents histograms from being deleted when the file is closed.

    // --- 1. Load the input file and all necessary histograms ---
    std::cout << "--> Loading histograms from fake_rates.root..." << std::endl;
    auto inputFile = std::make_unique<TFile>("outputs/fake_rates.root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open outputs/fake_rates.root" << std::endl;
        return;
    }

    // Nominal Histograms
    TH2D* h_fake_nominal = (TH2D*)inputFile->Get("fake_rate_2d_CR");
    TH2D* h_eff_nominal = (TH2D*)inputFile->Get("eff_2d_mc_CR");

    // Systematic Variation Histograms
    TH2D* h_fake_met_low = (TH2D*)inputFile->Get("fake_rate_2d_CR_MET_low");
    TH2D* h_fake_met_high = (TH2D*)inputFile->Get("fake_rate_2d_CR_MET_high");
    TH2D* h_eff_syst_sr = (TH2D*)inputFile->Get("eff_2d_mc_SR");

    if (!h_fake_nominal || !h_eff_nominal || !h_fake_met_low || !h_fake_met_high || !h_eff_syst_sr) {
        std::cerr << "Error: One or more required histograms not found in the file." << std::endl;
        return;
    }

    // --- 2. Calculate Relative Systematic Uncertainties ---
    std::cout << "--> Calculating relative systematic uncertainties..." << std::endl;

    // A) MET systematic on Fake Rate (epsilon_f)
    TH2D* h_fake_syst_met = (TH2D*)h_fake_nominal->Clone("h_fake_syst_met");
    h_fake_syst_met->SetTitle("Relative MET Syst. on Fake Rate (#epsilon_{f});Lepton p_{T} [GeV];Lepton |#eta|");

    for (int i = 1; i <= h_fake_nominal->GetNbinsX(); ++i) {
        for (int j = 1; j <= h_fake_nominal->GetNbinsY(); ++j) {
            double nominal_val = h_fake_nominal->GetBinContent(i, j);
            if (nominal_val <= 0) { // Avoid division by zero or nonsensical values
                h_fake_syst_met->SetBinContent(i, j, 0.0);
                continue;
            }
            double up_val = h_fake_met_high->GetBinContent(i, j);
            double down_val = h_fake_met_low->GetBinContent(i, j);

            // Take the larger of the two deviations (envelope method)
            double syst_abs = std::max(std::abs(up_val - nominal_val), std::abs(nominal_val - down_val));
            double syst_rel = syst_abs / nominal_val;
            h_fake_syst_met->SetBinContent(i, j, syst_rel);
        }
    }
    
    // B) Composition systematic on Real Efficiency (epsilon_r)
    TH2D* h_eff_syst_comp = (TH2D*)h_eff_nominal->Clone("h_eff_syst_comp");
    h_eff_syst_comp->SetTitle("Relative Composition Syst. on Real Efficiency (#epsilon_{r});Lepton p_{T} [GeV];Lepton |#eta|");
    
    // The systematic is the difference between the SR and CR efficiencies
    h_eff_syst_comp->Add(h_eff_syst_sr, -1.0); // h_eff_syst_comp = h_eff_nominal - h_eff_syst_sr
    h_eff_syst_comp->Divide(h_eff_nominal);    // Divide by nominal to get relative uncertainty
    
    // Take absolute value for visualization
    for (int i = 1; i <= h_eff_syst_comp->GetNbinsX(); ++i) {
        for (int j = 1; j <= h_eff_syst_comp->GetNbinsY(); ++j) {
            h_eff_syst_comp->SetBinContent(i, j, std::abs(h_eff_syst_comp->GetBinContent(i, j)));
        }
    }

    // --- 3. Generate and Save Plots ---
    std::cout << "--> Generating and saving plots..." << std::endl;
    
    // Plot for Fake Rate MET Systematic
    make_syst_plot(h_fake_nominal, h_fake_syst_met, "Fake Rate (#epsilon_{f}) with MET Systematic;Lepton p_{T} [GeV];Fake Rate", "outputs/syst_plot_fake_rate_met.png");

    // Plot for Real Efficiency Composition Systematic
    make_syst_plot(h_eff_nominal, h_eff_syst_comp, "Real Efficiency (#epsilon_{r}) with Composition Systematic;Lepton p_{T} [GeV];Real Efficiency", "outputs/syst_plot_real_eff_composition.png");
    
    std::cout << "--> Finished. Plots saved to current directory." << std::endl;
    
    // --- 4. (Optional) Write the calculated systematics to a new file for later use ---
    auto outputFile = std::make_unique<TFile>("outputs/systematics_processed.root", "RECREATE");
    h_fake_syst_met->Write();
    h_eff_syst_comp->Write();
    outputFile->Close();
    std::cout << "--> Relative systematics saved to systematics_processed.root" << std::endl;

}

// Helper function to create a 1D plot of a rate/efficiency with its systematic uncertainty band
void make_syst_plot(TH2D* h_nominal, TH2D* h_syst_relative, const std::string& title, const std::string& output_filename) {
    TCanvas* c = new TCanvas(output_filename.c_str(), "canvas", 800, 700);
    c->SetMargin(0.15, 0.05, 0.15, 0.08);

    // Project the 2D histograms to 1D (pT) to visualize
    TH1D* h_nominal_pt = h_nominal->ProjectionX((std::string(h_nominal->GetName()) + "_px").c_str());
    h_nominal_pt->SetTitle(title.c_str());
    h_nominal_pt->GetYaxis()->SetTitleOffset(1.5);
    h_nominal_pt->SetLineWidth(2);
    h_nominal_pt->SetLineColor(kBlack);

    // Create a TGraphAsymmErrors for the uncertainty band
    TGraphAsymmErrors* g_syst_band = new TGraphAsymmErrors(h_nominal_pt);

    for (int i = 1; i <= h_nominal_pt->GetNbinsX(); ++i) {
        double nominal_val = h_nominal_pt->GetBinContent(i);
        double rel_syst_avg = 0.0;
        int n_eta_bins = 0;

        // Average the relative uncertainty over the eta bins for this pT bin
        for (int j = 1; j <= h_syst_relative->GetNbinsY(); ++j) {
            rel_syst_avg += h_syst_relative->GetBinContent(i, j);
            n_eta_bins++;
        }
        if (n_eta_bins > 0) {
            rel_syst_avg /= n_eta_bins;
        }

        double abs_syst = nominal_val * rel_syst_avg;
        g_syst_band->SetPointError(i - 1, h_nominal_pt->GetBinWidth(i) / 2.0, h_nominal_pt->GetBinWidth(i) / 2.0, abs_syst, abs_syst);
    }
    
    g_syst_band->SetFillColorAlpha(kRed, 0.35);
    g_syst_band->SetFillStyle(1001);
    g_syst_band->SetMarkerSize(0);

    h_nominal_pt->Draw("AXIS");
    g_syst_band->Draw("2 SAME"); // "2" draws the rectangle for the band
    h_nominal_pt->Draw("HIST SAME"); // Draw the line on top

    TLegend* legend = new TLegend(0.6, 0.75, 0.9, 0.9);
    legend->AddEntry(h_nominal_pt, "Nominal Value", "l");
    legend->AddEntry(g_syst_band, "Systematic Uncertainty", "f");
    legend->SetBorderSize(0);
    legend->Draw();
    
    c->SaveAs(output_filename.c_str());
    delete c;
}

