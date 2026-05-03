//=============================================================================
// Filename: plot_fake_rates.C
// Description: ROOT macro to calculate and plot background-subtracted
//              fake lepton efficiencies using a custom RatePlotter class.
//
//              This macro operates in two phases:
//              Phase 1: Loads MC and Data histograms, performs real-subtraction
//                       (Data - MC) for tight (numerator) and loose (denominator)
//                       histograms, and saves these processed histograms
//                       to an intermediate ROOT file.
//              Phase 2: Initializes the RatePlotter class, directs it to the
//                       intermediate file, and calls its plotting functions
//                       to generate 1D and 2D fake efficiency plots.
//
// Author: Anara Research Workflow
// Date: February 25, 2026
//=============================================================================

// ROOT and C++ Standard Library Includes
#include <iostream>  // For standard input/output operations
#include <string>    // For std::string type
#include <vector>    // For std::vector container

// ROOT Specific Headers
#include "TROOT.h"       // For gROOT->LoadMacro
#include "TStyle.h"      // For gStyle (if needed, although AtlasStyle usually sets this)
#include "TCanvas.h"     // For drawing canvases
#include "TH1D.h"        // For 1D histograms
#include "TH2D.h"        // For 2D histograms
#include "TFile.h"       // For ROOT file operations
#include "TString.h"     // For ROOT's TString utility
#include "TLegend.h"     // For plot legends
#include "TEfficiency.h" // For efficiency calculations
#include "TSystem.h"     // For gSystem->mkdir (creating directories)

// Custom RatePlotter Class Inclusion
// IMPORTANT: The #include for RatePlotter.h MUST come AFTER gROOT->LoadMacro()
// when compiling with ACLiC (root -l -b -q your_macro.C).
// This ensures the class definition is available after the shared library is loaded.
// Adjust the path below if your RatePlotter.cxx/h are in a different subdirectory.
void plot_fake_rates() {
    // ------------------------------------------------------------------
    // 1. Compile and Load RatePlotter Class
    // ------------------------------------------------------------------
    std::cout << "INFO: Loading RatePlotter utility..." << std::endl;
    gROOT->LoadMacro("RatePlotter.cxx++");
    #include "RatePlotter.h" // Include after loading macro

    // ------------------------------------------------------------------
    // 2. Configuration Variables
    //    Define all file paths, histogram names, and plot settings.
    // ------------------------------------------------------------------

    // Output directory for generated plots
    std::string out_dir_name = "fake_rate_plots";
    // Intermediate ROOT file to store background-subtracted histograms
    std::string subtracted_hist_file_name = "subtracted_histograms.root";

    // Full paths to your original MC and Data ROOT files on EOS
    TString mc_file_path = "root://eoshome-s.cern.ch//eos/user/s/sabdelha/freshness_merge/fakelepton_BG_study/MC_eff/MC_histograms.root";
    TString data_file_path = "root://eoshome-s.cern.ch//eos/user/s/sabdelha/freshness_merge/fakelepton_BG_study/data_eff/Data_histograms.root";

    // List of 1D kinematic variables to plot
    std::vector<TString> variables = {"pt", "d0", "d0sig", "met", "eta"};
    // Corresponding LaTeX-formatted X-axis labels for 1D plots
    std::vector<TString> var_x_labels = {
        "Lepton $p_{T}$ [GeV]",
        "Lepton $d_{0}$ [mm]",
        "Lepton $d_{0}/\\sigma(d_{0})$",
        "$E_{T}^{miss}$ [GeV]",
        "Lepton $\\eta$"
    };
    
    // List of analysis regions for which histograms exist
    std::vector<TString> regions = {"CR", "SR", "CR_MET_low", "CR_MET_high", "CR_el_d0", "incl"};
    // Corresponding display names for regions (for potential use in titles/labels)
    // (Note: The RatePlotter version you have mostly infers titles from histogram names)
    // std::vector<TString> region_display_names = {
    //     "Control Region", "Signal Region", "Control Region (Low $E_{T}^{miss}$)",
    //     "Control Region (High $E_{T}^{miss}$)", "Control Region ($d_{0}$ Electron)", "Inclusive"
    // };

    // ------------------------------------------------------------------
    // 3. Phase 1: Create Background-Subtracted Histograms
    //    This phase calculates (Data_Tight - MC_Tight) / (Data_Loose - MC_Loose)
    //    and saves the numerator and denominator histograms to a new file.
    // ------------------------------------------------------------------

    std::cout << "\n--- Phase 1: Creating background-subtracted histograms ---" << std::endl;
    
    // Create the output directory for plots if it doesn't exist
    // Note: RatePlotter::setOutDir() exists but doesn't have a getOutDir()
    // so we manage directory creation manually and then pass the name.
    gSystem->mkdir(out_dir_name.c_str(), true);

    // Open input MC and Data files
    TFile* f_mc = TFile::Open(mc_file_path);
    TFile* f_data = TFile::Open(data_file_path);

    // Create the intermediate output file for subtracted histograms
    TFile* f_out = new TFile(subtracted_hist_file_name.c_str(), "RECREATE");

    // Basic error checking for file opening
    if (!f_mc || f_mc->IsZombie()) {
        std::cerr << "ERROR: Could not open MC file: " << mc_file_path << std::endl;
        return;
    }
    if (!f_data || f_data->IsZombie()) {
        std::cerr << "ERROR: Could not open Data file: " << data_file_path << std::endl;
        return;
    }
    if (!f_out || f_out->IsZombie()) {
        std::cerr << "ERROR: Could not create output file: " << subtracted_hist_file_name << std::endl;
        f_mc->Close(); delete f_mc;
        f_data->Close(); delete f_data;
        return;
    }
    std::cout << "INFO: Successfully opened input files and created intermediate output file." << std::endl;

    // --- Process 1D Histograms (Subtract MC from Data) ---
    for (const auto& region : regions) {
        for (size_t j = 0; j < variables.size(); ++j) {
            const auto& var = variables[j]; // e.g., "pt", "d0"
            const auto& x_label = var_x_labels[j]; // e.g., "Lepton $p_{T}$ [GeV]"

            TString tight_hist_name_in = Form("h_tight_%s_%s", var.Data(), region.Data());
            TString loose_hist_name_in = Form("h_loose_%s_%s", var.Data(), region.Data());

            // Retrieve histograms from input files
            TH1D* h_tight_data = (TH1D*)f_data->Get(tight_hist_name_in);
            TH1D* h_loose_data = (TH1D*)f_data->Get(loose_hist_name_in);
            TH1D* h_tight_mc = (TH1D*)f_mc->Get(tight_hist_name_in);
            TH1D* h_loose_mc = (TH1D*)f_mc->Get(loose_hist_name_in);

            // Skip if any required histogram is not found for this combination
            if (!h_tight_data ||!h_loose_data ||!h_tight_mc ||!h_loose_mc) {
                std::cerr << "WARNING: Missing 1D histogram(s) for '" << var << "' in region '" << region << "'. Skipping." << std::endl;
                continue;
            }

            // Define output names for the subtracted histograms
            TString num_name_out = Form("num_%s_%s", var.Data(), region.Data()); // Numerator
            TString den_name_out = Form("den_%s_%s", var.Data(), region.Data()); // Denominator

            // Clone Data histograms and perform MC subtraction
            TH1D* h_num = (TH1D*)h_tight_data->Clone(num_name_out); // Numerator: Data(tight) - MC(tight)
            h_num->Add(h_tight_mc, -1.0); // Subtract MC (prompt background)

            TH1D* h_den = (TH1D*)h_loose_data->Clone(den_name_out); // Denominator: Data(loose) - MC(loose)
            h_den->Add(h_loose_mc, -1.0); // Subtract MC (prompt background)

            // Ensure non-negative bin content for robust efficiency calculation
            // Efficiencies are typically [0,1], so negative counts are unphysical after subtraction.
            for (int bin_i = 1; bin_i <= h_num->GetNbinsX(); ++bin_i) {
                if (h_num->GetBinContent(bin_i) < 0) h_num->SetBinContent(bin_i, 0);
                if (h_den->GetBinContent(bin_i) < 0) h_den->SetBinContent(bin_i, 0);
            }

            // Set titles and axis labels directly on these histograms.
            // RatePlotter will automatically use these when drawing.
            h_num->SetTitle(Form("Fake Efficiency (%s);%s;Fake Efficiency", region.Data(), x_label.Data()));
            h_den->SetTitle(Form("Fake Efficiency (%s);%s;Fake Efficiency", region.Data(), x_label.Data()));

            // Write the processed histograms to the intermediate file
            f_out->cd();
            h_num->Write();
            h_den->Write();

            // Clean up cloned histograms from memory
            delete h_num;
            delete h_den;
            std::cout << "INFO: Prepared 1D histograms for '" << var << "' in region '" << region << "'." << std::endl;
        }
    }

    // --- Process 2D Histograms (Subtract MC from Data) ---
    for (const auto& region : regions) {
        TString tight_2d_name_in = Form("h2_tight_eta_pt_%s", region.Data());
        TString loose_2d_name_in = Form("h2_loose_eta_pt_%s", region.Data());

        // Retrieve 2D histograms
        TH2D* h2_tight_data = (TH2D*)f_data->Get(tight_2d_name_in);
        TH2D* h2_loose_data = (TH2D*)f_data->Get(loose_2d_name_in);
        TH2D* h2_tight_mc = (TH2D*)f_mc->Get(tight_2d_name_in);
        TH2D* h2_loose_mc = (TH2D*)f_mc->Get(loose_2d_name_in);

        // Skip if any required histogram is not found
        if (!h2_tight_data ||!h2_loose_data ||!h2_tight_mc ||!h2_loose_mc) {
            std::cerr << "WARNING: Missing 2D histogram(s) for region '" << region << "'. Skipping." << std::endl;
            continue;
        }

        // Define output names for the subtracted 2D histograms
        TString num_2d_name_out = Form("num_eta_pt_%s", region.Data()); // Numerator
        TString den_2d_name_out = Form("den_eta_pt_%s", region.Data()); // Denominator

        // Clone Data 2D histograms and perform MC subtraction
        TH2D* h2_num = (TH2D*)h2_tight_data->Clone(num_2d_name_out);
        h2_num->Add(h2_tight_mc, -1.0);

        TH2D* h2_den = (TH2D*)h2_loose_data->Clone(den_2d_name_out);
        h2_den->Add(h2_loose_mc, -1.0);

        // Ensure non-negative bin content for 2D histograms
        for (int bin_x = 1; bin_x <= h2_num->GetNbinsX(); ++bin_x) {
            for (int bin_y = 1; bin_y <= h2_num->GetNbinsY(); ++bin_y) {
                if (h2_num->GetBinContent(bin_x, bin_y) < 0) h2_num->SetBinContent(bin_x, bin_y, 0);
                if (h2_den->GetBinContent(bin_x, bin_y) < 0) h2_den->SetBinContent(bin_x, bin_y, 0);
            }
        }
        
        // Set titles and axis labels for 2D plots
        h2_num->SetTitle(Form("Fake Efficiency Map (%s);Lepton #eta;Lepton p_{T} [GeV]", region.Data()));
        h2_den->SetTitle(Form("Fake Efficiency Map (%s);Lepton #eta;Lepton p_{T} [GeV]", region.Data()));

        // Write the processed 2D histograms to the intermediate file
        f_out->cd();
        h2_num->Write();
        h2_den->Write();

        // Clean up cloned 2D histograms
        delete h2_num;
        delete h2_den;
        std::cout << "INFO: Prepared 2D histograms for region '" << region << "'." << std::endl;
    }

    // Close all opened files in Phase 1
    f_out->Close(); delete f_out;
    f_mc->Close(); delete f_mc;
    f_data->Close(); delete f_data;
    std::cout << "INFO: Phase 1 complete. Subtracted histograms saved to: " << subtracted_hist_file_name << std::endl;


    // ------------------------------------------------------------------
    // 4. Phase 2: Plot Fake Rates using RatePlotter
    //    This phase uses the RatePlotter to read the background-subtracted
    //    histograms and generate the final plots.
    // ------------------------------------------------------------------

    std::cout << "\n--- Phase 2: Plotting efficiencies with RatePlotter ---" << std::endl;

    // Initialize RatePlotter instance for plotting
    RatePlotter rp_plot;
    
    // Configure RatePlotter with common settings
    rp_plot.setOutDir(out_dir_name);         // Output directory (expects std::string)
    rp_plot.setPrint(true);                 // Enable saving plots to file
    rp_plot.setFigureFormat("pdf");         // Output format for plots
    rp_plot.setStylePath("Utils/AtlasStyle.C"); // Path to the ATLAS style file
    rp_plot.setStyle(true);                 // Apply ATLAS style
    rp_plot.setLabel("ATLAS Internal", "label"); // Set ATLAS label
    rp_plot.setLumi(140.0);                 // Set integrated luminosity in fb^-1 (example value)

    // Add the intermediate file (containing background-subtracted histograms)
    // as the "data" source for RatePlotter.
    // RatePlotter will then look for "num_..." and "den_..." histograms within this file.
    rp_plot.addDataFile(subtracted_hist_file_name.c_str());

    // --- Plot 1D Efficiencies ---
    std::cout << "\nINFO: Generating 1D fake rate plots..." << std::endl;
    for (const auto& region : regions) {
        for (const auto& var : variables) {
            TString num_name = Form("num_%s_%s", var.Data(), region.Data());
            TString den_name = Form("den_%s_%s", var.Data(), region.Data());

            // Check if the histograms exist in the intermediate file before attempting to plot
            // (Opening and closing file for each check is inefficient, but safe for robust check)
            TFile* temp_check_file = TFile::Open(subtracted_hist_file_name.c_str(), "READ");
            if (temp_check_file && temp_check_file->Get(num_name) && temp_check_file->Get(den_name)) {
                // Call makeRatePlot with the names of the numerator and denominator histograms
                rp_plot.makeRatePlot(num_name, den_name);
                std::cout << "INFO: Generated 1D plot for variable '" << var << "' in region '" << region << "'." << std::endl;
            } else {
                std::cerr << "WARNING: Skipping 1D plot for '" << var << "' in region '" << region << "' (histograms not found in intermediate file)." << std::endl;
            }
            if (temp_check_file) {temp_check_file->Close(); delete temp_check_file;}
        }
    }

    // --- Plot 2D Efficiencies ---
    std::cout << "\nINFO: Generating 2D fake efficiency maps..." << std::endl;
    for (const auto& region : regions) {
        TString num_2d_name = Form("num_eta_pt_%s", region.Data());
        TString den_2d_name = Form("den_eta_pt_%s", region.Data());

        // Check if the histograms exist in the intermediate file
        TFile* temp_check_file = TFile::Open(subtracted_hist_file_name.c_str(), "READ");
        if (temp_check_file && temp_check_file->Get(num_2d_name) && temp_check_file->Get(den_2d_name)) {
            // Call makeRatePlot2D with histogram names and specify "data" as the source
            rp_plot.makeRatePlot2D(num_2d_name, den_2d_name, "data");
            std::cout << "INFO: Generated 2D plot for region '" << region << "'." << std::endl;
        } else {
            std::cerr << "WARNING: Skipping 2D plot for region '" << region << "' (histograms not found in intermediate file)." << std::endl;
        }
        if (temp_check_file) {temp_check_file->Close(); delete temp_check_file;}
    }

    std::cout << "\nSUCCESS: Finished generating fake rate plots. Check '" << out_dir_name << "/' directory for results." << std::endl;
}


