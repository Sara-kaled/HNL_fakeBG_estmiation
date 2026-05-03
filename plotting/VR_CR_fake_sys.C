#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TLatex.h>
#include <TColor.h>
#include <TLine.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <tuple>
#include <algorithm> // For std::max
#include <memory>    // For std::unique_ptr
#include "../utils/SelectionUtils.h"

using namespace std;
using RVecF = ROOT::VecOps::RVec<float>;
using RVecI = ROOT::VecOps::RVec<int>;
using RVecC = ROOT::VecOps::RVec<char>;

// ==========================================================================================
// 1. ATLAS STYLE & PLOT CONFIGURATION
// ==========================================================================================

// --- Struct to hold info for each physics process ---
struct Process {
    string name;
    string legend;
    vector<string> files;
    int color;
    bool is_data = false;
};

// --- Struct to hold histograms for a single variable ---
struct HistSet {
    map<string, TH1F*> histograms;
    string xlabel;
    string title;
};

// --- Standard ATLAS Style ---
void SetATLASStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.16);
    gStyle->SetPadLeftMargin(0.16);
    gStyle->SetTitleSize(0.05, "xy");
    gStyle->SetLabelSize(0.045, "xy");
    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFillStyle(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
}

// --- Rewritten plotting function for stacked histograms with ratio plot ---
void DrawPlot(TCanvas* c, HistSet& hs,
              map<string, Process>& processes,
              bool logY = false)
 {
    c->cd();

    // Create two pads: upper for main plot, lower for ratio
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetTopMargin(0.07);
    pad1->Draw();
    pad1->cd();
    if (logY) pad1->SetLogy();

    // === Main Plot ===
    THStack *stack = new THStack("stack", "");
    TH1F* h_bkg_total = nullptr;

    // Changed stacking order: largest background at the bottom
    vector<string> stack_order = {"ttbar", "ztautau", "wjets", "singletop", "diboson", "higgs","fake"};

    for (const auto& proc_name : stack_order) {
        if (hs.histograms.count(proc_name)) {
            TH1F* h = hs.histograms[proc_name];
            h->SetFillColor(processes[proc_name].color);
            h->SetLineColor(kBlack);
            h->SetLineWidth(1);
            stack->Add(h);

            if (!h_bkg_total) {
                h_bkg_total = (TH1F*)h->Clone((string(h->GetName()) + "_total").c_str());
            } else {
                h_bkg_total->Add(h);
            }
        }
    }

    // Draw the stack
    stack->Draw("HIST");

    // === Total Background Statistical Uncertainty Band ===
    TH1F* h_bkg_total_stat = nullptr;
    if (h_bkg_total) {
        h_bkg_total_stat = (TH1F*)h_bkg_total->Clone("h_bkg_total_stat");
        h_bkg_total_stat->SetFillColor(kGray + 1);
        h_bkg_total_stat->SetFillStyle(3244);   // hatched pattern
        h_bkg_total_stat->SetMarkerSize(0);
        h_bkg_total_stat->SetLineWidth(0);
        h_bkg_total_stat->Draw("E2 SAME");      // Draw error band
    }

    // Draw Data on top
    TH1F* h_data = hs.histograms["data"];
    h_data->SetLineColor(kBlack);
    h_data->SetMarkerColor(kBlack);
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1.0);
    h_data->Draw("PE SAME");

    // Y-axis range
    double ymax = std::max(h_data->GetMaximum(), stack->GetMaximum());
    stack->SetMaximum(ymax * (logY ? 50.0 : 1.5));
    stack->SetMinimum(logY ? 0.1 : 0);

    stack->GetYaxis()->SetTitle("Events");
    stack->GetYaxis()->SetTitleOffset(1.5);
    stack->GetXaxis()->SetLabelSize(0);   // Hide x-labels on main plot

    // === Legend ===
    TLegend *leg = new TLegend(0.65, 0.65, 0.94, 0.92);
    leg->AddEntry(h_data, "Data", "PEL");
    if (h_bkg_total_stat) leg->AddEntry(h_bkg_total_stat, "Stat. unc.", "f");

    for (const auto& proc_name : stack_order) {
        if (hs.histograms.count(proc_name)) {
            leg->AddEntry(hs.histograms[proc_name], processes[proc_name].legend.c_str(), "F");
        }
    }
    leg->SetTextSize(0.04);
    leg->Draw();

    // === ATLAS Labels ===
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(72);
    latex.SetTextSize(0.06);
    latex.DrawLatex(0.2, 0.85, "ATLAS");
    latex.SetTextFont(42);
    latex.DrawLatex(0.34, 0.85, "Internal");
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.2, 0.78, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
    latex.DrawLatex(0.2, 0.72, hs.title.c_str());

    c->cd();  // Go back to main canvas

    // === Ratio Plot ===
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0.05);
    pad2->SetBottomMargin(0.35);
    pad2->Draw();
    pad2->cd();

    TH1F *h_ratio = (TH1F*)h_data->Clone("h_ratio");

    if (h_bkg_total && h_bkg_total->Integral() > 0) {
        // Proper error propagation for ratio
        for (int bin = 1; bin <= h_ratio->GetNbinsX(); ++bin) {
            double data_val = h_data->GetBinContent(bin);
            double bkg_val  = h_bkg_total->GetBinContent(bin);
            double bkg_err  = h_bkg_total->GetBinError(bin);

            if (bkg_val > 0) {
                double ratio_val = data_val / bkg_val;
                h_ratio->SetBinContent(bin, ratio_val);

                double rel_err = (data_val > 0) ?
                                 sqrt(1.0/data_val + pow(bkg_err/bkg_val, 2)) : 0.0;
                h_ratio->SetBinError(bin, rel_err * ratio_val);
            } else {
                h_ratio->SetBinContent(bin, 0);
                h_ratio->SetBinError(bin, 0);
            }
        }
    }

    h_ratio->GetYaxis()->SetRangeUser(0.5, 1.5);
    h_ratio->GetYaxis()->SetTitle("Data / Bkg");
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(0.12);
    h_ratio->GetYaxis()->SetLabelSize(0.1);
    h_ratio->GetYaxis()->SetTitleOffset(0.5);

    h_ratio->GetXaxis()->SetTitle(hs.xlabel.c_str());
    h_ratio->GetXaxis()->SetTitleSize(0.12);
    h_ratio->GetXaxis()->SetLabelSize(0.1);
    h_ratio->GetXaxis()->SetTitleOffset(1.2);

    // MC stat uncertainty band in ratio
    TH1F* h_ratio_band = nullptr;
    if (h_bkg_total) {
        h_ratio_band = (TH1F*)h_ratio->Clone("h_ratio_band");
        for (int bin = 1; bin <= h_ratio_band->GetNbinsX(); ++bin) {
            double bkg_val = h_bkg_total->GetBinContent(bin);
            double bkg_err = h_bkg_total->GetBinError(bin);
            if (bkg_val > 0) {
                h_ratio_band->SetBinError(bin, bkg_err / bkg_val);
            }
        }
        h_ratio_band->SetFillColor(kGray + 1);
        h_ratio_band->SetFillStyle(3244);
        h_ratio_band->Draw("E2 SAME");
    }

    h_ratio->Draw("PE SAME");

    // Line at y=1
    TLine *line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0,
                            h_ratio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kRed);
    line->SetLineStyle(2);
    line->Draw();

    c->Update();
}
// ==========================================================================================
// 2. MAIN ANALYSIS FUNCTION
// ==========================================================================================
void VR_CR_fake_sys() {
    SetATLASStyle();
    cout << "\n=== Initializing Analysis ===" << endl;

    // --- Define all processes, files, and styles ---
    map<string, Process> processes;
    string base_path = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    processes["data"] = {"data", "Data", {"data_with_fake_weights.root"}, kBlack, true};
    processes["fake"] = {"fake", "Fake (data-driven)", {}, kAzure - 9, false};
    processes["ttbar"] = {"ttbar", "t#bar{t}", {base_path + "ttbar/*.root"}, TColor::GetColor("#8E44AD")}; // Purple
    processes["ztautau"] = {"ztautau", "Z#rightarrow#tau#tau", {base_path + "ztautau/*.root"}, TColor::GetColor("#F39C12")}; // Orange
    processes["wjets"] = {"wjets", "W+jets", {base_path + "wjets/*.root"}, TColor::GetColor("#2ECC71")}; // Green
    processes["singletop"] = {"singletop", "Single Top", {base_path + "singletop/*.root"}, TColor::GetColor("#E74C3C")}; // Red
    processes["diboson"] = {"diboson", "Diboson", {base_path + "diboson/*.root"}, TColor::GetColor("#3498DB")}; // Blue
   processes["higgs"] = {"higgs", "Higgs", {base_path + "higgs/*.root"}, TColor::GetColor("#FF2D95")}; // Magenta/Pink
    vector<string> process_names = {"data", "ttbar", "ztautau", "wjets", "singletop", "diboson", "higgs","fake"};

    // --- Define variables to plot ---
    vector<tuple<string, int, double, double, string>> var_specs = {
        {"MET", 20, 0, 200, "Missing E_{T} [GeV]"},
        {"mll", 25, 0, 250, "m_{#mu e} [GeV]"},
        {"l1_pt", 25, 25, 150, "Muon p_{T} [GeV]"},
        {"l2_pt", 25, 25, 150, "Electron p_{T} [GeV]"},
        {"l2_eta", 10, 0, 2.5, "Electron |#eta|"}
    };

    // --- Create all histograms for all regions ---
    map<string, HistSet> histos_VR, histos_VR4;

  auto create_hists = [&](map<string, HistSet>& region_histos,
                        const string& tag, const string& title) {
    for (const auto& spec : var_specs) {
        string var_name = get<0>(spec);
        int    nbins    = get<1>(spec);
        double xmin     = get<2>(spec);
        double xmax     = get<3>(spec);

        // nominal histograms for each process (statistical band only)
        for (const auto& proc_name : process_names) {
            region_histos[var_name].histograms[proc_name] =
                new TH1F((proc_name + "_" + var_name + "_" + tag).c_str(),
                         "", nbins, xmin, xmax);
            region_histos[var_name].histograms[proc_name]->Sumw2();
        }
        region_histos[var_name].xlabel = get<4>(spec);
        region_histos[var_name].title  = title;
    }
};

create_hists(histos_VR,  "VR",  "Preselections");
create_hists(histos_VR4, "VR4", "VR: (m_{jj} < 400 GeV)");


    // --- Map for cutflow counters ---
    map<string, double> mc_yields;
    map<string, double> data_yields;
    vector<string> cut_names = {"Total", "Flavor", "Trigger", "pT_cuts", "Preselection", "OS_VR", "VR4"};
    for(const auto& cut : cut_names) {
        mc_yields[cut] = 0.0;
        data_yields[cut] = 0.0;
    }

    // ==========================================================================================
    // 3. MC EVENT LOOP
    // ==========================================================================================
    cout << "\n=== Processing MC Samples ===" << endl;
    for (const auto& proc_pair : processes) {
        string proc_name = proc_pair.first;
        if (proc_name == "data" || proc_name == "fake") continue;

        cout << "... processing " << proc_name << endl;
        TChain *tree = new TChain("tree");
        for (const auto& file : proc_pair.second.files) {
            tree->Add(file.c_str());
        }

        // --- Set branch addresses for MC via EventData ---
        EventData ev;
        SetupBranches(tree, ev);

        Long64_t nMC = tree->GetEntries();
        for (Long64_t i = 0; i < nMC; i++) {
            tree->GetEntry(i);
            ComputeDerived(ev);

            double total_weight_tight = mcWeightTight(ev);

            double mc_scale_factor = 1;
            if (proc_name == "ttbar") {
                mc_scale_factor = 1;
            } else if (proc_name == "ztautau") {
                mc_scale_factor = 1;
            }
            double w = mc_scale_factor * total_weight_tight;
            if (!isfinite(ev.total_weight)) continue;

            mc_yields["Total"] += w;

            // --- Event Selection ---
            if (!passTopology(ev)) continue;
            mc_yields["Flavor"] += w;

            if (!passTrigger(ev)) continue;
            mc_yields["Trigger"] += w;

            if (!passLeptonPt(ev)) continue;
            mc_yields["pT_cuts"] += w;

            if (!isMCPromptEvent(ev)) continue;

            if (!ev.mu_is_tight) continue;

            if (!passMuonIP(ev)) continue;
            if (!passElectronIP(ev)) continue;
            if (!passTriggerMatching(ev)) continue;
            if (!passCrackVeto(ev)) continue;

            mc_yields["Preselection"] += w;

            // --- Fill Regions ---
            if (ev.is_os) {
                mc_yields["OS_VR"] += w;
                for (auto const& [var_name, hist_set] : histos_VR) {
                    float val = 0;
                    if (var_name == "MET") val = ev.MET / 1000.0;
                    else if (var_name == "mll") val = ev.mll / 1000.0;
                    else if (var_name == "l1_pt") val = ev.l1_pt / 1000.0;
                    else if (var_name == "l2_pt") val = ev.l2_pt / 1000.0;
                    else if (var_name == "l2_eta") val = abs(ev.l2_eta);
                    hist_set.histograms.at(proc_name)->Fill(val, w);
                }
            }

            bool jet_cuts = passVBFJets(ev);   // j1>35 GeV, j2>20 GeV, |Δη|>3
            bool mjj_inverted = ev.mjj > 0 && ev.mjj < 400000.0;

            if (ev.is_os && jet_cuts && mjj_inverted && ev.no_bjets) {
                mc_yields["VR4"] += w;
                for (auto const& [var_name, hist_set] : histos_VR4) {
                    float val = 0;
                    if (var_name == "MET") val = ev.MET / 1000.0;
                    else if (var_name == "mll") val = ev.mll / 1000.0;
                    else if (var_name == "l1_pt") val = ev.l1_pt / 1000.0;
                    else if (var_name == "l2_pt") val = ev.l2_pt / 1000.0;
                    else if (var_name == "l2_eta") val = abs(ev.l2_eta);
                    hist_set.histograms.at(proc_name)->Fill(val, w);
                }
            }
        }
        delete tree;
    }

    // ==========================================================================================
    // 4. DATA & FAKE EVENT LOOP
    // ==========================================================================================
    cout << "\n=== Processing Data Sample ===" << endl;
    TChain *dataTree = new TChain("tree");
    dataTree->Add("../AMM/data_with_electron_fake_weights.root");

    // --- Set branch addresses for Data via EventData ---
    EventData evData;
    ResetTrigMatchPointers(evData);
    SetupBranches(dataTree, evData);

    Long64_t nData = dataTree->GetEntries();
    for (Long64_t i = 0; i < nData; i++) {
        dataTree->GetEntry(i);
        ComputeDerived(evData);

        data_yields["Total"] += 1.0;

        if (!passTopology(evData)) continue;
        data_yields["Flavor"] += 1.0;

        if (!passTrigger(evData)) continue;
        data_yields["Trigger"] += 1.0;

        if (!passLeptonPt(evData)) continue;
        data_yields["pT_cuts"] += 1.0;

        if (!evData.mu_is_tight) continue;
        if (!passMuonIP(evData)) continue;
        if (!passElectronIP(evData)) continue;
        if (!passTriggerMatching(evData)) continue;
        if (!passCrackVeto(evData)) continue;

        data_yields["Preselection"] += 1.0;

        double w_fake = evData.fw_nom;

        if (evData.is_os) {
    data_yields["OS_VR"] += 1.0;

    for (auto& [var_name, hist_set] : histos_VR) {
        float val = 0;
        if      (var_name == "MET")   val = evData.MET / 1000.0;
        else if (var_name == "mll")   val = evData.mll / 1000.0;
        else if (var_name == "l1_pt") val = evData.l1_pt / 1000.0;
        else if (var_name == "l2_pt") val = evData.l2_pt / 1000.0;
        else if (var_name == "l2_eta")val = fabs(evData.l2_eta);

        // data points
        hist_set.histograms.at("data")->Fill(val, 1.0);

        // fake prediction & its variations
        // (fake-syst variations removed — keep only statistical band)

        // If you also want to stack the nominal fake as its own process:
        hist_set.histograms.at("fake")->Fill(val, evData.fw_nom);
    }
}

        bool jet_cuts = passVBFJets(evData);   // j1>35 GeV, j2>20 GeV, |Δη|>3
        bool mjj_inverted = evData.mjj > 0 && evData.mjj < 400000.0;

      if (evData.is_os && jet_cuts && mjj_inverted && evData.no_bjets) {
    data_yields["VR4"] += 1.0;

    for (auto& [var_name, hist_set] : histos_VR4) {
        float val = 0;
        if      (var_name == "MET")   val = evData.MET / 1000.0;
        else if (var_name == "mll")   val = evData.mll / 1000.0;
        else if (var_name == "l1_pt") val = evData.l1_pt / 1000.0;
        else if (var_name == "l2_pt") val = evData.l2_pt / 1000.0;
        else if (var_name == "l2_eta")val = fabs(evData.l2_eta);

        hist_set.histograms.at("data")->Fill(val, 1.0);
        // (fake-syst variations removed — keep only statistical band)

        hist_set.histograms.at("fake")->Fill(val, evData.fw_nom);
    }
}
}

    // ==========================================================================================
    // 5. PRINT YIELDS & CREATE PLOTS
    // ==========================================================================================
// ==========================================================================================
// 5. PRINT YIELDS & CREATE PLOTS
// ==========================================================================================
cout << "\n=== Final Yields & Cutflow Summary ===\n";

// Compute totals
//double fake_vr   = histos_VR["MET"].histograms["fake"]->Integral();
double total_vr  = mc_yields["OS_VR"] ;
//double fake_vr4  = histos_VR4["MET"].histograms["fake"]->Integral();
double total_vr4 = mc_yields["VR4"] ;

// Table header
printf("+--------------------+-----------------------+-----------------------+\n");
printf("| %-18s | %-21s | %-21s |\n", "Cut", "MC (Weighted)", "Data (Events)");
printf("+--------------------+-----------------------+-----------------------+\n");

// Print cutflow
for (const auto& cut : cut_names) {
    printf("| %-18s | %21.1f | %21.0f |\n",
           cut.c_str(), mc_yields[cut], data_yields[cut]);
}

// Summary lines
printf("+--------------------+-----------------------+-----------------------+\n");
//printf("| %-18s | %21.1f | %21s |\n", "Fake (VR)",   fake_vr,   "-");
printf("| %-18s | %21.1f | %21.0f |\n", "Total Pred. (VR)", total_vr,  data_yields["OS_VR"]);
//printf("| %-18s | %21.1f | %21s |\n", "Fake (VR4)",  fake_vr4,  "-");
printf("| %-18s | %21.1f | %21.0f |\n", "Total Pred. (VR4)",total_vr4, data_yields["VR4"]);
printf("+--------------------+-----------------------+-----------------------+\n");




system("mkdir -p outputs/04_vr_cr_systematics");
TFile* outFile = new TFile("outputs/04_vr_cr_systematics/vr_cr_plots.root", "RECREATE");

// VR
for (auto& pair : histos_VR) {
    const string& var = pair.first;
    HistSet& hs = pair.second;

    TCanvas *c = new TCanvas(("c_" + var + "_VR").c_str(), "", 800, 800);
    DrawPlot(c, hs, processes, /*logY=*/true);
    c->SaveAs(("outputs/04_vr_cr_systematics/presel_" + var + ".pdf").c_str());
    c->Write();
}

// VR4
for (auto& pair : histos_VR4) {
    const string& var = pair.first;
    HistSet& hs = pair.second;

    TCanvas *c = new TCanvas(("c_" + var + "_VR4").c_str(), "", 800, 800);
    DrawPlot(c, hs, processes, /*logY=*/true);
    c->SaveAs(("outputs/04_vr_cr_systematics/VR_" + var + ".pdf").c_str());
    c->Write();
}


}
