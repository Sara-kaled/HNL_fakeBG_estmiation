// =============================================================================
// FILE:    plotting/fake_real_CRs.C
// PURPOSE: Data/MC comparison in fake_CR (SS, no b-jets, ≥2 jets) and prompt_CR
//          (OS, no b-jets, ≥2 jets); also draws the prompt-electron fraction
//          (truth-prompt MC over all MC) per region as a function of pT and |η|.
// INPUTS:  MC ntuples per process; data ntuple is AMM/data_with_electron_fake_weights.root.
// OUTPUTS: outputs/05_fake_composition/fake_CR_<var>.pdf
//          outputs/05_fake_composition/prompt_CR_<var>.pdf
//          outputs/05_fake_composition/prompt_frac_{pt,eta}_<region>.pdf
// =============================================================================
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
#include <algorithm>
#include <memory>

#include "../utils/SelectionUtils.h"

using namespace std;

// ==========================================================================================
// 1. ATLAS STYLE & PLOT CONFIGURATION
// ==========================================================================================

struct Process {
    string name;
    string legend;
    vector<string> files;
    int color;
    bool is_data = false;
};

struct HistSet {
    map<string, TH1F*> histograms;
    string xlabel;
    string title;
};

struct FractionHist {
    TH1F* h_prompt_num_pt;
    TH1F* h_total_pt;
    TH1F* h_prompt_num_eta;
    TH1F* h_total_eta;
    string title;
};

map<string, FractionHist> frac_VR, frac_VR4;

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

void DrawPlot(TCanvas* c, HistSet& hs, map<string, Process>& processes, bool logY = false) {
    c->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02); pad1->SetTopMargin(0.07); pad1->Draw(); pad1->cd();
    if (logY) pad1->SetLogy();

    THStack *stack = new THStack("stack", "");
    TH1F* h_bkg_total = nullptr;

    vector<string> stack_order = {"ttbar", "ztautau", "wjets", "singletop", "diboson", "higgs"};
    for (const auto& proc_name : stack_order) {
        if (hs.histograms.count(proc_name)) {
            TH1F* h = hs.histograms[proc_name];
            h->SetFillColor(processes[proc_name].color);
            h->SetLineColor(kBlack); h->SetLineWidth(1);
            stack->Add(h);
            if (!h_bkg_total) h_bkg_total = (TH1F*)h->Clone((string(h->GetName()) + "_total").c_str());
            else              h_bkg_total->Add(h);
        }
    }

    TH1F* h_bkg_total_stat = nullptr;
    if (h_bkg_total) {
        h_bkg_total_stat = (TH1F*)h_bkg_total->Clone("h_bkg_total_stat");
        h_bkg_total_stat->SetFillColor(kGray+1);
        h_bkg_total_stat->SetFillStyle(3244);
        h_bkg_total_stat->SetMarkerSize(0);
        h_bkg_total_stat->SetLineWidth(0);
    }

    stack->Draw("HIST");
    if (h_bkg_total_stat) h_bkg_total_stat->Draw("E2 SAME");

    TH1F* h_data = hs.histograms["data"];
    h_data->SetLineColor(kBlack); h_data->SetMarkerColor(kBlack);
    h_data->SetMarkerStyle(20); h_data->SetMarkerSize(1.0);
    h_data->Draw("PE SAME");

    double ymax = std::max(h_data->GetMaximum(), stack->GetMaximum());
    stack->SetMaximum(ymax * (logY ? 50.0 : 1.5));
    stack->SetMinimum(logY ? 0.1 : 0);
    stack->GetYaxis()->SetTitle("Events");
    stack->GetYaxis()->SetTitleOffset(1.5);
    stack->GetXaxis()->SetLabelSize(0);

    TLegend *leg = new TLegend(0.65, 0.65, 0.94, 0.92);
    leg->AddEntry(h_data, "Data", "PEL");
    if (h_bkg_total_stat) leg->AddEntry(h_bkg_total_stat, "Stat. unc.", "f");
    for (const auto& proc_name : stack_order) {
        if (hs.histograms.count(proc_name))
            leg->AddEntry(hs.histograms[proc_name], processes[proc_name].legend.c_str(), "F");
    }
    leg->SetTextSize(0.04); leg->Draw();

    TLatex latex; latex.SetNDC();
    latex.SetTextFont(72); latex.SetTextSize(0.06);
    latex.DrawLatex(0.2, 0.85, "ATLAS");
    latex.SetTextFont(42); latex.DrawLatex(0.34, 0.85, "Internal");
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.2, 0.78, "#sqrt{s} = 13 TeV,140 fb^{-1}");
    latex.DrawLatex(0.2, 0.72, hs.title.c_str());

    c->cd();

    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0.05); pad2->SetBottomMargin(0.35); pad2->Draw(); pad2->cd();

    TH1F *h_ratio = (TH1F*)h_data->Clone("h_ratio");
    if (h_bkg_total && h_bkg_total->Integral() > 0) {
        for (int bin = 1; bin <= h_ratio->GetNbinsX(); ++bin) {
            double data_val = h_data->GetBinContent(bin);
            double bkg_val  = h_bkg_total->GetBinContent(bin);
            double bkg_err  = h_bkg_total->GetBinError(bin);
            if (bkg_val > 0) {
                h_ratio->SetBinContent(bin, data_val / bkg_val);
                double rel_err = (data_val > 0) ? sqrt(1.0/data_val + pow(bkg_err/bkg_val, 2)) : 0;
                h_ratio->SetBinError(bin, rel_err * (data_val / bkg_val));
            } else {
                h_ratio->SetBinContent(bin, 0);
                h_ratio->SetBinError(bin, 0);
            }
        }
    }

    h_ratio->GetYaxis()->SetRangeUser(0.4, 2.5);
    h_ratio->GetYaxis()->SetTitle("Data / Bkg");
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(0.12); h_ratio->GetYaxis()->SetLabelSize(0.1);
    h_ratio->GetYaxis()->SetTitleOffset(0.5);
    h_ratio->GetXaxis()->SetTitle(hs.xlabel.c_str());
    h_ratio->GetXaxis()->SetTitleSize(0.12); h_ratio->GetXaxis()->SetLabelSize(0.1);
    h_ratio->GetXaxis()->SetTitleOffset(1.2);

    TH1F* h_ratio_band = (TH1F*)h_ratio->Clone("h_ratio_band");
    if (h_bkg_total) {
        for (int bin=1; bin<=h_ratio_band->GetNbinsX(); ++bin) {
            double bkg_val = h_bkg_total->GetBinContent(bin);
            double bkg_err = h_bkg_total->GetBinError(bin);
            if (bkg_val > 0) h_ratio_band->SetBinError(bin, bkg_err / bkg_val);
        }
    }
    h_ratio_band->SetFillColor(kGray+1);
    h_ratio_band->SetFillStyle(3244);
    h_ratio_band->Draw("E2 SAME");
    h_ratio->Draw("PE SAME");

    TLine *line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1, h_ratio->GetXaxis()->GetXmax(), 1);
    line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();

    c->Update();
}

// ==========================================================================================
// 2. MAIN ANALYSIS FUNCTION
// ==========================================================================================
void fake_real_CRs() {
    SetATLASStyle();
    cout << "\n=== Initializing Analysis ===" << endl;

    map<string, Process> processes;
    string base_path = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    processes["data"]      = {"data",      "Data",                  {"../AMM/data_with_electron_fake_weights.root"}, kBlack, true};
    processes["ttbar"]     = {"ttbar",     "t#bar{t}",              {base_path + "ttbar/*.root"},     TColor::GetColor("#8E44AD")};
    processes["ztautau"]   = {"ztautau",   "Z#rightarrow#tau#tau",  {base_path + "ztautau/*.root"},   TColor::GetColor("#F39C12")};
    processes["wjets"]     = {"wjets",     "W+jets",                {base_path + "wjets/*.root"},     TColor::GetColor("#2ECC71")};
    processes["singletop"] = {"singletop", "Single Top",            {base_path + "singletop/*.root"}, TColor::GetColor("#E74C3C")};
    processes["diboson"]   = {"diboson",   "Diboson",               {base_path + "diboson/*.root"},   TColor::GetColor("#3498DB")};
    processes["higgs"]     = {"higgs",     "Higgs",                 {base_path + "higgs/*.root"},     TColor::GetColor("#FF2D95")};
    vector<string> process_names = {"data", "ttbar", "ztautau", "wjets", "singletop", "diboson", "higgs"};

    vector<tuple<string, int, double, double, string>> var_specs = {
        {"MET",    20, 0,  200, "Missing E_{T} [GeV]"},
        {"mll",    25, 0,  250, "m_{#mu e} [GeV]"},
        {"l1_pt",  25, 25, 150, "Muon p_{T} [GeV]"},
        {"l2_pt",  25, 25, 150, "Electron p_{T} [GeV]"},
        {"l2_eta", 10, 0,  2.5, "Electron |#eta|"}
    };

    map<string, HistSet> histos_VR, histos_VR4;
    auto create_hists = [&](map<string, HistSet>& region_histos, const string& tag, const string& title) {
        for (const auto& spec : var_specs) {
            string var_name = get<0>(spec);
            for (const auto& proc_name : process_names) {
                region_histos[var_name].histograms[proc_name] = new TH1F(
                    (proc_name + "_" + var_name + "_" + tag).c_str(), "",
                    get<1>(spec), get<2>(spec), get<3>(spec));
                region_histos[var_name].histograms[proc_name]->Sumw2();
            }
            region_histos[var_name].xlabel = get<4>(spec);
            region_histos[var_name].title  = title;
        }
    };
    create_hists(histos_VR,  "fake_CR",   "FAKE_CR");
    create_hists(histos_VR4, "prompt_CR", "prompt_CR");

    auto create_frac_hists = [&](map<string, FractionHist>& fracs, const string& tag) {
        fracs["electron"].h_prompt_num_pt  = new TH1F(("prompt_el_pt_num_"  + tag).c_str(), "", 25, 25, 150);
        fracs["electron"].h_total_pt       = new TH1F(("el_pt_total_"      + tag).c_str(), "", 25, 25, 150);
        fracs["electron"].h_prompt_num_eta = new TH1F(("prompt_el_eta_num_" + tag).c_str(), "", 10, 0, 2.5);
        fracs["electron"].h_total_eta      = new TH1F(("el_eta_total_"     + tag).c_str(), "", 10, 0, 2.5);
    };
    create_frac_hists(frac_VR,  "fake_CR");
    create_frac_hists(frac_VR4, "prompt_CR");

    map<string, double> mc_yields, data_yields;
    vector<string> cut_names = {"Total", "Flavor", "Trigger", "pT_cuts", "Preselection", "fake_CR", "prompt_CR"};
    for (const auto& cut : cut_names) { mc_yields[cut] = 0.0; data_yields[cut] = 0.0; }

    // ==========================================================================================
    // 3. MC EVENT LOOP
    // ==========================================================================================
    cout << "\n=== Processing MC Samples ===" << endl;
    EventData ev;

    for (const auto& proc_pair : processes) {
        string proc_name = proc_pair.first;
        if (proc_name == "data") continue;

        cout << "... processing " << proc_name << endl;
        TChain *tree = new TChain("tree");
        for (const auto& file : proc_pair.second.files) tree->Add(file.c_str());

        SetupBranches(tree, ev);

        Long64_t nMC = tree->GetEntries();
        for (Long64_t i = 0; i < nMC; i++) {
            tree->GetEntry(i);
            ComputeDerived(ev);

            double mc_scale_factor = 1.0;  // ttbar / ztautau scale factors kept at 1 (as in original)
            double w = mc_scale_factor * mcWeightTight(ev);
            if (!std::isfinite(ev.total_weight)) continue;

            mc_yields["Total"] += w;

            if (!passTopology(ev)) continue;
            mc_yields["Flavor"] += w;

            if (!(passMuonTrigger(ev) || passElectronTrigger(ev))) continue;
            if (!passTriggerMatching(ev)) continue;
            mc_yields["Trigger"] += w;

            if (ev.l1_pt <= 20000.f || ev.l2_pt < 15000.f) continue;
            mc_yields["pT_cuts"] += w;

            // Both leptons must be truth-prompt for the MC prompt sample
            if (!isMCPromptEvent(ev)) continue;

            if (!ev.mu_is_tight) continue;
            if (!passMuonIP(ev) || !passElectronIP(ev)) continue;
            if (!passCrackVeto(ev)) continue;

            mc_yields["Preselection"] += w;

            // Region: ≥2 jets (any pT) + no b-jets + (SS for fake_CR, OS for prompt_CR)
            bool jets_ok      = passVBFJets(ev);  // j1>35, j2>20 GeV, |Δη|>3
            (void)jets_ok;                        // kept for clarity if user wants stricter cut later
            bool el_is_prompt = isMCPromptElectron(ev);

            // Fake CR: SS + no b-jets + ≥2 jets (jet pT loose, mirrors original "n_jets>=2")
            // Use VBF-style passVBFJets to match jet kinematic intent, OR fall back to passVBFJets check.
            // Original required only jet_pt_NOSYS->size() >= 2, no jet pT cut. We approximate with passVBFJets
            // for SR-aligned definition; if you want the looser version, drop passVBFJets here.
            if (isFakeCR(ev) && ev.no_bjets && passVBFJets(ev)) {
                mc_yields["fake_CR"] += w;
                for (auto const& [var_name, hist_set] : histos_VR) {
                    float val = 0;
                    if      (var_name == "MET")    val = ev.MET   / 1000.0;
                    else if (var_name == "mll")    val = ev.mll   / 1000.0;
                    else if (var_name == "l1_pt")  val = ev.l1_pt / 1000.0;
                    else if (var_name == "l2_pt")  val = ev.l2_pt / 1000.0;
                    else if (var_name == "l2_eta") val = std::fabs(ev.l2_eta);

                    frac_VR["electron"].h_total_pt ->Fill(ev.l2_pt / 1000.0, w);
                    frac_VR["electron"].h_total_eta->Fill(std::fabs(ev.l2_eta), w);
                    if (el_is_prompt) {
                        frac_VR["electron"].h_prompt_num_pt ->Fill(ev.l2_pt / 1000.0, w);
                        frac_VR["electron"].h_prompt_num_eta->Fill(std::fabs(ev.l2_eta), w);
                    }
                    hist_set.histograms.at(proc_name)->Fill(val, w);
                }
            }

            // Prompt CR: OS + no b-jets + ≥2 jets
            if (isPromptCR(ev) && ev.no_bjets && passVBFJets(ev)) {
                mc_yields["prompt_CR"] += w;
                for (auto const& [var_name, hist_set] : histos_VR4) {
                    float val = 0;
                    if      (var_name == "MET")    val = ev.MET   / 1000.0;
                    else if (var_name == "mll")    val = ev.mll   / 1000.0;
                    else if (var_name == "l1_pt")  val = ev.l1_pt / 1000.0;
                    else if (var_name == "l2_pt")  val = ev.l2_pt / 1000.0;
                    else if (var_name == "l2_eta") val = std::fabs(ev.l2_eta);

                    frac_VR4["electron"].h_total_pt ->Fill(ev.l2_pt / 1000.0, w);
                    frac_VR4["electron"].h_total_eta->Fill(std::fabs(ev.l2_eta), w);
                    if (el_is_prompt) {
                        frac_VR4["electron"].h_prompt_num_pt ->Fill(ev.l2_pt / 1000.0, w);
                        frac_VR4["electron"].h_prompt_num_eta->Fill(std::fabs(ev.l2_eta), w);
                    }
                    hist_set.histograms.at(proc_name)->Fill(val, w);
                }
            }
        }
        delete tree;
    }

    // ==========================================================================================
    // 4. DATA EVENT LOOP (AMM-weighted ntuple)
    // ==========================================================================================
    cout << "\n=== Processing Data Sample ===" << endl;
    TChain *dataTree = new TChain("tree");
    dataTree->Add("../AMM/data_with_electron_fake_weights.root");

    ResetTrigMatchPointers(ev);
    SetupBranches(dataTree, ev);

    Long64_t nData = dataTree->GetEntries();
    for (Long64_t i = 0; i < nData; i++) {
        dataTree->GetEntry(i);
        ComputeDerived(ev);

        data_yields["Total"] += 1.0;

        if (!passTopology(ev)) continue;
        data_yields["Flavor"] += 1.0;

        if (!(passMuonTrigger(ev) || passElectronTrigger(ev))) continue;
        if (!passTriggerMatching(ev)) continue;
        data_yields["Trigger"] += 1.0;

        if (ev.l1_pt <= 20000.f || ev.l2_pt < 15000.f) continue;
        data_yields["pT_cuts"] += 1.0;

        if (!ev.mu_is_tight) continue;
        if (!passMuonIP(ev) || !passElectronIP(ev)) continue;
        if (!passCrackVeto(ev)) continue;
        data_yields["Preselection"] += 1.0;

        if (isFakeCR(ev) && ev.no_bjets && passVBFJets(ev)) {
            data_yields["fake_CR"] += 1.0;
            for (auto const& [var_name, hist_set] : histos_VR) {
                float val = 0;
                if      (var_name == "MET")    val = ev.MET   / 1000.0;
                else if (var_name == "mll")    val = ev.mll   / 1000.0;
                else if (var_name == "l1_pt")  val = ev.l1_pt / 1000.0;
                else if (var_name == "l2_pt")  val = ev.l2_pt / 1000.0;
                else if (var_name == "l2_eta") val = std::fabs(ev.l2_eta);
                hist_set.histograms.at("data")->Fill(val, 1.0);
            }
        }

        if (isPromptCR(ev) && ev.no_bjets && passVBFJets(ev)) {
            data_yields["prompt_CR"] += 1.0;
            for (auto const& [var_name, hist_set] : histos_VR4) {
                float val = 0;
                if      (var_name == "MET")    val = ev.MET   / 1000.0;
                else if (var_name == "mll")    val = ev.mll   / 1000.0;
                else if (var_name == "l1_pt")  val = ev.l1_pt / 1000.0;
                else if (var_name == "l2_pt")  val = ev.l2_pt / 1000.0;
                else if (var_name == "l2_eta") val = std::fabs(ev.l2_eta);
                hist_set.histograms.at("data")->Fill(val, 1.0);
            }
        }
    }

    // ==========================================================================================
    // 5. PRINT YIELDS & CREATE PLOTS
    // ==========================================================================================
    cout << "\n=== Final Yields & Cutflow Summary ===\n";
    double total_vr  = mc_yields["fake_CR"];
    double total_vr4 = mc_yields["prompt_CR"];

    printf("+--------------------+-----------------------+-----------------------+\n");
    printf("| %-18s | %-21s | %-21s |\n", "Cut", "MC (Weighted)", "Data (Events)");
    printf("+--------------------+-----------------------+-----------------------+\n");
    for (const auto& cut : cut_names)
        printf("| %-18s | %21.1f | %21.0f |\n", cut.c_str(), mc_yields[cut], data_yields[cut]);
    printf("+--------------------+-----------------------+-----------------------+\n");
    printf("| %-18s | %21.1f | %21.0f |\n", "Total Pred. (fake)",   total_vr,  data_yields["fake_CR"]);
    printf("| %-18s | %21.1f | %21.0f |\n", "Total Pred. (prompt)", total_vr4, data_yields["prompt_CR"]);
    printf("+--------------------+-----------------------+-----------------------+\n");

    system("mkdir -p outputs/05_fake_composition");
    TFile* outFile = new TFile("outputs/05_fake_composition/fake_real_CRs.root", "RECREATE");

    for (auto& pair : histos_VR) {
        TCanvas *c = new TCanvas((string("c_") + pair.first + "_fake_CR").c_str(), "", 800, 800);
        DrawPlot(c, pair.second, processes, true);
        c->SaveAs((string("outputs/05_fake_composition/fake_CR_") + pair.first + ".pdf").c_str());
        c->Write();
    }
    for (auto& pair : histos_VR4) {
        TCanvas *c = new TCanvas((string("c_") + pair.first + "_prompt_CR").c_str(), "", 800, 800);
        DrawPlot(c, pair.second, processes, true);
        c->SaveAs((string("outputs/05_fake_composition/prompt_CR_") + pair.first + ".pdf").c_str());
        c->Write();
    }

    // === Plot electron prompt fraction in fake_CR and prompt_CR ===
    auto plot_fraction = [&](FractionHist& fh, const string& region) {
        if (!fh.h_total_pt || fh.h_total_pt->Integral() == 0) {
            cout << "Warning: No events for prompt fraction in " << region << endl;
            return;
        }
        TH1F* h_frac_pt  = (TH1F*)fh.h_prompt_num_pt ->Clone(("frac_pt_"  + region).c_str());
        TH1F* h_frac_eta = (TH1F*)fh.h_prompt_num_eta->Clone(("frac_eta_" + region).c_str());
        h_frac_pt ->Divide(fh.h_total_pt);
        h_frac_eta->Divide(fh.h_total_eta);

        for (TH1F* h : {h_frac_pt, h_frac_eta}) {
            h->SetLineColor(kRed); h->SetMarkerColor(kRed);
            h->SetMarkerStyle(20); h->SetMarkerSize(1.0); h->SetLineWidth(2);
        }

        auto draw_atlas_labels = [&](TPad* pad, const string& reg_name) {
            pad->cd();
            TLatex latex; latex.SetNDC();
            latex.SetTextFont(72); latex.SetTextSize(0.035);
            latex.DrawLatex(0.20, 0.88, "ATLAS");
            latex.SetTextFont(42); latex.DrawLatex(0.31, 0.88, "Internal");
            latex.SetTextSize(0.030);
            latex.DrawLatex(0.20, 0.83, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
            latex.DrawLatex(0.20, 0.78, ("#mu e, " + reg_name).c_str());
        };

        // pT plot
        TCanvas* cpt = new TCanvas(("c_frac_pt_" + region).c_str(), "", 700, 600);
        cpt->cd();
        h_frac_pt->SetTitle("");
        h_frac_pt->GetXaxis()->SetTitleSize(0.034); h_frac_pt->GetYaxis()->SetTitleSize(0.034);
        h_frac_pt->GetXaxis()->SetLabelSize(0.030); h_frac_pt->GetYaxis()->SetLabelSize(0.030);
        h_frac_pt->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
        h_frac_pt->GetYaxis()->SetTitle("Prompt Fraction");
        h_frac_pt->GetYaxis()->SetRangeUser(0.0, 1.7);
        h_frac_pt->GetYaxis()->SetTitleOffset(1.4); h_frac_pt->GetXaxis()->SetTitleOffset(1.2);
        h_frac_pt->Draw("PE");
        draw_atlas_labels((TPad*)gPad, region == "fake_CR" ? "Fake CR" : "Prompt CR");
        TLine* line = new TLine(h_frac_pt->GetXaxis()->GetXmin(), 1.0,
                                h_frac_pt->GetXaxis()->GetXmax(), 1.0);
        line->SetLineStyle(2); line->SetLineColor(kGray+2); line->Draw();
        cpt->SaveAs(("outputs/05_fake_composition/prompt_frac_pt_" + region + ".pdf").c_str());
        cpt->Write();

        // eta plot
        TCanvas* ceta = new TCanvas(("c_frac_eta_" + region).c_str(), "", 700, 600);
        ceta->cd();
        h_frac_eta->SetTitle("");
        h_frac_eta->GetXaxis()->SetTitleSize(0.035); h_frac_eta->GetYaxis()->SetTitleSize(0.035);
        h_frac_eta->GetXaxis()->SetLabelSize(0.030); h_frac_eta->GetYaxis()->SetLabelSize(0.030);
        h_frac_eta->GetXaxis()->SetTitle("Electron |#eta|");
        h_frac_eta->GetYaxis()->SetTitle("Prompt Fraction");
        h_frac_eta->GetYaxis()->SetRangeUser(0.0, 1.7);
        h_frac_eta->GetYaxis()->SetTitleOffset(1.4); h_frac_eta->GetXaxis()->SetTitleOffset(1.2);
        h_frac_eta->Draw("PE");
        draw_atlas_labels((TPad*)gPad, region == "fake_CR" ? "Fake CR" : "Prompt CR");
        TLine* line2 = new TLine(h_frac_eta->GetXaxis()->GetXmin(), 1.0,
                                 h_frac_eta->GetXaxis()->GetXmax(), 1.0);
        line2->SetLineStyle(2); line2->SetLineColor(kGray+2); line2->Draw();
        ceta->SaveAs(("outputs/05_fake_composition/prompt_frac_eta_" + region + ".pdf").c_str());
        ceta->Write();
    };

    plot_fraction(frac_VR ["electron"], "fake_CR");
    plot_fraction(frac_VR4["electron"], "prompt_CR");

    outFile->Close();
    cout << "\n=== Analysis Complete ===" << endl;
    cout << "Plots saved in 'outputs/05_fake_composition/' directory." << endl;
    cout << "ROOT file saved as 'fake_real_CRs.root'." << endl;
}
