// =============================================================================
// plotting/presel_VR_plots.C
// PURPOSE: Stacked MC + data distributions for the OS preselection (VR-loose)
//          and the kinematic VR (mjj < 400, no b-jets, VBF jets).
// PATTERN: Centralized EventData + SelectionUtils.h.
// INPUTS:  MC ntuples on EOS, ../AMM/data_with_electron_fake_weights.root
// OUTPUTS: outputs/10_presel_vr/
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
#include <TPad.h>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <map>
#include <tuple>
#include "../utils/SelectionUtils.h"

using namespace std;

// ---------------------------------------------------------------------------
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

// Stacked plot with statistical error band and ratio panel
void DrawPlot(TCanvas* c, HistSet& hs, map<string, Process>& processes, bool logY = false) {
    c->cd();
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.02);
    pad1->SetTopMargin(0.07);
    pad1->Draw(); pad1->cd();
    if (logY) pad1->SetLogy();

    THStack *stack = new THStack("stack", "");
    TH1F* h_bkg_total = nullptr;
    vector<string> stack_order = {"ttbar","ztautau","wjets","singletop","diboson","higgs"};
    for (const auto& proc_name : stack_order) {
        if (hs.histograms.count(proc_name)) {
            TH1F* h = hs.histograms[proc_name];
            h->SetFillColor(processes[proc_name].color);
            h->SetLineColor(kBlack); h->SetLineWidth(1);
            stack->Add(h);
            if (!h_bkg_total) h_bkg_total = (TH1F*)h->Clone((string(h->GetName())+"_total").c_str());
            else              h_bkg_total->Add(h);
        }
    }
    stack->Draw("HIST");

    TH1F* h_bkg_total_stat = nullptr;
    if (h_bkg_total) {
        h_bkg_total_stat = (TH1F*)h_bkg_total->Clone("h_bkg_total_stat");
        h_bkg_total_stat->SetFillColor(kGray+1);
        h_bkg_total_stat->SetFillStyle(3244);
        h_bkg_total_stat->SetMarkerSize(0);
        h_bkg_total_stat->SetLineWidth(0);
        h_bkg_total_stat->Draw("E2 SAME");
    }

    TH1F* h_data = hs.histograms["data"];
    h_data->SetLineColor(kBlack); h_data->SetMarkerColor(kBlack);
    h_data->SetMarkerStyle(20);   h_data->SetMarkerSize(1.0);
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
    for (const auto& proc_name : stack_order)
        if (hs.histograms.count(proc_name))
            leg->AddEntry(hs.histograms[proc_name], processes[proc_name].legend.c_str(), "F");
    leg->SetTextSize(0.04); leg->Draw();

    TLatex latex; latex.SetNDC(); latex.SetTextFont(72); latex.SetTextSize(0.06);
    latex.DrawLatex(0.2, 0.85, "ATLAS"); latex.SetTextFont(42);
    latex.DrawLatex(0.34, 0.85, "Internal");
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.2, 0.78, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
    latex.DrawLatex(0.2, 0.72, hs.title.c_str());

    c->cd();
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0.05); pad2->SetBottomMargin(0.35);
    pad2->Draw(); pad2->cd();

    TH1F *h_ratio = (TH1F*)h_data->Clone("h_ratio");
    if (h_bkg_total && h_bkg_total->Integral() > 0) {
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
            } else { h_ratio->SetBinContent(bin, 0); h_ratio->SetBinError(bin, 0); }
        }
    }
    h_ratio->GetYaxis()->SetRangeUser(0.4, 2.5);
    h_ratio->GetYaxis()->SetTitle("Data / Bkg");
    h_ratio->GetYaxis()->SetNdivisions(505);
    h_ratio->GetYaxis()->SetTitleSize(0.12);
    h_ratio->GetYaxis()->SetLabelSize(0.1);
    h_ratio->GetYaxis()->SetTitleOffset(0.5);
    h_ratio->GetXaxis()->SetTitle(hs.xlabel.c_str());
    h_ratio->GetXaxis()->SetTitleSize(0.12);
    h_ratio->GetXaxis()->SetLabelSize(0.1);
    h_ratio->GetXaxis()->SetTitleOffset(1.2);

    TH1F* h_ratio_band = nullptr;
    if (h_bkg_total) {
        h_ratio_band = (TH1F*)h_ratio->Clone("h_ratio_band");
        for (int bin = 1; bin <= h_ratio_band->GetNbinsX(); ++bin) {
            double bkg_val = h_bkg_total->GetBinContent(bin);
            double bkg_err = h_bkg_total->GetBinError(bin);
            if (bkg_val > 0) h_ratio_band->SetBinError(bin, bkg_err/bkg_val);
        }
        h_ratio_band->SetFillColor(kGray+1);
        h_ratio_band->SetFillStyle(3244);
        h_ratio_band->Draw("E2 SAME");
    }
    h_ratio->Draw("PE SAME");
    TLine *line = new TLine(h_ratio->GetXaxis()->GetXmin(), 1.0,
                            h_ratio->GetXaxis()->GetXmax(), 1.0);
    line->SetLineColor(kRed); line->SetLineStyle(2); line->Draw();
    c->Update();
}

// Helper: convert var-name to value (in GeV)
inline float varValue(const std::string& var, const EventData& ev) {
    if (var == "MET")    return ev.MET    / 1000.0f;
    if (var == "mll")    return ev.mll    / 1000.0f;
    if (var == "l1_pt")  return ev.l1_pt  / 1000.0f;
    if (var == "l2_pt")  return ev.l2_pt  / 1000.0f;
    if (var == "l2_eta") return std::fabs(ev.l2_eta);
    return 0.f;
}

void presel_VR_plots() {
    SetATLASStyle();
    cout << "\n=== presel_VR_plots ===" << endl;

    map<string, Process> processes;
    string base = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    processes["data"]      = {"data","Data",{},kBlack,true};
    processes["ttbar"]     = {"ttbar","t#bar{t}",         {base+"ttbar/*.root"},     TColor::GetColor("#8E44AD")};
    processes["ztautau"]   = {"ztautau","Z#rightarrow#tau#tau",{base+"ztautau/*.root"},   TColor::GetColor("#F39C12")};
    processes["wjets"]     = {"wjets","W+jets",           {base+"wjets/*.root"},     TColor::GetColor("#2ECC71")};
    processes["singletop"] = {"singletop","Single Top",   {base+"singletop/*.root"}, TColor::GetColor("#E74C3C")};
    processes["diboson"]   = {"diboson","Diboson",        {base+"diboson/*.root"},   TColor::GetColor("#3498DB")};
    processes["higgs"]     = {"higgs","Higgs",            {base+"higgs/*.root"},     TColor::GetColor("#FF2D95")};

    vector<string> proc_names = {"data","ttbar","ztautau","wjets","singletop","diboson","higgs"};

    vector<tuple<string,int,double,double,string>> var_specs = {
        {"MET",   20, 0,   200, "Missing E_{T} [GeV]"},
        {"mll",   25, 0,   250, "m_{#mu e} [GeV]"},
        {"l1_pt", 25, 25,  150, "Muon p_{T} [GeV]"},
        {"l2_pt", 25, 25,  150, "Electron p_{T} [GeV]"},
        {"l2_eta",10, 0,   2.5, "Electron |#eta|"}
    };

    map<string,HistSet> histos_VR, histos_VR4;
    auto create_hists = [&](map<string,HistSet>& region_histos, const string& tag, const string& title) {
        for (const auto& spec : var_specs) {
            string var = get<0>(spec);
            for (const auto& p : proc_names) {
                region_histos[var].histograms[p] = new TH1F(
                    (p + "_" + var + "_" + tag).c_str(), "",
                    get<1>(spec), get<2>(spec), get<3>(spec));
                region_histos[var].histograms[p]->Sumw2();
            }
            region_histos[var].xlabel = get<4>(spec);
            region_histos[var].title  = title;
        }
    };
    create_hists(histos_VR,  "VR", "Preselections");
    create_hists(histos_VR4, "VR", "VR:  (m_{jj} < 400 GeV)");

    map<string,double> mc_yields, data_yields;
    vector<string> cut_names = {"Total","Flavor","Trigger","pT_cuts","Preselection","OS_VR","VR4"};
    for (const auto& c : cut_names) { mc_yields[c]=0.; data_yields[c]=0.; }

    // ----- MC loop -----
    cout << "\n=== Processing MC Samples ===" << endl;
    for (const auto& proc_pair : processes) {
        const string& pname = proc_pair.first;
        if (pname == "data") continue;
        cout << "... " << pname << endl;
        TChain* tree = new TChain("tree");
        for (const auto& f : proc_pair.second.files) tree->Add(f.c_str());

        EventData ev; ResetTrigMatchPointers(ev);
        SetupBranches(tree, ev);

        Long64_t n = tree->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            tree->GetEntry(i); ComputeDerived(ev);
            double w = mcWeightTight(ev);
            if (!std::isfinite(w)) continue;
            mc_yields["Total"] += w;

            if (!passTopology(ev))                                          continue;
            mc_yields["Flavor"] += w;
            if (!(passTrigger(ev) && passTriggerMatching(ev)))              continue;
            mc_yields["Trigger"] += w;
            if (!passLeptonPt(ev))                                          continue;
            mc_yields["pT_cuts"] += w;

            if (!isMCPromptEvent(ev))                                       continue;
            if (!ev.mu_is_tight)                                            continue;
            if (!passMuonIP(ev) || !passElectronIP(ev))                     continue;
            if (!passCrackVeto(ev))                                         continue;
            mc_yields["Preselection"] += w;

            // OS preselection (presel) region
            if (ev.is_os) {
                mc_yields["OS_VR"] += w;
                for (const auto& [var, hs] : histos_VR)
                    hs.histograms.at(pname)->Fill(varValue(var, ev), w);
            }
            // VR: OS + VBF jets + mjj<400 GeV + no b-jets
            bool jet_cuts  = passVBFJets(ev);
            bool mjj_inv   = (ev.mjj > 0. && ev.mjj < 400000.);
            if (ev.is_os && jet_cuts && mjj_inv && ev.no_bjets) {
                mc_yields["VR4"] += w;
                for (const auto& [var, hs] : histos_VR4)
                    hs.histograms.at(pname)->Fill(varValue(var, ev), w);
            }
        }
        delete tree;
    }

    // ----- Data loop -----
    cout << "\n=== Processing Data Sample ===" << endl;
    TChain* dataTree = new TChain("tree");
    dataTree->Add("../AMM/data_with_electron_fake_weights.root");
    EventData ev; ResetTrigMatchPointers(ev);
    SetupBranches(dataTree, ev);

    Long64_t nD = dataTree->GetEntries();
    for (Long64_t i = 0; i < nD; ++i) {
        dataTree->GetEntry(i); ComputeDerived(ev);
        data_yields["Total"] += 1.0;
        if (!passTopology(ev))                              continue;
        data_yields["Flavor"] += 1.0;
        if (!(passTrigger(ev) && passTriggerMatching(ev)))  continue;
        data_yields["Trigger"] += 1.0;
        if (!passLeptonPt(ev))                              continue;
        data_yields["pT_cuts"] += 1.0;
        if (!ev.mu_is_tight)                                continue;
        if (!passMuonIP(ev) || !passElectronIP(ev))         continue;
        if (!passCrackVeto(ev))                             continue;
        data_yields["Preselection"] += 1.0;

        if (ev.is_os) {
            data_yields["OS_VR"] += 1.0;
            for (const auto& [var, hs] : histos_VR)
                hs.histograms.at("data")->Fill(varValue(var, ev), 1.0);
        }
        bool jet_cuts = passVBFJets(ev);
        bool mjj_inv  = (ev.mjj > 0. && ev.mjj < 400000.);
        if (ev.is_os && jet_cuts && mjj_inv && ev.no_bjets) {
            data_yields["VR4"] += 1.0;
            for (const auto& [var, hs] : histos_VR4)
                hs.histograms.at("data")->Fill(varValue(var, ev), 1.0);
        }
    }

    // ----- Yields table -----
    cout << "\n=== Final Yields & Cutflow Summary ===\n";
    printf("+--------------------+-----------------------+-----------------------+\n");
    printf("| %-18s | %-21s | %-21s |\n", "Cut", "MC (Weighted)", "Data (Events)");
    printf("+--------------------+-----------------------+-----------------------+\n");
    for (const auto& cut : cut_names)
        printf("| %-18s | %21.1f | %21.0f |\n",
               cut.c_str(), mc_yields[cut], data_yields[cut]);
    printf("+--------------------+-----------------------+-----------------------+\n");
    printf("| %-18s | %21.1f | %21.0f |\n", "Total Pred. (VR)",  mc_yields["OS_VR"], data_yields["OS_VR"]);
    printf("| %-18s | %21.1f | %21.0f |\n", "Total Pred. (VR4)", mc_yields["VR4"],   data_yields["VR4"]);
    printf("+--------------------+-----------------------+-----------------------+\n");

    // ----- Plots -----
    system("mkdir -p outputs/10_presel_vr");
    TFile* outFile = new TFile("outputs/10_presel_vr/presel_VR_plots.root", "RECREATE");
    for (auto& pair : histos_VR) {
        TCanvas *c = new TCanvas(("c_" + pair.first + "_VR").c_str(), "", 800, 800);
        DrawPlot(c, pair.second, processes, true);
        c->SaveAs(("outputs/10_presel_vr/presel_" + pair.first + ".pdf").c_str());
        c->Write();
    }
    for (auto& pair : histos_VR4) {
        TCanvas *c = new TCanvas(("c_" + pair.first + "_VR4").c_str(), "", 800, 800);
        DrawPlot(c, pair.second, processes, true);
        c->SaveAs(("outputs/10_presel_vr/VR_" + pair.first + ".pdf").c_str());
        c->Write();
    }
    outFile->Close();
    cout << "\n=== Analysis Complete ===" << endl;
    cout << "Plots saved in 'outputs/10_presel_vr/' directory." << endl;
}
