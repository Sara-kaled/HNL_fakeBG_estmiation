// =============================================================================
// FILE:    studies/bkg_composition.C
// PURPOSE: Show the data-driven fake-electron contribution as a fraction of
//          the total background in each analysis region:
//              fake-CR (SS μe)     — sanity (should be ≈100 % fake)
//              HNL_CR_tt           — ttbar-rich, OS, ≥1 b-jet, mjj>400
//              HNL_CR_ztautau      — Zττ-rich, OS, b-veto, 35<l1_pt<45
//              HNL_VR              — OS, b-veto, mjj<400
//              HNL_SR              — OS, b-veto, mjj>400, VBF
//
// HOW IT'S COMPUTED:
//   - Fake yield in a region = Σ_data fake_weight_nominal · 1[region cuts]
//     using the AMM-weighted ntuple from AMM/data_with_electron_fake_weights.root.
//   - Per-process prompt yield = run over each MC sample in
//     /eos/user/.../fresh_ntuples/<process>/*.root, apply the SAME region cuts,
//     filter to TRUTH-PROMPT events only (isMCPromptEvent) so we don't
//     double-count fakes that the data-driven method already covers, and sum
//     the standard MC weight (mcWeightTight).
//   - Total BG per region = Σ_processes prompt + fake.
//   - Output:
//       studies/outputs/bkg_composition.txt   — yield table
//       studies/outputs/bkg_composition.pdf   — bar chart per region
//
// RUNS AFTER: AMM_ElectronOnly.C  (needs both data + MC AMM ntuples to exist)
// =============================================================================
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include "ROOT/RVec.hxx"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <string>
#include <vector>
#include "../utils/SelectionUtils.h"
#include "../utils/AtlasStyle.h"

using namespace ROOT;
using namespace ROOT::VecOps;

// ---------------------------------------------------------------------------
// Region table — one row each, in the printout / plot order
// ---------------------------------------------------------------------------
enum Reg { kFakeCR=0, kCRtt=1, kCRztt=2, kVR=3, kSR=4, kNReg=5 };
static const char* kRegName[kNReg] = {
    "fake-CR (SS)", "HNL_CR_tt", "HNL_CR_ztautau", "HNL_VR", "HNL_SR"
};

// ---------------------------------------------------------------------------
// Decide which region(s) the event belongs to (events can fall into ≥1)
// ---------------------------------------------------------------------------
static void regionFlags(const EventData& ev, bool flags[kNReg])
{
    for (int r = 0; r < kNReg; ++r) flags[r] = false;
    if (!passBasePresel(ev)) return;
    flags[kFakeCR] = isFakeCR(ev) && ev.no_bjets
                    && (ev.jet1_pt > 0. && ev.jet2_pt > 0.);
    flags[kCRtt]   = isTtbarCR(ev);
    flags[kCRztt]  = isZtautauCR(ev);
    flags[kVR]     = isVR(ev);
    flags[kSR]     = isSR(ev);
}

// ---------------------------------------------------------------------------
// 1) Data fake yield per region from the AMM-weighted ntuple
// ---------------------------------------------------------------------------
static void sumFake(double yields[kNReg])
{
    for (int r = 0; r < kNReg; ++r) yields[r] = 0.;

    auto* fAMM = TFile::Open("../AMM/data_with_electron_fake_weights.root");
    if (!fAMM || fAMM->IsZombie()) {
        std::cerr << "ERROR: AMM/data_with_electron_fake_weights.root not found.\n"
                     "       Run AMM_ElectronOnly.C first.\n";
        return;
    }
    auto* tree = (TTree*)fAMM->Get("tree");
    if (!tree) { std::cerr << "ERROR: TTree 'tree' missing.\n"; return; }

    EventData ev;
    SetupBranches(tree, ev);
    double fake_w = 0.;
    tree->SetBranchAddress("fake_weight_nominal", &fake_w);

    Long64_t n = tree->GetEntries();
    std::cout << "  data ntuple: " << n << " events\n";
    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);
        ComputeDerived(ev);
        bool flags[kNReg];
        regionFlags(ev, flags);
        for (int r = 0; r < kNReg; ++r)
            if (flags[r]) yields[r] += fake_w;
    }
    fAMM->Close();
}

// ---------------------------------------------------------------------------
// 2) Per-process MC prompt yield (truth-prompt only — fakes are covered
//    by the data-driven side, not by MC, so we filter them out here).
// ---------------------------------------------------------------------------
static void sumMCProcess(const std::string& dir, double yields[kNReg])
{
    for (int r = 0; r < kNReg; ++r) yields[r] = 0.;
    const std::string base =
        "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    auto* tree = new TChain("tree");
    tree->Add((base + dir + "/*.root").c_str());

    EventData ev;
    SetupBranches(tree, ev);

    Long64_t n = tree->GetEntries();
    if (n == 0) {
        std::cerr << "  WARNING: " << dir << " is empty / not found.\n";
        delete tree;
        return;
    }
    std::cout << "  " << std::setw(12) << dir << ": " << n << " events\n";

    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);
        ComputeDerived(ev);
        if (!isMCPromptEvent(ev)) continue;       // FNP covered by data side
        bool flags[kNReg];
        regionFlags(ev, flags);
        if (!flags[kFakeCR] && !flags[kCRtt] && !flags[kCRztt]
            && !flags[kVR]   && !flags[kSR]) continue;
        const double w = mcWeightTight(ev);
        if (!std::isfinite(w) || w == 0.) continue;
        for (int r = 0; r < kNReg; ++r)
            if (flags[r]) yields[r] += w;
    }
    delete tree;
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
void bkg_composition()
{
    gStyle->SetOptStat(0);
    gSystem->mkdir("outputs", true);

    // -- per-process bookkeeping
    const std::vector<std::string> procs =
        {"ttbar", "wjets", "ztautau", "singletop", "diboson", "higgs"};
    std::map<std::string, std::vector<double>> y_proc;
    for (const auto& p : procs) y_proc[p].assign(kNReg, 0.);

    // -- 1) data fake yield
    std::cout << "\n[1/2] Reading AMM-weighted data ntuple…\n";
    double y_fake[kNReg];
    sumFake(y_fake);

    // -- 2) per-process MC prompt yields
    std::cout << "\n[2/2] Reading per-process MC ntuples…\n";
    for (const auto& p : procs) {
        std::cout << "  process: " << p << "\n";
        sumMCProcess(p, y_proc[p].data());
    }

    // ---------------------------------------------------------------------
    // Print + write yield table
    // ---------------------------------------------------------------------
    std::ofstream rep("outputs/bkg_composition.txt");
    auto out = [&](const std::string& s) { std::cout << s; rep << s; };

    char hdr[1024];
    snprintf(hdr, sizeof(hdr),
        "\n%-18s  %12s  %12s  %12s  %12s  %12s  %12s  %12s  %12s   %s\n",
        "region", "fake", "ttbar", "wjets", "ztautau", "singletop",
        "diboson", "higgs", "TOTAL", "fake/total");
    out(hdr);
    out(std::string(150, '-') + "\n");

    std::vector<double> frac(kNReg, 0.);
    for (int r = 0; r < kNReg; ++r) {
        double tot = y_fake[r];
        for (const auto& p : procs) tot += y_proc[p][r];
        frac[r] = (tot > 0.) ? y_fake[r] / tot : 0.;

        char line[1024];
        snprintf(line, sizeof(line),
            "%-18s  %12.2f  %12.2f  %12.2f  %12.2f  %12.2f  "
            "%12.2f  %12.2f  %12.2f   %5.1f %%\n",
            kRegName[r], y_fake[r],
            y_proc["ttbar"][r], y_proc["wjets"][r], y_proc["ztautau"][r],
            y_proc["singletop"][r], y_proc["diboson"][r], y_proc["higgs"][r],
            tot, 100. * frac[r]);
        out(line);
    }
    rep.close();
    std::cout << "\n✓ Wrote outputs/bkg_composition.txt\n";

    // ---------------------------------------------------------------------
    // Stacked-yield bar chart per region (fake first, then ttbar, …)
    // ---------------------------------------------------------------------
    auto* h_fake   = new TH1D("h_fake",   "", kNReg, 0, kNReg);
    auto* h_ttbar  = new TH1D("h_ttbar",  "", kNReg, 0, kNReg);
    auto* h_wjets  = new TH1D("h_wjets",  "", kNReg, 0, kNReg);
    auto* h_ztt    = new TH1D("h_ztt",    "", kNReg, 0, kNReg);
    auto* h_stop   = new TH1D("h_stop",   "", kNReg, 0, kNReg);
    auto* h_dibo   = new TH1D("h_dibo",   "", kNReg, 0, kNReg);
    auto* h_higgs  = new TH1D("h_higgs",  "", kNReg, 0, kNReg);

    // Brighter ATLAS-palette colors (defined in utils/AtlasStyle.h)
    auto styled = [](TH1D* h, int col){ h->SetFillColor(col); h->SetLineColor(kBlack); };
    styled(h_fake,  kATLAS_fake);
    styled(h_ttbar, kATLAS_ttbar);
    styled(h_wjets, kATLAS_Wjets);
    styled(h_ztt,   kATLAS_Ztautau);
    styled(h_stop,  kATLAS_singletop);
    styled(h_dibo,  kATLAS_diboson);
    styled(h_higgs, kATLAS_Higgs);

    for (int r = 0; r < kNReg; ++r) {
        h_fake ->SetBinContent(r+1, y_fake[r]);
        h_ttbar->SetBinContent(r+1, y_proc["ttbar"][r]);
        h_wjets->SetBinContent(r+1, y_proc["wjets"][r]);
        h_ztt  ->SetBinContent(r+1, y_proc["ztautau"][r]);
        h_stop ->SetBinContent(r+1, y_proc["singletop"][r]);
        h_dibo ->SetBinContent(r+1, y_proc["diboson"][r]);
        h_higgs->SetBinContent(r+1, y_proc["higgs"][r]);
        h_fake ->GetXaxis()->SetBinLabel(r+1, kRegName[r]);
    }

    auto* hs = new THStack("bkg_comp", ";region;weighted yield");
    hs->Add(h_fake);   // bottom: fake (the highlight)
    hs->Add(h_ttbar);
    hs->Add(h_wjets);
    hs->Add(h_ztt);
    hs->Add(h_stop);
    hs->Add(h_dibo);
    hs->Add(h_higgs);

    SetATLASStyle();
    auto* c = new TCanvas("c_bkg_comp", "", 1100, 650);
    c->SetLeftMargin(0.10); c->SetRightMargin(0.22); c->SetBottomMargin(0.14);
    c->SetTopMargin(0.10);
    hs->Draw("hist");
    hs->GetXaxis()->SetLabelSize(0.040);
    hs->GetYaxis()->SetTitleOffset(1.0);

    // Force y-axis headroom so the percentage labels above each bar
    // never get clipped or overlap with the bars themselves.
    double tot_max = 0.;
    for (int r = 0; r < kNReg; ++r) {
        double tot = y_fake[r];
        for (const auto& p : procs) tot += y_proc[p][r];
        if (tot > tot_max) tot_max = tot;
    }
    hs->SetMaximum(1.45 * tot_max);

    auto* leg = new TLegend(0.79, 0.40, 0.985, 0.85);
    leg->SetBorderSize(0); leg->SetFillStyle(0);
    leg->SetTextSize(0.038);
    leg->AddEntry(h_higgs, "Higgs",         "f");
    leg->AddEntry(h_dibo,  "Diboson",       "f");
    leg->AddEntry(h_stop,  "Single top",    "f");
    leg->AddEntry(h_ztt,   "Z#rightarrow#tau#tau","f");
    leg->AddEntry(h_wjets, "W+jets",        "f");
    leg->AddEntry(h_ttbar, "t#bar{t}",      "f");
    leg->AddEntry(h_fake,  "Fake (data-driven)", "f");
    leg->Draw();

    // Annotate fake fraction above each bar — placed at 1.10× the bar
    // top, with the y-axis already extended above (SetMaximum) to leave
    // room. Bigger text and bold so it's readable.
    TLatex t; t.SetTextSize(0.040); t.SetTextAlign(22); t.SetTextFont(62);
    for (int r = 0; r < kNReg; ++r) {
        double tot = y_fake[r];
        for (const auto& p : procs) tot += y_proc[p][r];
        if (tot <= 0.) continue;
        t.DrawLatex(r + 0.5, tot * 1.08,
            TString::Format("fake %.0f %%", 100. * frac[r]));
    }

    DrawAtlasLabel(0.12, 0.94, "Internal");

    TLatex hd; hd.SetNDC(); hd.SetTextFont(42); hd.SetTextSize(0.030);
    hd.DrawLatex(0.12, 0.84,
        "Fake / total background per region (post-AMM)");

    c->SaveAs("outputs/bkg_composition.pdf");
    std::cout << "✓ Wrote outputs/bkg_composition.pdf\n";
}
