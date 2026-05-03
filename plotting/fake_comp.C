// =============================================================================
// FILE:    plotting/fake_comp.C
// PURPOSE: Plot the truth composition of fake electrons in CR (SS preselection)
//          and SR (OS, no b-jets, VBF jet cuts) using IFFClass / truth_origin.
// INPUTS:  MC ntuples (ttbar, W+jets, Z→ττ, single-top, diboson, Higgs).
// OUTPUTS: outputs/05_fake_composition/fake_electron_<var>_<region>.{pdf,png}
// =============================================================================
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLatex.h>
#include <TLine.h>
#include <TSystem.h>

#include "ROOT/RVec.hxx"
#include "../utils/SelectionUtils.h"

#include <iostream>
#include <cmath>

using namespace ROOT;
using namespace ROOT::VecOps;

enum FakeSource {
  BottomHadron,
  CharmHadron,
  ConvertedPhoton,
  LightHadron,
  TauDecay,
  OtherSources,
  NumSources
};

enum Region {
  CR,
  SR,
  NumRegions
};

struct SourceInfo {
  int color;
  int marker;
  const char* label;
};

void fake_comp()
{
  // ------------------------------------------------------------------
  // 1. Input TChain
  // ------------------------------------------------------------------
  TChain* tree = new TChain("tree");
  tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/ztautau/*.root");
  tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/wjets/*.root");
  tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/ttbar/*.root");
  tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/singletop/*.root");
  tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/diboson/*.root");
  tree->Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/higgs/*.root");

  // ------------------------------------------------------------------
  // 2. Central EventData + shared branch setup
  // ------------------------------------------------------------------
  EventData ev;
  SetupBranches(tree, ev);

  // ---- Extra branches not in EventData ----
  RVecF* el_pt_NOSYS = nullptr;     // per-electron pT (MeV)
  RVecF* el_d0_NOSYS = nullptr;     // raw d0 value (mm)
  RVecI* el_IFFClass = nullptr;     // truth-classifier (electron)
  RVecI* mu_IFFClass = nullptr;     // truth-classifier (muon) — used for prompt muon ID
  RVecI* el_truth_origin     = nullptr;
  RVecI* el_truth_fromHadron = nullptr;
  RVecI* el_truth_fromTau    = nullptr;
  RVecI* el_truth_fromBSM    = nullptr;
  tree->SetBranchAddress("el_pt_NOSYS",         &el_pt_NOSYS);
  tree->SetBranchAddress("el_d0_NOSYS",         &el_d0_NOSYS);
  tree->SetBranchAddress("el_IFFClass",         &el_IFFClass);
  tree->SetBranchAddress("mu_IFFClass",         &mu_IFFClass);
  tree->SetBranchAddress("el_truth_origin",     &el_truth_origin);
  tree->SetBranchAddress("el_truth_fromHadron", &el_truth_fromHadron);
  tree->SetBranchAddress("el_truth_fromTau",    &el_truth_fromTau);
  tree->SetBranchAddress("el_truth_fromBSM",    &el_truth_fromBSM);

  // ------------------------------------------------------------------
  // 3. Histograms: fake composition vs variable, per region
  // ------------------------------------------------------------------
  const int NumVars = 5;
  const char* varNames[NumVars] = {"pt", "eta", "d0", "d0sig", "met"};

  double ptBins[]    = {10, 15, 20, 25, 30, 40, 50, 60, 70};
  int nPtBins        = sizeof(ptBins)/sizeof(double) - 1;
  double etaBins[]   = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5};
  int nEtaBins       = sizeof(etaBins)/sizeof(double) - 1;
  double d0Bins[]    = {0.0, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0};
  int nD0Bins        = sizeof(d0Bins)/sizeof(double) - 1;
  double d0sigBins[] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
  int nD0sigBins     = sizeof(d0sigBins)/sizeof(double) - 1;
  double metBins[]   = {0, 20, 40, 60, 80, 100, 150};
  int nMetBins       = sizeof(metBins)/sizeof(double) - 1;

  TH1D* hTotalFake[NumRegions][NumVars];
  TH1D* hFakeBySource[NumRegions][NumVars][NumSources];

  for (int r = 0; r < NumRegions; ++r) {
    hTotalFake[r][0] = new TH1D(Form("hTotal_pt_r%d", r),    "", nPtBins,    ptBins);
    hTotalFake[r][1] = new TH1D(Form("hTotal_eta_r%d", r),   "", nEtaBins,   etaBins);
    hTotalFake[r][2] = new TH1D(Form("hTotal_d0_r%d", r),    "", nD0Bins,    d0Bins);
    hTotalFake[r][3] = new TH1D(Form("hTotal_d0sig_r%d", r), "", nD0sigBins, d0sigBins);
    hTotalFake[r][4] = new TH1D(Form("hTotal_met_r%d", r),   "", nMetBins,   metBins);

    for (int v = 0; v < NumVars; ++v) {
      for (int s = 0; s < NumSources; ++s) {
        TString name = TString::Format("h_%s_r%d_s%d", varNames[v], r, s);
        if      (v == 0) hFakeBySource[r][v][s] = new TH1D(name, "", nPtBins,    ptBins);
        else if (v == 1) hFakeBySource[r][v][s] = new TH1D(name, "", nEtaBins,   etaBins);
        else if (v == 2) hFakeBySource[r][v][s] = new TH1D(name, "", nD0Bins,    d0Bins);
        else if (v == 3) hFakeBySource[r][v][s] = new TH1D(name, "", nD0sigBins, d0sigBins);
        else             hFakeBySource[r][v][s] = new TH1D(name, "", nMetBins,   metBins);
      }
    }
  }

  // ------------------------------------------------------------------
  // 4. Counters
  // ------------------------------------------------------------------
  double count_preselection = 0.0;
  double count_CR           = 0.0;
  double count_SR           = 0.0;

  Long64_t nentries = tree->GetEntries();
  std::cout << "Total number of entries in the tree: " << nentries << std::endl;

  // ------------------------------------------------------------------
  // 5. Event loop
  // ------------------------------------------------------------------
  for (Long64_t i = 0; i < nentries; ++i) {
    tree->GetEntry(i);
    ComputeDerived(ev);

    // Base preselection (topology + trigger + matching + mu tight + IPs + crack veto + lepton pT)
    if (!passBasePresel(ev)) continue;

    // Additional preselection: muon must be truth-prompt with IFFClass==4 (mirrors original cut)
    bool muon_tight_prompt = isMCPromptMuon(ev)
                          && mu_IFFClass && !mu_IFFClass->empty()
                          && (*mu_IFFClass)[0] == 4;
    if (!muon_tight_prompt) continue;
    count_preselection += 1.0;

    // SR: OS + no b-jets + VBF jet kinematics (no mll/mjj cuts here — original kept it loose)
    bool pass_SR = ev.is_os && ev.no_bjets && passVBFJets(ev);

    // CR: anything passing preselection that is NOT in SR (used for composition study)
    bool in_SR = pass_SR;
    bool in_CR = !pass_SR;
    if (!in_SR && !in_CR) continue;
    if (in_SR) count_SR += 1.0;
    if (in_CR) count_CR += 1.0;

    // Fake-electron requirement: NOT (truth-prompt && IFFClass==2)
    bool el_is_prompt = isMCPromptElectron(ev)
                     && el_IFFClass && !el_IFFClass->empty()
                     && (*el_IFFClass)[0] == 2;
    if (el_is_prompt) continue;

    // ---- Categorize fake source ----
    FakeSource source = OtherSources;
    int  truth_origin = (el_truth_origin     && !el_truth_origin->empty())     ? (*el_truth_origin)[0]      : 0;
    bool from_hadron  = (el_truth_fromHadron && !el_truth_fromHadron->empty()) ? (*el_truth_fromHadron)[0]  : false;
    bool from_tau     = (el_truth_fromTau    && !el_truth_fromTau->empty())    ? (*el_truth_fromTau)[0]     : false;
    bool from_bsm     = (el_truth_fromBSM    && !el_truth_fromBSM->empty())    ? (*el_truth_fromBSM)[0]     : false;

    if      (from_bsm)              source = OtherSources;
    else if (from_tau)              source = TauDecay;
    else if (truth_origin == 5)     source = ConvertedPhoton;
    else if (truth_origin == 45)    source = OtherSources;
    else if (from_hadron) {
      if      (truth_origin == 26)                              source = BottomHadron;
      else if (truth_origin == 25)                              source = CharmHadron;
      else if (truth_origin == 34 || truth_origin == 35)        source = LightHadron;
      else                                                       source = OtherSources;
    } else                          source = OtherSources;

    // ---- Region & weights ----
    int region = in_CR ? CR : SR;
    double fake_weight = mcWeightLoose(ev);

    double el_pt_val    = (el_pt_NOSYS && !el_pt_NOSYS->empty()) ? (*el_pt_NOSYS)[0] / 1000.0 : 0.0;
    double el_eta_val   = (ev.el_eta_vec && !ev.el_eta_vec->empty()) ? std::fabs((*ev.el_eta_vec)[0]) : 0.0;
    double el_d0_val    = (el_d0_NOSYS && !el_d0_NOSYS->empty()) ? std::fabs((*el_d0_NOSYS)[0]) : 0.0;
    double el_d0sig_val = (ev.el_d0sig && !ev.el_d0sig->empty()) ? std::fabs((*ev.el_d0sig)[0]) : 0.0;
    double met_val      = ev.MET / 1000.0;

    hTotalFake[region][0]->Fill(el_pt_val,    fake_weight);
    hTotalFake[region][1]->Fill(el_eta_val,   fake_weight);
    hTotalFake[region][2]->Fill(el_d0_val,    fake_weight);
    hTotalFake[region][3]->Fill(el_d0sig_val, fake_weight);
    hTotalFake[region][4]->Fill(met_val,      fake_weight);

    hFakeBySource[region][0][source]->Fill(el_pt_val,    fake_weight);
    hFakeBySource[region][1][source]->Fill(el_eta_val,   fake_weight);
    hFakeBySource[region][2][source]->Fill(el_d0_val,    fake_weight);
    hFakeBySource[region][3][source]->Fill(el_d0sig_val, fake_weight);
    hFakeBySource[region][4][source]->Fill(met_val,      fake_weight);
  }

  // ------------------------------------------------------------------
  // 6. ATLAS-style fake source composition plotting
  // ------------------------------------------------------------------
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadBottomMargin(0.14);
  gStyle->SetPadTopMargin(0.07);
  gStyle->SetLabelSize(0.045, "xy");
  gStyle->SetTitleSize(0.05, "xy");
  gStyle->SetTitleOffset(1.3, "y");
  gStyle->SetLegendBorderSize(0);

  SourceInfo sources[NumSources] = {
    {kBlue,        20, "bottom"},
    {kMagenta+1,   23, "charm"},
    {kCyan+1,      26, "#gamma conv."},
    {kGreen+2,     25, "light"},
    {kBlack,       22, "#tau"},
    {kOrange+7,    24, "other"}
  };

  auto plotCompositionVar = [&](Region r, int v, const char* xTitle)
  {
    const char* regionLabel = (r == CR) ? "CR" : "SR";
    TH1D*  hTot = hTotalFake[r][v];
    TH1D** hSrc = hFakeBySource[r][v];

    if (!hTot) {
      std::cerr << "plotCompositionVar: hTotalFake[" << r << "][" << v << "] is null\n";
      return;
    }

    TCanvas* c = new TCanvas(Form("c_fake_comp_%s_%s", xTitle, regionLabel), "", 800, 600);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.08);
    c->SetBottomMargin(0.12);
    c->SetLeftMargin(0.14);
    c->SetLogy(true);
    if (v == 1) c->SetLogx(false);
    else        c->SetLogx(true);

    TMultiGraph* mg  = new TMultiGraph();
    TLegend*     leg = new TLegend(0.62, 0.55, 0.88, 0.86);
    leg->SetFillStyle(0); leg->SetBorderSize(0);
    leg->SetTextFont(42); leg->SetTextSize(0.038);

    for (int s = 0; s < NumSources; ++s) {
      if (!hSrc[s]) continue;
      TH1D* h = (TH1D*) hSrc[s]->Clone(Form("%s_clone_%d_%d", hSrc[s]->GetName(), r, v));
      h->SetDirectory(nullptr);
      if (hTot->Integral() > 0.) { h->Divide(hTot); h->Scale(100.0); }
      else                       { h->Reset(); }

      TGraphErrors* g = new TGraphErrors(h);
      g->SetMarkerColor(sources[s].color);
      g->SetMarkerStyle(sources[s].marker);
      g->SetMarkerSize(1.2);
      g->SetLineColor(sources[s].color);
      g->SetLineWidth(2);
      mg->Add(g, "P");
      leg->AddEntry(g, sources[s].label, "P");
    }

    mg->Draw("AP");
    mg->GetYaxis()->SetTitle("Contribution [%]");
    mg->GetYaxis()->SetTitleOffset(1.4);
    mg->GetYaxis()->SetRangeUser(0.1, 290.0);
    mg->GetXaxis()->SetTitle(xTitle);
    mg->GetXaxis()->SetTitleOffset(1.0);

    switch (v) {
      case 0: mg->GetXaxis()->SetLimits(10., 200.); break;
      case 1: mg->GetXaxis()->SetLimits(0.,  2.5);  break;
      case 2: mg->GetXaxis()->SetLimits(0.,  2.0);  break;
      case 3: mg->GetXaxis()->SetLimits(0.,  5.0);  break;
      case 4: mg->GetXaxis()->SetLimits(0.,  200.); break;
      default: break;
    }

    if (v == 0) {
      for (int i = 1; i <= nPtBins; ++i) {
        double x = ptBins[i];
        TLine* l = new TLine(x, 0.1, x, 100.0);
        l->SetLineStyle(7); l->SetLineColor(kGray + 1); l->SetLineWidth(1);
        l->Draw("same");
      }
    }

    leg->Draw("same");

    TLatex lat; lat.SetNDC();
    lat.SetTextFont(72); lat.SetTextSize(0.035);
    lat.DrawLatex(0.18, 0.88, "ATLAS");
    lat.SetTextFont(42); lat.SetTextSize(0.035);
    lat.DrawLatex(0.28, 0.88, "Internal");
    lat.DrawLatex(0.18, 0.83, "#sqrt{s}=13 TeV,140 fb^{-1}");
    lat.DrawLatex(0.18, 0.79, Form("#mu e %s", regionLabel));

    c->RedrawAxis();

    gSystem->mkdir("outputs/05_fake_composition", true);
    c->SaveAs(Form("outputs/05_fake_composition/fake_electron_%s_%s.pdf", xTitle, regionLabel));
    c->SaveAs(Form("outputs/05_fake_composition/fake_electron_%s_%s.png", xTitle, regionLabel));
  };

  const char* xAxisTitle[NumVars] = {
    "p_{T}^{e} [GeV]",
    "|#eta^{e}|",
    "|d_{0}^{e}| [mm]",
    "|d_{0}^{e}| / #sigma_{d_{0}}",
    "E_{T}^{miss} [GeV]"
  };

  for (int v = 0; v < NumVars; ++v) {
    plotCompositionVar(CR, v, xAxisTitle[v]);
    plotCompositionVar(SR, v, xAxisTitle[v]);
  }

  // ------------------------------------------------------------------
  // 7. Summary
  // ------------------------------------------------------------------
  std::cout << "Events in preselection: " << count_preselection << std::endl;
  std::cout << "Events in CR:           " << count_CR           << std::endl;
  std::cout << "Events in SR:           " << count_SR           << std::endl;
}
