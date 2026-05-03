// =============================================================================
// FILE:    studies/sf_distributions.C
// PURPOSE: Visualize how el_ecids_effSF_tight_NOSYS and
//          el_charge_misid_effSF_tight_NOSYS depend on electron pT and |η|.
//          Also splits by truth-prompt vs not-prompt to show where the SFs
//          actually deviate from 1.
//
// INPUT:   MC ntuples (TChain)
// OUTPUT:  outputs/sf_distributions.root
//          outputs/sf_vs_pt_eta.pdf
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include "ROOT/RVec.hxx"
#include <cmath>
#include <iostream>

using ROOT::VecOps::RVec;

// ---------------------------------------------------------------------------
void sf_distributions(const char* outRoot = "outputs/sf_distributions.root",
                      const char* outPdf  = "outputs/sf_vs_pt_eta.pdf")
{
    gStyle->SetOptStat(0);

    // ----- Build chain -----
    const std::string base = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    TChain mc("tree");
    for (auto p : {"ttbar","ztautau","wjets","singletop","diboson","higgs"})
        mc.Add((base + p + "/*.root").c_str());

    // ----- Branch readers -----
    RVec<float>* el_pt          = nullptr;
    RVec<float>* el_eta         = nullptr;
    RVec<int>*   el_iff         = nullptr;
    RVec<int>*   el_isPrompt    = nullptr;
    RVec<float>* sf_ecids       = nullptr;
    RVec<float>* sf_qmis        = nullptr;
    mc.SetBranchAddress("el_pt_NOSYS",                            &el_pt);
    mc.SetBranchAddress("el_eta",                                 &el_eta);
    mc.SetBranchAddress("el_IFFClass",                            &el_iff);
    mc.SetBranchAddress("el_truth_isPrompt",                      &el_isPrompt);
    mc.SetBranchAddress("el_ecids_effSF_tight_NOSYS",             &sf_ecids);
    mc.SetBranchAddress("el_charge_misid_effSF_tight_NOSYS",      &sf_qmis);

    // ----- Histograms -----
    // 1D distributions of SF values (zoom in on deviation from 1)
    TH1D* h1_ecids_all     = new TH1D("h1_ecids_all",     ";ECIDS SF;electrons", 200, 0.95, 1.05);
    TH1D* h1_ecids_prompt  = new TH1D("h1_ecids_prompt",  ";ECIDS SF;electrons", 200, 0.95, 1.05);
    TH1D* h1_ecids_fake    = new TH1D("h1_ecids_fake",    ";ECIDS SF;electrons", 200, 0.95, 1.05);

    TH1D* h1_qmis_all      = new TH1D("h1_qmis_all",      ";charge-misID SF;electrons", 400, 0.9, 2.0);
    TH1D* h1_qmis_prompt   = new TH1D("h1_qmis_prompt",   ";charge-misID SF;electrons", 400, 0.9, 2.0);
    TH1D* h1_qmis_fake     = new TH1D("h1_qmis_fake",     ";charge-misID SF;electrons", 400, 0.9, 2.0);

    // SF vs pT (profile shows the mean per bin)
    TProfile* p_ecids_vs_pt   = new TProfile("p_ecids_vs_pt",
        "ECIDS SF vs p_{T};p_{T}^{e} [GeV];#LTECIDS SF#GT", 30, 0, 200, 0.9, 1.1);
    TProfile* p_qmis_vs_pt    = new TProfile("p_qmis_vs_pt",
        "charge-misID SF vs p_{T};p_{T}^{e} [GeV];#LTcharge-misID SF#GT", 30, 0, 200, 0.5, 3.0);

    // SF vs |eta|
    TProfile* p_ecids_vs_eta  = new TProfile("p_ecids_vs_eta",
        "ECIDS SF vs |#eta|;|#eta|;#LTECIDS SF#GT", 25, 0, 2.5, 0.9, 1.1);
    TProfile* p_qmis_vs_eta   = new TProfile("p_qmis_vs_eta",
        "charge-misID SF vs |#eta|;|#eta|;#LTcharge-misID SF#GT", 25, 0, 2.5, 0.5, 3.0);

    // 2D profile (pT, |η|) — shows the high-pT/high-|η| corner where charge flip lives
    TProfile2D* p2_qmis = new TProfile2D("p2_qmis_vs_pt_eta",
        "charge-misID SF;p_{T}^{e} [GeV];|#eta|", 20, 0, 200, 10, 0, 2.5);

    // ----- Loop -----
    Long64_t n = mc.GetEntries();
    std::cout << "Total entries: " << n << "\n";

    long n_used = 0, n_skipped_sentinel = 0;

    for (Long64_t i = 0; i < n; ++i) {
        if (i % 1000000 == 0) std::cout << "  " << i << " / " << n << "\n";
        mc.GetEntry(i);

        if (!el_pt || el_pt->empty()) continue;
        // Use the leading electron only — matches what AMM and ε_f measurement do.
        const double pt_gev = (*el_pt)[0] * 1e-3;
        const double aeta   = std::abs((*el_eta)[0]);
        if (pt_gev < 5.) continue;

        // SFs may be -1 sentinels (no electron / IFF failed).
        const double sf_e   = (sf_ecids && !sf_ecids->empty()) ? (double)(*sf_ecids)[0] : -1.;
        const double sf_q   = (sf_qmis  && !sf_qmis ->empty()) ? (double)(*sf_qmis )[0] : -1.;
        const bool   ecids_ok = (sf_e > 0. && std::isfinite(sf_e));
        const bool   qmis_ok  = (sf_q > 0. && std::isfinite(sf_q));

        if (!ecids_ok && !qmis_ok) { ++n_skipped_sentinel; continue; }

        const bool is_prompt = (el_isPrompt && !el_isPrompt->empty() && (*el_isPrompt)[0] == 1);

        if (ecids_ok) {
            h1_ecids_all->Fill(sf_e);
            (is_prompt ? h1_ecids_prompt : h1_ecids_fake)->Fill(sf_e);
            p_ecids_vs_pt ->Fill(pt_gev, sf_e);
            p_ecids_vs_eta->Fill(aeta,   sf_e);
        }
        if (qmis_ok) {
            h1_qmis_all->Fill(sf_q);
            (is_prompt ? h1_qmis_prompt : h1_qmis_fake)->Fill(sf_q);
            p_qmis_vs_pt ->Fill(pt_gev, sf_q);
            p_qmis_vs_eta->Fill(aeta,   sf_q);
            p2_qmis      ->Fill(pt_gev, aeta, sf_q);
        }
        ++n_used;
    }

    std::cout << "\nElectrons used: " << n_used
              << "   sentinels skipped: " << n_skipped_sentinel << "\n";

    // ----- Summary numbers -----
    std::cout << "\n=== Mean SFs ===\n";
    std::cout << "  ECIDS  : all=" << h1_ecids_all->GetMean()
              << ",  prompt=" << h1_ecids_prompt->GetMean()
              << ",  fake=" << h1_ecids_fake->GetMean() << "\n";
    std::cout << "  qmis   : all=" << h1_qmis_all->GetMean()
              << ",  prompt=" << h1_qmis_prompt->GetMean()
              << ",  fake=" << h1_qmis_fake->GetMean() << "\n";
    std::cout << "  qmis tail (>1.05): " << h1_qmis_all->Integral(h1_qmis_all->FindBin(1.05),
                                                                   h1_qmis_all->GetNbinsX())
              << " / " << h1_qmis_all->GetEntries() << " electrons\n";

    // ----- Plot -----
    TCanvas c("c","SF distributions", 1600, 1200);
    c.Divide(3, 2);

    auto styleP = [](TProfile* p, int color) {
        p->SetMarkerStyle(20);
        p->SetMarkerColor(color);
        p->SetLineColor(color);
        p->SetLineWidth(2);
    };
    auto styleH = [](TH1D* h, int color, int style=1) {
        h->SetLineColor(color);
        h->SetLineWidth(2);
        h->SetLineStyle(style);
    };

    // Pad 1: ECIDS SF distribution (linear, zoomed)
    c.cd(1); gPad->SetLogy();
    styleH(h1_ecids_all,    kBlack);
    styleH(h1_ecids_prompt, kBlue, 2);
    styleH(h1_ecids_fake,   kRed,  2);
    h1_ecids_all->Draw("HIST");
    h1_ecids_prompt->Draw("HIST SAME");
    h1_ecids_fake  ->Draw("HIST SAME");
    {
        TLegend* l = new TLegend(0.55,0.70,0.88,0.88);
        l->AddEntry(h1_ecids_all,    "all electrons", "l");
        l->AddEntry(h1_ecids_prompt, "truth-prompt",  "l");
        l->AddEntry(h1_ecids_fake,   "non-prompt",    "l");
        l->SetBorderSize(0); l->Draw();
    }

    // Pad 2: charge-misID SF distribution
    c.cd(2); gPad->SetLogy();
    styleH(h1_qmis_all,    kBlack);
    styleH(h1_qmis_prompt, kBlue, 2);
    styleH(h1_qmis_fake,   kRed,  2);
    h1_qmis_all->Draw("HIST");
    h1_qmis_prompt->Draw("HIST SAME");
    h1_qmis_fake  ->Draw("HIST SAME");
    {
        TLegend* l = new TLegend(0.55,0.70,0.88,0.88);
        l->AddEntry(h1_qmis_all,    "all electrons", "l");
        l->AddEntry(h1_qmis_prompt, "truth-prompt",  "l");
        l->AddEntry(h1_qmis_fake,   "non-prompt",    "l");
        l->SetBorderSize(0); l->Draw();
    }

    // Pad 3: 2D charge-misID SF heatmap
    c.cd(3);
    p2_qmis->Draw("COLZ");

    // Pad 4: ECIDS SF vs pT
    c.cd(4);
    styleP(p_ecids_vs_pt, kBlue);
    p_ecids_vs_pt->Draw();

    // Pad 5: charge-misID SF vs pT
    c.cd(5);
    styleP(p_qmis_vs_pt, kRed);
    p_qmis_vs_pt->Draw();

    // Pad 6: charge-misID SF vs |eta|
    c.cd(6);
    styleP(p_qmis_vs_eta, kRed);
    p_qmis_vs_eta->Draw();

    gSystem->mkdir("outputs", true);
    c.SaveAs(outPdf);

    TFile* fo = TFile::Open(outRoot, "RECREATE");
    h1_ecids_all->Write(); h1_ecids_prompt->Write(); h1_ecids_fake->Write();
    h1_qmis_all->Write();  h1_qmis_prompt->Write();  h1_qmis_fake->Write();
    p_ecids_vs_pt->Write(); p_qmis_vs_pt->Write();
    p_ecids_vs_eta->Write(); p_qmis_vs_eta->Write();
    p2_qmis->Write();
    fo->Close();

    std::cout << "\nWrote " << outRoot << "\nWrote " << outPdf << "\n";
}
