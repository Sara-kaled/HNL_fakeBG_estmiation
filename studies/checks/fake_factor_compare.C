// =============================================================================
// FILE:    studies/fake_factor_compare.C
// PURPOSE: Cross-check the matrix-method ε_f against a fake-factor formulation
//          f = ε_f / (1 − ε_f). The fake-factor avoids the (ε_r − ε_f)
//          denominator that blows up at high pT, so it's a robustness check
//          on the matrix method. If the two prescriptions diverge in their
//          per-event weights, the matrix method is sensitive to ε_r modeling
//          and the fake factor's prediction is the safer one.
//
//          Computes for each (pT, |η|) bin:
//             - ε_f   (matrix-method efficiency)
//             - ε_r   (real efficiency)
//             - f = ε_f/(1-ε_f)              (fake factor)
//             - w_L_AMM = ε_f·ε_r/(ε_r-ε_f)  (AMM loose weight)
//             - w_L_FF  = f                   (fake-factor loose weight)
//             - ratio AMM/FF
//
// OUTPUT:  outputs/fake_factor_compare.root
//          outputs/fake_factor_compare.pdf
// =============================================================================
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>
#include <cstdio>

#include "../../utils/SelectionUtils.h"

void fake_factor_compare()
{
    gStyle->SetOptStat(0);
    gSystem->mkdir("outputs", true);

    TFile* fF = TFile::Open("../../efficiency/outputs/fake_eff_with_systematics.root", "READ");
    TFile* fR = TFile::Open("../../efficiency/outputs/real_eff_histograms.root",      "READ");
    if (!fF || !fR || fF->IsZombie() || fR->IsZombie()) {
        std::cerr << "ERROR: cannot open input files\n"; return;
    }

    auto* hFake = (TH2D*)fF->Get("fake_rate_nominal");
    auto* hRT   = (TH2D*)fR->Get("h2_tight_eta_pt_CR");
    auto* hRL   = (TH2D*)fR->Get("h2_loose_eta_pt_CR");
    if (!hFake || !hRT || !hRL) {
        std::cerr << "ERROR: histograms missing\n"; return;
    }
    auto* hReal = (TH2D*)hRT->Clone("real_eff");
    hReal->Divide(hRL);

    const int nX = hFake->GetNbinsX();
    const int nY = hFake->GetNbinsY();

    // Per-bin computed quantities
    auto* hFF    = (TH2D*)hFake->Clone("fake_factor");
    auto* hWL_AMM = (TH2D*)hFake->Clone("wL_AMM");
    auto* hWL_FF  = (TH2D*)hFake->Clone("wL_FF");
    auto* hRatio = (TH2D*)hFake->Clone("ratio_AMM_over_FF");
    hFF->Reset(); hWL_AMM->Reset(); hWL_FF->Reset(); hRatio->Reset();

    printf("\n========================================================================\n");
    printf("  Fake-factor vs matrix-method per-bin comparison\n");
    printf("========================================================================\n");
    printf("  pT [GeV]      |η|       ε_f      ε_r      f=ε_f/(1-ε_f)   w_L^AMM   w_L^FF    AMM/FF\n");
    printf("  ----------------------------------------------------------------\n");

    for (int ix = 1; ix <= nX; ++ix) {
        for (int iy = 1; iy <= nY; ++iy) {
            double e   = hFake->GetBinContent(ix, iy);
            double r   = hReal->GetBinContent(ix, iy);
            // Clamp to physical range for safety
            e = std::max(1e-5, std::min(0.99,    e));
            r = std::max(1e-5, std::min(1.0-1e-5, r));
            double f   = e / (1.0 - e);
            double wAMM = (r > e + 1e-9) ? (e * r) / (r - e) : 0.;
            double wFF  = f;
            double rAMM_FF = (wFF > 0) ? wAMM / wFF : 0.;
            hFF    ->SetBinContent(ix, iy, f);
            hWL_AMM->SetBinContent(ix, iy, wAMM);
            hWL_FF ->SetBinContent(ix, iy, wFF);
            hRatio ->SetBinContent(ix, iy, rAMM_FF);
            printf("  [%4.0f,%4.0f]  [%.2f,%.2f]  %.4f   %.4f   %.4f          %.4f    %.4f    %.2f\n",
                   hFake->GetXaxis()->GetBinLowEdge(ix),
                   hFake->GetXaxis()->GetBinUpEdge(ix),
                   hFake->GetYaxis()->GetBinLowEdge(iy),
                   hFake->GetYaxis()->GetBinUpEdge(iy),
                   e, r, f, wAMM, wFF, rAMM_FF);
        }
    }
    printf("\n  Interpretation:\n");
    printf("    AMM/FF ~ 1: matrix method and fake factor agree → both prescriptions OK.\n");
    printf("    AMM/FF >> 1: matrix method amplifies the loose weight — sensitive to\n");
    printf("                 ε_r modeling. Bins with AMM/FF > 1.5 are matrix-method\n");
    printf("                 dominated and the fake-factor prediction is the safer\n");
    printf("                 one (independent of ε_r).\n");
    printf("    AMM < FF (ratio < 1): rare; check ε_r > ε_f condition.\n");

    // ----- Plot -----
    TCanvas c("c","FF vs AMM", 1500, 600);
    c.Divide(3, 1);

    // Per-pT averages
    auto avgPerPt = [&](TH2D* h) {
        TH1D* hp = new TH1D(Form("avg_%s", h->GetName()),
                            ";p_{T} [GeV];", nX, kPtBins);
        for (int ix = 1; ix <= nX; ++ix) {
            double s = 0., n = 0;
            for (int iy = 1; iy <= nY; ++iy) {
                double v = h->GetBinContent(ix, iy);
                if (v > 0) { s += v; ++n; }
            }
            if (n > 0) hp->SetBinContent(ix, s/n);
        }
        return hp;
    };

    c.cd(1);
    TH1D* a_eps  = avgPerPt(hFake);  a_eps->SetTitle(";p_{T} [GeV];#LT#varepsilon_{f}#GT");
    TH1D* a_real = avgPerPt(hReal);  a_real->SetTitle(";p_{T} [GeV];#LT#varepsilon_{r}#GT");
    a_eps->SetLineColor(kRed); a_eps->SetLineWidth(2); a_eps->SetMarkerStyle(20); a_eps->SetMarkerColor(kRed);
    a_real->SetLineColor(kBlue); a_real->SetLineWidth(2); a_real->SetMarkerStyle(21); a_real->SetMarkerColor(kBlue);
    a_real->Draw("PL"); a_real->GetYaxis()->SetRangeUser(0., 1.0);
    a_eps->Draw("PL SAME");
    {
        auto* l = new TLegend(0.6, 0.18, 0.85, 0.32);
        l->AddEntry(a_eps,  "#varepsilon_{f}", "lp");
        l->AddEntry(a_real, "#varepsilon_{r}", "lp");
        l->SetBorderSize(0); l->Draw();
    }

    c.cd(2);
    TH1D* a_AMM = avgPerPt(hWL_AMM); a_AMM->SetTitle(";p_{T} [GeV];#LTw_{L}#GT");
    TH1D* a_FF  = avgPerPt(hWL_FF);
    a_AMM->SetLineColor(kRed); a_AMM->SetLineWidth(2); a_AMM->SetMarkerStyle(20); a_AMM->SetMarkerColor(kRed);
    a_FF ->SetLineColor(kBlue); a_FF ->SetLineWidth(2); a_FF ->SetMarkerStyle(21); a_FF ->SetMarkerColor(kBlue);
    a_AMM->Draw("PL"); a_AMM->GetYaxis()->SetRangeUser(0., 1.5*a_AMM->GetMaximum());
    a_FF ->Draw("PL SAME");
    {
        auto* l = new TLegend(0.18, 0.70, 0.45, 0.85);
        l->AddEntry(a_AMM, "AMM w_{L}",         "lp");
        l->AddEntry(a_FF,  "fake-factor w_{L}", "lp");
        l->SetBorderSize(0); l->Draw();
    }

    c.cd(3);
    TH1D* a_rat = avgPerPt(hRatio); a_rat->SetTitle(";p_{T} [GeV];AMM / FF ratio");
    a_rat->SetLineColor(kBlack); a_rat->SetLineWidth(2); a_rat->SetMarkerStyle(20);
    a_rat->Draw("PL");

    c.SaveAs("outputs/fake_factor_compare.pdf");

    TFile* fo = TFile::Open("outputs/fake_factor_compare.root", "RECREATE");
    hFake->Write(); hReal->Write();
    hFF->Write(); hWL_AMM->Write(); hWL_FF->Write(); hRatio->Write();
    fo->Close();
    printf("\nWrote outputs/fake_factor_compare.{root,pdf}\n");
}
