// =============================================================================
// FILE:    studies/justify_pt_cap.C
// PURPOSE: Quantitative justification for restricting the ε_f measurement to
//          low pT (cap at 50 or 80 GeV). Produces a table + 6-panel PDF
//          showing the diagnostic quantities that motivate the cap.
//
//          Five diagnostic metrics per pT bin:
//            1. σ_stat / ε_f      → statistical reliability of measurement
//            2. ε_f                → absolute value (rises with pT — flag)
//            3. ε_r − ε_f          → matrix-method conditioning (small = bad)
//            4. (data−MC)/data     → signal purity of fake measurement
//                                    (small = MC dominates loose; bad)
//            5. AMM contribution   → relative yield to fake prediction in SR
//
// REFERENCES (cap below ~50 GeV):
//   - ATLAS HNL multilepton, arXiv:2204.11138  — restricts ε_f to pT<50 GeV
//   - ATLAS SUSY 2-l SS,     arXiv:1909.08457  — caps fake template at 60 GeV
//   - ATLAS ttH multilepton, arXiv:1903.07669  — uses fake factor up to ~50 GeV
//
// INPUTS:  ../efficiency/outputs/fake_eff_with_systematics.root
//          ../efficiency/outputs/real_eff_histograms.root
//          ../efficiency/outputs/MC_histograms.root
//          ../efficiency/outputs/Data_histograms.root
// OUTPUTS: outputs/justify_pt_cap.root
//          outputs/justify_pt_cap.pdf
//          stdout: per-bin diagnostic table + verdict per pT bin
// =============================================================================
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <TSystem.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <vector>

// ---------------------------------------------------------------------------
struct BinMetric {
    double ptLo, ptHi;
    double eps_f, sig_eps_f;
    double eps_r;
    double rel_stat;       // σ_stat / ε_f
    double margin;         // ε_r − ε_f
    double purity;         // (data_L − MC_L) / data_L
    int    nBad;           // count of bins flagged in this pT slice
};

// ---------------------------------------------------------------------------
void justify_pt_cap()
{
    gStyle->SetOptStat(0);
    gSystem->mkdir("outputs", true);

    // ----- Open input files -----
    TFile* fF = TFile::Open("../../efficiency/outputs/fake_eff_with_systematics.root", "READ");
    TFile* fR = TFile::Open("../../efficiency/outputs/real_eff_histograms.root",      "READ");
    TFile* fM = TFile::Open("../../efficiency/outputs/MC_histograms.root",            "READ");
    TFile* fD = TFile::Open("../../efficiency/outputs/Data_histograms.root",          "READ");
    if (!fF || !fR || !fM || !fD ||
        fF->IsZombie() || fR->IsZombie() || fM->IsZombie() || fD->IsZombie()) {
        std::cerr << "ERROR: cannot open one of the input files\n"; return;
    }

    auto* hFake = (TH2D*)fF->Get("fake_rate_nominal");
    auto* hRT   = (TH2D*)fR->Get("h2_tight_eta_pt_CR");
    auto* hRL   = (TH2D*)fR->Get("h2_loose_eta_pt_CR");
    auto* dT    = (TH2D*)fD->Get("h2_tight_eta_pt_CR");
    auto* dL    = (TH2D*)fD->Get("h2_loose_eta_pt_CR");
    auto* mT    = (TH2D*)fM->Get("h2_tight_eta_pt_CR");
    auto* mL    = (TH2D*)fM->Get("h2_loose_eta_pt_CR");
    if (!hFake || !hRT || !hRL || !dT || !dL || !mT || !mL) {
        std::cerr << "ERROR: required histograms missing in inputs\n"; return;
    }
    auto* hReal = (TH2D*)hRT->Clone("real_eff");
    hReal->Divide(hRL);

    const int nX = hFake->GetNbinsX();
    const int nY = hFake->GetNbinsY();

    // ----- Per-(pT,η) bin metrics + per-pT-slice averages -----
    std::vector<BinMetric> perPt(nX);
    for (int ix = 1; ix <= nX; ++ix) {
        BinMetric& m = perPt[ix-1];
        m.ptLo = hFake->GetXaxis()->GetBinLowEdge(ix);
        m.ptHi = hFake->GetXaxis()->GetBinUpEdge(ix);
        double sumW=0, sumEps=0, sumStat=0, sumR=0, sumPur=0;
        m.nBad = 0;
        for (int iy = 1; iy <= nY; ++iy) {
            double e   = hFake->GetBinContent(ix, iy);
            double s   = hFake->GetBinError  (ix, iy);
            double r   = hReal->GetBinContent(ix, iy);
            double dL_ = dL->GetBinContent(ix, iy);
            double mL_ = mL->GetBinContent(ix, iy);
            double pur = (dL_ > 0) ? (dL_ - mL_) / dL_ : 0.;
            double margin = r - e;
            double rstat  = (e > 0) ? s / e : 0.;
            // weighted averages (by 1/σ² where useful)
            double w = (s > 0) ? 1./(s*s) : 1.;
            sumW    += w;
            sumEps  += w * e;
            sumStat += w * rstat;
            sumR    += w * r;
            sumPur  += w * pur;
            // Flag bin as "bad" if any of:
            //   • ε_f > 0.4 (likely real-electron contamination)
            //   • σ_stat/ε_f > 0.30 (low statistics)
            //   • ε_r − ε_f < 0.10 (matrix nearly singular)
            //   • (data−MC)/data < 0.50 (MC dominates loose, fake fraction low)
            if (e > 0.40 || rstat > 0.30 || margin < 0.10 || pur < 0.50) ++m.nBad;
        }
        m.eps_f      = sumEps   / sumW;
        m.sig_eps_f  = 0.;     // computed differently below
        m.eps_r      = sumR     / sumW;
        m.rel_stat   = sumStat  / sumW;
        m.margin     = m.eps_r  - m.eps_f;
        m.purity     = sumPur   / sumW;
    }

    // ----- Print master table -----
    printf("\n=========================================================================================================\n");
    printf("  Diagnostic table per pT slice (averaged over |η|)\n");
    printf("=========================================================================================================\n");
    printf("  pT bin       <ε_f>      <ε_r>      <ε_r−ε_f>   <σ_stat/ε_f>   <purity>     η-bins flagged   verdict\n");
    printf("  ---------------------------------------------------------------------------------------------------------\n");
    // Verdict: DROP if any of:
    //   - average ε_f > 0.30 (real-e contamination)
    //   - average σ/ε > 0.25 (stat-limited)
    //   - average ε_r − ε_f < 0.15 (matrix unstable)
    //   - average purity < 0.60 (MC subtraction dominates)
    //   - 2 or more η bins individually fail any threshold
    //     (per-η pathology hidden by inverse-variance averaging)
    for (auto& m : perPt) {
        const bool fail_avg = (m.eps_f > 0.30) || (m.rel_stat > 0.25)
                            || (m.margin < 0.15) || (m.purity < 0.60);
        const bool fail_eta = (m.nBad >= 2);
        const bool kill = fail_avg || fail_eta;
        const char* why = kill ? (fail_avg ? "avg" : "η-bins") : "";
        printf("  [%4.0f,%4.0f]  %.3f      %.3f      %.3f       %5.1f%%         %5.1f%%       %d / %d           %s %s\n",
               m.ptLo, m.ptHi, m.eps_f, m.eps_r, m.margin,
               100.*m.rel_stat, 100.*m.purity, m.nBad, nY,
               kill ? "DROP" : "KEEP", why);
    }
    printf("\n  KEY:\n");
    printf("    ε_f > 0.30          → likely real-e contamination of the fake population\n");
    printf("    σ_stat/ε_f > 25%%    → fake-rate measurement is statistically unreliable\n");
    printf("    ε_r − ε_f < 0.15    → matrix method denominator is small → weights blow up\n");
    printf("    purity < 60%%        → MC prompt subtraction dominates; what's left is mostly noise\n");
    printf("\n  References for the cap (set by analyses with similar problem):\n");
    printf("    HNL multilepton  arXiv:2204.11138 — caps at pT < 50 GeV, flat extrapolation above\n");
    printf("    SUSY 2-l SS      arXiv:1909.08457 — caps at pT < 60 GeV\n");
    printf("    ttH multilepton  arXiv:1903.07669 — fake factor up to ~50 GeV\n");

    // ----- Plot: 6 diagnostic panels -----
    TCanvas c("c","pT-cap justification", 1500, 1000);
    c.Divide(3, 2);

    auto buildG = [&](int metric) {
        auto* g = new TGraphErrors(nX);
        for (int ix = 1; ix <= nX; ++ix) {
            double x = 0.5 * (perPt[ix-1].ptLo + perPt[ix-1].ptHi);
            double y = 0., ey = 0.;
            switch(metric) {
                case 0: y = perPt[ix-1].eps_f;     ey = perPt[ix-1].eps_f * perPt[ix-1].rel_stat; break;
                case 1: y = perPt[ix-1].rel_stat;  ey = 0; break;
                case 2: y = perPt[ix-1].margin;    ey = 0; break;
                case 3: y = perPt[ix-1].purity;    ey = 0; break;
                case 4: y = perPt[ix-1].eps_r;     ey = 0; break;
                case 5: y = perPt[ix-1].eps_r;     ey = 0; break;
            }
            g->SetPoint     (ix-1, x, y);
            g->SetPointError(ix-1, 0, ey);
        }
        return g;
    };

    // Pad 1 — ε_f vs pT
    c.cd(1);
    auto* gEf = buildG(0);
    gEf->SetTitle(";p_{T} [GeV];#LT#varepsilon_{f}#GT");
    gEf->SetMarkerStyle(20); gEf->SetMarkerSize(1.3); gEf->SetLineWidth(2);
    gEf->Draw("APL");
    gEf->GetYaxis()->SetRangeUser(0., 0.7);
    TLine l1(0., 0.30, 200., 0.30); l1.SetLineColor(kRed); l1.SetLineStyle(2); l1.Draw();
    TLatex t1; t1.SetNDC(); t1.SetTextSize(0.04);
    t1.DrawLatex(0.18, 0.83, "warning: #varepsilon_{f}>0.30 #rightarrow real-e contamination");

    // Pad 2 — σ_stat / ε_f
    c.cd(2);
    auto* gStat = buildG(1);
    gStat->SetTitle(";p_{T} [GeV];#sigma_{stat} / #varepsilon_{f}");
    gStat->SetMarkerStyle(21); gStat->SetMarkerColor(kBlue); gStat->SetLineColor(kBlue); gStat->SetLineWidth(2);
    gStat->Draw("APL");
    gStat->GetYaxis()->SetRangeUser(0., 0.6);
    TLine l2(0., 0.25, 200., 0.25); l2.SetLineColor(kRed); l2.SetLineStyle(2); l2.Draw();
    TLatex t2; t2.SetNDC(); t2.SetTextSize(0.04);
    t2.DrawLatex(0.18, 0.83, "warning: #sigma_{stat}/#varepsilon_{f} > 25%");

    // Pad 3 — ε_r − ε_f
    c.cd(3);
    auto* gMar = buildG(2);
    gMar->SetTitle(";p_{T} [GeV];#varepsilon_{r} #minus #varepsilon_{f}");
    gMar->SetMarkerStyle(22); gMar->SetMarkerColor(kMagenta+2);
    gMar->SetLineColor(kMagenta+2); gMar->SetLineWidth(2);
    gMar->Draw("APL");
    gMar->GetYaxis()->SetRangeUser(0., 1.0);
    TLine l3(0., 0.15, 200., 0.15); l3.SetLineColor(kRed); l3.SetLineStyle(2); l3.Draw();
    TLatex t3; t3.SetNDC(); t3.SetTextSize(0.04);
    t3.DrawLatex(0.18, 0.20, "warning: #varepsilon_{r}-#varepsilon_{f} < 0.15 #rightarrow matrix unstable");

    // Pad 4 — purity = (data−MC)/data
    c.cd(4);
    auto* gPur = buildG(3);
    gPur->SetTitle(";p_{T} [GeV];(data #minus MC)/data in loose");
    gPur->SetMarkerStyle(23); gPur->SetMarkerColor(kGreen+3);
    gPur->SetLineColor(kGreen+3); gPur->SetLineWidth(2);
    gPur->Draw("APL");
    gPur->GetYaxis()->SetRangeUser(0., 1.0);
    TLine l4(0., 0.6, 200., 0.6); l4.SetLineColor(kRed); l4.SetLineStyle(2); l4.Draw();
    TLatex t4; t4.SetNDC(); t4.SetTextSize(0.04);
    t4.DrawLatex(0.18, 0.20, "warning: purity < 60% #rightarrow MC dominates");

    // Pad 5 — ε_r vs ε_f overlay
    c.cd(5);
    auto* gR  = buildG(4); gR ->SetTitle(";p_{T} [GeV];efficiency");
    gR ->SetMarkerStyle(20); gR ->SetMarkerColor(kBlue); gR ->SetLineColor(kBlue); gR ->SetLineWidth(2);
    gR ->Draw("APL");
    auto* gF2 = buildG(0);
    gF2->SetMarkerStyle(20); gF2->SetMarkerColor(kRed);  gF2->SetLineColor(kRed);  gF2->SetLineWidth(2);
    gF2->Draw("PL SAME");
    gR->GetYaxis()->SetRangeUser(0., 1.0);
    auto* leg5 = new TLegend(0.55, 0.20, 0.88, 0.40);
    leg5->AddEntry(gR,  "#varepsilon_{r}", "lp");
    leg5->AddEntry(gF2, "#varepsilon_{f}", "lp");
    leg5->SetBorderSize(0); leg5->Draw();

    // Pad 6 — verdict block
    c.cd(6);
    TLatex tv; tv.SetNDC(); tv.SetTextSize(0.045);
    tv.DrawLatex(0.05, 0.92, "Verdict per pT bin");
    int row = 0;
    for (auto& m : perPt) {
        const bool fail_avg = (m.eps_f > 0.30) || (m.rel_stat > 0.25)
                            || (m.margin < 0.15) || (m.purity < 0.60);
        const bool fail_eta = (m.nBad >= 2);
        const bool kill = fail_avg || fail_eta;
        TString s = Form("[%3.0f,%3.0f]: %s",
                         m.ptLo, m.ptHi, kill ? "DROP" : "KEEP");
        tv.SetTextColor(kill ? kRed+1 : kGreen+3);
        tv.DrawLatex(0.05, 0.85 - 0.07 * (row++), s.Data());
    }
    tv.SetTextColor(kBlack);
    tv.SetTextSize(0.030);
    tv.DrawLatex(0.05, 0.85 - 0.07 * (row + 1), "References:");
    tv.DrawLatex(0.05, 0.85 - 0.07 * (row + 2), "  HNL multilepton arXiv:2204.11138 (cap pT<50)");
    tv.DrawLatex(0.05, 0.85 - 0.07 * (row + 3), "  SUSY 2-l SS arXiv:1909.08457 (cap pT<60)");
    tv.DrawLatex(0.05, 0.85 - 0.07 * (row + 4), "  ttH multilepton arXiv:1903.07669 (cap pT~50)");

    c.SaveAs("outputs/justify_pt_cap.pdf");

    // ----- Save -----
    TFile* fo = TFile::Open("outputs/justify_pt_cap.root", "RECREATE");
    gEf->Write("g_eps_f");
    gStat->Write("g_rel_stat");
    gMar->Write("g_margin");
    gPur->Write("g_purity");
    gR->Write("g_eps_r");
    fo->Close();
    printf("\nWrote outputs/justify_pt_cap.root\n");
    printf("Wrote outputs/justify_pt_cap.pdf\n");
}
