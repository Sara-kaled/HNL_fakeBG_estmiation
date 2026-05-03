// =============================================================================
// FILE:    studies/smoothing_test/smooth_fake_eff.C
// PURPOSE: Take the binned nominal fake-rate map (fake_rate_nominal) and
//          produce a SMOOTHED version using a cubic spline fit through bin
//          centres for each |η| slice. The output has a finer pT binning
//          (100 bins) so that the AMM lookup sees a continuous ε_f(pT)
//          instead of step jumps at bin edges.
//
// INPUT:   ../../efficiency/outputs/fake_eff_with_systematics.root
// OUTPUT:  outputs/fake_eff_smoothed.root          — smoothed map
//          outputs/smoothing_diagnostic.pdf         — overlay per |η| slice
// =============================================================================
#include <TFile.h>
#include <TH2F.h>
#include <TGraphErrors.h>
#include <TSpline.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include <iostream>
#include <vector>

void smooth_fake_eff()
{
    TFile* fIn = TFile::Open("../../efficiency/outputs/fake_eff_with_systematics.root", "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "ERROR: cannot open input fake_eff_with_systematics.root\n";
        return;
    }
    TH2F* hRaw = (TH2F*)fIn->Get("fake_rate_nominal");
    if (!hRaw) { std::cerr << "ERROR: fake_rate_nominal not found\n"; return; }

    const int nPt  = hRaw->GetNbinsX();
    const int nEta = hRaw->GetNbinsY();
    const double ptLo = hRaw->GetXaxis()->GetXmin();
    const double ptHi = hRaw->GetXaxis()->GetXmax();
    std::vector<double> etaEdges(nEta+1);
    for (int j = 0; j <= nEta; ++j)
        etaEdges[j] = hRaw->GetYaxis()->GetBinLowEdge(j+1);

    // ----- Output histogram: fine pT binning, same eta binning -----
    const int nPtFine = 100;
    TH2F* hSm = new TH2F("fake_rate_nominal", "Smoothed #varepsilon_{f}(p_{T},|#eta|);"
                         "p_{T}^{e} [GeV];|#eta^{e}|",
                         nPtFine, ptLo, ptHi,
                         nEta, etaEdges.data());

    // ----- Diagnostic canvas -----
    gSystem->mkdir("outputs", true);
    TCanvas c("c","Smoothing comparison", 1200, 800);
    c.Divide(3, 2);

    // ----- Per-|η| spline fit -----
    for (int j = 1; j <= nEta; ++j) {
        // Collect (bin centre pT, ε_f) for this |η| slice
        std::vector<double> xv, yv, ex, ey;
        for (int i = 1; i <= nPt; ++i) {
            xv.push_back(hRaw->GetXaxis()->GetBinCenter(i));
            yv.push_back(hRaw->GetBinContent(i, j));
            ex.push_back(0.);
            ey.push_back(hRaw->GetBinError(i, j));
        }

        // Build cubic spline through the bin centres
        TSpline3* spl = new TSpline3(Form("spl_etabin%d", j),
                                     xv.data(), yv.data(), (int)xv.size());

        // Fill the fine-binning histogram by evaluating the spline at every
        // bin centre, clamped to [first,last] bin centre to avoid wild
        // extrapolation outside the measured pT range.
        const double xMin = xv.front();
        const double xMax = xv.back();
        for (int i = 1; i <= nPtFine; ++i) {
            double x = hSm->GetXaxis()->GetBinCenter(i);
            double xClamp = std::max(xMin, std::min(xMax, x));
            double yEval = spl->Eval(xClamp);
            yEval = std::max(0.001, std::min(0.99, yEval));   // guard for AMM
            hSm->SetBinContent(i, j, yEval);
            // Approximate the per-bin error via the nearest measured bin's err
            int iRaw = hRaw->GetXaxis()->FindBin(xClamp);
            iRaw = std::max(1, std::min(iRaw, nPt));
            hSm->SetBinError  (i, j, hRaw->GetBinError(iRaw, j));
        }

        // ----- Overlay: raw points + smoothed line -----
        c.cd(j);
        TGraphErrors* g = new TGraphErrors((int)xv.size(),
                                           xv.data(), yv.data(),
                                           ex.data(), ey.data());
        g->SetTitle(Form("|#eta| in [%.2f, %.2f];p_{T} [GeV];#varepsilon_{f}",
                        etaEdges[j-1], etaEdges[j]));
        g->SetMarkerStyle(20);
        g->SetMarkerColor(kBlack);
        g->Draw("AP");

        // Smoothed curve at fine grid
        const int nDraw = 200;
        std::vector<double> xd(nDraw), yd(nDraw);
        for (int k = 0; k < nDraw; ++k) {
            double xx = ptLo + (ptHi - ptLo) * (k + 0.5) / nDraw;
            double xClamp = std::max(xv.front(), std::min(xv.back(), xx));
            xd[k] = xx;
            yd[k] = std::max(0.001, std::min(0.99, spl->Eval(xClamp)));
        }
        TGraph* gs = new TGraph(nDraw, xd.data(), yd.data());
        gs->SetLineColor(kRed);
        gs->SetLineWidth(2);
        gs->Draw("L SAME");

        TLegend* leg = new TLegend(0.55, 0.20, 0.85, 0.35);
        leg->AddEntry(g, "raw bins", "ep");
        leg->AddEntry(gs, "spline smooth", "l");
        leg->SetBorderSize(0);
        leg->Draw();
    }

    c.SaveAs("outputs/smoothing_diagnostic.pdf");

    TFile* fOut = TFile::Open("outputs/fake_eff_smoothed.root", "RECREATE");
    hSm->Write();
    fOut->Close();
    fIn->Close();
    std::cout << "\n✓ Wrote outputs/fake_eff_smoothed.root\n";
    std::cout << "  → " << nPtFine << " pT bins × " << nEta << " |η| bins\n";
    std::cout << "✓ Wrote outputs/smoothing_diagnostic.pdf\n";
}
