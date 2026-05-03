// =============================================================================
// FILE:    studies/met_split_test.C
// PURPOSE: Cleaner alternative to anti-iso-μ multijet test. Splits the
//          standard SS μe fake CR by event MET:
//             low-MET  (MET < 30 GeV)  → multijet-enriched
//             high-MET (MET > 50 GeV)  → multijet-depleted
//          Compares ε_f(pT, |η|) between the two. If they agree, no multijet
//          contamination. If low-MET ε_f is significantly higher, multijet
//          IS contaminating your standard CR and you can restrict the fake
//          measurement to high-MET to clean it up.
//
//          Both regions use the same MC prompt subtraction (no anti-iso muon
//          population issues), so statistics are decent in both.
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <cmath>
#include <cstdio>
#include <iostream>

#include "../../utils/SelectionUtils.h"

static inline void fillEle(TH2D* hT, TH2D* hL, const EventData& ev, double w)
{
    const double pt_gev = ev.l2_pt * 1e-3;
    if (ev.el_is_tight) hT->Fill(pt_gev, ev.el_aeta, w);
    if (ev.el_is_loose) hL->Fill(pt_gev, ev.el_aeta, w);
}

static TH2D* computeEpsF(TH2D* dT, TH2D* dL, TH2D* mT, TH2D* mL, const char* name)
{
    auto* h = (TH2D*)dT->Clone(name);
    h->Reset();
    for (int ix = 1; ix <= dT->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= dT->GetNbinsY(); ++iy) {
            double nT = dT->GetBinContent(ix, iy) - mT->GetBinContent(ix, iy);
            double nL = dL->GetBinContent(ix, iy) - mL->GetBinContent(ix, iy);
            double eT = std::sqrt(dT->GetBinError(ix,iy)*dT->GetBinError(ix,iy)
                                + mT->GetBinError(ix,iy)*mT->GetBinError(ix,iy));
            double eL = std::sqrt(dL->GetBinError(ix,iy)*dL->GetBinError(ix,iy)
                                + mL->GetBinError(ix,iy)*mL->GetBinError(ix,iy));
            if (nL <= 0 || nT < 0) {
                h->SetBinContent(ix, iy, 0.);
                h->SetBinError  (ix, iy, 0.);
                continue;
            }
            double r = nT / nL;
            double rel2 = (nT > 0 ? eT*eT/(nT*nT) : 0.) + eL*eL/(nL*nL);
            // Cap at 1.0 to flag unphysical values
            r = std::min(1.0, r);
            h->SetBinContent(ix, iy, r);
            h->SetBinError  (ix, iy, r * std::sqrt(rel2));
        }
    }
    return h;
}

// ---------------------------------------------------------------------------
void met_split_test(double met_low_max  = 30000.,   // < 30 GeV → multijet-rich
                    double met_high_min = 50000.)   // > 50 GeV → multijet-poor
{
    gStyle->SetOptStat(0);
    gSystem->mkdir("outputs", true);

    TChain dataCh("tree");
    dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data15.root");
    dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data16.root");
    dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_d/data17.root");
    dataCh.Add("/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/data18.root");

    TChain mcCh("tree");
    const std::string base = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    for (auto p : {"ttbar","ztautau","wjets","singletop","diboson","higgs"})
        mcCh.Add((base + p + "/*.root").c_str());

    auto book = [](const char* name) {
        return new TH2D(name, ";p_{T}^{e} [GeV];|#eta^{e}|",
                        kNPtBins, kPtBins, kNEtaBins, kEtaBins);
    };
    TH2D *dT_lo = book("dT_lo"), *dL_lo = book("dL_lo"),
         *mT_lo = book("mT_lo"), *mL_lo = book("mL_lo");
    TH2D *dT_hi = book("dT_hi"), *dL_hi = book("dL_hi"),
         *mT_hi = book("mT_hi"), *mL_hi = book("mL_hi");
    for (auto* h : {dT_lo,dL_lo,mT_lo,mL_lo,dT_hi,dL_hi,mT_hi,mL_hi}) h->Sumw2();

    auto loop = [&](TTree* t, bool isMC, const char* tag) {
        std::cout << "\n=== Loop " << tag << " ===\n";
        EventData ev;
        SetupBranches(t, ev);
        const Long64_t n = t->GetEntries();
        Long64_t nLo = 0, nHi = 0;
        for (Long64_t i = 0; i < n; ++i) {
            if (i % 1000000 == 0) std::cout << "  " << i << " / " << n << "\n";
            t->GetEntry(i);
            ComputeDerived(ev);

            if (!passBasePresel_FakeMeasurement(ev)) continue;
            if (!isFakeCR(ev)) continue;
            if (!(ev.jet1_pt > 0. && ev.jet2_pt > 0.)) continue;
            if (!ev.no_bjets) continue;

            if (isMC && !isMCPromptElectron(ev)) continue;

            const double w = isMC ? mcWeightLoose(ev) : 1.0;

            if (ev.MET < met_low_max) {
                fillEle(isMC ? mT_lo : dT_lo, isMC ? mL_lo : dL_lo, ev, w);
                ++nLo;
            } else if (ev.MET > met_high_min) {
                fillEle(isMC ? mT_hi : dT_hi, isMC ? mL_hi : dL_hi, ev, w);
                ++nHi;
            }
        }
        std::cout << "  filled  low-MET: " << nLo
                  << "  high-MET: " << nHi << "\n";
    };

    loop(&mcCh,   true,  "MC");
    loop(&dataCh, false, "DATA");

    TH2D* eps_lo = computeEpsF(dT_lo, dL_lo, mT_lo, mL_lo, "eps_lo");
    TH2D* eps_hi = computeEpsF(dT_hi, dL_hi, mT_hi, mL_hi, "eps_hi");

    printf("\n========================================================================\n");
    printf("  MET-split test: ε_f in low-MET (multijet-rich) vs high-MET (clean)\n");
    printf("    low-MET cut:  MET < %.0f GeV\n", met_low_max  * 1e-3);
    printf("    high-MET cut: MET > %.0f GeV\n", met_high_min * 1e-3);
    printf("========================================================================\n");
    printf("  pT [GeV]      |η|        ε_f (low-MET)    ε_f (high-MET)   ratio\n");
    printf("  ----------------------------------------------------------------\n");
    int nValid = 0;
    double sumRatio = 0., maxRatio = 0.;
    for (int ix = 1; ix <= eps_lo->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= eps_lo->GetNbinsY(); ++iy) {
            double e1 = eps_lo->GetBinContent(ix, iy);
            double s1 = eps_lo->GetBinError  (ix, iy);
            double e2 = eps_hi->GetBinContent(ix, iy);
            double s2 = eps_hi->GetBinError  (ix, iy);
            double r  = (e1 > 0 && e2 > 0) ? e1 / e2 : 0.;
            double ptLo = eps_lo->GetXaxis()->GetBinLowEdge(ix);
            double ptHi = eps_lo->GetXaxis()->GetBinUpEdge(ix);
            double etLo = eps_lo->GetYaxis()->GetBinLowEdge(iy);
            double etHi = eps_lo->GetYaxis()->GetBinUpEdge(iy);
            printf("  [%4.0f,%4.0f]  [%.2f,%.2f]  %.4f±%.4f   %.4f±%.4f   %.2f\n",
                   ptLo, ptHi, etLo, etHi, e1, s1, e2, s2, r);
            if (r > 0) { sumRatio += r; ++nValid; if (r > maxRatio) maxRatio = r; }
        }
    }
    printf("  ----------------------------------------------------------------\n");
    if (nValid > 0)
        printf("  Average ratio (low-MET / high-MET): %.2f   max bin: %.2f\n",
               sumRatio / nValid, maxRatio);
    printf("\n  Interpretation:\n");
    printf("    ratio ~ 1: no multijet contamination — your CR is clean.\n");
    printf("    ratio > 1.3: low-MET region has higher ε_f → multijet IS in CR.\n");
    printf("                 Apply MET > %.0f GeV to your fake-rate measurement.\n",
           met_high_min * 1e-3);
    printf("    ratio < 0.7: low-MET has LOWER ε_f than high-MET — unusual,\n");
    printf("                 might indicate over-subtraction at low-MET.\n");

    TCanvas c("c","MET split test", 1500, 900);
    c.Divide(3, 2);
    for (int iy = 1; iy <= kNEtaBins; ++iy) {
        c.cd(iy);
        TH1D* h1 = eps_lo->ProjectionX(Form("lo_etabin%d", iy), iy, iy);
        TH1D* h2 = eps_hi->ProjectionX(Form("hi_etabin%d", iy), iy, iy);
        for (int ix = 1; ix <= eps_lo->GetNbinsX(); ++ix) {
            h1->SetBinContent(ix, eps_lo->GetBinContent(ix, iy));
            h1->SetBinError  (ix, eps_lo->GetBinError  (ix, iy));
            h2->SetBinContent(ix, eps_hi->GetBinContent(ix, iy));
            h2->SetBinError  (ix, eps_hi->GetBinError  (ix, iy));
        }
        h1->SetTitle(Form("|#eta| in [%.2f, %.2f];p_{T} [GeV];#varepsilon_{f}",
                          eps_lo->GetYaxis()->GetBinLowEdge(iy),
                          eps_lo->GetYaxis()->GetBinUpEdge(iy)));
        h1->SetMarkerStyle(20); h1->SetMarkerColor(kRed);   h1->SetLineColor(kRed);
        h2->SetMarkerStyle(21); h2->SetMarkerColor(kBlue);  h2->SetLineColor(kBlue);
        h1->GetYaxis()->SetRangeUser(0., 0.6);
        h1->Draw("E1");
        h2->Draw("E1 SAME");
        TLegend* l = new TLegend(0.55, 0.75, 0.88, 0.88);
        l->AddEntry(h1, Form("MET < %.0f GeV", met_low_max  * 1e-3), "lp");
        l->AddEntry(h2, Form("MET > %.0f GeV", met_high_min * 1e-3), "lp");
        l->SetBorderSize(0); l->Draw();
    }
    c.SaveAs("outputs/met_split_test.pdf");

    TFile* fo = TFile::Open("outputs/met_split_test.root", "RECREATE");
    eps_lo->Write(); eps_hi->Write();
    fo->Close();
    printf("\nWrote outputs/met_split_test.root\nWrote outputs/met_split_test.pdf\n");
}
