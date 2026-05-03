// =============================================================================
// FILE:    studies/dR_jet_distribution.C
// PURPOSE: ΔR(electron, nearest jet) distribution for loose fake electrons in
//          the SS μe fake CR, per pT bin and split by truth source.
//          - Real jet→e fakes cluster at small ΔR (< 0.4): the electron IS
//            embedded in the parent jet.
//          - Heavy-flavor real e's from leptonic decays of b/c hadrons sit
//            at small-to-medium ΔR (electron close to but not on the jet axis).
//          - Real prompt electrons sit at large ΔR (well-separated).
//
//          A high-pT loose population peaking at large ΔR strongly suggests
//          those "fakes" are not really fakes but real isolated electrons.
//
// OUTPUT:  outputs/dR_jet_distribution.root
//          outputs/dR_jet_distribution.pdf
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TMath.h>
#include "ROOT/RVec.hxx"
#include <iostream>
#include <cstdio>
#include <cmath>

#include "../../utils/SelectionUtils.h"

using ROOT::VecOps::RVec;

enum Cat { kPrompt=0, kHF=1, kLF=2, kConv=3, kOther=4, kNCat=5 };
static const char* kCatName[kNCat] = {"prompt","HF (b/c)","LF","Conv","Other"};
static const int   kCatColor[kNCat] = {kBlack, kRed+1, kGreen+2, kBlue+1, kGray+1};

static int classify(int iff, bool isPrompt) {
    if (isPrompt) return kPrompt;
    if (iff == 8 || iff == 9) return kHF;
    if (iff == 10) return kLF;
    if (iff == 5)  return kConv;
    return kOther;
}

void dR_jet_distribution()
{
    gStyle->SetOptStat(0);
    gSystem->mkdir("outputs", true);

    TChain mc("tree");
    const std::string base = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    for (auto p : {"ttbar","ztautau","wjets","singletop","diboson","higgs"})
        mc.Add((base + p + "/*.root").c_str());

    EventData ev;
    SetupBranches(&mc, ev);
    RVecI* el_IFFClass = nullptr;
    RVecF* jet_eta_v   = nullptr;
    RVecF* jet_phi_v   = nullptr;
    mc.SetBranchAddress("el_IFFClass", &el_IFFClass);
    mc.SetBranchAddress("jet_eta",     &jet_eta_v);
    mc.SetBranchAddress("jet_phi",     &jet_phi_v);

    std::vector<std::vector<TH1D*>> h(kNPtBins, std::vector<TH1D*>(kNCat, nullptr));
    for (int ix = 0; ix < kNPtBins; ++ix) {
        for (int s = 0; s < kNCat; ++s) {
            h[ix][s] = new TH1D(Form("h_dR_pt%d_%s", ix, kCatName[s]),
                                Form(";#DeltaR(loose e, nearest jet) (p_{T} #in [%.0f,%.0f]);electrons",
                                     kPtBins[ix], kPtBins[ix+1]),
                                40, 0., 4.);
            h[ix][s]->Sumw2();
            h[ix][s]->SetLineColor(kCatColor[s]);
            h[ix][s]->SetLineWidth(2);
        }
    }

    auto deltaR = [](double e1, double p1, double e2, double p2) {
        double dphi = TMath::Abs(p1 - p2);
        if (dphi > M_PI) dphi = 2. * M_PI - dphi;
        double deta = e1 - e2;
        return std::sqrt(deta*deta + dphi*dphi);
    };

    const Long64_t n = mc.GetEntries();
    std::cout << "MC entries: " << n << "\n";
    for (Long64_t i = 0; i < n; ++i) {
        if (i % 1000000 == 0) std::cout << "  " << i << " / " << n << "\n";
        mc.GetEntry(i);
        ComputeDerived(ev);
        if (!passBasePresel_FakeMeasurement(ev)) continue;
        if (!isFakeCR(ev)) continue;
        if (!ev.no_bjets) continue;
        if (!ev.el_is_loose) continue;
        if (!el_IFFClass || el_IFFClass->empty()) continue;
        if (!ev.el_eta_vec || ev.el_eta_vec->empty()) continue;
        if (!jet_eta_v || jet_eta_v->empty()) continue;

        const double pt = ev.l2_pt * 1e-3;
        int ix = -1;
        for (int b = 0; b < kNPtBins; ++b)
            if (pt >= kPtBins[b] && pt < kPtBins[b+1]) { ix = b; break; }
        if (ix < 0) continue;

        // ΔR to nearest jet
        double el_eta = (*ev.el_eta_vec)[0];
        double el_phi = (ev.el_phi_vec && !ev.el_phi_vec->empty())
                        ? (*ev.el_phi_vec)[0] : 0.;
        double dRmin = 999.;
        for (size_t j = 0; j < jet_eta_v->size(); ++j) {
            double dR = deltaR(el_eta, el_phi, (*jet_eta_v)[j], (*jet_phi_v)[j]);
            if (dR < dRmin) dRmin = dR;
        }

        const int cat = classify((*el_IFFClass)[0], isMCPromptElectron(ev));
        h[ix][cat]->Fill(dRmin, mcWeightLoose(ev));
    }

    printf("\n========================================================================\n");
    printf("  Mean ΔR(loose e, nearest jet) per pT bin × truth category\n");
    printf("========================================================================\n");
    printf("  pT [GeV]      ");
    for (int s = 0; s < kNCat; ++s) printf("%9s   ", kCatName[s]);
    printf("\n  -------------------------------------------------------------------\n");
    for (int ix = 0; ix < kNPtBins; ++ix) {
        printf("  [%4.0f,%4.0f]  ", kPtBins[ix], kPtBins[ix+1]);
        for (int s = 0; s < kNCat; ++s)
            printf("  %.2f      ", h[ix][s]->GetMean());
        printf("\n");
    }
    printf("\n  Interpretation:\n");
    printf("    LF jet→e and HF should peak at small ΔR (< 0.4) — embedded in jet.\n");
    printf("    Prompt + Conv sit at larger ΔR (> 1) — separated from jets.\n");
    printf("    If high-pT 'fakes' have <ΔR> > 1, they're not real fakes.\n");

    TCanvas c("c","ΔR(e,jet) split", 1500, 1000);
    c.Divide(3, (kNPtBins+2)/3);
    for (int ix = 0; ix < kNPtBins; ++ix) {
        c.cd(ix+1); gPad->SetLogy();
        double maxY = 0.;
        for (int s = 0; s < kNCat; ++s) maxY = std::max(maxY, h[ix][s]->GetMaximum());
        for (int s = 0; s < kNCat; ++s) {
            h[ix][s]->GetYaxis()->SetRangeUser(0.5, 5.0 * maxY);
            h[ix][s]->Draw(s == 0 ? "HIST" : "HIST SAME");
        }
        if (ix == 0) {
            auto* leg = new TLegend(0.55, 0.55, 0.88, 0.88);
            for (int s = 0; s < kNCat; ++s) leg->AddEntry(h[ix][s], kCatName[s], "l");
            leg->SetBorderSize(0); leg->Draw();
        }
    }
    c.SaveAs("outputs/dR_jet_distribution.pdf");

    TFile* fo = TFile::Open("outputs/dR_jet_distribution.root", "RECREATE");
    for (int ix = 0; ix < kNPtBins; ++ix)
        for (int s = 0; s < kNCat; ++s) h[ix][s]->Write();
    fo->Close();
    printf("\nWrote outputs/dR_jet_distribution.{root,pdf}\n");
}
