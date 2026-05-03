// =============================================================================
// FILE:    studies/checks/fake_source_diagnostics.C
// PURPOSE: One-stop diagnostic of fake-electron truth-source breakdown vs.
//          relevant variables, in the SS μe fake CR. Replaces the older
//          truth_source_breakdown.C (kept its file but this one is more
//          comprehensive). Produces fake-source distributions (B / C / Conv /
//          LF / Tau / Other) as a function of:
//             pT                MET                |Δη_jj|
//             |η|               |d0/σ|             m_{ll}
//             m_{jj}            ΔR(e, j_lead)
//
//          All plots in ATLAS Internal style, full Run-2 luminosity, and
//          every panel has a legend.
//
// INPUT:   MC ntuples (TChain) — same list as the rest of the framework.
// OUTPUT:  outputs/fake_source_diagnostics.root
//          outputs/fake_source_diagnostics.pdf  (multi-panel)
// =============================================================================
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <THStack.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TSystem.h>
#include <TMath.h>
#include "ROOT/RVec.hxx"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

#include "../../utils/SelectionUtils.h"
#include "../../utils/AtlasStyle.h"

using ROOT::VecOps::RVec;

enum FakeSrc { kB=0, kC=1, kConv=2, kLF=3, kTau=4, kOther=5, kNSrc=6 };
static const char* kSrcName[kNSrc] = {"B","C","Conv","LF","Tau","Other"};

static int classifyIFF(int iff) {
    switch(iff) {
        case 8:  return kB;
        case 9:  return kC;
        case 5:  return kConv;
        case 10: return kLF;
        case 7:  return kTau;
        default: return kOther;
    }
}

struct VarDef {
    std::string name;
    std::string xLabel;
    int    nBins;
    double xLo, xHi;
};

static int colorOf(int s) {
    switch (s) {
        case kB:    return kATLAS_src_B;
        case kC:    return kATLAS_src_C;
        case kConv: return kATLAS_src_Conv;
        case kLF:   return kATLAS_src_LF;
        case kTau:  return kATLAS_src_Tau;
        default:    return kATLAS_src_Other;
    }
}

void fake_source_diagnostics()
{
    SetATLASStyle();
    gSystem->mkdir("outputs", true);

    // ----- MC chain -----
    const std::string base = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";
    TChain mc("tree");
    for (auto p : {"ttbar","ztautau","wjets","singletop","diboson","higgs"})
        mc.Add((base + p + "/*.root").c_str());

    EventData ev;
    SetupBranches(&mc, ev);
    RVecI* el_iff   = nullptr;
    RVecF* jet_eta_v= nullptr;
    RVecF* jet_phi_v= nullptr;
    mc.SetBranchAddress("el_IFFClass", &el_iff);
    mc.SetBranchAddress("jet_eta",     &jet_eta_v);
    mc.SetBranchAddress("jet_phi",     &jet_phi_v);

    // ----- Variable definitions -----
    const std::vector<VarDef> vars = {
        {"pt",     "Electron p_{T} [GeV]",          40,  0., 200.},
        {"eta",    "Electron |#eta|",                25,  0.,  2.5},
        {"d0sig",  "Electron |d_{0}/#sigma|",        40,  0.,  8. },
        {"met",    "E_{T}^{miss} [GeV]",             40,  0., 200.},
        {"mll",    "m_{ll} [GeV]",                   40,  0., 200.},
        {"mjj",    "m_{jj} [GeV]",                   40,  0., 1000.},
        {"detajj", "|#Delta#eta_{jj}|",              25,  0.,  6. },
        {"dRej1",  "#DeltaR(e, leading jet)",        40,  0.,  4. }
    };

    // ----- Histograms: h[varIndex][source] -----
    const int V = (int)vars.size();
    std::vector<std::vector<TH1D*>> h(V, std::vector<TH1D*>(kNSrc, nullptr));
    for (int v = 0; v < V; ++v) {
        for (int s = 0; s < kNSrc; ++s) {
            h[v][s] = new TH1D(Form("h_%s_%s", vars[v].name.c_str(), kSrcName[s]),
                               Form(";%s;Events", vars[v].xLabel.c_str()),
                               vars[v].nBins, vars[v].xLo, vars[v].xHi);
            h[v][s]->Sumw2();
            h[v][s]->SetLineColor (colorOf(s));
            h[v][s]->SetFillColor (colorOf(s));
            h[v][s]->SetMarkerColor(colorOf(s));
            h[v][s]->SetLineWidth(2);
        }
    }

    auto deltaR = [](double e1,double p1,double e2,double p2){
        double dphi = TMath::Abs(p1 - p2);
        if (dphi > M_PI) dphi = 2*M_PI - dphi;
        double deta = e1 - e2;
        return std::sqrt(deta*deta + dphi*dphi);
    };

    // ----- Loop -----
    const Long64_t n = mc.GetEntries();
    std::cout << "MC entries: " << n << "\n";
    long nUsed = 0;
    for (Long64_t i = 0; i < n; ++i) {
        if (i % 1000000 == 0) std::cout << "  " << i << " / " << n << "\n";
        mc.GetEntry(i);
        ComputeDerived(ev);
        if (!passBasePresel(ev))      continue;
        if (!isFakeCR(ev))            continue;
        if (!ev.no_bjets)             continue;
        if (!ev.el_is_loose)          continue;
        if (isMCPromptElectron(ev))   continue;        // truth-fakes only
        if (!el_iff || el_iff->empty()) continue;

        const int s = classifyIFF((*el_iff)[0]);
        const double w = mcWeightLoose(ev);

        const double pt   = ev.l2_pt * 1e-3;
        const double aeta = ev.el_aeta;
        const double d0s  = (ev.el_d0sig && !ev.el_d0sig->empty())
                            ? std::fabs((*ev.el_d0sig)[0]) : 0.;
        const double met  = ev.MET * 1e-3;
        const double mll  = ev.mll * 1e-3;
        const double mjj  = ev.mjj * 1e-3;
        const double dej  = std::fabs(ev.deta_jj);

        // ΔR(e, leading jet) — use first jet in the vector (leading by convention)
        double dRej = -1.;
        if (jet_eta_v && !jet_eta_v->empty() &&
            ev.el_eta_vec && !ev.el_eta_vec->empty()) {
            double el_eta = (*ev.el_eta_vec)[0];
            double el_phi = (ev.el_phi_vec && !ev.el_phi_vec->empty())
                            ? (*ev.el_phi_vec)[0] : 0.;
            dRej = deltaR(el_eta, el_phi, (*jet_eta_v)[0], (*jet_phi_v)[0]);
        }

        h[0][s]->Fill(pt,   w);
        h[1][s]->Fill(aeta, w);
        h[2][s]->Fill(d0s,  w);
        h[3][s]->Fill(met,  w);
        h[4][s]->Fill(mll,  w);
        h[5][s]->Fill(mjj,  w);
        h[6][s]->Fill(dej,  w);
        if (dRej >= 0) h[7][s]->Fill(dRej, w);
        ++nUsed;
    }
    std::cout << "Loose fake electrons used: " << nUsed << "\n";

    // ----- Plot: 4×2 grid of stacks, each panel with its own legend -----
    TCanvas c("c","Fake-source diagnostics", 1600, 900);
    c.Divide(4, 2);

    for (int v = 0; v < V; ++v) {
        c.cd(v+1);
        gPad->SetLogy();

        THStack* hs = new THStack(Form("hs_%s", vars[v].name.c_str()),
                                  Form(";%s;Events", vars[v].xLabel.c_str()));
        for (int s = 0; s < kNSrc; ++s) hs->Add(h[v][s]);
        hs->Draw("HIST");
        hs->GetYaxis()->SetTitleOffset(1.4);

        TLegend* leg = new TLegend(0.65, 0.55, 0.95, 0.92);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.040);
        for (int s = 0; s < kNSrc; ++s)
            leg->AddEntry(h[v][s], kSrcName[s], "f");
        leg->Draw();

        if (v == 0) DrawAtlasLabel(0.18, 0.88);
        TLatex lat; lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.04);
        lat.DrawLatex(0.18, 0.80, "Fake CR (SS #mue), truth fakes");
    }

    c.SaveAs("outputs/fake_source_diagnostics.pdf");

    // ----- Persist -----
    TFile* fo = TFile::Open("outputs/fake_source_diagnostics.root", "RECREATE");
    for (int v = 0; v < V; ++v)
        for (int s = 0; s < kNSrc; ++s) h[v][s]->Write();
    fo->Close();

    // ----- Per-pT integrated fractions table -----
    printf("\n========================================================================\n");
    printf("  Fake-source fractions vs pT (SS μe fake CR, truth-matched MC)\n");
    printf("========================================================================\n");
    printf("  pT bin       ");
    for (int s = 0; s < kNSrc; ++s) printf("  %-7s", kSrcName[s]);
    printf("    total\n");
    printf("  ----------------------------------------------------------------\n");
    auto* hPt = h[0][kB];   // any source, used only for the X axis
    for (int ib = 1; ib <= hPt->GetNbinsX(); ib += 4) {   // 4-bin grouping
        double tot = 0.;
        for (int s = 0; s < kNSrc; ++s)
            for (int j = 0; j < 4 && ib+j <= hPt->GetNbinsX(); ++j)
                tot += h[0][s]->GetBinContent(ib+j);
        if (tot < 1) continue;
        double lo = hPt->GetXaxis()->GetBinLowEdge(ib);
        double hi = hPt->GetXaxis()->GetBinUpEdge(std::min(ib+3, hPt->GetNbinsX()));
        printf("  [%5.1f,%5.1f]  ", lo, hi);
        for (int s = 0; s < kNSrc; ++s) {
            double f = 0;
            for (int j = 0; j < 4 && ib+j <= hPt->GetNbinsX(); ++j)
                f += h[0][s]->GetBinContent(ib+j);
            printf("  %5.1f%%", 100.*f/tot);
        }
        printf("    %.1f\n", tot);
    }

    std::cout << "\nWrote outputs/fake_source_diagnostics.{root,pdf}\n";
}
