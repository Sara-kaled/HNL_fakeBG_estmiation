// =============================================================================
// plotting/plot_prefit.C
// PURPOSE: Pre-fit stacked distributions in SR, VR, CR_tt, CR_ztautau.
//          MC with nominal weights (no fit normalization).
//          Fake contribution from AMM weights.
//          Data shown in VR, CR_tt, CR_ztautau. SR is blinded.
//          Stacking order: largest background at bottom (sorted by yield).
// INPUTS:  MC ntuples on EOS, AMM/data_with_electron_fake_weights.root
// OUTPUTS: outputs/06_prefit/
// =============================================================================
#include "PlottingUtils.h"
#include "../utils/SelectionUtils.h"
#include <TChain.h>

// ----- Region definitions (custom prefit regions; closely related to
//       SelectionUtils SR/VR/ttCR/ZttCR but with prefit-specific tweaks
//       on the obs binning and the ztautau mll window 30-80 GeV) ------------
struct Region {
    const char* name;
    const char* label;
    const char* obs;       // "mll" or "mjj"
    const char* xTitle;
    int nBins; double xLo, xHi;
    bool blind;
};

// Pass region using EventData fields and the el_is_tight / loose flags
inline bool passPrefitRegion(const EventData& ev, const Region& reg, bool el_id_pass)
{
    if (!el_id_pass) return false;
    bool vbf = passVBFJets(ev);
    double mll_g = ev.mll * 1e-3, mjj_g = ev.mjj * 1e-3, l1_g = ev.l1_pt * 1e-3;
    bool mjj_hi = (mjj_g > 400.);
    bool mjj_lo = (mjj_g > 0. && mjj_g < 400.);
    TString rn(reg.name);
    if (rn == "SR")
        return ev.is_os && ev.no_bjets  && vbf && mjj_hi;
    if (rn == "VR")
        return ev.is_os && ev.no_bjets  && vbf && mjj_lo;
    if (rn == "CR_tt")
        return ev.is_os && !ev.no_bjets && vbf && mjj_hi && l1_g > 45.;
    if (rn == "CR_ztautau")
        return ev.is_os && ev.no_bjets  && vbf && mjj_hi
            && mll_g > 30. && mll_g < 80. && l1_g > 35. && l1_g < 45.
            && ev.l2_pt > 15000.f;
    return false;
}

void plot_prefit()
{
    setStyle();
    makeDir("outputs/06_prefit");

    std::vector<Region> regions = {
        {"SR",         "Signal Region",             "mll", "m_{ll} [GeV]", 30, 0, 150, true },
        {"VR",         "Validation Region",         "mjj", "m_{jj} [GeV]", 20, 0, 400, false},
        {"CR_tt",      "CR (t#bar{t})",             "mll", "m_{ll} [GeV]", 30, 0, 150, false},
        {"CR_ztautau", "CR (Z#rightarrow#tau#tau)", "mll", "m_{ll} [GeV]", 6,  0, 150, false},
    };
    const int NR = (int)regions.size();
    std::vector<std::string> procNames = {"ttbar","ztautau","wjets","singletop","diboson","higgs"};
    auto meta = sampleMeta();
    std::string base = "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/";

    // Histogram allocation
    std::map<std::string,std::vector<TH1D*>> hMC;
    std::vector<TH1D*> hFake(NR), hData(NR), hTotal(NR);
    for (const auto& p : procNames) {
        hMC[p].resize(NR);
        for (int r=0;r<NR;++r) {
            hMC[p][r] = new TH1D(TString::Format("hMC_%s_%s",p.c_str(),regions[r].name),
                                 "", regions[r].nBins,regions[r].xLo,regions[r].xHi);
            hMC[p][r]->Sumw2();
        }
    }
    for (int r=0;r<NR;++r) {
        hFake[r]  = new TH1D(TString::Format("hFake_%s", regions[r].name),"",
                             regions[r].nBins,regions[r].xLo,regions[r].xHi);
        hData[r]  = new TH1D(TString::Format("hData_%s", regions[r].name),"",
                             regions[r].nBins,regions[r].xLo,regions[r].xHi);
        hTotal[r] = new TH1D(TString::Format("hTotal_%s",regions[r].name),"",
                             regions[r].nBins,regions[r].xLo,regions[r].xHi);
        hFake[r]->Sumw2(); hData[r]->Sumw2();
    }

    // ----- MC loop -----
    for (const auto& pname : procNames) {
        printf("Processing MC: %s\n", pname.c_str());
        TChain chain("tree");
        chain.Add((base + pname + "/*.root").c_str());
        EventData ev; ResetTrigMatchPointers(ev);
        SetupBranches(&chain, ev);

        Long64_t nE = chain.GetEntries();
        for (Long64_t i=0;i<nE;++i) {
            chain.GetEntry(i); ComputeDerived(ev);
            if (!passBasePresel(ev)) continue;
            if (!isMCPromptEvent(ev)) continue;          // prefit MC = prompt-truth-matched
            double w = mcWeightTight(ev);
            if (!std::isfinite(w)) continue;
            double mll_g = ev.mll*1e-3, mjj_g = ev.mjj*1e-3;
            for (int r=0;r<NR;++r) {
                double val = (TString(regions[r].obs)=="mll") ? mll_g : mjj_g;
                if (passPrefitRegion(ev, regions[r], ev.el_is_tight))
                    hMC[pname][r]->Fill(val, w);
            }
        }
    }

    // ----- Data + AMM fake-weight loop -----
    printf("Processing data (fake weights + observed)\n");
    TChain dTree("tree");
    dTree.Add("../AMM/data_with_electron_fake_weights.root");
    EventData ev; ResetTrigMatchPointers(ev);
    SetupBranches(&dTree, ev);

    Long64_t nD = dTree.GetEntries();
    for (Long64_t i=0;i<nD;++i) {
        dTree.GetEntry(i); ComputeDerived(ev);
        if (!passBasePresel(ev)) continue;
        double mll_g = ev.mll*1e-3, mjj_g = ev.mjj*1e-3;

        // Fake contribution: loose baseline + AMM weight
        if (ev.el_is_loose && std::isfinite(ev.fw_nom)) {
            for (int r=0;r<NR;++r) {
                double val=(TString(regions[r].obs)=="mll")?mll_g:mjj_g;
                if (passPrefitRegion(ev, regions[r], ev.el_is_loose))
                    hFake[r]->Fill(val, ev.fw_nom);
            }
        }
        // Observed data: tight electron, w=1
        if (ev.el_is_tight) {
            for (int r=0;r<NR;++r) {
                if (regions[r].blind) continue;
                double val=(TString(regions[r].obs)=="mll")?mll_g:mjj_g;
                if (passPrefitRegion(ev, regions[r], ev.el_is_tight))
                    hData[r]->Fill(val, 1.0);
            }
        }
    }

    // ----- Build total + plot -----
    TFile* fout = TFile::Open("outputs/06_prefit/prefit_histograms.root","RECREATE");

    for (int r=0;r<NR;++r) {
        for (const auto& p : procNames) hTotal[r]->Add(hMC[p][r]);
        hTotal[r]->Add(hFake[r]);

        std::vector<std::pair<double,std::string>> yields;
        for (const auto& p:procNames) yields.push_back({hMC[p][r]->Integral(),p});
        std::sort(yields.begin(),yields.end(),
                  [](auto&a,auto&b){return a.first>b.first;});
        std::vector<std::string> order;
        for (auto& y:yields) order.push_back(y.second);

        bool hasData = !regions[r].blind && hData[r]->Integral()>0;
        TCanvas* c = new TCanvas(TString::Format("c_%s",regions[r].name),"",
                                 700, hasData?700:500);
        TPad *upper=nullptr,*lower=nullptr;
        if (hasData) { auto [u,l]=makeRatioCanvas(c); upper=u; lower=l; upper->cd(); }

        TLegend* leg = new TLegend(0.55, hasData?0.48:0.55, 0.92, 0.92);
        leg->SetTextSize(0.034);
        std::map<std::string,TH1D*> hMCr;
        for (const auto& p:procNames) hMCr[p]=hMC[p][r];

        TH1D* hUnc = (TH1D*)hTotal[r]->Clone("hUnc");
        drawStack(order, hMCr, hFake[r],
                  hasData ? hData[r] : nullptr,
                  hUnc, leg, regions[r].xTitle);

        TLatex lat; lat.SetNDC(); lat.SetTextFont(72); lat.SetTextSize(0.045);
        lat.DrawLatex(0.15,0.88,"ATLAS");
        lat.SetTextFont(42); lat.DrawLatex(0.27,0.88,"Internal");
        lat.SetTextSize(0.038);
        lat.DrawLatex(0.15,0.82,regions[r].label);
        if (regions[r].blind) {
            lat.SetTextColor(kRed); lat.DrawLatex(0.15,0.76,"Blinded");
            lat.SetTextColor(kBlack);
        }
        if (hasData) drawRatio(lower, hData[r], hTotal[r], regions[r].xTitle);

        TString outPath = TString::Format("outputs/06_prefit/prefit_%s.pdf",
                                          regions[r].name);
        c->SaveAs(outPath); printf("Saved %s\n", outPath.Data());

        fout->cd();
        hTotal[r]->Write(TString::Format("total_%s",regions[r].name));
        hFake[r]->Write(TString::Format("fake_%s",regions[r].name));
        for (const auto& p:procNames)
            hMC[p][r]->Write(TString::Format("%s_%s",p.c_str(),regions[r].name));
        if (!regions[r].blind)
            hData[r]->Write(TString::Format("data_%s",regions[r].name));
        delete c;
    }

    fout->Close();
    printf("\n=== plot_prefit done. Output in outputs/06_prefit/ ===\n");
}
