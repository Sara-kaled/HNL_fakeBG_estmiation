// =============================================================================
// studies/05_postfit/plot_postfit.C
// PURPOSE: Post-fit stacked distributions — reads pre-fit histograms from
//          plots/04_prefit/prefit_histograms.root.
//          NOTE: ttbar / ztautau normalization factors removed (set to 1.0).
//                Re-enable by setting mu_tt / mu_ztautau if a future fit
//                produces them.
//          Shows statistical uncertainty only (no systematics).
// INPUTS:  plots/04_prefit/prefit_histograms.root
// OUTPUTS: outputs/07_postfit/
// =============================================================================
#include "PlottingUtils.h"

void plot_postfit()
{
    setStyle();
    makeDir("outputs/07_postfit");

    // Normalization factors: removed (kept at 1.0). The TRExFitter post-fit
    // will provide the proper normalization once it is wired up.
    const double mu_tt     = 1.0;
    const double mu_ztautau= 1.0;

    TFile* fin = TFile::Open("outputs/06_prefit/prefit_histograms.root");
    if (!fin || fin->IsZombie()) {
        printf("ERROR: prefit_histograms.root not found — run plot_prefit.C first\n");
        return;
    }

    std::vector<std::string> regNames   = {"SR","VR","CR_tt","CR_ztautau"};
    std::vector<const char*> regLabels  = {"Signal Region","Validation Region",
                                           "CR (t#bar{t})","CR (Z#rightarrow#tau#tau)"};
    std::vector<const char*> xTitles    = {"m_{ll} [GeV]","m_{jj} [GeV]",
                                           "m_{ll} [GeV]","m_{ll} [GeV]"};
    std::vector<bool>        blind       = {true, false, false, false};

    std::vector<std::string> procNames  = {"ttbar","ztautau","wjets","singletop","diboson","higgs"};
    auto meta = sampleMeta();

    for (int r = 0; r < (int)regNames.size(); ++r) {
        const std::string& rn = regNames[r];

        // Load histograms
        std::map<std::string, TH1D*> hMC;
        for (const auto& p : procNames) {
            TH1D* h = (TH1D*)fin->Get(TString::Format("%s_%s", p.c_str(), rn.c_str()));
            if (!h) continue;
            hMC[p] = (TH1D*)h->Clone(TString::Format("pf_%s_%s", p.c_str(), rn.c_str()));
            hMC[p]->SetDirectory(0);
        }
        TH1D* hFake = (TH1D*)fin->Get(TString::Format("fake_%s", rn.c_str()));
        TH1D* hDataRaw = blind[r] ? nullptr :
            (TH1D*)fin->Get(TString::Format("data_%s", rn.c_str()));

        if (hFake) { hFake = (TH1D*)hFake->Clone(TString::Format("pf_fake_%s", rn.c_str())); hFake->SetDirectory(0); }
        if (hDataRaw) { hDataRaw = (TH1D*)hDataRaw->Clone(TString::Format("pf_data_%s", rn.c_str())); hDataRaw->SetDirectory(0); }

        // Apply normalization factors
        if (hMC.count("ttbar")   && hMC["ttbar"])    hMC["ttbar"]->Scale(mu_tt);
        if (hMC.count("ztautau") && hMC["ztautau"])  hMC["ztautau"]->Scale(mu_ztautau);

        // Rebuild total (stat uncertainty only — already in bin errors)
        int nBins = hFake ? hFake->GetNbinsX() : hMC.begin()->second->GetNbinsX();
        double xLo = hFake ? hFake->GetXaxis()->GetXmin() : hMC.begin()->second->GetXaxis()->GetXmin();
        double xHi = hFake ? hFake->GetXaxis()->GetXmax() : hMC.begin()->second->GetXaxis()->GetXmax();
        TH1D* hTotal = new TH1D(TString::Format("pf_total_%s",rn.c_str()),"",nBins,xLo,xHi);
        hTotal->Sumw2();
        for (const auto& p:procNames) if(hMC.count(p)&&hMC[p]) hTotal->Add(hMC[p]);
        if (hFake) hTotal->Add(hFake);

        // Sort by yield (largest at bottom)
        std::vector<std::pair<double,std::string>> yields;
        for (const auto& p:procNames)
            if (hMC.count(p)&&hMC[p]) yields.push_back({hMC[p]->Integral(),p});
        std::sort(yields.begin(),yields.end(),[](auto&a,auto&b){return a.first>b.first;});
        std::vector<std::string> order;
        for (auto& y:yields) order.push_back(y.second);

        bool hasData = !blind[r] && hDataRaw && hDataRaw->Integral()>0;

        TCanvas* c = new TCanvas(TString::Format("cpf_%s",rn.c_str()),"",700,hasData?700:500);
        TPad *upper=nullptr, *lower=nullptr;
        if (hasData) {
            auto [u,l] = makeRatioCanvas(c);
            upper=u; lower=l; upper->cd();
        }

        TLegend* leg = new TLegend(0.55, hasData?0.48:0.55, 0.92, 0.92);
        leg->SetTextSize(0.034);

        TH1D* hUnc = (TH1D*)hTotal->Clone("hUncPF");

        drawStack(order, hMC, hFake,
                  hasData ? hDataRaw : nullptr,
                  hUnc, leg, xTitles[r]);

        TLatex lat; lat.SetNDC(); lat.SetTextFont(72); lat.SetTextSize(0.045);
        lat.DrawLatex(0.15,0.88,"ATLAS"); lat.SetTextFont(42);
        lat.DrawLatex(0.27,0.88,"Internal"); lat.SetTextSize(0.038);
        lat.DrawLatex(0.15,0.82,regLabels[r]);
        lat.SetTextSize(0.030);
        lat.DrawLatex(0.15,0.76,TString::Format(
            "#mu_{t#bar{t}}=%.2f,  #mu_{Z#tau#tau}=%.2f, stat. unc. only",
            mu_tt, mu_ztautau).Data());

        if (hasData) drawRatio(lower, hDataRaw, hTotal, xTitles[r]);

        TString outPath = TString::Format("outputs/07_postfit/postfit_%s.pdf", rn.c_str());
        c->SaveAs(outPath);
        printf("Saved %s\n", outPath.Data());
        delete c;
    }

    fin->Close();
    printf("\n=== plot_postfit done. Output in outputs/07_postfit/ ===\n");
    printf("  Applied: mu_tt=%.2f, mu_ztautau=%.2f, stat unc only\n", mu_tt, mu_ztautau);
}
