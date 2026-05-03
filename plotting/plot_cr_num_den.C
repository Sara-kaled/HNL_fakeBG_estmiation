// =============================================================================
// studies/03_cr_distributions/plot_cr_num_den.C
// PURPOSE: Stacked distributions of loose (denominator) and tight (numerator)
//          electrons in the SS fake CR, with data overlay.
//          Validates that data > MC in the SS CR and shows the fake-enriched
//          sample used to measure the fake rate.
// INPUTS:  MC_eff/MC_histograms.root, data_eff/Data_histograms.root
// OUTPUTS: outputs/02_1d_rates/cr_num_den/
// =============================================================================
#include "PlottingUtils.h"

void plotNumDen(TH2* hData, TH2* hMC,
                bool useTight,
                const char* label, const char* outPath)
{
    // Project onto pT axis (X axis of the 2D histogram)
    TH1D* hD = (TH1D*)hData->ProjectionX(
        TString::Format("hD_%s_%s", useTight?"tight":"loose", label).Data());
    TH1D* hM = (TH1D*)hMC->ProjectionX(
        TString::Format("hM_%s_%s", useTight?"tight":"loose", label).Data());

    TCanvas c("c","c",800,700);
    auto [upper, lower] = makeRatioCanvas(&c);
    upper->cd();

    double ymax = std::max(hD->GetMaximum(), hM->GetMaximum()) * 1.4;
    hM->SetFillColor(kAzure-4); hM->SetLineColor(kBlack);
    hM->GetYaxis()->SetTitle("Events");
    hM->GetYaxis()->SetRangeUser(0, ymax);
    hM->GetXaxis()->SetLabelSize(0);
    hM->Draw("HIST");

    hD->SetMarkerStyle(20); hD->SetMarkerSize(0.9); hD->SetLineColor(kBlack);
    hD->Draw("E SAME");

    TLegend leg(0.55,0.70,0.92,0.88);
    leg.SetTextSize(0.038);
    leg.AddEntry(hD,"Data","pe");
    leg.AddEntry(hM,"MC (prompt)","f");
    leg.Draw();

    TLatex lat; lat.SetNDC(); lat.SetTextFont(72); lat.SetTextSize(0.045);
    lat.DrawLatex(0.15,0.88,"ATLAS");
    lat.SetTextFont(42); lat.DrawLatex(0.27,0.88,"Internal");
    lat.SetTextSize(0.038);
    lat.DrawLatex(0.15,0.82, TString::Format("SS CR — %s electrons",
        useTight?"tight (numerator)":"loose (denominator)").Data());

    lower->cd();
    TH1D* ratio = (TH1D*)hD->Clone("ratio");
    ratio->Divide(hM);
    ratio->SetMarkerStyle(20); ratio->SetMarkerSize(0.8);
    ratio->GetYaxis()->SetTitle("Data / MC");
    ratio->GetYaxis()->SetNdivisions(505);
    ratio->GetYaxis()->SetTitleSize(0.13); ratio->GetYaxis()->SetLabelSize(0.11);
    ratio->GetYaxis()->SetTitleOffset(0.38);
    ratio->GetXaxis()->SetTitle("Electron p_{T} [MeV]");
    ratio->GetXaxis()->SetTitleSize(0.13); ratio->GetXaxis()->SetLabelSize(0.11);
    ratio->GetYaxis()->SetRangeUser(0., 3.0);
    ratio->Draw("E");
    TLine* line = new TLine(ratio->GetXaxis()->GetXmin(),1,
                             ratio->GetXaxis()->GetXmax(),1);
    line->SetLineStyle(2); line->SetLineColor(kRed); line->Draw();

    c.SaveAs(outPath);
    printf("Saved %s\n", outPath);
}

void plot_cr_num_den()
{
    setStyle();
    makeDir("outputs/02_1d_rates/cr_num_den");

    TFile* fData = TFile::Open("../efficiency/outputs/Data_histograms.root");
    TFile* fMC   = TFile::Open("../efficiency/outputs/MC_histograms.root");
    if (!fData || fData->IsZombie() || !fMC || fMC->IsZombie()) {
        printf("ERROR: histogram files not found\n"); return;
    }

    TH2F* hDL = (TH2F*)fData->Get("h2_loose_eta_pt_CR");
    TH2F* hDT = (TH2F*)fData->Get("h2_tight_eta_pt_CR");
    TH2F* hML = (TH2F*)fMC->Get("h2_loose_eta_pt_CR");
    TH2F* hMT = (TH2F*)fMC->Get("h2_tight_eta_pt_CR");

    if (!hDL||!hDT||!hML||!hMT) {
        printf("ERROR: histograms not found\n"); fData->ls(); return;
    }

    hDL->SetDirectory(0); hDT->SetDirectory(0);
    hML->SetDirectory(0); hMT->SetDirectory(0);
    fData->Close(); fMC->Close();

    // Denominator (loose)
    plotNumDen(hDL, hML, false,
               "loose",
               "outputs/02_1d_rates/cr_num_den/ss_cr_loose_denominator.pdf");

    // Numerator (tight)
    plotNumDen(hDT, hMT, true,
               "tight",
               "outputs/02_1d_rates/cr_num_den/ss_cr_tight_numerator.pdf");

    // -------------------------------------------------------------------------
    // Data minus MC (fake rate numerator = what we actually measure)
    // -------------------------------------------------------------------------
    TCanvas c2("c2","c2",700,500);
    c2.SetLeftMargin(0.13); c2.SetBottomMargin(0.13);

    TH1D* hDL1 = (TH1D*)hDL->ProjectionX("hDL1");
    TH1D* hML1 = (TH1D*)hML->ProjectionX("hML1");
    TH1D* hSub = (TH1D*)hDL1->Clone("hSub");
    hSub->Add(hML1, -1.);
    hSub->SetLineColor(kRed+1); hSub->SetLineWidth(2);
    hSub->GetXaxis()->SetTitle("Electron p_{T} [MeV]");
    hSub->GetYaxis()->SetTitle("Data - MC (fake numerator)");
    hSub->Draw("HIST E");

    TLatex lat; lat.SetNDC(); lat.SetTextFont(72); lat.SetTextSize(0.045);
    lat.DrawLatex(0.15,0.88,"ATLAS");
    lat.SetTextFont(42); lat.DrawLatex(0.27,0.88,"Internal");
    lat.SetTextSize(0.038);
    lat.DrawLatex(0.15,0.82,"SS CR — data #minus prompt MC");

    c2.SaveAs("outputs/02_1d_rates/cr_num_den/ss_cr_data_minus_mc.pdf");

    printf("\n=== plot_cr_num_den done. Output in outputs/02_1d_rates/cr_num_den/ ===\n");
}
