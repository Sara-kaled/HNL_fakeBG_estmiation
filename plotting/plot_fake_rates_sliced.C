#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <iostream>
#include <string>

using namespace std;

void SetATLASStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetPadTopMargin(0.07);
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetPadBottomMargin(0.16);
    gStyle->SetPadLeftMargin(0.16);
}

void DrawATLASLabels(TPad* pad, const string& extra = "") {
    pad->cd();
    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(72); latex.SetTextSize(0.055);
    latex.DrawLatex(0.20, 0.88, "ATLAS");
    latex.SetTextFont(42);
    latex.DrawLatex(0.34, 0.88, "Internal");
    latex.SetTextSize(0.045);
    latex.DrawLatex(0.20, 0.81, "#sqrt{s} = 13 TeV, 140 fb^{-1}");
    latex.DrawLatex(0.20, 0.74, ("#mu e" + extra).c_str());
}

void plot_fake_rates_sliced() {
    SetATLASStyle();

    TFile* f = TFile::Open("../efficiency/outputs/fake_rates.root", "READ");
    if (!f || f->IsZombie()) {
        cout << "ERROR: Cannot open file" << endl;
        return;
    }

    TH2F* h2d = (TH2F*)f->Get("fake_rate_2d_CR;1");

    TGraphAsymmErrors* g_pt  = (TGraphAsymmErrors*)f->Get("fake_rate_1d_pt_CR;1");
    TGraphAsymmErrors* g_eta = (TGraphAsymmErrors*)f->Get("fake_rate_1d_eta_CR;1");
    TGraphAsymmErrors* g_d0  = (TGraphAsymmErrors*)f->Get("fake_rate_1d_d0_CR;1");
    TGraphAsymmErrors* g_met = (TGraphAsymmErrors*)f->Get("fake_rate_1d_met_CR;1");

    system("mkdir -p outputs/02_1d_rates/sliced");

    // === 1. Individual 1D plots (most reliable) ===
    auto plotSingle = [&](TGraphAsymmErrors* gr, const string& name, const string& xTitle) {
        if (!gr) return;
        TCanvas* c = new TCanvas(("c_" + name).c_str(), "", 700, 600);
        c->cd();

        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.0);
        gr->SetMarkerColor(kBlue+1);
        gr->SetLineColor(kBlue+1);
        gr->SetLineWidth(2);

        gr->Draw("AP");
        gr->GetXaxis()->SetTitle(xTitle.c_str());
        gr->GetYaxis()->SetTitle("Fake Rate");
        gr->GetYaxis()->SetRangeUser(0.0, 0.5);   // <-- Changed to 0-0.5, common for fake rates

        DrawATLASLabels((TPad*)gPad);

        c->SaveAs(("outputs/02_1d_rates/sliced/fake_rate_" + name + ".pdf").c_str());
    };

    plotSingle(g_pt,  "pt",  "Electron p_{T} [GeV]");
    plotSingle(g_eta, "eta", "Electron |#eta|");
    plotSingle(g_d0,  "d0",  "Electron d_{0} [mm]");
    plotSingle(g_met, "met", "MET [GeV]");

    // === 2. Sliced from 2D (with safety) ===
    if (h2d) {
        cout << "2D histogram found. Making sliced plots..." << endl;

        // vs |eta| sliced in pT
        TCanvas* c1 = new TCanvas("c_eta_sliced", "", 700, 600);
        c1->cd();

        vector<pair<double,double>> pt_cuts = {{0,20}, {20,40}, {40,200}};
        vector<int> cols = {kRed, kBlue, kGreen+2};
        vector<string> labs = {"p_{T} < 20 GeV", "20 < p_{T} < 40 GeV", "p_{T} > 40 GeV"};

        TLegend* leg = new TLegend(0.55, 0.70, 0.92, 0.92);

        for (size_t i = 0; i < pt_cuts.size(); ++i) {
            int bx1 = h2d->GetXaxis()->FindBin(pt_cuts[i].first);
            int bx2 = h2d->GetXaxis()->FindBin(pt_cuts[i].second) - 1;

            TH1F* h = (TH1F*)h2d->ProjectionY(Form("py_%d",i), bx1, bx2);
            if (h->GetEntries() == 0) continue;

            h->SetMarkerStyle(20+i);
            h->SetMarkerColor(cols[i]);
            h->SetLineColor(cols[i]);

            if (i == 0) {
                h->Draw("PE");
                h->GetXaxis()->SetTitle("Electron |#eta|");
                h->GetYaxis()->SetTitle("Fake Rate");
                h->GetYaxis()->SetRangeUser(0.0, 0.5);
            } else {
                h->Draw("PE SAME");
            }
            leg->AddEntry(h, labs[i].c_str(), "P");
        }
        leg->Draw();
        DrawATLASLabels((TPad*)gPad, ", p_{T} slices");
        c1->SaveAs("outputs/02_1d_rates/sliced/fake_rate_eta_vs_pt_slices.pdf");

        // vs pT sliced in |eta|
        TCanvas* c2 = new TCanvas("c_pt_sliced", "", 700, 600);
        c2->cd();

        vector<pair<double,double>> eta_cuts = {{0,1.37}, {1.52,2.5}};
        vector<int> cols2 = {kRed, kBlue};
        vector<string> labs2 = {"|#eta| < 1.37", "|#eta| > 1.52"};

        TLegend* leg2 = new TLegend(0.60, 0.75, 0.92, 0.90);

        for (size_t i = 0; i < eta_cuts.size(); ++i) {
            int by1 = h2d->GetYaxis()->FindBin(eta_cuts[i].first);
            int by2 = h2d->GetYaxis()->FindBin(eta_cuts[i].second) - 1;

            TH1F* h = (TH1F*)h2d->ProjectionX(Form("px_%d",i), by1, by2);

            h->SetMarkerStyle(20+i);
            h->SetMarkerColor(cols2[i]);
            h->SetLineColor(cols2[i]);

            if (i == 0) {
                h->Draw("PE");
                h->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
                h->GetYaxis()->SetTitle("Fake Rate");
                h->GetYaxis()->SetRangeUser(0.0, 0.5);
            } else {
                h->Draw("PE SAME");
            }
            leg2->AddEntry(h, labs2[i].c_str(), "P");
        }
        leg2->Draw();
        DrawATLASLabels((TPad*)gPad, ", |#eta| slices");
        c2->SaveAs("outputs/02_1d_rates/sliced/fake_rate_pt_vs_eta_slices.pdf");
    }

    f->Close();
    cout << "\nPlots saved in outputs/02_1d_rates/sliced/" << endl;
}