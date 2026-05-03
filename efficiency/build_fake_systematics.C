// =============================================================================
// FILE:       build_fake_systematics.C
// PURPOSE:    Build systematic variations on the 2D fake rate map.
//             Three sources:
//             1. MET-dependence: fake rate measured in CR_MET_low vs CR_MET_high.
//                Captures the change in fake composition (QCD vs W+jets) with MET.
//             2. Composition: MC prompt efficiency compared between CR and SR.
//                Captures the shift in fake population between VBF-exclusive SR and
//                the inclusive SS CR where the rate is measured.
//             3. Statistical: per-bin Poisson errors from limited CR statistics.
//             Total systematic = sqrt(MET^2 + composition^2), stat kept separate.
// INPUTS:     fake_rates.root  (nominal fake rates + MC efficiencies)
// OUTPUTS:    fake_eff_with_systematics.root  (nominal + up/down variations)
// RUNS AFTER: calculate_fake_rates.C
// =============================================================================
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TString.h"
#include "TMath.h"
#include "TLatex.h"
#include "TLegend.h"
#include <iostream>
#include <iomanip>
#include <string>

//======================================================================
// ATLAS Style
//======================================================================
void SetAtlasStyle() {
    gStyle->SetOptStat(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetFrameFillColor(0);
    gStyle->SetFrameLineColor(1);
    gStyle->SetFrameLineWidth(1);
    gStyle->SetPadBorderMode(0);
    gStyle->SetPadColor(0);
    gStyle->SetPadGridX(false);
    gStyle->SetPadGridY(false);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadRightMargin(0.16);
    gStyle->SetPadTopMargin(0.05);
    gStyle->SetPadBottomMargin(0.13);
    gStyle->SetLabelFont(42, "xyz");
    gStyle->SetLabelSize(0.05, "xyz");
    gStyle->SetTitleFont(42, "xyz");
    gStyle->SetTitleSize(0.05, "xyz");
    gStyle->SetTitleOffset(1.4, "y");
    gStyle->SetTitleOffset(1.2, "x");
    gStyle->SetTextFont(42);
    gStyle->SetTextSize(0.05);
    gStyle->SetPalette(57); // viridis
}

//======================================================================
// ATLAS label
//======================================================================
void DrawAtlasLabel(Double_t x = 0.16, Double_t y = 0.93) {
    TLatex tex;
    tex.SetNDC();
    
    // ATLAS Internal - bigger, bold
    tex.SetTextFont(72);
    tex.SetTextSize(0.045);
    tex.DrawLatex(x, y, "ATLAS");
    tex.SetTextFont(62);
    tex.DrawLatex(x + 0.075, y, "Internal");
    
    // Data / lumi - smaller, normal font, moved down a bit
    tex.SetTextFont(42);
    tex.SetTextSize(0.035);
    tex.DrawLatex(x, y - 0.055, "#sqrt{s} = 13 TeV, 139 fb^{-1}");
}

//======================================================================
// Fractional variation map: (alt - nom)/nom
//======================================================================
TH2D* ComputeFractionalVariation(TH2D* nom, TH2D* alt, const TString& name) {
    TH2D* var = (TH2D*)nom->Clone(name);
    var->SetDirectory(nullptr);
    var->Reset();
    for (int ix = 1; ix <= nom->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= nom->GetNbinsY(); ++iy) {
            double n = nom->GetBinContent(ix, iy);
            double a = alt->GetBinContent(ix, iy);
            var->SetBinContent(ix, iy, (n > 0) ? (a - n)/n : 0.0);
        }
    }
    return var;
}

//======================================================================
// Build symmetric up/down from fractional map
//======================================================================
void BuildUpDown(TH2D* nom, TH2D* frac, TH2D*& up, TH2D*& down,
                 const TString& nameUp, const TString& nameDown) {
    up   = (TH2D*)nom->Clone(nameUp);
    down = (TH2D*)nom->Clone(nameDown);
    up->SetDirectory(nullptr);
    down->SetDirectory(nullptr);

    for (int ix = 1; ix <= nom->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= nom->GetNbinsY(); ++iy) {
            double n   = nom->GetBinContent(ix, iy);
            double err = nom->GetBinError(ix, iy);
            // Clip the applied fractional variation at 1.5 so a single
            // pathological bin (very different SR/CR composition or low stats)
            // can't drive the up-side fake rate to >2.5x nominal.
            double f   = TMath::Min(TMath::Abs(frac->GetBinContent(ix, iy)), 1.5);
            double u   = n * (1.0 + f);
            double d   = n * TMath::Max(0.0, 1.0 - f);
            up->SetBinContent(ix, iy, u);
            up->SetBinError(ix, iy, err * (1.0 + f));
            down->SetBinContent(ix, iy, d);
            down->SetBinError(ix, iy, err * TMath::Max(0.0, 1.0 - f));
        }
    }
}

//======================================================================
// Statistical up/down (asymmetric, using stat error only)
//======================================================================
void BuildStatUpDown(TH2D* nom, TH2D*& statUp, TH2D*& statDown) {
    statUp   = (TH2D*)nom->Clone("fake_rate_statUp");
    statDown = (TH2D*)nom->Clone("fake_rate_statDown");
    statUp->SetDirectory(nullptr);
    statDown->SetDirectory(nullptr);

    for (int ix = 1; ix <= nom->GetNbinsX(); ++ix) {
        for (int iy = 1; iy <= nom->GetNbinsY(); ++iy) {
            double c = nom->GetBinContent(ix, iy);
            double e = nom->GetBinError(ix, iy);
            statUp->SetBinContent(ix, iy, c + e);
            statDown->SetBinContent(ix, iy, TMath::Max(0.0, c - e));
        }
    }
}

//======================================================================
void Draw2DMap(TH2D* h, const TString& title, const TString& fname) {
    TCanvas* c = new TCanvas("c", "", 820, 720);
    c->SetLeftMargin(0.15);
    c->SetRightMargin(0.17);
    c->SetBottomMargin(0.14);
    c->SetTopMargin(0.10);

    h->SetTitle("");
    h->GetXaxis()->SetTitle("Electron p_{T} [GeV]");
    h->GetYaxis()->SetTitle("Electron |#eta|");
    h->GetZaxis()->SetTitle("Fake rate");
    h->GetZaxis()->SetTitleOffset(1.35);

    // Smaller font for axes titles and numbers
    h->GetXaxis()->SetTitleSize(0.045);
    h->GetYaxis()->SetTitleSize(0.045);
    h->GetZaxis()->SetTitleSize(0.045);
    
    h->GetXaxis()->SetLabelSize(0.038);
    h->GetYaxis()->SetLabelSize(0.038);
    h->GetZaxis()->SetLabelSize(0.038);

    h->Draw("colz text");

    gPad->SetLogx();

    // Draw bin contents with ± error, 2 decimal places, smaller font
    TLatex latex;
    latex.SetTextAlign(22);          // center
    latex.SetTextFont(42);
    latex.SetTextSize(0.022);        // smaller font for bin text

    int nx = h->GetNbinsX();
    int ny = h->GetNbinsY();

    for (int ix = 1; ix <= nx; ++ix) {
        for (int iy = 1; iy <= ny; ++iy) {
            double val   = h->GetBinContent(ix, iy);
            double err   = h->GetBinError(ix, iy);
            
            if (val <= 0.0001) continue;  // skip empty/very small bins

            double x = h->GetXaxis()->GetBinCenter(ix);
            double y = h->GetYaxis()->GetBinCenter(iy);

            TString label = Form("%.2f#pm%.2f", val, err);
            latex.DrawLatex(x, y, label);
        }
    }

    DrawAtlasLabel(0.16, 0.94);

    c->SaveAs("outputs/syst/" + fname + ".pdf");
    delete c;
}

//======================================================================
void DrawVariation1D(TH1D* nom, TH1D* up, TH1D* down, const TString& systName, const TString& fname) {
    TCanvas* c = new TCanvas("cvar", "", 900, 700);
    c->Divide(1, 2, 0, 0);

    // ────────────────────────────────────────────────
    // Top pad – absolute fake rate
    // ────────────────────────────────────────────────
    TPad* p1 = (TPad*)c->cd(1);
    p1->SetBottomMargin(0.10);
    p1->SetLeftMargin(0.13);
    p1->SetRightMargin(0.05);
    p1->SetTickx(1);
    p1->SetTicky(1);

    // Force x-range 0–70 GeV for all three histos
    nom->GetXaxis()->SetRangeUser(0, 70);
    up->GetXaxis()->SetRangeUser(0, 70);
    down->GetXaxis()->SetRangeUser(0, 70);

    nom->SetTitle("");
    nom->GetXaxis()->SetTitle("");
    nom->GetYaxis()->SetTitle("Fake rate");
    nom->GetYaxis()->SetTitleOffset(1.1);
    nom->SetLineColor(kBlack);
    nom->SetLineWidth(2);
    nom->SetMarkerStyle(20);
    nom->SetMarkerColor(kBlack);
    nom->Draw("E1");

    up->SetLineColor(kRed);
    up->SetLineWidth(2);
    up->SetMarkerStyle(21);
    up->SetMarkerColor(kRed);
    up->Draw("E1 same");

    down->SetLineColor(kAzure+2);
    down->SetLineWidth(2);
    down->SetMarkerStyle(22);
    down->SetMarkerColor(kAzure+2);
    down->Draw("E1 same");

    // (ymax computation left as-is; can be used if you want SetMaximum)

    TLegend* leg = new TLegend(0.60, 0.70, 0.92, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(nom,  "Nominal", "lep");
    leg->AddEntry(up,   systName + " Up",   "lep");
    leg->AddEntry(down, systName + " Down", "lep");
    leg->Draw();

    DrawAtlasLabel(0.17, 0.93);

    // ────────────────────────────────────────────────
    // Bottom pad – ratio
    // ────────────────────────────────────────────────
    TPad* p2 = (TPad*)c->cd(2);
    p2->SetTopMargin(0.05);
    p2->SetBottomMargin(0.42);
    p2->SetLeftMargin(0.13);
    p2->SetRightMargin(0.05);
    p2->SetGridy();

    TH1D* rUp   = (TH1D*)up->Clone("rUp");
    TH1D* rDown = (TH1D*)down->Clone("rDown");
    rUp->Divide(nom);
    rDown->Divide(nom);

    // Force x-range 0–70 GeV for ratios as well
    rUp->GetXaxis()->SetRangeUser(0, 70);
    rDown->GetXaxis()->SetRangeUser(0, 70);

    rUp->SetTitle("");
    rUp->GetXaxis()->SetTitle("Electron p_{T} [GeV]  (|#eta| < 1.8)");
    rUp->GetYaxis()->SetTitle("Variation / nominal");

    rUp->GetYaxis()->SetTitleSize(0.09);
    rUp->GetYaxis()->SetLabelSize(0.055);
    rUp->GetXaxis()->SetTitleSize(0.09);
    rUp->GetXaxis()->SetLabelSize(0.06);

    rDown->GetYaxis()->SetTitleSize(0.09);
    rDown->GetYaxis()->SetLabelSize(0.055);
    rDown->GetXaxis()->SetTitleSize(0.09);
    rDown->GetXaxis()->SetLabelSize(0.06);

    rUp->GetYaxis()->SetRangeUser(0.6, 1.4);

    rUp->SetLineColor(kRed);
    rUp->SetMarkerColor(kRed);
    rUp->SetMarkerStyle(21);
    rUp->Draw("E1");

    rDown->SetLineColor(kAzure+2);
    rDown->SetMarkerColor(kAzure+2);
    rDown->SetMarkerStyle(22);
    rDown->Draw("E1 same");

    TF1* one = new TF1("one","1", 0, 70);
    one->SetLineColor(kBlack);
    one->SetLineStyle(7);
    one->Draw("same");

    c->SaveAs("outputs/syst/" + fname + ".pdf");
    delete c;
}

//======================================================================
// Print table (only for console check)
//======================================================================
void PrintTable(TH2D* h, const TString& title) {
    std::cout << "\n" << std::string(70,'=') << "\n" << title << "\n" << std::string(70,'=') << "\n";
    int nx = h->GetNbinsX(), ny = h->GetNbinsY();
    std::cout << std::setw(12) << "pT bin";
    for (int iy=1; iy<=ny; ++iy) {
        std::cout << std::setw(14) << Form("[%.2f,%.2f]", h->GetYaxis()->GetBinLowEdge(iy), h->GetYaxis()->GetBinUpEdge(iy));
    }
    std::cout << "\n";
    for (int ix=1; ix<=nx; ++ix) {
        std::cout << Form("[%.0f,%.0f] ", h->GetXaxis()->GetBinLowEdge(ix), h->GetXaxis()->GetBinUpEdge(ix));
        for (int iy=1; iy<=ny; ++iy) {
            std::cout << std::setw(14) << Form("%.4f±%.4f", h->GetBinContent(ix,iy), h->GetBinError(ix,iy));
        }
        std::cout << "\n";
    }
}

//======================================================================
// Main function
//======================================================================
void build_fake_systematics() {
    SetAtlasStyle();
    gROOT->ForceStyle();

    gSystem->mkdir("outputs/syst", kTRUE);

    TFile* fIn = TFile::Open("outputs/fake_rates.root", "READ");
    if (!fIn || fIn->IsZombie()) {
        std::cerr << "Cannot open outputs/fake_rates.root\n";
        return;
    }

    TFile* fOut = TFile::Open("outputs/fake_eff_with_systematics.root", "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        std::cerr << "Cannot create output file\n";
        fIn->Close();
        return;
    }

    // Retrieve histograms (nominal + MET-split rates)
    TH2D* nom      = (TH2D*)fIn->Get("fake_rate_2d_CR;2");
    TH2D* met_low  = (TH2D*)fIn->Get("fake_rate_2d_CR_MET_low;2");
    TH2D* met_high = (TH2D*)fIn->Get("fake_rate_2d_CR_MET_high;2");

    if (!nom || !met_low || !met_high) {
        std::cerr << "Missing required histograms in fake_rates.root\n";
        fIn->Close();
        return;
    }

    nom->Sumw2(); met_low->Sumw2(); met_high->Sumw2();

    // MET systematic — selection-dependence test (composition proxy via MET split)
    TH2D* frac_met_low  = ComputeFractionalVariation(nom, met_low,  "frac_MET_low");
    TH2D* frac_met_high = ComputeFractionalVariation(nom, met_high, "frac_MET_high");

    TH2D* frac_met_max = (TH2D*)frac_met_low->Clone("frac_MET_max");
    for (int ix=1; ix<=nom->GetNbinsX(); ++ix)
        for (int iy=1; iy<=nom->GetNbinsY(); ++iy)
            frac_met_max->SetBinContent(ix,iy,
                TMath::Max(TMath::Abs(frac_met_low->GetBinContent(ix,iy)),
                           TMath::Abs(frac_met_high->GetBinContent(ix,iy))));

    TH2D *met_up = nullptr, *met_down = nullptr;
    BuildUpDown(nom, frac_met_max, met_up, met_down, "fake_rate_MET_up", "fake_rate_MET_down");

    // -------------------------------------------------------------------
    // Composition systematic — proper truth-based fake-source reweighting
    // (CR→SR composition shift, ATLAS JINST 18 T11004).
    // Reads frac_comp from outputs/composition_systematic.root
    // (produced by build_composition_systematic.C — must run first).
    //
    // Old proxy (mc_cr/mc_sr prompt-eff ratio) is removed because that
    // captures kinematic VBF-cut effects on REAL electrons, not the change
    // in the FAKE-source mixture, which is what the systematic should test.
    // -------------------------------------------------------------------
    TFile* fComp = TFile::Open("outputs/composition_systematic.root", "READ");
    TH2D* frac_comp = nullptr;
    if (fComp && !fComp->IsZombie()) {
        frac_comp = (TH2D*)fComp->Get("frac_comp");
        if (frac_comp) frac_comp->SetDirectory(nullptr);
    }
    if (!frac_comp) {
        std::cerr << "WARNING: outputs/composition_systematic.root or 'frac_comp' "
                     "not found.\n         Run build_composition_systematic.C first. "
                     "Falling back to zero composition variation.\n";
        frac_comp = (TH2D*)nom->Clone("frac_comp_FALLBACK_ZERO");
        frac_comp->SetDirectory(nullptr);
        frac_comp->Reset();   // all-zero ⇒ no variation
    }

    TH2D *comp_up = nullptr, *comp_down = nullptr;
    BuildUpDown(nom, frac_comp, comp_up, comp_down, "fake_rate_comp_up", "fake_rate_comp_down");

    // -------------------------------------------------------------------
    // Real-lepton (prompt-MC subtraction) systematic
    // -------------------------------------------------------------------
    // ε_f = (data_T − k·MC_T) / (data_L − k·MC_L), with k = 1.0 nominal.
    // We vary k by ±20 % (literature norm — ATLAS HNL multilepton
    // arXiv:2204.11138 quotes 20–30 %) and re-derive ε_f.  The maximum of
    // |Δ_real_up|, |Δ_real_down| in fractional units is then turned into
    // symmetric up/down rate maps, the same way MET and Composition are.
    //
    // Direction convention:
    //   k=0.80 → under-subtract prompt → larger numerator/denominator
    //              residuals → ε_f drifts (sign depends on which factor
    //              the prompt component dominates).
    //   k=1.20 → over-subtract → opposite drift.
    // We take the max of the two as a single fractional variation.
    // -------------------------------------------------------------------
    auto* fData = TFile::Open("outputs/Data_histograms.root", "READ");
    auto* fMC   = TFile::Open("outputs/MC_histograms.root",   "READ");

    TH2D* frac_real = (TH2D*)nom->Clone("frac_real");
    frac_real->SetDirectory(nullptr);
    frac_real->Reset();

    if (!fData || fData->IsZombie() || !fMC || fMC->IsZombie()) {
        std::cerr << "WARNING: Data_histograms.root or MC_histograms.root "
                     "missing — real-lepton systematic skipped (set to 0).\n";
    } else {
        auto* dT = (TH2D*)fData->Get("h2_tight_eta_pt_CR");
        auto* dL = (TH2D*)fData->Get("h2_loose_eta_pt_CR");
        auto* mT = (TH2D*)fMC  ->Get("h2_tight_eta_pt_CR");
        auto* mL = (TH2D*)fMC  ->Get("h2_loose_eta_pt_CR");

        if (!dT || !dL || !mT || !mL) {
            std::cerr << "WARNING: raw CR histograms not found "
                         "(h2_{tight,loose}_eta_pt_CR) — real-lepton "
                         "systematic skipped.\n";
        } else {
            const double kRealNorm = 0.20;   // ±20 %
            for (int ix = 1; ix <= nom->GetNbinsX(); ++ix) {
                for (int iy = 1; iy <= nom->GetNbinsY(); ++iy) {
                    double dT_ = dT->GetBinContent(ix, iy);
                    double dL_ = dL->GetBinContent(ix, iy);
                    double mT_ = mT->GetBinContent(ix, iy);
                    double mL_ = mL->GetBinContent(ix, iy);
                    double eN  = nom->GetBinContent(ix, iy);

                    auto rate = [&](double k) {
                        double num = dT_ - k * mT_;
                        double den = dL_ - k * mL_;
                        return (den > 0.) ? num / den : 0.;
                    };
                    double eUp = rate(1.0 - kRealNorm);   // under-subtract
                    double eDn = rate(1.0 + kRealNorm);   // over-subtract
                    double f   = (eN > 0.)
                                 ? std::max(std::fabs(eUp - eN),
                                            std::fabs(eDn - eN)) / eN
                                 : 0.;
                    if (!std::isfinite(f)) f = 0.;
                    frac_real->SetBinContent(ix, iy, f);
                }
            }
        }
    }
    if (fData) fData->Close();
    if (fMC)   fMC->Close();

    TH2D *real_up = nullptr, *real_down = nullptr;
    BuildUpDown(nom, frac_real, real_up, real_down,
                "fake_rate_real_up", "fake_rate_real_down");

    // Statistical variation
    TH2D *stat_up = nullptr, *stat_down = nullptr;
    BuildStatUpDown(nom, stat_up, stat_down);

    // Total systematic (MET + composition in quadrature)
    TH2D* tot_up   = (TH2D*)nom->Clone("fake_rate_total_systUp");
    TH2D* tot_down = (TH2D*)nom->Clone("fake_rate_total_systDown");
    tot_up->SetDirectory(nullptr);
    tot_down->SetDirectory(nullptr);

    // Total systematic = quadrature sum of MET-split + composition fractional
    // variations.  No special-case override on the last pT bin: the merged
    // [50,200] GeV bin (kPtBins update in SelectionUtils.h) now has enough
    // statistics for a normal MET-split test, so the previous flat 100 %
    // top-bin treatment is no longer needed.  If a future bin redefinition
    // produces a low-stat tail, fall back to a stat-driven cap (e.g., raise
    // ftot ⇒ 1.0 only when σ_stat / ε_f > 50 %) rather than hard-coding by
    // bin index.
    for (int ix=1; ix<=nom->GetNbinsX(); ++ix) {
        for (int iy=1; iy<=nom->GetNbinsY(); ++iy) {
            double n   = nom->GetBinContent(ix,iy);
            double err = nom->GetBinError(ix,iy);

            double fm  = TMath::Abs(frac_met_max->GetBinContent(ix,iy));
            double fc  = TMath::Abs(frac_comp   ->GetBinContent(ix,iy));
            double fr  = TMath::Abs(frac_real   ->GetBinContent(ix,iy));
            double ftot = TMath::Sqrt(fm*fm + fc*fc + fr*fr);

            double u = n * (1 + ftot);
            double d = n * TMath::Max(0., 1 - ftot);
            tot_up->SetBinContent(ix,iy, u);
            tot_up->SetBinError(ix,iy, err * (1 + ftot));
            tot_down->SetBinContent(ix,iy, d);
            tot_down->SetBinError(ix,iy, err * TMath::Max(0., 1 - ftot));
        }
    }

    // 2D maps
    Draw2DMap(nom,      "Nominal fake rate (CR)",               "fake_rate_nominal");
    Draw2DMap(met_up,   "Fake rate  MET syst Up",               "fake_rate_MET_up");
    Draw2DMap(met_down, "Fake rate  MET syst Down",             "fake_rate_MET_down");
    Draw2DMap(comp_up,  "Fake rate  Composition syst Up",       "fake_rate_comp_up");
    Draw2DMap(comp_down,"Fake rate  Composition syst Down",     "fake_rate_comp_down");
    Draw2DMap(real_up,  "Fake rate  Real-lep syst Up (#pm20%)", "fake_rate_real_up");
    Draw2DMap(real_down,"Fake rate  Real-lep syst Down (#pm20%)","fake_rate_real_down");
    Draw2DMap(tot_up,   "Fake rate  Total syst Up (MET+Comp+Real)",  "fake_rate_total_up");
    Draw2DMap(tot_down, "Fake rate  Total syst Down (MET+Comp+Real)","fake_rate_total_down");

    // Central η only (|η| < 1.8) – project and then show only 0–70 GeV
    TH1D* central_nom   = nom->ProjectionX("central_nom", 1, 1);
    TH1D* central_met_u = met_up->ProjectionX("central_met_u", 1, 1);
    TH1D* central_met_d = met_down->ProjectionX("central_met_d", 1, 1);

    central_nom->GetXaxis()->SetRangeUser(0, 70);
    central_met_u->GetXaxis()->SetRangeUser(0, 70);
    central_met_d->GetXaxis()->SetRangeUser(0, 70);

    central_nom->SetBinContent(central_nom->GetNbinsX()+1, 0);
    central_met_u->SetBinContent(central_met_u->GetNbinsX()+1, 0);
    central_met_d->SetBinContent(central_met_d->GetNbinsX()+1, 0);

    DrawVariation1D(central_nom, central_met_u, central_met_d, "MET", "variation_pt_MET_central");

    TH1D* central_comp_u = comp_up->ProjectionX("central_comp_u", 1, 1);
    TH1D* central_comp_d = comp_down->ProjectionX("central_comp_d", 1, 1);
    central_comp_u->GetXaxis()->SetRangeUser(0, 70);
    central_comp_d->GetXaxis()->SetRangeUser(0, 70);
    DrawVariation1D(central_nom, central_comp_u, central_comp_d, "Composition", "variation_pt_comp_central");

    TH1D* central_tot_u = tot_up->ProjectionX("central_tot_u", 1, 1);
    TH1D* central_tot_d = tot_down->ProjectionX("central_tot_d", 1, 1);
    central_tot_u->GetXaxis()->SetRangeUser(0, 70);
    central_tot_d->GetXaxis()->SetRangeUser(0, 70);
    DrawVariation1D(central_nom, central_tot_u, central_tot_d, "Total Syst", "variation_pt_total_syst_central");

    // Console summary
    PrintTable(nom,      "NOMINAL");
    PrintTable(stat_up,  "STAT UP");
    PrintTable(stat_down,"STAT DOWN");
    PrintTable(met_up,   "MET UP");
    PrintTable(met_down, "MET DOWN");
    PrintTable(comp_up,  "COMPOSITION UP");
    PrintTable(comp_down,"COMPOSITION DOWN");
    PrintTable(real_up,  "REAL-LEPTON UP (prompt-MC norm +20%)");
    PrintTable(real_down,"REAL-LEPTON DOWN (prompt-MC norm -20%)");
    PrintTable(tot_up,   "TOTAL SYST UP   (MET (+) Comp (+) Real)");
    PrintTable(tot_down, "TOTAL SYST DOWN (MET (+) Comp (+) Real)");

    // Save to ROOT file
    fOut->cd();
    nom->Write("fake_rate_nominal");
    met_up->Write();
    met_down->Write();
    comp_up->Write();
    comp_down->Write();
    real_up->Write();
    real_down->Write();
    stat_up->Write();
    stat_down->Write();
    tot_up->Write();
    tot_down->Write();

    frac_met_max->Write("frac_MET_max");
    frac_comp   ->Write("frac_comp");
    frac_real   ->Write("frac_real");

    fOut->Close();
    fIn->Close();

    std::cout << "\nDone.\n";
    std::cout << "   Plots   →  ./syst/*.pdf\n";
    std::cout << "   ROOT file →  fake_eff_with_systematics.root\n";
}
