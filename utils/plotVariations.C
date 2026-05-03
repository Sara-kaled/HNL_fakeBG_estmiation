#include <iostream>
#include <vector>
#include <string>
#include "TStyle.h"
#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TMath.h"

static bool print = false;

void drawATLASLabel(TH1* h){
  const char *flavor(""), *type("");
  TString name = Form("%s",h->GetName());
  if( name.Contains("Fake") ) type   = "fake";
  if( name.Contains("Real") ) type   = "real";
  if( name.Contains("_mu")  ) flavor = "muon";
  if( name.Contains("_el")  ) flavor = "electron";

  TLatex n;
  n.SetNDC();
  n.SetTextColor(kBlack);
  n.SetTextSize(h->GetDimension()==1 ? 0.05 : 0.04);
  n.SetTextFont(72); 
  n.DrawLatex(0.12,0.92,"ATLAS");
  n.SetTextFont(42); 
  n.DrawLatex(0.22,0.92,"Internal");
  n.SetTextFont(42);
  n.DrawLatex(0.35, 0.92, Form("Rel. uncertainties of %s %s efficiencies [%s]", flavor, type, "%"));
  return;
}

void drawGridLines(TH1D* h){
  if(!h) return;
  for(int i(2); i<=h->GetNbinsX(); i++){
    float x = h->GetXaxis()->GetBinLowEdge(i);

    TLine xLine;
    xLine.SetLineWidth(1);
    xLine.SetLineColor(kGray+1);
    xLine.SetLineStyle(7);
    xLine.DrawLine(x, -109, x, 109.);
  }
  TLine yLine;
  yLine.SetLineWidth(2);
  yLine.SetLineColor(kBlue-7);
  yLine.SetLineStyle(9);
  yLine.DrawLine(h->GetXaxis()->GetXmin(), 0.0, h->GetXaxis()->GetXmax(), 0.0);
}

void draw2DGridLines(TH2F* h){
  if(!h) return;
  TLine line;
  line.SetLineWidth(1);
  line.SetLineColor(kWhite);
  line.SetLineStyle(1);

  for(int i(2); i<=h->GetNbinsX(); i++){
    float x = h->GetXaxis()->GetBinLowEdge(i);
    line.DrawLine(x, 0, x, h->GetYaxis()->GetXmax());
  }
  for(int i(2); i<=h->GetNbinsY(); i++){
    float y = h->GetYaxis()->GetBinLowEdge(i);
    line.DrawLine(0, y, h->GetXaxis()->GetXmax(), y);
  }
}

bool isTotalVar(TH1* h){
  if(!h) return false;
  TString name = h->GetName();
  return (bool)name.Contains("TOTAL");
}

float averageEta(int bin, TH2F *h){
  bool isMu = ((TString)h->GetName()).Contains("_mu") ? true : false;

  std::vector<float> wMU = {0.225, 0.237, 0.209, 0.191, 0.137};
  std::vector<float> wEL = {0.294, 0.312, 0.080, 0.190, 0.123};
  std::vector<float> weights = isMu ? wMU : wEL;

  float av(0);
  for(int i(1); i<=h->GetNbinsY(); i++) 
    av += h->GetBinContent(bin, i)*weights[i-1];
  return av;
}

float averagePt(int bin, TH2F *h){
  bool isMu = ((TString)h->GetName()).Contains("_mu") ? true : false;

  std::vector<float> wMU = {0.508, 0.307, 0.149, 0.027, 0.009};
  std::vector<float> wEL = {0.456, 0.243, 0.157, 0.056, 0.088};
  std::vector<float> weights = isMu ? wMU : wEL;

  float av(0);
  for(int i(1); i<=h->GetNbinsX(); i++)
    av += h->GetBinContent(i, bin)*weights[i-1];
  return av;
}


void checkTitleUnit(TH1 *h){
  if(!h) return;
  TString xtitle(""), ytitle("");
  
  switch( (int)h->GetDimension() ){
  case 1:
    xtitle = h->GetXaxis()->GetTitle();
    if( xtitle.Contains("p_{T}") && !xtitle.Contains("GeV")){ xtitle += " [GeV]"; h->GetXaxis()->SetTitle(xtitle);} 
    break;
  case 2:
    xtitle = h->GetXaxis()->GetTitle();
    ytitle = h->GetYaxis()->GetTitle();
    if( xtitle.Contains("p_{T}") && !xtitle.Contains("GeV")){ xtitle += " [GeV]"; h->GetXaxis()->SetTitle(xtitle);}
    if( ytitle.Contains("p_{T}") && !ytitle.Contains("GeV")){ ytitle += " [GeV]"; h->GetXaxis()->SetTitle(ytitle);}
    break;
  }
  return;
}

void setProjStyle(TH1D* h){
  if(!h) return;
  h->SetFillColor(kBlue-10);
  h->SetLineColor(h->GetFillColor());
  h->SetLineWidth(1);
  h->SetMarkerColor(h->GetLineColor());
  h->SetMarkerStyle(1);
  h->SetMarkerSize(0);
  h->SetFillStyle(1111);
  h->GetYaxis()->SetRangeUser(-109.,109.);
  h->GetYaxis()->SetNdivisions(20);
  h->GetYaxis()->SetTitle( Form("Rel. uncertainty [%s]","%"));
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.);
  h->GetXaxis()->SetTitleOffset(1.3);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelOffset(0.01);

}

void plot(TH1F *h){
  if(!h) return;
  TCanvas *c = new TCanvas(Form("c_%s",h->GetName()), Form("c_%s",h->GetName()));
  c->cd();
  c->SetTicks(1,1);
  if( ((TString)h->GetXaxis()->GetTitle()).Contains("p_{T}") ) c->SetLogx();
  h->GetXaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->SetLabelOffset(0.001);
  checkTitleUnit(h);

  h->Draw("HIST");
  gPad->RedrawAxis();

  drawGridLines((TH1D*)h);
  drawATLASLabel(h);
  if(print) 
    c->Print(Form("%s.eps", c->GetName()));
  return;
}

void plot(TH2F *h){
  if(!h) return;
  TCanvas *c = new TCanvas(Form("c_%s",h->GetName()), Form("c_%s",h->GetName()));
  c->cd();
  c->SetTicks(1,1);

  gStyle->SetOptStat(0);
  gStyle->SetPaintTextFormat(".1f %%");
  gStyle->SetNumberContours(40);
  gStyle->SetPalette(104);
  if( ((TString)h->GetXaxis()->GetTitle()).Contains("p_{T}") ) c->SetLogx();
  if( ((TString)h->GetYaxis()->GetTitle()).Contains("p_{T}") ) c->SetLogy();
  h->GetXaxis()->SetTitleOffset(1.15);
  h->GetXaxis()->SetLabelOffset(0.001);
  checkTitleUnit(h);

  h->Draw("COLZ TEXT45");
  gPad->RedrawAxis();

  draw2DGridLines(h);
  drawATLASLabel(h);
  if(print) 
    c->Print(Form("%s.eps",c->GetName()));
  return;
}

void plotProjections(TH2F *h){
  if(!h) return;
  TCanvas *cx = new TCanvas(Form("ProjX_%s",h->GetName()),Form("ProjX_%s",h->GetName()), 1, 10, 770, 460);
  TCanvas *cy = new TCanvas(Form("ProjY_%s",h->GetName()),Form("ProjY_%s",h->GetName()), 1, 10, 770, 460);

  TPad *px = new TPad(Form("%s_p",cx->GetName()),Form("%s_p",cx->GetName()), 0.00, 0.00, 1.00, 0.95, -1, 0, 0);
  TPad *py = new TPad(Form("%s_p",cy->GetName()),Form("%s_p",cy->GetName()), 0.00, 0.00, 1.00, 0.95, -1, 0, 0);

  px->SetLeftMargin(0.10);
  py->SetLeftMargin(0.10);
  px->SetBottomMargin(0.15);
  py->SetBottomMargin(0.15);
  px->SetTicks(1,1);
  py->SetTicks(1,1);  

  TH1D *hx = h->ProjectionX( Form("%s_pX",h->GetName()) );
  TH1D *hy = h->ProjectionY( Form("%s_pY",h->GetName()) );
  hx->Reset();
  hy->Reset();

  for(int i(1); i<=hx->GetNbinsX(); i++) hx->SetBinError(i, averageEta(i,h) );
  for(int i(1); i<=hy->GetNbinsX(); i++) hy->SetBinError(i, averagePt(i,h) );

  cx->cd();
  px->Draw();
  px->SetLogx();
  px->cd();

  setProjStyle(hx);
  hx->Draw("E2");

  drawGridLines(hx);
  drawATLASLabel(hx);
  gPad->RedrawAxis();

  cy->cd();
  py->Draw();
  py->cd();

  setProjStyle(hy);
  hy->Draw("E2");

  drawGridLines(hy);
  drawATLASLabel(hy);
  gPad->RedrawAxis();

  if(print){
    cx->Print(Form("%s.eps",cx->GetName()));
    cy->Print(Form("%s.eps",cy->GetName()));
  }
}

void makeDiffPlot(TH1F* nom, TH1F* var){
  if(!nom || !var) return;
  std::cout << Form("Plotting relative variations between %s and %s", nom->GetName(), var->GetName()) << std::endl;

  TH1F *hTemp = (TH1F*)nom->Clone(Form("Diff_%s", var->GetName()));
  hTemp->Reset();
  hTemp->GetYaxis()->SetRangeUser(0,101);
  
  if(!isTotalVar(var)){
    for(int i(1); i<=hTemp->GetNbinsX(); i++){
      float ratio = var->GetBinContent(i);
      ratio = 100 * TMath::Min(ratio, (float)1.);
      
      std::cout << Form("Bin (%i) Rel. uncertainty = %.1f %s", i, ratio, "%") << std::endl;
      hTemp->SetBinContent(i,ratio);
    }
    plot(hTemp);
    return;
  }

  for(int i(1); i<=hTemp->GetNbinsX();i++){
    float ratio = nom->GetBinContent(i)>0 ? (var->GetBinContent(i) / nom->GetBinContent(i)) : 1.;
    ratio = 100 * TMath::Min(ratio, (float)1.);

    std::cout << Form("Bin (%i) : Rel. uncertainty = %.1f %s", i, ratio, "%") << std::endl;
    hTemp->SetBinContent(i,ratio);
  }
  plot(hTemp);
  return;
}

void makeDiffPlot(TH2F* nom, TH2F* var){
  if(!nom || !var) return;
  std::cout << Form("Plotting relative variations between %s and %s", nom->GetName(), var->GetName()) << std::endl;

  TH2F *hTemp = (TH2F*)nom->Clone(Form("Diff_%s", var->GetName()));
  hTemp->Reset();
  hTemp->GetZaxis()->SetRangeUser(0,101);

  if(!isTotalVar(var)){
    for(int i(1); i<=hTemp->GetNbinsX(); i++){
      for(int j(1); j<=hTemp->GetNbinsY(); j++){
	float ratio = var->GetBinContent(i,j);
	ratio = 100 * TMath::Min(ratio, (float)1.);

	std::cout << Form("Bin (%i|%i) Rel. uncertainty = %.1f %s", i, j, ratio, "%") << std::endl;
	hTemp->SetBinContent(i,j,ratio);
      }
    }
    plot(hTemp);
    plotProjections(hTemp);
    return;
  }
  
  for(int i(1); i<=hTemp->GetNbinsX(); i++){
    for(int j(1); j<=hTemp->GetNbinsY(); j++){
      float ratio = nom->GetBinContent(i,j)>0 ? (var->GetBinContent(i,j) / nom->GetBinContent(i,j)) : 1.;
      ratio = 100 * TMath::Min(ratio, (float)1.);
      
      std::cout << Form("Bin (%i|%i) [var=%.3f|nom=%.3f]: Rel. uncertainty = %.1f %s", i, j, var->GetBinContent(i,j), nom->GetBinContent(i,j), ratio, "%") << std::endl;
      hTemp->SetBinContent(i,j,ratio);
    }
  }
  plot(hTemp);
  plotProjections(hTemp);
  return;
}

void getUnc1D(TFile *f, std::vector<TString> names, TString var){
  if( !names.size() ) return;
  if( !var.Length() ){ std::cout << "No variation selected" << std::endl; return; } 

  TH1F* Nom_FakeMu(0), *Nom_FakeEl(0), *Nom_RealEl(0), *Nom_RealMu(0);
  TH1F* Var_FakeMu(0), *Var_FakeEl(0), *Var_RealEl(0), *Var_RealMu(0);

  for(auto name : names){
    if(!Nom_FakeMu && name.Contains("Fake") && name.Contains("mu") && !name.Contains(var)) Nom_FakeMu = (TH1F*)f->Get(name);
    if(!Nom_FakeEl && name.Contains("Fake") && name.Contains("el") && !name.Contains(var)) Nom_FakeEl = (TH1F*)f->Get(name);
    if(!Nom_RealMu && name.Contains("Real") && name.Contains("mu") && !name.Contains(var)) Nom_RealMu = (TH1F*)f->Get(name);
    if(!Nom_RealEl && name.Contains("Real") && name.Contains("el") && !name.Contains(var)) Nom_RealEl = (TH1F*)f->Get(name);

    if(!Var_FakeMu && name.Contains("Fake") && name.Contains("mu") &&  name.Contains(var)) Var_FakeMu = (TH1F*)f->Get(name);
    if(!Var_FakeEl && name.Contains("Fake") && name.Contains("el") &&  name.Contains(var)) Var_FakeEl = (TH1F*)f->Get(name);
    if(!Var_RealMu && name.Contains("Real") && name.Contains("mu") &&  name.Contains(var)) Var_RealMu = (TH1F*)f->Get(name);
    if(!Var_RealEl && name.Contains("Real") && name.Contains("el") &&  name.Contains(var)) Var_RealEl = (TH1F*)f->Get(name);
  }
  makeDiffPlot(Nom_FakeMu, Var_FakeMu);
  makeDiffPlot(Nom_FakeEl, Var_FakeEl);
  makeDiffPlot(Nom_RealMu, Var_RealMu);
  makeDiffPlot(Nom_RealEl, Var_RealEl);
  return;
}

void getUnc2D(TFile *f, std::vector<TString> names, TString var){
  if( !names.size() ) return;
  if( !var.Length() ){std::cout << "No variation selected" <<std::endl; return; }
  
  TH2F* Nom_FakeMu(0), *Nom_FakeEl(0), *Nom_RealEl(0), *Nom_RealMu(0);
  TH2F* Var_FakeMu(0), *Var_FakeEl(0), *Var_RealEl(0), *Var_RealMu(0);

  for(auto name : names){
    if(!Nom_FakeMu && name.Contains("Fake") && name.Contains("mu") && !name.Contains(var)) Nom_FakeMu = (TH2F*)f->Get(name);
    if(!Nom_FakeEl && name.Contains("Fake") && name.Contains("el") && !name.Contains(var)) Nom_FakeEl = (TH2F*)f->Get(name);
    if(!Nom_RealMu && name.Contains("Real") && name.Contains("mu") && !name.Contains(var)) Nom_RealMu = (TH2F*)f->Get(name);
    if(!Nom_RealEl && name.Contains("Real") && name.Contains("el") && !name.Contains(var)) Nom_RealEl = (TH2F*)f->Get(name);

    if(!Var_FakeMu && name.Contains("Fake") && name.Contains("mu") &&  name.Contains(var)) Var_FakeMu = (TH2F*)f->Get(name);
    if(!Var_FakeEl && name.Contains("Fake") && name.Contains("el") &&  name.Contains(var)) Var_FakeEl = (TH2F*)f->Get(name);
    if(!Var_RealMu && name.Contains("Real") && name.Contains("mu") &&  name.Contains(var)) Var_RealMu = (TH2F*)f->Get(name);
    if(!Var_RealEl && name.Contains("Real") && name.Contains("el") &&  name.Contains(var)) Var_RealEl = (TH2F*)f->Get(name);
  }
  makeDiffPlot(Nom_FakeMu, Var_FakeMu);
  makeDiffPlot(Nom_FakeEl, Var_FakeEl);
  makeDiffPlot(Nom_RealMu, Var_RealMu);
  makeDiffPlot(Nom_RealEl, Var_RealEl);
  return;
}

void plotVariations(const char *filename, const char *variation="", bool doPrint=false){

  TFile *f = new TFile(filename); 
  if(!f){ std::cout << Form("File %s not found", filename) << std::endl; return; }
  std::cout << Form("Found %s", f->GetName()) << std::endl;

  std::vector<TString> names1D(0),names2D(0);
  TKey *key(0);
  TList* Objects = f->GetListOfKeys();
  Objects->Sort();  

  TIter next(Objects);
  while(( key = (TKey*)next() )){
    TObject *obj = key->ReadObj();
    TH1 *h = (TH1*)obj;
    switch( h->GetDimension() ){
    case 1: names1D.push_back( (TString)h->GetName() ); break;
    case 2: names2D.push_back( (TString)h->GetName() ); break;
    default: break;
    }
  }
  std::cout << Form("Histograms found (1D|2D) : (%i|%i)", (int)names1D.size(),(int)names2D.size()) << std::endl;

  print = doPrint;
  getUnc1D(f, names1D, Form("%s",variation));
  getUnc2D(f, names2D, Form("%s",variation));
  return;
}

