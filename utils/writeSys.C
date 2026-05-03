#include <iostream>
#include <vector>
#include "TFile.h"
#include "TString.h"
#include "TKey.h"
#include "TList.h"
#include "TObject.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TMath.h"

void INFO(const char* msg){std::cout << Form("writeSys::INFO\t %s",msg) << std::endl;}

void ERROR(const char *msg){std::cout << Form("writeSys::ERROR\t %s",msg) << std::endl; exit(1);}

bool isPtHist(TH1F* h){
  if(!h) return false;
  if( ((TString)h->GetName()).Contains("pt")) return true;
  return false;
}

bool errorCheck(int i, TH1F* h){
  if(!h) return false;
  if(h->GetBinContent(i) - h->GetBinError(i) < 0.0001){
    INFO( Form("Histogram %s bin(%i) = %.4f +/- %.4f setting uncertainty to 1.0",h->GetName(),i,h->GetBinContent(i),h->GetBinError(i)) );
    return false;
  }
  return true;
}

bool errorCheck(int i, int j, TH2F *h){
  if(!h) return false;
  if(h->GetBinContent(i,j) - h->GetBinError(i,j) < 0.0001){
    INFO( Form("Histogram %s bin(%i|%i) = %.4f +/- %.4f setting uncertainty to 1.0",h->GetName(),i,j,h->GetBinContent(i,j),h->GetBinError(i,j)) );
    return false;
  }
  return true;
}

const char* getXTitle(TH1F* h){
  TString name = Form("%s",h->GetName());
  if( name.Contains("el_pt") )  return "Electron p_{T} [GeV]";
  if( name.Contains("mu_pt") )  return "Muon p_{T} [GeV]";
  if( name.Contains("el_eta") ) return "Electron |#eta|";
  if( name.Contains("mu_eta") ) return "Muon |#eta|";
  return "";
}

const char* makeLegEntry(const char *suf){
  TString s = Form("%s",suf);
  if(s == "SS_3L")  return "SS/3L";
  if(s == "CF")     return "charge-flip";
  if(s == "nBJets") return "N_{b-jets}";
  if(s == "prompt") return "prompt";
  return "";
}

void setStyle1D(TH1F* h){
  if(!h) return;
  h->SetFillColor(kBlue-10);
  h->SetLineColor(h->GetFillColor());
  h->SetLineWidth(1);
  h->SetMarkerColor(h->GetLineColor());
  h->SetMarkerStyle(1);
  h->SetMarkerSize(0);
  h->SetFillStyle(1111);
  h->GetYaxis()->SetRangeUser(-109.,149.);
  h->GetYaxis()->SetNdivisions(20);
  h->GetYaxis()->SetTitle( Form("Rel. uncertainty [%s]","%"));
  h->GetXaxis()->SetTitle( getXTitle(h) );
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelOffset(0.01);
}

void setStyle2D(TH2F* h){
  if(!h) return;
  TString name = Form("%s",h->GetName());
  
  gStyle->SetPaintTextFormat(".1f %%");
  gStyle->SetNumberContours(40);
  gStyle->SetPalette(104);

  h->GetZaxis()->SetRangeUser(0.01, 99.);
  h->GetZaxis()->SetLabelSize(0.04);

  if( name.Contains("mu_pt_eta") ){ h->GetXaxis()->SetTitle("Muon p_{T} [GeV]");     h->GetYaxis()->SetTitle("Muon |#eta|"); }
  if( name.Contains("el_pt_eta") ){ h->GetXaxis()->SetTitle("Electron p_{T} [GeV]"); h->GetYaxis()->SetTitle("Electron |#eta|"); }
}

void drawGridLines(TH1F* h){
  if(!h) return;
  for(int i(2); i<=h->GetNbinsX(); i++){
    float x = h->GetXaxis()->GetBinLowEdge(i);
        
    TLine xLine;
    xLine.SetLineWidth(1);
    xLine.SetLineColor(kGray+1);
    xLine.SetLineStyle(7);
    xLine.DrawLine(x, -109., x, 149.);
  }

  TLine yLine;
  yLine.SetLineWidth(2);
  yLine.SetLineColor(kBlue-7);
  yLine.SetLineStyle(9);
  yLine.DrawLine(h->GetXaxis()->GetXmin(), 0.0, h->GetXaxis()->GetXmax(), 0.0);

}

void drawGridLines2D(TH2F* h){
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

void drawLeptonFlavor(TH1F *h){
  if(!h) return;
  std::string flavor("");
  TString name = Form("%s",h->GetName());
  if( name.Contains("_mu") ) flavor = "Muons";
  if( name.Contains("_el") ) flavor = "Electrons";
  if(!flavor.length()) return;

  TLegend *l = new TLegend(0.14, 0.66, 0.40, 0.84);
  l->SetFillStyle(1111);
  l->SetFillColor(kWhite);
  l->SetBorderSize(1);
  l->Draw("SAME");

  TLatex n;
  n.SetNDC();
  n.SetTextColor(kBlack);;
  n.SetTextFont(62);
  n.SetTextSize(0.06);
  n.DrawLatex(0.16, 0.69, flavor.c_str());

  n.SetTextFont(72);
  n.DrawLatex(0.16, 0.77, "ATLAS");
  n.SetTextFont(42);
  n.DrawLatex(0.27, 0.77, "Internal");
  return;
}

void draw2DplotLabel(const char* entry){

  TLatex n;
  n.SetNDC();
  n.SetTextColor(kBlack);;
  n.SetTextSize(0.04);
  n.SetTextFont(72); n.DrawLatex(0.11,0.93,"ATLAS");
  n.SetTextFont(42); n.DrawLatex(0.20,0.93,"Internal");
  n.SetTextSize(0.04);
  n.SetTextFont(42);
  n.DrawLatex(0.31, 0.93, Form("Relative efficiency uncertainty %s [%s]", makeLegEntry(entry), "%"));
}

std::vector<TH1F*> getHistos1D(TFile *f, TString name=""){
  std::vector<TH1F*> H;
  if(!f) return H;

  TKey *key(0);
  TList *Objects = f->GetListOfKeys();
  Objects->Sort();

  TIter next(Objects);
  while( (key = (TKey*)next()) ){
    TObject *obj = key->ReadObj();
    if(!obj->InheritsFrom("TH1")) continue;
    TH1F* h = (TH1F*)obj;
    TString hName = h->GetName();
    if(name.Length() && !hName.Contains(name)) continue;
   
    INFO( Form("Found %s/%s",f->GetName(),hName.Data()) );
    H.push_back(h);
  }
  return H;
}

std::vector<TH2F*> getHistos2D(TFile *f, TString name=""){
  std::vector<TH2F*> H;
  if(!f) return H;

  TKey *key(0);
  TList *Objects = f->GetListOfKeys();
  Objects->Sort();

  TIter next(Objects);
  while( (key = (TKey*)next()) ){
    TObject *obj = key->ReadObj();
    if(!obj->InheritsFrom("TH2")) continue;
    TH2F* h = (TH2F*)obj;
    TString hName = h->GetName();
    if(name.Length() && !hName.Contains(name)) continue;

    INFO( Form("Found %s/%s",f->GetName(),hName.Data()) );
    H.push_back(h);
  }
  return H;
}

TH1F *hDiff1D(TH1F* h1, TH1F* h2){
  if(!h1 || !h2) return nullptr;
  TH1F *h = (TH1F*)h1->Clone( Form("Diff_%s_%s",h1->GetName(),h2->GetName()) );
  h->Reset();
  
  for(int i(1); i<=h->GetNbinsX(); i++){
    float mean = h1->GetBinContent(i) + h2->GetBinContent(i);
    float err  = TMath::Abs(h1->GetBinContent(i) - h2->GetBinContent(i));
   
    float ratio = (errorCheck(i,h1) && errorCheck(i,h2)) ? err/mean : 1;
    h->SetBinContent(i, ratio);
  }
  INFO( Form("Created %s",h->GetName()) );
  return h;
}

TH2F *hDiff2D(TH2F* h1, TH2F* h2){
  if(!h1 || !h2) return nullptr;
  TH2F *h = (TH2F*)h2->Clone( Form("Diff_%s_%s",h1->GetName(),h2->GetName()) );
  h->Reset();

  for(int i(1); i<=h->GetNbinsX(); i++){
    for(int j(1); j<=h->GetNbinsY(); j++){
      float mean = h1->GetBinContent(i,j) + h2->GetBinContent(i,j);
      float err  = TMath::Abs(h1->GetBinContent(i,j) - h2->GetBinContent(i,j));
    
      float ratio = (errorCheck(i,j,h1) && errorCheck(i,j,h2)) ? err/mean : 1;
      h->SetBinContent(i,j, ratio);
    }
  }
  INFO( Form("Created %s",h->GetName()) );
  return h;
}

TH1F * getErrorHist(TH1F *hTemp){
  if(!hTemp) return nullptr;
  TH1F *h = (TH1F*)hTemp->Clone( Form("Err_%s",hTemp->GetName()));
  h->Reset();
  for(int i(1); i<=h->GetNbinsX(); i++)  h->SetBinContent(i, hTemp->GetBinError(i));
  return h;
}

void writeToFile(TH1 *h, const char* suffix){

  TFile *f(0);
  const char *filename("");
  TString name = h->GetName();
  const char *type(""), *flavor(""), *var("");

  if(name.Contains("Fake")) type = "Fake";
  if(name.Contains("Real")) type = "Real";

  if(name.Contains("_mu")) flavor = "mu";
  if(name.Contains("_el")) flavor = "el";
  
  switch( (int)h->GetDimension() ){
  case 1:
    filename = "Efficiency1D.root";
    f = TFile::Open(filename, "UPDATE");
    f->cd();
    if(name.Contains("_pt"))  var = "pt";
    if(name.Contains("_eta")) var = "eta";
    
    h->SetName( Form("%sEfficiency1D_%s_%s__%s",type,flavor,var,suffix) );
    h->Write();
    INFO( Form("Wrote %s/%s",filename,h->GetName()) );
    break;
  case 2:
    filename = "Efficiency2D.root";
    f = TFile::Open(filename, "UPDATE");
    f->cd();
    if(name.Contains("_pt_eta"))  var = "pt_eta";
    
    h->SetName( Form("%sEfficiency2D_%s_%s__%s",type,flavor,var,suffix) );
    h->Write();
    INFO( Form("Wrote %s/%s",filename,h->GetName()) );
    break;
  default: break;
  }
  if(f->IsOpen()) f->Close();
  return;
}

void write1D(TFile* f1, TFile *f2,  const char *sysName="", bool write=false){
  
  std::vector<TH1F*> H1 = getHistos1D(f1);
  std::vector<TH1F*> H2 = getHistos1D(f2);
  
  if(H1.size() != H2.size())
    ERROR( Form("Inconsistent file content: %i != %i",(int)H1.size(),(int)H2.size()) );

  for(unsigned int i(0); i<H1.size(); i++){
    
    TH1F *h = hDiff1D(H1.at(i),H2.at(i));
    setStyle1D(h);

    TString cname = Form("C_%s_%s",sysName,h->GetName());
    TCanvas * c = new TCanvas(cname, cname, 1, 10, 770, 460);
    TPad *p = new TPad(cname+"_p",cname+"_p", 0.00, 0.00, 1.00, 0.95, -1, 0, 0);

    p->SetLeftMargin(0.10);
    p->SetBottomMargin(0.15);
    p->Draw();
    p->SetTicks(1,1);
    if( isPtHist(h) ) p->SetLogx(); 
    p->cd();

    TH1F *hRel = (TH1F*)h->Clone( Form("Rel_%s",h->GetName()) );
    hRel->Reset();
    for(int i(1); i<=hRel->GetNbinsX(); i++){
      hRel->SetBinContent(i,0.);
      hRel->SetBinError(i, h->GetBinContent(i)*100);;
    }
    hRel->Draw("E2");

    drawGridLines(h);
    drawLeptonFlavor(h);

    TLegend *leg = new TLegend(0.45, 0.77, 0.85, 0.84);
    leg->SetNColumns(1);
    leg->SetBorderSize(1);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->SetFillStyle(1111);
    leg->SetFillColor(kWhite);
    leg->AddEntry(h, Form(" Rel. uncertainty %s", makeLegEntry(sysName)) );
    leg->Draw("SAME");
    gPad->RedrawAxis();

    if(write){
      c->Print(Form("%s.eps",cname.Data()));
      writeToFile(h, sysName);
    }
  }
  return;  
}

void write2D(TFile *f1, TFile *f2,  const char *sysName="", bool write=false){

  std::vector<TH2F*> H1 = getHistos2D(f1);
  std::vector<TH2F*> H2 = getHistos2D(f2);

  if(H1.size() != H2.size())
    ERROR( Form("Inconsistent file content: %i != %i",(int)H1.size(),(int)H2.size()) );

  for(unsigned int i(0); i<H1.size(); i++){

    TH2F *h = hDiff2D(H1.at(i),H2.at(i));
    setStyle2D(h);

    TString cname = Form("C_%s_%s",sysName,h->GetName());
    TCanvas *c  = new TCanvas(cname, cname, 1, 10, 770, 560);
    c->cd();

    TPad *p = new TPad(cname+"_p",cname+"_p", 0.00, 0.00, 1.00, 0.95, -1, 0, 0);
    p->SetLeftMargin(0.10);
    p->SetRightMargin(0.15);
    p->SetBottomMargin(0.15);
    p->SetTicks(1,1);
    p->Draw();
    p->cd();

    if( ((TString)h->GetName()).Contains("pt_eta") ){
      h->GetXaxis()->SetRangeUser(1., h->GetXaxis()->GetXmax());
      p->SetLogx();
    }
    
    TH2F* hRel = (TH2F*)h->Clone( Form("Rel_%s",h->GetName()) );
    hRel->Scale(100.);
    setStyle2D(hRel);
    hRel->Draw("COLZ TEXT45");
   
    drawGridLines2D(h);
    draw2DplotLabel(sysName);

    if(write){
      c->Print(Form("%s.eps",cname.Data()));
      writeToFile(h, sysName);
    }
  }
  return;
}

void writeSys(const char* file1, const char *file2, const char *sysName="", bool write=false){
  
  TFile *f1 = new TFile(file1);
  TFile *f2 = new TFile(file2);
  if(f1->IsZombie() || f2->IsZombie())
    ERROR( Form("Cannot find file %s|%s ",file1,file2) );

  TString fname1 = f1->GetName(), fname2 = f1->GetName();
  INFO( Form("Found files\n- %s\n- %s ", fname1.Data(), fname2.Data()) );

  if(fname1.Contains("1D") && fname2.Contains("1D")){
    INFO("Writing 1D histograms");
    write1D(f1, f2, sysName, write);
  }
  else if(fname1.Contains("2D") && fname2.Contains("2D")){
    INFO("Writing 2D histograms");
    write2D(f1, f2, sysName, write);
  }
  else
    ERROR("Mixing 1D and 2D histograms ...exiting");
}

std::vector<TH1F*> getTotalVar1D(std::vector<TH1F*> nominals, std::vector<TH1F*> variations, const char* suffix, bool showBins=false){
  std::vector<TH1F*> totalVars(0);
  if(!nominals.size()) return totalVars;

  for(auto nom : nominals){
    TString nameNom = nom->GetName();

    TH1F *hOut = (TH1F*)nom->Clone( Form("%s__%s", nameNom.Data(), suffix) );
    hOut->Reset();

    //Statistical errors 
    for(int i(1); i<=hOut->GetNbinsX(); i++) hOut->SetBinContent(i, nom->GetBinError(i));

    //Systematic errors
    for(auto var : variations){
      TString nameVar = var->GetName();

      if(!nameVar.Contains(nameNom)) continue;
      INFO( Form("Nominal(%s): Variation: %s",nameNom.Data(),nameVar.Data()) );

      for(int i(1); i<=hOut->GetNbinsX(); i++){
	if(showBins) INFO( Form("%s bin(%i) = %.3f \t variation %s = %.3f",nameNom.Data(),i,nom->GetBinContent(i),nameVar.Data(),var->GetBinContent(i)) );

	float errVar = var->GetBinContent(i)*nom->GetBinContent(i);
	float errTot = TMath::Power(errVar,2) + TMath::Power(hOut->GetBinContent(i),2);
	
	hOut->SetBinContent(i, TMath::Sqrt(errTot));
      }
    }
    //Total error <= nominal value
    for(int i(1); i<=hOut->GetNbinsX(); i++){ if(hOut->GetBinContent(i) > nom->GetBinContent(i)) hOut->SetBinContent(i,nom->GetBinContent(i)); } 

    if(!hOut->GetSumOfWeights()) continue;
    if(showBins){for(int i(1); i<=hOut->GetNbinsX(); i++){ INFO( Form("%s bin(%i) = %.3f",hOut->GetName(),i,hOut->GetBinContent(i)) ); } } 
    totalVars.push_back(hOut);
  }
  return totalVars;
}

std::vector<TH2F*> getTotalVar2D(std::vector<TH2F*> nominals, std::vector<TH2F*> variations, const char* suffix, bool showBins=false){
  std::vector<TH2F*> totalVars(0);
  if(!nominals.size()) return totalVars;

  for(auto nom : nominals){
    TString nameNom = nom->GetName();
    
    TH2F *hOut = (TH2F*)nom->Clone( Form("%s__%s", nameNom.Data(), suffix) );
    hOut->Reset();

    //Statistical errors
    for(int i(1); i<=hOut->GetNbinsX(); i++){ for(int j(1); j<=hOut->GetNbinsY(); j++){ hOut->SetBinContent(i, j, nom->GetBinError(i,j)); } }

    // Systematic errors
    for(auto var : variations){
      TString nameVar = var->GetName();

      if(!nameVar.Contains(nameNom)) continue;
      INFO( Form("Nominal(%s): Variation: %s",nameNom.Data(),nameVar.Data()) );

      for(int i(1); i<=hOut->GetNbinsX(); i++){
	for(int j(1); j<=hOut->GetNbinsY(); j++){
	  if(showBins) INFO( Form("%s bin(%i|%i) = %.3f \t variation %s = %.3f",nameNom.Data(),i,j,nom->GetBinContent(i,j),nameVar.Data(),var->GetBinContent(i,j)) );

	  float errVar = var->GetBinContent(i,j)*nom->GetBinContent(i,j);
	  float errTot = TMath::Power(errVar,2) + TMath::Power(hOut->GetBinContent(i,j),2);
	  
	  hOut->SetBinContent(i, j, TMath::Sqrt(errTot));
	}
      }
    }
    //Total error <= nominal value
    for(int i(1); i<=hOut->GetNbinsX(); i++){ for(int j(1); j<=hOut->GetNbinsY(); j++){ if(hOut->GetBinContent(i,j) > nom->GetBinContent(i,j)) hOut->SetBinContent(i,j,nom->GetBinContent(i,j)); } }

    if(!hOut->GetSumOfWeights()) continue;
    if(showBins){for(int i(1); i<=hOut->GetNbinsX(); i++){for(int j(1); j<=hOut->GetNbinsY(); j++){ INFO( Form("%s bin(%i|%i) = %.3f",hOut->GetName(),i,j,hOut->GetBinContent(i,j)) ); } } }
    totalVars.push_back(hOut);
  }
  return totalVars;
}

void makeTotalError1D(TFile *f, bool write, const char *suffix){

  std::vector<TH1F*> H = getHistos1D(f);
  std::vector<TH1F*> nominals(0), variations(0);

  for(auto h : H){
    TString name = h->GetName();
    name.Contains("__") ? variations.push_back(h) : nominals.push_back(h);
  }
  std::vector<TH1F*> hVarTotal = getTotalVar1D(nominals, variations, suffix);

  if(!write) return;
  for(auto var : hVarTotal){
    var->Write();
    INFO( Form("Wrote %s/%s", f->GetName(), var->GetName()) );
  }
}

void makeTotalError2D(TFile *f, bool write, const char *suffix){
  
  std::vector<TH2F*> H = getHistos2D(f);
  std::vector<TH2F*> nominals(0), variations(0);

  for(auto h : H){
    TString name = h->GetName();
    name.Contains("__") ? variations.push_back(h) : nominals.push_back(h);  
  }
  std::vector<TH2F*> hVarTotal = getTotalVar2D(nominals, variations, suffix);

  if(!write) return;
  for(auto var : hVarTotal){ 
    var->Write();
    INFO( Form("Wrote %s/%s", f->GetName(), var->GetName()) );
  }
}

void makeTotalError(const char* filename, bool write=false, const char *suffix="TOTAL"){
 
  TFile *f = TFile::Open(filename, "UPDATE");
  if(f->IsZombie())
    ERROR( Form("Cannot find file %s", filename) );

  if( ((TString)f->GetName()).Contains("1D") ){
    INFO("Writing 1D histograms");
    makeTotalError1D(f, write, suffix);
  }
  else if( ((TString)f->GetName()).Contains("2D") ){
    INFO("Writing 2D histograms");
    makeTotalError2D(f, write, suffix);
  }
  else
    ERROR("Mixing 1D and 2D histograms ...exiting");
}

