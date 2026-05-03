#include <iostream>
#include <vector>
#include <string>
#include "TStyle.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TString.h"
#include "TH1F.h"
#include "TKey.h"
#include "TCanvas.h"
#include "TLegend.h"

std::vector<TH1F*> getClassificationHists(TFile *f, TString directory, TString type){
  
  std::vector<TH1F*> hVec(0);
  if(!f){std::cout << "ERROR: No files selected" << std::endl; return hVec;}
  if(f->IsZombie()){std::cout << "ERROR: No files selected" << std::endl; return hVec;}

  TDirectory *d = (TDirectory*)f->Get(directory);
  if(!d){std::cout << "ERROR: No directory found" << std::endl; return hVec;}

  TH1F* hAll(0);
  if(type=="el") hAll = (TH1F*)d->Get("histoLoose_el0");
  if(type=="mu") hAll = (TH1F*)d->Get("histoLoose_mu0");

  TKey *key(0);
  TList* Objects = d->GetListOfKeys();
  Objects->Sort();

  TIter next(Objects);
  while(( key = (TKey*)next() )){
    TObject *obj = key->ReadObj();
    if (!obj->IsA()->InheritsFrom("TH1F")) continue;

    TH1F *hist = (TH1F*)obj;
    TString histname(hist->GetName());

    bool passClass = (histname.Contains("histoLoose_")
		      && !histname.Contains("all_")
		      && histname.Contains("0")
		      && histname != "histoLoose_el0"
		      && histname != "histoLoose_mu0");

    bool passType(false);
    if(type=="mu") passType = histname.Contains(type);
    if(type=="el") passType = histname.Contains(type) || histname.Contains("conversion") || histname.Contains("charge_flip"); 

    if(passClass && passType){ 
      hist->Scale(1./hAll->Integral());
      hVec.push_back(hist);
    }
  }
  return hVec;
}

std::vector<float> getYields(std::vector<TH1F*> vec){
  std::vector<float> yields(0);
  for(auto h : vec) yields.push_back(h->Integral());
  return yields;
}

std::vector<TString> getNames(std::vector<TH1F*> vec){
  std::vector<TString> names(0);
  for(auto h : vec){
    TString name = h->GetName();
    if(name.Contains("_HF"))            names.push_back("heavy flavor");
    if(name.Contains("_LF"))            names.push_back("light flavor");
    if(name.Contains("Muon_electron"))  names.push_back("muon-electron");
    if(name.Contains("Tau"))            names.push_back("tau decays");
    if(name.Contains("not_classified")) names.push_back("unclassified");
    if(name.Contains("prompt"))         names.push_back("prompt lepton");
    if(name.Contains("conversion"))     names.push_back("conversion");
    if(name.Contains("charge_flip"))    names.push_back("charge-flip");
  }
  return names;
}

std::vector<TH1F*> runFile(TString filename){
  std::vector<TH1F*> out(0);

  TFile *f = new TFile(filename);
  TString dirname  = "Efficiencies_Selection_1";
  
  std::vector<TH1F*> hClass_el = getClassificationHists(f, dirname, "el");
  std::vector<TH1F*> hClass_mu = getClassificationHists(f, dirname, "mu");

  std::vector<float> events_elClasses = getYields(hClass_el);
  std::vector<float> events_muClasses = getYields(hClass_mu);

  std::vector<TString> names_elClasses = getNames(hClass_el);
  std::vector<TString> names_muClasses = getNames(hClass_mu);

  for(unsigned int i(0); i<hClass_el.size(); i++) std::cout << Form("--> El Class :: %s \t fraction %.1f %s", names_elClasses.at(i).Data(), events_elClasses.at(i)*100, "%") << std::endl;
  for(unsigned int i(0); i<hClass_mu.size(); i++) std::cout << Form("--> Mu Class :: %s \t fraction %.1f %s", names_muClasses.at(i).Data(), events_muClasses.at(i)*100, "%") << std::endl;
  std::cout << std::endl;

  TH1F * hMu = new TH1F("hMu", "hMu", hClass_mu.size(), 0, hClass_mu.size());
  TH1F * hEl = new TH1F("hEl", "hEl", hClass_el.size(), 0, hClass_el.size());

  for(int bin(1); bin<=hMu->GetNbinsX(); bin++){
    hMu->GetXaxis()->SetBinLabel(bin, names_muClasses.at(bin-1));
    hMu->SetBinContent(bin, events_muClasses.at(bin-1));
  }
  for(int bin(1); bin<=hEl->GetNbinsX(); bin++){
    hEl->GetXaxis()->SetBinLabel(bin, names_elClasses.at(bin-1));
    hEl->SetBinContent(bin, events_elClasses.at(bin-1));
  }
  out.push_back(hMu);
  out.push_back(hEl);

  for(auto h : out){
    h->SetLineWidth(2);
    if(((TString)h->GetName()).Contains("Mu")) h->SetLineColor(kGreen-6);
    if(((TString)h->GetName()).Contains("El")) h->SetLineColor(kViolet-8);
    h->GetYaxis()->SetRangeUser(0.001,1);
    h->GetYaxis()->SetTitle("Fraction of events");
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetLabelOffset(0.01);
    h->GetYaxis()->SetLabelOffset(0.01);
  }
  return out;
}

void compareFiles(TString File1, TString File2, TString labelFile1, TString labelFile2, int print=0){

  std::cout << "Comparing files: \n-- " << File1 << "\n-- " << File2 << std::endl;
  std::cout << "Compare: [ " << labelFile1 << " || " << labelFile2 << " ]" << std::endl;
  std::vector<TH1F*> hFile1 = runFile(File1);
  std::vector<TH1F*> hFile2 = runFile(File2);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  for(unsigned int i(0); i<hFile1.size(); i++){

    TString cname("");
    if( ((TString)(hFile1.at(i))->GetName()).Contains("Mu") ) cname = "Muons";
    if( ((TString)(hFile1.at(i))->GetName()).Contains("El") ) cname = "Electrons";    

    TCanvas *c = new TCanvas(cname, cname);
    c->cd();
    c->SetGridx();
    c->SetTicks(1,1);
    (hFile1.at(i))->Draw("HIST");

    (hFile2.at(i))->SetLineStyle(7);
    (hFile2.at(i))->Draw("HIST SAME");

    TLegend *leg = new TLegend(0.52, 0.72, 0.82, 0.88);
    leg->SetHeader(cname);    
    leg->AddEntry(hFile1.at(i), labelFile1);
    leg->AddEntry(hFile2.at(i), labelFile2);
    leg->SetFillStyle(1111);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.04);
    leg->Draw("SAME");
  
    if(print) 
      c->Print(Form("%s.pdf",c->GetName())); 
  }
  return;
}

