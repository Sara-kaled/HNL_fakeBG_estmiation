#include "RatePlotter.h"

void RatePlotter::setDebug(bool debug){
  Debug = debug;
  INFO("setDebug", Form("Debug %i",Debug));
  return;
}

void RatePlotter::setPrint(bool print){
  Print = print;
  INFO("setPrint", Form("Print %i",Print));
}

void RatePlotter::setOutDir(std::string dir){
  outDir = dir;
  INFO("setOutDir", Form("Output dir: %s",outDir.c_str()));
}

void RatePlotter::setEffDirectory(const char* dir){
  effDir = dir;
  INFO("setEffDirectory", Form("Efficiency dir: %s",effDir));
}

void RatePlotter::writeHistFile(const char* outname, bool writeH){
  writeHist = writeH;
  outFile   = outname;
  INFO("writeHistFile", Form("Write=%i, Outfile name: %s",writeHist, outFile));
}

void RatePlotter::setStylePath(const char* path){
  stylePath = path;
  INFO("setStylePath", Form("Style path: %s",path));
}

void RatePlotter::setStyle(bool setAtlas){
  Style = true;
  INFO("setStyle", Form("AtlasStyle %i",setAtlas));
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(000000);
  if(!setAtlas) return;
  gROOT->LoadMacro(stylePath);
  gROOT->ProcessLine("SetAtlasStyle()");
  return;
}

TFile *RatePlotter::findFile(const char* name){
  for(auto f : InFiles){
    if( (TString)name == (TString)f->GetName() ) return f;
  }
  return nullptr;
}

void RatePlotter::setRateType(const char *type){
  RateType = type;
  if(!(RateType=="Fake" || RateType=="Real")){    
    ERROR("setRateType", "No valid rate type selected. Please set rate type [Fake|Real]");
  }
  INFO("setRateType", Form("Rate type: %s", RateType.Data()));
}

void RatePlotter::setLabel(std::string label, std::string option){
  if(option=="MC"){
    INFO("setMCLabel", Form("MC label: %s",label.c_str()));
    mcLabel = label;
    return;
  }
  else if(option=="Data"){
    INFO("setDataLabel", Form("DATA label: %s",label.c_str()));
    dataLabel = label;
    return;
  }
  else{ INFO("setLabel", "no label set"); }
  return;
}

const char* RatePlotter::addSuffix(const char* name){
  if(!sysSuffix.length()) return name;
  return Form("%s__%s", name, sysSuffix.c_str());
}

void RatePlotter::setHistRange(float min, float max){
  yMin = min; yMax = max;
  INFO("setHistRange", Form("Set y range to [%.2f|%.2f]",yMin,yMax));
}

void RatePlotter::setRatioRange(float min, float max){
  yRMin = min; yRMax = max;
  INFO("setRatioRange", Form("Set y(ratio) range to [%.2f|%.2f]",yRMin,yRMax));
}

void RatePlotter::getFiles(const char* dirname, const char* key1, const char* key2){

  TString command = Form("ls -1 %s/*.root | sort &> filelist.txt",dirname);
  if(strlen(key1) && !strlen(key2)) command = Form("ls -1 %s/*.root | sort | grep %s &> filelist.txt",dirname,key1);
  if(strlen(key1) &&  strlen(key2)) command = Form("ls -1 %s/*.root | sort | grep %s | grep %s &> filelist.txt",dirname,key1,key2);
  gSystem->Exec(command);

  std::ifstream File;
  File.open("filelist.txt");
  if( !File.is_open() ){ INFO("addMCFiles", "No input list found"); return;}
  while(true){
    if( File.eof() ) break;
    TString n;
    File >> n;
    if(!n.Contains(".root")) continue;
    bool isDataFile = n.Contains("data") || n.Contains("Data");

    if(isDataFile) this->addDataFile(n.Data());
    else this->addMCFile(n.Data());
  }
  gSystem->Exec("rm filelist.txt");
}

void RatePlotter::addMCFile(const char* filename){
  DEBUG("addFile", Form("Adding MC file %s",filename));
  MCFiles.push_back(filename);
  MCRates = true;
  return;
}

void RatePlotter::addDataFile(const char* filename){
  DEBUG("addFile", Form("Adding data file %s",filename));
  DataFiles.push_back(filename);
  DataRates = true;
  return;
}

void RatePlotter::addPromptFile(const char* filename){
  DEBUG("addFile", Form("Adding prompt MC file %s",filename));
  PromptMCFiles.push_back(filename);
  Prompt = true;
  return;
}

void RatePlotter::setProcessSubtraction(TString proc, float procSF){
  INFO("setProcSubtr", Form("MC process to subtract: %s",proc.Data()));
  subtractedProc.push_back(proc);
  subtractedProcSF.push_back(procSF);
}

void RatePlotter::setLumi(float lumi){
  Lumi = lumi;
  INFO("setLumiScale", Form("Set lumi scale to %.1f", Lumi));
}

std::vector<TGraphAsymmErrors*> RatePlotter::vec(TGraphAsymmErrors *g1, TGraphAsymmErrors *g2){
  std::vector<TGraphAsymmErrors*> v(0);
  if(g1) v.push_back(g1);
  if(g2) v.push_back(g2);
  return v;
}

void RatePlotter::lumiScale(std::vector<TH1F*> histos){
  if(histos.empty()) return; 
  for(auto h : histos) h->Scale(Lumi);
  INFO("luminosityScale",Form("Scale histograms by %.1f",Lumi));
}

std::vector<TH1F*> RatePlotter::getHistosFromList(std::vector<std::string> filelist, const char* dirname){ 
  if(filelist.empty()){ ERROR("getHistos", "No files selected"); }
  
  INFO("getHistos", Form("Retrieving histograms from %i files", (int)filelist.size()));
  const char* filename = filelist.at(0).c_str();

  std::vector<TH1F*> histos = this->getHistos(filename, dirname);
  if(filelist.size()==1) return histos; 
  
  for(unsigned int i(1); i<filelist.size(); i++){

    const char* filename = filelist.at(i).c_str();
    std::vector<TH1F*> histVec = this->getHistos(filename, dirname);

    for(unsigned int j(0); j<histVec.size(); j++) histos[j]->Add(histVec[j]); 
  }
  return histos;
}

std::vector<TH1F*> RatePlotter::getHistos(const char* filename, const char* dirname){
  TDirectory *d(0);
  TFile *file = this->findFile(filename);
  
  if(!file){ 
    file = new TFile(filename); 
    if(file->IsZombie()){ ERROR("getHistos", Form("Failed to open: %s", filename));}

    d = (TDirectory*)file->Get(dirname);
    if(d->IsZombie()){ ERROR("getHistos", Form("Failed to open: %s/%s", filename,dirname));}

    INFO("getHistos", Form("Opened: %s/%s", file->GetName(), d->GetName()));
    InFiles.push_back(file);
  }
  if(!d) d = (TDirectory*)file->Get(dirname);

  float scale = getMCNorm(file);
  std::vector<TH1F*> hVec(0);

  TKey *key(0);
  TList* Objects = d->GetListOfKeys();
  Objects->Sort();

  TIter next(Objects);
  while(( key = (TKey*)next() )){
    TObject *obj = key->ReadObj();                                                
    if (!obj->IsA()->InheritsFrom("TH1F")) continue;
    TH1F *hist = (TH1F*)obj;
    hist->Scale(scale);
    hVec.push_back(hist);
  }
  if(!FakeSourcesEl.empty()) this->addFakeHist(hVec, "El");
  if(!FakeSourcesMu.empty()) this->addFakeHist(hVec, "Mu");

  DEBUG("getHistos",Form("Retrieved %i histograms from file",(int)hVec.size()));
  return hVec;
}

TH1F* RatePlotter::findHisto(TString name, std::vector<TH1F*> histos){
  for (auto h : histos)
    if(h->GetName()==name) return h;
  INFO("findHisto", Form("Histogram %s not found", name.Data()));
  return nullptr;
}

float RatePlotter::getProcessSF(TString proc){
  for(unsigned int i(0); i<subtractedProc.size(); i++){ if(subtractedProc[i] == proc) return subtractedProcSF[i]; }
  return 1.;
}

float RatePlotter::getMCNorm(TFile *f){

  TH1F *hNorm = (TH1F*)f->Get("MCLumiHist");
  if(!hNorm){ ERROR("getMCNorm", Form("No normalization histogram found in file %s", f->GetName()));}

  float mcLumi = hNorm->GetBinContent(1);
  if(mcLumi > 0.){
    DEBUG("getMCNorm", Form("File %s : MC virtual lumi %.1f", f->GetName(), mcLumi)); 
    return 1./mcLumi;
  } 
  return 1.;
}

TString RatePlotter::GetXTitle(TH1F* h){
  if(!h) return "";
  TString name = Form("%s",h->GetName());
  if( name.Contains("el0") || name.Contains("electron0") ) return "Electron #it{p}_{T} [GeV]";
  if( name.Contains("el1") || name.Contains("electron1") ) return "Electron |#it{#eta}|";
  if( name.Contains("mu0") || name.Contains("muon0")     ) return "Muon #it{p}_{T} [GeV]";
  if( name.Contains("mu1") || name.Contains("muon1")     ) return "Muon |#it{#eta}|";
  
  if( name.Contains("charge_flip0")|| name.Contains("conversion0") ) return "Electron #it{p}_{T} [GeV]";
  if( name.Contains("charge_flip1")|| name.Contains("conversion1") ) return "Electron |#it{#eta}|";
  return "";  
}

const char* RatePlotter::getAxisPar(TString name){
  if( name.Contains("{p}_{T}") )  return "pt";
  if( name.Contains("{#eta}")   ) return "eta";
  return "";
}

TGraphAsymmErrors* RatePlotter::getRateGraph(TH1F *hPass, TH1F *hTotal, TString source){
  TGraphAsymmErrors *g(0);
  if( !checkEntries(hPass,hTotal) ) return g;

  g = new TGraphAsymmErrors(hPass, hTotal, "n");
  DEBUG("getRateGraph", Form("Created graph (%s): %i points (pass=%s, tot=%s)",source.Data(),(int)g->GetN(),hPass->GetName(),hTotal->GetName()));

  g->SetLineWidth(2);
  g->SetMarkerStyle(20);
  g->SetMarkerSize(0.8);
  if(source=="Data")
    g->SetLineColor(kBlack);
  else
    this->getMCcolor(g,source);
  g->SetMarkerColor(g->GetLineColor());

  return g;
}

TH1F* RatePlotter::divideTH1(TH1F* hPass, TH1F *hTotal){
  TH1F *h(0);
  if( !checkEntries(hPass,hTotal) ) return h;

  h = (TH1F*)hTotal->Clone(Form("%s_%s",hPass->GetName(),hTotal->GetName()));
  h->Reset();

  for(int i(1); i<=h->GetNbinsX(); i++){
    float pass  = hPass->GetBinContent(i);
    float total = hTotal->GetBinContent(i);

    float val = total > 0. ? pass/total : 0.;
    float err = 0.5*(TEfficiency::Normal(total,pass,0.68,1) - TEfficiency::Normal(total,pass,0.68,0));
    if(val<=0.){ 
      val = hPass->Integral() / hTotal->Integral(); 
      err = val; 
    }
    DEBUG("divideTH1", Form("Bin (%i): N(pass)=%.3f, N(tot)=%.3f \t Rate=%.5f (err=%.5f)", i, pass, total, val, err));
    h->SetBinContent(i, val);
    h->SetBinError(i, err);
  }
  return h;
}

TH2F* RatePlotter::divideTH2(TH2F* hPass, TH2F *hTotal){
  TH2F *h(0);
  if( !checkEntries(hPass,hTotal) ) return h;

  h = (TH2F*)hTotal->Clone(Form("%s_over_%s",hPass->GetName(),hTotal->GetName()));
  h->Reset();

  for(int x(1); x<=h->GetNbinsX(); x++){
    for(int y(1); y<=h->GetNbinsY(); y++){

      float pass  = hPass->GetBinContent(x,y);
      float total = hTotal->GetBinContent(x,y); 

      float val = total > 0. ? pass/total : 0.;
      float err = 0.5*(TEfficiency::Normal(total,pass,0.68,1) - TEfficiency::Normal(total,pass,0.68,0));
      if(val<=0.){
	val = hPass->Integral(x-1,x,0,hPass->GetNbinsY()) /  hTotal->Integral(x-1,x,0,hTotal->GetNbinsY());
	err = val;
      }
      DEBUG("divideTH2", Form("Bin (%i|%i): N(pass)=%.3f, N(tot)=%.3f \t Rate=%.5f (err=%.5f)", x, y, pass, total, val, err));
      h->SetBinContent(x,y,val);
      h->SetBinError(x,y,err);
    }
  }
  return h;
}

void RatePlotter::setHistStyle(TH1F* h){
  if(!h) return;
  h->GetYaxis()->SetRangeUser(yMin, yMax);
  h->GetYaxis()->SetTitle("Efficiency");
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetYaxis()->SetTitleOffset(1.);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetLabelOffset(0.01);

  h->GetXaxis()->SetTitle(GetXTitle(h));
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetTitleOffset(1.50);
  h->GetXaxis()->SetLabelOffset(0.02);
  if(isAllHist(h))
    h->GetXaxis()->SetLabelSize(0.00);
  else
    h->GetXaxis()->SetLabelSize(0.04);
  return;
}

void RatePlotter::setHistStyle(TH2F *h, float zMax){
  if(!h) return;
  gStyle->SetPaintTextFormat(".1f");
  gStyle->SetNumberContours(40);
  gStyle->SetPalette(104);

  h->GetZaxis()->SetRangeUser(0, zMax);
  h->GetZaxis()->SetLabelSize(0.04);
  h->GetZaxis()->SetTitle("Efficiency [%]");
  h->GetZaxis()->SetTitleSize(0.04);
  h->GetZaxis()->SetTitleOffset(1.0);
  
  h->GetXaxis()->SetTitle("Lepton #it{p}_{T} [GeV]");
  if( ((TString)h->GetName()).Contains("_mu")) h->GetXaxis()->SetTitle("Muon #it{p}_{T} [GeV]");
  if( ((TString)h->GetName()).Contains("_el")) h->GetXaxis()->SetTitle("Electron #it{p}_{T} [GeV]");

  h->GetYaxis()->SetTitle("Lepton |#eta|");
  if( ((TString)h->GetName()).Contains("_mu")) h->GetYaxis()->SetTitle("Muon |#it{#eta}|");
  if( ((TString)h->GetName()).Contains("_el")) h->GetYaxis()->SetTitle("Electron |#it{#eta}|");

  h->GetYaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleOffset(1.15);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelOffset(0.01);

  h->GetXaxis()->SetTitleSize(0.04);
  h->GetXaxis()->SetTitleOffset(1.55);
  h->GetXaxis()->SetLabelOffset(0.01);
  h->GetXaxis()->SetLabelSize(0.04);

  h->SetMarkerSize(1.4);
  h->SetMarkerColor(kWhite);
  return;
}

void RatePlotter::setRatioHistStyle(TH1F* h){
  if(!h) return;
  h->GetYaxis()->SetRangeUser(yRMin, yRMax);
  h->GetYaxis()->SetTitle("Ratio");
  h->GetYaxis()->SetTitleSize(0.12);
  h->GetYaxis()->SetTitleOffset(0.4);
  h->GetYaxis()->SetLabelSize(0.12);
  h->GetYaxis()->SetLabelOffset(0.011);
  h->GetYaxis()->SetNdivisions(6);

  h->GetXaxis()->SetTitleSize(0.12);
  h->GetXaxis()->SetTitleOffset(1.45);
  h->GetXaxis()->SetLabelOffset(0.02);
  if(isAllHist(h))
    h->GetXaxis()->SetLabelSize(0.00);
  else
    h->GetXaxis()->SetLabelSize(0.12);
  return;
}

void RatePlotter::setHistStyleNoRatio(TH1F *h){
  if(!h) return;
  h->SetLineWidth(0);
  h->SetLineColor(kWhite);
  h->GetYaxis()->SetTitle("Events");
  h->GetXaxis()->SetTitle(GetXTitle(h));

  h->GetYaxis()->SetTitleOffset(1.35);
  h->GetXaxis()->SetTitleOffset(1.6);
  
  h->GetXaxis()->SetLabelSize(0.04);
  h->GetYaxis()->SetLabelSize(0.04);
  h->GetXaxis()->SetTitleSize(0.04);
  h->GetYaxis()->SetTitleSize(0.04);

  h->GetXaxis()->SetLabelOffset(0.010);
  h->GetYaxis()->SetLabelOffset(0.010);
  return;
}

void RatePlotter::getMCcolor(TGraphAsymmErrors *g, TString col){
  
  if(col == "MC_blue"){   g->SetLineColor(kBlue-7);    return; }
  if(col == "MC_green"){  g->SetLineColor(kGreen-6);   return; }
  if(col == "MC_red"){    g->SetLineColor(kRed+3);     return; }
  if(col == "MC_violet"){ g->SetLineColor(kMagenta-7); return; }
  if(col == "MC_orange"){ g->SetLineColor(kOrange-9);  return; }
  if(col == "MC_cyan"){   g->SetLineColor(kCyan-3);    return; }

  g->SetLineColor(kBlue-7);
  return;
}

void RatePlotter::drawGridLines(TH1F* h, float yEnd){
  if(!h) return;
  for(int i(2); i<=h->GetNbinsX(); i++){
    float x = h->GetXaxis()->GetBinLowEdge(i);
        
    TLine xLine;
    xLine.SetLineWidth(1);
    xLine.SetLineColor(kGray+1);
    xLine.SetLineStyle(7);
    xLine.DrawLine(x, yEnd, x, h->GetMaximum());
  }
}

void RatePlotter::draw2DGridLines(TH2F* h){
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

void RatePlotter::drawEtaRegions(TH1F* h, float yEnd, bool binLabels, int etaBins){
  if(!h) return;
  for(int i(1); i<=h->GetNbinsX(); i++){
    float x = h->GetXaxis()->GetBinUpEdge(i);
    if(i % etaBins) continue;

    TLine xLine;
    xLine.SetLineWidth(1);
    xLine.SetLineColor(kGray+1);
    xLine.SetLineStyle(7);
    xLine.DrawLine(x, yEnd, x, h->GetMaximum());
    
    TLatex p;
    p.SetNDC();
    p.SetTextColor(kGray+3);;
    p.SetTextFont(62);
    p.SetTextSize(0.12);
    if(binLabels) 
      p.DrawLatex(0.165*i/etaBins, 0.18, Form("#eta bin %i",i/etaBins) );
  }
}

void RatePlotter::drawLeptonFlavor(TH1F *h, bool no2pad){
  if(!h) return;
  std::string flavor("");

  TString name = Form("%s",h->GetName());
  if( name.Contains("_mu") ) flavor = "Muons";
  if( name.Contains("_el") ) flavor = "Electrons";
  if( name.Contains("charge_flip")|| name.Contains("conversion") ) flavor = "Electrons"; 
  if(!flavor.length()) return;

  float x1 = 0.18;
  float x2 = AtlasLabel ? 0.44 : 0.32;
  float y1 = AtlasLabel ? 0.71 : 0.76;
  float y2 = AtlasLabel ? 0.89 : 0.88;
  if(no2pad){y1 += 0.025; y2 -= 0.015; x2 -= 0.02;}

  if(!AtlasLabel && flavor == "Electrons") x2 += 0.03;

  TLegend *l = new TLegend(x1, y1, x2, y2);
  l->SetFillStyle(1111);
  l->SetFillColor(kWhite);
  l->SetBorderSize(1);
  l->Draw("SAME");

  TLatex n;
  n.SetNDC();
  n.SetTextColor(kBlack);;
  n.SetTextFont(62);
  n.SetTextSize(no2pad ? 0.04 : 0.06);
  if(no2pad) n.DrawLatex(0.20, AtlasLabel ? 0.76 : 0.81, flavor.c_str());
  else n.DrawLatex(0.20, AtlasLabel ? 0.73 : 0.80, flavor.c_str());

  if(AtlasLabel){
    n.SetTextFont(72);
    n.DrawLatex(0.20, 0.82, "ATLAS");
    n.SetTextFont(42); 
    n.DrawLatex(0.31, 0.82, "Internal"); 
  }
  return;
}

void RatePlotter::draw2DplotLabel(TH2F* h, TString type, TString source){
  
  const char *flavor("");
  TString name = Form("%s",h->GetName());
  if( name.Contains("_mu") ) flavor = "Muon";
  if( name.Contains("_el") ) flavor = "Electron";
  if( name.Contains("charge_flip")|| name.Contains("conversion") ) flavor = "Electron";
  type.ToLower();

  TLatex n;
  n.SetNDC();
  n.SetTextColor(kBlack);;
  n.SetTextSize(0.04);
  if(AtlasLabel){
    n.SetTextFont(72); n.DrawLatex(0.12,0.97,"ATLAS");
    n.SetTextFont(42); n.DrawLatex(0.22,0.97,"Internal");
  }
  n.SetTextSize(0.04);
  n.SetTextFont(42);
  n.DrawLatex(AtlasLabel ? 0.33 : 0.12, 0.97, Form("%s %s efficiencies (%s)", flavor, type.Data(), source.Data()) );
}

TLegend* RatePlotter::makeLegend(TGraphAsymmErrors *g1, TGraphAsymmErrors *g2, TString type1, TString type2){
  int entries(1);
  if(g1 && g2) entries++;

  TLegend *leg(0);
  switch(entries){
  case 1: leg = new TLegend(0.58, 0.80, 0.88, 0.89); break;
  case 2: leg = new TLegend(0.58, 0.74, 0.88, 0.89); break;
  default: return leg;
  }
  leg->SetNColumns(1); 
  leg->SetBorderSize(1);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetFillStyle(1111);
  leg->SetFillColor(kWhite);
  if(g1) leg->AddEntry(g1, Form(" Rates %s", type1.Data()) );
  if(g2) leg->AddEntry(g2, Form(" Rates %s", type2.Data()) );  
  return leg;
}

void RatePlotter::drawRatio(std::vector< TGraphAsymmErrors* > graphs, TH1F* h){

  if(graphs.empty()){ INFO("drawRatio", "No input graphs provided"); return;} 
  TGraphAsymmErrors *gNom = graphs.at(0); 

  std::vector<TGraphAsymmErrors*> VR(0);
  for(auto g : graphs){
    if(graphs.size()>1 && g==gNom) continue;

    TGraphAsymmErrors *gR = new TGraphAsymmErrors();
    double x(0),xN(0);
    double y(0),yN(0);

    for(int i(0); i<g->GetN(); i++){
      gNom->GetPoint(i,xN,yN);
      g->GetPoint(i,x,y);
      float ratio = (y>0 && yN>0) ? y/yN : 0;

      DEBUG("drawRatio", Form("Point %i, x=%.1f, ratio=%.2f", i, x, ratio));
      gR->SetPoint(i,x,ratio);
      
      gR->SetPointEXhigh(i, g->GetErrorXhigh(i));
      gR->SetPointEXlow(i,  g->GetErrorXlow(i));
      gR->SetPointEYhigh(i, g->GetErrorYhigh(i)/yN);
      gR->SetPointEYlow(i,  g->GetErrorYlow(i)/yN);

      gR->SetLineWidth(2);
      gR->SetMarkerStyle(20);
      gR->SetMarkerSize(0.8);
      gR->SetLineColor(g->GetLineColor());
      gR->SetMarkerColor(gR->GetLineColor());
    }
    VR.push_back(gR);
  }

  TH1F *temp = (TH1F*)h->Clone(Form("Template_%s",h->GetName()));
  setRatioHistStyle(temp);
  temp->Draw("AXIS");

  if(!isAllHist(temp))
    this->drawGridLines(temp, yRMin);
  else
    this->drawEtaRegions(temp, yRMin, true);

  TLine *l = new TLine();
  l->SetLineColor(kGray+1);
  l->SetLineStyle(9);
  l->SetLineWidth(1);
  l->DrawLine(temp->GetXaxis()->GetXmin(), 1., temp->GetXaxis()->GetXmax(), 1.);

  for(auto g : VR) g->Draw("SAME P0");
}

void RatePlotter::makeRatePlot(TString namePass, TString nameTot){
  
  std::cout << std::endl;
  if(!MCRates && !DataRates){ INFO("makeRatePlot", "No input (MC or Data) provided"); return;}

  INFO("makeRatePlot", Form("[MC|Data] = [%i|%i]",(int)MCRates,(int)DataRates));
  
  if(MCRates)   histosMC   = getHistosFromList(MCFiles, effDir);
  if(DataRates) histosData = getHistosFromList(DataFiles, effDir);
  this->lumiScale(histosMC);

  if(Prompt && DataRates){
    INFO("makeRatePlot", Form("Subtracting prompt processes from %i files", (int)PromptMCFiles.size()));
    std::vector<TH1F*> histosPrompt = getHistosFromList(PromptMCFiles, effDir);
   
    this->lumiScale(histosPrompt);
    this->subtractPrompt(histosData, histosPrompt);
  }
  if(!namePass.Length() || !nameTot.Length()){ INFO("makeRatePlot", "No input names provided"); return;}

  INFO("makeRatePlot", Form("Calculating rates from %s over %s",namePass.Data(),nameTot.Data()));
  
  TGraphAsymmErrors *gMC(0), *gData(0);
  TH1F *hMC(0), *hData(0), *h1_MC(0), *h2_MC(0), *h1_Data(0), *h2_Data(0);

  if(MCRates){
    h1_MC = findHisto(nameTot, histosMC);
    h2_MC = findHisto(namePass,histosMC);
    this->subtractMCProcess(h1_MC, h2_MC, histosMC);

    hMC = divideTH1(h2_MC, h1_MC);
    gMC = getRateGraph(h2_MC, h1_MC, "MC_blue");
  }

  if(DataRates){
    h1_Data = findHisto(nameTot, histosData);
    h2_Data = findHisto(namePass,histosData);
    this->subtractMCProcess(h1_Data, h2_Data, histosMC);

    hData = divideTH1(h2_Data, h1_Data);
    gData = getRateGraph(h2_Data, h1_Data, "Data");
  }
  
  if(!Style) this->setStyle(1);  
  TString cname = Form("%s_over_%s",namePass.Data(),nameTot.Data());
  cname = addSuffix(cname.Data());
  TCanvas *c  = new TCanvas(cname, cname, 1, 10, 770, 560);

  TPad *p1 = new TPad(cname+"_p1",cname+"_p1", 0.00, 0.30, 1.00, 1.00, -1, 0, 0);
  TPad *p2 = new TPad(cname+"_p2",cname+"_p1", 0.00, 0.00, 1.00, 0.30, -1, 0, 0);
  c->cd();
  p2->SetLeftMargin(0.10);
  p1->SetLeftMargin(0.10);
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.4);  
  p1->SetBottomMargin(0.01);
  p2->Draw();
  p1->Draw();

  p1->cd();
  TH1F* hTemp = h1_Data ? h1_Data : h1_MC;
  hTemp->Reset();
  setHistStyle(hTemp);
  if( ((TString)hTemp->GetXaxis()->GetTitle()).Contains("{p}_{T}") ){
    hTemp->GetXaxis()->SetRangeUser(1., hTemp->GetXaxis()->GetXmax());
    p1->SetLogx();
    p2->SetLogx();
  };
  hTemp->Draw("AXIS");
  
  if(gMC) gMC->Draw("SAME P0");
  if(gData) gData->Draw("SAME P0");

  if(!isAllHist(hTemp))
    this->drawGridLines(hTemp, yMin);
  else
    this->drawEtaRegions(hTemp, yMin);
  
  TLegend *leg = this->makeLegend(gMC, gData, mcLabel.c_str(), dataLabel.c_str());
  leg->Draw("SAME");
  this->drawLeptonFlavor(hTemp);

  p2->cd();
  this->drawRatio(vec(gMC,gData), hTemp);

  if(hMC && writeHist)
    this->writeToFile(hMC, RateType, "MC", outFile);

  if(hData && writeHist)
    this->writeToFile(hData, RateType, "Data", outFile);


  if(Print){
    if(!(bool)outDir.length()) outDir = ".";
    this->checkDir(outDir);
    c->Print(Form("%s/%s.%s",outDir.c_str(),c->GetName(),figType));
    INFO("makeRatePlot", Form("Created %s/%s.%s",outDir.c_str(),c->GetName(),figType));
  }
}

void RatePlotter::makeRatePlot2D(TString namePass, TString nameTot, TString source){
  
  std::cout << std::endl;
  if(!source.Length() || !RateType.Length()){ INFO("makeRatePlot2D", "Please select rate type [Fake|Real] and source [Data|MC]"); }

  if( (source=="MC" && MCFiles.empty()) || (source=="Data" && DataFiles.empty()) ){ INFO("makeRatePlot2D", "No data or MC files found"); return; }

  std::vector<std::string> files(0);
  if(source=="MC")   files = MCFiles;
  if(source=="Data") files = DataFiles;
  if(files.empty()){ INFO("makeRatePlot2D", "No source [Data|MC] selected"); return; }

  TFile *f0 = new TFile(files.front().c_str());
  TH2F *hPass = (TH2F*)f0->Get( Form("%s/%s",effDir,namePass.Data()) );
  TH2F *hTot  = (TH2F*)f0->Get( Form("%s/%s",effDir,nameTot.Data())  );

  if(!hPass || !hTot){ INFO("makeRatePlot2D", "No matching histogram found"); return; }
  hPass->Reset();
  hTot->Reset();

  TH2F *hRate(0), *hRatePerc(0);
  for(auto file : files){
    TFile *f = new TFile(file.c_str());
    float scale = (source=="MC") ? getMCNorm(f) : 1.;
    
    TH2F *hp = (TH2F*)f->Get( Form("%s/%s",effDir,namePass.Data()) );
    TH2F *ht = (TH2F*)f->Get( Form("%s/%s",effDir,nameTot.Data()) );
    hp->Scale(scale); 
    ht->Scale(scale);
    
    hPass->Add(hp);
    hTot->Add(ht);

    INFO("makeRatePlot2D", Form("File %s/%s : Calculating rates (%s) from %s over %s",f->GetName(),effDir,source.Data(),namePass.Data(),nameTot.Data()));
  }
  if(source=="MC"){ hPass->Scale(Lumi); hTot->Scale(Lumi); }

  if(Prompt && source=="Data"){
    INFO("makeRatePlot2D", Form("Subtracting prompt processes from %i files", (int)PromptMCFiles.size()));

    for(auto promptFile : PromptMCFiles){
      TFile *f = new TFile(promptFile.c_str());
      float scale = getMCNorm(f);
      INFO("makeRatePlot2D", Form("Subtract prompt rates from file %s (MC norm: 1/%.1f)",f->GetName(),1./scale));

      TH2F *hp_prompt = (TH2F*)f->Get( Form("%s/%s",effDir,namePass.Data()) );
      TH2F *ht_prompt = (TH2F*)f->Get( Form("%s/%s",effDir,nameTot.Data()) );
      hp_prompt->Scale(Lumi*scale);
      ht_prompt->Scale(Lumi*scale);

      this->subtract(hPass, hp_prompt);
      this->subtract(hTot,  ht_prompt);
    }
  }
  this->subtractMCProcess2D(hTot, hPass, MCFiles);
  hRate = divideTH2(hPass, hTot);

  //Same plot in % for drawing
  hRatePerc = divideTH2(hPass, hTot);
  hRatePerc->Scale(100.);

  if(!Style) this->setStyle(1);
  setHistStyle(hRate, 0.99);
  setHistStyle(hRatePerc, 99.);

  TString cname = Form("%s_over_%s_%s_%s",namePass.Data(),nameTot.Data(),RateType.Data(),source.Data());
  cname = addSuffix(cname.Data());
  TCanvas *c  = new TCanvas(cname, cname, 1, 10, 770, 560);
  c->cd();

  TPad *p = new TPad(cname+"_p",cname+"_p", 0.00, 0.00, 1.00, 0.95, -1, 0, 0);
  p->SetLeftMargin(0.10);
  p->SetRightMargin(0.15);
  p->SetBottomMargin(0.15);
  p->Draw();
  p->cd();
  if( ((TString)hRatePerc->GetXaxis()->GetTitle()).Contains("{p}_{T}") ){
    hRate->GetXaxis()->SetRangeUser(1., hRate->GetXaxis()->GetXmax());
    hRatePerc->GetXaxis()->SetRangeUser(1., hRatePerc->GetXaxis()->GetXmax());
    p->SetLogx();
  };
  hRatePerc->Draw("COLZ TEXT45 E");
  
  this->draw2DGridLines(hRatePerc);
  this->draw2DplotLabel(hRatePerc, RateType, source);

  if(writeHist) 
    this->writeToFile(hRate, RateType, source, outFile);

  if(Print){
    if(!(bool)outDir.length()) outDir = ".";
    this->checkDir(outDir);
    c->Print(Form("%s/%s.%s",outDir.c_str(),c->GetName(),figType));
    INFO("makeRatePlot2D", Form("Created %s/%s.%s",outDir.c_str(),c->GetName(),figType));
  }
}

TString RatePlotter::getOriginLabel(TString name){
  if(name.Contains("Fakes"))              return "All fakes";
  else if(name.Contains("HF"))            return "Heavy flavor";
  else if(name.Contains("LF"))            return "Light flavor";
  else if(name.Contains("conversion"))    return "Conversion";
  else if(name.Contains("prompt"))        return "Prompt";
  else if(name.Contains("charge_flip"))   return "Charge flip";
  else if(name.Contains("Muon_electron")) return "Muon elec.";
  else if(name.Contains("Tau_electron"))  return "Tau elec.";
  else if(name.Contains("Tau_muon"))      return "Tau muon";
  else if(name.Contains("not_classified"))return "Unclassified";
  else return "";
}

void RatePlotter::compareMCRates(TString namePass1, TString nameTot1, TString namePass2, TString nameTot2, TString col1, TString col2){
  
  std::cout << std::endl;
  if(!MCRates){ INFO("compareMCRates", "No MC input provided"); return;}

  histosMC = getHistosFromList(MCFiles, effDir);
  this->lumiScale(histosMC);  

  if(!namePass1.Length() || !namePass2.Length() || !nameTot1.Length()  || !nameTot2.Length()){ INFO("compareMCRates", "No input names provided"); return; }

  INFO("compareMCRates", Form("Calculating rates from %s over %s",namePass1.Data(),nameTot1.Data()));
  INFO("compareMCRates", Form("Calculating rates from %s over %s",namePass2.Data(),nameTot2.Data()));
  
  TGraphAsymmErrors *gMC1(0), *gMC2(0);
  TH1F *h1_MC1(0), *h2_MC1(0), *h1_MC2(0), *h2_MC2(0);

  h1_MC1 = findHisto(nameTot1, histosMC);
  h2_MC1 = findHisto(namePass1,histosMC);

  h1_MC2 = findHisto(nameTot2, histosMC);
  h2_MC2 = findHisto(namePass2,histosMC);
 
  gMC1 = getRateGraph(h2_MC1, h1_MC1, Form("MC_%s",col1.Data()));
  gMC2 = getRateGraph(h2_MC2, h1_MC2, Form("MC_%s",col2.Data()));
  
  if(!Style) this->setStyle(1);  
  TString cname = Form("%s_over_%s_AND_%s_over_%s",namePass1.Data(),nameTot1.Data(),namePass2.Data(),nameTot2.Data());
  cname = addSuffix(cname.Data());  
  TCanvas *c  = new TCanvas(cname, cname, 1, 10, 770, 560);

  TPad *p1 = new TPad(cname+"_p1",cname+"_p1", 0.00, 0.30, 1.00, 1.00, -1, 0, 0);
  TPad *p2 = new TPad(cname+"_p2",cname+"_p1", 0.00, 0.00, 1.00, 0.30, -1, 0, 0);
  c->cd();
  p2->SetLeftMargin(0.10);
  p1->SetLeftMargin(0.10);
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.4);
  p1->SetBottomMargin(0.01);
  p2->Draw();
  p1->Draw();

  p1->cd();
  TH1F* hTemp = h1_MC1;
  hTemp->Reset();
  setHistStyle(hTemp);
  if( ((TString)hTemp->GetXaxis()->GetTitle()).Contains("{p}_{T}") ){
    hTemp->GetXaxis()->SetRangeUser(1., hTemp->GetXaxis()->GetXmax());
    p1->SetLogx();
    p2->SetLogx();
  };
  hTemp->Draw("AXIS");
  
  if(gMC1) gMC1->Draw("SAME P0");
  if(gMC2) gMC2->Draw("SAME P0");

  if(!isAllHist(hTemp))
    this->drawGridLines(hTemp, yMin);
  else
    this->drawEtaRegions(hTemp, yMin);
  
  TLegend *leg = this->makeLegend(gMC1, gMC2, getOriginLabel(namePass1), getOriginLabel(namePass2));
  leg->Draw("SAME");
  this->drawLeptonFlavor(hTemp);

  p2->cd();
  this->drawRatio(vec(gMC1,gMC2), hTemp);

  if(Print){
    if(!(bool)outDir.length()) outDir = ".";
    this->checkDir(outDir);
    c->Print(Form("%s/%s.%s",outDir.c_str(),c->GetName(),figType));
    INFO("makeRatePlot2D", Form("Created %s/%s.%s",outDir.c_str(),c->GetName(),figType));
  }
}

void RatePlotter::compareSelections(TString namePass, TString nameTot, std::vector< std::pair<std::string, std::string> > selections, TString source){
  
  std::cout << std::endl;
  if(!source.Length()){   INFO("compareSelec", "No source [Data|MC] selected" ); return; }
  if(selections.empty()){ INFO("compareSelec", "No list of selection provided"); return; } 

  if( (source=="MC" && MCFiles.empty()) || (source=="Data" && DataFiles.empty()) ){ INFO("compareSelec", "No data or MC files found"); return; }

  std::vector<TGraphAsymmErrors*> RateGraphs(0); 
  TH1F* hTemp(0);

  for(auto selection : selections){
    const char *dir = selection.first.c_str();
    INFO("compareSelec", Form("Selection dir: %s (%s)", dir, selection.second.c_str()));

    std::vector<TH1F*> histos(0);
    if(source=="Data"){ 
      histos = getHistosFromList(DataFiles, dir);
      if(Prompt){
	INFO("compareSelec", Form("Subtracting prompt processes from %i files", (int)PromptMCFiles.size()));
	std::vector<TH1F*> histosPrompt = getHistosFromList(PromptMCFiles, dir);

	this->lumiScale(histosPrompt);
	this->subtractPrompt(histos, histosPrompt);
      }
    }
    if(source=="MC"){   
      histos = getHistosFromList(MCFiles, dir);
      this->lumiScale(histos);
    }
    if(histos.empty()){ INFO("compareSelec", "No source [Data|MC] selected"); return; }

    TH1F *h1(0), *h2(0);
    h1 = findHisto(nameTot,  histos);
    h2 = findHisto(namePass, histos);
    
    if(!subtractedProc.empty()){
      std::vector<TH1F*> subtractionHistos = getHistosFromList(MCFiles, dir);
      this->lumiScale(subtractionHistos);
      this->subtractMCProcess(h1, h2, subtractionHistos);
    }        
    TGraphAsymmErrors *g = getRateGraph(h2, h1);
    hTemp = h1;
    RateGraphs.push_back(g);
  }
  INFO("compareSelec", Form("Created rate plots (%s) for %i selections", source.Data(), (int)RateGraphs.size()));

  if(!Style) this->setStyle(1);  
  TString cname = Form("%s_over_%s_Selections_%s",namePass.Data(),nameTot.Data(),source.Data());
  cname = addSuffix(cname.Data());  
  TCanvas *c  = new TCanvas(cname, cname, 1, 10, 770, 560);

  TPad *p1 = new TPad(cname+"_p1",cname+"_p1", 0.00, 0.30, 1.00, 1.00, -1, 0, 0);
  TPad *p2 = new TPad(cname+"_p2",cname+"_p1", 0.00, 0.00, 1.00, 0.30, -1, 0, 0);
  c->cd();
  p2->SetLeftMargin(0.10);
  p1->SetLeftMargin(0.10);
  p2->SetTopMargin(0.02);
  p2->SetBottomMargin(0.4);
  p1->SetBottomMargin(0.01);
  p2->Draw();
  p1->Draw();

  p1->cd();
  hTemp->Reset();
  setHistStyle(hTemp);
  if( ((TString)hTemp->GetXaxis()->GetTitle()).Contains("{p}_{T}") ){
    hTemp->GetXaxis()->SetRangeUser(1., hTemp->GetXaxis()->GetXmax());
    p1->SetLogx();
    p2->SetLogx();
  }
  hTemp->Draw("AXIS");

  for(unsigned int i(0); i<RateGraphs.size(); i++){
    switch( (int)i ){
    case 0: RateGraphs[i]->SetLineColor(kBlue);    break;
    case 1: RateGraphs[i]->SetLineColor(kCyan);    break;
    case 2: RateGraphs[i]->SetLineColor(kViolet);  break;
    case 3: RateGraphs[i]->SetLineColor(kRed);     break;
    case 4: RateGraphs[i]->SetLineColor(kGreen);   break;
    case 5: RateGraphs[i]->SetLineColor(kMagenta); break;
    default: break;
    }
    RateGraphs[i]->SetMarkerColor(RateGraphs[i]->GetLineColor());
    RateGraphs[i]->Draw("SAME P0");
  }

  if(!isAllHist(hTemp))
    this->drawGridLines(hTemp, yMin);
  else
    this->drawEtaRegions(hTemp, yMin);

  TLegend *leg(0);
  switch( (int)RateGraphs.size() ){
  case 1: leg = new TLegend(0.58, 0.83, 0.91, 0.91); break;
  case 2: leg = new TLegend(0.58, 0.78, 0.91, 0.91); break;
  case 3: leg = new TLegend(0.58, 0.73, 0.91, 0.91); break;
  case 4: leg = new TLegend(0.58, 0.68, 0.91, 0.91); break;
  case 5: leg = new TLegend(0.58, 0.63, 0.91, 0.91); break;
  case 6: leg = new TLegend(0.58, 0.58, 0.91, 0.91); break;
  default: break;
  }
  leg->SetNColumns(1); 
  leg->SetBorderSize(1);
  leg->SetTextFont(42);
  leg->SetTextSize(0.05);
  leg->SetFillStyle(1111);
  leg->SetFillColor(kWhite);
  for(unsigned int i(0); i<RateGraphs.size(); i++) 
    leg->AddEntry(RateGraphs[i], Form(" Rates %s: %s", source.Data(), selections.at(i).second.c_str())); 
  leg->Draw("SAME");
  
  this->drawLeptonFlavor(hTemp);

  p2->cd();
  this->drawRatio(RateGraphs, hTemp);

  if(Print){
    if(!(bool)outDir.length()) outDir = ".";
    this->checkDir(outDir);
    c->Print(Form("%s/%s.%s",outDir.c_str(),c->GetName(),figType));
    INFO("compareSelec", Form("Created %s/%s.%s",outDir.c_str(),c->GetName(),figType));
  }
}

void RatePlotter::subtractPrompt(std::vector<TH1F*> data, std::vector<TH1F*> prompt){
  if(data.empty() || prompt.empty() || (data.size() != prompt.size())){ 
    INFO("subtractPrompt", "Lists are empty or N(data) != N(prompt)");
    return;
  }
  for(unsigned int i(0); i<data.size(); i++) this->subtract(data[i], prompt[i]);
}

void RatePlotter::subtractMCProcess(TH1F* histInputTot, TH1F *histInputPass, std::vector<TH1F*> histProc){
  
  if(subtractedProc.empty() || !histProc.size() || !histInputPass || !histInputTot){ 
    INFO("subtractMCProc", "No MC processes subtracted"); return; 
  }
  TString name = histInputTot->GetName();

  for(auto proc : subtractedProc){ 
    TH1F* hProcTot(0), *hProcPass(0);
    TString nameProcPass(""), nameProcTot("");
    float sf = getProcessSF(proc);

    if(name.Contains("el0")){
      nameProcPass = Form("histoTight_%s_electron0",proc.Data());
      nameProcTot  = Form("histoLoose_%s_electron0",proc.Data());

      if(proc.Contains("charge_flip") || proc.Contains("conversion")){
	nameProcPass = Form("histoTight_%s0",proc.Data()); 
	nameProcTot  = Form("histoLoose_%s0",proc.Data());
      }
    }
    if(name.Contains("el1")){
      nameProcPass = Form("histoTight_%s_electron1",proc.Data());
      nameProcTot  = Form("histoLoose_%s_electron1",proc.Data());

      if(proc.Contains("charge_flip") || proc.Contains("conversion")){
	nameProcPass = Form("histoTight_%s1",proc.Data());
	nameProcTot  = Form("histoLoose_%s1",proc.Data());
      }
    }
    if(name.Contains("mu0")){
      nameProcPass = Form("histoTight_%s_muon0",proc.Data());
      nameProcTot  = Form("histoLoose_%s_muon0",proc.Data());

      if(proc.Contains("charge_flip") || proc.Contains("conversion")) continue;
    }
    if(name.Contains("mu1")){
      nameProcPass = Form("histoTight_%s_muon1",proc.Data());
      nameProcTot  = Form("histoLoose_%s_muon1",proc.Data());

      if(proc.Contains("charge_flip") || proc.Contains("conversion")) continue;
    }
    hProcTot  = findHisto(nameProcTot, histProc);
    hProcPass = findHisto(nameProcPass,histProc);

    if(!hProcTot || !hProcPass){ 
      INFO("subtractMCProc", Form("No histograms [%s|%s] found", nameProcTot.Data(), nameProcPass.Data())); continue; 
    }
    INFO("subtractMCProc", Form("Subtracting histograms [%s|%s] (SF=%.1f) from [%s|%s]", hProcPass->GetName(), hProcTot->GetName(), sf, histInputPass->GetName(), histInputTot->GetName()));

    this->subtract(histInputTot,  hProcTot,  sf);
    this->subtract(histInputPass, hProcPass, sf);
  }
  return;
}

void RatePlotter::subtractMCProcess2D(TH2F* histInputTot, TH2F *histInputPass, std::vector<std::string> files){

  if(subtractedProc.empty() || !files.size() || !histInputPass || !histInputTot){
    INFO("subtractMCProc", "No MC processes subtracted"); return;
  }
  TString name = histInputTot->GetName();

  for(auto proc : subtractedProc){
    TH2F* histTempTot  = (TH2F*)histInputTot->Clone(Form("Template_%s_%s",proc.Data(),histInputTot->GetName()));
    TH2F* histTempPass = (TH2F*)histInputPass->Clone(Form("Template_%s_%s",proc.Data(),histInputPass->GetName()));
    histTempTot->Reset();
    histTempPass->Reset();
    float sf = getProcessSF(proc);
  
    TString nameProcPass(""), nameProcTot("");
    if(name.Contains("_el")){
      nameProcPass = Form("histo2D_Tight_%s_electron",proc.Data());
      nameProcTot  = Form("histo2D_Loose_%s_electron",proc.Data());

      if(proc.Contains("charge_flip") || proc.Contains("conversion")){
	nameProcPass = Form("histo2D_Tight_%s",proc.Data());
	nameProcTot  = Form("histo2D_Loose_%s",proc.Data());
      }
    }
    if(name.Contains("_mu")){
      nameProcPass = Form("histo2D_Tight_%s_muon",proc.Data());
      nameProcTot  = Form("histo2D_Loose_%s_muon",proc.Data());

      if(proc.Contains("charge_flip") || proc.Contains("conversion")) continue;
    }
    for(auto file : files){
      TFile *f = new TFile(file.c_str());
      float scale = getMCNorm(f);

      TH2F *ht = (TH2F*)f->Get( Form("%s/%s",effDir,nameProcTot.Data()) );
      TH2F *hp = (TH2F*)f->Get( Form("%s/%s",effDir,nameProcPass.Data()) );
      ht->Scale(Lumi*scale);
      hp->Scale(Lumi*scale);

      histTempTot->Add(ht);
      histTempPass->Add(hp);
    }
    INFO("subtractMCProc", Form("Subtracting histograms [%s|%s] (SF=%.1f) from [%s|%s]", histTempPass->GetName(), histTempTot->GetName(), sf, histInputPass->GetName(), histInputTot->GetName()));

    this->subtract(histInputTot,  histTempTot,  sf);
    this->subtract(histInputPass, histTempPass, sf);
  }
  return;
}

bool RatePlotter::checkEntries(TH1F* pass, TH1F* total){
  if(!pass || !total) return false;

  int bins1(total->GetNbinsX()), bins2(pass->GetNbinsX()), bins(0);
  if(bins1 != bins2){ INFO("checkEntries", "Bin numbers for h(pass) != h(tot)"); return false; } 
  bins = bins1;

  for(int i(0); i<=bins+1; i++){ if(pass->GetBinContent(i) < 0.) pass->SetBinContent(i,0); }
  for(int i(0); i<=bins+1; i++){ if(total->GetBinContent(i) < 0.) total->SetBinContent(i,0); }
  for(int i(0); i<=bins+1; i++){ if(pass->GetBinContent(i) > total->GetBinContent(i)) pass->SetBinContent(i, total->GetBinContent(i)); }
  return true;
}

bool RatePlotter::checkEntries(TH2F* pass, TH2F* total){
  if(!pass || !total) return false;

  int bins1x(total->GetNbinsX()), bins1y(total->GetNbinsY()), bins2x(pass->GetNbinsX()), bins2y(total->GetNbinsY()), binsx(0), binsy(0);
  if( (bins1x != bins2x) || (bins1y != bins2y) ){ INFO("checkEntries", "Bin numbers for h(pass) != h(tot)"); return false; } 
  binsx = bins1x; 
  binsy = bins1y;

  for(int i(0); i<=binsx+1; i++){
    for(int j(0); j<=binsy+1; j++){
      if(pass->GetBinContent(i,j) < 0.) pass->SetBinContent(i,j,0);
    }
  }
  for(int i(0); i<=binsx+1; i++){
    for(int j(0); j<=binsy+1; j++){
      if(total->GetBinContent(i,j) < 0.) total->SetBinContent(i,j,0);
    }
  }
  for(int i(0); i<=binsx+1; i++){
    for(int j(0); j<=binsy+1; j++){
      if(pass->GetBinContent(i,j) > total->GetBinContent(i,j) ) pass->SetBinContent(i,j, total->GetBinContent(i,j));
    }
  }
  return true;
}

void RatePlotter::writeToFile(TH1* hTemp, TString type, TString source, const char* outname){
  if(!hTemp){ INFO("writeToFile", "Nothing to write"); return; }
  this->checkDir(outDir);
  
  TString filename = Form("%s/%s%s_%s.root", outDir.c_str(), outname, Form("%iD",(int)hTemp->GetDimension()), source.Data());
  TFile *f = TFile::Open(filename, "UPDATE");
  f->cd();

  TString hName = hTemp->GetName();
  const char *flavor("");
  if( hName.Contains("_mu") ) flavor = "mu";
  if( hName.Contains("_el") ) flavor = "el";

  TH1F *hOut1D(0);
  TH2F *hOut2D(0);
  const char *newName("");
  const char *paraX(""), *paraY("");
  
  switch( (int)hTemp->GetDimension()){
  case 1:
    paraX = getAxisPar( this->GetXTitle( (TH1F*)hTemp) );

    newName = Form("%s%s%iD_%s_%s", type.Data(), outname, (int)hTemp->GetDimension(), flavor, paraX);
    newName = addSuffix(newName);
    hOut1D = (TH1F*) hTemp->Clone(newName);
    break;
  case 2:
    paraX = getAxisPar( (TString)hTemp->GetXaxis()->GetTitle() );
    paraY = getAxisPar( (TString)hTemp->GetYaxis()->GetTitle() );

    newName = Form("%s%s%iD_%s_%s_%s",type.Data(), outname, (int)hTemp->GetDimension(), flavor, paraX, paraY);
    newName = addSuffix(newName);
    hOut2D = (TH2F*) hTemp->Clone(newName);
    break;
  default: break;
  }
  if(hOut1D){
    subtractNominal(f,hOut1D);
    hOut1D->Write();
    INFO("writeToFile", Form("Wrote histogram %s/%s",f->GetName(),hOut1D->GetName()));
  }
  if(hOut2D){
    subtractNominal(f,hOut2D);
    hOut2D->Write();
    INFO("writeToFile", Form("Wrote histogram %s/%s",f->GetName(),hOut2D->GetName()));
  }
  f->Close();
  return;
}

void RatePlotter::addFakeHist(std::vector<TH1F*> &histos, TString opt){
  if(histos.empty() || !opt.Length()){ INFO("addFakeHist", "Empty input or no lepton type selected.. returning"); return; }

  TH1F* h0_loose(0), *h1_loose(0), *hA_loose(0);
  TH1F* h0_tight(0), *h1_tight(0), *hA_tight(0);

  if( opt=="Mu" ){
    for(auto suf : FakeSourcesMu) DEBUG("addFakeHist", Form("Merging muon histograms with suffix : %s", suf.Data()));

    h0_loose = (TH1F*) (this->findHisto("histoLoose_mu0",    histos))->Clone("histoLoose_Fakes_muon0");
    h1_loose = (TH1F*) (this->findHisto("histoLoose_mu1",    histos))->Clone("histoLoose_Fakes_muon1");
    hA_loose = (TH1F*) (this->findHisto("all_histoLoose_mu", histos))->Clone("all_histoLoose_Fakes_muon");
    
    h0_tight = (TH1F*) (this->findHisto("histoTight_mu0",    histos))->Clone("histoTight_Fakes_muon0");
    h1_tight = (TH1F*) (this->findHisto("histoTight_mu1",    histos))->Clone("histoTight_Fakes_muon1");
    hA_tight = (TH1F*) (this->findHisto("all_histoTight_mu", histos))->Clone("all_histoTight_Fakes_muon");
    
    h0_loose->Reset(); h0_tight->Reset();
    h1_loose->Reset(); h1_tight->Reset();
    hA_loose->Reset(); hA_tight->Reset();

    for(auto h : histos){
      TString hname = h->GetName();

      for(auto suf : FakeSourcesMu){

	if(hname.Contains(suf) && hname.Contains("histoLoose")     && hname.Contains("muon0")) h0_loose->Add(h); 
	if(hname.Contains(suf) && hname.Contains("histoLoose")     && hname.Contains("muon1")) h1_loose->Add(h);
	if(hname.Contains(suf) && hname.Contains("all_histoLoose") && hname.Contains("muon") ) hA_loose->Add(h);

	if(hname.Contains(suf) && hname.Contains("histoTight")     && hname.Contains("muon0")) h0_tight->Add(h);
	if(hname.Contains(suf) && hname.Contains("histoTight")     && hname.Contains("muon1")) h1_tight->Add(h);
	if(hname.Contains(suf) && hname.Contains("all_histoTight") && hname.Contains("muon") ) hA_tight->Add(h);  
      }
    }
  }

  if( opt=="El" ){
    for(auto suf : FakeSourcesEl) DEBUG("addFakeHist", Form("Merging electron histograms with suffix : %s", suf.Data()));

    h0_loose = (TH1F*) (this->findHisto("histoLoose_el0",    histos))->Clone("histoLoose_Fakes_electron0");
    h1_loose = (TH1F*) (this->findHisto("histoLoose_el1",    histos))->Clone("histoLoose_Fakes_electron1");
    hA_loose = (TH1F*) (this->findHisto("all_histoLoose_el", histos))->Clone("all_histoLoose_Fakes_electron");
    
    h0_tight = (TH1F*) (this->findHisto("histoTight_el0",    histos))->Clone("histoTight_Fakes_electron0");
    h1_tight = (TH1F*) (this->findHisto("histoTight_el1",    histos))->Clone("histoTight_Fakes_electron1");
    hA_tight = (TH1F*) (this->findHisto("all_histoTight_el", histos))->Clone("all_histoTight_Fakes_electron");
    
    h0_loose->Reset(); h0_tight->Reset();
    h1_loose->Reset(); h1_tight->Reset();
    hA_loose->Reset(); hA_tight->Reset();

    for(auto h : histos){
      TString hname = h->GetName();

      for(auto suf : FakeSourcesEl){

	if(hname.Contains(suf) && hname.Contains("histoLoose")     && hname.Contains("electron0")) h0_loose->Add(h);
	if(hname.Contains(suf) && hname.Contains("histoLoose")     && hname.Contains("electron1")) h1_loose->Add(h);
	if(hname.Contains(suf) && hname.Contains("all_histoLoose") && hname.Contains("electron") ) hA_loose->Add(h);
	
	if( (suf=="charge_flip" || suf=="conversion") ){ 
	  if( hname.Contains(suf) && hname.Contains("histoLoose")  && hname.Contains("0") ) h0_loose->Add(h);
	  if( hname.Contains(suf) && hname.Contains("histoLoose")  && hname.Contains("1") ) h1_loose->Add(h);
	  if( hname.Contains(suf) && hname.Contains("all_histoLoose")                     ) hA_loose->Add(h);
	}

	if(hname.Contains(suf) && hname.Contains("histoTight")     && hname.Contains("electron0")) h0_tight->Add(h);
	if(hname.Contains(suf) && hname.Contains("histoTight")     && hname.Contains("electron1")) h1_tight->Add(h);
	if(hname.Contains(suf) && hname.Contains("all_histoTight") && hname.Contains("electron") ) hA_tight->Add(h);
	
	if( (suf=="charge_flip" || suf=="conversion") ){  
	  if( hname.Contains(suf) && hname.Contains("histoTight")  && hname.Contains("0") ) h0_tight->Add(h);
	  if( hname.Contains(suf) && hname.Contains("histoTight")  && hname.Contains("1") ) h1_tight->Add(h);
	  if( hname.Contains(suf) && hname.Contains("all_histoTight")                     ) hA_tight->Add(h);
	}
      }
    }
  } 
  histos.push_back(h0_loose); histos.push_back(h0_tight);
  histos.push_back(h1_loose); histos.push_back(h1_tight);
  histos.push_back(hA_loose); histos.push_back(hA_tight);
  return;
}


void RatePlotter::getMCSources(TString flavor, TString quality, bool log){

  if(flavor!="el" && flavor!="mu"){ INFO("getMCSources", "No lepton type selected. Please set [el|mu]"); return; }
  if(quality!="Loose" && quality!="Tight"){ INFO("getMCSources", "No lepton quality selected. Please set [Tight|Loose]"); return; }
  
  std::cout << std::endl;
  if(!MCRates){ INFO("getMCSources", "No MC input provided"); return;}

  histosMC = getHistosFromList(MCFiles, effDir);
  this->lumiScale(histosMC);

  std::vector<TString> sources(0);
  if(flavor=="el") sources = FakeSourcesEl;
  if(flavor=="mu") sources = FakeSourcesMu;
  if(sources.empty()){ INFO("getMCSources", "No MC sources selected. Please call setFakeSources(Electron/Muon) before"); return; }

  sources.insert(sources.begin(), "prompt");
  sources.push_back("Fakes");

  std::vector<TH1F*> hSources0(0), hSources1(0);
  for(auto source : sources){
    for(auto h : histosMC){
      if(isAllHist(h)) continue;
      TString name = h->GetName();
    
      if(name.Contains(flavor) && name.Contains(quality) && name.Contains(source) && name.Contains("0")) hSources0.push_back(h);
      if(name.Contains(flavor) && name.Contains(quality) && name.Contains(source) && name.Contains("1")) hSources1.push_back(h);

      if(flavor=="el" && name.Contains(quality) && name.Contains(source) && source=="charge_flip" && name.Contains("0")) hSources0.push_back(h);
      if(flavor=="el" && name.Contains(quality) && name.Contains(source) && source=="conversion"  && name.Contains("0")) hSources0.push_back(h);
      
      if(flavor=="el" && name.Contains(quality) && name.Contains(source) && source=="charge_flip" && name.Contains("1")) hSources1.push_back(h);
      if(flavor=="el" && name.Contains(quality) && name.Contains(source) && source=="conversion"  && name.Contains("1")) hSources1.push_back(h);
    }
  }
  for(auto h : hSources0) DEBUG("getMCSources", Form("Looking at %s \t :: Nevents=%.8f",h->GetName(),h->Integral()));
  for(auto h : hSources1) DEBUG("getMCSources", Form("Looking at %s \t :: Nevents=%.8f",h->GetName(),h->Integral()));

  TCanvas *c[2];
  TString cname[2];
  for(unsigned int i(0); i<2; i++){
    cname[i] = Form("Sources_%s%s_p%i",flavor.Data(),quality.Data(),i);
    cname[i] = addSuffix(cname[i].Data());

    c[i] =  new TCanvas(cname[i], cname[i], 1, 10, 770, 560);
    c[i]->SetLeftMargin(0.10);
    if(log) c[i]->SetLogy();
  }

  TH1F* hTemp0 = (TH1F*)(hSources0.back())->Clone(Form("Template_%s",hSources0.back()->GetName()));
  TH1F* hTemp1 = (TH1F*)(hSources1.back())->Clone(Form("Template_%s",hSources1.back()->GetName()));

  hTemp0->GetYaxis()->SetRangeUser(1, (log ? hTemp0->GetMaximum()*700 : hTemp0->GetMaximum()*2));
  hTemp1->GetYaxis()->SetRangeUser(1, (log ? hTemp1->GetMaximum()*700 : hTemp1->GetMaximum()*2));
  hTemp0->Reset();
  hTemp1->Reset();
  setHistStyleNoRatio(hTemp0);
  setHistStyleNoRatio(hTemp1);
  
  if( ((TString)hTemp0->GetXaxis()->GetTitle()).Contains("{p}_{T}") ){
    hTemp0->GetXaxis()->SetRangeUser(1., hTemp0->GetXaxis()->GetXmax());
    c[0]->SetLogx();  
  };
  if( ((TString)hTemp1->GetXaxis()->GetTitle()).Contains("{p}_{T}") ){
    hTemp1->GetXaxis()->SetRangeUser(1., hTemp1->GetXaxis()->GetXmax());
    c[1]->SetLogx();
  };
  c[0]->cd();
  hTemp0->Draw("AXIS");
  this->drawGridLines(hTemp0, 1);

  float off0(0.05);
  for(auto h : hSources0){
    this->setSourceStyle(h);
    h->SetBarWidth(0.11);
    h->SetBarOffset(off0);
    h->Draw("BAR SAME");
    off0 += 0.11;
  }
  c[1]->cd();
  hTemp1->Draw("AXIS");
  this->drawGridLines(hTemp1, 1);

  float off1(0.05);
  for(auto h : hSources1){
    this->setSourceStyle(h);
    h->SetBarWidth(0.11);
    h->SetBarOffset(off1);
    h->Draw("BAR SAME");
    off1 += 0.11;
  }

  c[0]->cd();
  TLegend *leg(0);
  switch( (int)hSources0.size() ){
  case 3: leg = new TLegend(0.63, 0.83, 0.85, 0.91); break;
  case 4: leg = new TLegend(0.63, 0.78, 0.85, 0.91); break;
  case 5: leg = new TLegend(0.63, 0.73, 0.85, 0.91); break;
  case 6: leg = new TLegend(0.63, 0.68, 0.85, 0.91); break;
  case 7: leg = new TLegend(0.63, 0.63, 0.85, 0.91); break;
  case 8: leg = new TLegend(0.63, 0.58, 0.85, 0.91); break;
  default: break;
  }
  leg->SetNColumns(1);
  leg->SetBorderSize(1);
  leg->SetTextFont(42);
  leg->SetTextSize(0.03);
  leg->SetFillStyle(1111);
  leg->SetFillColor(kWhite);
  for(unsigned int i(0); i<hSources0.size(); i++)
    leg->AddEntry(hSources0[i], getOriginLabel(hSources0[i]->GetName()) );
  leg->Draw("SAME");
  this->drawLeptonFlavor(hTemp0, true);
  c[0]->RedrawAxis();

  c[1]->cd();
  leg->Draw("SAME");
  this->drawLeptonFlavor(hTemp1, true);
  c[1]->RedrawAxis();

  if(Print){
    if(!(bool)outDir.length()) outDir = ".";
    this->checkDir(outDir);

    for(unsigned int i(0); i<2; i++){
      c[i]->Print(Form("%s/%s.%s",outDir.c_str(),c[i]->GetName(),figType));
      INFO("getMCSources", Form("Created %s/%s.%s",outDir.c_str(),c[i]->GetName(),figType));
    }
  }
  return;
}

void::RatePlotter::setSourceStyle(TH1F *h){
  if(!h) return;
  TString name = h->GetName();
  if(name.Contains("_HF"))              h->SetFillColor(kGreen-6);
  if(name.Contains("_LF"))              h->SetFillColor(kBlue-7);
  if(name.Contains("_charge_flip"))     h->SetFillColor(kMagenta-7);
  if(name.Contains("_conversion"))      h->SetFillColor(kCyan-3);  
  if(name.Contains("_Tau"))             h->SetFillColor(kViolet-8);
  if(name.Contains("_not_classified"))  h->SetFillColor(kOrange-9);
  if(name.Contains("_Fakes"))           h->SetFillColor(kRed+3);
  if(name.Contains("_prompt"))          h->SetFillColor(kGray);

  h->SetLineColor(kWhite);
  h->SetMarkerColor(kWhite);
  h->SetLineWidth(0);
  h->SetMarkerSize(0);
  return;
}

void::RatePlotter::subtractNominal(TFile *f, TH1 *hVar){
  if(!subNomRate || !sysSuffix.length() || !hVar) return;
  TString nameVar(hVar->GetName());

  TH1* h(0); TKey *key(0);
  TIter next(f->GetListOfKeys());
  while(( key = (TKey*)next() )){
    TObject *obj = key->ReadObj();
    if(!obj->IsA()->InheritsFrom("TH1")) continue;
    h = (TH1*)obj;
    TString name = h->GetName();
    if( nameVar.Contains(name) && !name.Contains("__") ) break;
  }
  if(!h){ INFO("subtractNominal", Form("No nominal histogram for variation %s is found. Check root file", hVar->GetName())); return; }
  INFO("subtractNominal", Form("Subtract %s (nominal) from %s (syst. variation)",h->GetName(),hVar->GetName()));

  switch( (int)hVar->GetDimension() ){
  case 1:
    for(int i(1); i<=h->GetNbinsX(); i++){
      float nom = h->GetBinContent(i), var = hVar->GetBinContent(i);

      hVar->SetBinContent(i, TMath::Abs(var-nom));
      DEBUG("subtractNominal", Form("Bin (%i) Subtract %.3f vom %.3f --> %.3f",i,nom,var,hVar->GetBinContent(i)));
    }
    break;
  case 2:
    for(int i(1); i<=h->GetNbinsX(); i++){
      for(int j(1); j<=h->GetNbinsY(); j++){
	float nom = h->GetBinContent(i,j), var = hVar->GetBinContent(i,j);

	hVar->SetBinContent(i, j, TMath::Abs(var-nom));
	DEBUG("subtractNominal", Form("Bin (%i|%i) Subtract %.3f vom %.3f --> %.3f",i,j,nom,var,hVar->GetBinContent(i,j)));
      }
    }
    break;
  default: break;
  }
  return;
}

void RatePlotter::subtract(TH1 *h1, TH1*h2, float sf){
  if(!h1 || !h2) return;
  if(h1->GetDimension() != h2->GetDimension()){ INFO("subtractHist", "Cannot subtract TH1s with different dimensions"); return; }

  switch( (int)h1->GetDimension() ){
  case 1:
    for(int i(0); i<=h1->GetNbinsX()+1; i++){
      float val = TMath::Max(h1->GetBinContent(i) - h2->GetBinContent(i)*sf,0.);
      float err = TMath::Sqrt(TMath::Power(h1->GetBinError(i),2) + TMath::Power(h2->GetBinError(i),2));
      h1->SetBinContent(i,val);
      h1->SetBinError(i,err);
    }
    break;
  case 2:
    for(int i(0); i<=h1->GetNbinsX()+1; i++){
      for(int j(0); j<=h1->GetNbinsY()+1; j++){
	float val = TMath::Max(h1->GetBinContent(i,j) - h2->GetBinContent(i,j)*sf,0.);
	float err = TMath::Sqrt(TMath::Power(h1->GetBinError(i,j),2)+TMath::Power(h2->GetBinError(i,j),2));
	h1->SetBinContent(i,j,val);
	h1->SetBinError(i,j,err);
      }
    }
    break;
  default: break;
  }
  return;
}

