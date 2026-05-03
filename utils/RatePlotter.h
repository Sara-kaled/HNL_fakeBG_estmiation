#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <math.h>
#include "TSystem.h"
#include "TStyle.h"
#include "TFile.h"
#include "TKey.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TEfficiency.h"
#include "TLegend.h"
#include "TLine.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"

class RatePlotter
{
 public:
  RatePlotter(std::string name = "RatePlotter"){
    CNAME = name;
    std::cout << "Initialize Class " << CNAME << std::endl;
    Debug = 0;
    Print = 0;
    Style = 0;
    MCRates    = 0;
    DataRates  = 0;
    Prompt     = 0;
    writeHist  = 0;
    AtlasLabel = 0;
    subNomRate = 0;
    yMin      = 0.0;
    yMax      = 1.0;
    yRMin     = 0.0;
    yRMax     = 1.0;
    Lumi      = 1.0;
    InFiles.clear();
    MCFiles.clear();
    DataFiles.clear();
    PromptMCFiles.clear();
    outDir    = "";
    effDir    = "";
    RateType  = "";
    mcLabel   = "";
    dataLabel = "";
    stylePath = "";
    outFile   = "";
    sysSuffix = "";
    figType   = "pdf";
    histosMC.clear();
    histosData.clear();
    FakeSourcesEl.clear();
    FakeSourcesMu.clear();
    subtractedProc.clear();
    subtractedProcSF.clear();
  };
  ~RatePlotter(){};

 public:
  void setDebug(bool debug);
  void setPrint(bool print);
  void setStyle(bool setAtlas);
  void setLumi(float lumi);
  void setOutDir(std::string dir);
  void setLabel(std::string label, std::string option);
  void setHistRange(float min=0.01, float max=1.1);
  void setRatioRange(float min=0.1, float max=2.3);
  void getMCcolor(TGraphAsymmErrors *g, TString col); 

  void getFiles(const char* dirname, const char* key1="", const char* key2="");
  void addMCFile(const char* filename);
  void addDataFile(const char* filename);
  void addPromptFile(const char* filename);
  void setEffDirectory(const char* dir);
  void setStylePath(const char* path);
  void setRateType(const char *type);
  void setFigureFormat(const char* format){figType = format;}
  void writeHistFile(const char* outname, bool writeH=false);
  void writeToFile(TH1* hTemp, TString type, TString source, const char* outname);

  void unsetData(){ DataFiles.clear(); DataRates = false; }
  void unsetMC(){ MCFiles.clear(); MCRates = false; }
  void unsetPrompt(){ PromptMCFiles.clear(); Prompt = false; }

  void addFakeHist(std::vector<TH1F*> &histos, TString opt);
  void setFakeSourcesMuon(std::vector<TString> s){ FakeSourcesMu = s; }
  void setFakeSourcesElectron(std::vector<TString> s){ FakeSourcesEl = s; }
  void setSysSuffix(std::string suf){ sysSuffix = suf; }

  void setHistStyle(TH1F* h);
  void setHistStyle(TH2F *h, float zMax);
  void setSourceStyle(TH1F *h); 
  void setHistStyleNoRatio(TH1F *h);
  void setRatioHistStyle(TH1F* h);
  void drawGridLines(TH1F* h, float yEnd);
  void draw2DGridLines(TH2F* h);
  void drawLeptonFlavor(TH1F* h, bool no2pad=false);
  void draw2DplotLabel(TH2F* h, TString type, TString source);
  void drawEtaRegions(TH1F* h, float yEnd, bool binLabels=false, int etaBins=5);
  void drawRatio(std::vector<TGraphAsymmErrors*> graphs, TH1F* h);
  void drawAtlasLabel(bool draw){AtlasLabel = draw;}
  void lumiScale(std::vector<TH1F*> histos);
  
  void subtract(TH1 *h1, TH1* h2, float sf=1.);
  void subtractPrompt(std::vector<TH1F*> data, std::vector<TH1F*> prompt);
  void setProcessSubtraction(TString proc, float sf=1.);
  void subtractMCProcess(TH1F* histInputTot, TH1F *histInputPass, std::vector<TH1F*> histProc);
  void subtractMCProcess2D(TH2F* histInputTot, TH2F *histInputPass, std::vector<std::string> files);

  void subtractNominal(TFile *f, TH1 *hVar);
  void subtractNominalRates(bool sub){ subNomRate = sub; }

  bool checkEntries(TH1F *pass, TH1F *total);
  bool checkEntries(TH2F *pass, TH2F *total);

  void checkDir(std::string dir){ if(!gSystem->OpenDirectory(outDir.c_str())){ throw std::invalid_argument("Directory "+ outDir + " does not exist");}}

  void getMCSources(TString flavor, TString quality, bool log=true);
  void makeRatePlot(TString namePass, TString nameTot);
  void makeRatePlot2D(TString namePass, TString nameTot, TString source);
  void compareSelections(TString namePass, TString nameTot, 
			 std::vector< std::pair<std::string, std::string> > selections, TString source="");
  
  void compareMCRates(TString namePass1, TString nameTot1,
		      TString namePass2, TString nameTot2, TString col1="blue", TString col2="green");

  bool isAllHist(TH1F* h){return ((TString)(h->GetName())).Contains("all_"); }
  
  TString MSG(const char* fname){return Form("%s::%s() \t",CNAME.c_str(),fname);}

  void INFO(const char* app,  const char* msg){std::cout << Form("%s::%s() \t\t INFO \t %s",CNAME.c_str(),app,msg) << std::endl;}
  void DEBUG(const char* app, const char* msg){if(Debug) std::cout << Form("%s::%s() \t\t DEBUG \t %s",CNAME.c_str(),app,msg) << std::endl;}
  void ERROR(const char* app, const char* msg){std::cout << Form("%s::%s() \t\t ERROR \t %s",CNAME.c_str(),app,msg) << std::endl; exit(1);}

  TString GetXTitle(TH1F *h);
  TString getOriginLabel(TString name);

  float getMCNorm(TFile *f);
  float getProcessSF(TString proc);

  const char* getAxisPar(TString name);
  const char* addSuffix(const char* name);  

  TLegend* makeLegend(TGraphAsymmErrors *g1, TGraphAsymmErrors *g2, TString type1, TString type2);

  TFile* findFile(const char* name);
  TH1F*  findHisto(TString name, std::vector<TH1F*> histos);
  
  TH1F* divideTH1(TH1F* hPass, TH1F *hTotal);
  TH2F* divideTH2(TH2F* hPass, TH2F *hTotal);

  TGraphAsymmErrors* getRateGraph(TH1F *hPass, TH1F *hTotal, TString source="");

  std::vector<TH1F*> getHistos(const char* filename, const char* dirname);
  std::vector<TH1F*> getHistosFromList(std::vector<std::string> filelist, const char* dirname);

  std::vector<TGraphAsymmErrors*> vec(TGraphAsymmErrors *g1, TGraphAsymmErrors *g2);
  
 private:
  std::string CNAME;
  bool  Debug;
  bool  Print;
  bool  Style;

  bool MCRates;
  bool DataRates;
  bool Prompt;
  bool writeHist;
  bool AtlasLabel;
  bool subNomRate;

  float Lumi;
  float yMin;
  float yMax;
  float yRMin;
  float yRMax;

  TString RateType;

  const char *effDir;
  const char *stylePath;
  const char *outFile;
  const char *figType;

  std::vector<TString> FakeSourcesEl;
  std::vector<TString> FakeSourcesMu;
  std::vector<TString> subtractedProc;
  std::vector<float>   subtractedProcSF;

  std::string outDir;
  std::string mcLabel;
  std::string dataLabel;
  std::string sysSuffix;

  std::vector<TFile*> InFiles;
  
  std::vector<std::string> MCFiles;
  std::vector<std::string> DataFiles;
  std::vector<std::string> PromptMCFiles;

  std::vector<TH1F*> histosMC;
  std::vector<TH1F*> histosData;
  
};

