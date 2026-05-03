{
  gROOT->LoadMacro("RatePlotter.cxx++");
  gROOT->ProcessLine("RatePlotter Plotter");
  gErrorIgnoreLevel = kFatal;

  Plotter.setDebug(0);
  Plotter.setPrint(0);
  Plotter.setFigureFormat("eps");
  Plotter.setEffDirectory("Efficiencies_Selection_2");

  Plotter.setStylePath("~/Atlas/atlasrootstyle/AtlasStyle.C"); //location of AtlasStyle.C (for plotting)
  Plotter.setStyle(1);
  Plotter.drawAtlasLabel(1);
  
  Plotter.writeHistFile("Efficiency", false);
  Plotter.setSysSuffix("");
  Plotter.subtractNominalRates(1);

  TString type = "mu";

  float lumi = 138965.2;
  Plotter.setLumi(lumi);

  TString dataType = Form("Data (%.0f fb^{-1})",lumi/1000.);
  TString mcProc   = "MC"; 

  Plotter.setRateType("Fake");
  Plotter.setLabel(mcProc.Data(),   "MC");
  Plotter.setLabel(dataType.Data(), "Data");

  std::vector<TString> sourcesMuon     = {"LF", "HF", "Tau","not_classified"};
  std::vector<TString> sourcesElectron = {"LF", "HF", "Tau", "charge_flip", "conversion", "not_classified"};
  Plotter.setFakeSourcesMuon(sourcesMuon);
  Plotter.setFakeSourcesElectron(sourcesElectron);

  Plotter.setHistRange(0.01, 1.29);
  Plotter.setRatioRange(0.1, 2.30);
  Plotter.setOutDir("/home/cardillo/FakeRates_ttZ/Plots/TEST/");
  
  Plotter.getFiles(""); //location of input files (2nd argument is wildcard-selection

  //Truth composition plots
  if(type=="el"){
    Plotter.getMCSources("el", "Loose",1);
    Plotter.getMCSources("el", "Tight",1);
  }
  if(type=="mu"){
    Plotter.getMCSources("mu", "Loose",1);
    Plotter.getMCSources("mu", "Tight",1);
  }

  //Lepton source to be subtracted from the rates
  Plotter.setProcessSubtraction("charge_flip", 1.0);
  Plotter.setProcessSubtraction("prompt",      1.0);
  */
  //Standard efficiency plots
  if(type=="el"){
    Plotter.makeRatePlot("histoTight_el0", "histoLoose_el0");
    Plotter.makeRatePlot("histoTight_el1", "histoLoose_el1");
  
    Plotter.makeRatePlot2D("histo2D_Tight_el", "histo2D_Loose_el", "Data");
  }
  if(type=="mu"){
    Plotter.makeRatePlot("histoTight_mu0", "histoLoose_mu0");
    Plotter.makeRatePlot("histoTight_mu1", "histoLoose_mu1");

    Plotter.makeRatePlot2D("histo2D_Tight_mu", "histo2D_Loose_mu", "Data");
  }
  
  //Rates for different origins
  bool originPlots = false;
  if(originPlots){
    Plotter.compareMCRates("histoTight_Fakes_electron0",    "histoLoose_Fakes_electron0",    "histoTight_HF_electron0",    "histoLoose_HF_electron0",   "red", "green");
    Plotter.compareMCRates("histoTight_Fakes_electron1",    "histoLoose_Fakes_electron1",    "histoTight_HF_electron1",    "histoLoose_HF_electron1",   "red", "green");
    Plotter.compareMCRates("all_histoTight_Fakes_electron", "all_histoLoose_Fakes_electron", "all_histoTight_HF_electron", "all_histoLoose_HF_electron","red", "green");

    Plotter.compareMCRates("histoTight_HF_electron0",    "histoLoose_HF_electron0",    "histoTight_LF_electron0",    "histoLoose_LF_electron0",   "green", "blue");
    Plotter.compareMCRates("histoTight_HF_electron1",    "histoLoose_HF_electron1",    "histoTight_LF_electron1",    "histoLoose_LF_electron1",   "green", "blue");
    Plotter.compareMCRates("all_histoTight_HF_electron", "all_histoLoose_HF_electron", "all_histoTight_LF_electron", "all_histoLoose_LF_electron","green", "blue");

    Plotter.compareMCRates("histoTight_Fakes_muon0",    "histoLoose_Fakes_muon0",    "histoTight_HF_muon0",    "histoLoose_HF_muon0",   "red", "green");
    Plotter.compareMCRates("histoTight_Fakes_muon1",    "histoLoose_Fakes_muon1",    "histoTight_HF_muon1",    "histoLoose_HF_muon1",   "red", "green");
    Plotter.compareMCRates("all_histoTight_Fakes_muon", "all_histoLoose_Fakes_muon", "all_histoTight_HF_muon", "all_histoLoose_HF_muon","red", "green");

    Plotter.compareMCRates("histoTight_HF_muon0",    "histoLoose_HF_muon0",    "histoTight_LF_muon0",    "histoLoose_LF_muon0",    "green",  "blue");
    Plotter.compareMCRates("histoTight_HF_muon1",    "histoLoose_HF_muon1",    "histoTight_LF_muon1",    "histoLoose_LF_muon1",    "green",  "blue");
    Plotter.compareMCRates("all_histoTight_HF_muon", "all_histoLoose_HF_muon", "all_histoTight_LF_muon", "all_histoLoose_LF_muon", "green",  "blue");
  
    Plotter.compareMCRates("histoTight_HF_electron0",    "histoLoose_HF_electron0",    "histoTight_charge_flip0",    "histoLoose_charge_flip0",    "green", "violet");
    Plotter.compareMCRates("histoTight_HF_electron1",    "histoLoose_HF_electron1",    "histoTight_charge_flip1",    "histoLoose_charge_flip1",    "green", "violet");
    Plotter.compareMCRates("all_histoTight_HF_electron", "all_histoLoose_HF_electron", "all_histoTight_charge_flip", "all_histoLoose_charge_flip", "green", "violet");

    Plotter.compareMCRates("histoTight_HF_muon0",    "histoLoose_HF_muon0",    "histoTight_not_classified_muon0",    "histoLoose_not_classified_muon0",    "green", "orange");
    Plotter.compareMCRates("histoTight_HF_muon1",    "histoLoose_HF_muon1",    "histoTight_not_classified_muon1",    "histoLoose_not_classified_muon1",    "green", "orange");
    Plotter.compareMCRates("all_histoTight_HF_muon", "all_histoLoose_HF_muon", "all_histoTight_not_classified_muon", "all_histoLoose_not_classified_muon", "green", "orange");
  }

  //Rates for different preselections
  std::vector< std::pair<std::string, std::string> > selections = { {"Efficiencies_Selection_1", "0 b-jets" },
								    {"Efficiencies_Selection_2", "#geq 1 b-jets" },
								    {"Efficiencies_Selection_3", "#geq 2 b-jets" } };

  if(type=="el"){
    Plotter.compareSelections("histoTight_el0", "histoLoose_el0", selections, "MC");
    Plotter.compareSelections("histoTight_el1", "histoLoose_el1", selections, "MC");

    Plotter.compareSelections("histoTight_el0", "histoLoose_el0", selections, "Data");
    Plotter.compareSelections("histoTight_el1", "histoLoose_el1", selections, "Data");
  }
  if(type=="mu"){
    Plotter.compareSelections("histoTight_mu0", "histoLoose_mu0", selections, "MC");
    Plotter.compareSelections("histoTight_mu1", "histoLoose_mu1", selections, "MC");
    
    Plotter.compareSelections("histoTight_mu0", "histoLoose_mu0", selections, "Data");
    Plotter.compareSelections("histoTight_mu1", "histoLoose_mu1", selections, "Data");
  }
}

