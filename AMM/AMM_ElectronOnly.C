// =============================================================================
// FILE:    AMM/AMM_ElectronOnly.C
// PURPOSE: Apply the Asymptotic Matrix Method weight per event for the electron
//          side. Loops over input ntuple, computes per-event AMM weight from
//          ε_f and ε_r maps, writes a clone of the input tree with new branches:
//             fake_weight_nominal
//             fake_weight_MET_up / MET_down
//             fake_weight_comp_up / comp_down
//             fake_weight_statUp / statDown
//
// PRODUCES TWO OUTPUTS:
//   - data_with_electron_fake_weights.root   (Run-2 data)
//   - mc_with_electron_fake_weights.root     (MC for closure / plotting)
//
// Default entry point AMM_ElectronOnly() runs both sequentially.
// =============================================================================
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH2F.h>
#include <TTree.h>
#include "ROOT/RVec.hxx"
#include <iostream>
#include <vector>
#include <string>
#include <cmath>

using namespace ROOT;
using namespace ROOT::VecOps;

struct FakeEfficiencySet { TH2F* el_eff = nullptr; };
struct EventWeightSet    { double fake_weight = 0.0; };

// ---------------------------------------------------------------------------
// Per-electron AMM weight: w_T (tight) and w_L (loose-only).
// Sum over all loose electrons in the event = predicted fake yield contribution.
// ---------------------------------------------------------------------------
static double computeLeptonWeight(bool isTight, double real_eff, double fake_eff)
{
    real_eff = std::max(1e-5, std::min(1.0 - 1e-5, real_eff));
    fake_eff = std::max(1e-5, std::min(real_eff - 1e-5, fake_eff));
    double denom = real_eff - fake_eff;
    if (std::abs(denom) < 1e-9) return 0.0;
    if (isTight) return (fake_eff / denom) * (real_eff - 1.0);  // w_T < 0
    else         return (fake_eff / denom) *  real_eff;          // w_L > 0
}

static double getEfficiency(TH2F* hist, double pt, double eta)
{
    if (!hist) return 0.0;
    double abs_eta = std::abs(eta);
    int binX = std::max(1, std::min(hist->GetXaxis()->FindBin(pt),
                                    hist->GetNbinsX()));
    int binY = std::max(1, std::min(hist->GetYaxis()->FindBin(abs_eta),
                                    hist->GetNbinsY()));
    return hist->GetBinContent(binX, binY);
}

// ---------------------------------------------------------------------------
// Core processing routine: read input ntuples, apply AMM, write output.
// ---------------------------------------------------------------------------
static int AMM_runOne(const std::vector<std::string>& inputs,
                     const std::string& outputFile,
                     const std::string& label)
{
    std::cout << "\n==============================================================\n"
              << "  AMM pass: " << label << "\n"
              << "  → " << outputFile << "\n"
              << "==============================================================\n";

    // ----- Efficiency histograms -----
    TFile* fakeFile = TFile::Open("../efficiency/outputs/fake_eff_with_systematics.root");
    TFile* realFile = TFile::Open("../efficiency/outputs/real_eff_histograms.root");
    if (!fakeFile || fakeFile->IsZombie()) {
        std::cerr << "ERROR: cannot open fake_eff_with_systematics.root\n"; return 1;
    }
    if (!realFile || realFile->IsZombie()) {
        std::cerr << "ERROR: cannot open real_eff_histograms.root\n"; return 1;
    }

    TH2F* real_tight = (TH2F*)realFile->Get("h2_tight_eta_pt_CR");
    TH2F* real_loose = (TH2F*)realFile->Get("h2_loose_eta_pt_CR");
    if (!real_tight || !real_loose) {
        std::cerr << "ERROR: real-eff histograms missing\n"; return 1;
    }
    TH2F* real_eff_el = (TH2F*)real_tight->Clone("real_efficiency_el");
    real_eff_el->Divide(real_loose);

    FakeEfficiencySet nominal_effs   = {(TH2F*)fakeFile->Get("fake_rate_nominal")};
    FakeEfficiencySet met_up_effs    = {(TH2F*)fakeFile->Get("fake_rate_MET_up")};
    FakeEfficiencySet met_down_effs  = {(TH2F*)fakeFile->Get("fake_rate_MET_down")};
    FakeEfficiencySet comp_up_effs   = {(TH2F*)fakeFile->Get("fake_rate_comp_up")};
    FakeEfficiencySet comp_down_effs = {(TH2F*)fakeFile->Get("fake_rate_comp_down")};
    FakeEfficiencySet real_up_effs   = {(TH2F*)fakeFile->Get("fake_rate_real_up")};
    FakeEfficiencySet real_down_effs = {(TH2F*)fakeFile->Get("fake_rate_real_down")};
    FakeEfficiencySet stat_up_effs   = {(TH2F*)fakeFile->Get("fake_rate_statUp")};
    FakeEfficiencySet stat_down_effs = {(TH2F*)fakeFile->Get("fake_rate_statDown")};

    // ----- Input chain -----
    TChain* tree = new TChain("tree");
    for (const auto& f : inputs) tree->Add(f.c_str());

    Long64_t n = tree->GetEntries();
    if (n == 0) {
        std::cerr << "ERROR: input chain is empty for " << label << "\n";
        return 1;
    }
    std::cout << "Total entries: " << n << "\n";

    RVecF*      el_pt_NOSYS            = nullptr;
    RVecF*      el_eta                 = nullptr;
    RVec<char>* el_select_loose_NOSYS  = nullptr;
    RVec<char>* el_select_tight_NOSYS  = nullptr;
    tree->SetBranchAddress("el_pt_NOSYS",            &el_pt_NOSYS);
    tree->SetBranchAddress("el_eta",                 &el_eta);
    tree->SetBranchAddress("el_select_loose_NOSYS",  &el_select_loose_NOSYS);
    tree->SetBranchAddress("el_select_tight_NOSYS",  &el_select_tight_NOSYS);

    // ----- Output -----
    TFile* outFile = new TFile(outputFile.c_str(), "RECREATE");
    TTree* outTree = tree->CloneTree(0);
    outTree->SetDirectory(outFile);

    EventWeightSet nominal_w, met_up_w, met_down_w,
                   comp_up_w, comp_down_w, real_up_w, real_down_w,
                   stat_up_w, stat_down_w;
    // data15 lacks per-trigger matching branches and uses globalTriggerMatch
    // instead.  We persist the flag per event so downstream macros can read
    // it directly from the AMM-cloned ntuple (which mixes all data years
    // into one file).
    Bool_t isData15_out = false;
    outTree->Branch("isData15",              &isData15_out,           "isData15/O");
    outTree->Branch("fake_weight_nominal",   &nominal_w.fake_weight);
    outTree->Branch("fake_weight_MET_up",    &met_up_w.fake_weight);
    outTree->Branch("fake_weight_MET_down",  &met_down_w.fake_weight);
    outTree->Branch("fake_weight_comp_up",   &comp_up_w.fake_weight);
    outTree->Branch("fake_weight_comp_down", &comp_down_w.fake_weight);
    outTree->Branch("fake_weight_real_up",   &real_up_w.fake_weight);
    outTree->Branch("fake_weight_real_down", &real_down_w.fake_weight);
    outTree->Branch("fake_weight_statUp",    &stat_up_w.fake_weight);
    outTree->Branch("fake_weight_statDown",  &stat_down_w.fake_weight);

    for (Long64_t i = 0; i < n; ++i) {
        if (i % 1000000 == 0)
            std::cout << "  entry " << i << " / " << n << "\n";
        tree->GetEntry(i);

        // Detect data15 from the filename of the currently-open input file.
        // For MC inputs this stays false (paths contain mc16*, fresh_ntuples/<proc>).
        TFile* curr = tree->GetCurrentFile();
        std::string cf = curr ? std::string(curr->GetName()) : std::string();
        isData15_out  = (cf.find("data15") != std::string::npos);

        nominal_w.fake_weight   = 0.0;
        met_up_w.fake_weight    = 0.0;
        met_down_w.fake_weight  = 0.0;
        comp_up_w.fake_weight   = 0.0;
        comp_down_w.fake_weight = 0.0;
        real_up_w.fake_weight   = 0.0;
        real_down_w.fake_weight = 0.0;
        stat_up_w.fake_weight   = 0.0;
        stat_down_w.fake_weight = 0.0;

        if (el_pt_NOSYS && el_select_loose_NOSYS) {
            for (size_t j = 0; j < el_pt_NOSYS->size(); ++j) {
                bool isLoose = (j < el_select_loose_NOSYS->size())
                               && (*el_select_loose_NOSYS)[j];
                if (!isLoose) continue;
                bool isTight = (j < el_select_tight_NOSYS->size())
                               && (*el_select_tight_NOSYS)[j];
                double pt  = (*el_pt_NOSYS)[j];
                double eta = (*el_eta)[j];
                double r   = getEfficiency(real_eff_el, pt, eta);

                nominal_w.fake_weight   += computeLeptonWeight(isTight, r,
                    getEfficiency(nominal_effs.el_eff,   pt, eta));
                met_up_w.fake_weight    += computeLeptonWeight(isTight, r,
                    getEfficiency(met_up_effs.el_eff,    pt, eta));
                met_down_w.fake_weight  += computeLeptonWeight(isTight, r,
                    getEfficiency(met_down_effs.el_eff,  pt, eta));
                comp_up_w.fake_weight   += computeLeptonWeight(isTight, r,
                    getEfficiency(comp_up_effs.el_eff,   pt, eta));
                comp_down_w.fake_weight += computeLeptonWeight(isTight, r,
                    getEfficiency(comp_down_effs.el_eff, pt, eta));
                real_up_w.fake_weight   += computeLeptonWeight(isTight, r,
                    getEfficiency(real_up_effs.el_eff,   pt, eta));
                real_down_w.fake_weight += computeLeptonWeight(isTight, r,
                    getEfficiency(real_down_effs.el_eff, pt, eta));
                stat_up_w.fake_weight   += computeLeptonWeight(isTight, r,
                    getEfficiency(stat_up_effs.el_eff,   pt, eta));
                stat_down_w.fake_weight += computeLeptonWeight(isTight, r,
                    getEfficiency(stat_down_effs.el_eff, pt, eta));
            }
        }
        outTree->Fill();
    }

    outFile->cd();
    outTree->Write();
    outFile->Close();
    std::cout << "✓ Wrote " << outputFile << "\n";

    delete fakeFile;
    delete realFile;
    return 0;
}

// ---------------------------------------------------------------------------
// Default entry point — DATA only.
// MC is no longer cloned here: the closure test computes AMM weights for MC
// on-the-fly inside studies/mc_closure_fake.C.
// ---------------------------------------------------------------------------
void AMM_ElectronOnly()
{
    const std::vector<std::string> dataInputs = {
        "/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data15.root",
        "/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_a/data16.root",
        "/eos/user/s/sabdelha/freshness/hnl-recon/freshntuples_d/data17.root",
        "/eos/user/s/sabdelha/freshness/hnl-recon/fresh_ntuples/data18.root",
    };
    AMM_runOne(dataInputs, "data_with_electron_fake_weights.root", "DATA");
}
