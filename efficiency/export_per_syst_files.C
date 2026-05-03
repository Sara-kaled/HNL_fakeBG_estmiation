// =============================================================================
// FILE:    efficiency/export_per_syst_files.C
// PURPOSE: Produce one default-format ROOT file per systematic variation,
//          where the syst-varied fake efficiency is written as the NOMINAL
//          histogram. Lets you run N independent FakeBkgCalculator blocks
//          (one per syst) and get one named branch per syst in the output
//          ntuple — no shadow file or runSystematics needed.
//
// Inputs (relative to the current dir):
//   outputs/fake_eff_with_systematics.root
//   outputs/real_eff_histograms.root
//
// Outputs (one file per syst, named by setupName suffix):
//   outputs/fake_eff_input_nominal.root     ε_f = fake_rate_nominal
//   outputs/fake_eff_input_METup.root       ε_f = fake_rate_MET_up
//   outputs/fake_eff_input_METdown.root     ε_f = fake_rate_MET_down
//   outputs/fake_eff_input_compup.root      ε_f = fake_rate_comp_up
//   outputs/fake_eff_input_compdown.root    ε_f = fake_rate_comp_down
//   outputs/fake_eff_input_realup.root      ε_f = fake_rate_real_up
//   outputs/fake_eff_input_realdown.root    ε_f = fake_rate_real_down
//   outputs/fake_eff_input_statup.root      ε_f = fake_rate_statUp
//   outputs/fake_eff_input_statdown.root    ε_f = fake_rate_statDown
//
// Each file contains exactly the four histograms the FakeBkgTool needs:
//   RealEfficiency2D_el_pt_eta   (nominal ε_r — same in every file)
//   FakeEfficiency2D_el_pt_eta   (the syst-varied ε_f — different per file)
//   RealEfficiency_mu_pt         (1-bin dummy ≈ 1)
//   FakeEfficiency_mu_pt         (1-bin dummy ≈ 0)
// =============================================================================
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TSystem.h>
#include <iostream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
static void writeOne(TH2D* hRealEff, TH2D* hFakeEff, const std::string& outName)
{
    TFile* fOut = TFile::Open(outName.c_str(), "RECREATE");
    if (!fOut || fOut->IsZombie()) {
        std::cerr << "ERROR: cannot create " << outName << "\n";
        return;
    }

    hRealEff->Write("RealEfficiency2D_el_pt_eta");
    hFakeEff->Write("FakeEfficiency2D_el_pt_eta");

    TH1D* hMuReal = new TH1D("RealEfficiency_mu_pt",
                             ";p_{T}^{#mu} [GeV];#varepsilon_{r}^{#mu}",
                             1, 0., 1.e6);
    hMuReal->SetBinContent(1, 0.999);
    hMuReal->SetBinError  (1, 0.);
    hMuReal->Write();

    TH1D* hMuFake = new TH1D("FakeEfficiency_mu_pt",
                             ";p_{T}^{#mu} [GeV];#varepsilon_{f}^{#mu}",
                             1, 0., 1.e6);
    hMuFake->SetBinContent(1, 0.001);
    hMuFake->SetBinError  (1, 0.);
    hMuFake->Write();

    fOut->Close();
    std::cout << "  wrote " << outName << "\n";
}

// ---------------------------------------------------------------------------
void export_per_syst_files()
{
    TFile* fFake = TFile::Open("outputs/fake_eff_with_systematics.root", "READ");
    TFile* fReal = TFile::Open("outputs/real_eff_histograms.root", "READ");
    if (!fFake || fFake->IsZombie() || !fReal || fReal->IsZombie()) {
        std::cerr << "ERROR: missing input file(s) under outputs/\n";
        return;
    }

    auto* hRealTight = (TH2D*)fReal->Get("h2_tight_eta_pt_CR");
    auto* hRealLoose = (TH2D*)fReal->Get("h2_loose_eta_pt_CR");
    if (!hRealTight || !hRealLoose) {
        std::cerr << "ERROR: real-eff histograms not found\n"; return;
    }
    auto* hRealEff = (TH2D*)hRealTight->Clone("hRealEff_tmp");
    hRealEff->Divide(hRealLoose);

    // (variant key, source histogram name in fake_eff_with_systematics.root)
    const std::vector<std::pair<std::string,std::string>> variants = {
        {"nominal",  "fake_rate_nominal"},
        {"METup",    "fake_rate_MET_up"},
        {"METdown",  "fake_rate_MET_down"},
        {"compup",   "fake_rate_comp_up"},
        {"compdown", "fake_rate_comp_down"},
        {"realup",   "fake_rate_real_up"},
        {"realdown", "fake_rate_real_down"},
        {"statup",   "fake_rate_statUp"},
        {"statdown", "fake_rate_statDown"},
    };

    gSystem->mkdir("outputs", true);
    std::cout << "\n=== Writing per-systematic FakeBkgTool input files ===\n";

    for (auto& v : variants) {
        const std::string& key = v.first;
        const std::string& src = v.second;

        auto* hFake = (TH2D*)fFake->Get(src.c_str());
        if (!hFake) {
            std::cerr << "  skip '" << key << "': histogram '" << src
                      << "' not found in fake_eff_with_systematics.root\n";
            continue;
        }
        const std::string out = "outputs/fake_eff_input_" + key + ".root";
        writeOne(hRealEff, hFake, out);
    }

    fFake->Close();
    fReal->Close();
    std::cout << "\n✓ Done. Stage outputs/fake_eff_input_*.root into MyAnalysis/data/\n";
}
