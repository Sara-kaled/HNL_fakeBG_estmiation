// =============================================================================
// FILE:    AMM/filter_only_emu.C
// PURPOSE: Produce a copy of data_with_electron_fake_weights.root that only
//          contains events passing the YAML 'only_emu' selection:
//             JET_N baselineJvtLoose 20000 >= 2
//             EL_N  4500  == 1
//             MU_N  15000 == 1
//          No weight recalculation — every branch is preserved as-is.
//
// USAGE:
//   root -b -q -l 'filter_only_emu.C("data_with_electron_fake_weights.root",
//                                    "data_only_emu.root")'
// =============================================================================
#include <TFile.h>
#include <TTree.h>
#include "ROOT/RVec.hxx"
#include <iostream>

using ROOT::VecOps::RVec;

template <class V>
static int count_above(const V* v, double thresh_mev) {
    if (!v) return 0;
    int n = 0;
    for (size_t i = 0; i < v->size(); ++i) if ((*v)[i] > thresh_mev) ++n;
    return n;
}

static int count_good_jets(const RVec<float>* pt, const RVec<char>* jvt,
                           double pt_thresh_mev) {
    if (!pt || !jvt) return 0;
    int n = 0;
    const size_t N = std::min(pt->size(), jvt->size());
    for (size_t i = 0; i < N; ++i)
        if ((*pt)[i] > pt_thresh_mev && (*jvt)[i]) ++n;
    return n;
}

void filter_only_emu(const char* inFile  = "data_with_electron_fake_weights.root",
                     const char* outFile = "data_only_emu.root")
{
    TFile* fi = TFile::Open(inFile, "READ");
    if (!fi || fi->IsZombie()) { std::cerr << "ERROR: cannot open " << inFile << "\n"; return; }
    TTree* tin = (TTree*)fi->Get("tree");
    if (!tin) { std::cerr << "ERROR: 'tree' not found in " << inFile << "\n"; return; }

    // Branches needed only for cut evaluation. Output tree clones everything.
    RVec<float>* el_pt          = nullptr;
    RVec<float>* mu_pt          = nullptr;
    RVec<float>* jet_pt         = nullptr;
    RVec<char>*  jet_jvt_loose  = nullptr;
    tin->SetBranchAddress("el_pt_NOSYS",                       &el_pt);
    tin->SetBranchAddress("mu_pt_NOSYS",                       &mu_pt);
    tin->SetBranchAddress("jet_pt_NOSYS",                      &jet_pt);
    tin->SetBranchAddress("jet_select_baselineJvtLoose_NOSYS", &jet_jvt_loose);

    TFile* fo = TFile::Open(outFile, "RECREATE");
    fo->SetCompressionLevel(1);
    TTree* tout = tin->CloneTree(0);
    tout->SetDirectory(fo);

    const Long64_t n = tin->GetEntries();
    std::cout << "Total entries: " << n << "\n";

    Long64_t nPass = 0;
    for (Long64_t i = 0; i < n; ++i) {
        tin->GetEntry(i);
        if (i % 100000 == 0) std::cout << "  entry " << i << " / " << n << "\n";

        if (count_good_jets(jet_pt, jet_jvt_loose, 20000.) < 2) continue;
        if (count_above(el_pt, 4500.) != 1)                     continue;
        if (count_above(mu_pt, 15000.) != 1)                    continue;

        tout->Fill();
        ++nPass;
    }

    fo->cd();
    tout->Write();
    fo->Close();
    fi->Close();

    std::cout << "\n✓ Wrote " << outFile << "\n";
    std::cout << "  Events passing 'only_emu':  " << nPass << " / " << n
              << "  (" << 100.*nPass/double(n ? n : 1) << "%)\n";
}
