// =============================================================================
// FILE:    AMM/rename_syst_branches.C
// PURPOSE: Rename the generic FakeBkgTools systematic branches
//          (FAKEBKG_SYST_VAR0, ..._VAR1, ...) into human-readable names
//          (METdep, Composition, ...) in a copy of the framework ntuple.
//
//          The mapping is the hardcoded `kSystMap` table below — edit it once
//          to match the output of dump_syst_mapping.C, then run.
//
// USAGE:
//   root -b -q -l 'rename_syst_branches.C("framework_out.root","framework_out_renamed.root")'
//
// HOW IT WORKS:
//   - Opens the input tree
//   - For each branch whose name appears as a key in kSystMap, sets its alias
//     to the desired new name AND writes a renamed copy in the output tree
//   - All other branches are cloned untouched
// =============================================================================
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <iostream>
#include <map>
#include <string>

// ---------------------------------------------------------------------------
// MAPPING TABLE — edit to match dump_syst_mapping.C output
//   key   = original branch suffix as written by FakeBkgTool
//   value = desired suffix in the output tree
// ---------------------------------------------------------------------------
// Map matches the histograms currently in fake_eff_input.root:
//   __METdep, __Composition, __PromptSubtr
// Add VAR3 → StatExtra below if you regenerate with stat_up/down histograms.
static const std::map<std::string, std::string> kSystMap = {
    {"FAKEBKG_SYST_VAR0__1up",   "METdep__1up"},
    {"FAKEBKG_SYST_VAR0__1down", "METdep__1down"},
    {"FAKEBKG_SYST_VAR1__1up",   "Composition__1up"},
    {"FAKEBKG_SYST_VAR1__1down", "Composition__1down"},
    {"FAKEBKG_SYST_VAR2__1up",   "PromptSubtr__1up"},
    {"FAKEBKG_SYST_VAR2__1down", "PromptSubtr__1down"},
};

// Branch prefix that holds the fake weight (everything else is left alone)
static const std::string kPrefix = "weight_fake_electronAMM_";

// ---------------------------------------------------------------------------
static std::string renamed(const std::string& origName)
{
    if (origName.compare(0, kPrefix.size(), kPrefix) != 0) return origName;
    const std::string suffix = origName.substr(kPrefix.size());
    auto it = kSystMap.find(suffix);
    if (it == kSystMap.end()) return origName;
    return kPrefix + it->second;
}

void rename_syst_branches(const char* inFile  = "framework_out.root",
                          const char* outFile = "framework_out_renamed.root",
                          const char* treeName = "reco")
{
    TFile* fi = TFile::Open(inFile, "READ");
    if (!fi || fi->IsZombie()) { std::cerr << "Cannot open " << inFile << "\n"; return; }
    TTree* tin = (TTree*)fi->Get(treeName);
    if (!tin) { std::cerr << "Tree '" << treeName << "' not in " << inFile << "\n"; return; }

    TFile* fo = TFile::Open(outFile, "RECREATE");
    fo->SetCompressionLevel(1);
    TTree* tout = tin->CloneTree(-1, "fast");   // clone all events + branches
    tout->SetDirectory(fo);

    // Walk the cloned branches and rename the ones in the map.
    int nRenamed = 0;
    TIter next(tout->GetListOfBranches());
    while (auto* b = (TBranch*)next()) {
        const std::string oldName = b->GetName();
        const std::string newName = renamed(oldName);
        if (newName == oldName) continue;

        b->SetName(newName.c_str());
        b->SetTitle(newName.c_str());
        // Rename the leaf too so TTree::Draw / Scan use the new name
        if (auto* leaf = (TLeaf*)b->GetListOfLeaves()->First()) {
            leaf->SetName(newName.c_str());
            leaf->SetTitle(newName.c_str());
        }
        std::cout << "  " << oldName << "  →  " << newName << "\n";
        ++nRenamed;
    }

    fo->cd();
    tout->Write("", TObject::kOverwrite);
    fo->Close();
    fi->Close();

    std::cout << "\n✓ Wrote " << outFile << "  (" << nRenamed << " branches renamed)\n";
}
