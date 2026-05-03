// =============================================================================
// FILE:    AMM/dump_syst_mapping.C
// PURPOSE: Initialize CP::AsymptMatrixTool with fake_eff_input.root and print
//          the mapping between the generic FAKEBKG_*_VARn names that appear
//          in the output ntuple branches and the human-readable systematic
//          labels (METdep, Composition, PromptSubtr, StatExtra, ...).
//
//          Run once after each change to the efficiency input file. Pipe the
//          output to a text file and use it to rename branches downstream.
//
// USAGE (with AnalysisBase set up):
//   root -b -q -l 'dump_syst_mapping.C("MyAnalysis/fake_eff_input.root")'
// =============================================================================
#include "AsgTools/AnaToolHandle.h"
#include "FakeBkgTools/IFakeBkgTool.h"
#include "FakeBkgTools/IFakeBkgSystDescriptor.h"
#include "PATInterfaces/SystematicVariation.h"
#include "PATInterfaces/SystematicSet.h"

#include <iostream>
#include <vector>
#include <string>

void dump_syst_mapping(const char* effFile = "MyAnalysis/fake_eff_input.root")
{
    asg::AnaToolHandle<CP::IFakeBkgTool> tool("CP::AsymptMatrixTool/AMM");
    tool.setProperty("InputFiles",   std::vector<std::string>{effFile}).ignore();
    tool.setProperty("Process",      ">=1F").ignore();
    tool.setProperty("Selection",    "=1T").ignore();
    tool.setProperty("EnergyUnit",   "GeV").ignore();
    tool.setProperty("ConvertWhenMissing", true).ignore();
    tool.setProperty("TightDecoration", "TightForFakeBkgCalculation,as_char").ignore();
    if (tool.initialize().isFailure()) {
        std::cerr << "Tool initialize failed\n";
        return;
    }

    auto sysvars = tool->affectingSystematics();
    std::cout << "\n========================================================\n";
    std::cout << "  FakeBkgTools systematic-variation mapping\n";
    std::cout << "  (input: " << effFile << ")\n";
    std::cout << "========================================================\n";
    std::cout << "  variation name              | description\n";
    std::cout << "--------------------------------------------------------\n";
    for (auto& sv : sysvars) {
        std::string desc = tool->getSystDescriptor()->getUncertaintyDescription(sv);
        std::cout << "  " << sv.name();
        for (size_t k = sv.name().size(); k < 28; ++k) std::cout << ' ';
        std::cout << "| " << desc << "\n";
    }
    std::cout << "========================================================\n\n";
}
