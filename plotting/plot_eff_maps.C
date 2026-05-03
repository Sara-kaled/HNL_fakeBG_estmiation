// =============================================================================
// studies/01_efficiency_maps/plot_eff_maps.C
// PURPOSE: Wrapper that calls the existing plotSubtractedFakeElectronEfficiencies
//          (data - MC prompt subtraction) for fake eff maps, and plot_real for
//          real efficiency maps. Runs for all available regions.
// RUNS FROM: fake_study_final/studies/01_efficiency_maps/
// =============================================================================

#include <TSystem.h>
#include <vector>
// Load the existing macros
#include "plotSubtractedFakeElectronEfficiencies.c"
#include "plot_real.C"

void plot_eff_maps()
{
    // Correct relative paths from this directory
    const char* mc_fake  = "../efficiency/outputs/MC_histograms.root";
    const char* mc_real  = "../efficiency/outputs/real_eff_histograms.root";
    const char* data     = "../efficiency/outputs/Data_histograms.root";

    // Output directories
    const char* outFake = "outputs/01_eff_maps/fake";
    const char* outReal = "outputs/01_eff_maps/real";

    gSystem->mkdir("outputs/01_eff_maps/fake", true);
    gSystem->mkdir("outputs/01_eff_maps/real",  true);

    // Regions to plot
    std::vector<const char*> regions = {"CR", "CR_MET_high", "CR_MET_low"};

    printf("\n=== Fake efficiency maps (data - MC subtraction) ===\n");
    for (const char* reg : regions) {
        printf("  Region: %s\n", reg);
        plotSubtractedFakeElectronEfficiencies(
            mc_fake, data,
            reg,
            true,    // plot2D
            false,   // logx
            0.0,     // yMin
            0.5,     // yMax — fake rate typically < 0.3
            outFake,
            "pdf"
        );
    }

    printf("\n=== Real efficiency maps (MC tight/loose ratio) ===\n");
    for (const char* reg : regions) {
        printf("  Region: %s\n", reg);
        plot_real(
            mc_real,
            reg,
            true,    // plot2D
            false,   // logx
            0.0,     // yMin
            1.0,     // yMax
            outReal,
            "pdf"
        );
    }

    printf("\n=== plot_eff_maps done ===\n");
    printf("  Fake maps: outputs/01_eff_maps/fake/\n");
    printf("  Real maps: outputs/01_eff_maps/real/\n");
}

