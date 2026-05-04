#!/bin/bash
# =============================================================================
# slides/collect_figures.sh
# Copies every PDF the framework produced into slides/figures/, renaming the
# few that need disambiguation so the LaTeX \includegraphics references match.
# =============================================================================
set -u
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
ROOT_DIR="${SCRIPT_DIR}/.."
FIGS="${SCRIPT_DIR}/figures"
mkdir -p "$FIGS"

cp_if_exists() {
    local src="$1"; local dst="$2"
    if [ -f "$src" ]; then
        cp -f "$src" "$dst" && echo "  + $(basename "$src") -> $(basename "$dst")"
    else
        echo "  - missing: $src"
    fi
}

cd "$ROOT_DIR"

# ----- ε_f map and systematic variations -----
cp_if_exists "efficiency/outputs/syst/fake_rate_nominal.pdf"            "$FIGS/fake_rate_nominal.pdf"
cp_if_exists "efficiency/outputs/syst/fake_rate_MET_up.pdf"             "$FIGS/fake_rate_MET_up.pdf"
cp_if_exists "efficiency/outputs/syst/fake_rate_MET_down.pdf"           "$FIGS/fake_rate_MET_down.pdf"
cp_if_exists "efficiency/outputs/syst/fake_rate_comp_up.pdf"            "$FIGS/fake_rate_comp_up.pdf"
cp_if_exists "efficiency/outputs/syst/fake_rate_comp_down.pdf"          "$FIGS/fake_rate_comp_down.pdf"
cp_if_exists "efficiency/outputs/syst/fake_rate_real_up.pdf"            "$FIGS/fake_rate_real_up.pdf"
cp_if_exists "efficiency/outputs/syst/fake_rate_real_down.pdf"          "$FIGS/fake_rate_real_down.pdf"
cp_if_exists "efficiency/outputs/syst/fake_rate_total_up.pdf"           "$FIGS/fake_rate_total_up.pdf"
cp_if_exists "efficiency/outputs/syst/fake_rate_total_down.pdf"         "$FIGS/fake_rate_total_down.pdf"
cp_if_exists "efficiency/outputs/syst/variation_pt_MET_central.pdf"     "$FIGS/variation_pt_MET_central.pdf"
cp_if_exists "efficiency/outputs/syst/variation_pt_comp_central.pdf"    "$FIGS/variation_pt_comp_central.pdf"
cp_if_exists "efficiency/outputs/syst/variation_pt_total_syst_central.pdf" "$FIGS/variation_pt_total_syst_central.pdf"

# ----- 1D rates from plotting/plot_rates_1d.C -----
# Fake (epsilon_f) vs each observable
cp_if_exists "plotting/outputs/02_1d_rates/fake/CR_pt.pdf"              "$FIGS/CR_pt_fake.pdf"
cp_if_exists "plotting/outputs/02_1d_rates/fake/CR_eta.pdf"             "$FIGS/CR_eta_fake.pdf"
cp_if_exists "plotting/outputs/02_1d_rates/fake/CR_d0sig.pdf"           "$FIGS/CR_d0sig_fake.pdf"
cp_if_exists "plotting/outputs/02_1d_rates/fake/CR_d0.pdf"              "$FIGS/CR_d0_fake.pdf"
cp_if_exists "plotting/outputs/02_1d_rates/fake/CR_met.pdf"             "$FIGS/CR_met_fake.pdf"
# MET-split versions (for syst motivation)
cp_if_exists "plotting/outputs/02_1d_rates/fake/CR_MET_low_pt.pdf"      "$FIGS/CR_MET_low_pt_fake.pdf"
cp_if_exists "plotting/outputs/02_1d_rates/fake/CR_MET_high_pt.pdf"     "$FIGS/CR_MET_high_pt_fake.pdf"
# Real efficiency
cp_if_exists "plotting/outputs/02_1d_rates/real/CR_pt.pdf"              "$FIGS/CR_pt_real.pdf"
cp_if_exists "plotting/outputs/02_1d_rates/real/CR_eta.pdf"             "$FIGS/CR_eta_real.pdf"
cp_if_exists "plotting/outputs/02_1d_rates/real/CR_d0sig.pdf"           "$FIGS/CR_d0sig_real.pdf"
cp_if_exists "plotting/outputs/02_1d_rates/real/CR_met.pdf"             "$FIGS/CR_met_real.pdf"

# ----- Eff-maps from plot_eff_maps.C and plot_real.C -----
cp_if_exists "plotting/outputs/eff_maps/fake_eff_2D.pdf"                "$FIGS/fake_eff_2D.pdf"
cp_if_exists "plotting/outputs/eff_maps/real_eff_2D.pdf"                "$FIGS/real_eff_2D.pdf"

# ----- CR vs SR fake-composition comparison (efficiency/plot_composition_cr_sr.C) -----
cp_if_exists "efficiency/outputs/composition_cr.pdf"                    "$FIGS/comp_cr.pdf"
cp_if_exists "efficiency/outputs/composition_sr.pdf"                    "$FIGS/comp_sr.pdf"

# ----- AMM validation in the VR (plotting/amm_validation_plots.C) -----
cp_if_exists "plotting/outputs/amm_validation/VR.pdf"                   "$FIGS/amm_validation_VR.pdf"

# ----- BG composition + studies -----
cp_if_exists "studies/outputs/bkg_composition.pdf"                      "$FIGS/bkg_composition.pdf"

# ----- Cross-checks -----
cp_if_exists "studies/checks/outputs/sf_vs_pt_eta.pdf"                  "$FIGS/sf_vs_pt_eta.pdf"
cp_if_exists "studies/checks/outputs/met_split_test.pdf"                "$FIGS/met_split_test.pdf"
cp_if_exists "studies/checks/outputs/eps_f_vs_met_cut.pdf"              "$FIGS/eps_f_vs_met_cut.pdf"
cp_if_exists "studies/checks/outputs/truth_source_breakdown.pdf"        "$FIGS/truth_source_breakdown.pdf"
cp_if_exists "studies/checks/outputs/data_mc_loose_ratio.pdf"           "$FIGS/data_mc_loose_ratio.pdf"
cp_if_exists "studies/checks/outputs/d0sig_distribution.pdf"            "$FIGS/d0sig_distribution.pdf"
cp_if_exists "studies/checks/outputs/dR_jet_distribution.pdf"           "$FIGS/dR_jet_distribution.pdf"
cp_if_exists "studies/checks/outputs/fake_factor_compare.pdf"           "$FIGS/fake_factor_compare.pdf"
cp_if_exists "studies/checks/outputs/justify_pt_cap.pdf"                "$FIGS/justify_pt_cap.pdf"
cp_if_exists "studies/checks/outputs/multijet_test.pdf"                 "$FIGS/multijet_test.pdf"
cp_if_exists "studies/checks/outputs/fake_source_diagnostics.pdf"       "$FIGS/fake_source_diagnostics.pdf"

# ----- Preselection plots (data/MC, no fake estimate) — used in backup -----
cp_if_exists "plotting/outputs/presel/mll.pdf"                          "$FIGS/presel_mll.pdf"
cp_if_exists "plotting/outputs/presel/mjj.pdf"                          "$FIGS/presel_mjj.pdf"
cp_if_exists "plotting/outputs/presel/l1_pt.pdf"                        "$FIGS/presel_l1pt.pdf"
cp_if_exists "plotting/outputs/presel/l2_pt.pdf"                        "$FIGS/presel_l2pt.pdf"
cp_if_exists "plotting/outputs/presel/met.pdf"                          "$FIGS/presel_met.pdf"

# ----- Smoothing -----
cp_if_exists "studies/smoothing_test/outputs/smoothing_diagnostic.pdf"  "$FIGS/smoothing_diagnostic.pdf"
cp_if_exists "studies/smoothing_test/outputs/compare_smoothing.pdf"     "$FIGS/compare_smoothing.pdf"

# ----- Pre/post-fit (only matter once a fit is run) -----
cp_if_exists "plotting/outputs/06_prefit/prefit_SR.pdf"                 "$FIGS/prefit_SR.pdf"
cp_if_exists "plotting/outputs/07_postfit/postfit_SR.pdf"               "$FIGS/postfit_SR.pdf"
cp_if_exists "plotting/outputs/06_prefit/prefit_VR.pdf"                 "$FIGS/prefit_VR.pdf"
cp_if_exists "plotting/outputs/07_postfit/postfit_VR.pdf"               "$FIGS/postfit_VR.pdf"

echo
echo "Done.  Figures in: $FIGS"
ls -1 "$FIGS" 2>/dev/null
