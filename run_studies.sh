#!/bin/bash
# =============================================================================
# FILE:    run_studies.sh
# PURPOSE: Master driver — runs the full fake-electron background study from
#          scratch in six phases:
#
#          PHASE 1 — Efficiency chain (ε_f, ε_r maps + systematics)
#          PHASE 2 — AMM per-event weights (data ntuple)
#          PHASE 3 — Main studies (cutflows, closure, yields, LaTeX tables)
#          PHASE 4 — Cross-check diagnostics (in studies/checks/)
#          PHASE 5 — Smoothing study (in studies/smoothing_test/)
#          PHASE 6 — Plotting (1D/2D rates, syst variations, BG composition,
#                              pre/post-fit, AMM validation, region plots)
#
#          Each macro logs to <subdir>/outputs/logs/<name>.log.  Failures
#          DO NOT abort the script (`set -e` deliberately omitted) so a
#          single bad macro doesn't block the rest.
#
# USAGE:   cd /eos/.../fake_study_final
#          chmod +x run_studies.sh
#          ./run_studies.sh                                     # full chain
#          ./run_studies.sh --skip-eff                          # skip phase 1
#          ./run_studies.sh --only-plots                        # only phase 6
#          nohup ./run_studies.sh > run_studies.log 2>&1 &      # background
# =============================================================================
set -u

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

# ---------- options ----------
SKIP_EFF=0; SKIP_AMM=0; SKIP_STUDIES=0; SKIP_CHECKS=0; SKIP_SMOOTH=0; SKIP_PLOTS=0
ONLY=""
for arg in "$@"; do
    case "$arg" in
        --skip-eff)     SKIP_EFF=1 ;;
        --skip-amm)     SKIP_AMM=1 ;;
        --skip-studies) SKIP_STUDIES=1 ;;
        --skip-checks)  SKIP_CHECKS=1 ;;
        --skip-smooth)  SKIP_SMOOTH=1 ;;
        --skip-plots)   SKIP_PLOTS=1 ;;
        --only-eff)     ONLY=eff      ;;
        --only-amm)     ONLY=amm      ;;
        --only-studies) ONLY=studies  ;;
        --only-checks)  ONLY=checks   ;;
        --only-smooth)  ONLY=smooth   ;;
        --only-plots)   ONLY=plots    ;;
        *) echo "Unknown option: $arg" >&2; exit 2 ;;
    esac
done
should_run() {
    [ -n "$ONLY" ] && [ "$ONLY" != "$1" ] && return 1
    return 0
}

START_TS=$(date +%s)

banner() {
    echo
    echo "=============================================================="
    echo "  $1"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "=============================================================="
}

# Run a macro from a sub-directory; log to <subdir>/outputs/logs/<name>.log.
# Args: <subdir>  <macro.C>  [pretty-name]
run_one() {
    local subdir="$1"
    local macro="$2"
    local label="${3:-$macro}"
    local name="${macro%.C}"
    local logdir="${SCRIPT_DIR}/${subdir}/outputs/logs"
    mkdir -p "$logdir"
    local log="${logdir}/${name}.log"

    echo "----- [RUN ] ${label}  (${subdir}/${macro})"
    local t0=$(date +%s)
    pushd "${SCRIPT_DIR}/${subdir}" >/dev/null
    if root -b -q -l "${macro}" >"${log}" 2>&1; then
        local rc=0
    else
        local rc=$?
    fi
    popd >/dev/null
    local t1=$(date +%s)
    local dt=$((t1 - t0))
    if [ "$rc" -eq 0 ]; then
        echo "      [DONE] ${label}  (${dt}s)"
    else
        echo "      [FAIL] ${label}  (rc=${rc}, ${dt}s)  log: ${log}"
    fi
}

# =============================================================================
# PHASE 1 — Efficiency chain
# =============================================================================
if [ $SKIP_EFF -eq 0 ] && should_run eff; then
    banner "PHASE 1 — Efficiency chain"
#    run_one "efficiency" "MC_eff.C"                       "MC prompt histograms"
#    run_one "efficiency" "data_eff.C"                     "Data tight/loose histograms"
#    run_one "efficiency" "real_eff.C"                     "Real efficiency epsilon_r (MC truth-prompt)"
#    run_one "efficiency" "calculate_fake_rates.C"         "Compute epsilon_f map"
#    run_one "efficiency" "build_composition_systematic.C" "Composition NP (truth-source reweighting)"
 #   run_one "efficiency" "build_fake_systematics.C"       "Total systematic variations"
  #  run_one "efficiency" "build_uncertainties.C"          "Uncertainty breakdown table"
  #  run_one "efficiency" "export_for_FakeBkgTools.C"      "Export epsilon_f / epsilon_r for FakeBkgTools"
fi

# =============================================================================
# PHASE 2 — Per-event AMM weights
# =============================================================================
if [ $SKIP_AMM -eq 0 ] && should_run amm; then
    banner "PHASE 2 — AMM per-event weights"
  #  run_one "AMM" "AMM_ElectronOnly.C" "Apply AMM to data (writes data_with_electron_fake_weights.root)"
fi

# =============================================================================
# PHASE 3 — Main studies (cutflows, closure, yields)
# =============================================================================
if [ $SKIP_STUDIES -eq 0 ] && should_run studies; then
    banner "PHASE 3 — Main studies"
 #   run_one "studies" "MC_cutflow.C"                  "MC cutflow"
 #   run_one "studies" "data_cutflow.C"                "Data cutflow"
 #   run_one "studies" "fake_cutflow.C"                "Fake cutflow (AMM-weighted)"
 #   run_one "studies" "fake_yields_correct.C"         "Fake yield table"
  #  run_one "studies" "mc_closure_fake.C"             "MC closure with diagnostics"
  #  run_one "studies" "met_systematic_investigation.C" "MET-split significance test"
  #  run_one "studies" "bkg_composition.C"             "BG composition per region"
  #  run_one "studies" "print_yields_latex.C"          "Yields LaTeX table"
fi

# =============================================================================
# PHASE 4 — Cross-checks
# =============================================================================
if [ $SKIP_CHECKS -eq 0 ] && should_run checks; then
    banner "PHASE 4 — Cross-checks"
  #  run_one "studies/checks" "fake_factor_compare.C"      "Fake factor vs Matrix method"
  #  run_one "studies/checks" "justify_pt_cap.C"           "pT-cap diagnostic"
  #  run_one "studies/checks" "sf_distributions.C"         "Charge-misID & ECIDS SF"
  #  run_one "studies/checks" "truth_source_breakdown.C"   "Truth source breakdown"
  #  run_one "studies/checks" "data_mc_loose_ratio.C"      "data / MC_prompt loose ratio"
  #  run_one "studies/checks" "d0sig_distribution.C"       "|d0sig| of loose fakes"
  #  run_one "studies/checks" "dR_jet_distribution.C"      "DeltaR(e, jet) of loose fakes"
  #  run_one "studies/checks" "fake_source_diagnostics.C"  "Fake source vs all variables"
  #  run_one "studies/checks" "met_split_test.C"           "MET-split multijet test"
  #  run_one "studies/checks" "multijet_test.C"            "Anti-iso muon multijet test"
  #  run_one "studies/checks" "eps_f_vs_met_cut.C"         "epsilon_f vs MET cut scan"
fi

# =============================================================================
# PHASE 5 — Smoothing study
# =============================================================================
if [ $SKIP_SMOOTH -eq 0 ] && should_run smooth; then
    banner "PHASE 5 — Smoothing study"
 #   run_one "studies/smoothing_test" "smooth_fake_eff.C"      "Smooth epsilon_f map (cubic spline)"
 #   run_one "studies/smoothing_test" "amm_compare_smoothed.C" "AMM raw vs smoothed (data16)"
  #  run_one "studies/smoothing_test" "plot_compare.C"         "Comparison plot"
fi

# =============================================================================
# PHASE 6 — Plotting (final ATLAS-Internal plots)
# =============================================================================
if [ $SKIP_PLOTS -eq 0 ] && should_run plots; then
    banner "PHASE 6 — Plotting"
  #  run_one "plotting" "plot_eff_maps.C"            "2D epsilon_f and epsilon_r maps"
  #  run_one "plotting" "plot_real.C"                "Real efficiency plots"
  #  run_one "plotting" "plot_rates_1d.C"            "1D fake & real efficiencies"
  #  run_one "plotting" "plot_fake_rates_sliced.C"   "Fake rate sliced views"
    run_one "plotting" "fake_real_CRs.C"            "Fake/real CR data-MC plots"
    run_one "plotting" "fake_comp.C"                "Fake composition plot"
    run_one "plotting" "plot_cr_num_den.C"          "CR numerator/denominator"
    run_one "plotting" "presel_VR_plots.C"          "Preselection / VR plots"
    run_one "plotting" "amm_validation_plots.C"     "AMM validation (closure shapes)"
    run_one "plotting" "VR_CR_fake_sys.C"           "Fake systematic in VR/CR"
    run_one "plotting" "plot_prefit.C"              "Pre-fit stacked distributions"
    run_one "plotting" "plot_postfit.C"             "Post-fit (mu_tt=mu_Zttau=1.0)"
fi

# =============================================================================
# Summary
# =============================================================================
END_TS=$(date +%s)
TOTAL=$((END_TS - START_TS))
banner "All phases finished — total wall time: ${TOTAL}s"

echo "Output PDFs across all phases:"
find efficiency/outputs studies/outputs studies/checks/outputs \
     studies/smoothing_test/outputs plotting/outputs \
     -maxdepth 3 -name '*.pdf' 2>/dev/null | sort

echo
echo "Logs:"
echo "  efficiency/outputs/logs/"
echo "  AMM/outputs/logs/"
echo "  studies/outputs/logs/"
echo "  studies/checks/outputs/logs/"
echo "  studies/smoothing_test/outputs/logs/"
echo "  plotting/outputs/logs/"
echo
echo "To copy all plots into the slides figures folder:"
echo "  bash slides/collect_figures.sh"
