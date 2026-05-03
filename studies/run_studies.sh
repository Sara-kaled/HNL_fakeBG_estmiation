#!/bin/bash
# =============================================================================
# FILE:    studies/run_studies.sh
# PURPOSE: Run all remaining studies and cross-checks AFTER the efficiency
#          chain (MC_eff, data_eff, real_eff, calculate_fake_rates,
#          build_composition_systematic, build_fake_systematics,
#          export_for_FakeBkgTools) has been run, AND after
#          met_systematic_investigation.C has been run separately.
#
#          Runs three groups in order:
#            1. main studies (closure, yields, cutflows, latex tables)
#            2. checks/ diagnostic macros
#            3. smoothing_test/ smoothing study
#
#          Continues past failures (no `set -e`), logs each macro's stdout/
#          stderr to outputs/logs/<macro>.log, and prints a summary at the end.
#
# USAGE:   cd /eos/.../fake_study_final/studies
#          ./run_studies.sh
#          # background + log everything to one file:
#          # nohup ./run_studies.sh >run_studies.log 2>&1 &
# =============================================================================
set -u

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

mkdir -p outputs/logs checks/outputs/logs smoothing_test/outputs/logs
START_TS=$(date +%s)

banner() {
    echo ""
    echo "=================================================================="
    echo "  $1"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "=================================================================="
}

# Run a macro from a sub-directory; log to <subdir>/outputs/logs/<name>.log.
# Args: <subdir>  <macro.C>  [optional: pretty-name]
run_one() {
    local subdir="$1"
    local macro="$2"
    local label="${3:-$macro}"
    local name="${macro%.C}"
    local logdir="${subdir}/outputs/logs"
    mkdir -p "$logdir"
    local log="${logdir}/${name}.log"

    banner "[RUN] ${label} (in ${subdir}/)"
    local t0=$(date +%s)
    pushd "$subdir" >/dev/null
    if root -b -q -l "${macro}" >"../${log}" 2>&1; then
        local rc=0
    else
        local rc=$?
    fi
    popd >/dev/null
    local t1=$(date +%s)
    local dt=$((t1 - t0))

    if [ "$rc" -eq 0 ]; then
        echo "[DONE] ${label}  (${dt}s)"
    else
        echo "[FAIL] ${label}  (rc=${rc}, ${dt}s)  see ${log}"
    fi
}

# -----------------------------------------------------------------------------
# 1. Main studies (in studies/)
# -----------------------------------------------------------------------------
banner "PHASE 1 — main studies"

# closure first (it's slow; nice to have it running early so logs are ready)
run_one "."  "mc_closure_fake.C"        "MC closure"

# yields and cutflows
run_one "."  "MC_cutflow.C"             "MC cutflow"
run_one "."  "data_cutflow.C"           "data cutflow"
run_one "."  "fake_cutflow.C"           "fake cutflow"
run_one "."  "fake_yields_correct.C"    "fake yields"
run_one "."  "bkg_composition.C"        "bkg composition"
run_one "."  "print_yields_latex.C"     "yields LaTeX table"

# -----------------------------------------------------------------------------
# 2. Cross-check diagnostics (in studies/checks/)
# -----------------------------------------------------------------------------
banner "PHASE 2 — cross-checks"

# Fast first, slow last
run_one "checks" "fake_factor_compare.C"        "AMM vs fake-factor"
run_one "checks" "justify_pt_cap.C"             "pT-cap justification"
run_one "checks" "sf_distributions.C"           "charge-misID SF distrib."
run_one "checks" "truth_source_breakdown.C"     "truth-source breakdown"
run_one "checks" "data_mc_loose_ratio.C"        "data/MC loose ratio"
run_one "checks" "d0sig_distribution.C"         "|d0sig| of fakes"
run_one "checks" "dR_jet_distribution.C"        "ΔR(e, jet) of fakes"
run_one "checks" "fake_source_diagnostics.C"    "fake-source vs all vars"
run_one "checks" "met_split_test.C"             "MET-split test"
run_one "checks" "multijet_test.C"              "multijet (anti-iso μ) test"
run_one "checks" "eps_f_vs_met_cut.C"           "ε_f vs MET cut scan"

# -----------------------------------------------------------------------------
# 3. Smoothing study (in studies/smoothing_test/)
# -----------------------------------------------------------------------------
banner "PHASE 3 — smoothing study"

run_one "smoothing_test" "smooth_fake_eff.C"      "smooth ε_f"
run_one "smoothing_test" "amm_compare_smoothed.C" "AMM raw vs smoothed (data16)"
run_one "smoothing_test" "plot_compare.C"         "plot comparison"

# -----------------------------------------------------------------------------
END_TS=$(date +%s)
TOTAL=$((END_TS - START_TS))
banner "All studies finished — total wall time: ${TOTAL}s"

echo ""
echo "Outputs:"
ls -la outputs/*.pdf outputs/*.root              2>/dev/null
ls -la checks/outputs/*.pdf checks/outputs/*.root 2>/dev/null
ls -la smoothing_test/outputs/*.pdf smoothing_test/outputs/*.root 2>/dev/null

echo ""
echo "Logs:"
echo "  outputs/logs/                  — main studies"
echo "  checks/outputs/logs/           — cross-checks"
echo "  smoothing_test/outputs/logs/   — smoothing study"
