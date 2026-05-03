#!/bin/bash
# =============================================================================
# FILE:    studies/run_diagnostics.sh
# PURPOSE: Run the five diagnostic macros sequentially, log each one's output
#          to a separate file, and continue past failures so a crash in one
#          doesn't block the rest.
#
# USAGE:   cd /eos/.../fake_study_final/studies
#          ./run_diagnostics.sh
#          # or to log everything to a single file as well:
#          ./run_diagnostics.sh 2>&1 | tee run_diagnostics.log
# =============================================================================

set -u   # fail on undefined vars; intentionally NOT using `set -e` so a single
         # macro failure doesn't stop the whole batch

# Always run from the directory where this script lives
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

mkdir -p outputs/logs
LOGDIR="outputs/logs"

# Macros in order: fastest first, so failures show up early
MACROS=(
    "fake_factor_compare.C"      # ~30s, just reads existing maps
    "truth_source_breakdown.C"   # ~10 min, MC chain
    "data_mc_loose_ratio.C"      # ~15 min, data + MC
    "d0sig_distribution.C"       # ~10 min, MC only
    "dR_jet_distribution.C"      # ~12 min, MC only
)

banner() {
    echo ""
    echo "=========================================================="
    echo "  $1"
    echo "  $(date '+%Y-%m-%d %H:%M:%S')"
    echo "=========================================================="
}

START_TS=$(date +%s)
banner "Diagnostic batch start"

for m in "${MACROS[@]}"; do
    name="${m%.C}"
    log="${LOGDIR}/${name}.log"
    banner "[RUN]  ${m}    log: ${log}"

    t0=$(date +%s)
    if root -b -q -l "${m}" >"${log}" 2>&1; then
        rc=0
    else
        rc=$?
    fi
    t1=$(date +%s)
    dt=$((t1 - t0))

    if [ "${rc}" -eq 0 ]; then
        echo "[DONE]  ${m}   (${dt}s)"
    else
        echo "[FAIL]  ${m}   (rc=${rc}, ${dt}s)  — see ${log}"
    fi
done

END_TS=$(date +%s)
TOTAL=$((END_TS - START_TS))
banner "All diagnostics finished — total wall time: ${TOTAL}s"

echo ""
echo "Outputs:"
ls -la outputs/*.pdf 2>/dev/null || true
ls -la outputs/*.root 2>/dev/null || true
echo ""
echo "Per-macro logs in ${LOGDIR}/"
