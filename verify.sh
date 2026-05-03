#!/bin/bash
# =============================================================================
# verify.sh — Compile-check every macro in the framework without running it.
#
# Runs `root -b -l -q -e '.L <macro>+'` for each .C / .c file in
# efficiency/, studies/, plotting/, AMM/, plus the root-level
# convert_fake_ntuple.C. Prints PASS/FAIL per file and exits non-zero if any
# macro fails to compile.
#
# Usage:    bash verify.sh
# Requires: ROOT in $PATH (e.g. after `lsetup root` on lxplus).
# =============================================================================
set -u
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "$SCRIPT_DIR"

if ! command -v root >/dev/null 2>&1; then
    echo "ERROR: root not in PATH. Run 'lsetup root' (or equivalent) first."
    exit 2
fi

PASS=0; FAIL=0; FAIL_LIST=()
LOGDIR="$SCRIPT_DIR/.verify_logs"
rm -rf "$LOGDIR"; mkdir -p "$LOGDIR"

# Each entry: "<subdir> <basename>".  Empty subdir = project root.
FILES=(
    " convert_fake_ntuple.C"
    "efficiency MC_eff.C"
    "efficiency data_eff.C"
    "efficiency real_eff.C"
    "efficiency calculate_fake_rates.C"
    "efficiency build_fake_systematics.C"
    "efficiency build_uncertainties.C"
    "efficiency process_systematics.C"
    "efficiency build_composition_systematic.C"
    "efficiency export_for_FakeBkgTools.C"
    "studies MC_cutflow.C"
    "studies data_cutflow.C"
    "studies fake_cutflow.C"
    "studies fake_yields_correct.C"
    "studies print_yields_latex.C"
    "studies mc_closure_fake.C"
    "studies met_systematic_investigation.C"
    "studies systematics_evidence.C"
    "studies bkg_composition.C"
    "plotting plot_eff_maps.C"
    "plotting plot_rates_1d.C"
    "plotting plot_fake_rates_sliced.C"
    "plotting plotSubtractedFakeElectronEfficiencies.c"
    "plotting plot_real.C"
    "plotting plot_cr_num_den.C"
    "plotting amm_validation_plots.C"
    "plotting VR_CR_fake_sys.C"
    "plotting fake_comp.C"
    "plotting fake_real_CRs.C"
    "plotting presel_VR_plots.C"
    "plotting plot_prefit.C"
    "plotting plot_postfit.C"
    "AMM AMM_ElectronOnly.C"
)

printf "%-55s  %s\n" "FILE" "RESULT"
printf -- "------------------------------------------------------------\n"

for entry in "${FILES[@]}"; do
    subdir="${entry%% *}"
    macro="${entry#* }"
    target_dir="$SCRIPT_DIR"
    [ -n "$subdir" ] && target_dir="$SCRIPT_DIR/$subdir"
    rel_path="${subdir:+$subdir/}$macro"

    if [ ! -f "$SCRIPT_DIR/$rel_path" ]; then
        printf "%-55s  %s\n" "$rel_path" "MISSING"
        FAIL=$((FAIL+1)); FAIL_LIST+=("$rel_path (missing)")
        continue
    fi

    log="$LOGDIR/$(echo "$rel_path" | tr '/' '_').log"
    # ACLiC compile check: load the macro with the + (or ++ for clean rebuild)
    # in the macro's own directory so its #include "../utils/..." resolves.
    ( cd "$target_dir" && root -b -l -q -e ".L $macro+" ) >"$log" 2>&1
    rc=$?
    # ROOT often returns 0 even on errors; also grep the log
    if [ $rc -eq 0 ] && ! grep -qE "Error|error:|fatal error" "$log"; then
        printf "%-55s  \033[32mPASS\033[0m\n" "$rel_path"
        PASS=$((PASS+1))
    else
        printf "%-55s  \033[31mFAIL\033[0m  (see %s)\n" "$rel_path" "$log"
        FAIL=$((FAIL+1)); FAIL_LIST+=("$rel_path")
    fi
done

printf -- "------------------------------------------------------------\n"
echo "Summary: $PASS passed, $FAIL failed"

if [ $FAIL -gt 0 ]; then
    echo ""
    echo "Failed files:"
    for f in "${FAIL_LIST[@]}"; do echo "  - $f"; done
    echo ""
    echo "Per-file logs in: $LOGDIR/"
    exit 1
fi

# Cleanup ACLiC build artifacts
find . -maxdepth 3 \( -name "*_C.so" -o -name "*_C.d" -o -name "*_C.pcm" \
    -o -name "*_c.so" -o -name "*_c.d" -o -name "*_c.pcm" \
    -o -name "*_C_ACLiC_dict_rdict.pcm" \) -delete 2>/dev/null

echo "All macros compile cleanly."
