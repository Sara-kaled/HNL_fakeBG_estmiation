# Fake-background slides

Beamer presentation (`fake_background_study.tex`) summarizing the
HNL μe VBF fake-electron measurement.

## Structure

**Main body** — high-level walk-through (one idea per slide, plots/tables only):
1. Introduction (HNL signal, backgrounds, what is a fake)
2. Object definitions (electrons, muons, jets, MET)
3. Region table + triggers
4. The Asymptotic Matrix Method
5. ε_f and ε_r measurements (2D maps + 1D vs pT, |η|)
6. Systematic uncertainties (table + variation plots)
7. AMM application + per-region yields + bkg composition
8. Summary

**Backup A — cross-checks** (each on one slide with a small conclusion box):
- Charge-misID SF distributions
- MET-split test
- ε_f vs MET cut scan
- Truth-source breakdown
- Data / MC_prompt loose ratio
- |d0σ| of loose fakes
- ΔR(e, jet) of loose fakes
- Fake factor vs Matrix method
- pT-cap justification (per-bin diagnostic)
- Multijet test (anti-iso μ)
- 8-panel fake-source diagnostics (vs pT, η, MET, mll, mjj, |Δηjj|, ΔR, |d0σ|)
- MC closure
- Smoothing study
- MET systematic significance

**Backup B — code reference**: directory map and step-by-step run sequence.

## How to build (Overleaf or local)

1. Upload `fake_background_study.tex` and the `figures/` folder to Overleaf.
2. Build with `pdflatex` (twice, for cross-references).

## Filling in the figures

Each `\includegraphics` line lists the source PDF written by the framework,
e.g.:

```
\includegraphics[width=0.55\textwidth]{figures/fake_rate_nominal.pdf}
% ↑ produced by efficiency/build_fake_systematics.C  →  outputs/syst/fake_rate_nominal.pdf
```

Copy each output PDF into `slides/figures/` keeping the basename the same.
A simple shell helper:

```bash
# from fake_study_final/
SLIDES_FIG=slides/figures
mkdir -p $SLIDES_FIG
cp efficiency/outputs/syst/fake_rate_nominal.pdf            $SLIDES_FIG/
cp efficiency/outputs/syst/fake_rate_MET_up.pdf             $SLIDES_FIG/
cp efficiency/outputs/syst/fake_rate_MET_down.pdf           $SLIDES_FIG/
cp efficiency/outputs/syst/fake_rate_comp_up.pdf            $SLIDES_FIG/
cp efficiency/outputs/syst/fake_rate_comp_down.pdf          $SLIDES_FIG/
cp efficiency/outputs/syst/variation_pt_total_syst_central.pdf $SLIDES_FIG/
cp plotting/outputs/02_1d_rates/fake/CR_pt.pdf              $SLIDES_FIG/CR_pt_fake.pdf
cp plotting/outputs/02_1d_rates/fake/CR_eta.pdf             $SLIDES_FIG/CR_eta_fake.pdf
cp plotting/outputs/02_1d_rates/real/CR_pt.pdf              $SLIDES_FIG/CR_pt_real.pdf
cp plotting/outputs/02_1d_rates/real/CR_eta.pdf             $SLIDES_FIG/CR_eta_real.pdf
cp studies/outputs/bkg_composition.pdf                      $SLIDES_FIG/
cp studies/checks/outputs/sf_vs_pt_eta.pdf                  $SLIDES_FIG/
cp studies/checks/outputs/met_split_test.pdf                $SLIDES_FIG/
cp studies/checks/outputs/eps_f_vs_met_cut.pdf              $SLIDES_FIG/
cp studies/checks/outputs/truth_source_breakdown.pdf        $SLIDES_FIG/
cp studies/checks/outputs/data_mc_loose_ratio.pdf           $SLIDES_FIG/
cp studies/checks/outputs/d0sig_distribution.pdf            $SLIDES_FIG/
cp studies/checks/outputs/dR_jet_distribution.pdf           $SLIDES_FIG/
cp studies/checks/outputs/fake_factor_compare.pdf           $SLIDES_FIG/
cp studies/checks/outputs/justify_pt_cap.pdf                $SLIDES_FIG/
cp studies/checks/outputs/multijet_test.pdf                 $SLIDES_FIG/
cp studies/checks/outputs/fake_source_diagnostics.pdf       $SLIDES_FIG/
cp studies/smoothing_test/outputs/smoothing_diagnostic.pdf  $SLIDES_FIG/
```

If a figure file is missing on first compile, Beamer will complain — just
re-run the corresponding macro and copy the PDF in. The slide numbering
itself doesn't change.

## ATLAS style notes

- All plots produced by the framework already include `ATLAS Internal`,
  `√s = 13 TeV, 140 fb⁻¹` via `utils/AtlasStyle.h`.
- The Beamer theme is set to ATLAS-blue / red accents.
- Conclusion callouts use the `tcolorbox` package — small rounded boxes
  with a blue border, exactly as you asked for.
