# Fake-Electron Background Framework — HNL μe VBF

Data-driven fake-electron background estimation for the HNL μe VBF analysis,
based on the Asymptotic Matrix Method (AMM). Full Run-2 (140 fb⁻¹) at
√s = 13 TeV.

## Structure

```
fake_study_final/
├── utils/
│   ├── SelectionUtils.h       Central EventData + every region/cut helper
│   └── AtlasStyle.h            Plotting style + ATLAS palette
├── efficiency/                 ε_f and ε_r measurements + systematics
│   ├── MC_eff.C, data_eff.C, real_eff.C, calculate_fake_rates.C
│   ├── build_composition_systematic.C
│   ├── build_fake_systematics.C, build_uncertainties.C
│   └── export_for_FakeBkgTools.C
├── AMM/
│   └── AMM_ElectronOnly.C      Per-event AMM weight ntuple producer
├── studies/                    Cutflows, closure, yields
│   ├── MC_cutflow.C, data_cutflow.C, fake_cutflow.C
│   ├── mc_closure_fake.C, met_systematic_investigation.C
│   ├── bkg_composition.C, fake_yields_correct.C
│   ├── checks/                 Cross-check macros (one per question)
│   └── smoothing_test/         Optional ε_f smoothing study
├── plotting/                   ATLAS-style plot producers
├── slides/                     Beamer presentation + collect_figures.sh
└── run_studies.sh              Master driver — runs everything in 6 phases
```

## Run

Full chain (≈ 1–2 h on lxplus):
```bash
./run_studies.sh
```

Per-phase:
```bash
./run_studies.sh --only-eff       # ε_f, ε_r, systematics
./run_studies.sh --only-amm       # AMM weights ntuple
./run_studies.sh --only-studies   # cutflows + closure + yields
./run_studies.sh --only-checks    # cross-checks
./run_studies.sh --only-smooth    # smoothing study
./run_studies.sh --only-plots     # final plots
```

Skip phases:
```bash
./run_studies.sh --skip-eff --skip-amm   # plotting / checks only
```

Each macro logs to `<subdir>/outputs/logs/<name>.log`. Failures don't abort
the script, so a single bad macro doesn't block the rest.

## Slides

Beamer presentation in `slides/fake_background_study.tex` with the main
narrative + 14 backup cross-check slides (each with motivation, plot,
conclusion).

After `run_studies.sh` finishes:
```bash
bash slides/collect_figures.sh   # copy all output PDFs into slides/figures/
```

To deploy to Overleaf:
- Either upload `slides/` as a zip
- Or connect the GitHub repo to Overleaf via "Import from GitHub"

## Method (one-line)

The AMM weight per loose electron is

> w_T = ε_f (ε_r − 1) / (ε_r − ε_f) for tight events,
> w_L = ε_f ε_r / (ε_r − ε_f) for loose-only events

with ε_f measured in the SS μe fake CR and ε_r in the OS prompt CR.
Both are measured under the **muon-only trigger** to keep the electron
unbiased; the AMM is then **applied** in SR/VR with the full mu‖e trigger.

## Systematic uncertainties

| Source | Method | Output histogram |
|---|---|---|
| Statistical | Bin-wise Poisson | `fake_rate_statUp/Down` |
| MET dependence | Split CR at MET = 50 GeV | `fake_rate_MET_up/down` |
| Composition | Truth-source reweighting (CR → SR) | `fake_rate_comp_up/down` |
| Prompt subtraction | ±20% on prompt-MC norm | `fake_rate_real_up/down` |
| **Total** | quadrature sum, capped at 150% | `fake_rate_total_systUp/Down` |

## Inputs (not in this repo — on EOS)

- Data: `freshntuples_a/data{15,16}.root`, `freshntuples_d/data17.root`,
        `fresh_ntuples/data18.root`
- MC:   `fresh_ntuples/{ttbar,wjets,ztautau,singletop,diboson,higgs}/*.root`

The framework expects these paths. Adjust at the top of each input macro
if your storage layout differs.

## License & contact

Internal ATLAS analysis code. Contact: Sara Abdelhameed.
