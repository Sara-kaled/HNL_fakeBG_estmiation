# HNL Fake Lepton Background Framework

ATLAS Run 2 analysis framework for estimating fake/non-prompt lepton backgrounds in a
Heavy Neutral Lepton (HNL) search in the μe VBF final state.

---

## 1. Physics Context

### What is an HNL?

Heavy Neutral Leptons are proposed right-handed (sterile) neutrinos that mix with the
active SM neutrinos via a mixing parameter V. In the Type-I seesaw model, they are produced
at the LHC via Vector Boson Fusion (VBF): two quarks exchange a W/Z boson, radiating two
forward jets with large pseudorapidity separation (Δη_jj > 3) and large invariant mass
(m_jj > 400 GeV). The HNL decays via its mixing to produce a lepton pair.

The signal topology for this analysis:
- **Leading lepton (l1)**: muon from W decay, pT > 20 GeV
- **Subleading lepton (l2)**: electron from HNL decay, pT > 3.5 GeV
- **Two VBF jets**: high Δη, high m_jj, no b-jets
- **Charge**: opposite-sign (OS) μe pair

### Why does the μe channel have backgrounds?

Several SM processes produce a real or fake μe pair with VBF-like jets:

| Process | Contribution | Why it enters |
|---------|-------------|---------------|
| **W+jets** | Fake dominant | W→μν + b/c-quark jet → fake e |
| **ttbar** | Real + fake | t→Wb→μνb, second b→fake e, OR real OS pair |
| **Z→ττ** | Real | τ→μνν̄, τ→eνν̄ |
| **WW, WZ, ZZ** | Real | Prompt μe from diboson |
| **Single top** | Mixed | t→Wb, associated W→μ or e |
| **QCD multijet** | Fake | Both leptons are fakes from jet misidentification |

---

## 2. What is a "Fake" Lepton?

A **fake** (or **FNP — Fake/Non-Prompt**) electron is a reconstructed electron object that
does NOT originate from a prompt hard-scatter electroweak decay (W, Z, H, or signal).

### Sources of fake electrons in this analysis

| Source | Mechanism | Dominant process | Typical |d0| |
|--------|-----------|-----------------|---------|
| **Heavy-flavor (b-decay)** | b→eνX semi-leptonic | W+jets, ttbar | Large (~mm) — displaced vertex |
| **Charm-decay** | c→eνX semi-leptonic | W+jets, Drell-Yan | Medium |
| **Photon conversions** | γ→e⁺e⁻ in tracker material | Z+jets | Small (prompt-like IP) |
| **Light hadron misID** | π⁰→γγ, one γ converts; or π± shower fluctuation | QCD multijet | Small |
| **Tau decays** | τ→eνν̄ | W→τν, Z→ττ | Medium |

The dominant source in a W+jets-enriched same-sign CR (used here) is **b/c-hadron decays**.
They produce soft electrons with a large transverse impact parameter (d0), which is why
the electron d0sig < 5.0 cut is important for suppressing them in the SR while still being
loose enough to maintain statistics for fake rate measurement.

### Why do fakes survive the selection?

At low pT (3.5–15 GeV), electron isolation criteria are less efficient at rejecting
semi-leptonic hadron decays because the electron is close to the parent jet. The loose
electron working point (`el_select_loose_NOSYS`) admits more fakes than the tight WP
(`el_select_tight_NOSYS`). The ratio of these two rates is the fake factor / fake rate
that drives the Matrix Method.

---

## 3. Why the Same-Sign CR Works

SM processes producing **real** OS μe pairs (Z→ℓℓ, WW, etc.) are charge-conjugation
symmetric at leading order: they produce μ⁺e⁻ and μ⁻e⁺ with equal probability, but
they cannot produce μ⁺e⁺ or μ⁻e⁻ at LO without one lepton being fake.

Therefore a **same-sign (SS) μe event** is almost certainly:
- A fake electron (non-prompt, the dominant case), OR
- A charge-misidentified real electron (small, subtracted via MC)

After subtracting the prompt MC contribution (charge flip, real prompt SS from WW with
one leptonic/hadronic confusion, etc.), the SS data becomes a **pure fake-enriched sample**
ideal for measuring the fake rate as a function of electron kinematics.

The b-jet veto (`no_bjets`) in the CR removes ttbar events (which have real b-jets and
can produce real OS leptons), keeping the CR enriched in W+jets-like fakes.

---

## 4. The Matrix Method

### Setup

Define two lepton working points:
- **Loose**: `el_select_loose_NOSYS == 1` — lower ID requirement, higher fake acceptance
- **Tight**: `el_select_tight_NOSYS == 1` — higher ID requirement, lower fake acceptance

And two true populations of electrons:
- **Real** (prompt): ε_r = P(tight | real) — real efficiency
- **Fake** (non-prompt): ε_f = P(tight | fake) — fake rate (what we measure)

### The matrix equation

```
N_tight = ε_r · N_real + ε_f · N_fake
N_loose = (1 − ε_r) · N_real + (1 − ε_f) · N_fake
```

Inverting:

```
N_fake = (N_loose · ε_r − N_tight) / (ε_r − ε_f)
N_real = (N_tight − N_loose · ε_f) / (ε_r − ε_f)
```

### Per-event weights (Applied Matrix Method, AMM)

Instead of solving for totals, each event gets a weight:

- **Tight event**: `w_T = ε_f · (ε_r − 1) / (ε_r − ε_f)` ← negative contribution
- **Loose-only event**: `w_L = ε_f · ε_r / (ε_r − ε_f)` ← positive contribution

The sum Σ(w_T + w_L) over all data events in the SS CR gives the total fake prediction.
Weights depend on ε_f(pT, |η|) looked up from the 2D fake rate map and ε_r(pT, |η|) from
the 2D prompt efficiency map (derived from MC).

### Why data minus MC for the fake rate?

The SS CR contains not only fake electrons but also a small fraction of **real prompt**
electrons (from charge-flip, W+γ→Weν, etc.). These are subtracted using the MC prediction
of prompt-only events, isolating the truly fake population:

```
N_fake_tight = N_tight,data − N_tight,MC_prompt
N_fake_loose = N_loose,data − N_loose,MC_prompt
ε_f = N_fake_tight / N_fake_loose
```

---

## 5. Parameterization in pT × |η|

The fake rate is measured as a 2D map in **electron pT × |η|** because:

- **pT dependence**: Fake probability is high at low pT (soft b-decay electrons fail
  isolation less consistently) and falls at high pT (harder electrons are more isolated).
  The bins used: [0, 15, 25, 35, 50, 80, 200] GeV.

- **|η| dependence**: The ATLAS detector geometry changes between barrel (|η| < 1.37),
  transition region (1.37–1.52, excluded), and endcap (|η| > 1.52). The material budget
  increases with |η|, producing more photon conversions in the endcap.
  The bins used: [0.0, 0.7, 1.4, 1.6, 2.0, 2.5].

---

## 6. Systematic Uncertainties

### MET-dependence systematic

The fake composition in the SS CR changes with MET:
- **Low MET (< 30 GeV)**: QCD multijet-enriched — both leptons may be fake, or the W
  is off-shell. Fake electrons tend to be lighter-hadron misidentifications.
- **High MET (> 30 GeV)**: W+jets-enriched — genuine neutrino from W decay produces real
  MET. Fake electrons are more likely from b/c-hadron decays.

The fake rate is measured separately in these two sub-regions. The fractional difference
becomes the MET-dependence systematic:
```
σ_MET = max(|f_MET_low − f_nom|, |f_MET_high − f_nom|) / f_nom
```

### Composition systematic

The fake population in the SS CR (no VBF cuts, no m_jj requirement) may differ from
the fake population in the OS SR (VBF cuts, m_jj > 400 GeV). In the SR, the VBF topology
selects events where the W recoils against forward jets — the b/c-quark fraction and its
pT spectrum may shift relative to the inclusive CR.

This is estimated by comparing the MC prompt efficiency ε_r between the CR and SR:
```
σ_comp = |ε_r(SR) − ε_r(CR)| / ε_r(CR)
```
A change in prompt efficiency between regions signals that the fake composition has shifted.

### Statistical systematic

Per-bin Poisson uncertainty from limited SS CR data statistics, propagated through the
division. In VBF-enriched selections with full Run 2 (139 fb⁻¹), individual pT-|η| bins
can have 10–50% statistical uncertainty.

### Total

```
σ_total = sqrt(σ²_MET + σ²_comp)   [systematic, in quadrature]
σ_stat   — kept separate, uncorrelated
```

---

## 7. MC Closure Test

The MC closure test validates the Matrix Method before applying it to real data.

1. **Truth loop**: Count MC events where `el_truth_isPrompt != 1` — these are the
   actual fake electrons according to generator truth.
2. **MM loop**: Apply the AMM fake weights to the same MC events (without truth filtering).
3. **Compare**: If N_MM / N_truth ≈ 1 (within ~10–20%) in each region, the method
   correctly predicts the fake background yield.

Deviations in the closure indicate either:
- Genuine non-closure of the MM in this topology (→ add closure uncertainty)
- A bug in the IP cut consistency or region definition
- Insufficient statistics to test

---

## 8. Analysis Chain — Execution Order

Run these macros in sequence (see `run_all.sh`):

```
1. MC_eff/fake_bg_estimation.C      → MC_histograms.root
   (prompt electron tight/loose histograms in pT × |η| for each region)

2. data_eff/fake_bg_estimation.C    → Data_histograms.root
   (all electron tight/loose histograms from SS data)

3. fake_eff_calculation/calculate_fake_rates.C   → fake_rates.root
   (fake rate = (data−MC)_tight / (data−MC)_loose per bin)

4. fake_eff_calculation/build_fake_systematics.C → fake_eff_with_systematics.root
   (nominal + MET up/down + comp up/down + stat up/down variations)

5. AMM/AMM_ElectronOnly.C           → data_with_fake_weights.root
   (attaches fake_weight_nominal and systematic variations to each data event)

6. mc_closure_fake.C                → mc_closure_fake_SR_shapes.root
   (validates MM: N_MM vs N_truth in CR/VR/SR)

7. fake_yields_simple.C             → printed yield table
   (final fake + prompt MC yields vs data in each region)

8. studies/VR_CR_fake_sys.C         → validation stacked plots
   (data vs MC+fake with systematic band in VR)
```

---

## 9. Region Definitions

| Region | Charge | m_jj | Purpose |
|--------|--------|------|---------|
| **Preselection** | — | — | Common base: trigger + pT + IP + muon tight |
| **CR** | SS | no cut | Measure fake rate from data |
| **VR** | OS | < 400 GeV | Validate fake prediction (signal-like but low m_jj) |
| **SR** | OS | > 400 GeV | Search region for HNL signal |

All regions except preselection require:
- ≥ 2 jets: jet1_pT ≥ 30 GeV, jet2_pT ≥ 20 GeV
- |Δη_jj| > 3.0
- No b-tagged jets (GN2v01 at 85% WP)

---

## 10. Literature Benchmarks

Compare your results against these published analyses:

| Analysis | Method | Electron fake rate | Total fake syst |
|----------|--------|--------------------|-----------------|
| ATLAS HNL multilepton (arXiv:2204.11138) | Matrix Method | 5–15% | ~20–40% of yield |
| ATLAS SS WW (arXiv:2212.08009) | Matrix Method | 5–20% | ~25–40% |
| CMS HNL (arXiv:1806.10905) | Tight-to-loose ratio | 5–15% | ~20–30% |

**Sanity checks for your fake rates:**
- All 2D bins should be in [0, 1] — if negative, data−MC subtraction went wrong
- Rates above 30% in a high-statistics bin indicate a selection inconsistency
- MC closure ratio N_MM/N_truth in the CR should be within ~15%
- VR data/prediction ratio should sit near 1.0 within the uncertainty band

---

## 11. Known Fixed Bugs (history)

| Bug | File | Fix applied |
|-----|------|-------------|
| b-jet veto inverted (`if no_bjets continue`) | `fake_yields_correct.C:464` | `!no_bjets` |
| pT preselection AND instead of OR | `MC_eff/fake_bg_estimation.C:594`, `data_eff/fake_bg_estimation.C:491` | `&&` → `\|\|` |
| `safe()` returning 1.0 for bad weights | 10 files | `1.0` → `0.0` |
| Electron IP cuts commented out | `MC_eff/fake_bg_estimation.C`, `data_eff/fake_bg_estimation.C` | Uncommented |
| Bitwise `\|` instead of logical `\|\|` in trigger | `fake_comp.C`, `AMM/`, `MC_eff/`, `data_eff/`, `fake_composition/` | `\|` → `\|\|` |
| Division by zero in fake rate | `calculate_fake_rates.C` | Added bin guard with warning |
