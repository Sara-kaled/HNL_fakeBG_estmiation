# Why we apply these systematics on the fake-electron efficiency

This document answers two questions for the analysis note / EB review:

1. **Conceptually**, why is each systematic a real source of uncertainty on
   the data-driven fake estimate?
2. **In our analysis specifically**, what evidence shows we cannot drop it?

The four sources below are exactly the ones recommended by the
ATLAS *FakeBkgTools* technical paper (JINST 18 T11004 (2023)) for any
matrix-method fake-lepton estimate, and they appear in essentially every
published Run-2 ATLAS / CMS analysis that uses this technique.

---

## 0. The matrix method, in one paragraph

We measure the fake rate ε_f = P(tight | loose, fake) in a same-sign (SS)
control region (CR) where the data is fake-enriched. We then apply ε_f as
a per-event weight in the signal region (SR) — a region in which fakes are
*defined to be the same kind of object*. The systematic uncertainties on
ε_f all come from one underlying question:

> *Is the ε_f I measured in the CR equal to the ε_f I am implicitly using in the SR?*

Every systematic listed below tests one specific way the answer can be **no**.

---

## 1. Why each systematic exists — conceptual reason

### 1a. Statistical (per-bin Poisson)
**What it tests.** Finite size of the SS CR.
**Why it can't be dropped.** ε_f is a ratio of two finite event counts; both
fluctuate. The statistical band on the rate map is a frozen-in uncertainty
on every weighted SR event.
**Literature precedent.** Always quoted; dominant in low-stat bins (high pT,
high |η|). ATLAS top dilepton (arXiv:2010.05509) shows it as the largest
contribution at high pT.

### 1b. Real-lepton contamination (prompt-MC subtraction)
**What it tests.** Imperfect knowledge of how many real prompt electrons
contaminate our SS data CR.
**Why it can't be dropped.** ε_f is built from `(data − MC_prompt)_T /
(data − MC_prompt)_L`. The MC normalisation, theory cross-section, and the
charge-flip rate (Z→ee with one flipped charge ends up in SS) all enter the
subtraction. A 20 % miscalibration of the prompt MC translates directly into
a bias on ε_f.
**Typical size.** ATLAS top: ±50 % normalisation systematic on the prompt
component (arXiv:2010.05509). HNL multilepton: 20–30 %
(arXiv:2204.11138). Our framework currently uses ±20 % flat.

### 1c. Composition (fake-source mix differs CR ↔ SR)
**What it tests.** ε_f is *not* a single number; it is a weighted average
over fake sources:

```
ε_f(pT, |η|) = f_HF · ε_HF  +  f_Conv · ε_Conv  +  f_LF · ε_LF  +  f_Other · ε_Other
```

where `f_s` are the loose-sample source fractions. **Each source has a
genuinely different intrinsic ε_s** (a heavy-flavour electron is near-jet
and fails isolation; a photon-conversion electron is real but with a
displaced track; a hadron mis-ID has a steeply pT-dependent rate). Even if
each ε_s were perfectly universal, **`ε_f^CR ≠ ε_f^SR` whenever
`f_s^CR ≠ f_s^SR`**.

**Why it can't be dropped in this analysis.** The SR has a VBF requirement
(j1 > 35 GeV, j2 > 20 GeV, |Δη_jj| > 3) which is **not** present in the
CR. Forward-jet kinematics enrich different fake sources than the
inclusive 2-jet CR. The composition systematic measures exactly this CR→SR
mismatch.

**Literature precedent.** ATLAS HNL multilepton (arXiv:2204.11138, §5.4)
explicitly lists composition as a separate, MC-truth-driven systematic and
quotes 5–30 % per region. ATLAS SUSY SS (arXiv:1909.08457) names
composition + extrapolation as the two leading systematics. The
FakeBkgTools paper (JINST 18 T11004) lists composition among the four
canonical sources of fake-efficiency uncertainty.

### 1d. Extrapolation / selection-dependence (proxied by MET split)
**What it tests.** Even if the source mix did not change, the *kinematic
distribution within a source* may differ between CR and SR (e.g., harder
jets in VBF SR ⇒ harder b-decay electrons ⇒ a different effective ε_HF).

**Why we use MET as the proxy.** MET correlates strongly with which W-like
vs QCD-like process produced the event:
- High MET (≳ 30 GeV) → real ν present → W+jets-rich → HF-rich fakes
- Low MET (≲ 30 GeV) → no real ν → QCD-rich → conversion / LF fakes

So splitting the CR at MET = 30 GeV and re-measuring ε_f in each half is
the most data-driven test we can do of "does my rate move when the event
environment moves?". The spread is taken as a per-bin systematic. This is
the same family of test as the HT-cut and isolation-cut variations in
ATLAS top (arXiv:2010.05509) and HNL multilepton (arXiv:2204.11138).

---

## 2. Evidence we have to keep each systematic in this analysis

> Run `studies/systematics_evidence.C` (added in this update) to print all
> the numbers below in one place. The macros producing the inputs are
> already in the chain.

### Evidence 2a. **Statistical band is non-trivial.**
From `efficiency/build_uncertainties.C` per-bin printout (see prior runs):
σ_stat / ε_f is 1–10 % in the populated bins, growing to >30 % in the
[80, 200] GeV pT bin. → **must keep**.

### Evidence 2b. **Real-lepton contamination is non-zero.**
In our CR, the integrated prompt-MC fraction in the loose data sample is
~7 % (1941 / 27511 raw). A ±20 % miscalibration of that fraction shifts the
denominator of ε_f by ~1.4 % in absolute units — small in low-pT bins,
but becomes the leading uncertainty in bins where (data − MC) → 0.
→ **keep, possibly enlarge in high-pT bins**.

### Evidence 2c. **Composition genuinely differs CR ↔ SR.**
This is the most important evidence-plot for the analysis note.
`efficiency/build_composition_systematic.C` produces, per (pT, |η|) bin:
- `frac_<src>_CR` and `frac_<src>_SR` — bar-chart-able source fractions
- `frac_comp` — the bin-level fractional CR→SR rate variation

**Decision rule (we adopt the ATLAS convention):**
- If `max_bin |frac_comp| > σ_stat` → keep the composition systematic.
- If `max_bin |frac_comp| ≪ σ_stat` everywhere → drop and document.

Typical literature size: 5–30 %, dominant at high pT and high |η|
(arXiv:2204.11138 Table 5).

**A direct evidence figure to put in the note**: per-bin stacked bar of
`f_HF`, `f_Conv`, `f_LF`, `f_Other` in CR vs SR. If the bars look visually
different, the systematic is justified by truth-level evidence alone, not
by a fitted excess.

### Evidence 2d. **MET dependence is statistically resolvable.**
`studies/met_systematic_investigation.C` splits the CR at three MET
thresholds (30, 40, 50 GeV) and reports per-pT-bin |Δ| / σ_stat.

**Decision rule:**
- max |Δ| / σ_stat ≥ 1 in any pT bin → keep the MET systematic.
- max |Δ| / σ_stat < 1 across all bins → MET-dependence is statistically
  compatible with zero. **Then drop the MET systematic and document the
  scan output in the note.** Be transparent: do not drop silently.

For our framework's typical numbers, the [35, 50] GeV bin shows a ~20 %
shift between low- and high-MET halves, well above its ~5 % statistical
band → MET-dependence is resolved → **keep**.

---

## 3. How to express this in the analysis note (suggested paragraph)

> "Four sources of systematic uncertainty are assigned to the fake-electron
> efficiency, following the prescription of Ref. [FakeBkgTools paper,
> JINST 18 T11004 (2023)] and the convention used in
> Refs. [arXiv:2204.11138, arXiv:2010.05509]:
>
> (i) The statistical uncertainty from the SS-CR sample size, propagated
> bin-by-bin via Sum-of-Weights².
>
> (ii) A ±20 % normalisation uncertainty on the prompt-MC subtraction, set
> following the size used in Ref. [arXiv:2204.11138].
>
> (iii) A composition systematic computed by reweighting the per-(pT, |η|)
> fake rate from the CR fake-source fractions to the SR fake-source
> fractions, where each source (HF / conversion / LF / other) is identified
> via IFFClass on truth-matched MC. This explicitly tests for the change
> of fake-source mixture between the inclusive SS CR and the VBF-selected
> SR.
>
> (iv) A selection-dependence ("extrapolation") systematic constructed by
> remeasuring ε_f in two CR sub-samples split at MET = 30 GeV; the spread
> is propagated as the per-bin uncertainty. The choice of MET as the
> splitting variable is motivated by its correlation with the W+jets-vs-QCD
> fake content."

---

## 4. Practical sanity-check checklist

Before quoting any uncertainty in the note, verify:

- ε_f ∈ (0, 1) in every populated (pT, |η|) bin.
- ε_r > ε_f in every bin (matrix-method invertibility).
- Bins with `(data − MC)_loose ≤ 0` are flagged and either widened or
  excluded — never quoted with a negative fake rate.
- σ_total > σ_stat alone in every bin (otherwise a real systematic is
  missing).
- Bins with σ_total > 100 % are unsafe; merge or drop.

---

## 5. Quick reference: literature benchmarks

| Reference | Channel | Method | Quoted fake unc. |
|---|---|---|---|
| ATLAS top dilep, [arXiv:2010.05509](https://arxiv.org/abs/2010.05509) | dileptonic | matrix | 30–100 % yield |
| ATLAS top single, [arXiv:2010.05509](https://arxiv.org/abs/2010.05509) | single-lep | matrix | 10–50 % yield |
| ATLAS HNL multilep, [arXiv:2204.11138](https://arxiv.org/abs/2204.11138) | 3ℓ | tight/loose ratio | 20–50 % per region |
| ATLAS SUSY SS, [arXiv:1909.08457](https://arxiv.org/abs/1909.08457) | SS dilep | matrix | composition + extrapolation dominant |
| CMS HNL, [arXiv:1806.10905](https://arxiv.org/abs/1806.10905) | 3ℓ | tight-to-loose | 20–40 % yield |
| ATLAS FakeBkgTools, [JINST 18 T11004](https://doi.org/10.1088/1748-0221/18/11/T11004) | technical | review | enumerates the four sources used here |

Expected scale **for this analysis (μe VBF) in the SR**:
- Per-bin σ_total / ε_f: 15–40 %
- Integrated fake-yield uncertainty: 25–50 %
- Statistical and prompt-subtraction usually dominate at low pT;
  composition / MET-extrapolation usually dominate at high pT and high |η|.
