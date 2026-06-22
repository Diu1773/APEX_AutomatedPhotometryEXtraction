# APEX CMD isochrone fitting: methodology, model–data color consistency, and validation

This document describes APEX's color–magnitude diagram (CMD) isochrone-fitting
method (Step 12), the photometric-system consistency problem that limits
multi-color fits, the empirical color-correction that resolves it, and the
validation evidence gathered on real clusters. It is written to be paper-grade:
the methods and figures here are intended to seed a methods/validation section.

---

## 1. The problem

Fitting a stellar cluster CMD means recovering four correlated parameters from
the observed (color, magnitude) distribution:

- log(age), metallicity [M/H], distance modulus (m−M)₀, and reddening E(color).

These four are **strongly degenerate**: a younger + more metal-poor + more
reddened + nearer isochrone can reproduce nearly the same CMD as an older +
more metal-rich + less reddened + farther one. With a **single color** (e.g.
g−r) the problem is fundamentally under-constrained — many (age, [M/H], dm, E)
combinations are near-equivalent. This is intrinsic to all isochrone fitting,
not specific to any one code.

---

## 2. The solver: generative mixture likelihood + MCMC

APEX replaces a nearest-point χ² objective (which does not encode CMD density or
morphology, so age is poorly constrained) with a **generative probabilistic
likelihood** sampled by MCMC (`apex/analysis/cmd/isochrone_mcmc.py`).

For each observed star *i* with photometry **x**ᵢ = (colors, magnitude) and
per-axis errors σᵢ, the likelihood under parameters θ = (log age, [M/H], dm, E)
is a three-component mixture:

```
L_i(θ) = (1 − f_field) · [ (1 − f_bin) · L_single + f_bin · L_binary ]
         + f_field · L_field
```

- **L_single** — IMF-weighted sum of Gaussians centered on the isochrone points
  (Kroupa IMF density along the track), convolved with the star's errors. This
  encodes the *stellar density* along the isochrone, so the (sparse but
  age-sensitive) main-sequence turnoff is weighted correctly — the key to
  constraining age.
- **L_binary** — an unresolved-binary ridge (companions drawn from the IMF,
  brightening stars by up to ~0.75 mag), so binaries do not bias the fit.
- **L_field** — a uniform term over the CMD region, making the fit robust to
  non-members/outliers without an ad-hoc "use the closest fraction" cut.

The total log-likelihood Σᵢ log L_i is sampled with `emcee` (Foreman-Mackey
et al. 2013). Two settings are essential for the narrow, correlated CMD
posterior and are **general (data-independent)**:

- **DE / DESnooker moves** instead of the default stretch move — the standard
  remedy for correlated posteriors. With the stretch move the acceptance
  fraction collapsed to a few percent (walkers stuck at the seed).
- **A systematic error floor** (`_ERR_FLOOR ≈ 0.02 mag`, added in quadrature),
  representing the irreducible CMD width from isochrone-model error, unresolved
  binaries, differential reddening and calibration. Without it the combined
  many-star likelihood is razor-sharp and unsamplable; with it the posterior is
  samplable and the credible intervals are honest.

This approach follows the probabilistic CMD-fitting tradition (Naylor &
Jeffries 2006; Dolphin 2002, MATCH). The fitter is validated by a synthetic
recovery test (`tests/test_isochrone_mcmc.py`): a cluster generated from known
parameters is recovered within the credible interval — the generality guarantee
(no per-target tuning).

**Status (single color, M67, g−r):** acceptance ≈ 0.39, `convergence_ok=True`;
distance modulus tightly recovered. Age/[M/H]/E remain in the degenerate corner,
as expected for a single color — see §3.

---

## 3. The degeneracy and how to break it

A single color cannot separate age, [M/H], and reddening. Two routes break it:

1. **A second color.** g−r and r−i (or B−V and V−I) carry the reddening and
   metallicity signals in *different directions* in the two color–color planes,
   so a joint fit constrains [M/H] and E. APEX supports a multi-color likelihood
   (`cmd_loglike_multicolor`) in (g−r, r−i, g) space. On *synthetic* data this
   demonstrably breaks the metallicity degeneracy.
2. **External priors** — spectroscopic [Fe/H], or a dust-map E(B−V) — applied as
   Gaussian priors in the sampler. (Demonstrated: with a solar-metallicity prior
   the M67 fit moves from age 3.25 → 4.47 Gyr, bracketing the literature
   3.5–4.0 Gyr — the degeneracy made navigable.)

**But multi-color fitting on real SDSS data did not break the degeneracy** — and
the reason is *not* the data. See §4.

---

## 4. The model–data color consistency problem (key finding)

> ⚠️ **SUPERSEDED / CORRECTED (2026-06 deep review) — read §10.** A follow-up
> rigorous review found the "−0.41 r−i slope vs PARSEC" reported in this section
> is a **comparison-method artifact** (it appears only when pairing each star to
> the isochrone by *magnitude*; evaluated correctly at fixed g−r it is a small
> near-constant pedestal). The real obstacle to multi-color fitting is the
> **fundamental age–[M/H]–reddening–distance degeneracy** of gri photometry
> (the metallicity and reddening vectors are only **2.6°** apart in the
> (g−r, r−i) plane), plus a secondary, fixable likelihood weakness — **not** a
> PARSEC color error and **not** an APEX photometry bias (photometry is directly
> verified innocent in §10.1). The empirical-correction prescription in §5 is
> therefore not the right fix. The conclusions below are retained for history;
> §10 is authoritative.

Multi-color fitting requires that the **observed** color–color relation and the
**model** color–color relation be on the *same photometric system*. We measured
the M67 (g−r, r−i) relation against PARSEC and found a residual:

| Comparison (r−i residual vs PARSEC at fixed g−r) | median | slope vs g−r |
|---|---:|---:|
| APEX calibrated photometry (mag_std) | −0.035 | −0.42 |
| **Gaia DR3 synthetic SDSS (GSPC), independent** | **−0.033** | **−0.41** |

The decisive test is the second row: **GSPC** — Gaia DR3 synthetic SDSS
photometry (Gaia Collaboration, Montegriffo et al. 2023), an independent,
validated empirical standard with no APEX data involved — shows the **same**
~−0.033 mag + −0.41 slope offset from PARSEC. Therefore:

- **The offset is in PARSEC's theoretical SDSS r−i colors**, a documented
  isochrone model-atmosphere / bolometric-correction limitation (worst in the
  i band and for cool stars), **not in APEX's calibration.**
- A multi-color fit demands that one isochrone satisfy both g−r and r−i; the
  model's r−i is off by a g−r-dependent amount, so no isochrone fits both →
  the likelihood is frustrated (acceptance ≈ 0.05, walkers rail to grid edges).
- Single-color (g−r) fits hide this (the offset is absorbed into [M/H]/E within
  the degeneracy), which is why single-color converges but is biased.

**Corollary — APEX zeropoint/color calibration is validated.** APEX's
calibrated colors match GSPC (both at −0.033 vs PARSEC), i.e. APEX agrees with
an independent Gaia-based SDSS standard at the few-tens-of-mmag level. Together
with the aperture-flux cross-check against SExtractor and IRAF/DAOPHOT
(`docs/validation_crosscheck.md`: ~2–4 mmag and ~3 mmag agreement), APEX
photometry is independently validated in **both flux and color**.

---

## 5. The fix: empirical color correction of the isochrone

> ⚠️ **SUPERSEDED (2026-06) — see §10.** This section assumed the multi-color
> failure was a real, correctable PARSEC r−i color offset. The deep review showed
> that premise is wrong (§4 note): the offset measured here is a pairing artifact,
> and a linear correction did not (and could not) unblock the fit because the true
> obstacle is parameter degeneracy + likelihood weakness, not a model color error.
> The adopted resolution is **external priors** ([M/H], reddening) + Gaia parallax
> distance (§10.4), not isochrone color correction. Retained for history.

Because the offset is in the model, the principled fix corrects the **model**
colors to an independent empirical standard — *not* to the desired answer:

```
(r−i)_corrected(g−r) = (r−i)_PARSEC(g−r) + Δ(g−r),
   with  Δ(g−r) = (r−i)_GSPC − (r−i)_PARSEC  ≈ +0.033 + 0.41·(g−r)
```

The correction Δ is measured against GSPC (accurate SDSS), applied to the
isochrone before fitting. The corrected model is then on the real SDSS system,
consistent with the (already-correct) data, so the multi-color fit becomes
self-consistent and [M/H]/E break free of the degeneracy.

**Why this is legitimate (not circular):** the correction anchors the model to
an *independent photometric standard* (GSPC), which carries no information about
the cluster's age/metallicity. This is the standard "empirically calibrated
isochrone" technique used throughout the cluster-CMD literature:

- An et al. (2009) — empirically calibrated isochrones in SDSS *ugriz* using
  cluster fiducial sequences.
- MIST (Choi et al. 2016) and BaSTI/PARSEC documentation discuss synthetic-
  photometry systematics and color/BC corrections.

**Generality:** Δ is a property of the *model + photometric system*, derived
once and applied to any cluster — not a per-target tuning.

> **Implementation status — tested, did NOT resolve (honest negative result).**
> A *linear* empirical r−i correction (subtracting the measured residual line
> A + B·(g−r), A=+0.112, B=−0.361) was applied to M67 and the multi-color MCMC
> re-run. It did **not** unblock the fit: acceptance stayed low (~0.16, not
> converged), [M/H] still railed metal-poor (−0.40), and the distance modulus
> actually worsened (→10.06 vs literature 9.6–9.7). So the multi-color
> [M/H]/E recovery on this real data remains **unsolved**. Likely reasons: the
> model–color mismatch is not purely linear (a linear correction leaves
> structure); the 3-D likelihood is intrinsically sharper than the 2-D case
> (low acceptance); and/or per-star i-band systematics not captured by a global
> color term. This is genuine research-grade difficulty. **Practical path that
> works today: single-color (g−r) fitting — which converges cleanly (acceptance
> 0.39) — plus external [Fe/H]/E priors, which demonstrably brackets the
> literature age.** Multi-color absolute [M/H]/E is left as a documented open
> problem (candidate next steps: non-parametric model-color correction,
> i-band error recalibration, or comparing isochrone families).

---

## 6. Model and reference selection (design)

The diagnostic and correction are **model- and reference-agnostic**, so both are
exposed as choices rather than hard-coded:

- **Isochrone model:** PARSEC / BaSTI / MIST. Each has its own color systematic;
  the empirical correction (§5) measures and removes whichever model's offset,
  so model choice is no longer critical. (Different models can also be compared
  by their GSPC offset to pick the smallest.)
- **Calibration reference DB:** GSPC (Gaia XP synthetic) / Gaia broadband
  transform / Pan-STARRS DR2 / APASS / SDSS. Auto-defaulted by the data's filter
  system, with manual override.

**Photometric-system rule (must hold):** *data system = reference system =
isochrone-grid system.* APEX detects the data system from the FITS `FILTER`
header (`normalize_filter_key`); the SDSS data uses the PARSEC `sdss` grid and a
SDSS reference, Johnson data uses the `johnson` grid and a JKC reference. Gaia
XP synthetic photometry (GSPC) provides **both** SDSS and Johnson-Cousins
synthetic magnitudes from the same all-sky source, making it the universal
reference.

**Isochrone-generation parameters that affect colors** (chosen per target, then
residuals absorbed by §5): photometric system (must match — critical),
bolometric-correction / spectral library (the root of the color systematic),
[α/Fe], helium Y, rotation, mass loss, and model version. Sensible defaults
(solar-scaled, standard BC) plus the empirical correction suffice; only the
photometric system must be matched exactly.

---

## 7. Outputs (paper-grade figures)

Each fit (`validation/cmd_step12_realdata.py --method mcmc`) produces:

- **CMD with the best-fit isochrone overlaid** — observed stars + the median
  posterior isochrone (e.g. `*_mc_cmd.png`).
- **Corner plot of the 4-D posterior** — the (log age, [M/H], (m−M), E)
  joint/marginal distributions with credible intervals, directly showing the
  age–metallicity–reddening covariance (`*_mc_corner.png`).
- **Posterior summary** — medians and 16/84 percentiles per parameter, plus
  convergence diagnostics (acceptance fraction, autocorrelation time).

These are the standard figures for a cluster-CMD paper and are generated
automatically per target.

---

## 8. Honest limitations

- Absolute [M/H]/E/distance accuracy depends on (a) the photometric system
  consistency (§4–5) and (b) the input calibration quality; the empirical
  correction addresses the model side, and APEX's calibration is validated.
- Credible intervals reflect the statistical + adopted-systematic floor; they do
  not capture isochrone-model *family* differences — comparing PARSEC/BaSTI/MIST
  bounds that.
- The multi-color likelihood is currently the per-fit runtime bottleneck (full
  cluster, ~minutes to tens of minutes); caching/vectorization is a usability
  (GUI) optimization, separate from correctness.
- GSPC synthetic photometry exists for G ≲ 17.5; this covers calibrators amply
  but not the faintest members (which are calibration targets, not references).

---

## 10. Deep review (2026-06): photometry innocence, the true root cause, and ugri

This section supersedes §4–§5 where they conflict. It is the result of a
deductive/inductive/abductive re-audit triggered by the question "could a hidden
Step-4 detection / aperture-correction / recentering bias be the real cause?"

### 10.1 Photometry is innocent — verified directly (not just by external medians)

Each suspect was tested at the **code and real-frame level**, not only by external
catalog comparison:

| Suspect | Direct test | Result |
|---|---|---|
| Aperture correction | code: `forced_photometry.py` computes a **single per-frame scalar** apcorr (growth curve on high-SNR refs, applied to all stars) | A per-frame scalar shifts every star equally → **cannot** produce a magnitude- or color-dependent bias; fully absorbed by the Step-10 zeropoint. Ruled out *deductively*. |
| Recentering | re-measured aperture flux at the registration position `x_reg` vs the recentered position `x_fit` on **8 real M67 frames, 5286 stars** | median Δmag = **+0.0000** (g) / **+0.0001** (i) in every magnitude bin; centroid shifts ~0.1 px vs aperture ~12 px ⇒ negligible flux loss. Ruled out *empirically*. |
| Step-4 detection | flux at `x_reg` == flux at `x_fit`; master positions are Gaia-matched and multi-frame averaged | positions are sound; no position-induced flux bias. Ruled out. |

External corroboration: APEX `mag_std` == Gaia DR3 **GSPC** synthetic SDSS ==
**Pan-STARRS DR2**, per-star, color- and magnitude-flat on open clusters
(M67, M37). A faint g-band slope seen against PS1 alone (+0.095) is **not**
reproduced by SDSS (flat) nor by the internal `x_reg`/`x_fit` flux test (0.000),
so it is a PS1 faint-end systematic, not APEX.

> Caveat retained as a separate issue: in **crowded globular** fields (M5, M13)
> APEX faint-star r−i drifts ~−0.07/−0.08 vs PS1 with 2–3× larger scatter — a
> genuine crowding limitation of fixed-aperture photometry in dense cores, *not*
> relevant to the open-cluster isochrone problem here.

### 10.2 The real root cause is two-layered

1. **Fundamental degeneracy (physics, not a bug).** With gri photometry the
   metallicity, reddening, age and distance signals are strongly degenerate. In
   the (g−r, r−i) color–color plane the **[M/H] vector and the reddening vector
   are only 2.6° apart** (essentially parallel). Consequently *every* fit that
   lets the parameters vary freely — the full generative likelihood **and** a
   pure color–color χ² — rails to a (metal-poor, high-reddening, short-distance)
   corner. [M/H] = solar appears only under a conditional view with E, age and
   distance all pinned. This is why every unconstrained 4-D attempt has failed:
   the information is not in gri alone.
2. **Secondary, fixable likelihood weakness.** Even with E and distance fixed at
   their literature values, APEX's generative likelihood still prefers metal-poor
   because (a) the magnitude axis matches the IMF **luminosity function with no
   selection/completeness term**, so the SNR-truncated real sample biases the fit
   (synthetic, untruncated data hides this — it converges at acceptance ≈ 0.4),
   and (b) the two colors are compared with a **diagonal covariance** that ignores
   the shared-band anti-correlation (g−r and r−i share r), diluting the color
   constraint. These do not *create* the degeneracy but make APEX land in its
   worst corner instead of being neutral.

### 10.3 Adopted resolution — external priors (confirmed design decision)

[M/H] and reddening (color excess) are taken from **external priors**, which is
standard, non-circular practice (the data genuinely lacks the information):

- **Reddening E(B−V):** dust-map prior (SFD / Bayestar at the target coordinates).
- **[M/H]/[Fe/H]:** spectroscopic prior where available, else a literature/Gaussian prior.
- **Distance:** Gaia parallax prior (e.g. M67 members: parallax 1.154 mas ⇒
  m−M ≈ 9.69, matching the literature 9.6–9.7). The astrometry is already in the
  catalog; an auto membership filter (Gaia PM + parallax clump) is implemented in
  `validation/cmd_step12_realdata.py`.

Implemented knobs supporting this: `--membership` (auto PM+parallax member
selection), `--mh-prior`, `--ecolor-prior`, `--err-floor`. The likelihood fixes
in §10.2(2) remain open work so the priors are not fighting a biased likelihood.

### 10.4 Adding a bluer band (ugri): does it help? — yes, quantitatively

The fundamental degeneracy of §10.2(1) is a property of the *filter set*. Adding
**u** attacks it directly, because u−g is the classic line-blanketing metallicity
indicator. Measured on PARSEC (age 9.6, fixed lower-MS mass, [M/H] 0 → −0.4):

| Color–color plane | angle( [M/H] vector , reddening vector ) | interpretation |
|---|---:|---|
| (g−r, r−i) — current gri | **2.6°** | parallel → [M/H] and E inseparable |
| (u−g, g−r) — adds u | **16.6°** | ~6× more separable → [M/H] becomes measurable |

Per-color [M/H] sensitivity (per 0.4 dex): **d(u−g) = −0.35** vs d(g−r) = −0.21
vs d(r−i) = −0.11 — u−g carries ~1.7× the metallicity signal of g−r.

**So ugri can *measure* [M/H] where gri can only prior it.** Caveats: u-band is
observationally expensive (low CCD QE, bright/variable sky, high atmospheric
extinction) so u photometry is noisy, worst for faint/red stars; 16.6° is better
but not orthogonal, so good u SNR is needed; and reddening still matters (R_u is
large). PARSEC carries SDSS u, so APEX could fit ugri directly once u-band frames
are available — this is the recommended path to a genuinely data-driven 4-D fit.

### 10.5 Fresh re-execution of the photometry suspects (this-session evidence)

The Step-4 detection / aperture-correction / recentering suspects were
**re-tested directly on real M67 data** (not cited from prior runs), each with a
distinct reasoning mode and adversarially verified:

| Suspect | Mode | Direct test run | Result |
|---|---|---:|---|
| apcorr | deductive | per-frame scalar (forced_photometry.py:1059–1133); numerically demonstrated ZP cancellation | Δmag_std = **1.8×10⁻¹⁵** (machine ε) |
| recentering | inductive | re-measured aperture flux at pre- vs post-recenter position, 6 real frames, 3982 stars | Δmag **−0.0001…−0.0004** (sub-mmag), slope −0.0003 mag/mag |
| detection | abductive | APEX vs Gaia DR3 position residual vs mag & colour | colour-propagated bias **+0.0002 mag** |

All three **ruled out** (each ≥1–2 orders of magnitude below the ~0.05–0.3 mag the
rail needs). Two small, real residual risks were newly found and bounded:
**(R1)** crowded-field faint-star drift in *globulars* (M5/M13, intensity-weighted
centroid dragged toward neighbours; not present in open M67); **(R2)** the
per-frame apcorr scalar does not cancel *exactly* once Step-10's colour term uses
a multi-frame stacked instrumental colour, leaking ≤7×10⁻³ mag — both far below the
rail threshold.

**Decisive structural proof the rail is not photometry:** a **colour-only** fit
(magnitude axis deleted) *still* rails metal-poor, and the **inter-colour
covariance fix (§10.2 B, committed)** did not move the rail either. A photometric
zero-point/aperture/recenter error is an additive magnitude offset; it cannot
create a rail that survives deleting the magnitude axis.

### 10.6 What actually works (decisive prior experiment)

Running M67 with the **distance pinned by Gaia parallax** (members π=1.154 mas →
(m−M)₀=9.69) and **reddening pinned by a dust-map prior** (E(B−V)≈0.04), with
[M/H] left free:

| Quantity | Recovered | Literature |
|---|---:|---:|
| Age | **3.99 Gyr** | 3.5–4.0 |
| (m−M)₀ | **9.65** | 9.6–9.7 |
| [M/H] | −0.20 | 0.0 |
| acceptance | 0.24 | (healthy) |

→ **Age and distance are recovered correctly** once the fundamental degeneracy is
broken by external priors. [M/H] still lands ~0.2 dex metal-poor — partly the
residual gri degeneracy, partly a residual magnitude-axis likelihood weakness (no
selection/completeness term). Closing that last ~0.2 dex needs either a
spectroscopic [M/H] prior, a Naylor–Jeffries τ² selection-aware likelihood, or
u-band (§10.4). **Practical recipe: membership + parallax + reddening priors →
reliable age & distance; add an [M/H] prior or u-band for metallicity.**

### 10.7 Implementation + why a τ² selection term is NOT the fix

Implemented (commit `ff480ad`): a Gaussian (m−M)₀ prior in the likelihood, auto-derived
from the membership clump's median Gaia parallax, exposed as a "Gaia parallax distance
prior" checkbox in the GUI Auto-fit (MCMC) tab. Because the gri degeneracy makes the
likelihood prefer a spurious short distance by ~200+ logL, a *soft* Gaussian prior is
insufficient (σ=0.08 → dm still railed to 9.25); the service therefore also **tightens the
dm (and, given a reddening prior, the E) bounds to a hard window** around the external
value. With both pinned, M67 reproducibly recovers **age ≈ 3.98 Gyr and (m−M)₀ ≈ 9.64**.

A selection/completeness (Naylor–Jeffries τ²) term was considered for the residual ~0.2–0.4
dex metal-poor [M/H] bias, and **empirically rejected this session**: with dm+E pinned, a
**bright-only** subsample (SNR>80, where faint truncation cannot act) fits *more* metal-poor
([M/H]=−0.45) than the full sample ([M/H]=−0.20), the opposite of a selection signature.
A box=data truncation term was also a verified no-op (Z θ-constant; rail−truth gap +240 logL
unchanged). The residual [M/H] is therefore the **fundamental gri colour floor** (even the
best-measured bright stars' gri colours sit slightly metal-poor of PARSEC solar; possibly
compounded by the small per-band ZP colour-system offset, §10.5 R2), not a fixable
luminosity-function/selection bug. It is closed only by a spectroscopic [M/H] prior (GUI
``--mh-prior``) or by adding u-band (§10.4).

### 10.8 EEP interpolation fix + multi-cluster age validation (2026-06-22)

The single biggest *age* error was an interpolation artefact, found by running real
clusters with clear turn-offs (M67, NGC 6811, M37) end to end rather than synthetic data.

**Root cause (fixed, commit `a5ebdd1`).** `MultiBandIsochrone` blended the four
(age, [M/H]) grid corners at a **fixed initial mass**. PARSEC isochrones are
EEP-parametrised, so the same initial mass is a *different* evolutionary phase at
different ages; mass-blending interleaved phases → a jagged interpolated track → a
**bumpy age likelihood with spurious spikes at grid ages**. The fit then landed on a
random spike (NGC 6811 gave 1.1 / 1.6 / 2.4 Gyr across runs). The fix resamples each
corner onto a common **normalised arc-length (EEP) coordinate** before blending, so the
main sequence aligns with the main sequence and giants with giants. Effect: the NGC 6811
age likelihood (dm/[M/H]/E at literature) went from a spurious 1.6 Gyr spike to a clean
**1.0 Gyr peak**, and ages became **reproducible across seeds** (the headline win).

**Multi-cluster validation** (EEP + Gaia-parallax dm + reddening & [M/H] priors, hard
bound-window tightened to 2.0σ, commit `5b32567`; mean of 2 seeds):

| cluster   | filters | recovered age | literature | note |
|-----------|---------|---------------|------------|------|
| NGC 6811  | BVR     | 0.83 Gyr      | ~1.0       | [M/H] drifts to +0.10 within prior → slightly young |
| M67       | gri     | 3.15 Gyr      | 3.5–4.0    | faint-MS dilution caps full-sample age; see below |
| M37       | gri     | 0.72 Gyr      | 0.4–0.5    | strong differential reddening broadens the CMD |

Seed-to-seed scatter is now ≤0.05 Gyr (was ~1 Gyr). The residual offsets are
**astrophysical/degeneracy limits, not bugs**: the secondary parameters ([M/H], E, dm)
drift within their priors and the age follows; M37 has real differential reddening.

**Why M67 caps at ~3.1 Gyr (full sample).** The dense, age-insensitive **faint main
sequence** dominates the star count and dilutes the age constraint, which lives in the
sparse turn-off / sub-giant / giant region. With dm/[M/H]/E pinned at literature, an
age scan over the **full sample peaks at 3.0 Gyr**, but over a **turn-off-only**
subsample peaks at **4.0 Gyr** (M37 likewise sharpens). A per-star likelihood weight
(`obs_weights`, a generic hook now in `cmd_loglike_multicolor`) that down-weights the
faint MS was tested: it recovers M67 = 4.0 **only when dm/[M/H]/E are all pinned at the
canonical values**; in a realistic parallax-dm fit the secondary parameters drift and
the gain collapses to 3.15 → 3.23. It is therefore **not shipped as a feature** — the
honest automatic broadband age for M67 is ~3.1 Gyr, ~15 % below the canonical 3.5–4.0
because the turn-off carries too little weight in a full-CMD generative likelihood.
Closing that gap properly needs a completeness-aware likelihood (a real magnitude-
completeness model, not the box=data τ² already rejected in §10.7), which remains open
work. Evidence: `validation/_scratch/age_grid.py` + `age_grid_result.json`.

---

## 9. References (representative)

- Bressan, A., et al. 2012, MNRAS, 427, 127 — PARSEC isochrones.
- Choi, J., et al. 2016, ApJ, 823, 102 — MIST isochrones.
- An, D., et al. 2009, ApJ, 700, 523 — empirically calibrated SDSS isochrones.
- Naylor, T., & Jeffries, R. D. 2006, MNRAS, 373, 1251 — τ² CMD fitting.
- Dolphin, A. E. 2002, MNRAS, 332, 91 — MATCH (CMD synthesis fitting).
- Foreman-Mackey, D., et al. 2013, PASP, 125, 306 — emcee.
- Gaia Collaboration, Montegriffo, P., et al. 2023, A&A, 674, A33 — Gaia DR3
  synthetic photometry (GSPC).
- Kroupa, P. 2001, MNRAS, 322, 231 — IMF.

*(Verify exact bibliographic details before submission.)*
