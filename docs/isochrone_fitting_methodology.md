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

> Implementation status: the diagnostic (model vs GSPC) is built and quantified;
> the correction module + multi-color re-validation on M67 (expected: [M/H]→~0,
> E(B−V)→~0.04, age→3.5–4 Gyr) is the next step.

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
