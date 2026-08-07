# Figure 6 — Frame-QC validation: injected frame defects vs automatic decisions

**Figure 6.** Validation of the automatic frame-quality-control module
(`apex/analysis/frame_qc.py`, `evaluate_frame_qc`) on a synthetic "night" of 44
frames with defects injected by construction. Each frame is an independent
realization of the self-contained synthetic reference field (170 Gaussian-PSF
stars, uniform true magnitudes in $[13.5, 19.0]$, $640^2$ px, Poisson photon +
Gaussian read noise in the electron domain, gain $g = 1.5$
e$^-$/ADU, zeropoint 25.0, FWHM 3.5 px, sky 150 ADU, RN 5 e$^-$; star
separation fixed at 21 px for *all* classes so every frame carries the same
170-star geometric truth). Twenty defective frames alter exactly one physical
property: **cloud** (zeropoint $-0.7$ mag; grey transparency loss, sky
unchanged), **bad seeing** (FWHM $\times 1.8 \to 6.3$ px), **bright sky**
(background $\times 5 \to 750$ ADU), and **noisy readout** (realized with RN
$= 25$ e$^-$ while the header — and therefore the QC
input — still claims 5 e$^-$: an undocumented-electronics
"header lie" that only a measured-vs-expected noise check can expose). Every
frame was processed by the **production Step-4 detection service** (the same
Qt-free `run_detection` code path the GUI executes; SEP engine, $\sigma =
3.2$, repository `parameters.toml`), and the per-frame metrics were assembled
exactly as the production Step-4 window assembles them before calling
`evaluate_frame_qc`: `fwhm_med` is the detector's median radial-profile FWHM of
the brightest $\leq 25$ sources, and `sky_med` / `sky_sigma` are the
detector-reported global background level and RMS (SEP `Background`
`globalback` / `globalrms`) — no auxiliary pixel-level estimate was
substituted. All frames enter QC with `filter = "r"`, `airmass = 1.0`,
`gain_e_per_adu = 1.5`, `rdnoise_e = 5.0`;
the per-source-quality and candidate-count columns are omitted
(`evaluate_frame_qc` skips absent metrics NaN-safely). Because the airmass is
constant, the Kolmogorov seeing model falls back to the night-median FWHM, so
`fwhm_model_ratio` is the frame FWHM over the night median. Decisions use the
shipped conservative `FrameQCThresholds()` defaults (review/fail): FWHM $z$
3.0/4.5, FWHM model ratio
1.3/1.6, sky $z$
4.0/5.5, sky-noise ratio
2.5/6.0, $N_{\rm src}$ $z$
$-3.0$/$-4.5$, elongation
1.22/1.35, depth cost
0.5/1.5 mag (estimated limiting-depth
loss vs the night median from $F_{\rm lim} \propto \sigma_{\rm sky,e}
\cdot \mathrm{FWHM}$). **(a)** Decision matrix, truth class vs
QC decision; green shading marks correct outcomes (clean $\to$ PASS, defect
$\to$ REVIEW or FAIL — the policy is deliberately conservative, reserving FAIL
for frames likely to poison downstream steps), red shading marks wrong ones.
**(b)** The two headline diagnostics added by `evaluate_frame_qc`: the seeing
ratio (`fwhm_model_ratio`, abscissa) against the measured-to-expected sky-noise
ratio (`sky_noise_ratio` $= \sigma_{\rm sky} g \, / \sqrt{B g +
\mathrm{RN}^2}$, log ordinate), colored by truth class with marker shape
encoding the decision (circle PASS, triangle REVIEW, cross FAIL); grey lines
are the REVIEW (dashed) and FAIL (solid) thresholds of the two ratios.

## Per-class outcome (measured, defaults untouched)

| Truth class | N | PASS | REVIEW | FAIL | QC reasons fired (frames) |
|---|---|---|---|---|---|
| clean | 24 | 24 | 0 | 0 | — |
| cloud (ZP $-0.7$ mag) | 5 | 5 | 0 | 0 | — |
| bad seeing (FWHM $\times 1.8$) | 5 | 0 | 0 | 5 | `high_fwhm` (5), `depth_warning` (5) |
| bright sky ($\times 5$) | 5 | 0 | 5 | 0 | `depth_warning` (5), `sky_warning` (5) |
| noisy readout (RN 25, hdr 5 e$^-$) | 5 | 0 | 5 | 0 | `depth_warning` (5), `sky_warning` (5) |

Headline: **24/24 clean frames PASS** and
**15/20 defective frames are flagged** (REVIEW or
FAIL) at the shipped default thresholds.

## Median per-frame metrics by class

| Truth class | $N_{\rm src}$ | FWHM (px) | sky (ADU) | $\sigma_{\rm sky}$ (ADU) | FWHM ratio | noise ratio | depth cost (mag) |
|---|---|---|---|---|---|---|---|
| clean | 120 | 3.54 | 150 | 10.46 | 1.00 | 0.99 | -0.01 |
| cloud (ZP $-0.7$ mag) | 100 | 3.52 | 150 | 10.44 | 0.99 | 0.99 | -0.01 |
| bad seeing (FWHM $\times 1.8$) | 88 | 6.36 | 150 | 10.56 | 1.79 | 1.00 | +0.64 |
| bright sky ($\times 5$) | 91 | 3.52 | 750 | 22.40 | 0.99 | 0.99 | +0.82 |
| noisy readout (RN 25, hdr 5 e$^-$) | 96 | 3.58 | 150 | 19.25 | 1.01 | 1.83 | +0.67 |

## Which check caught what

* **bad seeing** — the only class driven to **FAIL**, by the FWHM checks:
  night-relative robust $z$ (`fwhm_z`, median
  29.8 $\gg$ 4.5) and the
  seeing-model ratio (median 1.79
  $>$ 1.6 FAIL cut); the depth-cost check
  (median +0.64 mag) fired in
  support. Reasons: `high_fwhm` (5), `depth_warning` (5).
* **bright sky** — caught by the sky-level outlier $z$-score (`sky_z`, median
  7303 $>$ 4.0; the enormous
  $z$ reflects the near-zero clean-frame scatter of the synthetic sky) plus
  the depth-cost check (median
  +0.82 mag). Its
  `sky_noise_ratio` stays at 0.99
  because a genuinely brighter sky raises the CCD-equation expectation in
  step — the absolute-consistency and night-relative sky checks are
  complementary by design.
* **noisy readout (header lie)** — the measured-vs-expected sky-noise ratio rose to a median of 1.83$\times$ (clean frames: 0.99) but stayed below the conservative 2.5 review threshold at this sky level; the frames were flagged anyway, by the *night-relative* sky-noise outlier check (`sky_sigma_z` $>$ 4.0 on 5/5 frames — the same measured-noise-anomaly family, referenced to the night's own frames rather than to the CCD equation) and by the depth-cost check (median 0.67 mag of estimated depth loss).
* **cloud** — All five cloud frames **PASS** — the expected, and measured, blind spot. A grey 0.7 mag transparency loss leaves FWHM, sky level, sky noise, and star shapes untouched (median metrics are indistinguishable from clean, see table), so the shape, sky, noise-consistency, and depth checks are *structurally* blind to it: even the depth-cost proxy is unchanged because extinction costs *source* flux, not sky noise. The only metric that responds at all is the detected source count (median $N_{\rm src}$ 100 vs 120 clean — an LF-dependent deficit from the ~0.7 mag slice of stars pushed below the detection limit), and its night-relative robust $z$ (median -0.6) stayed far above the $-3.0$ review cut for two honest reasons: (i) the deficit itself depends on the luminosity function near the limit, and (ii) with 20/44 frames defective, the MAD-based night scale is inflated by the defect population itself (clean-only scale 4.4 counts vs 18.5 whole-night; referenced to clean frames alone the same cloud frames would sit at $z = -6.6$ to $-3.0$, i.e. mostly flaggable). Grey transparency loss is thus only robustly observable relative to a photometric reference; this measured limitation is precisely the motivation for the planned post-photometry transparency-QC stage (frame zeropoint / comparison-star flux monitoring).

The 44 production detections took 22 s of wall time
(single worker). Full per-frame metrics and decisions:
`validation/paper/data/frame_qc/frame_qc_night.csv`.
