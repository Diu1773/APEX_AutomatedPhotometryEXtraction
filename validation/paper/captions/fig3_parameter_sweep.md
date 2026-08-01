# Figure — parameter & observing-conditions sensitivity (rebuilt 2026-08-02)

**Data: synthetic only.** APEX synthetic generator
(`apex.benchmark.synthetic_frame`), seed 20260702;
1024x1024 px,
435 field stars (star density matches the old 640 px
/ 170-star frame), gain 1.5 e-/ADU, read noise 5 e-, zero point 25.0, r band.
Every point is an
artificial-star benchmark through the production code path (empirical-PSF
injection -> Step-4 detection -> forced aperture photometry;
`apex.benchmark.runner.run_benchmark`) with **30 trials x
40 injected stars = 1200 injections per
point** (the old figure used 4 trials; its sawtooth was Monte-Carlo noise, see
verdict below). Error bars are 95% CIs from a cluster bootstrap over trials
(1000 resamples for MAD/RMSE; the depth CI comes from the production
completeness fit's own 500-sample cluster bootstrap). Reproduce:
`validation/paper/fig3_parameter_sweep.py` (distilled numbers in
`validation/paper/data_parameter_sweep/`).

**(a) Aperture radius — scatter minimum at 0.80xFWHM.** One fixed frame
with FWHM 6.0 px (detected 6.01 px), sky 150
ADU, injections 14.5-19.0 mag. The forced-aperture radius is
`max(min_r_ap_px = 4 px, scale x FWHM)`. **Why FWHM 6 px:** on the old 3.5 px
frame every requested scale <= 1.1xFWHM collapsed onto the same 4 px floor
radius, and the old panel drew those duplicate measurements as a "flat"
segment — a clamp artifact, not physics. On this frame the floor sits at
0.67xFWHM, so 0.7-2.5xFWHM are all distinct apertures.
The one deliberately clamped request (0.5xFWHM -> 4 px, open marker) lands on
the floor. Scatter (MAD of the forced-photometry magnitude error) is minimized
at 0.80xFWHM with MAD 0.0513 mag
[0.0464, 0.0557]; larger apertures admit sky noise
(0.181 mag at 2.5xFWHM). The grey band marks radii the
production pipeline cannot reach (below the 4 px floor); the dashed line is
the default scale 0.8.

**(b) Position RMSE is aperture-independent.** Median 0.26 px across
all unclamped apertures; the centroid comes from Step-4 detection, which does
not use the photometry aperture.

**(c,d) Sky background** (frame regenerated per level, injections 14.0-20.0
mag, log x): sky 50 -> 1000 ADU makes the 50% completeness depth
18.38 -> 17.06 mag
(1.32 mag shallower) and the
scatter 0.036 -> 0.088 mag: a brighter sky
raises the shot-noise floor.

**(e,f) Seeing** (frame regenerated per FWHM, sky 150 ADU): FWHM
3 -> 6 px makes the depth
18.07 -> 16.96 mag
(1.11 mag shallower) and the
scatter 0.042 -> 0.060 mag: the same flux
spreads over more pixels, lowering peak SNR. The grid stays inside the
production FWHM-QC window (`[fwhm] px_min = 3.0`): below it the frame-FWHM
estimator is outside its envelope (a true 2.5 px frame was measured as 8.5 px
from 5 stars), so a 2.5 px point would not be a pipeline measurement.

**Sawtooth verdict (old R5 defect).** With 4 trials the scatter curves wobbled
non-monotonically. At 30 trials the largest remaining downward step between
adjacent scatter points is 0.0 mmag (sky) / 1.6
mmag (seeing) against median CI half-widths of 8.8 /
5.4 mmag — the old sawtooth was Monte-Carlo noise, not
structure. What survives at 30 trials is consistent with monotone trends
within the CIs.

**Detection-threshold companion sweep (not drawn) — the old "threshold-robust"
claim was a no-op artifact.** The old run overrode only `detect_sigma`, which
the per-filter `detect_sigma_r` in parameters.toml silently masks, so all five
"swept" thresholds ran at 3.2-sigma and returned bit-identical results; the
"depth constant to < 0.01 mag over 2-6 sigma" sentence measured nothing. With
the per-filter key overridden too (each point's `sigma_used` is verified equal
to the requested value), the threshold sets the depth directly: m50
18.36 at 2 sigma ->
17.32 at 6 sigma
(1.04 mag), while false detections stay
at <= 2 per 30-trial point on this clean low-crowding
field. Lowering the threshold is therefore free depth *on this field*; the
real-frame cost side (false-detection contamination vs threshold) is measured
on real frames in the detection-threshold section (its own figure).
Production default: 3.2 sigma.

## Sweep 1 — aperture (fixed frame, FWHM 6.0 px)

| requested_scale | r_ap_px | mad | mad_lo | mad_hi | rmse |
|---|---|---|---|---|---|
| 0.50 | 4.00 | 0.0541 | 0.0486 | 0.0579 | 0.255 |
| 0.70 | 4.21 | 0.0521 | 0.0473 | 0.0586 | 0.255 |
| 0.80 | 4.81 | 0.0513 | 0.0464 | 0.0557 | 0.255 |
| 0.90 | 5.41 | 0.0550 | 0.0504 | 0.0605 | 0.255 |
| 1.00 | 6.01 | 0.0601 | 0.0541 | 0.0683 | 0.255 |
| 1.10 | 6.62 | 0.0648 | 0.0585 | 0.0723 | 0.255 |
| 1.20 | 7.22 | 0.0736 | 0.0653 | 0.0817 | 0.255 |
| 1.40 | 8.42 | 0.0877 | 0.0812 | 0.0946 | 0.255 |
| 1.70 | 10.22 | 0.1124 | 0.1004 | 0.1239 | 0.255 |
| 2.00 | 12.03 | 0.1236 | 0.1135 | 0.1369 | 0.255 |
| 2.50 | 15.04 | 0.1810 | 0.1619 | 0.2042 | 0.255 |

## Sweep 2 — sky background (frame per value)

| background | m50 | m50_lo | m50_hi | mad | mad_lo | mad_hi |
|---|---|---|---|---|---|---|
| 50 | 18.380 | 18.273 | 18.482 | 0.0360 | 0.0325 | 0.0408 |
| 100 | 18.055 | 17.933 | 18.168 | 0.0411 | 0.0358 | 0.0471 |
| 150 | 17.908 | 17.788 | 18.012 | 0.0450 | 0.0382 | 0.0503 |
| 225 | 17.681 | 17.581 | 17.761 | 0.0521 | 0.0468 | 0.0590 |
| 300 | 17.583 | 17.489 | 17.660 | 0.0652 | 0.0558 | 0.0749 |
| 450 | 17.366 | 17.268 | 17.452 | 0.0681 | 0.0563 | 0.0785 |
| 600 | 17.236 | 17.156 | 17.307 | 0.0776 | 0.0670 | 0.0846 |
| 800 | 17.175 | 17.095 | 17.240 | 0.0806 | 0.0695 | 0.0936 |
| 1000 | 17.060 | 17.009 | 17.111 | 0.0883 | 0.0769 | 0.1010 |

## Sweep 3 — seeing / PSF FWHM (frame per value)

| fwhm_px | m50 | m50_lo | m50_hi | mad | mad_lo | mad_hi |
|---|---|---|---|---|---|---|
| 3.0 | 18.069 | 17.959 | 18.158 | 0.0415 | 0.0367 | 0.0471 |
| 3.5 | 17.908 | 17.788 | 18.012 | 0.0450 | 0.0382 | 0.0503 |
| 4.0 | 17.705 | 17.622 | 17.797 | 0.0486 | 0.0419 | 0.0552 |
| 4.5 | 17.556 | 17.474 | 17.634 | 0.0470 | 0.0410 | 0.0519 |
| 5.0 | 17.366 | 17.304 | 17.425 | 0.0473 | 0.0429 | 0.0537 |
| 5.5 | 17.148 | 17.073 | 17.221 | 0.0570 | 0.0520 | 0.0624 |
| 6.0 | 16.959 | 16.925 | 16.989 | 0.0596 | 0.0531 | 0.0654 |

## Sweep 4 — detection threshold (companion, not drawn)

| detect_sigma | sigma_used | m50 | m50_lo | m50_hi | false_detections |
|---|---|---|---|---|---|
| 2.0 | 2.0 | 18.364 | 18.266 | 18.475 | 1 |
| 3.2 | 3.2 | 17.908 | 17.788 | 18.012 | 1 |
| 4.0 | 4.0 | 17.675 | 17.577 | 17.756 | 2 |
| 5.0 | 5.0 | 17.433 | 17.330 | 17.541 | 2 |
| 6.0 | 6.0 | 17.323 | 17.229 | 17.412 | 2 |
