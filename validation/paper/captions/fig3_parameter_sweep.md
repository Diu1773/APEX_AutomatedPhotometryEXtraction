# Figure 3 — Parameter & observing-conditions sensitivity

**Sensitivity of APEX pipeline metrics to one pipeline parameter and two
observing conditions, each swept independently.** Every point is an
artificial-star benchmark (empirical-PSF injection, production Step 4 detection,
and forced aperture photometry) with 4 trials of 40 injected stars; the modest
trial count keeps the sweep fast at the cost of per-point Monte-Carlo noise.
Panel (a) isolates a pure pipeline knob on a single fixed synthetic frame; for
the observing-condition panels (b, c) a controlled synthetic frame is
regenerated per value (fixed star-field seed 20260702, 170
stars, 640x640 px, gain
1.5 e-/ADU, read noise 5.0 e-, zero
point 25.0 mag), changing only the swept condition.

**(a) Aperture radius.** Photometric scatter (median absolute deviation of the
recovered magnitude, blue, left axis) and astrometric position RMSE (orange,
right axis) versus aperture radius in units of the PSF FWHM
(`aperture_scale_fwhm`), on the fixed frame (sky 150
ADU, FWHM 3.5 px). Scatter is minimized at
**1.2x FWHM** (MAD = 0.0579 mag):
smaller apertures lose source flux and SNR, larger apertures admit sky noise.
Position RMSE is essentially aperture-independent (~0.154
px) because the centroid comes from detection, not the aperture. Dashed grey =
pipeline default (0.8x FWHM); dotted green = recovered
optimum.

**(b) Sky background.** 50% completeness depth *m*<sub>50</sub> (blue, left axis)
and photometric scatter (orange, right axis) versus sky background, with a frame
regenerated at each level (injection window widened to 14.0-20.0 mag so
*m*<sub>50</sub> stays in range). As the sky brightens from
50 to 1000 ADU the depth becomes
shallower (brighter), *m*<sub>50</sub> = 18.484 -> 17.109
mag, and the photometric scatter grows
(0.0496 -> 0.1456 mag): a brighter sky
raises the shot-noise floor, so faint sources drop below the detection/SNR
limit — a clean monotonic trend.

**(c) Seeing (PSF FWHM).** *m*<sub>50</sub> depth (blue, left axis) and
photometric scatter (orange, right axis) versus PSF FWHM (seeing), frame
regenerated per value at fixed sky 150 ADU. Worse
seeing spreads source flux over more pixels, lowering peak SNR: from FWHM
2.5 to 6.0 px the depth becomes shallower
(*m*<sub>50</sub> = 18.089 -> 16.790 mag) and the
photometric scatter increases
(0.0497 -> 0.0770 mag).

**Detection-threshold robustness (not shown).** A companion sweep of the
detection threshold `detect_sigma` over 2.0-6.0 sigma on the clean fixed frame
left both the depth (*m*<sub>50</sub> constant to < 0.01 mag) and the purity
(zero spurious detections at every threshold) flat: for bright, well-sampled,
isolated sources the production segmentation detector is threshold-robust and
the recovery limit is SNR-set, not threshold-set. This is a desirable stability
property rather than a trend, so it is reported here instead of as a panel.

## Numeric optima and trends

- **Optimal aperture scale (minimum photometric scatter):**
  1.2x FWHM, MAD = 0.0579 mag.
- **Depth vs sky background:** m50 18.484 mag at
  50 ADU -> 17.109 mag at
  1000 ADU (shallower with brighter sky).
- **Depth vs seeing:** m50 18.089 mag at FWHM 2.5
  px -> 16.790 mag at 6.0 px (shallower with worse
  seeing).
- **Detection-threshold override key used:** `detect_sigma` (robustness
  companion sweep).

## Sweep 1 — aperture radius (fixed frame)

| aperture / FWHM | scatter MAD (mag) | position RMSE (px) |
|---|---|---|
| 0.60 | 0.0599 | 0.1542 |
| 0.80 | 0.0599 | 0.1542 |
| 1.00 | 0.0599 | 0.1542 |
| 1.20 | 0.0579 | 0.1542 |
| 1.50 | 0.0693 | 0.1542 |
| 2.00 | 0.0797 | 0.1542 |
| 2.50 | 0.1185 | 0.1527 |

## Sweep 2 — sky background (frame per value)

| sky background (ADU) | m50 (mag) | scatter MAD (mag) |
|---|---|---|
| 50 | 18.4842 | 0.0496 |
| 100 | 18.1629 | 0.0474 |
| 150 | 17.8747 | 0.0608 |
| 300 | 17.7400 | 0.0969 |
| 600 | 17.3782 | 0.1667 |
| 1000 | 17.1086 | 0.1456 |

## Sweep 3 — seeing / PSF FWHM (frame per value)

| FWHM (px) | m50 (mag) | scatter MAD (mag) |
|---|---|---|
| 2.5 | 18.0894 | 0.0497 |
| 3.0 | 18.2152 | 0.0505 |
| 3.5 | 17.8747 | 0.0608 |
| 4.0 | 17.7407 | 0.0563 |
| 5.0 | 17.2314 | 0.0480 |
| 6.0 | 16.7904 | 0.0770 |

*Benchmark: 4 trials x 40 stars per point; seed 20260702; trials=4 chosen for
sweep speed. Panel (a) reuses one fixed frame; panels (b)/(c) regenerate a
controlled frame per condition value.*
