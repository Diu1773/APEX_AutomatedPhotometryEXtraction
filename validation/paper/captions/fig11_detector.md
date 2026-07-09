# Figure 11 — Detector characterisation from the data

**Figure 11.** Gain, read noise and dark current of the Moravian C3-61000
(Sony IMX455, 2×2 binned) measured directly from APEX's own calibration frames.
**(a)** Photon-transfer relation: the variance of same-level flat-pair
differences (which cancel fixed-pattern noise) versus signal, over
12 clean pairs. The slope is 1/gain, giving
gain = 0.681 ± 0.014 e⁻/ADU with read noise 2.35 e⁻. The
measured value is consistent with the IMX455 laboratory value
(Alarcón et al. 2023, 0.763 e⁻/ADU, native resolution) and the
vendor full-well specification (>50 ke⁻ over the 16-bit range ⇒ ≈0.76
e⁻/ADU); the small difference is the 2×2 binning. The gain is *measured*, not
taken from the FITS header: for this camera the MaxIm/ASCOM `EGAIN` keyword is a
factor of ≈16 too small (a documented 12-bit→16-bit ADC left-shift), so it is
not used. **(b)** Dark current from the source-free background versus exposure
across a 10–480 s ladder: linear (R² = 0.9978, residuals in the lower panel),
slope 0.0077 e⁻/s at +5 °C. These measured values anchor APEX's photometric
error model to the detector's real physics.
