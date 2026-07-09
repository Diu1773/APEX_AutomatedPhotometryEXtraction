# Figure 11 — Detector characterisation from the data

**Figure 11.** Gain, read noise, and dark current of the Moravian C3-61000
(Sony IMX455, 2×2 binned) measured directly from APEX's own calibration frames.
**(a)** Photon-transfer curve: the variance of same-level flat-pair differences
(which cancel PRNU and vignette) versus signal, over 12 clean pairs
spanning 19120–36510 ADU. The slope gives
gain = 0.681 ± 0.014 e⁻/ADU (stored pixel); the read noise from a
bias-pair difference is 2.35 e⁻ (stored) / 1.18 e⁻
(native), consistent with the IMX455 laboratory value (Alarcón et al. 2023). The
FITS-header EGAIN (0.0495 e⁻/ADU, the nominal MaxIm max-gain value) implies
a 14× steeper line and is ruled out at 46σ — so the gain must be measured, not
read from the header. **(b)** Dark current from the source-free background versus
exposure across a 10–480 s ladder: linear (R² = 0.9978), slope
0.0077 e⁻/s at +5 °C. These measured values anchor APEX's photometric error
model to the detector's real physics.
