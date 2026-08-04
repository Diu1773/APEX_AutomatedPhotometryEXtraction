# Figure 3 — detector constants from three reductions

**Figure 3.** Detector constants obtained from the same Moravian C3-61000
(Sony IMX455, 2×2) calibration set by three reduction paths: APEX, the Python
`ccdproc` package, and the IRAF `ccdproc` task. **(a)** The table reports the
values used by the error model: flat-pair photon-transfer gain
(0.681 ± 0.014 e⁻/ADU), bias-pair read noise (2.35 e⁻), and the
linear dark-ladder slope (0.0077 e⁻/s). **(b)** The grouped bars show one
physical scale per quantity, so the agreement is visible without putting
incommensurate units on a single axis. The three values are identical at the
shown precision because bias/dark subtraction is additive and cancels in the
flat-pair variance, while all reductions use the same median dark ladder; this
is an agreement result, not a claim that the pixel arrays are byte-identical.
The pixel-level residuals, including the independent IRAF flat-normalisation
and full-chain differences, are reported in Figure 4. No FITS-header, vendor,
or laboratory gain is used in this comparison. Inputs: 8 bias frames, 5 B
flats, and the 10–480 s dark ladder from 2026-06-11.
