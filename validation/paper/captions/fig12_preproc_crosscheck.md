# Figure — per-step preprocessing cross-check vs ccdproc

**(a)** Maximum absolute APEX−ccdproc difference at every calibration stage.
Master bias/dark/flat construction and the three individual corrections are
bit-identical; open circles are labelled measured zeros and placed on a display
floor only so that all stages remain visible. The full chain differs by at most
8.6e-04 DN (robust σ =
2.6e-05 DN) because of float32 rounding.
**(b)** The end-to-end difference beside the read-noise (3.45 DN) and
sky-shot-noise (41 DN) scales. Inputs: 8 bias,
8 darks (60 s), 5 flats, one 60 s NGC 6811 $B$
light (Moravian C3-61000, 2×2, night 2026-06-11); reference
astropy ccdproc 2.5.1. The cosmetic (L.A.Cosmic + hot-pixel) stage is disabled
here because it repairs pixels by design and is validated separately by
injection. Generator: `calib_crosscheck_ngc6811.py`.
