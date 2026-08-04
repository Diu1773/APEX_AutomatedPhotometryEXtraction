# Figure — per-step preprocessing cross-check vs Python ccdproc

**(a)** Numeric audit against the independent Python `ccdproc` package
(`ccdproc` is not the IRAF task of the same name). The table reports the
maximum absolute difference and robust scatter for each stage; six rows are
bit-identical at every pixel, while the full chain leaves only
8.6e-04 DN (robust σ =
2.6e-05 DN) from float32 rounding.
**(b)** The end-to-end difference beside the read-noise (3.45 DN) and
sky-shot-noise (41 DN) scales. Inputs: 8 bias,
8 darks (60 s), 5 flats, one 60 s NGC 6811 $B$
light (Moravian C3-61000, 2×2, night 2026-06-11); reference
astropy ccdproc 2.5.1. The cosmetic (L.A.Cosmic + hot-pixel) stage is disabled
here because it repairs pixels by design and is validated separately by
injection. Generator: `calib_crosscheck_ngc6811.py`.
