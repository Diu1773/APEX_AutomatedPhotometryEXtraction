# Figure — detector preprocessing comparison

**(a)** Pixel-level comparison of the APEX reduction with two independent
implementations applied to the same NGC 6811 raw set: the Python `ccdproc`
package and the IRAF `ccdproc` task (run through PyRAF). Each entry is
`max |APEX − reference| / robust σ` in DN. The bias and dark masters are bit
identical for both references. The flat comparison uses the same unit-median
normalisation; the remaining 2.78e-03 DN maximum
is below the displayed detector scales. **(b)** The maximum differences are
shown as grouped bars; exact zeros are placed at a display floor of
1e-06 DN, while the table retains the measured zero. Dotted and dashed
lines mark the 3.5 DN read noise and 41 DN sky-shot-noise scales. The full
pipeline difference is 8.61e-04 DN for Python `ccdproc` and
34.39 DN for IRAF `ccdproc`; the latter reflects the independent
IRAF combination/flat-correction path, not a photometry comparison. Inputs:
8 bias, 8 darks (60 s), 5 flats and one 60 s B-band light, Moravian C3-61000,
night 2026-06-11. Cosmetic repair was disabled because it intentionally
changes pixels and is validated separately.
