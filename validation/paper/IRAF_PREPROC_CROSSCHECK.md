# IRAF detector-preprocessing cross-check

Figure 4 uses the same NGC 6811 raw set for all three reductions: eight bias
frames, eight 60-s dark frames, five B flats, and one 60-s B science frame from
2026-06-11 (Moravian C3-61000, 2×2). APEX and Python `ccdproc` use the existing
generator `calib_crosscheck_ngc6811.py`. The independent IRAF reduction was run
in WSL with PyRAF 2.2.3.dev9 and IRAF `ccdred`:

1. `zerocombine` and `darkcombine`, median combine with rejection disabled;
2. `ccdproc` bias subtraction on the dark library, followed by `darkcombine`;
3. `ccdproc` bias/dark correction on each flat, `flatcombine` with median
   scaling, and unit-median normalisation;
4. `ccdproc` bias + dark + flat correction on the science frame.

Overscan, trim, cosmic-ray, and hot-pixel repair were disabled so the comparison
tests detector-calibration arithmetic only. The compact measured summary is
`data/iraf_preproc_stats.json`; the large IRAF FITS intermediates are ignored by
the repository. Figure generation is performed by
`fig12_preproc_crosscheck.py`.
