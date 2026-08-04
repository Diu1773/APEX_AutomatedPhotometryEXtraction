# IRAF detector-preprocessing cross-check

Figure 4 uses the same NGC 6811 raw set for all three reductions: eight bias
frames, eight 60-s dark frames, five B flats, and one 60-s B science frame from
2026-06-11 (Moravian C3-61000, 2x2). APEX and Python `ccdproc` use the existing
generator `calib_crosscheck_ngc6811.py`. The Python path was run in the repository
Python 3.12.3 environment with Astropy `ccdproc` 2.5.1 (Craig et al. 2015,
ASCL:1510.007). FITS arrays are loaded by the shared APEX loader and then wrapped
as Astropy `CCDData`; the calibration arithmetic itself uses `ccdproc.combine`,
`subtract_bias`, `subtract_dark`, and `flat_correct`. The independent IRAF
reduction was run in WSL with PyRAF 2.2.3.dev9 and IRAF `ccdred` (Tody 1986,
1993). The retained run log does not contain the IRAF distribution build
number, so no IRAF version is inferred here:

1. `zerocombine` and `darkcombine`, median combine with rejection disabled;
2. `ccdproc` bias subtraction on the dark library, followed by `darkcombine`;
3. `ccdproc` bias/dark correction on each flat, `flatcombine` with median
   scaling, and unit-median normalisation;
4. `ccdproc` bias + dark + flat correction on the science frame.

Overscan, trim, cosmic-ray, and hot-pixel repair were disabled so the comparison
tests detector-calibration arithmetic only. The Python package citation documents
the software; it is not an external validation of equivalence to IRAF. The
equivalence claim in the manuscript is the same-raw cross-check reported here.
The compact pixel-residual audit is
`data/iraf_preproc_stats.json`; the absolute product values used in Figure 4 are
recorded separately in `data/preproc_absolute_summary.json`. The large IRAF FITS
intermediates are ignored by the repository. Figure generation is performed by
`fig12_preproc_crosscheck.py`; the plotted quantity is the absolute product
level, not an APEX-minus-reference residual.
