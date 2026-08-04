# Figure — absolute detector-preprocessing products

The table lists the absolute output level of each calibration product as the
full-frame median ± robust σ (1.4826 × MAD), rather than subtracting APEX from
either reference.  Columns are the production APEX reduction, an independent
Python `ccdproc` reduction, and an independent IRAF `ccdproc` task run through
PyRAF.  Rows are the master bias, the 60-s master dark, the unit-median
normalised master flat, and the science frame after the complete bias–dark–flat
chain.  The displayed Python values are 512 ± 1.483 DN (bias), 1 ± 2.224 DN
(60-s dark), 1 ± 0.04428 DN (flat), and 633.5 ± 31.39 DN (full chain), the
same printed values shown in the APEX column.  The IRAF values are 512 ± 1.483,
1 ± 2.224, 1.001 ± 0.04431, and 633 ± 31.37 DN, respectively.  The separate
pixel-residual audit remains in `data/iraf_preproc_stats.json` and is not used
as the plotted quantity.  Inputs are eight bias frames, eight 60-s
darks, five B flats, and one 60-s NGC 6811 B science frame from 2026-06-11
(Moravian C3-61000, 2×2).  Cosmetic repair was disabled because this table
isolates detector calibration arithmetic.
