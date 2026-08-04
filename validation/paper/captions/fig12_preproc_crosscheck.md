# Figure — absolute detector-preprocessing products

The table lists the absolute output level of each calibration product as the
full-frame median ± robust σ (1.4826 × MAD), rather than subtracting APEX from
either reference.  Columns are the production APEX reduction, an independent
Python `ccdproc` reduction, and an independent IRAF `ccdproc` task run through
PyRAF.  Rows are the master bias, the 60-s master dark, the unit-median
normalised master flat, and the science frame after the complete bias–dark–flat
chain.  The Python column agrees with APEX at the displayed precision; the
separate pixel-residual audit remains in `data/iraf_preproc_stats.json` and is
not used as the plotted quantity.  Inputs are eight bias frames, eight 60-s
darks, five B flats, and one 60-s NGC 6811 B science frame from 2026-06-11
(Moravian C3-61000, 2×2).  Cosmetic repair was disabled because this table
isolates detector calibration arithmetic.
