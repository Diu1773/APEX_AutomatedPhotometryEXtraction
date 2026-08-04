# Figure — detector calibration (Step 0)

**(a)** Master bias image and **(b)** its pixel-value distribution, including
the median and within-frame scatter. **(c)** Master dark image and **(d)** its
distribution; the long positive tail is retained in the displayed percentile
stretch. **(e)** Master flat image and **(f)** its horizontal profile, showing a
12 per cent response
gradient. **(g)** Raw NGC 6811 $B$ science frame and **(h)** the calibrated
frame after bias/dark subtraction and flat division. The small inset in (h)
shows the sky profile before and after flat-fielding, normalised to each
median; the peak-to-peak gradient falls from 12.9 to 1.6
per cent. All frames are real
(Moravian C3-61000, 2x2, night 2026-06-11; 8 bias,
8 darks, 5 flats), shown at 1/12
scale. Numerical equivalence to the independent Python `ccdproc` package is
tested separately in Fig. 4.
