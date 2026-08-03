# Figure — detector calibration (Step 0)

**(a)-(c)** The three master frames APEX builds from the night's own
calibration exposures: bias (median 512 DN),
dark (1.0 DN in 60 s) and flat
(unit median). **(d)** The flat's horizontal profile — a
12 per cent
response gradient across the field, which is what the division removes.
**(e)-(g)** The same NGC 6811 $B$ science frame after each operation: raw,
then bias- and dark-subtracted, then flat-fielded. Panels (f) and (g) share a
greyscale so the change is the flat-fielding alone. **(h)** The sky profile
before and after flat-fielding, each normalised to its own median: the
peak-to-peak gradient falls from 12.9 to 1.6 per
cent. All frames are real (Moravian C3-61000, 2x2, night 2026-06-11;
8 bias, 8 darks, 5 flats), shown at
1/12 scale. Numerical validation of these operations is
separate: recovery of injected truth (text), the detector constants (Fig. 3),
pixel-for-pixel agreement with ccdproc (Fig. 4) and reproduction on two more
cameras (Fig. 5).
