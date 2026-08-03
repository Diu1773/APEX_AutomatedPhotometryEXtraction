# Figure — detector calibration (Step 0)

**(a)** A real M13 $V$ frame before and after APEX Step 0: the optical vignette
and the fixed-pattern structure are removed. **(b)** Controlled truth recovery —
known bias, dark and flat are injected into a synthetic frame and inverted
through Step 0; the residual against the true science frame has systematic
offset -0.004 DN and scatter MAD 3.13 DN, the latter equal
to the injected read-noise floor (3.0 DN). The master bias is
recovered to RMS 0.97 DN and the unit-median flat reduces the
corner vignette residual from 10.6 to 0.7
DN. **(c)** The optional cosmetic stage (L.A.Cosmic via astroscrappy) removes
100 per cent of injected cosmic-ray pixels and
94 per cent of hot pixels while touching
0.000 per cent of star-core pixels; aperture
fluxes shift by 0.00 mmag (median).
**(d)** The same raw frames reduced by the independent AstralImage/AIPPI engine:
across 19 datasets in 9 bands the calibrated frames are bit-identical
(difference RMS and maximum both 0 DN), and the master frames agree to
$\leq$0.30 DN RMS.
Generator: `fig_calibration_step0.py`; every number is recomputed at figure time.
