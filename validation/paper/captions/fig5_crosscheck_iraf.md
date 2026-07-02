# Figure 5 — External cross-check on real data: APEX vs IRAF/DAOPHOT

**Figure 5.** Agreement between APEX and IRAF `phot` (DAOPHOT, run via PyRAF)
measuring the **same 278 stars** at the **same fixed sky coordinates** on a real
observed frame of the open cluster **NGC 457** in the **g band**
(`pp_-0016-gfilter_20240907.fit`, observed 2024-09-07; FWHM $= 5.31$ px, gain
$= 0.049524$ e$^-$/ADU, read noise
$= 2.78$ e$^-$; IRAF zmag $= 25$). APEX measures with
its production forced-aperture photometry; IRAF measures the identical
coordinates with `phot`. **(a)** APEX vs IRAF instrumental magnitude after a
single median-zeropoint alignment (the grey line is $y=x$); the constant
$\approx-0.060$ mag offset absorbs the aperture-correction / gain /
exposure difference between the two flux scales. **(b)** Zeropoint-aligned
residual $\Delta$ versus magnitude, with the median (solid) and a shaded
$\pm$MAD band.

**Agreement metrics** (from `fixed_comparison.csv`, matching the committed
`fixed_summary.json`):

| Quantity | Value |
|---|---|
| $N$ (matched, SNR $>$ 20) | 278 |
| Zeropoint offset | -0.0602 mag |
| median $\Delta$ | +0.0000 mag |
| MAD | 0.0043 mag |
| RMS$(\Delta)$ | 0.0083 mag |
| Pearson $r$ | 0.99998 |
| position RMS | 0.0004 px |
| fraction $|\Delta| < 0.01$ mag | 0.950 |
| fraction $|\Delta| < 0.02$ mag | 0.996 |

**Verdict.** MAD $= 0.0043$ mag and Pearson $r = 0.99998$ over 278
real stars: APEX and IRAF **agree at photometric grade**.
Because IRAF/DAOPHOT shares no code with APEX and this is a real observed frame
(not a synthetic model), this is the most stringent single accuracy test in the
suite; the milli-magnitude agreement confirms APEX's aperture photometry is
correct on-sky, not just self-consistent. Complements Fig 4 (independent `sep`
engine on a synthetic frame with a known truth catalog).
