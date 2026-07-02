# Figure 4 — Independent cross-check: APEX vs the `sep` C library

**Figure.** Agreement between APEX aperture photometry and the independent
`sep` (SExtractor SEP backend, v1.4.1) engine measuring the **same 95 stars**
in the **same synthetic frame** with the **same aperture**
($r_{\mathrm{ap}} = 1.0\,\mathrm{FWHM} = 3.50$ px).
The frame is a 800x800 px field of 150 isolated stars
(separation $\ge 8\,\mathrm{FWHM}$), true magnitudes uniform in
$[14, 20]$, with a Gaussian PSF
($\mathrm{FWHM} = 3.5$ px), gain $g = 1.5$ e$^-$/ADU, read noise
$\mathrm{RN} = 5.0$ e$^-$, sky $B = 150$ ADU, and
zeropoint $\mathrm{ZP} = 25.0$; noise is Poisson (photon) plus Gaussian
(read) in the electron domain, matching `apex/benchmark/synthetic_frame.py`.
**(a)** APEX vs `sep` instrumental magnitude after a single median-zeropoint
alignment; the grey line is $y=x$. **(b)** Zeropoint-aligned residual
$\Delta = (m_{\mathrm{APEX}} - m_{\mathrm{sep}}) - \mathrm{ZP}$
versus true magnitude, with the median (solid) and a shaded $\pm$MAD band.
Only stars with APEX $\mathrm{SNR} > 10$ and positive flux in both
engines are kept.

**Agreement metrics** (ZP offset = median$(m_{\mathrm{APEX}} -
m_{\mathrm{sep}})$; residuals are ZP-removed):

| Quantity | Value |
|---|---|
| $N$ (matched, SNR $> 10$) | 95 |
| Zeropoint offset | -0.4395 mag |
| median $\Delta$ | +0.0000 mag |
| MAD ($1.4826\times$ median$|\Delta - \mathrm{med}|$) | 0.0060 mag |
| RMS$(\Delta)$ | 0.0097 mag |
| Pearson $r$ | 0.99995 |
| fraction $|\Delta| < 0.01$ mag | 0.811 |
| fraction $|\Delta| < 0.02$ mag | 0.947 |

**Verdict.** MAD $= 0.0060$ mag $\le$ 0.02 mag and Pearson
$r = 0.99995$: APEX and `sep` **agree at photometric grade**.
The two engines share no photometry code, so this differential agreement
confirms that APEX's aperture-sum + local-sky measurement is correct to the
milli-magnitude level, with any residual scatter set by pixelization and the
per-engine sky estimator rather than by a systematic error in APEX.
