# Figure 7 — Cross-catalog reference validation: a Gaia BP faint bias

**Figure 7.** Cross-matching 1928 NGC 6811 stars to Pan-STARRS 1 (PS1) — a
deep ($g \sim 22$), Gaia-independent survey — separates a *reference*
systematic from a *measurement* systematic, band by band. Each series is a
colour-term fit ($\Delta = {\rm zp} + {\rm ct}\,(g-r)$) against PS1 $g$,
defined on bright well-measured stars, with the
14.5-15.5 mag anchor pinned to zero; points are
binned medians (error bars $= {\rm MAD}/\sqrt N$). **(a)** In $B$ the
Gaia-transformed reference (Pancino et al. 2022, from BP-RP) drifts
$\Delta_{\rm faint} = +0.022$ mag faint-ward against PS1 —
larger than APEX's own $B$ ($\Delta_{\rm faint} = +0.010$
mag): the $B$ faint drift is dominated by the *reference*. **(b)** In $V$ the
pattern reverses — the Gaia reference (derived from the robust $G$ band, not BP)
is flat (+0.007 mag) while APEX carries a small
$+0.024$ mag residual, the last trace of aperture sky
over-subtraction.

**Interpretation.** The dominant systematic *swaps catalog between bands*,
proving the two are independent. In $B$ the culprit is Gaia's BP channel: BP is
slitless-prism spectrophotometry, and for faint (intrinsically blue-faint)
sources the BP flux is background-dominated and systematically mismeasured
(Riello et al. 2021), biasing every magnitude transformed from BP-RP. APEX
aperture $B$, compared to the independent PS1 scale, is flat to $\sim0.01$ mag
— so the "B-band faint drift" originally seen against Gaia is chiefly a
*reference* artifact, not an APEX error. In $V$, where the reference is clean,
only a residual $\sim0.02$ mag APEX sky effect remains — at the systematic
floor. Practical consequence: bright stars anchor the zeropoint safely with
Gaia, but faint-end *validation* should use a $G$-based or PS1 reference, never
BP-RP-transformed magnitudes.
