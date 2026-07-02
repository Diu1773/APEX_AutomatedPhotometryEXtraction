# Figure 1 — Detection completeness of the APEX pipeline

**Figure 1.** Detection completeness of the APEX forced-photometry pipeline
measured with artificial-star tests. Blue points show the fraction of injected
artificial stars recovered as a function of injected magnitude, binned in
0.25 mag intervals; error bars are exact Wilson 95% binomial confidence
intervals on each bin's recovered fraction. The vermillion curve is the fitted
completeness model, a logistic function
$C(m) = \mathrm{expit}\!\left[(m_{50}-m)/w\right]$, where $m_{50}$ is the
50% completeness depth and $w$ the transition width; its parameters and their
uncertainties are obtained from a cluster bootstrap over the 12 injection trials
(500 resamples) and are *not* re-fit to the binned points. The dashed grey line
and shaded band mark $m_{50}$ and its bootstrap 95% CI; the horizontal dotted
line marks $C=0.5$. Dash-dotted orange lines indicate the 90% and 10%
completeness depths $m_{90}$ and $m_{10}$. The lower strip shows the number of
eligible injections $N_{\mathrm{elig}}$ per magnitude bin, confirming roughly
uniform sampling across the depth range.

The pipeline reaches a 50% completeness depth of
$m_{50} = 17.56\,^{+0.08}_{-0.07}$ mag (bootstrap 95% CI
$17.48$–$17.63$ mag), with a logistic transition width of
$w = 0.33$ mag ($0.28$–$0.38$ mag, 95% CI). The completeness stays above 90%
brighter than $m_{90} = 16.83$ mag and falls below 10% fainter than
$m_{10} = 18.28$ mag — a sharp $\sim$1.4 mag transition. The measurement uses
$N = 840$ injected stars (837 eligible) across 12 independent trials. The
logistic model tracks the recovered fraction well through the transition; the
few-percent shortfall of the brightest bins from unit completeness reflects a
small magnitude-independent loss (injections landing on defective pixels or
close blends) that the two-parameter logistic form does not attempt to absorb,
and does not affect the recovered $m_{50}$ or width.
