# Figure — Detection threshold: contamination measured from the frame itself

**Figure X.** Spurious-detection contamination as a function of the detection
threshold, measured without any external catalogue. Re-running the *same*
detector on the sign-flipped, background-subtracted frame counts the
noise-origin spurious detections directly, because the noise is symmetric about
zero while real sources are positive (Serra, Jurek & Flöer 2012; applied to
optical images by Molino et al. 2014). The measurement needs no WCS solution and
no injection, so it runs inside the detection step, before plate solving, at
3.6 per cent of that step's cost (0.21 s on a 4788×3194 frame).

**(a)** Estimated contamination $N_-/N_+$ for five real single exposures — two
globular clusters (M13, M3) and two open clusters (M67, NGC 6811), in $B$, $R$
and $g'$ — against the threshold at which they were extracted. Open carets are
upper limits: frames in which the sign-flipped image yielded *no* detection at
all, plotted at $1/N_+$ and excluded from the connecting line, since joining them
would draw a slope set by $N_+$ rather than by contamination. Large open symbols
mark each frame's own floor, the lowest threshold at which contamination stays
below 5 per cent. **Those floors are 1.5, 1.5, 1.8, 2.0 and 2.2** — they differ
by frame, and the two M13 curves show the split occurring between two filters of
the *same cluster on the same night*, so it is not a property of the target.
At the shipped default of $3.2\sigma$ (vertical line) every frame sits at or
below 2 per cent contamination. The shaded band is where all five collapse
together, between $1.5\sigma$ and $1.2\sigma$: the probability that noise clears
the threshold over the minimum connected area depends on the threshold alone,
not on how many stars the frame contains.

**(b)** The estimate validated against Gaia DR3. Detections matched to a Gaia
source within 2.0 arcsec are counted real and the remainder is the realized
spurious count (abscissa); the ordinate is the catalogue-free estimate from
panel (a). The 44 points, spanning four decades, follow the one-to-one line
inside a factor of two (shaded). In the collapse region where the number matters
most the agreement is 2.4 per cent (1309 estimated against 1341 realized for
M13 $R$ at $1.2\sigma$). Above the collapse the estimate runs *higher* than the
Gaia-realized count, which is the conservative direction for a gate and is partly
expected: real stars fainter than the Gaia catalogue limit are counted as
spurious on the abscissa.

## Data

- Frames: `pp_messier13-0004-R`, `pp_messier13-0001-B`, `pp_messier3-0001-B`,
  `pp_Messier67-0001-g`, `pp_NGC6811-0001-B` — Moravian C3-61000, 4788×3194,
  60 s each, from the raw→science reprocessing tree.
- Detector settings held fixed across the scan at the pipeline defaults
  (`minarea = 3`, `deblend_nthresh = 64`, `deblend_cont = 0.004`); only the
  threshold varies.
- Truth for panel (b): Gaia DR3 field catalogue from the frame's own step-5
  product, restricted to sources landing inside the image.

## Caveats stated in the text

- The quantity measured is **noise-origin** contamination. Cosmic rays, hot
  pixels and satellite trails are positive-only and have no negative
  counterpart, so they are invisible to this test and are handled by the
  separate cosmic-ray rejection stage.
- The symmetry assumption was checked rather than assumed: after background
  subtraction the pixel distribution has skewness $+0.058$ (M13 $R$) and
  $+0.051$ (NGC 6811 $B$), matching the Poisson expectation $+0.045$ and
  $+0.048$ for those sky levels, with a median offset of $0.02\sigma$.
- Five frames, one camera. The floors quoted are not a recommendation
  transferable to another instrument.
