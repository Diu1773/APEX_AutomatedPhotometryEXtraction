# Figure 9 — Crowded-field validation on a real globular cluster (M5)

**Figure 9.** Figs 1-8 validate APEX on open clusters and uncrowded synthetic
fields; this figure extends into a genuinely denser real field. M5 (NGC 5904)
was re-reduced end-to-end with the *current* codebase (Step-7 forced
photometry with the fixed sky annulus, Step-8 PSF photometry, Step-10
calibration with the colour-solve/quadratic/Gaia-quality fixes), 34
$r$-band frames, 1406 master sources.

**(a)** Radial source density around the field's true density peak (found
with APEX's own `psf_core` auto-density-center routine, not the field
centroid): M5's core is enhanced **34$\times$** over the field
background, compared to **11$\times$** for the open cluster
NGC 6811 used in Figs 1-8 — by this same metric, M5 is a genuinely denser
field. **(b)** Aperture-vs-PSF magnitude residual (APEX's own two
independent measurement methods, zeropoint-aligned on
219 bright/isolated stars) versus each star's nearest-neighbour
separation (`neighbor_dist_px`, an APEX master-catalog product), for
683 $r$-band stars matched positionally between the two methods.
The grey band marks the dataset's detection/deduplication floor
(~10 px, set by the master-catalog minimum separation).

**Honest result.** Two probes were run before drawing this figure: (1)
residual against the Gaia-transformed reference, binned by neighbour
separation at fixed magnitude, and (2) the aperture-vs-PSF comparison shown
here. Neither shows a clean crowding-driven degradation — panel (b)'s
binned medians are flat and consistent with zero (within $\pm$0.01-0.02 mag)
from the resolution floor out to isolated separations. Probe (1) is
confounded by Gaia's own RUWE/$C^*$ quality cuts, which preferentially
reject the worst blends before a star ever reaches the calibrator table
(a survivorship bias, not evidence of APEX robustness); probe (2), shown
here, is not subject to that bias since it compares two APEX-internal
methods on the same detected/matched star list. **Within the $\sim$10 px
($\sim$4$''$) separation this ground-based, 2$\times$2-binned dataset
resolves, APEX detects, forced-photometers, and PSF-fits M5's core
correctly, with no detected crowding-dependent bias between the two
methods.** This is a genuine positive finding about the domain of validity
established here — not a claim that crowding never degrades aperture
photometry. Sub-resolution blending (separations below this dataset's own
detection floor, as would be probed by space-based or lucky imaging) is not
and cannot be tested by this dataset, and remains open.
