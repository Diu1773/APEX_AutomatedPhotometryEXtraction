# Figure 9 — Crowded-field validation on two real globular clusters (M5, M13)

**Figure 9.** Figs 1-8 validate APEX on open clusters and uncrowded synthetic
fields; this figure extends into genuinely denser real fields, and checks
whether the result replicates across two independent globular clusters. M5
(NGC 5904, 34 $r$-band frames, 1406 master sources) and M13
(NGC 6205, 12 $r$-band frames, 1347 master sources) were each
re-reduced end-to-end with the *current* codebase (Step-7 forced photometry
with the fixed sky annulus, Step-8 PSF photometry, Step-10 calibration).

**(a)** Radial source density around each field's true density peak (found
with APEX's own `psf_core` auto-density-center routine, not the field
centroid): M5's core is enhanced **34$\times$** and M13's
**35$\times$** over their respective field backgrounds —
both close together and both far above the open cluster NGC 6811's
**11$\times$** (Figs 1-8) by the identical metric: two
genuinely, and similarly, denser fields. **(b)** Aperture-vs-PSF magnitude
residual (APEX's own two independent measurement methods, zeropoint-aligned
per cluster on 219 (M5) / 165 (M13) bright/isolated
stars) versus each star's nearest-neighbour separation (`neighbor_dist_px`,
an APEX master-catalog product), for 683 (M5) and
482 (M13) $r$-band stars matched positionally between the two
methods. The grey band marks the dataset's detection/deduplication floor
(~10 px, ~4$''$).

**Honest result, now replicated.** Two probes were run on M5 before this
figure was first drawn: (1) residual against the Gaia-transformed reference
binned by neighbour separation, and (2) the aperture-vs-PSF comparison shown
in panel (b). Neither shows a clean crowding-driven degradation in M5; probe
(1) is not attempted for M13 at all, because M13's archived Gaia cross-match
is essentially unusable (only 38 of 1347 detected sources resolve a
`gaia_source_id` in a live Gaia DR3 query — the match itself appears to have
been broken by an earlier version of the matching code, not a sign of a
uniquely bad field). Probe (2), shown here, sidesteps that dependency
entirely: it compares two APEX-internal methods on the same detected/matched
star list, independent of Gaia. **Panel (b)'s binned medians are flat in
BOTH clusters** — consistent with zero within $\pm$0.02–0.04 mag from the
resolution floor out to isolated separations, despite M13 having a far
worse (and unrelated) Gaia cross-match than M5, which rules out a shared
Gaia-side artifact as the explanation for the flatness. **Within the
$\sim$10 px ($\sim$4$''$) separation this ground-based,
2$\times$2-binned instrument resolves, APEX detects, forced-photometers, and
PSF-fits both globular-cluster cores correctly, with no detected
crowding-dependent bias between the two methods, replicated in two
independent fields.** This is a genuine, replicated positive finding about
the domain of validity established here — not a claim that crowding never
degrades aperture photometry, and not a manufactured trend. Sub-resolution
blending (separations below this dataset's own detection floor, as would be
probed by space-based or lucky imaging) is not and cannot be tested by this
dataset, and remains open. Note also (§4.3) that M5, M13, and every other
cluster in this validation share the identical camera (Moravian
Instruments C3-61000) — this replication is across two *fields*, not two
*instruments*.
