# Figure 12 — Per-step preprocessing cross-check vs ccdproc

**Figure 12.** Each APEX detector-calibration stage compared, pixel-for-pixel,
against the equivalent operation in **astropy ccdproc** — the community-standard
Python CCD-reduction package — on the same real Moravian C3-61000 frames
(NGC 6811, B, 2x2; 8 bias, 8 darks, 5
flats, one science frame). **(a)** The maximum per-pixel disagreement for master
bias/dark/flat construction, bias/dark/flat application, and the full pipeline.
Master bias, master dark, and all three application steps are **bit-identical**
(delta = 0); master-flat construction and the full pipeline agree to
6e-08 and 5e-04
DN (float32 rounding). All stages sit four to nine orders of magnitude below the
detector read noise (3.45 DN) and the sky shot noise
(31 DN). **(b)** The same three quantities as a noise budget.
The two independently-written pipelines are numerically identical, so APEX's
reduction implements the standard bias/dark/flat arithmetic correctly. This is a
cross-implementation check (analogous to the sep and IRAF photometry
cross-checks), not a ground-truth validation — the latter is the synthetic
inject->recover test. Cosmetic correction uses astroscrappy, the L.A.Cosmic
reference implementation.
