# Figure 10 — Detector calibration (Step 0) validation

**Figure 10.** Three independent lines of evidence that APEX's built-in detector
calibration (bias/dark/flat, Step 0) is correct. **(a,b)** A real M13 *V* frame
before and after APEX calibration (identical asinh stretch): the optical
vignette and fixed-pattern structure are removed. **(c)** Equivalence to the
reference pipeline: APEX vs the AstralImage/AIPPI engine (documented to match
PixInsight WBPP within tolerance) on 19 target/filter combinations
across several nights, filters (B, H, R, V, g, g', i, i', r, r') and exposures. Master
bias/dark frames agree to $\sim$0.3 DN (0.06 %) and the fully calibrated
science frames are **bit-identical** (max$|\Delta| = 0.000$ DN). **(d)**
Cosmetic correction (cosmic rays + hot pixels) via L.A.Cosmic
(astroscrappy; van Dokkum 2001): injected artefacts are removed
(100 % of cosmic-ray pixels,
94 % of hot pixels) while real stars are
preserved (star-core false-positive rate 0.000 %,
aperture-flux change 0.00 mmag).

**Verdict.** APEX reproduces the reference calibration to the bit on real data
across targets/filters/exposures, and its optional cosmic-ray/hot-pixel step
follows the standard L.A.Cosmic algorithm without harming photometry — so APEX
performs detector calibration end-to-end (raw$\to$science), not just aperture
photometry on externally-calibrated frames.
