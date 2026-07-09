# Figure 13 — Cross-instrument + cross-pipeline validation (LCO)

**Figure 13.** APEX reduces public raw frames from two different Las Cumbres
Observatory cameras and is compared pixel-for-pixel against the archive's
BANZAI-processed product — an independent, published pipeline — testing the
calibration on foreign detectors. **Top:** a QHY600 CMOS camera (single
amplifier, 0.4 m; Proxima Cen field). The whole frame agrees to a uniform
+0.063 e⁻ offset (robust σ 0.080 e⁻); the
difference image is featureless. **Bottom:** a Sinistro CCD (four amplifiers,
1 m; NGC 5985 field). The sky and sources agree to ≈0.3 %, but the difference
shows a four-quadrant pattern (Δmedian +1.95 e⁻,
σ 0.74 e⁻): the per-amplifier assembly — gain, overscan and
cross-talk — that BANZAI performs with dedicated Sinistro handling and a generic
reduction does not. The bias/dark/flat calibration arithmetic generalises across
cameras; multi-amplifier detector assembly is instrument-specific (APEX targets
single-CCD detectors). ZScale stretch; raw data from archive.lco.global.
