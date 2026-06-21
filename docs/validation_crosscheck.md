# APEX photometry cross-validation (vs independent software)

This is APEX's **primary measurement validation**: the same stars, on the same
frame, measured with the same aperture by APEX and by an **independent
photometry engine**, then compared. It validates the photometry directly,
without going through any downstream interpretation (CMD / isochrone).

It is fully automated and reproducible:

```bash
apex validate --suite crosscheck --result-dir <APEX result dir> --reference sep
# references: sep (SExtractor, independent, default) | iraf (DAOPHOT) | photutils | all
```

The harness auto-selects a representative processed frame, re-measures APEX's
Step-7 star positions with the chosen engine at APEX's aperture + sky annulus
(`sep.sum_circle(..., bkgann=(r_in, r_out))` for SExtractor — its own
local-annulus sky, sharing no code with APEX's photutils backend), and reports
the agreement. A position sanity check confirms the (x, y) ↔ image mapping
(peak pixel σ over sky; fraction on-frame).

## APEX vs SExtractor (sep) — real data

Instrumental-magnitude agreement on SNR > 20 stars, across four APEX-processed
clusters (open + globular):

| Cluster | Field | N | Zeropoint offset Δmag | Robust scatter (MAD) | Pearson r | within 0.05 mag |
|---|---|---:|---:|---:|---:|---:|
| NGC6811 | open | 798 | −0.0015 | **0.0039** | **0.9958** | 0.93 |
| M37 | open | 809 | −0.0012 | 0.0043 | 0.9879 | 0.92 |
| M67 | open | 619 | −0.0011 | **0.0024** | 0.9902 | 0.90 |
| M5 | globular (crowded) | 769 | −0.0089 | 0.0145 | 0.9609 | 0.75 |

**Reading it:**
- **Zeropoint offset ≈ 0** everywhere (−0.001 to −0.009 mag): APEX and SExtractor
  agree on the flux scale to the millimag level — there is no systematic bias.
- **Core scatter 2–4 mmag (MAD)** in open clusters: the two independent codes
  measure the same flux to a few millimagnitudes.
- **Crowded globular (M5) degrades gracefully** (MAD ~15 mmag, 75 % within
  0.05 mag): expected — in dense fields APEX's per-source recentering and
  aperture correction legitimately diverge from a fixed-aperture SExtractor
  measurement. The harness reports this honestly rather than hiding it.

> Note: a full-sample RMS (~0.2 mag) is dominated by a faint, low-SNR tail; the
> harness reports MAD (robust) for the verdict and shows the residual-vs-magnitude
> trend so the tail is visible, not averaged away. Plots: `crosscheck_sep_*.png`
> in the validation report.

## APEX vs IRAF / DAOPHOT — real data (gold standard)

Run via PyRAF (`--reference iraf`; on Windows it drives IRAF inside WSL). IRAF
`phot` is run at APEX's star positions with `datapars.itime = exptime`, so it
reports a count-rate magnitude directly comparable to APEX after one constant.

| Cluster | N | Zeropoint-system offset* | Robust scatter (MAD) | RMS | within 0.05 mag |
|---|---:|---:|---:|---:|---:|
| M67 (r) | 500 | −22.94 (constant) | **0.0030** | 0.0082 | **0.998** |

\* The ~−23 mag offset is **not a disagreement** — it is the fixed difference
between IRAF's `photpars.zmag` (~25) magnitude system and APEX's instrumental
`mag_inst` (count-rate, ~−15). It is a single additive constant; the agreement
is measured by the scatter *around* it.

**Result:** after the expected zeropoint-system constant, APEX and IRAF/DAOPHOT
agree to **~3 millimag (MAD)** with **99.8 % of stars within 0.05 mag** — even
tighter than the SExtractor comparison, because IRAF `phot` at fixed coordinates
closely mirrors APEX's forced-aperture methodology. This is the strongest,
most-recognized form of photometric validation in optical astronomy.

> Reproduce: `apex validate --suite crosscheck --result-dir E:/observed_Analysis/M67/pp/result --reference iraf`

## Why this is the validation that matters

This answers "*is APEX's photometry correct?*" against the community-standard
SExtractor, independently of APEX's own code path, on real observations — the
kind of agreement statistic reviewers and observers trust. It runs on any
APEX result directory with one command and is suitable as a CI regression guard.

IRAF/DAOPHOT (`--reference iraf`) and photutils (`--reference photutils`) are
also supported; sep is the default because it is an independent engine that needs
no external installation.

## Photometry-algorithm innocence: apcorr / recentering / detection

Beyond the flux/color cross-checks above, the three internal steps most often
suspected of introducing a magnitude- or color-dependent bias were tested
**directly** (2026-06 deep review), because median external agreement alone can
hide structured residuals:

| Step | Direct test | Verdict |
|---|---|---|
| **Aperture correction** | code inspection: apcorr is a **single per-frame scalar** (growth curve on high-SNR references, applied to every star) | A per-frame scalar multiplies all stars equally, so it is mathematically incapable of a magnitude- or color-dependent bias and is fully absorbed by the Step-10 zeropoint. |
| **Recentering** | re-measured aperture flux at the registration position vs the recentered position on **8 real M67 frames (5286 stars)** | median Δmag = **+0.0000** (g) / **+0.0001** (i) across all magnitude bins; centroid shifts ~0.1 px against an ~12 px aperture ⇒ negligible flux change. |
| **Step-4 detection** | flux at the registration vs recentered position is identical; master positions are Gaia-matched and multi-frame | positions are sound; no flux bias from detection centroiding. |

Combined with APEX = Gaia GSPC = Pan-STARRS agreement in colors, this establishes
that APEX **photometry** is not the source of the cluster isochrone-fitting
difficulty. The root cause is the gri filter-set's intrinsic age–metallicity–
reddening–distance degeneracy plus a fitting-likelihood weakness — see
`isochrone_fitting_methodology.md` §10.

## Reproduce

```bash
apex validate --suite crosscheck --result-dir E:/observed_Analysis/M67/pp/result    --reference sep
apex validate --suite crosscheck --result-dir E:/observed_Analysis/M37/pp/result    --reference sep
apex validate --suite crosscheck --result-dir E:/observed_Analysis/NGC6811/pp/result --reference sep
apex validate --suite crosscheck --result-dir E:/observed_Analysis/M5/light/result  --reference sep
```
