# APEX A&A Section 15 submission-readiness audit

Target: **Astronomy & Astrophysics, Section 15 — Numerical methods and codes**. The official author guide lists Section 15 under the numerical-methods/code sections; the target is not the instrumentation section. See the [A&A author guide](https://www.aanda.org/doc_journal/instructions/aadoc.pdf).

## Readiness snapshot

| Area | Status | Evidence / blocker |
|---|---|---|
| Scientific scope | **Promising, not frozen** | APEX is best framed as an integrated, parameterised photometry workflow with scoped validation, not as a new universal detector/photometry algorithm. |
| Target fit | **Good pending reframing** | Section 15 is compatible with a code/method paper if the implementation, reproducibility contract and numerical validation are explicit. AutoPhOT is a useful precedent but is not evidence of APEX equivalence. |
| Code path | **Amber** | Steps 0–7 shared core is clear; PSF Step 8 and zeropoint Step 10 remain GUI-worker centric. |
| Reproducibility | **Amber/red** | Benchmarks use seeds in many places, but the normal pipeline manifest omits resolved config, package versions, input hashes, external solver/database versions and network/cache state. |
| Validation | **Amber** | Strong synthetic and cross-check infrastructure exists; coverage by boundary and negative controls is incomplete. Optional external suites can skip. |
| Documentation | **Amber** | Paper, README/manual and code disagree on WCS priority and headless scope. `docs/audit` now records the discrepancies; source docs still need a controlled update. |
| Release hygiene | **Amber** | Tracked documentation still contains local absolute paths. **The API-key blocker was a false positive — cleared 2026-08-12** (below). |
| Tests | **Green for regression baseline** | Full oracle run on 2026-08-11: **1,039 passed** in 710 s (2026-08-07 baseline was 905). This is software regression evidence, not a substitute for paper validation. |

## A&A Section 15 package checklist

### Must be completed

1. **Methods boundary:** one canonical diagram and table separating R/N-A/N-W/N-O/N-S/M components.
2. **Execution contract:** documented commands for the shared headless runner, GUI workflow, Step 8/10 off-screen scripts, and exact config resolution.
3. **Numerical validation:** retained input manifests and machine-readable results for artificial stars, SEP, IRAF/DAOPHOT, Python `ccdproc`, WCS engine comparison, PSF, time-series and CMD tests.
4. **Reproducibility manifest:** git commit, Python/package versions, camera/filter/frame counts, parameter snapshot, seed, solver engine/fallback, Gaia/cache status, worker count and output hashes.
5. **Scope limits:** detector type, crowding/FWHM, magnitude/SNR, sky structure, calibration assumptions, Gaia availability and external-runtime requirements.
6. **Terminology repair:** Python `ccdproc` versus IRAF/PyRAF task; internal WCS default versus optional ASTAP/astrometry.net; same-backend photutils consistency versus independent references.
7. **Native rationale:** describe only APEX-specific algorithms and special policies in the dedicated methods subsection; state the alternative package, unmet contract, validation and performance boundary for each.
8. **Bottleneck contract:** add wrapper parity tests and record `HAS_BOTTLENECK` plus package versions in benchmark manifests; do not present selected reduction acceleration as a pipeline-wide speed claim.
9. **Public release:** remove private paths and internal tracking/working material from the public release surface; add a clean example config, citation metadata, license, version tag and archival DOI plan.

### Do not claim yet

* universal instrument portability;
* fully Qt-free end-to-end processing;
* automatic, prior-free recovery of cluster physical parameters;
* linear speed-up or “optimized” performance without a baseline;
* independent validation from photutils remeasurement;
* equivalence of Python and IRAF `ccdproc` merely because a plot looks similar;
* successful optional suites when the validator reports `skipped`.

## Public/private repository decision

### The API-key finding was a false positive (verified 2026-08-12)

The earlier audit reported "several tracked HTML files match a Google-style
API-key pattern" and made that a release blocker. It is a regex collision.

Exactly one tracked file matches today, `validation/paper/STUDY_GUIDE.html`, and
the match sits inside a **204,332-character base64 run** that begins
`<img src="data:image/jpeg;base64,` — the embedded JPEG of Figure 8. The four
characters `AIza` occur by chance in the image bytes. A Google API key appears
as a standalone token (a `key=` parameter or a config value), never buried
mid-blob.

History was checked as well, not just the working tree. Six commits changed the
occurrence count (`git log --all -S"AIza"`), all of them paper/figure/UI commits
— which is what regenerating an embedded image looks like. Extracting every
matching blob at each of those commits and testing its context gives **4 matches,
all inside the same base64 image run, 0 outside base64.**

So there is no credential to rotate, and **the history rewrite this document
warned against must not be performed on these grounds.** If a scanner flags it
again, the test is whether the match sits inside a `data:` URI.

Local absolute paths in tracked documentation are a separate item and still
stand.

### Public/private sequence

The current repository is a useful private research workspace but is not yet a
clean public release. Recommended sequence:

1. ~~Rotate/revoke any credential-like value found in tracked HTML~~ →
   **2026-08-12: no credential exists.** Instead, add a scanner-suppression note
   (or strip embedded data URIs from the published copy of `STUDY_GUIDE.html`)
   so the same false positive does not block release again.
2. Keep the full research repository private, including tracking notes, local paths, raw-data references and scratch provenance.
3. Create a release branch/repository containing source, tests, reproducible benchmark fixtures, paper scripts, non-sensitive manifests and the seven audit documents (or a public summary of them).
4. If history cleaning is required, make a verified backup and obtain explicit approval before any filter-repo/force-push operation.

## Go/no-go gate

**No-go today** for a public A&A submission package, because claim scope is not
yet controlled. **Go after:** the seven audit actions are addressed, the
WCS/headless/manual contradictions are repaired, optional validation status is
explicit, and a clean release is rerun from a pinned environment.

The credential half of the release-hygiene blocker is cleared (2026-08-12, above);
local absolute paths remain. The largest remaining *claim* item is that the
abstract's "IRAF/DAOPHOT" cross-check is IRAF's **aperture** task `phot` — the
PSF engine's first external comparison was run on 2026-08-12 and is preliminary
(`validation/psf_engines/`). See the claim matrix.
