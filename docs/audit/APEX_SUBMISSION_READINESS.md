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
| Release hygiene | **Red until cleaned** | Tracked documentation contains local absolute paths, and several tracked HTML files match a Google-style API-key pattern. Do not publish or rewrite history until the credential is revoked/rotated and the public/private boundary is decided. |
| Tests | **Green for regression baseline** | Full oracle run on 2026-08-07: **905 passed, 24 warnings** in 362.66 s. This is software regression evidence, not a substitute for paper validation. |

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

The current repository is a useful private research workspace but is not yet a clean public release. Recommended sequence:

1. Rotate/revoke any credential-like value found in tracked HTML before opening or publishing anything.
2. Keep the full research repository private, including tracking notes, local paths, raw-data references and scratch provenance.
3. Create a release branch/repository containing source, tests, reproducible benchmark fixtures, paper scripts, non-sensitive manifests and the seven audit documents (or a public summary of them).
4. If history cleaning is required, make a verified backup and obtain explicit approval before any filter-repo/force-push operation.

## Go/no-go gate

**No-go today** for a public A&A submission package, because claim scope and release hygiene are not yet controlled. **Go after:** the seven audit actions are addressed, the WCS/headless/manual contradictions are repaired, optional validation status is explicit, and a clean release is rerun from a pinned environment.
