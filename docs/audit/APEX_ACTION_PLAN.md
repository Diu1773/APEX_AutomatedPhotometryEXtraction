# APEX audit action plan

This plan is deliberately ordered before manuscript layout work. It separates evidence cleanup from optional architecture improvements.

## P0 — block public release and paper claim freeze

| Action | Files/commands | Acceptance criterion |
|---|---|---|
| ~~Rotate/revoke credential-like values~~ **Cleared 2026-08-12 — false positive.** The only match sits inside the base64 JPEG of Fig. 8 in `validation/paper/STUDY_GUIDE.html`; history has 4 matches, all inside that image run, 0 outside base64. | see `APEX_SUBMISSION_READINESS.md` | Replaced by: strip embedded `data:` URIs from the published copy, or add a scanner-suppression note, so the same collision does not block release again. **Do not rewrite history on these grounds.** |
| **Fix the "IRAF/DAOPHOT" wording** (new 2026-08-12) | abstract, §1, §3.11, Table 2 caption, §5, acknowledgements | The IRAF cross-check is IRAF's *aperture* task `phot`; the harness never calls `psf`/`allstar`. Wherever "DAOPHOT" stands alone a referee reads PSF-fitting photometry — and the paper has a separate PSF section. Say "IRAF's aperture task `phot`". |
| Decide private research repo vs clean release repo | repository owner decision | A written boundary exists; no history rewrite or force-push is performed without approval. |
| Correct WCS terminology and default | `apex/analysis/wcs_solve.py`, `parameters.example.toml`, `docs/manual/02-shared-steps.md`, paper §2/§3.6 | All surfaces agree: internal default; ASTAP/astrometry.net explicit; conditional fallback only. |
| Correct ccdproc terminology | `validation/paper/IRAF_PREPROC_CROSSCHECK.md`, paper §3.2/figure 4 caption, references | Python Astropy `ccdproc` and IRAF/PyRAF task are identified as separate implementations with versions/runtime. |
| Mark optional suites honestly | `apex/benchmark/validate.py`, validation report template, paper methods | `skipped`, `failed`, and `passed` have distinct status; no skipped external engine contributes to a headline metric. |

## P1 — make the current evidence paper-grade

| Action | Evidence | Acceptance criterion |
|---|---|---|
| Build experiment manifest | `validation/paper/run_all.py`, `fig*.py`, retained `data/*.json` | Every headline figure has input IDs, camera/filter, engine, parameters, seed, frame/star count and software versions. |
| Add boundary matrices | `apex/benchmark/runner.py`, `artificial_stars.py`, `photometry_crosscheck.py`, WCS validators | Completeness, error pull, WCS residual, cross-check scatter and PSF agreement are reported by magnitude/SNR/crowding/FWHM/engine where applicable. |
| Add negative controls | tests/validation for detection, WCS, apcorr, PSF iteration, SYSREM, PDM/FAP | Known bad frames/coordinates/sky or shuffled labels are rejected at a documented rate. |
| Separate independent and same-backend checks | paper figures/table and `docs/validation_crosscheck.md` | SEP/IRAF are labelled independent; photutils is labelled same-backend consistency; truth injections remain the absolute reference. |
| Reconcile config defaults | `parameters_cmd.py`, `parameters_lc.py`, `schema.py`, `parameters.example.toml`, tests | One source of truth per mode; tests assert the documented default detector/WCS/worker settings. |
| Review Bottleneck boundary | `docs/audit/APEX_BOTTLENECK_REVIEW.md`, `apex/utils/fast_stats.py`, `tests/` | Add NumPy/Bottleneck parity tests; record `HAS_BOTTLENECK` and versions; benchmark only large-array candidates before routing more reductions through the wrapper. |
| Profile LC loading before adding more Bottleneck calls | `docs/audit/APEX_LC_PERFORMANCE_REVIEW.md`, Step 8--11 and LC tools | Measure cold/warm load, source-map reads, PSF reads, cache hit, build/QC/plot time and RSS; then fix shared map/cache/index paths before changing reductions. |

## P2 — strengthen reproducibility and architecture

| Action | Acceptance criterion |
|---|---|
| Enrich `pipeline_run.json` | Includes resolved parameter snapshot, package versions, input signatures, seed, worker counts, WCS engine/fallback, Gaia/cache status and output hashes. |
| Extract PSF service | `Step6PSFWorker` becomes a thin Qt adapter around a Qt-free request/result service; `run_step8_headless.py` no longer imports PyQt5. |
| Extract zeropoint service | Step 10 GUI/CLI share a Qt-free orchestration path; headless script is genuinely Qt-free. |
| Add LC pipeline runner | A single reproducible command covers the documented LC steps, or the paper narrows the end-to-end claim to the services that are actually runnable. |
| Add performance baseline | Fixed fixtures, workers 1/2/4, wall time, peak RSS, cache hit/miss, output count and variability are committed as machine-readable results. |

## P3 — manuscript and release polish

* Rewrite the abstract and contribution list from the claim matrix, keeping only scoped evidence.
* Add an implementation table with exact backend/package/function roles and citations.
* Add a hardware/software table, including Python, Astropy, photutils, SEP, ccdproc, PyRAF/IRAF, ASTAP/astrometry.net status and camera metadata.
* Keep cover/TOC/figures/layout work after the scientific and release contracts are frozen; otherwise figures will be regenerated against moving claims.
* Run the full oracle and the paper figure harness in a pinned environment; archive logs without local paths or credentials.

## Immediate next three tasks

Revised 2026-08-12. The previous item 1 asked the user to rotate a credential
that does not exist.

1. **Finish the PSF A-axis comparison** — the harness runs
   (`validation/psf_engines/`, IRAF ALLSTAR) but rests on 25 implanted stars.
   Re-score at ~400 so each magnitude × crowding cell carries statistics, then
   settle the remaining `varorder` confounder. Until this lands, the paper has
   **no external comparison for its PSF engine at all**, while its abstract says
   "IRAF/DAOPHOT".
2. Apply the WCS / `ccdproc` / headless / **`phot`-not-DAOPHOT** wording
   corrections to the manual and manuscript, using
   `APEX_MANUSCRIPT_CLAIM_MATRIX.md` as the gate.
3. User approves the private-vs-public repository boundary (the credential half
   of that decision is now moot).
