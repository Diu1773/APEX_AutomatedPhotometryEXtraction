# APEX manuscript claim matrix

Source checked: `validation/paper/논문작업/MANUSCRIPT_ko.md` and nearby paper/README documents. The matrix is intentionally conservative: a number can be correct for one retained run while the surrounding generalisation is still unsafe.

| Manuscript claim (location) | Code/evidence checked | Status | Safer wording |
|---|---|---|---|
| APEX integrates raw frames to CMD or multi-night light curves without scripting (§1/abstract) | GUI workflows, shared Steps 1–7, CMD/LC modules | **Partial** | “Provides a GUI workflow for CMD and LC analysis; Steps 0–10 all run from a Qt-free runner (2026-08-17), while Steps 11–12 remain GUI-only.” |
| GUI and headless run the same computation (§1–2) | `analysis/detection.py`, `forced_photometry.py`, `wcs_solve.py`; `analysis/cmd/psf_photometry_runner.py`, `zeropoint_runner.py` | **Supported for Steps 0–10** (2026-08-17) | Say it as identity, not similarity: the window subclasses the same runner, so `Step6PSFWorker.run is PsfPhotometryRunner.run`. Measured on M13 with PyQt5 absent — Steps 8 and 10 both reproduced all 48 stored tables to 0.0e+00. Requires equal config, equal input **and equal package versions**: a first attempt differed on 5 of 22,305 crowded measurements by ≤5e-5 mag purely because the fresh install pulled scipy 1.18.0 instead of 1.17.1, and pinning it back made the difference vanish. That is why the run manifest now records the environment. See `docs/audit/HEADLESS_WITHOUT_QT.md`. |
| APEX does not propose new algorithms (§1/§2) | Native quad solver, PDM, SYSREM, aperture-correction policy, PSF iteration | **Defensible with qualification** | “Uses established primitives and contributes an integrated workflow plus project-specific policies/implementations; algorithmic novelty is not claimed.” |
| Four native components are quad solver, aperture correction, PDM, SYSREM (§2) | `astrometry/*`, `forced_photometry.py`, `period_analysis_service.py`, `sysrem.py` | **Incomplete** | Also describe native acceptance/QC, ID propagation, PSF iteration and calibration workflow as policies, not silently as primitives. |
| Built-in WCS is the default; ASTAP/astrometry.net are optional (§2/§3.6) | `wcs_solve.py:resolve_wcs_engine` and `run_wcs_solve` | **Supported in current path** | Explicitly state internal is default; ASTAP/astrometry.net are user-selected, and ASTAP→astrometry fallback is conditional. |
| WCS accepted with “20 matches, RMS <2.5 px, P99 <5 px” (§3.6) | Internal worker QC and configuration; paper text | **Needs exact code reconciliation** | Quote the thresholds actually used by the selected engine/config and report per-frame distributions; do not present hard thresholds as universal accuracy. |
| Detector has SEP, segmentation and DAO engines (§3.1/§3.2) | `detection.py` branches; CMD default `segm`, LC default `dao` | **Supported but mode-dependent** | Report engine, preset and sigma per experiment. Avoid “the detector” when results came from only one engine. |
| Completeness m50 ±0.08 mag (§abstract/§3.5) | `apex/benchmark/artificial_stars.py`, completeness fit and paper figure scripts | **Potentially supported for the named synthetic setup** | Give seed, PSF/crowding/sky, placement exclusions, number injected, fit model and confidence interval; do not generalise to all cameras/fields. |
| Error pull standard deviation 1.014 (§abstract/§3.7) | Benchmark metrics and forced-photometry error model | **Single-run evidence** | “For the reported benchmark configuration, pull scatter was …”; add residual-vs-SNR and coverage before calling the error model calibrated. |
| SEP agrees at 0.006 mag; IRAF/DAOPHOT at 0.0097 mag (§abstract) | `photometry_crosscheck.py`, `iraf_crosscheck.py` | **Conditional** | State robust scatter after frame zeropoint alignment, matched coordinates, selection cuts, frame/camera, and whether IRAF suite ran or skipped. |
| Python `ccdproc` and IRAF `ccdproc` are both comparison packages (§3.2/figure 4) | `calibration_crosscheck.py` imports Astropy `ccdproc`; IRAF path is PyRAF/IRAF `ccdred` | **Terminology correction required** | “Astropy-affiliated Python `ccdproc`” versus “IRAF/PyRAF `ccdproc` task”; they are separate implementations. |
| PSF and aperture agree across three cameras (§abstract/§3.10–3.11) | `fig_psf_validation.py`, PSF tests, `psf_policy.py` | **Evidence exists, scope limited** | Report cameras, frames, source-selection/crowding cuts, and that aperture agreement is a consistency check, not truth. |
| No density-dependent bias in globular clusters (§abstract) | Paper figures and retained cluster products | **Needs boundary language** | “No bias was detected within the tested density/FWHM/magnitude range”; quantify limits and exclusions. |
| APEX performs detector calibration from raw to science (§abstract) | Optional Step 0 and `calibration_run.py` | **Supported as a feature** | State Step 0 is optional/off-chain and list accepted frame types/overscan limitations. |
| APEX restores isochrone physical parameters (§§3.12–3.13) | `isochrone_mcmc.py`, synthetic validator, real-cluster scripts | **Overstated without priors** | “Fits a generative isochrone mixture under stated priors; synthetic recovery and real-data consistency were tested.” |
| “Automatic”/“reproducible” GUI processing | GUI state, run manifest, caches | **Partial** | Use “parameterised and resumable”; reproducibility requires captured config, package versions, external solver/catalog state and seed. |
| All science code is GUI-independent (§2) | `apex/pipeline/**` imports; fresh base install without PyQt5 | **Supported** (2026-08-17) | The Step 8/10 workers were lifted into `apex.analysis`; importing the whole pipeline registry now loads no `PyQt5` or `apex.gui` module. Two claims to keep separate: photometry Steps 0–10 need no Qt, while Steps 11–12 (isochrone, viewer) are still GUI-only. |
| “IRAF/DAOPHOT” as the independent photometry engine (abstract, §1, §3.11, Table 2 caption, §5, acknowledgements) | `apex/benchmark/iraf_crosscheck.py` calls `iraf.phot` only; it loads the daophot *package* but never `psf`/`allstar` (verified 2026-08-12) | **Ambiguous where it matters most** | The §3.11 body is already precise — “IRAF's `phot` (DAOPHOT, via PyRAF)”. The abstract and §1 compress this to “IRAF/DAOPHOT”, which a referee reads as PSF-fitting photometry — and the paper has a separate PSF section, so the two get connected. Say “IRAF's aperture task `phot`” wherever the word DAOPHOT stands alone. **The PSF engine had no external comparison until 2026-08-12** (`validation/psf_engines/`, ALLSTAR, preliminary). |
| “No crowding-dependent bias between aperture and PSF photometry” in two globulars (abstract) | `validation/psf_archive/README.md` (PSF−aperture scatter 0.0241 → 0.1421 mag from ≥6 FWHM to <1.5 FWHM); `validation/crowding_aperture_vs_psf.py` (2026-08-12, M13: aperture crowded/isolated zero-point residual 3.30× in V, PSF 2.05×) | **Needs the bias/scatter distinction stated, or it reads as refuted** | The claim is defensible only about the *median offset*. The aperture−PSF **scatter** grows ~6× into the core, and against an independent Gaia reference the aperture residual is 3.3× worse in the crowded quartile while PSF is 2.05×. Say “no crowding-dependent **offset** was detected; the scatter does grow with crowding, by X”, and quote both numbers. As written, today's own measurement is the counter-example a referee would find. |

## Native implementation and performance wording

The methods section should answer two separate questions: (1) why an established
package was not used as the complete solution, and (2) which local optimisations are
actually present. The defensible summary is:

> APEX reuses established primitives whenever they match the operation: SEP or
> photutils for detection, photutils for aperture/ePSF fitting, Astropy for WCS
> transforms and Lomb--Scargle/BLS, and SciPy for spatial searches. Native code is
> limited to missing operations (the Stellingwerf PDM statistic, the Tamuz SYSREM
> iteration, and an in-process blind WCS path) or to workflow contracts that the
> primitives do not define (automatic aperture-correction selection, QC, iteration,
> provenance and fallbacks). The implementation uses batched NumPy/cKDTree work,
> bounded trial grids, chunked calibration stacks and fixed-sky reuse. These are
> engineering choices; no package-wide speed superiority is claimed without a
> controlled benchmark.

The Bottleneck statement must be equally narrow: `apex.utils.fast_stats` routes
selected NaN reductions to Bottleneck when importable and otherwise to NumPy. It is
used by calibration, detection, cosmetic robust-stat helpers and the Step 4 helper;
many other modules call NumPy reductions directly. Thus the manuscript should say
“Bottleneck accelerates selected wrapper-routed reductions” rather than “Bottleneck
accelerates the pipeline” or “NumPy was too slow”.

The whole-code review found 555 direct `numpy.nan*` calls in 62 Python files. This
does not make every call a performance problem: independent benchmark calculations
should not share the acceleration path, and many calls are on small per-source or
per-frame-summary arrays. It does mean that the paper should identify Bottleneck as
a selected reduction backend and that a release benchmark must record
`HAS_BOTTLENECK`, package versions, array shape and worker/BLAS settings.

| Avoid | Use instead |
|---|---|
| “The native solver is faster than ASTAP/astrometry.net.” | “The internal solver removes executable/index/WSL dependencies and uses batched matching with bounded RANSAC; speed was not ranked without a benchmark.” |
| “NumPy was too slow, so PDM/SYSREM was rewritten.” | “PDM was absent from the selected Astropy API; SYSREM required the stated weighting, missing-data and target-exclusion contract. Their inner operations are vectorised.” |
| “APEX uses Bottleneck for its statistics.” | “Selected reductions are routed through an optional Bottleneck/NumPy compatibility wrapper; direct NumPy reductions remain elsewhere.” |
| “APEX replaces photutils.” | “APEX retains photutils primitives and adds selection, QC, iteration, correction and output policy.” |

## Required manuscript edits before submission

* Add version-pinned software references for Bottleneck/NumPy/SciPy (and the exact photutils/ccdproc releases) rather than citing an unnamed “Python calculation”.
* Add a one-row-per-experiment table with camera, filter, engine, config, frame/star count, seed, reference and metric.
* Replace “identical” with “same named function path under the same resolved configuration” where parity was actually checked.
* Separate software-library validation from APEX workflow validation. A package citation is not a result-level validation.
* State skipped optional suites and environmental requirements (PyRAF/WSL, ASTAP databases, astrometry.net indices, Gaia network/cache).
* Remove any “first”, “novel”, “universal”, “automatic recovery”, or “optimized” wording unless backed by a scoped comparison or benchmark.
