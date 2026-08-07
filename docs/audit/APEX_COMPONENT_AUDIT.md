# APEX component audit (source-of-truth map)

Audit date: 2026-08-07. This is a code-trace audit, not a claim that every branch has been run. Paths and function names refer to the current repository; the function name is the stable locator when line numbers move.

## Audited environment snapshot

The following versions were queried from `.venv-deploy` on the audit date. They are
the versions for which the current verification run (905 passed tests) is meaningful;
they are not a substitute for pinning the release environment in the paper.

| Package | Version | Role in this audit |
|---|---:|---|
| Python | 3.12.3 (repository supports 3.10+) | execution |
| NumPy | 2.4.4 | arrays, vectorised native operations, fallback reductions |
| SciPy | 1.17.1 | `cKDTree`, optimisation, signal/linear algebra |
| Astropy | 7.2.0 | FITS, coordinates/WCS, time-series LS/BLS |
| photutils | 2.3.0 | detection/segmentation, aperture and PSF primitives |
| SEP | 1.4.1 | independent source-extraction engine |
| Bottleneck | 1.6.0 | selected `fast_stats` reductions |
| astroscrappy | 1.3.0 | L.A.Cosmic cosmic-ray primitive |
| pandas | 2.3.3 | tables, SYSREM pivot and workflow data |
| ccdproc | 2.5.1 | validation-only Python CCD reduction cross-check |

## Execution spine

| Boundary | Actual path | Audit result |
|---|---|---|
| GUI | `apex/gui/workflow/step*_*.py` → Qt worker | UI owns selection, parameters, progress, and review decisions. |
| Shared headless | `apex/cli.py:_cmd_run` → `apex/pipeline/registry.py:get_steps` → `PipelineRunner.run` | Steps 1–7 are registered; detector calibration is a separate optional Step 0 (`get_calibration_step`). |
| Runtime context | `apex/pipeline/context.py:RunContext.build` | Reads mode-specific parameters and synchronises CLI path overrides into `params.P`. Context is Qt-free. |
| Step contract | `apex/pipeline/base.py:PipelineStep` / `StepResult` | Prerequisites, completion markers, status, outputs, and exceptions are represented. |
| Manifest | `apex/pipeline/runner.py:_write_manifest` | Records mode, start/end, and per-step status/duration/output; it does **not** snapshot resolved parameters or package versions. |

## Component map

| Component / role | Actual implementation and backend | Class | GUI/headless status | Safe interpretation / missing evidence |
|---|---|---|---|---|
| FITS scan and frame classification | `apex/analysis/calibration_scan.py:classify_type`, `read_frame_info`, `scan_folder`; `apex/analysis/calibration.py:_load_fits` | N-W + R (Astropy FITS) | Shared by Step 0 GUI and `CalibrationStep` | Header/filename heuristics, night/temperature/exposure grouping are workflow policy. Need a mislabeled-frame confusion matrix and multi-extension FITS tests. |
| Overscan | `apex/analysis/overscan.py:correct_overscan`, `correct_overscan_from_header`; called by `calibration.py:_apply_overscan` | N-S/N-W | Qt-free; Step 0 uses it | Row/column median subtraction and optional trim are conventional arithmetic. Disabled by default (`CalibrationOptions.overscan_enable=False`); no multi-amplifier path is claimed. Validate on real overscan geometry before advertising detector generality. |
| Master bias/dark/flat | `calibration.py:combine_frames`, `build_master_bias`, `build_master_dark`, `build_master_flat` | N-A/N-W | Qt-free, GUI and headless Step 0 share `calibration_run.run_calibration` | Median or per-pixel MAD sigma-clipped mean; dark exposure scaling and temperature matching are explicit policies. Independent Python `ccdproc` cross-check is in `apex/benchmark/calibration_crosscheck.py:run_crosscheck`; IRAF is a separate PyRAF tool, not the Python package. |
| Light calibration | `calibration.py:calibrate_light`, `calibrate_light_file`; formula is `(raw-bias)-k*dark`, divided by flat, then pedestal | N-S/N-W | Qt-free | This is standard CCD arithmetic with project-specific dark/pedestal safeguards. The evidence does not establish universal detector linearity or a calibrated uncertainty model. |
| Cosmetic correction | `apex/analysis/cosmetic.py:clean_frame`, `hot_pixel_mask`, `star_protect_mask` | M: astroscrappy R + native mask/interpolation N-A | Qt-free function; optional Step 0 integration | Cosmic-ray detection is delegated to `astroscrappy.detect_cosmics`; persistent hot pixels use a native master-dark threshold and 3×3 median interpolation. Do not call the whole stage “native L.A.Cosmic.” Validate recovery of injected cosmic rays and star-core preservation. |
| Crop/mask | `apex/analysis/crop.py:run_crop`, `_is_crop_cache_valid`; `pipeline/steps/crop.py` resolves GUI rectangle/config | N-W/N-O | Shared Qt-free path | Crop is a workflow/cache policy, not a new science algorithm. FITS metadata/WCS propagation and stale-cache behavior need explicit release tests. |
| Sky and frame-level QC | `apex/pipeline/steps/sky_qc.py`; `apex/analysis/frame_qc.py:evaluate_frame_qc`, `summarize_frame_qc` | N-W/N-S | QC core shared; UI presents review | Robust z-scores, FWHM/elongation/sky/excess-detection gates, and PASS/REVIEW/FAIL are project policy. Threshold sensitivity and false-rejection rates are not yet a complete validation. |
| Detection, deblend, centroid | `apex/analysis/detection.py:run_detection`; SEP branch (`sep.extract`), segmentation branch (`photutils.detect_sources/deblend_sources`), DAO refinement (`DAOStarFinder`) | M: SEP/photutils R + native orchestration/QC N-W | Headless core is called by `step4_source_detection.DetectionWorker`; callbacks are Qt-only | Engine and presets are configurable. `parameters_cmd.py` defaults to `segm`, while `parameters_lc.py` defaults to `dao`; a manuscript-wide “the detector” statement is unsafe without mode/parameter specification. |
| Detection cache | `apex/utils/cache_utils.py:build_detection_cache_signature`, `detection_cache_signature_matches`; writes in `detection.py` | N-W/N-O | Shared | Cache checks schema, engine, path/crop and size/mtime with a relaxed WCS-header rule. It is not content hashing; manual edits with unchanged size/mtime can survive. |
| Gaia query/cache and WCS refinement | `apex/utils/gaia_catalog_service.py:GaiaCatalogService`; used by `wcs_solve.py` and `refbuild.py` | R + N-W | Qt-free service, GUI workers call it | Network/cache availability, magnitude limits, radius, retries and no-cache behavior are configurable. Gaia absence can skip refinement/residual statistics; this is a declared fallback, not a guarantee of identical results. |
| WCS default solver | `apex/analysis/wcs_solve.py:resolve_wcs_engine`, `run_wcs_solve`; `apex/analysis/astrometry/quad_matcher.py` and `solver.py` | N-A + N-W | Qt-free `InternalWcsWorkerBase`; GUI wraps it with `QThread` | **Actual default is the internal Python quad/RANSAC/TAN-SIP path.** ASTAP and astrometry.net are explicit external engines; astrometry.net is attempted after ASTAP only when the ASTAP path and fallback option are selected. This corrects stale docs that describe ASTAP→astrometry.net→Internal priority. |
| External WCS engines | `wcs_solve.py:WcsWorkerBase` (ASTAP subprocess), `AstrometryNetWorkerBase` (solve-field/WSL subprocess) | R + N-W wrapper | GUI and headless dispatcher support explicit engine | Availability, database/index files, WSL, timeout, subprocess return code and artifact ingestion are environmental. Compare accepted WCS and residuals, not raw solver exit codes alone. |
| Master catalogue / IDs | `apex/analysis/refbuild.py:run_refbuild`, `build_master_catalog`, `match_detections_to_gaia`; `apex/analysis/merge/id_match.py` | M: Astropy matching R + native ID/QC policy N-W | Shared Step 6 core; GUI worker wraps it | WCS-QC, match radii, duplicate resolution, frame coverage and Gaia ID propagation are policy. Need an injected-coordinate ID recovery matrix across density, WCS error and missing-Gaia cases. |
| Forced aperture photometry | `apex/analysis/forced_photometry.py:run_forced_photometry`, `_phot_frame`; `apex.utils.photometry_utils:phot_vectorized` | M: photutils R + native orchestration/N-S noise model | Shared Qt-free Step 7 core and GUI wrapper | Circular aperture/annulus is delegated to photutils; recentering, registration, saturation/nonlinearity flags, depth QC and output policy are native. The growth-curve implementation reuses one annulus-stat pass (`_growth_curve_fixed_sky`) but this is an optimization, not accuracy evidence. |
| Aperture correction | `forced_photometry.py:_growth_curve_fixed_sky`, `_simple_apcorr` | N-W/N-S built on photutils | Shared | Automatic reference-star selection, cap, minimum SNR/count and fallback median ratio are project policy. It is not a new photutils algorithm; validate against injected total flux and crowded/structured sky. |
| PSF/ePSF photometry | Main fit is in `apex/gui/workflow/cmd/step8_psf_photometry.py:Step6PSFWorker`; uses `photutils.psf` (`EPSFBuilder`, `PSFPhotometry`/fitters) | M: photutils R + native policy/iteration N-W | **Not a pure `apex.analysis` core.** `scripts/run_step8_headless.py` imports PyQt5 and this GUI worker, so “headless” means off-screen Qt, not Qt-free. | ePSF star selection, fit window, residual re-detection, caps, quality flags and aperture flux scaling are native (`psf_policy.py`, `psf_iteration.py`, `psf_flux_scale.py`). A claim that PSF computation is fully GUI-independent is currently too strong. |
| Photometric QC/error labels | `apex/analysis/photometric_qc.py`, `utils/source_quality.py`, `frame_qc.py` | N-W/N-S | Mixed shared/UI | QC labels are decision policies, not calibrated probabilities. The reported error model needs pull/coverage tests under varied gain, read noise, sky structure and saturation. |
| Zeropoint/color calibration | `apex/gui/workflow/cmd/step10_zeropoint_calibration.py:ZeropointCalibrationWorker`; `apex/analysis/cmd/standard_anchor.py` | M: Astropy/catalog services R + native robust fit/policy | GUI worker; `run_step10_headless.py` imports the GUI worker with `QCoreApplication` | Standard-star discovery, anchor selection, color terms and Gaia/PS1 fallbacks are workflow decisions. Instrument-to-instrument generality is not established by one camera/night. |
| CMD/isochrone | `apex/analysis/cmd/isochrone_fit_service.py`, `isochrone_mcmc.py:fit_isochrone_mcmc`; GUI `step12_isochrone_model.py` and `scripts/run_step12_headless.py` | M: emcee/NumPy/SciPy R + native likelihood/prior model N-A/N-W | Step 12 has a genuine Qt-free service/runner | MCMC is stochastic but a seed is accepted and global NumPy RNG is explicitly seeded in `isochrone_mcmc.py` near the sampler. Degeneracy and prior sensitivity remain scientific limitations; “automatic physical recovery” is unsafe. |
| Time-series/detrending | `apex/analysis/light_curve/global_ensemble.py`, `sysrem.py:sysrem`, `detrend_output_service.py` | N-A/N-W + Astropy where used | GUI workers call analysis services; some LC orchestration remains GUI-centric | SYSREM/global ensemble are native implementations and need injection tests, comparison-star leakage checks, and multi-night missing-data tests before “robust detrending” is claimed. |
| Period analysis | `period_analysis_service.py:compute_ls`, `compute_pdm`, `compute_bls`, `bootstrap_fap`; `period_alias_service.py` | M: Astropy Lomb–Scargle/BLS R + native PDM/alias policy N-A/N-W | Service is Qt-free; GUI worker wraps it | PDM and alias diagnostics are real native code; bootstrap FAP has seeded permutations and thread parallelism. Validate false-alarm calibration under the actual cadence and red/systematic noise, not only sine injections. |

## Classification rule used

`R` means an external library/service is called for the scientific primitive; `N-A` is a project algorithm, `N-W` a project workflow/acceptance policy, `N-O` a computational optimization, `N-S` standard arithmetic, and `M` combines these. Imports alone were never treated as native work. A component is not “novel” merely because it is assembled in APEX.

## Native-vs-reused rationale (engineering evidence)

The word "native" is used narrowly here. It does not mean that APEX reimplemented a
library primitive for speed. In most stages APEX calls the established primitive and
implements only the surrounding contract: parameter presets, reference selection,
fallbacks, QC, provenance, or output schema. A native algorithm is retained only
when the required operation is absent from the selected in-process stack or when an
external executable would break the intended offline, per-frame workflow. The table
below records the reason that can be supported by source inspection. It does **not**
claim a universal speed advantage; no hardware-normalised APEX-vs-package benchmark
has been committed.

| Native or mixed component | Why a package primitive was not used as the whole solution | Where the performance decision is real | What is still reused | Safe manuscript wording |
|---|---|---|---|---|
| Internal WCS/quad solver (`astrometry/quad_matcher.py`, `solver.py`) | `astropy.wcs` transforms coordinates but is not a blind plate solver. ASTAP and astrometry.net are external executables with separate databases/indexes (and, for solve-field, WSL/network/environment requirements). The internal path supplies the same per-frame Gaia refinement, acceptance gates, callbacks and sidecar/header handling without those runtime requirements. | Quad construction and catalogue matching are vectorised; `cKDTree.query` batches all source quads instead of a Python per-quad query, and RANSAC is bounded by `max_trials`. These are implementation choices, not a measured claim that the solver is faster than ASTAP or astrometry.net. | Astropy coordinates/WCS, SciPy `cKDTree`/linear algebra, and Gaia service/cache. | “The internal solver is an offline in-process default; ASTAP and astrometry.net remain selectable external engines. The implementation uses batched vector operations and bounded RANSAC, but solver speed was not claimed without a benchmark.” |
| PDM (`period_analysis_service.py:_pdm_theta_vectorized`) | Astropy supplies Lomb–Scargle and BLS, but not the Stellingwerf PDM statistic used here. A subprocess/tool replacement would require a different input/output and diagnostic contract. | Trial periods are processed in NumPy batches with a memory bound (about 50 MB per batch); the requested grid is capped at 50,000 trials. The source comment reports a scalar-loop speedup estimate, but it is not a release benchmark. | NumPy, SciPy peak finding, Astropy Lomb–Scargle for the complementary method. | “PDM is a project implementation of the Stellingwerf statistic because the selected Astropy time-series API does not provide PDM; it is vectorised and memory-bounded, not asserted to be universally faster.” |
| SYSREM (`analysis/light_curve/sysrem.py`) | The core dependency set has no SYSREM primitive with this missing-data/error model and explicit exclusion of the target from component estimation. The implementation is based on Tamuz et al. (2005), not presented as a new algorithm. | The alternating coefficient updates are matrix/vector operations; the DataFrame adapter uses a pivot rather than `iterrows`. This reduces Python overhead but is not an accuracy or universal throughput claim. | NumPy/pandas for the numerical and table operations; the literature algorithm and citation define the method. | “SYSREM follows Tamuz et al. with an APEX-specific missing-data, weighting and target-preservation contract; the implementation is vectorised.” |
| Aperture-correction workflow (`_growth_curve_fixed_sky`, `_simple_apcorr`) | photutils provides aperture primitives and growth-curve ingredients, but not APEX's automatic reference-star selection, S/N and count cuts, correction cap, fallback median ratio, and application to every frame. Replacing the workflow with DAOPHOT/IRAF would add an external format/process boundary. | The fixed-sky growth curve computes annulus statistics once and reuses them for the radius grid; this is an I/O/arithmetic reduction, not a new aperture algorithm. | photutils aperture masks/statistics and Astropy units/FITS. | “APEX adds an automatic, QC-controlled aperture-correction workflow around photutils; it does not replace photutils aperture photometry.” |
| Calibration combination (`calibration.py`) | `ccdproc` is retained as an independent cross-check, while APEX controls frame grouping, exposure/temperature matching, pedestal and output metadata. Using a package call alone would not express the full APEX policy or manifest. | Stack dtype/chunking and bounded intermediate arrays address peak memory; the code does not establish that its median/sigma-clipped combine is faster than `ccdproc`. | Astropy FITS I/O; `ccdproc` only in the validation cross-check, not the production implementation. | “APEX implements the declared calibration policy and verifies it against Python ccdproc; no speed superiority is implied.” |
| Detection/QC and PSF iteration | SEP, photutils segmentation/deblending, DAOStarFinder, and photutils PSF fitting remain the scientific primitives. Native code chooses engines/presets, merges/refines detections, selects PSF stars, iterates residuals, and assigns QC flags. | `cKDTree` proximity filtering, vectorised image operations, label caps and configurable worker pools bound cost; there is no evidence that APEX's policy layer is a faster detector than the underlying engines. | SEP, photutils, SciPy and Astropy. | “APEX contributes orchestration and acceptance policy around established detection/PSF primitives.” |
| `fast_stats` / Bottleneck | This is not a native statistics algorithm. `apex.utils.fast_stats` is a compatibility wrapper: it calls Bottleneck for standard reductions when importable and falls back to NumPy with the same operation. `bottleneck` is declared as a core/required dependency in `pyproject.toml`, `requirements.txt`, the Windows spec and `apex/cli.py`, while the fallback keeps module imports and reduced environments usable. | Bottleneck can accelerate large NaN-aware reductions; a local spot-check showed a gain on calibration-sized arrays but no release benchmark establishes a fixed percentage. It is used by calibration, detection, cosmetic robust-stat helpers, and the Step 4 GUI helper. | NumPy is the scientific fallback and remains the direct implementation for many small or algorithm-specific reductions. | “Bottleneck is a declared core dependency used through a compatibility wrapper for selected reductions; NumPy fallback preserves semantics. It is not used throughout the pipeline and no universal speed claim is made.” |

### Bottleneck boundary (important for reproducibility)

The audited call graph imports `fast_stats` in `calibration.py`, `detection.py`,
`cosmetic.py`, and `gui/workflow/step4_source_detection.py`. The wrapper exposes
`nanmedian`, `nanmean`, `nanstd`, `nansum`, and `nanmax`, but only calls routed through
that wrapper can use Bottleneck. Many other modules still call `numpy.nanmedian`,
`numpy.nanmean`, `numpy.nanstd`, or `numpy.nansum` directly (for example forced
photometry, WCS residual summaries, light-curve services, overscan, photometry
helpers, and benchmark scripts). The source does not document a deliberate
“Bottleneck-off” scientific policy for each of those calls. Therefore the paper must
not say that Bottleneck accelerates APEX globally; if this distinction matters to a
performance result, record the exact module path, package versions, array shape and
benchmark condition.

## Immediate documentation corrections

1. State `ccdproc` as the Astropy-affiliated Python package and IRAF `ccdred/ccdproc` as a separate PyRAF/IRAF execution; they are not the same package.
2. State the internal WCS solver as the default. ASTAP and astrometry.net are optional, explicit engines; ASTAP→astrometry fallback is conditional.
3. Split “Qt-free headless” into shared Steps 0–7 core, genuine Step 12 service, and off-screen Qt wrappers for Steps 8/10 (and the legacy Step 7 script).
4. Treat photutils comparison as same-backend consistency, not an independent validation; SEP and IRAF/DAOPHOT are the stronger independent comparisons.
