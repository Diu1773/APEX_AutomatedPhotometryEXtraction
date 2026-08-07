# APEX performance and parallelism audit

## Evidence boundary

The repository contains timing hooks and benchmark runners, but this audit did not find a committed, hardware-normalised baseline that supports a universal speed-up percentage. Therefore “parallel”, “vectorized”, or “optimized” below describes code structure only; it is not a measured performance claim.

## Why some code is native and some is not

The performance rationale is local, not a blanket claim that NumPy or an existing
package is too slow. The production path deliberately keeps established compiled
primitives where they match the scientific operation (SEP extraction, photutils
segmentation/apertures/PSF fitting, Astropy WCS transforms, Lomb--Scargle/BLS,
SciPy spatial searches). Native code surrounds those calls with workflow policy, or
implements a missing operation (PDM, SYSREM, internal blind WCS). The following
optimisations are visible in source and should be described as bounded engineering
choices rather than as benchmark results:

* WCS quad matching uses a batched `scipy.spatial.cKDTree.query` and bounded RANSAC
  trials instead of a Python query for every quad.
* PDM evaluates trial periods in NumPy batches and caps the grid at 50,000 trials;
  the batch size is reduced to bound temporary memory.
* SYSREM uses matrix/vector updates and a pandas pivot rather than row-wise table
  iteration, while preserving the explicit target-exclusion contract.
* Calibration controls stack dtype/chunking to limit peak memory; forced-photometry
  growth curves reuse one fixed-sky annulus calculation across radii.
* Frame-level concurrency is implemented with thread pools because the underlying
  Astropy/photutils/SciPy kernels may release the GIL, but scaling is not assumed.

These choices answer different bottlenecks: Python-loop overhead, temporary-array
memory, repeated I/O/statistics, or external-process latency. They should not be
collapsed into the sentence “NumPy was too slow”.

## Bottleneck is selective, not global

`apex/utils/fast_stats.py` is a small compatibility layer, not a second statistical
engine. For `nanmedian`, `nanmean`, `nanstd`, `nansum`, and `nanmax`, it calls
Bottleneck when available and otherwise calls the corresponding NumPy function.
The audited imports are `analysis/calibration.py`, `analysis/detection.py`,
`analysis/cosmetic.py`, and the Step 4 source-detection GUI helper. The package is
declared in the core dependency set, but the defensive fallback allows tests and
partial installations to retain NumPy semantics.

Many other paths call NumPy reductions directly: forced photometry, WCS residual
summaries, overscan, photometry helpers, light-curve services, and benchmark
scripts. The source does not give a module-by-module reason for bypassing the
wrapper. Some calls operate on small per-source arrays or are embedded in a
specialised formula, but that is an inference, not an audited performance result.
Consequently, report Bottleneck only for the wrapper-routed operations and include
the package version and benchmark condition if a speed number is added. Do not write
that “APEX uses Bottleneck throughout” or that it makes the pipeline faster without
a controlled comparison.

## Whole-code review result (2026-08-07)

The repository-wide search found `bottleneck` itself in the dependency/doctor
surfaces and `fast_stats.py`; only four consumers import the wrapper:
`analysis/calibration.py`, `analysis/detection.py`, `analysis/cosmetic.py`, and the
Step 4 source-detection helper. A separate search found 555 direct
`numpy.nan*` calls in 62 Python files (170 in `analysis`, 302 in `gui`, 49 in
`benchmark`, 32 in `utils`, and 2 in `io`; the wrapper's fallback calls are included
in these counts). This is not automatically a defect: benchmark/reference code
should remain independent, and many GUI or per-source calls operate on small arrays.
The defect is documentation/test coverage: the source does not state a module-level
policy for when a direct NumPy reduction is intentional.

### Findings and release actions

| Priority | Finding | Action before a performance claim |
|---|---|---|
| R1 | No test currently targets `fast_stats` parity with NumPy/Bottleneck, including all-NaN, empty, float32/float64, axis and `ddof` cases. | Add a small parity test; keep the NumPy fallback and compare values/NaN masks, not warning text. |
| R1 | `pyproject.toml`, `requirements.txt`, `deploy/apex_windows.spec`, and `apex/cli.py` treat Bottleneck as a core/required dependency, while `fast_stats.py` silently falls back when it cannot import it. | Keep this as an intentional compatibility fallback and state “declared core dependency; NumPy fallback” in the software/reproducibility section, or change the packaging policy consistently. Do not call it an optional install extra while the doctor marks it required. |
| R1 | A run manifest does not record whether the Bottleneck path was active. | Record `HAS_BOTTLENECK`, package version, NumPy version and worker/BLAS settings in the benchmark manifest. |
| R2 | High-volume direct NumPy reductions remain in `overscan.py`, light-curve matrix/ensemble paths, WCS summaries and other core modules. | Benchmark only candidate large arrays (calibration stack, overscan region, ensemble matrix, WCS catalogue); route a function through `fast_stats` only when the measured wall/RSS gain outweighs conversion and review cost. |
| R2 | Direct reductions in benchmark scripts are mixed with production reductions. | Keep independent baselines on NumPy; label them “reference calculation”, never “APEX Bottleneck result”. |

### Local spot-check (not a manuscript result)

On the audit workstation (Python 3.12.3, NumPy 2.4.4, Bottleneck 1.6.0), a single
microbenchmark with float32 arrays and `axis=0` gave Bottleneck/NumPy time ratios
of roughly 3.5–4.2 for `nanmedian`, `nanmean` and `nanstd` on a 20×2048×2048
calibration-like stack; `nansum` was about 2.2×. For a 14×500 growth-curve-like
array the absolute times were below a millisecond despite ratios up to about 5×,
so the result is unlikely to dominate a frame's photometry runtime. These numbers
are machine- and shape-specific and are retained only to motivate the candidate
benchmark above, not as a paper speed claim.

## Cost map

| Stage | Dominant work and current implementation | Parallel unit / default | Main risk | Required measurement |
|---|---|---|---|---|
| Step 0 scan/calibration | FITS reads are `memmap=False`; raw arrays are converted to float64 in `calibration.py:_load_fits`; master combine uses float32 stack but float64 output and may hold stack/mask/clipped arrays | Calibration run is sequential by night/group in `calibration_run.py` | Peak RAM scales with frame count × image size; repeated reads for master and science | peak RSS and wall time versus frames, image size, combine method; report dtype and I/O throughput |
| Overscan/crop | `overscan.py` makes copies; `crop.py:run_crop` can use `ThreadPoolExecutor` per file | `get_parallel_workers(params)`; thread pool | Large FITS copies and simultaneous decompression can saturate memory/disk | single-frame and N-frame scaling, with/without crop, peak RSS |
| Detection | `detection.py:run_detection` loads each frame, builds background/segmentation/deblend/DAO products, writes JSON/CSV caches; SEP branch casts contiguous float32, other branches use float arrays; `cKDTree` is used for proximity filtering | ThreadPool per frame (`detection.py` around `max_workers` and executor creation) | NumPy/photutils kernels may release GIL but each worker holds full image plus labels; deblend label caps hide cost; unordered future completion affects log order | frame throughput and peak RSS across workers 1,2,4,8; separate SEP/segmentation/DAO; record cache-hit time |
| WCS | Internal solver builds quads/KD trees and performs RANSAC/lstsq; ASTAP/solve-field are subprocesses; Gaia catalog loaded once per worker family but per-frame FITS/detection reads remain | Internal and ASTAP paths use `ThreadPoolExecutor`; `wcs_max_workers` defaults to 1 in parameters models | Oversubscription with external solver threads and simultaneous FITS/header writes; external timeout dominates; results are emitted as completed, not input order | per-engine solve time distribution, timeout rate, Gaia query time, workers×RAM; record accepted/failed frame ordering |
| Refbuild | Pandas/FITS reloads, Astropy `SkyCoord.match_to_catalog_sky`, optional KDTree; writes master TSVs | No single clear pool in `analysis/refbuild.py` | Catalog copies and coordinate object allocation dominate for large catalogs | N sources × N frames, memory profile, match radius and duplicate-policy sensitivity |
| Forced photometry | `forced_photometry.py:run_forced_photometry` reads/caches headers/images/sky/detection tables; per-frame aperture work is vectorized; growth curve uses 14 radii and one annulus-stat pass; outputs TSV/JSON | ThreadPool per frame; WCS/header caches protected by `Lock` | Each task retains an image and pandas frames; thread contention on disk and `photutils` calls; failed tasks can leave partial outputs | wall/RSS across frame count, workers, aperture mode, apcorr on/off; compare one-pass vs fallback apcorr |
| PSF | `step8_psf_photometry.py:Step6PSFWorker` builds per-frame ePSF, groups/fits sources, residual re-detects, and may do a second pass; several policies cap stars/iterations | Qt worker plus an internal thread pool | Per-frame ePSF and fit tables are memory-heavy; nested BLAS/OpenMP threads are not controlled by APEX | timing by ePSF build, fit, residual pass; workers and BLAS thread settings; peak RSS |
| CMD MCMC | `isochrone_mcmc.py:fit_isochrone_mcmc`; grid scan/DE initialisation plus emcee burn/production | Sampler is CPU-bound; no project-level process pool | Runtime scales with walkers×steps×stars; GUI progress does not imply cancellation at every likelihood call | fixed seed, stars, walkers, burn, steps; effective samples/sec and convergence diagnostics |
| LC periods | LS/PDM/BLS in `period_analysis_service.py`; bootstrap FAP permutes data and computes LS repeatedly | Bootstrap uses a thread pool up to `get_parallel_workers()`; PDM is vectorized over trials | Memory for permutations; BLAS/thread oversubscription; FAP runtime can dominate | n points × trials × bootstrap × workers, with deterministic seed and RSS |

## Parallelism map and correctness concerns

* The configuration exposes `parallel_mode` and `max_workers` (`parameters_cmd.py`/`parameters_lc.py`), but the audited science paths consistently instantiate `ThreadPoolExecutor`; there is no `ProcessPoolExecutor` path in the main pipeline. The mode field is therefore a compatibility/configuration surface, not evidence of process parallelism.
* `get_parallel_workers` (`apex/utils/constants.py`) supplies an automatic worker count. WCS has its own `wcs_max_workers` and PSF its own `psf_parallel_workers`, so a future GUI orchestration can oversubscribe if an outer task pool is added.
* WCS, detection, forced photometry and PSF collect futures as they complete. Output files are named by frame, which is safe for independent products, but logs/summary row order is not a deterministic ordering guarantee.
* Cancellation is cooperative in WCS (`should_stop` watcher and worker stop), and callback-based in detection/forced photometry. External subprocesses have timeouts, but a killed subprocess and partially copied WSL artifacts require a post-cancel cleanup test.
* Exceptions are converted to per-step failure by `PipelineRunner`, and several per-frame workers keep going after a frame exception. The audit found no single documented policy that distinguishes an incomplete run from a scientifically valid partial output; downstream consumers must inspect summary/QC files.
* Astropy/photutils/scipy may release the GIL in compiled kernels, but this is package/version dependent. No claim of linear scaling is supportable without the benchmark above.

## Caching and I/O

* Detection cache signatures (`apex/utils/cache_utils.py:build_file_signature`, `detection_cache_signature_matches`) use path/crop/schema/engine/mtime/size, with a relaxed rule to tolerate WCS header edits. This is practical for Step 5 in-place header writes but is weaker than a pixel-content hash.
* `variable_analysis_cache.py` uses a SHA-256 request key and atomic `.tmp` replace; this is the strongest cache contract in the tree. It is not representative of every step.
* `PipelineRunner` skips a complete step based on output presence/step markers and writes a simple manifest. The manifest does not include the full parameter snapshot, source hashes, environment variables, BLAS thread settings, or external executable versions.
* FITS files are repeatedly opened with `memmap=False` in core calibration and WCS paths. That avoids stale mapped headers and Windows file locking but increases copies and I/O. A safe optimisation must preserve header/write semantics and be benchmarked on local and external drives.

## Bottleneck ranking (current evidence)

1. PSF per-frame ePSF + iterative fitting (CPU and memory), especially with residual passes.
2. Detection deblending/segmentation and label arrays on crowded full-resolution frames.
3. Forced-photometry image/table reloads and growth-curve/apcorr work across many frames.
4. External WCS subprocess latency and Gaia/network/cache waits.
5. Bootstrap FAP / MCMC loops, whose cost is intentional but not bounded by a release-level budget.

## Safe next benchmark

Add a machine-readable benchmark manifest (without changing science code) containing commit, package versions, CPU/RAM, worker counts, input dimensions, cache state, per-stage wall time, peak RSS, output count, and failure/cancellation status. Run each workload at workers 1/2/4 and report medians plus variability. Until then, manuscript wording should say “supports frame-level threading” or “uses vectorized aperture calculations,” not “is faster” or “scales linearly.”
