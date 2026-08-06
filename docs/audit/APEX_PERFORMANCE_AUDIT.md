# APEX performance and parallelism audit

## Evidence boundary

The repository contains timing hooks and benchmark runners, but this audit did not find a committed, hardware-normalised baseline that supports a universal speed-up percentage. Therefore “parallel”, “vectorized”, or “optimized” below describes code structure only; it is not a measured performance claim.

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
