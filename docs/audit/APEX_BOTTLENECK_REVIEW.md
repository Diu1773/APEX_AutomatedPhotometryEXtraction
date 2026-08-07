# APEX Bottleneck whole-code review

Audit date: 2026-08-07. This is a source and environment review of the Bottleneck
integration, not a claim that every GUI branch was executed.

## Scope and inventory

`bottleneck` appears in four dependency/runtime surfaces (`pyproject.toml`,
`requirements.txt`, `deploy/apex_windows.spec`, and `apex/cli.py`) and in the
compatibility module `apex/utils/fast_stats.py`. The wrapper is imported by exactly
four production consumers:

| Consumer | Wrapper functions actually used | Data role |
|---|---|---|
| `apex/analysis/calibration.py` | `nanmedian`, `nanmean`, `finite_nanmedian` | master bias/dark/flat stacks and summary statistics |
| `apex/analysis/detection.py` | `finite_nanmedian`, `finite_nanstd`, `robust_median_mad` | sky/shape/QC summaries |
| `apex/analysis/cosmetic.py` | `robust_median_mad` | hot-pixel/star-protection statistics |
| `apex/gui/workflow/step4_source_detection.py` | `robust_median_mad`, finite helpers | GUI detection review statistics |

The wrapper exposes `nanmedian`, `median`, `nanmean`, `nanstd`, `nansum` and
`nanmax`. It calls Bottleneck when importable and NumPy otherwise. It does not
change the scientific formula, and it deliberately normalises `finite_values` to
float64 before robust summaries.

A repository-wide search found 555 direct `numpy.nan*` calls in 62 Python files:

| Tree | Calls | Interpretation |
|---|---:|---|
| `apex/analysis` | 170 | mixed: some large arrays, many residual/QC summaries |
| `apex/gui` | 302 | mostly interactive diagnostics and worker summaries |
| `apex/benchmark` | 49 | should remain independent NumPy reference calculations |
| `apex/utils` | 32 | small helpers and wrapper fallback |
| `apex/io` | 2 | export summaries |

Direct NumPy use is therefore not evidence of a bug or of an “unoptimised” code
path. It does mean that the manuscript must describe Bottleneck as a selected
backend, not as a package-wide acceleration layer.

## Numerical parity spot-check

On `.venv-deploy` (Python 3.12.3, NumPy 2.4.4, Bottleneck 1.6.0), random and edge
arrays were compared for all wrapper operations at axes `None`, `0` and `-1`, with
float32/float64, integer, empty and all-NaN cases; `nanstd` was also checked for
`ddof=0,1,2`. No value or NaN-mask mismatch was observed in the spot-check. This
is not a substitute for a committed regression test; the absence of a
`test_fast_stats.py` test is a release finding.

## Performance spot-check (not a paper result)

The same environment was used for a three-repeat minimum-time comparison at
`axis=0`:

| Array | `nanmedian` NumPy/Bottleneck | `nanmean` | `nanstd` | `nansum` |
|---|---:|---:|---:|---:|
| `float32`, `(20,2048,2048)` | 4.23× | 3.85× | 3.50× | 2.20× |
| `float32`, `(14,500)` | 4.88× | 15.1× | 5.27× | 2.80× |
| `float32`, `(1_000_000,)` | 0.96× | 3.92× | 3.99× | 3.17× |

Ratios are hardware-, dtype- and shape-dependent. The second row completed in
sub-millisecond absolute time, so its ratio is not evidence that Bottleneck is a
meaningful end-to-end photometry speedup. These measurements justify using the
wrapper in calibration-like reductions, but not replacing every direct NumPy call.

## Findings before submission

1. Add a parity regression test covering the edge cases above. Compare values and
   NaN masks; do not make warning text part of the contract.
2. Keep the packaging statement consistent: Bottleneck is a declared core/required
   dependency, while the wrapper has a defensive NumPy fallback. The paper should
   say exactly that, rather than calling Bottleneck an optional extra.
3. Record `fast_stats.HAS_BOTTLENECK`, package versions, array shape, worker count,
   BLAS thread settings and cache state in any published benchmark manifest.
4. Benchmark only high-volume candidates before routing more calls through the
   wrapper: calibration stacks, enabled overscan regions, ensemble/detrend matrices
   and large WCS catalogues. Keep benchmark/reference scripts on direct NumPy so
   the comparison remains independent.
5. Do not put the local ratios above in the paper. If a speed claim is desired,
   rerun a hardware-normalised benchmark with cold/warm cache, RSS, worker count,
   package versions and repeated runs.

## Why and when to rerun the pipeline

The Step 0 run currently in progress in another session is the pre-optimisation
baseline. It should finish and be preserved rather than being overwritten. Its role is
to provide the numerical reference needed to distinguish a genuine performance change
from a silent change in calibration, masking or row selection.

Bottleneck is not being introduced as a package-wide speed switch. The first changes
target memory and data movement: streaming calibration avoids image-count-scaled OOM,
LC caches avoid repeated CSV/pandas work, and bounded workers avoid CPU/RAM/I/O
oversubscription. Only after those changes are measured should large reductions in a
calibration or compact ensemble matrix be routed through `fast_stats`.

The acceptance sequence is:

1. Record the baseline manifest (raw input, output directory, commit, parameters,
   versions, worker count, wall time, peak RSS and Bottleneck active/fallback state).
2. Apply one optimisation class at a time and compare a fixed subset against the
   baseline, including values and NaN masks—not merely elapsed time.
3. Freeze the accepted code and perform a clean Step 0-to-downstream run from raw
   input. This is the only run eligible for final manuscript numbers.

The reason for the full rerun is scientific provenance: all downstream tables and
figures must come from one coherent code/configuration state, with no stale cache or
mixed pre/post-optimisation products. Git history cleanup is useful documentation but
is not a prerequisite for this numerical acceptance gate.

## Manuscript-safe sentence

> Bottleneck is a declared core dependency used through `apex.utils.fast_stats` for
> selected NaN-aware reductions in calibration and detection-related paths; the
> wrapper falls back to NumPy if the compiled extension cannot be imported. Other
> modules retain direct NumPy reductions, and no package-wide speed claim is made.
