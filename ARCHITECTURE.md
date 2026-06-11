# APEX Architecture

이 문서는 현재 APEX 코드의 모듈 경계, workflow dispatch, 데이터 흐름을
정리합니다. 실제 step 번호와 출력 경로의 최종 기준은
`apex/gui/main_window.py`와 `apex/utils/step_paths*.py`입니다.

## Runtime Entry Points

| Entry point | 역할 |
| --- | --- |
| `main.py` | CMD/LC launcher, packaged-app setup, `--smoke` |
| `apex/cmd/main.py` | CMD mode direct launch |
| `apex/lightcurve/main.py` | LC mode direct launch |
| `apex/gui/main_window.py` | mode별 step 목록과 window dispatch |

두 모드는 `MainWindowWorkflow(mode="cmd" | "lc")`를 공유합니다.
Step 1은 mode-specific wrapper를 사용하고, Steps 2-7은 동일한 window
class를 사용합니다.

## Package Layout

```text
apex/
  analysis/
    astrometry/          internal quad matcher and WCS fitter
    cmd/                 isochrone fitting and CMD science services
    light_curve/         ensemble, detrend, period, SysRem, output services
    merge/               multi-workspace scanning, ID reconciliation, build
  benchmark/             benchmark implementations used by CLI wrappers
  config/
    parameters_cmd.py    CMD runtime parameter mapping
    parameters_lc.py     LC runtime parameter mapping
    schema.py            configuration validation helpers
  core/
    project_state.py     workflow completion and project state
    file_manager.py      selected FITS and result-path ownership
    instrument.py        instrument model
  gui/
    main_window.py       unified CMD/LC shell
    workflow/            shared and mode-specific step windows
    tools/               analysis and diagnostic tools
    widgets/             reusable PyQt5 widgets and FITS viewers
  utils/
    step_paths*.py       canonical result-directory paths
    cache*.py            cache and manifest helpers
    photometry*.py       aperture photometry and table loading
    astro_utils.py       time, airmass, filter, coordinate helpers
    gaia_*.py            Gaia query, matching, and derived fields
```

GUI modules should coordinate user interaction and delegate reusable numerical
logic to `apex/analysis/` or `apex/utils/`.

## Workflow Dispatch

`MainWindowWorkflow._open_step_window()` maps a zero-based UI index to a step
window.

### Shared Steps

| UI Step | Module | Output |
| ---: | --- | --- |
| 1 | `workflow/{cmd,lc}/step1_file_selection.py` | `step1_file_selection/` |
| 2 | `workflow/step2_crop_selector.py` | `step2_crop/` |
| 3 | `workflow/step3_sky_preview.py` | `step3_sky_preview/` |
| 4 | `workflow/step4_source_detection.py` | `step4_detection/` |
| 5 | `workflow/step5_wcs_plate_solving.py` | `step5_wcs/` |
| 6 | `workflow/step6_ref_build.py` | `step6_refbuild/` |
| 7 | `workflow/step7_forced_aperture_phot.py` | `step7_forced_phot/` |

### CMD Branch

| UI Step | Module | Output |
| ---: | --- | --- |
| 8 | `workflow/cmd/step8_psf_photometry.py` | `cmd_psf/` |
| 9 | `workflow/cmd/step9_master_id_editor.py` | `cmd_selection/` |
| 10 | `workflow/cmd/step10_zeropoint_calibration.py` | `cmd_zeropoint/` |
| 11 | `workflow/cmd/step11_cmd_plot.py` | `cmd_plot/` |
| 12 | `workflow/cmd/step12_isochrone_model.py` | `cmd_isochrone/` |

### LC Branch

| UI Step | Module | Output |
| ---: | --- | --- |
| 8 | `workflow/lc/step8_target_selection.py` | `lc_selection/` |
| 9 | `workflow/lc/step9_lightcurve_builder.py` | `lc_lightcurve/` |
| 10 | `workflow/lc/step10_detrend_merge.py` | `lc_detrend/` |
| 11 | `workflow/lc/step11_period_analysis.py` | `lc_period/` |

Legacy path helper aliases remain for old result folders, but new writes must
use the canonical names above.

## Data Flow

```text
raw FITS
  -> Step 1 header/frame manifest
  -> Step 2 optional cropped FITS
  -> Step 3 frame quality preview
  -> Step 4 detections + frame_quality.csv
  -> Step 5 solved WCS + Gaia catalog/QC
  -> Step 6 master catalog + stable source IDs
  -> Step 7 per-frame forced photometry TSV + index/statistics
       -> CMD: PSF -> selection -> zeropoint -> CMD -> isochrone
       -> LC: target/comps -> raw LC -> detrend/merge -> period products
```

Step 7 owns the shared science photometry table. Its principal outputs are:

- `photometry_index.csv`
- `photometry_<frame>.tsv`
- `master_sources.csv`
- `frame_stats.csv`
- `apcorr_summary.csv`
- `centering_stats.csv`
- `filter_frames.json`

## Internal WCS Solver

`apex.analysis.astrometry.solver.solve()` implements the built-in solver:

1. Load Step 4 source coordinates and a Gaia catalog.
2. Select bright, usable source and catalog subsets.
3. Build translation/rotation-invariant quad codes.
4. Match codes with a 4-D `cKDTree`.
5. Verify candidate transforms with RANSAC and unique inlier counting.
6. Fit a TAN WCS and optionally a SIP model when enough pairs exist.
7. Sigma-clip, rematch, refit, and independently verify on the full catalog.

The Step 5 worker applies a per-frame hint chain:

```text
FITS header pointing -> Step 1 target -> local blind retry
```

Each hint uses a catalog-size ladder to avoid both sparse-pattern and
crowded-field failures. Solver provenance and QC values are written to the
same Step 5 summary path used by external solvers.

## Configuration Flow

`parameters.example.toml` is the committed default. Runtime code resolves a
local `parameters.toml`, then maps TOML paths to mode-specific `params.P`
attributes through `parameters_cmd.py` or `parameters_lc.py`.

Rules:

- Add new canonical defaults to `parameters.example.toml`.
- Map the key in every mode that consumes it.
- Keep path and unit names explicit, for example `_px`, `_arcsec`, `_adu`.
- Treat `docs/parameter-inventory.md` as generated diagnostic output, not a
  hand-maintained schema.

## Cache and Output Ownership

Persistent outputs belong to their workflow step directory. In-memory FITS,
plot, or table caches do not imply persisted step reuse.

Complete-output reuse currently matters most for:

- Step 4 detection
- Step 5 local astrometry.net outputs
- Step 6 master/reference build
- Step 7 forced photometry
- CMD Step 8 PSF photometry

See `docs/cache-manager-design.md` for invalidation and migration rules.

## Release Architecture

The Windows release pipeline is:

```text
source preflight
  -> compileall + pytest
  -> PyInstaller onedir bundle
  -> APEX.exe --smoke
  -> portable ZIP
  -> Inno Setup installer
  -> release verification
```

`deploy/apex_windows.spec` must collect package data required at runtime,
including scientific-package resources that are not discovered through normal
imports. The smoke test imports critical workflow, Gaia/SIMBAD, and Tools-menu
modules without opening the GUI.

## Extension Guidelines

### Add a shared step

1. Add a window under `apex/gui/workflow/`.
2. Add a canonical helper in `apex/utils/step_paths.py`.
3. Add dispatch and completion state in `main_window.py`.
4. Keep reusable science logic outside the window module.
5. Add focused tests for path, cache, and numerical behavior.

### Add a mode-specific step

1. Add the window under `workflow/cmd/` or `workflow/lc/`.
2. Add its path helper to `step_paths_cmd.py` or `step_paths_lc.py`.
3. Update the mode step list and dispatch.
4. Document inputs, outputs, and invalidation dependencies.
