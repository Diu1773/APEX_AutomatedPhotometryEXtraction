# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What APEX Does

APEX is a PyQt5 desktop app for astronomical aperture and PSF photometry. It has two operational modes:

- **CMD mode** (`apex/cmd/`): Cluster photometry — source detection through CMD diagram and isochrone fitting, 12 steps.
- **LC mode** (`apex/lightcurve/`): Light curve analysis — multi-night photometry, detrending, and period analysis, 11 steps.

Both modes share a common pipeline through Step 7: file selection, crop, sky
preview, source detection, WCS plate solving, master catalog build, and forced
aperture photometry. CMD and LC branch after Step 7.

## Commands

```bash
# Run the launcher (choose mode interactively)
python main.py

# Run a specific mode directly
python apex/cmd/main.py
python apex/lightcurve/main.py

# Syntax-check the whole package after edits
python -m compileall apex main.py

# Run tests
python -m pytest tests

# Smoke-test all module imports
python scripts/smoke_steps.py
```

## Configuration

Runtime parameters are read from `parameters.toml` in the working directory (not committed). Copy from `parameters.example.toml` to create it. The TOML is parsed by `apex/config/parameters_cmd.py` and `apex/config/parameters_lc.py`. It covers I/O paths, target coordinates, telescope/camera specs, detection thresholds, WCS solving, and airmass settings.

The app caches the last-used parameter file path in `~/.apex/last_param.txt`.

## Architecture

### Package layout

```
apex/
  core/        — ProjectState (step tracking + JSON persistence), FileManager, Instrument
  config/      — TOML-backed parameter models for each mode + schema validators
  utils/       — step_paths*.py, photometry_utils, astro_utils, io_utils, cache_utils, logging_utils, …
  analysis/    — Science services: isochrone fitting (cmd/), light curve / detrend / period (light_curve/), multi-night merge (merge/)
  gui/
    main_window.py          — Unified main window; dispatches step windows by mode
    workflow/               — Step windows (step_window_base.py + stepN_*.py)
      cmd/                  — CMD-specific steps (8–12)
      lc/                   — LC-specific steps (8–11)
    widgets/                — Reusable widgets (e.g. image_viewer.py for zoomable FITS display)
    tools/                  — Standalone analysis dialogs (extinction fit, Gaia 3D viewer, transit fitting, …)
  resources/   — SVG assets (logo_base.svg, logo_cmd.svg, logo_lc.svg)
  cmd/main.py  — CMD entry point
  lightcurve/main.py — LC entry point
main.py        — Root launcher; spawns subprocess for chosen mode
```

### Key relationships

- `main_window.py` owns `ProjectState` and routes `_open_step_window(step_index)` to the correct step window class depending on `mode`.
- All step windows inherit from `step_window_base.py`.
- Step results land under the path helpers in `apex/utils/step_paths*.py`.
  Shared steps write `step1_file_selection/` through `step7_forced_phot/`;
  mode-specific steps use `cmd_*/` or `lc_*/` directories.
- Path helpers in `apex/utils/step_paths.py` (shared Step 1-7),
  `step_paths_cmd.py` (CMD Step 8-12), and `step_paths_lc.py` (LC Step 8-11)
  are the canonical source for output paths — always use them instead of
  constructing paths manually.
- Caches (header scan, detection, WCS) live under `result_dir/cache/` and are managed by `cache_utils.py` and `header_cache.py`.

## Coding Style

- 4-space indentation; `snake_case` for functions/modules, `PascalCase` for classes, `UPPER_CASE` for constants.
- Use `pathlib.Path` for all filesystem work.
- Step files are named by current UI step index and purpose, for example
  `step7_forced_aperture_phot.py` and `lc/step9_lightcurve_builder.py`.
- GUI changes must follow existing PyQt5 patterns in `main_window.py` and `step_window_base.py`.

## Commits and PRs

Use concise lowercase prefixes: `feat:`, `fix:`, `refactor:`, `remove:`, `chore:`. Keep commits scoped and imperative. PRs should note which mode is affected, list validation commands run, and include screenshots for visible GUI changes.
