# Repository Guidelines

## Project Structure & Module Organization

APEX is a Python 3.10+ PyQt5 desktop application for astronomical photometry. The root `main.py` launches either mode. Mode entry points live in `apex/cmd/main.py` for cluster CMD work and `apex/lightcurve/main.py` for light-curve analysis. Shared application code is under `apex/`: `core/` contains state and file-management primitives, `config/` holds TOML-backed parameter models, `utils/` contains path, cache, I/O, logging, and astronomy helpers, `analysis/` contains science services, and `gui/` contains windows, workflow steps, tools, and widgets. SVG assets are in `apex/resources/`. Runtime settings are read from `parameters.toml`.

## Build, Test, and Development Commands

- `python3 main.py`: start the launcher and choose CMD or LC mode.
- `python3 apex/cmd/main.py`: run CMD mode directly.
- `python3 apex/lightcurve/main.py`: run LC mode directly.
- `python3 -m compileall apex main.py`: syntax-check the package after edits.
- `python3 -m pytest tests`: run tests once a `tests/` suite is added.

The README references `requirements.txt`, but this checkout does not currently include it. Keep dependency changes documented in `README.md` and include PyQt5, astropy, photutils, numpy, pandas, scipy, matplotlib, astroquery, and local astrometry.net requirements.

## Coding Style & Naming Conventions

Use 4-space indentation and standard Python naming: `snake_case` for functions, methods, and modules; `PascalCase` for classes; uppercase for constants such as `_CMD_MAIN`. Prefer `pathlib.Path` for filesystem work and keep workflow step files named by step and purpose, for example `step7_forced_aperture_phot.py`. Keep GUI changes consistent with the existing PyQt5 patterns in `apex/gui/main_window.py` and `apex/gui/workflow/step_window_base.py`.

## Testing Guidelines

There is no committed test suite yet. For new logic, add focused pytest tests under `tests/` using names like `test_step_paths_lc.py` or `test_period_analysis_service.py`. Prefer unit tests for `apex/utils/`, `apex/core/`, and `apex/analysis/`; GUI-heavy changes should at least pass `compileall` and include a short manual validation note.

## Commit & Pull Request Guidelines

Recent history uses concise lowercase prefixes such as `feat:`, `fix:`, and `remove:`. Keep commits scoped and imperative, for example `fix: preserve lc parameter cache path`. Pull requests should summarize the mode affected, list validation commands, note dependency or `parameters.toml` changes, and include screenshots for visible GUI updates.
