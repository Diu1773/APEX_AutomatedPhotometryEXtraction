# Changelog

All notable changes to APEX are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims
to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **MIT license** (`LICENSE`) — the project is now openly licensed.
- **pip-installable packaging** (`pyproject.toml`) with console entry points
  `apex` (headless CLI) and `apex-gui` (desktop launcher), and a dynamic version
  sourced from `apex.__version__`.
- **Headless command-line interface** (`apex` / `python -m apex`):
  - `apex doctor [--network]` — diagnose Python, dependencies, external WCS
    solvers (ASTAP / astrometry.net), and optional Gaia/SIMBAD reachability.
  - `apex config init|path|show` — manage `parameters.toml`.
  - `apex version`, `apex gui [--mode cmd|lc]`.
- **Headless pipeline orchestrator** (`apex/pipeline/`, Qt-free): step contract,
  run context, runner with prerequisite checks / idempotent reuse / JSON run
  manifest, and a step registry. Exposed via `apex run --mode {cmd,lc}
  [--steps 1-7] [--force] [--dry-run] [--result-dir] [--data-dir]`.
- **Step 1 (scan)** ported to headless execution: scans FITS headers, writes
  `step1_file_selection/{headers.csv,selection.json}`, and resolves the target
  from config or headers. Validated end-to-end on real data.
- Unit tests for the orchestrator (`tests/test_pipeline_runner.py`).

### Changed
- **Dependencies reorganized for maintainability.** The core install is now
  headless (scientific stack + CLI, no PyQt5); the desktop GUI moved to a `gui`
  extra. New extras: `gui`, `build`, `test`, `docs`, `dev`, `all`. Install the
  desktop app with `pip install apex-photometry[gui]`.
- Single-sourced the version via `apex/__init__.py`.

### Notes
- Steps 2–7 are recognised by the orchestrator (planning + completion detection
  of GUI-produced outputs) but not yet executable headless; they are being
  ported incrementally in the order 4 → 7 → 3 → 6 → 5 → 2.
