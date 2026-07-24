# Changelog

All notable changes to APEX are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims
to follow [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Per-frame predicted detection limit + depth QC gate (Step 7).** New
  Qt-free module `apex/analysis/detection_limit.py`: the 50%-completeness
  magnitude of a frame is predicted from its background noise and PSF peak
  fraction alone via `m50 = ZP − 2.5·log10(S/N₅₀·σ_e/p_peak)`, with
  `S/N₅₀ = 4.05 ± 0.18` calibrated by artificial-star injection into 7 real
  cluster frames (residual RMS ≈ 0.05 mag; constant in
  `apex/utils/constants.py`, gain from the runtime config). Step 7 now writes
  `sky_sigma_e`, `p_peak_frame`, `predicted_m50`, `observed_m50` (empirical
  detection-fraction rolloff over the master catalog), `depth_delta_mag`, and
  `depth_qc_flag` to `frame_stats.csv`, flagging frames whose observed depth
  deviates from the prediction by more than `depth_qc_tolerance_mag`
  (default 0.5 mag — focus / clouds / tracking / defect suspects). The same
  injection-calibrated completeness also predicts the number of detected
  real catalog stars to ~6%, grounding the predicted-vs-observed comparison
  (see `validation/paper/논문작업/COMPLETENESS_REALFRAME_INVESTIGATION.md`).
  Tests: `tests/test_detection_limit.py` (unit + 7-run calibration
  reproduction, residual RMS 0.048 mag). The gate is configurable via
  `[photometry.depth_qc] tolerance_mag / min_snr` (both modes), included in
  the Step 7 cache signature so changing them triggers recomputation.
  Validated on real data (M67 9 frames g/r/i: |predicted − observed| ≤ 0.073
  mag; NGC 6811 3 frames) and end-to-end in the CMD GUI (identical numbers,
  flag + log warning fire when the tolerance is exceeded).
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
