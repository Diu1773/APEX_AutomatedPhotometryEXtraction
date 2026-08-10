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

**JSON is the single authority** (2026-08-06, "TOML 완전 제거"): each workspace
owns an `apex_config.json` (uncommitted) covering I/O paths, target
coordinates, instrument specs, detection thresholds, WCS solving, and airmass
settings. All read/write goes through `apex/config/config_io.py` —
`load_config_data()` / `save_config_data()` — and is parsed into runtime
params by `apex/config/parameters_cmd.py` / `parameters_lc.py`.

- **Legacy TOML**: a `parameters*.toml` is migrated to its JSON sibling on
  first load (`parameters.toml` → `apex_config.json`,
  `parameters_X.toml` → `apex_config_X.json`) and never re-read afterwards —
  a newer TOML mtime only produces an "IGNORED" warning. Never write TOML.
- **Choosing a workspace**: GUI entry points accept
  `--params <apex_config.json|legacy.toml|dir>`; without it the repo-root
  config is used. (`~/.apex/last_param.txt` was documented before but never
  implemented — do not rely on it.)
- **Identity guard**: `check_workspace_identity()` warns when `target.name`
  does not match the io paths (the "M13 paths + NGC 6811 target" mixing
  accident). Warn-only by design.
- `parameters.example.toml` remains only as the commented template; first-run
  bootstrap converts it to `apex_config.json`.

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

### GUI layout rules (`apex/gui/layout_rules.py`)

Window sizing/anti-clipping is centralized — do not re-solve it per window.

- **Window sizing is automatic.** Every step/tool window auto-fits to its content
  and clamps to the monitor on first show. Step/tool windows get this free via
  `WindowChromeMixin.showEvent`; raw `QMainWindow` windows must mix in
  `AutoFitMixin` *before* `QMainWindow` (sip routes the `showEvent` virtual to a
  Python method only when it's defined on the class, never an instance attribute).
- **Dialogs.** Use `FittedDialog(parent)` instead of `QDialog(parent)` for any
  modal/popup dialog so it auto-fits and clamps like the windows do. Parameter
  dialogs built via `run_param_dialog`/`build_scroll_param_dialog` are already
  clamped (`configure_parameter_dialog`). A dialog that stacks many groups must
  put them in a `QScrollArea` with the button row *outside* the scroll, so a
  short screen scrolls instead of clipping Save/Cancel.
- **Never `setMinimumSize`/`resize` larger than the screen.** A minimum larger
  than the monitor permanently clips the bottom row. Pass desired sizes through
  `clamp_to_screen(w, h, self)`.
- **Embedded matplotlib canvases must be tamed.** A bare `FigureCanvas` reports
  `minimumSizeHint() == 10×10` and collapses to a sliver next to a table/controls.
  Wrap with `tame_canvas(canvas)` (min size + Expanding) and add it with `stretch=1`.
- **Splitters with a plot pane** must call `prevent_collapse(splitter)` (or
  `setChildrenCollapsible(False)`) so a pane can't be dragged/laid out to 0 px.

### Button hierarchy & color (`apex/gui/theme.py`)

The global stylesheet is the single source of look; `apply_theme(app)` is called
in all three entry points (`main.py`, `apex/cmd/main.py`, `apex/lightcurve/main.py`).

- **Never hand-paint a button** with `setStyleSheet("background-color: …")`. That
  is exactly the inconsistency to avoid. Assign a *role* instead with
  `style_button(btn, variant, height=Tokens.H_*)`.
- **One hierarchy:** `primary` (filled accent — the single main action: Run /
  Next / Save) · `danger` (filled red — Stop / destructive) · `ghost` (accent
  text — tertiary: Log / 가이드) · *no variant* = neutral default (everything
  else: Parameters / Browse / Export / Previous).
- **One size scale:** `Tokens.H_ACTION` (38, bottom action row) ·
  `Tokens.H_BUTTON` (32, standard) · `Tokens.H_COMPACT` (28, header cluster).
- Disabled primary/danger are themed automatically (muted fill) — just call
  `setEnabled(...)`, don't swap stylesheets per state.

### Icons & spacing grid (the "keyline" discipline)

- **Button glyphs come from `theme.ICON`**, never pasted emoji literals. Bare
  emoji (⚙ 📜 🔒 📂 💾) render as multicolor OS emoji on Windows and break the
  flat look; `ICON` appends U+FE0E (text presentation) to force the monochrome
  symbol. Use `ICON["params"|"log"|"guide"|"locked"|"input"|"output"|…]`.
- **Layout snaps to an 8px grid.** Margins/spacing use `Tokens.MARGIN` (16),
  `Tokens.S3` (12), `Tokens.GAP`/`Tokens.S2` (8) — never 5/6/10. The shared
  window bases already set the outer rhythm; new panels should follow it.
- **Placement is fixed by convention:** bottom action row = `[Previous] … (stretch) … [Next]`;
  run bar = `[Run][Stop] … (stretch) … [Log]`; header cluster (right of title) =
  subclass actions, then Parameters, then Log, then 가이드.

## Commits and PRs

Use concise lowercase prefixes: `feat:`, `fix:`, `refactor:`, `remove:`, `chore:`. Keep commits scoped and imperative. PRs should note which mode is affected, list validation commands run, and include screenshots for visible GUI changes.

## Review Domain Notes

Domain facts for code review (the generic `/review-math`, `/review-deps`,
`/review-perf` commands read this section to build project context).

### Math / numerics

- **Magnitude system**: instrumental (`mag_inst`) → per-reference → absolute. Larger value = fainter source. Magnitude errors scale as `MAG_ERR_COEFF / SNR`.
- **Airmass**: `X ≈ sec(z)`; `X = 1` at zenith (alt 90°), diverges near the horizon.
- **Extinction model**: `m_ij = s_i + z_j + k1·X_j` (`s_i` = star brightness, `z_j` = frame offset, `k1` = extinction coefficient). `z_j` and `k1` are degenerate, so the frame-offset basis must be SVD-projected to remove any constant and airmass-linear component — otherwise `k1` is not identifiable.
- **SYSREM** (Tamuz+ 2005, MNRAS 356, 1466): iterate `r_ij -= a_i·c_j` to convergence. Each `a_i`/`c_j` update needs its denominator `Σ_j w_ij·c_j²` (resp. `Σ_i w_ij·a_i²`) `> 0`; missing data carry weight 0.
- **PDM** (Stellingwerf 1978): `θ = (Σ_k SS_k / Σ_k DOF_k) / σ²_total`. Each bin needs ≥ 2 points (DOF ≥ 1); `σ²_total` uses sample variance (`ddof=1`).
- **BJD_TDB**: from `JD_UTC` including light-travel + relativistic corrections (~±8 min); easy to get the sign wrong.
- **WCS TAN**: `CDELT[0] < 0` — RA decreases toward increasing pixel x (east is −x).
- **Weights**: photometric weight is `w = 1/σ²`. When feeding `np.linalg.lstsq`, rows/values are multiplied by `√w = 1/σ` — do not conflate the two forms.

### Architecture / dependencies

- **Layers**: `gui/` (presentation: Qt, workflow steps, tools) → `analysis/` (pure science calc), `core/` (state/config/files), `utils/` (shared), `config/` (TOML param models). Allowed: gui→analysis/core/utils, analysis→utils, core→utils. Forbidden: analysis/utils/config/core → gui.
- **Path helpers**: `step_paths.py` (shared Step 1–7), `step_paths_cmd.py` (CMD 8–12), `step_paths_lc.py` (LC 8–11). Never build output paths by string concat.
- **Filter keys**: always via `normalize_filter_key()`. Johnson = uppercase (B,V,R,I), SDSS = lowercase (g,r,i,z), narrowband = title case (Ha, OIII).
- **source_id**: int64; convert via `coerce_int64_source_id()` (direct casts risk sign errors).
- **ProjectState**: `store_step_data(key, dict)` / `get_step_data(key)`; a mistyped key silently returns None.
- **Cache invalidation**: `StepCacheManager`; a parameter missing from the cache signature means stale results are reused after that parameter changes.
- **QThread safety**: mutate GUI widgets only on the main thread; workers emit signals, main-thread slots touch widgets.

### Performance

- **Typical scale**: N_frames 50–500 (single night), up to ~3000 (multi-night); N_stars 100–5000 (master catalog), 20–200 (references); ~16M px/frame (4000×4000); up to 50,000 PDM/LS trial periods.
- **Frame loops**: Step 7 forced phot runs in **processes** headless (`APEX_FORCEDPHOT_PROCESSES=0` forces threads) and threads in the GUI; LC step9 `_build_star_mag_series` file iteration.
- **Preload cache**: LC frames go through `FramePhotometryCache` (`apex/utils/photometry_loader.py`), keyed by `(result_dir, frame)` and bounded by bytes — flag code that bypasses it and re-reads per star/frame.
- **Worker count**: use `get_parallel_workers(params, stage=…, frame_bytes=…)` (`apex/utils/constants.py`), never hardcode. Stage ceilings and the RAM budget are measured, not assumed (`benchmark/perf/20260807/RESULTS.md`).
- **Threads do not always help numeric work.** Measured on Step 7: 83 % CPU at one worker, 109 % at twelve — it never leaves one core, because photutils' per-aperture statistics hold the GIL between many small numpy calls. Twelve threads cost 4.2× wall time. Ask what a stage *waits on*: frame reads gain from threads (34.1 s → 18.7 s), per-aperture computation needs processes (244.2 s → 69.6 s, byte-identical output).
- **Batched numerics**: `_pdm_theta_vectorized` caps each batch at ~50 MB by design — preserve such memory bounds.
- **Qt tables**: bulk `setItem` is faster with `setSortingEnabled(False)` around the batch.

## Research OS 세션 브리지

연구 질문이나 실험 방향을 새로 잡을 때는 제어층 하네스를 먼저 동기화한다.
컨텍스트에 `SESSION_SYNC_V1`이 없으면 작업 전에 한 번 실행한다.

```text
python -X utf8 C:/Users/bmffr/Desktop/Main/scripts/session_sync.py --cwd .
```

원시 아이디어는 `C:/Users/bmffr/Desktop/Main/RESEARCH_INBOX.md`에 남길 수 있다.
하네스는 먼저 질문 후보만 확장하며, 논문 검색·novelty 판정·승인·프로젝트 파일 변경은
사람이 선택한 다음 단계에서만 한다.
ResearchCandidate 검토용 또는 승인된 theory/experiment handoff는 이 프로젝트의
`.research-os/handoff/`에 JSON으로 도착한다. 새 작업 전 최신 파일과 이 문서의 기존 테스트를 읽고, 자동 실행하지 않는다.
