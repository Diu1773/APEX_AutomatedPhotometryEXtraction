"""
Lightcurve-mode step path helpers.
Re-exports shared paths and retains the pre-PSF LC helper names as stable APIs.

LC pipeline layout:
  step1_file_selection/   File selection (multi-night)
  step2_crop/             Image crop
  step3_sky_preview/      Sky preview QC
  step4_detection/        Source detection + frame QC
  step5_wcs/              WCS plate solving
  step6_refbuild/         Reference catalog build
  step7_forced_phot/      Forced aperture photometry
  cmd_psf/                Step 8 optional PSF photometry
  lc_selection/           Step 9 target/comparison selection
  lc_lightcurve/          Step 10 light curve builder
  lc_detrend/             Step 11 detrend & night merge
  lc_period/              Step 12 period analysis

Functions such as step8_selection_dir() keep their historical names so existing
projects and plugins do not break when the optional PSF step is inserted.
"""

from __future__ import annotations
import csv
from pathlib import Path
from typing import Union

# Re-export everything from the shared base
from apex.utils.step_paths import *  # noqa: F401, F403
from apex.utils.step_paths import step_dir, _as_path

PathLike = Union[str, Path]

# ── LC-specific step directories ──────────────────────────────────────────────

LC_SELECTION_DIRNAME = "lc_selection"
LC_LC_DIRNAME        = "lc_lightcurve"
LC_DETREND_DIRNAME   = "lc_detrend"
LC_PERIOD_DIRNAME    = "lc_period"


def _is_usable_lightcurve_csv(path: Path) -> bool:
    try:
        if not path.is_file() or path.stat().st_size <= 0:
            return False
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            reader = csv.reader(handle)
            header = next(reader, [])
            first_row = next(reader, [])
        return len(header) >= 2 and len(first_row) >= 2
    except (OSError, UnicodeError, csv.Error):
        return False


def step8_selection_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, LC_SELECTION_DIRNAME)


def step9_lc_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, LC_LC_DIRNAME)


def step10_detrend_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, LC_DETREND_DIRNAME)


def step11_period_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, LC_PERIOD_DIRNAME)


# ── Legacy-aware input resolution (read-only; for old RESULT_* workspaces) ─────

def selection_input_dir(result_dir: PathLike) -> Path:
    """Selection/master-catalog dir for READING a (possibly legacy) workspace."""
    return first_existing_dir(result_dir, LC_SELECTION_DIRNAME, "step9_selection")


def selection_input_dirs(result_dir: PathLike) -> list[Path]:
    """Return all existing selection input dirs, newest schema first."""
    root = _as_path(result_dir)
    candidates = [root / LC_SELECTION_DIRNAME, root / "step9_selection"]
    existing: list[Path] = []
    for path in candidates:
        if path.exists() and path not in existing:
            existing.append(path)
    return existing or [candidates[0]]


def lightcurve_input_dir(result_dir: PathLike) -> Path:
    """Light-curve dir for READING a (possibly legacy) input workspace."""
    return first_existing_dir(result_dir, LC_LC_DIRNAME, "step10_lightcurve")


def step9_selection_dir(result_dir: PathLike) -> Path:
    """Deprecated alias; use the stable step8_selection_dir() API."""
    import warnings
    warnings.warn(
        "step9_selection_dir is deprecated; use step8_selection_dir()",
        DeprecationWarning, stacklevel=2,
    )
    return step8_selection_dir(result_dir)


def step10_lc_dir(result_dir: PathLike) -> Path:
    """Deprecated alias; use the stable step9_lc_dir() API."""
    import warnings
    warnings.warn(
        "step10_lc_dir is deprecated; use step9_lc_dir()",
        DeprecationWarning, stacklevel=2,
    )
    return step9_lc_dir(result_dir)


def step11_detrend_dir(result_dir: PathLike) -> Path:
    """Deprecated alias; use the stable step10_detrend_dir() API."""
    import warnings
    warnings.warn(
        "step11_detrend_dir is deprecated; use step10_detrend_dir()",
        DeprecationWarning, stacklevel=2,
    )
    return step10_detrend_dir(result_dir)


def step12_period_dir(result_dir: PathLike) -> Path:
    """Deprecated alias; use the stable step11_period_dir() API."""
    import warnings
    warnings.warn(
        "step12_period_dir is deprecated; use step11_period_dir()",
        DeprecationWarning, stacklevel=2,
    )
    return step11_period_dir(result_dir)


def find_best_lightcurve_csv(result_dir: PathLike, star_id: int) -> Path | None:
    """Find the preferred usable light-curve CSV for a target star."""
    candidates = list_lightcurve_csvs(result_dir, target_id=star_id)
    return candidates[0] if candidates else None


def step5_photometry_dir(result_dir: PathLike) -> Path:
    """Legacy alias for Step 7 forced aperture photometry output."""
    from apex.utils.step_paths import step7_forced_phot_dir
    return step7_forced_phot_dir(result_dir)


# ── LC short-name aliases used by analysis modules ───────────────────────────

def step10_dir(result_dir: PathLike) -> Path:
    """DEPRECATED — misleading name; resolves to lc_lightcurve/ (step9). Use step9_lc_dir()."""
    import warnings
    warnings.warn(
        "step10_dir is deprecated; use step9_lc_dir()",
        DeprecationWarning, stacklevel=2,
    )
    return step9_lc_dir(result_dir)


def step11_dir(result_dir: PathLike) -> Path:
    """DEPRECATED — misleading name; resolves to lc_detrend/ (step10). Use step10_detrend_dir()."""
    import warnings
    warnings.warn(
        "step11_dir is deprecated; use step10_detrend_dir()",
        DeprecationWarning, stacklevel=2,
    )
    return step10_detrend_dir(result_dir)


def step10_current_meta_path(result_dir: PathLike, target_id: int) -> Path:
    """Path to current Step 11 detrend metadata; name retained for compatibility."""
    return step10_detrend_dir(result_dir) / f"result_ID{int(target_id)}_current.json"


def step10_current_lc_path(result_dir: PathLike, target_id: int) -> Path:
    return step10_detrend_dir(result_dir) / f"lightcurve_ID{int(target_id)}_current.csv"


def _step10_current_lc_path(result_dir: PathLike, target_id: int) -> Path:
    return step10_current_lc_path(result_dir, target_id)


def step10_current_params_path(result_dir: PathLike, target_id: int) -> Path:
    return step10_detrend_dir(result_dir) / f"params_ID{int(target_id)}_current.csv"


def step10_current_summary_path(result_dir: PathLike, target_id: int) -> Path:
    return step10_detrend_dir(result_dir) / f"summary_ID{int(target_id)}_current.txt"


def step10_current_plot_path(result_dir: PathLike, target_id: int) -> Path:
    return step10_detrend_dir(result_dir) / f"plot_ID{int(target_id)}_current.png"


def step10_current_global_zp_path(result_dir: PathLike, target_id: int) -> Path:
    return step10_detrend_dir(result_dir) / f"global_zp_ID{int(target_id)}_current.csv"


def step10_current_global_mean_path(result_dir: PathLike, target_id: int) -> Path:
    return step10_detrend_dir(result_dir) / f"global_mean_ID{int(target_id)}_current.csv"


def step10_current_global_diag_path(result_dir: PathLike, target_id: int) -> Path:
    return step10_detrend_dir(result_dir) / f"global_diagnostics_ID{int(target_id)}_current.json"


def step10_history_dir(result_dir: PathLike) -> Path:
    return step10_detrend_dir(result_dir) / "_history"


def step11_current_meta_path(result_dir: PathLike, target_id: int) -> Path:
    """Legacy alias for Step 10 detrend metadata."""
    return step10_current_meta_path(result_dir, target_id)


def step11_current_lc_path(result_dir: PathLike, target_id: int) -> Path:
    """Legacy alias for Step 10 detrend light curve."""
    return step10_current_lc_path(result_dir, target_id)


def _step11_current_lc_path(result_dir: PathLike, target_id: int) -> Path:
    return step10_current_lc_path(result_dir, target_id)


def step11_current_params_path(result_dir: PathLike, target_id: int) -> Path:
    """Legacy alias for Step 10 detrend parameters."""
    return step10_current_params_path(result_dir, target_id)


def step11_current_summary_path(result_dir: PathLike, target_id: int) -> Path:
    """Legacy alias for Step 10 detrend summary."""
    return step10_current_summary_path(result_dir, target_id)


def step11_current_plot_path(result_dir: PathLike, target_id: int) -> Path:
    """Legacy alias for Step 10 detrend plot."""
    return step10_current_plot_path(result_dir, target_id)


def step11_current_global_zp_path(result_dir: PathLike, target_id: int) -> Path:
    """Legacy alias for Step 10 global zeropoints."""
    return step10_current_global_zp_path(result_dir, target_id)


def step11_current_global_mean_path(result_dir: PathLike, target_id: int) -> Path:
    """Legacy alias for Step 10 global mean table."""
    return step10_current_global_mean_path(result_dir, target_id)


def step11_current_global_diag_path(result_dir: PathLike, target_id: int) -> Path:
    """Legacy alias for Step 10 global diagnostics."""
    return step10_current_global_diag_path(result_dir, target_id)


def step11_history_dir(result_dir: PathLike) -> Path:
    """Legacy alias for Step 10 detrend history."""
    return step10_history_dir(result_dir)


def load_detrend_preference(result_dir: PathLike, target_id: int | None = None) -> str | None:
    """Read the adopted correction mode from Step 10 detrend metadata."""
    import json as _json
    d = _as_path(result_dir)
    if target_id is not None:
        meta = step10_current_meta_path(d, target_id)
        if meta.exists():
            try:
                data = _json.loads(meta.read_text(encoding="utf-8"))
                return data.get("mode", "").lower() or None
            except Exception:
                pass
    step10_out = step10_detrend_dir(d)
    if step10_out.exists():
        for mp in sorted(step10_out.glob("result_ID*_current.json")):
            try:
                data = _json.loads(mp.read_text(encoding="utf-8"))
                return data.get("mode", "").lower() or None
            except Exception:
                continue
    return None


def list_lightcurve_csvs(result_dir: PathLike, target_id: int | None = None) -> list[Path]:
    """Return candidate light curve CSVs ordered by preferred analysis priority."""
    d = _as_path(result_dir)
    step9_out = step9_lc_dir(d)
    step10_out = step10_detrend_dir(d)
    candidates: list[Path] = []
    if target_id is not None:
        candidates.append(_step10_current_lc_path(d, target_id))
        for mode in ("global", "sysrem", "color", "offset"):
            candidates.append(step10_out / f"lightcurve_ID{int(target_id)}_{mode}.csv")
        candidates.append(step9_out / f"lightcurve_combined_ID{int(target_id)}_raw.csv")
        candidates.append(step9_out / f"lightcurve_ID{int(target_id)}_raw.csv")
        candidates.extend(sorted(step10_out.glob(f"lc_{int(target_id)}_*.csv"), reverse=True))
        candidates.append(step10_out / f"lc_{int(target_id)}.csv")
        candidates.extend(sorted(step9_out.glob(f"lc_{int(target_id)}_*.csv"), reverse=True))
        candidates.append(step9_out / f"lc_{int(target_id)}.csv")
    else:
        if step10_out.exists():
            candidates.extend(sorted(step10_out.glob("lightcurve_ID*_current.csv"), reverse=True))
            for mode in ("global", "sysrem", "color", "offset"):
                candidates.extend(sorted(step10_out.glob(f"lightcurve_ID*_{mode}.csv"), reverse=True))
        if step9_out.exists():
            candidates.extend(sorted(step9_out.glob("lightcurve_combined_ID*_raw.csv"), reverse=True))
            candidates.extend(sorted(step9_out.glob("lightcurve_ID*_raw.csv"), reverse=True))
    for base_dir in (step10_out, step9_out):
        if base_dir.exists():
            for f in sorted(base_dir.glob("lightcurve_*.csv"), reverse=True):
                if f not in candidates:
                    candidates.append(f)

    existing: list[Path] = []
    seen: set[Path] = set()
    for path in candidates:
        try:
            key = path.resolve()
        except OSError:
            continue
        if _is_usable_lightcurve_csv(path) and key not in seen:
            seen.add(key)
            existing.append(path)
    return existing
