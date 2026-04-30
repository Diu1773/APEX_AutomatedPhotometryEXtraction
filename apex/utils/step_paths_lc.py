"""
Lightcurve-mode step path helpers.
Re-exports shared paths (step1-8) and adds LC-specific paths (step9-12).

LC pipeline layout:
  step1_file_selection/   File selection (multi-night)
  step2_crop/             Image crop
  step3_sky_preview/      Sky preview QC
  step4_detection/        Source detection + frame QC
  step5_aperture/         Aperture photometry (forced)
  step6_wcs/              WCS plate solving       } shared
  step7_refbuild/         Reference catalog build  }
  step8_idmatch/          Star ID matching         }
  lc_selection/           Target/comparison selection
  lc_lightcurve/          Light curve builder
  lc_detrend/             Detrend & night merge
  lc_period/              Period analysis
"""

from __future__ import annotations
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


def step9_selection_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, LC_SELECTION_DIRNAME)


def step10_lc_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, LC_LC_DIRNAME)


def step11_detrend_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, LC_DETREND_DIRNAME)


def step12_period_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, LC_PERIOD_DIRNAME)


def find_best_lightcurve_csv(result_dir: PathLike, star_id: int) -> Path | None:
    """Find the best detrended lightcurve CSV for a given star ID."""
    lc_dir = step10_lc_dir(result_dir)
    detrend_dir = step11_detrend_dir(result_dir)
    for base_dir in (detrend_dir, lc_dir):
        for pattern in (f"lc_{star_id}_*.csv", f"lc_{star_id}.csv"):
            matches = sorted(base_dir.glob(pattern))
            if matches:
                return matches[-1]
    return None


def step5_photometry_dir(result_dir: PathLike) -> Path:
    """LC aperture photometry directory alias."""
    return step5_aperture_dir(result_dir)


# ── LC short-name aliases used by analysis modules ───────────────────────────

def step10_dir(result_dir: PathLike) -> Path:
    """Alias for step10_lc_dir (lc_lightcurve/)."""
    return step10_lc_dir(result_dir)


def step11_dir(result_dir: PathLike) -> Path:
    """Alias for step11_detrend_dir (lc_detrend/)."""
    return step11_detrend_dir(result_dir)


def step11_current_meta_path(result_dir: PathLike, target_id: int) -> Path:
    """Path to the current detrend result metadata JSON for a target star."""
    return step11_dir(result_dir) / f"result_ID{int(target_id)}_current.json"


def step11_current_lc_path(result_dir: PathLike, target_id: int) -> Path:
    return step11_dir(result_dir) / f"lightcurve_ID{int(target_id)}_current.csv"


def _step11_current_lc_path(result_dir: PathLike, target_id: int) -> Path:
    return step11_current_lc_path(result_dir, target_id)


def step11_current_params_path(result_dir: PathLike, target_id: int) -> Path:
    return step11_dir(result_dir) / f"params_ID{int(target_id)}_current.csv"


def step11_current_summary_path(result_dir: PathLike, target_id: int) -> Path:
    return step11_dir(result_dir) / f"summary_ID{int(target_id)}_current.txt"


def step11_current_plot_path(result_dir: PathLike, target_id: int) -> Path:
    return step11_dir(result_dir) / f"plot_ID{int(target_id)}_current.png"


def step11_current_global_zp_path(result_dir: PathLike, target_id: int) -> Path:
    return step11_dir(result_dir) / f"global_zp_ID{int(target_id)}_current.csv"


def step11_current_global_mean_path(result_dir: PathLike, target_id: int) -> Path:
    return step11_dir(result_dir) / f"global_mean_ID{int(target_id)}_current.csv"


def step11_current_global_diag_path(result_dir: PathLike, target_id: int) -> Path:
    return step11_dir(result_dir) / f"global_diagnostics_ID{int(target_id)}_current.json"


def step11_history_dir(result_dir: PathLike) -> Path:
    return step11_dir(result_dir) / "_history"


def load_detrend_preference(result_dir: PathLike, target_id: int | None = None) -> str | None:
    """Read the adopted correction mode from step11 meta JSON."""
    import json as _json
    d = _as_path(result_dir)
    if target_id is not None:
        meta = step11_current_meta_path(d, target_id)
        if meta.exists():
            try:
                data = _json.loads(meta.read_text(encoding="utf-8"))
                return data.get("mode", "").lower() or None
            except Exception:
                pass
    s11 = step11_dir(d)
    if s11.exists():
        for mp in sorted(s11.glob("result_ID*_current.json")):
            try:
                data = _json.loads(mp.read_text(encoding="utf-8"))
                return data.get("mode", "").lower() or None
            except Exception:
                continue
    return None


def list_lightcurve_csvs(result_dir: PathLike, target_id: int | None = None) -> list[Path]:
    """Return candidate light curve CSVs ordered by preferred analysis priority."""
    d = _as_path(result_dir)
    step10_out = step10_dir(d)
    step11_out = step11_dir(d)
    candidates: list[Path] = []
    if target_id is not None:
        candidates.append(_step11_current_lc_path(d, target_id))
        for mode in ("global", "color", "offset"):
            candidates.append(step11_out / f"lightcurve_ID{int(target_id)}_{mode}.csv")
        candidates.append(step10_out / f"lightcurve_combined_ID{int(target_id)}_raw.csv")
        candidates.append(step10_out / f"lightcurve_ID{int(target_id)}_raw.csv")
    else:
        if step11_out.exists():
            candidates.extend(sorted(step11_out.glob("lightcurve_ID*_current.csv"), reverse=True))
            for mode in ("global", "color", "offset"):
                candidates.extend(sorted(step11_out.glob(f"lightcurve_ID*_{mode}.csv"), reverse=True))
        if step10_out.exists():
            candidates.extend(sorted(step10_out.glob("lightcurve_combined_ID*_raw.csv"), reverse=True))
            candidates.extend(sorted(step10_out.glob("lightcurve_ID*_raw.csv"), reverse=True))
    for base_dir in (step11_out, step10_out):
        if base_dir.exists():
            for f in sorted(base_dir.glob("lightcurve_*.csv"), reverse=True):
                if f not in candidates:
                    candidates.append(f)
    return candidates
