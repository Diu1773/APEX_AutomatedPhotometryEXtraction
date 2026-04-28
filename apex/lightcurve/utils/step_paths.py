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
from apex.common.utils.step_paths_base import *  # noqa: F401, F403
from apex.common.utils.step_paths_base import step_dir, _as_path

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


# ── Legacy AAPKL compatibility aliases ───────────────────────────────────────

def legacy_step5_photometry_dir(result_dir: PathLike) -> Path:
    """AAPKL legacy: step5_photometry/ was aperture photometry output."""
    d = _as_path(result_dir)
    for name in ("step5_photometry", "step5_aperture"):
        c = d / name
        if c.exists():
            return c
    return d / "step5_aperture"


def legacy_step9_selection_dir(result_dir: PathLike) -> Path:
    d = _as_path(result_dir)
    for name in ("step9_selection", "step8_selection", "lc_selection"):
        c = d / name
        if c.exists():
            return c
    return d / "lc_selection"


def legacy_step10_lc_dir(result_dir: PathLike) -> Path:
    d = _as_path(result_dir)
    for name in ("step10_lightcurve", "lc_lightcurve"):
        c = d / name
        if c.exists():
            return c
    return d / "lc_lightcurve"


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


# step5_photometry_dir: AAPKL canonical name alias
def step5_photometry_dir(result_dir: PathLike) -> Path:
    """AAPKL canonical aperture photometry directory alias."""
    return step5_aperture_dir(result_dir)
