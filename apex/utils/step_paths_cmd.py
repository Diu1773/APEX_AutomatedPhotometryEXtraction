"""
CMD-mode step path helpers.
Re-exports shared paths (step1-7) and adds CMD-specific paths (step8-12).

CMD pipeline layout:
  step1_file_selection/   File selection
  step2_crop/             Image crop
  step3_sky_preview/      Sky preview QC
  step4_detection/        Source detection + frame QC
  step5_wcs/              WCS plate solving       } shared
  step6_refbuild/         Reference catalog build  }
  step7_forced_phot/       Forced aperture photometry }
  cmd_psf/                Step 8 PSF photometry (optional)
  cmd_selection/          Step 9 Master ID editor
  cmd_zeropoint/          Step 10 zeropoint calibration
  cmd_plot/               Step 11 CMD diagram
  cmd_isochrone/          Step 12 isochrone fitting
"""

from __future__ import annotations
from pathlib import Path
from typing import Union

# Re-export everything from the shared base
from apex.utils.step_paths import *  # noqa: F401, F403
from apex.utils.step_paths import step_dir, _as_path

PathLike = Union[str, Path]

# ── CMD-specific step directories ─────────────────────────────────────────────

CMD_PSF_DIRNAME        = "cmd_psf"
CMD_SELECTION_DIRNAME  = "cmd_selection"
CMD_ZP_DIRNAME         = "cmd_zeropoint"
CMD_PLOT_DIRNAME       = "cmd_plot"
CMD_ISO_DIRNAME        = "cmd_isochrone"


def step8_psf_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, CMD_PSF_DIRNAME)


def step9_selection_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, CMD_SELECTION_DIRNAME)


def step10_zp_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, CMD_ZP_DIRNAME)


def step11_cmd_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, CMD_PLOT_DIRNAME)


def step12_iso_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, CMD_ISO_DIRNAME)


def step6_psf_dir(result_dir: PathLike) -> Path:
    """Legacy alias for Step 8 PSF output."""
    return step8_psf_dir(result_dir)


def step10_selection_dir(result_dir: PathLike) -> Path:
    """Legacy alias for Step 9 Master ID Editor output."""
    return step9_selection_dir(result_dir)


def step11_zp_dir(result_dir: PathLike) -> Path:
    """Legacy alias for Step 10 zeropoint output."""
    return step10_zp_dir(result_dir)


def step12_cmd_dir(result_dir: PathLike) -> Path:
    """Legacy alias for Step 11 CMD plot output."""
    return step11_cmd_dir(result_dir)


def step13_iso_dir(result_dir: PathLike) -> Path:
    """Legacy alias for Step 12 isochrone output."""
    return step12_iso_dir(result_dir)

