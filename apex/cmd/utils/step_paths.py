"""
CMD-mode step path helpers.
Re-exports shared paths (step1-8) and adds CMD-specific paths (step9-13).

CMD pipeline layout:
  step1_file_selection/   File selection
  step2_crop/             Image crop
  step3_sky_preview/      Sky preview QC
  step4_detection/        Source detection + frame QC
  step5_aperture/         Aperture photometry
  cmd_psf/                PSF photometry (optional)
  step6_wcs/              WCS plate solving       } shared
  step7_refbuild/         Reference catalog build  }
  step8_idmatch/          Star ID matching         }
  cmd_selection/          Master ID editor
  cmd_zeropoint/          Zeropoint calibration
  cmd_plot/               CMD diagram
  cmd_isochrone/          Isochrone fitting
"""

from __future__ import annotations
from pathlib import Path
from typing import Union

# Re-export everything from the shared base
from apex.common.utils.step_paths_base import *  # noqa: F401, F403
from apex.common.utils.step_paths_base import step_dir, _as_path

PathLike = Union[str, Path]

# ── CMD-specific step directories ─────────────────────────────────────────────

CMD_PSF_DIRNAME        = "cmd_psf"
CMD_SELECTION_DIRNAME  = "cmd_selection"
CMD_ZP_DIRNAME         = "cmd_zeropoint"
CMD_PLOT_DIRNAME       = "cmd_plot"
CMD_ISO_DIRNAME        = "cmd_isochrone"


def step6_psf_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, CMD_PSF_DIRNAME)


def step10_selection_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, CMD_SELECTION_DIRNAME)


def step11_zp_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, CMD_ZP_DIRNAME)


def step12_cmd_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, CMD_PLOT_DIRNAME)


def step13_iso_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, CMD_ISO_DIRNAME)


# ── Legacy AAPKC compatibility aliases ───────────────────────────────────────
# Old AAPKC used different directory names; these aliases allow existing
# result directories to remain readable during migration.

def legacy_step5_aperture_dir(result_dir: PathLike) -> Path:
    """AAPKC legacy: step5_aperture/ was the aperture photometry dir."""
    return step_dir(result_dir, "step5_aperture")


def legacy_step6_psf_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, "step6_psf")


def legacy_step7_wcs_dir(result_dir: PathLike) -> Path:
    """AAPKC legacy WCS directory (step5_wcs/ or step7_wcs/)."""
    d = _as_path(result_dir)
    for name in ("step7_wcs", "step5_wcs", "step6_wcs"):
        c = d / name
        if c.exists():
            return c
    return d / "step6_wcs"  # canonical APEX name


def legacy_step8_refbuild_dir(result_dir: PathLike) -> Path:
    d = _as_path(result_dir)
    for name in ("step8_refbuild", "step6_refbuild", "step7_refbuild"):
        c = d / name
        if c.exists():
            return c
    return d / "step7_refbuild"


def legacy_step9_idmatch_dir(result_dir: PathLike) -> Path:
    d = _as_path(result_dir)
    for name in ("step9_idmatch", "step7_idmatch", "step8_idmatch"):
        c = d / name
        if c.exists():
            return c
    return d / "step8_idmatch"


def legacy_step10_selection_dir(result_dir: PathLike) -> Path:
    d = _as_path(result_dir)
    for name in ("step10_selection", "step8_selection", "cmd_selection"):
        c = d / name
        if c.exists():
            return c
    return d / "cmd_selection"


# step9_dir: legacy AAPKC forced-phot / step9_photometry directory alias
def step9_dir(result_dir: PathLike) -> Path:
    """Legacy AAPKC step9 photometry directory (step9_photometry/)."""
    d = _as_path(result_dir)
    for name in ("step9_photometry", "step5_aperture"):
        c = d / name
        if c.exists():
            return c
    return d / "step5_aperture"
