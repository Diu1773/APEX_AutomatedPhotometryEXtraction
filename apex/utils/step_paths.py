"""
Step path helpers shared across all APEX modes.

Steps 1-8 use the same canonical directory names in every pipeline variant.
Steps 9+ diverge between cmd/ and lightcurve/ and are defined in each mode's
own step_paths.py.

Canonical layout:
  step1_file_selection/   File selection + FITS header scan
  step2_crop/             Crop region + cropped images
  step3_sky_preview/      Sky preview QC metadata
  step4_detection/        Source detection + frame QC
  step5_aperture/         Aperture photometry
  step6_wcs/              WCS plate solving
  step7_refbuild/         Reference catalog build
  step8_idmatch/          Star ID matching
"""

from __future__ import annotations
from pathlib import Path
from typing import Union

PathLike = Union[str, Path]

STEP1_DIRNAME = "step1_file_selection"
STEP2_DIRNAME = "step2_crop"
STEP3_DIRNAME = "step3_sky_preview"
STEP4_DIRNAME = "step4_detection"


def _as_path(p: PathLike) -> Path:
    return Path(p) if not isinstance(p, Path) else p


def step_dir(result_dir: PathLike, dirname: str) -> Path:
    return _as_path(result_dir) / dirname


def step1_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, STEP1_DIRNAME)


def step2_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, STEP2_DIRNAME)


def step2_cropped_dir(result_dir: PathLike) -> Path:
    return step2_dir(result_dir) / "cropped"


def crop_rect_path(result_dir: PathLike) -> Path:
    return step2_dir(result_dir) / "crop_rect.json"


def crop_is_active(result_dir: PathLike) -> bool:
    return crop_rect_path(result_dir).exists()


def step3_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, STEP3_DIRNAME)


def step4_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, STEP4_DIRNAME)


# ── Steps 5-8: shared pipeline steps ─────────────────────────────────────────

STEP5_APERTURE_DIRNAME = "step5_aperture"
STEP6_WCS_DIRNAME      = "step6_wcs"
STEP7_REFBUILD_DIRNAME = "step7_refbuild"
STEP8_IDMATCH_DIRNAME  = "step8_idmatch"
TOOL_EXTINCTION_DIRNAME = "tool_extinction"


def step5_aperture_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, STEP5_APERTURE_DIRNAME)


def step6_wcs_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, STEP6_WCS_DIRNAME)


def step7_refbuild_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, STEP7_REFBUILD_DIRNAME)


def step8_idmatch_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, STEP8_IDMATCH_DIRNAME)


def tool_extinction_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, TOOL_EXTINCTION_DIRNAME)
