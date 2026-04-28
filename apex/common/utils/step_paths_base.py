"""
Step path helpers shared across all APEX modes.

Only steps 1-4 (file selection, crop, sky preview, detection) have canonical
directory names that are the same in every pipeline variant.  Steps 5+ diverge
between cmd/ and lightcurve/ and are defined in each mode's own step_paths.py.
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
