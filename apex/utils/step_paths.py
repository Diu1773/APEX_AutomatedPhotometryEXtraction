"""
Step path helpers shared across all APEX modes.

Canonical layout:
  step0_calibration/      Detector calibration (optional): masters + calibrated frames
  step1_file_selection/   File selection + FITS header scan
  step2_crop/             Crop region + cropped images
  step3_sky_preview/      Sky preview QC metadata
  step4_detection/        Source detection + frame QC
  step5_wcs/              WCS plate solving
  step6_refbuild/         Master catalog build (MasterBuild)
  step7_forced_phot/       Forced aperture photometry (master-driven)
"""

from __future__ import annotations
from pathlib import Path
from typing import Union

PathLike = Union[str, Path]

# Detector calibration is an optional off-chain pre-stage ("Step 0"); it is not
# part of the numbered step chain, but its outputs live alongside the others.
STEP0_CALIBRATION_DIRNAME = "step0_calibration"
STEP1_DIRNAME = "step1_file_selection"
STEP2_DIRNAME = "step2_crop"
STEP3_DIRNAME = "step3_sky_preview"
STEP4_DIRNAME = "step4_detection"


def _as_path(p: PathLike) -> Path:
    return Path(p) if not isinstance(p, Path) else p


def step_dir(result_dir: PathLike, dirname: str) -> Path:
    return _as_path(result_dir) / dirname


def step0_calibration_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, STEP0_CALIBRATION_DIRNAME)


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


# ── Shared pipeline steps ─────────────────────────────────────────────────────

STEP5_WCS_DIRNAME         = "step5_wcs"
STEP6_REFBUILD_DIRNAME    = "step6_refbuild"
STEP7_FORCED_PHOT_DIRNAME  = "step7_forced_phot"
TOOL_EXTINCTION_DIRNAME   = "tool_extinction"


def step5_wcs_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, STEP5_WCS_DIRNAME)


def step6_refbuild_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, STEP6_REFBUILD_DIRNAME)


def step7_forced_phot_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, STEP7_FORCED_PHOT_DIRNAME)


def tool_extinction_dir(result_dir: PathLike) -> Path:
    return step_dir(result_dir, TOOL_EXTINCTION_DIRNAME)


# ── Legacy-aware input resolution ─────────────────────────────────────────────
# The step directories were renumbered/renamed over time. Code that READS a
# pre-existing input workspace (e.g. the multi-night merger) must tolerate the
# older on-disk names so old RESULT_* folders stay consumable. Writers keep
# using the canonical helpers above — a freshly created workspace has none of
# the legacy dirs, so first_existing_dir() returns the canonical name.

def first_existing_dir(result_dir: PathLike, *dirnames: str) -> Path:
    """Return the first of `dirnames` that exists under result_dir; if none
    exist, return the first (canonical) one. Read-only callers only."""
    paths = [step_dir(result_dir, d) for d in dirnames]
    for p in paths:
        if p.exists():
            return p
    return paths[0]


def forced_phot_input_dir(result_dir: PathLike) -> Path:
    """Forced-photometry dir for READING a (possibly legacy) input workspace."""
    return first_existing_dir(result_dir, STEP7_FORCED_PHOT_DIRNAME, "step5_photometry")
