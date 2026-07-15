"""Quality metrics and iteration policy for CPU PSF photometry."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import IntFlag
from typing import Callable

import numpy as np


class PSFFitFlag(IntFlag):
    """PSF fit flags aligned with ``photutils.psf.PSFPhotometry``."""

    NONE = 0
    INCOMPLETE_REGION = 1
    OUTSIDE_IMAGE = 2
    NONPOSITIVE_FLUX = 4
    NONCONVERGENCE = 8
    MISSING_COVARIANCE = 16
    NEAR_BOUND = 32
    NO_OVERLAP = 64
    FULLY_MASKED = 128
    TOO_FEW_PIXELS = 256
    NONFINITE_POSITION = 512
    NONFINITE_FLUX = 1024
    NONFINITE_LOCAL_BACKGROUND = 2048
    SATURATED = 4096
    CROWDING_UNRELIABLE = 8192


@dataclass(frozen=True)
class IterationDecision:
    """Decision made after detecting candidates in an actual residual image."""

    stop_now: bool
    stop_after_refit: bool
    reason: str
    candidate_fraction: float


@dataclass(frozen=True)
class IterationSnapshot:
    """Serializable summary of one completed fit/subtraction pass."""

    iteration: int
    n_fit: int
    n_candidates_raw: int
    n_candidates_unique: int
    n_candidates_accepted: int
    residual_std: float
    median_qfit: float
    median_reduced_chi2: float
    elapsed_s: float
    stop_reason: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class PSFQualityMetrics:
    qfit: np.ndarray
    cfit: np.ndarray
    reduced_chi2: np.ndarray
    flux_err: np.ndarray
    flags: np.ndarray
    n_pixels_fit: np.ndarray


@dataclass(frozen=True)
class PSFFrameAssessment:
    """Conservative post-fit quality grade for one PSF frame."""

    status: str
    reasons: tuple[str, ...]
    score: float
    clean_fraction: float
    fit_failure_fraction: float
    crowding_unreliable_fraction: float

    def to_dict(self) -> dict:
        payload = asdict(self)
        payload["reasons"] = ",".join(self.reasons)
        return payload


def assess_psf_frame_quality(
    *,
    n_sources: int,
    n_good: int,
    n_crowding_unreliable: int,
    median_qfit_noise_ratio: float,
    epsf_n_selected: int,
    epsf_median_contamination: float,
    frame_fwhm_px: float,
    frame_fwhm_max_px: float,
) -> PSFFrameAssessment:
    """Grade a completed PSF fit without treating aperture photometry as truth.

    The grade combines fit survival, noise-normalized residual shape, ePSF
    reference quality, and the configured seeing envelope.  ``REVIEW`` is
    intentionally common for difficult crowded frames; only clearly unusable
    outputs are marked ``FAIL``.
    """

    total = max(0, int(n_sources))
    good = min(total, max(0, int(n_good)))
    crowding = min(total, max(0, int(n_crowding_unreliable)))
    if total > 0:
        clean_fraction = good / total
        failure_fraction = 1.0 - clean_fraction
        crowding_fraction = crowding / total
    else:
        clean_fraction = 0.0
        failure_fraction = 1.0
        crowding_fraction = 0.0

    fail: list[str] = []
    review: list[str] = []
    score = 100.0

    if total <= 0 or good <= 0:
        fail.append("no_clean_psf_sources")
        score -= 70.0
    elif clean_fraction < 0.70:
        fail.append("low_clean_fraction")
        score -= 40.0
    elif clean_fraction < 0.90:
        review.append("clean_fraction_warning")
        score -= 15.0

    qratio = float(median_qfit_noise_ratio)
    if np.isfinite(qratio):
        if qratio > 3.0:
            fail.append("high_qfit_noise_ratio")
            score -= 35.0
        elif qratio > 1.5:
            review.append("qfit_noise_warning")
            score -= 15.0
    else:
        review.append("qfit_noise_missing")
        score -= 8.0

    if int(epsf_n_selected) < 3:
        fail.append("too_few_epsf_stars")
        score -= 35.0
    elif int(epsf_n_selected) < 8:
        review.append("few_epsf_stars")
        score -= 12.0

    contamination = float(epsf_median_contamination)
    if np.isfinite(contamination):
        if contamination > 0.50:
            fail.append("high_epsf_contamination")
            score -= 35.0
        elif contamination > 0.25:
            review.append("epsf_contamination_warning")
            score -= 15.0

    if crowding_fraction > 0.15:
        fail.append("many_crowding_unreliable")
        score -= 30.0
    elif crowding_fraction > 0.03:
        review.append("crowding_warning")
        score -= 12.0

    fwhm = float(frame_fwhm_px)
    fwhm_max = float(frame_fwhm_max_px)
    if np.isfinite(fwhm) and np.isfinite(fwhm_max) and fwhm > 0 and fwhm_max > 0:
        seeing_ratio = fwhm / fwhm_max
        if seeing_ratio > 1.25:
            fail.append("fwhm_far_above_config")
            score -= 35.0
        elif seeing_ratio > 1.0:
            review.append("fwhm_above_config")
            score -= 15.0

    reasons = tuple(dict.fromkeys(fail + review if fail else review))
    status = "FAIL" if fail else ("REVIEW" if review else "PASS")
    return PSFFrameAssessment(
        status=status,
        reasons=reasons,
        score=float(np.clip(score, 0.0, 100.0)),
        clean_fraction=float(clean_fraction),
        fit_failure_fraction=float(failure_fraction),
        crowding_unreliable_fraction=float(crowding_fraction),
    )


def qfit_noise_diagnostics(
    qfit: np.ndarray | float,
    n_pixels_fit: np.ndarray | float,
    snr: np.ndarray | float,
    psf_nea_px: np.ndarray | float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return expected noise-only qfit and the observed/expected ratio.

    Raw qfit grows with the number of fitted pixels and falls with source SNR.
    This normalization makes a single threshold meaningful across fit-window
    sizes and source brightnesses.
    """

    qfit_arr, n_pixels_arr, snr_arr, nea_arr = np.broadcast_arrays(
        np.asarray(qfit, dtype=float),
        np.asarray(n_pixels_fit, dtype=float),
        np.asarray(snr, dtype=float),
        np.asarray(psf_nea_px, dtype=float),
    )
    expected = np.full(qfit_arr.shape, np.nan, dtype=float)
    valid_model = (
        np.isfinite(n_pixels_arr)
        & (n_pixels_arr > 0)
        & np.isfinite(snr_arr)
        & (snr_arr > 0)
        & np.isfinite(nea_arr)
        & (nea_arr > 0)
    )
    np.divide(
        np.sqrt(2.0 / np.pi) * n_pixels_arr,
        snr_arr * np.sqrt(nea_arr),
        out=expected,
        where=valid_model,
    )
    ratio = np.full(qfit_arr.shape, np.nan, dtype=float)
    np.divide(
        qfit_arr,
        expected,
        out=ratio,
        where=np.isfinite(qfit_arr) & np.isfinite(expected) & (expected > 0),
    )
    return expected, ratio


def fit_parameters_changed(
    x_old: float,
    y_old: float,
    flux_old: float,
    x_new: float,
    y_new: float,
    flux_new: float,
    *,
    flux_fraction: float,
    position_pixels: float = 0.01,
) -> bool:
    """Return whether any fitted parameter remains outside its tolerance."""

    dflux = abs(float(flux_new) - float(flux_old)) / max(abs(float(flux_old)), 1e-10)
    return bool(
        dflux > max(0.0, float(flux_fraction))
        or abs(float(x_new) - float(x_old)) > max(0.0, float(position_pixels))
        or abs(float(y_new) - float(y_old)) > max(0.0, float(position_pixels))
    )


def decide_residual_iteration(
    *,
    n_candidates_raw: int,
    n_candidates_unique: int,
    n_candidates_accepted: int,
    n_current: int,
    convergence_fraction: float,
) -> IterationDecision:
    """Choose a stop condition using pre-cap candidate counts.

    The convergence fraction is intentionally evaluated before the CPU cap is
    applied. Otherwise a small cap can force every residual pass to converge.
    """

    raw = max(0, int(n_candidates_raw))
    unique = max(0, int(n_candidates_unique))
    accepted = max(0, int(n_candidates_accepted))
    denominator = max(1, int(n_current))
    candidate_fraction = unique / denominator

    if raw == 0:
        return IterationDecision(True, False, "no_candidates", candidate_fraction)
    if unique == 0:
        return IterationDecision(True, False, "duplicates_only", candidate_fraction)
    if accepted == 0:
        return IterationDecision(True, False, "no_accepted_candidates", candidate_fraction)
    if convergence_fraction > 0 and candidate_fraction <= convergence_fraction:
        return IterationDecision(False, True, "candidate_fraction", candidate_fraction)
    return IterationDecision(False, False, "continue", candidate_fraction)


def measure_psf_fit_quality(
    data: np.ndarray,
    model: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    flux: np.ndarray,
    eval_psf: Callable,
    *,
    fit_shape: int,
    background_rms: float,
    gain: float,
    fit_ok: np.ndarray | None = None,
    initial_xy: np.ndarray | None = None,
    xy_bound: float | None = None,
) -> PSFQualityMetrics:
    """Measure Photutils-compatible fit diagnostics from the final residual.

    ``qfit`` is ``sum(abs(data-model)) / flux`` over the fit footprint and
    ``cfit`` is the central residual divided by flux. Reduced chi-square and
    the formal flux error use a source-plus-background variance model.
    """

    # Keep full-frame arrays in their native dtype.  Quality metrics only need
    # small source footprints, where conversion to float64 is inexpensive.
    image = np.asarray(data)
    fitted_model = np.asarray(model)
    x_arr = np.asarray(x, dtype=float)
    y_arr = np.asarray(y, dtype=float)
    flux_arr = np.asarray(flux, dtype=float)
    n_sources = len(flux_arr)
    if image.shape != fitted_model.shape:
        raise ValueError("data and model must have the same shape")
    if len(x_arr) != n_sources or len(y_arr) != n_sources:
        raise ValueError("x, y, and flux must have the same length")

    qfit = np.full(n_sources, np.nan, dtype=float)
    cfit = np.full(n_sources, np.nan, dtype=float)
    reduced_chi2 = np.full(n_sources, np.nan, dtype=float)
    flux_err = np.full(n_sources, np.nan, dtype=float)
    flags = np.zeros(n_sources, dtype=np.int32)
    n_pixels_fit = np.zeros(n_sources, dtype=np.int32)
    half = max(1, int(fit_shape) // 2)
    h, w = image.shape
    bkg_var = max(float(background_rms), 1e-6) ** 2
    gain_safe = max(float(gain), 1e-6)
    ok_arr = (
        np.ones(n_sources, dtype=bool)
        if fit_ok is None
        else np.asarray(fit_ok, dtype=bool)
    )
    anchor = None if initial_xy is None else np.asarray(initial_xy, dtype=float)
    anchor_valid = anchor is not None and anchor.shape == (n_sources, 2)

    for index in range(n_sources):
        xi = x_arr[index]
        yi = y_arr[index]
        fi = flux_arr[index]
        flag = PSFFitFlag.NONE
        if not np.isfinite(xi) or not np.isfinite(yi):
            flags[index] = int(PSFFitFlag.NONFINITE_POSITION)
            continue
        if not np.isfinite(fi):
            flags[index] = int(PSFFitFlag.NONFINITE_FLUX)
            continue
        if fi <= 0:
            flag |= PSFFitFlag.NONPOSITIVE_FLUX
        if index >= len(ok_arr) or not ok_arr[index]:
            flag |= PSFFitFlag.NONCONVERGENCE

        x_center = int(round(float(xi)))
        y_center = int(round(float(yi)))
        if x_center < 0 or x_center >= w or y_center < 0 or y_center >= h:
            flags[index] = int(flag | PSFFitFlag.OUTSIDE_IMAGE | PSFFitFlag.NO_OVERLAP)
            continue

        y0 = max(0, y_center - half)
        y1 = min(h, y_center + half + 1)
        x0 = max(0, x_center - half)
        x1 = min(w, x_center + half + 1)
        if (y1 - y0) < (2 * half + 1) or (x1 - x0) < (2 * half + 1):
            flag |= PSFFitFlag.INCOMPLETE_REGION

        yy, xx = np.mgrid[y0:y1, x0:x1]
        data_patch = np.asarray(image[y0:y1, x0:x1], dtype=float)
        model_patch = np.asarray(fitted_model[y0:y1, x0:x1], dtype=float)
        res_patch = data_patch - model_patch
        psf_patch = np.asarray(eval_psf(xx - xi, yy - yi), dtype=float)
        valid = np.isfinite(res_patch) & np.isfinite(model_patch) & np.isfinite(psf_patch)
        n_valid = int(np.sum(valid))
        n_pixels_fit[index] = n_valid
        if n_valid == 0:
            flags[index] = int(flag | PSFFitFlag.FULLY_MASKED | PSFFitFlag.NO_OVERLAP)
            continue
        if n_valid <= 4:
            flag |= PSFFitFlag.TOO_FEW_PIXELS

        flux_norm = max(abs(float(fi)), 1e-20)
        qfit[index] = float(np.sum(np.abs(res_patch[valid]))) / flux_norm
        center_residual = float(image[y_center, x_center]) - float(
            fitted_model[y_center, x_center]
        )
        if np.isfinite(center_residual):
            cfit[index] = float(center_residual) / flux_norm

        variance = bkg_var + np.clip(model_patch, 0.0, None) / gain_safe
        valid_var = valid & np.isfinite(variance) & (variance > 0)
        dof = int(np.sum(valid_var)) - 3
        if dof > 0:
            reduced_chi2[index] = float(
                np.sum((res_patch[valid_var] ** 2) / variance[valid_var]) / dof
            )
        information = float(np.sum((psf_patch[valid_var] ** 2) / variance[valid_var]))
        if information > 0 and np.isfinite(information):
            flux_err[index] = 1.0 / np.sqrt(information)
        else:
            flag |= PSFFitFlag.MISSING_COVARIANCE

        if anchor_valid and xy_bound is not None and xy_bound > 0:
            dx = abs(float(xi) - float(anchor[index, 0]))
            dy = abs(float(yi) - float(anchor[index, 1]))
            if max(dx, dy) >= 0.95 * float(xy_bound):
                flag |= PSFFitFlag.NEAR_BOUND
        flags[index] = int(flag)

    return PSFQualityMetrics(
        qfit=qfit,
        cfit=cfit,
        reduced_chi2=reduced_chi2,
        flux_err=flux_err,
        flags=flags,
        n_pixels_fit=n_pixels_fit,
    )
