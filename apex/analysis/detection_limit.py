"""Per-frame predicted detection limit (m50) and depth QC.

Physical model
--------------
Detection of a point source is governed by its peak-pixel S/N: injection of
artificial stars into 7 real cluster frames (sky sigma 5-58 e-/px, FWHM
5.2-9.0 px) shows every frame's completeness curve collapses onto a single
erf law in peak S/N with S/N_50 = PEAK_SN50_DETECTION (4.05 +/- 0.18,
consistent with the 3.2-sigma matched-filter threshold). The 50% depth of a
frame is therefore predictable from two numbers measured on the frame itself:

    F50_e = S/N_50 * sigma_e / p_peak            [total electrons]
    m50   = ZMAG - 2.5 log10(F50_e / exptime)    [instrumental mag scale]

where sigma_e is the background noise per pixel in electrons (sep.Background
rms x gain) and p_peak is the PSF peak fraction (peak-pixel flux / total
flux). Over the 7 calibration runs this predicts the injection-measured m50
with residual RMS ~0.05 mag.

The same injection-derived completeness curve also predicts the *number* of
real catalog stars detected in a frame to ~6% (NGC 6811 R 480 s: predicted
1,629 vs 1,728 actually flagged detected), so the injection-calibrated m50
and the frame's empirical detection rolloff are directly comparable — the
basis of the QC gate in :func:`frame_depth_qc`.

Derivation and calibration data:
validation/paper/논문작업/COMPLETENESS_REALFRAME_INVESTIGATION.md
(runs: validation/paper/data_realframe_*/artificial_star/benchmark_run/).

This module is Qt-free and depends only on numpy (analysis layer).
"""
from __future__ import annotations

from typing import Optional

import numpy as np

from apex.utils.constants import (
    DEPTH_QC_TOLERANCE_MAG,
    INSTRUMENTAL_ZMAG,
    PEAK_SN50_DETECTION,
)

__all__ = [
    "predict_frame_m50",
    "peak_fraction_from_psf",
    "estimate_peak_fraction_from_stars",
    "detection_fraction_rolloff",
    "frame_depth_qc",
]


def predict_frame_m50(
    sky_sigma_e: float,
    p_peak: float,
    *,
    exptime_s: float = 1.0,
    sn50: float = PEAK_SN50_DETECTION,
    zmag: float = INSTRUMENTAL_ZMAG,
) -> float:
    """Predicted 50%-completeness magnitude of a frame.

    Parameters
    ----------
    sky_sigma_e : background noise per pixel in electrons
        (median sep.Background rms in ADU x gain).
    p_peak : PSF peak fraction — flux in the brightest pixel divided by the
        total flux of a point source (0 < p_peak <= 1).
    exptime_s : exposure time. With the default 1.0 the result is on the
        total-electron scale (zmag - 2.5 log10(F50_e)), the scale of the
        injection benchmark; pass the frame exposure time to get the
        pipeline's count-rate ``mag_inst`` scale.
    sn50 : peak-pixel S/N at 50% completeness (calibrated constant).
    zmag : instrumental zeropoint.

    Returns NaN when any input is missing or non-physical.
    """
    sky_sigma_e = float(sky_sigma_e)
    p_peak = float(p_peak)
    exptime_s = float(exptime_s)
    if not (
        np.isfinite(sky_sigma_e) and sky_sigma_e > 0
        and np.isfinite(p_peak) and 0 < p_peak <= 1
        and np.isfinite(exptime_s) and exptime_s > 0
    ):
        return float("nan")
    f50_e = float(sn50) * sky_sigma_e / p_peak
    return float(zmag) - 2.5 * float(np.log10(f50_e / exptime_s))


def peak_fraction_from_psf(psf: np.ndarray) -> float:
    """Peak fraction of a PSF stamp: max pixel / total (after normalization).

    NaN when the stamp is empty, non-finite, or has non-positive total flux.
    """
    psf = np.asarray(psf, dtype=float)
    if psf.size == 0 or not np.all(np.isfinite(psf)):
        return float("nan")
    total = float(psf.sum())
    if total <= 0:
        return float("nan")
    return float(psf.max() / total)


def estimate_peak_fraction_from_stars(
    peak_e: np.ndarray,
    flux_e: np.ndarray,
    mask: Optional[np.ndarray] = None,
) -> tuple[float, int]:
    """Empirical PSF peak fraction from measured stars.

    ``peak_e`` is the background-subtracted peak-pixel signal in electrons and
    ``flux_e`` the aperture-corrected total flux in electrons, per star.
    The estimate is the median of ``peak_e / flux_e`` over the stars selected
    by ``mask`` (all stars when None) with finite, physical ratios in (0, 1].

    Returns ``(p_peak, n_used)``; ``(nan, 0)`` when no star qualifies.
    """
    peak_e = np.asarray(peak_e, dtype=float)
    flux_e = np.asarray(flux_e, dtype=float)
    sel = np.ones(peak_e.shape, dtype=bool) if mask is None else np.asarray(mask, dtype=bool)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = peak_e / flux_e
    good = sel & np.isfinite(ratio) & (ratio > 0) & (ratio <= 1.0) & (flux_e > 0)
    n_used = int(good.sum())
    if n_used == 0:
        return float("nan"), 0
    return float(np.median(ratio[good])), n_used


def detection_fraction_rolloff(
    mag: np.ndarray,
    detected: np.ndarray,
    *,
    bin_width: float = 0.25,
    min_per_bin: int = 15,
    level: float = 0.5,
) -> float:
    """Magnitude at which the detection fraction falls through ``level``.

    Bins ``mag`` (any magnitude scale), takes the detected fraction per bin
    (bins with fewer than ``min_per_bin`` entries are skipped), and linearly
    interpolates the first falling crossing of ``level`` scanning faintward
    **from the bin of maximum completeness** — saturated or shape-rejected
    stars can depress the detected fraction at the bright end, and starting
    at the maximum keeps such a dip from faking an early rolloff. This is
    the same estimator applied to the injection catalogs in the calibration
    set, so its result is directly comparable to :func:`predict_frame_m50`.

    Returns NaN when there is no crossing (e.g. the catalog does not reach
    the frame's limit, or too few stars per bin).
    """
    mag = np.asarray(mag, dtype=float)
    detected = np.asarray(detected, dtype=bool)
    ok = np.isfinite(mag)
    mag, detected = mag[ok], detected[ok]
    if mag.size == 0:
        return float("nan")
    edges = np.arange(
        np.floor(mag.min() / bin_width) * bin_width,
        mag.max() + bin_width,
        bin_width,
    )
    centers, fractions = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (mag >= lo) & (mag < hi)
        if int(m.sum()) >= int(min_per_bin):
            centers.append(0.5 * (lo + hi))
            fractions.append(float(detected[m].mean()))
    if not centers:
        return float("nan")
    start = int(np.argmax(fractions))
    for i in range(start, len(centers) - 1):
        c0, c1 = fractions[i], fractions[i + 1]
        if c0 >= level >= c1:
            den = c1 - c0
            f = (level - c0) / den if den else 0.0
            return float(centers[i] + f * (centers[i + 1] - centers[i]))
    return float("nan")


def frame_depth_qc(
    *,
    sky_sigma_adu: float,
    gain_e_per_adu: float,
    exptime_s: float,
    peak_e: np.ndarray,
    flux_e: np.ndarray,
    peak_star_mask: np.ndarray,
    mag_inst: np.ndarray,
    detected: np.ndarray,
    rolloff_mask: Optional[np.ndarray] = None,
    tolerance_mag: float = DEPTH_QC_TOLERANCE_MAG,
    sn50: float = PEAK_SN50_DETECTION,
    zmag: float = INSTRUMENTAL_ZMAG,
) -> dict:
    """Frame depth QC: predicted m50 vs the observed detection rolloff.

    ``peak_e``/``flux_e``/``peak_star_mask`` feed the per-frame PSF peak
    fraction (bright clean stars only); ``mag_inst``/``detected`` (optionally
    restricted by ``rolloff_mask``) feed the observed 50% detection rolloff
    over the master catalog. A depth anomaly — |predicted - observed| above
    ``tolerance_mag`` — points at focus drift, clouds, tracking, or
    calibration defects (flag ``depth_shallow``), or at a censored/mismatched
    comparison (flag ``depth_deep``).

    Caveats (documented, acceptable at a 0.5-mag gate): the observed rolloff
    bins *measured* magnitudes, so it carries a small Eddington bias near the
    limit, and stars whose forced flux is non-positive cannot enter the
    denominator; both effects are well below the tolerance in the
    calibration set.

    Returns a flat dict of QC columns (NaN/"" when not computable):
    sky_sigma_adu, sky_sigma_e, p_peak_frame, n_peak_stars, predicted_m50,
    observed_m50, depth_delta_mag, depth_qc_flag.
    """
    gain = float(gain_e_per_adu)
    sky_sigma_adu = float(sky_sigma_adu)
    sky_sigma_e = sky_sigma_adu * gain if (
        np.isfinite(sky_sigma_adu) and np.isfinite(gain) and gain > 0
    ) else float("nan")

    peak_e = np.asarray(peak_e, dtype=float)
    flux_e = np.asarray(flux_e, dtype=float)
    p_peak, n_peak = estimate_peak_fraction_from_stars(peak_e, flux_e, peak_star_mask)

    predicted = predict_frame_m50(
        sky_sigma_e, p_peak, exptime_s=exptime_s, sn50=sn50, zmag=zmag
    )

    mag_inst = np.asarray(mag_inst, dtype=float)
    detected = np.asarray(detected, dtype=bool)
    if rolloff_mask is not None:
        keep = np.asarray(rolloff_mask, dtype=bool)
        mag_inst, detected = mag_inst[keep], detected[keep]
    observed = detection_fraction_rolloff(mag_inst, detected)

    delta = float("nan")
    flag = ""
    if np.isfinite(predicted) and np.isfinite(observed):
        delta = observed - predicted
        tol = float(tolerance_mag)
        if delta < -tol:
            flag = "depth_shallow"
        elif delta > tol:
            flag = "depth_deep"
        else:
            flag = "ok"

    return {
        "sky_sigma_adu": sky_sigma_adu,
        "sky_sigma_e": sky_sigma_e,
        "p_peak_frame": p_peak,
        "n_peak_stars": n_peak,
        "predicted_m50": predicted,
        "observed_m50": observed,
        "depth_delta_mag": delta,
        "depth_qc_flag": flag,
    }
