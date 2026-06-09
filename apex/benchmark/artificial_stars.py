"""Empirical-PSF extraction and artificial-star injection."""

from __future__ import annotations

from dataclasses import dataclass
import math

import numpy as np
import pandas as pd
from scipy.ndimage import shift as ndi_shift
from scipy.spatial import cKDTree


@dataclass(frozen=True)
class EmpiricalPSFResult:
    data: np.ndarray
    n_candidates: int
    n_used: int
    stamp_size: int
    selection_mode: str


def _odd_size(value: float, minimum: int = 9, maximum: int = 101) -> int:
    size = int(math.ceil(float(value)))
    size = max(minimum, min(maximum, size))
    return size if size % 2 == 1 else size + 1


def _positive_stamp(
    image: np.ndarray,
    x: float,
    y: float,
    size: int,
) -> np.ndarray | None:
    half = size // 2
    xi = int(round(x))
    yi = int(round(y))
    y0, y1 = yi - half, yi + half + 1
    x0, x1 = xi - half, xi + half + 1
    if y0 < 0 or x0 < 0 or y1 > image.shape[0] or x1 > image.shape[1]:
        return None

    stamp = np.asarray(image[y0:y1, x0:x1], dtype=float)
    if stamp.shape != (size, size) or not np.isfinite(stamp).any():
        return None

    yy, xx = np.mgrid[0:size, 0:size]
    rr = np.hypot(xx - half, yy - half)
    border = stamp[rr >= 0.72 * half]
    border = border[np.isfinite(border)]
    background = float(np.nanmedian(border)) if border.size else float(np.nanmedian(stamp))
    signal = np.where(np.isfinite(stamp), stamp - background, 0.0)
    signal[signal < 0] = 0.0
    total = float(signal.sum())
    if not np.isfinite(total) or total <= 0:
        return None

    cy = float((yy * signal).sum() / total)
    cx = float((xx * signal).sum() / total)
    aligned = ndi_shift(
        signal,
        shift=(half - cy, half - cx),
        order=3,
        mode="constant",
        cval=0.0,
        prefilter=True,
    )
    aligned = np.where(np.isfinite(aligned), aligned, 0.0)
    aligned[aligned < 0] = 0.0
    aligned_sum = float(aligned.sum())
    if aligned_sum <= 0:
        return None
    return aligned / aligned_sum


def extract_empirical_psf(
    image: np.ndarray,
    detections: pd.DataFrame,
    fwhm_px: float,
    *,
    saturation_adu: float,
    radius_fwhm: float = 4.0,
    isolation_fwhm: float = 4.0,
    max_stars: int = 40,
    min_stars: int = 3,
    allow_nonisolated_fallback: bool = False,
) -> EmpiricalPSFResult:
    """Build a normalized empirical PSF from isolated real stars.

    Detection is used only to locate candidate stars. Pixel values come from
    the unmodified observation, and the final PSF is a median stack.
    """
    if not {"x", "y"} <= set(detections.columns):
        raise ValueError("Detection catalog must contain x and y columns")
    if not np.isfinite(fwhm_px) or fwhm_px <= 0:
        raise ValueError("fwhm_px must be positive")

    xy = detections[["x", "y"]].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    finite = np.isfinite(xy).all(axis=1)
    xy = xy[finite]
    if len(xy) < min_stars:
        raise RuntimeError(f"Too few finite detections for empirical PSF: {len(xy)}")

    size = _odd_size(2.0 * radius_fwhm * fwhm_px + 1.0)
    half = size // 2
    h, w = image.shape
    edge_ok = (
        (xy[:, 0] >= half)
        & (xy[:, 0] < w - half)
        & (xy[:, 1] >= half)
        & (xy[:, 1] < h - half)
    )

    if len(xy) > 1:
        nearest, _ = cKDTree(xy).query(xy, k=2)
        isolated = nearest[:, 1] >= isolation_fwhm * fwhm_px
    else:
        isolated = np.ones(len(xy), dtype=bool)

    ix = np.clip(np.rint(xy[:, 0]).astype(int), 0, w - 1)
    iy = np.clip(np.rint(xy[:, 1]).astype(int), 0, h - 1)
    peaks = np.asarray(image[iy, ix], dtype=float)
    valid_peaks = peaks[np.isfinite(peaks) & (peaks < 0.9 * saturation_adu)]
    if valid_peaks.size:
        peak_lo = float(np.nanpercentile(valid_peaks, 55))
        peak_hi = float(np.nanpercentile(valid_peaks, 98))
        brightness_ok = np.isfinite(peaks) & (peaks >= peak_lo) & (peaks <= peak_hi)
    else:
        brightness_ok = np.zeros(len(xy), dtype=bool)

    candidate_mask = edge_ok & isolated & brightness_ok
    candidate_indices = np.flatnonzero(candidate_mask)
    selection_mode = "isolated"
    if len(candidate_indices) < min_stars:
        candidate_indices = np.flatnonzero(
            edge_ok
            & isolated
            & np.isfinite(peaks)
            & (peaks < 0.9 * saturation_adu)
        )
        selection_mode = "isolated_brightness_relaxed"
    if len(candidate_indices) < min_stars and allow_nonisolated_fallback:
        candidate_indices = np.flatnonzero(edge_ok & brightness_ok)
        selection_mode = "nonisolated_brightness_fallback"
    if len(candidate_indices) < min_stars and allow_nonisolated_fallback:
        candidate_indices = np.flatnonzero(edge_ok & np.isfinite(peaks))
        selection_mode = "nonisolated_edge_fallback"
    if len(candidate_indices) < min_stars:
        raise RuntimeError(
            f"Empirical PSF requires {min_stars} isolated stars; found "
            f"{len(candidate_indices)}. Reduce psf_isolation_fwhm only with "
            "a documented reason, or explicitly enable nonisolated fallback."
        )

    order = candidate_indices[np.argsort(peaks[candidate_indices])[::-1]]
    stamps: list[np.ndarray] = []
    for idx in order[: max(max_stars * 2, min_stars)]:
        stamp = _positive_stamp(image, xy[idx, 0], xy[idx, 1], size)
        if stamp is not None:
            stamps.append(stamp)
        if len(stamps) >= max_stars:
            break

    if len(stamps) < min_stars:
        raise RuntimeError(
            f"Empirical PSF requires at least {min_stars} usable stars; found {len(stamps)}"
        )

    psf = np.nanmedian(np.stack(stamps, axis=0), axis=0)
    psf = np.where(np.isfinite(psf), psf, 0.0)
    psf[psf < 0] = 0.0
    total = float(psf.sum())
    if total <= 0:
        raise RuntimeError("Empirical PSF stack has zero flux")
    psf /= total
    return EmpiricalPSFResult(
        data=psf.astype(np.float64),
        n_candidates=int(len(candidate_indices)),
        n_used=int(len(stamps)),
        stamp_size=int(size),
        selection_mode=selection_mode,
    )


def sample_injection_catalog(
    image_shape: tuple[int, int],
    real_positions: np.ndarray,
    *,
    count: int,
    magnitude_min: float,
    magnitude_max: float,
    fwhm_px: float,
    psf_size: int,
    rng: np.random.Generator,
    isolated_fraction: float = 0.5,
    isolated_min_real_sep_fwhm: float = 3.0,
    min_injected_sep_fwhm: float = 5.0,
    extra_margin_px: float = 0.0,
) -> pd.DataFrame:
    """Sample low-density artificial-star positions and crowding strata."""
    if count <= 0:
        raise ValueError("count must be positive")
    h, w = image_shape
    margin = int(math.ceil(psf_size / 2 + extra_margin_px))
    if w <= 2 * margin or h <= 2 * margin:
        raise ValueError("Image is too small for the configured PSF and margin")

    real_xy = np.asarray(real_positions, dtype=float)
    real_xy = real_xy[np.isfinite(real_xy).all(axis=1)] if real_xy.size else np.empty((0, 2))
    real_tree = cKDTree(real_xy) if len(real_xy) else None
    target_isolated = int(round(count * float(np.clip(isolated_fraction, 0.0, 1.0))))
    selected: list[tuple[float, float, str, float]] = []
    min_injected_sep = max(1.0, min_injected_sep_fwhm * fwhm_px)
    max_attempts = max(5000, count * 1000)

    def nearest_real(x: float, y: float) -> float:
        if real_tree is None:
            return float("inf")
        distance, _ = real_tree.query([x, y], k=1)
        return float(distance)

    def accept(x: float, y: float) -> bool:
        if not selected:
            return True
        current = np.asarray([(row[0], row[1]) for row in selected], dtype=float)
        return bool(np.min(np.hypot(current[:, 0] - x, current[:, 1] - y)) >= min_injected_sep)

    for desired_stratum, desired_count in (
        ("isolated", target_isolated),
        ("field", count - target_isolated),
    ):
        made = 0
        attempts = 0
        while made < desired_count and attempts < max_attempts:
            attempts += 1
            x = float(rng.uniform(margin, w - margin))
            y = float(rng.uniform(margin, h - margin))
            distance = nearest_real(x, y)
            if desired_stratum == "isolated" and distance < isolated_min_real_sep_fwhm * fwhm_px:
                continue
            if not accept(x, y):
                continue
            selected.append((x, y, desired_stratum, distance))
            made += 1
        if made < desired_count:
            raise RuntimeError(
                f"Could place only {made}/{desired_count} {desired_stratum} stars. "
                "Reduce stars_per_trial or separation constraints."
            )

    magnitudes = rng.uniform(magnitude_min, magnitude_max, size=count)
    rows = []
    for idx, ((x, y, stratum, nearest), mag) in enumerate(zip(selected, magnitudes), start=1):
        rows.append(
            {
                "injection_id": idx,
                "x_true": x,
                "y_true": y,
                "magnitude_true": float(mag),
                "placement_stratum": stratum,
                "nearest_real_sep_px": nearest,
                "nearest_real_sep_fwhm": nearest / fwhm_px,
            }
        )
    return pd.DataFrame(rows)


def inject_catalog(
    image: np.ndarray,
    psf: np.ndarray,
    catalog: pd.DataFrame,
    *,
    gain_e_per_adu: float,
    zeropoint_mag: float,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, pd.DataFrame]:
    """Inject source Poisson realizations into an existing observation.

    Existing image noise is retained unchanged. Only the additional source
    photons are sampled, so read noise and sky noise are not double-counted.
    """
    if gain_e_per_adu <= 0:
        raise ValueError("gain_e_per_adu must be positive")
    psf_arr = np.asarray(psf, dtype=float)
    if psf_arr.ndim != 2 or psf_arr.shape[0] % 2 != 1 or psf_arr.shape[1] % 2 != 1:
        raise ValueError("PSF must be a 2D odd-sized array")
    psf_sum = float(np.nansum(psf_arr))
    if not np.isfinite(psf_sum) or psf_sum <= 0:
        raise ValueError("PSF must have positive finite flux")
    psf_arr = np.where(np.isfinite(psf_arr), psf_arr / psf_sum, 0.0)

    injected = np.asarray(image, dtype=float).copy()
    expected_signal_adu = np.zeros_like(injected, dtype=float)
    realized_signal_adu = np.zeros_like(injected, dtype=float)
    half_y = psf_arr.shape[0] // 2
    half_x = psf_arr.shape[1] // 2
    out = catalog.copy()
    expected_fluxes = []
    realized_fluxes = []

    for row in out.itertuples(index=False):
        x = float(row.x_true)
        y = float(row.y_true)
        mag = float(row.magnitude_true)
        # Injected total source electrons for a star of magnitude ``mag`` under
        # the injection zeropoint. APEX Step 7 recovers an IRAF-style count-rate
        # magnitude mag_inst = INSTRUMENTAL_ZMAG - 2.5*log10(flux_e/exptime); for
        # a unit-exposure synthetic frame with zeropoint_mag == INSTRUMENTAL_ZMAG
        # the recovered mag_inst equals ``mag`` directly. Step 10 then adds zp_frame.
        expected_flux_e = float(10.0 ** (-0.4 * (mag - zeropoint_mag)))
        xi, yi = int(round(x)), int(round(y))
        x0, x1 = xi - half_x, xi + half_x + 1
        y0, y1 = yi - half_y, yi + half_y + 1
        if x0 < 0 or y0 < 0 or x1 > injected.shape[1] or y1 > injected.shape[0]:
            raise ValueError(f"Injection {row.injection_id} falls outside the image")

        shifted = ndi_shift(
            psf_arr,
            shift=(y - yi, x - xi),
            order=3,
            mode="constant",
            cval=0.0,
            prefilter=True,
        )
        shifted = np.where(np.isfinite(shifted), shifted, 0.0)
        shifted[shifted < 0] = 0.0
        shifted_sum = float(shifted.sum())
        if shifted_sum <= 0:
            raise RuntimeError(f"Shifted PSF has zero flux for injection {row.injection_id}")
        shifted /= shifted_sum

        expected_e = shifted * expected_flux_e
        realized_e = rng.poisson(expected_e).astype(float)
        expected_adu = expected_e / gain_e_per_adu
        realized_adu = realized_e / gain_e_per_adu
        expected_signal_adu[y0:y1, x0:x1] += expected_adu
        realized_signal_adu[y0:y1, x0:x1] += realized_adu
        injected[y0:y1, x0:x1] += realized_adu
        expected_fluxes.append(expected_flux_e)
        realized_fluxes.append(float(realized_e.sum()))

    out["flux_expected_e"] = expected_fluxes
    out["flux_realized_e"] = realized_fluxes
    return injected, expected_signal_adu, realized_signal_adu, out
