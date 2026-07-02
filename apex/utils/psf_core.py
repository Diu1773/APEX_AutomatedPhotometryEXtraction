"""Automatic crowded-core masking for PSF photometry."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PSFCoreCut:
    enabled: bool
    center_x: float = float("nan")
    center_y: float = float("nan")
    radius_px: float = float("nan")
    method: str = "off"
    n_total: int = 0
    n_excluded: int = 0
    n_kept: int = 0
    density_ratio: float = float("nan")
    reason: str = ""


def _finite_xy(xy: np.ndarray) -> np.ndarray:
    arr = np.asarray(xy, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2:
        return np.zeros((0, 2), dtype=float)
    arr = arr[:, :2]
    keep = np.isfinite(arr[:, 0]) & np.isfinite(arr[:, 1])
    return arr[keep]


def _auto_density_center(
    xy: np.ndarray,
    image_shape: tuple[int, int],
    fwhm_px: float,
    cell_fwhm_mult: float,
) -> tuple[float, float, float]:
    h, w = int(image_shape[0]), int(image_shape[1])
    cell = max(16.0, min(float(fwhm_px) * float(cell_fwhm_mult), max(16.0, min(h, w) / 4.0)))
    nx = max(2, int(np.ceil(w / cell)))
    ny = max(2, int(np.ceil(h / cell)))
    x_edges = np.linspace(0.0, float(w), nx + 1)
    y_edges = np.linspace(0.0, float(h), ny + 1)
    hist, _, _ = np.histogram2d(xy[:, 1], xy[:, 0], bins=[y_edges, x_edges])
    iy, ix = np.unravel_index(int(np.nanargmax(hist)), hist.shape)
    peak = float(hist[iy, ix])
    mean = float(len(xy) / max(1, hist.size))
    ratio = peak / max(mean, 1e-6)

    x0, x1 = x_edges[ix], x_edges[ix + 1]
    y0, y1 = y_edges[iy], y_edges[iy + 1]
    cx0 = 0.5 * (x0 + x1)
    cy0 = 0.5 * (y0 + y1)
    near = (
        (np.abs(xy[:, 0] - cx0) <= cell)
        & (np.abs(xy[:, 1] - cy0) <= cell)
    )
    if int(np.sum(near)) >= 5:
        return float(np.nanmedian(xy[near, 0])), float(np.nanmedian(xy[near, 1])), ratio
    return float(cx0), float(cy0), ratio


def _auto_radius(
    xy: np.ndarray,
    center_x: float,
    center_y: float,
    image_shape: tuple[int, int],
    fwhm_px: float,
    fallback_fwhm_mult: float,
    density_ratio: float,
) -> float:
    h, w = int(image_shape[0]), int(image_shape[1])
    r = np.hypot(xy[:, 0] - float(center_x), xy[:, 1] - float(center_y))
    r = r[np.isfinite(r)]
    if len(r) < 20:
        return float(fallback_fwhm_mult) * float(fwhm_px)

    max_r_edge = min(center_x, center_y, float(w) - center_x, float(h) - center_y)
    max_r_data = float(np.nanpercentile(r, 90.0))
    max_r = max(3.0 * float(fwhm_px), min(max_r_edge, max_r_data))
    bin_w = max(8.0, 5.0 * float(fwhm_px))
    if max_r <= 2.0 * bin_w:
        return float(fallback_fwhm_mult) * float(fwhm_px)

    edges = np.arange(0.0, max_r + bin_w, bin_w)
    if len(edges) < 5:
        return float(fallback_fwhm_mult) * float(fwhm_px)
    counts, _ = np.histogram(r, bins=edges)
    area = np.pi * (edges[1:] ** 2 - edges[:-1] ** 2)
    dens = counts / np.maximum(area, 1e-6)
    finite = np.isfinite(dens) & (area > 0)
    if int(np.sum(finite)) < 4:
        return float(fallback_fwhm_mult) * float(fwhm_px)

    outer = dens[max(1, len(dens) // 2):]
    outer = outer[np.isfinite(outer)]
    if len(outer) == 0:
        return float(fallback_fwhm_mult) * float(fwhm_px)
    bg = float(np.nanmedian(outer))
    if not np.isfinite(bg) or bg <= 0:
        return float(fallback_fwhm_mult) * float(fwhm_px)
    threshold = bg * max(1.0, float(density_ratio))
    min_r = max(3.0 * float(fwhm_px), edges[2] if len(edges) > 2 else 0.0)
    for i, d in enumerate(dens):
        rout = float(edges[i + 1])
        if rout < min_r:
            continue
        if np.isfinite(d) and d <= threshold:
            return rout
    return float(fallback_fwhm_mult) * float(fwhm_px)


def estimate_psf_core_cut(
    xy: np.ndarray,
    image_shape: tuple[int, int],
    fwhm_px: float,
    *,
    enabled: bool = True,
    center_mode: str = "auto",
    manual_center: tuple[float, float] | None = None,
    radius_px: float = 0.0,
    radius_fwhm_mult: float = 25.0,
    auto_cell_fwhm_mult: float = 8.0,
    auto_min_density_ratio: float = 1.5,
    auto_min_sources: int = 50,
    max_exclude_frac: float = 0.70,
) -> PSFCoreCut:
    """Estimate a central crowded-core exclusion from detection coordinates."""

    arr = _finite_xy(xy)
    n_total = int(len(arr))
    if not enabled:
        return PSFCoreCut(False, n_total=n_total, n_kept=n_total, reason="disabled")
    if n_total < int(auto_min_sources):
        return PSFCoreCut(False, n_total=n_total, n_kept=n_total, reason="too_few_sources")

    h, w = int(image_shape[0]), int(image_shape[1])
    fwhm = max(1.0, float(fwhm_px))
    mode = str(center_mode or "auto").strip().lower()
    density_ratio = float("nan")

    if mode == "manual":
        if manual_center is None:
            return PSFCoreCut(False, n_total=n_total, n_kept=n_total, reason="manual_center_missing")
        cx, cy = manual_center
        if not (np.isfinite(cx) and np.isfinite(cy)):
            return PSFCoreCut(False, n_total=n_total, n_kept=n_total, reason="manual_center_missing")
        method = "manual"
    elif mode == "image":
        cx, cy = 0.5 * (w - 1), 0.5 * (h - 1)
        method = "image_center"
    else:
        cx, cy, density_ratio = _auto_density_center(arr, image_shape, fwhm, auto_cell_fwhm_mult)
        if np.isfinite(density_ratio) and density_ratio < float(auto_min_density_ratio):
            return PSFCoreCut(
                False,
                center_x=cx,
                center_y=cy,
                n_total=n_total,
                n_kept=n_total,
                density_ratio=density_ratio,
                reason="density_peak_below_threshold",
            )
        method = "auto_density"

    try:
        radius = float(radius_px)
    except Exception:
        radius = 0.0
    if not np.isfinite(radius) or radius <= 0:
        radius = _auto_radius(
            arr,
            float(cx),
            float(cy),
            image_shape,
            fwhm,
            float(radius_fwhm_mult),
            float(auto_min_density_ratio),
        )
        method = f"{method}+auto_radius"
    else:
        method = f"{method}+manual_radius"

    radius = max(0.0, float(radius))
    if radius <= 0:
        return PSFCoreCut(False, n_total=n_total, n_kept=n_total, reason="invalid_radius")

    r = np.hypot(arr[:, 0] - float(cx), arr[:, 1] - float(cy))
    max_frac = min(max(float(max_exclude_frac), 0.05), 0.95)
    if np.mean(r < radius) > max_frac:
        radius = float(np.nanpercentile(r, 100.0 * max_frac))
        method = f"{method}+clamped"

    excluded = r < radius
    n_excl = int(np.sum(excluded))
    n_kept = int(n_total - n_excl)
    if n_excl <= 0:
        return PSFCoreCut(
            False,
            center_x=float(cx),
            center_y=float(cy),
            radius_px=float(radius),
            method=method,
            n_total=n_total,
            n_excluded=0,
            n_kept=n_total,
            density_ratio=density_ratio,
            reason="no_sources_inside_cut",
        )
    if n_kept <= 0:
        return PSFCoreCut(
            False,
            center_x=float(cx),
            center_y=float(cy),
            radius_px=float(radius),
            method=method,
            n_total=n_total,
            n_excluded=n_excl,
            n_kept=0,
            density_ratio=density_ratio,
            reason="no_sources_left",
        )

    return PSFCoreCut(
        True,
        center_x=float(cx),
        center_y=float(cy),
        radius_px=float(radius),
        method=method,
        n_total=n_total,
        n_excluded=n_excl,
        n_kept=n_kept,
        density_ratio=density_ratio,
    )


def psf_core_keep_mask(xy: np.ndarray, cut: PSFCoreCut) -> np.ndarray:
    arr = np.asarray(xy, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 2:
        return np.zeros((0,), dtype=bool)
    valid = np.isfinite(arr[:, 0]) & np.isfinite(arr[:, 1])
    if not cut.enabled or not np.isfinite(cut.radius_px) or cut.radius_px <= 0:
        return valid
    r = np.hypot(arr[:, 0] - float(cut.center_x), arr[:, 1] - float(cut.center_y))
    return valid & np.isfinite(r) & (r >= float(cut.radius_px))
