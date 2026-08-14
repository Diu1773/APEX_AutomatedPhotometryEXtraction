"""Real-image artificial-star validation for APEX Step 4 and Step 8.

The module intentionally does not import the Step 8 GUI worker.  It mirrors
only the public numerical contract needed to make an injection frame and to
compare production tables without using ``det_uid`` or ``seed_uid``.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Any, Callable, Iterable

import numpy as np
import pandas as pd
from scipy.ndimage import map_coordinates, shift as ndi_shift
from scipy.spatial import cKDTree

from apex.analysis.psf_iteration import PSFFitFlag, qfit_noise_diagnostics


DEFAULT_TARGET_SNRS = (5.0, 10.0, 20.0, 50.0, 100.0)
DEFAULT_RADIUS_EDGES = (0.0, 0.25, 0.50, 0.75, 1.0)
DEFAULT_CROWDING_EDGES_FWHM = (0.75, 1.5, 3.0, 6.0, math.inf)


def _number(value: Any, default: float = np.nan) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if np.isfinite(result) else default


def _numeric(df: pd.DataFrame, name: str, default: float = np.nan) -> np.ndarray:
    if name not in df:
        return np.full(len(df), default, dtype=float)
    return pd.to_numeric(df[name], errors="coerce").to_numpy(float)


def _xy(df: pd.DataFrame, *, fitted: bool = False) -> tuple[np.ndarray, np.ndarray]:
    names = (("x_fit", "y_fit"), ("x", "y")) if fitted else (("x", "y"), ("x_fit", "y_fit"))
    for x_name, y_name in names:
        if x_name in df.columns and y_name in df.columns:
            return _numeric(df, x_name), _numeric(df, y_name)
    raise ValueError("table must contain x/y or x_fit/y_fit")


def _model_data_header(model: Any, header: Any | None) -> tuple[np.ndarray, Any | None]:
    if hasattr(model, "data"):
        data = np.asarray(model.data, dtype=float)
        if header is None and hasattr(model, "header"):
            header = model.header
    else:
        data = np.asarray(model, dtype=float)
    return data, header


def oversampled_epsf_to_native_kernel(
    model: Any,
    *,
    header: Any | None = None,
    oversampling: int | None = None,
    epsf_size: int | None = None,
    phase_x: float = 0.0,
    phase_y: float = 0.0,
) -> np.ndarray:
    """Convert a Step 8 oversampled ePSF into a normalized native kernel.

    Step 8 evaluates ``data[dy * os + cy, dx * os + cx]`` by linear
    interpolation and divides by ``data.sum() / os**2``.  This function uses
    that same coordinate mapping at native pixel offsets before normalizing
    the finite injection kernel. ``phase_x`` and ``phase_y`` are the source
    offsets from the nearest native-pixel centre. Directly normalizing the
    FITS array, or shifting an already downsampled kernel, would not preserve
    the Step 8 model for an undersampled PSF.
    """
    data, header = _model_data_header(model, header)
    if data.ndim != 2 or min(data.shape) < 3:
        raise ValueError("ePSF model must be a two-dimensional array")
    if oversampling is None:
        oversampling = int(_number(header.get("OVERSAMPL", 1), 1)) if header is not None else 1
    os = max(1, int(oversampling))
    if epsf_size is None:
        epsf_size = int(_number(header.get("EPSFSIZE", (min(data.shape) - 1) // os), 1)) if header is not None else (min(data.shape) - 1) // os
    size = int(epsf_size)
    if size < 3:
        raise ValueError("EPSFSIZE must be at least 3")
    if size % 2 == 0:
        size += 1

    cy, cx = data.shape[0] // 2, data.shape[1] // 2
    half = size // 2
    offsets = np.arange(-half, half + 1, dtype=float)
    yy, xx = np.meshgrid(offsets, offsets, indexing="ij")
    phase_x = float(phase_x)
    phase_y = float(phase_y)
    if not (np.isfinite(phase_x) and np.isfinite(phase_y)):
        raise ValueError("subpixel phases must be finite")
    coords = np.vstack((
        ((yy.ravel() - phase_y) * os + cy),
        ((xx.ravel() - phase_x) * os + cx),
    ))
    values = map_coordinates(data, coords, order=1, mode="constant", cval=0.0).reshape(size, size)

    # Preserve the evaluator's scale, then normalize the finite native kernel
    # so a requested source flux is the total flux actually injected.
    evaluator_norm = float(np.nansum(data)) / float(os * os)
    if not np.isfinite(evaluator_norm) or evaluator_norm <= 0:
        raise ValueError("ePSF model has no positive finite normalization")
    values = np.where(np.isfinite(values), values / evaluator_norm, 0.0)
    values = np.maximum(values, 0.0)
    total = float(values.sum())
    if not np.isfinite(total) or total <= 0:
        raise ValueError("native ePSF kernel has no positive finite flux")
    return (values / total).astype(np.float64)


def psf_noise_equivalent_area(kernel: np.ndarray) -> float:
    """Return the optimal-PSF noise-equivalent area in native pixels."""
    p = np.asarray(kernel, dtype=float)
    total = float(np.nansum(p))
    if p.ndim != 2 or total <= 0:
        raise ValueError("kernel must be a positive two-dimensional array")
    p = np.where(np.isfinite(p), p / total, 0.0)
    return float(1.0 / np.sum(p * p))


def optimal_psf_flux_for_snr(
    target_snr: float,
    *,
    gain_e_per_adu: float,
    background_rms_adu: float,
    psf_nea_px: float,
) -> tuple[float, float]:
    """Invert ``SNR = F / sqrt(F + NEA * sigma_background_e**2)``.

    Returns ``(source_electrons, source_adu)``.  The background/read noise in
    the real image is retained; only the source Poisson realization is added.
    """
    snr = float(target_snr)
    gain = float(gain_e_per_adu)
    rms = float(background_rms_adu)
    nea = float(psf_nea_px)
    if not (np.isfinite(snr) and snr > 0 and np.isfinite(gain) and gain > 0 and np.isfinite(rms) and rms >= 0 and np.isfinite(nea) and nea > 0):
        raise ValueError("target_snr, gain, background RMS, and NEA must be valid")
    background_variance_e = (rms * gain) ** 2 * nea
    snr2 = snr * snr
    flux_e = 0.5 * (snr2 + math.sqrt(snr2 * snr2 + 4.0 * snr2 * background_variance_e))
    return float(flux_e), float(flux_e / gain)


def _bin_label(value: float, edges: Iterable[float], *, suffix: str = "") -> str:
    edges = tuple(float(v) for v in edges)
    for lo, hi in zip(edges[:-1], edges[1:]):
        if value >= lo and value < hi:
            hi_label = "inf" if not np.isfinite(hi) else f"{hi:g}"
            return f"{lo:g}-{hi_label}{suffix}"
    return "out-of-range"


def sample_stratified_injections(
    image_shape: tuple[int, int],
    real_positions: np.ndarray,
    *,
    count: int,
    fwhm_px: float,
    rng: np.random.Generator,
    center_xy: tuple[float, float] | None = None,
    pixel_scale_arcsec: float | None = None,
    radius_edges: tuple[float, ...] = DEFAULT_RADIUS_EDGES,
    crowding_edges_fwhm: tuple[float, ...] = DEFAULT_CROWDING_EDGES_FWHM,
    min_real_sep_fwhm: float = 0.75,
    min_injected_sep_fwhm: float = 3.0,
    psf_size: int = 25,
    max_attempts: int | None = None,
    pair_fraction: float = 0.0,
    pair_separations_fwhm: tuple[float, ...] = (0.8, 1.2, 2.0, 3.0),
) -> pd.DataFrame:
    """Sample sparse positions across normalized radius and crowding bins.

    Candidate positions below ``min_real_sep_fwhm`` from any real source are
    rejected.  A global rejection check also enforces hard separation between
    injected stars, independent of their assigned stratum.

    ``pair_fraction`` puts the blend regime under the experimenter's control
    rather than the field's. Note what it is *not* needed for: a field whose
    stars are far apart on average can still host tight blends, because an
    injection only has to sit beside *one* real star and the stratified sampler
    goes looking for such spots. The LCO 1 m frame has a median nearest
    neighbour of 13 FWHM and nothing inside 1.5 FWHM, yet 400 injections filled
    the tightest bin with 102 stars — more than the M13 frame's 86 (2026-08-14).
    A predicted failure that did not happen.

    What it does buy is control and reach: a chosen separation instead of
    whatever the field offers, blends in frames with too few stars to sit beside,
    and separations closer than ``min_real_sep_fwhm`` allows against real stars.

    With ``pair_fraction`` at 0 the behaviour is unchanged, including the
    meaning of ``crowding_bin``; companions are only created when it is
    positive, and only then does the bin come from the nearest neighbour of
    *either* kind rather than from real stars alone.
    """
    if count <= 0 or not np.isfinite(fwhm_px) or fwhm_px <= 0:
        raise ValueError("count and fwhm_px must be positive")
    h, w = map(int, image_shape)
    half = int(psf_size) // 2
    if h <= 2 * half + 2 or w <= 2 * half + 2:
        raise ValueError("image is too small for the requested PSF kernel")
    real = np.asarray(real_positions, dtype=float)
    real = real.reshape((-1, 2)) if real.size else np.empty((0, 2), dtype=float)
    real = real[np.isfinite(real).all(axis=1)]
    tree = cKDTree(real) if len(real) else None
    cx, cy = center_xy if center_xy is not None else ((w - 1) / 2.0, (h - 1) / 2.0)
    max_radius = max(math.hypot(x - cx, y - cy) for x, y in ((0, 0), (w - 1, 0), (0, h - 1), (w - 1, h - 1)))
    min_real_sep = max(0.0, float(min_real_sep_fwhm)) * float(fwhm_px)
    min_injected_sep = max(1.0, float(min_injected_sep_fwhm) * float(fwhm_px))
    attempts_limit = max_attempts or max(10000, count * 3000)

    candidates: list[dict[str, Any]] = []
    attempts = 0
    while attempts < attempts_limit and len(candidates) < count * 120:
        attempts += 1
        x = float(rng.uniform(half + 1, w - half - 1))
        y = float(rng.uniform(half + 1, h - half - 1))
        real_sep = float(tree.query([x, y], k=1)[0]) if tree is not None else math.inf
        if real_sep < min_real_sep:
            continue
        radius_px = float(math.hypot(x - cx, y - cy))
        radius_fraction = min(1.0, radius_px / max_radius)
        crowding = real_sep / float(fwhm_px)
        candidates.append({
            "x_true": x,
            "y_true": y,
            "radius_px": radius_px,
            "radius_fraction": radius_fraction,
            "nearest_real_sep_px": real_sep,
            "nearest_real_sep_fwhm": crowding,
            "radius_bin": _bin_label(radius_fraction, radius_edges),
            "crowding_bin": _bin_label(crowding, crowding_edges_fwhm, suffix=" FWHM"),
        })
    if not candidates:
        raise RuntimeError("could not find any valid artificial-star positions")

    # Round-robin through the two-dimensional strata, then fill from the
    # remaining valid pool when a real image has an empty stratum.
    cells: dict[tuple[str, str], list[int]] = {}
    for i, candidate in enumerate(candidates):
        cells.setdefault((candidate["radius_bin"], candidate["crowding_bin"]), []).append(i)
    ordered_cells = [(r, c) for r in (_bin_label(v, radius_edges) for v in radius_edges[:-1]) for c in (_bin_label(v, crowding_edges_fwhm, suffix=" FWHM") for v in crowding_edges_fwhm[:-1])]
    chosen: list[int] = []
    chosen_xy: list[tuple[float, float]] = []

    def take(index: int) -> bool:
        point = candidates[index]
        if chosen_xy:
            xy = np.asarray(chosen_xy, dtype=float)
            if np.any(np.hypot(xy[:, 0] - point["x_true"], xy[:, 1] - point["y_true"]) < min_injected_sep):
                return False
        chosen.append(index)
        chosen_xy.append((point["x_true"], point["y_true"]))
        return True

    quota = max(1, int(math.ceil(count / max(1, len(ordered_cells)))))
    for cell in ordered_cells:
        for index in cells.get(cell, ()):
            if len(chosen) >= count:
                break
            if take(index) and sum(1 for i in chosen if (candidates[i]["radius_bin"], candidates[i]["crowding_bin"]) == cell) >= quota:
                break
        if len(chosen) >= count:
            break
    if len(chosen) < count:
        for index in range(len(candidates)):
            if index in chosen:
                continue
            if take(index) and len(chosen) >= count:
                break
    if len(chosen) < count:
        raise RuntimeError(f"could place only {len(chosen)}/{count} injections with the configured separations")

    rows: list[dict[str, Any]] = []
    for injection_id, index in enumerate(chosen, start=1):
        row = dict(candidates[index])
        row["injection_id"] = injection_id
        row["pair_id"] = 0
        row["is_pair_companion"] = False
        rows.append(row)

    fraction = min(1.0, max(0.0, float(pair_fraction)))
    if fraction > 0.0 and rows:
        separations = [s for s in pair_separations_fwhm
                       if np.isfinite(s) and float(s) > 0.0]
        if not separations:
            raise ValueError("pair_separations_fwhm must contain a positive value")
        # Companions replace primaries rather than adding to the total, so the
        # requested count stays the count.
        n_pairs = int(round(len(rows) * fraction / 2.0))
        n_pairs = max(0, min(n_pairs, len(rows) // 2))
        for pair_index in range(n_pairs):
            primary = rows[pair_index]
            companion = rows[len(rows) - 1 - pair_index]
            separation = float(rng.choice(separations)) * float(fwhm_px)
            for _ in range(200):
                angle = float(rng.uniform(0.0, 2.0 * math.pi))
                cx_new = primary["x_true"] + separation * math.cos(angle)
                cy_new = primary["y_true"] + separation * math.sin(angle)
                if not (half + 1 <= cx_new <= w - half - 1
                        and half + 1 <= cy_new <= h - half - 1):
                    continue
                if tree is not None and float(tree.query([cx_new, cy_new], k=1)[0]) < min_real_sep:
                    continue
                companion["x_true"], companion["y_true"] = cx_new, cy_new
                companion["radius_px"] = float(math.hypot(cx_new - cx, cy_new - cy))
                companion["radius_fraction"] = min(1.0, companion["radius_px"] / max_radius)
                companion["radius_bin"] = _bin_label(companion["radius_fraction"], radius_edges)
                real_sep_new = (float(tree.query([cx_new, cy_new], k=1)[0])
                                if tree is not None else math.inf)
                companion["nearest_real_sep_px"] = real_sep_new
                companion["nearest_real_sep_fwhm"] = real_sep_new / float(fwhm_px)
                companion["is_pair_companion"] = True
                companion["pair_id"] = pair_index + 1
                primary["pair_id"] = pair_index + 1
                break

        # With companions present the field's own crowding no longer describes
        # what each star sees, so the bin comes from the nearest neighbour of
        # either kind.
        injected = np.array([[r["x_true"], r["y_true"]] for r in rows], dtype=float)
        combined = np.vstack([real, injected]) if len(real) else injected
        combined_tree = cKDTree(combined)
        for row in rows:
            distances, _ = combined_tree.query([row["x_true"], row["y_true"]], k=2)
            nearest = float(distances[1])
            row["nearest_any_sep_px"] = nearest
            row["nearest_any_sep_fwhm"] = nearest / float(fwhm_px)
            row["crowding_bin"] = _bin_label(
                row["nearest_any_sep_fwhm"], crowding_edges_fwhm, suffix=" FWHM")
    else:
        for row in rows:
            row["nearest_any_sep_px"] = row["nearest_real_sep_px"]
            row["nearest_any_sep_fwhm"] = row["nearest_real_sep_fwhm"]

    for row in rows:
        if pixel_scale_arcsec is not None and np.isfinite(pixel_scale_arcsec):
            row["radius_arcmin"] = row["radius_px"] * float(pixel_scale_arcsec) / 60.0
        else:
            row["radius_arcmin"] = np.nan
    return pd.DataFrame(rows)


def inject_flux_catalog(
    image: np.ndarray,
    kernel: np.ndarray,
    catalog: pd.DataFrame,
    *,
    gain_e_per_adu: float,
    rng: np.random.Generator,
    kernel_sampler: Callable[[float, float], np.ndarray] | None = None,
    return_layers: bool = True,
) -> tuple[np.ndarray, np.ndarray | None, np.ndarray | None, pd.DataFrame]:
    """Inject source Poisson noise into a real image while retaining its noise."""
    if return_layers:
        image_arr = np.asarray(image, dtype=float)
    else:
        image_arr = np.asarray(image)
        if not np.issubdtype(image_arr.dtype, np.floating):
            image_arr = image_arr.astype(np.float32)
    p = np.asarray(kernel, dtype=float)
    if image_arr.ndim != 2 or p.ndim != 2 or p.shape[0] % 2 != 1 or p.shape[1] % 2 != 1:
        raise ValueError("image and kernel must be two-dimensional; kernel dimensions must be odd")
    if gain_e_per_adu <= 0:
        raise ValueError("gain_e_per_adu must be positive")
    p = np.where(np.isfinite(p), p, 0.0)
    p[p < 0] = 0.0
    total = float(p.sum())
    if total <= 0:
        raise ValueError("kernel must have positive flux")
    p /= total
    required = "x_true" if "x_true" in catalog else "x"
    x_name = required
    y_name = "y_true" if "y_true" in catalog else "y"
    flux_name = "true_flux_e" if "true_flux_e" in catalog else "flux_true_e"
    if flux_name not in catalog:
        raise ValueError(f"catalog must contain {flux_name}")
    injected = image_arr.copy()
    expected = np.zeros_like(image_arr, dtype=float) if return_layers else None
    realized = np.zeros_like(image_arr, dtype=float) if return_layers else None
    hy, hx = p.shape[0] // 2, p.shape[1] // 2
    out = catalog.copy()
    realized_fluxes: list[float] = []
    for row in out.itertuples(index=False):
        x, y = float(getattr(row, x_name)), float(getattr(row, y_name))
        flux_e = float(getattr(row, flux_name))
        xi, yi = int(round(x)), int(round(y))
        x0, x1, y0, y1 = xi - hx, xi + hx + 1, yi - hy, yi + hy + 1
        if x0 < 0 or y0 < 0 or x1 > image_arr.shape[1] or y1 > image_arr.shape[0]:
            raise ValueError("injection falls outside the image")
        phase_x = x - xi
        phase_y = y - yi
        if kernel_sampler is None:
            shifted = ndi_shift(
                p,
                shift=(phase_y, phase_x),
                order=1,
                mode="constant",
                cval=0.0,
                prefilter=False,
            )
        else:
            shifted = np.asarray(kernel_sampler(phase_x, phase_y), dtype=float)
            if shifted.shape != p.shape:
                raise ValueError("phase-aware kernel shape does not match the reference kernel")
        shifted = np.maximum(np.where(np.isfinite(shifted), shifted, 0.0), 0.0)
        shifted /= max(float(shifted.sum()), 1e-30)
        expected_e = shifted * flux_e
        realized_e = rng.poisson(np.maximum(expected_e, 0.0)).astype(float)
        expected_adu = expected_e / float(gain_e_per_adu)
        realized_adu = realized_e / float(gain_e_per_adu)
        if expected is not None and realized is not None:
            expected[y0:y1, x0:x1] += expected_adu
            realized[y0:y1, x0:x1] += realized_adu
        injected[y0:y1, x0:x1] += realized_adu
        realized_fluxes.append(float(realized_e.sum()))
    out["flux_realized_e"] = realized_fluxes
    out["flux_realized_adu"] = np.asarray(realized_fluxes) / float(gain_e_per_adu)
    return injected, expected, realized, out


def measure_preinjection_psf_residual(
    residual_image: np.ndarray,
    catalog: pd.DataFrame,
    *,
    kernel_sampler: Callable[[float, float], np.ndarray],
) -> pd.DataFrame:
    """Measure PSF-like signal already present before artificial injection.

    A PSF amplitude plus a local tilted plane is fit at every truth position.
    The result separates confusion/background structure from recovery errors
    introduced by source matching and PSF fitting.
    """

    image = np.asarray(residual_image)
    if image.ndim != 2:
        raise ValueError("residual_image must be two-dimensional")
    if not {"x_true", "y_true"} <= set(catalog.columns):
        raise ValueError("catalog must contain x_true and y_true")

    out = catalog.copy()
    amplitudes: list[float] = []
    local_rms: list[float] = []
    for row in out.itertuples(index=False):
        x = float(row.x_true)
        y = float(row.y_true)
        xi, yi = int(round(x)), int(round(y))
        kernel = np.asarray(kernel_sampler(x - xi, y - yi), dtype=float)
        if kernel.ndim != 2 or kernel.shape[0] % 2 != 1 or kernel.shape[1] % 2 != 1:
            raise ValueError("phase-aware kernel must have odd two-dimensional shape")
        hy, hx = kernel.shape[0] // 2, kernel.shape[1] // 2
        x0, x1 = xi - hx, xi + hx + 1
        y0, y1 = yi - hy, yi + hy + 1
        if x0 < 0 or y0 < 0 or x1 > image.shape[1] or y1 > image.shape[0]:
            amplitudes.append(np.nan)
            local_rms.append(np.nan)
            continue
        stamp = np.asarray(image[y0:y1, x0:x1], dtype=float)
        yy, xx = np.indices(kernel.shape, dtype=float)
        xx -= hx
        yy -= hy
        valid = np.isfinite(stamp) & np.isfinite(kernel)
        design = np.column_stack([
            kernel[valid],
            np.ones(int(valid.sum()), dtype=float),
            xx[valid],
            yy[valid],
        ])
        if int(valid.sum()) <= design.shape[1]:
            amplitudes.append(np.nan)
            local_rms.append(np.nan)
            continue
        coefficients, _, rank, _ = np.linalg.lstsq(design, stamp[valid], rcond=None)
        if rank < design.shape[1] or not np.all(np.isfinite(coefficients)):
            amplitudes.append(np.nan)
            local_rms.append(np.nan)
            continue
        residual = stamp[valid] - design @ coefficients
        amplitudes.append(float(coefficients[0]))
        local_rms.append(_robust_scatter(residual))

    out["preinjection_psf_residual_adu"] = amplitudes
    out["preinjection_local_rms_adu"] = local_rms
    if "true_flux_adu" in out:
        true_flux = _numeric(out, "true_flux_adu")
        fraction = np.divide(
            np.asarray(amplitudes, dtype=float),
            true_flux,
            out=np.full(len(out), np.nan),
            where=np.isfinite(true_flux) & (true_flux > 0),
        )
        out["preinjection_psf_residual_frac"] = fraction
        target_snr = _numeric(out, "target_snr")
        clean_limit = np.maximum(
            0.05,
            np.divide(
                3.0,
                target_snr,
                out=np.full(len(out), np.inf),
                where=np.isfinite(target_snr) & (target_snr > 0),
            ),
        )
        out["preinjection_clean_limit_frac"] = clean_limit
        out["confusion_clean"] = np.isfinite(fraction) & (np.abs(fraction) <= clean_limit)
    return out


def add_forced_truth_to_step7(
    step7: pd.DataFrame,
    truth: pd.DataFrame,
    *,
    filename: str,
) -> pd.DataFrame:
    """Append artificial-star truth positions as Step 7 forced PSF seeds.

    Existing real-star rows and their aperture measurements are preserved for
    flux-scale estimation. Artificial rows provide only the position and
    positive flux seed fields consumed by Step 8.
    """
    required = {"x_true", "y_true", "true_flux_e", "true_flux_adu"}
    if not required <= set(truth.columns):
        raise ValueError(f"truth table is missing {sorted(required - set(truth.columns))}")

    base = step7.copy()
    context: dict[str, Any] = {}
    if not base.empty:
        first = base.iloc[0]
        for name in (
            "filter",
            "FILTER",
            "gain_e_per_adu",
            "exptime",
            "rdnoise_e",
            "binning_x",
            "binning_y",
            "gain_source",
            "rdnoise_source",
        ):
            if name in base.columns:
                context[name] = first[name]

    master = pd.to_numeric(base.get("master_id", pd.Series(dtype=float)), errors="coerce")
    finite_master = master[np.isfinite(master)]
    first_master = int(finite_master.max()) + 1 if len(finite_master) else 1
    rows: list[dict[str, Any]] = []
    for offset, item in enumerate(truth.itertuples(index=False)):
        x = float(item.x_true)
        y = float(item.y_true)
        flux_e = float(item.true_flux_e)
        flux_adu = float(item.true_flux_adu)
        row = {name: np.nan for name in base.columns}
        row.update(context)
        row.update({
            "file": filename,
            "master_id": first_master + offset,
            "source_id": -1,
            "ID": first_master + offset,
            "det_uid": -1,
            "x": x,
            "y": y,
            "x_pred": x,
            "y_pred": y,
            "x_reg": x,
            "y_reg": y,
            "x_fit": x,
            "y_fit": y,
            "registration_method": "artificial_truth",
            "registration_anchor": False,
            "detected_flag": False,
            "forced_flag": True,
            "recentered_flag": False,
            "recenter_method": "artificial_truth",
            "step4_quality_ok": False,
            "step4_quality_used": False,
            "step4_anchor_candidate": False,
            "step4_apcorr_candidate": False,
            "centroid_outlier": False,
            "centering_quality": "artificial_forced",
            "off_frame_flag": False,
            "flux": flux_e,
            "flux_e": flux_e,
            "flux_net_adu": flux_adu,
            "is_saturated": False,
            "is_nonlinear": False,
            "bad_phot_flag": False,
        })
        rows.append(row)

    artificial = pd.DataFrame(rows)
    if base.empty:
        return artificial
    return pd.concat([base, artificial.reindex(columns=base.columns)], ignore_index=True)


def _global_match(truth_xy: np.ndarray, product_xy: np.ndarray, radius_px: float) -> tuple[np.ndarray, np.ndarray]:
    matches = np.full(len(truth_xy), -1, dtype=int)
    distances = np.full(len(truth_xy), np.nan, dtype=float)
    valid_truth = np.isfinite(truth_xy).all(axis=1)
    valid_product = np.isfinite(product_xy).all(axis=1)
    if not valid_truth.any() or not valid_product.any():
        return matches, distances
    tree = cKDTree(product_xy[valid_product])
    product_indices = np.flatnonzero(valid_product)
    candidates: list[tuple[float, int, int]] = []
    for i in np.flatnonzero(valid_truth):
        for local_j in tree.query_ball_point(truth_xy[i], float(radius_px)):
            j = int(product_indices[local_j])
            candidates.append((float(np.hypot(*(truth_xy[i] - product_xy[j]))), int(i), j))
    used_truth: set[int] = set()
    used_product: set[int] = set()
    for distance, i, j in sorted(candidates):
        if i in used_truth or j in used_product:
            continue
        used_truth.add(i)
        used_product.add(j)
        matches[i] = j
        distances[i] = distance
    return matches, distances


def apply_recovery_quality_policy(recovery: pd.DataFrame) -> pd.DataFrame:
    """Recompute the science-quality recovery mask from persisted diagnostics."""
    result = recovery.copy()
    if "psf_positive_recovered" in result:
        quality = result["psf_positive_recovered"].astype(bool).to_numpy(copy=True)
    elif "flux_recovered_e" in result:
        flux = _numeric(result, "flux_recovered_e")
        quality = np.isfinite(flux) & (flux > 0)
        result["psf_positive_recovered"] = quality
    elif "psf_recovered" in result:
        quality = result["psf_recovered"].astype(bool).to_numpy(copy=True)
        result["psf_positive_recovered"] = quality
    else:
        quality = np.zeros(len(result), dtype=bool)
        result["psf_positive_recovered"] = quality
    qfit = _numeric(result, "qfit")
    n_pixels = _numeric(result, "n_pixels_fit")
    snr = _numeric(result, "snr_psf")
    nea = _numeric(result, "psf_nea_px")
    expected_qfit, qfit_noise_ratio = qfit_noise_diagnostics(
        qfit,
        n_pixels,
        snr,
        nea,
    )
    result["qfit_noise_expected"] = expected_qfit
    result["qfit_noise_ratio"] = qfit_noise_ratio
    qfit_ok = ~np.isfinite(qfit_noise_ratio) | (qfit_noise_ratio <= 3.0)
    result["qfit_quality_ok"] = qfit_ok
    quality &= qfit_ok
    if "reduced_chi2" in result:
        values = _numeric(result, "reduced_chi2")
        if np.isfinite(values).any():
            quality &= np.isfinite(values) & (values <= 25.0)
    if "cfit" in result:
        values = _numeric(result, "cfit")
        if np.isfinite(values).any():
            quality &= np.isfinite(values) & (np.abs(values) <= 0.1)
    if "flags_psf" in result:
        flags = _numeric(result, "flags_psf")
        if np.isfinite(flags).any():
            severe_flags = int(
                PSFFitFlag.NONPOSITIVE_FLUX
                | PSFFitFlag.NONCONVERGENCE
                | PSFFitFlag.NO_OVERLAP
                | PSFFitFlag.NONFINITE_POSITION
                | PSFFitFlag.NONFINITE_FLUX
            )
            flags_int = np.where(np.isfinite(flags), flags, severe_flags).astype(np.int64)
            flags_ok = np.isfinite(flags) & ((flags_int & severe_flags) == 0)
            result["fit_flags_ok"] = flags_ok
            quality &= flags_ok
    if "nearest_real_sep_fwhm" in result:
        nearest = _numeric(result, "nearest_real_sep_fwhm")
        crowding_ok = np.isfinite(nearest) & (nearest >= 1.5)
        result["crowding_ok"] = crowding_ok
        quality &= crowding_ok
    result["psf_recovered"] = quality
    detected = (
        result["detection_recovered"].astype(bool).to_numpy()
        if "detection_recovered" in result
        else np.zeros(len(result), dtype=bool)
    )
    result["blind_psf_recovered"] = quality & detected
    return result


def match_injections_to_products(
    truth: pd.DataFrame,
    detections: pd.DataFrame,
    psf_photometry: pd.DataFrame,
    *,
    radius_px: float = 1.5,
) -> pd.DataFrame:
    """Match truth to Step4 and Step8 products independently, one-to-one."""
    if not {"x_true", "y_true"} <= set(truth.columns):
        raise ValueError("truth table must contain x_true and y_true")
    if radius_px <= 0:
        raise ValueError("radius_px must be positive")
    truth_xy = truth[["x_true", "y_true"]].to_numpy(float)
    det_xy = np.column_stack(_xy(detections, fitted=False)) if len(detections) else np.empty((0, 2))
    psf_xy = np.column_stack(_xy(psf_photometry, fitted=True)) if len(psf_photometry) else np.empty((0, 2))
    det_match, det_distance = _global_match(truth_xy, det_xy, radius_px)
    psf_match, psf_distance = _global_match(truth_xy, psf_xy, radius_px)
    result = truth.copy().reset_index(drop=True)
    result["detection_recovered"] = det_match >= 0
    result["detection_row"] = np.where(det_match >= 0, det_match, np.nan)
    result["detection_positional_error_px"] = det_distance
    result["psf_row"] = np.where(psf_match >= 0, psf_match, np.nan)
    result["psf_positional_error_px"] = psf_distance

    flux = np.full(len(result), np.nan, dtype=float)
    for name in ("flux_psf_e", "flux_e", "flux"):
        if name in psf_photometry.columns:
            values = _numeric(psf_photometry, name)
            good = (psf_match >= 0) & np.isfinite(psf_match)
            flux[good] = values[psf_match[good]]
            break
    result["flux_recovered_e"] = flux
    result["flux_recovered_adu"] = flux / _numeric(result, "gain_e_per_adu", 1.0) if "gain_e_per_adu" in result else np.nan
    flux_error = np.full(len(result), np.nan, dtype=float)
    if "flux_psf_err_e" in psf_photometry.columns:
        product_flux_error = _numeric(psf_photometry, "flux_psf_err_e")
        good = psf_match >= 0
        flux_error[good] = product_flux_error[psf_match[good]]
    result["flux_recovered_err_e"] = flux_error
    positive_recovered = (psf_match >= 0) & np.isfinite(flux) & (flux > 0)
    result["psf_positive_recovered"] = positive_recovered
    true_flux = _numeric(result, "true_flux_e")
    ratio = np.divide(flux, true_flux, out=np.full(len(result), np.nan), where=np.isfinite(flux) & (true_flux > 0))
    result["flux_frac_error"] = ratio - 1.0
    result["flux_frac_formal_err"] = np.divide(
        flux_error,
        true_flux,
        out=np.full(len(result), np.nan),
        where=np.isfinite(flux_error) & (flux_error > 0) & (true_flux > 0),
    )
    result["flux_pull"] = np.divide(
        flux - true_flux,
        flux_error,
        out=np.full(len(result), np.nan),
        where=np.isfinite(flux_error) & (flux_error > 0),
    )
    delta_mag = np.full(len(result), np.nan, dtype=float)
    positive_ratio = ratio > 0
    delta_mag[positive_ratio] = -2.5 * np.log10(ratio[positive_ratio])
    result["delta_mag"] = delta_mag

    raw_flux = np.full(len(result), np.nan, dtype=float)
    if "flux_psf_raw_e" in psf_photometry.columns:
        product_raw_flux = _numeric(psf_photometry, "flux_psf_raw_e")
        good = psf_match >= 0
        raw_flux[good] = product_raw_flux[psf_match[good]]
    result["flux_recovered_raw_e"] = raw_flux
    raw_flux_error = np.full(len(result), np.nan, dtype=float)
    if "flux_psf_err_raw_e" in psf_photometry.columns:
        product_raw_error = _numeric(psf_photometry, "flux_psf_err_raw_e")
        good = psf_match >= 0
        raw_flux_error[good] = product_raw_error[psf_match[good]]
    result["flux_recovered_raw_err_e"] = raw_flux_error
    raw_ratio = np.divide(
        raw_flux,
        true_flux,
        out=np.full(len(result), np.nan),
        where=np.isfinite(raw_flux) & (true_flux > 0),
    )
    result["raw_flux_frac_error"] = raw_ratio - 1.0
    raw_delta_mag = np.full(len(result), np.nan, dtype=float)
    positive_raw_ratio = raw_ratio > 0
    raw_delta_mag[positive_raw_ratio] = -2.5 * np.log10(
        raw_ratio[positive_raw_ratio]
    )
    result["raw_delta_mag"] = raw_delta_mag
    result["positional_error_px"] = np.where(np.isfinite(psf_distance), psf_distance, det_distance)

    for name in (
        "qfit",
        "cfit",
        "reduced_chi2",
        "flags_psf",
        "snr_psf",
        "iter_found",
        "n_pixels_fit",
    ):
        values = np.full(len(result), np.nan if name != "flags_psf" else np.nan, dtype=float)
        if name in psf_photometry.columns:
            product_values = _numeric(psf_photometry, name)
            good = psf_match >= 0
            values[good] = product_values[psf_match[good]]
        result[name] = values

    forced = np.zeros(len(result), dtype=bool)
    if "forced_psf" in psf_photometry.columns:
        raw_forced = psf_photometry["forced_psf"]
        if raw_forced.dtype == bool:
            product_forced = raw_forced.fillna(False).to_numpy(bool)
        else:
            product_forced = raw_forced.astype(str).str.strip().str.lower().isin(
                {"1", "true", "t", "yes", "y"}
            ).to_numpy(bool)
        good = psf_match >= 0
        forced[good] = product_forced[psf_match[good]]
    result["forced_psf"] = forced
    return apply_recovery_quality_policy(result)


def _robust_scatter(values: np.ndarray) -> float:
    values = values[np.isfinite(values)]
    return float(1.4826 * np.median(np.abs(values - np.median(values)))) if len(values) else np.nan


def _metric_row(group: pd.DataFrame, scope: str, label: str) -> dict[str, Any]:
    frac = _numeric(group, "flux_frac_error")
    mag = _numeric(group, "delta_mag")
    raw_frac = _numeric(group, "raw_flux_frac_error")
    raw_mag = _numeric(group, "raw_delta_mag")
    pull = _numeric(group, "flux_pull")
    formal_frac = _numeric(group, "flux_frac_formal_err")
    recovered = group["psf_recovered"].astype(bool).to_numpy()
    valid = recovered & np.isfinite(frac)
    raw_valid = recovered & np.isfinite(raw_frac)
    frac_valid = frac[valid]
    mag_valid = mag[recovered & np.isfinite(mag)]
    raw_frac_valid = raw_frac[raw_valid]
    raw_mag_valid = raw_mag[recovered & np.isfinite(raw_mag)]
    pull_valid = pull[recovered & np.isfinite(pull)]
    formal_frac_valid = formal_frac[recovered & np.isfinite(formal_frac)]
    flux_bias = float(np.median(frac_valid)) if len(frac_valid) else np.nan
    flux_scatter = _robust_scatter(frac_valid)
    flux_rmse = float(np.sqrt(np.mean(frac_valid * frac_valid))) if len(frac_valid) else np.nan
    pull_scatter = _robust_scatter(pull_valid)
    median_formal_frac = (
        float(np.median(formal_frac_valid)) if len(formal_frac_valid) else np.nan
    )
    systematic_floor_frac = (
        float(np.sqrt(max(flux_scatter * flux_scatter - median_formal_frac * median_formal_frac, 0.0)))
        if np.isfinite(flux_scatter) and np.isfinite(median_formal_frac)
        else np.nan
    )
    psf_completeness = float(group["psf_recovered"].astype(bool).mean()) if len(group) else np.nan
    blind_recovered = (
        group["blind_psf_recovered"].astype(bool)
        if "blind_psf_recovered" in group
        else group["psf_recovered"].astype(bool) & group["detection_recovered"].astype(bool)
    )
    return {
        "scope": scope,
        "label": label,
        "n_truth": int(len(group)),
        "n_detection_recovered": int(group["detection_recovered"].astype(bool).sum()),
        "n_psf_recovered": int(group["psf_recovered"].astype(bool).sum()),
        "n_blind_psf_recovered": int(blind_recovered.sum()),
        "n_valid_flux": int(len(frac_valid)),
        "detection_completeness": float(group["detection_recovered"].astype(bool).mean()) if len(group) else np.nan,
        "psf_completeness": psf_completeness,
        "blind_psf_completeness": float(blind_recovered.mean()) if len(group) else np.nan,
        "flux_bias": flux_bias,
        "flux_scatter_robust": flux_scatter,
        "flux_rmse": flux_rmse,
        "bias": flux_bias,
        "scatter": flux_scatter,
        "rmse": flux_rmse,
        "completeness": psf_completeness,
        "delta_mag_bias": float(np.median(mag_valid)) if len(mag_valid) else np.nan,
        "delta_mag_scatter_robust": _robust_scatter(mag_valid),
        "delta_mag_rmse": float(np.sqrt(np.mean(mag_valid * mag_valid))) if len(mag_valid) else np.nan,
        "raw_flux_bias": float(np.median(raw_frac_valid)) if len(raw_frac_valid) else np.nan,
        "raw_flux_scatter_robust": _robust_scatter(raw_frac_valid),
        "raw_flux_rmse": (
            float(np.sqrt(np.mean(raw_frac_valid * raw_frac_valid)))
            if len(raw_frac_valid)
            else np.nan
        ),
        "raw_delta_mag_bias": float(np.median(raw_mag_valid)) if len(raw_mag_valid) else np.nan,
        "raw_delta_mag_scatter_robust": _robust_scatter(raw_mag_valid),
        "raw_delta_mag_rmse": (
            float(np.sqrt(np.mean(raw_mag_valid * raw_mag_valid)))
            if len(raw_mag_valid)
            else np.nan
        ),
        "flux_pull_median": float(np.median(pull_valid)) if len(pull_valid) else np.nan,
        "flux_pull_scatter_robust": pull_scatter,
        "recommended_error_scale": (
            float(max(1.0, pull_scatter)) if np.isfinite(pull_scatter) else np.nan
        ),
        "median_formal_flux_frac_err": median_formal_frac,
        "empirical_systematic_floor_frac": systematic_floor_frac,
        "empirical_systematic_floor_mag": (
            float(2.5 / np.log(10.0) * systematic_floor_frac)
            if np.isfinite(systematic_floor_frac)
            else np.nan
        ),
    }


def aggregate_recovery_metrics(recovery: pd.DataFrame) -> pd.DataFrame:
    """Aggregate robust bias, scatter, RMSE, and completeness by requested strata."""
    recovery = apply_recovery_quality_policy(recovery)
    rows = [_metric_row(recovery, "overall", "all")]
    specifications = (("target_snr", "target_snr"), ("radius_bin", "radius"), ("crowding_bin", "crowding"))
    for field, scope in specifications:
        if field not in recovery:
            continue
        for label, group in recovery.groupby(
            field,
            dropna=False,
            sort=field == "target_snr",
        ):
            rows.append(_metric_row(group, scope, str(label)))
    if "target_snr" in recovery:
        for field, scope in (
            ("radius_bin", "snr_x_radius"),
            ("crowding_bin", "snr_x_crowding"),
        ):
            if field not in recovery:
                continue
            for (snr, stratum), group in recovery.groupby(
                ["target_snr", field],
                dropna=False,
                sort=True,
            ):
                rows.append(
                    _metric_row(
                        group,
                        scope,
                        f"SNR={float(snr):g} | {stratum}",
                    )
                )
    if "confusion_clean" in recovery:
        clean_mask = recovery["confusion_clean"].astype(bool)
        clean = recovery.loc[clean_mask]
        rows.append(_metric_row(clean, "confusion_clean", "all"))
        if "target_snr" in clean:
            for snr, group in clean.groupby("target_snr", dropna=False, sort=True):
                rows.append(
                    _metric_row(group, "target_snr_clean", str(snr))
                )
    return pd.DataFrame(rows)


def plot_recovery_summary(summary: pd.DataFrame, output_path: str | Path) -> None:
    """Write a readable four-panel benchmark summary figure."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5), constrained_layout=True)
    for ax, scope, title in zip(axes.flat, ("target_snr", "radius", "crowding", "target_snr"), ("By target SNR", "By normalized cluster radius", "By crowding", "Flux bias and scatter by SNR")):
        group = summary[summary["scope"] == scope].copy()
        x = np.arange(len(group))
        if scope == "target_snr" and title.startswith("Flux"):
            ax.axhline(0.0, color="0.3", lw=0.8)
            enough = pd.to_numeric(group["n_valid_flux"], errors="coerce") >= 5
            ax.plot(x, group["raw_flux_bias"].where(enough), "o-", label="raw PSF bias")
            ax.plot(x, group["flux_bias"].where(enough), "^-", label="scaled PSF bias")
            ax.plot(x, group["flux_scatter_robust"].where(enough), "s--", label="robust scatter")
            clean = summary[summary["scope"] == "target_snr_clean"].copy()
            if len(clean):
                clean = clean.set_index("label").reindex(group["label"].astype(str))
                clean_enough = pd.to_numeric(clean["n_valid_flux"], errors="coerce") >= 3
                ax.plot(
                    x,
                    clean["flux_scatter_robust"].where(clean_enough),
                    "d:",
                    label="clean-background scatter",
                )
            ax.set_ylabel("fractional flux")
        else:
            ax.plot(x, group["psf_completeness"], "o-", label="usable forced PSF")
            ax.plot(x, group["detection_completeness"], "s--", label="Step4 completeness")
            ax.set_ylabel("completeness")
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(group["label"].astype(str), rotation=35, ha="right")
        ax.grid(axis="y", alpha=0.25)
        ax.legend(fontsize=8)
    fig.suptitle("APEX real-image PSF artificial-star recovery")
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def write_benchmark_outputs(
    truth: pd.DataFrame,
    recovery: pd.DataFrame,
    summary: pd.DataFrame,
    output_dir: str | Path,
    *,
    metadata: dict[str, Any] | None = None,
) -> dict[str, Path]:
    """Write the required CSV, JSON, and PNG products."""
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    paths = {
        "truth": output / "truth.csv",
        "recovery": output / "recovery.csv",
        "summary": output / "summary.csv",
        "summary_json": output / "summary.json",
        "plot": output / "psf_artificial_stars.png",
    }
    truth.to_csv(paths["truth"], index=False)
    recovery.to_csv(paths["recovery"], index=False)
    summary.to_csv(paths["summary"], index=False)
    payload = {"metadata": metadata or {}, "summary": summary.to_dict(orient="records")}
    paths["summary_json"].write_text(json.dumps(payload, indent=2, allow_nan=True, default=str), encoding="utf-8")
    plot_recovery_summary(summary, paths["plot"])
    return paths


__all__ = [
    "DEFAULT_TARGET_SNRS",
    "add_forced_truth_to_step7",
    "aggregate_recovery_metrics",
    "apply_recovery_quality_policy",
    "inject_flux_catalog",
    "match_injections_to_products",
    "measure_preinjection_psf_residual",
    "optimal_psf_flux_for_snr",
    "oversampled_epsf_to_native_kernel",
    "plot_recovery_summary",
    "psf_noise_equivalent_area",
    "sample_stratified_injections",
    "write_benchmark_outputs",
]
