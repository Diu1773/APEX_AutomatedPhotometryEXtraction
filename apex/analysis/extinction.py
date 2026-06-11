"""Numerical routines for atmospheric extinction fitting."""

from __future__ import annotations

import numpy as np
import pandas as pd

from apex.utils.io_utils import coerce_int64_source_id
from apex.utils.photometry_utils import MAG_ERR_COEFF, MAD_TO_SIGMA, get_numeric_array


def robust_linfit(x, y, clip_sigma: float = 3.0, iters: int = 5, min_n: int = 10):
    """Robust line fit y = zp + k*x with MAD clipping."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    m0 = np.isfinite(x) & np.isfinite(y)
    x = x[m0]
    y = y[m0]
    if len(x) < min_n:
        return (np.nan, np.nan, np.nan, 0, np.nan)
    k, zp = np.polyfit(x, y, 1)
    base_n = len(x)
    for _ in range(int(iters)):
        yhat = zp + k * x
        r = y - yhat
        med = np.nanmedian(r)
        mad = np.nanmedian(np.abs(r - med)) + 1e-12
        sig = MAD_TO_SIGMA * mad
        keep = np.abs(r - med) <= float(clip_sigma) * sig
        if keep.sum() < min_n:
            break
        if keep.sum() == len(x):
            break
        x, y = x[keep], y[keep]
        k, zp = np.polyfit(x, y, 1)
    yhat = zp + k * x
    scatter = float(np.nanstd(y - yhat)) if len(x) else np.nan
    outlier_frac = float(1.0 - (len(x) / max(base_n, 1)))
    return (float(k), float(zp), scatter, int(len(x)), outlier_frac)


def robust_weighted_linfit(
    x,
    y,
    w=None,
    clip_sigma: float = 3.0,
    iters: int = 5,
    min_n: int = 10,
):
    """Robust weighted line fit y = zp + k*x with MAD clipping."""
    x = np.asarray(x, float)
    y = np.asarray(y, float)
    if w is not None:
        w = np.asarray(w, float)
        m0 = np.isfinite(x) & np.isfinite(y) & np.isfinite(w) & (w > 0)
    else:
        m0 = np.isfinite(x) & np.isfinite(y)
    x = x[m0]
    y = y[m0]
    w = w[m0] if w is not None else None
    if len(x) < min_n:
        return (np.nan, np.nan, np.nan, 0, np.nan)

    def _fit(xx, yy, ww):
        if ww is None:
            k_fit, zp_fit = np.polyfit(xx, yy, 1)
        else:
            k_fit, zp_fit = np.polyfit(xx, yy, 1, w=ww)
        return k_fit, zp_fit

    k, zp = _fit(x, y, w)
    base_n = len(x)
    for _ in range(int(iters)):
        yhat = zp + k * x
        r = y - yhat
        med = np.nanmedian(r)
        mad = np.nanmedian(np.abs(r - med)) + 1e-12
        sig = MAD_TO_SIGMA * mad
        keep = np.abs(r - med) <= float(clip_sigma) * sig
        if keep.sum() < min_n:
            break
        if keep.sum() == len(x):
            break
        x = x[keep]
        y = y[keep]
        if w is not None:
            w = w[keep]
        k, zp = _fit(x, y, w)

    yhat = zp + k * x
    scatter = float(np.nanstd(y - yhat)) if len(x) else np.nan
    outlier_frac = float(1.0 - (len(x) / max(base_n, 1)))
    return (float(k), float(zp), scatter, int(len(x)), outlier_frac)


def _frame_offset_basis(frame_airmass: np.ndarray) -> np.ndarray:
    """Basis for frame offsets with no constant or airmass-linear component."""
    x = np.asarray(frame_airmass, float)
    n_frame = int(len(x))
    if n_frame <= 0:
        return np.zeros((0, 0), dtype=float)
    x_centered = x - np.nanmean(x)
    constraints = [np.ones(n_frame, dtype=float)]
    if np.isfinite(x_centered).any() and np.nanmax(np.abs(x_centered)) > 0:
        constraints.append(np.nan_to_num(x_centered, nan=0.0))
    cmat = np.vstack(constraints)
    _, svals, vh = np.linalg.svd(cmat, full_matrices=True)
    tol = np.finfo(float).eps * max(cmat.shape) * (svals[0] if len(svals) else 1.0)
    rank = int((svals > tol).sum())
    basis = vh[rank:].T
    if basis.size == 0:
        return np.zeros((n_frame, 0), dtype=float)
    return np.asarray(basis, dtype=float)


def fit_ensemble_extinction_once(
    sub: pd.DataFrame,
    clip_sigma: float,
    iters: int,
    min_n: int,
):
    """Fit m_ij = s_i + z_j + k1*X_j for one date/filter group.

    The frame offset term z_j is constrained to contain no constant or
    airmass-linear component. Without that constraint, z_j and k1 are
    degenerate and the extinction coefficient is not identifiable.
    """
    if sub.empty:
        return None

    sub = sub.reset_index(drop=True)
    x = get_numeric_array(sub, "airmass")
    y = get_numeric_array(sub, "mag_inst")
    source_ids = coerce_int64_source_id(sub["source_id"])
    source_valid = source_ids.notna().to_numpy()
    star_ids = np.full(len(sub), -1, dtype=np.int64)
    if source_valid.any():
        star_ids[source_valid] = source_ids[source_valid].astype("int64").to_numpy()
    frame_ids = sub["file"].astype(str).to_numpy()

    base_mask = np.isfinite(x) & np.isfinite(y) & source_valid
    if base_mask.sum() < min_n:
        return None

    unique_stars = np.unique(star_ids[base_mask])
    unique_frames = np.unique(frame_ids[base_mask])
    n_star = len(unique_stars)
    n_frame = len(unique_frames)
    if n_star == 0 or n_frame == 0:
        return None

    star_index = {int(sid): i for i, sid in enumerate(unique_stars)}
    frame_index = {str(fid): i for i, fid in enumerate(unique_frames)}
    ref_frame = unique_frames[0]

    star_idx = np.full(len(sub), -1, dtype=int)
    for sid, idx_star in star_index.items():
        star_idx[star_ids == sid] = idx_star
    frame_idx = np.array([frame_index.get(str(fid), -1) for fid in frame_ids], dtype=int)
    base_mask &= (star_idx >= 0) & (frame_idx >= 0)
    base_n = int(base_mask.sum())
    if base_n < min_n:
        return None

    frame_airmass = np.full(n_frame, np.nan, dtype=float)
    for fid, idx_frame in frame_index.items():
        vals = x[base_mask & (frame_ids == fid)]
        if np.isfinite(vals).any():
            frame_airmass[idx_frame] = float(np.nanmedian(vals))
    z_basis = _frame_offset_basis(frame_airmass)
    n_z = int(z_basis.shape[1])

    if "snr" in sub.columns:
        snr_arr = get_numeric_array(sub, "snr", default=1.0)
        sig_mag = MAG_ERR_COEFF / np.clip(snr_arr, 1e-6, None)
        w_all = 1.0 / np.clip(sig_mag ** 2, 1e-12, None)
    else:
        w_all = np.ones(len(sub), dtype=float)

    def build_matrix(mask):
        idx = np.where(mask)[0]
        n_obs = len(idx)
        n_params = n_star + n_z + 1
        amat = np.zeros((n_obs, n_params), dtype=float)
        rows = np.arange(n_obs)
        amat[rows, star_idx[idx]] = 1.0
        if n_z:
            amat[:, n_star:n_star + n_z] = z_basis[frame_idx[idx], :]
        amat[:, -1] = x[idx]
        return amat, y[idx], idx

    mask_fit = base_mask.copy()
    for _ in range(int(iters)):
        amat, yv, idx = build_matrix(mask_fit)
        if len(yv) < min_n:
            break
        weights = np.sqrt(w_all[idx])
        coef, _, _, _ = np.linalg.lstsq(weights[:, None] * amat, weights * yv, rcond=None)
        resid_fit = yv - (amat @ coef)
        med = np.nanmedian(resid_fit)
        mad = np.nanmedian(np.abs(resid_fit - med)) + 1e-12
        sig = MAD_TO_SIGMA * mad
        keep = np.abs(resid_fit - med) <= float(clip_sigma) * sig
        if keep.sum() < min_n:
            break
        if keep.sum() == len(yv):
            break
        new_mask = mask_fit.copy()
        new_mask[idx] = keep
        mask_fit = new_mask

    amat, yv, idx = build_matrix(mask_fit)
    if len(yv) < min_n:
        return None
    weights = np.sqrt(w_all[idx])
    coef, _, _, _ = np.linalg.lstsq(weights[:, None] * amat, weights * yv, rcond=None)

    s_vals = coef[:n_star]
    if n_z:
        z_vals = z_basis @ coef[n_star:n_star + n_z]
    else:
        z_vals = np.zeros(n_frame, dtype=float)
    k1 = float(coef[-1])

    yhat = np.full(len(sub), np.nan, dtype=float)
    valid_all = base_mask.copy()
    yhat[valid_all] = (
        s_vals[star_idx[valid_all]]
        + z_vals[frame_idx[valid_all]]
        + k1 * x[valid_all]
    )
    resid = y - yhat
    delta_m = np.full(len(sub), np.nan, dtype=float)
    delta_m[valid_all] = y[valid_all] - (
        s_vals[star_idx[valid_all]] + z_vals[frame_idx[valid_all]]
    )
    scatter = float(np.nanstd(resid[base_mask])) if base_n else np.nan
    outlier_frac = float(1.0 - (mask_fit.sum() / max(base_n, 1)))

    return {
        "k1": k1,
        "s_vals": s_vals,
        "z_vals": z_vals,
        "star_index": star_index,
        "frame_index": frame_index,
        "ref_frame": ref_frame,
        "resid": resid,
        "delta_m": delta_m,
        "scatter": scatter,
        "n_used": int(mask_fit.sum()),
        "outlier_frac": outlier_frac,
        "n_star": n_star,
        "n_frame": n_frame,
    }


__all__ = [
    "fit_ensemble_extinction_once",
    "robust_linfit",
    "robust_weighted_linfit",
]
