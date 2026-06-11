"""SYSREM systematic-noise removal for differential light curves.

Implements the algorithm of Tamuz, Mazeh & Zucker (2005), MNRAS 356, 1466.

The method iteratively identifies and removes common systematic trends shared
across comparison stars.  Each iteration extracts one linear systematic
component (c_j, the per-frame trend vector) and the corresponding per-star
amplitude (a_i).  The residuals are then passed to the next iteration.

Typical usage
-------------
    from apex.analysis.light_curve.sysrem import sysrem, apply_sysrem

    # mag_matrix: shape (n_stars, n_frames)  — differential magnitudes
    # err_matrix: shape (n_stars, n_frames)  — photometric errors
    components, corrections = sysrem(mag_matrix, err_matrix, n_iter=3)

    # Apply to a single target light curve
    corrected = apply_sysrem(target_mag, target_err, corrections)
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Optional


@dataclass
class SysremComponent:
    """One extracted systematic component."""
    c: np.ndarray   # per-frame trend vector, shape (n_frames,)
    a: np.ndarray   # per-star amplitude, shape (n_stars,)
    rms_before: float
    rms_after: float
    iteration: int


@dataclass
class SysremResult:
    """Full SYSREM output."""
    components: list[SysremComponent] = field(default_factory=list)
    residuals: Optional[np.ndarray] = None   # final residuals, (n_stars, n_frames)


def sysrem(
    mag_matrix: np.ndarray,
    err_matrix: np.ndarray,
    n_iter: int = 5,
    max_inner_iter: int = 50,
    tol: float = 1e-6,
) -> SysremResult:
    """Extract *n_iter* systematic components from a differential light curve matrix.

    Parameters
    ----------
    mag_matrix : (n_stars, n_frames) array of differential magnitudes.
        NaN marks missing observations.
    err_matrix : (n_stars, n_frames) array of photometric errors (> 0).
        NaN or non-positive values → treated as large errors (low weight).
    n_iter : number of systematic components to remove.
    max_inner_iter : maximum iterations for the inner convergence loop.
    tol : convergence tolerance on the relative change in ``c``.

    Returns
    -------
    SysremResult with extracted components and final residual matrix.
    """
    mag = np.array(mag_matrix, dtype=float)
    err = np.array(err_matrix, dtype=float)

    n_stars, n_frames = mag.shape

    # Weights: w_ij = 1 / sigma_ij^2; zero only for missing magnitudes.
    # Missing/invalid errors are retained with a low fallback weight so callers
    # without formal errors still get a usable SYSREM solution.
    valid_mag = np.isfinite(mag)
    valid_err = np.isfinite(err) & (err > 0)
    if np.any(valid_err):
        fallback_err = float(np.nanmedian(err[valid_err])) * 10.0
        if not np.isfinite(fallback_err) or fallback_err <= 0:
            fallback_err = 1.0
    else:
        fallback_err = 1.0
    safe_err = np.where(valid_err, err, fallback_err)
    w = np.where(valid_mag, 1.0 / safe_err**2, 0.0)
    observed = w > 0  # mask of valid measurements

    # Subtract per-star weighted mean (vectorised)
    r = np.where(observed, mag, 0.0)
    sw = w.sum(axis=1, keepdims=True)                # (n_stars, 1)
    weighted_sum = (w * r).sum(axis=1, keepdims=True)
    mean_w = np.divide(
        weighted_sum,
        sw,
        out=np.zeros_like(weighted_sum),
        where=sw > 0,
    )
    r -= mean_w
    r = np.where(observed, r, 0.0)

    def _rms(arr: np.ndarray) -> float:
        # RMS over observed cells only; avoids phantom contributions from
        # the missing-data positions (which are zeroed for matrix ops but
        # accumulate -a·c after each component subtraction).
        if not observed.any():
            return float("nan")
        return float(np.sqrt(np.mean(arr[observed] ** 2)))

    result = SysremResult(residuals=np.where(observed, r, np.nan))

    for it in range(n_iter):
        rms_before = _rms(r)

        # Initialise c_j: weighted mean of residuals across stars per frame
        sw_f = w.sum(axis=0)                          # (n_frames,)
        c = np.where(sw_f > 0, (w * r).sum(axis=0) / sw_f, 0.0)

        # Iterative (a_i, c_j) estimation — fully vectorised
        for _ in range(max_inner_iter):
            with np.errstate(invalid="ignore", divide="ignore"):
                # a_i = Σ_j w_ij r_ij c_j / Σ_j w_ij c_j²
                denom_a = (w * c[None, :] ** 2).sum(axis=1)   # (n_stars,)
                numer_a = (w * r * c[None, :]).sum(axis=1)
                a = np.where(denom_a > 0, numer_a / denom_a, 0.0)

                # c_j = Σ_i w_ij r_ij a_i / Σ_i w_ij a_i²
                denom_c = (w * a[:, None] ** 2).sum(axis=0)   # (n_frames,)
                numer_c = (w * r * a[:, None]).sum(axis=0)
                c_new = np.where(denom_c > 0, numer_c / denom_c, 0.0)

            norm_c = np.linalg.norm(c_new)
            delta = np.linalg.norm(c_new - c) / norm_c if norm_c > 0 else 0.0
            c = c_new
            if delta < tol:
                break

        # Final a with converged c
        with np.errstate(invalid="ignore", divide="ignore"):
            denom_a = (w * c[None, :] ** 2).sum(axis=1)
            numer_a = (w * r * c[None, :]).sum(axis=1)
            a = np.where(denom_a > 0, numer_a / denom_a, 0.0)

        # Subtract this component, then re-mask so the missing-data positions
        # stay at zero (the next iteration's diagnostics and any exposed
        # residuals reflect only the observed cells).
        r -= np.outer(a, c)
        r = np.where(observed, r, 0.0)
        rms_after = _rms(r)

        result.components.append(SysremComponent(
            c=c.copy(), a=a.copy(),
            rms_before=rms_before,
            rms_after=rms_after,
            iteration=it + 1,
        ))

    result.residuals = np.where(observed, r, np.nan)
    return result


def apply_sysrem(
    target_mag: np.ndarray,
    target_err: Optional[np.ndarray],
    result: SysremResult,
    n_components: Optional[int] = None,
) -> np.ndarray:
    """Apply pre-computed SYSREM components to a target light curve.

    The amplitude ``a_t`` for the target is estimated from the target data
    and each component ``c_j``.  This avoids using the target in the
    component extraction (which would bias towards removing intrinsic signal).

    Parameters
    ----------
    target_mag : (n_frames,) differential magnitudes of the target.
    target_err : (n_frames,) photometric errors, or None.
    result : SysremResult from ``sysrem()`` on the comparison stars.
    n_components : how many components to apply (default: all).

    Returns
    -------
    (n_frames,) corrected differential magnitudes.
    """
    mag = np.array(target_mag, dtype=float)
    if target_err is not None:
        err = np.array(target_err, dtype=float)
        w = np.where(np.isfinite(err) & (err > 0), 1.0 / err**2, 0.0)
    else:
        w = np.where(np.isfinite(mag), 1.0, 0.0)
    w = np.where(np.isfinite(mag), w, 0.0)

    # Remove per-target weighted mean
    sw = w.sum()
    if sw > 0:
        mean_t = float(np.nansum(w * mag)) / sw
        r = np.where(np.isfinite(mag), mag - mean_t, 0.0)
    else:
        r = np.where(np.isfinite(mag), mag, 0.0)
        mean_t = 0.0

    components = result.components
    if n_components is not None:
        components = components[:n_components]

    for comp in components:
        c = comp.c
        denom = float(np.sum(w * c**2))
        if denom > 0:
            a_t = float(np.sum(w * r * c)) / denom
            r -= a_t * c

    # Restore mean
    corrected = np.where(np.isfinite(target_mag), r + mean_t, np.nan)
    return corrected


def sysrem_dataframe(
    comp_df: "pd.DataFrame",
    target_series: "pd.Series",
    time_col: str = "time_id",
    mag_col: str = "mag_inst",
    err_col: str = "err",
    n_iter: int = 5,
    n_apply: Optional[int] = None,
) -> tuple["np.ndarray", "SysremResult"]:
    """High-level interface: run SYSREM on a pandas comparison DataFrame.

    Parameters
    ----------
    comp_df : long-format DataFrame with columns [time_id, star_id, mag, err].
    target_series : Series indexed by time_id with target differential magnitudes.
    time_col / mag_col / err_col : column names.
    n_iter : number of SYSREM iterations.
    n_apply : number of components to apply to target (default: n_iter).

    Returns
    -------
    (corrected_target, SysremResult)
    """
    import pandas as pd

    frames = sorted(comp_df[time_col].unique())
    stars  = sorted(comp_df["star_id"].unique())
    n_f = len(frames)
    n_s = len(stars)

    frame_idx = {f: i for i, f in enumerate(frames)}
    star_idx  = {s: i for i, s in enumerate(stars)}

    # Build matrices via pivot (much faster than iterrows for large DataFrames)
    mag_wide = comp_df.pivot_table(
        index="star_id", columns=time_col, values=mag_col, aggfunc="first"
    ).reindex(index=stars, columns=frames)
    mag_mat = mag_wide.to_numpy(float)

    if err_col in comp_df.columns:
        err_wide = comp_df.pivot_table(
            index="star_id", columns=time_col, values=err_col, aggfunc="first"
        ).reindex(index=stars, columns=frames)
        err_mat = err_wide.to_numpy(float)
    else:
        err_mat = np.full_like(mag_mat, np.nan)

    result = sysrem(mag_mat, err_mat, n_iter=n_iter)

    # Align target series to the common frame order
    tgt = np.array([float(target_series.get(f, np.nan)) for f in frames])

    corrected = apply_sysrem(tgt, None, result, n_components=n_apply)
    return corrected, result
