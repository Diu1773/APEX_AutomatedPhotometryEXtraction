"""Global ensemble photometry (inhomogeneous ensemble) solver.

This implements method C: solve per-frame zeropoints (Z_t) and per-star
means (M_i) simultaneously with weighted least squares:
    mag_inst(i,t,f) = M_i,f + Z_t,f + eps

The solution is per-filter by default. It supports iterative sigma-clipping
and comp-star rejection, then produces:
  - zp_df: per-frame zeropoints
  - mean_df: per-star mean instrumental mags
  - lc_df: target light curve corrected by Z_t
  - diagnostics: residual stats, removed comps, outliers
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsqr

from apex.utils.common_helpers import normalize_filter_key
from apex.utils.constants import MAD_TO_SIGMA

REQUIRED_COLS = {"time_id", "jd", "filter", "star_id", "mag_inst", "err"}


def comparison_network_diagnostics(df: pd.DataFrame) -> dict:
    """Describe connected components in the comparison-star/frame graph."""
    if df is None or df.empty:
        return {"connected": False, "n_components": 0, "components": []}

    work = df.dropna(subset=["time_id", "star_id"]).copy()
    if work.empty:
        return {"connected": False, "n_components": 0, "components": []}
    work["time_id"] = work["time_id"].astype(str)
    work["star_id"] = pd.to_numeric(work["star_id"], errors="coerce")
    work = work.dropna(subset=["star_id"])
    work["star_id"] = work["star_id"].astype(int)

    frame_to_stars = {
        str(frame): set(group["star_id"].astype(int).tolist())
        for frame, group in work.groupby("time_id")
    }
    star_to_frames = {
        int(star): set(group["time_id"].astype(str).tolist())
        for star, group in work.groupby("star_id")
    }
    unseen_frames = set(frame_to_stars)
    components = []
    while unseen_frames:
        pending_frames = [unseen_frames.pop()]
        component_frames: set[str] = set()
        component_stars: set[int] = set()
        while pending_frames:
            frame = pending_frames.pop()
            if frame in component_frames:
                continue
            component_frames.add(frame)
            for star in frame_to_stars.get(frame, set()):
                if star in component_stars:
                    continue
                component_stars.add(star)
                for linked_frame in star_to_frames.get(star, set()):
                    if linked_frame not in component_frames:
                        unseen_frames.discard(linked_frame)
                        pending_frames.append(linked_frame)
        n_obs = int(work[work["time_id"].isin(component_frames)].shape[0])
        components.append(
            {
                "n_frames": int(len(component_frames)),
                "n_stars": int(len(component_stars)),
                "n_observations": n_obs,
                "frames": sorted(component_frames),
                "star_ids": sorted(component_stars),
            }
        )
    components.sort(key=lambda item: item["n_frames"], reverse=True)
    return {
        "connected": len(components) == 1,
        "n_components": int(len(components)),
        "components": components,
    }


def _validate_comparison_design(df: pd.DataFrame) -> dict:
    network = comparison_network_diagnostics(df)
    if not network["connected"]:
        sizes = ", ".join(
            f"{item['n_frames']} frames/{item['n_stars']} stars"
            for item in network["components"]
        )
        raise ValueError(
            "Comparison network is disconnected; multi-night zeropoints are not "
            f"identifiable ({sizes}). Include at least one chain of shared comparison stars."
        )
    n_stars = int(df["star_id"].nunique())
    n_frames = int(df["time_id"].nunique())
    n_params = n_stars + n_frames - 1
    n_obs = int(len(df))
    if n_obs <= n_params:
        raise ValueError(
            f"Underconstrained ensemble fit: n_obs={n_obs}, n_params={n_params}."
        )
    network["n_parameters"] = n_params
    network["n_observations"] = n_obs
    return network


def select_comparisons_from_qc(
    qc_df: pd.DataFrame,
    rms_max: float,
    outlier_frac_max: float,
    min_points: int,
    min_coverage: float = 0.8,
    min_count: int = 3,
) -> dict:
    """Select stable comparisons, with a relative fallback for noisy fields."""
    if qc_df is None or qc_df.empty or "comp_id" not in qc_df.columns:
        return {"selected_ids": [], "method": "none", "scores": []}
    work = qc_df.copy()
    for column in (
        "comp_id",
        "n",
        "rms",
        "sigma_nights",
        "outlier_frac",
        "coverage_fraction",
    ):
        if column not in work.columns:
            work[column] = np.nan
        work[column] = pd.to_numeric(work[column], errors="coerce")
    work = work[np.isfinite(work["comp_id"])].copy()
    if work.empty:
        return {"selected_ids": [], "method": "none", "scores": []}
    work["comp_id"] = work["comp_id"].astype(int)
    work["coverage_fraction"] = work["coverage_fraction"].fillna(1.0)

    eligible = (
        (work["n"] >= int(min_points))
        & (work["coverage_fraction"] >= float(min_coverage))
        & np.isfinite(work["rms"])
    )
    limit = float(rms_max) if float(rms_max) > 0 else np.inf
    night_stable = ~np.isfinite(work["sigma_nights"]) | (
        work["sigma_nights"] <= limit
    )
    absolute = (
        eligible
        & (work["rms"] <= limit)
        & (work["outlier_frac"] <= float(outlier_frac_max))
        & night_stable
    )
    absolute_ids = work.loc[absolute, "comp_id"].astype(int).tolist()
    if len(absolute_ids) >= min(int(min_count), int(np.sum(eligible))):
        selected = sorted(set(absolute_ids))
        return {
            "selected_ids": selected,
            "method": "absolute",
            "scores": work[["comp_id", "rms", "sigma_nights", "outlier_frac", "coverage_fraction"]].to_dict("records"),
        }

    pool = work[eligible].copy()
    if pool.empty:
        return {"selected_ids": [], "method": "none", "scores": []}

    def _robust_upper(values: pd.Series, floor: float | None = None) -> float:
        numeric = pd.to_numeric(values, errors="coerce").to_numpy(float)
        numeric = numeric[np.isfinite(numeric)]
        if not len(numeric):
            return np.inf
        center = float(np.nanmedian(numeric))
        spread = _robust_sigma(numeric)
        upper = center + 3.0 * spread if np.isfinite(spread) else np.inf
        return max(float(floor), upper) if floor is not None else upper

    guard_specs = [
        (
            "outlier_frac",
            _robust_upper(pool["outlier_frac"], 2.0 * float(outlier_frac_max)),
            False,
        ),
        ("sigma_nights", _robust_upper(pool["sigma_nights"]), True),
        ("rms", _robust_upper(pool["rms"]), False),
    ]
    for column, upper, allow_missing in guard_specs:
        values = pd.to_numeric(pool[column], errors="coerce")
        keep = values <= upper
        if allow_missing:
            keep |= ~np.isfinite(values)
        guarded = pool[keep].copy()
        if len(guarded) >= min(int(min_count), len(pool)):
            pool = guarded

    def _rank_metric(values: pd.Series) -> pd.Series:
        numeric = pd.to_numeric(values, errors="coerce")
        worst = float(np.nanmax(numeric.to_numpy(float))) if np.any(np.isfinite(numeric)) else 1.0
        numeric = numeric.fillna(worst + max(abs(worst) * 0.1, 1e-6))
        return numeric.rank(method="average", pct=True)

    pool["qc_score"] = (
        0.45 * _rank_metric(pool["rms"])
        + 0.30 * _rank_metric(pool["sigma_nights"])
        + 0.20 * _rank_metric(pool["outlier_frac"])
        + 0.05 * (1.0 - pool["coverage_fraction"].clip(0.0, 1.0))
    )
    pool.loc[pool["comp_id"].isin(absolute_ids), "qc_score"] -= 1.0
    target_count = min(
        len(pool),
        max(int(min_count), int(np.ceil(len(pool) / 2.0))),
    )
    selected = (
        pool.sort_values(["qc_score", "rms", "comp_id"])
        .head(target_count)["comp_id"]
        .astype(int)
        .tolist()
    )
    return {
        "selected_ids": sorted(set(selected)),
        "method": "relative_fallback",
        "scores": pool[
            [
                "comp_id",
                "qc_score",
                "rms",
                "sigma_nights",
                "outlier_frac",
                "coverage_fraction",
            ]
        ].to_dict("records"),
    }


def solve_global_ensemble(
    df: pd.DataFrame,
    target_id: int,
    comp_ids: Iterable[int],
    min_comps: int = 3,
    sigma: float = 3.0,
    n_iter: int = 3,
    gauge: str = "meanZ0",
    per_filter: bool = True,
    robust: bool = True,
    rms_clip_pct: float = 10.0,
    rms_clip_threshold: float | None = None,
    frame_sigma: float = 3.0,
    interp_missing: bool = False,
    normalize_target: bool = False,
    max_dense_params: int = 2000,
    rescale_errors: bool = True,
    log: Optional[Callable[[str], None]] = None,
) -> Dict[str, pd.DataFrame | dict]:
    """Solve global ensemble photometry (method C).

    Args:
        df: DataFrame with REQUIRED_COLS.
        target_id: Target star ID (excluded from solving).
        comp_ids: Comparison star IDs.
        min_comps: Minimum comps required per frame to solve Z_t.
        sigma: Sigma clip threshold.
        n_iter: Iteration count for comp/outlier rejection.
        gauge: "meanZ0" (default) or "ref".
        per_filter: Solve independently per filter.
        robust: Use MAD-based sigma for clipping.
        rms_clip_pct: Maximum fraction of statistically high-RMS comps removed
            per iteration. Stable comps are not removed to fill this quota.
        rms_clip_threshold: Absolute RMS cutoff for comps (optional).
        frame_sigma: Frame-level outlier sigma threshold.
        interp_missing: Interpolate Z_t for frames with too few comps.
        normalize_target: Subtract median from corrected target curve.
        max_dense_params: Dense covariance threshold.
        rescale_errors: Inflate per-frame errors based on chi2_red (floor at 1.0).
        log: Optional logger callback.
    """

    _log = log or (lambda _: None)
    if df is None or df.empty:
        raise ValueError("Input DataFrame is empty")
    missing = REQUIRED_COLS - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    comp_ids = [int(c) for c in comp_ids if str(c).strip()]
    comp_ids = [c for c in comp_ids if c != int(target_id)]
    if not comp_ids:
        raise ValueError("comp_ids is empty after removing target_id")

    data = df.copy()
    data["time_id"] = data["time_id"].astype(str)
    # Canonical case: Johnson uppercase (B/V/R/I), SDSS lowercase (g/r/i/z) — preserved.
    data["filter"] = data["filter"].astype(str).map(normalize_filter_key)
    data["star_id"] = pd.to_numeric(data["star_id"], errors="coerce").astype("Int64")
    data["mag_inst"] = pd.to_numeric(data["mag_inst"], errors="coerce")
    data["err"] = pd.to_numeric(data["err"], errors="coerce")
    data = data.dropna(subset=["time_id", "filter", "star_id", "mag_inst"])
    data = data[data["star_id"].notna()].copy()
    data["star_id"] = data["star_id"].astype(int)
    data["error_imputed"] = ~(np.isfinite(data["err"]) & (data["err"] > 0))
    for fkey, indices in data.groupby("filter").groups.items():
        values = data.loc[indices, "err"].to_numpy(float)
        valid = values[np.isfinite(values) & (values > 0)]
        fallback = float(np.nanmedian(valid)) if len(valid) else 1.0
        bad = data.loc[indices, "error_imputed"]
        data.loc[np.asarray(indices)[bad.to_numpy(bool)], "err"] = fallback
        n_bad = int(np.sum(bad))
        if n_bad:
            _log(
                f"[GLOBAL] filter={fkey}: imputed {n_bad} invalid uncertainties "
                f"with {fallback:.6g} mag."
            )

    if per_filter:
        filters = sorted(set(data["filter"].unique()))
    else:
        filters = [""]
        data["filter"] = ""

    zp_frames = []
    mean_frames = []
    lc_frames = []
    diagnostics: dict = {
        "filters": {},
        "removed_comps": [],
        "removed_comp_details": [],
        "outliers": [],
    }

    for fkey in filters:
        sub = data[data["filter"] == fkey].copy()
        if sub.empty:
            continue
        _log(f"[GLOBAL] filter={fkey} rows={len(sub)} comps={len(comp_ids)}")

        result = _solve_one_filter(
            sub,
            target_id=target_id,
            comp_ids=comp_ids,
            min_comps=min_comps,
            sigma=sigma,
            n_iter=n_iter,
            gauge=gauge,
            robust=robust,
            rms_clip_pct=rms_clip_pct,
            rms_clip_threshold=rms_clip_threshold,
            frame_sigma=frame_sigma,
            interp_missing=interp_missing,
            normalize_target=normalize_target,
            max_dense_params=max_dense_params,
            rescale_errors=rescale_errors,
            log=_log,
        )

        zp = result["zp_df"]
        mean = result["mean_df"]
        lc = result["lc_df"]
        diag = result["diagnostics"]

        zp["filter"] = fkey
        mean["filter"] = fkey
        lc["filter"] = fkey

        zp_frames.append(zp)
        mean_frames.append(mean)
        lc_frames.append(lc)
        diagnostics["filters"][fkey] = diag
        diagnostics["removed_comps"].extend(diag.get("removed_comps", []))
        diagnostics["removed_comp_details"].extend(
            diag.get("removed_comp_details", [])
        )
        diagnostics["outliers"].extend(diag.get("outliers", []))

    zp_df = pd.concat(zp_frames, ignore_index=True) if zp_frames else pd.DataFrame()
    mean_df = pd.concat(mean_frames, ignore_index=True) if mean_frames else pd.DataFrame()
    lc_df = pd.concat(lc_frames, ignore_index=True) if lc_frames else pd.DataFrame()

    return {
        "zp_df": zp_df,
        "mean_df": mean_df,
        "lc_df": lc_df,
        "diagnostics": diagnostics,
    }


def _solve_one_filter(
    df: pd.DataFrame,
    target_id: int,
    comp_ids: List[int],
    min_comps: int,
    sigma: float,
    n_iter: int,
    gauge: str,
    robust: bool,
    rms_clip_pct: float,
    rms_clip_threshold: float | None,
    frame_sigma: float,
    interp_missing: bool,
    normalize_target: bool,
    max_dense_params: int,
    rescale_errors: bool,
    log: Callable[[str], None],
) -> Dict[str, pd.DataFrame | dict]:
    comp_active = [c for c in comp_ids if c != int(target_id)]
    removed_comps: List[int] = []
    removed_comp_details: List[dict] = []
    outliers: List[dict] = []
    suspect_frames: List[dict] = []
    diag_iters = []

    target_df = df[df["star_id"] == int(target_id)].copy()
    comp_df_all = df[df["star_id"].isin(comp_active)].copy()

    if comp_df_all.empty:
        raise ValueError("No comparison star measurements found")

    # Track measurement-level removals
    comp_df = comp_df_all.copy()
    comp_df["keep"] = True

    for it in range(max(1, n_iter)):
        comp_df = comp_df[comp_df["star_id"].isin(comp_active)].copy()
        if comp_df.empty:
            break

        comp_df = _drop_frames_with_few_comps(comp_df, min_comps)
        if comp_df.empty:
            break

        fit = _solve_wls(comp_df, gauge=gauge, max_dense_params=max_dense_params, log=log)
        resid = comp_df["mag_inst"].to_numpy(float) - fit["model"]
        comp_df["resid"] = resid

        # Chi²_red error rescaling (conservative: inflate only, floor at 1.0)
        if rescale_errors and it < max(1, n_iter) - 1:
            if "err_orig" not in comp_df.columns:
                comp_df["err_orig"] = comp_df["err"].copy()
            err_cur = comp_df["err_orig"].to_numpy(float)
            _ok = np.isfinite(resid) & np.isfinite(err_cur) & (err_cur > 0)
            if np.any(_ok):
                def _chi2_red_frame(g):
                    r = g["resid"].to_numpy(float)
                    e = g["err_orig"].to_numpy(float)
                    v = np.isfinite(r) & np.isfinite(e) & (e > 0)
                    if np.sum(v) < 2:
                        return 1.0
                    return float(np.sum((r[v] / e[v]) ** 2) / max(int(np.sum(v)) - 1, 1))
                chi2_map = comp_df.groupby("time_id")[["resid", "err_orig"]].apply(
                    _chi2_red_frame
                )
                scale = np.sqrt(np.clip(
                    comp_df["time_id"].map(chi2_map).to_numpy(float), 1.0, np.inf
                ))
                comp_df["err"] = err_cur * scale
                n_inflated = len(
                    set(comp_df.loc[scale > 1.0, "time_id"].astype(str).tolist())
                )
                if n_inflated > 0:
                    log(f"[RESCALE] Iter {it+1}: {n_inflated}/{len(chi2_map)} frames inflated "
                        f"(max scale={float(np.nanmax(scale)):.2f})")

        # Outlier rejection (measurement-level)
        sigma_global = _robust_sigma(resid) if robust else np.nanstd(resid)
        if not np.isfinite(sigma_global) or sigma_global <= 0:
            sigma_global = np.nanstd(resid) if np.isfinite(np.nanstd(resid)) else 0.0

        frame_sigmas = _frame_sigmas(comp_df, sigma_global, robust)
        comp_df["resid_scale_frame"] = (
            comp_df["time_id"].astype(str).map(frame_sigmas).fillna(sigma_global)
        )
        local_scale = np.maximum(
            comp_df["resid_scale_frame"].to_numpy(float),
            comp_df["err"].to_numpy(float),
        )
        outlier_mask = (
            np.isfinite(local_scale)
            & (local_scale > 0)
            & (np.abs(comp_df["resid"].to_numpy(float)) > float(sigma) * local_scale)
        )
        comp_df["outlier"] = outlier_mask
        out_rows = comp_df[comp_df["outlier"]].copy()
        if not out_rows.empty:
            trial = comp_df[~comp_df["outlier"]].copy()
            try:
                _validate_comparison_design(trial)
            except ValueError:
                log(
                    f"[GLOBAL] Iter {it + 1}: retained {len(out_rows)} clipped "
                    "measurements because removing them disconnects the ensemble."
                )
                out_rows = out_rows.iloc[0:0]
            for _, r in out_rows.iterrows():
                outliers.append(
                    dict(
                        time_id=str(r["time_id"]),
                        star_id=int(r["star_id"]),
                        resid=float(r["resid"]),
                        reason="measurement_sigma",
                    )
                )
            if not out_rows.empty:
                comp_df = trial

        frame_score = pd.Series(
            {
                str(frame): float(
                    np.sqrt(
                        np.mean(
                            np.square(
                                group["resid"].to_numpy(float)
                                / group["err"].to_numpy(float)
                            )
                        )
                    )
                )
                for frame, group in comp_df.groupby("time_id")
            },
            dtype=float,
        )
        bad_frames: list[str] = []
        if len(frame_score) >= 5:
            frame_center = float(np.nanmedian(frame_score.to_numpy(float)))
            frame_spread = _robust_sigma(frame_score.to_numpy(float))
            if np.isfinite(frame_spread) and frame_spread > 0:
                threshold = frame_center + float(frame_sigma) * frame_spread
                bad_frames = frame_score[frame_score > threshold].index.tolist()
        if bad_frames:
            suspect_frames.extend(
                {"time_id": frame, "reason": "frame_scatter"}
                for frame in bad_frames
            )

        # Comparison RMS and rejection
        rms_by_star = comp_df.groupby("star_id")["resid"].apply(
            lambda values: float(np.sqrt(np.mean(np.square(values))))
        )
        rms_by_star = rms_by_star.dropna()
        drop_reasons: dict[int, set[str]] = {}
        if rms_clip_threshold is not None:
            for star_id in rms_by_star[
                rms_by_star > float(rms_clip_threshold)
            ].index:
                drop_reasons.setdefault(int(star_id), set()).add("absolute_rms")
        if rms_clip_pct and rms_clip_pct > 0 and len(rms_by_star) > min_comps:
            rms_values = rms_by_star.to_numpy(float)
            rms_center = float(np.nanmedian(rms_values))
            rms_spread = _robust_sigma(rms_values)
            if np.isfinite(rms_spread) and rms_spread > 0:
                statistical = rms_by_star[
                    rms_by_star > rms_center + float(sigma) * rms_spread
                ].sort_values(ascending=False)
                max_relative = int(
                    np.ceil(len(rms_by_star) * float(rms_clip_pct) / 100.0)
                )
                max_relative = max(
                    0, min(len(rms_by_star) - min_comps, max_relative)
                )
                for star_id in statistical.head(max_relative).index:
                    drop_reasons.setdefault(int(star_id), set()).add("relative_rms")
        drop_ids = sorted(
            drop_reasons,
            key=lambda star_id: float(rms_by_star.get(star_id, -np.inf)),
            reverse=True,
        )[: max(len(comp_active) - min_comps, 0)]

        if drop_ids:
            trial_active = [c for c in comp_active if c not in drop_ids]
            trial = comp_df[comp_df["star_id"].isin(trial_active)].copy()
            try:
                _validate_comparison_design(trial)
            except ValueError:
                log(
                    f"[GLOBAL] Iter {it + 1}: retained {len(drop_ids)} noisy "
                    "comparisons because removing them disconnects the ensemble."
                )
                drop_ids = []
            if drop_ids:
                removed_comps.extend(drop_ids)
                for star_id in drop_ids:
                    removed_comp_details.append(
                        {
                            "star_id": int(star_id),
                            "iteration": int(it + 1),
                            "rms": float(rms_by_star.get(star_id, np.nan)),
                            "reasons": sorted(drop_reasons.get(star_id, set())),
                        }
                    )
                comp_active = trial_active

        diag_iters.append(
            dict(
                iter=it + 1,
                n_comp=len(comp_active),
                n_obs=len(comp_df),
                n_outliers=int(len(out_rows)),
                n_frames_flagged=int(len(bad_frames)),
                n_comps_removed=int(len(drop_ids)),
            )
        )

        if len(comp_active) <= min_comps:
            break

    # Final solve
    comp_df = comp_df[comp_df["star_id"].isin(comp_active)].copy()
    comp_df = _drop_frames_with_few_comps(comp_df, min_comps)
    if comp_df.empty:
        raise ValueError("No usable comp measurements after clipping")

    fit = _solve_wls(comp_df, gauge=gauge, max_dense_params=max_dense_params, log=log)
    comp_df["resid"] = comp_df["mag_inst"].to_numpy(float) - fit["model"]

    # ZP table
    zp_df = _build_zp_df(comp_df, fit, min_comps=min_comps)
    if interp_missing:
        zp_df = _interp_missing_zp(zp_df)

    # Mean table
    mean_df = _build_mean_df(comp_df, fit)

    # Target light curve
    lc_df = _build_target_lc(
        target_df,
        comp_df,
        zp_df,
        mean_df,
        min_comps=min_comps,
        normalize_target=normalize_target,
    )

    diagnostics = dict(
        removed_comps=removed_comps,
        removed_comp_details=removed_comp_details,
        outliers=outliers,
        suspect_frames=suspect_frames,
        iterations=diag_iters,
        n_comp_final=len(comp_active),
        n_obs_final=len(comp_df),
        n_errors_imputed=int(
            comp_df_all.get("error_imputed", pd.Series(dtype=bool)).sum()
        ),
        network=fit["network"],
        comparison_stats=_build_comparison_stats(
            comp_df_all, comp_df, removed_comp_details
        ),
    )

    return {"zp_df": zp_df, "mean_df": mean_df, "lc_df": lc_df, "diagnostics": diagnostics}


def _build_comparison_stats(
    original: pd.DataFrame,
    final: pd.DataFrame,
    removed_details: List[dict],
) -> List[dict]:
    total_frames = max(int(original["time_id"].nunique()), 1)
    removed_map = {
        int(item["star_id"]): item for item in removed_details if "star_id" in item
    }
    rows = []
    for star_id, group in original.groupby("star_id"):
        final_group = final[final["star_id"] == int(star_id)]
        residual = pd.to_numeric(
            final_group.get("resid", pd.Series(dtype=float)), errors="coerce"
        ).to_numpy(float)
        residual = residual[np.isfinite(residual)]
        rms = float(np.sqrt(np.mean(np.square(residual)))) if len(residual) else None
        removed = removed_map.get(int(star_id), {})
        rows.append(
            {
                "star_id": int(star_id),
                "active": not final_group.empty,
                "n_observations": int(len(group)),
                "n_frames": int(group["time_id"].nunique()),
                "coverage_fraction": float(group["time_id"].nunique() / total_frames),
                "final_rms": rms,
                "error_imputed_count": int(
                    group.get("error_imputed", pd.Series(dtype=bool)).sum()
                ),
                "removed_reasons": list(removed.get("reasons", [])),
            }
        )
    return sorted(rows, key=lambda item: item["star_id"])


def _solve_wls(
    df: pd.DataFrame,
    gauge: str,
    max_dense_params: int,
    log: Callable[[str], None],
) -> dict:
    network = _validate_comparison_design(df)
    # Map IDs
    star_ids = sorted(df["star_id"].unique().tolist())
    time_ids = sorted(df["time_id"].unique().tolist())
    star_map = {sid: i for i, sid in enumerate(star_ids)}
    time_map = {tid: i for i, tid in enumerate(time_ids)}

    ref_time = time_ids[0]
    n_star = len(star_ids)
    n_time = len(time_ids)
    n_time_params = n_time - 1  # ref Z fixed to 0

    n_params = n_star + n_time_params
    n_obs = len(df)

    # Build sparse A
    row = np.arange(n_obs, dtype=int)
    cols = np.empty(n_obs * 2, dtype=int)
    data = np.ones(n_obs * 2, dtype=float)

    star_idx = df["star_id"].map(star_map).to_numpy(int)
    time_idx = df["time_id"].map(time_map).to_numpy(int)
    time_col = np.where(time_idx == 0, -1, time_idx - 1)  # ref time removed

    cols[:n_obs] = star_idx
    cols[n_obs:] = n_star + time_col
    keep = cols[n_obs:] >= n_star
    cols = np.concatenate([cols[:n_obs], cols[n_obs:][keep]])
    data = np.concatenate([data[:n_obs], data[n_obs:][keep]])
    row = np.concatenate([row, row[keep]])

    A = csr_matrix((data, (row, cols)), shape=(n_obs, n_params))

    y = df["mag_inst"].to_numpy(float)
    err = df["err"].to_numpy(float)
    w = np.where(np.isfinite(err) & (err > 0), 1.0 / (err * err), 1.0)

    # Weighted LS via lsqr
    Aw = A.multiply(np.sqrt(w)[:, None])
    yw = y * np.sqrt(w)
    sol = lsqr(Aw, yw, atol=1e-10, btol=1e-10, iter_lim=2000)
    if int(sol[1]) not in {1, 2}:
        raise ValueError(
            f"Global ensemble LSQR did not converge (istop={int(sol[1])})."
        )
    x = sol[0]

    # Build M and Z
    M = x[:n_star].copy()
    Z = np.zeros(n_time, dtype=float)
    Z[0] = 0.0
    if n_time_params > 0:
        Z[1:] = x[n_star:]

    gauge_coeff = np.zeros(n_params, dtype=float)
    if gauge.lower() == "meanz0":
        weights = df.groupby("time_id")["mag_inst"].count().reindex(time_ids).fillna(1).to_numpy(float)
        wsum = np.sum(weights)
        if wsum > 0:
            if n_time_params > 0:
                gauge_coeff[n_star:] = weights[1:] / wsum
            z_mean = float(np.sum(Z * weights) / wsum)
            Z = Z - z_mean
            M = M + z_mean

    model = M[star_idx] + Z[time_idx]

    # Error estimates
    approx = False
    M_err = np.full(n_star, np.nan)
    Z_err = np.full(n_time, np.nan)
    try:
        if n_params <= max_dense_params:
            AtW = A.T.multiply(w)
            AtWA = (AtW @ A).toarray()
            cov = np.linalg.pinv(AtWA)
            diag = np.diag(cov)
            if gauge.lower() == "meanz0":
                cov_with_mean = cov @ gauge_coeff
                var_mean = float(gauge_coeff @ cov_with_mean)
                m_diag = (
                    diag[:n_star]
                    + var_mean
                    + 2.0 * cov_with_mean[:n_star]
                )
                z_diag = np.full(n_time, var_mean, dtype=float)
                if n_time_params > 0:
                    z_diag[1:] = (
                        diag[n_star:]
                        + var_mean
                        - 2.0 * cov_with_mean[n_star:]
                    )
            else:
                m_diag = diag[:n_star]
                z_diag = np.zeros(n_time, dtype=float)
                if n_time_params > 0:
                    z_diag[1:] = diag[n_star:]
            M_err = np.sqrt(np.clip(m_diag, 0, np.inf))
            Z_err = np.sqrt(np.clip(z_diag, 0, np.inf))
        else:
            AtW = A.T.multiply(w)
            AtWA = (AtW @ A).diagonal()
            diag = np.where(AtWA > 0, 1.0 / AtWA, np.nan)
            if gauge.lower() == "meanz0":
                var_mean = float(np.nansum(np.square(gauge_coeff) * diag))
                m_diag = diag[:n_star] + var_mean
                z_diag = np.full(n_time, var_mean, dtype=float)
                if n_time_params > 0:
                    z_diag[1:] = diag[n_star:] + var_mean
            else:
                m_diag = diag[:n_star]
                z_diag = np.zeros(n_time, dtype=float)
                if n_time_params > 0:
                    z_diag[1:] = diag[n_star:]
            M_err = np.sqrt(np.clip(m_diag, 0, np.inf))
            Z_err = np.sqrt(np.clip(z_diag, 0, np.inf))
            approx = True
    except Exception as e:
        approx = True
        log(f"[GLOBAL] Warning: covariance failed ({e})")

    return dict(
        M=M,
        Z=Z,
        M_err=M_err,
        Z_err=Z_err,
        model=model,
        star_ids=star_ids,
        time_ids=time_ids,
        approx_errors=approx,
        network=network,
    )


def _build_zp_df(comp_df: pd.DataFrame, fit: dict, min_comps: int) -> pd.DataFrame:
    time_ids = fit["time_ids"]
    Z = fit["Z"]
    Z_err = fit["Z_err"]

    stats = comp_df.groupby("time_id")[["resid", "err"]].apply(_frame_stats)
    stats = stats.reindex(time_ids)
    n_used = stats.get("n_used", pd.Series([0] * len(time_ids))).to_numpy(int)
    chi2_red = stats.get("chi2_red", pd.Series([np.nan] * len(time_ids))).to_numpy(float)

    zp_df = pd.DataFrame(
        dict(
            time_id=time_ids,
            Z=Z,
            Z_err=Z_err,
            n_used=n_used,
            chi2_red=chi2_red,
        )
    )
    zp_df.loc[zp_df["n_used"] < min_comps, ["Z", "Z_err"]] = np.nan
    return zp_df


def _build_mean_df(comp_df: pd.DataFrame, fit: dict) -> pd.DataFrame:
    star_ids = fit["star_ids"]
    M = fit["M"]
    M_err = fit["M_err"]
    counts = comp_df.groupby("star_id")["mag_inst"].count().reindex(star_ids).fillna(0).to_numpy(int)
    return pd.DataFrame(
        dict(
            star_id=star_ids,
            M=M,
            M_err=M_err,
            n_used=counts,
        )
    )


def _build_target_lc(
    target_df: pd.DataFrame,
    comp_df: pd.DataFrame,
    zp_df: pd.DataFrame,
    mean_df: pd.DataFrame,
    min_comps: int,
    normalize_target: bool,
) -> pd.DataFrame:
    if target_df.empty:
        return pd.DataFrame()

    def _weighted_mean(g):
        m = g["mag_inst"].to_numpy(float)
        e = g["err"].to_numpy(float)
        ok = np.isfinite(m) & np.isfinite(e) & (e > 0)
        if np.any(ok):
            w = 1.0 / (e[ok] ** 2)
            return float(np.sum(m[ok] * w) / np.sum(w))
        return float(np.nanmean(m))
    comp_mean = comp_df.groupby("time_id")[["mag_inst", "err"]].apply(_weighted_mean)
    comp_n = comp_df.groupby("time_id")["mag_inst"].count()

    comp_ref = np.nan
    if mean_df is not None and not mean_df.empty:
        ref = mean_df.copy()
        ref["M"] = pd.to_numeric(ref.get("M"), errors="coerce")
        ref["M_err"] = pd.to_numeric(ref.get("M_err"), errors="coerce")
        ref["n_used"] = pd.to_numeric(ref.get("n_used"), errors="coerce")
        ref = ref[np.isfinite(ref["M"])].copy()
        if not ref.empty:
            w = np.where(
                np.isfinite(ref["M_err"]) & (ref["M_err"] > 0),
                1.0 / np.square(ref["M_err"].to_numpy(float)),
                ref["n_used"].fillna(1.0).to_numpy(float),
            )
            if np.any(np.isfinite(w) & (w > 0)):
                comp_ref = float(np.average(ref["M"].to_numpy(float), weights=np.where(np.isfinite(w) & (w > 0), w, 1.0)))
            else:
                comp_ref = float(np.nanmean(ref["M"].to_numpy(float)))

    lc = target_df.copy()
    lc = lc.rename(columns={"mag_inst": "mag"})
    lc["diff_mag_raw"] = lc["mag"] - lc["time_id"].map(comp_mean)
    lc.loc[lc["time_id"].map(comp_n) < min_comps, "diff_mag_raw"] = np.nan

    zp_map = zp_df.set_index("time_id")["Z"]
    zp_err_map = zp_df.set_index("time_id")["Z_err"]
    lc["mag_ensemble_corr"] = lc["mag"] - lc["time_id"].map(zp_map)
    lc["comp_ref_mean"] = comp_ref
    if np.isfinite(comp_ref):
        lc["diff_mag_corr"] = lc["mag_ensemble_corr"] - comp_ref
    else:
        lc["diff_mag_corr"] = lc["mag_ensemble_corr"]
    lc["diff_err_corr"] = np.sqrt(
        np.square(lc["err"].to_numpy(float)) + np.square(lc["time_id"].map(zp_err_map).to_numpy(float))
    )

    if normalize_target:
        med = np.nanmedian(lc["diff_mag_corr"].to_numpy(float))
        if np.isfinite(med):
            lc["diff_mag_corr"] = lc["diff_mag_corr"] - med

    lc = lc.rename(columns={"err": "diff_err"})
    core_cols = [
        "jd",
        "time_id",
        "filter",
        "star_id",
        "mag",
        "mag_ensemble_corr",
        "comp_ref_mean",
        "diff_mag_raw",
        "diff_mag_corr",
        "diff_err",
        "diff_err_corr",
    ]
    extra_cols = [c for c in lc.columns if c not in core_cols]
    return lc[[c for c in core_cols if c in lc.columns] + extra_cols]


def _drop_frames_with_few_comps(df: pd.DataFrame, min_comps: int) -> pd.DataFrame:
    counts = df.groupby("time_id")["star_id"].nunique()
    keep = counts[counts >= min_comps].index
    return df[df["time_id"].isin(keep)].copy()


def _frame_stats(group: pd.DataFrame) -> pd.Series:
    resid = group["resid"].to_numpy(float)
    err = group["err"].to_numpy(float)
    ok = np.isfinite(resid) & np.isfinite(err) & (err > 0)
    chi2 = np.sum((resid[ok] / err[ok]) ** 2) if np.any(ok) else np.nan
    dof = max(int(np.sum(ok)) - 1, 1)
    return pd.Series(
        dict(
            n_used=int(np.sum(ok)),   # measurements with finite resid AND positive err
            chi2_red=float(chi2 / dof) if np.isfinite(chi2) else np.nan,
        )
    )


def _frame_sigmas(df: pd.DataFrame, sigma_global: float, robust: bool) -> Dict[str, float]:
    sigmas: Dict[str, float] = {}
    for tid, sub in df.groupby("time_id"):
        vals = sub["resid"].to_numpy(float)
        if len(vals) >= 5:
            sig = _robust_sigma(vals) if robust else np.nanstd(vals)
            sig = sig if np.isfinite(sig) and sig > 0 else sigma_global
        else:
            sig = sigma_global
        sigmas[str(tid)] = float(sig)
    return sigmas


def _robust_sigma(x: np.ndarray) -> float:
    x = np.asarray(x, float)
    med = np.nanmedian(x)
    mad = np.nanmedian(np.abs(x - med))
    if np.isfinite(mad) and mad > 0:
        return float(MAD_TO_SIGMA * mad)
    return float(np.nanstd(x))


def _interp_missing_zp(zp_df: pd.DataFrame) -> pd.DataFrame:
    if zp_df.empty or "Z" not in zp_df.columns:
        return zp_df
    zp_df = zp_df.copy()
    z = pd.to_numeric(zp_df["Z"], errors="coerce")
    zp_df["Z"] = z.interpolate(limit_direction="both")
    return zp_df


def generate_synthetic_data(
    n_frames: int = 200,
    n_comps: int = 20,
    period: float = 0.35,
    noise: float = 0.01,
    seed: int = 42,
) -> Tuple[pd.DataFrame, int, List[int], dict]:
    """Create synthetic dataset for validation."""
    rng = np.random.default_rng(seed)
    time_ids = [f"f{i:04d}.fit" for i in range(n_frames)]
    jd = np.linspace(0.0, 2.0, n_frames) + 2450000.0
    filters = ["g"] * n_frames
    target_id = 1
    comp_ids = list(range(2, 2 + n_comps))

    # True per-frame ZP offsets (nightly jumps)
    Z_true = np.zeros(n_frames)
    for k in range(0, n_frames, 50):
        Z_true[k:k + 50] = rng.normal(0, 0.05)

    # Star means
    M_comp = rng.normal(12.0, 0.5, size=n_comps)
    M_target = 12.3
    target_signal = 0.2 * np.sin(2 * np.pi * (jd - jd.min()) / period)

    rows = []
    for t, tid in enumerate(time_ids):
        for i, sid in enumerate(comp_ids):
            mag = M_comp[i] + Z_true[t] + rng.normal(0, noise)
            rows.append(dict(time_id=tid, jd=jd[t], filter="g", star_id=sid, mag_inst=mag, err=noise))
        mag_t = M_target + Z_true[t] + target_signal[t] + rng.normal(0, noise)
        rows.append(dict(time_id=tid, jd=jd[t], filter="g", star_id=target_id, mag_inst=mag_t, err=noise))

    df = pd.DataFrame(rows)
    truth = dict(Z_true=Z_true, jd=jd)
    return df, target_id, comp_ids, truth


def run_synthetic_test() -> dict:
    """Basic verification: recover injected frame offsets."""
    df, target_id, comp_ids, truth = generate_synthetic_data()
    result = solve_global_ensemble(df, target_id, comp_ids, min_comps=5, n_iter=2)
    zp = result["zp_df"]
    z = pd.to_numeric(zp["Z"], errors="coerce").to_numpy(float)
    ok = np.isfinite(z)
    if np.any(ok):
        corr = np.corrcoef(z[ok], truth["Z_true"][ok])[0, 1]
    else:
        corr = np.nan
    return {"corr_Z": float(corr), "n_points": int(np.sum(ok))}


def example_usage() -> None:
    """Example usage snippet."""
    df, target_id, comp_ids, _ = generate_synthetic_data()
    result = solve_global_ensemble(df, target_id, comp_ids)
    print(result["zp_df"].head())
