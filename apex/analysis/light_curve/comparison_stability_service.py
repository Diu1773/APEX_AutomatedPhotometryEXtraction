"""Time-domain stability diagnostics for comparison-star selection."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from apex.utils.constants import MAD_TO_SIGMA


@dataclass(frozen=True)
class ComparisonSelectionConfig:
    min_coverage: float = 0.90
    min_points: int = 20
    min_ensemble: int = 2
    min_comparisons: int = 3
    target_count: int = 5
    suspect_score: float = 2.0
    reject_score: float = 4.0
    max_iterations: int = 20
    magnitude_weight: float = 0.40
    color_weight: float = 0.10


def _robust_sigma(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size < 2:
        return float("nan")
    center = float(np.nanmedian(finite))
    return float(MAD_TO_SIGMA * np.nanmedian(np.abs(finite - center)))


def _safe_positive_scale(values: np.ndarray, floor: float = 1e-4) -> float:
    scale = _robust_sigma(values)
    if np.isfinite(scale) and scale > floor:
        return float(scale)
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size >= 2:
        fallback = float(np.nanstd(finite, ddof=1))
        if np.isfinite(fallback) and fallback > floor:
            return fallback
    return float(floor)


def _von_neumann_eta(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size < 3:
        return float("nan")
    variance = float(np.var(finite, ddof=1))
    if not np.isfinite(variance) or variance <= 0:
        return float("nan")
    return float(np.sum(np.diff(finite) ** 2) / ((finite.size - 1) * variance))


def _max_run_fraction(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size < 2:
        return float("nan")
    centered = finite - float(np.nanmedian(finite))
    signs = np.sign(centered)
    signs = signs[signs != 0]
    if signs.size == 0:
        return 0.0
    longest = 1
    current = 1
    for previous, value in zip(signs[:-1], signs[1:]):
        if value == previous:
            current += 1
            longest = max(longest, current)
        else:
            current = 1
    return float(longest / signs.size)


def _condition_correlation(group: pd.DataFrame) -> tuple[float, str]:
    best = 0.0
    best_name = ""
    residual = pd.to_numeric(group["residual"], errors="coerce").to_numpy(float)
    for column in ("airmass", "fwhm", "sky", "x", "y"):
        if column not in group.columns:
            continue
        values = pd.to_numeric(group[column], errors="coerce").to_numpy(float)
        valid = np.isfinite(values) & np.isfinite(residual)
        if int(valid.sum()) < 8 or np.nanstd(values[valid]) <= 0:
            continue
        corr = spearmanr(values[valid], residual[valid], nan_policy="omit").statistic
        if np.isfinite(corr) and abs(float(corr)) > best:
            best = abs(float(corr))
            best_name = column
    return float(best), best_name


def compute_leave_one_out_residuals(
    measurements: pd.DataFrame,
    candidate_ids: Iterable[int],
    *,
    min_ensemble: int = 2,
) -> pd.DataFrame:
    """Build candidate-minus-ensemble residuals without self dilution."""
    required = {"frame", "star_id", "mag"}
    if measurements is None or measurements.empty or not required <= set(measurements.columns):
        return pd.DataFrame()

    ids = sorted({int(value) for value in candidate_ids})
    if not ids:
        return pd.DataFrame()
    work = measurements.copy()
    work["star_id"] = pd.to_numeric(work["star_id"], errors="coerce").astype("Int64")
    work["mag"] = pd.to_numeric(work["mag"], errors="coerce")
    work = work[work["star_id"].isin(ids) & np.isfinite(work["mag"])].copy()
    if work.empty:
        return pd.DataFrame()
    work["star_id"] = work["star_id"].astype("int64")
    work = work.sort_values(["frame", "star_id"]).drop_duplicates(
        subset=["frame", "star_id"], keep="first"
    )

    baselines = work.groupby("star_id")["mag"].median()
    work["centered_mag"] = work["mag"] - work["star_id"].map(baselines)
    centered = work.pivot(index="frame", columns="star_id", values="centered_mag")
    centered = centered.reindex(columns=ids)

    if "mag_err" in work.columns:
        err_work = work.copy()
        err_work["mag_err"] = pd.to_numeric(err_work["mag_err"], errors="coerce")
        errors = err_work.pivot(index="frame", columns="star_id", values="mag_err")
        errors = errors.reindex(index=centered.index, columns=ids)
    else:
        errors = pd.DataFrame(np.nan, index=centered.index, columns=ids)

    frame_meta_columns = [
        column
        for column in ("time", "night_id", "airmass", "fwhm", "sky")
        if column in work.columns
    ]
    frame_meta = (
        work.groupby("frame", as_index=True)[frame_meta_columns].first()
        if frame_meta_columns
        else pd.DataFrame(index=centered.index)
    )
    xy = work.set_index(["frame", "star_id"])[
        [column for column in ("x", "y") if column in work.columns]
    ]

    rows: list[pd.DataFrame] = []
    for star_id in ids:
        if star_id not in centered.columns:
            continue
        other_columns = [value for value in ids if value != star_id]
        if len(other_columns) < int(min_ensemble):
            continue
        others = centered[other_columns]
        ensemble_n = others.notna().sum(axis=1)
        reference = others.median(axis=1, skipna=True)
        reference = reference.where(ensemble_n >= int(min_ensemble))

        star_centered = centered[star_id]
        residual = star_centered - reference
        star_error = pd.to_numeric(errors[star_id], errors="coerce")
        other_errors = errors[other_columns]
        reference_error = np.sqrt((other_errors**2).sum(axis=1, skipna=True)) / ensemble_n
        reference_error = 1.253 * reference_error.where(ensemble_n >= int(min_ensemble))
        residual_error = np.sqrt(star_error**2 + reference_error**2)

        result = pd.DataFrame(
            {
                "frame": centered.index.astype(str),
                "star_id": int(star_id),
                "mag": (
                    work.loc[work["star_id"] == star_id]
                    .set_index("frame")["mag"]
                    .reindex(centered.index)
                    .to_numpy(float)
                ),
                "centered_mag": star_centered.to_numpy(float),
                "ensemble_mag": reference.to_numpy(float),
                "residual": residual.to_numpy(float),
                "residual_err": residual_error.to_numpy(float),
                "ensemble_n": ensemble_n.to_numpy(int),
            }
        )
        for column in frame_meta_columns:
            result[column] = frame_meta[column].reindex(centered.index).to_numpy()
        if not xy.empty and star_id in xy.index.get_level_values("star_id"):
            star_xy = xy.xs(star_id, level="star_id")
            for column in xy.columns:
                result[column] = pd.to_numeric(
                    star_xy[column].reindex(centered.index), errors="coerce"
                ).to_numpy(float)
        rows.append(result)
    if not rows:
        return pd.DataFrame()
    output = pd.concat(rows, ignore_index=True)
    if "time" in output.columns:
        output["time"] = pd.to_numeric(output["time"], errors="coerce")
        output = output.sort_values(["star_id", "time", "frame"], na_position="last")
    else:
        output = output.sort_values(["star_id", "frame"])
    return output.reset_index(drop=True)


def compute_stability_metrics(
    residuals: pd.DataFrame,
    *,
    total_frames: int | None = None,
) -> pd.DataFrame:
    """Summarize scatter, temporal coherence, night offsets, and systematics."""
    if residuals is None or residuals.empty or "star_id" not in residuals.columns:
        return pd.DataFrame()
    n_total = int(total_frames or residuals["frame"].nunique())
    rows: list[dict] = []
    for star_id, group in residuals.groupby("star_id"):
        ordered = group.copy()
        if "time" in ordered.columns and np.any(
            np.isfinite(pd.to_numeric(ordered["time"], errors="coerce"))
        ):
            ordered = ordered.sort_values("time", na_position="last")
        values = pd.to_numeric(ordered["residual"], errors="coerce").to_numpy(float)
        valid = np.isfinite(values)
        values = values[valid]
        n_points = int(values.size)
        if n_points:
            center = float(np.nanmedian(values))
            centered_values = values - center
            rms = float(np.sqrt(np.nanmean(centered_values**2)))
            robust_sigma = _robust_sigma(values)
        else:
            center = rms = robust_sigma = float("nan")
            centered_values = np.array([], dtype=float)

        if "residual_err" in ordered.columns:
            error = pd.to_numeric(
                ordered["residual_err"], errors="coerce"
            ).to_numpy(float)
            error = error[valid]
        else:
            error = np.full(n_points, np.nan)
        chi_mask = np.isfinite(error) & (error > 0) & np.isfinite(centered_values)
        chi2_red = (
            float(np.sum((centered_values[chi_mask] / error[chi_mask]) ** 2) / max(int(chi_mask.sum()) - 1, 1))
            if int(chi_mask.sum()) >= 2
            else float("nan")
        )

        night_scatter = float("nan")
        n_nights = 0
        if "night_id" in ordered.columns:
            night_frame = ordered.loc[valid, ["night_id"]].copy()
            night_frame["residual"] = values
            night_frame["night_id"] = pd.to_numeric(night_frame["night_id"], errors="coerce")
            night_medians = (
                night_frame[np.isfinite(night_frame["night_id"]) & (night_frame["night_id"] > 0)]
                .groupby("night_id")["residual"]
                .median()
                .to_numpy(float)
            )
            n_nights = int(len(night_medians))
            if n_nights >= 2:
                night_scatter = float(np.std(night_medians, ddof=1))

        time_span = float("nan")
        slope_span = float("nan")
        if "time" in ordered.columns:
            times = pd.to_numeric(ordered["time"], errors="coerce").to_numpy(float)
            times = times[valid]
            fit_mask = np.isfinite(times) & np.isfinite(values)
            if int(fit_mask.sum()) >= 3:
                x = times[fit_mask]
                y = values[fit_mask]
                time_span = float(np.nanmax(x) - np.nanmin(x))
                if time_span > 0:
                    slope = float(np.polyfit(x - np.nanmedian(x), y, 1)[0])
                    slope_span = abs(slope) * time_span

        condition_corr, condition_name = _condition_correlation(ordered.loc[valid])
        outlier_fraction = float("nan")
        if n_points >= 3 and np.isfinite(robust_sigma) and robust_sigma > 0:
            outlier_fraction = float(
                np.mean(np.abs(values - center) > 4.0 * robust_sigma)
            )
        rows.append(
            {
                "star_id": int(star_id),
                "n": n_points,
                "coverage": float(n_points / max(n_total, 1)),
                "mean_mag": float(np.nanmedian(pd.to_numeric(ordered["mag"], errors="coerce"))),
                "rms": rms,
                "robust_sigma": robust_sigma,
                "chi2_red": chi2_red,
                "eta": _von_neumann_eta(values),
                "max_run_fraction": _max_run_fraction(values),
                "night_scatter": night_scatter,
                "n_nights": n_nights,
                "slope_span": slope_span,
                "time_span": time_span,
                "outlier_fraction": outlier_fraction,
                "condition_corr": condition_corr,
                "condition_name": condition_name,
            }
        )
    return pd.DataFrame(rows).sort_values("star_id").reset_index(drop=True)


def score_stability_metrics(
    metrics: pd.DataFrame,
    *,
    target_mag: float = float("nan"),
    target_color: float = float("nan"),
    color_by_id: dict[int, float] | None = None,
    config: ComparisonSelectionConfig | None = None,
) -> pd.DataFrame:
    """Assign adaptive variability scores using a local magnitude noise floor."""
    cfg = config or ComparisonSelectionConfig()
    if metrics is None or metrics.empty:
        return pd.DataFrame()
    work = metrics.copy().reset_index(drop=True)
    for column in (
        "mean_mag",
        "rms",
        "robust_sigma",
        "chi2_red",
        "eta",
        "max_run_fraction",
        "night_scatter",
        "slope_span",
        "condition_corr",
    ):
        work[column] = pd.to_numeric(work.get(column), errors="coerce")

    scatter_z: list[float] = []
    for _, row in work.iterrows():
        distances = np.abs(work["mean_mag"] - float(row["mean_mag"]))
        local_count = min(len(work), max(5, int(np.ceil(np.sqrt(len(work)) * 2))))
        local = work.loc[distances.sort_values().index[:local_count], "robust_sigma"].to_numpy(float)
        local = local[np.isfinite(local)]
        value = float(row["robust_sigma"])
        if not len(local) or not np.isfinite(value):
            scatter_z.append(float("inf"))
            continue
        center = float(np.nanmedian(local))
        scale = max(_safe_positive_scale(local), 0.15 * max(abs(center), 1e-4), 1e-4)
        scatter_z.append(max(0.0, (value - center) / scale))
    work["scatter_z"] = scatter_z

    night_values = work["night_scatter"].to_numpy(float)
    finite_night = night_values[np.isfinite(night_values)]
    night_center = float(np.nanmedian(finite_night)) if finite_night.size else 0.0
    night_scale = _safe_positive_scale(finite_night) if finite_night.size else 1e-4
    work["night_z"] = np.where(
        np.isfinite(work["night_scatter"]),
        np.maximum(0.0, (work["night_scatter"] - night_center) / night_scale),
        0.0,
    )
    work["eta_score"] = np.where(
        np.isfinite(work["eta"]), np.maximum(0.0, (1.5 - work["eta"]) / 0.30), 0.0
    )
    expected_run = np.log2(np.maximum(work["n"].astype(float), 2.0)) / np.maximum(
        work["n"].astype(float), 2.0
    )
    work["run_score"] = np.where(
        np.isfinite(work["max_run_fraction"]),
        np.maximum(0.0, (work["max_run_fraction"] - expected_run - 0.05) / 0.05),
        0.0,
    )
    work["slope_score"] = np.where(
        np.isfinite(work["slope_span"]) & np.isfinite(work["robust_sigma"]),
        np.maximum(
            0.0,
            work["slope_span"] / np.maximum(work["robust_sigma"], 1e-4) - 2.0,
        ),
        0.0,
    )
    work["condition_score"] = np.where(
        np.isfinite(work["condition_corr"]),
        np.maximum(0.0, (work["condition_corr"] - 0.55) / 0.15),
        0.0,
    )
    work["variability_score"] = (
        work["scatter_z"]
        + 0.45 * work["eta_score"]
        + 0.35 * work["night_z"]
        + 0.20 * work["run_score"]
        + 0.20 * work["slope_score"]
    )
    work["quality_score"] = work["variability_score"] + 0.25 * work["condition_score"]

    colors = color_by_id or {}
    work["color"] = work["star_id"].map(
        lambda value: float(colors.get(int(value), np.nan))
    )
    work["d_mag"] = (
        np.abs(work["mean_mag"] - float(target_mag)) if np.isfinite(target_mag) else 0.0
    )
    work["d_color"] = (
        np.abs(work["color"] - float(target_color)) if np.isfinite(target_color) else np.nan
    )
    work["selection_score"] = (
        work["quality_score"]
        + float(cfg.magnitude_weight) * work["d_mag"].fillna(3.0)
        + float(cfg.color_weight) * np.minimum(work["d_color"].fillna(1.0), 2.0)
    )

    max_points = max(int(work["n"].max()), 1)
    required_points = min(int(cfg.min_points), max_points)
    hard_fail = (work["coverage"] < float(cfg.min_coverage)) | (
        work["n"] < required_points
    )
    work["status"] = "stable"
    work.loc[work["quality_score"] >= float(cfg.suspect_score), "status"] = "suspect"
    work.loc[work["quality_score"] >= float(cfg.reject_score), "status"] = "reject"
    work.loc[hard_fail, "status"] = "reject"

    reasons: list[str] = []
    for _, row in work.iterrows():
        row_reasons: list[str] = []
        if float(row["coverage"]) < float(cfg.min_coverage):
            row_reasons.append("low_coverage")
        if float(row["scatter_z"]) >= 3.0:
            row_reasons.append("excess_scatter")
        if float(row["eta_score"]) >= 1.0:
            row_reasons.append("time_correlated")
        if float(row["night_z"]) >= 2.0:
            row_reasons.append("night_shift")
        if float(row["condition_corr"]) >= 0.70:
            name = str(row.get("condition_name", "condition") or "condition")
            row_reasons.append(f"correlated_{name}")
        if not row_reasons:
            row_reasons.append("stable")
        reasons.append(",".join(row_reasons))
    work["reasons"] = reasons
    return work.sort_values(["selection_score", "star_id"]).reset_index(drop=True)


def recommend_check_candidate(
    metrics: pd.DataFrame,
    *,
    excluded_ids: Iterable[int] = (),
    min_coverage: float = 0.90,
) -> dict | None:
    """Return the best stable source that remains independent of the ensemble."""
    required = {"star_id", "status", "coverage", "selection_score"}
    if metrics is None or metrics.empty or not required <= set(metrics.columns):
        return None

    work = metrics.copy()
    work["star_id"] = pd.to_numeric(work["star_id"], errors="coerce").astype("Int64")
    work["coverage"] = pd.to_numeric(work["coverage"], errors="coerce")
    work["selection_score"] = pd.to_numeric(
        work["selection_score"], errors="coerce"
    )
    excluded = {int(value) for value in excluded_ids}
    eligible = work[
        work["star_id"].notna()
        & (work["status"].astype(str).str.lower() == "stable")
        & (work["coverage"] >= float(min_coverage))
        & np.isfinite(work["selection_score"])
        & ~work["star_id"].isin(excluded)
    ].copy()
    if eligible.empty:
        return None

    eligible["star_id"] = eligible["star_id"].astype("int64")
    sort_columns = ["selection_score"]
    ascending = [True]
    if "robust_sigma" in eligible.columns:
        eligible["robust_sigma"] = pd.to_numeric(
            eligible["robust_sigma"], errors="coerce"
        )
        sort_columns.append("robust_sigma")
        ascending.append(True)
    sort_columns.extend(["coverage", "star_id"])
    ascending.extend([False, True])
    return eligible.sort_values(
        sort_columns, ascending=ascending, na_position="last"
    ).iloc[0].to_dict()


def select_adaptive_ensemble(
    measurements: pd.DataFrame,
    metrics: pd.DataFrame,
    *,
    min_comparisons: int = 3,
    max_comparisons: int = 12,
    near_best_fraction: float = 0.03,
) -> dict:
    """Choose ensemble size by internal LOO scatter and reserve one check star."""
    required = {"star_id", "status", "selection_score"}
    empty = {
        "selected_ids": [],
        "check_id": None,
        "check_metrics": {},
        "ensemble_trials": pd.DataFrame(),
    }
    if metrics is None or metrics.empty or not required <= set(metrics.columns):
        return {**empty, "reason": "no stability metrics"}

    ranked = metrics.copy()
    ranked["star_id"] = pd.to_numeric(ranked["star_id"], errors="coerce").astype("Int64")
    ranked["selection_score"] = pd.to_numeric(
        ranked["selection_score"], errors="coerce"
    )
    ranked = ranked[
        ranked["star_id"].notna()
        & (ranked["status"].astype(str).str.lower() == "stable")
        & np.isfinite(ranked["selection_score"])
    ].copy()
    ranked["star_id"] = ranked["star_id"].astype("int64")
    ranked = ranked.sort_values(
        ["selection_score", "star_id"], na_position="last"
    ).drop_duplicates("star_id", keep="first")

    minimum = max(3, int(min_comparisons))
    if len(ranked) < minimum + 1:
        return {
            **empty,
            "reason": (
                f"only {len(ranked)} stable candidates; {minimum} comparisons "
                "plus one check star are required"
            ),
        }

    ordered_ids = ranked["star_id"].astype("int64").tolist()
    largest = min(int(max_comparisons), len(ordered_ids) - 1)
    total_frames = int(measurements["frame"].nunique()) if not measurements.empty else 0
    trials: list[dict] = []
    for count in range(minimum, largest + 1):
        comparison_ids = ordered_ids[:count]
        residuals = compute_leave_one_out_residuals(
            measurements, comparison_ids, min_ensemble=2
        )
        trial_metrics = compute_stability_metrics(
            residuals, total_frames=total_frames
        )
        scatter = pd.to_numeric(
            trial_metrics.get("robust_sigma"), errors="coerce"
        ).to_numpy(float)
        scatter = scatter[np.isfinite(scatter)]
        if scatter.size == 0:
            continue
        trials.append(
            {
                "comparison_count": int(count),
                "median_loo_sigma": float(np.nanmedian(scatter)),
                "max_loo_sigma": float(np.nanmax(scatter)),
                "comparison_source_ids": comparison_ids,
            }
        )
    trial_frame = pd.DataFrame(trials)
    if trial_frame.empty:
        return {**empty, "reason": "ensemble scatter could not be evaluated"}

    best_sigma = float(trial_frame["median_loo_sigma"].min())
    threshold = best_sigma * (1.0 + max(float(near_best_fraction), 0.0))
    chosen = trial_frame[
        trial_frame["median_loo_sigma"] <= threshold
    ].sort_values("comparison_count").iloc[0]
    selected_ids = [int(value) for value in chosen["comparison_source_ids"]]
    check_id = next(
        int(value) for value in ordered_ids if int(value) not in set(selected_ids)
    )

    check_residuals = compute_leave_one_out_residuals(
        measurements, [check_id, *selected_ids], min_ensemble=2
    )
    check_only = check_residuals[check_residuals["star_id"] == check_id].copy()
    check_metrics_frame = compute_stability_metrics(
        check_only, total_frames=total_frames
    )
    check_metrics = (
        check_metrics_frame.iloc[0].to_dict()
        if not check_metrics_frame.empty
        else {}
    )
    return {
        "selected_ids": selected_ids,
        "check_id": check_id,
        "check_metrics": check_metrics,
        "ensemble_trials": trial_frame,
        "reason": "adaptive LOO scatter",
    }


def select_stable_comparisons(
    measurements: pd.DataFrame,
    candidate_ids: Iterable[int],
    *,
    target_mag: float = float("nan"),
    target_color: float = float("nan"),
    color_by_id: dict[int, float] | None = None,
    config: ComparisonSelectionConfig | None = None,
) -> dict:
    """Iteratively remove unstable candidates and return a final comparison set."""
    cfg = config or ComparisonSelectionConfig()
    active = sorted({int(value) for value in candidate_ids})
    removed_rows: list[dict] = []
    total_frames = int(measurements["frame"].nunique()) if not measurements.empty else 0

    for iteration in range(max(1, int(cfg.max_iterations))):
        if len(active) <= int(cfg.min_comparisons):
            break
        residuals = compute_leave_one_out_residuals(
            measurements, active, min_ensemble=cfg.min_ensemble
        )
        metrics = compute_stability_metrics(residuals, total_frames=total_frames)
        metrics = score_stability_metrics(
            metrics,
            target_mag=target_mag,
            target_color=target_color,
            color_by_id=color_by_id,
            config=cfg,
        )
        rejected = metrics[metrics["status"] == "reject"].copy()
        if rejected.empty:
            break
        worst = rejected.sort_values(
            ["quality_score", "coverage", "star_id"],
            ascending=[False, True, True],
        ).iloc[0]
        star_id = int(worst["star_id"])
        if len(active) - 1 < int(cfg.min_comparisons):
            break
        removed_rows.append(
            {
                **worst.to_dict(),
                "iteration_removed": int(iteration + 1),
            }
        )
        active.remove(star_id)

    final_residuals = compute_leave_one_out_residuals(
        measurements, active, min_ensemble=cfg.min_ensemble
    )
    final_metrics = compute_stability_metrics(final_residuals, total_frames=total_frames)
    final_metrics = score_stability_metrics(
        final_metrics,
        target_mag=target_mag,
        target_color=target_color,
        color_by_id=color_by_id,
        config=cfg,
    )
    if not final_metrics.empty:
        usable = final_metrics[final_metrics["status"] != "reject"].copy()
        selected = (
            usable.sort_values(["selection_score", "star_id"])
            .head(int(cfg.target_count))["star_id"]
            .astype("int64")
            .tolist()
        )
    else:
        selected = []

    removed = pd.DataFrame(removed_rows)
    if not removed.empty:
        removed["status"] = "reject"
    report = pd.concat([final_metrics, removed], ignore_index=True, sort=False)
    if not report.empty:
        report = report.sort_values(["status", "selection_score", "star_id"], na_position="last")
    return {
        "selected_ids": sorted(set(int(value) for value in selected)),
        "active_ids": active,
        "metrics": report.reset_index(drop=True),
        "residuals": final_residuals.reset_index(drop=True),
        "removed_ids": sorted(
            set(int(value) for value in removed.get("star_id", pd.Series(dtype=int)).dropna())
        ),
    }


def build_target_difference(
    measurements: pd.DataFrame,
    target_id: int,
    comparison_id: int,
) -> pd.DataFrame:
    """Return a centered target-minus-comparison differential series."""
    if measurements is None or measurements.empty:
        return pd.DataFrame()
    subset = measurements[measurements["star_id"].isin([int(target_id), int(comparison_id)])]
    pivot = subset.pivot_table(index="frame", columns="star_id", values="mag", aggfunc="first")
    if int(target_id) not in pivot.columns or int(comparison_id) not in pivot.columns:
        return pd.DataFrame()
    difference = pivot[int(target_id)] - pivot[int(comparison_id)]
    difference = difference - float(np.nanmedian(difference.to_numpy(float)))
    output = pd.DataFrame({"frame": pivot.index.astype(str), "value": difference.to_numpy(float)})
    meta_columns = [column for column in ("time", "night_id") if column in measurements.columns]
    if meta_columns:
        meta = measurements.groupby("frame")[meta_columns].first()
        for column in meta_columns:
            output[column] = meta[column].reindex(pivot.index).to_numpy()
    return output
