"""Matching and summary metrics for artificial-star tests."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import minimize
from scipy.special import expit, logit
from scipy.spatial import cKDTree


def match_truth_to_detections(
    truth: pd.DataFrame,
    detections: pd.DataFrame,
    *,
    radius_px: float,
) -> pd.DataFrame:
    """Greedily create one-to-one truth/detection matches by distance."""
    out = truth.copy().reset_index(drop=True)
    out["recovered"] = False
    out["detection_index"] = -1
    out["x_detected"] = np.nan
    out["y_detected"] = np.nan
    out["position_error_px"] = np.nan
    if len(out) == 0 or len(detections) == 0:
        return out

    truth_xy = out[["x_true", "y_true"]].to_numpy(float)
    det_xy = detections[["x", "y"]].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    finite_det = np.isfinite(det_xy).all(axis=1)
    valid_indices = np.flatnonzero(finite_det)
    det_xy_valid = det_xy[finite_det]
    if len(det_xy_valid) == 0:
        return out

    tree = cKDTree(det_xy_valid)
    neighbor_lists = tree.query_ball_point(truth_xy, r=float(radius_px))
    candidates: list[tuple[float, int, int]] = []
    for truth_idx, neighbors in enumerate(neighbor_lists):
        for local_det_idx in neighbors:
            distance = float(np.hypot(*(truth_xy[truth_idx] - det_xy_valid[local_det_idx])))
            candidates.append((distance, truth_idx, int(local_det_idx)))
    candidates.sort()

    used_truth: set[int] = set()
    used_det: set[int] = set()
    for distance, truth_idx, local_det_idx in candidates:
        if truth_idx in used_truth or local_det_idx in used_det:
            continue
        used_truth.add(truth_idx)
        used_det.add(local_det_idx)
        det_idx = int(valid_indices[local_det_idx])
        out.loc[truth_idx, "recovered"] = True
        out.loc[truth_idx, "detection_index"] = det_idx
        out.loc[truth_idx, "x_detected"] = det_xy[det_idx, 0]
        out.loc[truth_idx, "y_detected"] = det_xy[det_idx, 1]
        out.loc[truth_idx, "position_error_px"] = distance
    return out


def mark_baseline_confounding(
    matched: pd.DataFrame,
    baseline_detections: pd.DataFrame,
    *,
    radius_px: float,
) -> pd.DataFrame:
    """Flag truth positions already occupied by a baseline detection."""
    out = matched.copy()
    out["baseline_confounded"] = False
    out["nearest_baseline_detection_px"] = np.nan
    if len(out) == 0 or len(baseline_detections) == 0:
        return out
    base_xy = baseline_detections[["x", "y"]].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    base_xy = base_xy[np.isfinite(base_xy).all(axis=1)]
    if len(base_xy) == 0:
        return out
    distances, _ = cKDTree(base_xy).query(out[["x_true", "y_true"]].to_numpy(float), k=1)
    out["nearest_baseline_detection_px"] = distances
    out["baseline_confounded"] = distances <= float(radius_px)
    return out


def count_new_false_detections(
    matched: pd.DataFrame,
    injected_detections: pd.DataFrame,
    baseline_detections: pd.DataFrame,
    *,
    baseline_radius_px: float,
) -> tuple[int, int]:
    """Return (new detections, new detections not matched to injected truth)."""
    if len(injected_detections) == 0:
        return 0, 0
    inj_xy = injected_detections[["x", "y"]].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    finite = np.isfinite(inj_xy).all(axis=1)
    new_mask = finite.copy()
    if len(baseline_detections):
        base_xy = baseline_detections[["x", "y"]].apply(pd.to_numeric, errors="coerce").to_numpy(float)
        base_xy = base_xy[np.isfinite(base_xy).all(axis=1)]
        if len(base_xy):
            distances, _ = cKDTree(base_xy).query(inj_xy[finite], k=1)
            new_mask[np.flatnonzero(finite)] = distances > float(baseline_radius_px)

    matched_indices = set(
        pd.to_numeric(
            matched.loc[matched["recovered"], "detection_index"], errors="coerce"
        ).dropna().astype(int)
    )
    new_indices = set(np.flatnonzero(new_mask).tolist())
    false_indices = new_indices - matched_indices
    return len(new_indices), len(false_indices)


def summarize_results(stars: pd.DataFrame) -> dict:
    eligible = ~stars["baseline_confounded"].astype(bool)
    recovered = stars["recovered"].astype(bool)
    good_mag = np.isfinite(pd.to_numeric(stars.get("forced_mag_error"), errors="coerce"))
    pos = pd.to_numeric(stars.get("position_error_px"), errors="coerce")
    mag_error = pd.to_numeric(stars.get("forced_mag_error"), errors="coerce")

    n_eligible = int(eligible.sum())
    n_recovered = int((eligible & recovered).sum())
    pos_values = pos[eligible & recovered & np.isfinite(pos)].to_numpy(float)
    mag_values = mag_error[eligible & good_mag].to_numpy(float)
    return {
        "n_injected": int(len(stars)),
        "n_baseline_confounded": int((~eligible).sum()),
        "n_completeness_eligible": n_eligible,
        "n_recovered": n_recovered,
        "completeness": float(n_recovered / n_eligible) if n_eligible else float("nan"),
        "position_rmse_px": (
            float(np.sqrt(np.mean(pos_values**2))) if len(pos_values) else float("nan")
        ),
        "forced_mag_bias_median": (
            float(np.median(mag_values)) if len(mag_values) else float("nan")
        ),
        "forced_mag_scatter_mad": (
            float(1.4826 * np.median(np.abs(mag_values - np.median(mag_values))))
            if len(mag_values)
            else float("nan")
        ),
    }


def summarize_by_placement(stars: pd.DataFrame) -> dict[str, dict]:
    if "placement_stratum" not in stars.columns:
        return {}
    return {
        str(name): summarize_results(group.reset_index(drop=True))
        for name, group in stars.groupby("placement_stratum", dropna=False)
    }


def magnitude_bin_summary(stars: pd.DataFrame, bin_width: float = 1.0) -> pd.DataFrame:
    if len(stars) == 0:
        return pd.DataFrame()
    mags = pd.to_numeric(stars["magnitude_true"], errors="coerce")
    finite = mags[np.isfinite(mags)]
    if finite.empty:
        return pd.DataFrame()
    lo = np.floor(finite.min() / bin_width) * bin_width
    hi = np.ceil(finite.max() / bin_width) * bin_width + bin_width
    edges = np.arange(lo, hi + 0.5 * bin_width, bin_width)
    labels = 0.5 * (edges[:-1] + edges[1:])
    bin_ids = pd.cut(mags, bins=edges, labels=False, include_lowest=True, right=False)
    rows = []
    for idx, center in enumerate(labels):
        group = stars[bin_ids == idx]
        if group.empty:
            continue
        eligible = ~group["baseline_confounded"].astype(bool)
        recovered = group["recovered"].astype(bool)
        mag_error = pd.to_numeric(group.get("forced_mag_error"), errors="coerce")
        good = eligible & np.isfinite(mag_error)
        rows.append(
            {
                "magnitude_center": float(center),
                "n_injected": int(len(group)),
                "n_eligible": int(eligible.sum()),
                "n_recovered": int((eligible & recovered).sum()),
                "completeness": (
                    float((eligible & recovered).sum() / eligible.sum())
                    if eligible.sum()
                    else np.nan
                ),
                "forced_mag_bias_median": (
                    float(np.nanmedian(mag_error[good])) if good.any() else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def magnitude_point_summary(stars: pd.DataFrame) -> pd.DataFrame:
    """Summarize exact injected magnitude points with Wilson 95% intervals."""
    if len(stars) == 0:
        return pd.DataFrame()
    z = 1.959963984540054
    rows = []
    for magnitude, group in stars.groupby("magnitude_true", sort=True):
        eligible = ~group["baseline_confounded"].astype(bool)
        recovered = group["recovered"].astype(bool)
        n = int(eligible.sum())
        k = int((eligible & recovered).sum())
        p = float(k / n) if n else np.nan
        if n:
            denominator = 1.0 + z * z / n
            center = (p + z * z / (2.0 * n)) / denominator
            half = z * np.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / denominator
            ci_low = max(0.0, float(center - half))
            ci_high = min(1.0, float(center + half))
        else:
            ci_low = np.nan
            ci_high = np.nan
        mag_error = pd.to_numeric(group.get("forced_mag_error"), errors="coerce")
        good_mag = eligible & np.isfinite(mag_error)
        rows.append(
            {
                "magnitude": float(magnitude),
                "n_injected": int(len(group)),
                "n_eligible": n,
                "n_recovered": k,
                "completeness": p,
                "ci95_low": ci_low,
                "ci95_high": ci_high,
                "forced_mag_bias_median": (
                    float(np.nanmedian(mag_error[good_mag])) if good_mag.any() else np.nan
                ),
                "forced_mag_scatter_mad": (
                    float(
                        1.4826
                        * np.nanmedian(
                            np.abs(mag_error[good_mag] - np.nanmedian(mag_error[good_mag]))
                        )
                    )
                    if good_mag.any()
                    else np.nan
                ),
            }
        )
    return pd.DataFrame(rows)


def _fit_logistic_once(magnitude: np.ndarray, recovered: np.ndarray) -> tuple[float, float]:
    magnitude = np.asarray(magnitude, dtype=float)
    recovered = np.asarray(recovered, dtype=float)
    finite = np.isfinite(magnitude) & np.isfinite(recovered)
    magnitude = magnitude[finite]
    recovered = recovered[finite]
    if len(magnitude) < 10 or recovered.min() == recovered.max():
        raise RuntimeError("Completeness fit requires both recovered and missed stars")

    def objective(theta):
        m50 = float(theta[0])
        width = float(np.exp(theta[1]))
        probability = np.clip(expit((m50 - magnitude) / width), 1e-12, 1.0 - 1e-12)
        return float(
            -np.sum(
                recovered * np.log(probability)
                + (1.0 - recovered) * np.log1p(-probability)
            )
        )

    initial_m50 = float(np.nanmedian(magnitude))
    result = minimize(
        objective,
        x0=np.array([initial_m50, np.log(0.25)], dtype=float),
        method="L-BFGS-B",
        bounds=[
            (float(np.nanmin(magnitude) - 2.0), float(np.nanmax(magnitude) + 2.0)),
            (np.log(0.01), np.log(5.0)),
        ],
    )
    if not result.success:
        raise RuntimeError(f"Completeness logistic fit failed: {result.message}")
    return float(result.x[0]), float(np.exp(result.x[1]))


def fit_completeness_logistic(
    stars: pd.DataFrame,
    *,
    bootstrap_samples: int = 1000,
    seed: int = 20260607,
) -> dict:
    """Fit a binomial logistic completeness curve and cluster-bootstrap trials."""
    eligible = ~stars["baseline_confounded"].astype(bool)
    data = stars.loc[eligible, ["magnitude_true", "recovered", "trial"]].copy()
    data["magnitude_true"] = pd.to_numeric(data["magnitude_true"], errors="coerce")
    data = data[np.isfinite(data["magnitude_true"])].reset_index(drop=True)
    magnitude = data["magnitude_true"].to_numpy(float)
    recovered = data["recovered"].astype(bool).to_numpy(float)
    m50, width = _fit_logistic_once(magnitude, recovered)

    trial_values = data["trial"].dropna().unique()
    rng = np.random.default_rng(seed)
    estimates = []
    requested = max(0, int(bootstrap_samples))
    for _ in range(requested):
        if len(trial_values) >= 2:
            sampled_trials = rng.choice(trial_values, size=len(trial_values), replace=True)
            pieces = []
            for bootstrap_trial, original_trial in enumerate(sampled_trials):
                piece = data[data["trial"] == original_trial].copy()
                piece["trial"] = bootstrap_trial
                pieces.append(piece)
            sample = pd.concat(pieces, ignore_index=True)
        else:
            indices = rng.integers(0, len(data), size=len(data))
            sample = data.iloc[indices]
        try:
            estimates.append(
                _fit_logistic_once(
                    sample["magnitude_true"].to_numpy(float),
                    sample["recovered"].astype(bool).to_numpy(float),
                )
            )
        except RuntimeError:
            continue

    m90 = float(m50 - width * logit(0.9))
    m10 = float(m50 - width * logit(0.1))
    result = {
        "model": "p(recovered)=expit((m50-magnitude)/width)",
        "n_stars": int(len(data)),
        "n_trials": int(len(trial_values)),
        "m50": m50,
        "width_mag": width,
        "m90": m90,
        "m10": m10,
        "bootstrap_requested": requested,
        "bootstrap_successful": int(len(estimates)),
    }
    if estimates:
        estimate_array = np.asarray(estimates, dtype=float)
        m50_boot = estimate_array[:, 0]
        width_boot = estimate_array[:, 1]
        m90_boot = m50_boot - width_boot * logit(0.9)
        m10_boot = m50_boot - width_boot * logit(0.1)
        for name, values in (
            ("m50", m50_boot),
            ("width_mag", width_boot),
            ("m90", m90_boot),
            ("m10", m10_boot),
        ):
            result[f"{name}_ci95_low"] = float(np.percentile(values, 2.5))
            result[f"{name}_ci95_high"] = float(np.percentile(values, 97.5))
    return result
