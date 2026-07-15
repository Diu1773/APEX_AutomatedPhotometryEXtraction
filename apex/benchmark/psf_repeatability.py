"""CPU-only exact cross-validated repeatability for Step 8 PSF photometry."""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.optimize import nnls
from itertools import combinations


def _num(df: pd.DataFrame, name: str, default=np.nan) -> np.ndarray:
    if name not in df:
        return np.full(len(df), default, dtype=float)
    return pd.to_numeric(df[name], errors="coerce").to_numpy(float)


def _xy(df: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    for x, y in (("x_fit", "y_fit"), ("x", "y")):
        if x in df and y in df:
            return _num(df, x), _num(df, y)
    raise ValueError("photometry table needs x_fit/y_fit or x/y")


def _greedy_match(psf: pd.DataFrame, forced: pd.DataFrame, radius: float) -> tuple[np.ndarray, np.ndarray]:
    """Generate local candidates with cKDTree, then resolve globally by distance."""
    px, py = _xy(psf); fx, fy = _xy(forced)
    ids = _num(forced, "master_id")
    out = np.full(len(psf), np.nan); distances = np.full(len(psf), np.nan)
    pgood = np.isfinite(px) & np.isfinite(py); fgood = np.isfinite(fx) & np.isfinite(fy)
    if not pgood.any() or not fgood.any():
        return out, distances
    tree = cKDTree(np.column_stack((fx[fgood], fy[fgood])))
    forced_rows = np.flatnonzero(fgood)
    candidates = []
    for i in np.flatnonzero(pgood):
        point = np.array([px[i], py[i]])
        for local_j in tree.query_ball_point(point, float(radius)):
            j = int(forced_rows[local_j])
            candidates.append((float(np.hypot(px[i] - fx[j], py[i] - fy[j])), int(i), j))
    used_p: set[int] = set(); used_f: set[int] = set()
    for distance, i, j in sorted(candidates):
        if i in used_p or j in used_f:
            continue
        used_p.add(i); used_f.add(j); out[i] = ids[j]; distances[i] = distance
    return out, distances


def match_psf_to_step7(psf: pd.DataFrame, forced: pd.DataFrame, radius_px: float = 1.5) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Restore stable ``master_id`` by same-frame one-to-one positional matching."""
    if "master_id" not in forced:
        raise ValueError("Step7 table has no master_id column")
    result = psf.copy(); ids, distances = _greedy_match(result, forced, radius_px)
    result["master_id"] = ids; result["step7_match_distance_px"] = distances
    matched = np.isfinite(ids); finite = distances[matched]
    return result, {
        "method": "same_frame_step7_position_cKDTree_global_greedy",
        "radius_px": float(radius_px), "psf_rows": int(len(result)),
        "matched_rows": int(matched.sum()), "unmatched_rows": int((~matched).sum()),
        "match_fraction": float(matched.mean()) if len(result) else 0.0,
        "match_distance_median_px": float(np.median(finite)) if len(finite) else None,
        "match_distance_p95_px": float(np.percentile(finite, 95)) if len(finite) else None,
        "id_policy": "det_uid/seed_uid ignored; master_id restored from Step7",
    }


def _huber_weights(residual: np.ndarray, scale: float) -> np.ndarray:
    u = np.abs(residual) / max(float(scale), 1e-6)
    return np.minimum(1.0, 1.345 / np.maximum(u, 1e-12))


def fit_additive_model(data: pd.DataFrame, max_iter: int = 80) -> dict[str, Any]:
    """Fit ``mag_psf = star_mag + frame_zp`` and return named residual products.

    ``raw_in_sample_residual`` is deliberately separate from
    ``exact_leave_one_out_cv_residual``.  The latter compares each corrected
    observation with the weighted mean of the *other* observations of that
    star, and is the only residual used for repeatability statistics.
    """
    y = _num(data, "mag_psf"); err = _num(data, "mag_psf_err")
    stars = data["master_id"].to_numpy(); frames = data["frame"].to_numpy()
    unique_s = list(pd.unique(stars)); unique_f = list(pd.unique(frames))
    si = {v: i for i, v in enumerate(unique_s)}; fi = {v: i for i, v in enumerate(unique_f)}
    sidx = np.array([si[v] for v in stars]); fidx = np.array([fi[v] for v in frames])
    base_w = 1.0 / np.maximum(err, 1e-4) ** 2
    star = np.zeros(len(unique_s)); zp = np.zeros(len(unique_f))
    counts = np.bincount(sidx, weights=base_w, minlength=len(star))
    star[:] = np.bincount(sidx, weights=base_w * y, minlength=len(star)) / np.maximum(counts, 1e-12)
    for iteration in range(max_iter):
        old = np.r_[star, zp]
        raw = y - star[sidx] - zp[fidx]
        scale = 1.4826 * np.median(np.abs(raw - np.median(raw)))
        weights = base_w * _huber_weights(raw, max(scale, np.nanmedian(err)))
        for k in range(len(unique_s)):
            mask = sidx == k
            star[k] = np.sum(weights[mask] * (y[mask] - zp[fidx[mask]])) / max(np.sum(weights[mask]), 1e-12)
        for k in range(1, len(unique_f)):
            mask = fidx == k
            zp[k] = np.sum(weights[mask] * (y[mask] - star[sidx[mask]])) / max(np.sum(weights[mask]), 1e-12)
        shift = zp[0]; zp -= shift; star += shift
        if np.max(np.abs(np.r_[star, zp] - old)) < 1e-7:
            break
    raw = y - star[sidx] - zp[fidx]
    corrected = y - zp[fidx]
    cv = np.full(len(data), np.nan); cv_pull = np.full(len(data), np.nan)
    for k in range(len(unique_s)):
        mask = sidx == k; others = np.flatnonzero(mask)
        for i in others:
            rest = others[others != i]
            w = base_w[rest]
            if len(rest):
                other_mean = float(np.sum(w * corrected[rest]) / np.sum(w))
                denominator = np.sqrt(err[i] ** 2 + 1.0 / np.sum(w))
                cv[i] = corrected[i] - other_mean
                cv_pull[i] = cv[i] / denominator
    return {
        "star_mag": dict(zip(unique_s, star)), "frame_zp": dict(zip(unique_f, zp)),
        "corrected_mag": corrected, "raw_in_sample_residual": raw,
        "exact_leave_one_out_cv_residual": cv, "exact_leave_one_out_cv_pull": cv_pull,
        "iterations": iteration + 1,
        "residual_note": "exact leave-one-out cross-validated residual; current observation excluded from same-star weighted mean",
    }


def _mad(x: np.ndarray) -> float:
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if len(x) else np.nan


def _metrics(data: pd.DataFrame, scope: str, label: str, lower: float | None, upper: float | None, values: np.ndarray) -> dict[str, Any]:
    r = values[np.isfinite(values)]; err = _num(data, "mag_psf_err")
    pull = values / np.maximum(np.sqrt(err * err + _num(data, "other_mean_variance", 0.0)), 1e-12)
    pull = pull[np.isfinite(pull)]
    return {"scope": scope, "label": label, "lower": lower, "upper": upper, "n": int(len(r)),
            "median_residual": float(np.median(r)) if len(r) else np.nan, "cv_robust_scatter_mad": _mad(r),
            "cv_rmse": float(np.sqrt(np.mean(r * r))) if len(r) else np.nan,
            "cv_abs_residual_gt_0_2_frac": float(np.mean(np.abs(r) > .2)) if len(r) else np.nan,
            "cv_normalized_pull_scatter": _mad(pull)}


def _bool_column(df: pd.DataFrame, name: str) -> np.ndarray:
    if name not in df: return np.zeros(len(df), dtype=bool)
    value = df[name]
    if pd.api.types.is_bool_dtype(value): return value.fillna(False).to_numpy(bool)
    return value.astype(str).str.strip().str.lower().isin(("1", "true", "yes", "y")).to_numpy()


def _pairwise_products(data: pd.DataFrame, snr_min: float = 0.0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build corrected frame pairs and solve robust/RMSE variance components."""
    eligible = data.loc[data["star_min_snr"] >= float(snr_min), "master_id"].unique()
    data = data[data.master_id.isin(eligible)]
    pair_rows = []
    for star_id, group in data.groupby("master_id"):
        by_frame = {row.frame: row for row in group.itertuples()}
        for left, right in combinations(sorted(by_frame), 2):
            a, b = by_frame[left], by_frame[right]
            delta = float(a.corrected_mag - b.corrected_mag)
            formal = float(np.sqrt(a.mag_psf_err ** 2 + b.mag_psf_err ** 2))
            pair_rows.append({"master_id": star_id, "frame_i": left, "frame_j": right,
                              "pair_label": f"{left}__{right}", "corrected_mag_i": a.corrected_mag,
                              "corrected_mag_j": b.corrected_mag, "pair_delta": delta,
                              "formal_error": formal, "formal_pull": delta / formal if formal > 0 else np.nan})
    observations = pd.DataFrame(pair_rows)
    if observations.empty:
        return observations, pd.DataFrame()
    pair_summary = []
    for label, group in observations.groupby("pair_label", sort=True):
        values = group.pair_delta.to_numpy(float); pulls = group.formal_pull.to_numpy(float)
        pair_summary.append({"frame_i": group.frame_i.iloc[0], "frame_j": group.frame_j.iloc[0], "pair_label": label,
                             "n": len(values), "n_stars": int(group.master_id.nunique()), "snr_min": float(snr_min),
                             "median_delta": np.median(values), "robust_mad": _mad(values),
                             "rmse": np.sqrt(np.mean(values ** 2)), "abs_delta_gt_0_2_frac": np.mean(np.abs(values) > .2),
                             "formal_pull_scatter": _mad(pulls)})
    pairs = pd.DataFrame(pair_summary)
    frames = sorted(set(pairs.frame_i) | set(pairs.frame_j)); frame_index = {name: i for i, name in enumerate(frames)}
    noise_rows = []
    for metric, variance in (("robust", pairs.robust_mad.to_numpy(float) ** 2), ("rmse", pairs.rmse.to_numpy(float) ** 2)):
        valid = np.isfinite(variance); matrix = np.zeros((int(valid.sum()), len(frames)))
        for row, (_, pair) in zip(matrix, pairs.loc[valid].iterrows()):
            row[frame_index[pair.frame_i]] = 1.; row[frame_index[pair.frame_j]] = 1.
        reason = "ok"; estimate = np.full(len(frames), np.nan)
        if len(frames) < 3:
            reason = "underconstrained: at least 3 frames required"
        elif len(matrix) < len(frames) or np.linalg.matrix_rank(matrix) < len(frames):
            reason = "underconstrained: pair design matrix is rank deficient"
        elif len(matrix):
            estimate = np.sqrt(nnls(matrix, variance[valid])[0])
        for name, value in zip(frames, estimate):
            info = data.loc[data.frame == name].iloc[0]
            noise_rows.append({"frame": name, "frame_label": str(name).removeprefix("photometry_").removesuffix(".tsv"),
                               "fwhm_px": info.get("frame_fwhm_px", np.nan), "fwhm_arcsec": info.get("frame_fwhm_arcsec", np.nan),
                               "noise_metric": metric, "snr_min": float(snr_min), "n_stars": int(data.master_id.nunique()),
                               "frame_scatter": value, "status": reason})
    return pairs, pd.DataFrame(noise_rows)


def _geometry(root: Path, frame_name: str) -> dict[str, float]:
    token = frame_name.removeprefix("photometry_").removesuffix(".tsv")
    residual_path = root / "cmd_psf" / f"residual_meta_{token}.json"
    detect_path = root / "cache" / f"detect_{token}.json"
    try: residual = json.loads(residual_path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError): residual = {}
    try: detect = json.loads(detect_path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError): detect = {}
    core = residual.get("core_cut", {})
    fwhm_px = float(detect.get("fwhm_px", np.nan)); fwhm_arcsec = float(detect.get("fwhm_arcsec", np.nan))
    return {"center_x": float(core.get("center_x", np.nan)), "center_y": float(core.get("center_y", np.nan)),
            "fwhm_px": fwhm_px, "fwhm_arcsec": fwhm_arcsec,
            "pixel_scale_arcsec": fwhm_arcsec / fwhm_px if fwhm_px > 0 else np.nan}


def analyze_psf_repeatability(result_dir: str | Path, output_dir: str | Path | None = None, match_radius_px: float = 1.5,
                              min_frames: int = 3, snr_min: float | None = None, qfit_max: float | None = None,
                              crowding_min_px: float | None = None) -> dict[str, Any]:
    root = Path(result_dir); psf_dir = root / "cmd_psf"; step7_dir = root / "step7_forced_phot"
    rows = []; matching = []
    for path in sorted(psf_dir.glob("photometry_*.tsv")):
        forced_path = step7_dir / path.name
        if not forced_path.exists(): continue
        psf, meta = match_psf_to_step7(pd.read_csv(path, sep="\t"), pd.read_csv(forced_path, sep="\t"), match_radius_px)
        psf["frame"] = path.name; geo = _geometry(root, path.name)
        psf["radius_px"] = np.hypot(_num(psf, "x_fit") - geo["center_x"], _num(psf, "y_fit") - geo["center_y"])
        psf["radius_fwhm"] = psf["radius_px"] / geo["fwhm_px"]
        psf["radius_arcmin"] = psf["radius_px"] * geo["pixel_scale_arcsec"] / 60.0
        psf["frame_fwhm_px"] = geo["fwhm_px"]
        psf["frame_fwhm_arcsec"] = geo["fwhm_arcsec"]
        matching.append({"frame": path.name, **meta, **geo}); good = np.isfinite(_num(psf, "master_id")) & (_num(psf, "flags_psf", 1) == 0)
        good &= np.isfinite(_num(psf, "mag_psf")) & np.isfinite(_num(psf, "mag_psf_err")) & (_num(psf, "mag_psf_err") > 0)
        if snr_min is not None: good &= _num(psf, "snr_psf") >= snr_min
        if qfit_max is not None: good &= _num(psf, "qfit") <= qfit_max
        if crowding_min_px is not None: good &= _num(psf, "neighbor_dist_px") >= crowding_min_px
        good &= ~_bool_column(psf, "crowding_unreliable_psf"); rows.append(psf.loc[good].copy())
    if not rows: raise FileNotFoundError("no matching cmd_psf and step7_forced_phot TSV pairs")
    data = pd.concat(rows, ignore_index=True)
    counts = data.groupby("master_id")["frame"].nunique(); keep = counts[counts >= min_frames].index
    data = data[data.master_id.isin(keep)].reset_index(drop=True)
    if data.empty: raise ValueError(f"no stars have at least min_frames={min_frames} observations")
    if "snr_psf" in data:
        star_min_snr = data.groupby("master_id")["snr_psf"].min()
        data["star_min_snr"] = data["master_id"].map(star_min_snr)
    else:
        data["star_min_snr"] = 0.0
    fit = fit_additive_model(data); data["corrected_mag"] = fit["corrected_mag"]
    data["raw_in_sample_residual"] = fit["raw_in_sample_residual"]; data["cv_residual"] = fit["exact_leave_one_out_cv_residual"]
    data["cv_pull"] = fit["exact_leave_one_out_cv_pull"]
    corrected = data["corrected_mag"].to_numpy(); err = _num(data, "mag_psf_err"); other_var = np.zeros(len(data))
    for i in range(len(data)):
        rest = (data.master_id.to_numpy() == data.master_id.iloc[i]) & (np.arange(len(data)) != i)
        other_var[i] = 1.0 / np.sum(1.0 / np.maximum(err[rest], 1e-4) ** 2)
    data["other_mean_variance"] = other_var
    cv = data["cv_residual"].to_numpy(); summary = [_metrics(data, "overall", "all", None, None, cv)]
    specs = [("snr_psf", "snr", [0, 20, 50, 100, np.inf]),
             ("radius_arcmin", "radius", [0, 1, 2, 4, 8, np.inf]),
             ("neighbor_dist_fwhm", "neighbor", [0, 1, 2, 4, 8, np.inf])]
    for field, scope, edges in specs:
        values = _num(data, field)
        for lo, hi in zip(edges[:-1], edges[1:]):
            mask = (values >= lo) & (values < hi)
            if mask.any(): summary.append(_metrics(data.loc[mask], scope, f"{lo:g}-{hi:g}", float(lo), float(hi) if np.isfinite(hi) else None, cv[mask]))
    pair_tables = []; noise_tables = []
    for cut in (0.0, 20.0, 50.0, 100.0):
        pair_table, noise_table = _pairwise_products(data, cut)
        pair_tables.append(pair_table); noise_tables.append(noise_table)
    pair_summary = pd.concat(pair_tables, ignore_index=True) if pair_tables else pd.DataFrame()
    frame_noise = pd.concat(noise_tables, ignore_index=True) if noise_tables else pd.DataFrame()
    star_rows = []
    for star_id, group in data.groupby("master_id"):
        value = group.cv_residual.to_numpy(); star_rows.append({"master_id": star_id, "n_frames": len(group), "median_cv_residual": np.nanmedian(value), "cv_robust_scatter_mad": _mad(value), "cv_rmse": np.sqrt(np.nanmean(value * value))})
    result = {"summary": summary, "frame_offsets": {str(k): float(v) for k, v in fit["frame_zp"].items()}, "matching": matching,
              "fit": {"iterations": fit["iterations"], "residual_note": fit["residual_note"]}, "min_frames": min_frames,
              "n_input_good": int(len(data)), "n_stars": int(data.master_id.nunique()), "n_frames": int(data.frame.nunique()),
              "design": "PSF mag_psf only; det_uid/seed_uid ignored; neighbor_dist_fwhm preferred",
              "interpretation": {"overall_snr_radius_neighbor": "exact leave-one-out predictive disagreement",
                                 "frame_noise": "pairwise variance decomposition; individual frame performance estimate"}}
    if output_dir is not None:
        out = Path(output_dir); out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(summary).to_csv(out / "psf_repeatability_summary.csv", index=False)
        pair_summary.to_csv(out / "psf_repeatability_pairs.csv", index=False)
        frame_noise.to_csv(out / "psf_repeatability_frame_noise.csv", index=False)
        data.to_csv(out / "psf_repeatability_observations.csv", index=False)
        pd.DataFrame(star_rows).to_csv(out / "psf_repeatability_by_star.csv", index=False)
        (out / "psf_repeatability_summary.json").write_text(json.dumps(result, indent=2, allow_nan=True), encoding="utf-8")
    return result


__all__ = ["match_psf_to_step7", "fit_additive_model", "analyze_psf_repeatability"]
