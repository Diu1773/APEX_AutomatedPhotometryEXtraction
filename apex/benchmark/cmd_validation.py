"""Validation reports for CMD artificial-star benchmark campaigns."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass
class CmdValidationConfig:
    batch_root: str = ""
    project_root: str = ""
    output_root: str = "benchmark/runs/cmd_validation"
    repeatability_min_frames: int = 3
    repeatability_min_snr: float = 5.0
    magnitude_bin_width: float = 0.5


def _normal_filter(value: Any) -> str:
    return str(value).strip().lower().rstrip("'")


def _to_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False).astype(bool)
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def _num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _finite(values: Any) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    return arr[np.isfinite(arr)]


def _mad_sigma(values: Any) -> float:
    arr = _finite(values)
    if arr.size == 0:
        return float("nan")
    med = float(np.median(arr))
    return float(1.4826 * np.median(np.abs(arr - med)))


def _rms(values: Any) -> float:
    arr = _finite(values)
    if arr.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(arr * arr)))


def _robust_inlier_mask(values: Any, *, sigma: float = 3.0, iters: int = 5) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    keep = np.isfinite(arr)
    for _ in range(int(iters)):
        subset = arr[keep]
        if subset.size < 3:
            break
        med = float(np.median(subset))
        scale = _mad_sigma(subset)
        if not np.isfinite(scale) or scale <= 0:
            break
        new_keep = np.isfinite(arr) & (np.abs(arr - med) <= float(sigma) * scale)
        if np.array_equal(new_keep, keep):
            break
        keep = new_keep
    return keep


def _read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _condition_rank(value: Any) -> int:
    return {"best": 0, "median": 1, "worst": 2}.get(str(value), 99)


def _load_precision_summary(batch_root: Path) -> pd.DataFrame:
    precision_path = batch_root / "cmd_precision_summary.csv"
    if precision_path.is_file():
        precision = pd.read_csv(precision_path)
    else:
        summary_path = batch_root / "cmd_batch_summary.csv"
        if not summary_path.is_file():
            raise FileNotFoundError(f"CMD batch summary is missing: {summary_path}")
        all_rows = pd.read_csv(summary_path)
        precision = all_rows[all_rows["stage"].astype(str) == "precision"].copy()
    if precision.empty:
        raise RuntimeError("CMD batch does not contain precision rows")
    if "run_dir" not in precision.columns:
        raise RuntimeError("CMD precision summary must include run_dir")
    precision["filter"] = precision["filter"].map(_normal_filter)
    precision["condition_rank"] = precision["condition"].map(_condition_rank)
    return precision.sort_values(["filter", "condition_rank"]).reset_index(drop=True)


def collect_ast_tables(batch_root: str | Path) -> dict[str, pd.DataFrame]:
    """Collect completeness, bias, crowding, and false-positive AST tables."""
    root = Path(batch_root).expanduser().resolve()
    precision = _load_precision_summary(root)
    point_rows: list[pd.DataFrame] = []
    false_rows: list[dict[str, Any]] = []
    crowding_rows: list[dict[str, Any]] = []

    for _, row in precision.iterrows():
        run_dir = Path(str(row["run_dir"]))
        if not run_dir.is_absolute():
            run_dir = (root / run_dir).resolve()
        if not run_dir.is_dir():
            raise FileNotFoundError(f"Precision run directory is missing: {run_dir}")

        points = pd.read_csv(run_dir / "magnitude_points.csv")
        points.insert(0, "condition", row["condition"])
        points.insert(0, "filter", row["filter"])
        points["file"] = row.get("file", "")
        points["fwhm_px"] = row.get("fwhm_px_benchmark", row.get("fwhm_px_step7", np.nan))
        points["run_dir"] = str(run_dir)
        point_rows.append(points)

        trials = pd.read_csv(run_dir / "trials.csv")
        summary = _read_json(run_dir / "summary.json")
        false_rows.append(
            {
                "filter": row["filter"],
                "condition": row["condition"],
                "file": row.get("file", ""),
                "fwhm_px": row.get("fwhm_px_benchmark", summary.get("fwhm_px")),
                "n_trials": int(len(trials)),
                "n_injected": int(trials["n_injected"].sum()),
                "new_detections": int(trials["new_detections"].sum()),
                "new_false_detections": int(trials["new_false_detections"].sum()),
                "false_per_trial": float(trials["new_false_detections"].mean()),
                "false_per_1000_injected": float(
                    1000.0
                    * trials["new_false_detections"].sum()
                    / max(1, trials["n_injected"].sum())
                ),
                "false_fraction_of_new": float(
                    trials["new_false_detections"].sum()
                    / max(1, trials["new_detections"].sum())
                ),
                "run_dir": str(run_dir),
            }
        )

        fit_by_placement = summary.get("completeness_fit_by_placement") or {}
        by_placement = summary.get("by_placement") or {}
        for placement, metrics in by_placement.items():
            fit = fit_by_placement.get(str(placement), {})
            crowding_rows.append(
                {
                    "filter": row["filter"],
                    "condition": row["condition"],
                    "file": row.get("file", ""),
                    "fwhm_px": row.get("fwhm_px_benchmark", summary.get("fwhm_px")),
                    "placement": placement,
                    "n_injected": metrics.get("n_injected"),
                    "n_recovered": metrics.get("n_recovered"),
                    "completeness": metrics.get("completeness"),
                    "m50": fit.get("m50"),
                    "m50_ci95_low": fit.get("m50_ci95_low"),
                    "m50_ci95_high": fit.get("m50_ci95_high"),
                    "forced_mag_bias_median": metrics.get("forced_mag_bias_median"),
                    "forced_mag_scatter_mad": metrics.get("forced_mag_scatter_mad"),
                    "position_rmse_px": metrics.get("position_rmse_px"),
                    "fit_error": fit.get("error"),
                    "run_dir": str(run_dir),
                }
            )

    return {
        "precision": precision,
        "magnitude_points": pd.concat(point_rows, ignore_index=True),
        "false_positive_summary": pd.DataFrame(false_rows),
        "crowding_summary": pd.DataFrame(crowding_rows),
    }


def _write_ast_plots(tables: dict[str, pd.DataFrame], output_dir: Path) -> None:
    points = tables["magnitude_points"].copy()
    precision = tables["precision"].copy()
    false_positive = tables["false_positive_summary"].copy()
    crowding = tables["crowding_summary"].copy()

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), sharey=True)
    for ax, filter_name in zip(axes, sorted(points["filter"].unique())):
        sub = points[points["filter"] == filter_name]
        for condition, group in sub.groupby("condition", sort=False):
            group = group.sort_values("magnitude")
            ax.plot(
                group["magnitude"],
                group["completeness"],
                marker="o",
                linewidth=1.5,
                label=str(condition),
            )
        ax.set_title(str(filter_name))
        ax.set_xlabel("Injected magnitude")
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    axes[0].set_ylabel("Completeness")
    fig.tight_layout()
    fig.savefig(output_dir / "ast_completeness_curves.png", dpi=180)
    plt.close(fig)

    fig, axes = plt.subplots(2, 1, figsize=(8.5, 7.0), sharex=False)
    markers = {"best": "o", "median": "s", "worst": "^"}
    for (filter_name, condition), group in points.groupby(["filter", "condition"], sort=True):
        group = group.sort_values("magnitude")
        label = f"{filter_name} {condition}"
        axes[0].plot(
            group["magnitude"],
            group["forced_mag_bias_median"],
            marker=markers.get(str(condition), "o"),
            linewidth=1.2,
            label=label,
        )
        axes[1].plot(
            group["magnitude"],
            group["forced_mag_scatter_mad"],
            marker=markers.get(str(condition), "o"),
            linewidth=1.2,
            label=label,
        )
    axes[0].axhline(0.0, color="black", linewidth=1)
    axes[0].set_ylabel("Median mag bias")
    axes[1].set_ylabel("MAD scatter (mag)")
    axes[1].set_xlabel("Injected magnitude")
    for ax in axes:
        ax.grid(alpha=0.25)
    axes[1].legend(ncol=3, fontsize=8)
    fig.tight_layout()
    fig.savefig(output_dir / "ast_photometric_bias_scatter.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for filter_name, group in precision.groupby("filter", sort=True):
        group = group.sort_values("fwhm_px_benchmark")
        y = group["m50"].astype(float)
        yerr = np.vstack(
            [
                np.maximum(0.0, y - group["m50_ci95_low"].astype(float)),
                np.maximum(0.0, group["m50_ci95_high"].astype(float) - y),
            ]
        )
        ax.errorbar(
            group["fwhm_px_benchmark"],
            y,
            yerr=yerr,
            marker="o",
            linewidth=1.7,
            capsize=3,
            label=str(filter_name),
        )
    ax.set(xlabel="Benchmark FWHM (px)", ylabel="50% completeness magnitude")
    ax.grid(alpha=0.25)
    ax.legend(title="Filter")
    fig.tight_layout()
    fig.savefig(output_dir / "ast_m50_vs_fwhm.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.8, 4.8))
    for filter_name, group in false_positive.groupby("filter", sort=True):
        group = group.sort_values("fwhm_px")
        ax.plot(
            group["fwhm_px"],
            group["false_per_1000_injected"],
            marker="o",
            linewidth=1.5,
            label=str(filter_name),
        )
    ax.set(xlabel="Benchmark FWHM (px)", ylabel="False detections per 1000 injected")
    ax.grid(alpha=0.25)
    ax.legend(title="Filter")
    fig.tight_layout()
    fig.savefig(output_dir / "ast_false_positive_rate.png", dpi=180)
    plt.close(fig)

    finite_crowding = crowding[np.isfinite(_num(crowding["m50"]))].copy()
    if not finite_crowding.empty:
        fig, ax = plt.subplots(figsize=(8.0, 5.0))
        for placement, group in finite_crowding.groupby("placement", sort=True):
            ax.scatter(
                group["fwhm_px"],
                group["m50"],
                s=55,
                label=str(placement),
                alpha=0.85,
            )
        ax.set(xlabel="Benchmark FWHM (px)", ylabel="Placement-specific m50")
        ax.grid(alpha=0.25)
        ax.legend(title="Placement")
        fig.tight_layout()
        fig.savefig(output_dir / "ast_crowding_m50.png", dpi=180)
        plt.close(fig)


def _resolve_step7_path(row: pd.Series, project_root: Path) -> Path:
    raw = row.get("path")
    if isinstance(raw, str) and raw.strip():
        path = Path(raw)
        if path.is_file():
            return path
    filename = str(row["file"])
    candidate = project_root / "result" / "step7_forced_phot" / f"photometry_{filename}.tsv"
    if candidate.is_file():
        return candidate
    matches = list((project_root / "result" / "step7_forced_phot").glob(f"*{filename}*.tsv"))
    if len(matches) == 1:
        return matches[0]
    raise FileNotFoundError(f"Step 7 photometry table is missing for {filename}")


def load_repeatability_observations(
    project_root: str | Path,
    *,
    min_snr: float = 5.0,
) -> pd.DataFrame:
    """Load calibrated per-frame aperture photometry for repeatability tests."""
    root = Path(project_root).expanduser().resolve()
    index_path = root / "result" / "step7_forced_phot" / "photometry_index.csv"
    zp_path = root / "result" / "cmd_zeropoint" / "frame_zeropoint.csv"
    if not index_path.is_file():
        raise FileNotFoundError(f"Step 7 index is missing: {index_path}")
    if not zp_path.is_file():
        raise FileNotFoundError(f"Step 10 frame zeropoints are missing: {zp_path}")
    index = pd.read_csv(index_path)
    zp = pd.read_csv(zp_path)
    index["filter_key"] = index["filter"].map(_normal_filter)
    zp["filter_key"] = zp["filter"].map(_normal_filter)
    zp["zp_frame"] = _num(zp["zp_frame"])
    merged = index.merge(
        zp[["file", "filter_key", "zp_frame"]],
        on=["file", "filter_key"],
        how="inner",
    )
    if "status" in merged.columns:
        merged = merged[merged["status"].astype(str).str.lower().eq("ok")]
    if "wcs_ok" in merged.columns:
        merged = merged[_to_bool(merged["wcs_ok"])]

    usecols = {
        "ID",
        "master_id",
        "mag_inst",
        "mag_err",
        "snr",
        "bad_phot_flag",
        "is_saturated",
        "is_nonlinear",
        "off_frame_flag",
    }
    frames: list[pd.DataFrame] = []
    for _, row in merged.iterrows():
        path = _resolve_step7_path(row, root)
        phot = pd.read_csv(path, sep="\t", usecols=lambda col: col in usecols)
        id_col = "ID" if "ID" in phot.columns else "master_id"
        if id_col not in phot.columns:
            continue
        phot = phot.rename(columns={id_col: "ID"})
        phot["file"] = row["file"]
        phot["filter"] = row["filter_key"]
        phot["zp_frame"] = float(row["zp_frame"])
        phot["mag_inst"] = _num(phot["mag_inst"])
        phot["mag_err"] = _num(phot["mag_err"]) if "mag_err" in phot.columns else np.nan
        phot["snr"] = _num(phot["snr"]) if "snr" in phot.columns else np.nan
        phot["mag_cal_frame"] = phot["mag_inst"] + phot["zp_frame"]
        good = np.isfinite(phot["mag_cal_frame"].to_numpy(float))
        if "snr" in phot.columns and np.isfinite(float(min_snr)):
            good &= phot["snr"].to_numpy(float) >= float(min_snr)
        for flag in ("bad_phot_flag", "is_saturated", "is_nonlinear", "off_frame_flag"):
            if flag in phot.columns:
                good &= ~_to_bool(phot[flag]).to_numpy(bool)
        frames.append(phot.loc[good, ["file", "filter", "ID", "mag_cal_frame", "mag_err", "snr"]])

    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def summarize_repeatability(
    observations: pd.DataFrame,
    *,
    min_frames: int = 3,
    bin_width: float = 0.5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if observations.empty:
        return pd.DataFrame(), pd.DataFrame()
    rows = []
    for (filter_name, star_id), group in observations.groupby(["filter", "ID"], sort=False):
        mags = _finite(group["mag_cal_frame"])
        if mags.size < int(min_frames):
            continue
        rows.append(
            {
                "filter": filter_name,
                "ID": star_id,
                "n_frames": int(mags.size),
                "mag_median": float(np.median(mags)),
                "mag_std": float(np.std(mags, ddof=1)) if mags.size > 1 else float("nan"),
                "mag_mad_sigma": _mad_sigma(mags),
                "mag_err_median": float(np.nanmedian(group["mag_err"])),
                "snr_median": float(np.nanmedian(group["snr"])),
            }
        )
    stars = pd.DataFrame(rows)
    if stars.empty:
        return stars, pd.DataFrame()
    stars["mag_bin"] = (
        np.floor(stars["mag_median"].astype(float) / float(bin_width)) * float(bin_width)
        + 0.5 * float(bin_width)
    ).round(3)
    bins = (
        stars.groupby(["filter", "mag_bin"], as_index=False)
        .agg(
            n_stars=("ID", "count"),
            mag_median=("mag_median", "median"),
            repeatability_std_median=("mag_std", "median"),
            repeatability_mad_median=("mag_mad_sigma", "median"),
            mag_err_median=("mag_err_median", "median"),
            snr_median=("snr_median", "median"),
            n_frames_median=("n_frames", "median"),
        )
        .sort_values(["filter", "mag_bin"])
    )
    return stars, bins


def _write_repeatability_plot(bins: pd.DataFrame, output_dir: Path) -> None:
    if bins.empty:
        return
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for filter_name, group in bins.groupby("filter", sort=True):
        group = group.sort_values("mag_bin")
        ax.plot(
            group["mag_median"],
            group["repeatability_mad_median"],
            marker="o",
            linewidth=1.5,
            label=f"{filter_name} robust RMS",
        )
        ax.plot(
            group["mag_median"],
            group["mag_err_median"],
            linestyle="--",
            linewidth=1.2,
            alpha=0.7,
            label=f"{filter_name} median formal err",
        )
    ax.set(xlabel="Median calibrated magnitude", ylabel="Magnitude scatter")
    ax.set_yscale("log")
    ax.grid(alpha=0.25, which="both")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(output_dir / "repeatability_rms_vs_magnitude.png", dpi=180)
    plt.close(fig)


def compute_zeropoint_residuals(
    calibrators: pd.DataFrame,
    coefficients: pd.DataFrame,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for _, coeff in coefficients.iterrows():
        filter_name = _normal_filter(coeff["filter"])
        delta_col = f"delta_{filter_name}"
        if delta_col not in calibrators.columns:
            continue
        zp = float(coeff["zp"])
        ct = float(coeff.get("ct", 0.0))
        color_name = str(coeff.get("color_col", "none"))
        if color_name == "none":
            color = np.zeros(len(calibrators), dtype=float)
        else:
            color_col = f"color_{color_name}"
            if color_col not in calibrators.columns:
                continue
            color = _num(calibrators[color_col]).to_numpy(float)
        delta = _num(calibrators[delta_col]).to_numpy(float)
        residual = delta - (zp + ct * color)
        ref_col = f"ref_{filter_name}"
        ref_mag = (
            _num(calibrators[ref_col]).to_numpy(float)
            if ref_col in calibrators.columns
            else np.full(len(calibrators), np.nan)
        )
        snr_col = f"snr_{filter_name}"
        snr = (
            _num(calibrators[snr_col]).to_numpy(float)
            if snr_col in calibrators.columns
            else np.full(len(calibrators), np.nan)
        )
        table = pd.DataFrame(
            {
                "filter": filter_name,
                "ID": calibrators["ID"] if "ID" in calibrators.columns else np.arange(len(calibrators)),
                "ref_mag": ref_mag,
                "color": color,
                "delta": delta,
                "residual_mag": residual,
                "snr": snr,
                "zp": zp,
                "ct": ct,
                "color_col": color_name,
                "fit_scatter_rms": coeff.get("scatter_rms", np.nan),
            }
        )
        table = table[np.isfinite(table["residual_mag"])]
        rows.append(table)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def mark_zeropoint_inliers(residuals: pd.DataFrame) -> pd.DataFrame:
    if residuals.empty:
        return residuals.copy()
    out = residuals.copy()
    out["robust_inlier"] = False
    for filter_name, group in out.groupby("filter", sort=False):
        mask = _robust_inlier_mask(group["residual_mag"])
        out.loc[group.index, "robust_inlier"] = mask
    return out


def summarize_zeropoint_residuals(residuals: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for filter_name, group in residuals.groupby("filter", sort=True):
        raw = _finite(group["residual_mag"])
        if "robust_inlier" in group.columns:
            clipped = _finite(group.loc[group["robust_inlier"].astype(bool), "residual_mag"])
        else:
            clipped = raw[_robust_inlier_mask(raw)]
        rows.append(
            {
                "filter": filter_name,
                "n_raw": int(raw.size),
                "n_clipped": int(clipped.size),
                "residual_median_raw": float(np.median(raw)) if raw.size else np.nan,
                "residual_rms_raw": _rms(raw),
                "residual_mad_sigma_raw": _mad_sigma(raw),
                "residual_p95_abs_raw": float(np.percentile(np.abs(raw), 95)) if raw.size else np.nan,
                "residual_median_clipped": float(np.median(clipped)) if clipped.size else np.nan,
                "residual_rms_clipped": _rms(clipped),
                "residual_mad_sigma_clipped": _mad_sigma(clipped),
                "residual_p95_abs_clipped": (
                    float(np.percentile(np.abs(clipped), 95)) if clipped.size else np.nan
                ),
                "zp": float(group["zp"].iloc[0]),
                "ct": float(group["ct"].iloc[0]),
                "color_col": group["color_col"].iloc[0],
                "fit_scatter_rms": float(group["fit_scatter_rms"].iloc[0])
                if "fit_scatter_rms" in group.columns
                else np.nan,
            }
        )
    return pd.DataFrame(rows)


def load_zeropoint_tables(project_root: str | Path) -> dict[str, pd.DataFrame]:
    root = Path(project_root).expanduser().resolve()
    zp_dir = root / "result" / "cmd_zeropoint"
    calibrator_path = zp_dir / "gaia_sdss_calibrator_by_ID.csv"
    coeff_path = zp_dir / "zp_fit_coefficients.csv"
    frame_path = zp_dir / "frame_zeropoint.csv"
    index_path = root / "result" / "step7_forced_phot" / "photometry_index.csv"
    if not calibrator_path.is_file():
        raise FileNotFoundError(f"Calibrator table is missing: {calibrator_path}")
    if not coeff_path.is_file():
        raise FileNotFoundError(f"ZP coefficient table is missing: {coeff_path}")
    calibrators = pd.read_csv(calibrator_path)
    coefficients = pd.read_csv(coeff_path)
    residuals = mark_zeropoint_inliers(compute_zeropoint_residuals(calibrators, coefficients))
    residual_summary = summarize_zeropoint_residuals(residuals)

    frame_summary = pd.DataFrame()
    if frame_path.is_file():
        frame = pd.read_csv(frame_path)
        frame["filter"] = frame["filter"].map(_normal_filter)
        for column in ("zp_frame", "zp_scatter", "n_ref", "snr_med"):
            if column in frame.columns:
                frame[column] = _num(frame[column])
        if index_path.is_file():
            index = pd.read_csv(index_path)
            index["filter"] = index["filter"].map(_normal_filter)
            frame = frame.merge(
                index[["file", "filter", "fwhm_px"]] if "fwhm_px" in index.columns else index[["file", "filter"]],
                on=["file", "filter"],
                how="left",
            )
        rows = []
        for filter_name, group in frame.groupby("filter", sort=True):
            rows.append(
                {
                    "filter": filter_name,
                    "n_frames": int(len(group)),
                    "zp_frame_median": float(np.nanmedian(group["zp_frame"])),
                    "zp_frame_std": float(np.nanstd(group["zp_frame"], ddof=1)),
                    "zp_scatter_median": float(np.nanmedian(group["zp_scatter"])),
                    "n_ref_median": float(np.nanmedian(group["n_ref"])),
                    "snr_med_median": float(np.nanmedian(group["snr_med"])),
                }
            )
        frame_summary = pd.DataFrame(rows)
    else:
        frame = pd.DataFrame()

    return {
        "zeropoint_residuals": residuals,
        "zeropoint_residual_summary": residual_summary,
        "frame_zeropoints": frame,
        "frame_zeropoint_summary": frame_summary,
    }


def _write_zeropoint_plots(tables: dict[str, pd.DataFrame], output_dir: Path) -> None:
    residuals = tables["zeropoint_residuals"]
    frame = tables["frame_zeropoints"]
    if not residuals.empty:
        filters = sorted(residuals["filter"].unique())
        fig, axes = plt.subplots(1, len(filters), figsize=(5.0 * len(filters), 4.0), sharey=True)
        axes = np.atleast_1d(axes)
        for ax, filter_name in zip(axes, filters):
            group = residuals[residuals["filter"] == filter_name]
            inlier = (
                group["robust_inlier"].astype(bool).to_numpy()
                if "robust_inlier" in group.columns
                else np.ones(len(group), dtype=bool)
            )
            ax.scatter(
                group.loc[~inlier, "ref_mag"],
                group.loc[~inlier, "residual_mag"],
                s=8,
                alpha=0.25,
                color="0.65",
                label="clipped",
            )
            ax.scatter(
                group.loc[inlier, "ref_mag"],
                group.loc[inlier, "residual_mag"],
                s=8,
                alpha=0.45,
                label="inlier",
            )
            ax.axhline(0.0, color="black", linewidth=1)
            ax.set_title(str(filter_name))
            ax.set_xlabel("Reference magnitude")
            ax.grid(alpha=0.25)
            ax.legend(fontsize=8)
        axes[0].set_ylabel("ZP fit residual (mag)")
        fig.tight_layout()
        fig.savefig(output_dir / "zeropoint_residuals_vs_mag.png", dpi=180)
        plt.close(fig)

    if not frame.empty and "fwhm_px" in frame.columns:
        fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))
        for filter_name, group in frame.groupby("filter", sort=True):
            median_zp = float(np.nanmedian(group["zp_frame"]))
            axes[0].scatter(group["fwhm_px"], group["zp_frame"] - median_zp, label=str(filter_name), alpha=0.8)
            axes[1].scatter(group["fwhm_px"], group["zp_scatter"], label=str(filter_name), alpha=0.8)
        axes[0].axhline(0.0, color="black", linewidth=1)
        axes[0].set(xlabel="Step 7 FWHM (px)", ylabel="Frame ZP - filter median")
        axes[1].set(xlabel="Step 7 FWHM (px)", ylabel="Per-frame ZP scatter")
        for ax in axes:
            ax.grid(alpha=0.25)
            ax.legend(title="Filter")
        fig.tight_layout()
        fig.savefig(output_dir / "zeropoint_frame_stability.png", dpi=180)
        plt.close(fig)


def _ast_loss_table(precision: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for filter_name, group in precision.groupby("filter", sort=True):
        by_condition = group.set_index("condition")
        if {"best", "worst"} <= set(by_condition.index):
            best = float(by_condition.loc["best", "m50"])
            worst = float(by_condition.loc["worst", "m50"])
            rows.append(
                {
                    "filter": filter_name,
                    "m50_best": best,
                    "m50_worst": worst,
                    "m50_loss_best_minus_worst": best - worst,
                }
            )
    return pd.DataFrame(rows)


def _write_report(
    output_dir: Path,
    *,
    config: CmdValidationConfig,
    ast_tables: dict[str, pd.DataFrame],
    repeatability_bins: pd.DataFrame,
    zp_tables: dict[str, pd.DataFrame],
) -> None:
    loss = _ast_loss_table(ast_tables["precision"])
    false_positive = ast_tables["false_positive_summary"]
    zp_summary = zp_tables["zeropoint_residual_summary"]
    lines = [
        "# CMD Benchmark Validation Report",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Scope",
        "",
        "This report summarizes CMD benchmark evidence from artificial-star recovery,",
        "forced-photometry repeatability, and Step 10 calibration residuals.",
        "",
        "## Completeness Loss",
        "",
    ]
    if not loss.empty:
        for _, row in loss.iterrows():
            lines.append(
                f"- {row['filter']}: m50 loss best-to-worst seeing = "
                f"{row['m50_loss_best_minus_worst']:.3f} mag"
            )
    else:
        lines.append("- Not available")

    lines.extend(["", "## False Detections", ""])
    if not false_positive.empty:
        total_false = int(false_positive["new_false_detections"].sum())
        total_injected = int(false_positive["n_injected"].sum())
        rate = 1000.0 * total_false / max(1, total_injected)
        lines.append(
            f"- New unmatched detections: {total_false} across {total_injected} injected stars "
            f"({rate:.2f} per 1000 injected)"
        )
    else:
        lines.append("- Not available")

    lines.extend(["", "## Repeatability", ""])
    if not repeatability_bins.empty:
        for filter_name, group in repeatability_bins.groupby("filter", sort=True):
            lines.append(
                f"- {filter_name}: median binned robust scatter = "
                f"{float(np.nanmedian(group['repeatability_mad_median'])):.4f} mag"
            )
    else:
        lines.append("- Not available")

    lines.extend(["", "## Zeropoint Residuals", ""])
    if not zp_summary.empty:
        for _, row in zp_summary.iterrows():
            lines.append(
                f"- {row['filter']}: clipped residual RMS = "
                f"{row['residual_rms_clipped']:.4f} mag, "
                f"MAD sigma = {row['residual_mad_sigma_clipped']:.4f} mag, "
                f"N={int(row['n_clipped'])}/{int(row['n_raw'])}, "
                f"Step10 fit scatter = {row['fit_scatter_rms']:.4f} mag"
            )
    else:
        lines.append("- Not available")

    lines.extend(
        [
            "",
            "## Outputs",
            "",
            "- ast_completeness_curves.png",
            "- ast_photometric_bias_scatter.png",
            "- ast_m50_vs_fwhm.png",
            "- ast_false_positive_rate.png",
            "- ast_crowding_m50.png",
            "- repeatability_rms_vs_magnitude.png",
            "- zeropoint_residuals_vs_mag.png",
            "- zeropoint_frame_stability.png",
            "",
            "## Configuration",
            "",
            "```json",
            json.dumps(asdict(config), indent=2, sort_keys=True),
            "```",
            "",
        ]
    )
    (output_dir / "README.md").write_text("\n".join(lines), encoding="utf-8")


def run_cmd_validation(
    config: CmdValidationConfig,
    *,
    batch_override: str | Path | None = None,
    project_override: str | Path | None = None,
    output_override: str | Path | None = None,
) -> Path:
    repo_root = Path(__file__).resolve().parents[2]
    batch_root = Path(batch_override or config.batch_root).expanduser().resolve()
    project_root = Path(project_override or config.project_root).expanduser().resolve()
    output_root = Path(output_override or config.output_root).expanduser()
    if not output_root.is_absolute():
        output_root = repo_root / output_root
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    ast_tables = collect_ast_tables(batch_root)
    for name, table in ast_tables.items():
        table.to_csv(output_root / f"{name}.csv", index=False)
    _ast_loss_table(ast_tables["precision"]).to_csv(output_root / "ast_m50_loss.csv", index=False)
    _write_ast_plots(ast_tables, output_root)

    observations = load_repeatability_observations(
        project_root,
        min_snr=float(config.repeatability_min_snr),
    )
    observations.to_csv(output_root / "repeatability_observations.csv", index=False)
    repeatability_stars, repeatability_bins = summarize_repeatability(
        observations,
        min_frames=int(config.repeatability_min_frames),
        bin_width=float(config.magnitude_bin_width),
    )
    repeatability_stars.to_csv(output_root / "repeatability_star_summary.csv", index=False)
    repeatability_bins.to_csv(output_root / "repeatability_by_magnitude.csv", index=False)
    _write_repeatability_plot(repeatability_bins, output_root)

    zp_tables = load_zeropoint_tables(project_root)
    for name, table in zp_tables.items():
        table.to_csv(output_root / f"{name}.csv", index=False)
    _write_zeropoint_plots(zp_tables, output_root)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "batch_root": str(batch_root),
        "project_root": str(project_root),
        "effective_config": asdict(config),
    }
    (output_root / "cmd_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    _write_report(
        output_root,
        config=config,
        ast_tables=ast_tables,
        repeatability_bins=repeatability_bins,
        zp_tables=zp_tables,
    )
    return output_root
