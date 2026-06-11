"""Batch orchestration for IRAF/PyRAF daofind cross-checks."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from apex.benchmark.iraf_crosscheck import (
    IRAFCrosscheckConfig,
    add_iraf_calibrated_equivalent_columns,
    run_iraf_crosscheck,
)


@dataclass
class IRAFCrosscheckBatchConfig:
    project_root: str = ""
    output_root: str = "benchmark/runs/iraf_daofind_batch"
    filters: list[str] = field(default_factory=lambda: ["g", "r", "i"])
    threshold_grid: list[float] = field(default_factory=lambda: [12.0, 9.0, 7.0, 5.0])
    mode: str = "daofind_phot"
    min_snr: float = 20.0
    zmag: float = 25.0
    daofind_max_sources: int = 2500
    daofind_max_ratio_to_apex: float = 2.0
    match_radius_px: float = 2.0
    limit: int | None = None
    overwrite: bool = False
    resume: bool = True
    runtime_cmd: list[str] = field(default_factory=list)


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _normal_filter(value: Any) -> str:
    return str(value).strip().lower().rstrip("'")


def _to_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series.fillna(False).astype(bool)
    return series.astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y"})


def _num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _resolve_frame_path(project_root: Path, filename: str) -> Path:
    direct = project_root / filename
    if direct.is_file():
        return direct.resolve()
    matches = [path for path in project_root.rglob(filename) if path.is_file()]
    if len(matches) == 1:
        return matches[0].resolve()
    if not matches:
        raise FileNotFoundError(f"FITS frame is missing: {project_root / filename}")
    raise RuntimeError(f"Frame name is ambiguous under {project_root}: {filename}")


def _resolve_step7_path(project_root: Path, row: pd.Series) -> Path:
    raw = row.get("path")
    if isinstance(raw, str) and raw.strip():
        path = Path(raw)
        if path.is_file():
            return path.resolve()
    filename = str(row["file"])
    candidate = project_root / "result" / "step7_forced_phot" / f"photometry_{filename}.tsv"
    if candidate.is_file():
        return candidate.resolve()
    matches = list((project_root / "result" / "step7_forced_phot").glob(f"*{filename}*.tsv"))
    if len(matches) == 1:
        return matches[0].resolve()
    raise FileNotFoundError(f"Step 7 photometry table is missing for {filename}")


def select_iraf_batch_frames(
    project_root: str | Path,
    *,
    filters: list[str] | tuple[str, ...] = ("g", "r", "i"),
) -> pd.DataFrame:
    """Select all calibrated OK/WCS frames for IRAF daofind comparison."""
    root = Path(project_root).expanduser().resolve()
    index_path = root / "result" / "step7_forced_phot" / "photometry_index.csv"
    zp_path = root / "result" / "cmd_zeropoint" / "frame_zeropoint.csv"
    if not index_path.is_file():
        raise FileNotFoundError(f"Step 7 photometry index is missing: {index_path}")
    if not zp_path.is_file():
        raise FileNotFoundError(f"Step 10 frame zeropoints are missing: {zp_path}")

    index = pd.read_csv(index_path)
    zeropoints = pd.read_csv(zp_path)
    required_index = {"file", "filter", "status", "fwhm_px"}
    required_zp = {"file", "filter", "zp_frame"}
    if not required_index <= set(index.columns):
        missing = sorted(required_index - set(index.columns))
        raise ValueError(f"Step 7 index is missing columns: {', '.join(missing)}")
    if not required_zp <= set(zeropoints.columns):
        missing = sorted(required_zp - set(zeropoints.columns))
        raise ValueError(f"Step 10 table is missing columns: {', '.join(missing)}")

    index = index.copy()
    zeropoints = zeropoints.copy()
    index["filter_key"] = index["filter"].map(_normal_filter)
    zeropoints["filter_key"] = zeropoints["filter"].map(_normal_filter)
    wanted = {_normal_filter(value) for value in filters}
    index = index[index["filter_key"].isin(wanted)]
    index = index[index["status"].astype(str).str.strip().str.lower().eq("ok")]
    if "wcs_ok" in index.columns:
        index = index[_to_bool(index["wcs_ok"])]
    index["fwhm_px"] = _num(index["fwhm_px"])
    index = index[np.isfinite(index["fwhm_px"])]
    zeropoints["zp_frame"] = _num(zeropoints["zp_frame"])

    zp_cols = ["file", "filter_key", "zp_frame"]
    for optional in ("zp_scatter", "n_ref", "snr_med"):
        if optional in zeropoints.columns:
            zp_cols.append(optional)
    merged = index.merge(
        zeropoints[zp_cols],
        on=["file", "filter_key"],
        how="inner",
        validate="one_to_one",
    )
    merged = merged[np.isfinite(merged["zp_frame"])].copy()
    rows: list[dict[str, Any]] = []
    for _, row in merged.iterrows():
        out = row.to_dict()
        out["input_fits"] = str(_resolve_frame_path(root, str(row["file"])))
        out["step7_tsv"] = str(_resolve_step7_path(root, row))
        rows.append(out)
    selected = pd.DataFrame(rows)
    if selected.empty:
        raise RuntimeError("No calibrated OK/WCS frames selected for IRAF cross-check")
    return selected.sort_values(["filter_key", "file"], kind="stable").reset_index(drop=True)


def _output_root(path: str | Path) -> Path:
    out = Path(path).expanduser()
    if not out.is_absolute():
        out = _repo_root() / out
    return out.resolve()


def _frame_output_dir(output_root: Path, row: pd.Series) -> Path:
    return output_root / "frames" / str(row["filter_key"]) / Path(str(row["file"])).stem


def _read_frame_summary(frame_dir: Path) -> pd.DataFrame:
    path = frame_dir / "daofind_phot" / "daofind_threshold_summary.csv"
    if not path.is_file():
        return pd.DataFrame()
    return pd.read_csv(path)


def _threshold_tag(value: float) -> str:
    return str(float(value)).replace(".", "p")


def _mad_sigma(values: Any) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    med = float(np.median(arr))
    return float(1.4826 * np.median(np.abs(arr - med)))


def _median_or_nan(values: Any) -> float:
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    return float(np.median(arr)) if arr.size else float("nan")


def collect_batch_summary(output_root: str | Path, selected: pd.DataFrame) -> pd.DataFrame:
    root = Path(output_root).expanduser().resolve()
    rows: list[pd.DataFrame] = []
    for _, row in selected.iterrows():
        frame_dir = _frame_output_dir(root, row)
        frame_summary = _read_frame_summary(frame_dir)
        if frame_summary.empty:
            continue
        frame_summary = frame_summary.copy()
        frame_summary.insert(0, "run_dir", str(frame_dir))
        frame_summary.insert(0, "step7_tsv", row["step7_tsv"])
        frame_summary.insert(0, "input_fits", row["input_fits"])
        frame_summary.insert(0, "zp_frame", row["zp_frame"])
        frame_summary.insert(0, "fwhm_px", row["fwhm_px"])
        frame_summary.insert(0, "file", row["file"])
        frame_summary.insert(0, "filter", row["filter_key"])
        for source_col, output_col in (
            ("n_master", "apex_n_master"),
            ("n_detected", "apex_n_detected"),
            ("n_forced", "apex_n_forced"),
            ("n_valid_phot", "apex_n_valid_phot"),
            ("detected_rate", "apex_detected_rate"),
            ("forced_rate", "apex_forced_rate"),
        ):
            if source_col in row.index:
                frame_summary[output_col] = row[source_col]
        rows.append(frame_summary)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def collect_apex_reference_observations(
    output_root: str | Path,
    selected: pd.DataFrame,
) -> pd.DataFrame:
    root = Path(output_root).expanduser().resolve()
    rows: list[pd.DataFrame] = []
    for _, row in selected.iterrows():
        path = _frame_output_dir(root, row) / "apex_reference_all.csv"
        if not path.is_file():
            continue
        table = pd.read_csv(path)
        if "master_id" not in table.columns or "mag_inst" not in table.columns:
            continue
        table = table.copy()
        table["filter"] = row["filter_key"]
        table["file"] = row["file"]
        table["zp_frame"] = float(row["zp_frame"])
        table["apex_mag_cal_frame"] = pd.to_numeric(table["mag_inst"], errors="coerce") + float(row["zp_frame"])
        rows.append(table)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def collect_paired_repeatability_observations(
    output_root: str | Path,
    selected: pd.DataFrame,
    *,
    thresholds: list[float] | tuple[float, ...],
    zmag: float = 25.0,
) -> pd.DataFrame:
    root = Path(output_root).expanduser().resolve()
    rows: list[pd.DataFrame] = []
    for _, row in selected.iterrows():
        frame_dir = _frame_output_dir(root, row)
        for threshold in thresholds:
            match_path = frame_dir / "daofind_phot" / f"matches_thr_{_threshold_tag(float(threshold))}.csv"
            if not match_path.is_file():
                continue
            matches = pd.read_csv(match_path)
            required = {"master_id", "mag_inst", "iraf_mag"}
            if not required <= set(matches.columns):
                continue
            matches = matches.copy()
            matches["filter"] = row["filter_key"]
            matches["file"] = row["file"]
            matches["threshold"] = float(threshold)
            matches["zp_frame"] = float(row["zp_frame"])
            matches["apex_mag_cal_frame"] = pd.to_numeric(matches["mag_inst"], errors="coerce") + float(row["zp_frame"])
            matches["iraf_mag"] = pd.to_numeric(matches["iraf_mag"], errors="coerce")
            matches = add_iraf_calibrated_equivalent_columns(
                matches,
                zmag=float(zmag),
                frame_zeropoint=float(row["zp_frame"]),
            )
            finite = np.isfinite(matches["apex_mag_cal_frame"].to_numpy(float)) & np.isfinite(
                matches["iraf_mag"].to_numpy(float)
            )
            offset = _median_or_nan(matches.loc[finite, "apex_mag_cal_frame"] - matches.loc[finite, "iraf_mag"])
            matches["iraf_to_apex_frame_offset"] = offset
            matches["iraf_mag_aligned_frame"] = matches["iraf_mag"] + offset
            rows.append(matches)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def summarize_apex_reference_repeatability(
    observations: pd.DataFrame,
    *,
    min_frames: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if observations.empty:
        return pd.DataFrame(), pd.DataFrame()
    work = observations.copy()
    work["master_id"] = work["master_id"].astype(str)
    work["apex_mag_cal_frame"] = pd.to_numeric(work["apex_mag_cal_frame"], errors="coerce")
    rows = []
    for (filter_name, master_id), group in work.groupby(["filter", "master_id"], sort=False):
        mags = group["apex_mag_cal_frame"].to_numpy(float)
        mags = mags[np.isfinite(mags)]
        if mags.size < int(min_frames):
            continue
        rows.append(
            {
                "filter": filter_name,
                "master_id": master_id,
                "n_apex_frames": int(mags.size),
                "apex_mag_median": float(np.median(mags)),
                "apex_repeat_mad_sigma": _mad_sigma(mags),
                "apex_repeat_std": float(np.std(mags, ddof=1)) if mags.size > 1 else float("nan"),
            }
        )
    by_star = pd.DataFrame(rows)
    if by_star.empty:
        return by_star, pd.DataFrame()
    summary = (
        by_star.groupby("filter", as_index=False)
        .agg(
            n_apex_reference_stars=("master_id", "count"),
            n_apex_frames_median=("n_apex_frames", "median"),
            apex_reference_repeat_mad_median=("apex_repeat_mad_sigma", "median"),
            apex_reference_repeat_std_median=("apex_repeat_std", "median"),
        )
        .sort_values("filter")
    )
    return by_star, summary


def summarize_paired_repeatability(
    observations: pd.DataFrame,
    *,
    min_frames: int = 3,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if observations.empty:
        return pd.DataFrame(), pd.DataFrame()
    work = observations.copy()
    work["master_id"] = work["master_id"].astype(str)
    for column in (
        "apex_mag_cal_frame",
        "iraf_mag",
        "iraf_mag_cal_apcorr_zp",
        "iraf_mag_aligned_frame",
        "mag_err",
        "iraf_merr",
    ):
        if column in work.columns:
            work[column] = pd.to_numeric(work[column], errors="coerce")
    rows = []
    for (filter_name, threshold, master_id), group in work.groupby(["filter", "threshold", "master_id"], sort=False):
        apex = group["apex_mag_cal_frame"].to_numpy(float)
        iraf_raw = group["iraf_mag"].to_numpy(float)
        iraf = group["iraf_mag_aligned_frame"].to_numpy(float)
        if "iraf_mag_cal_apcorr_zp" in group.columns:
            iraf_cal = group["iraf_mag_cal_apcorr_zp"].to_numpy(float)
        else:
            iraf_cal = np.full(len(group), np.nan, dtype=float)
        good = np.isfinite(apex) & np.isfinite(iraf_raw) & np.isfinite(iraf)
        apex = apex[good]
        iraf_raw = iraf_raw[good]
        iraf = iraf[good]
        iraf_cal = iraf_cal[good]
        if apex.size < int(min_frames):
            continue
        cal_good = np.isfinite(iraf_cal)
        iraf_cal_finite = iraf_cal[cal_good]
        apex_mad = _mad_sigma(apex)
        iraf_raw_mad = _mad_sigma(iraf_raw)
        iraf_mad = _mad_sigma(iraf)
        iraf_cal_mad = _mad_sigma(iraf_cal_finite) if iraf_cal_finite.size >= int(min_frames) else float("nan")
        rows.append(
            {
                "filter": filter_name,
                "threshold": float(threshold),
                "master_id": master_id,
                "n_paired_frames": int(apex.size),
                "n_iraf_calibrated_frames": int(iraf_cal_finite.size),
                "apex_mag_median": float(np.median(apex)),
                "iraf_mag_median": float(np.median(iraf)),
                "iraf_calibrated_mag_median": (
                    float(np.median(iraf_cal_finite)) if iraf_cal_finite.size else float("nan")
                ),
                "apex_repeat_mad_sigma": apex_mad,
                "iraf_raw_repeat_mad_sigma": iraf_raw_mad,
                "iraf_calibrated_repeat_mad_sigma": iraf_cal_mad,
                "iraf_repeat_mad_sigma": iraf_mad,
                "apex_repeat_std": float(np.std(apex, ddof=1)) if apex.size > 1 else float("nan"),
                "iraf_raw_repeat_std": float(np.std(iraf_raw, ddof=1)) if iraf_raw.size > 1 else float("nan"),
                "iraf_calibrated_repeat_std": (
                    float(np.std(iraf_cal_finite, ddof=1)) if iraf_cal_finite.size > 1 else float("nan")
                ),
                "iraf_repeat_std": float(np.std(iraf, ddof=1)) if iraf.size > 1 else float("nan"),
                "iraf_raw_over_apex_repeat_mad": (
                    float(iraf_raw_mad / apex_mad) if np.isfinite(apex_mad) and apex_mad > 0 else float("nan")
                ),
                "iraf_calibrated_over_apex_repeat_mad": (
                    float(iraf_cal_mad / apex_mad)
                    if np.isfinite(iraf_cal_mad) and np.isfinite(apex_mad) and apex_mad > 0
                    else float("nan")
                ),
                "iraf_over_apex_repeat_mad": (
                    float(iraf_mad / apex_mad) if np.isfinite(apex_mad) and apex_mad > 0 else float("nan")
                ),
                "apex_mag_err_median": _median_or_nan(group.get("mag_err", pd.Series(dtype=float))),
                "iraf_merr_median": _median_or_nan(group.get("iraf_merr", pd.Series(dtype=float))),
            }
        )
    by_star = pd.DataFrame(rows)
    if by_star.empty:
        return by_star, pd.DataFrame()
    summary = (
        by_star.groupby(["filter", "threshold"], as_index=False)
        .agg(
            n_paired_stars=("master_id", "count"),
            n_paired_frames_median=("n_paired_frames", "median"),
            n_iraf_calibrated_frames_median=("n_iraf_calibrated_frames", "median"),
            apex_repeat_mad_median=("apex_repeat_mad_sigma", "median"),
            iraf_raw_repeat_mad_median=("iraf_raw_repeat_mad_sigma", "median"),
            iraf_calibrated_repeat_mad_median=("iraf_calibrated_repeat_mad_sigma", "median"),
            iraf_repeat_mad_median=("iraf_repeat_mad_sigma", "median"),
            apex_repeat_std_median=("apex_repeat_std", "median"),
            iraf_raw_repeat_std_median=("iraf_raw_repeat_std", "median"),
            iraf_calibrated_repeat_std_median=("iraf_calibrated_repeat_std", "median"),
            iraf_repeat_std_median=("iraf_repeat_std", "median"),
            iraf_raw_over_apex_repeat_mad_median=("iraf_raw_over_apex_repeat_mad", "median"),
            iraf_calibrated_over_apex_repeat_mad_median=("iraf_calibrated_over_apex_repeat_mad", "median"),
            iraf_over_apex_repeat_mad_median=("iraf_over_apex_repeat_mad", "median"),
            apex_formal_error_median=("apex_mag_err_median", "median"),
            iraf_formal_error_median=("iraf_merr_median", "median"),
        )
        .sort_values(["filter", "threshold"])
    )
    return by_star, summary


def summarize_batch_thresholds(summary: pd.DataFrame) -> pd.DataFrame:
    if summary.empty:
        return pd.DataFrame()
    work = summary.copy()
    for column in (
        "threshold",
        "n_detected",
        "n_matched",
        "recall_vs_apex_reference",
        "matched_fraction_of_iraf",
        "mad_sigma_delta_mag_zp_aligned",
        "p95_abs_delta_mag_zp_aligned",
        "n_apex_reference",
    ):
        if column in work.columns:
            work[column] = pd.to_numeric(work[column], errors="coerce")
    ok = work[work["status"].astype(str).eq("ok")].copy()
    if ok.empty:
        return pd.DataFrame()
    grouped = (
        ok.groupby(["filter", "threshold"], as_index=False)
        .agg(
            n_frames=("file", "nunique"),
            n_detected_median=("n_detected", "median"),
            n_detected_p90=("n_detected", lambda s: float(np.nanpercentile(s, 90))),
            n_matched_median=("n_matched", "median"),
            recall_median=("recall_vs_apex_reference", "median"),
            matched_fraction_median=("matched_fraction_of_iraf", "median"),
            residual_mad_median=("mad_sigma_delta_mag_zp_aligned", "median"),
            residual_p95_median=("p95_abs_delta_mag_zp_aligned", "median"),
            apex_reference_median=("n_apex_reference", "median"),
        )
        .sort_values(["filter", "threshold"], ascending=[True, False])
    )
    return grouped


def recommend_thresholds(summary_by_threshold: pd.DataFrame) -> pd.DataFrame:
    if summary_by_threshold.empty:
        return pd.DataFrame()
    rows = []
    for filter_name, group in summary_by_threshold.groupby("filter", sort=True):
        candidates = group.copy()
        if len(candidates) == 1:
            chosen = candidates.iloc[0]
            reason = "only tested threshold; not an independently optimized recommendation"
            row = chosen.to_dict()
            row["reason"] = reason
            rows.append(row)
            continue
        high_recall = candidates[candidates["recall_median"] >= 0.95]
        if not high_recall.empty:
            chosen = high_recall.sort_values(
                ["matched_fraction_median", "threshold"],
                ascending=[False, False],
            ).iloc[0]
            reason = "highest purity among thresholds with median recall >= 0.95"
        else:
            chosen = candidates.sort_values(
                ["recall_median", "matched_fraction_median", "threshold"],
                ascending=[False, False, False],
            ).iloc[0]
            reason = "highest median recall; no threshold reached 0.95"
        row = chosen.to_dict()
        row["reason"] = reason
        rows.append(row)
    return pd.DataFrame(rows)


def summarize_batch_filters(summary: pd.DataFrame) -> pd.DataFrame:
    if summary.empty:
        return pd.DataFrame()
    work = summary[summary["status"].astype(str).eq("ok")].copy()
    numeric_columns = (
        "n_detected",
        "n_matched",
        "n_apex_reference",
        "recall_vs_apex_reference",
        "matched_fraction_of_iraf",
        "mad_sigma_delta_mag_zp_aligned",
        "p95_abs_delta_mag_zp_aligned",
        "median_apex_mag_err",
        "median_iraf_merr",
        "apex_n_master",
        "apex_n_detected",
        "apex_n_forced",
        "apex_n_valid_phot",
        "apex_detected_rate",
        "apex_forced_rate",
    )
    for column in numeric_columns:
        if column in work.columns:
            work[column] = pd.to_numeric(work[column], errors="coerce")
    if work.empty:
        return pd.DataFrame()

    rows = []
    for filter_name, group in work.groupby("filter", sort=True):
        row = {
            "filter": filter_name,
            "n_frames": int(group["file"].nunique()),
            "iraf_detected_median": float(np.nanmedian(group["n_detected"])),
            "matched_median": float(np.nanmedian(group["n_matched"])),
            "apex_reference_median": float(np.nanmedian(group["n_apex_reference"])),
            "recall_median": float(np.nanmedian(group["recall_vs_apex_reference"])),
            "matched_fraction_median": float(np.nanmedian(group["matched_fraction_of_iraf"])),
            "residual_mad_median": float(np.nanmedian(group["mad_sigma_delta_mag_zp_aligned"])),
            "residual_p95_median": float(np.nanmedian(group["p95_abs_delta_mag_zp_aligned"])),
            "apex_formal_error_median": float(np.nanmedian(group["median_apex_mag_err"])),
            "iraf_formal_error_median": float(np.nanmedian(group["median_iraf_merr"])),
        }
        matched_fraction = row["matched_fraction_median"]
        if matched_fraction < 0.80:
            row["detection_qc"] = "likely over-detection; raise threshold"
        elif matched_fraction < 0.90:
            row["detection_qc"] = "review threshold purity"
        else:
            row["detection_qc"] = "balanced at tested threshold"
        for column in (
            "apex_n_master",
            "apex_n_detected",
            "apex_n_forced",
            "apex_n_valid_phot",
            "apex_detected_rate",
            "apex_forced_rate",
        ):
            row[f"{column}_median"] = (
                float(np.nanmedian(group[column])) if column in group.columns else float("nan")
            )
        rows.append(row)
    return pd.DataFrame(rows)


def _write_batch_plots(
    summary: pd.DataFrame,
    by_threshold: pd.DataFrame,
    by_filter: pd.DataFrame,
    paired_repeatability: pd.DataFrame,
    output_root: Path,
) -> None:
    if summary.empty:
        return
    ok = summary[summary["status"].astype(str).eq("ok")].copy()
    if ok.empty:
        return
    for column in (
        "threshold",
        "n_detected",
        "recall_vs_apex_reference",
        "matched_fraction_of_iraf",
        "mad_sigma_delta_mag_zp_aligned",
        "fwhm_px",
    ):
        if column in ok.columns:
            ok[column] = pd.to_numeric(ok[column], errors="coerce")

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for filter_name, group in by_threshold.groupby("filter", sort=True):
        ordered = group.sort_values("threshold", ascending=False)
        ax.plot(
            ordered["threshold"],
            ordered["n_detected_median"],
            marker="o",
            linewidth=1.6,
            label=str(filter_name),
        )
    ax.invert_xaxis()
    ax.set(xlabel="IRAF findpars.threshold", ylabel="Median detected sources per frame")
    ax.grid(alpha=0.25)
    ax.legend(title="Filter")
    fig.tight_layout()
    fig.savefig(output_root / "batch_detected_sources_vs_threshold.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for filter_name, group in by_threshold.groupby("filter", sort=True):
        ordered = group.sort_values("threshold", ascending=False)
        ax.plot(
            ordered["threshold"],
            ordered["recall_median"],
            marker="o",
            linewidth=1.6,
            label=f"{filter_name} recall",
        )
        ax.plot(
            ordered["threshold"],
            ordered["matched_fraction_median"],
            marker="s",
            linestyle="--",
            linewidth=1.2,
            label=f"{filter_name} matched fraction",
        )
    ax.invert_xaxis()
    ax.set(xlabel="IRAF findpars.threshold", ylabel="Median fraction")
    ax.set_ylim(0, 1.05)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(output_root / "batch_match_fractions_vs_threshold.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for filter_name, group in by_threshold.groupby("filter", sort=True):
        ordered = group.sort_values("threshold", ascending=False)
        ax.plot(
            ordered["threshold"],
            ordered["residual_mad_median"],
            marker="o",
            linewidth=1.6,
            label=str(filter_name),
        )
    ax.invert_xaxis()
    ax.set(xlabel="IRAF findpars.threshold", ylabel="Median ZP-aligned residual MAD (mag)")
    ax.grid(alpha=0.25)
    ax.legend(title="Filter")
    fig.tight_layout()
    fig.savefig(output_root / "batch_residual_mad_vs_threshold.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for filter_name, group in ok.groupby("filter", sort=True):
        ax.scatter(
            group["fwhm_px"],
            group["mad_sigma_delta_mag_zp_aligned"],
            s=28,
            alpha=0.7,
            label=str(filter_name),
        )
    ax.set(xlabel="APEX Step 7 FWHM (px)", ylabel="ZP-aligned residual MAD (mag)")
    ax.grid(alpha=0.25)
    ax.legend(title="Filter")
    fig.tight_layout()
    fig.savefig(output_root / "batch_residual_mad_vs_fwhm.png", dpi=180)
    plt.close(fig)

    if not by_filter.empty:
        labels = by_filter["filter"].astype(str).tolist()
        x = np.arange(len(labels), dtype=float)
        width = 0.2
        fig, ax = plt.subplots(figsize=(8.5, 5.2))
        ax.bar(
            x - 1.5 * width,
            by_filter["apex_n_valid_phot_median"],
            width,
            label="APEX valid photometry",
        )
        ax.bar(
            x - 0.5 * width,
            by_filter["apex_n_detected_median"],
            width,
            label="APEX direct detections",
        )
        ax.bar(
            x + 0.5 * width,
            by_filter["iraf_detected_median"],
            width,
            label="IRAF daofind+phot",
        )
        ax.bar(
            x + 1.5 * width,
            by_filter["matched_median"],
            width,
            label="IRAF matched to APEX",
        )
        ax.set_xticks(x, labels)
        ax.set(xlabel="Filter", ylabel="Median sources per frame")
        ax.grid(alpha=0.2, axis="y")
        ax.legend(fontsize=8, ncol=2)
        fig.tight_layout()
        fig.savefig(output_root / "batch_photometry_coverage_by_filter.png", dpi=180)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8.0, 5.2))
        ax.bar(
            x - width / 2,
            by_filter["recall_median"],
            width,
            label="IRAF recall vs APEX reference",
        )
        ax.bar(
            x + width / 2,
            by_filter["matched_fraction_median"],
            width,
            label="IRAF matched fraction",
        )
        ax.set_xticks(x, labels)
        ax.set(xlabel="Filter", ylabel="Median fraction", ylim=(0, 1.05))
        ax.grid(alpha=0.2, axis="y")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(output_root / "batch_recall_purity_by_filter.png", dpi=180)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8.0, 5.2))
        ax.bar(
            x - width / 2,
            by_filter["apex_formal_error_median"],
            width,
            label="APEX median formal error",
        )
        ax.bar(
            x + width / 2,
            by_filter["iraf_formal_error_median"],
            width,
            label="IRAF median MERR",
        )
        ax.scatter(
            x,
            by_filter["residual_mad_median"],
            color="black",
            marker="D",
            s=45,
            label="APEX-IRAF residual MAD",
            zorder=3,
        )
        ax.set_xticks(x, labels)
        ax.set(xlabel="Filter", ylabel="Magnitude error / scatter")
        ax.set_yscale("log")
        ax.grid(alpha=0.2, axis="y", which="both")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(output_root / "batch_formal_error_by_filter.png", dpi=180)
        plt.close(fig)

    if not paired_repeatability.empty:
        labels = paired_repeatability["filter"].astype(str).tolist()
        x = np.arange(len(labels), dtype=float)
        width = 0.22
        fig, ax = plt.subplots(figsize=(8.0, 5.2))
        ax.bar(
            x - 1.5 * width,
            paired_repeatability["apex_repeat_mad_median"],
            width,
            label="APEX repeatability MAD",
        )
        ax.bar(
            x - 0.5 * width,
            paired_repeatability["iraf_raw_repeat_mad_median"],
            width,
            label="IRAF raw repeatability MAD",
        )
        if "iraf_calibrated_repeat_mad_median" in paired_repeatability.columns:
            ax.bar(
                x + 0.5 * width,
                paired_repeatability["iraf_calibrated_repeat_mad_median"],
                width,
                label="IRAF apcorr+ZP-equivalent MAD",
            )
        ax.bar(
            x + 1.5 * width,
            paired_repeatability["iraf_repeat_mad_median"],
            width,
            label="IRAF frame-aligned repeatability MAD",
        )
        ax.set_xticks(x, labels)
        ax.set(xlabel="Filter", ylabel="Median per-star repeatability MAD (mag)")
        ax.set_yscale("log")
        ax.grid(alpha=0.2, axis="y", which="both")
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_root / "batch_repeatability_by_filter.png", dpi=180)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8.0, 5.2))
        ax.bar(
            x - width,
            paired_repeatability["iraf_raw_over_apex_repeat_mad_median"],
            width=width,
            label="IRAF raw / APEX",
        )
        if "iraf_calibrated_over_apex_repeat_mad_median" in paired_repeatability.columns:
            ax.bar(
                x,
                paired_repeatability["iraf_calibrated_over_apex_repeat_mad_median"],
                width=width,
                label="IRAF apcorr+ZP / APEX",
            )
        ax.bar(
            x + width,
            paired_repeatability["iraf_over_apex_repeat_mad_median"],
            width=width,
            label="IRAF aligned / APEX",
        )
        ax.axhline(1.0, color="black", linewidth=1)
        ax.set_xticks(x, labels)
        ax.set(xlabel="Filter", ylabel="Median scatter ratio")
        ax.grid(alpha=0.2, axis="y")
        ax.legend()
        fig.tight_layout()
        fig.savefig(output_root / "batch_repeatability_ratio_by_filter.png", dpi=180)
        plt.close(fig)


def _write_report(
    output_root: Path,
    *,
    config: IRAFCrosscheckBatchConfig,
    selected: pd.DataFrame,
    summary: pd.DataFrame,
    by_threshold: pd.DataFrame,
    by_filter: pd.DataFrame,
    paired_repeatability: pd.DataFrame,
    apex_reference_repeatability: pd.DataFrame,
    recommendations: pd.DataFrame,
) -> None:
    lines = [
        "# IRAF daofind Batch Cross-Check",
        "",
        f"Generated: {datetime.now(timezone.utc).isoformat()}",
        "",
        "## Scope",
        "",
        "This report runs IRAF/DAOPHOT `daofind+phot` for every selected APEX CMD frame",
        "and matches the results against APEX Step 7 sources.",
        "",
        "## Frames",
        "",
        f"- Selected frames: {len(selected)}",
    ]
    for filter_name, count in selected.groupby("filter_key").size().items():
        lines.append(f"- {filter_name}: {int(count)} frames")

    lines.extend(["", "## Threshold Summary", ""])
    if by_threshold.empty:
        lines.append("- Not available")
    else:
        for _, row in by_threshold.iterrows():
            lines.append(
                f"- {row['filter']} threshold={row['threshold']:.2f}: "
                f"frames={int(row['n_frames'])}, "
                f"detected median={row['n_detected_median']:.0f}, "
                f"recall median={row['recall_median']:.3f}, "
                f"matched fraction median={row['matched_fraction_median']:.3f}, "
                f"residual MAD median={row['residual_mad_median']:.4f} mag"
            )

    lines.extend(["", "## Threshold Selection", ""])
    if recommendations.empty:
        lines.append("- Not available")
    else:
        for _, row in recommendations.iterrows():
            lines.append(
                f"- {row['filter']}: threshold={row['threshold']:.2f} "
                f"(recall={row['recall_median']:.3f}, matched fraction={row['matched_fraction_median']:.3f}); "
                f"{row['reason']}"
            )

    lines.extend(["", "## APEX vs IRAF Coverage", ""])
    if by_filter.empty:
        lines.append("- Not available")
    else:
        for _, row in by_filter.iterrows():
            lines.append(
                f"- {row['filter']}: APEX valid median={row['apex_n_valid_phot_median']:.0f}, "
                f"APEX direct median={row['apex_n_detected_median']:.0f}, "
                f"APEX forced median={row['apex_n_forced_median']:.0f}, "
                f"IRAF detected median={row['iraf_detected_median']:.0f}, "
                f"IRAF matched median={row['matched_median']:.0f}; "
                f"QC={row['detection_qc']}"
            )

    lines.extend(["", "## Photometric Agreement", ""])
    if by_filter.empty:
        lines.append("- Not available")
    else:
        for _, row in by_filter.iterrows():
            lines.append(
                f"- {row['filter']}: residual MAD={row['residual_mad_median']:.4f} mag, "
                f"APEX formal error median={row['apex_formal_error_median']:.4f} mag, "
                f"IRAF MERR median={row['iraf_formal_error_median']:.4f} mag"
            )

    lines.extend(["", "## Empirical Repeatability", ""])
    if paired_repeatability.empty:
        lines.append("- Not available")
    else:
        for _, row in paired_repeatability.iterrows():
            lines.append(
                f"- {row['filter']}: paired stars={int(row['n_paired_stars'])}, "
                f"APEX repeat MAD={row['apex_repeat_mad_median']:.4f} mag, "
                f"IRAF raw repeat MAD={row['iraf_raw_repeat_mad_median']:.4f} mag "
                f"(raw/APEX={row['iraf_raw_over_apex_repeat_mad_median']:.2f}), "
                f"IRAF apcorr+ZP-equivalent repeat MAD={row['iraf_calibrated_repeat_mad_median']:.4f} mag "
                f"(apcorr+ZP/APEX={row['iraf_calibrated_over_apex_repeat_mad_median']:.2f}), "
                f"IRAF frame-aligned repeat MAD={row['iraf_repeat_mad_median']:.4f} mag "
                f"(aligned/APEX={row['iraf_over_apex_repeat_mad_median']:.2f})"
            )
    if not apex_reference_repeatability.empty:
        lines.append("")
        lines.append("APEX reference coverage:")
        for _, row in apex_reference_repeatability.iterrows():
            lines.append(
                f"- {row['filter']}: APEX high-SNR stars with >=3 frames={int(row['n_apex_reference_stars'])}, "
                f"repeat MAD={row['apex_reference_repeat_mad_median']:.4f} mag"
            )

    lines.extend(
        [
            "",
            "## Plots",
            "",
            "- batch_detected_sources_vs_threshold.png",
            "- batch_match_fractions_vs_threshold.png",
            "- batch_residual_mad_vs_threshold.png",
            "- batch_residual_mad_vs_fwhm.png",
            "- batch_photometry_coverage_by_filter.png",
            "- batch_recall_purity_by_filter.png",
            "- batch_formal_error_by_filter.png",
            "- batch_repeatability_by_filter.png",
            "- batch_repeatability_ratio_by_filter.png",
            "",
            "## Interpretation Notes",
            "",
            "- High threshold generally increases purity but may miss faint APEX sources.",
            "- Low threshold can inflate IRAF detections; source-count guards prevent runaway lower-threshold processing.",
            "- Residual scatter is measured after Step 10 frame-ZP alignment and median IRAF-to-APEX offset fitting.",
            "- Raw IRAF repeatability includes uncorrected frame-to-frame photometric offsets; frame-aligned IRAF isolates random measurement scatter.",
            "- IRAF apcorr+ZP-equivalent repeatability applies APEX's gain, aperture correction, and Step 10 frame ZP scale to IRAF `phot` measurements.",
            "- This is a detector-plus-photometry agreement benchmark, not an absolute truth catalog.",
            "",
            "## Configuration",
            "",
            "```json",
            json.dumps(asdict(config), indent=2, sort_keys=True),
            "```",
            "",
        ]
    )
    (output_root / "README.md").write_text("\n".join(lines), encoding="utf-8")


def run_iraf_crosscheck_batch(config: IRAFCrosscheckBatchConfig) -> Path:
    project_root = Path(config.project_root).expanduser().resolve()
    if not project_root.is_dir():
        raise FileNotFoundError(f"Project root does not exist: {project_root}")
    output_root = _output_root(config.output_root)
    if output_root.exists() and any(output_root.iterdir()) and not (config.overwrite or config.resume):
        raise FileExistsError(f"Batch output already exists: {output_root}")
    output_root.mkdir(parents=True, exist_ok=True)

    selected = select_iraf_batch_frames(project_root, filters=config.filters)
    if config.limit is not None:
        selected = selected.head(int(config.limit)).copy()
    selected.to_csv(output_root / "selected_frames.csv", index=False)

    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": str(project_root),
        "effective_config": asdict(config),
    }
    (output_root / "iraf_daofind_batch_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )

    summary = pd.DataFrame()
    for index, row in selected.reset_index(drop=True).iterrows():
        frame_dir = _frame_output_dir(output_root, row)
        completed_summary = frame_dir / "daofind_phot" / "daofind_threshold_summary.csv"
        if config.resume and completed_summary.is_file() and not config.overwrite:
            pass
        else:
            frame_config = IRAFCrosscheckConfig(
                input_fits=str(row["input_fits"]),
                step7_tsv=str(row["step7_tsv"]),
                output_root=str(frame_dir),
                mode=str(config.mode),
                filter_name=str(row["filter_key"]),
                fwhm_px=float(row["fwhm_px"]),
                min_snr=float(config.min_snr),
                zmag=float(config.zmag),
                daofind_threshold_grid=[float(value) for value in config.threshold_grid],
                daofind_max_sources=int(config.daofind_max_sources),
                daofind_max_ratio_to_apex=float(config.daofind_max_ratio_to_apex),
                daofind_match_radius_px=float(config.match_radius_px),
                overwrite=True,
                runtime_cmd=list(config.runtime_cmd),
            )
            run_iraf_crosscheck(frame_config)
        summary = collect_batch_summary(output_root, selected)
        summary.to_csv(output_root / "iraf_daofind_batch_summary.csv", index=False)
        by_threshold = summarize_batch_thresholds(summary)
        by_threshold.to_csv(output_root / "iraf_daofind_threshold_summary.csv", index=False)
        by_filter = summarize_batch_filters(summary)
        by_filter.to_csv(output_root / "iraf_daofind_filter_summary.csv", index=False)
        apex_obs = collect_apex_reference_observations(output_root, selected)
        apex_obs.to_csv(output_root / "apex_reference_repeatability_observations.csv", index=False)
        apex_reference_by_star, apex_reference_repeatability = summarize_apex_reference_repeatability(apex_obs)
        apex_reference_by_star.to_csv(output_root / "apex_reference_repeatability_by_star.csv", index=False)
        apex_reference_repeatability.to_csv(output_root / "apex_reference_repeatability_summary.csv", index=False)
        paired_obs = collect_paired_repeatability_observations(
            output_root,
            selected,
            thresholds=[float(value) for value in config.threshold_grid],
            zmag=float(config.zmag),
        )
        paired_obs.to_csv(output_root / "iraf_apex_paired_repeatability_observations.csv", index=False)
        paired_by_star, paired_repeatability = summarize_paired_repeatability(paired_obs)
        paired_by_star.to_csv(output_root / "iraf_apex_paired_repeatability_by_star.csv", index=False)
        paired_repeatability.to_csv(output_root / "iraf_apex_paired_repeatability_summary.csv", index=False)
        recommendations = recommend_thresholds(by_threshold)
        recommendations.to_csv(output_root / "iraf_daofind_recommended_thresholds.csv", index=False)
        _write_batch_plots(summary, by_threshold, by_filter, paired_repeatability, output_root)
        _write_report(
            output_root,
            config=config,
            selected=selected,
            summary=summary,
            by_threshold=by_threshold,
            by_filter=by_filter,
            paired_repeatability=paired_repeatability,
            apex_reference_repeatability=apex_reference_repeatability,
            recommendations=recommendations,
        )
        print(f"[IRAF batch] {index + 1}/{len(selected)} {row['file']} complete", flush=True)

    return output_root
