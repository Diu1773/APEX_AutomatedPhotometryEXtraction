"""Step 4 QC products, drawn once for both the window and the batch.

The window has built `step4_qc_summary.csv`, `step4_qc_overview.png` and
`step4_detection_overlay_examples.png` since it existed. A headless run left
that directory with no figures at all, because the export was a method on the
dialog — the same shape as the Step 8 comparison figure and the Step 12
isochrone plots, and the same fix: the drawing moves next to the calculation and
the window keeps the widgets.

Everything needed is already on disk. Step 4 writes one CSV and one metadata
JSON per frame either way, so the frame table can be rebuilt without a window
having assembled it.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from astropy.io import fits

from apex.analysis.frame_qc import (
    FAIL, PASS, REVIEW, FrameQCThresholds, evaluate_frame_qc, summarize_frame_qc,
)
from apex.utils.step_paths import step2_cropped_dir, step4_dir


def _detect_source_table_path(params, fname: str):
    """Where Step 4 left this frame's source table — cache first, then output."""
    cache_dir = Path(getattr(params.P, "cache_dir", "") or "")
    for path in (cache_dir / f"detect_{fname}.csv",
                 Path(step4_dir(params.P.result_dir)) / f"detect_{fname}.csv"):
        if path.exists():
            return path
    return None


def _resolve_fits_path(params, fname: str, use_cropped: bool):
    if use_cropped:
        candidate = Path(step2_cropped_dir(params.P.result_dir)) / fname
        if candidate.exists():
            return candidate
    try:
        return Path(params.get_file_path(fname))
    except Exception:                               # noqa: BLE001
        return None


def build_frame_table(result_dir, params=None) -> pd.DataFrame:
    """Rebuild the per-frame QC table from what Step 4 left on disk.

    The window assembles this while it scans; a batch run has no window, so read
    the metadata JSONs back. Fields absent from an older run come back as NaN
    rather than raising — a figure from a partial run is better than no figure.
    """
    out_dir = Path(step4_dir(result_dir))
    rows = []
    for meta_path in sorted(out_dir.glob("detect_*.json")):
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except Exception:                       # noqa: BLE001
            continue
        name = meta_path.name[len("detect_"):-len(".json")]
        rows.append({
            "file": name,
            "filter": meta.get("filter", ""),
            "sky_med": _safe_float(meta.get("bkg_median")),
            "sky_sigma": _safe_float(meta.get("bkg_rms")),
            "sky_e": _safe_float(meta.get("sky_e")),
            "sky_sigma_e": _safe_float(meta.get("sky_sigma_e")),
            "sky_sigma_expected_e": _safe_float(meta.get("sky_sigma_expected_e")),
            "fwhm_med": _safe_float(meta.get("fwhm_px")),
            "n_sources": int(meta.get("n_sources", 0) or 0),
            "n_raw_detections": int(meta.get("n_raw_detections", 0) or 0),
            "n_after_shape_filter": int(meta.get("n_after_shape_filter", 0) or 0),
            "detect_capped": bool(meta.get("detect_capped", False)),
            "elong_med": _safe_float(meta.get("median_elongation")),
            "round_med": _safe_float(meta.get("median_roundness")),
            "airmass": _safe_float(meta.get("airmass")),
            "time_val": _safe_float(meta.get("time_val")),
            "qc_status": str(meta.get("qc_status", "") or ""),
            "quality_score_median": _safe_float(meta.get("quality_score_median")),
        })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["time_index"] = np.arange(len(df), dtype=float)
    df = _join_header_scan(result_dir, df)
    if params is not None:
        # Without this every frame carries an empty verdict, every mask in the
        # overview is empty, and the figure renders blank under a title reading
        # "PASS=0 REVIEW=0 FAIL=0". Measured on the first headless render.
        df = evaluate_frame_qc(df, params.P, auto_qc_thresholds(params))
    return df


def _join_header_scan(result_dir, df: pd.DataFrame) -> pd.DataFrame:
    """Bring airmass and time across from Step 1's `headers.csv`.

    Without it the QC overview plots against a bare frame index while the window
    plots against airmass — the same numbers on a different axis, which is the
    kind of difference that makes two figures of the same run look unrelated.
    """
    from apex.utils.step_paths import step1_dir

    scan = Path(step1_dir(result_dir)) / "headers.csv"
    if not scan.exists():
        return df
    try:
        head = pd.read_csv(scan)
    except Exception:                               # noqa: BLE001
        return df
    if "Filename" not in head.columns:
        return df

    head = head.rename(columns={"Filename": "file"})
    keep = ["file"]
    if "AIRMASS" in head.columns:
        head["airmass_scan"] = pd.to_numeric(head["AIRMASS"], errors="coerce")
        keep.append("airmass_scan")
    if "JD" in head.columns:
        head["time_scan"] = pd.to_numeric(head["JD"], errors="coerce")
        keep.append("time_scan")

    merged = df.merge(head[keep], on="file", how="left")
    if "airmass_scan" in merged.columns:
        merged["airmass"] = merged["airmass"].fillna(merged.pop("airmass_scan"))
    if "time_scan" in merged.columns:
        merged["time_val"] = merged["time_val"].fillna(merged.pop("time_scan"))
    return merged


def auto_qc_thresholds(params) -> FrameQCThresholds:
    """The window's thresholds, minus the spin boxes it reads them from.

    Only one is configurable outside the dialog — `fwhm_elong_max` — so the rest
    fall back to the same dataclass defaults the dialog starts at.
    """
    base = FrameQCThresholds()
    elong_fail = max(1.0, _safe_float(
        getattr(params.P, "fwhm_elong_max", base.elong_fail), base.elong_fail))
    return FrameQCThresholds(
        fwhm_z_review=base.fwhm_z_review,
        fwhm_z_fail=base.fwhm_z_review + 1.5,
        fwhm_model_ratio_review=base.fwhm_model_ratio_review,
        fwhm_model_ratio_fail=base.fwhm_model_ratio_fail,
        sky_z_review=base.sky_z_review,
        sky_z_fail=base.sky_z_review + 1.5,
        nsrc_z_review=base.nsrc_z_review,
        nsrc_z_fail=base.nsrc_z_review + 1.5,
        elong_review=min(base.elong_review, max(1.0, elong_fail - 0.08)),
        elong_fail=elong_fail,
    )


def export_qc_products(result_dir, df=None, *, exclude_reasons=None,
                       params=None) -> list[Path]:
    """Write the three products. Returns what was written.

    `df` is the window's table when a window has one, and rebuilt from disk when
    it does not — so both routes draw the same figure from the same numbers.
    """
    out_dir = Path(step4_dir(result_dir))
    out_dir.mkdir(parents=True, exist_ok=True)
    if df is None or getattr(df, "empty", True):
        df = build_frame_table(result_dir, params)
    if df.empty:
        return []

    saved = [write_qc_summary_csv(df, out_dir)]

    fig = Figure(figsize=(10.5, 7.2), dpi=120)
    draw_qc_overview(fig, df, exclude_reasons=exclude_reasons, params=params)
    path = out_dir / "step4_qc_overview.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    saved.append(path)

    examples = select_overlay_example_rows(df)
    if examples:
        fig = Figure(figsize=(11.0, 5.8), dpi=120)
        if draw_detection_overlay_examples(fig, examples, params=params):
            path = out_dir / "step4_detection_overlay_examples.png"
            fig.savefig(path, dpi=160, bbox_inches="tight")
            saved.append(path)
    return [p for p in saved if p is not None]


def _safe_float(value, default=np.nan):
    try:
        if value is None:
            return default
        return float(value)
    except Exception:
        return default


def _qc_plot_x(df: pd.DataFrame) -> tuple[np.ndarray, str]:
    airmass = pd.to_numeric(df.get("airmass", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
    time_vals = pd.to_numeric(df.get("time_val", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
    if int(np.isfinite(airmass).sum()) >= 3:
        return airmass, "Airmass"
    if int(np.isfinite(time_vals).sum()) >= 3:
        return time_vals, "Time"
    return pd.to_numeric(df.get("time_index", pd.Series(np.arange(len(df)), index=df.index)), errors="coerce").to_numpy(float), "Index"


def draw_qc_overview(fig, df, exclude_reasons=None, params=None) -> None:
    exclude_reasons = exclude_reasons or {}
    fig.clear()
    ax_sky = fig.add_subplot(2, 2, 1)
    ax_fwhm = fig.add_subplot(2, 2, 2)
    ax_nsrc = fig.add_subplot(2, 2, 3)
    ax_elong = fig.add_subplot(2, 2, 4)

    x_vals, x_label = _qc_plot_x(df)
    files = df["file"].astype(str).tolist()
    excluded = np.array([len(exclude_reasons.get(f, set())) > 0 for f in files])
    status = df.get("qc_status", pd.Series([""] * len(df))).fillna("").astype(str).str.upper().to_numpy()
    masks = [
        ((status == PASS) & ~excluded, "#212121", "o", "PASS"),
        ((status == REVIEW) & ~excluded, "#F9A825", "^", "REVIEW"),
        ((status == FAIL) & ~excluded, "#D32F2F", "s", "FAIL"),
        (excluded, "#9E9E9E", "x", "excluded"),
    ]

    def _scatter(ax, x, y):
        handles = []
        for mask, color, marker, label in masks:
            finite = mask & np.isfinite(x) & np.isfinite(y)
            if not np.any(finite):
                continue
            handles.append(
                ax.scatter(x[finite], y[finite], s=22, c=color, marker=marker, alpha=0.9, label=label)
            )
        return handles

    sky_x = pd.to_numeric(df.get("sky_e", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
    sky_y = pd.to_numeric(df.get("sky_sigma_e", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
    sky_x_label = "sky_e"
    sky_y_label = "sky_sigma_e"
    if int(np.isfinite(sky_x).sum()) < 3 or int(np.isfinite(sky_y).sum()) < 3:
        sky_x = pd.to_numeric(df.get("sky_med", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
        sky_y = pd.to_numeric(df.get("sky_sigma", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
        sky_x_label = "sky_med"
        sky_y_label = "sky_sigma"
    _scatter(ax_sky, sky_x, sky_y)
    expected = pd.to_numeric(df.get("sky_sigma_expected_e", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
    model_ok = sky_x_label == "sky_e" and np.isfinite(sky_x) & np.isfinite(expected) & (expected > 0)
    if np.any(model_ok):
        order = np.argsort(sky_x[model_ok])
        ax_sky.plot(sky_x[model_ok][order], expected[model_ok][order], color="#1565C0", lw=1.2, label="sqrt(sky + RN^2)")
    ax_sky.set_title("Sky Noise")
    ax_sky.set_xlabel(sky_x_label)
    ax_sky.set_ylabel(sky_y_label)

    fwhm_y = pd.to_numeric(df.get("fwhm_med", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
    _scatter(ax_fwhm, x_vals, fwhm_y)
    fwhm_model = pd.to_numeric(df.get("fwhm_model_px", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
    ok = np.isfinite(x_vals) & np.isfinite(fwhm_model) & (fwhm_model > 0)
    if np.any(ok):
        order = np.argsort(x_vals[ok])
        label = "FWHM0 * X^(3/5)" if x_label == "Airmass" else "robust FWHM model"
        ax_fwhm.plot(x_vals[ok][order], fwhm_model[ok][order], color="#1565C0", lw=1.2, label=label)
    fwhm_cut = pd.to_numeric(df.get("fwhm_high_cut_px", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
    ok = np.isfinite(x_vals) & np.isfinite(fwhm_cut)
    if np.any(ok):
        order = np.argsort(x_vals[ok])
        ax_fwhm.plot(x_vals[ok][order], fwhm_cut[ok][order], color="#E53935", lw=1.0, ls="--", label="review cut")
    ax_fwhm.set_title("Seeing / FWHM")
    ax_fwhm.set_xlabel(x_label)
    ax_fwhm.set_ylabel("fwhm_px")

    nsrc_y = pd.to_numeric(df.get("n_sources", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
    _scatter(ax_nsrc, x_vals, nsrc_y)
    nsrc_trend = pd.to_numeric(df.get("n_sources_trend", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
    ok = np.isfinite(x_vals) & np.isfinite(nsrc_trend)
    if np.any(ok):
        order = np.argsort(x_vals[ok])
        ax_nsrc.plot(x_vals[ok][order], nsrc_trend[ok][order], color="#1565C0", lw=1.2, label="robust median")
    nsrc_low = pd.to_numeric(df.get("n_sources_low_cut", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
    ok = np.isfinite(x_vals) & np.isfinite(nsrc_low)
    if np.any(ok):
        order = np.argsort(x_vals[ok])
        ax_nsrc.plot(x_vals[ok][order], nsrc_low[ok][order], color="#E53935", lw=1.0, ls="--", label="review cut")
    ax_nsrc.set_title("Detected Sources")
    ax_nsrc.set_xlabel(x_label)
    ax_nsrc.set_ylabel("n_sources")

    elong_y = pd.to_numeric(df.get("elong_med", pd.Series(np.nan, index=df.index)), errors="coerce").to_numpy(float)
    _scatter(ax_elong, x_vals, elong_y)
    elong_cut = _safe_float(getattr(params.P, "fwhm_elong_max", 1.3), 1.3)
    if np.isfinite(elong_cut) and elong_cut > 0:
        ax_elong.axhline(elong_cut, color="#E53935", lw=1.0, ls="--", label=f"elong cut={elong_cut:.2f}")
    ax_elong.set_title("Shape")
    ax_elong.set_xlabel(x_label)
    ax_elong.set_ylabel("median elongation")

    counts = summarize_frame_qc(df)
    fig.suptitle(
        f"Step 4 Auto QC | PASS={counts.get(PASS, 0)} REVIEW={counts.get(REVIEW, 0)} FAIL={counts.get(FAIL, 0)}",
        fontsize=12,
    )
    for ax in (ax_sky, ax_fwhm, ax_nsrc, ax_elong):
        ax.grid(True, alpha=0.2)
        if ax.get_legend_handles_labels()[0]:
            ax.legend(loc="best", fontsize=7, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.96))


def write_qc_summary_csv(df, out_dir) -> Path:
    rows = []
    df_work = df.copy()
    df_work["filter"] = df_work.get("filter", "").fillna("").astype(str)
    groups = [("ALL", df_work)]
    groups.extend((str(filt) or "(none)", grp) for filt, grp in df_work.groupby("filter", sort=True))
    for label, grp in groups:
        if grp.empty:
            continue
        passed = grp.get("passed", pd.Series([True] * len(grp), index=grp.index)).astype(bool)
        counts = summarize_frame_qc(grp)
        rows.append({
            "filter": label,
            "n_frames": int(len(grp)),
            "n_passed_pipeline": int(passed.sum()),
            "n_excluded_pipeline": int((~passed).sum()),
            "qc_pass": counts.get(PASS, 0),
            "qc_review": counts.get(REVIEW, 0),
            "qc_fail": counts.get(FAIL, 0),
            "median_fwhm_px": float(pd.to_numeric(grp.get("fwhm_med"), errors="coerce").median()),
            "median_sky": float(pd.to_numeric(grp.get("sky_med"), errors="coerce").median()),
            "median_n_sources": float(pd.to_numeric(grp.get("n_sources"), errors="coerce").median()),
            "median_elongation": float(pd.to_numeric(grp.get("elong_med"), errors="coerce").median()),
        })
    out_path = out_dir / "frame_quality_summary.csv"
    pd.DataFrame(rows).to_csv(out_path, index=False)
    return out_path


def _load_detect_sources_for_overlay(params, fname: str, max_sources: int = 2500) -> pd.DataFrame:
    path = _detect_source_table_path(params, fname)
    if path is None:
        return pd.DataFrame()
    try:
        df = pd.read_csv(path)
    except Exception:
        return pd.DataFrame()
    if df.empty or "x" not in df.columns or "y" not in df.columns:
        return pd.DataFrame()
    df = df.copy()
    df["x"] = pd.to_numeric(df["x"], errors="coerce")
    df["y"] = pd.to_numeric(df["y"], errors="coerce")
    df = df[df["x"].notna() & df["y"].notna()]
    if len(df) > max_sources:
        if "quality_score" in df.columns:
            df["_sort_quality"] = pd.to_numeric(df["quality_score"], errors="coerce").fillna(-np.inf)
            df = df.sort_values("_sort_quality", ascending=False).head(max_sources)
            df = df.drop(columns=["_sort_quality"])
        else:
            df = df.head(max_sources)
    return df


def select_overlay_example_rows(df) -> list[tuple[str, pd.Series]]:
    if df.empty or "file" not in df.columns:
        return []
    work = df.copy()
    work["qc_score"] = pd.to_numeric(work.get("qc_score", pd.Series(np.nan, index=work.index)), errors="coerce")
    work["qc_status"] = work.get("qc_status", pd.Series("", index=work.index)).fillna("").astype(str).str.upper()
    work["passed"] = work.get("passed", pd.Series([True] * len(work), index=work.index)).astype(bool)

    examples: list[tuple[str, pd.Series]] = []
    pass_df = work[(work["qc_status"] == PASS) & work["passed"]]
    if not pass_df.empty:
        examples.append(("Best PASS", pass_df.sort_values("qc_score", ascending=False).iloc[0]))

    fail_df = work[work["qc_status"] == FAIL]
    if fail_df.empty:
        fail_df = work[work["qc_status"] == REVIEW]
    if fail_df.empty:
        fail_df = work[~work["passed"]]
    if not fail_df.empty:
        worst = fail_df.sort_values("qc_score", ascending=True).iloc[0]
        if not examples or str(worst.get("file")) != str(examples[0][1].get("file")):
            examples.append(("Worst QC", worst))

    if not examples and not work.empty:
        examples.append(("Example", work.sort_values("qc_score", ascending=False).iloc[0]))
    return examples[:2]


def _read_overlay_image(params, fname: str, use_cropped: bool = False):
    path = _resolve_fits_path(params, fname, use_cropped)
    if path is None or not path.exists():
        return None
    try:
        with fits.open(path, memmap=True) as hdul:
            data = np.asarray(hdul[0].data, dtype=float)
    except Exception:
        return None
    if data.ndim > 2:
        data = np.squeeze(data)
    if data.ndim != 2:
        return None
    return data


def draw_detection_overlay_examples(fig, examples, params=None,
                                    use_cropped: bool = False) -> bool:
    fig.clear()
    drawable: list[tuple[str, pd.Series, np.ndarray, pd.DataFrame]] = []
    for label, row in examples:
        fname = str(row.get("file", "") or "")
        if not fname:
            continue
        image = _read_overlay_image(params, fname, use_cropped)
        sources = _load_detect_sources_for_overlay(params, fname)
        if image is None or image.size == 0 or sources.empty:
            continue
        drawable.append((label, row, image, sources))
    if not drawable:
        return False

    axes = fig.subplots(1, len(drawable), squeeze=False)[0]
    for ax, (label, row, image, sources) in zip(axes, drawable):
        h, w = image.shape
        stride = max(1, int(np.ceil(max(h, w) / 1800.0)))
        disp = image[::stride, ::stride]
        finite = disp[np.isfinite(disp)]
        if finite.size:
            vmin, vmax = np.nanpercentile(finite, [1.0, 99.5])
            if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
                vmin, vmax = np.nanmedian(finite), np.nanmax(finite)
        else:
            vmin, vmax = 0.0, 1.0
        ax.imshow(disp, origin="lower", cmap="gray", vmin=vmin, vmax=vmax)
        ax.scatter(
            sources["x"].to_numpy(float) / stride,
            sources["y"].to_numpy(float) / stride,
            s=5,
            facecolors="none",
            edgecolors="#00E5FF",
            linewidths=0.35,
            alpha=0.75,
        )
        fname = str(row.get("file", "") or "")
        status = str(row.get("qc_status", "") or "")
        reasons = str(row.get("qc_reasons", "") or "").strip()
        score = _safe_float(row.get("qc_score"), np.nan)
        title = f"{label}: {status} score={score:.1f}\n{fname}"
        ax.set_title(title, fontsize=9)
        ax.set_xlim(0, w / stride)
        ax.set_ylim(0, h / stride)
        ax.set_xlabel("x / downsample")
        ax.set_ylabel("y / downsample")
        if reasons:
            ax.text(
                0.01,
                0.01,
                reasons,
                transform=ax.transAxes,
                fontsize=7,
                color="white",
                va="bottom",
                ha="left",
                bbox={"facecolor": "black", "alpha": 0.55, "edgecolor": "none", "pad": 3},
            )
    fig.suptitle("Step 4 Detection Overlay Examples", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return True
