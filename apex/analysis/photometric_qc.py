"""Post-photometry (second-stage) frame QC from matched-star magnitudes.

Step-4 frame QC sees only image-level statistics (FWHM, sky, N_sources), which
are blind to *transparency* loss: thin cirrus dims every star by tenths of a
magnitude while leaving FWHM, sky and the source count nearly unchanged (see
validation/paper fig6, where injected cloud frames sail through image-level
QC). After Step 7 forced photometry the cure is nearly free: the same stars
are measured on every frame, so a per-frame magnitude offset against each
star's night median is a direct transparency measurement — the same frame
offset (z_j) the LC detrend and extinction fits solve for.

This module is Qt-free (analysis layer) and consumes Step 7 photometry tables.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd

PASS = "PASS"
REVIEW = "REVIEW"
FAIL = "FAIL"
SKIP = "SKIP"


@dataclass(frozen=True)
class PhotometricQCThresholds:
    """Conservative gates for the transparency / stability checks."""

    offset_review_mag: float = 0.3   # frame dimmer than night median by this -> REVIEW
    offset_fail_mag: float = 0.8     # -> FAIL (thick cloud)
    scatter_ratio_review: float = 3.0  # frame scatter vs night median scatter
    scatter_floor_mag: float = 0.02    # ignore scatter blowups below this absolute level
    min_snr: float = 20.0            # QC stars must be this bright
    min_frame_fraction: float = 0.7  # ... and present in this fraction of frames
    min_stars: int = 8               # below this the verdict is SKIP, never FAIL


def _star_frame_matrix(
    frames: Mapping[str, pd.DataFrame],
    filter_name: str,
    thr: PhotometricQCThresholds,
) -> pd.DataFrame:
    """Long table -> (star x frame) instrumental-mag matrix for one filter."""
    pieces: list[pd.DataFrame] = []
    for fname, df in frames.items():
        if df is None or df.empty:
            continue
        cols = df.columns
        if "mag_inst" not in cols or "source_id" not in cols:
            continue
        sub = pd.DataFrame(
            {
                "source_id": pd.to_numeric(df["source_id"], errors="coerce"),
                "mag_inst": pd.to_numeric(df["mag_inst"], errors="coerce"),
                "snr": pd.to_numeric(df.get("snr"), errors="coerce"),
            }
        )
        if "filter" in cols:
            filt = df["filter"].fillna("").astype(str)
            sub = sub[filt == filter_name]
        sub = sub[np.isfinite(sub["mag_inst"]) & np.isfinite(sub["source_id"])]
        if sub.empty:
            continue
        sub = sub.drop_duplicates(subset="source_id", keep="first")
        sub["file"] = fname
        pieces.append(sub)
    if not pieces:
        return pd.DataFrame()
    long = pd.concat(pieces, ignore_index=True)

    # Bright, well-covered QC stars only.
    per_star = long.groupby("source_id").agg(
        snr_med=("snr", "median"), n_frames=("file", "nunique")
    )
    n_frames_total = long["file"].nunique()
    good = per_star[
        (per_star["snr_med"] >= thr.min_snr)
        & (per_star["n_frames"] >= thr.min_frame_fraction * n_frames_total)
    ].index
    long = long[long["source_id"].isin(good)]
    if long.empty:
        return pd.DataFrame()
    return long.pivot_table(index="source_id", columns="file", values="mag_inst", aggfunc="first")


def evaluate_photometric_qc(
    frames: Mapping[str, pd.DataFrame],
    thresholds: PhotometricQCThresholds | None = None,
) -> pd.DataFrame:
    """Per-frame transparency offset + scatter with PASS/REVIEW/FAIL/SKIP.

    ``frames`` maps frame filename -> Step 7 photometry table (needs
    ``source_id``, ``mag_inst``; ``snr`` and ``filter`` are used when present).
    Returns one row per (frame, filter) with:

    - ``transparency_offset_mag`` — median over QC stars of (m_ij − star night
      median); positive = frame is DIMMER (cloud). One refinement pass excludes
      strongly offset frames from the star medians so a few cloudy frames do
      not bias the reference.
    - ``frame_scatter_mag`` — robust (MAD) scatter of those residuals.
    - ``phot_qc_status`` / ``phot_qc_reasons`` — conservative gates; frames in
      filters with too few QC stars are SKIP (never punished for sparse data).
    """
    thr = thresholds or PhotometricQCThresholds()

    filters: set[str] = set()
    for df in frames.values():
        if df is None or df.empty:
            continue
        if "filter" in df.columns:
            filters.update(df["filter"].fillna("").astype(str).unique())
        else:
            filters.add("")
    if not filters:
        filters = {""}

    rows: list[dict] = []
    for filt in sorted(filters):
        matrix = _star_frame_matrix(frames, filt, thr)
        frame_names = [
            f for f, df in frames.items()
            if df is not None and not df.empty
            and ("filter" not in df.columns or (df["filter"].fillna("").astype(str) == filt).any())
        ]
        if matrix.empty or len(matrix.index) < thr.min_stars:
            for fname in frame_names:
                rows.append(
                    {
                        "file": fname,
                        "filter": filt,
                        "n_qc_stars": int(len(matrix.index)) if not matrix.empty else 0,
                        "transparency_offset_mag": np.nan,
                        "frame_scatter_mag": np.nan,
                        "phot_qc_status": SKIP,
                        "phot_qc_reasons": "too_few_qc_stars",
                    }
                )
            continue

        values = matrix.to_numpy(float)  # (n_stars, n_frames)

        # Pass 1: star medians over all frames -> provisional frame offsets.
        star_med = np.nanmedian(values, axis=1, keepdims=True)
        offsets = np.nanmedian(values - star_med, axis=0)
        # Refinement: rebuild star medians from frames that look clear, so a
        # handful of cloudy frames cannot drag the reference faint.
        clear = np.abs(offsets - np.nanmedian(offsets)) <= thr.offset_review_mag
        if int(np.sum(clear)) >= max(3, values.shape[1] // 3):
            star_med = np.nanmedian(values[:, clear], axis=1, keepdims=True)
        resid = values - star_med
        offsets = np.nanmedian(resid, axis=0)
        # Center on the night so offsets read as "vs typical clear frame".
        offsets = offsets - np.nanmedian(offsets)
        scatter = 1.4826 * np.nanmedian(np.abs(resid - offsets[None, :]), axis=0)
        night_scatter = float(np.nanmedian(scatter))

        n_stars_per_frame = np.sum(np.isfinite(values), axis=0)
        for j, fname in enumerate(matrix.columns):
            off = float(offsets[j]) if np.isfinite(offsets[j]) else np.nan
            sc = float(scatter[j]) if np.isfinite(scatter[j]) else np.nan
            reasons: list[str] = []
            status = PASS
            if np.isfinite(off):
                if off > thr.offset_fail_mag:
                    status = FAIL
                    reasons.append("transparency_loss")
                elif off > thr.offset_review_mag:
                    status = REVIEW
                    reasons.append("transparency_warning")
                elif off < -thr.offset_review_mag:
                    status = REVIEW
                    reasons.append("negative_offset_anomaly")
            if (
                np.isfinite(sc)
                and np.isfinite(night_scatter)
                and night_scatter > 0
                and sc > max(thr.scatter_ratio_review * night_scatter, thr.scatter_floor_mag)
            ):
                if status == PASS:
                    status = REVIEW
                reasons.append("frame_scatter")
            rows.append(
                {
                    "file": str(fname),
                    "filter": filt,
                    "n_qc_stars": int(n_stars_per_frame[j]),
                    "transparency_offset_mag": off,
                    "frame_scatter_mag": sc,
                    "phot_qc_status": status,
                    "phot_qc_reasons": ",".join(reasons),
                }
            )
        # Frames of this filter that contributed no QC-star measurements.
        seen = {str(c) for c in matrix.columns}
        for fname in frame_names:
            if fname not in seen:
                rows.append(
                    {
                        "file": fname,
                        "filter": filt,
                        "n_qc_stars": 0,
                        "transparency_offset_mag": np.nan,
                        "frame_scatter_mag": np.nan,
                        "phot_qc_status": SKIP,
                        "phot_qc_reasons": "no_qc_star_measurements",
                    }
                )
    return pd.DataFrame(rows)


def summarize_photometric_qc(qc_df: pd.DataFrame) -> dict[str, int]:
    counts = {PASS: 0, REVIEW: 0, FAIL: 0, SKIP: 0}
    if qc_df is None or qc_df.empty or "phot_qc_status" not in qc_df.columns:
        return counts
    series = qc_df["phot_qc_status"].fillna("").astype(str).str.upper()
    for key in counts:
        counts[key] = int((series == key).sum())
    return counts


def run_photometric_qc(
    result_dir: Path | str,
    filenames: list[str] | None = None,
    thresholds: PhotometricQCThresholds | None = None,
    write_csv: bool = True,
) -> pd.DataFrame:
    """Load Step 7 tables for ``filenames`` (default: photometry_index.csv),
    evaluate, and (optionally) write ``phot_quality.csv`` next to them."""
    from apex.utils.photometry_loader import load_frame_photometry
    from apex.utils.step_paths import step7_forced_phot_dir

    result_dir = Path(result_dir)
    phot_dir = step7_forced_phot_dir(result_dir)
    if filenames is None:
        index_path = phot_dir / "photometry_index.csv"
        if not index_path.exists():
            return pd.DataFrame()
        idx = pd.read_csv(index_path)
        name_col = "file" if "file" in idx.columns else idx.columns[0]
        filenames = idx[name_col].dropna().astype(str).tolist()

    frames: dict[str, pd.DataFrame] = {}
    for fname in filenames:
        try:
            df = load_frame_photometry(result_dir, fname, sid_map={})
        except Exception:
            df = None
        if df is not None and not df.empty:
            frames[fname] = df

    qc = evaluate_photometric_qc(frames, thresholds)
    if write_csv and not qc.empty:
        out_path = phot_dir / "phot_quality.csv"
        try:
            qc.to_csv(out_path, index=False)
        except OSError:
            pass
    return qc
