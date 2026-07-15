"""Figure 16: repeatability of legacy PSF and aperture photometry on NGC 6811.

The comparison uses only same-frame sources matched within one pixel. Per-frame
median offsets are removed independently for each method before calculating
the standard deviation of each star across repeated frames.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "validation" / "paper"))

from apex_paper_style import C, DOUBLE_COL, PALETTE, apply_paper_style, save_fig  # noqa: E402


DEFAULT_RESULT_DIR = Path(r"E:\observed_Analysis\NGC6811\pp\result")
OUTDIR = REPO / "validation" / "paper" / "figures"
DATADIR = REPO / "validation" / "paper" / "data" / "ngc6811_psf_aperture"


def _numeric(frame: pd.DataFrame, columns: tuple[str, ...]) -> None:
    for column in columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")


def _true_mask(values: pd.Series) -> np.ndarray:
    return values.astype(str).str.lower().isin(("true", "1", "yes")).to_numpy()


def _load_matches(result_dir: Path) -> pd.DataFrame:
    aperture_dir = result_dir / "step7_forced_phot"
    psf_dir = result_dir / "cmd_psf"
    rows: list[pd.DataFrame] = []
    for psf_path in sorted(psf_dir.glob("photometry_*.tsv")):
        aperture_path = aperture_dir / psf_path.name
        if not aperture_path.exists():
            continue
        psf = pd.read_csv(psf_path, sep="\t")
        aperture = pd.read_csv(aperture_path, sep="\t", low_memory=False)
        _numeric(psf, ("x_fit", "y_fit", "mag_psf", "snr_psf", "flags_psf"))
        _numeric(aperture, ("x_fit", "y_fit", "mag_inst", "snr", "master_id"))

        psf_good = (
            np.isfinite(psf["x_fit"])
            & np.isfinite(psf["y_fit"])
            & np.isfinite(psf["mag_psf"])
            & np.isfinite(psf["snr_psf"])
            & (psf["flags_psf"].fillna(0) == 0)
        )
        aperture_good = (
            np.isfinite(aperture["x_fit"])
            & np.isfinite(aperture["y_fit"])
            & np.isfinite(aperture["mag_inst"])
            & np.isfinite(aperture["snr"])
            & np.isfinite(aperture["master_id"])
        )
        if "bad_phot_flag" in aperture:
            aperture_good &= ~_true_mask(aperture["bad_phot_flag"])
        psf = psf.loc[psf_good].copy()
        aperture = aperture.loc[aperture_good].copy()
        if psf.empty or aperture.empty:
            continue

        distance, index = cKDTree(aperture[["x_fit", "y_fit"]]).query(
            psf[["x_fit", "y_fit"]], k=1, workers=1
        )
        keep = distance <= 1.0
        psf = psf.loc[keep].reset_index(drop=True)
        matched_aperture = aperture.iloc[index[keep]].reset_index(drop=True)
        matched = pd.DataFrame(
            {
                "file": psf_path.name.removeprefix("photometry_").removesuffix(".tsv"),
                "filter": psf["FILTER"].astype(str).str.upper(),
                "master_id": matched_aperture["master_id"].astype("int64"),
                "mag_psf": psf["mag_psf"],
                "mag_aperture": matched_aperture["mag_inst"],
                "snr_psf": psf["snr_psf"],
                "snr_aperture": matched_aperture["snr"],
                "match_distance_px": distance[keep],
            }
        )
        rows.append(
            matched.sort_values("match_distance_px").drop_duplicates("master_id")
        )
    if not rows:
        raise RuntimeError(f"No matched Step 7/Step 8 tables under {result_dir}")
    return pd.concat(rows, ignore_index=True)


def _remove_frame_offsets(frame: pd.DataFrame, column: str) -> pd.DataFrame:
    work = frame.copy()
    star_level = work.groupby("master_id")[column].median()
    for _ in range(8):
        frame_offset = (
            work[column] - work["master_id"].map(star_level)
        ).groupby(work["file"]).median()
        corrected = work[column] - work["file"].map(frame_offset)
        star_level = corrected.groupby(work["master_id"]).median()
    work[f"{column}_corrected"] = work[column] - work["file"].map(frame_offset)
    work[f"{column}_residual"] = (
        work[f"{column}_corrected"] - work["master_id"].map(star_level)
    )
    return work


def _repeatability(matches: pd.DataFrame) -> pd.DataFrame:
    tables = []
    for filter_name, subset in matches.groupby("filter"):
        subset = subset[
            (subset["snr_psf"] >= 20.0) & (subset["snr_aperture"] >= 20.0)
        ].copy()
        psf = _remove_frame_offsets(subset, "mag_psf")
        aperture = _remove_frame_offsets(subset, "mag_aperture")
        merged = psf[
            ["file", "master_id", "mag_psf_corrected", "mag_psf_residual"]
        ].merge(
            aperture[
                [
                    "file",
                    "master_id",
                    "mag_aperture_corrected",
                    "mag_aperture_residual",
                ]
            ],
            on=["file", "master_id"],
        )
        minimum_frames = max(3, int(np.ceil(subset["file"].nunique() / 2)))
        counts = merged.groupby("master_id").size()
        merged = merged[
            merged["master_id"].isin(counts[counts >= minimum_frames].index)
        ]
        stars = merged.groupby("master_id").agg(
            n_frames=("file", "size"),
            magnitude=("mag_aperture_corrected", "median"),
            scatter_psf=("mag_psf_residual", "std"),
            scatter_aperture=("mag_aperture_residual", "std"),
        )
        stars = stars.replace([np.inf, -np.inf], np.nan).dropna().reset_index()
        stars.insert(0, "filter", filter_name)
        tables.append(stars)
    return pd.concat(tables, ignore_index=True)


def _residual_yield(result_dir: Path) -> pd.DataFrame:
    rows = []
    for path in sorted((result_dir / "cmd_psf").glob("residual_meta_*.json")):
        metadata = json.loads(path.read_text(encoding="utf-8"))
        iterations = metadata.get("iters") or []
        if not iterations:
            continue
        initial = float(iterations[0].get("n_fit", np.nan))
        added = float(metadata.get("n_new_raw", iterations[-1].get("n_new_kept", np.nan)))
        rows.append(
            {
                "file": metadata.get("file", path.stem),
                "filter": str(metadata.get("filter", "?")).upper(),
                "initial_sources": initial,
                "residual_sources": added,
                "residual_fraction_pct": 100.0 * added / initial if initial > 0 else np.nan,
            }
        )
    return pd.DataFrame(rows)


def _binned_curve(frame: pd.DataFrame, column: str) -> tuple[np.ndarray, np.ndarray]:
    quantiles = np.linspace(0.0, 1.0, 7)
    edges = np.unique(frame["magnitude"].quantile(quantiles).to_numpy(float))
    x_values, y_values = [], []
    for low, high in zip(edges[:-1], edges[1:]):
        keep = (frame["magnitude"] >= low) & (frame["magnitude"] <= high)
        if np.sum(keep) < 10:
            continue
        x_values.append(float(frame.loc[keep, "magnitude"].median()))
        y_values.append(float(frame.loc[keep, column].median()))
    return np.asarray(x_values), np.asarray(y_values)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--result-dir", type=Path, default=DEFAULT_RESULT_DIR)
    args = parser.parse_args()

    apply_paper_style()
    matches = _load_matches(args.result_dir)
    repeatability = _repeatability(matches)
    residual_yield = _residual_yield(args.result_dir)
    DATADIR.mkdir(parents=True, exist_ok=True)
    matches.to_csv(DATADIR / "matched_measurements.csv", index=False)
    repeatability.to_csv(DATADIR / "repeatability.csv", index=False)
    residual_yield.to_csv(DATADIR / "residual_yield.csv", index=False)

    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL, 5.8))
    for ax, filter_name in zip(axes.ravel()[:3], ("B", "V", "R")):
        subset = repeatability[repeatability["filter"] == filter_name]
        for column, label, color, marker in (
            ("scatter_psf", "legacy PSF", C["data"], "o"),
            ("scatter_aperture", "aperture", C["accent"], "s"),
        ):
            ax.scatter(
                subset["magnitude"], subset[column], s=5, alpha=0.08,
                color=color, edgecolors="none"
            )
            x_values, y_values = _binned_curve(subset, column)
            ax.plot(x_values, y_values, marker=marker, color=color, label=label)
        med_psf = float(subset["scatter_psf"].median())
        med_aperture = float(subset["scatter_aperture"].median())
        ax.set_yscale("log")
        ax.set_ylim(0.003, 0.2)
        ax.set_xlabel("median instrumental magnitude")
        ax.set_ylabel("repeatability scatter (mag)")
        ax.set_title(
            f"({chr(97 + list(('B', 'V', 'R')).index(filter_name))}) {filter_name}: "
            f"{med_psf:.3f} vs {med_aperture:.3f} mag",
            loc="left",
        )
        ax.legend()

    ax_yield = axes.ravel()[3]
    rng = np.random.default_rng(20260714)
    for index, filter_name in enumerate(("B", "V", "R")):
        values = residual_yield.loc[
            residual_yield["filter"] == filter_name, "residual_fraction_pct"
        ].dropna().to_numpy(float)
        ax_yield.scatter(
            index + rng.uniform(-0.08, 0.08, len(values)), values,
            s=22, color=PALETTE["green"], alpha=0.75
        )
        if len(values):
            ax_yield.plot(
                [index - 0.18, index + 0.18], [np.median(values)] * 2,
                color="black", lw=1.2
            )
    ax_yield.set_xticks(range(3), ("B", "V", "R"))
    ax_yield.set_xlabel("filter")
    ax_yield.set_ylabel("new residual detections / initial (%)")
    ax_yield.set_title("(d) residual-source yield", loc="left")

    fig.suptitle("NGC 6811: legacy PSF versus aperture repeatability", y=1.01)
    fig.tight_layout()
    paths = save_fig(fig, "fig16_ngc6811_psf_aperture", OUTDIR)
    plt.close(fig)
    print("saved:", ", ".join(str(path) for path in paths.values()))
    for filter_name, subset in repeatability.groupby("filter"):
        print(
            filter_name,
            "N=", len(subset),
            "PSF=", float(subset["scatter_psf"].median()),
            "aperture=", float(subset["scatter_aperture"].median()),
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
