"""Summarize the one-CPU M5 PSF seeing sweep."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from apex.analysis.psf_diagnostics import epsf_shape_metrics


def _robust_scatter(values: pd.Series | np.ndarray) -> float:
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if data.size == 0:
        return np.nan
    median = float(np.median(data))
    return float(1.4826 * np.median(np.abs(data - median)))


def _rmse(values: pd.Series | np.ndarray) -> float:
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    return float(np.sqrt(np.mean(np.square(data)))) if data.size else np.nan


def _frame_summary(result_dir: Path, index: pd.Series) -> dict[str, object]:
    filename = str(index["file"])
    psf_dir = result_dir / "cmd_psf"
    meta = json.loads(
        (psf_dir / f"residual_meta_{filename}.json").read_text(encoding="utf-8")
    )
    detect = json.loads(
        (result_dir / "cache" / f"detect_{filename}.json").read_text(
            encoding="utf-8"
        )
    )
    frame_stem = Path(filename).stem
    epsf_path = next(psf_dir.glob(f"epsf_model_*_{frame_stem}.fits"))
    shape = epsf_shape_metrics(np.asarray(fits.getdata(epsf_path), dtype=float))
    fit_elapsed = sum(float(item.get("elapsed_s", 0.0)) for item in meta["iters"])
    fit_failure_fraction = float(index.get(
        "psf_fit_failure_fraction",
        float(index["n_fail"]) / max(float(index["n"]), 1.0),
    ))

    summary: dict[str, object] = {
        "file": filename,
        "fwhm_px": float(detect["fwhm_px"]),
        "fwhm_arcsec": float(detect["fwhm_arcsec"]),
        "step4_detected": int(detect["n_sources"]),
        "psf_fit": int(index["n"]),
        "psf_good": int(index["n_goodmag"]),
        "psf_fail": int(index["n_fail"]),
        "psf_clean_fraction": float(index.get(
            "psf_clean_fraction", 1.0 - fit_failure_fraction
        )),
        "psf_fit_failure_fraction": fit_failure_fraction,
        "psf_qc_status": str(index.get("psf_qc_status", "")),
        "psf_qc_score": float(index.get("psf_qc_score", np.nan)),
        "psf_qc_reasons": str(index.get("psf_qc_reasons", "")),
        "forced": int(index["n_forced"]),
        "crowding_unreliable": int(index["n_crowding_unreliable"]),
        "median_qfit": float(index["median_qfit"]),
        "median_qfit_noise_ratio": float(index.get("median_qfit_noise_ratio", np.nan)),
        "median_reduced_chi2": float(index["median_reduced_chi2"]),
        "epsf_candidates": int(index["epsf_n_candidates"]),
        "epsf_target": int(index["epsf_target"]),
        "epsf_selected": int(index["epsf_n_selected"]),
        "epsf_fallback": int(index["epsf_n_fallback_selected"]),
        "epsf_median_contamination": float(index["epsf_median_contamination"]),
        "epsf_ellipticity": float(shape["ellipticity"]),
        "epsf_rotation_asymmetry": float(shape["rotation_asymmetry"]),
        "fit_window_px": int(index.get("fit_window_px", 0)),
        "fit_window_energy": float(index.get("fit_window_energy", np.nan)),
        "psf_nea_px": float(index.get("psf_nea_px", np.nan)),
        "residual_std_final": float(meta["iters"][-1].get("residual_std", np.nan)),
        "fit_elapsed_s": fit_elapsed,
        "core_center_x": float(meta["core_cut"]["center_x"]),
        "core_center_y": float(meta["core_cut"]["center_y"]),
        "core_radius_px": float(meta["core_cut"]["radius_px"]),
        "core_method": str(meta["core_cut"]["method"]),
    }
    return summary


def _plot(summary: pd.DataFrame, output: Path) -> None:
    summary = summary.sort_values("fwhm_arcsec").reset_index(drop=True)
    x = summary["fwhm_arcsec"].to_numpy(dtype=float)
    labels = [f"{value:.2f}" for value in x]
    colors = ("#0072B2", "#E69F00", "#D55E00")
    positions = np.arange(len(summary))

    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.6), constrained_layout=True)

    axes[0, 0].plot(x, summary["step4_detected"], "o-", label="Step4 detected")
    axes[0, 0].plot(x, summary["psf_good"], "^-", label="clean PSF flux")
    axes[0, 0].set_xlabel("FWHM (arcsec)")
    axes[0, 0].set_ylabel("sources")
    axes[0, 0].set_title("(a) usable source yield", loc="left")
    axes[0, 0].legend(fontsize=8)
    for _, row in summary.iterrows():
        axes[0, 0].annotate(
            f"clean={100 * row['psf_clean_fraction']:.1f}%",
            (row["fwhm_arcsec"], row["psf_good"]),
            xytext=(0, 8), textcoords="offset points", ha="center", fontsize=7,
        )

    width = 0.36
    axes[0, 1].bar(
        positions - width / 2,
        summary["epsf_target"],
        width,
        color="0.78",
        label="candidate target",
    )
    axes[0, 1].bar(
        positions + width / 2,
        summary["epsf_selected"],
        width,
        color=colors,
        label="selected",
    )
    for pos, row in summary.iterrows():
        axes[0, 1].text(
            pos + width / 2,
            row["epsf_selected"] + 1.5,
            f"fallback={int(row['epsf_fallback'])}\nC={row['epsf_median_contamination']:.2f}",
            ha="center",
            va="bottom",
            fontsize=7,
        )
    axes[0, 1].set_ylim(0, max(20.0, 1.28 * float(summary["epsf_target"].max())))
    axes[0, 1].set_xticks(positions, labels)
    axes[0, 1].set_xlabel("FWHM (arcsec)")
    axes[0, 1].set_ylabel("ePSF reference stars")
    axes[0, 1].set_title("(b) ePSF reference quality", loc="left")
    axes[0, 1].legend(fontsize=8)

    axes[0, 2].plot(x, summary["epsf_ellipticity"], "o-", label="ellipticity")
    axes[0, 2].plot(
        x,
        summary["epsf_rotation_asymmetry"],
        "s-",
        label="180-degree asymmetry",
    )
    axes[0, 2].set_xlabel("FWHM (arcsec)")
    axes[0, 2].set_ylabel("shape metric")
    axes[0, 2].set_title("(c) ePSF shape", loc="left")
    axes[0, 2].legend(fontsize=8)

    axes[1, 0].plot(x, summary["median_qfit_noise_ratio"], "o-", label="qfit / noise")
    axes[1, 0].plot(x, summary["median_reduced_chi2"], "s-", label="reduced chi2")
    axes[1, 0].axhline(1.0, color="0.35", linewidth=0.8, linestyle=":")
    axes[1, 0].set_xlabel("FWHM (arcsec)")
    axes[1, 0].set_ylabel("dimensionless fit metric")
    axes[1, 0].set_title("(d) PSF residual consistency", loc="left")
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].plot(x, summary["fit_window_px"], "o-", color="#0072B2", label="fit window")
    axes[1, 1].set_xlabel("FWHM (arcsec)")
    axes[1, 1].set_ylabel("fit window (px)", color="#0072B2")
    axes[1, 1].tick_params(axis="y", labelcolor="#0072B2")
    nea_axis = axes[1, 1].twinx()
    nea_axis.plot(x, summary["psf_nea_px"], "s--", color="#D55E00", label="PSF NEA")
    nea_axis.set_ylabel("noise-equivalent area (px)", color="#D55E00")
    nea_axis.tick_params(axis="y", labelcolor="#D55E00")
    axes[1, 1].set_title("(e) adaptive fit footprint", loc="left")

    axes[1, 2].bar(positions, summary["fit_elapsed_s"], color=colors)
    for pos, row in summary.iterrows():
        axes[1, 2].text(
            pos,
            row["fit_elapsed_s"] + 2,
            f"{row['psf_qc_status']} ({row['psf_qc_score']:.0f})\n"
            f"res={row['residual_std_final']:.2f}",
            ha="center",
            va="bottom",
            fontsize=8,
        )
    axes[1, 2].set_ylim(0, 1.22 * float(summary["fit_elapsed_s"].max()))
    axes[1, 2].set_xticks(positions, labels)
    axes[1, 2].set_xlabel("FWHM (arcsec)")
    axes[1, 2].set_ylabel("iterative fit CPU time (s)")
    axes[1, 2].set_title("(f) one-CPU cost", loc="left")

    for axis in axes.flat:
        axis.grid(alpha=0.18)
    fig.suptitle(
        "M5 CPU-only PSF validation across seeing: same night, same i filter",
        fontsize=15,
    )
    fig.savefig(output, dpi=200, facecolor="white")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "result_dir",
        type=Path,
        nargs="?",
        default=Path("validation/real_gui_run/m5_seeing_sweep_v2/result"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation/real_gui_run/m5_seeing_sweep_v2"),
    )
    args = parser.parse_args()

    index = pd.read_csv(args.result_dir / "cmd_psf" / "photometry_index.csv")
    summary = pd.DataFrame(
        [_frame_summary(args.result_dir, row) for _, row in index.iterrows()]
    ).sort_values("fwhm_arcsec")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "m5_psf_seeing_summary.csv"
    figure_path = args.output_dir / "m5_psf_seeing_comparison.png"
    summary.to_csv(csv_path, index=False)
    _plot(summary, figure_path)
    print(summary.to_string(index=False))
    print(figure_path.resolve())
    print(csv_path.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
