"""Run the CPU-only PSF repeatability benchmark and write a four-panel figure."""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from apex.benchmark.psf_repeatability import analyze_psf_repeatability


def _plot(summary: pd.DataFrame, frame_noise: pd.DataFrame, output: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.4), constrained_layout=True)
    frame = frame_noise[frame_noise.noise_metric == "robust"].copy()
    ax = axes.flat[0]
    for cut, group in frame.groupby("snr_min", sort=True):
        group = group.sort_values("fwhm_arcsec")
        x = group.fwhm_arcsec.to_numpy(float)
        if not np.isfinite(x).any():
            x = group.fwhm_px.to_numpy(float); xlabel = "FWHM (px)"
        else:
            xlabel = "FWHM (arcsec)"
        ax.plot(x, group.frame_scatter, "o-", label=f"SNR >= {cut:g}")
    ax.set_xlabel(xlabel)
    ax.set_title("Frame quality (NNLS intrinsic robust scatter)"); ax.set_ylabel("mag")
    ax.grid(axis="y", alpha=.25); ax.legend(fontsize=8)
    for ax, scope in zip(axes.flat[1:], ("snr", "radius", "neighbor")):
        group = summary[summary.scope == scope].copy()
        x = np.arange(len(group)); labels = group.label.astype(str).tolist()
        ax.plot(x, group.cv_robust_scatter_mad, "o-", label="CV robust scatter (MAD)")
        ax.plot(x, group.cv_rmse, "s--", label="CV RMSE")
        ax.set_title(scope.capitalize()); ax.set_ylabel("mag"); ax.set_xticks(x); ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.grid(axis="y", alpha=.25); ax.legend(fontsize=8)
    fig.suptitle("PSF repeatability: frame noise decomposition and exact CV predictive disagreement")
    fig.savefig(output, dpi=220); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("result_dir", type=Path); parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--match-radius-px", type=float, default=1.5); parser.add_argument("--min-frames", type=int, default=3)
    parser.add_argument("--snr-min", type=float); parser.add_argument("--qfit-max", type=float); parser.add_argument("--crowding-min-px", type=float)
    args = parser.parse_args(); output = args.output_dir or args.result_dir / "psf_repeatability"
    result = analyze_psf_repeatability(args.result_dir, output, args.match_radius_px, args.min_frames, args.snr_min, args.qfit_max, args.crowding_min_px)
    frame_noise = pd.read_csv(output / "psf_repeatability_frame_noise.csv")
    _plot(pd.DataFrame(result["summary"]), frame_noise, output / "psf_repeatability.png")
    print(f"Wrote {output / 'psf_repeatability_summary.csv'}")
    print(pd.DataFrame(result["summary"]).query("scope == 'overall'").to_string(index=False))


if __name__ == "__main__":
    main()
