"""Compare M5 artificial-star recovery for good and poor seeing frames."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


SNR_ORDER = (5.0, 10.0, 20.0, 50.0, 100.0)
CROWDING_ORDER = (
    "0.75-1.5 FWHM",
    "1.5-3 FWHM",
    "3-6 FWHM",
    "6-inf FWHM",
)


def _rows(summary: pd.DataFrame, scope: str, labels) -> pd.DataFrame:
    frame = summary[summary["scope"] == scope].copy()
    frame["label"] = frame["label"].astype(str)
    order = {str(label): index for index, label in enumerate(labels)}
    frame["_order"] = frame["label"].map(order)
    return frame[frame["_order"].notna()].sort_values("_order")


def _load(path: Path, label: str, fwhm_arcsec: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = pd.read_csv(path / "summary.csv")
    recovery = pd.read_csv(path / "recovery.csv")
    summary["seeing_label"] = label
    summary["fwhm_arcsec"] = float(fwhm_arcsec)
    recovery["seeing_label"] = label
    recovery["fwhm_arcsec"] = float(fwhm_arcsec)
    return summary, recovery


def _plot(summary: pd.DataFrame, recovery: pd.DataFrame, output: Path) -> None:
    styles = {
        "good": {"color": "#0072B2", "marker": "o"},
        "poor": {"color": "#D55E00", "marker": "s"},
    }
    fig, axes = plt.subplots(2, 3, figsize=(15.5, 8.8), constrained_layout=True)

    for name, group in summary.groupby("seeing_label", sort=False):
        style = styles[name]
        snr = _rows(group, "target_snr", SNR_ORDER)
        x = pd.to_numeric(snr["label"], errors="coerce").to_numpy(float)
        legend = f"{name} ({snr['fwhm_arcsec'].iloc[0]:.2f} arcsec)"
        axes[0, 0].plot(
            x, snr["psf_completeness"],
            color=style["color"], marker=style["marker"], label=legend,
        )
        axes[0, 1].plot(
            x, snr["delta_mag_bias"],
            color=style["color"], marker=style["marker"], label=legend,
        )
        axes[0, 2].plot(
            x, snr["delta_mag_scatter_robust"],
            color=style["color"], marker=style["marker"], label=legend,
        )
        axes[1, 1].plot(
            x, snr["recommended_error_scale"],
            color=style["color"], marker=style["marker"], label=legend,
        )

        crowd = _rows(group, "crowding", CROWDING_ORDER)
        cx = np.arange(len(CROWDING_ORDER), dtype=float)
        offset = -0.18 if name == "good" else 0.18
        axes[1, 0].bar(
            cx + offset,
            crowd["psf_completeness"],
            width=0.34,
            color=style["color"],
            alpha=0.85,
            label=legend,
        )

    snr_grid = np.logspace(np.log10(5.0), np.log10(100.0), 200)
    expected_sigma = 1.085736 / snr_grid
    axes[0, 2].plot(
        snr_grid, expected_sigma, color="0.25", linestyle=":", linewidth=1.2,
        label="photon limit: 1.086 / SNR",
    )

    axes[0, 0].set_title("(a) forced-PSF completeness", loc="left")
    axes[0, 0].set_ylabel("flags=0 recovery fraction")
    axes[0, 0].set_ylim(-0.05, 1.05)

    axes[0, 1].axhline(0.0, color="0.35", linestyle=":", linewidth=0.9)
    axes[0, 1].set_title("(b) magnitude bias", loc="left")
    axes[0, 1].set_ylabel("median recovered - true (mag)")

    axes[0, 2].set_title("(c) empirical precision", loc="left")
    axes[0, 2].set_ylabel("robust scatter (mag)")
    axes[0, 2].set_yscale("log")

    axes[1, 0].set_title("(d) crowding limit", loc="left")
    axes[1, 0].set_xticks(np.arange(len(CROWDING_ORDER)), CROWDING_ORDER, rotation=18)
    axes[1, 0].set_ylabel("flags=0 recovery fraction")
    axes[1, 0].set_ylim(-0.05, 1.05)

    axes[1, 1].axhline(1.0, color="0.35", linestyle=":", linewidth=0.9)
    axes[1, 1].set_title("(e) formal-error calibration", loc="left")
    axes[1, 1].set_ylabel("recommended error scale")

    ax = axes[1, 2]
    for name, group in recovery.groupby("seeing_label", sort=False):
        style = styles[name]
        valid = group[
            group["psf_recovered"].astype(str).str.lower().isin({"true", "1"})
            & np.isfinite(pd.to_numeric(group["delta_mag"], errors="coerce"))
        ].copy()
        rng = np.random.default_rng(420 if name == "good" else 421)
        jitter = np.exp(rng.normal(0.0, 0.025, len(valid)))
        ax.scatter(
            pd.to_numeric(valid["target_snr"], errors="coerce") * jitter,
            pd.to_numeric(valid["delta_mag"], errors="coerce"),
            s=18, alpha=0.65, color=style["color"], marker=style["marker"],
            label=f"{name} recovered",
        )
    ax.plot(snr_grid, expected_sigma, color="0.25", linestyle=":", linewidth=1.0)
    ax.plot(snr_grid, -expected_sigma, color="0.25", linestyle=":", linewidth=1.0)
    ax.axhline(0.0, color="0.45", linewidth=0.8)
    ax.set_title("(f) recovered artificial stars", loc="left")
    ax.set_ylabel("recovered - true (mag)")
    ax.set_ylim(-0.75, 0.75)

    for axis in axes.flat:
        axis.set_xscale("log") if axis not in (axes[1, 0],) else None
        if axis not in (axes[1, 0],):
            axis.set_xlabel("target SNR")
        axis.grid(alpha=0.18)
        handles, labels = axis.get_legend_handles_labels()
        if handles:
            axis.legend(fontsize=8, frameon=False)

    fig.suptitle(
        "M5 CPU-only PSF artificial-star validation (2 trials, 50 injections/frame)",
        fontsize=15,
    )
    fig.savefig(output, dpi=210, facecolor="white")
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--good", type=Path, required=True)
    parser.add_argument("--poor", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--good-fwhm", type=float, default=2.223461)
    parser.add_argument("--poor-fwhm", type=float, default=4.055816)
    args = parser.parse_args()

    good_summary, good_recovery = _load(args.good, "good", args.good_fwhm)
    poor_summary, poor_recovery = _load(args.poor, "poor", args.poor_fwhm)
    summary = pd.concat([good_summary, poor_summary], ignore_index=True)
    recovery = pd.concat([good_recovery, poor_recovery], ignore_index=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_dir / "m5_artificial_star_comparison.csv", index=False)
    output = args.output_dir / "m5_artificial_star_comparison.png"
    _plot(summary, recovery, output)
    print(output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
