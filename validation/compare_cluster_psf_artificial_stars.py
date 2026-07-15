"""Build a cross-cluster PSF artificial-star validation figure."""

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
DATASETS = (
    ("M3 3.41 arcsec", "m3", "#009E73", "o"),
    ("M13 2.73 arcsec", "m13", "#CC79A7", "D"),
    ("M5 good 2.22 arcsec", "m5_good", "#0072B2", "^"),
    ("M5 poor 4.06 arcsec", "m5_poor", "#D55E00", "s"),
)


def _ordered(frame: pd.DataFrame, scope: str, labels) -> pd.DataFrame:
    out = frame[frame["scope"] == scope].copy()
    order = {str(label): index for index, label in enumerate(labels)}
    out["_order"] = out["label"].astype(str).map(order)
    return out[out["_order"].notna()].sort_values("_order")


def _robust_scatter(values: np.ndarray) -> float:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return np.nan
    median = float(np.median(finite))
    return float(1.4826 * np.median(np.abs(finite - median)))


def _load(path: Path, key: str, label: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary = pd.read_csv(path / "summary.csv")
    recovery = pd.read_csv(path / "recovery.csv")
    summary["dataset"] = key
    summary["dataset_label"] = label
    recovery["dataset"] = key
    recovery["dataset_label"] = label
    return summary, recovery


def _recovered_mask(frame: pd.DataFrame) -> np.ndarray:
    recovered = frame["psf_recovered"]
    if pd.api.types.is_bool_dtype(recovered):
        return recovered.fillna(False).to_numpy(bool)
    return recovered.astype(str).str.lower().isin({"true", "1"}).to_numpy(bool)


def _plot(summary: pd.DataFrame, recovery: pd.DataFrame, output: Path) -> pd.DataFrame:
    fig, axes = plt.subplots(2, 3, figsize=(16.0, 9.0), constrained_layout=True)
    crowd_x = np.arange(len(CROWDING_ORDER), dtype=float)
    bar_offsets = np.linspace(-0.30, 0.30, len(DATASETS))
    high_snr_rows = []

    for dataset_index, (label, key, color, marker) in enumerate(DATASETS):
        group = summary[summary["dataset"] == key]
        snr = _ordered(group, "target_snr", SNR_ORDER)
        x = pd.to_numeric(snr["label"], errors="coerce").to_numpy(float)
        axes[0, 0].plot(x, snr["psf_completeness"], marker=marker, color=color, label=label)
        axes[0, 1].plot(x, snr["delta_mag_bias"], marker=marker, color=color, label=label)
        axes[0, 2].plot(
            x, snr["delta_mag_scatter_robust"], marker=marker, color=color, label=label
        )

        crowd = _ordered(group, "crowding", CROWDING_ORDER)
        axes[1, 0].bar(
            crowd_x + bar_offsets[dataset_index],
            crowd["psf_completeness"],
            width=0.19,
            color=color,
            alpha=0.85,
            label=label,
        )
        axes[1, 1].plot(
            crowd_x,
            crowd["delta_mag_bias"],
            marker=marker,
            color=color,
            label=label,
        )

        rec = recovery[recovery["dataset"] == key].copy()
        rec = rec[
            _recovered_mask(rec)
            & (pd.to_numeric(rec["target_snr"], errors="coerce") >= 50)
            & (pd.to_numeric(rec["nearest_real_sep_fwhm"], errors="coerce") >= 6)
        ]
        values = pd.to_numeric(rec["delta_mag"], errors="coerce").to_numpy(float)
        high_snr_rows.append({
            "dataset": key,
            "dataset_label": label,
            "n": int(np.isfinite(values).sum()),
            "bias_mag": float(np.nanmedian(values)) if np.isfinite(values).any() else np.nan,
            "scatter_mag": _robust_scatter(values),
        })

    grid = np.logspace(np.log10(5.0), np.log10(100.0), 200)
    axes[0, 2].plot(
        grid, 1.085736 / grid, color="0.2", linestyle=":", linewidth=1.2,
        label="photon limit: 1.086 / SNR",
    )

    axes[0, 0].set_title("(a) forced-PSF completeness", loc="left")
    axes[0, 0].set_ylabel("flags=0 recovery fraction")
    axes[0, 0].set_ylim(-0.05, 1.05)

    axes[0, 1].axhline(0.0, color="0.35", linestyle=":", linewidth=0.9)
    axes[0, 1].set_title("(b) magnitude bias", loc="left")
    axes[0, 1].set_ylabel("median recovered - true (mag)")

    axes[0, 2].set_title("(c) conditional precision", loc="left")
    axes[0, 2].set_ylabel("robust scatter for flags=0 (mag)")
    axes[0, 2].set_yscale("log")

    axes[1, 0].set_title("(d) common crowding boundary", loc="left")
    axes[1, 0].set_ylabel("flags=0 recovery fraction")
    axes[1, 0].set_ylim(-0.05, 1.05)
    axes[1, 0].set_xticks(crowd_x, CROWDING_ORDER, rotation=18)

    axes[1, 1].axhline(0.0, color="0.35", linestyle=":", linewidth=0.9)
    axes[1, 1].set_title("(e) crowding-dependent bias", loc="left")
    axes[1, 1].set_ylabel("median recovered - true (mag)")
    axes[1, 1].set_xticks(crowd_x, CROWDING_ORDER, rotation=18)

    high = pd.DataFrame(high_snr_rows)
    hx = np.arange(len(high), dtype=float)
    colors = [item[2] for item in DATASETS]
    axes[1, 2].bar(hx, high["bias_mag"], color=colors, alpha=0.85)
    axes[1, 2].errorbar(
        hx,
        high["bias_mag"],
        yerr=high["scatter_mag"],
        fmt="none",
        ecolor="0.15",
        capsize=4,
        linewidth=1.0,
    )
    axes[1, 2].axhline(0.0, color="0.35", linestyle=":", linewidth=0.9)
    axes[1, 2].set_xticks(hx, [f"{row.dataset}\nN={row.n}" for row in high.itertuples()])
    axes[1, 2].set_ylabel("bias +/- robust scatter (mag)")
    axes[1, 2].set_title("(f) SNR >= 50 and neighbor >= 6 FWHM", loc="left")

    for axis in axes[0, :]:
        axis.set_xscale("log")
        axis.set_xlabel("target SNR")
    for axis in axes.flat:
        axis.grid(alpha=0.18)
    for axis in (axes[0, 0], axes[0, 1], axes[0, 2]):
        axis.legend(fontsize=7.5, frameon=False)

    fig.suptitle(
        "APEX CPU-only PSF artificial-star validation across globular clusters",
        fontsize=15,
    )
    fig.savefig(output, dpi=210, facecolor="white")
    plt.close(fig)
    return high


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--m3", type=Path, required=True)
    parser.add_argument("--m13", type=Path, required=True)
    parser.add_argument("--m5-good", type=Path, required=True)
    parser.add_argument("--m5-poor", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    paths = {"m3": args.m3, "m13": args.m13, "m5_good": args.m5_good, "m5_poor": args.m5_poor}
    summaries = []
    recoveries = []
    for label, key, _, _ in DATASETS:
        summary, recovery = _load(paths[key], key, label)
        summaries.append(summary)
        recoveries.append(recovery)
    summary_all = pd.concat(summaries, ignore_index=True)
    recovery_all = pd.concat(recoveries, ignore_index=True)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary_all.to_csv(args.output_dir / "cluster_psf_artificial_summary.csv", index=False)
    output = args.output_dir / "cluster_psf_artificial_comparison.png"
    high = _plot(summary_all, recovery_all, output)
    high.to_csv(args.output_dir / "cluster_psf_high_snr_isolated.csv", index=False)
    print(high.to_string(index=False))
    print(output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
