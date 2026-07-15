"""Validate contamination-aware ePSF reference selection on the M13 run."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from validation.compare_m13_psf_variants import _load_variant, _robust_scatter


SNR_BINS = (
    ("5-10", 5.0, 10.0),
    ("10-20", 10.0, 20.0),
    ("20-50", 20.0, 50.0),
    ("50-100", 50.0, 100.0),
    (">=100", 100.0, np.inf),
)


def _bin_metrics(variant: dict) -> pd.DataFrame:
    matched = variant["matched"]
    rows = []
    for label, low, high in SNR_BINS:
        keep = (matched["snr_aperture"] >= low) & (matched["snr_aperture"] < high)
        values = matched.loc[keep, "delta_centered_mag"].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        rows.append({
            "variant": variant["summary"]["variant"],
            "snr_bin": label,
            "n": int(values.size),
            "bias_mag": float(np.median(values)) if values.size else np.nan,
            "robust_scatter_mag": _robust_scatter(values),
            "rmse_mag": float(np.sqrt(np.mean(values**2))) if values.size else np.nan,
            "outlier_fraction_gt_0p2mag": (
                float(np.mean(np.abs(values) > 0.2)) if values.size else np.nan
            ),
        })
    return pd.DataFrame(rows)


def _reference_catalog(result_dir: Path) -> tuple[pd.DataFrame, dict]:
    psf_dir = result_dir / "cmd_psf"
    meta_path = next(psf_dir.glob("residual_meta_*.json"))
    meta = json.loads(meta_path.read_text(encoding="utf-8"))
    catalog_name = meta["epsf_reference"]["catalog_path"]
    reference = pd.read_csv(psf_dir / catalog_name)
    for column in ("selected", "isolated", "low_contamination", "core_safe"):
        reference[column] = reference[column].astype(str).str.lower().eq("true")
    return reference, meta


def _plot_metric(ax, table: pd.DataFrame, column: str, ylabel: str, title: str) -> None:
    colors = ("#777777", "#0072B2")
    x = np.arange(len(SNR_BINS))
    for (label, group), color in zip(table.groupby("variant", sort=False), colors):
        ordered = group.set_index("snr_bin").loc[[item[0] for item in SNR_BINS]]
        ax.plot(x, ordered[column], marker="o", linewidth=1.8, color=color, label=label)
    ax.set_xticks(x, [item[0] for item in SNR_BINS])
    ax.set_xlabel("Step7 aperture SNR bin")
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left")
    ax.grid(alpha=0.2)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--variant-root",
        type=Path,
        default=Path("validation/real_gui_run/m13_variants"),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("validation/real_gui_run"),
    )
    args = parser.parse_args()

    baseline_dir = args.variant_root / "epsf57_no_core_forced_all" / "result"
    aware_dir = args.variant_root / "epsf_contam" / "result"
    baseline = _load_variant("candidate-only isolation", baseline_dir)
    aware = _load_variant("all-source + local contamination", aware_dir)
    variants = (baseline, aware)
    binned = pd.concat([_bin_metrics(item) for item in variants], ignore_index=True)
    reference, meta = _reference_catalog(aware_dir)

    core = meta["core_cut"]
    radius = max(float(core["radius_px"]), 1.0)
    finite = np.isfinite(reference["contamination_score"])
    rho, p_value = spearmanr(
        reference.loc[finite, "core_distance_px"],
        reference.loc[finite, "contamination_score"],
    )

    fig, axes = plt.subplots(2, 3, figsize=(16.2, 9.2), constrained_layout=True)
    rejected = ~reference["selected"]
    points = axes[0, 0].scatter(
        reference.loc[rejected, "x"],
        reference.loc[rejected, "y"],
        c=np.log10(1.0 + reference.loc[rejected, "contamination_score"].clip(lower=0.0)),
        s=18,
        cmap="inferno",
        alpha=0.75,
        linewidths=0,
        label="rejected candidate",
    )
    selected = reference["selected"]
    axes[0, 0].scatter(
        reference.loc[selected, "x"],
        reference.loc[selected, "y"],
        facecolors="none",
        edgecolors="#0072B2",
        marker="s",
        s=48,
        linewidths=1.25,
        label="selected ePSF reference",
    )
    axes[0, 0].add_patch(Circle(
        (float(core["center_x"]), float(core["center_y"])),
        radius,
        facecolor="none",
        edgecolor="#D55E00",
        linewidth=1.5,
        label="reference exclusion radius",
    ))
    axes[0, 0].set_aspect("equal", adjustable="box")
    axes[0, 0].set_xlabel("detector x (px)")
    axes[0, 0].set_ylabel("detector y (px)")
    axes[0, 0].set_title("(a) ePSF reference selection", loc="left")
    axes[0, 0].legend(fontsize=8, loc="upper left")
    colorbar = fig.colorbar(points, ax=axes[0, 0], fraction=0.046, pad=0.03)
    colorbar.set_label("log10(1 + contamination score)")

    axes[0, 1].scatter(
        reference.loc[rejected, "core_distance_px"] / radius,
        reference.loc[rejected, "contamination_score"],
        color="#D55E00",
        s=18,
        alpha=0.55,
        label="rejected",
    )
    axes[0, 1].scatter(
        reference.loc[selected, "core_distance_px"] / radius,
        reference.loc[selected, "contamination_score"],
        facecolors="none",
        edgecolors="#0072B2",
        marker="s",
        s=42,
        linewidths=1.1,
        label="selected",
    )
    axes[0, 1].axvline(1.0, color="0.35", linestyle="--", linewidth=1.0)
    axes[0, 1].set_yscale("symlog", linthresh=0.05)
    axes[0, 1].set_xlabel("distance from cluster center / reference radius")
    axes[0, 1].set_ylabel("local annulus contamination score")
    axes[0, 1].set_title(f"(b) radial prior is not enough: rho={rho:.2f}, p={p_value:.2g}", loc="left")
    axes[0, 1].grid(alpha=0.2)
    axes[0, 1].legend(fontsize=8)

    _plot_metric(
        axes[0, 2], binned, "robust_scatter_mag", "robust scatter (mag)",
        "(c) PSF-aperture robust scatter",
    )
    axes[0, 2].legend(fontsize=8)
    _plot_metric(
        axes[1, 0], binned, "rmse_mag", "RMSE (mag)",
        "(d) PSF-aperture RMSE",
    )
    _plot_metric(
        axes[1, 1], binned, "bias_mag", "median centered delta (mag)",
        "(e) SNR-dependent bias",
    )
    axes[1, 1].axhline(0.0, color="0.35", linestyle="--", linewidth=1.0)

    hist_colors = ("#777777", "#0072B2")
    for variant, color in zip(variants, hist_colors):
        matched = variant["matched"]
        values = matched.loc[
            matched["snr_aperture"] >= 50,
            "delta_centered_mag",
        ].to_numpy(dtype=float)
        values = values[np.isfinite(values)]
        axes[1, 2].hist(
            values,
            bins=np.linspace(-0.35, 0.35, 36),
            histtype="step",
            density=True,
            linewidth=1.8,
            color=color,
            label=f"{variant['summary']['variant']} (N={values.size})",
        )
    axes[1, 2].axvline(0.0, color="0.35", linestyle="--", linewidth=1.0)
    axes[1, 2].set_xlabel("centered PSF - aperture (mag), SNR >= 50")
    axes[1, 2].set_ylabel("density")
    axes[1, 2].set_title("(f) high-SNR error floor", loc="left")
    axes[1, 2].legend(fontsize=8)

    epsf_meta = meta["epsf_reference"]
    fig.suptitle(
        "M13 contamination-aware ePSF validation: "
        f"{epsf_meta['n_candidates']} candidates, {epsf_meta['n_isolated']} isolated, "
        f"{epsf_meta['n_selected']} selected",
        fontsize=15,
    )

    summary_rows = []
    for variant in variants:
        row = dict(variant["summary"])
        row.update({
            "reference_n_candidates": (
                epsf_meta["n_candidates"] if variant is aware else np.nan
            ),
            "reference_n_isolated": (
                epsf_meta["n_isolated"] if variant is aware else np.nan
            ),
            "reference_n_core_rejected": (
                epsf_meta["n_core_rejected"] if variant is aware else np.nan
            ),
            "reference_median_contamination": (
                epsf_meta["selected_median_contamination"] if variant is aware else np.nan
            ),
        })
        summary_rows.append(row)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure_path = args.output_dir / "m13_epsf_contamination_comparison.png"
    binned_path = args.output_dir / "m13_epsf_contamination_snr_bins.csv"
    summary_path = args.output_dir / "m13_epsf_contamination_summary.csv"
    fig.savefig(figure_path, dpi=200, facecolor="white")
    plt.close(fig)
    binned.to_csv(binned_path, index=False)
    pd.DataFrame(summary_rows).to_csv(summary_path, index=False)
    print(figure_path.resolve())
    print(binned_path.resolve())
    print(summary_path.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
