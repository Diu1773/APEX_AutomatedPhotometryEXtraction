"""Compare controlled M13 PSF configuration variants."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO))

from validation.analyze_m13_psf_asymmetry import (
    _epsf_metrics,
    _match_catalogues,
    _quantile_stats,
)


SNR_CUTS = (5, 10, 20, 50, 100)
SNR_BINS = (
    ("5-10", 5, 10),
    ("10-20", 10, 20),
    ("20-50", 20, 50),
    ("50-100", 50, 100),
    (">=100", 100, np.inf),
)


def _robust_scatter(values: pd.Series | np.ndarray) -> float:
    data = np.asarray(values, dtype=float)
    data = data[np.isfinite(data)]
    if data.size == 0:
        return np.nan
    median = float(np.median(data))
    return float(1.4826 * np.median(np.abs(data - median)))


def _load_variant(label: str, result_dir: Path) -> dict:
    psf_dir = result_dir / "cmd_psf"
    meta = json.loads(next(psf_dir.glob("residual_meta_*.json")).read_text(encoding="utf-8"))
    index = pd.read_csv(psf_dir / "photometry_index.csv").iloc[0]
    run = json.loads((result_dir / "real_gui_run.json").read_text(encoding="utf-8"))
    matched = _match_catalogues(result_dir, meta["file"])
    reference = matched.loc[matched["snr_aperture"] >= 100, "delta_mag"]
    reference_median = float(np.median(reference))
    matched["delta_centered_mag"] = matched["delta_mag"] - reference_median
    epsf = _epsf_metrics(next(psf_dir.glob("epsf_model_*.fits")))

    summary = {
        "variant": label,
        "n_fit": int(index["n"]),
        "n_good": int(index["n_goodmag"]),
        "n_fail": int(index["n_fail"]),
        "median_qfit": float(index["median_qfit"]),
        "median_reduced_chi2": float(index["median_reduced_chi2"]),
        "epsf_selected": int(index["epsf_n_selected"]),
        "step8_elapsed_s": float(run["step8_elapsed_s"]),
        "residual_std": float(meta["iters"][-1]["residual_std"]),
        "epsf_rotation_asymmetry": float(epsf["rotation_asymmetry"]),
        "epsf_ellipticity": float(epsf["ellipticity"]),
        "high_snr_offset_mag": reference_median,
        "matched": int(len(matched)),
    }
    for cut in SNR_CUTS:
        values = matched.loc[matched["snr_aperture"] >= cut, "delta_mag"]
        summary[f"scatter_snr_ge_{cut}"] = _robust_scatter(values)
        summary[f"n_snr_ge_{cut}"] = int(len(values))
    for bin_label, low, high in SNR_BINS:
        selected = (matched["snr_aperture"] >= low) & (matched["snr_aperture"] < high)
        values = matched.loc[selected, "delta_centered_mag"]
        summary[f"bias_{bin_label}"] = float(np.median(values)) if len(values) else np.nan
        summary[f"asym_{bin_label}"] = float(_quantile_stats(values)["quantile_asymmetry"])
    for detected, class_label in ((True, "detected"), (False, "forced")):
        values = matched.loc[matched["detected"] == detected, "delta_centered_mag"]
        summary[f"scatter_{class_label}"] = _robust_scatter(values)
        summary[f"median_{class_label}"] = float(np.median(values)) if len(values) else np.nan
        summary[f"n_{class_label}"] = int(len(values))
    return {"summary": summary, "matched": matched}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--validation-root",
        type=Path,
        default=Path(r"E:\APEX_validation\reprocess"),
    )
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

    specs = (
        ("30 stars, 8 iter", args.validation_root / "M13" / "result"),
        ("57 stars, 8 iter", args.variant_root / "epsf57" / "result"),
        ("57 + groups", args.variant_root / "epsf57_grouped" / "result"),
        ("57 stars, 16 iter", args.variant_root / "epsf57_iter16" / "result"),
    )
    variants = [_load_variant(label, path) for label, path in specs]
    summary = pd.DataFrame([variant["summary"] for variant in variants])

    colors = ("#999999", "#0072B2", "#D55E00", "#009E73")
    fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.2), constrained_layout=True)
    for row, color in zip(summary.itertuples(index=False), colors):
        axes[0, 0].plot(
            SNR_CUTS,
            [getattr(row, f"scatter_snr_ge_{cut}") for cut in SNR_CUTS],
            marker="o",
            label=row.variant,
            color=color,
        )
        axes[0, 1].plot(
            np.arange(len(SNR_BINS)),
            [getattr(row, f"bias_{label.replace('>=', '_')}", np.nan) for label, _, _ in SNR_BINS],
            marker="o",
            label=row.variant,
            color=color,
        )

    axes[0, 0].set_xscale("log")
    axes[0, 0].set_xlabel("minimum aperture SNR")
    axes[0, 0].set_ylabel("robust scatter (mag)")
    axes[0, 0].set_title("(a) PSF-aperture scatter", loc="left")
    axes[0, 0].grid(alpha=0.2)
    axes[0, 0].legend(fontsize=8)

    # Pandas sanitizes namedtuple fields, so draw low-SNR bias directly from the frame.
    axes[0, 1].cla()
    for (_, row), color in zip(summary.iterrows(), colors):
        axes[0, 1].plot(
            np.arange(len(SNR_BINS)),
            [row[f"bias_{label}"] for label, _, _ in SNR_BINS],
            marker="o",
            label=row["variant"],
            color=color,
        )
    axes[0, 1].axhline(0.0, color="0.4", linewidth=0.8, linestyle="--")
    axes[0, 1].set_xticks(np.arange(len(SNR_BINS)), [row[0] for row in SNR_BINS])
    axes[0, 1].set_xlabel("aperture SNR bin")
    axes[0, 1].set_ylabel("median bias after high-SNR centering (mag)")
    axes[0, 1].set_title("(b) low-SNR flux bias", loc="left")
    axes[0, 1].grid(alpha=0.2)

    x = np.arange(len(summary))
    axes[0, 2].bar(x, summary["n_fail"], color=colors)
    axes[0, 2].set_xticks(x, summary["variant"], rotation=18, ha="right", fontsize=8)
    axes[0, 2].set_ylabel("flagged / failed sources")
    axes[0, 2].set_title("(c) numerical failures", loc="left")

    width = 0.36
    axes[1, 0].bar(x - width / 2, summary["scatter_detected"], width, label="Step4 detected", color="#56B4E9")
    axes[1, 0].bar(x + width / 2, summary["scatter_forced"], width, label="forced only", color="#E69F00")
    axes[1, 0].set_xticks(x, summary["variant"], rotation=18, ha="right", fontsize=8)
    axes[1, 0].set_ylabel("robust scatter (mag)")
    axes[1, 0].set_title("(d) catalogue class", loc="left")
    axes[1, 0].legend(fontsize=8)

    axes[1, 1].bar(x, summary["epsf_rotation_asymmetry"], color=colors)
    axes[1, 1].set_xticks(x, summary["variant"], rotation=18, ha="right", fontsize=8)
    axes[1, 1].set_ylabel("ePSF 180-degree asymmetry")
    axes[1, 1].set_title("(e) ePSF shape", loc="left")

    axes[1, 2].bar(x, summary["step8_elapsed_s"], color=colors)
    axes[1, 2].set_xticks(x, summary["variant"], rotation=18, ha="right", fontsize=8)
    axes[1, 2].set_ylabel("Step8 elapsed time (s)")
    axes[1, 2].set_title("(f) one-CPU runtime", loc="left")
    fig.suptitle("Controlled M13 PSF configuration comparison", fontsize=15)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    figure_path = args.output_dir / "m13_psf_variant_comparison.png"
    csv_path = args.output_dir / "m13_psf_variant_comparison.csv"
    fig.savefig(figure_path, dpi=200, facecolor="white")
    plt.close(fig)
    summary.to_csv(csv_path, index=False)
    print(figure_path.resolve())
    print(csv_path.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
