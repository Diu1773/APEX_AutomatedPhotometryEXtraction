"""Compare hard-core masking with full-field forced PSF policies on M13."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO))

from validation.compare_m13_psf_variants import SNR_BINS, SNR_CUTS, _load_variant


ROOT = Path("validation/real_gui_run/m13_variants")
OUTPUT_DIR = Path("validation/real_gui_run")
FILENAME = "pp_messier13-0003-B.fit"
SPECS = (
    ("masked, free/positive", ROOT / "epsf57" / "result"),
    ("full field, free/positive", ROOT / "epsf57_no_core" / "result"),
    ("full field, fixed/signed", ROOT / "epsf57_no_core_forced_all" / "result"),
)
COLORS = ("#777777", "#D55E00", "#0072B2")


def _catalogue_metrics(result_dir: Path, center_x: float, center_y: float, radius_px: float) -> dict:
    path = result_dir / "cmd_psf" / f"photometry_{FILENAME}.tsv"
    frame = pd.read_csv(path, sep="\t")
    for column in (
        "x_fit",
        "y_fit",
        "flux_psf_e",
        "mag_psf",
        "flags_psf",
        "forced_psf",
    ):
        if column in frame:
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
    radius = np.hypot(frame["x_fit"] - center_x, frame["y_fit"] - center_y)
    good = np.isfinite(frame["mag_psf"]) & (frame["flags_psf"] == 0)
    forced = frame.get("forced_psf", pd.Series(False, index=frame.index)).fillna(False).astype(bool)
    flux = frame["flux_psf_e"]
    flags = frame["flags_psf"].fillna(0).astype(int)
    return {
        "core_rows": int(np.sum(radius < radius_px)),
        "core_good": int(np.sum((radius < radius_px) & good)),
        "outside_good": int(np.sum((radius >= radius_px) & good)),
        "forced_rows": int(np.sum(forced)),
        "forced_negative": int(np.sum(forced & np.isfinite(flux) & (flux <= 0))),
        "crowding_flagged": int(np.sum((flags & 8192) != 0)),
    }


def main() -> int:
    masked_meta = json.loads(
        next((SPECS[0][1] / "cmd_psf").glob("residual_meta_*.json")).read_text(encoding="utf-8")
    )
    core = masked_meta["core_cut"]
    center_x = float(core["center_x"])
    center_y = float(core["center_y"])
    radius_px = float(core["radius_px"])

    loaded = [_load_variant(label, path) for label, path in SPECS]
    rows = []
    for (_, result_dir), variant in zip(SPECS, loaded):
        row = dict(variant["summary"])
        row.update(_catalogue_metrics(result_dir, center_x, center_y, radius_px))
        rows.append(row)
    summary = pd.DataFrame(rows)

    fig, axes = plt.subplots(2, 3, figsize=(15.0, 8.4), constrained_layout=True)
    x_bins = np.arange(len(SNR_BINS))
    for (_, row), color in zip(summary.iterrows(), COLORS):
        axes[0, 0].plot(
            x_bins,
            [row[f"bias_{label}"] for label, _, _ in SNR_BINS],
            marker="o",
            color=color,
            label=row["variant"],
        )
        axes[0, 1].plot(
            SNR_CUTS,
            [row[f"scatter_snr_ge_{cut}"] for cut in SNR_CUTS],
            marker="o",
            color=color,
            label=row["variant"],
        )

    axes[0, 0].axhline(0.0, color="0.4", linestyle="--", linewidth=0.8)
    axes[0, 0].set_xticks(x_bins, [label for label, _, _ in SNR_BINS])
    axes[0, 0].set_xlabel("aperture SNR bin")
    axes[0, 0].set_ylabel("median bias after high-SNR centering (mag)")
    axes[0, 0].set_title("(a) Low-SNR flux bias", loc="left")
    axes[0, 0].grid(alpha=0.2)
    axes[0, 0].legend(fontsize=8)

    axes[0, 1].set_xscale("log")
    axes[0, 1].set_xlabel("minimum aperture SNR")
    axes[0, 1].set_ylabel("robust PSF-aperture scatter (mag)")
    axes[0, 1].set_title("(b) Scatter by SNR cut", loc="left")
    axes[0, 1].grid(alpha=0.2)

    x = np.arange(len(summary))
    width = 0.36
    axes[0, 2].bar(x - width / 2, summary["n_good"], width, color="#009E73", label="good mag")
    axes[0, 2].bar(x + width / 2, summary["n_fail"], width, color="#CC79A7", label="failed/flagged")
    axes[0, 2].set_xticks(x, summary["variant"], rotation=16, ha="right", fontsize=8)
    axes[0, 2].set_ylabel("sources")
    axes[0, 2].set_title("(c) Output quality", loc="left")
    axes[0, 2].legend(fontsize=8)

    axes[1, 0].bar(x, summary["step8_elapsed_s"], color=COLORS)
    axes[1, 0].set_xticks(x, summary["variant"], rotation=16, ha="right", fontsize=8)
    axes[1, 0].set_ylabel("seconds, one CPU core")
    axes[1, 0].set_title("(d) Step 8 runtime", loc="left")

    axes[1, 1].bar(x - width / 2, summary["core_good"], width, color="#56B4E9", label="core good")
    axes[1, 1].bar(x + width / 2, summary["outside_good"], width, color="#009E73", label="outside good")
    axes[1, 1].set_xticks(x, summary["variant"], rotation=16, ha="right", fontsize=8)
    axes[1, 1].set_ylabel("good magnitude sources")
    axes[1, 1].set_title(f"(e) Retained at r={radius_px:.1f}px", loc="left")
    axes[1, 1].legend(fontsize=8)

    axes[1, 2].bar(x - width / 2, summary["forced_negative"], width, color="#0072B2", label="signed negative")
    axes[1, 2].bar(x + width / 2, summary["crowding_flagged"], width, color="#E69F00", label="unresolved flag")
    axes[1, 2].set_xticks(x, summary["variant"], rotation=16, ha="right", fontsize=8)
    axes[1, 2].set_ylabel("sources")
    axes[1, 2].set_title("(f) Honest forced/QC outcomes", loc="left")
    axes[1, 2].legend(fontsize=8)

    fig.suptitle("M13 PSF core policy: radial deletion vs fixed signed forced fitting", fontsize=15)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    figure_path = OUTPUT_DIR / "m13_psf_core_policy_comparison.png"
    csv_path = OUTPUT_DIR / "m13_psf_core_policy_comparison.csv"
    fig.savefig(figure_path, dpi=200, facecolor="white")
    plt.close(fig)
    summary.to_csv(csv_path, index=False)
    print(figure_path.resolve())
    print(csv_path.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
