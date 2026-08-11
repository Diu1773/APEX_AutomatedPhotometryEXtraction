"""What PSF photometry buys in a globular cluster, measured on the same frames.

APEX's batch reprocessing runs steps 1-7 and then step 10, skipping the PSF
step, because step 10 reads star IDs straight from the forced-aperture tables
and does not need it. That shortcut was never re-examined for crowded fields.
M13 and M67 were taken in almost the same seeing (2.72" vs 2.67") yet M13's
zero point scatters three times as much, and swapping the reference system for
Gaia synthetic photometry did not help. What separates them is how close the
stars sit: median FWHM / nearest-neighbour distance is 0.307 against 0.119.

This figure runs the identical pipeline twice on the identical frames, changing
only the photometry source, and shows the three things that move: the CMD, the
dependence of the zero-point residual on neighbour distance, and the number of
stars that survive calibration.

Both panels come from the same raw frames, the same detection, the same WCS,
the same master catalogue and the same step-10 code. The only difference is
whether step 8 ran.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO / "validation" / "paper"))
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

import matplotlib.pyplot as plt  # noqa: E402

from apex_paper_style import (  # noqa: E402
    C, DOUBLE_COL, PALETTE, apply_paper_style, save_fig,
)

PROVENANCE = ("Moravian C3-61000  ·  M13  ·  2026-05-15  ·  "
              "15 frames (5xB, 5xV, 5xR), 60 s  ·  seeing 2.72\"")
CMD_FILE = "median_by_ID_filter_wide_cmd.csv"


def cmd_points(path: Path) -> tuple[np.ndarray, np.ndarray]:
    d = pd.read_csv(path)
    b = pd.to_numeric(d["mag_std_B"], errors="coerce").to_numpy(float)
    v = pd.to_numeric(d["mag_std_V"], errors="coerce").to_numpy(float)
    m = np.isfinite(b) & np.isfinite(v)
    return (b - v)[m], v[m]


def ridge_width(colour: np.ndarray, mag: np.ndarray,
                lo: float, hi: float, deg: int = 3) -> tuple[float, int]:
    """Robust width of the giant branch: MAD about a fitted ridge line.

    A narrower branch means less photometric scatter, since the cluster's
    intrinsic sequence is the same in both panels. The MAD, not the standard
    deviation, so that the binaries and field stars sitting off the sequence
    do not set the number.
    """
    m = np.isfinite(colour) & np.isfinite(mag) & (mag > lo) & (mag < hi)
    if int(m.sum()) < deg + 10:
        return float("nan"), int(m.sum())
    # Fit colour as a function of magnitude: the giant branch is steep in the
    # CMD, so colour(mag) is single-valued where mag(colour) is not.
    coeffs = np.polyfit(mag[m], colour[m], deg)
    for _ in range(5):
        resid = colour - np.polyval(coeffs, mag)
        s = 1.4826 * np.nanmedian(np.abs(resid[m] - np.nanmedian(resid[m])))
        keep = m & (np.abs(resid - np.nanmedian(resid[m])) < 3.0 * s)
        if keep.sum() < deg + 10 or keep.sum() == m.sum():
            break
        m = keep
        coeffs = np.polyfit(mag[m], colour[m], deg)
    resid = colour[m] - np.polyval(coeffs, mag[m])
    return float(1.4826 * np.median(np.abs(resid - np.median(resid)))), int(m.sum())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--result-dir",
                    default=r"E:\APEX_validation\phase3\M13\result")
    ap.add_argument("--aperture-zp", default="cmd_zeropoint_APERTURE")
    ap.add_argument("--psf-zp", default="cmd_zeropoint")
    ap.add_argument("--crowding-json",
                    default=str(REPO / "validation"
                                / "crowding_aperture_vs_psf_M13.json"))
    ap.add_argument("--outdir", default=str(REPO / "validation" / "paper" / "figures"))
    args = ap.parse_args()

    root = Path(args.result_dir)
    crowding = json.loads(Path(args.crowding_json).read_text(encoding="utf-8"))
    bands = {r["label"]: r["bands"] for r in crowding["runs"]}

    apply_paper_style()
    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COL, 2.9))

    # ── (a),(b) the two CMDs ───────────────────────────────────────────────
    widths = {}
    for ax, (label, sub, title) in zip(
        axes[:2],
        (("aperture", args.aperture_zp, "(a) forced aperture"),
         ("psf", args.psf_zp, "(b) PSF (step 8)")),
    ):
        colour, mag = cmd_points(root / sub / CMD_FILE)
        # Identical ink in both panels. A lighter colour on one side would make
        # the comparison look decided before the reader measures anything.
        ax.plot(colour, mag, ".", ms=1.6, mew=0, alpha=0.55,
                color=PALETTE["black"])
        width, n_ridge = ridge_width(colour, mag, 13.0, 16.5)
        widths[label] = (width, n_ridge)
        ax.set_title(title, fontsize=8)
        ax.set_xlabel("$B-V$ (mag)")
        ax.set_xlim(-0.2, 1.8)
        ax.set_ylim(18.2, 11.5)
        ax.text(0.04, 0.05,
                f"N = {colour.size}\nRGB width = {width * 1000:.0f} mmag",
                transform=ax.transAxes, fontsize=6.5, va="bottom")
    axes[0].set_ylabel("$V$ (mag)")
    axes[1].tick_params(labelleft=False)

    # ── (c) the crowding dependence that motivated the change ──────────────
    ax = axes[2]
    width = 0.34
    order = ("B", "V", "R")
    for offset, (label, colour, hatch) in enumerate(
        (("aperture", PALETTE["grey"], ""), ("psf", C["data"], "//"))
    ):
        ratios = [bands[label][b]["ratio"] for b in order]
        ax.bar(np.arange(len(order)) + (offset - 0.5) * width, ratios,
               width, color=colour, hatch=hatch, edgecolor="white",
               linewidth=0.4,
               label="forced aperture" if label == "aperture" else "PSF")
    # Label the reference level in the legend; every position inside the axes
    # collides with one of the bars.
    ax.axhline(1.0, color=PALETTE["black"], lw=0.7, ls=":",
               label="no crowding penalty")
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order)
    ax.set_xlabel("band")
    ax.set_ylabel("crowded / isolated\n|zero-point residual|")
    ax.set_title("(c) crowding penalty", fontsize=8)
    ax.set_ylim(0, 3.9)
    ax.legend(fontsize=6, frameon=False, loc="upper left")

    fig.text(0.5, 0.005, PROVENANCE, ha="center", fontsize=5.6,
             color=PALETTE["grey"])
    fig.tight_layout(rect=(0, 0.035, 1, 1))

    paths = save_fig(fig, "fig_psf_vs_aperture_m13", args.outdir)

    print(f"{'band':>5}{'scatter ap':>12}{'scatter psf':>13}{'N ap':>7}"
          f"{'N psf':>7}{'ratio ap':>10}{'ratio psf':>11}")
    for b in order:
        a, p = bands["aperture"][b], bands["psf"][b]
        print(f"{b:>5}{a['scatter_mmag']:>11.1f}m{p['scatter_mmag']:>12.1f}m"
              f"{a['n_inliers']:>7}{p['n_inliers']:>7}"
              f"{a['ratio']:>10.2f}{p['ratio']:>11.2f}")
    print(f"\nRGB ridge width (13 < V < 16.5): "
          f"aperture {widths['aperture'][0] * 1000:.0f} mmag "
          f"(n={widths['aperture'][1]})  ->  "
          f"PSF {widths['psf'][0] * 1000:.0f} mmag (n={widths['psf'][1]})")
    for kind, path in paths.items():
        print(f"[{kind}] {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
