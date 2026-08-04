"""Report figures for the U-B degeneracy / zero-point re-anchor result.

Produces (next to REPORT_UB_DEGENERACY.md):
  fig_ub_zeropoint_diagnosis.png/pdf — per-band offsets vs the MMJ93 standards
  fig_ub_rail_resolution.png/pdf     — cluster parameters across the four fits

Run:
    .venv-deploy/Scripts/python.exe -X utf8 validation/psf_crossinstrument/fig_ub_rail.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "paper"))

from apex_paper_style import C, DOUBLE_COL, apply_paper_style, save_fig  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

BASE = Path(r"E:/APEX_validation/psf_crossinstrument/m67_ubv")
PROV = ("M67 | LCO 1m elp/kb74 SBIG STL-6303, 2015-04-01 (U 2x300 s, B 3x70 s, V 3x30 s) | "
        "standards: Montgomery+1993 (VizieR J/AJ/106/181, Landolt-tied)")

MAD = 1.4826


def fig_zeropoint_diagnosis() -> None:
    res = pd.read_csv(BASE / "external_ref" / "mmj93_match_residuals.csv")
    bands = [("dU", "U"), ("dB", "B"), ("dV", "V")]

    fig, axes = plt.subplots(3, 1, figsize=(DOUBLE_COL, 6.2), sharex=True)
    for ax, (col, band) in zip(axes, bands):
        d = res[["BV_mmj", col]].dropna()
        med = float(np.median(d[col]))
        sig = float(MAD * np.median(np.abs(d[col] - med)))
        ax.axhline(0.0, color="0.25", lw=0.8)
        ax.axhspan(med - sig, med + sig, color=C["data"], alpha=0.12, lw=0)
        ax.scatter(d["BV_mmj"], d[col], s=9, color=C["data"], alpha=0.55,
                   edgecolors="none", label="matched star (measured)")
        ax.axhline(med, color=C["model"], lw=1.4,
                   label=f"median {med:+.3f} mag")
        gaia_sig = {"U": 0.200, "B": 0.025, "V": 0.030}[band]
        ax.annotate(
            f"$\\Delta${band} = {med:+.3f} mag  (robust $\\sigma$ {sig:.3f}, "
            f"N={len(d)});  Gaia ref quoted $\\sigma$ = {gaia_sig:.3f}",
            xy=(0.02, 0.92), xycoords="axes fraction", va="top", fontsize=8)
        ax.set_ylabel(f"{band}(std) $-$ {band}(APEX)  [mag]")
        ax.set_ylim(med - 0.45, med + 0.45)
        if band == "U":
            ax.legend(loc="lower left", fontsize=7, framealpha=0.9)
    axes[-1].set_xlabel("(B$-$V) standard  [mag]")
    axes[0].set_title("Zero-point offsets vs external standards — the Gaia 'approx' "
                      "U reference is off by $-$0.13 mag", fontsize=9)
    fig.text(0.01, 0.005, PROV, fontsize=6, color="0.35")
    fig.tight_layout(rect=(0, 0.02, 1, 1))
    save_fig(fig, "fig_ub_zeropoint_diagnosis", HERE)


FITS = [
    ("fit_BV_only_ustd.json", "result_ustd", "B$-$V", "floor"),
    ("fit_UB_BV.json", "result", "+U$-$B\n(Gaia U)", "bad"),
    ("fit_UB_BV_ustd.json", "result_ustd", "+U$-$B\n(std U)", "data"),
    ("fit_UB_BV_psf_ustd.json", "result_psf_ustd", "+U$-$B\n(std U, PSF)", "reference"),
]

PANELS = [
    ("age_gyr", "age  [Gyr]", (3.5, 4.5)),          # VandenBerg+ / Bellini+: 3.5-4.5
    ("metallicity", "[M/H]  [dex]", (-0.05, 0.10)),  # [Fe/H] ~ 0.0..+0.05
    ("e_bv", "E(B$-$V)  [mag]", (0.03, 0.05)),
    ("distance_mod", "(m$-$M)$_0$  [mag]", (9.6, 9.7)),
]


def fig_rail_resolution() -> None:
    rows = []
    for fname, resdir, label, ckey in FITS:
        d = json.loads((BASE / resdir / "cmd_isochrone" / fname).read_text(encoding="utf-8"))
        rows.append((label, ckey, d["summary"]))

    fig, axes = plt.subplots(1, 4, figsize=(DOUBLE_COL, 2.9))
    x = np.arange(len(rows))
    for ax, (key, ylabel, lit) in zip(axes, PANELS):
        ax.axhspan(lit[0], lit[1], color="0.82", alpha=0.6, lw=0,
                   label="literature")
        for i, (label, ckey, s) in enumerate(rows):
            p16, p50, p84 = s[key]
            ax.errorbar(i, p50, yerr=[[p50 - p16], [p84 - p50]],
                        fmt="o", ms=4.5, color=C[ckey], capsize=2.5, lw=1.2)
        ax.set_xticks(x)
        ax.set_xticklabels([r[0] for r in rows], fontsize=5.6,
                           rotation=28, ha="right")
        ax.set_ylabel(ylabel, fontsize=8)
        ax.tick_params(axis="y", labelsize=7)
        ax.set_xlim(-0.6, len(rows) - 0.4)
    axes[3].legend(loc="lower left", fontsize=6, framealpha=0.9)
    fig.suptitle("Same data, same stars, same MCMC — only the zero-point anchor changes",
                 fontsize=9, y=1.02)
    fig.text(0.01, 0.005,
             PROV + " | membership + Gaia-parallax dm prior, 32 walkers x 2000 steps, seed 2024",
             fontsize=5.4, color="0.35")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    save_fig(fig, "fig_ub_rail_resolution", HERE)


if __name__ == "__main__":
    apply_paper_style()
    fig_zeropoint_diagnosis()
    fig_rail_resolution()
    print("done")
