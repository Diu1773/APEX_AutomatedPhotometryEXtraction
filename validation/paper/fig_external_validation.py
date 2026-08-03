# -*- coding: utf-8 -*-
"""External-catalog residuals and the resulting Johnson CMD in one figure."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).absolute().parents[2]
PAPER = REPO / "validation" / "paper"
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(PAPER))

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL
from fig7_reference_crosscheck import (
    BRIGHT_ANCHOR, FAINT, binned, faint_drift, load_matched, residual_vs_mag,
)
from fig8_cmd_reproduction import CAL_CSV, CMD_CSV, ridgeline

apply_paper_style()

OUTDIR = PAPER / "figures"
CAPDIR = PAPER / "captions"


def main() -> int:
    matched = load_matched()
    gr = (matched["ps1_g"] - matched["ps1_r"]).to_numpy(float)
    edges = np.arange(14.0, 18.75, 0.5)
    residual_data = {}
    drifts = {}
    for band in ("B", "V"):
        refmag = matched[f"ref_{band}"]
        snr = matched[f"snr_{band}"].to_numpy(float)
        apex_r, apex_good = residual_vs_mag(
            matched[f"mag_inst_{band}"], matched["ps1_g"], gr, snr, refmag
        )
        gaia_r, gaia_good = residual_vs_mag(
            matched[f"ref_{band}"], matched["ps1_g"], gr, None, refmag
        )
        residual_data[band] = (refmag, apex_r, apex_good, gaia_r, gaia_good)
        drifts[f"{band}_apex"] = faint_drift(apex_r, apex_good, refmag)
        drifts[f"{band}_gaia"] = faint_drift(gaia_r, gaia_good, refmag)

    cmd = pd.read_csv(CMD_CSV)
    cal = pd.read_csv(CAL_CSV)[["ID", "ref_B", "ref_V"]]
    cmd = cmd.merge(cal, on="ID", how="left")
    apex_v = pd.to_numeric(cmd["mag_std_V"], errors="coerce").to_numpy(float)
    apex_bv = (
        pd.to_numeric(cmd["mag_std_B"], errors="coerce")
        - pd.to_numeric(cmd["mag_std_V"], errors="coerce")
    ).to_numpy(float)
    gaia_v = pd.to_numeric(cmd["ref_V"], errors="coerce").to_numpy(float)
    gaia_bv = (
        pd.to_numeric(cmd["ref_B"], errors="coerce")
        - pd.to_numeric(cmd["ref_V"], errors="coerce")
    ).to_numpy(float)
    ridge_edges = np.arange(12.8, 17.8, 0.4)
    apex_rl, apex_rm = ridgeline(apex_bv, apex_v, ridge_edges)
    gaia_rl, gaia_rm = ridgeline(gaia_bv, gaia_v, ridge_edges)
    common = np.intersect1d(apex_rm, gaia_rm)
    ai, gi = np.isin(apex_rm, common), np.isin(gaia_rm, common)
    ridge_rms = float(np.sqrt(np.nanmean((apex_rl[ai] - gaia_rl[gi]) ** 2)))
    n_cmd = int(np.isfinite(apex_v + apex_bv + gaia_v + gaia_bv).sum())

    fig, axes = plt.subplots(
        1, 3, figsize=(DOUBLE_COL, 3.05),
        gridspec_kw={"width_ratios": [1.0, 1.0, 1.08]},
    )
    for ax, band, letter in zip(axes[:2], ("B", "V"), ("a", "b")):
        refmag, apex_r, apex_good, gaia_r, gaia_good = residual_data[band]
        ax.axhline(0, color=PALETTE["grey"], ls=":", lw=0.8)
        ax.axvspan(FAINT, edges[-1], color=PALETTE["skyblue"], alpha=0.35, lw=0)
        for values, good, marker, linestyle, label, key in (
            (apex_r, apex_good, "o", "-", "APEX − PS1", f"{band}_apex"),
            (gaia_r, gaia_good, "s", "--", "Gaia ref − PS1", f"{band}_gaia"),
        ):
            xx, yy, ee = binned(values, good, refmag, edges)
            ax.errorbar(
                xx, yy, yerr=ee, fmt=marker, ls=linestyle,
                mfc="white" if marker == "s" else C["data"],
                mec=C["data"], color=C["data"], ms=4.0, lw=1.2, capsize=2,
                label=f"{label} ({drifts[key]:+.3f})",
            )
        ax.set_xlabel(f"${band}$ magnitude")
        ax.set_title(f"({letter}) {band}: faint-end drift", loc="left", fontsize=8.4)
        ax.set_ylim(-0.045, 0.055)
        ax.legend(loc="upper left", fontsize=6.2)
    axes[0].set_ylabel(r"residual vs PS1 (mag; bright anchor = 0)")

    axc = axes[2]
    axc.scatter(apex_bv, apex_v, s=5, color=C["data"], edgecolors="none",
                alpha=0.30, label="APEX")
    axc.scatter(gaia_bv, gaia_v, s=7, facecolors="none", edgecolors=C["model"],
                linewidths=0.35, alpha=0.32, label="Gaia-transformed")
    axc.plot(apex_rl, apex_rm, color=C["data"], ls="-", lw=1.7,
             label="APEX ridge")
    axc.plot(gaia_rl, gaia_rm, color=C["model"], ls="--", lw=1.7,
             label="Gaia ridge")
    axc.set_xlim(-0.1, 1.8)
    axc.set_ylim(18.4, 11.0)
    axc.set_xlabel(r"$B-V$")
    axc.set_ylabel(r"$V$")
    axc.set_title(f"(c) CMD ridge RMS {ridge_rms * 1000:.0f} mmag", loc="left",
                  fontsize=8.4)
    axc.legend(loc="upper right", fontsize=5.9, ncol=1)
    axc.text(0.04, 0.04, f"N = {n_cmd}", transform=axc.transAxes, fontsize=6.5)

    provenance = (
        f"NGC 6811 · (a,b) {len(matched)} APEX/Gaia/PS1 matches; PS1 $g,r$; "
        f"colour-term fit on bright stars; {BRIGHT_ANCHOR[0]}–{BRIGHT_ANCHOR[1]} mag anchor. "
        f"(c) {n_cmd} common Johnson $B,V$ stars; APEX ground-based vs Gaia-transformed.\n"
        "APEX source: median_by_ID_filter_wide_cmd.csv · Gaia reference: "
        "gaia_sdss_calibrator_by_ID.csv · PS1 cache: ps1_match_ngc6811.csv."
    )
    fig.tight_layout(rect=(0, 0.16, 1, 1), w_pad=0.8)
    fig.text(0.005, 0.015, provenance, fontsize=5.4, color=PALETTE["grey"],
             va="bottom")
    paths = save_fig(fig, "fig_external_validation", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    (CAPDIR / "fig_external_validation.md").write_text(
        f"""# Figure — external-catalog and CMD validation

NGC 6811 measurements compared with Pan-STARRS 1 (PS1) and a
Gaia-transformed Johnson reference. **(a,b)** Binned residuals after a
colour-term fit and a {BRIGHT_ANCHOR[0]}–{BRIGHT_ANCHOR[1]} mag anchor. In
$B$, the faint-end changes are {drifts['B_apex']:+.3f} mag for APEX−PS1 and
{drifts['B_gaia']:+.3f} mag for Gaia-reference−PS1; in $V$ they are
{drifts['V_apex']:+.3f} and {drifts['V_gaia']:+.3f} mag, respectively.
**(c)** Johnson CMDs from APEX and the Gaia-transformed reference for {n_cmd}
common stars. Their clipped main-sequence ridgelines differ by an RMS of
{ridge_rms * 1000:.0f} mmag. The CMD comparison tests the plotted photometric
product; it does not validate isochrone parameters.
""",
        encoding="utf-8",
    )
    print(f"drifts: {drifts}; CMD N={n_cmd}; ridge RMS={ridge_rms*1000:.1f} mmag")
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
