# -*- coding: utf-8 -*-
"""Combine the synthetic SEP and real-frame IRAF photometry cross-checks."""
from __future__ import annotations

import json
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
from fig4_crosscheck_sep import (
    FWHM_PX, FRAME_SIZE, GAIN, MAG_MAX, MAG_MIN, N_STARS, READ_NOISE_E,
    R_AP, SEED, SKY_ADU, SNR_MIN, _mad, _rms, build_synthetic_frame,
    measure_apex, measure_sep,
)

apply_paper_style()

OUTDIR = PAPER / "figures"
CAPDIR = PAPER / "captions"
_RUN_REL = Path("benchmark") / "runs" / "ngc6811_iraf_allapex_v1" / "phot_fixed_coords"
_RUN_CANDIDATES = [
    REPO / _RUN_REL,
    Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction") / _RUN_REL,
]
RUN_DIR = next(
    (path for path in _RUN_CANDIDATES if (path / "fixed_comparison.csv").exists()),
    _RUN_CANDIDATES[0],
)


def _running_median(x: np.ndarray, y: np.ndarray, width: float = 0.6):
    edges = np.arange(np.floor(x.min()), np.ceil(x.max()) + width, width)
    bx, by = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        keep = (x >= lo) & (x < hi)
        if keep.sum() >= 5:
            bx.append(0.5 * (lo + hi))
            by.append(float(np.median(y[keep])))
    return np.asarray(bx), np.asarray(by)


def main() -> int:
    image, x, y, true_mag = build_synthetic_frame()
    apex_flux, apex_snr = measure_apex(image, x, y)
    sep_flux, sep_flag = measure_sep(image, x, y)
    with np.errstate(divide="ignore", invalid="ignore"):
        mag_apex = -2.5 * np.log10(np.where(apex_flux > 0, apex_flux, np.nan))
        mag_sep = -2.5 * np.log10(np.where(sep_flux > 0, sep_flux, np.nan))
    keep = (
        np.isfinite(mag_apex) & np.isfinite(mag_sep) & np.isfinite(apex_snr)
        & (apex_snr > SNR_MIN) & (sep_flag == 0)
    )
    ma, ms, mt = mag_apex[keep], mag_sep[keep], true_mag[keep]
    sep_offset = float(np.median(ma - ms))
    sep_delta = (ma - ms) - sep_offset
    sep_mad = _mad(sep_delta)
    sep_rms = _rms(sep_delta)
    sep_r = float(np.corrcoef(ma, ms)[0, 1])

    iraf = pd.read_csv(RUN_DIR / "fixed_comparison.csv")
    summary = json.loads((RUN_DIR / "fixed_summary.json").read_text(encoding="utf-8"))
    iraf_mag = pd.to_numeric(iraf["iraf_mag"], errors="coerce").to_numpy(float)
    apex_iraf = pd.to_numeric(iraf["apex_mag_iraf_units"], errors="coerce").to_numpy(float)
    iraf_delta = pd.to_numeric(
        iraf["delta_mag_units_centered"], errors="coerce"
    ).to_numpy(float)
    finite = np.isfinite(iraf_mag) & np.isfinite(apex_iraf) & np.isfinite(iraf_delta)
    iraf_mag, apex_iraf, iraf_delta = (
        iraf_mag[finite], apex_iraf[finite], iraf_delta[finite]
    )
    iraf_mad = _mad(iraf_delta)
    iraf_rms = _rms(iraf_delta)
    iraf_r = float(np.corrcoef(apex_iraf, iraf_mag)[0, 1])

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.05))

    axa.axhspan(-sep_mad, sep_mad, color=PALETTE["skyblue"], alpha=0.45, lw=0)
    axa.axhline(0, color=PALETTE["grey"], ls="--", lw=0.9)
    axa.scatter(mt, sep_delta, s=13, facecolors="white", edgecolors=C["data"],
                linewidths=0.6, alpha=0.8)
    bx, by = _running_median(mt, sep_delta)
    axa.plot(bx, by, color=C["model"], ls="-", lw=1.6, label="binned median")
    axa.set_xlabel("true magnitude")
    axa.set_ylabel(r"$(m_{\rm APEX}-m_{\rm SEP})-\mathrm{median}$ (mag)")
    axa.set_title(f"(a) SEP: MAD {sep_mad * 1000:.1f} mmag", loc="left")
    axa.text(
        0.03, 0.95,
        f"N = {len(mt)}\nRMS = {sep_rms * 1000:.1f} mmag\nr = {sep_r:.5f}",
        transform=axa.transAxes, va="top", fontsize=7.0,
        bbox=dict(fc="white", ec=PALETTE["grey"], alpha=0.9, pad=2),
    )

    shown = np.abs(iraf_delta) <= 0.16
    axb.axhline(0, color=PALETTE["grey"], ls="--", lw=0.9)
    axb.scatter(iraf_mag[shown], iraf_delta[shown], s=11, color=C["data"],
                edgecolors="none", alpha=0.45)
    bx, by = _running_median(iraf_mag, iraf_delta)
    axb.plot(bx, by, color=C["model"], ls="-", lw=1.6, label="binned median")
    axb.set_xlabel("IRAF magnitude")
    axb.set_ylabel(r"$m_{\rm APEX}-m_{\rm IRAF}$ (mag)")
    axb.set_ylim(-0.16, 0.16)
    axb.set_title(f"(b) IRAF: MAD {iraf_mad * 1000:.1f} mmag", loc="left")
    axb.text(
        0.03, 0.95,
        f"N = {len(iraf_mag)}\nRMS = {iraf_rms * 1000:.1f} mmag\nr = {iraf_r:.5f}",
        transform=axb.transAxes, va="top", fontsize=7.0,
        bbox=dict(fc="white", ec=PALETTE["grey"], alpha=0.9, pad=2),
    )
    n_hidden = int((~shown).sum())
    if n_hidden:
        axb.text(0.98, 0.04, f"{n_hidden} points outside shown range",
                 transform=axb.transAxes, ha="right", fontsize=6.0,
                 color=PALETTE["grey"])
    for axis in (axa, axb):
        axis.legend(loc="lower left", fontsize=6.8)

    provenance = (
        f"(a) synthetic frame: seed {SEED}, {FRAME_SIZE}$^2$ px, {N_STARS} isolated stars, "
        f"FWHM {FWHM_PX:g} px, sky {SKY_ADU:.0f} ADU, gain {GAIN:g} e$^-$/ADU, "
        f"RN {READ_NOISE_E:g} e$^-$; same known centres and $r_{{ap}}={R_AP:g}$ px.\n"
        "(b) real NGC 6811 $V$, 30 s, Moravian C3-61000 (2$\\times$2), "
        "2026-06-11; 499 fixed centres; APEX-reduced frame; matched apertures; "
        "one zeropoint alignment."
    )
    fig.tight_layout(rect=(0, 0.16, 1, 1), w_pad=1.0)
    fig.text(0.005, 0.015, provenance, fontsize=5.5, color=PALETTE["grey"],
             va="bottom")
    paths = save_fig(fig, "fig_photometry_crosschecks", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    formal = float(summary["median_quadrature_formal_err"])
    (CAPDIR / "fig_photometry_crosschecks.md").write_text(
        f"""# Figure — independent photometry-engine cross-checks

**(a)** APEX and SEP aperture photometry at the same known centres in a
synthetic frame. After one zeropoint alignment, the {len(mt)} retained stars
have MAD {sep_mad * 1000:.1f} mmag, RMS {sep_rms * 1000:.1f} mmag, and
$r={sep_r:.5f}$. **(b)** APEX and IRAF/DAOPHOT measurements at 499 fixed
centres in one APEX-reduced NGC 6811 $V$ exposure. The residual has MAD
{iraf_mad * 1000:.1f} mmag and RMS {iraf_rms * 1000:.1f} mmag; the median
quadrature formal error is {formal * 1000:.1f} mmag. The lines show binned
medians, not fitted corrections. Synthetic settings and real-frame provenance
are printed below the panels.
""",
        encoding="utf-8",
    )
    print(f"SEP: N={len(mt)} MAD={sep_mad*1000:.1f} mmag RMS={sep_rms*1000:.1f}")
    print(f"IRAF: N={len(iraf_mag)} MAD={iraf_mad*1000:.1f} mmag RMS={iraf_rms*1000:.1f}")
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
