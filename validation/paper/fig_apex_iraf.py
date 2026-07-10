"""Figure (all-APEX) — Independent photometry validation: APEX vs IRAF/DAOPHOT.

The honest, non-circular accuracy check for the fully-APEX (raw->science) NGC 6811
reduction. IRAF/DAOPHOT (via PyRAF) re-measures the SAME stars at the SAME fixed
sky coordinates on the SAME all-APEX-calibrated frame. IRAF shares no code with
APEX, so a milli-magnitude differential agreement confirms APEX's photometry is
correct on-sky — not merely self-consistent, and not a comparison against the
author's own preprocessing.

Data produced by:
    apex.benchmark.iraf_crosscheck_cli --input <all-APEX NGC6811 V frame>
        --step7 <its Step-7 TSV> --mode phot_fixed_coords --runtime-cmd wsl python3

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_apex_iraf.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

RUN = REPO / "benchmark" / "runs" / "ngc6811_iraf_matched_v1"
CSV = RUN / "phot_fixed_coords" / "fixed_comparison.csv"
MANIFEST = RUN / "iraf_crosscheck_manifest.json"
OUTDIR = REPO / "validation" / "paper" / "figures"

OBJ = "NGC 6811"
DATE = "2026-06-11"  # all-APEX reprocess night (E:\observe_raw_Analysis\20260611)


def _mad(v):
    v = v[np.isfinite(v)]
    return float(1.4826 * np.median(np.abs(v - np.median(v)))) if v.size else float("nan")


def _rms(v):
    v = v[np.isfinite(v)]
    return float(np.sqrt(np.mean(v * v))) if v.size else float("nan")


def main() -> int:
    df = pd.read_csv(CSV)
    man = json.loads(MANIFEST.read_text(encoding="utf-8"))
    meta = man.get("frame_metadata", {})
    cfg = man.get("effective_config", {})
    frame_name = Path(str(man.get("input_fits", cfg.get("input_fits", "")))).name
    band = str(meta.get("filter_name", cfg.get("filter_name", "V"))).strip() or "V"
    fwhm = float(meta.get("fwhm_px", float("nan")))
    zmag = float(cfg.get("zmag", 25.0))

    apex_mag = pd.to_numeric(df["apex_mag_iraf_units"], errors="coerce").to_numpy(float)
    iraf_mag = pd.to_numeric(df["iraf_mag"], errors="coerce").to_numpy(float)
    delta = pd.to_numeric(df["delta_mag_units_centered"], errors="coerce").to_numpy(float)
    dr_px = pd.to_numeric(df.get("iraf_dr_px"), errors="coerce").to_numpy(float)

    keep = np.isfinite(apex_mag) & np.isfinite(iraf_mag) & np.isfinite(delta)
    apex_mag, iraf_mag, delta = apex_mag[keep], iraf_mag[keep], delta[keep]
    dr_px = dr_px[keep] if dr_px.size == keep.size else dr_px
    n = int(keep.sum())

    offset = float(np.median(apex_mag - iraf_mag))
    apex_aligned = apex_mag - offset
    med_delta, mad, rms = float(np.median(delta)), _mad(delta), _rms(delta)
    pearson_r = float(np.corrcoef(apex_aligned, iraf_mag)[0, 1])
    pos_rms = _rms(dr_px) if np.isfinite(dr_px).any() else float("nan")

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.35))

    lo = float(min(apex_aligned.min(), iraf_mag.min()))
    hi = float(max(apex_aligned.max(), iraf_mag.max()))
    pad = 0.04 * (hi - lo)
    ax_a.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=C["floor"], lw=1.0,
              zorder=1, label="$y=x$")
    ax_a.scatter(iraf_mag, apex_aligned, s=12, alpha=0.7, color=C["reference"],
                 edgecolors="none", zorder=2)
    ax_a.set_xlim(lo - pad, hi + pad); ax_a.set_ylim(lo - pad, hi + pad)
    ax_a.set_aspect("equal", adjustable="box")
    ax_a.set_xlabel(r"$m_{\mathrm{IRAF}}$  (instrumental, zmag$=%.0f$)" % zmag)
    ax_a.set_ylabel(r"$m_{\mathrm{APEX}}$  (ZP-aligned to IRAF)")
    ax_a.text(0.05, 0.95, f"$r = {pearson_r:.5f}$\n$N = {n}$", transform=ax_a.transAxes,
              va="top", ha="left", fontsize=8.0,
              bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85,
                    "edgecolor": PALETTE["grey"]})
    ax_a.legend(loc="lower right", fontsize=8.0)
    ax_a.set_title("(a) all-APEX vs IRAF", loc="left")

    ax_b.axhspan(med_delta - mad, med_delta + mad, color=C["model"], alpha=0.18,
                 zorder=0, label=r"$\pm$MAD")
    ax_b.axhline(med_delta, color=C["model"], lw=1.2, zorder=2, label="median")
    ax_b.axhline(0.0, color=C["floor"], lw=0.8, ls=":", zorder=1)
    ax_b.scatter(iraf_mag, delta, s=12, alpha=0.7, color=C["reference"],
                 edgecolors="none", zorder=3)
    ylim = max(0.02, 1.15 * float(np.nanpercentile(np.abs(delta), 99)))
    n_clip = int(np.sum(np.abs(delta) > ylim))
    ax_b.set_ylim(-ylim, ylim)
    if n_clip:
        ax_b.text(0.95, 0.05, f"{n_clip} outlier(s) beyond axis", transform=ax_b.transAxes,
                  va="bottom", ha="right", fontsize=7.0, color=PALETTE["grey"])
    ax_b.set_xlabel(r"$m_{\mathrm{IRAF}}$  (instrumental)")
    ax_b.set_ylabel(r"$\Delta = (m_{\mathrm{APEX}} - m_{\mathrm{IRAF}}) - \mathrm{ZP}$  (mag)")
    ax_b.text(0.05, 0.95, f"offset $= {offset:+.4f}$\nMAD $= {mad:.4f}$\nRMS $= {rms:.4f}$",
              transform=ax_b.transAxes, va="top", ha="left", fontsize=8.0,
              bbox={"boxstyle": "round,pad=0.3", "facecolor": "white", "alpha": 0.85,
                    "edgecolor": PALETTE["grey"]})
    ax_b.legend(loc="lower left", fontsize=7.5, ncol=2)
    ax_b.set_title("(b) ZP-aligned residuals", loc="left")

    provenance = (
        f"{OBJ}  ·  {band}-band  ·  {frame_name}  ({DATE})  ·  fully-APEX reduction  ·  "
        f"N = {n} stars  ·  FWHM = {fwhm:.1f}px  ·  "
        f"APEX forced-aperture vs IRAF phot at fixed coordinates (PyRAF)"
    )
    fig.suptitle(provenance, fontsize=7.2, y=1.005, color="#333333")
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    paths = save_fig(fig, "fig_apex_iraf", OUTDIR)
    plt.close(fig)
    print(f"N={n}  offset={offset:+.4f}  MAD={mad:.4f}  RMS={rms:.4f}  r={pearson_r:.5f}  posRMS={pos_rms:.4f}px")
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
