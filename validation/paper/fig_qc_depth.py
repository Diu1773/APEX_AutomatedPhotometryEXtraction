"""Figure — per-frame detection depth as an operational QC product:
predicted (from each frame's noise + PSF, via the S/N50 law) versus realized
(from real master-catalogue stars), across all reprocessed frames.

Presentation follows the fake/known-source monitoring practice of DES-SN
(Kessler+2015, per-epoch depth distributions) and ATLAS (known-source recovery);
the prediction side — a frame's depth computed BEFORE looking at its detections,
m50 = ZP - 2.5 log10(SN50 * sigma_e / p_peak) with SN50 = 4.05 calibrated by
injection — is the addition of this work, implemented in the pipeline as the
step-7 depth QC gate (apex/analysis/detection_limit.py).

Inputs per frame (no image reads):
  sigma_e  = median per-star local background std (TSV bkg_std_adu) x gain
  p_peak   = median (peak_adu - local bkg)*gain / flux_e over bright clean
             detected stars (TSV joined to the step-4 detect list on det_uid),
             via the canonical estimate_peak_fraction_from_stars
  realized = 50% crossing of detected_flag vs median-frame magnitude
             (data_qc_depth/realized_m50.csv; Eddington-safe, circularity-guarded)

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_qc_depth.py
"""
from __future__ import annotations
import sys
from pathlib import Path
REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from apex.analysis.detection_limit import (
    predict_frame_m50, estimate_peak_fraction_from_stars)
from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

REPRO = Path(r"E:\APEX_validation\reprocess")
REAL_CSV = REPO / "validation" / "paper" / "data_qc_depth" / "realized_m50.csv"
OUTDIR = REPO / "validation" / "paper" / "figures"
TCOL = {"M67": C["data"], "NGC6811": C["reference"], "M13": C["accent"]}
TLBL = {"M67": "M67", "NGC6811": "NGC 6811", "M13": "M13"}
INJECTED = {"pp_Messier67-0008-i.fit", "pp_Messier67-0003-r.fit",
            "pp_Messier67-0004-g.fit", "pp_messier13-0004-R.fit",
            "pp_messier13-0001-V.fit", "pp_NGC6811-0005-R.fit",
            "pp_NGC6811-0008-R.fit"}


def predict_for_frame(tgt: str, fname: str) -> float:
    base = REPRO / tgt / "result"
    tsv = base / "step7_forced_phot" / f"photometry_{fname}.tsv"
    det = base / "step4_detection" / f"detect_{fname}.csv"
    if not tsv.exists() or not det.exists():
        return float("nan")
    d = pd.read_csv(tsv, sep="\t",
                    usecols=["det_uid", "flux_e", "bkg_median_adu", "bkg_std_adu",
                             "gain_e_per_adu", "snr", "detected_flag",
                             "bad_phot_flag", "is_saturated", "off_frame_flag"])
    gain = float(pd.to_numeric(d.gain_e_per_adu, errors="coerce").median())
    clean = d[(d.off_frame_flag == False) & (d.bad_phot_flag == False)]
    sigma_e = float(pd.to_numeric(clean.bkg_std_adu, errors="coerce").median()) * gain
    # join detected stars to the detect list for their peak pixel
    dd = pd.read_csv(det, usecols=["det_uid", "peak_adu"])
    j = clean[clean.detected_flag == True].merge(dd, on="det_uid", how="inner")
    peak_e = (pd.to_numeric(j.peak_adu, errors="coerce")
              - pd.to_numeric(j.bkg_median_adu, errors="coerce")).to_numpy(float) * gain
    flux_e = pd.to_numeric(j.flux_e, errors="coerce").to_numpy(float)
    bright = ((pd.to_numeric(j.snr, errors="coerce") > 20)
              & (j.is_saturated == False)).to_numpy(bool)
    p_peak, n_used = estimate_peak_fraction_from_stars(peak_e, flux_e, mask=bright)
    if not np.isfinite(p_peak) or n_used < 10:
        return float("nan")
    return predict_frame_m50(sigma_e, p_peak)


def main() -> int:
    r = pd.read_csv(REAL_CSV)
    r["predicted_m50"] = [predict_for_frame(t, f) for t, f in zip(r.target, r.file)]
    ok = np.isfinite(r.predicted_m50) & np.isfinite(r.realized_m50)
    val = r[ok & r.depth_valid]
    inv = r[ok & ~r.depth_valid]
    resid = (val.realized_m50 - val.predicted_m50).to_numpy(float)
    rms = float(np.sqrt(np.mean(resid ** 2)))
    print(f"[qc-depth] frames={ok.sum()} valid={len(val)}  "
          f"realized-predicted: mean {resid.mean():+.3f}  RMS {rms:.3f} mag")

    fig = plt.figure(figsize=(DOUBLE_COL, 3.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.25], wspace=0.27)

    # ── (a) distribution of realized depths (Kessler Fig.7 presentation) ──
    axa = fig.add_subplot(gs[0])
    bins = np.arange(13.9, 18.2, 0.35)
    bottom = np.zeros(len(bins) - 1)
    for tgt in ["M67", "NGC6811", "M13"]:
        h, _ = np.histogram(val[val.target == tgt].realized_m50, bins=bins)
        axa.bar(0.5 * (bins[:-1] + bins[1:]), h, width=0.33, bottom=bottom,
                color=TCOL[tgt], alpha=0.85, label=TLBL[tgt])
        bottom += h
    axa.set_xlabel("realized $m_{50}$ (total-$e^-$ scale, ZP = 25)")
    axa.set_ylabel("frames / bin")
    axa.legend(loc="upper left", fontsize=7.0)
    axa.set_title(f"(a) every frame's depth, monitored", loc="left", fontsize=9.5)
    axa.text(0.97, 0.95, "each entry =\none single exposure",
             transform=axa.transAxes, ha="right", va="top",
             fontsize=6.8, color="0.3")

    # ── (b) predicted (before detection) vs realized (real stars) ──
    axb = fig.add_subplot(gs[1])
    lo, hi = 13.8, 18.3
    axb.fill_between([lo, hi], [lo - 0.5, hi - 0.5], [lo + 0.5, hi + 0.5],
                     color=PALETTE["grey"], alpha=0.18, zorder=1, lw=0)
    axb.plot([lo, hi], [lo, hi], "-", color=PALETTE["grey"], lw=1.0, zorder=2)
    for tgt in ["M67", "NGC6811", "M13"]:
        s = val[val.target == tgt]
        injm = s.file.isin(INJECTED)
        axb.plot(s[~injm].predicted_m50, s[~injm].realized_m50, "o",
                 color=TCOL[tgt], ms=4.2, mew=0, alpha=0.85, zorder=4,
                 label=TLBL[tgt])
        axb.plot(s[injm].predicted_m50, s[injm].realized_m50, "o",
                 color=TCOL[tgt], ms=5.2, mec="k", mew=1.0, zorder=5)
    axb.plot(inv.predicted_m50, inv.realized_m50, "o", color="none",
             mec="0.6", ms=4.2, mew=0.9, zorder=3,
             label="excluded (master-limit circularity)")
    axb.plot([], [], "o", color="0.4", ms=5.2, mec="k", mew=1.0,
             label="injection-calibrated frame")
    axb.set_xlim(lo, hi); axb.set_ylim(lo, hi)
    axb.set_aspect("equal", adjustable="box")
    axb.set_xlabel("predicted $m_{50}$ (noise + PSF only)")
    axb.set_ylabel("realized $m_{50}$ (real stars)")
    axb.legend(loc="lower right", fontsize=6.6, framealpha=0.92)
    axb.set_title("(b) depth is predicted before detection", loc="left",
                  fontsize=9.5)
    axb.text(0.04, 0.96,
             f"RMS = {rms*1000:.0f} mmag ({len(val)} frames)\n"
             "band: QC-gate tolerance ±0.5 mag",
             transform=axb.transAxes, va="top", ha="left", fontsize=6.9,
             color="0.25")

    paths = save_fig(fig, "fig_qc_depth", OUTDIR)
    plt.close(fig)
    r.to_csv(REPO / "validation/paper/data_qc_depth/qc_depth_summary.csv", index=False)
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
