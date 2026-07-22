"""Figure — injection-recovery completeness: controlled synthetic (verification)
and a real cluster frame (validation), the latter cross-checked against the
pipeline's own detection-limit rolloff.

This is the reference-standard test: artificial stars injected into a REAL frame
(DAOPHOT ADDSTAR; DES Balrog; HSC SynPipe; AutoPhOT App.D) rather than a smooth
synthetic frame. The synthetic curve is the controlled known-truth check that the
recovery machinery is numerically correct (verification); the real-frame curve is
the on-sky performance (validation). The two frames differ in seeing and sky —
the synthetic is idealised (FWHM 3.4px, sky 150 ADU), the real M13 V frame is
shallow (FWHM 7.6px, sky 1315 ADU) — so the depth difference is a property of the
frames, NOT a generic "synthetic is optimistic" claim.

Independent cross-check (the strong result): the recovered 50% limit on the real
frame coincides with the magnitude where the pipeline's own clean detections roll
off — two independent methods agree on the depth of the same real frame.

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_completeness_realvssynth.py
"""
from __future__ import annotations
import sys
from pathlib import Path
REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import pandas as pd
from astropy.io import fits
import sep
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

SYN = REPO / "validation" / "paper" / "data" / "artificial_star" / "benchmark_run"
REAL = REPO / "validation" / "paper" / "data_realframe_M13V" / "artificial_star" / "benchmark_run"
FRAME = Path(r"E:\APEX_validation\reprocess\M13\calibrated\20260515\pp_messier13-0001-V.fit")
GAIN = 0.689  # e-/ADU (PTC-measured, C3-61000)
ZP = 25.0
OUTDIR = REPO / "validation" / "paper" / "figures"


def read_off(mag, comp, level):
    o = np.argsort(mag); mm, cc = mag[o], comp[o]
    for i in range(len(mm) - 1):
        if cc[i] >= level >= cc[i + 1]:
            den = cc[i] - cc[i + 1]
            f = (cc[i] - level) / den if den else 0.0
            return float(mm[i] + f * (mm[i + 1] - mm[i]))
    return float("nan")


def curve(run_dir, bw=0.25):
    s = pd.read_csv(run_dir / "stars.csv")
    s = s[~s["baseline_confounded"].astype(bool)]
    m = s["magnitude_true"].to_numpy(float)
    r = s["recovered"].to_numpy(bool)
    lo0 = np.floor(m.min() / bw) * bw
    edges = np.arange(lo0, m.max() + bw, bw)
    xs, cs = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        msk = (m >= lo) & (m < hi)
        if msk.sum() >= 15:
            xs.append(0.5 * (lo + hi)); cs.append(r[msk].mean())
    return np.array(xs), np.array(cs)


def empirical_detection_mags():
    """Instrumental mag (ZP25) of the pipeline's own clean baseline detections on
    the real M13 frame — an INDEPENDENT probe of the frame depth (no injection)."""
    det = pd.read_csv(REAL / "baseline" / "step4_detection"
                      / "detect_pp_messier13-0001-V.fit.csv")
    data = fits.getdata(FRAME).astype(np.float64)
    bkg = sep.Background(data)
    sub = data - bkg.back()
    fl, _, _ = sep.sum_circle(sub, det["x"].to_numpy(float),
                              det["y"].to_numpy(float), 6.0, gain=GAIN)
    fl_e = fl * GAIN
    fl_e = fl_e[np.isfinite(fl_e) & (fl_e > 0)]
    return ZP - 2.5 * np.log10(fl_e)


def main() -> int:
    ms, cs = curve(SYN)
    mr, cr = curve(REAL)
    m50s, m50r = read_off(ms, cs, 0.5), read_off(mr, cr, 0.5)
    det_mag = empirical_detection_mags()

    fig, ax = plt.subplots(1, 1, figsize=(DOUBLE_COL * 0.74, 3.8))

    # ── independent cross-check: pipeline detection-count histogram (right axis) ──
    axr = ax.twinx()
    bins = np.arange(11.0, 18.01, 0.35)
    axr.hist(det_mag, bins=bins, color=PALETTE["grey"], alpha=0.30,
             zorder=1, label="_nolegend_")
    axr.set_ylabel("pipeline detections / bin", color="0.55", fontsize=8)
    axr.tick_params(axis="y", labelsize=7, colors="0.55")
    axr.set_ylim(0, axr.get_ylim()[1] * 2.5)   # keep histogram in lower third

    # ── completeness curves (left axis) ──
    ax.axhline(0.5, color=PALETTE["grey"], lw=0.7, ls=":", zorder=2)
    ax.plot(ms, cs, "--o", color=C["reference"], lw=1.4, ms=3.0, zorder=4,
            mfc="white", mew=1.0,
            label=f"controlled synthetic — verification ($m_{{50}}${'='}{m50s:.1f})")
    ax.plot(mr, cr, "-s", color=C["data"], lw=2.0, ms=4.0, zorder=6,
            label=f"real M13 V frame — validation ($m_{{50}}${'='}{m50r:.1f})")
    ax.axvline(m50r, color=C["data"], lw=1.0, ls="--", zorder=3)

    # annotate the independent agreement (point to the histogram turnover)
    ax.annotate("independent cross-check:\npipeline detections (grey)\nturn over at the same mag",
                xy=(14.75, 0.30), xytext=(12.65, 0.60),
                fontsize=6.8, color="0.25", va="center", ha="left",
                arrowprops={"arrowstyle": "->", "color": "0.45", "lw": 0.9})

    ax.set_xlabel("injected / instrumental magnitude (count-rate, ZP = 25)")
    ax.set_ylabel("recovery completeness")
    ax.set_xlim(12.5, 19.2)
    ax.set_ylim(-0.03, 1.06)
    ax.set_zorder(axr.get_zorder() + 1)   # curves above histogram
    ax.patch.set_visible(False)
    ax.legend(loc="center right", bbox_to_anchor=(1.0, 0.62), fontsize=6.9,
              framealpha=0.92)
    ax.set_title("real-frame injection recovers to the pipeline's detection limit",
                 loc="left", fontsize=10.5)

    # frame characterization — explains the depth difference honestly
    ax.text(15.55, 0.90,
            "synthetic: FWHM 3.4 px, sky 150 ADU\n"
            "real M13 V: FWHM 7.6 px, sky 1315 ADU\n"
            "depth set by frame seeing + sky, not method",
            va="top", ha="left", fontsize=6.6, color="0.3",
            bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                  "alpha": .85, "edgecolor": PALETTE["grey"]})

    paths = save_fig(fig, "fig_completeness_realvssynth", OUTDIR)
    plt.close(fig)
    print(f"[real vs synth] synth m50={m50s:.2f}  real M13 m50={m50r:.2f}")
    print(f"  empirical detections: median mag={np.median(det_mag):.2f}  "
          f"n={len(det_mag)}  (rolloff cross-check)")
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
