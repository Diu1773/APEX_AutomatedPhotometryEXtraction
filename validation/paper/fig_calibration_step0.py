# -*- coding: utf-8 -*-
"""Figure — detector calibration (Step 0): the frames applied, and their effect.

2026-08-05 paired design. The first attempt put three abstract bar charts next
to one before/after pair; the user's objections were all correct:
  * the before/after cutouts looked identical, so the panel proved nothing
  * "recovered - true" was an opaque axis label
  * a "star cores touched 0.000%" bar sat on the same axis as two ~100% bars
    and read as noise
  * what a reader actually wants is the calibration frames themselves, and the
    light frame before and after each operation

So this figure pairs each master image with a distribution or profile, then
shows the science frame before and after calibration. The quantitative
pixel-level cross-check remains Fig. 4; this figure explains what was applied
and makes the effect legible at a larger panel size.

Data: data/calib_stages.npz, written by calib_crosscheck_ngc6811.py from the
real 2026-06-11 night (8 bias, 8x60 s dark, 5 B flats, one 60 s NGC 6811 B).

Run: .venv-deploy\\Scripts\\python -X utf8 validation\\paper\\fig_calibration_step0.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

DATA = REPO / "validation" / "paper" / "data"
OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"


def show(ax, img, title, lo=1.0, hi=99.0, vmin=None, vmax=None, cmap="gray"):
    v = img[np.isfinite(img)]
    if vmin is None:
        vmin, vmax = np.percentile(v, lo), np.percentile(v, hi)
    ax.imshow(img, cmap=cmap, origin="lower", aspect="equal",
              vmin=vmin, vmax=vmax, interpolation="nearest")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, loc="left", fontsize=7.6, pad=2.5)
    return vmin, vmax


def note(ax, text, color=None):
    """이미지 패널 아래 한 줄. 눈금이 없는 패널 전용."""
    ax.text(0.5, -0.045, text, transform=ax.transAxes, fontsize=6.6,
            ha="center", va="top", color=color or PALETTE["grey"])


def inset_note(ax, text, x=0.5, y=0.055, ha="center", color=None):
    """눈금이 있는 선 그래프는 축 안에 적는다 — 축 밖에 두면 눈금과 겹친다."""
    ax.text(x, y, text, transform=ax.transAxes, fontsize=6.6,
            ha=ha, va="bottom", color=color or PALETTE["grey"],
             bbox=dict(fc="white", ec="none", alpha=0.8, pad=1.0))


def histogram(ax, image, title, color, stat_text, lo=0.2, hi=99.8):
    """Show a calibration-frame value distribution beside its image."""
    v = image[np.isfinite(image)]
    low, high = np.percentile(v, [lo, hi])
    if high <= low:
        high = low + 1.0
    ax.hist(v, bins=42, range=(low, high), color=color, edgecolor="white",
            linewidth=0.25, density=True)
    med = float(np.median(v))
    ax.axvline(med, color="k", lw=0.85, ls="--")
    ax.set_title(title, loc="left", fontsize=7.6, pad=2.5)
    ax.set_xlabel("DN", fontsize=7.0, labelpad=1)
    ax.set_ylabel("density", fontsize=7.0, labelpad=1)
    ax.tick_params(labelsize=6.2)
    ax.text(0.03, 0.95, stat_text, transform=ax.transAxes, va="top",
            fontsize=6.2, color=PALETTE["grey"],
            bbox=dict(fc="white", ec="none", alpha=0.82, pad=1.0))


def profile(ax, values, title, ylabel, color, labels=True):
    x = np.linspace(0, 1, values.size)
    ax.plot(x, values, color=color, lw=1.25)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1.0])
    if labels:
        ax.set_xticklabels(["left", "centre", "right"], fontsize=6.2)
    else:
        ax.set_xticklabels([])
    ax.set_ylabel(ylabel, fontsize=7.0, labelpad=1)
    ax.tick_params(labelsize=6.2)
    ax.set_title(title, loc="left", fontsize=7.6, pad=2.5)


def main() -> int:
    d = np.load(DATA / "calib_stages.npz")
    res = json.loads((DATA / "calib_crosscheck_ngc6811.json").read_text(encoding="utf-8"))

    mb, md, mf = d["master_bias"], d["master_dark"], d["master_flat"]
    raw, bd, cal = d["light_raw"], d["light_bd"], d["light_cal"]
    p_bd, p_cal, p_flat = d["prof_bd"], d["prof_cal"], d["prof_flat"]
    fig = plt.figure(figsize=(DOUBLE_COL, 4.75))
    gs = fig.add_gridspec(2, 4, wspace=0.28, hspace=0.62,
                          left=0.040, right=0.985, top=0.920, bottom=0.155)

    # ── row 1: each master is paired with a quantitative view ──
    ax = fig.add_subplot(gs[0, 0]); show(ax, mb, "(a) master bias")
    note(ax, f"median {float(d['bias_med']):.0f} DN · 8 frames")
    v = mb[np.isfinite(mb)]
    histogram(fig.add_subplot(gs[0, 1]), mb, "(b) bias distribution", C["data"],
              f"median {np.median(v):.2f}\nσ {np.std(v):.2f} DN")

    ax = fig.add_subplot(gs[0, 2]); show(ax, md, "(c) master dark", lo=1.0, hi=99.5)
    note(ax, f"median {float(d['dark_med']):.1f} DN · {float(d['dark_exp']):.0f} s · 8 frames")
    v = md[np.isfinite(md)]
    histogram(fig.add_subplot(gs[0, 3]), md, "(d) dark distribution", C["reference"],
              f"median {np.median(v):.2f}\n95th pct {np.percentile(v, 95):.2f} DN",
              lo=0.2, hi=99.5)

    # ── row 2: flat pair and science before/after ──
    ax = fig.add_subplot(gs[1, 0]); show(ax, mf, "(e) master flat", lo=0.5, hi=99.5)
    note(ax, f"median {float(d['flat_med']):.3f} · 5 frames")
    axf = fig.add_subplot(gs[1, 1]); profile(axf, p_flat, "(f) flat profile", "flat", C["reference"])
    axf.axhline(1.0, color=PALETTE["grey"], lw=0.8, ls=":")
    lo_f, hi_f = float(np.min(p_flat)), float(np.max(p_flat))
    axf.set_ylim(lo_f - 0.015, hi_f + 0.02)
    inset_note(axf, f"{(hi_f - lo_f) * 100:.0f}% across field", x=0.03, y=0.80, ha="left")

    v = bd[np.isfinite(bd)]
    vmin, vmax = np.percentile(v, 1.0), np.percentile(v, 99.0)
    ax = fig.add_subplot(gs[1, 2]); show(ax, raw, "(g) raw science frame")
    note(ax, f"NGC 6811 $B$, {float(d['light_exp']):.0f} s")
    axh = fig.add_subplot(gs[1, 3]); show(axh, cal, "(h) calibrated science", vmin=vmin, vmax=vmax)
    note(axh, "raw $\\rightarrow$ $-$bias$-$dark $\\rightarrow$ $\\div$flat")

    # Compact sky-profile inset: it is subordinate to the calibrated image,
    # rather than consuming a full additional panel.
    xb = np.linspace(0, 1, p_bd.size)
    n_bd = p_bd / float(np.median(p_bd))
    n_cal = p_cal / float(np.median(p_cal))
    inset = axh.inset_axes([0.05, 0.05, 0.54, 0.28], facecolor="white")
    inset.plot(xb, n_bd, color=C["bad"], lw=1.0, label="before")
    inset.plot(xb, n_cal, color=C["data"], lw=1.0, label="after")
    inset.axhline(1.0, color=PALETTE["grey"], lw=0.6, ls=":")
    inset.set_xlim(0, 1); inset.set_ylim(0.90, 1.10)
    inset.set_xticks([]); inset.set_yticks([])
    span_before = (float(np.max(n_bd)) - float(np.min(n_bd))) * 100
    span_after = (float(np.max(n_cal)) - float(np.min(n_cal))) * 100
    inset.text(0.02, 0.92, "before", transform=inset.transAxes,
               fontsize=5.2, color=C["bad"], va="top",
               bbox=dict(fc="white", ec="none", alpha=0.82, pad=0.2))
    inset.text(0.02, 0.72, "after", transform=inset.transAxes,
               fontsize=5.2, color=C["data"], va="top",
               bbox=dict(fc="white", ec="none", alpha=0.82, pad=0.2))
    inset.text(0.98, 0.10, f"{span_before:.1f}% $\\to$ {span_after:.1f}%",
               transform=inset.transAxes, fontsize=5.0, ha="right", va="bottom",
               color=PALETTE["grey"],
               bbox=dict(fc="white", ec="none", alpha=0.82, pad=0.2))

    # 한 줄로 쓰면 그림 폭을 넘어 tight-bbox 가 캔버스를 왼쪽으로 늘린다
    # (2026-08-03 실측: 그림 왼쪽에 빈 여백 2인치). 두 줄로 나눈다.
    prov1 = (f"All frames real: Moravian C3-61000 (2$\\times$2), night {res['night']} · "
             f"{res['n_bias']} bias, {res['n_dark']} dark ({float(d['dark_exp']):.0f} s), "
              f"{res['n_flat']} flat, one {float(d['light_exp']):.0f} s NGC 6811 $B$ light")
    prov2 = (f"displayed at 1/{int(d['shrink'])} scale (block-mean) · "
             "gen: calib_crosscheck_ngc6811.py $\\to$ fig_calibration_step0.py")
    fig.text(0.985, 0.048, prov1, ha="right", va="bottom", fontsize=5.9,
             color=PALETTE["grey"])
    fig.text(0.985, 0.010, prov2, ha="right", va="bottom", fontsize=5.9,
             color=PALETTE["grey"])

    paths = save_fig(fig, "fig_calibration_step0", OUTDIR)
    plt.close(fig)

    CAPDIR.mkdir(parents=True, exist_ok=True)
    (CAPDIR / "fig_calibration_step0.md").write_text(
        f"""# Figure — detector calibration (Step 0)

**(a)** Master bias image and **(b)** its pixel-value distribution, including
the median and within-frame scatter. **(c)** Master dark image and **(d)** its
distribution; the long positive tail is retained in the displayed percentile
stretch. **(e)** Master flat image and **(f)** its horizontal profile, showing a
{(float(np.max(p_flat)) - float(np.min(p_flat))) * 100:.0f} per cent response
gradient. **(g)** Raw NGC 6811 $B$ science frame and **(h)** the calibrated
frame after bias/dark subtraction and flat division. The small inset in (h)
shows the sky profile before and after flat-fielding, normalised to each
median; the peak-to-peak gradient falls from {span_before:.1f} to {span_after:.1f}
per cent. All frames are real
(Moravian C3-61000, 2x2, night {res['night']}; {res['n_bias']} bias,
{res['n_dark']} darks, {res['n_flat']} flats), shown at 1/{int(d['shrink'])}
scale. Numerical equivalence to the independent Python `ccdproc` package is
tested separately in Fig. 4.
""", encoding="utf-8")

    print("=== fig_calibration_step0 (v3) ===")
    print(f"  flat gradient {(float(np.max(p_flat)) - float(np.min(p_flat))) * 100:.1f}%")
    print(f"  sky peak-to-peak {span_before:.2f}% -> {span_after:.2f}%")
    print(f"  bias {float(d['bias_med']):.1f} DN · dark {float(d['dark_med']):.2f} DN"
          f" · flat {float(d['flat_med']):.4f}")
    for ext, p in paths.items():
        print(f"wrote {ext}: {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
