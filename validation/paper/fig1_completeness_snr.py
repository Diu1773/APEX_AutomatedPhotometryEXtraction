"""Figure 1 (v2) — Detection completeness, two panels.

(a) Completeness vs injected magnitude: measured recovery fractions (Wilson 95%
    CIs) with an empirical logistic *summary* (not a theoretical prediction);
    m50/m90/m10 marked.
(b) Completeness vs source S/N: the same recoveries collapse onto a single clean
    threshold near S/N ~ 7.4 — the physical reason the magnitude completeness has
    its sigmoid shape (detection is S/N-governed; the logistic is the smooth
    approximation to this noise-threshold crossing).

Per-star S/N is derived from the injected flux and a noise term calibrated on the
recovered stars' own reported magnitude errors, so it applies to detected and
undetected injections alike.

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig1_completeness_snr.py
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
from scipy.special import expit, erf
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

RUN = REPO / "validation" / "paper" / "data" / "artificial_star" / "benchmark_run"
BINS = RUN / "magnitude_bins.csv"
FIT = RUN / "completeness_fit.json"
STARS = RUN / "stars.csv"
OUTDIR = REPO / "validation" / "paper" / "figures"


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def erfmod(x, x50, w):
    return 0.5 * (1 + erf((x - x50) / (np.sqrt(2) * w)))


def main() -> int:
    fit = json.loads(FIT.read_text())
    m50, width = fit["m50"], fit["width_mag"]
    m90, m10 = fit["m90"], fit["m10"]
    m50lo, m50hi = fit["m50_ci95_low"], fit["m50_ci95_high"]

    # per-star data (eligible = reasonably isolated injections)
    s = pd.read_csv(STARS)
    s = s[s["nearest_real_sep_fwhm"] > 2].copy()
    flux = s["flux_realized_e"].to_numpy(float)
    rec = s["recovered"].to_numpy(bool)
    magt = s["magnitude_true"].to_numpy(float)

    # (a) completeness vs magnitude, binned at 0.25 mag from the stars themselves
    bw = 0.25
    lo0 = np.floor(magt.min() / bw) * bw
    medges = np.arange(lo0, magt.max() + bw, bw)
    mag, comp, ci = [], [], []
    for lo, hi in zip(medges[:-1], medges[1:]):
        m = (magt >= lo) & (magt < hi)
        if m.sum() >= 5:
            k, n = int(rec[m].sum()), int(m.sum())
            mag.append(0.5 * (lo + hi)); comp.append(k / n); ci.append(wilson(k, n))
    mag = np.array(mag); comp = np.array(comp); ci = np.array(ci)
    me = pd.to_numeric(s["forced_mag_error"], errors="coerce").to_numpy(float)
    ok = rec & np.isfinite(me) & (me > 0)
    snr_meas = 1.0857 / me[ok]
    Carr = (flux[ok] / snr_meas) ** 2 - flux[ok]
    Cnoise = float(np.median(Carr[np.isfinite(Carr) & (Carr > 0)]))
    snr = flux / np.sqrt(flux + Cnoise)
    p_snr, _ = curve_fit(erfmod, snr, rec.astype(float), p0=[7, 1], maxfev=20000)

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.5))

    # ── (a) completeness vs magnitude ──
    axa.axvspan(m50lo, m50hi, color=C["floor"], alpha=0.16, lw=0, zorder=1)
    axa.axvline(m50, color=C["floor"], lw=1.1, ls="--", zorder=2)
    axa.axhline(0.5, color=C["floor"], lw=0.8, ls=":", zorder=1)
    for mm, lo, hi in [(m90, None, None), (m10, None, None)]:
        axa.axvline(mm, color=C["accent"], lw=0.9, ls="-.", zorder=1, alpha=.8)
    grid = np.linspace(mag.min() - 0.3, mag.max() + 0.3, 300)
    axa.plot(grid, expit((m50 - grid) / width), color=C["model"], lw=1.8,
             zorder=4, label="empirical summary (logistic)")
    axa.errorbar(mag, comp,
                 yerr=[np.clip(comp - ci[:, 0], 0, None), np.clip(ci[:, 1] - comp, 0, None)],
                 fmt="o", color=C["data"], ms=4, lw=1, capsize=2, zorder=5,
                 label="injections (Wilson 95%)")
    axa.set_xlabel("injected magnitude"); axa.set_ylabel("completeness")
    axa.set_ylim(-0.03, 1.05)
    axa.text(0.05, 0.30,
             rf"$m_{{50}}={m50:.2f}^{{+{m50hi-m50:.2f}}}_{{-{m50-m50lo:.2f}}}$"
             + f"\nwidth $={width:.2f}$\n$m_{{90}}={m90:.2f}$, $m_{{10}}={m10:.2f}$",
             transform=axa.transAxes, va="top", fontsize=7.6,
             bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                   "alpha": .85, "edgecolor": PALETTE["grey"]})
    axa.legend(loc="upper right", fontsize=7.2)
    axa.set_title("(a) vs magnitude — measured + summary", loc="left")

    # ── (b) completeness vs S/N ──
    edges = np.array([2, 4, 5, 6, 7, 8, 9, 11, 14, 20, 40])
    cx, cy, clo, chi = [], [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (snr >= lo) & (snr < hi)
        if m.sum() >= 8:
            k, n = int(rec[m].sum()), int(m.sum())
            cx.append(np.median(snr[m])); cy.append(k / n)
            l, h = wilson(k, n); clo.append(max(0, k/n - l)); chi.append(max(0, h - k/n))
    cx = np.array(cx)
    gg = np.linspace(2, 30, 300)
    axb.axvline(p_snr[0], color=C["floor"], lw=1.1, ls="--", zorder=2)
    axb.axhline(0.5, color=C["floor"], lw=0.8, ls=":", zorder=1)
    axb.plot(gg, erfmod(gg, *p_snr), color=C["reference"], lw=1.8, zorder=4,
             label=r"S/N-threshold (erf)")
    axb.errorbar(cx, cy, yerr=[clo, chi], fmt="s", color=C["data"], ms=4, lw=1,
                 capsize=2, zorder=5, label="recovered fraction")
    axb.set_xscale("log")
    axb.set_xlabel("source S/N"); axb.set_ylabel("completeness")
    axb.set_ylim(-0.03, 1.05); axb.set_xlim(2.5, 30)
    axb.text(0.05, 0.92, rf"50% at S/N $\approx {p_snr[0]:.1f}$",
             transform=axb.transAxes, va="top", fontsize=8,
             bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                   "alpha": .85, "edgecolor": PALETTE["grey"]})
    axb.legend(loc="lower right", fontsize=7.2)
    axb.set_title("(b) vs S/N — detection is S/N-governed", loc="left")

    fig.tight_layout()
    paths = save_fig(fig, "fig1_completeness", OUTDIR)
    plt.close(fig)
    print(f"m50={m50:.2f} width={width:.2f} | SNR50={p_snr[0]:.2f} width={p_snr[1]:.2f} Cnoise={Cnoise:.0f}")
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
