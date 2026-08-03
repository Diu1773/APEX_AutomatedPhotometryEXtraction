"""LC-mode reproduction figure: YZ Boötis (HADS delta Scuti).

Literature: YZ Boo P = 0.104091579 d, HADS, V pk-pk ~0.42 mag
(Yang et al. 2018, RAA 18, 2 = arXiv:1709.08798).

Panel (a): a single good night (2026-03-28, r band, 5.2 h) folded on the
LITERATURE period reproduces the published high-amplitude sawtooth.
Panel (b): Lomb-Scargle periodograms -- the single night peaks at the true
period, while a 2-night (1-day-gap) merge's top peak jumps to the +1 c/d
alias (0.095 d); the true period survives only as a secondary peak. This is
the standard ground-based multi-night spectral-window aliasing, reported
honestly as a domain limitation, not hidden.

Data SSD required. Run:
    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_lc_yzboo.py
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.timeseries import LombScargle

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

P_LIT = 0.104091579  # Yang+2018
GOOD = REPO.parent  # unused placeholder
GOOD_NIGHT = r"E:/observed_Analysis/YZbootis/RESULT_YZbootis_20260328/step10_lightcurve/lightcurve_ID1_raw.csv"
MERGED = r"E:/observed_Analysis/YZbootis/MERGED_YZ bootis_20250429_20250430/step10_lightcurve/lightcurve_ID1_raw.csv"
OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"


def load(path, band="r"):
    """Load a single photometric band. The raw light-curve tables interleave
    g/r/i frames with different zero-points, so a band cut is required before
    folding or period analysis (mixing bands creates spurious parallel
    sequences and inflates the amplitude)."""
    d = pd.read_csv(path)
    if band is not None and "filter" in d.columns:
        d = d[d["filter"] == band]
    t = d["BJD_TDB"].to_numpy(float)
    y = d["diff_mag_raw"].to_numpy(float) if "diff_mag_raw" in d.columns else d["mag"].to_numpy(float)
    dy = d["mag_err"].to_numpy(float) if "mag_err" in d.columns else None
    m = np.isfinite(t) & np.isfinite(y)
    return t[m], y[m], (dy[m] if dy is not None else None)


def periodogram(t, y, pmin=0.03, pmax=1.5):
    f, p = LombScargle(t, y).autopower(minimum_frequency=1/pmax, maximum_frequency=1/pmin,
                                       samples_per_peak=20)
    return 1.0 / f, p


def main() -> int:
    tg, yg, dyg = load(GOOD_NIGHT)
    tm, ym, dym = load(MERGED)

    # best period on the good night (fine)
    fg, pg = LombScargle(tg, yg).autopower(minimum_frequency=1/0.30, maximum_frequency=1/0.03,
                                           samples_per_peak=50)
    p_good = 1.0 / fg[np.argmax(pg)]

    per_g, pow_g = periodogram(tg, yg)
    per_m, pow_m = periodogram(tm, ym)
    p_merged = per_m[np.argmax(pow_m)]

    # sigma-clip outliers for a clean fold display
    phase = ((tg - tg.min()) / P_LIT) % 1.0
    med = np.median(yg)
    sig = 1.4826 * np.median(np.abs(yg - med))
    keep = np.abs(yg - med) < 4 * sig
    amp = float(np.ptp(np.percentile(yg[keep], [2, 98])))

    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.2))

    # (a) phase fold on literature period, two cycles
    for shift in (0.0, 1.0):
        ax_a.scatter(phase[keep] + shift, yg[keep], s=9, alpha=0.55,
                     color=C["data"], edgecolors="none")
    ax_a.invert_yaxis()
    ax_a.set_xlim(0, 2)
    ax_a.set_xlabel(r"phase  ($P_{\rm lit} = %.5f$ d)" % P_LIT)
    ax_a.set_ylabel(r"$\Delta r$  (target $-$ comparison, mag)")
    ax_a.text(0.04, 0.06,
              f"YZ Boo, single night (2026-03-28, r)\n"
              f"N={keep.sum()}, span {np.ptp(tg)*24:.1f} h\n"
              f"APEX best P = {p_good:.5f} d\n"
              f"pk-pk $\\approx$ {amp:.2f} mag",
              transform=ax_a.transAxes, va="bottom", ha="left", fontsize=7.0,
              bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
                    "alpha": 0.85, "edgecolor": PALETTE["grey"]})
    ax_a.set_title("(a) folded on literature period", loc="left")

    # (b) periodograms single vs merged
    ax_b.plot(per_g, pow_g / pow_g.max(), color=C["data"], lw=1.0,
              label=f"single night (peak {p_good:.4f} d)")
    ax_b.plot(per_m, pow_m / pow_m.max(), color=C["accent"], lw=1.0, ls="--", alpha=0.85,
              label=f"2-night merge (peak {p_merged:.4f} d)")
    ax_b.axvline(P_LIT, color=C["model"], ls="--", lw=1.2,
                 label=f"literature {P_LIT:.4f} d")
    ax_b.set_xscale("log")
    ax_b.set_xlim(0.03, 1.5)
    ax_b.set_xlabel("period (d, log)")
    ax_b.set_ylabel("normalized LS power")
    ax_b.legend(loc="upper left", fontsize=6.6)
    ax_b.set_title("(b) multi-night 1-day aliasing", loc="left")

    fig.suptitle("YZ Boötis (HADS $\\delta$ Sct) — LC-mode period reproduction & multi-night alias",
                 fontsize=8.0, y=1.02, color="#333333")
    fig.tight_layout(rect=(0, 0, 0.99, 1))
    paths = save_fig(fig, "fig_lc_yzboo", OUTDIR)
    plt.close(fig)

    print(f"good-night best P = {p_good:.6f} d (lit {P_LIT:.6f}, {abs(p_good-P_LIT)/P_LIT*100:.1f}% off)")
    print(f"merged 2-night top P = {p_merged:.6f} d  (1-day alias of true period)")
    print(f"fold amplitude pk-pk ~ {amp:.3f} mag (r band)")
    for ext, p in paths.items():
        print(f"wrote {ext}: {p}  exists={p.exists()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
