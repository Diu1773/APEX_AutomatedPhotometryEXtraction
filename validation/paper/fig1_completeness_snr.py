"""Figure — Detection completeness (dense injection-recovery, 21,000 sources).

Top panel (DES-Balrog style): completeness vs injected magnitude from 300
Monte-Carlo trials x 70 injections. With ~21k injections the binned recovery
fraction is itself a smooth curve, so — following DES Balrog, AutoPhOT, and
the artificial-star literature — NO functional form is fitted in magnitude
space; the depths m90/m50/m10 are read off where the measured completeness
crosses each fraction. The isolated / field split quantifies crowding: field
placements (no exclusion around real stars) lose depth relative to isolated
placements. The erf detection-probability form (Masci 2011; Kashyap+2010;
the AutoPhOT appendix) lives in S/N space and is cited in the text.

Bottom row (AutoPhOT Fig-12 style): cutouts of sources injected with the SAME
production injection code (empirical PSF, Poisson photons, gain 0.689) at
representative magnitudes, so the reader sees what a star at the measured
50 per cent depth actually looks like against the noise.

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
from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Circle

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

RUN = REPO / "validation" / "paper" / "data" / "artificial_star" / "benchmark_run"
SUITE = RUN.parent
OUTDIR = REPO / "validation" / "paper" / "figures"

GAIN = 0.689          # e-/ADU, from run manifest (instrument.gain_e_per_adu)
ZP = 25.0             # injection zeropoint, from run summary
CUTOUT_MAGS = [16.0, 17.0, 17.5, 18.3]
CUT = 21              # cutout half-size in px


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def binned_completeness(mags, rec, bw):
    lo0 = np.floor(mags.min() / bw) * bw
    edges = np.arange(lo0, mags.max() + bw, bw)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (mags >= lo) & (mags < hi)
        if m.sum() >= 20:
            k, n = int(rec[m].sum()), int(m.sum())
            w = wilson(k, n)
            out.append((0.5 * (lo + hi), k / n, w[0], w[1]))
    a = np.array(out)
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3]


def read_off(mag, comp, level):
    """Magnitude where the measured completeness crosses `level` (linear interp)."""
    o = np.argsort(mag)
    mm, cc = mag[o], comp[o]
    for i in range(len(mm) - 1):
        if cc[i] >= level >= cc[i + 1]:
            den = cc[i] - cc[i + 1]
            f = (cc[i] - level) / den if den else 0.0
            return float(mm[i] + f * (mm[i + 1] - mm[i]))
    return float("nan")


def main() -> int:
    s = pd.read_csv(RUN / "stars.csv")
    eligible = ~s["baseline_confounded"].astype(bool)
    s = s[eligible].copy()
    mags = s["magnitude_true"].to_numpy(float)
    rec = s["recovered"].to_numpy(bool)
    iso = s["placement_stratum"].astype(str).eq("isolated").to_numpy()

    bw = 0.15
    m_all, c_all, lo_all, hi_all = binned_completeness(mags, rec, bw)
    m_iso, c_iso, *_ = binned_completeness(mags[iso], rec[iso], bw)
    m_fld, c_fld, *_ = binned_completeness(mags[~iso], rec[~iso], bw)

    m50 = read_off(m_all, c_all, 0.5)
    m90 = read_off(m_all, c_all, 0.9)
    m10 = read_off(m_all, c_all, 0.1)
    m50_iso = read_off(m_iso, c_iso, 0.5)
    m50_fld = read_off(m_fld, c_fld, 0.5)

    # completeness at the cutout magnitudes — read off the ISOLATED curve,
    # because the demo cutouts are injected at isolated positions
    comp_at = {m: float(np.interp(m, m_iso, c_iso)) for m in CUTOUT_MAGS}

    # ── demo cutouts: inject with the production code at isolated positions ──
    from apex.benchmark.artificial_stars import inject_catalog
    frame = fits.getdata(SUITE / "synthetic_reference.fits").astype(float)
    psf = fits.getdata(RUN / "empirical_psf.fits").astype(float)
    pos_pool = s[iso & (s["nearest_real_sep_fwhm"] > 4)][["x_true", "y_true"]]
    rng = np.random.default_rng(20260712)
    picks = pos_pool.sample(len(CUTOUT_MAGS), random_state=7).to_numpy()
    cat = pd.DataFrame({
        "injection_id": np.arange(len(CUTOUT_MAGS)),
        "x_true": picks[:, 0], "y_true": picks[:, 1],
        "magnitude_true": CUTOUT_MAGS,
    })
    injected, _, _, _ = inject_catalog(
        frame, psf, cat, gain_e_per_adu=GAIN, zeropoint_mag=ZP, rng=rng)

    # ── figure ──
    fig = plt.figure(figsize=(DOUBLE_COL, 5.4))
    gs = GridSpec(2, len(CUTOUT_MAGS), height_ratios=[2.35, 1.0],
                  hspace=0.34, wspace=0.06, figure=fig)
    ax = fig.add_subplot(gs[0, :])

    # AutoPhOT-style scatter cloud: one dot = one trial's recovery fraction
    # in one magnitude bin (300 MC trials -> a dense cloud), x-jittered
    bws = 0.25
    s["mbin"] = np.floor(s["magnitude_true"] / bws) * bws + bws / 2
    per_trial = (s.groupby(["trial", "mbin"])["recovered"]
                   .agg(["mean", "size"]).reset_index())
    per_trial = per_trial[per_trial["size"] >= 3]
    rng_j = np.random.default_rng(11)
    jit = rng_j.uniform(-0.085, 0.085, len(per_trial))
    ax.plot(per_trial["mbin"] + jit, per_trial["mean"], "o",
            color=C["data"], ms=2.0, alpha=0.10, mew=0, zorder=2,
            label=f"per-trial fractions (300 trials, N={len(s):,})")
    # pooled completeness over all injections
    ax.plot(m_all, c_all, "-", color="k", lw=1.8, zorder=5,
            label="pooled completeness")
    # depth read-offs
    ax.axhline(0.5, color=PALETTE["grey"], lw=0.7, ls=":", zorder=1)
    ax.axhline(0.9, color=PALETTE["grey"], lw=0.5, ls=":", zorder=1, alpha=0.6)
    ax.axvline(m50, color=C["data"], lw=1.2, ls="--", zorder=3)
    ax.annotate(rf"$m_{{50}}={m50:.2f}$", xy=(m50, 0.52), xytext=(-64, 10),
                textcoords="offset points", fontsize=9, fontweight="bold",
                color=C["data"])
    # cutout magnitude markers on the curve
    for mtag in CUTOUT_MAGS:
        ax.plot(mtag, comp_at[mtag], "v", color="k", ms=5, zorder=6)
    ax.set_xlabel("injected magnitude")
    ax.set_ylabel("completeness")
    ax.set_ylim(-0.03, 1.06)
    ax.set_xlim(m_all.min() - 0.2, m_all.max() + 0.2)
    ax.text(0.03, 0.42,
            rf"$m_{{90}}={m90:.2f}$" + "\n" + rf"$m_{{50}}={m50:.2f}$"
            + "\n" + rf"$m_{{10}}={m10:.2f}$"
            + "\n" + rf"$\Delta m_{{50}}$(crowding)$={m50_iso - m50_fld:.2f}$",
            transform=ax.transAxes, va="top", fontsize=8,
            bbox={"boxstyle": "round,pad=0.32", "facecolor": "white",
                  "alpha": .9, "edgecolor": PALETTE["grey"]})
    ax.legend(loc="lower left", fontsize=7.4)
    # right axis in percent
    ax2 = ax.twinx()
    ax2.set_ylim(-3, 106)
    ax2.set_ylabel("recovered [%]")

    # ── cutout row ──
    bg = np.median(frame)
    mad = 1.4826 * np.median(np.abs(frame - bg))
    vmin, vmax = bg - 2.0 * mad, bg + 9.0 * mad
    for j, mtag in enumerate(CUTOUT_MAGS):
        axc = fig.add_subplot(gs[1, j])
        x, y = picks[j]
        xi, yi = int(round(x)), int(round(y))
        cut = injected[yi - CUT:yi + CUT + 1, xi - CUT:xi + CUT + 1]
        axc.imshow(cut, origin="lower", cmap="gray_r", vmin=vmin, vmax=vmax,
                   interpolation="nearest")
        axc.add_patch(Circle((CUT + (x - xi), CUT + (y - yi)), 6.5, fill=False,
                             ec=C["accent"], lw=1.1))
        axc.set_xticks([]); axc.set_yticks([])
        axc.set_title(rf"$m={mtag:.1f}$  ($C\approx{100*comp_at[mtag]:.0f}$%)",
                      fontsize=8)
    fig.text(0.5, 0.015,
             "sources injected with the production pipeline (empirical PSF, "
             f"Poisson photons, gain {GAIN} e$^-$/ADU)",
             ha="center", fontsize=7, color=PALETTE["grey"])

    paths = save_fig(fig, "fig1_completeness", OUTDIR)
    plt.close(fig)
    print(f"[dense] N={len(s)} m50={m50:.2f} m90={m90:.2f} m10={m10:.2f} "
          f"| iso={m50_iso:.2f} field={m50_fld:.2f} dCrowd={m50_iso-m50_fld:.2f}")
    print(f"saved: {paths['png']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
