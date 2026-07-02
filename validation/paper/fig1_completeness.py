"""Figure 1 — APEX detection completeness from artificial-star tests.

Question: how deep does APEX's forced-photometry pipeline reliably *recover*
sources? The standard answer is the completeness function C(m) = fraction of
injected artificial stars recovered as a function of injected magnitude, and its
characteristic depth m50 (the magnitude at which C = 0.5).

This script is a pure CONSUMER of an already-generated artificial-star benchmark
(``validation/paper/data/artificial_star/benchmark_run/``). Nothing is
re-simulated here:

  * ``magnitude_points.csv`` holds the per-injection outcome for every one of the
    injected artificial stars: its injected magnitude, whether it was eligible
    (fell on a usable part of a frame), and whether it was recovered. We bin
    these raw Bernoulli outcomes into magnitude bins and attach an *exact*
    Wilson 95% binomial interval to each bin's recovered fraction — the honest
    small-sample interval for a proportion.
  * ``completeness_fit.json`` holds the logistic completeness model fitted with a
    cluster bootstrap over the 12 injection trials:
        p(m) = expit((m50 - m) / width)
    together with m50, m90, m10, the logistic width, and each quantity's
    bootstrap 95% CI. We overlay that model and mark m50 / m90 / m10; we do NOT
    re-fit it.

One clean panel: binned Wilson points + logistic model + m50 (with its 95% CI
band) + m90 / m10 depth markers, with a slim lower strip showing how many
eligible injections fall in each magnitude bin (sampling context).
"""

from __future__ import annotations

import sys
import json
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.special import expit

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE

apply_paper_style()

# ── Input data (already generated — never re-run the benchmark) ──────────────
DATA_DIR = REPO / "validation" / "paper" / "data" / "artificial_star" / "benchmark_run"
POINTS_CSV = DATA_DIR / "magnitude_points.csv"
FIT_JSON = DATA_DIR / "completeness_fit.json"

# Magnitude binning for the plotted Wilson points. The raw CSV is one row per
# injected star (n_injected = 1), so we aggregate into fixed-width bins to form
# real binomial samples; the logistic fit itself is untouched (comes from JSON).
BIN_WIDTH = 0.25  # mag


def wilson_interval(k: int, n: int, z: float = 1.959963984540054):
    """Wilson score 95% interval for k successes out of n Bernoulli trials.

    This is the same small-sample binomial interval the benchmark reports on the
    per-magnitude rows; recomputing it here on the *binned* counts keeps the
    plotted error bars internally consistent with that convention.
    """
    if n == 0:
        return np.nan, np.nan
    phat = k / n
    denom = 1.0 + z * z / n
    centre = (phat + z * z / (2.0 * n)) / denom
    half = (z / denom) * np.sqrt(phat * (1.0 - phat) / n + z * z / (4.0 * n * n))
    return centre - half, centre + half


def load_data():
    df = pd.read_csv(POINTS_CSV)
    with open(FIT_JSON, "r", encoding="utf-8") as fh:
        fit = json.load(fh)
    return df, fit


def bin_completeness(df: pd.DataFrame, bin_width: float = BIN_WIDTH):
    """Aggregate per-injection rows into magnitude bins.

    Only *eligible* injections (n_eligible == 1) count toward the denominator —
    an injection that landed off a usable region can neither be recovered nor
    missed, so including it would bias completeness low.
    """
    elig = df[df["n_eligible"] == 1].copy()
    mmin = float(np.floor(elig["magnitude"].min() / bin_width) * bin_width)
    mmax = float(np.ceil(elig["magnitude"].max() / bin_width) * bin_width)
    edges = np.arange(mmin, mmax + 0.5 * bin_width, bin_width)
    centres = 0.5 * (edges[:-1] + edges[1:])

    idx = np.clip(np.digitize(elig["magnitude"].to_numpy(), edges) - 1, 0, len(centres) - 1)
    rows = []
    for b in range(len(centres)):
        sel = idx == b
        n = int(sel.sum())
        if n == 0:
            continue
        k = int(elig["n_recovered"].to_numpy()[sel].sum())
        phat = k / n
        lo, hi = wilson_interval(k, n)
        rows.append(
            dict(mag=centres[b], n_eligible=n, n_recovered=k,
                 completeness=phat, ci_low=lo, ci_high=hi)
        )
    return pd.DataFrame(rows)


def make_figure(binned: pd.DataFrame, fit: dict, n_injected: int):
    m50 = fit["m50"]
    width = fit["width_mag"]
    m90, m10 = fit["m90"], fit["m10"]
    m50_lo, m50_hi = fit["m50_ci95_low"], fit["m50_ci95_high"]

    grid = np.linspace(binned["mag"].min() - 0.15, binned["mag"].max() + 0.15, 600)
    model = expit((m50 - grid) / width)

    mags = binned["mag"].to_numpy()
    comp = binned["completeness"].to_numpy()
    ci_lo = binned["ci_low"].to_numpy()
    ci_hi = binned["ci_high"].to_numpy()
    # Clip tiny negative values from float rounding (e.g. at C=1 bins where the
    # Wilson upper edge coincides with the point); errorbar rejects negatives.
    yerr = np.clip(np.vstack([comp - ci_lo, ci_hi - comp]), 0.0, None)

    fig = plt.figure(figsize=(4.6, 3.6))
    gs = fig.add_gridspec(
        2, 1, height_ratios=[4.2, 1.0], hspace=0.06, left=0.135,
        right=0.975, top=0.975, bottom=0.115,
    )
    ax = fig.add_subplot(gs[0])
    ax_n = fig.add_subplot(gs[1], sharex=ax)

    # ── m50 95% CI band + guide lines ────────────────────────────────────────
    ax.axvspan(m50_lo, m50_hi, color=C["floor"], alpha=0.16, lw=0,
               zorder=1, label=r"$m_{50}$ 95% CI")
    ax.axhline(0.5, color=C["floor"], lw=0.9, ls=":", zorder=1)
    ax.axvline(m50, color=C["floor"], lw=1.1, ls="--", zorder=2)

    # Faint m90 / m10 depth markers.
    for md, lab in ((m90, r"$m_{90}$"), (m10, r"$m_{10}$")):
        ax.axvline(md, color=C["accent"], lw=0.9, ls="-.", alpha=0.7, zorder=2)
        yv = 0.9 if lab.endswith("90}$") else 0.1
        ax.annotate(lab, xy=(md, yv), xytext=(3, 0), textcoords="offset points",
                    va="center", ha="left", fontsize=7.5, color=C["accent"])

    # ── Logistic model curve ─────────────────────────────────────────────────
    ax.plot(grid, model, color=C["model"], lw=1.8, zorder=4,
            label=r"logistic fit")

    # ── Binned Wilson injection points ───────────────────────────────────────
    ax.errorbar(mags, comp, yerr=yerr, fmt="o", ms=4.0, mfc=C["data"],
                mec="white", mew=0.5, ecolor=C["data"], elinewidth=0.9,
                capsize=1.6, capthick=0.9, zorder=5,
                label="injections (Wilson 95%)")

    ax.set_ylim(-0.03, 1.03)
    ax.set_ylabel("Completeness")
    ax.set_xlim(grid.min(), grid.max())
    plt.setp(ax.get_xticklabels(), visible=False)

    # Compact info box (m50, width, m90/m10, N).
    txt = (
        rf"$m_{{50}} = {m50:.2f}\,^{{+{m50_hi - m50:.2f}}}_{{-{m50 - m50_lo:.2f}}}$"
        + "\n"
        + rf"width $= {width:.2f}$ mag"
        + "\n"
        + rf"$m_{{90}} = {m90:.2f}$, $m_{{10}} = {m10:.2f}$"
        + "\n"
        + rf"$N = {n_injected}$ stars / {fit['n_trials']} trials"
    )
    ax.text(0.035, 0.10, txt, transform=ax.transAxes, va="bottom", ha="left",
            fontsize=7.6, linespacing=1.35,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec=C["floor"], lw=0.6))

    ax.legend(loc="upper right", frameon=False, fontsize=7.6,
              handlelength=1.5, borderaxespad=0.3)

    # ── Lower strip: eligible injections per magnitude bin ───────────────────
    ax_n.bar(mags, binned["n_eligible"].to_numpy(), width=BIN_WIDTH * 0.9,
             color=C["data"], alpha=0.30, edgecolor="none", zorder=2)
    ax_n.axvspan(m50_lo, m50_hi, color=C["floor"], alpha=0.16, lw=0, zorder=1)
    ax_n.axvline(m50, color=C["floor"], lw=1.1, ls="--", zorder=2)
    ax_n.set_ylabel(r"$N_{\rm elig}$", fontsize=8)
    ax_n.set_xlabel("Injected magnitude")
    ax_n.set_ylim(0, binned["n_eligible"].max() * 1.25)
    ax_n.tick_params(labelsize=7)

    stats = dict(m50=m50, m50_lo=m50_lo, m50_hi=m50_hi, width=width,
                 m90=m90, m10=m10, n_injected=n_injected, n_eligible=fit["n_stars"],
                 n_trials=fit["n_trials"], n_bins=int(len(binned)))
    return fig, stats


def _model_vs_points(binned: pd.DataFrame, fit: dict) -> dict:
    """Quick goodness check: do binned Wilson points bracket the logistic model?"""
    p = expit((fit["m50"] - binned["mag"].to_numpy()) / fit["width_mag"])
    inside = (p >= binned["ci_low"].to_numpy()) & (p <= binned["ci_high"].to_numpy())
    return dict(n_bins=int(len(binned)),
                n_inside=int(np.sum(inside)),
                frac_inside=float(np.mean(inside)))


def main():
    print("[fig1] loading artificial-star benchmark ...")
    df, fit = load_data()
    n_injected = int(len(df))
    binned = bin_completeness(df)
    print(f"[fig1] {n_injected} injections ({int((df['n_eligible'] == 1).sum())} eligible) "
          f"-> {len(binned)} magnitude bins (width {BIN_WIDTH} mag)")

    fig, stats = make_figure(binned, fit, n_injected)
    outdir = REPO / "validation" / "paper" / "figures"
    paths = save_fig(fig, "fig1_completeness", outdir)
    print(f"[fig1] wrote {paths['pdf']}")
    print(f"[fig1] wrote {paths['png']}")

    check = _model_vs_points(binned, fit)
    print("\n=== KEY NUMBERS ===")
    print(f"m50   = {stats['m50']:.3f}  (95% CI {stats['m50_lo']:.3f} .. {stats['m50_hi']:.3f})")
    print(f"m90   = {stats['m90']:.3f}")
    print(f"m10   = {stats['m10']:.3f}")
    print(f"width = {stats['width']:.3f} mag")
    print(f"N     = {stats['n_injected']} injected ({stats['n_eligible']} eligible) / {stats['n_trials']} trials")
    print(f"model-vs-Wilson: {check['n_inside']}/{check['n_bins']} bins "
          f"({100 * check['frac_inside']:.0f}%) contain the logistic curve inside their 95% CI")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
