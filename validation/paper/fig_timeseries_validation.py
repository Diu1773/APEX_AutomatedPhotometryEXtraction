"""Figure — the two self-implemented time-series modules, checked against injected truth.

Table 1 lists PDM and SYSREM as written for APEX rather than taken from a
library, so both are validated directly here. Everything runs through the same
entry points the LC-mode step windows call:

    apex.analysis.light_curve.period_analysis_service.run_period_analysis
    apex.analysis.light_curve.sysrem.sysrem

(a) PDM period recovery. A high-amplitude delta-Scuti-like signal (fundamental
    plus one harmonic, so the curve is asymmetric like a real sawtooth) is
    injected on the *same* observing pattern as the science application of
    Section 4 — one night, 5.2 h, 80 points — over a grid of periods and noise
    levels. Recovery is scored as |P_rec - P_true| / P_true.

(b) PDM against Lomb-Scargle on the identical series. LS comes from
    astropy.timeseries (a library component); agreement between the two is the
    cross-check that the hand-written PDM is not solving a different problem.

(c) SYSREM. A star x frame matrix carries a common systematic (per-frame
    transparency times per-star sensitivity) on top of constant stars, with one
    genuine variable hidden among them. What matters is not only that the
    systematic goes away but that the injected variability survives: a detrender
    that flattens everything would score perfectly on the first test alone.

Fixed seed, no external data.

    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig_timeseries_validation.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import matplotlib.pyplot as plt

from apex.analysis.light_curve.period_analysis_service import run_period_analysis
from apex.analysis.light_curve.sysrem import sysrem, apply_sysrem
from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, DOUBLE_COL

apply_paper_style()

OUTDIR = REPO / "validation" / "paper" / "figures"
DATADIR = REPO / "validation" / "paper" / "data_timeseries"
SEED = 20260801

# Section 4's observing pattern: one night, 5.2 h, 80 points.
N_POINTS = 80
SPAN_HOURS = 5.2
PERIODS = [0.06, 0.08, 0.104092, 0.13, 0.16, 0.20]   # days; 0.104092 = YZ Boo
NOISE = [0.005, 0.010, 0.020, 0.040]                  # mag rms
N_TRIAL = 12                                          # noise realizations per cell


def make_series(rng, period, noise, n=N_POINTS, span_h=SPAN_HOURS):
    """delta-Scuti-like: fundamental + harmonic, so rise and fall differ."""
    t = np.sort(rng.uniform(0, span_h / 24.0, n))
    ph = 2 * np.pi * t / period
    mag = 0.21 * np.sin(ph) + 0.06 * np.sin(2 * ph + 0.7)
    mag = mag - mag.mean() + 12.0
    mag = mag + rng.normal(0, noise, n)
    err = np.full(n, noise)
    return t, mag, err


def recover(t, mag, err):
    """Run the production period analysis; return (P_pdm, P_ls)."""
    res = run_period_analysis(
        time=t, mag_raw=mag, mag_corr=None, mag_err=err,
        min_period=0.04, max_period=0.30, samples_per_peak=12,
        methods=["ls", "pdm"], pdm_n_bins=10,
    )
    # keys are prefixed by the series label the service assigns ("raw")
    def pick(key):
        d = res.get(key) or {}
        v = d.get("best_period")
        return float(v) if v is not None and np.isfinite(v) else float("nan")
    return pick("raw_pdm"), pick("raw_ls")


def run_pdm_grid():
    rng = np.random.default_rng(SEED)
    rows = []
    for p in PERIODS:
        for s in NOISE:
            errs_pdm, errs_ls, pairs = [], [], []
            for _ in range(N_TRIAL):
                t, m, e = make_series(rng, p, s)
                p_pdm, p_ls = recover(t, m, e)
                if np.isfinite(p_pdm):
                    errs_pdm.append(abs(p_pdm - p) / p)
                if np.isfinite(p_ls):
                    errs_ls.append(abs(p_ls - p) / p)
                if np.isfinite(p_pdm) and np.isfinite(p_ls):
                    pairs.append((p_pdm, p_ls))
            rows.append({
                "period": p, "noise": s,
                "pdm_med_relerr": float(np.median(errs_pdm)) if errs_pdm else float("nan"),
                "ls_med_relerr": float(np.median(errs_ls)) if errs_ls else float("nan"),
                "n_ok_pdm": len(errs_pdm), "n_trial": N_TRIAL,
                "pairs": pairs,
            })
            print(f"  P={p:.6f} noise={s:.3f}  PDM {rows[-1]['pdm_med_relerr']*100:6.3f}%"
                  f"  LS {rows[-1]['ls_med_relerr']*100:6.3f}%", flush=True)
    return rows


def run_sysrem():
    """Common systematic + one real variable, run the way the pipeline runs it.

    APEX solves the components on the *comparison* stars only and then applies
    them to the target (`global_ensemble` drops the target id before solving;
    `apply_sysrem` fits only the target's own amplitude for each component).
    We also run the naive alternative — target included in the solve — to show
    what that design decision buys.
    """
    rng = np.random.default_rng(SEED + 1)
    n_star, n_frame = 60, 120
    t = np.sort(rng.uniform(0, 0.22, n_frame))

    sens = rng.normal(1.0, 0.35, n_star)[:, None]        # per-star sensitivity a_i
    transp = 0.08 * np.sin(2 * np.pi * t / 0.17) + 0.05 * (t - t.mean()) / np.ptp(t)
    systematic = sens * transp[None, :]                   # a_i * c_j

    var_idx = 7
    var_period, var_amp = 0.104092, 0.21
    signal = np.zeros(n_frame)
    signal = var_amp * np.sin(2 * np.pi * t / var_period)

    noise_rms = 0.008
    noise = rng.normal(0, noise_rms, (n_star, n_frame))
    err = np.full((n_star, n_frame), noise_rms)

    const = np.ones(n_star, bool)
    const[var_idx] = False

    # differential magnitudes (base level already removed by the ensemble step)
    comp_mag = systematic[const] + noise[const]
    tgt_mag = systematic[var_idx] + signal + noise[var_idx]

    # --- pipeline path: solve on comparisons, apply to target -------------
    res = sysrem(comp_mag, err[const], n_iter=3)
    comp_resid = np.asarray(res.residuals)
    tgt_corr = np.asarray(apply_sysrem(tgt_mag, err[var_idx], res))

    # --- naive path: target inside the solve ------------------------------
    all_mag = np.vstack([comp_mag, tgt_mag[None, :]])
    all_err = np.vstack([err[const], err[var_idx][None, :]])
    res_bad = sysrem(all_mag, all_err, n_iter=3)
    tgt_bad = np.asarray(res_bad.residuals)[-1]

    def amp(x):
        return float(np.ptp(np.asarray(x) - np.nanmedian(x)))

    return {
        "n_star": n_star, "n_comp": int(const.sum()), "n_frame": n_frame,
        "noise_rms": noise_rms,
        "rms_before": float(np.median(np.nanstd(comp_mag, axis=1))),
        "rms_after": float(np.median(np.nanstd(comp_resid, axis=1))),
        "amp_injected": amp(signal),
        "amp_pipeline": amp(tgt_corr),
        "amp_naive": amp(tgt_bad),
        "ratio_pipeline": amp(tgt_corr) / amp(signal),
        "ratio_naive": amp(tgt_bad) / amp(signal),
        "t": t.tolist(),
        "var_before": tgt_mag.tolist(),
        "var_after": (tgt_corr - float(np.nanmedian(tgt_corr))).tolist(),
        "var_naive": (tgt_bad - float(np.nanmedian(tgt_bad))).tolist(),
        "var_truth": signal.tolist(),
        "const_rms_before": np.nanstd(comp_mag, axis=1).tolist(),
        "const_rms_after": np.nanstd(comp_resid, axis=1).tolist(),
    }


def main() -> int:
    DATADIR.mkdir(parents=True, exist_ok=True)
    print("PDM grid ...")
    pdm_rows = run_pdm_grid()
    print("SYSREM ...")
    sr = run_sysrem()
    print(f"  비교성 RMS {sr['rms_before']*1000:.1f} -> {sr['rms_after']*1000:.1f} mmag")
    print(f"  변광 진폭 유지: 파이프라인 경로 {sr['ratio_pipeline']*100:.1f}% · "
          f"대상 포함(잘못된 경로) {sr['ratio_naive']*100:.1f}%")

    (DATADIR / "pdm_grid.json").write_text(
        json.dumps([{k: v for k, v in r.items() if k != "pairs"} for r in pdm_rows],
                   indent=1), encoding="utf-8")
    (DATADIR / "sysrem.json").write_text(json.dumps(sr, indent=1), encoding="utf-8")

    fig, axes = plt.subplots(1, 3, figsize=(DOUBLE_COL, 2.8))

    # (a) PDM relative error vs noise, one line per period
    ax = axes[0]
    line_styles = ["-", "--", ":", "-.", (0, (5, 1, 1, 1)), (0, (2, 1))]
    markers = ["o", "s", "^", "D", "v", "P"]
    greys = ["#111111", "#333333", "#555555", "#777777", "#999999", "#BBBBBB"]
    for i, p in enumerate(PERIODS):
        sub = [r for r in pdm_rows if r["period"] == p]
        ax.plot([r["noise"] * 1000 for r in sub],
                [r["pdm_med_relerr"] * 100 for r in sub],
                marker=markers[i], ls=line_styles[i], ms=3, color=greys[i],
                label=f"{p:g} d")
    ax.set_xscale("log"); ax.set_yscale("log")
    # default log ticks collide at 5-40; label the four sampled values instead
    ax.set_xticks([5, 10, 20, 40])
    ax.set_xticklabels(["5", "10", "20", "40"])
    ax.set_xticks([], minor=True)
    ax.set_xlim(4.2, 48)
    ax.set_ylim(0.28, 13)
    ax.set_yticks([0.5, 1, 2, 5, 10])
    ax.set_yticklabels(["0.5", "1", "2", "5", "10"])
    ax.set_yticks([], minor=True)
    ax.set_xlabel("injected noise (mmag)")
    ax.set_ylabel("PDM period error (per cent)")
    ax.set_title("(a) PDM recovers injected period", loc="left", fontsize=9)
    ax.legend(fontsize=5.6, ncol=3, title="injected $P$ (d)", title_fontsize=5.8,
              loc="upper center", handlelength=1.2, columnspacing=0.7,
              handletextpad=0.35, borderaxespad=0.25)

    # (b) PDM vs LS on the same series
    ax = axes[1]
    xs, ys = [], []
    for r in pdm_rows:
        for a, b in r["pairs"]:
            xs.append(b); ys.append(a)
    ax.plot([0.04, 0.30], [0.04, 0.30], "-", color=PALETTE["black"], lw=0.8, zorder=1)
    ax.plot(xs, ys, "o", ms=2.4, color=C["data"], alpha=0.45, zorder=2)
    ax.set_xlabel("Lomb-Scargle period (d)")
    ax.set_ylabel("PDM period (d)")
    ax.set_title(f"(b) against library LS, {len(xs)} series", loc="left", fontsize=9)

    # (c) SYSREM: systematic removed, variability kept — and what happens if
    #     the target is left inside the solve
    ax = axes[2]
    t = np.asarray(sr["t"])
    order = np.argsort(t)
    ax.plot(t * 24, sr["var_before"], "o", ms=2.0, mfc="white",
            mec=PALETTE["grey"], color=PALETTE["grey"],
            alpha=0.6, label="before")
    ax.plot(t[order] * 24, np.asarray(sr["var_truth"])[order], "-",
            color=C["model"], lw=1.3, zorder=4, label="injected")
    ax.plot(t * 24, sr["var_after"], "o", ms=2.4, color=C["data"],
            alpha=0.9, zorder=3,
            label=f"comparisons only ({sr['ratio_pipeline']*100:.0f} per cent)")
    ax.plot(t * 24, sr["var_naive"], "x", ms=2.6, color=PALETTE["purple"],
            alpha=0.8, mew=0.7,
            label=f"target in solve ({sr['ratio_naive']*100:.1f} per cent)")
    ax.set_xlabel("time (h)")
    ax.set_ylabel("differential magnitude")
    ax.invert_yaxis()
    ax.legend(fontsize=6.0, loc="lower right", handletextpad=0.4,
              frameon=True, framealpha=0.85, edgecolor="none",
              facecolor="white")
    ax.set_title("(c) SYSREM keeps the variable", loc="left", fontsize=9)

    prov = (
        f"Synthetic only, fixed seed {SEED}, no external data. Series on the Section-4 "
        f"observing pattern (one night, {SPAN_HOURS:g} h, {N_POINTS} points), analysed "
        "through the production LC entry points (run_period_analysis: PDM + astropy "
        "Lomb-Scargle; sysrem).\n"
        f"(a,b): P = {PERIODS[0]:g}-{PERIODS[-1]:g} d $\\times$ noise "
        f"{NOISE[0] * 1000:g}-{NOISE[-1] * 1000:g} mmag $\\times$ {N_TRIAL} realizations "
        f"= {len(PERIODS) * len(NOISE) * N_TRIAL} series; (c): 60 stars $\\times$ 120 "
        "frames, one injected variable (P = 0.104092 d, amplitude 0.21 mag)."
    )
    fig.tight_layout(pad=0.5, w_pad=1.6, rect=(0, 0.115, 1, 1))
    fig.text(0.005, 0.005, prov, fontsize=5.6, color=PALETTE["grey"], va="bottom")
    for ext, p in save_fig(fig, "fig_timeseries_validation", OUTDIR).items():
        print(f"[saved] {p}")

    ok = [r for r in pdm_rows if np.isfinite(r["pdm_med_relerr"])]
    worst = max(ok, key=lambda r: r["pdm_med_relerr"])
    print(f"\nPDM 최악 셀: P={worst['period']:g} d, noise={worst['noise']*1000:.0f} mmag "
          f"-> {worst['pdm_med_relerr']*100:.3f}%")
    print(f"PDM 20 mmag 이하 전체 중앙값: "
          f"{np.median([r['pdm_med_relerr'] for r in ok if r['noise'] <= 0.02])*100:.3f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
