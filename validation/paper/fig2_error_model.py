"""Figure 2 — validation of APEX's photometric error model.

Question: are APEX's *reported* magnitude uncertainties honest? A correctly
calibrated error model means the pull ((m_meas - m_true) / sigma_m) is drawn
from a unit Gaussian N(0, 1).

Design — a controlled Monte-Carlo experiment that mirrors the exact noise
realization used by ``apex/benchmark/synthetic_frame.py`` (gain, read noise,
sky, Poisson + Gaussian read noise, ZP, Gaussian PSF), then MEASURES every
injected star with APEX's REAL photometry routine
``apex.utils.photometry_utils.phot_vectorized``. Nothing about the recovered
flux, its error, or the SNR is invented here — those come straight out of APEX,
so the pull directly tests APEX's reported sigma.

We inject a magnitude ladder (14.5 -> 20.0, 0.5 steps), many independent seeded
realizations per magnitude, stars on a well-separated grid (separation >=
2*r_out) so blending never contaminates the aperture.

Panels (2x2):
  (a) bias  (m_meas - m_true) vs m_true
  (b) empirical scatter vs magnitude, with the CCD-equation prediction
  (c) pull distribution vs unit Gaussian
  (d) reliability: reported sigma vs empirical RMS, with y=x
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, SINGLE_COL, DOUBLE_COL
from apex.utils.photometry_utils import phot_vectorized

apply_paper_style()

# ── Detector / measurement constants (match synthetic_frame.py defaults) ─────
GAIN = 1.5          # e-/ADU
RN = 5.0            # read noise, e-
SKY_B = 150.0       # sky background, ADU
ZP = 25.0           # magnitude zero point
FWHM = 3.5          # Gaussian PSF FWHM, px

# Aperture geometry (state exact values in the caption).
R_AP = 1.0 * FWHM   # = 3.5 px
R_IN = 4.0 * FWHM   # = 14.0 px
R_OUT = 6.0 * FWHM  # = 21.0 px

SAT_ADU = 65000.0
DATAMAX_ADU = 60000.0

SIGMA_PSF = FWHM / 2.3548200450309493

# ── Experiment design ────────────────────────────────────────────────────────
MAG_LADDER = np.arange(14.5, 20.0 + 1e-9, 0.5)   # 14.5 .. 20.0
STARS_PER_FRAME = 36                              # 6x6 grid of isolated stars
# realizations chosen so recovered count per magnitude >> 300
N_REALIZATIONS = 12                               # -> 12*36 = 432 per magnitude
BASE_SEED = 20260702

# Grid geometry — separation >= 2*r_out so apertures/annuli never overlap.
GRID_N = 6                                         # 6x6 = 36 stars
SEP = 2.0 * R_OUT + 6.0                             # comfortably > 2*r_out
MARGIN = R_OUT + 6.0                                # keep annuli inside frame
FRAME_SIZE = int(np.ceil(MARGIN * 2 + (GRID_N - 1) * SEP)) + 2


def _grid_positions() -> np.ndarray:
    """(N,2) col-major (x,y) grid of well-separated star centres."""
    xs = MARGIN + np.arange(GRID_N) * SEP
    ys = MARGIN + np.arange(GRID_N) * SEP
    gx, gy = np.meshgrid(xs, ys)
    return np.column_stack([gx.ravel(), gy.ravel()])


def _gaussian_stamp(sigma: float, dx: float, dy: float, half: int) -> np.ndarray:
    """Unit-sum Gaussian PSF stamp shifted by (dx, dy), matching synthetic_frame."""
    yy, xx = np.mgrid[-half : half + 1, -half : half + 1]
    stamp = np.exp(-((xx - dx) ** 2 + (yy - dy) ** 2) / (2.0 * sigma * sigma))
    total = float(stamp.sum())
    return stamp / total if total > 0 else stamp


def _measure_aperture_correction() -> tuple[float, float]:
    """Fraction of PSF flux inside r_ap, MEASURED with APEX's own aperture.

    A finite aperture of r_ap = 1.0*FWHM captures only ~94% of a Gaussian PSF;
    the missing wings make m_meas fainter by a fixed, deterministic offset (the
    classic aperture correction). This is NOT part of the random noise budget,
    so APEX's per-star sigma correctly excludes it — but left uncorrected it
    would contaminate the pull with a constant bias. We remove it by measuring
    the flux fraction on a *noiseless, high-flux* frame with the exact same
    aperture + sky annulus APEX uses, then defining the truth as the
    aperture-corrected magnitude. Returns (flux_fraction, apcorr_mag).
    """
    Fe = 1.0e7  # bright + noiseless -> pure geometric flux fraction
    half = int(np.ceil(4.0 * FWHM))
    size = 100
    yy, xx = np.mgrid[-half : half + 1, -half : half + 1]
    st = np.exp(-((xx) ** 2 + (yy) ** 2) / (2.0 * SIGMA_PSF ** 2))
    st = st / st.sum()
    img_e = np.zeros((size, size))
    img_e[50 - half : 50 + half + 1, 50 - half : 50 + half + 1] += st * Fe
    img_e += SKY_B * GAIN
    flux_e, *_ = phot_vectorized(
        img_e / GAIN, np.array([[50.0, 50.0]]), R_AP, R_IN, R_OUT,
        gain=GAIN, rn_param_e=RN, sky_frame_e=np.nan, sky_sigma_mode="local",
        sat_adu=1e12, datamax_adu=1e12,
    )
    frac = float(flux_e[0]) / Fe
    apcorr = -2.5 * np.log10(frac)
    return frac, apcorr


def _inject_frame(mag_true: float, positions: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Build one realized ADU frame with all stars at ``mag_true``.

    Mirrors apex/benchmark/synthetic_frame.py exactly:
      signal_e   = sum unit-sum PSF * F_e
      expected_e = signal_e + B*gain
      realized_e = Poisson(expected_e) + Normal(0, RN)
      image_adu  = realized_e / gain
    """
    flux_e = 10.0 ** (-0.4 * (mag_true - ZP))
    half = int(np.ceil(4.0 * FWHM))
    signal_e = np.zeros((FRAME_SIZE, FRAME_SIZE), dtype=np.float64)
    for (x, y) in positions:
        xi = int(round(x))
        yi = int(round(y))
        dx = x - xi
        dy = y - yi
        stamp = _gaussian_stamp(SIGMA_PSF, dx, dy, half) * flux_e
        signal_e[yi - half : yi + half + 1, xi - half : xi + half + 1] += stamp

    background_e = SKY_B * GAIN
    expected_e = signal_e + background_e
    realized_e = rng.poisson(np.clip(expected_e, 0.0, None)).astype(np.float64)
    realized_e += rng.normal(0.0, RN, size=realized_e.shape)
    image_adu = realized_e / GAIN
    # Do NOT clip to saturation here: at these magnitudes peaks stay well below
    # SAT, and clipping would only bias bright pixels. (synthetic_frame clips at
    # 65000 ADU; our brightest star peak << that — verified in run.)
    return image_adu


def run_experiment():
    positions = _grid_positions()
    n_pix_ap = np.pi * R_AP ** 2

    flux_frac, apcorr = _measure_aperture_correction()
    print(f"[fig2] measured aperture flux fraction = {flux_frac:.5f} "
          f"(r_ap = {R_AP:.2f} px = 1.0*FWHM) -> apcorr = {apcorr:+.5f} mag")

    records = {
        "m_true": [], "m_meas": [], "sigma_m": [], "snr": [],
    }

    for mi, mag_true in enumerate(MAG_LADDER):
        for r in range(N_REALIZATIONS):
            seed = BASE_SEED + mi * 1000 + r
            rng = np.random.default_rng(seed)
            image_adu = _inject_frame(mag_true, positions, rng)

            flux_e, flux_err_e, snr, sky_med, is_sat, is_nonlin = phot_vectorized(
                image_adu, positions, R_AP, R_IN, R_OUT,
                gain=GAIN, rn_param_e=RN, sky_frame_e=np.nan,
                sky_sigma_mode="local", sat_adu=SAT_ADU, datamax_adu=DATAMAX_ADU,
            )

            good = (
                np.isfinite(flux_e) & (flux_e > 0)
                & np.isfinite(flux_err_e) & (flux_err_e > 0)
                & np.isfinite(snr) & ~is_sat & ~is_nonlin
            )
            if not good.any():
                continue
            fe = flux_e[good]
            fee = flux_err_e[good]
            # Aperture-corrected measured magnitude (removes the deterministic
            # finite-aperture flux loss; noise budget is unaffected).
            m_meas = ZP - 2.5 * np.log10(fe) - apcorr
            sigma_m = (2.5 / np.log(10.0)) * fee / fe   # = 1.0857 / snr

            records["m_true"].append(np.full(good.sum(), mag_true))
            records["m_meas"].append(m_meas)
            records["sigma_m"].append(sigma_m)
            records["snr"].append(snr[good])

    out = {k: np.concatenate(v) for k, v in records.items()}
    out["dmag"] = out["m_meas"] - out["m_true"]
    out["pull"] = out["dmag"] / out["sigma_m"]
    out["n_pix_ap"] = n_pix_ap
    out["flux_frac"] = flux_frac
    out["apcorr"] = apcorr
    return out


def _running_median(x, y, nbins=40):
    order = np.argsort(x)
    xs, ys = x[order], y[order]
    edges = np.linspace(xs.min(), xs.max(), nbins + 1)
    idx = np.clip(np.digitize(xs, edges) - 1, 0, nbins - 1)
    mx, my = [], []
    for b in range(nbins):
        sel = idx == b
        if sel.sum() >= 5:
            mx.append(np.median(xs[sel]))
            my.append(np.median(ys[sel]))
    return np.array(mx), np.array(my)


def _per_mag_stats(res):
    mags = np.unique(res["m_true"])
    emp_rms, rep_sigma, mean_bias = [], [], []
    for m in mags:
        sel = res["m_true"] == m
        emp_rms.append(np.sqrt(np.mean(res["dmag"][sel] ** 2)))
        rep_sigma.append(np.median(res["sigma_m"][sel]))
        mean_bias.append(np.mean(res["dmag"][sel]))
    return mags, np.array(emp_rms), np.array(rep_sigma), np.array(mean_bias)


def _ccd_sigma_theory(mags, n_pix_ap):
    """CCD-equation magnitude uncertainty for a source of magnitude m."""
    F_e = 10.0 ** (-0.4 * (mags - ZP))
    snr = F_e / np.sqrt(F_e + n_pix_ap * (SKY_B * GAIN + RN ** 2))
    return 1.0857 / snr, snr


def make_figure(res):
    mags, emp_rms, rep_sigma, mean_bias = _per_mag_stats(res)
    fine_mags = np.linspace(MAG_LADDER.min(), MAG_LADDER.max(), 200)
    sigma_theory, _ = _ccd_sigma_theory(fine_mags, res["n_pix_ap"])

    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL, 5.4))
    (ax_a, ax_b), (ax_c, ax_d) = axes

    # ── (a) Bias ─────────────────────────────────────────────────────────────
    # faintest ladder step whose mean bias is still below 5 mmag
    small = np.abs(mean_bias) <= 0.005
    m_zero = float(mags[np.max(np.nonzero(small))]) if small.any() else float("nan")
    ax_a.axhline(0.0, color=C["floor"], lw=1.0, ls="--", zorder=1)
    ax_a.scatter(res["m_true"], res["dmag"], s=3, alpha=0.06,
                 color=C["data"], edgecolors="none", rasterized=True, zorder=2)
    rmx, rmy = _running_median(res["m_true"], res["dmag"])
    ax_a.plot(rmx, rmy, color=C["model"], lw=1.6, zorder=4, label="running median")
    ax_a.set_xlabel(r"$m_{\rm true}$ (mag)")
    ax_a.set_ylabel(r"$m_{\rm meas}-m_{\rm true}$ (mag)")
    lim = np.percentile(np.abs(res["dmag"]), 99.0)
    n_out = int((np.abs(res["dmag"]) > lim).sum())
    ax_a.set_ylim(-lim, lim)
    # honesty: the axis range hides the extreme faint-end tail — say so
    ax_a.text(0.02, 0.03,
              f"{n_out:,}/{res['dmag'].size:,} points beyond $\\pm${lim:.2f} "
              "mag not shown",
              transform=ax_a.transAxes, fontsize=6.0, color=PALETTE["grey"])
    ax_a.annotate(f"median {rmy[-1] * 1000:+.0f} mmag\nat m = {rmx[-1]:.1f}",
                  (rmx[-1], rmy[-1]), textcoords="offset points",
                  xytext=(-64, 26), fontsize=6.4, color=C["model"],
                  arrowprops=dict(arrowstyle="-", lw=0.7, color=C["model"]))
    ax_a.set_title(f"(a) bias consistent with zero to m $\\approx$ {m_zero:.1f}",
                   loc="left", fontsize=9)
    ax_a.legend(loc="upper left", fontsize=6.8)

    # ── (b) Scatter vs magnitude ─────────────────────────────────────────────
    decades = np.log10(emp_rms.max() / emp_rms.min())
    ax_b.plot(fine_mags, sigma_theory, color=C["model"], lw=1.6,
              label=r"CCD eq. $\sigma_m = 1.0857/\mathrm{SNR}$", zorder=3)
    ax_b.scatter(mags, emp_rms, s=22, color=C["data"], zorder=4,
                 label="empirical RMS", edgecolors="white", linewidths=0.5)
    ax_b.set_yscale("log")
    ax_b.set_xlabel(r"$m_{\rm true}$ (mag)")
    ax_b.set_ylabel(r"$\sigma_m$ (mag)")
    ax_b.annotate(f"{emp_rms[0] * 1000:.1f} mmag", (mags[0], emp_rms[0]),
                  textcoords="offset points", xytext=(6, -3), fontsize=6.4,
                  color=C["data"])
    ax_b.annotate(f"{emp_rms[-1]:.2f} mag", (mags[-1], emp_rms[-1]),
                  textcoords="offset points", xytext=(-42, 2), fontsize=6.4,
                  color=C["data"])
    ax_b.set_title(f"(b) scatter follows the CCD equation over "
                   f"{decades:.1f} dex", loc="left", fontsize=9)
    ax_b.legend(loc="upper left", fontsize=6.8)

    # ── (c) Pull distribution (SNR > 5) ──────────────────────────────────────
    sel = res["snr"] > 5.0
    pull = res["pull"][sel]
    p_mean, p_std, p_n = float(np.mean(pull)), float(np.std(pull, ddof=1)), int(pull.size)
    n_pull_out = int((np.abs(pull) > 5).sum())
    bins = np.linspace(-5, 5, 61)
    ax_c.hist(pull, bins=bins, density=True, color=C["data"], alpha=0.55,
              edgecolor="white", linewidth=0.3, label="APEX pull")
    xg = np.linspace(-5, 5, 400)
    ax_c.plot(xg, np.exp(-xg ** 2 / 2) / np.sqrt(2 * np.pi),
              color=C["model"], lw=1.6, label=r"$\mathcal{N}(0,1)$")
    ax_c.set_xlabel(r"pull $=(m_{\rm meas}-m_{\rm true})/\sigma_m$  (SNR $>$ 5)")
    ax_c.set_ylabel("density")
    ax_c.set_xlim(-5, 5)
    ax_c.text(0.03, 0.97,
              f"mean = {p_mean:+.3f}\nstd = {p_std:.3f}\nN = {p_n:,}"
              + (f"\n{n_pull_out} beyond $\\pm$5" if n_pull_out else ""),
              transform=ax_c.transAxes, va="top", ha="left", fontsize=7.6,
              bbox=dict(boxstyle="round,pad=0.35", fc="white", ec=C["floor"], lw=0.6))
    ax_c.set_title(f"(c) reported $\\sigma$ is honest: pull std "
                   f"{p_std:.3f}", loc="left", fontsize=9)
    ax_c.legend(loc="upper right", fontsize=6.8)

    # ── (d) Reliability: reported sigma vs empirical RMS ─────────────────────
    ratio = rep_sigma / emp_rms
    med_ratio = float(np.median(ratio))
    worst_i = int(np.argmin(ratio))
    lo = min(rep_sigma.min(), emp_rms.min()) * 0.8
    hi = max(rep_sigma.max(), emp_rms.max()) * 1.2
    ax_d.plot([lo, hi], [lo, hi], color=C["floor"], lw=1.0, ls="--",
              label=r"$y=x$", zorder=1)
    sc = ax_d.scatter(rep_sigma, emp_rms, c=mags, cmap="viridis", s=28,
                      zorder=3, edgecolors="white", linewidths=0.5)
    cb = fig.colorbar(sc, ax=ax_d, pad=0.02, fraction=0.05)
    cb.set_label(r"$m_{\rm true}$ (mag)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ax_d.annotate(
        f"m = {mags[worst_i]:.1f}: reported\n{(1 - ratio[worst_i]) * 100:.0f}% "
        "below empirical",
        (rep_sigma[worst_i], emp_rms[worst_i]), textcoords="offset points",
        xytext=(-104, -30), fontsize=6.4, color=C["model"],
        arrowprops=dict(arrowstyle="-", lw=0.7, color=C["model"]))
    ax_d.set_xscale("log")
    ax_d.set_yscale("log")
    ax_d.set_xlim(lo, hi)
    ax_d.set_ylim(lo, hi)
    ax_d.set_aspect("equal", adjustable="box")
    ax_d.set_xlabel(r"APEX-reported median $\sigma_m$ (mag)")
    ax_d.set_ylabel(r"empirical RMS (mag)")
    ax_d.set_title(f"(d) reported vs empirical: median ratio "
                   f"{med_ratio:.2f}", loc="left", fontsize=9)
    ax_d.legend(loc="upper left", fontsize=6.8)

    # ── in-figure data provenance ────────────────────────────────────────────
    prov = (
        "Synthetic Monte-Carlo, no external data: noise model of "
        "apex.benchmark.synthetic_frame (gain 1.5 e$^-$/ADU · RN 5 e$^-$ · sky 150 ADU · "
        f"ZP 25 · Gaussian PSF FWHM {FWHM:g} px), measured with APEX's production "
        "photometry (apex.utils.photometry_utils.phot_vectorized; "
        f"$r_{{\\rm ap}}$ {R_AP:g} px = 1.0$\\times$FWHM, annulus {R_IN:g}-{R_OUT:g} px).\n"
        f"Ladder m = {MAG_LADDER.min():.1f}-{MAG_LADDER.max():.1f} in 0.5-mag steps "
        f"$\\times$ {STARS_PER_FRAME} isolated grid stars $\\times$ {N_REALIZATIONS} "
        f"noise realizations, seed {BASE_SEED}; aperture correction "
        f"{res['apcorr']:+.4f} mag measured on a noiseless frame with the same "
        "aperture."
    )
    fig.tight_layout(pad=0.6, rect=(0, 0.055, 1, 1))
    fig.text(0.005, 0.005, prov, fontsize=5.6, color=PALETTE["grey"], va="bottom")

    stats = {
        "pull_mean": p_mean, "pull_std": p_std, "pull_n": p_n,
        "mag_min": float(MAG_LADDER.min()), "mag_max": float(MAG_LADDER.max()),
        "mags": mags, "emp_rms": emp_rms, "rep_sigma": rep_sigma,
        "mean_bias": mean_bias,
        "n_total": int(res["m_true"].size),
        "flux_frac": float(res["flux_frac"]), "apcorr": float(res["apcorr"]),
        "m_zero_bias": m_zero, "median_ratio": med_ratio,
        "worst_mag": float(mags[worst_i]),
        "worst_deficit_pct": float((1 - ratio[worst_i]) * 100),
    }
    return fig, stats


def main():
    print("[fig2] running Monte-Carlo error-model experiment ...")
    res = run_experiment()
    print(f"[fig2] total recovered measurements: {res['m_true'].size}")
    for m in np.unique(res["m_true"]):
        n = int((res["m_true"] == m).sum())
        print(f"        m_true={m:.1f}  N={n}")

    fig, stats = make_figure(res)
    outdir = REPO / "validation" / "paper" / "figures"
    paths = save_fig(fig, "fig2_error_model", outdir)
    print(f"[fig2] wrote {paths['pdf']}")
    print(f"[fig2] wrote {paths['png']}")

    # Report the key numbers.
    print("\n=== KEY NUMBERS ===")
    print(f"aperture flux fraction = {stats['flux_frac']:.5f}  "
          f"(apcorr = {stats['apcorr']:+.5f} mag)")
    print(f"pull mean (SNR>5) = {stats['pull_mean']:+.4f}")
    print(f"pull std  (SNR>5) = {stats['pull_std']:.4f}")
    print(f"pull N    (SNR>5) = {stats['pull_n']}")
    print(f"mag range = {stats['mag_min']:.1f} .. {stats['mag_max']:.1f}  "
          f"(N_total = {stats['n_total']})")
    print(f"bias consistent with zero (|bias| <= 5 mmag) to m = {stats['m_zero_bias']:.1f}")
    print(f"median reported/empirical = {stats['median_ratio']:.3f}; worst bin "
          f"m={stats['worst_mag']:.1f} reported {stats['worst_deficit_pct']:.0f}% low")
    ratio = stats["rep_sigma"] / stats["emp_rms"]
    print("per-mag reported/empirical sigma ratio (panel d):")
    for m, rr in zip(stats["mags"], ratio):
        print(f"    m={m:.1f}  reported/empirical = {rr:.3f}")
    print(f"median reported/empirical ratio = {np.median(ratio):.3f}")

    # Pull std stratified by SNR (the honest test: std should track 1).
    print("pull std by SNR band:")
    for lo, hi in [(5, 10), (10, 20), (20, 50), (50, 1e9)]:
        m = (res["snr"] >= lo) & (res["snr"] < hi)
        if m.sum() > 20:
            print(f"    {lo:>3}<=SNR<{hi:<5.0f}  mean={res['pull'][m].mean():+.3f}"
                  f"  std={res['pull'][m].std(ddof=1):.3f}  N={int(m.sum())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
