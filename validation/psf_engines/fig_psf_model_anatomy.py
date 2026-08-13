"""Why DAOPHOT's PSF is not "an ePSF by another name".

Both engines end up with a grid of numbers, and both use that grid to correct a
star's fit, so it is fair to ask whether the hybrid is just an ePSF with extra
steps. It is not, and the difference is where the *noise* sits.

An ePSF is empirical all the way down: every pixel of the model, including the
bright core that dominates a fit, is an average of real star pixels and
therefore carries those stars' noise. DAOPHOT splits the model in two. An
analytic Moffat — five numbers fitted to sixty stars — carries the bulk of the
light and is smooth by construction, so it injects essentially no noise. Only
the leftover goes into a grid, and that grid is a small correction.

The panels show the real arrays from one M13 frame: APEX's ePSF as written by
step 8, and DAOPHOT's analytic parameters and residual look-up table as written
by `psf`. The analytic panel is reconstructed from the header parameters
(PAR1, PAR2, PAR3, beta) using DAOPHOT's Moffat form; it is a rendering of what
those numbers mean, not a file read off disk, and it is labelled as such.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
from astropy.io import fits

HERE = Path(__file__).absolute().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "validation" / "paper"))

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import SymLogNorm  # noqa: E402

try:
    from apex_paper_style import DOUBLE_COL, PALETTE, apply_paper_style, save_fig
    _STYLE = True
except Exception:  # pragma: no cover
    _STYLE = False
    PALETTE = {"black": "#111", "grey": "#777", "blue": "#0072B2",
               "vermillion": "#D55E00"}
    DOUBLE_COL = 7.0


def daophot_moffat(shape: tuple[int, int], par1: float, par2: float,
                   par3: float, beta: float, oversample: float) -> np.ndarray:
    """Render DAOPHOT's analytic Moffat from its stored parameters.

    DAOPHOT writes the widths in native pixels while the look-up table is
    stored oversampled, so the grid coordinates are divided down before the
    profile is evaluated.
    """
    ny, nx = shape
    y, x = np.mgrid[:ny, :nx].astype(float)
    x = (x - (nx - 1) / 2.0) / oversample
    y = (y - (ny - 1) / 2.0) / oversample
    z = (x / par1) ** 2 + (y / par2) ** 2 + x * y * par3
    return np.power(1.0 + np.clip(z, 0.0, None), -beta)


def radial_profile(a: np.ndarray, oversample: float = 1.0) -> tuple[np.ndarray, np.ndarray]:
    ny, nx = a.shape
    y, x = np.mgrid[:ny, :nx].astype(float)
    r = np.hypot(x - (nx - 1) / 2.0, y - (ny - 1) / 2.0) / oversample
    edges = np.arange(0, r.max(), 0.5)
    idx = np.digitize(r.ravel(), edges) - 1
    prof = np.array([np.nanmean(a.ravel()[idx == k]) if np.any(idx == k) else np.nan
                     for k in range(len(edges))])
    return edges + 0.25, prof


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--daophot-psf", default=(
        r"E:\APEX_validation\psf_engines\M13_ast2x200_moffat"
        r"\daophot_work_trial1\psf.fits"))
    ap.add_argument("--apex-epsf", default=(
        r"E:\APEX_validation\psf_engines\M13_ast2x200_moffat\trial_0001"
        r"\result\cmd_psf\epsf_model_B_pp_messier13-0005-B.fits"))
    ap.add_argument("--outdir", default=str(HERE))
    args = ap.parse_args()

    with fits.open(args.daophot_psf) as h:
        dh, look = h[0].header, h[0].data.astype(float)
    with fits.open(args.apex_epsf) as h:
        eh, epsf = h[0].header, h[0].data.astype(float)

    beta = float(dh.get("PAR4", 1.5))
    par1, par2, par3 = (float(dh["PAR1"]), float(dh["PAR2"]),
                        float(dh.get("PAR3", 0.0)))
    peak = float(dh["PSFHEIGH"])
    # The table spans 2*psfrad+1 native pixels; recover its sampling.
    os_dao = (look.shape[0] - 1) / (2.0 * float(dh["PSFRAD"]))
    os_apex = float(eh.get("OVERSAMPL", 2))

    analytic = daophot_moffat(look.shape, par1, par2, par3, beta, os_dao) * peak
    total = analytic + look

    # What fraction of the model's light comes from the grid rather than the
    # smooth part — the quantity that decides how much grid noise a fit inherits.
    frac_grid = abs(look.sum()) / (abs(analytic.sum()) + abs(look.sum()))

    if _STYLE:
        apply_paper_style()
    fig = plt.figure(figsize=(DOUBLE_COL, 4.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.75], hspace=0.45, wspace=0.28)

    def show(ax, a, title, *, symlog=False):
        v = np.nanpercentile(np.abs(a), 99.5)
        norm = SymLogNorm(linthresh=max(v * 1e-3, 1e-6), vmin=-v, vmax=v) if symlog else None
        ax.imshow(a, origin="lower", cmap="gray" if not symlog else "coolwarm",
                  norm=norm, vmin=None if symlog else 0, vmax=None if symlog else v)
        ax.set_title(title, fontsize=7.5)
        ax.set_xticks([]); ax.set_yticks([])

    show(fig.add_subplot(gs[0, 0]), analytic,
         f"DAOPHOT analytic (Moffat $\\beta$={beta:g})\n"
         "rebuilt from 5 header numbers — smooth")
    show(fig.add_subplot(gs[0, 1]), look,
         "DAOPHOT residual look-up table (file)\n"
         "what the analytic part missed", symlog=True)
    show(fig.add_subplot(gs[0, 2]), epsf,
         "APEX ePSF (file)\nempirical everywhere, core included")

    ax = fig.add_subplot(gs[1, :])
    r_a, p_a = radial_profile(analytic, os_dao)
    r_l, p_l = radial_profile(np.abs(look), os_dao)
    r_e, p_e = radial_profile(epsf, os_apex)
    ax.plot(r_a, p_a / np.nanmax(p_a), lw=1.6, color=PALETTE["black"],
            label="DAOPHOT analytic (noise-free)")
    ax.plot(r_l, p_l / np.nanmax(p_a), lw=1.4, color=PALETTE["vermillion"],
            label="DAOPHOT residual |value|")
    ax.plot(r_e, p_e / np.nanmax(p_e), lw=1.4, ls="--", color=PALETTE["blue"],
            label="APEX ePSF (all empirical)")
    ax.set_yscale("log"); ax.set_ylim(1e-5, 2)
    ax.set_xlim(0, 30)
    ax.set_xlabel("radius (native pixel)")
    ax.set_ylabel("normalised value")
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    ax.set_title(
        f"grid-carried light: DAOPHOT {frac_grid*100:.1f} %  vs  APEX ePSF 100 %  "
        f"— how much grid noise a fit inherits", fontsize=7.5)

    fig.text(0.5, 0.005,
             "Moravian C3-61000 · M13 · pp_messier13-0005-B (injected frame) · "
             "DAOPHOT psf.fits / APEX epsf_model — real products",
             ha="center", fontsize=6, color=PALETTE["grey"])
    fig.tight_layout(rect=(0, 0.03, 1, 1))

    outdir = Path(args.outdir)
    if _STYLE:
        for k, p in save_fig(fig, "fig_psf_model_anatomy", outdir).items():
            print(f"[{k}] {p}")
    else:
        p = outdir / "fig_psf_model_anatomy.png"
        fig.savefig(p, dpi=150); print(f"[png] {p}")

    print(f"\n해석부 적분 {analytic.sum():,.0f} · 잔차표 순합 {look.sum():+,.0f}")
    print(f"잔차표가 차지하는 비율 {frac_grid*100:.1f} %  (ePSF 는 정의상 100 %)")
    print(f"오버샘플링: DAOPHOT {os_dao:.2f}x · APEX {os_apex:.0f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
