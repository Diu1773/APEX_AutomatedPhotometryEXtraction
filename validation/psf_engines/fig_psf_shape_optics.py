"""What the two optical systems' PSFs actually look like.

A Corrected Dall-Kirkham exists to do one thing: keep the image round and the
field flat, which a plain Dall-Kirkham does not — its ellipsoidal primary and
spherical secondary leave strong coma a short way off axis. So the question
"is the CDK's signature reproduced in the data?" has a testable form: does the
PSF stay round all the way to the corner, and is any elongation *not* aligned
with the radius from field centre? Coma stretches a star radially; seeing,
tracking and focus do not.

What cannot be seen here is worth stating first. A 508 mm aperture is
diffraction-limited near 0.27 arcsec at 550 nm, while these stars have a FWHM
of 2.85 arcsec. The image is ~10x wider than the diffraction core, so the
Airy pattern, the central-obstruction ring structure and the spider spikes are
all buried under the atmosphere. Nothing in this figure is a diffraction
measurement; it is the delivered PSF, which is what photometry actually fits.

Two systems at nearly the same sampling make the comparison fair: 0.393
arcsec/px on the CDK and 0.390 on the LCO 1 m, so the pixel grid is not doing
the talking.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.optimize import curve_fit

HERE = Path(__file__).absolute().parent
REPO = HERE.parents[1]
sys.path.insert(0, str(REPO / "validation" / "paper"))

import matplotlib.pyplot as plt  # noqa: E402
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401,E402

try:
    from apex_paper_style import C, DOUBLE_COL, PALETTE, apply_paper_style, save_fig
    _STYLE = True
except Exception:  # pragma: no cover
    _STYLE = False
    PALETTE = {"black": "#111", "grey": "#777", "blue": "#0072B2",
               "vermillion": "#D55E00", "green": "#009E73"}
    C = {"data": PALETTE["blue"], "model": PALETTE["vermillion"],
         "reference": PALETTE["green"], "floor": PALETTE["grey"]}
    DOUBLE_COL = 7.0


def load_epsf(path: Path) -> tuple[np.ndarray, int]:
    with fits.open(path) as h:
        return np.asarray(h[0].data, dtype=float), int(h[0].header.get("OVERSAMPL", 1))


def radial_profile(psf: np.ndarray, oversampling: int) -> tuple[np.ndarray, np.ndarray]:
    ny, nx = psf.shape
    yy, xx = np.mgrid[:ny, :nx].astype(float)
    r = np.hypot(xx - (nx - 1) / 2.0, yy - (ny - 1) / 2.0) / oversampling
    edges = np.arange(0.0, r.max(), 0.25)
    prof = np.array([np.nanmean(psf[(r >= a) & (r < a + 0.25)])
                     if np.any((r >= a) & (r < a + 0.25)) else np.nan
                     for a in edges])
    return edges + 0.125, prof


def moffat_beta(radius: np.ndarray, profile: np.ndarray) -> tuple[float, float]:
    ok = np.isfinite(profile) & (radius < 12.0)

    def model(r, amplitude, gamma, beta):
        return amplitude * (1.0 + (r / gamma) ** 2) ** (-beta)

    popt, _ = curve_fit(model, radius[ok], profile[ok],
                        p0=[np.nanmax(profile), 5.0, 3.0], maxfev=40000)
    return float(popt[1]), float(popt[2])


def profile_fwhm(radius: np.ndarray, profile: np.ndarray) -> float:
    peak = np.nanmax(profile)
    below = np.where(profile < peak / 2.0)[0]
    return float(2.0 * radius[below[0]]) if len(below) else float("nan")


def surface(ax, psf: np.ndarray, oversampling: int, scale: float, title: str,
            colour: str) -> None:
    """One PSF as a 3-D surface, trimmed to the part that carries the light."""
    ny, nx = psf.shape
    half = min(ny, nx) // 2
    keep = int(min(half, round(3.0 * oversampling * 3)))
    cy, cx = ny // 2, nx // 2
    cut = psf[cy - keep:cy + keep + 1, cx - keep:cx + keep + 1]
    cut = cut / cut.max()
    n = cut.shape[0]
    axis = (np.arange(n) - (n - 1) / 2.0) / oversampling * scale
    xx, yy = np.meshgrid(axis, axis)
    ax.plot_surface(xx, yy, cut, cmap="magma", linewidth=0, antialiased=True,
                    rcount=60, ccount=60, alpha=0.97)
    ax.contour(xx, yy, cut, levels=[0.05, 0.25, 0.5], colors=colour,
               linewidths=0.6, offset=0.0)
    ax.set_title(title, fontsize=7.5, pad=1)
    ax.set_xlabel("arcsec", fontsize=6, labelpad=-6)
    ax.set_ylabel("arcsec", fontsize=6, labelpad=-6)
    ax.set_zlabel("peak-normalised", fontsize=6, labelpad=-6)
    ax.tick_params(labelsize=5.5, pad=-2)
    ax.set_zlim(0, 1.0)
    ax.view_init(elev=32, azim=-125)
    ax.set_box_aspect((1, 1, 0.62))


def whisker_panel(ax, table: pd.DataFrame, shape: tuple[int, int], title: str,
                  colour: str) -> None:
    """Ellipticity as short sticks: length is elongation, angle is its axis.

    Coma from an uncorrected Dall-Kirkham would point every stick at the field
    centre and lengthen them outwards. A flat, randomly oriented field means
    the corrector is doing its job.
    """
    x = table["x"].to_numpy(float)
    y = table["y"].to_numpy(float)
    e = table["ell"].to_numpy(float)
    pa = np.radians(table["pa"].to_numpy(float))
    length = e / max(0.08, np.nanpercentile(e, 95)) * (shape[1] * 0.045)
    ax.quiver(x, y, length * np.cos(pa), length * np.sin(pa),
              angles="xy", scale_units="xy", scale=1, width=0.0035,
              headwidth=0, headlength=0, headaxislength=0,
              color=colour, alpha=0.85)
    ax.quiver(x, y, -length * np.cos(pa), -length * np.sin(pa),
              angles="xy", scale_units="xy", scale=1, width=0.0035,
              headwidth=0, headlength=0, headaxislength=0,
              color=colour, alpha=0.85)
    ax.set_xlim(0, shape[1]); ax.set_ylim(0, shape[0])
    ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title(title, fontsize=7.5, pad=2)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cdk-epsf", default=(
        r"E:\APEX_validation\psf_engines\M13_ast2x200_moffat\apexfix_trial1"
        r"\result\cmd_psf\epsf_model_B_pp_messier13-0005-B.fits"))
    ap.add_argument("--rc-epsf", default=(
        r"E:\APEX_validation\psf_engines\NGC5985_sinistro\trial_0001"
        r"\result\cmd_psf\epsf_model_rp_pp_ngc5985-0001-rp.fits"))
    ap.add_argument("--cdk-field", default=str(REPO / "scratchpad" / "psf_field_cdk.csv"))
    ap.add_argument("--rc-field", default=str(REPO / "scratchpad" / "psf_field_rc.csv"))
    ap.add_argument("--cdk-shape", type=int, nargs=2, default=[3194, 4788])
    ap.add_argument("--rc-shape", type=int, nargs=2, default=[4096, 4096])
    ap.add_argument("--outdir", default=str(HERE))
    args = ap.parse_args()

    cdk, cdk_os = load_epsf(Path(args.cdk_epsf))
    rc, rc_os = load_epsf(Path(args.rc_epsf))
    cdk_scale, rc_scale = 0.393, 0.390

    if _STYLE:
        apply_paper_style()
    # The paper style uses a serif face with no Hangul; this figure's labels are
    # Korean, so fall back to a face that has the glyphs rather than emit boxes.
    for face in ("Malgun Gothic", "Batang", "Gulim", "NanumGothic"):
        try:
            import matplotlib.font_manager as fm
            fm.findfont(fm.FontProperties(family=face), fallback_to_default=False)
        except Exception:
            continue
        plt.rcParams["font.family"] = face
        plt.rcParams["axes.unicode_minus"] = False
        break
    fig = plt.figure(figsize=(DOUBLE_COL, 6.4))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.15, 1.0, 0.92],
                          hspace=0.42, wspace=0.24)

    surface(fig.add_subplot(gs[0, 0], projection="3d"), cdk, cdk_os, cdk_scale,
            "CDK 508 mm f/7.8 — M13, B", C["data"])
    surface(fig.add_subplot(gs[0, 1], projection="3d"), rc, rc_os, rc_scale,
            "LCO 1 m (RC 계열) — NGC 5985, rp", C["reference"])

    ax = fig.add_subplot(gs[1, :])
    for psf, os_, scale, colour, label in (
            (cdk, cdk_os, cdk_scale, C["data"], "CDK 508 mm"),
            (rc, rc_os, rc_scale, C["reference"], "LCO 1 m")):
        radius, prof = radial_profile(psf / psf.sum(), os_)
        gamma, beta = moffat_beta(radius, prof)
        fwhm = profile_fwhm(radius, prof)
        ax.plot(radius * scale, prof / np.nanmax(prof), lw=1.5, color=colour,
                label=f"{label} · FWHM {fwhm * scale:.2f}″ · Moffat β={beta:.2f}")
    ax.set_yscale("log"); ax.set_ylim(3e-5, 1.6); ax.set_xlim(0, 6.0)
    # The Hangul face has no mathtext minus, so the default 10^-n labels come out
    # as boxes; plain decimals say the same thing in glyphs it has.
    ax.set_yticks([1.0, 1e-1, 1e-2, 1e-3, 1e-4])
    ax.set_yticklabels(["1", "0.1", "0.01", "0.001", "0.0001"])
    ax.minorticks_off()
    ax.set_xlabel("반경 (arcsec)"); ax.set_ylabel("정규화 세기")
    ax.legend(fontsize=6.5, frameon=False, loc="upper right")
    ax.set_title("반경 프로파일 — 같은 각크기로 겹쳐 보면 날개 두께가 갈린다",
                 fontsize=7.5)

    for column, (path, shape, colour, name) in enumerate((
            (args.cdk_field, tuple(args.cdk_shape), C["data"], "CDK 508 mm"),
            (args.rc_field, tuple(args.rc_shape), C["reference"], "LCO 1 m"))):
        table = pd.read_csv(path)
        median_e = float(np.nanmedian(table["ell"]))
        radial = float(np.nanmedian(table["radial_align"]))
        whisker_panel(fig.add_subplot(gs[2, column]), table, shape,
                      f"{name} · 타원율 중앙 {median_e:.3f} · 장축–방사각 {radial:.0f}°",
                      colour)

    fig.text(0.5, 0.055,
             "막대 방향 = 장축, 길이 ∝ 타원율. 코마라면 모두 시야중심을 향하고 "
             "가장자리에서 길어진다 — 둘 다 그렇지 않다(방사각 약 45° = 무작위 배향).",
             ha="center", fontsize=6.2, color=PALETTE["grey"])
    fig.text(0.5, 0.012,
             "실측 ePSF · Moravian C3-61000 @ CDK 508 mm (M13, B, 0.393″/px) · "
             "LCO 1m0-08 Sinistro (NGC 5985, rp, 0.390″/px) · "
             "회절한계 0.27″ 는 FWHM 2.85″ 의 1/10 — 시상 지배, 회절무늬 관측 불가",
             ha="center", fontsize=5.8, color=PALETTE["grey"])
    fig.tight_layout(rect=(0, 0.075, 1, 1))

    outdir = Path(args.outdir)
    if _STYLE:
        for kind, path in save_fig(fig, "fig_psf_shape_optics", outdir).items():
            print(f"[{kind}] {path}")
    else:
        path = outdir / "fig_psf_shape_optics.png"
        fig.savefig(path, dpi=170)
        print(f"[png] {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
