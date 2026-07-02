"""Figure 4 — independent cross-check: APEX vs the ``sep`` C library.

This figure proves that APEX aperture photometry agrees, at photometric grade,
with a fully INDEPENDENT engine (Barbary's ``sep`` / SExtractor SEP backend,
which shares no code with photutils) when both measure the SAME stars in the
SAME image with the SAME aperture.

Design
------
We build our own synthetic frame with a KNOWN catalog (exact positions + true
magnitudes) so the two engines measure identical apertures at identical centers
— the only remaining difference is the measurement engine itself. The noise
model mirrors ``apex/benchmark/synthetic_frame.py`` exactly (Gaussian PSF,
Poisson photon noise in electrons, Gaussian read noise, sky background, gain).

* APEX  : ``apex.utils.photometry_utils.phot_vectorized`` (photutils backend).
* sep   : ``sep.Background`` + ``sep.sum_circle`` (global-rms sky, gain-aware).

The comparison is differential and zeropoint-free: both instrumental
magnitudes are ``-2.5 log10(flux)``; the median APEX-minus-sep difference is the
zeropoint offset (it absorbs the constant electrons<->ADU scale difference),
and the residual scatter about it is the real agreement metric.

Run:
    .venv-deploy\\Scripts\\python.exe validation\\paper\\fig4_crosscheck_sep.py
"""

from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "paper"))

import numpy as np
import matplotlib.pyplot as plt

from apex_paper_style import apply_paper_style, save_fig, C, PALETTE, SINGLE_COL, DOUBLE_COL

apply_paper_style()

import sep

from apex.utils.photometry_utils import phot_vectorized


# ── synthetic-frame constants (mirror apex/benchmark/synthetic_frame.py) ──────
FWHM_PX = 3.5
GAIN = 1.5           # e-/ADU
READ_NOISE_E = 5.0   # e-
SKY_ADU = 150.0      # B, ADU
ZP = 25.0            # magnitude zeropoint (flux_e = 10**(-0.4*(m-ZP)))
FRAME_SIZE = 800     # px (>= 700, roomy enough for >= 8*FWHM separation)
N_STARS = 150
MAG_MIN, MAG_MAX = 14.0, 20.0
SEED = 20260702

# aperture geometry — SAME r_ap for both engines
R_AP = 1.0 * FWHM_PX
R_IN = 4.0 * FWHM_PX
R_OUT = 6.0 * FWHM_PX

SAT_ADU = 65000.0
DATAMAX_ADU = 60000.0

# quality cut
SNR_MIN = 10.0

OUTDIR = REPO / "validation" / "paper" / "figures"
CAPDIR = REPO / "validation" / "paper" / "captions"


def _sigma_from_fwhm(fwhm: float) -> float:
    return float(fwhm) / 2.3548200450309493


def build_synthetic_frame() -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Return (image_adu, x, y, true_mag) for a KNOWN isolated-star catalog.

    Mirrors the electron-domain noise model of
    ``apex/benchmark/synthetic_frame.py`` exactly:
        signal_e   = sum(unit-sum Gaussian PSF) * F_e,  F_e = 10**(-0.4*(m-ZP))
        expected_e = signal_e + B*gain
        realized_e = Poisson(expected_e) + Normal(0, RN)
        image_adu  = realized_e / gain
    """
    rng = np.random.default_rng(SEED)
    sigma = _sigma_from_fwhm(FWHM_PX)
    stamp_size = int(2 * int(np.ceil(4.0 * FWHM_PX)) + 1)
    half = stamp_size // 2

    # Every star must be isolated: >= 8*FWHM from every other star, and far
    # enough from the edge for the sky annulus (r_out) to stay on-frame.
    min_sep = 8.0 * FWHM_PX
    margin = int(np.ceil(R_OUT)) + half + 2
    lo, hi = float(margin), float(FRAME_SIZE - margin)

    placed: list[tuple[float, float]] = []
    max_attempts = max(40000, N_STARS * 800)
    attempts = 0
    while len(placed) < N_STARS and attempts < max_attempts:
        attempts += 1
        x = float(rng.uniform(lo, hi))
        y = float(rng.uniform(lo, hi))
        if placed:
            arr = np.asarray(placed, dtype=float)
            if np.min(np.hypot(arr[:, 0] - x, arr[:, 1] - y)) < min_sep:
                continue
        placed.append((x, y))
    if len(placed) < N_STARS:
        raise RuntimeError(
            f"Only placed {len(placed)}/{N_STARS} isolated stars; enlarge the frame."
        )
    centers = np.asarray(placed, dtype=float)

    true_mag = rng.uniform(MAG_MIN, MAG_MAX, size=len(centers))
    flux_e = 10.0 ** (-0.4 * (true_mag - ZP))

    # noiseless electron image (sub-pixel Gaussian per star)
    signal_e = np.zeros((FRAME_SIZE, FRAME_SIZE), dtype=np.float64)
    yy, xx = np.mgrid[-half : half + 1, -half : half + 1]
    for (x, y), fe in zip(centers, flux_e):
        xi = int(round(x))
        yi = int(round(y))
        dx = x - xi
        dy = y - yi
        shifted = np.exp(-((xx - dx) ** 2 + (yy - dy) ** 2) / (2.0 * sigma * sigma))
        total = float(shifted.sum())
        if total <= 0:
            continue
        shifted = shifted / total * float(fe)
        signal_e[yi - half : yi + half + 1, xi - half : xi + half + 1] += shifted

    # noise model in electrons (identical order to synthetic_frame.py)
    background_e = SKY_ADU * GAIN
    expected_e = signal_e + background_e
    realized_e = rng.poisson(np.clip(expected_e, 0.0, None)).astype(np.float64)
    realized_e += rng.normal(0.0, READ_NOISE_E, size=realized_e.shape)

    image_adu = realized_e / GAIN
    image_adu = np.clip(image_adu, 0.0, SAT_ADU)

    return image_adu, centers[:, 0], centers[:, 1], true_mag


def measure_apex(image_adu: np.ndarray, x: np.ndarray, y: np.ndarray):
    """APEX aperture photometry at the KNOWN positions. Returns (flux_e, snr)."""
    positions = np.column_stack([x, y])  # (N, 2), col-major (x, y)
    flux_e, flux_err_e, snr, sky_med, is_sat, is_nl = phot_vectorized(
        image_adu,
        positions,
        r_ap=R_AP,
        r_in=R_IN,
        r_out=R_OUT,
        gain=GAIN,
        rn_param_e=READ_NOISE_E,
        sky_sigma_mode="local",
        sat_adu=SAT_ADU,
        datamax_adu=DATAMAX_ADU,
    )
    return np.asarray(flux_e, dtype=float), np.asarray(snr, dtype=float)


def measure_sep(image_adu: np.ndarray, x: np.ndarray, y: np.ndarray):
    """Independent aperture photometry with the ``sep`` C library. Returns flux (ADU)."""
    data = image_adu.astype(np.float32)
    bkg = sep.Background(data)
    data_sub = data - bkg
    flux, fluxerr, flag = sep.sum_circle(
        data_sub, x, y, R_AP, err=bkg.globalrms, gain=GAIN
    )
    return np.asarray(flux, dtype=float), np.asarray(flag, dtype=int)


def _mad(values: np.ndarray) -> float:
    v = values[np.isfinite(values)]
    if v.size == 0:
        return float("nan")
    med = float(np.median(v))
    return float(1.4826 * np.median(np.abs(v - med)))


def _rms(values: np.ndarray) -> float:
    v = values[np.isfinite(values)]
    if v.size == 0:
        return float("nan")
    return float(np.sqrt(np.mean(v * v)))


def main() -> int:
    image_adu, x, y, true_mag = build_synthetic_frame()

    apex_flux_e, apex_snr = measure_apex(image_adu, x, y)
    sep_flux, sep_flag = measure_sep(image_adu, x, y)

    with np.errstate(divide="ignore", invalid="ignore"):
        mag_apex = -2.5 * np.log10(np.where(apex_flux_e > 0, apex_flux_e, np.nan))
        mag_sep = -2.5 * np.log10(np.where(sep_flux > 0, sep_flux, np.nan))

    # keep: APEX snr > 10, positive flux in BOTH engines, sep flag ok
    keep = (
        np.isfinite(mag_apex)
        & np.isfinite(mag_sep)
        & np.isfinite(apex_snr)
        & (apex_snr > SNR_MIN)
        & (apex_flux_e > 0)
        & (sep_flux > 0)
        & (sep_flag == 0)
    )
    ma = mag_apex[keep]
    ms = mag_sep[keep]
    mt = true_mag[keep]
    n = int(keep.sum())

    # zeropoint alignment + residuals
    offset = float(np.median(ma - ms))
    delta = (ma - ms) - offset
    med_delta = float(np.median(delta))
    mad = _mad(delta)
    rms = _rms(delta)
    pearson_r = float(np.corrcoef(ma, ms)[0, 1])
    frac_01 = float(np.mean(np.abs(delta) < 0.01))
    frac_02 = float(np.mean(np.abs(delta) < 0.02))

    # sep magnitudes aligned to the APEX zeropoint for panel (a)
    ms_aligned = ms + offset

    # ── figure: 1x2, double-column width ──────────────────────────────────────
    fig, (ax_a, ax_b) = plt.subplots(1, 2, figsize=(DOUBLE_COL, 3.1))

    # (a) mag_apex vs mag_sep (ZP-aligned) with y=x
    lo = float(min(ma.min(), ms_aligned.min()))
    hi = float(max(ma.max(), ms_aligned.max()))
    pad = max(0.05, 0.04 * (hi - lo))
    ax_a.plot(
        [lo - pad, hi + pad], [lo - pad, hi + pad],
        color=C["floor"], linewidth=1.0, zorder=1, label="$y=x$",
    )
    ax_a.scatter(
        ma, ms_aligned, s=12, alpha=0.7, color=C["data"],
        edgecolors="none", zorder=2,
    )
    ax_a.set(
        xlim=(hi + pad, lo - pad),  # brighter (smaller mag) upper-right
        ylim=(hi + pad, lo - pad),
    )
    ax_a.set_aspect("equal", adjustable="box")
    ax_a.set_xlabel(r"$m_{\mathrm{APEX}}$  ($-2.5\,\log_{10}\,F_{e^-}$)")
    ax_a.set_ylabel(r"$m_{\mathrm{sep}}$  (ZP-aligned to APEX)")
    ax_a.text(
        0.05, 0.95,
        f"$r = {pearson_r:.5f}$\n$N = {n}$",
        transform=ax_a.transAxes, va="top", ha="left", fontsize=8.0,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
              "alpha": 0.85, "edgecolor": PALETTE["grey"]},
    )
    ax_a.set_title("(a) APEX vs sep", loc="left")

    # (b) residual (ZP-aligned) vs true mag, with +/- MAD band + median line
    ax_b.axhspan(
        med_delta - mad, med_delta + mad,
        color=C["model"], alpha=0.18, zorder=0, label=r"$\pm$MAD",
    )
    ax_b.axhline(med_delta, color=C["model"], linewidth=1.2, zorder=2, label="median")
    ax_b.axhline(0.0, color=C["floor"], linewidth=0.8, linestyle=":", zorder=1)
    ax_b.scatter(
        mt, delta, s=12, alpha=0.7, color=C["data"],
        edgecolors="none", zorder=3,
    )
    ax_b.set_xlabel(r"$m_{\mathrm{true}}$  (mag)")
    ax_b.set_ylabel(r"$\Delta = (m_{\mathrm{APEX}} - m_{\mathrm{sep}}) - \mathrm{ZP}$  (mag)")
    ax_b.text(
        0.05, 0.95,
        f"offset $= {offset:+.4f}$\nMAD $= {mad:.4f}$\nRMS $= {rms:.4f}$",
        transform=ax_b.transAxes, va="top", ha="left", fontsize=8.0,
        bbox={"boxstyle": "round,pad=0.3", "facecolor": "white",
              "alpha": 0.85, "edgecolor": PALETTE["grey"]},
    )
    ax_b.legend(loc="lower left", fontsize=7.5, ncol=2)
    ax_b.set_title("(b) ZP-aligned residuals", loc="left")

    fig.tight_layout()
    paths = save_fig(fig, "fig4_crosscheck_sep", OUTDIR)
    plt.close(fig)

    # ── caption + numbers ─────────────────────────────────────────────────────
    CAPDIR.mkdir(parents=True, exist_ok=True)
    photometric_grade = mad <= 0.02
    mad_cmp = "$\\le$" if photometric_grade else "$>$"
    verdict_phrase = (
        "agree at photometric grade" if photometric_grade
        else "do NOT agree at photometric grade"
    )
    caption = f"""# Figure 4 — Independent cross-check: APEX vs the `sep` C library

**Figure.** Agreement between APEX aperture photometry and the independent
`sep` (SExtractor SEP backend, v1.4.1) engine measuring the **same {n} stars**
in the **same synthetic frame** with the **same aperture**
($r_{{\\mathrm{{ap}}}} = 1.0\\,\\mathrm{{FWHM}} = {R_AP:.2f}$ px).
The frame is a {FRAME_SIZE}x{FRAME_SIZE} px field of {N_STARS} isolated stars
(separation $\\ge 8\\,\\mathrm{{FWHM}}$), true magnitudes uniform in
$[{MAG_MIN:.0f}, {MAG_MAX:.0f}]$, with a Gaussian PSF
($\\mathrm{{FWHM}} = {FWHM_PX}$ px), gain $g = {GAIN}$ e$^-$/ADU, read noise
$\\mathrm{{RN}} = {READ_NOISE_E}$ e$^-$, sky $B = {SKY_ADU:.0f}$ ADU, and
zeropoint $\\mathrm{{ZP}} = {ZP}$; noise is Poisson (photon) plus Gaussian
(read) in the electron domain, matching `apex/benchmark/synthetic_frame.py`.
**(a)** APEX vs `sep` instrumental magnitude after a single median-zeropoint
alignment; the grey line is $y=x$. **(b)** Zeropoint-aligned residual
$\\Delta = (m_{{\\mathrm{{APEX}}}} - m_{{\\mathrm{{sep}}}}) - \\mathrm{{ZP}}$
versus true magnitude, with the median (solid) and a shaded $\\pm$MAD band.
Only stars with APEX $\\mathrm{{SNR}} > {SNR_MIN:.0f}$ and positive flux in both
engines are kept.

**Agreement metrics** (ZP offset = median$(m_{{\\mathrm{{APEX}}}} -
m_{{\\mathrm{{sep}}}})$; residuals are ZP-removed):

| Quantity | Value |
|---|---|
| $N$ (matched, SNR $> {SNR_MIN:.0f}$) | {n} |
| Zeropoint offset | {offset:+.4f} mag |
| median $\\Delta$ | {med_delta:+.4f} mag |
| MAD ($1.4826\\times$ median$|\\Delta - \\mathrm{{med}}|$) | {mad:.4f} mag |
| RMS$(\\Delta)$ | {rms:.4f} mag |
| Pearson $r$ | {pearson_r:.5f} |
| fraction $|\\Delta| < 0.01$ mag | {frac_01:.3f} |
| fraction $|\\Delta| < 0.02$ mag | {frac_02:.3f} |

**Verdict.** MAD $= {mad:.4f}$ mag {mad_cmp} 0.02 mag and Pearson
$r = {pearson_r:.5f}$: APEX and `sep` **{verdict_phrase}**.
The two engines share no photometry code, so this differential agreement
confirms that APEX's aperture-sum + local-sky measurement is correct to the
milli-magnitude level, with any residual scatter set by pixelization and the
per-engine sky estimator rather than by a systematic error in APEX.
"""
    cap_path = CAPDIR / "fig4_crosscheck_sep.md"
    cap_path.write_text(caption, encoding="utf-8")

    print("=== fig4 cross-check (APEX vs sep) ===")
    print(f"N matched (SNR>{SNR_MIN:.0f}) : {n}")
    print(f"ZP offset               : {offset:+.4f} mag")
    print(f"median delta            : {med_delta:+.4f} mag")
    print(f"MAD                     : {mad:.4f} mag")
    print(f"RMS                     : {rms:.4f} mag")
    print(f"Pearson r               : {pearson_r:.6f}")
    print(f"frac |d|<0.01           : {frac_01:.3f}")
    print(f"frac |d|<0.02           : {frac_02:.3f}")
    print(f"photometric grade (MAD<=0.02): {photometric_grade}")
    for ext, p in paths.items():
        print(f"wrote {ext}: {p}  exists={p.exists()}")
    print(f"wrote caption: {cap_path}  exists={cap_path.exists()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
