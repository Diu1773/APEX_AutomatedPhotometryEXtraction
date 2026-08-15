"""Which sky estimator is right — asked where the answer is known.

The IRAF cross-check ended with the two engines' aperture corrections differing
by exactly the amount their sky estimators differ (M13 −14.3 mmag, M67 −2.0,
both predicted to within 1 mmag from the measured sky offset). That settles
*what* the difference is and leaves the question a cross-check cannot answer:
IRAF's `mode` or APEX's sigma-clipped median — which one is closer to the sky
that is actually there?

Injection answers it. Plant stars of known flux into a real frame, measure the
annulus with each estimator, and see whose recovered magnitude lands on truth.
The frame supplies the thing that makes this hard and that a synthetic sky
would not: the annulus of a real star in a real cluster contains other stars,
and how an estimator handles that contamination is the whole question.

Split by crowding, because that is where the two diverge. In an empty field a
symmetric sky distribution makes mode and median agree by construction; the
gap opens when faint neighbours put a positive tail on it.

    python validation/apcorr/sky_estimator_truth.py --workspace <phase3/M13>
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))

MAD_TO_SIGMA = 1.4826


def sky_sigma_clipped_median(values: np.ndarray, sigma: float = 3.0,
                             maxiters: int = 5) -> float:
    """APEX's estimator: clip outliers, take the median of what is left."""
    from astropy.stats import sigma_clipped_stats

    if values.size == 0:
        return float("nan")
    _mean, median, _std = sigma_clipped_stats(values, sigma=sigma, maxiters=maxiters)
    return float(median)


def sky_mode(values: np.ndarray) -> float:
    """IRAF's `mode`: 3*median - 2*mean, on sigma-clipped pixels.

    DAOPHOT's `fitskypars.salgorithm = "mode"` uses the classic empirical
    relation for a mildly skewed unimodal distribution. It is not the histogram
    peak; it is an estimate of it that costs one pass. Contamination pushes the
    mean further than the median, so this reads *below* the median exactly when
    neighbours are present — which is the sign seen against APEX.
    """
    from astropy.stats import sigma_clipped_stats

    if values.size == 0:
        return float("nan")
    mean, median, _std = sigma_clipped_stats(values, sigma=3.0, maxiters=5)
    return float(3.0 * median - 2.0 * mean)


def sky_plain_median(values: np.ndarray) -> float:
    """No clipping at all — the control that shows what clipping buys."""
    return float(np.median(values)) if values.size else float("nan")


ESTIMATORS = {
    "APEX (시그마클립 중앙값)": sky_sigma_clipped_median,
    "IRAF (mode)": sky_mode,
    "무보정 중앙값": sky_plain_median,
}


def annulus_pixels(image: np.ndarray, x: float, y: float,
                   r_in: float, r_out: float) -> np.ndarray:
    """Pixel values in the sky annulus, as the estimators see them."""
    lo_x, hi_x = int(np.floor(x - r_out)), int(np.ceil(x + r_out)) + 1
    lo_y, hi_y = int(np.floor(y - r_out)), int(np.ceil(y + r_out)) + 1
    lo_x, lo_y = max(lo_x, 0), max(lo_y, 0)
    hi_x, hi_y = min(hi_x, image.shape[1]), min(hi_y, image.shape[0])
    if hi_x <= lo_x or hi_y <= lo_y:
        return np.empty(0)
    patch = image[lo_y:hi_y, lo_x:hi_x]
    yy, xx = np.mgrid[lo_y:hi_y, lo_x:hi_x]
    radius = np.hypot(xx - x, yy - y)
    inside = (radius >= r_in) & (radius < r_out) & np.isfinite(patch)
    return patch[inside]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workspace", required=True, type=Path)
    ap.add_argument("--frame", default=None,
                    help="프레임 이름 (기본: apcorr_summary.csv 첫 줄)")
    ap.add_argument("--n-inject", type=int, default=400)
    ap.add_argument("--snr", type=float, default=60.0,
                    help="주입 별의 목표 SNR — 하늘 오차가 등급에 미치는 크기를 "
                         "밝기와 분리하려고 한 값으로 고정한다")
    ap.add_argument("--seed", type=int, default=20260816)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    from apex.config.parameters_cmd import read_params
    from apex.analysis.forced_photometry import _to_float
    from photutils.aperture import CircularAperture, aperture_photometry

    params = read_params(args.workspace / "apex_config.json")
    P = params.P
    step7 = Path(P.result_dir) / "step7_forced_phot"
    recorded = pd.read_csv(step7 / "apcorr_summary.csv")
    frame = args.frame or str(recorded.iloc[0]["file"])
    stats = pd.read_csv(step7 / "frame_stats.csv")
    row = stats[stats["file"].astype(str) == frame]
    if row.empty:
        print(f"[error] frame_stats 에 {frame} 없음")
        return 1

    fwhm = float(pd.to_numeric(row.iloc[0]["fwhm_px"], errors="coerce"))
    r_ap = max(_to_float(getattr(P, "min_r_ap_px", 4.0), 4.0),
               _to_float(getattr(P, "forced_r_ap_scale", 0.8), 0.8) * fwhm)
    r_ref = max(r_ap + 2.0, _to_float(getattr(P, "forced_ref_ap_scale", 2.4), 2.4) * fwhm)
    gap = _to_float(getattr(P, "annulus_min_gap_px", 6.0), 6.0)
    r_in = max(r_ref + gap, _to_float(getattr(P, "fitsky_annulus_scale", 4.0), 4.0) * fwhm)
    r_out = r_in + max(gap, _to_float(getattr(P, "fitsky_dannulus_scale", 2.0), 2.0) * fwhm)

    with fits.open(Path(P.data_dir) / frame, memmap=False) as hdul:
        image = np.asarray(hdul[0].data, dtype=float)

    phot = pd.read_csv(step7 / f"photometry_{frame}.tsv", sep="\t")
    sky_level = float(pd.to_numeric(phot["sky"], errors="coerce").median())
    sky_noise = float(pd.to_numeric(phot["sky_std"], errors="coerce").median())
    gain = float(pd.to_numeric(phot["gain_e_per_adu"], errors="coerce").median())

    # Existing sources, so an injected star can be placed at a known distance
    # from the nearest one — crowding is the variable under test, not a nuisance.
    real = np.column_stack([
        pd.to_numeric(phot["x_fit"], errors="coerce").to_numpy(float),
        pd.to_numeric(phot["y_fit"], errors="coerce").to_numpy(float)])
    real = real[np.isfinite(real).all(axis=1)]
    from scipy.spatial import cKDTree
    tree = cKDTree(real)

    rng = np.random.default_rng(args.seed)
    height, width = image.shape
    margin = int(np.ceil(r_out)) + 6
    xs = rng.uniform(margin, width - margin, args.n_inject)
    ys = rng.uniform(margin, height - margin, args.n_inject)
    separation = tree.query(np.column_stack([xs, ys]), k=1)[0]

    # Flux for the requested SNR, using the frame's own noise. A Gaussian of
    # this FWHM has a noise-equivalent area of 4*pi*sigma^2.
    sigma_px = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    nea = 4.0 * np.pi * sigma_px ** 2
    total_flux = args.snr * sky_noise * np.sqrt(nea)

    truth = image.copy()
    yy, xx = np.mgrid[0:height, 0:width]
    for x, y in zip(xs, ys):
        lo_x, hi_x = int(x - 5 * sigma_px), int(x + 5 * sigma_px) + 1
        lo_y, hi_y = int(y - 5 * sigma_px), int(y + 5 * sigma_px) + 1
        lo_x, lo_y = max(lo_x, 0), max(lo_y, 0)
        hi_x, hi_y = min(hi_x, width), min(hi_y, height)
        gx = xx[lo_y:hi_y, lo_x:hi_x] - x
        gy = yy[lo_y:hi_y, lo_x:hi_x] - y
        profile = np.exp(-(gx ** 2 + gy ** 2) / (2.0 * sigma_px ** 2))
        truth[lo_y:hi_y, lo_x:hi_x] += total_flux * profile / (2.0 * np.pi * sigma_px ** 2)

    # The quantity under test is the aperture *correction*, not one aperture.
    # A sky error moves the small and the large aperture by delta*pi*r^2, so it
    # is the ratio between them that is sensitive to it — the r_out aperture
    # here has ~12x the area of r_ap. Measuring a single aperture would report
    # a few mmag and hide the effect that produced the 14 mmag engine gap.
    r_outer = r_ref * 1.15                       # APEX's outermost curve radius
    small = CircularAperture(np.column_stack([xs, ys]), r=r_ap)
    large = CircularAperture(np.column_stack([xs, ys]), r=r_outer)
    raw_small = aperture_photometry(truth, small)["aperture_sum"].data.astype(float)
    raw_large = aperture_photometry(truth, large)["aperture_sum"].data.astype(float)

    # A Gaussian's enclosed fraction is known, so the correct apcorr is too.
    frac_small = 1.0 - np.exp(-(r_ap ** 2) / (2.0 * sigma_px ** 2))
    frac_large = 1.0 - np.exp(-(r_outer ** 2) / (2.0 * sigma_px ** 2))
    truth_apcorr = frac_large / frac_small

    rows = []
    for index, (x, y) in enumerate(zip(xs, ys)):
        pixels = annulus_pixels(truth, x, y, r_in, r_out)
        entry = {"x": x, "y": y, "sep_px": separation[index],
                 "sep_fwhm": separation[index] / fwhm,
                 "truth_apcorr": truth_apcorr}
        for label, estimator in ESTIMATORS.items():
            sky = estimator(pixels)
            net_small = raw_small[index] - sky * small.area
            net_large = raw_large[index] - sky * large.area
            measured = net_large / net_small if net_small > 0 and net_large > 0 else np.nan
            # Positive = this estimator's apcorr is too large.
            entry[label] = (2.5 * np.log10(measured / truth_apcorr)
                            if np.isfinite(measured) and measured > 0 else np.nan)
            entry[f"sky::{label}"] = sky
        rows.append(entry)

    frame_table = pd.DataFrame(rows)
    bins = [(0.0, 2.0, "0-2 FWHM"), (2.0, 4.0, "2-4"), (4.0, 8.0, "4-8"),
            (8.0, np.inf, "8+")]
    print(f"{frame} · FWHM {fwhm:.2f} px · 하늘 {sky_level:.1f}±{sky_noise:.2f} ADU · "
          f"주입 {args.n_inject} · 목표 SNR {args.snr:.0f}")
    print(f"고리 {r_in:.1f}~{r_out:.1f} px · r_ap {r_ap:.2f}\n")
    print(f'{"최근접 이웃":<12}{"n":>5}' + "".join(f"{name:>26}" for name in ESTIMATORS))
    print("-" * (17 + 26 * len(ESTIMATORS)))
    summary = []
    for lo, hi, label in bins:
        mask = (frame_table["sep_fwhm"] >= lo) & (frame_table["sep_fwhm"] < hi)
        if mask.sum() < 5:
            continue
        cells = []
        for name in ESTIMATORS:
            values = frame_table.loc[mask, name].to_numpy(float)
            values = values[np.isfinite(values)]
            median = np.median(values) * 1000 if values.size else np.nan
            scatter = MAD_TO_SIGMA * np.median(np.abs(values - np.median(values))) * 1000
            cells.append(f"{median:+8.1f} ± {scatter:5.1f} mmag")
            summary.append({"bin": label, "estimator": name, "n": int(mask.sum()),
                            "bias_mmag": median, "scatter_mmag": scatter})
        print(f'{label:<12}{int(mask.sum()):>5}' + "".join(f"{c:>26}" for c in cells))

    print("\n편향이 0 에 가까운 쪽이 진리값에 가깝다 (양수 = 어둡게 잼).")
    if args.output:
        pd.DataFrame(summary).to_csv(args.output, index=False)
        args.output.with_suffix(".inputs.json").write_text(json.dumps({
            "workspace": str(args.workspace), "frame": frame, "fwhm_px": fwhm,
            "n_inject": args.n_inject, "target_snr": args.snr, "seed": args.seed,
            "r_ap_px": r_ap, "r_in_px": r_in, "r_out_px": r_out,
            "sky_adu": sky_level, "sky_noise_adu": sky_noise, "gain": gain,
            "injected_profile": "Gaussian at the frame's own FWHM",
            "measured_quantity": "apcorr = F(r_outer)/F(r_ap), the aperture correction",
            "truth": "ratio of Gaussian enclosed fractions, known analytically",
            "r_outer_px": float(r_ref * 1.15),
        }, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"표 -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
