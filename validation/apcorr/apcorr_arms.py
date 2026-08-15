"""Two more arms on the same stars: an independent kernel, and independent choices.

The control arm (`reproduce_apcorr.py`) showed the reimplementation lands on
APEX's number to the recorded precision, which means the orchestration is
understood. It could not say whether the number is *right*, because the
aperture summation was photutils on both sides.

This adds the two arms that can:

**Kernel** — swap photutils' `aperture_photometry` for `sep.sum_circle`
(Bertin's SExtractor backend, a genuinely separate C implementation) and change
nothing else. A difference here is the pixel-summation itself: sub-pixel
weighting, edge handling, how a partially covered pixel is counted.

**Choices** — keep the kernel and vary one orchestration decision at a time.
APEX makes six that are not forced by physics, and the question each variant
asks is "how much does apcorr move if a reasonable person chose differently".
A choice that moves it by less than the measurement noise is not worth
defending; one that moves it a lot is a documented sensitivity, not a bug.

Variables held fixed across every arm (V6/V9): the same frames, the same
reference stars reconstructed from Step 7's own record, the same per-star sky,
the same radius grid — unless the variant under test is precisely that.

    python validation/apcorr/apcorr_arms.py --workspace <phase3/M67>
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
sys.path.insert(0, str(Path(__file__).absolute().parent))

from reproduce_apcorr import (  # noqa: E402
    GC_N_STEPS, _num, apcorr_from_curve, growth_curve, reference_mask,
)


def growth_curve_sep(image: np.ndarray, positions: np.ndarray,
                     radii: np.ndarray, sky_median: np.ndarray) -> np.ndarray:
    """The same growth curve with SExtractor's summation instead of photutils.

    `sep` wants a native-byte-order contiguous array and sums exactly the
    aperture area, so the sky is removed the same way — the only thing that
    differs is which C code adds the pixels up.
    """
    import sep

    data = np.ascontiguousarray(image, dtype=np.float64)
    if data.dtype.byteorder not in ("=", "|"):
        data = data.byteswap().view(data.dtype.newbyteorder("="))
    x, y = positions[:, 0], positions[:, 1]
    flux = np.full((len(radii), len(positions)), np.nan, dtype=float)
    for index, radius in enumerate(radii):
        total, _err, _flag = sep.sum_circle(data, x, y, float(radius),
                                            subpix=0)  # exact, matches photutils
        flux[index] = np.asarray(total, float) - sky_median * np.pi * radius ** 2
    return flux


def _curve_from_flux(flux: np.ndarray, min_n: int = 12):
    outer = flux[-1, :]
    good = np.isfinite(outer) & (outer > 0)
    if int(good.sum()) < min_n:
        return None
    return np.nanmedian(flux[:, good] / outer[good], axis=1), good


def variant_apcorr(flux: np.ndarray, radii: np.ndarray, r_ap: float,
                   *, stack: str = "median", lookup: str = "nearest") -> float:
    """apcorr under one changed choice, everything else as APEX has it.

    `stack`  — median (APEX) vs mean: does one odd star matter?
    `lookup` — nearest grid point (APEX) vs linear interpolation: is the
               14-point grid coarse enough to move the answer?
    """
    result = _curve_from_flux(flux)
    if result is None:
        return float("nan")
    per_star = flux[:, result[1]] / flux[-1, result[1]]
    curve = (np.nanmedian(per_star, axis=1) if stack == "median"
             else np.nanmean(per_star, axis=1))
    if lookup == "nearest":
        enclosed = float(curve[int(np.argmin(np.abs(radii - r_ap)))])
    else:
        enclosed = float(np.interp(r_ap, radii, curve))
    return float(1.0 / enclosed) if enclosed > 0.05 else float("nan")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workspace", required=True, type=Path)
    ap.add_argument("--frames", type=int, default=0)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    from apex.config.parameters_cmd import read_params
    from apex.analysis.forced_photometry import _to_float

    params = read_params(args.workspace / "apex_config.json")
    P = params.P
    step7 = Path(P.result_dir) / "step7_forced_phot"
    data_dir = Path(P.data_dir)

    recorded = pd.read_csv(step7 / "apcorr_summary.csv")
    stats = pd.read_csv(step7 / "frame_stats.csv")
    master_path = step7 / "master_sources.csv"
    master = pd.read_csv(master_path) if master_path.exists() else None

    r_ap_scale = _to_float(getattr(P, "forced_r_ap_scale", 0.8), 0.8)
    ref_scale = _to_float(getattr(P, "forced_ref_ap_scale", 2.4), 2.4)
    min_r_ap = _to_float(getattr(P, "min_r_ap_px", 4.0), 4.0)
    fwhm_column = next(c for c in ("fwhm_px", "fwhm", "fwhm_median_px")
                       if c in stats.columns)

    rows = []
    frames = recorded if args.frames <= 0 else recorded.head(args.frames)
    for entry in frames.itertuples(index=False):
        name = str(entry.file)
        table, image_path = step7 / f"photometry_{name}.tsv", data_dir / name
        stat_row = stats[stats["file"].astype(str) == name]
        if not table.exists() or not image_path.exists() or stat_row.empty:
            continue

        fwhm = float(pd.to_numeric(stat_row.iloc[0][fwhm_column], errors="coerce"))
        r_ap = max(min_r_ap, r_ap_scale * fwhm)
        r_ref = max(r_ap + 2.0, ref_scale * fwhm)
        radii = np.linspace(max(2.0, r_ap * 0.4), r_ref * 1.15, GC_N_STEPS)

        phot = pd.read_csv(table, sep="\t")
        chosen = phot[reference_mask(phot, master)]
        positions = np.column_stack([_num(chosen, "x_fit"), _num(chosen, "y_fit")])
        finite = np.isfinite(positions).all(axis=1)
        positions, chosen = positions[finite], chosen[finite]
        sky = _num(chosen, "sky")

        with fits.open(image_path, memmap=False) as hdul:
            image = np.asarray(hdul[0].data, dtype=float)

        photutils_flux = growth_curve(image, positions, radii, sky)
        sep_flux = growth_curve_sep(image, positions, radii, sky)

        base, _curve, n_used = apcorr_from_curve(photutils_flux, radii, r_ap)

        # A denser grid over the same span: the 14 points are APEX's choice,
        # and the nearest-point lookup makes that choice part of the answer.
        dense = np.linspace(radii[0], radii[-1], 141)
        dense_flux = growth_curve(image, positions, dense, sky)

        rows.append({
            "frame": name, "filter": entry.filter,
            "fwhm_px": fwhm, "r_ap_px": r_ap, "n_stars": int(len(positions)),
            "apex": float(entry.apcorr),
            "photutils": base,
            "sep": apcorr_from_curve(sep_flux, radii, r_ap)[0],
            "mean_stack": variant_apcorr(photutils_flux, radii, r_ap, stack="mean"),
            "interpolated": variant_apcorr(photutils_flux, radii, r_ap, lookup="linear"),
            "dense_grid": variant_apcorr(dense_flux, dense, r_ap),
        })
        print(f"  {name}", flush=True)

    if not rows:
        print("[error] 처리한 프레임이 없다")
        return 1

    frame = pd.DataFrame(rows)
    arms = [("sep 커널", "sep"), ("평균 쌓기", "mean_stack"),
            ("선형 보간", "interpolated"), ("촘촘한 격자 141", "dense_grid")]
    print(f"\n기준: photutils 재구현 (APEX 재현 확인됨) · {len(frame)} 프레임\n")
    print(f'{"팔":<18}{"중앙 Δapcorr":>14}{"최대 |Δ|":>12}{"등급 환산 중앙":>16}{"최대":>10}')
    print("-" * 72)
    for label, column in arms:
        delta = (frame[column] - frame["photutils"]).to_numpy(float)
        delta = delta[np.isfinite(delta)]
        if delta.size == 0:
            print(f"{label:<18}(값 없음)"); continue
        # apcorr is a flux multiplier; what matters downstream is the magnitude.
        mag = 2.5 * np.log10(frame[column] / frame["photutils"])
        mag = mag.to_numpy(float)[np.isfinite(mag)]
        print(f"{label:<18}{np.median(delta):>+14.2e}{np.max(np.abs(delta)):>12.2e}"
              f"{np.median(mag) * 1000:>+15.2f}m{np.max(np.abs(mag)) * 1000:>9.1f}m")

    if args.output:
        frame.to_csv(args.output, index=False)
        args.output.with_suffix(".inputs.json").write_text(json.dumps({
            "workspace": str(args.workspace), "frames": len(frame),
            "held_fixed": ["frames", "reference stars (from Step 7 record)",
                           "per-star sky", "radius grid (except dense_grid arm)"],
            "arms": {"sep": "SExtractor kernel, subpix=0 (exact)",
                     "mean_stack": "mean instead of median across stars",
                     "interpolated": "linear interpolation instead of nearest grid point",
                     "dense_grid": "141 radii over the same span"},
            "gc_n_steps": GC_N_STEPS,
        }, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"\n표 -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
