"""Rebuild APEX's aperture correction from photutils primitives — the control arm.

This is step one of the aperture-correction cross-check (`①′` in
`docs/audit/APEX_CROSSCHECK_PLAN.md`), and it deliberately compares APEX
against *itself*. Before asking "is APEX's aperture correction right", the
question has to be "do I know what it computes" — a reimplementation that
lands on a different number for an unknown reason cannot then be used to
judge anything.

What is shared and what is not (checklist item V8, and the honest limit of
this arm): the aperture *summation* is photutils in both cases — APEX's
`_growth_curve_fixed_sky` imports `CircularAperture`/`aperture_photometry`
itself. So this arm tests the orchestration only:

    star selection · radius grid · fixed-sky handling · per-star
    normalisation · median stacking · which radius apcorr is read at

An independent *primitive* needs `sep` (Bertin's SExtractor backend), and an
independent *orchestration* needs different choices at each of those six
points. Both are separate arms, not this one.

The reference-star mask is reconstructed from what Step 7 recorded rather than
guessed, so the two arms see the same stars (V6). Every condition in
`forced_photometry.py`'s `apcorr_mask` maps to a persisted column; the one
that lives in the master catalogue (`crowding_flag`) is joined back by ID.

    python validation/apcorr/reproduce_apcorr.py --workspace <phase3/M67>
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

REPO = Path(__file__).absolute().parents[2]

# Mirrors apex/analysis/forced_photometry.py. Kept as literals on purpose: if
# the engine changes one of these, this arm must fail loudly rather than
# silently follow along.
GC_N_STEPS = 14
APCORR_MIN_SNR = 40.0
APCORR_MIN_N = 12
APCORR_MAX_SOURCES = 250


def _bool(frame: pd.DataFrame, column: str, default: bool = False) -> np.ndarray:
    if column not in frame.columns:
        return np.full(len(frame), default, dtype=bool)
    values = frame[column]
    if values.dtype == object:
        text = values.astype(str).str.strip().str.lower()
        return text.isin({"true", "1", "yes", "t"}).to_numpy()
    return values.fillna(default).astype(bool).to_numpy()


def _num(frame: pd.DataFrame, column: str) -> np.ndarray:
    if column not in frame.columns:
        return np.full(len(frame), np.nan)
    return pd.to_numeric(frame[column], errors="coerce").to_numpy(float)


def reference_mask(phot: pd.DataFrame, master: pd.DataFrame | None,
                   min_snr: float = APCORR_MIN_SNR,
                   max_sources: int = APCORR_MAX_SOURCES) -> np.ndarray:
    """The same stars the engine used, rebuilt from its own output.

    `crowding_flag` is the only condition Step 7 does not repeat into the
    per-frame table — it belongs to the master catalogue, so it is joined on
    the ID the two tables share.
    """
    crowding = np.zeros(len(phot), dtype=bool)
    if master is not None and "crowding_flag" in master.columns:
        key = next((c for c in ("master_id", "ID", "source_id")
                    if c in phot.columns and c in master.columns), None)
        if key is not None:
            lookup = dict(zip(master[key], _bool(master, "crowding_flag")))
            crowding = phot[key].map(lookup).fillna(False).astype(bool).to_numpy()

    snr = _num(phot, "snr")
    flux = _num(phot, "flux")
    mask = (
        _bool(phot, "detected_flag")
        & ~crowding
        & _bool(phot, "step4_quality_ok", True)
        & _bool(phot, "step4_apcorr_candidate", True)
        & ~_bool(phot, "centroid_outlier")
        & np.isfinite(flux) & (flux > 0)
        & np.isfinite(snr) & (snr >= min_snr)
        & ~_bool(phot, "is_saturated")
        & ~_bool(phot, "is_nonlinear")
    )
    if max_sources > 0 and mask.sum() > max_sources:
        index = np.flatnonzero(mask)
        order = np.argsort(np.nan_to_num(snr[index], nan=-np.inf))[::-1]
        capped = np.zeros(len(phot), dtype=bool)
        capped[index[order[:max_sources]]] = True
        mask = capped
    return mask


def growth_curve(image: np.ndarray, positions: np.ndarray, radii: np.ndarray,
                 sky_median: np.ndarray) -> np.ndarray:
    """Aperture sums at each radius, with the sky held fixed per star.

    Fixed sky is APEX's choice and is reproduced here rather than improved on:
    the annulus does not move while the aperture grows, so re-measuring it per
    radius would change the number for a reason that is not the one under test.
    """
    from photutils.aperture import CircularAperture, aperture_photometry

    flux = np.full((len(radii), len(positions)), np.nan, dtype=float)
    for index, radius in enumerate(radii):
        aperture = CircularAperture(positions, r=float(radius))
        raw = aperture_photometry(image, aperture)["aperture_sum"].data.astype(float)
        flux[index] = raw - sky_median * aperture.area
    return flux


def apcorr_from_curve(flux: np.ndarray, radii: np.ndarray, r_ap: float,
                      min_n: int = APCORR_MIN_N) -> tuple[float, np.ndarray, int]:
    """1 / enclosed fraction at the science radius — APEX's definition.

    Normalised per star by its own outermost flux, then the *median* across
    stars: an outlier star moves nothing. The radius is picked by nearest
    grid point, not by interpolation, so the grid itself is part of the answer.
    """
    outer = flux[-1, :]
    good = np.isfinite(outer) & (outer > 0)
    if int(good.sum()) < min_n:
        return float("nan"), np.full(len(radii), np.nan), int(good.sum())
    curve = np.nanmedian(flux[:, good] / outer[good], axis=1)
    # Interpolated at r_ap, mirroring the engine since 2026-08-16. It used to
    # read the nearest grid point, and any arm that kept doing so would carry
    # the same 0.1 mag bias — two engines sharing a bias look like agreement.
    finite = np.isfinite(curve)
    enclosed = (float(np.interp(r_ap, radii[finite], curve[finite]))
                if finite.sum() >= 2 else 0.0)
    value = float(1.0 / enclosed) if enclosed > 0.05 else 1.0
    return value, curve, int(good.sum())


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workspace", required=True, type=Path,
                    help="phase3 워크스페이스 (apex_config.json 이 있는 곳)")
    ap.add_argument("--frames", type=int, default=6,
                    help="확인할 프레임 수 (0 = 전부)")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    import sys
    sys.path.insert(0, str(REPO))
    from apex.config.parameters_cmd import read_params
    from apex.analysis.forced_photometry import _to_float

    params = read_params(args.workspace / "apex_config.json")
    P = params.P
    result_dir = Path(P.result_dir)
    data_dir = Path(P.data_dir)
    step7 = result_dir / "step7_forced_phot"

    recorded = pd.read_csv(step7 / "apcorr_summary.csv")
    stats = pd.read_csv(step7 / "frame_stats.csv")
    master_path = step7 / "master_sources.csv"
    master = pd.read_csv(master_path) if master_path.exists() else None

    r_ap_scale = _to_float(getattr(P, "forced_r_ap_scale", 0.8), 0.8)
    ref_scale = _to_float(getattr(P, "forced_ref_ap_scale", 2.4), 2.4)
    min_r_ap = _to_float(getattr(P, "min_r_ap_px", 4.0), 4.0)

    fwhm_column = next((c for c in ("fwhm_px", "fwhm", "fwhm_median_px")
                        if c in stats.columns), None)
    if fwhm_column is None:
        print(f"[error] frame_stats.csv 에 FWHM 열이 없다: {list(stats.columns)}")
        return 1

    rows = []
    frames = recorded if args.frames <= 0 else recorded.head(args.frames)
    print(f"{'프레임':<30}{'APEX':>10}{'재현':>10}{'차이':>12}{'선정':>8}{'사용':>8}")
    print("-" * 80)
    for entry in frames.itertuples(index=False):
        name = str(entry.file)
        table = step7 / f"photometry_{name}.tsv"
        image_path = data_dir / name
        if not table.exists() or not image_path.exists():
            print(f"{name:<30}(입력 없음)")
            continue

        stat_row = stats[stats["file"].astype(str) == name]
        if stat_row.empty:
            print(f"{name:<30}(frame_stats 없음)")
            continue
        fwhm = float(pd.to_numeric(stat_row.iloc[0][fwhm_column], errors="coerce"))
        r_ap = max(min_r_ap, r_ap_scale * fwhm)
        r_ref = max(r_ap + 2.0, ref_scale * fwhm)
        radii = np.linspace(max(2.0, r_ap * 0.4), r_ref * 1.15, GC_N_STEPS)

        phot = pd.read_csv(table, sep="\t")
        mask = reference_mask(phot, master)
        chosen = phot[mask]
        positions = np.column_stack([_num(chosen, "x_fit"), _num(chosen, "y_fit")])
        finite = np.isfinite(positions).all(axis=1)
        positions, chosen = positions[finite], chosen[finite]

        with fits.open(image_path, memmap=False) as hdul:
            image = np.asarray(hdul[0].data, dtype=float)

        flux = growth_curve(image, positions, radii, _num(chosen, "sky"))
        value, _curve, n_used = apcorr_from_curve(flux, radii, r_ap)

        # Two different counts, and they must not be confused: the engine
        # reports the *mask* size, while the curve is built from the stars
        # whose outermost flux is usable. Comparing the second against the
        # first read as a -1 mismatch on every 250-capped frame.
        n_selected = int(len(positions))
        difference = value - float(entry.apcorr)
        print(f"{name:<30}{entry.apcorr:>10.6f}{value:>10.6f}"
              f"{difference:>+12.2e}{n_selected:>8}{n_used:>8}")
        rows.append({
            "frame": name, "filter": entry.filter,
            "apex_apcorr": float(entry.apcorr), "reproduced_apcorr": value,
            "difference": difference,
            "apex_n_stars": int(entry.n_apcorr_stars),
            "reproduced_n_selected": n_selected, "reproduced_n_used": n_used,
            "fwhm_px": fwhm, "r_ap_px": r_ap, "r_ref_px": r_ref,
        })

    if rows:
        frame = pd.DataFrame(rows)
        worst = frame["difference"].abs().max()
        star_match = int((frame["apex_n_stars"] == frame["reproduced_n_selected"]).sum())
        print(f"\n최대 차이 {worst:.3e} · 별 수 일치 {star_match}/{len(frame)}")
        print("판정:", "재현됨" if worst < 1e-6 and star_match == len(frame)
              else "불일치 — 원인 규명 전에는 이 대조를 쓰지 않는다")
        if args.output:
            frame.to_csv(args.output, index=False)
            args.output.with_suffix(".inputs.json").write_text(json.dumps({
                "workspace": str(args.workspace),
                "shared_with_apex": ["photutils aperture summation", "star positions (x_fit)",
                                     "per-star sky", "radius grid", "fixed-sky policy"],
                "independent_here": ["reimplemented selection, normalisation, "
                                     "median stacking, radius lookup"],
                "r_ap_scale": r_ap_scale, "ref_ap_scale": ref_scale,
                "min_r_ap_px": min_r_ap, "gc_n_steps": GC_N_STEPS,
                "apcorr_min_snr": APCORR_MIN_SNR, "apcorr_max_sources": APCORR_MAX_SOURCES,
            }, indent=2, ensure_ascii=False), encoding="utf-8")
            print(f"표 -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
