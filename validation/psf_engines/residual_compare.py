"""What is left in the pixels after each engine subtracts its model.

The statistical decomposition ran out of room at 274 stars. Model, interpolation
order, group size, fit window, background level, per-star sky, pixel phase,
iteration count and position error were each tested and closed; the one effect
that survives multiple-testing correction is confined to the tightest blends,
where ALLSTAR's magnitudes are 41 mmag closer to truth. Slicing that sample
further only finds noise.

Both engines write the frame with their fitted stars removed — ALLSTAR to
`sub.fits`, APEX to `residual_iter*.fits` — so the question can be asked of the
pixels instead. If APEX is leaving a neighbour's light behind, it shows up as
structure at the implanted star's position that ALLSTAR's residual does not
have.

Two things make the images not directly comparable, and both are handled here:
ALLSTAR keeps the sky in its residual while APEX subtracts it first, so each
image is offset to a zero median measured in blank field; and the two engines
fit different star lists, so only positions where both actually removed
something are compared.

Reported per crowding bin, in a fixed aperture at each implanted star:
the residual sum (light left behind, in units of the star's own flux) and the
residual scatter (how structured what is left is).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits

HERE = Path(__file__).absolute().parent
sys.path.insert(0, str(HERE))

from compare_recovery import robust_scatter  # noqa: E402


def sky_zeroed(path: Path, boxes: int = 400, half: int = 24,
               seed: int = 0) -> np.ndarray:
    """Load a residual image and put its blank-field level at zero.

    ALLSTAR leaves the sky in; APEX removes it beforehand. The offset is taken
    as the median of many small boxes' medians, which is insensitive to the
    stars still present in either image.
    """
    data = fits.getdata(path).astype(float)
    rng = np.random.default_rng(seed)
    ny, nx = data.shape
    ys = rng.integers(half, ny - half, boxes)
    xs = rng.integers(half, nx - half, boxes)
    levels = [np.median(data[y - half:y + half, x - half:x + half])
              for y, x in zip(ys, xs)]
    return data - float(np.median(levels))


def aperture_stats(image: np.ndarray, x: float, y: float,
                   radius: float) -> tuple[float, float]:
    """Sum and robust scatter of the residual inside one aperture."""
    xi, yi = int(round(x)), int(round(y))
    r = int(np.ceil(radius))
    if xi - r < 0 or yi - r < 0 or xi + r + 1 > image.shape[1] or yi + r + 1 > image.shape[0]:
        return float("nan"), float("nan")
    cut = image[yi - r:yi + r + 1, xi - r:xi + r + 1]
    yy, xx = np.mgrid[:cut.shape[0], :cut.shape[1]]
    inside = np.hypot(xx - (x - (xi - r)), yy - (y - (yi - r))) <= radius
    values = cut[inside]
    return float(values.sum()), robust_scatter(values)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--work", default=r"E:\APEX_validation\psf_engines\M13_ast2x200_moffat")
    ap.add_argument("--frame", default="pp_messier13-0005-B.fit")
    ap.add_argument("--apex-run", default="apexhybrid")
    ap.add_argument("--trials", type=int, default=2)
    ap.add_argument("--radius-fwhm", type=float, default=1.0)
    ap.add_argument("--fwhm-px", type=float, default=7.052391)
    ap.add_argument("--output", default=str(HERE / "residual_compare.csv"))
    args = ap.parse_args()

    work = Path(args.work)
    radius = args.radius_fwhm * args.fwhm_px
    truth_all = pd.read_csv(work / "truth.csv")

    rows = []
    for trial in range(1, args.trials + 1):
        truth = truth_all[truth_all["trial"] == trial].reset_index(drop=True)
        dao_path = work / f"daophot_work_trial{trial}" / "sub.fits"
        apex_dir = work / f"{args.apex_run}_trial{trial}" / "result" / "cmd_psf"
        # `residual_<frame>` is written after the final pass; the numbered
        # `residual_iter*` files are snapshots of the residual passes and do not
        # include what the last fit removed. Comparing an intermediate APEX
        # residual against ALLSTAR's final one would blame APEX for light it
        # subtracts a step later.
        final = apex_dir / f"residual_{args.frame}"
        apex_candidates = ([final] if final.exists()
                           else sorted(apex_dir.glob(f"residual_iter*_{args.frame}")))
        if truth.empty or not dao_path.exists() or not apex_candidates:
            print(f"[건너뜀] trial {trial}")
            continue
        apex_path = apex_candidates[-1]
        print(f"  trial {trial}: APEX {apex_path.name} · ALLSTAR sub.fits")
        dao = sky_zeroed(dao_path)
        apex = sky_zeroed(apex_path)

        for _, star in truth.iterrows():
            x, y = float(star["x_true"]), float(star["y_true"])
            flux = float(star["flux_realized_adu"])
            a_sum, a_sct = aperture_stats(apex, x, y, radius)
            d_sum, d_sct = aperture_stats(dao, x, y, radius)
            rows.append({
                "trial": trial, "crowding_bin": star["crowding_bin"],
                "target_snr": star["target_snr"], "flux_adu": flux,
                "apex_left_frac": a_sum / flux, "dao_left_frac": d_sum / flux,
                "apex_scatter": a_sct, "dao_scatter": d_sct,
            })

    if not rows:
        raise SystemExit("비교할 회차가 없다")
    d = pd.DataFrame(rows)
    d.to_csv(args.output, index=False)

    print(f"\n인공별 자리에 남은 빛 (별 자신의 밝기 대비) · 반경 {args.radius_fwhm:g}·FWHM")
    print(f"{'혼잡 구간':>16}{'N':>5}{'APEX 남은빛':>13}{'ALLSTAR 남은빛':>16}"
          f"{'APEX 잔차산포':>15}{'ALLSTAR 잔차산포':>17}")
    print("-" * 84)
    for label in list(d["crowding_bin"].unique()) + ["전체"]:
        part = d if label == "전체" else d[d["crowding_bin"] == label]
        part = part[np.isfinite(part["apex_left_frac"]) & np.isfinite(part["dao_left_frac"])]
        if len(part) < 8:
            continue
        print(f"{label:>16}{len(part):>5}"
              f"{np.median(part['apex_left_frac']):>+13.4f}"
              f"{np.median(part['dao_left_frac']):>+16.4f}"
              f"{np.median(part['apex_scatter']):>15.2f}"
              f"{np.median(part['dao_scatter']):>17.2f}")
    print("\n남은 빛이 양수면 모형이 별을 덜 뺀 것이고, 음수면 더 뺀 것이다.\n"
          "잔차 산포가 크면 뺀 자리에 구조가 남아 있다는 뜻이다.")
    print(f"saved -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
