"""Union of every filter's picks, or one filter's screened set?

The Step 8 window screens each filter separately and writes
`selection_<filter>.json` for each. What it does *downstream* is not per-filter:
`_load_selection_ids` takes the **union** when the sets differ, and Step 9 builds
one curve from that. On YZ Boo that is eleven stars against the six the pipeline
carries from the busiest filter — and the union necessarily contains stars that
some other filter's screening threw out.

(An earlier version of this script compared "each filter's own set" against
"g's set". Neither program does that: the window never uses a per-filter set to
build the curve. The union is what it actually uses, so that is the comparison.)

The target cannot answer it — YZ Boo swings 1.2 mag, so the scatter of its
differential curve is its own pulsation, not the ensemble's noise. A star that
does *not* vary can: its differential curve should be flat, and whatever scatter
is left is what the ensemble put there. So every star in the union of the three
ensembles is used as a stand-in check star, held out of the ensemble while it is
being measured.

Scatter is measured after removing a per-night offset, because that is what
Step 10 does and the question is what survives it.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).absolute().parent.parent))

import numpy as np
import pandas as pd

from apex.analysis.light_curve.photometry_source_service import (
    load_filter_photometry_timeseries,
)
from apex.utils.step_paths import master_catalog_path

RESULT = Path(r"E:\APEX_validation\reprocess\YZBoo_2n\result")
SELECTION = RESULT / "lc_selection"
MAD = 1.4826
NIGHT_GAP_DAYS = 0.3


def window_ensembles() -> dict[str, list[int]]:
    out = {}
    for path in sorted(SELECTION.glob("selection_*.json")):
        body = json.loads(path.read_text(encoding="utf-8"))
        flt = str(body.get("filter") or path.stem.split("_", 1)[-1])
        out[flt] = [int(v) for v in (body.get("comparison_ids") or [])]
    return out


def id_to_source_id() -> dict[int, int]:
    catalog = pd.read_csv(master_catalog_path(RESULT), sep="\t")
    frame = catalog[["ID", "source_id"]].apply(pd.to_numeric, errors="coerce").dropna()
    return {int(r.ID): int(r.source_id) for r in frame.itertuples()}


def nightly_flattened_scatter(frame: pd.DataFrame) -> float:
    """Robust scatter of a differential curve after a per-night offset."""
    if frame.empty:
        return float("nan")
    times = frame["time"].to_numpy(float)
    order = np.argsort(times)
    times, values = times[order], frame["diff"].to_numpy(float)[order]
    night = np.concatenate([[0], np.cumsum(np.diff(times) > NIGHT_GAP_DAYS)])
    flat = np.empty_like(values)
    for label in np.unique(night):
        m = night == label
        flat[m] = values[m] - np.median(values[m])
    good = flat[np.isfinite(flat)]
    if good.size < 10:
        return float("nan")
    return float(MAD * np.median(np.abs(good - np.median(good))))


def differential(measurements: pd.DataFrame, star: int,
                 ensemble: list[int]) -> pd.DataFrame:
    """`star` minus the inverse-variance weighted mean of `ensemble`, per frame.

    The same weighting `raw_lightcurve` uses, so the number means what the
    pipeline's `comp_avg` means.
    """
    members = [s for s in ensemble if s != star]
    if len(members) < 2:
        return pd.DataFrame(columns=["time", "diff"])
    wanted = measurements[measurements["star_id"].isin(members + [star])]
    rows = []
    for frame_name, group in wanted.groupby("frame"):
        target = group[group["star_id"] == star]
        comps = group[group["star_id"].isin(members)]
        if target.empty or comps.empty:
            continue
        mag = float(target["mag"].iloc[0])
        cm = comps["mag"].to_numpy(float)
        ce = comps["mag_err"].to_numpy(float)
        ok = np.isfinite(cm) & np.isfinite(ce) & (ce > 0)
        if not ok.any() or not np.isfinite(mag):
            continue
        w = 1.0 / ce[ok] ** 2
        rows.append({"time": float(group["time"].iloc[0]) if "time" in group else np.nan,
                     "frame": frame_name,
                     "diff": mag - float(np.sum(cm[ok] * w) / np.sum(w))})
    out = pd.DataFrame(rows)
    if not out.empty and out["time"].isna().all():
        out["time"] = np.arange(len(out), dtype=float)
    return out


def main() -> None:
    ensembles = window_ensembles()
    sid = id_to_source_id()
    donor = "g"                       # the filter the pipeline would pick
    union = sorted({v for ids in ensembles.values() for v in ids})
    print("  창 스텝 8 기록: " + " · ".join(f"{f}={len(v)}" for f, v in ensembles.items()))
    print(f"  창 스텝 9 가 실제로 쓰는 것: 합집합 {len(union)} 개 {union}")
    print(f"  파이프라인이 쓰는 것: {donor} 의 {len(ensembles[donor])} 개 "
          f"{sorted(ensembles[donor])}")
    print(f"  검사: 그 {len(union)} 개를 각각 체크성처럼 쓴다\n")
    print(f"  {'필터':<5}{'검사 별':>8}{'합집합 11':>11}{'선별 6':>10}"
          f"{'차이(mmag)':>12}{'판정':>16}")
    print("  " + "-" * 66)

    for flt in ("g", "r", "i"):
        measurements, _info = load_filter_photometry_timeseries(RESULT, flt, None)
        if measurements.empty:
            print(f"  {flt:<5} 측광 없음")
            continue
        own = [sid[i] for i in union if i in sid]                  # 창 스텝 9
        borrowed = [sid[i] for i in ensembles[donor] if i in sid]  # 파이프라인
        own_s, bor_s = [], []
        for display_id in union:
            star = sid.get(display_id)
            if star is None:
                continue
            a = nightly_flattened_scatter(differential(measurements, star, own))
            b = nightly_flattened_scatter(differential(measurements, star, borrowed))
            if np.isfinite(a) and np.isfinite(b):
                own_s.append(a)
                bor_s.append(b)
        if not own_s:
            print(f"  {flt:<5} 비교 가능한 별 없음")
            continue
        a, b = float(np.median(own_s)), float(np.median(bor_s))
        delta = (b - a) * 1000.0
        # delta = (6 개 산포) − (합집합 산포). 음수면 6 개 쪽이 덜 흩어진다.
        verdict = ("선별 6 개가 낫다" if delta < -1.0 else
                   "합집합이 낫다" if delta > 1.0 else "차이 없음")
        print(f"  {flt:<5}{len(own_s):>8}{a:>11.5f}{b:>10.5f}"
              f"{delta:>+12.2f}{verdict:>16}")

    print("\n  숫자는 밤별 영점을 뺀 뒤의 robust 산포 중앙값(작을수록 좋다). "
          "\n  차이 = 선별 6 개 − 합집합 11 개. 음수 = 선별한 6 개가 덜 흩어진다.")


if __name__ == "__main__":
    main()
