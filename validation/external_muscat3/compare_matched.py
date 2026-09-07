"""같은 별 하나의 집합으로 세 쌍을 견준다 — 변인통제.

**왜 다시 재나.** 앞선 비교(`compare_magnitudes.py`)는 쌍마다 표본이 달랐다 —
MuSCAT3−Moravian 은 211 별, MuSCAT3−PS1 은 244 별, Moravian−PS1 은 409 별.
그러면 중앙값끼리 더하고 빼도 맞지 않는다. 실제로 r 밴드에서
MuSCAT3−Moravian 이 −0.019 인데 PS1 을 거쳐 가면 −0.057 이 나와 세 배 어긋났다.

**세 곳 모두에 있는 별만 남긴다.** 그러면 `A−B = (A−PS1) − (B−PS1)` 가 항등식으로
성립하고, 어긋남이 표본 차이인지 진짜 차이인지 헷갈릴 일이 없다.

## 여기서 내는 것

    밴드마다 · 같은 별 · 세 쌍
      중앙 차이와 MAD
      등급 구간별 차이      — 어두운 쪽에서 휘는가
      SNR 구간별 차이       — 휘는 것이 밝기 때문인가 신호대잡음 때문인가
      색 (g−i) 기울기       — 색항이 남았는가

**절대값과 흩어짐을 함께 낸다.** 「상관 제거」는 전부 균일하게 나쁘게 만들어도
달성되므로 절대값을 늘 같이 본다.

실행:
    python -X utf8 validation/external_muscat3/compare_matched.py
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).absolute().parents[2]))

from compare_magnitudes import (  # noqa: E402
    BANDS, join, load_workspace, mad, query_ps1,
)


def _bin_stats(delta, axis, edges, label):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = np.isfinite(delta) & np.isfinite(axis) & (axis >= lo) & (axis < hi)
        if m.sum() < 5:
            continue
        out.append({label: f"{lo:g}~{hi:g}", "n": int(m.sum()),
                    "median": float(np.median(delta[m])), "mad": mad(delta[m])})
    return out


def main(argv: list[str] | None = None) -> int:
    here = Path(__file__).absolute().parent
    ap = argparse.ArgumentParser(description="같은 별로 세 쌍을 견준다")
    ap.add_argument("--muscat3", default=str(here / "results"))
    ap.add_argument("--moravian",
                    default="E:/APEX_validation/reprocess/M67/result")
    ap.add_argument("--out", default=str(here / "magnitude_matched.json"))
    a = ap.parse_args(argv)

    m3 = load_workspace(Path(a.muscat3), "MuSCAT3")
    mv = load_workspace(Path(a.moravian), "Moravian")

    ra0, dec0 = float(np.nanmedian(m3["ra_deg"])), float(np.nanmedian(m3["dec_deg"]))
    span = float(np.nanmax(np.hypot(
        (m3["ra_deg"] - ra0) * np.cos(np.deg2rad(dec0)), m3["dec_deg"] - dec0)))
    ps1 = query_ps1(ra0, dec0, span + 0.02)
    if ps1 is None:
        raise SystemExit("PS1 을 못 받아서 변인통제 비교를 할 수 없다.")
    ps1["gaia_source_id"] = pd.NA

    # MuSCAT3 를 축으로 두 상대를 차례로 붙인다 — 그러면 남는 것이 교집합이다.
    step1 = join(m3, mv)
    step1["ra_deg"] = step1["ra_deg_A"]
    step1["dec_deg"] = step1["dec_deg_A"]
    step1["gaia_source_id"] = step1["gaia_source_id_A"]
    both = join(step1, ps1)
    print(f"\n세 곳 모두에 있는 별: {len(both)}")

    cols = {
        "MuSCAT3": "mag_cal_{b}_A_A",
        "Moravian": "mag_cal_{b}_B_A",
        "PS1": "mag_cal_{b}_B",
    }
    missing = [n for n, t in cols.items() if t.format(b="g") not in both.columns]
    if missing:
        raise SystemExit(f"열을 못 찾았다: {missing}\n있는 열: {list(both.columns)[:40]}")

    summary = {"generated": time.strftime("%Y-%m-%dT%H:%M:%S"),
               "n_common": int(len(both)), "bands": {}}

    print()
    print("같은 별 · 같은 표본으로 다시 잰다")
    print(f"{'밴드':<5}{'N':>5}"
          f"{'M3−Mor':>10}{'MAD':>8}"
          f"{'M3−PS1':>10}{'MAD':>8}"
          f"{'Mor−PS1':>10}{'MAD':>8}{'항등식':>9}")
    print("-" * 76)
    for b in BANDS:
        m3v = pd.to_numeric(both[cols["MuSCAT3"].format(b=b)], errors="coerce").to_numpy(float)
        mvv = pd.to_numeric(both[cols["Moravian"].format(b=b)], errors="coerce").to_numpy(float)
        p1v = pd.to_numeric(both[cols["PS1"].format(b=b)], errors="coerce").to_numpy(float)
        m = np.isfinite(m3v) & np.isfinite(mvv) & np.isfinite(p1v)
        if m.sum() < 10:
            print(f"{b:<5}{int(m.sum()):>5}   별이 모자라 건너뛴다")
            continue
        d_mm = m3v[m] - mvv[m]
        d_m1 = m3v[m] - p1v[m]
        d_v1 = mvv[m] - p1v[m]
        ident = float(np.median(d_mm) - (np.median(d_m1) - np.median(d_v1)))
        print(f"{b:<5}{int(m.sum()):>5}"
              f"{np.median(d_mm):>+10.4f}{mad(d_mm):>8.4f}"
              f"{np.median(d_m1):>+10.4f}{mad(d_m1):>8.4f}"
              f"{np.median(d_v1):>+10.4f}{mad(d_v1):>8.4f}{ident:>+9.4f}")

        gi = None
        gA = pd.to_numeric(both[cols["MuSCAT3"].format(b="g")], errors="coerce").to_numpy(float)
        iA = pd.to_numeric(both[cols["MuSCAT3"].format(b="i")], errors="coerce").to_numpy(float)
        if np.isfinite(gA[m]).any() and np.isfinite(iA[m]).any():
            gi = (gA - iA)[m]
        snr_col = f"snr_{b}_A_A"
        snr = (pd.to_numeric(both[snr_col], errors="coerce").to_numpy(float)[m]
               if snr_col in both.columns else None)

        entry = {"n": int(m.sum())}
        for name, d in (("m3_minus_mor", d_mm), ("m3_minus_ps1", d_m1),
                        ("mor_minus_ps1", d_v1)):
            e = {"median": float(np.median(d)), "mad": mad(d)}
            e["by_mag"] = _bin_stats(d, m3v[m], np.arange(12, 21.5, 1.0), "mag")
            if snr is not None:
                e["by_snr"] = _bin_stats(d, snr, np.array([0, 20, 50, 100, 300, 1e9]), "snr")
            if gi is not None and np.isfinite(gi).sum() > 10:
                gm = np.isfinite(gi) & np.isfinite(d)
                e["colour_slope"] = float(np.polyfit(gi[gm], d[gm], 1)[0])
            entry[name] = e
        summary["bands"][b] = entry

    print()
    print("항등식 = (M3−Mor) − [(M3−PS1) − (Mor−PS1)].  같은 별을 썼으면 0 이어야 한다.")

    # 밝기·SNR 구간표
    for b in BANDS:
        e = summary["bands"].get(b)
        if not e:
            continue
        print()
        print(f"--- {b} 밴드 · MuSCAT3 − Moravian 이 밝기에 따라 어떻게 변하나 ---")
        print(f"{'등급':>10}{'N':>6}{'중앙':>10}{'MAD':>9}")
        for row in e["m3_minus_mor"]["by_mag"]:
            print(f"{row['mag']:>10}{row['n']:>6}{row['median']:>+10.4f}{row['mad']:>9.4f}")

    Path(a.out).write_text(json.dumps(summary, ensure_ascii=False, indent=2),
                           encoding="utf-8")
    print(f"\n썼다: {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
