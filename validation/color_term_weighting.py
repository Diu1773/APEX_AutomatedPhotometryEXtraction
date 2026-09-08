"""2 차 색항이 시야마다 흔들리는 것이 가중치 탓인가, 절단 탓인가, 자료 탓인가.

**D-012 의 근거를 다시 세우는 측정이다** (2026-08-21). 그 결정은 한 카메라의
다섯 시야에서 R 의 ct2 가 부호까지 뒤집히는 것을 보고 「시야마다 적합해 쓰지
말자」로 닫았다. 그 값들은 전부 가중 적합(1/σ²)에서 나온 것이라, 흔들림이 하늘의
성질인지 적합기의 성질인지 갈라 놓지 않았다.

## 왜 가중치를 의심하게 됐나

`robust_weighted_polyfit` 은 **가중으로 곡선을 맞추면서 절단은 무가중 MAD 로**
한다(126~137 줄). 곡선이 밝은 별 쪽으로 끌려가면 어두운 별에 큰 잔차가 남고,
그 잔차를 무가중 잣대로 자르니 어두운 별이 떨어져 나간다. 떨어져 나가면 색의
지렛대가 짧아지고 곡선은 더 끌려간다. **되먹임이다.**

MuSCAT3 의 rp 에서 이것이 실제로 보였다 — 별이 184 개에서 105 개로 줄고 ct2 가
+0.019 에서 −0.985 로 뒤집혔다.

## 네 가지를 같은 별로 잰다

앞선 판에는 결함이 있었다. 무가중 쪽은 오차 열이 비어도 별을 세고 가중 쪽은
못 세니 **두 N 이 애초에 다른 별을 봤다.** 여기서는 색·차이·오차가 모두 있는
별만 남기고 그 하나의 집합에 네 적합을 돌린다.

    A. 무가중 · 절단 있음      기준선
    B. 가중  · 절단 있음      APEX 가 실제로 쓰는 것
    C. 가중  · 절단 없음      가중치만의 효과
    D. 무가중 · 절단 없음      아무것도 안 한 날것

그러면 **D→C 가 가중치의 몫이고 C→B 가 절단 되먹임의 몫이다.**

**결론을 미리 정하지 않는다.** 넷이 다 비슷하게 흔들리면 흔들림은 적합기가 아니라
자료에서 온 것이고, D-012 는 그대로 선다.

실행:
    python -X utf8 validation/color_term_weighting.py
"""
from __future__ import annotations

import glob
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).absolute().parents[1]))

from apex.analysis.cmd.zeropoint_runner import (  # noqa: E402
    _quad_coefficient_sigma,
    robust_weighted_polyfit,
)

ROOTS = ("E:/APEX_validation/reprocess/*/result",
         "E:/observed_Analysis/*/*/result")
EXTRA = {"MuSCAT3": Path("validation/external_muscat3/results"),
         "kb26": Path("validation/external_kb26/results")}

KIND = {"M13": "구상", "M3": "구상", "M5": "구상",
        "M67": "산개", "NGC6811": "산개", "M37": "산개", "NGC457": "산개",
        "MuSCAT3": "산개", "kb26": "산개"}   # 둘 다 M67 을 다른 기기로 찍은 것

BAND_COLOUR = {"g": "g_r", "r": "g_r", "i": "r_i",
               "B": "B_V", "V": "B_V", "R": "V_R", "I": "V_R"}

#: 네 적합. (이름, 가중치 쓰나, 절단 횟수)
VARIANTS = (("A_무가중_절단", False, 5),
            ("B_가중_절단", True, 5),
            ("C_가중_생", True, 1),
            ("D_무가중_생", False, 1))


def field_of(path: str) -> str:
    for name in KIND:
        if re.search(rf"[\\/]{name}[\\/]", path):
            return name
    return "?"


def one_fit(x, y, w, use_w: bool, iters: int):
    """(ct2, 1σ, 남은 별 수). x·y·w 는 이미 같은 집합으로 걸러진 것."""
    ww = w if use_w else np.ones_like(w)
    c, n, _ = robust_weighted_polyfit(x, y, w=ww, degree=2, clip_sigma=3.0,
                                      iters=iters, min_n=10)
    if c is None or not np.all(np.isfinite(c)):
        return np.nan, np.nan, 0
    return float(c[0]), float(_quad_coefficient_sigma(x, y, ww, c)), int(n)


def rows_for(field: str, result_dir: Path) -> list[dict]:
    p = result_dir / "cmd_zeropoint" / "gaia_sdss_calibrator_by_ID.csv"
    if not p.exists():
        return []
    try:
        d = pd.read_csv(p)
    except Exception:  # noqa: BLE001
        return []
    out = []
    for band, colour in BAND_COLOUR.items():
        dc, cc, ec = f"delta_{band}", f"color_{colour}", f"mag_inst_err_{band}"
        if dc not in d.columns or cc not in d.columns or ec not in d.columns:
            continue
        x = pd.to_numeric(d[cc], errors="coerce").to_numpy(float)
        y = pd.to_numeric(d[dc], errors="coerce").to_numpy(float)
        e = pd.to_numeric(d[ec], errors="coerce").to_numpy(float)

        # **같은 별 하나의 집합.** 오차가 없는 별은 가중 적합이 못 쓰므로
        # 무가중 쪽에서도 뺀다 — 안 그러면 두 적합이 다른 자료를 본다.
        m = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0)
        if int(m.sum()) < 20:
            continue
        xs, ys = x[m], y[m]
        w = 1.0 / np.maximum(e[m], 1e-6) ** 2

        row = dict(field=field, kind=KIND.get(field, "?"), band=band,
                   colour=colour, n_input=int(m.sum()),
                   colour_span=float(np.percentile(xs, 95) - np.percentile(xs, 5)),
                   err_ratio=float(np.percentile(e[m], 95) / np.percentile(e[m], 5)))
        for name, use_w, iters in VARIANTS:
            ct2, sig, n = one_fit(xs, ys, w, use_w, iters)
            row[f"ct2_{name}"] = ct2
            row[f"sig_{name}"] = sig
            row[f"n_{name}"] = n
        out.append(row)
    return out


def main() -> int:
    rows: list[dict] = []
    seen: set[str] = set()
    for pat in ROOTS:
        for rd in sorted(glob.glob(pat)):
            fld = field_of(str(Path(rd) / "x"))
            if fld == "?" or fld in seen:
                continue
            r = rows_for(fld, Path(rd))
            if r:
                seen.add(fld)
                rows += r
    for name, rd in EXTRA.items():
        rows += rows_for(name, rd)

    d = pd.DataFrame(rows)
    if d.empty:
        print("적합할 자료를 못 찾았다.")
        return 1

    print("=== 같은 별을 네 가지로 적합 (괄호는 절단 뒤 남은 별) ===")
    print("{:<9}{:<4}{:<6}{:>6}{:>17}{:>17}{:>11}{:>11}".format(
        "시야", "밴드", "색", "별수", "A 무가중·절단", "B 가중·절단(APEX)",
        "C 가중·생", "D 무가중·생"))
    print("-" * 82)
    for _, r in d.sort_values(["band", "field"]).iterrows():
        print("{:<9}{:<4}{:<6}{:>6}{:>+11.3f}({:>3}){:>+11.3f}({:>3}){:>+11.3f}{:>11.3f}".format(
            r["field"], r["band"], r["colour"], r["n_input"],
            r["ct2_A_무가중_절단"], r["n_A_무가중_절단"],
            r["ct2_B_가중_절단"], r["n_B_가중_절단"],
            r["ct2_C_가중_생"], r["ct2_D_무가중_생"]))

    print()
    print("=== 어느 것이 흔들림을 만드나 — 밴드마다 시야 사이의 폭 ===")
    print("{:<5}{:<6}{:>5}{:>11}{:>11}{:>11}{:>11}{:>10}".format(
        "밴드", "색", "시야", "A 폭", "B 폭", "C 폭", "D 폭", "부호뒤집힘"))
    print("-" * 70)
    for (band, col), g in d.groupby(["band", "colour"]):
        vals, ptp, flip = {}, {}, {}
        for name, _, _ in VARIANTS:
            v = g[f"ct2_{name}"].to_numpy(float)
            v = v[np.isfinite(v)]
            vals[name] = v
            ptp[name] = float(v.max() - v.min()) if len(v) >= 2 else np.nan
            flip[name] = bool(len(v) >= 2 and v.max() > 0 and v.min() < 0)
        tag = "".join("ABCD"[i] for i, (n, _, _) in enumerate(VARIANTS) if flip[n]) or "없음"
        print("{:<5}{:<6}{:>5}{:>11.3f}{:>11.3f}{:>11.3f}{:>11.3f}{:>10}".format(
            band, col, len(vals["A_무가중_절단"]), ptp["A_무가중_절단"],
            ptp["B_가중_절단"], ptp["C_가중_생"], ptp["D_무가중_생"], tag))

    print()
    print("=== 절단 되먹임이 실제로 도나 — 가중 적합이 버리는 별 ===")
    keep_b = d["n_B_가중_절단"] / d["n_input"]
    keep_a = d["n_A_무가중_절단"] / d["n_input"]
    print(f"  가중·절단이 남기는 비율   중앙 {keep_b.median():.2f} · 최소 {keep_b.min():.2f}")
    print(f"  무가중·절단이 남기는 비율 중앙 {keep_a.median():.2f} · 최소 {keep_a.min():.2f}")
    worst = d.loc[keep_b.nsmallest(4).index]
    print("  가중 절단이 가장 많이 버린 넷:")
    for _, r in worst.iterrows():
        print(f"    {r['field']:<9}{r['band']:<3}"
              f" {r['n_input']:>4} -> {r['n_B_가중_절단']:>4}"
              f"   ct2 {r['ct2_D_무가중_생']:+.3f} -> {r['ct2_B_가중_절단']:+.3f}")

    print()
    print("=== 몫 나누기 — 무가중 날것에서 APEX 값까지 얼마씩 움직였나 ===")
    dw = (d["ct2_C_가중_생"] - d["ct2_D_무가중_생"]).abs()
    dc = (d["ct2_B_가중_절단"] - d["ct2_C_가중_생"]).abs()
    print(f"  가중치가 옮긴 몫    중앙 {dw.median():.3f} · 최대 {dw.max():.3f}")
    print(f"  절단이 옮긴 몫      중앙 {dc.median():.3f} · 최대 {dc.max():.3f}")

    out = Path("validation/color_term_weighting.json")
    out.write_text(json.dumps(d.to_dict("records"), ensure_ascii=False,
                              indent=2, default=str), encoding="utf-8")
    print()
    print(f"썼다: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
