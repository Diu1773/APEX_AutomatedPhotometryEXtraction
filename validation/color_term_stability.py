"""1 차 색항이 적합에 안 쓴 보정성에서도 맞는가 — 밖을 안 보고 하는 검사.

**MuSCAT3 의 r 밴드에서 이것이 실제로 틀어져 있었다** (2026-09-08). 두 기기의
r 등급이 SNR 50~100 에서 0.27 어긋나는 원인을 좇다가 나왔다.

    적합이 쓴 별로 (SNR>=20 · Gaia 품질 · 3σ 절단, 153 별)   ct = -0.204
    보정성 표 전체로 (241 별)                                 ct = +0.021

**부호가 반대다.** 그리고 −0.204 를 적용한 결과, 보정된 등급에 색 기울기 −0.26 이
남는다 — 없애려던 것과 비슷한 크기를 반대로 넣은 셈이다. Moravian 의 같은 밴드는
어느 절단에서도 +0.03 언저리로 흔들리지 않는다.

## 왜 이 검사가 필요한가

적합은 자기가 고른 별에서는 언제나 잔차 기울기가 0 이다. 최소제곱이 그렇게
만든다. 그러니 **적합 안의 어떤 통계도 이 고장을 못 잡는다** — 산포도, 유의도도,
결정계수도 다 좋게 나온다. 실제로 MuSCAT3 r 의 적합 산포는 0.029 로 세 밴드 중
제일 나쁘긴 했지만 「나쁘다」와 「부호가 틀렸다」는 다른 말이다.

**밖을 볼 필요는 없다.** 적합이 버린 별들이 그 자리에 있다. 버린 별에서 기울기가
남으면, 적합이 고른 부분집합이 전체를 대표하지 못한 것이다.

## 재는 것

시야마다 · 밴드마다 셋을 낸다.

    적합 ct       적합이 쓴 별로 다시 맞춘 값 (기록된 것과 맞는지 확인용)
    전체 ct       보정성 표 전체로 맞춘 값
    차이          |전체 − 적합|. **이것이 크면 적합이 전체를 대표하지 않는다.**

**판정선을 미리 정하지 않는다.** 여덟 시야에서 이 값이 어떻게 분포하는지 보고
MuSCAT3 r 이 실제로 튀는지부터 확인한다.

실행:
    python -X utf8 validation/color_term_stability.py
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

from apex.analysis.cmd.zeropoint_runner import robust_weighted_polyfit  # noqa: E402
from apex.utils.gaia_quality import gaia_quality_report  # noqa: E402

ROOTS = ("E:/APEX_validation/reprocess/*/result",
         "E:/observed_Analysis/*/*/result")
EXTRA = {"MuSCAT3": Path("validation/external_muscat3/results"),
         "kb26": Path("validation/external_kb26/results")}

KIND = {"M13": "구상", "M3": "구상", "M5": "구상",
        "M67": "산개", "NGC6811": "산개", "M37": "산개", "NGC457": "산개",
        "MuSCAT3": "산개", "kb26": "산개"}   # 둘 다 M67 을 다른 기기로 찍은 것

#: APEX 가 밴드마다 쓰는 색지수 (`_FILTER_COLOR_PREF` 의 첫 후보).
BAND_COLOUR = {"g": "g_r", "r": "g_r", "i": "r_i",
               "B": "B_V", "V": "B_V", "R": "V_R", "I": "V_R"}

SNR_CUT = 20.0   # gaia_snr_calib_min 의 기본값


def field_of(path: str) -> str:
    for name in KIND:
        if re.search(rf"[\\/]{name}[\\/]", path):
            return name
    return "?"


def linear_ct(x, y, w, m):
    if int(m.sum()) < 20:
        return float("nan"), 0
    c, n, _ = robust_weighted_polyfit(x[m], y[m], w=w[m], degree=1,
                                      clip_sigma=3.0, iters=5, min_n=10)
    if c is None or not np.all(np.isfinite(c)):
        return float("nan"), 0
    return float(c[0]), int(n)


def rows_for(field: str, result_dir: Path) -> list[dict]:
    zp_dir = result_dir / "cmd_zeropoint"
    cal_p = zp_dir / "gaia_sdss_calibrator_by_ID.csv"
    if not cal_p.exists():
        return []
    try:
        d = pd.read_csv(cal_p)
    except Exception:  # noqa: BLE001
        return []
    # **색축은 산출물에서 읽는다.** APEX 는 밴드마다 후보 목록(`FILTER_COLOR_PREF`)
    # 을 순서대로 보고 자료에 있는 첫 색을 쓴다. 그래서 같은 r 밴드라도 g 가 있으면
    # g−r, 없으면 r−i 로 맞춘다 — kb26 이 그 경우다. 여기에 표를 박아 두면 그런
    # 워크스페이스를 조용히 건너뛴다(실제로 한 번 그렇게 건너뛰었다).
    recorded: dict[str, float] = {}
    axis: dict[str, str] = {}
    co_p = zp_dir / "zp_fit_coefficients.csv"
    if co_p.exists():
        try:
            co = pd.read_csv(co_p).set_index("filter")
            recorded = {str(k): float(v) for k, v in co["ct"].items()}
            if "color_col" in co.columns:
                axis = {str(k): str(v) for k, v in co["color_col"].items()
                        if isinstance(v, str) or pd.notna(v)}
        except Exception:  # noqa: BLE001
            pass
    try:
        qual = np.asarray(gaia_quality_report(d, cstar_nsigma=None)[0], bool)
    except Exception:  # noqa: BLE001
        qual = np.ones(len(d), bool)

    out = []
    for band, fallback in BAND_COLOUR.items():
        # 산출물이 적어 둔 축이 정본. 없을 때만 표의 첫 후보로 간다.
        colour = axis.get(band, fallback)
        dc, cc, ec, sc = (f"delta_{band}", f"color_{colour}",
                          f"mag_inst_err_{band}", f"snr_{band}")
        if dc not in d.columns or ec not in d.columns:
            continue
        if cc not in d.columns:      # 적어 둔 축의 열이 표에 없으면 물러선다
            colour = fallback
            cc = f"color_{colour}"
            if cc not in d.columns:
                continue
        x = pd.to_numeric(d[cc], errors="coerce").to_numpy(float)
        y = pd.to_numeric(d[dc], errors="coerce").to_numpy(float)
        e = pd.to_numeric(d[ec], errors="coerce").to_numpy(float)
        s = (pd.to_numeric(d[sc], errors="coerce").to_numpy(float)
             if sc in d.columns else np.full(len(d), np.inf))
        w = 1.0 / np.maximum(e, 1e-6) ** 2

        base = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0)
        m_fit = base & np.isfinite(s) & (s >= SNR_CUT) & qual
        if int(base.sum()) < 40 or int(m_fit.sum()) < 20:
            continue

        ct_fit, n_fit = linear_ct(x, y, w, m_fit)
        ct_all, n_all = linear_ct(x, y, w, base)
        out.append(dict(field=field, kind=KIND.get(field, "?"), band=band,
                        colour=colour, ct_recorded=recorded.get(band, float("nan")),
                        ct_fit=ct_fit, n_fit=n_fit, ct_all=ct_all, n_all=n_all,
                        n_input=int(base.sum()),
                        gap=abs(ct_all - ct_fit)))
    return out


def main() -> int:
    rows: list[dict] = []
    seen: set[str] = set()
    for pat in ROOTS:
        for rd in sorted(glob.glob(pat)):
            f = field_of(str(Path(rd) / "x"))
            if f == "?" or f in seen:
                continue
            r = rows_for(f, Path(rd))
            if r:
                seen.add(f)
                rows += r
    for name, rd in EXTRA.items():
        rows += rows_for(name, rd)

    d = pd.DataFrame(rows)
    if d.empty:
        print("적합할 자료를 못 찾았다.")
        return 1

    print("=== 적합이 쓴 별로 잰 색항 vs 보정성 전체로 잰 색항 ===")
    print("{:<9}{:<4}{:<6}{:>10}{:>10}{:>6}{:>10}{:>6}{:>9}".format(
        "시야", "밴드", "색", "기록된 ct", "적합 ct", "N", "전체 ct", "N", "차이"))
    print("-" * 70)
    for _, r in d.sort_values("gap", ascending=False).iterrows():
        print("{:<9}{:<4}{:<6}{:>+10.4f}{:>+10.4f}{:>6}{:>+10.4f}{:>6}{:>9.4f}".format(
            r["field"], r["band"], r["colour"], r["ct_recorded"], r["ct_fit"],
            r["n_fit"], r["ct_all"], r["n_all"], r["gap"]))

    g = d["gap"].to_numpy(float)
    g = g[np.isfinite(g)]
    print()
    print("=== 차이의 분포 ===")
    print(f"  조합 {len(g)} 개 · 중앙 {np.median(g):.4f} · 90 % 지점 "
          f"{np.percentile(g, 90):.4f} · 최대 {g.max():.4f}")
    big = d[d["gap"] > 0.10]
    print(f"  0.10 을 넘는 조합 {len(big)} 개"
          + ("" if big.empty else ": "
             + ", ".join(f"{r['field']} {r['band']} ({r['gap']:.3f})"
                         for _, r in big.iterrows())))
    print()
    print("  ※ 판정선을 정하려는 표가 아니다. 이 차이가 어디에 몰려 있는지 보고,")
    print("     MuSCAT3 r 이 실제로 튀는지부터 확인하는 것이 목적이다.")

    out = Path("validation/color_term_stability.json")
    out.write_text(json.dumps(d.to_dict("records"), ensure_ascii=False,
                              indent=2, default=str), encoding="utf-8")
    print()
    print(f"썼다: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
