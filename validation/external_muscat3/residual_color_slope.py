"""보정한 뒤에도 남는 색 기울기가 어디서 오나 — 두 조각으로 나눈다.

**풀리지 않은 것이 이것이다** (`MAGNITUDE_ERROR_ANALYSIS.md`). 두 기기의 r 등급이
SNR 50~100 에서 0.27 어긋나고, 같은 별 243 개에서 보정 뒤에도 색 기울기가 남는다.

    밴드   MuSCAT3 − PS1    Moravian − PS1
    g        +0.2045          +0.0908
    r        −0.1156          +0.0735      <- 부호가 반대다
    i        +0.0497          +0.0640

**이상한 점.** 두 워크스페이스는 **같은 기준**을 쓴다 — Gaia 등급을 SDSS 로 옮긴
`ref_g`·`ref_r`·`ref_i`. 그리고 APEX 는 1 차 색항을 적합해 적용하므로, 보정 뒤
등급에서 기준을 빼면 그 적합의 잔차가 되고 그 잔차는 적합에 쓴 색에 대해 기울기가
0 이다. 그러니 두 기기의 남은 기울기는 **같아야 한다.** 다른 것이 물음이다.

## 나누는 방법

항등식 하나로 나뉜다.

    (보정 등급 − PS1) = (보정 등급 − 기준) + (기준 − PS1)
      남은 기울기          적합의 잔차          기준 자체의 색 오차

**뒤쪽은 두 기기가 공유한다** — 같은 별의 같은 Gaia 등급을 같은 식으로 옮긴 값이다.
그러니 두 기기의 차이는 **앞쪽에서만** 올 수 있다.

앞쪽이 0 이 아닐 수 있는 이유가 있다. 적합은 그 워크스페이스의 **보정성 전체**로
하는데(MuSCAT3 241 별, Moravian 903 별), 여기서 재는 것은 **둘이 공유하는
부분집합**이다. 전체에서 기울기가 0 이어도 부분집합에서는 0 이 아닐 수 있다.

## 세우는 예측

MuSCAT3 는 보정성이 241 개이고 공유 부분집합이 그와 거의 같으므로 **적합의 잔차
기울기가 0 에 가까워야** 하고, 그러면 MuSCAT3 의 남은 기울기가 곧 **기준 자체의
색 오차**다. Moravian 은 903 개로 맞췄으므로 부분집합에서 잔차 기울기가 남을 수 있다.

**예측이 맞으면** g 의 `기준 − PS1` 기울기가 +0.20 근처로 나오고, 이것은 기기
문제가 아니라 **Gaia→SDSS 변환의 색 의존 오차**다.

**예측이 틀리면** 두 워크스페이스의 `ref` 가 같은 별에서 다르다는 뜻이고, 그것은
그것대로 찾아야 할 것이다. 그 검사도 같이 한다.

실행:
    python -X utf8 validation/external_muscat3/residual_color_slope.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).absolute().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))

from apex.utils.io_utils import read_csv_int64_source_id  # noqa: E402
from compare_magnitudes import BANDS, mad, query_ps1  # noqa: E402

M3 = HERE / "results"
MOR = Path("E:/APEX_validation/reprocess/M67/result")


def load(result_dir: Path, label: str) -> pd.DataFrame:
    """보정성 표 + 보정된 등급을 `ID` 로 붙인다."""
    cal = read_csv_int64_source_id(
        result_dir / "cmd_zeropoint" / "gaia_sdss_calibrator_by_ID.csv")
    mags = read_csv_int64_source_id(
        result_dir / "cmd_zeropoint" / "median_by_ID_filter_wide.csv")
    keep = ["ID"] + [f"mag_cal_{b}" for b in BANDS if f"mag_cal_{b}" in mags.columns]
    out = cal.merge(mags[keep], on="ID", how="left")
    n = {b: int(out[f"ref_{b}"].notna().sum()) for b in BANDS if f"ref_{b}" in out}
    print(f"[{label}] 보정성 {len(out)} 별 · 기준값 있는 별 {n}")
    return out


def match(a: pd.DataFrame, b: pd.DataFrame, tol_arcsec=1.0):
    """a 의 어느 줄이 b 의 어느 줄과 짝인지 **자리 번호로** 돌려준다.

    표를 바로 잘라 돌려주면 그다음에 세 번째 표를 붙일 때 「어느 줄이 살아
    남았나」를 다시 찾아야 해서, 좌표로 되짚는 위험한 코드가 생긴다.
    """
    from astropy.coordinates import SkyCoord

    ca = SkyCoord(a["ra_deg"].to_numpy(float), a["dec_deg"].to_numpy(float), unit="deg")
    cb = SkyCoord(b["ra_deg"].to_numpy(float), b["dec_deg"].to_numpy(float), unit="deg")
    k, sep, _ = ca.match_to_catalog_sky(cb)
    ok = np.flatnonzero(sep.arcsec <= tol_arcsec)
    return ok, np.asarray(k)[ok]


def slope(x, y, m):
    """**주어진 하나의 마스크로만** 잰다.

    조각마다 자기 유한값으로 재면 세 조각이 다른 별을 보게 되고, 그러면
    항등식이 성립하지 않아 분해가 뜻을 잃는다. 한 번 그렇게 재서 합이
    안 맞는 표를 냈다.
    """
    if int(m.sum()) < 15:
        return float("nan")
    return float(np.polyfit(x[m], y[m], 1)[0])


def main() -> int:
    m3, mv = load(M3, "MuSCAT3"), load(MOR, "Moravian")

    ra0 = float(np.nanmedian(m3["ra_deg"]))
    dec0 = float(np.nanmedian(m3["dec_deg"]))
    span = float(np.nanmax(np.hypot(
        (m3["ra_deg"] - ra0) * np.cos(np.deg2rad(dec0)), m3["dec_deg"] - dec0)))
    ps1 = query_ps1(ra0, dec0, span + 0.02)
    if ps1 is None:
        raise SystemExit("PS1 을 못 받아서 분해할 수 없다.")

    i_m3, i_mv = match(m3, mv)
    a = m3.iloc[i_m3].reset_index(drop=True)
    j_a, j_ps1 = match(a, ps1)
    # 세 표를 **같은 순서의 같은 별**로 맞춘다.
    a2 = a.iloc[j_a].reset_index(drop=True)
    b2 = mv.iloc[i_mv[j_a]].reset_index(drop=True)
    p = ps1.iloc[j_ps1].reset_index(drop=True)
    print(f"\n세 곳 모두에 있는 별 {len(a2)}")

    colour = (pd.to_numeric(p["mag_cal_g"], errors="coerce").to_numpy(float)
              - pd.to_numeric(p["mag_cal_i"], errors="coerce").to_numpy(float))

    print()
    print("=== 0. 두 워크스페이스의 기준값이 같은 별에서 같은가 ===")
    print("같은 Gaia 등급을 같은 식으로 옮긴 값이므로 0 이어야 한다.")
    print(f"{'밴드':<5}{'중앙 차이':>12}{'MAD':>10}{'|차이|>0.001 인 별':>18}")
    print("-" * 45)
    ref_same = True
    for band in BANDS:
        col = f"ref_{band}"
        if col not in a2.columns or col not in b2.columns:
            continue
        d = (pd.to_numeric(a2[col], errors="coerce").to_numpy(float)
             - pd.to_numeric(b2[col], errors="coerce").to_numpy(float))
        d = d[np.isfinite(d)]
        if d.size == 0:
            continue
        n_diff = int((np.abs(d) > 0.001).sum())
        if n_diff:
            ref_same = False
        print(f"{band:<5}{np.median(d):>+12.5f}{mad(d):>10.5f}{n_diff:>12} / {d.size}")
    print("→ 기준값은 " + ("같다." if ref_same else "**다르다. 이것부터 봐야 한다.**"))

    print()
    print("=== 1. 남은 기울기를 두 조각으로 나눈다 (PS1 의 g−i 에 대한 기울기) ===")
    print("항등식: (보정 − PS1) = (보정 − 기준) + (기준 − PS1)")
    print()
    print(f"{'밴드':<4}{'기기':<10}{'N':>5}{'보정−PS1':>11}{'보정−기준':>11}"
          f"{'기준−PS1':>11}{'합':>9}{'맞나':>7}")
    print("-" * 68)
    rows = []
    for band in BANDS:
        pc = f"mag_cal_{band}"
        if pc not in p.columns:
            continue
        ps = pd.to_numeric(p[pc], errors="coerce").to_numpy(float)
        for tag, t in (("MuSCAT3", a2), ("Moravian", b2)):
            cc, rc = f"mag_cal_{band}", f"ref_{band}"
            if cc not in t.columns or rc not in t.columns:
                continue
            cal = pd.to_numeric(t[cc], errors="coerce").to_numpy(float)
            ref = pd.to_numeric(t[rc], errors="coerce").to_numpy(float)
            m = (np.isfinite(colour) & np.isfinite(ps) & np.isfinite(cal)
                 & np.isfinite(ref))
            s_tot = slope(colour, cal - ps, m)
            s_fit = slope(colour, cal - ref, m)
            s_ref = slope(colour, ref - ps, m)
            ok = "예" if abs(s_tot - (s_fit + s_ref)) < 0.002 else "아니오"
            print(f"{band:<4}{tag:<10}{int(m.sum()):>5}{s_tot:>+11.4f}{s_fit:>+11.4f}"
                  f"{s_ref:>+11.4f}{s_fit + s_ref:>+9.4f}{ok:>7}")
            rows.append(dict(band=band, instrument=tag, n=int(m.sum()), total=s_tot,
                             fit_residual=s_fit, reference=s_ref))
        print()

    print("=== 1b. 잰 축을 바꿔 본다 — APEX 가 실제로 쓴 색으로 ===")
    print("APEX 는 g·r 을 g−r 로, i 를 r−i 로 없앴다. 그 축에서 잔차가 0 이면")
    print("적합은 제 할 일을 한 것이고, 남은 것은 색축 사이의 관계 문제다.")
    print()
    print(f"{'밴드':<4}{'기기':<10}{'적합 축':<7}{'N':>5}"
          f"{'그 축에서의 잔차 기울기':>22}{'PS1 g−i 에서':>14}")
    print("-" * 62)
    axis_rows = []
    for band in BANDS:
        ax_col = "color_r_i" if band == "i" else "color_g_r"
        for tag, t in (("MuSCAT3", a2), ("Moravian", b2)):
            cc, rc = f"mag_cal_{band}", f"ref_{band}"
            if cc not in t.columns or ax_col not in t.columns:
                continue
            cal = pd.to_numeric(t[cc], errors="coerce").to_numpy(float)
            ref = pd.to_numeric(t[rc], errors="coerce").to_numpy(float)
            ax = pd.to_numeric(t[ax_col], errors="coerce").to_numpy(float)
            m = (np.isfinite(ax) & np.isfinite(cal) & np.isfinite(ref)
                 & np.isfinite(colour))
            s_own = slope(ax, cal - ref, m)
            s_ps1 = slope(colour, cal - ref, m)
            print(f"{band:<4}{tag:<10}{ax_col[6:]:<7}{int(m.sum()):>5}"
                  f"{s_own:>+22.4f}{s_ps1:>+14.4f}")
            axis_rows.append(dict(band=band, instrument=tag, axis=ax_col,
                                  slope_own_axis=s_own, slope_ps1_axis=s_ps1))
        print()

    print("=== 2. 보정성 집합이 얼마나 다른가 ===")
    print(f"{'기기':<10}{'보정성':>8}{'공유':>8}{'공유 비율':>10}"
          f"{'g−r 5~95%':>14}{'중앙 밝기(r)':>14}")
    print("-" * 64)
    for tag, full, sub in (("MuSCAT3", m3, a2), ("Moravian", mv, b2)):
        cg = pd.to_numeric(full.get("color_g_r"), errors="coerce").to_numpy(float)
        cg = cg[np.isfinite(cg)]
        rr = pd.to_numeric(full.get("ref_r"), errors="coerce").to_numpy(float)
        rr = rr[np.isfinite(rr)]
        lo, hi = (np.percentile(cg, 5), np.percentile(cg, 95)) if cg.size else (np.nan,) * 2
        print(f"{tag:<10}{len(full):>8}{len(sub):>8}{len(sub)/len(full):>10.2f}"
              f"{lo:>+7.2f}~{hi:<6.2f}{np.median(rr) if rr.size else np.nan:>14.2f}")

    out = HERE / "residual_color_slope.json"
    out.write_text(json.dumps({"n_common": int(len(a2)), "slopes": rows,
                              "axes": axis_rows},
                              ensure_ascii=False, indent=2), encoding="utf-8")
    print()
    print(f"썼다: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
