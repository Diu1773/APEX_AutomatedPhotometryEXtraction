"""남은 몫은 흩어짐이 아니라 밝기를 따라 미끄러지는 계통이다.

**오차 예산의 산포 항이 잘못 세워져 있었다** (`ERROR_BUDGET.md` 3 절). 두 기기로
잰 같은 별의 차이를 등급 구간으로 나눠 보면, 차이가 밝은 쪽에서 +0.18 이고 어두운
쪽에서 −0.21 로 **한 방향으로 미끄러진다.** 구간 안의 산포는 0.03 밖에 안 된다.
전부를 한 통에 넣고 MAD 를 재면 이 미끄러짐이 「별마다의 흩어짐」으로 둔갑한다.

    MuSCAT3 − kb26 (i)     10~12 등급  +0.180
                           13~14       +0.052
                           16~17       −0.214      구간 안 MAD 0.024~0.099

## 무엇을 가르나

다섯 가지를 잰다. **거의 다 적합을 쓰지 않는다** — 구간 안에서만 재므로 매끄러운
모형을 맞췄다가 산포가 줄어드는 것을 진짜 무늬로 착각할 여지가 없다. 적합을 쓰는
것은 넷째 하나뿐이고, 거기서는 두 모형의 미지수 개수를 맞춰 견준다.

    1. 밝기 구간별     구간마다 중앙 차이와 MAD. 중앙값이 미끄러지면 계통이고,
                      구간 안 MAD 가 전체 MAD 보다 작으면 그만큼이 계통이었다.

    2. 밝기 × 반지름   두 축을 동시에 자른다. **반지름 무늬가 밝기의 그림자일 수
                      있기 때문이다** — 성단은 가운데가 밝고 바깥이 어두우므로
                      밝기 계통 하나가 반지름 추세로 보일 수 있다.

    3. 기기 하나씩     보정성 별의 잔차 `delta − (zp + ct·색)` 를 밝기로 나눈다.
                      기기 사이 비교가 아니라 **그 기기 혼자의 문제**를 본다.
                      한 기기만 미끄러지면 범인이 그 기기다.

    4. 하늘 항인가     하늘을 과대추정하면 별마다 같은 **양**을 잃으므로 곡선의
                      모양이 정해진다. 그 곡선과 직선을 같은 미지수 개수로 맞춰
                      어느 쪽이 자료에 맞는지 본다.

    5. 색이냐 밝기냐   성단에서는 둘이 얽히므로 한쪽을 구간으로 묶어 고정한 채
                      다른 쪽의 폭을 잰다. 그리고 밝기 구간의 중앙값을 뺐을 때
                      잔차가 얼마나 주는지를 **밝기를 섞은 대조군과 함께** 낸다.

실행:
    python -X utf8 validation/residual_magnitude.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "validation" / "external_muscat3"))

from apex.utils.gaia_quality import gaia_quality_report  # noqa: E402
from compare_magnitudes import join, load_workspace, mad  # noqa: E402

TRIO = {
    "Moravian": Path("E:/APEX_validation/reprocess/M67/result"),
    "MuSCAT3": REPO / "validation/external_muscat3/results",
    "kb26": REPO / "validation/external_kb26/results",
}
BANDS = ("g", "r", "i")
MAG_EDGES = [10.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 20.0]
RAD_EDGES = [0.0, 4.0, 9.0, 1e9]          # 분각
SNR_CUT = 20.0
MIN_CELL = 8


def radius_arcmin(ra, dec):
    """시야 중심(별들의 중앙값)에서의 각거리, 분각."""
    ra0, dec0 = float(np.nanmedian(ra)), float(np.nanmedian(dec))
    dx = (ra - ra0) * np.cos(np.deg2rad(dec0))
    return np.hypot(dx, dec - dec0) * 60.0


def _binned(d, key, edges):
    """구간마다 (중앙값, MAD, 별수). 별이 적은 구간은 건너뛴다."""
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        k = np.isfinite(d) & np.isfinite(key) & (key >= lo) & (key < hi)
        if k.sum() < MIN_CELL:
            rows.append((lo, hi, np.nan, np.nan, int(k.sum())))
            continue
        rows.append((lo, hi, float(np.median(d[k])), mad(d[k]), int(k.sum())))
    return rows


def pairs_table(ws: dict) -> list[dict]:
    names = list(ws)
    out: list[dict] = []
    print()
    print("=== 1. 두 기기의 차이를 밝기 구간으로 나누면 ===")
    print("구간 안 MAD 가 전체 MAD 보다 작으면, 그 차이는 흩어짐이 아니라 계통이었다.")
    print()
    hdr = "{:<22}{:<4}{:>7}{:>10}{:>11}{:>7}".format(
        "두 기기", "밴드", "N", "전체 MAD", "구간안 MAD", "줄어든 배")
    hdr += "   " + " ".join(f"{lo:g}~{hi:g}".rjust(8)
                            for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]))
    print(hdr)
    print("-" * len(hdr))

    for i, a in enumerate(names):
        for b in names[i + 1:]:
            j = join(ws[a], ws[b])
            num = lambda c: pd.to_numeric(j[c], errors="coerce").to_numpy(float)  # noqa: E731
            for band in BANDS:
                ca, cb = f"mag_cal_{band}_A", f"mag_cal_{band}_B"
                if ca not in j.columns or cb not in j.columns:
                    continue
                d, mag = num(ca) - num(cb), num(ca)
                base = np.isfinite(d) & np.isfinite(mag)
                if base.sum() < 40:
                    continue
                rows = _binned(d, mag, MAG_EDGES)
                mads = [r[3] for r in rows if np.isfinite(r[3])]
                if not mads:
                    continue
                pooled, within = mad(d[base]), float(np.median(mads))
                cells = " ".join(
                    (f"{r[2]:+.3f}".rjust(8) if np.isfinite(r[2]) else "·".rjust(8))
                    for r in rows)
                print("{:<22}{:<4}{:>7}{:>10.4f}{:>11.4f}{:>7.1f}   {}".format(
                    f"{a} vs {b}", band, int(base.sum()), pooled, within,
                    pooled / within if within else np.nan, cells))
                out.append(dict(kind="magnitude", pair=f"{a} vs {b}", band=band,
                                n=int(base.sum()), mad_pooled=pooled,
                                mad_within_bin=within,
                                bins=[dict(lo=r[0], hi=r[1], median=r[2],
                                           mad=r[3], n=r[4]) for r in rows]))
    return out


def two_way(ws: dict) -> list[dict]:
    """밝기와 반지름을 동시에 자른다 — 어느 쪽이 진짜 축인가."""
    names = list(ws)
    out: list[dict] = []
    print()
    print("=== 2. 밝기와 반지름을 동시에 자르면 ===")
    print("반지름 무늬가 밝기의 그림자인지 본다. 같은 밝기 구간 안에서도 반지름을 따라")
    print("중앙값이 움직이면 반지름이 따로 사는 축이고, 안 움직이면 밝기 하나였다.")
    print()
    hdr = "{:<22}{:<4}{:>10}".format("두 기기", "밴드", "밝기 구간")
    hdr += "".join(f"{f'{lo:g}~{hi:g}′':>18}" for lo, hi in
                   zip(RAD_EDGES[:-1], RAD_EDGES[1:-1]))
    hdr += f"{'9′ 이상':>18}"
    print(hdr)
    print("-" * len(hdr))

    for i, a in enumerate(names):
        for b in names[i + 1:]:
            j = join(ws[a], ws[b])
            num = lambda c: pd.to_numeric(j[c], errors="coerce").to_numpy(float)  # noqa: E731
            rad = radius_arcmin(num("ra_deg_A"), num("dec_deg_A"))
            for band in BANDS:
                ca, cb = f"mag_cal_{band}_A", f"mag_cal_{band}_B"
                if ca not in j.columns or cb not in j.columns:
                    continue
                d, mag = num(ca) - num(cb), num(ca)
                if (np.isfinite(d) & np.isfinite(mag)).sum() < 80:
                    continue
                printed = False
                for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]):
                    inmag = np.isfinite(d) & (mag >= lo) & (mag < hi)
                    if inmag.sum() < 3 * MIN_CELL:
                        continue
                    cells, vals = [], []
                    for rlo, rhi in zip(RAD_EDGES[:-1], RAD_EDGES[1:]):
                        k = inmag & (rad >= rlo) & (rad < rhi)
                        if k.sum() < MIN_CELL:
                            cells.append("·".rjust(18))
                            continue
                        v = float(np.median(d[k]))
                        vals.append(v)
                        cells.append(f"{v:+.3f}({int(k.sum())})".rjust(18))
                    lab = "{:<22}{:<4}".format(f"{a} vs {b}" if not printed else "",
                                               band if not printed else "")
                    printed = True
                    print(lab + "{:>10}".format(f"{lo:g}~{hi:g}") + "".join(cells))
                    out.append(dict(kind="magnitude_x_radius", pair=f"{a} vs {b}",
                                    band=band, mag_lo=lo, mag_hi=hi,
                                    radial_ptp=(float(max(vals) - min(vals))
                                                if len(vals) >= 2 else np.nan)))
                if printed:
                    print()
    ptp = np.array([r["radial_ptp"] for r in out if np.isfinite(r.get("radial_ptp", np.nan))])
    if ptp.size:
        print(f"  같은 밝기 구간 안에서 반지름이 만든 폭: 중앙 {np.median(ptp):.4f} 등급 "
              f"· 최대 {ptp.max():.4f}")
    return out


def per_instrument() -> list[dict]:
    """기기 하나씩 — 보정성 별의 잔차가 자기 밝기를 따라 미끄러지나.

    두 기기를 견주면 **누가 미끄러지는지** 알 수 없다. 한 기기 안에서 잔차
    `delta − (zp + ct·색)` 를 밝기로 나누면 그 기기 혼자의 문제가 드러난다.
    """
    out: list[dict] = []
    print()
    print("=== 3. 기기 하나씩 — 누가 미끄러지나 ===")
    print("보정성 별의 잔차 = delta − (zp + ct·색). 밝기 축은 셋이 공유하는 기준 등급이다.")
    print()
    hdr = "{:<12}{:<4}{:>7}{:>10}".format("기기", "밴드", "N", "잔차 MAD")
    hdr += "   " + " ".join(f"{lo:g}~{hi:g}".rjust(8)
                            for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]))
    hdr += "{:>9}".format("폭")
    print(hdr)
    print("-" * len(hdr))

    for label, rd in TRIO.items():
        zp_dir = Path(rd) / "cmd_zeropoint"
        cal_p, co_p = (zp_dir / "gaia_sdss_calibrator_by_ID.csv",
                       zp_dir / "zp_fit_coefficients.csv")
        if not (cal_p.exists() and co_p.exists()):
            continue
        cal, co = pd.read_csv(cal_p), pd.read_csv(co_p).set_index("filter")
        try:
            qual = np.asarray(gaia_quality_report(cal, cstar_nsigma=None)[0], bool)
        except Exception:  # noqa: BLE001
            qual = np.ones(len(cal), bool)
        for band in [str(b) for b in co.index]:
            colour = str(co.loc[band, "color_col"])
            dc, cc, sc = f"delta_{band}", f"color_{colour}", f"snr_{band}"
            if dc not in cal.columns or cc not in cal.columns:
                continue
            d = pd.to_numeric(cal[dc], errors="coerce").to_numpy(float)
            c = pd.to_numeric(cal[cc], errors="coerce").to_numpy(float)
            s = (pd.to_numeric(cal[sc], errors="coerce").to_numpy(float)
                 if sc in cal.columns else np.full(len(cal), np.inf))
            # 밝기 축은 **기준 등급** `ref_{band}` 다 — Gaia 를 표준계로 옮긴
            # 값이라 세 기기가 똑같은 축을 쓴다. 자기 등급을 축으로 쓰면 기기마다
            # 축이 달라져 견줄 수 없다. 이 기준 자체에 어두운 쪽 치우침이 있어도
            # 세 기기에 똑같이 실리므로, **기기 사이의 차이**는 그대로 남는다.
            rc = f"ref_{band}"
            if rc not in cal.columns:
                continue
            mag = pd.to_numeric(cal[rc], errors="coerce").to_numpy(float)
            resid = d - (float(co.loc[band, "zp"]) + float(co.loc[band, "ct"]) * c)
            m = (np.isfinite(resid) & np.isfinite(mag) & np.isfinite(s)
                 & (s >= SNR_CUT) & qual)
            if m.sum() < 40:
                continue
            r = resid[m] - np.median(resid[m])
            rows = _binned(r, mag[m], MAG_EDGES)
            meds = [x[2] for x in rows if np.isfinite(x[2])]
            span = float(max(meds) - min(meds)) if len(meds) >= 2 else np.nan
            cells = " ".join(
                (f"{x[2]:+.3f}".rjust(8) if np.isfinite(x[2]) else "·".rjust(8))
                for x in rows)
            print("{:<12}{:<4}{:>7}{:>10.4f}   {}{:>9.4f}".format(
                label, band, int(m.sum()), mad(r), cells, span))
            out.append(dict(kind="per_instrument", instrument=label, band=band,
                            n=int(m.sum()), resid_mad=mad(r), span=span,
                            bins=[dict(lo=x[0], hi=x[1], median=x[2],
                                       mad=x[3], n=x[4]) for x in rows]))
    print()
    print("  「폭」은 가장 높은 구간과 가장 낮은 구간의 차이다. 한 기기만 크면 범인이 그 기기다.")
    return out


def flux_offset_model() -> list[dict]:
    """미끄러짐이 **일정한 밝기 값을 빼먹은 것**으로 설명되나.

    방향이 후보를 좁힌다. kb26 은 밝은 별을 더 밝게, 어두운 별을 더 어둡게 낸다.

        혼잡·블렌딩    어두운 별이 이웃과 겹쳐 **더 밝아진다** → 방향이 반대다.
        검출 편향      잡음이 어두운 별을 밝은 쪽으로 흩뜨린다 → 역시 반대다.
        조리개 손실    별마다 같은 **비율**을 잃으므로 밝기에 안 움직인다.
        **하늘 과대추정**  별마다 같은 **양**을 잃는다. 밝은 별에는 티가 안 나고
                       어두운 별에서 커진다 → **방향과 모양이 둘 다 맞는다.**

    그래서 하늘을 얼마 더 뺐다고 하면 곡선이 맞는지 본다. 별의 밝기 값을 f,
    더 뺀 양을 Δf 라 하면

        잰 등급 − 참 등급 = −2.5 log10(1 − Δf/f)

    이고, Δf 와 같은 밝기 값을 갖는 등급을 m_c 라 하면 Δf/f = 10^(0.4(m − m_c)) 다.
    영점이 상수를 흡수하므로 맞출 것은 **둘뿐이다** — m_c 와 상수.

        잔차(m) = 2.5 log10(1 − 10^(0.4(m − m_c))) + c

    **직선과 견준다.** 두 모형 다 두 개짜리라 자유도가 같으므로, 남는 잔차를
    그대로 견주어도 공정하다.
    """
    from scipy.optimize import least_squares

    rows = [r for r in _LAST_PER_INSTRUMENT if r.get("kind") == "per_instrument"]
    out: list[dict] = []
    print()
    print("=== 4. 하늘을 더 뺀 것으로 설명되나 ===")
    print("잔차(m) = 2.5 log10(1 − 10^(0.4(m − m_c))) + c 를 맞춘다. m_c 는 더 뺀 양과")
    print("같은 밝기의 등급이다 — 작을수록 많이 뺐다는 뜻이다.")
    print()
    print("{:<12}{:<4}{:>6}{:>9}{:>11}{:>11}{:>9}".format(
        "기기", "밴드", "점수", "m_c", "하늘모형 RMS", "직선 RMS", "이긴 쪽"))
    print("-" * 62)
    for r in rows:
        pts = [(0.5 * (b["lo"] + b["hi"]), b["median"]) for b in r["bins"]
               if np.isfinite(b.get("median", np.nan))]
        if len(pts) < 4:
            continue
        m = np.array([p[0] for p in pts], float)
        y = np.array([p[1] for p in pts], float)

        def res(par):
            mc, c = par
            x = 10.0 ** (0.4 * (m - mc))
            # 1 − x <= 0 이면 그 등급의 별은 하늘에 묻힌다. 큰 벌점으로 막는다.
            safe = np.clip(1.0 - x, 1e-6, None)
            return 2.5 * np.log10(safe) + c - y

        best = None
        for mc0 in (m.max() + 0.5, m.max() + 1.5, m.max() + 3.0):
            try:
                f = least_squares(res, [mc0, 0.0],
                                  bounds=([m.max() + 0.01, -5.0], [m.max() + 12.0, 5.0]))
            except Exception:  # noqa: BLE001
                continue
            if best is None or f.cost < best.cost:
                best = f
        if best is None:
            continue
        rms_sky = float(np.sqrt(np.mean(best.fun ** 2)))
        lin = np.polyfit(m, y, 1)
        rms_lin = float(np.sqrt(np.mean((np.polyval(lin, m) - y) ** 2)))
        winner = "하늘" if rms_sky < rms_lin else "직선"
        print("{:<12}{:<4}{:>6}{:>9.2f}{:>11.4f}{:>11.4f}{:>9}".format(
            r["instrument"], r["band"], len(m), float(best.x[0]),
            rms_sky, rms_lin, winner))
        out.append(dict(kind="flux_offset", instrument=r["instrument"],
                        band=r["band"], n_points=len(m), m_c=float(best.x[0]),
                        offset_const=float(best.x[1]), rms_sky_model=rms_sky,
                        rms_linear=rms_lin, winner=winner))
    print()
    print("  두 모형 다 미지수가 둘이라 자유도가 같다. RMS 가 작은 쪽이 자료에 맞는 모형이다.")
    return out


def colour_vs_magnitude() -> list[dict]:
    """색이냐 밝기냐 — 둘을 서로 붙잡아 놓고 가른다.

    **성단에서는 이 둘이 얽힌다.** 밝은 별이 곧 특정한 색이므로, 색항이 틀린
    것 하나가 밝기 추세로 보일 수 있고 그 반대도 된다. 그래서 한쪽을 구간으로
    묶어 고정한 채 다른 쪽의 폭을 잰다.

    그리고 **밝기 구간의 중앙값을 빼면 잔차가 얼마나 주는지**를 같이 낸다.
    구간 중앙값을 빼는 것은 적합이 아니지만 그래도 자유도를 쓰므로, 밝기를 섞은
    대조군을 20 번 돌려 거기서 줄어드는 몫을 뺀다.
    """
    out: list[dict] = []
    print()
    print("=== 5. 색이냐 밝기냐 — 한쪽을 고정하고 다른 쪽을 잰다 ===")
    print("성단에서는 밝기와 색이 얽힌다. 고정한 뒤에도 남는 쪽이 진짜 축이다.")
    print()
    print("{:<10}{:<4}{:>6}{:>8}{:>9}{:>9}{:>11}{:>11}{:>10}".format(
        "기기", "밴드", "N", "색↔등급", "색 폭", "등급 폭",
        "등급고정 색", "색고정 등급", "밝기 뺀 뒤"))
    print("-" * 78)
    for label, rd in TRIO.items():
        zp_dir = Path(rd) / "cmd_zeropoint"
        cal_p, co_p = (zp_dir / "gaia_sdss_calibrator_by_ID.csv",
                       zp_dir / "zp_fit_coefficients.csv")
        if not (cal_p.exists() and co_p.exists()):
            continue
        cal, co = pd.read_csv(cal_p), pd.read_csv(co_p).set_index("filter")
        try:
            qual = np.asarray(gaia_quality_report(cal, cstar_nsigma=None)[0], bool)
        except Exception:  # noqa: BLE001
            qual = np.ones(len(cal), bool)
        for band in [str(b) for b in co.index]:
            colname = str(co.loc[band, "color_col"])
            cc, dc, rc = f"color_{colname}", f"delta_{band}", f"ref_{band}"
            if dc not in cal.columns or cc not in cal.columns or rc not in cal.columns:
                continue
            d = pd.to_numeric(cal[dc], errors="coerce").to_numpy(float)
            c = pd.to_numeric(cal[cc], errors="coerce").to_numpy(float)
            m0 = pd.to_numeric(cal[rc], errors="coerce").to_numpy(float)
            s = (pd.to_numeric(cal[f"snr_{band}"], errors="coerce").to_numpy(float)
                 if f"snr_{band}" in cal.columns else np.full(len(cal), np.inf))
            r = d - (float(co.loc[band, "zp"]) + float(co.loc[band, "ct"]) * c)
            m = (np.isfinite(r) & np.isfinite(c) & np.isfinite(m0)
                 & np.isfinite(s) & (s >= SNR_CUT) & qual)
            if m.sum() < 60:
                continue
            r, c, m0 = r[m] - np.median(r[m]), c[m], m0[m]
            rho = float(np.corrcoef(c, m0)[0, 1])
            ce = np.percentile(c, [0, 20, 40, 60, 80, 100])
            me = np.percentile(m0, [0, 20, 40, 60, 80, 100])

            def _span(key, edges):
                v = [np.median(r[(key >= lo) & (key <= hi)])
                     for lo, hi in zip(edges[:-1], edges[1:])
                     if ((key >= lo) & (key <= hi)).sum() >= MIN_CELL]
                return float(max(v) - min(v)) if len(v) >= 2 else np.nan

            def _span_held(primary, pe, other, oe):
                """다른 축을 구간으로 묶어 고정한 채 이 축의 폭을 잰다."""
                got = []
                for lo, hi in zip(oe[:-1], oe[1:]):
                    k = (other >= lo) & (other <= hi)
                    if k.sum() < 25:
                        continue
                    v = [np.median(r[k & (primary >= a) & (primary <= b)])
                         for a, b in zip(pe[:-1], pe[1:])
                         if (k & (primary >= a) & (primary <= b)).sum() >= MIN_CELL]
                    if len(v) >= 2:
                        got.append(max(v) - min(v))
                return float(np.median(got)) if got else np.nan

            def _after_removing(key, edges):
                """구간 중앙값을 뺀 뒤의 MAD."""
                z = r.copy()
                for lo, hi in zip(edges[:-1], edges[1:]):
                    k = (key >= lo) & (key <= hi)
                    if k.sum() >= MIN_CELL:
                        z[k] = z[k] - np.median(r[k])
                return mad(z)

            after = _after_removing(m0, me)
            rng = np.random.default_rng(20260910)
            null = [_after_removing(m0[rng.permutation(len(m0))], me)
                    for _ in range(20)]
            null_med = float(np.median(null))
            print("{:<10}{:<4}{:>6}{:>+8.2f}{:>9.4f}{:>9.4f}{:>11.4f}{:>11.4f}"
                  "{:>10.4f}".format(
                      label, band, int(m.sum()), rho, _span(c, ce), _span(m0, me),
                      _span_held(c, ce, m0, me), _span_held(m0, me, c, ce), after))
            out.append(dict(kind="colour_vs_magnitude", instrument=label, band=band,
                            n=int(m.sum()), rho_colour_mag=rho,
                            span_colour=_span(c, ce), span_mag=_span(m0, me),
                            span_colour_held=_span_held(c, ce, m0, me),
                            span_mag_held=_span_held(m0, me, c, ce),
                            resid_mad=mad(r), resid_mad_after_mag=after,
                            resid_mad_after_shuffled=null_med))
    print()
    print("  「등급고정 색」이 크면 색항 문제이고, 「색고정 등급」이 크면 밝기 문제다.")
    print("  「밝기 뺀 뒤」는 밝기 구간의 중앙값을 뺀 잔차 MAD 다 — 원래 MAD 와 견주면")
    print("  그 기기의 잔차 가운데 얼마가 밝기 계통이었는지 나온다.")
    return out


_LAST_PER_INSTRUMENT: list[dict] = []


def main() -> int:
    ws = {}
    for name, rd in TRIO.items():
        try:
            ws[name] = load_workspace(Path(rd), name)
        except SystemExit:
            print(f"[{name}] 등급 표가 없어 건너뛴다")
    global _LAST_PER_INSTRUMENT
    _LAST_PER_INSTRUMENT = per_instrument()
    rows = pairs_table(ws) + two_way(ws) + _LAST_PER_INSTRUMENT + flux_offset_model() + colour_vs_magnitude()
    p = REPO / "validation/residual_magnitude.json"
    p.write_text(json.dumps(rows, ensure_ascii=False, indent=2, default=str),
                 encoding="utf-8")
    print()
    print(f"썼다: {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
