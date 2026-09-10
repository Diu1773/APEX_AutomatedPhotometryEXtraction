"""남은 몫이 시야 바깥쪽에서 오는가 — 같은 별을 반지름으로 갈라 본다.

**오차 예산이 남긴 마지막 물음이다** (`ERROR_BUDGET.md` 6 절). 재현성과 위치
무늬를 넣고도 **Moravian 대 kb26 만 0.11~0.13 등급이 남는다.** 같은 kb26 인데
MuSCAT3 와 견주면 0 이다.

두 짝의 눈에 띄는 차이는 **겹치는 시야의 크기**다.

    MuSCAT3   2048 화소 × 0.266 초각 = 9 분각      겹친 별 173
    Moravian  4758 화소 × 0.393 초각 = 31 분각     겹친 별 541
    kb26      3100 화소 × 0.58  초각 = 30 분각

MuSCAT3 는 중심만 본다. Moravian 은 kb26 과 시야를 온전히 겹치므로 **바깥쪽까지
본다.** 위치 무늬가 반지름 방향이었으니(같은 문서 6 절), 바깥쪽이 더 나쁘면 남은
몫도 같은 이야기다.

## 두 가지로 가른다

    반지름 구간별 산포   중심에서의 각거리로 나눠 별마다의 차이의 MAD.
                        바깥으로 갈수록 커지면 위치 이야기다.

    **안쪽만 잘라 보기**   Moravian 짝을 MuSCAT3 가 덮는 반지름 안으로 자른다.
                        MuSCAT3 짝 수준으로 내려가면 두 짝의 차이는 시야 크기가
                        전부이고, 안 내려가면 다른 것이 있다.

둘째가 결정적이다 — 같은 kb26, 같은 반지름 범위로 맞춰 놓고 견주는 것이라
시야 크기 말고는 달라지는 것이 없다.

실행:
    python -X utf8 validation/residual_radius.py
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

from compare_magnitudes import join, load_workspace, mad  # noqa: E402

TRIO = {
    "Moravian": Path("E:/APEX_validation/reprocess/M67/result"),
    "MuSCAT3": REPO / "validation/external_muscat3/results",
    "kb26": REPO / "validation/external_kb26/results",
}
BANDS = ("g", "r", "i")
EDGES = np.array([0.0, 2.0, 4.0, 6.0, 9.0, 15.0, 1e9])   # 분각
MAG_EDGES = [10.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 20.0]


def radius_arcmin(ra, dec):
    """시야 중심(별들의 중앙값)에서의 각거리, 분각."""
    ra0, dec0 = float(np.nanmedian(ra)), float(np.nanmedian(dec))
    dx = (ra - ra0) * np.cos(np.deg2rad(dec0))
    return np.hypot(dx, dec - dec0) * 60.0


def main() -> int:
    ws = {}
    for name, rd in TRIO.items():
        try:
            ws[name] = load_workspace(Path(rd), name)
        except SystemExit:
            print(f"[{name}] 등급 표가 없어 건너뛴다")
    names = list(ws)
    out: list[dict] = []

    print()
    print("=== 별마다의 차이가 시야 바깥으로 갈수록 커지나 (M67) ===")
    print("중심에서의 각거리로 나눈 MAD (등급). 괄호는 별 수.")
    print()
    header = "{:<22}{:<3}".format("두 기기", "밴드")
    header += "".join(f"{f'{lo:g}~{hi:g}′':>13}" for lo, hi in zip(EDGES[:-1], EDGES[1:-1]))
    header += f"{'15′ 이상':>13}"
    print(header)
    print("-" * (25 + 13 * (len(EDGES) - 1)))

    pairs = [(a, b) for i, a in enumerate(names) for b in names[i + 1:]]
    joined = {}
    for a, b in pairs:
        j = join(ws[a], ws[b])
        joined[(a, b)] = j
        num = lambda c: pd.to_numeric(j[c], errors="coerce").to_numpy(float)  # noqa: E731
        rad = radius_arcmin(num("ra_deg_A"), num("dec_deg_A"))
        for band in BANDS:
            ca, cb = f"mag_cal_{band}_A", f"mag_cal_{band}_B"
            if ca not in j.columns or cb not in j.columns:
                continue
            d = num(ca) - num(cb)
            base = np.isfinite(d) & np.isfinite(rad)
            if base.sum() < 40:
                continue
            cells, row = [], dict(pair=f"{a} vs {b}", band=band, bins={})
            for lo, hi in zip(EDGES[:-1], EDGES[1:]):
                m = base & (rad >= lo) & (rad < hi)
                if m.sum() < 8:
                    cells.append("            -")
                    continue
                v = mad(d[m])
                cells.append(f"{v:.4f}({int(m.sum())})".rjust(13))
                row["bins"][f"{lo:g}~{hi:g}"] = dict(n=int(m.sum()), mad=v)
            print("{:<22}{:<3}".format(f"{a} vs {b}", band) + "".join(cells))
            out.append(row)

    # ── 결정적 대조 ────────────────────────────────────────────────────────
    # 세 짝을 **같은 반지름 범위**로 잘라 놓고 견준다. MuSCAT3 가 덮는 데까지만
    # 남기면 시야 크기라는 변인이 사라지므로, 그래도 남는 차이는 시야 크기가
    # 아닌 다른 것이다.
    print()
    print("=== 결정적 대조 — 세 짝을 같은 반지름 안으로 자르면 ===")
    r_cut = float("nan")
    if ("MuSCAT3", "kb26") in joined:
        jm = joined[("MuSCAT3", "kb26")]
        nm = lambda c: pd.to_numeric(jm[c], errors="coerce").to_numpy(float)  # noqa: E731
        r_cut = float(np.nanpercentile(
            radius_arcmin(nm("ra_deg_A"), nm("dec_deg_A")), 95))
        print(f"  MuSCAT3 가 덮는 반지름 (95 % 지점): {r_cut:.1f} 분각")
        print()
        print("{:<22}{:<4}{:>12}{:>8}{:>16}{:>8}".format(
            "두 기기", "밴드", "전체 MAD", "별수", f"{r_cut:.1f}′ 안쪽 MAD", "별수"))
        print("-" * 70)
        for a, b in pairs:
            j = joined[(a, b)]
            nj = lambda c: pd.to_numeric(j[c], errors="coerce").to_numpy(float)  # noqa: E731
            rad = radius_arcmin(nj("ra_deg_A"), nj("dec_deg_A"))
            for band in BANDS:
                ca, cb = f"mag_cal_{band}_A", f"mag_cal_{band}_B"
                if ca not in j.columns or cb not in j.columns:
                    continue
                d = nj(ca) - nj(cb)
                m_all = np.isfinite(d)
                m_in = m_all & (rad <= r_cut)
                if m_all.sum() < 40 or m_in.sum() < 20:
                    continue
                print("{:<22}{:<4}{:>12.4f}{:>8}{:>16.4f}{:>8}".format(
                    f"{a} vs {b}", band, mad(d[m_all]), int(m_all.sum()),
                    mad(d[m_in]), int(m_in.sum())))
                out.append(dict(pair=f"{a} vs {b}", band=band, kind="radius_cut",
                                r_cut_arcmin=r_cut,
                                mad_all=mad(d[m_all]), n_all=int(m_all.sum()),
                                mad_inner=mad(d[m_in]), n_inner=int(m_in.sum())))
        print()
        print("  같은 반지름에서 kb26 짝 둘이 비슷해지면, 두 짝의 차이는 시야 크기가")
        print("  전부다. 그래도 남으면 다른 것이 있다.")

    # ── 잘랐더니 준 것이 자리 때문인가, 밝기 구성 때문인가 ────────────────
    # **안쪽만 남기면 밝기 구성도 같이 바뀐다.** 성단 바깥은 어두운 별의 몫이 커서,
    # 위치를 자른 것이 사실은 밝기를 자른 것일 수 있다. 그래서 **안쪽과 똑같은
    # 밝기 구성을 갖도록 전체에서 다시 뽑은 대조군**을 200 번 만들어 견준다.
    # 이 대조군은 위치를 안 자르므로, 여기까지 내려온 몫은 밝기 구성의 것이고
    # 거기서 더 내려간 몫만 자리의 것이다.
    print()
    print("=== 잘랐더니 준 것이 자리 때문인가, 밝기 구성 때문인가 ===")
    if ("Moravian", "kb26") in joined:
        jv = joined[("Moravian", "kb26")]
        nv = lambda c: pd.to_numeric(jv[c], errors="coerce").to_numpy(float)  # noqa: E731
        r_mv = radius_arcmin(nv("ra_deg_A"), nv("dec_deg_A"))
        rng = np.random.default_rng(20260910)
        print("{:<6}{:>10}{:>12}{:>16}{:>12}{:>12}".format(
            "밴드", "전체 MAD", "안쪽 MAD", "밝기만 맞춘 대조군",
            "밝기의 몫", "자리의 몫"))
        print("-" * 68)
        for band in BANDS:
            ca, cb = f"mag_cal_{band}_A", f"mag_cal_{band}_B"
            if ca not in jv.columns or cb not in jv.columns:
                continue
            d, m0 = nv(ca) - nv(cb), nv(ca)
            allm = np.isfinite(d) & np.isfinite(m0)
            inn = allm & (r_mv <= r_cut)
            if allm.sum() < 80 or inn.sum() < 40:
                continue
            draws = []
            for _ in range(200):
                pick: list[int] = []
                for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]):
                    need = int((inn & (m0 >= lo) & (m0 < hi)).sum())
                    pool = np.flatnonzero(allm & (m0 >= lo) & (m0 < hi))
                    if need and pool.size:
                        pick += list(rng.choice(pool, size=min(need, pool.size),
                                                replace=False))
                if len(pick) >= 40:
                    draws.append(mad(d[np.array(pick)]))
            if not draws:
                continue
            a, b, c = mad(d[allm]), float(np.median(draws)), mad(d[inn])
            mag_part = float(np.sqrt(max(a ** 2 - b ** 2, 0.0)))
            pos_part = float(np.sqrt(max(b ** 2 - c ** 2, 0.0)))
            print("{:<6}{:>10.4f}{:>12.4f}{:>16.4f}{:>12.4f}{:>12.4f}".format(
                band, a, c, b, mag_part, pos_part))
            out.append(dict(pair="Moravian vs kb26", band=band,
                            kind="cut_decomposition", mad_all=a, mad_inner=c,
                            mad_mag_matched=b, magnitude_part=mag_part,
                            position_part=pos_part))
        print()
        print("  대조군은 안쪽과 밝기 구성만 같고 위치는 안 잘랐다. 거기까지 내려온")
        print("  몫이 밝기 구성의 것이고, 거기서 더 내려간 몫만 자리의 것이다.")

    p = REPO / "validation/residual_radius.json"
    p.write_text(json.dumps(out, ensure_ascii=False, indent=2, default=str),
                 encoding="utf-8")
    print()
    print(f"썼다: {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
