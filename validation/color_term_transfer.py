"""2 차 색항이 쓸모가 있는가 — 자기 자료 밖에서도 맞는지 본다.

**앞선 측정이 남긴 물음이다** (`color_term_weighting.py`). 거기서 밝힌 것은
두 가지였다. 시야마다 ct2 가 부호까지 뒤집히는 것이 적합 방식 탓이 아니라
자료에 있다는 것, 그리고 가중치가 값을 중앙 0.147 만큼 옮긴다는 것.

**그러면 그 값이 무엇을 뜻하는지는 아직 모른다.** 시야마다 다른 것이 필터의
성질이 시야마다 다르게 보이는 것인지, 아니면 애초에 잴 수 없는 것을 재고 있어서
매번 다른 수가 나오는 것인지 갈리지 않았다.

## 가르는 방법

**적합에 안 쓴 별로 시험한다.** 계수를 반쪽으로 맞추고 나머지 반쪽에서 잔차를
잰다. 2 차항이 진짜라면 안 쓴 별에서도 산포를 줄인다. 잡음을 외운 것이라면
자기 자료에서만 줄고 밖에서는 오히려 는다.

세 가지를 나란히 낸다. 셋 다 **안 쓴 별에서 잰 산포**다.

    같은 시야 안       반쪽으로 맞추고 나머지 반쪽에서 잰다 (40 번 섞어 중앙값)
    다른 시야로        시야 X 의 계수를 시야 Y 에 씌운다
    시야를 넘어 공유    모든 시야를 합쳐 하나로 맞추고 각 시야에서 잰다

**영점은 늘 대상 시야에서 다시 잡는다.** 밤도 기기도 다르니 상수까지 옮기면
그건 색항 시험이 아니다. 옮기는 것은 **색에 대한 모양뿐이다.**

## 읽는 법

낸 수는 전부 **1 차만 썼을 때의 산포 − 2 차까지 썼을 때의 산포**(밀리등급)다.

    양수   2 차항이 산포를 줄였다 = 쓸모가 있다
    음수   2 차항이 산포를 늘렸다 = 잡음을 외운 것이다

실행:
    python -X utf8 validation/color_term_transfer.py
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

ROOTS = ("E:/APEX_validation/reprocess/*/result",
         "E:/observed_Analysis/*/*/result")
EXTRA = {"MuSCAT3": Path("validation/external_muscat3/results")}

KIND = {"M13": "구상", "M3": "구상", "M5": "구상",
        "M67": "산개", "NGC6811": "산개", "M37": "산개", "NGC457": "산개"}
BAND_COLOUR = {"g": "g_r", "r": "g_r", "i": "r_i",
               "B": "B_V", "V": "B_V", "R": "V_R", "I": "V_R"}

N_SHUFFLE = 40
RNG = np.random.default_rng(20260908)


def field_of(path: str) -> str:
    for name in KIND:
        if re.search(rf"[\\/]{name}[\\/]", path):
            return name
    return "?"


def fit(x, y, w, degree):
    c, _, _ = robust_weighted_polyfit(x, y, w=w, degree=degree,
                                      clip_sigma=3.0, iters=5, min_n=10)
    return c


def scatter_with(coef, x, y):
    """계수의 색 모양만 씌우고 영점은 그 자리에서 다시 잡은 뒤의 산포(MAD)."""
    if coef is None or not np.all(np.isfinite(coef)):
        return float("nan")
    shape = np.polyval(np.append(coef[:-1], 0.0), x)   # 상수항을 뺀 모양
    r = y - shape
    r = r - np.median(r)
    return float(1.4826 * np.median(np.abs(r)))


def load(field: str, result_dir: Path) -> dict:
    p = result_dir / "cmd_zeropoint" / "gaia_sdss_calibrator_by_ID.csv"
    if not p.exists():
        return {}
    try:
        d = pd.read_csv(p)
    except Exception:  # noqa: BLE001
        return {}
    out = {}
    for band, colour in BAND_COLOUR.items():
        dc, cc, ec = f"delta_{band}", f"color_{colour}", f"mag_inst_err_{band}"
        if dc not in d.columns or cc not in d.columns or ec not in d.columns:
            continue
        x = pd.to_numeric(d[cc], errors="coerce").to_numpy(float)
        y = pd.to_numeric(d[dc], errors="coerce").to_numpy(float)
        e = pd.to_numeric(d[ec], errors="coerce").to_numpy(float)
        m = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0)
        if int(m.sum()) < 60:
            continue
        out[(band, colour)] = (x[m], y[m], 1.0 / np.maximum(e[m], 1e-6) ** 2)
    return out


def within_field(x, y, w, use_w: bool):
    """반쪽으로 맞추고 나머지 반쪽에서 잰다. 밀리등급 이득(1차 − 2차)."""
    gains = []
    n = len(x)
    ww = w if use_w else np.ones_like(w)
    for _ in range(N_SHUFFLE):
        idx = RNG.permutation(n)
        a, b = idx[: n // 2], idx[n // 2:]
        c1, c2 = fit(x[a], y[a], ww[a], 1), fit(x[a], y[a], ww[a], 2)
        s1, s2 = scatter_with(c1, x[b], y[b]), scatter_with(c2, x[b], y[b])
        if np.isfinite(s1) and np.isfinite(s2):
            gains.append((s1 - s2) * 1000)
    return (float(np.median(gains)) if gains else float("nan"),
            float(np.percentile(gains, 16)) if gains else float("nan"),
            float(np.percentile(gains, 84)) if gains else float("nan"))


def main() -> int:
    data: dict[str, dict] = {}
    for pat in ROOTS:
        for rd in sorted(glob.glob(pat)):
            f = field_of(str(Path(rd) / "x"))
            if f == "?" or f in data:
                continue
            got = load(f, Path(rd))
            if got:
                data[f] = got
    for name, rd in EXTRA.items():
        got = load(name, rd)
        if got:
            data[name] = got

    if not data:
        print("자료를 못 찾았다.")
        return 1

    keys = sorted({k for v in data.values() for k in v})
    summary: list[dict] = []

    print("=== 1. 같은 시야 안 — 반쪽으로 맞추고 나머지 반쪽에서 잰다 ===")
    print("2 차항이 준 이득 (밀리등급, 양수면 산포가 줄었다). 괄호는 40 번 섞은 16~84 %.")
    print()
    print("{:<9}{:<4}{:<6}{:>6}{:>22}{:>22}".format(
        "시야", "밴드", "색", "별수", "무가중 적합", "가중 적합(APEX)"))
    print("-" * 69)
    for band, colour in keys:
        for f, v in data.items():
            if (band, colour) not in v:
                continue
            x, y, w = v[(band, colour)]
            gu, lu, hu = within_field(x, y, w, False)
            gw, lw, hw = within_field(x, y, w, True)
            print("{:<9}{:<4}{:<6}{:>6}{:>+13.2f} ({:+.1f}~{:+.1f}){:>+13.2f} ({:+.1f}~{:+.1f})"
                  .format(f, band, colour, len(x), gu, lu, hu, gw, lw, hw))
            summary.append(dict(field=f, band=band, colour=colour, n=len(x),
                                within_unweighted=gu, within_weighted=gw,
                                within_unw_lo=lu, within_unw_hi=hu))

    print()
    print("=== 2. 다른 시야로 — 시야 X 의 계수를 시야 Y 에 씌운다 ===")
    print("Y 에서 잰 이득 (밀리등급). 대각선은 자기 자신이라 비교용으로만 본다.")
    cross: list[dict] = []
    for band, colour in keys:
        fields = [f for f, v in data.items() if (band, colour) in v]
        if len(fields) < 2:
            continue
        coefs = {}
        for f in fields:
            x, y, w = data[f][(band, colour)]
            coefs[f] = (fit(x, y, np.ones_like(w), 1), fit(x, y, np.ones_like(w), 2))
        print()
        print(f"--- {band} 밴드 ({colour}) ---")
        print("  줄 = 계수를 맞춘 시야 · 칸 = 그 계수를 씌워서 재 본 시야")
        print("{:<12}".format("") + "".join(f"{f:>10}" for f in fields))
        for src in fields:
            c1, c2 = coefs[src]
            cells = []
            for dst in fields:
                x, y, _ = data[dst][(band, colour)]
                g = (scatter_with(c1, x, y) - scatter_with(c2, x, y)) * 1000
                cells.append(g)
                if src != dst:
                    cross.append(dict(band=band, colour=colour, src=src,
                                      dst=dst, gain_mmag=float(g)))
            print("{:<12}".format(src) + "".join(f"{c:>+10.2f}" for c in cells))

    print()
    print("=== 3. 시야를 넘어 하나로 — 다 합쳐 맞춘 계수를 각 시야에서 잰다 ===")
    print("{:<5}{:<6}{:>6}{:>10}{:>34}".format("밴드", "색", "시야", "공유 ct2", "시야마다의 이득(밀리등급)"))
    print("-" * 61)
    shared: list[dict] = []
    for band, colour in keys:
        fields = [f for f, v in data.items() if (band, colour) in v]
        if len(fields) < 2:
            continue
        # 시야마다 영점이 다르므로 각자 중앙값을 뺀 뒤 합친다.
        xs, ys = [], []
        for f in fields:
            x, y, _ = data[f][(band, colour)]
            xs.append(x)
            ys.append(y - np.median(y))
        X, Y = np.concatenate(xs), np.concatenate(ys)
        W = np.ones_like(X)
        c1, c2 = fit(X, Y, W, 1), fit(X, Y, W, 2)
        gains = []
        for f in fields:
            x, y, _ = data[f][(band, colour)]
            gains.append((scatter_with(c1, x, y) - scatter_with(c2, x, y)) * 1000)
        print("{:<5}{:<6}{:>6}{:>+10.3f}   ".format(band, colour, len(fields), float(c2[0]))
              + " ".join(f"{f}{g:+.1f}" for f, g in zip(fields, gains)))
        shared.append(dict(band=band, colour=colour, ct2_shared=float(c2[0]),
                           gains={f: float(g) for f, g in zip(fields, gains)}))

    print()
    print("=== 요약 ===")
    w_un = np.array([s["within_unweighted"] for s in summary], float)
    w_w = np.array([s["within_weighted"] for s in summary], float)
    cg = np.array([c["gain_mmag"] for c in cross], float)
    print(f"  같은 시야 안 · 무가중 : 중앙 {np.nanmedian(w_un):+.2f} mmag · "
          f"양수인 조합 {int((w_un > 0).sum())}/{len(w_un)}")
    print(f"  같은 시야 안 · 가중   : 중앙 {np.nanmedian(w_w):+.2f} mmag · "
          f"양수인 조합 {int((w_w > 0).sum())}/{len(w_w)}")
    print(f"  다른 시야로           : 중앙 {np.nanmedian(cg):+.2f} mmag · "
          f"양수인 짝 {int((cg > 0).sum())}/{len(cg)}")

    out = Path("validation/color_term_transfer.json")
    out.write_text(json.dumps({"within": summary, "cross": cross, "shared": shared},
                              ensure_ascii=False, indent=2), encoding="utf-8")
    print()
    print(f"썼다: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
