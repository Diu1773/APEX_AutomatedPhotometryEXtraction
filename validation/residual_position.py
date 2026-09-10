"""남은 0.028 등급이 화면 위치에서 오는가 — 잔차를 검출기 좌표에 펼친다.

**오차 예산이 남긴 하나가 이것이다** (`ERROR_BUDGET.md` 6 절). 계통은 σ_cal 이
설명하고(비 1.2) 산포는 실제 재현성이 설명하는데(비 1.3), **0.028 등급이 남는다.**
혼잡은 아니다 — 이웃이 12 초각 넘게 떨어진 별에서도 바닥이 안 내려간다.

남은 후보 가운데 가장 먼저 잴 것이 **위치 의존**이다. 플랫 잔차든 위치에 따라
조리개가 담는 몫이 달라지는 것이든, 둘 다 검출기 위에서 매끄러운 무늬로 나타나고
**별마다 붙박이라 프레임을 쌓아도 안 줄어든다** — 우리가 쫓는 항의 성질과 같다.

## 재는 법

한 워크스페이스 안에서 끝난다. 보정성 별의 잔차를 영점 모형에서 뽑고

    잔차 = delta − (zp + ct·색)

그것을 검출기 좌표에 펼친다. 위치가 원인이면 **매끄러운 무늬**가 보이고, 아니면
잔차가 위치와 무관하게 흩어진다. 두 가지로 가른다.

    구역별 중앙값   화면을 3×3 으로 나눠 구역마다 잔차의 중앙값.
                    구역 사이 폭이 별 하나의 산포보다 크면 무늬가 있는 것이다.
    반지름 추세     중심에서의 거리에 대한 잔차의 기울기. 플랫·비네팅 잔차는
                    보통 반지름 방향이라 여기서 가장 먼저 걸린다.

**무늬가 있어도 그것이 남은 0.028 을 다 설명한다는 뜻은 아니다.** 설명하는 몫을
등급으로 적고, 남으면 남는다고 적는다.

실행:
    python -X utf8 validation/residual_position.py
"""
from __future__ import annotations

import glob
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO))

from apex.utils.gaia_quality import gaia_quality_report  # noqa: E402

SNR_CUT = 20.0
FIELDS = ("M13", "M3", "M5", "M67", "M37", "NGC6811", "NGC457")
ROOTS = ("E:/APEX_validation/reprocess/*/result", "E:/observed_Analysis/*/*/result")
EXTRA = {"MuSCAT3": REPO / "validation/external_muscat3/results",
         "kb26": REPO / "validation/external_kb26/results"}


def field_of(path: str) -> str:
    for name in FIELDS:
        if re.search(rf"[\\/]{name}[\\/]", path):
            return name
    return "?"


def mad(v) -> float:
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if v.size == 0:
        return float("nan")
    return float(1.4826 * np.median(np.abs(v - np.median(v))))


def rows_for(label: str, result_dir: Path) -> list[dict]:
    zp_dir = Path(result_dir) / "cmd_zeropoint"
    cal_p = zp_dir / "gaia_sdss_calibrator_by_ID.csv"
    co_p = zp_dir / "zp_fit_coefficients.csv"
    if not (cal_p.exists() and co_p.exists()):
        return []
    try:
        cal = pd.read_csv(cal_p)
        co = pd.read_csv(co_p).set_index("filter")
    except Exception:  # noqa: BLE001
        return []
    if "x_pix" not in cal.columns or "y_pix" not in cal.columns:
        return []
    try:
        qual = np.asarray(gaia_quality_report(cal, cstar_nsigma=None)[0], bool)
    except Exception:  # noqa: BLE001
        qual = np.ones(len(cal), bool)

    xp = pd.to_numeric(cal["x_pix"], errors="coerce").to_numpy(float)
    yp = pd.to_numeric(cal["y_pix"], errors="coerce").to_numpy(float)

    out = []
    for band in co.index:
        band = str(band)
        colour = str(co.loc[band, "color_col"])
        dc, cc, ec, sc = (f"delta_{band}", f"color_{colour}",
                          f"mag_inst_err_{band}", f"snr_{band}")
        if dc not in cal.columns or cc not in cal.columns:
            continue
        d = pd.to_numeric(cal[dc], errors="coerce").to_numpy(float)
        c = pd.to_numeric(cal[cc], errors="coerce").to_numpy(float)
        s = (pd.to_numeric(cal[sc], errors="coerce").to_numpy(float)
             if sc in cal.columns else np.full(len(cal), np.inf))
        zp = float(co.loc[band, "zp"])
        ct = float(co.loc[band, "ct"])
        resid = d - (zp + ct * c)

        m = (np.isfinite(resid) & np.isfinite(xp) & np.isfinite(yp)
             & np.isfinite(s) & (s >= SNR_CUT) & qual)
        if int(m.sum()) < 60:
            continue
        r, x, y = resid[m], xp[m], yp[m]
        r = r - np.median(r)

        # 구역별 중앙값 — 화면을 3×3 으로.
        xe = np.percentile(x, [0, 100 / 3, 200 / 3, 100])
        ye = np.percentile(y, [0, 100 / 3, 200 / 3, 100])
        cells = []
        for i in range(3):
            for j in range(3):
                k = ((x >= xe[i]) & (x <= xe[i + 1])
                     & (y >= ye[j]) & (y <= ye[j + 1]))
                if k.sum() >= 8:
                    cells.append(float(np.median(r[k])))
        cell_ptp = float(max(cells) - min(cells)) if len(cells) >= 4 else np.nan

        # 반지름 추세 — 중심에서의 거리에 대한 기울기.
        cx, cy = float(np.median(x)), float(np.median(y))
        rad = np.hypot(x - cx, y - cy)
        rad_n = rad / (np.percentile(rad, 95) or 1.0)
        slope = float(np.polyfit(rad_n, r, 1)[0]) if np.ptp(rad_n) > 0 else np.nan

        # 위치로 설명되는 몫 — 매끄러운 2 차 면을 맞춰 잔차가 줄어드는 정도.
        #
        # **면만 맞춰서는 안 된다.** 여섯 개짜리 모형은 순수 잡음에서도 산포를
        # 줄인다("중첩 모형을 잔차 산포로 비교하면 안 된다", TRACK 함정).
        # 그래서 **위치를 섞은 대조군**을 같이 돌린다 — 섞으면 별과 자리의
        # 짝이 끊기므로 거기서 줄어드는 몫은 전부 과적합이다. 진짜 무늬의
        # 크기는 둘의 차이다.
        xs = (x - cx) / (np.ptp(x) or 1.0)
        ys = (y - cy) / (np.ptp(y) or 1.0)

        def _explained(vals, ax, ay):
            A = np.column_stack([np.ones_like(ax), ax, ay,
                                 ax * ax, ay * ay, ax * ay])
            try:
                coef, *_ = np.linalg.lstsq(A, vals, rcond=None)
                model = A @ coef
            except np.linalg.LinAlgError:
                return 0.0, mad(vals)
            b, a = mad(vals), mad(vals - model)
            return float(np.sqrt(max(b ** 2 - a ** 2, 0.0))), a

        explained, after = _explained(r, xs, ys)
        before = mad(r)

        rng = np.random.default_rng(20260910)
        null = []
        for _ in range(20):
            k = rng.permutation(len(r))
            null.append(_explained(r, xs[k], ys[k])[0])
        null_med = float(np.median(null))
        genuine = float(np.sqrt(max(explained ** 2 - null_med ** 2, 0.0)))

        out.append(dict(field=label, band=band, n=int(m.sum()),
                        resid_mad=before, resid_mad_after_surface=after,
                        explained_mag=explained, shuffled_mag=null_med,
                        genuine_mag=genuine, cell_ptp=cell_ptp,
                        radial_slope=slope))
    return out


def main() -> int:
    rows: list[dict] = []
    seen: set[str] = set()
    for pat in ROOTS:
        for rd in sorted(glob.glob(pat)):
            f = field_of(str(Path(rd) / "x"))
            if f == "?" or f in seen:
                continue
            got = rows_for(f, Path(rd))
            if got:
                seen.add(f)
                rows += got
    for name, rd in EXTRA.items():
        rows += rows_for(name, rd)

    d = pd.DataFrame(rows)
    if d.empty:
        print("자료를 못 찾았다.")
        return 1

    print("=== 보정성 잔차가 검출기 위치와 상관되나 ===")
    print("잔차 = delta − (zp + ct·색). 매끄러운 2 차 면을 빼면 얼마나 줄어드나.")
    print()
    print("{:<9}{:<4}{:>6}{:>10}{:>11}{:>11}{:>11}{:>11}".format(
        "시야", "밴드", "별수", "잔차 MAD", "면이 준 몫", "섞었을 때",
        "진짜 무늬", "반지름 기울기"))
    print("-" * 72)
    for _, r in d.sort_values("genuine_mag", ascending=False).iterrows():
        print("{:<9}{:<4}{:>6}{:>10.4f}{:>11.4f}{:>11.4f}{:>11.4f}{:>+11.4f}".format(
            r["field"], r["band"], r["n"], r["resid_mad"],
            r["explained_mag"], r["shuffled_mag"], r["genuine_mag"],
            r["radial_slope"]))

    g = d["genuine_mag"].to_numpy(float)
    g = g[np.isfinite(g)]
    n_zero = int((g <= 1e-6).sum())
    print()
    print(f"=== 섞은 대조군을 뺀 진짜 위치 무늬 ===")
    print(f"  중앙 {np.median(g):.4f} 등급 · 최대 {g.max():.4f} · "
          f"무늬가 아예 없는 조합 {n_zero}/{len(g)}")
    print("  오차 예산이 남긴 몫은 0.028 등급이다. 이것이 그만 하면 위치가 그 정체이고,")
    print("  훨씬 작으면 위치는 아니다.")

    out = REPO / "validation/residual_position.json"
    out.write_text(json.dumps(d.to_dict("records"), ensure_ascii=False,
                              indent=2, default=str), encoding="utf-8")
    print()
    print(f"썼다: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
