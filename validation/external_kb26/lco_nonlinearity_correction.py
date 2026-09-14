"""LCO 가 공개한 역보정을 화소에 걸면 밝기 치우침이 사라지는가 (2026-09-14).

## 무엇을 확인하나

kb26 의 밝기 치우침은 **우리가 찾아낸 새 현상이 아니라 Las Cumbres Observatory 가
이름 붙여 문서로 올려 둔 결함**이다 — Harbeck & Volgenau, *"Photometric nonlinearity
in SBIG 6303 images with low light levels"* (2022-03-10, `lco.global/documentation/`,
사본은 이 폴더의 `LCO_Photometric_nonlinearity_SBIG_6303_20220310.pdf`). 보고서의
첫 절이 **"faint sources appear too faint"** 로 시작하고, 심각도가 **하늘 배경에
반비례**하며, 배경 20 e- 언저리에서 영점 폭이 0.4 등급이라고 적혀 있다.

우리 프레임이 정확히 그 조건이다.

    kb26 (lsc, 05-03)   달 위상 0.01 · 고도 −27°   r 하늘 21.9 e-   치우침 0.28~0.40
    kb27 (ogg, 04-15)   달 위상 0.86 · 고도 +66°   r 하늘 89.0 e-   치우침 0.12

**이 스크립트는 그 설명이 맞는지 되돌려서 확인한다.** 보고서가 준 검출기 응답 모형

    s = (x^k + z^k)^(1/k) − z        x 는 참값, s 는 카메라가 읽은 값

의 역함수를 화소마다 걸고

    x = ((s + z)^k − z^k)^(1/k)

다시 측광해서 치우침이 줄어드는지 본다. **줄면 원인이 확정된다.**

## z 와 k 를 어디서 가져오나

**카메라마다 다르고 kb26 값은 공개돼 있지 않다.** 보고서에 있는 것은 실험실에서 잰
kb96(z=238.5, k=1.046)과 과학 자료로 맞춘 kb95(z=297.5, k=1.37)뿐이고, 저자들도
*"different cameras exhibit the photometric nonlinearity with different severity"* 라고
적었다. 그래서 보고서 4 절이 하는 대로 **(z, k) 를 격자로 훑어 영점-등급 기울기를
0 에 가장 가깝게 만드는 값을 찾는다.** 보고서의 두 값은 격자 안에 들어 있으므로
자릿수가 맞는지도 같이 보인다.

단위 주의 — 보고서의 z 는 ADU 이고 **BANZAI 프레임은 이미 전자로 환산돼 있다**
(GAIN 1.0). 그래서 여기의 z 는 전자 단위이며 보고서 값과 직접 비교하면 안 되고,
SBIG 의 변환이득만큼 커야 맞다.

실행:
    python -X utf8 validation/external_kb26/lco_nonlinearity_correction.py
    python -X utf8 validation/external_kb26/lco_nonlinearity_correction.py --band i
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))

from apex.analysis.cmd.zeropoint_runner import ZeropointCalibrationRunner  # noqa: E402
from apex.utils.gaia_quality import gaia_quality_report  # noqa: E402

BANZAI = Path("E:/APEX_validation_output/external_kb26/raw/banzai")
RESULTS = REPO / "validation/external_kb26/results/cmd_zeropoint"
OUT = REPO / "validation/external_kb26/lco_nonlinearity_correction.json"

BAND_FOLDER = {"r": "rp", "i": "ip"}
MAG_EDGES = [10.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 20.0]
MIN_CELL = 8

#: 측광 설정. 반지름은 초각이고, 하늘은 고리의 중앙값이다 — APEX 의 기본과 같다.
APER_ARCSEC = 2.0
SKY_IN_ARCSEC, SKY_OUT_ARCSEC = 6.0, 9.0
STAMP = 24                      # 고리 바깥(9″ ≈ 15.5 px)을 덮는 반쪽 크기

#: 훑을 격자. z=0 은 아무것도 안 고친 것과 같아서 **대조군**이다.
Z_GRID = [0.0, 100.0, 200.0, 300.0, 450.0, 600.0, 800.0, 1100.0, 1500.0]
K_GRID = [1.05, 1.15, 1.30, 1.45, 1.60]


class _Silent:
    def _log(self, *a, **k) -> None:
        return None


def _robust_fit(x: np.ndarray, y: np.ndarray):
    return ZeropointCalibrationRunner._robust_linfit(
        _Silent(), x, y, w=np.ones(len(x)), clip_sigma=3.0, iters=5,
        slope_absmax=0.8, min_n=10)


def undo_nonlinearity(s: np.ndarray, z: float, k: float) -> np.ndarray:
    """카메라가 읽은 값 `s` 에서 참값을 되돌린다 — 보고서 4 쪽의 역함수.

    `z` 가 0 이면 항등이므로 **대조군이 같은 코드 경로를 탄다.** 읽은 값이 음수인
    화소는 0 으로 둔다(읽기 잡음·산탄 잡음 때문에 생기며, 보고서도 같게 한다).
    """
    if z <= 0:
        return s
    v = np.where(np.isfinite(s), np.maximum(s, 0.0), 0.0)
    return np.power(np.maximum(np.power(v + z, k) - z ** k, 0.0), 1.0 / k)


def _stamps(folder: str, ra: np.ndarray, dec: np.ndarray):
    """별마다 프레임마다 조각 하나. 화소 좌표는 조각 안 기준으로 돌려준다.

    조각을 한 번만 잘라 두면 (z, k) 를 바꿀 때마다 3054×2042 를 다시 건드리지
    않아도 되므로, 격자 훑기가 몇 초 단위로 끝난다.
    """
    out = []
    for path in sorted((BANZAI / folder).glob("*.fits.fz")):
        with fits.open(path, memmap=False) as hdul:
            head = hdul["SCI"].header
            data = np.asarray(hdul["SCI"].data, dtype=np.float64)
            exptime = float(head.get("EXPTIME") or 1.0)
            scale = float(head.get("PIXSCALE") or 0.58)
            wcs = WCS(head)
        x, y = wcs.all_world2pix(ra, dec, 0)
        cut, cx, cy, keep = [], [], [], []
        ny, nx = data.shape
        for i, (xi, yi) in enumerate(zip(x, y)):
            if not (np.isfinite(xi) and np.isfinite(yi)):
                continue
            x0, y0 = int(round(xi)) - STAMP, int(round(yi)) - STAMP
            if x0 < 0 or y0 < 0 or x0 + 2 * STAMP + 1 > nx or y0 + 2 * STAMP + 1 > ny:
                continue
            cut.append(data[y0:y0 + 2 * STAMP + 1, x0:x0 + 2 * STAMP + 1])
            cx.append(xi - x0)
            cy.append(yi - y0)
            keep.append(i)
        out.append(dict(name=path.name, exptime=exptime, scale=scale,
                        stamps=np.asarray(cut), cx=np.asarray(cx),
                        cy=np.asarray(cy), keep=np.asarray(keep, int)))
    return out


def _masks(scale: float):
    """조각 안에서 조리개와 하늘 고리에 드는 화소. 조각 크기가 같아 한 번만 만든다."""
    yy, xx = np.mgrid[0:2 * STAMP + 1, 0:2 * STAMP + 1]
    return yy, xx, (APER_ARCSEC / scale, SKY_IN_ARCSEC / scale, SKY_OUT_ARCSEC / scale)


def measure(frames, z: float, k: float, n_star: int) -> np.ndarray:
    """(z, k) 로 되돌린 화소에서 별마다의 기기 등급 중앙값.

    **길이는 보정성 별 표 그대로 돌려준다.** 가장자리에 걸려 조각을 못 자른 별이
    있으면 그만큼 짧아져서, 부르는 쪽의 기준 등급과 줄이 어긋나기 때문이다.
    """
    acc = np.full((len(frames), n_star), np.nan)
    for fi, f in enumerate(frames):
        yy, xx, (r_ap, r_in, r_out) = _masks(f["scale"])
        cube = undo_nonlinearity(f["stamps"], z, k)
        for j, (cx, cy, idx) in enumerate(zip(f["cx"], f["cy"], f["keep"])):
            rr = np.hypot(xx - cx, yy - cy)
            img = cube[j]
            ring = img[(rr >= r_in) & (rr < r_out)]
            ring = ring[np.isfinite(ring)]
            if ring.size < 20:
                continue
            sky = float(np.median(ring))
            inside = (rr <= r_ap)
            flux = float(np.sum(img[inside]) - sky * float(inside.sum()))
            if flux > 0:
                acc[fi, idx] = -2.5 * np.log10(flux / f["exptime"])
    with np.errstate(invalid="ignore"):
        return np.nanmedian(acc, axis=0)


def drift_of(inst: np.ndarray, ref: np.ndarray, colour: np.ndarray):
    """되돌린 등급으로 영점을 다시 맞추고 잔차가 밝기를 따라 얼마나 미끄러지나."""
    delta = ref - inst
    g = np.isfinite(delta) & np.isfinite(colour) & np.isfinite(ref)
    if int(g.sum()) < 60:
        return float("nan"), float("nan"), int(g.sum())
    zp, ct, _, _ = _robust_fit(colour[g], delta[g])
    if not np.isfinite(zp):
        return float("nan"), float("nan"), int(g.sum())
    resid = delta[g] - (zp + ct * colour[g])
    resid -= np.median(resid)
    mag = ref[g]
    meds, cent = [], []
    for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]):
        m = (mag >= lo) & (mag < hi)
        if int(m.sum()) >= MIN_CELL:
            meds.append(float(np.median(resid[m])))
            cent.append(0.5 * (lo + hi))
    if len(meds) < 3:
        return float("nan"), float("nan"), int(g.sum())
    span = float(max(meds) - min(meds))
    slope = float(np.polyfit(cent, meds, 1)[0])
    return span, slope, int(g.sum())


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="LCO 역보정을 걸어 치우침이 주는지 본다")
    ap.add_argument("--band", default="r", choices=sorted(BAND_FOLDER))
    a = ap.parse_args(argv)
    band, folder = a.band, BAND_FOLDER[a.band]

    cal = pd.read_csv(RESULTS / "gaia_sdss_calibrator_by_ID.csv")
    coeff = pd.read_csv(RESULTS / "zp_fit_coefficients.csv").set_index("filter")
    try:
        good = np.asarray(gaia_quality_report(cal, cstar_nsigma=None)[0], bool)
    except (KeyError, ValueError, TypeError) as exc:
        print(f"  [주의] Gaia 품질 표를 못 만들어 별을 안 거르고 전부 쓴다 — {exc}")
        good = np.ones(len(cal), bool)
    colour_col = f"color_{coeff.loc[band, 'color_col']}"
    sub = cal[good
              & np.isfinite(pd.to_numeric(cal[f"ref_{band}"], errors="coerce"))
              & np.isfinite(pd.to_numeric(cal[colour_col], errors="coerce"))].reset_index(drop=True)
    ref = pd.to_numeric(sub[f"ref_{band}"], errors="coerce").to_numpy(float)
    colour = pd.to_numeric(sub[colour_col], errors="coerce").to_numpy(float)

    print(f"보정성 별 {len(sub)} 개 · 밴드 {band} · 색축 {coeff.loc[band, 'color_col']}")
    frames = _stamps(folder, sub["ra_deg"].to_numpy(float),
                     sub["dec_deg"].to_numpy(float))
    print(f"프레임 {len(frames)} 장 · 조각 {2*STAMP+1}×{2*STAMP+1} · "
          f"조리개 {APER_ARCSEC}″ · 하늘 고리 {SKY_IN_ARCSEC}~{SKY_OUT_ARCSEC}″")
    print()

    rows, best = [], None
    print(f"{'z (e-)':>8}{'k':>7}{'별수':>7}{'폭':>9}{'기울기(등급/등급)':>20}")
    print("-" * 52)
    for z in Z_GRID:
        for k in (K_GRID if z > 0 else [1.0]):
            inst = measure(frames, z, k, len(sub))
            span, slope, n = drift_of(inst, ref, colour)
            if not np.isfinite(span):
                continue
            mark = ""
            if best is None or abs(slope) < abs(best["slope"]):
                best = dict(z=z, k=k, span=span, slope=slope, n=n)
                mark = "  ←"
            print(f"{z:>8.0f}{k:>7.2f}{n:>7}{span:>9.4f}{slope:>20.4f}{mark}")
            rows.append(dict(band=band, z=z, k=k, n=n, span=span, slope=slope))
    print()
    zero = next((r for r in rows if r["z"] == 0.0), None)
    if zero and best:
        print(f"안 고친 것(대조군)  폭 {zero['span']:.4f} · 기울기 {zero['slope']:+.4f}")
        print(f"가장 평평한 것      폭 {best['span']:.4f} · 기울기 {best['slope']:+.4f}"
              f"   (z={best['z']:.0f} e-, k={best['k']:.2f})")
        if zero["span"] > 0:
            print(f"폭이 {best['span']/zero['span']*100:.0f} % 로 줄었다")
        print()
        print("  보고서의 값은 kb96 z=238.5 ADU·k=1.046, kb95 z=297.5 ADU·k=1.37 인데")
        print("  둘 다 ADU 라 전자로 옮기면 SBIG 변환이득만큼 커진다. 자릿수만 본다.")

    OUT.write_text(json.dumps(dict(band=band, aper_arcsec=APER_ARCSEC,
                                   sky_annulus=[SKY_IN_ARCSEC, SKY_OUT_ARCSEC],
                                   rows=rows, best=best, uncorrected=zero),
                              ensure_ascii=False, indent=2, default=str),
                   encoding="utf-8")
    print()
    print(f"썼다: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
