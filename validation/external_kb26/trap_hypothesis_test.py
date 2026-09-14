"""전하 트랩 가설을 프레임 안에서 시험한다 — 읽기 방향과 국소 하늘 (2026-09-14).

## 왜 다시 보나

LCO 보고서는 이 현상을 **경험적 맞춤식으로만** 적어 두고 하드웨어 원인은 모른다고
했다(*"we cannot explain the flaw in the camera hardware"*). 그런데 보고서가 관측한
두 가지 — **어두운 별이 고정된 수의 전자를 잃는다**와 **하늘 배경이 밝을수록 덜하다**
— 는 CCD 물리에서 **전하 트랩의 교과서적 지문**이다.

    트랩이 각 전하 덩어리에서 일정 수의 전자를 붙잡는다
      → 어두운 별일수록 잃는 **비율**이 크다              (어두운 별이 더 어둡게)
    배경 전하가 먼저 트랩을 채워 두면 별의 전자를 못 잡는다
      → 하늘이 밝으면 손실이 줄어든다                      (보고서의 배경 의존)

뒤쪽이 옛날 CCD 에서 **fat zero**(작은 배경 전하를 일부러 넣어 트랩을 미리 채우는 것)
라고 불리며 선형성을 되찾는 표준 수법이다. 그러니까 **LCO 가 「원인 불명」이라고 적은
것에 이름이 있다** — 다만 정말 그것인지는 우리 자료로 확인해야 한다.

## 무엇으로 가르나

7 절에서 「전하 전송 손실이 아니다」로 판정했는데, 그 근거는 **사분위로 나눈 폭이
평평하다**는 것이었다. 그 시험을 다시, 더 예민하게 한다 — 폭이 아니라 **등급 기울기
자체가 자리에 따라 달라지는지**를 한 번에 적합하고 부트스트랩으로 오차를 붙인다.

    resid = a + b·M + (c + d·M)·X̂ + (e + f·M)·Ŷ + (g + h·M)·Ŝ

        M  = 기준 등급 − 14        (밝기)
        X̂  = 가로 위치            (읽어 내는 장치가 높은 x 쪽에 있다)
        Ŷ  = 세로 위치            (평행 전송 방향)
        Ŝ  = 별마다의 국소 하늘    (SEP 이 별마다 잰 배경)

**`d`·`f` 가 0 이 아니면 전송 횟수를 타는 것**이라 트랩이 전송 경로에 있다는 뜻이고,
**`h` 가 0 이 아니면 하늘이 트랩을 채우는 것**이라 프레임 **안에서도** 배경 의존이
보인다는 뜻이다. 뒤쪽은 밤을 견주는 것보다 훨씬 강한 증거다 — 같은 프레임, 같은
카메라, 같은 순간이라 다른 변인이 전부 고정돼 있기 때문이다.

실행:
    python -X utf8 validation/external_kb26/trap_hypothesis_test.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))

from apex.analysis.cmd.zeropoint_runner import ZeropointCalibrationRunner  # noqa: E402
from apex.utils.gaia_quality import gaia_quality_report  # noqa: E402

BANZAI = {"kb26": Path("E:/APEX_validation_output/external_kb26/raw/banzai"),
          "kb27": Path("E:/APEX_validation_output/external_kb27/raw/banzai")}
RESULTS = {c: REPO / f"validation/external_{c}/results/cmd_zeropoint" for c in BANZAI}
OUT = REPO / "validation/external_kb26/trap_hypothesis_test.json"

BAND_FOLDER = {"r": "rp", "i": "ip"}
MATCH_ARCSEC = 1.5
MAG_PIVOT = 14.0
N_BOOT = 400
RNG = np.random.default_rng(20260914)


class _Silent:
    def _log(self, *a, **k) -> None:
        return None


def _robust_fit(x, y):
    return ZeropointCalibrationRunner._robust_linfit(
        _Silent(), x, y, w=np.ones(len(x)), clip_sigma=3.0, iters=5,
        slope_absmax=0.8, min_n=10)


def _gather(camera: str, band: str) -> pd.DataFrame | None:
    """보정성 별마다 잔차·등급·검출기 위치·국소 하늘을 한 표로.

    프레임 여섯 장을 다 쓰되 별마다 중앙값으로 접는다 — 한 장만 쓰면 잡음이
    기울기를 덮고, 다 이어 붙이면 같은 별이 여섯 번 세어져 오차가 작아 보인다.
    """
    cal = pd.read_csv(RESULTS[camera] / "gaia_sdss_calibrator_by_ID.csv")
    coeff = pd.read_csv(RESULTS[camera] / "zp_fit_coefficients.csv").set_index("filter")
    if band not in coeff.index:
        return None
    try:
        good = np.asarray(gaia_quality_report(cal, cstar_nsigma=None)[0], bool)
    except (KeyError, ValueError, TypeError) as exc:
        print(f"  [주의] Gaia 품질 표를 못 만들어 전부 쓴다 — {exc}")
        good = np.ones(len(cal), bool)
    colour_col = f"color_{coeff.loc[band, 'color_col']}"
    sub = cal[good
              & np.isfinite(pd.to_numeric(cal[f"ref_{band}"], errors="coerce"))
              & np.isfinite(pd.to_numeric(cal[colour_col], errors="coerce"))
              ].reset_index(drop=True)
    ref = pd.to_numeric(sub[f"ref_{band}"], errors="coerce").to_numpy(float)
    colour = pd.to_numeric(sub[colour_col], errors="coerce").to_numpy(float)
    ca = SkyCoord(sub["ra_deg"].to_numpy(float), sub["dec_deg"].to_numpy(float),
                  unit="deg")

    cols = {name: [] for name in ("inst", "x", "y", "bkg")}
    for path in sorted((BANZAI[camera] / BAND_FOLDER[band]).glob("*.fits.fz")):
        with fits.open(path, memmap=False) as hdul:
            cat = hdul["CAT"].data
            exptime = float(hdul["SCI"].header.get("EXPTIME") or 1.0)
        flag = np.asarray(cat["flag"], int)
        flux = np.asarray(cat["flux"], float)
        keep = (flag == 0) & np.isfinite(flux) & (flux > 0)
        cb = SkyCoord(np.asarray(cat["ra"], float)[keep],
                      np.asarray(cat["dec"], float)[keep], unit="deg")
        k, sep, _ = ca.match_to_catalog_sky(cb)
        ok = np.asarray(sep.arcsec) <= MATCH_ARCSEC
        idx = np.asarray(k)
        for name, src in (("inst", -2.5 * np.log10(flux[keep] / exptime)),
                          ("x", np.asarray(cat["x"], float)[keep]),
                          ("y", np.asarray(cat["y"], float)[keep]),
                          ("bkg", np.asarray(cat["background"], float)[keep])):
            v = np.full(len(sub), np.nan)
            v[ok] = src[idx[ok]]
            cols[name].append(v)

    with np.errstate(invalid="ignore"):
        stack = {n: np.nanmedian(np.vstack(v), axis=0) for n, v in cols.items()}
    delta = ref - stack["inst"]
    g = np.isfinite(delta) & np.isfinite(colour) & np.isfinite(ref)
    if int(g.sum()) < 100:
        return None
    zp, ct, _, _ = _robust_fit(colour[g], delta[g])
    resid = np.full(len(sub), np.nan)
    resid[g] = delta[g] - (zp + ct * colour[g])
    out = pd.DataFrame(dict(resid=resid, mag=ref, x=stack["x"], y=stack["y"],
                            bkg=stack["bkg"]))
    out = out[np.isfinite(out).all(axis=1)].reset_index(drop=True)
    out["resid"] -= out["resid"].median()
    return out


def _design(t: pd.DataFrame, scale: dict) -> tuple[np.ndarray, list[str]]:
    """설계행렬. 위치와 하늘은 중앙 0 · 폭 1 로 맞춰 계수를 서로 견줄 수 있게 한다."""
    M = t["mag"].to_numpy(float) - MAG_PIVOT
    cols, names = [np.ones(len(t)), M], ["상수", "M(등급 기울기)"]
    for key, label in (("x", "X̂(가로)"), ("y", "Ŷ(세로)"), ("bkg", "Ŝ(국소 하늘)")):
        v = (t[key].to_numpy(float) - scale[key][0]) / scale[key][1]
        cols += [v, v * M]
        names += [label, f"M×{label}"]
    return np.vstack(cols).T, names


def fit_with_errors(t: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, list[str]]:
    scale = {k: (float(np.median(t[k])),
                 float(np.percentile(t[k], 84) - np.percentile(t[k], 16)) or 1.0)
             for k in ("x", "y", "bkg")}
    A, names = _design(t, scale)
    y = t["resid"].to_numpy(float)
    beta, *_ = np.linalg.lstsq(A, y, rcond=None)
    boot = np.empty((N_BOOT, A.shape[1]))
    n = len(t)
    for b in range(N_BOOT):
        pick = RNG.integers(0, n, n)
        bb, *_ = np.linalg.lstsq(A[pick], y[pick], rcond=None)
        boot[b] = bb
    return beta, boot.std(axis=0), names


def main() -> int:
    rows: list[dict] = []
    print("=== 전하 트랩 가설 — 등급 기울기가 자리와 하늘을 타는가 ===")
    print("계수는 「그 변수가 1 단위(16~84 백분위 폭) 달라질 때 등급 기울기가")
    print("얼마나 바뀌나」다. 트랩이 전송 경로에 있으면 M×X̂ 또는 M×Ŷ 가 0 이 아니고,")
    print("하늘이 트랩을 채우면 M×Ŝ 가 0 이 아니다.")
    print()
    for camera in ("kb26", "kb27"):
        for band in ("r", "i"):
            t = _gather(camera, band)
            if t is None or len(t) < 100:
                print(f"[{camera} {band}] 별이 모자라 건너뛴다")
                continue
            beta, err, names = fit_with_errors(t)
            print(f"-- {camera} {band}   별 {len(t)} 개 · "
                  f"x {t['x'].min():.0f}~{t['x'].max():.0f} · "
                  f"y {t['y'].min():.0f}~{t['y'].max():.0f} · "
                  f"국소 하늘 {t['bkg'].min():.0f}~{t['bkg'].max():.0f} e-")
            for nm, b, e in zip(names, beta, err):
                sig = abs(b) / e if e > 0 else 0.0
                mark = "  ←" if sig >= 3 and nm.startswith("M×") else ""
                print(f"     {nm:<16}{b:>9.4f} ± {e:.4f}   {sig:>5.1f}σ{mark}")
            rows.append(dict(camera=camera, band=band, n=len(t),
                             names=names, beta=beta.tolist(), err=err.tolist(),
                             bkg_lo=float(t["bkg"].min()),
                             bkg_hi=float(t["bkg"].max())))
            print()
    print("  M×X̂ · M×Ŷ 가 0 이면 전송 횟수를 안 탄다 = 트랩이 경로에 있지 않다.")
    print("  M×Ŝ 가 0 이 아니면 **같은 프레임 안에서도** 하늘이 밝은 자리의 별이")
    print("  덜 미끄러진다는 뜻이고, 그게 fat zero 가 하는 일이다.")

    OUT.write_text(json.dumps(rows, ensure_ascii=False, indent=2, default=str),
                   encoding="utf-8")
    print()
    print(f"썼다: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
