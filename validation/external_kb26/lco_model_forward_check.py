"""LCO 의 응답 모형을 앞으로 밀어서 APEX 가 잰 것과 대 본다 (2026-09-14).

## 왜 또 재나

앞선 `lco_nonlinearity_correction.py` 는 **맞춤**이다. (z, k) 를 훑어 치우침이 가장
작아지는 값을 골랐으니 치우침이 작아지는 것은 당연하고, 그것만으로는 모형이 맞는지
알 수 없다. **맞춤에 안 쓴 것을 예측해서 맞아야 모형이다.**

여기서 예측하는 것 셋은 전부 (z, k) 를 고를 때 안 본 것이다.

    1 조리개를 바꾸면 치우침이 얼마가 되나   맞춤은 2″ 조리개 하나만 썼다
    2 하늘이 밝으면 치우침이 얼마가 되나     맞춤은 kb26 의 한 밤만 썼다
    3 APEX 가 밴드마다 낸 기울기             맞춤은 r·i 만, 그것도 내 측광으로 했다

## 모형이 무엇인가

**물리 모형이 아니라 경험적 맞춤식이다.** 보고서도 그렇게 적었다 —
*"Our description of the nonlinearity is empirical; we cannot explain the flaw in
the camera hardware that causes the apparent loss of charge."*

    s = (x^k + z^k)^(1/k) − z        x 참값 · s 카메라가 읽은 값

**「z 를 직각으로 더한다」를 일반화한 것**이다. k=2 면 `s = √(x²+z²) − z` 이고, k 를
바꾸면 휘는 정도가 달라진다. 성질 셋만 알면 된다.

    k = 1      s = x 로 딱 떨어진다 — 아무것도 안 고치는 것과 같다
    x → 0      s → 0 (빛이 없으면 읽는 값도 0)
    x ≫ z      s → x − z 에 가까워진다 (밝은 쪽은 기울기 1 의 직선)

그래서 **밝은 쪽 직선을 뒤로 늘인 것에 견주면 어두운 쪽이 위로 뜬다.** 이것이 보고서
그림 4 의 모양이고, 어두운 화소를 실제보다 높게 읽는다는 뜻이다. 하늘이 그 어두운
화소이므로 **하늘이 과대추정되고, 그만큼 과대차감돼 어두운 별이 더 많이 잃는다.**

실행:
    python -X utf8 validation/external_kb26/lco_model_forward_check.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from scipy.special import erf

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))

OUT = REPO / "validation/external_kb26/lco_model_forward_check.json"
BANZAI = {"kb26": Path("E:/APEX_validation_output/external_kb26/raw/banzai"),
          "kb27": Path("E:/APEX_validation_output/external_kb27/raw/banzai")}
RESULTS = {k: REPO / f"validation/external_{k}/results/cmd_zeropoint"
           for k in BANZAI}
BAND_FOLDER = {"B": "B", "V": "V", "i": "ip", "r": "rp"}

#: r 밴드 맞춤에서 나온 값. **이 둘만 쓰고 나머지는 전부 예측이다.**
Z_FIT, K_FIT = 800.0, 1.15

#: 앞선 스크립트가 쓴 조리개와 구간. 예측을 같은 자로 재려고 그대로 가져온다.
MAG_EDGES = [10.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]
APER_ARCSEC = 2.0

#: `drift_what_it_is_not.py` 가 BANZAI 조리개 여섯으로 잰 폭. 맞춤에 안 썼다.
MEASURED_APERTURE = {
    "r": dict(enclosed=[0.109, 0.360, 0.607, 0.776, 0.879, 0.941],
              span=[0.2333, 0.2780, 0.3110, 0.3358, 0.3458, 0.3376]),
    "i": dict(enclosed=[0.110, 0.365, 0.617, 0.787, 0.888, 0.947],
              span=[0.2760, 0.2630, 0.2706, 0.2906, 0.2907, 0.2988]),
}


def read_signal(x, z: float = Z_FIT, k: float = K_FIT):
    """참값 `x` 를 카메라가 얼마로 읽나 — 보고서의 정방향 식."""
    return np.power(np.power(x, k) + z ** k, 1.0 / k) - z


def _frame_facts(camera: str, band: str) -> dict:
    """그 밴드의 하늘 배경(e-/화소)·FWHM·화소크기·노출. 여섯 장의 중앙값."""
    folder = BANZAI[camera] / BAND_FOLDER[band]
    sky, fwhm, scale, exp = [], [], [], []
    for path in sorted(folder.glob("*.fits.fz")):
        with fits.open(path, memmap=False) as hdul:
            h = hdul["SCI"].header
        sky.append(float(h.get("L1MEDIAN") or np.nan))
        fwhm.append(float(h.get("L1FWHM") or np.nan))
        scale.append(float(h.get("PIXSCALE") or 0.58))
        exp.append(float(h.get("EXPTIME") or 1.0))
    if not sky:
        return {}
    return dict(sky=float(np.nanmedian(sky)), fwhm=float(np.nanmedian(fwhm)),
                scale=float(np.nanmedian(scale)), exptime=float(np.nanmedian(exp)),
                n=len(sky))


def _flux_scale(camera: str, band: str, facts: dict) -> float:
    """등급을 전자로 옮기는 눈금 — `F = 10^((zp_e − mag)/2.5)`.

    **자료에서 꺼낸다.** 밝은 별은 비선형의 영향이 거의 없으므로, 기준 등급이
    있는 밝은 별들의 Kron 밝기로 눈금을 맞춘다.
    """
    cal = pd.read_csv(RESULTS[camera] / "gaia_sdss_calibrator_by_ID.csv")
    ref = pd.to_numeric(cal.get(f"ref_{band}"), errors="coerce").to_numpy(float)
    if not np.any(np.isfinite(ref)):
        return float("nan")
    folder = BANZAI[camera] / BAND_FOLDER[band]
    path = sorted(folder.glob("*.fits.fz"))[0]
    with fits.open(path, memmap=False) as hdul:
        cat = hdul["CAT"].data
    from astropy.coordinates import SkyCoord
    ca = SkyCoord(cal["ra_deg"].to_numpy(float), cal["dec_deg"].to_numpy(float),
                  unit="deg")
    cb = SkyCoord(np.asarray(cat["ra"], float), np.asarray(cat["dec"], float),
                  unit="deg")
    k, sep, _ = ca.match_to_catalog_sky(cb)
    ok = np.asarray(sep.arcsec) <= 1.5
    flux = np.asarray(cat["flux"], float)[np.asarray(k)]
    peak = np.asarray(cat["peak"], float)[np.asarray(k)]
    # 밝지만 포화하지 않은 별만 — 비선형이 거의 안 닿는 구간이다.
    use = ok & np.isfinite(ref) & (flux > 0) & (peak > 2000) & (peak < 50000)
    if int(use.sum()) < 10:
        use = ok & np.isfinite(ref) & (flux > 0) & (peak < 50000)
    return float(np.median(ref[use] + 2.5 * np.log10(flux[use])))


def magnitude_error(mag: float, facts: dict, zp_e: float, r_ap_px: float,
                    z: float = Z_FIT, k: float = K_FIT) -> float:
    """그 등급의 별이 이 모형 아래서 얼마나 어둡게 측정되나 (등급, 양수=더 어둡게).

    가우스 별을 화소에 깔고, 화소마다 읽은 값을 만들고, 하늘 고리에서 읽은
    하늘을 빼는 — 실제 조리개 측광이 하는 일을 그대로 한다. 하늘 고리는 별빛이
    거의 없으므로 **하늘만 읽은 값**이 된다.
    """
    sigma = facts["fwhm"] / 2.3548 / facts["scale"]
    total = 10.0 ** ((zp_e - mag) / 2.5)
    half = int(np.ceil(max(4 * sigma, r_ap_px + 2)))
    ax = np.arange(-half, half + 1, dtype=float)
    # 화소 하나에 드는 가우스의 몫 — 가장자리를 오차함수로 정확히 적분한다.
    edge = (erf((ax + 0.5) / (sigma * np.sqrt(2)))
            - erf((ax - 0.5) / (sigma * np.sqrt(2)))) / 2.0
    star = total * np.outer(edge, edge)
    yy, xx = np.mgrid[-half:half + 1, -half:half + 1]
    inside = (np.hypot(xx, yy) <= r_ap_px)
    sky = facts["sky"]
    measured = (read_signal(star[inside] + sky, z, k)
                - read_signal(sky, z, k)).sum()
    true = star[inside].sum()
    if measured <= 0 or true <= 0:
        return float("nan")
    return float(-2.5 * np.log10(measured / true))


def predicted_span(facts, zp_e, r_ap_px, edges=MAG_EDGES, **kw) -> float:
    """구간 중앙에서 잰 측정오차의 최대−최소 — 실제로 잰 「폭」과 같은 자."""
    cent = [0.5 * (a + b) for a, b in zip(edges[:-1], edges[1:])]
    v = [magnitude_error(m, facts, zp_e, r_ap_px, **kw) for m in cent]
    v = [x for x in v if np.isfinite(x)]
    return float(max(v) - min(v)) if len(v) >= 2 else float("nan")


def main() -> int:
    rows: list[dict] = []

    print("=== 0. 모형이 화소에 무엇을 하나 ===")
    print(f"z={Z_FIT:.0f} e- · k={K_FIT:.2f} 에서, 참값을 카메라가 얼마로 읽나.")
    print("밝은 쪽 직선(기울기 1)에 견준 뜬 양이 「어두운 화소를 높게 읽는」 정도다.")
    print()
    hi = 20000.0
    offset = hi - read_signal(hi)          # 밝은 쪽 직선의 절편
    print(f"{'참값 x (e-)':>13}{'읽은 값 s':>12}{'s − (x−절편)':>15}")
    for x in (5.0, 20.0, 50.0, 100.0, 300.0, 1000.0, 5000.0, 20000.0):
        s = float(read_signal(x))
        print(f"{x:>13.0f}{s:>12.1f}{s - (x - offset):>15.1f}")
    print(f"  (밝은 쪽 직선의 절편 {offset:.1f} e-)")
    print()

    print("=== 1. 조리개를 바꾸면 — 맞춤에 안 쓴 것 ===")
    print("맞춤은 2″ 조리개 하나만 봤다. 모형이 맞다면 BANZAI 의 조리개 여섯에서")
    print("잰 폭도 맞춰야 한다. 반지름은 담긴 비율에서 되짚는다.")
    print()
    for band in ("r", "i"):
        facts = _frame_facts("kb26", band)
        zp_e = _flux_scale("kb26", band, facts)
        sigma = facts["fwhm"] / 2.3548 / facts["scale"]
        print(f"-- kb26 {band}  하늘 {facts['sky']:.1f} e- · FWHM {facts['fwhm']:.2f}″ "
              f"· σ {sigma:.2f} px · 눈금 zp {zp_e:.3f}")
        print(f"{'담긴 빛':>9}{'반지름(px)':>12}{'잰 폭':>9}{'모형 예측':>11}{'차이':>9}")
        m = MEASURED_APERTURE[band]
        for phi, got in zip(m["enclosed"], m["span"]):
            r_px = sigma * np.sqrt(-2.0 * np.log(1.0 - phi))
            pred = predicted_span(facts, zp_e, r_px)
            print(f"{phi:>9.0%}{r_px:>12.2f}{got:>9.4f}{pred:>11.4f}{pred-got:>9.4f}")
            rows.append(dict(part=1, camera="kb26", band=band, enclosed=phi,
                             r_px=float(r_px), measured=got, predicted=pred))
        print()

    print("=== 2. 하늘이 밝으면 — 맞춤에 안 쓴 것 ===")
    print("APEX 가 밴드마다 낸 기울기(등급당)와 견준다. 밴드마다 등급 범위가 달라")
    print("폭이 아니라 **등급당 기울기**로 봐야 한다.")
    print()
    print(f"{'기기':<6}{'밴드':<5}{'하늘(e-)':>10}{'APEX 기울기':>13}{'모형 예측':>11}{'비':>8}")
    print("-" * 54)
    for camera in ("kb26", "kb27"):
        qc = pd.read_csv(RESULTS[camera] / "zp_qc_summary.csv")
        facts_band = {b: _frame_facts(camera, b) for b in BAND_FOLDER}
        for _, q in qc.iterrows():
            band = str(q["filter"])
            facts = facts_band.get(band) or {}
            if not facts:
                continue
            lever = float(q["drift_faint_mag"]) - float(q["drift_bright_mag"])
            if lever <= 0:
                continue
            got = float(q["bright_to_faint_drift"]) / lever
            zp_e = _flux_scale(camera, band, facts)
            if not np.isfinite(zp_e):
                print(f"{camera:<6}{band:<5}{facts['sky']:>10.1f}"
                      f"{got:>13.4f}{'기준없음':>11}{'':>8}")
                continue
            sigma = facts["fwhm"] / 2.3548 / facts["scale"]
            r_px = APER_ARCSEC / facts["scale"]
            e1 = magnitude_error(float(q["drift_bright_mag"]), facts, zp_e, r_px)
            e2 = magnitude_error(float(q["drift_faint_mag"]), facts, zp_e, r_px)
            # APEX 의 부호는 (기준 − 기기) 잔차라 측정오차와 반대다.
            pred = -(e2 - e1) / lever
            ratio = pred / got if got else float("nan")
            print(f"{camera:<6}{band:<5}{facts['sky']:>10.1f}"
                  f"{got:>13.4f}{pred:>11.4f}{ratio:>8.2f}")
            rows.append(dict(part=2, camera=camera, band=band, sky=facts["sky"],
                             apex_slope=got, predicted_slope=pred,
                             lever=lever, sigma_px=float(sigma)))
    print()
    print("  비가 1 에 가까우면 모형이 그 칸을 맞춘 것이다. z·k 는 kb26 의 r 하나로")
    print("  정했고 나머지 일곱 칸은 전부 예측이며, kb27 은 카메라가 달라서 제 값이")
    print("  따로 있어야 맞다 — 보고서도 카메라마다 심각도가 다르다고 적었다.")

    print()
    print("=== 3. 골짜기 안의 다른 (z, k) 로는 맞나 ===")
    print("맞춤은 z 와 k 가 서로 축퇴돼서 여러 쌍이 비슷하게 평평해진다. 그 쌍들이")
    print("조리개 여섯을 맞추는지 보면 **어느 쌍이 진짜인지**를 가릴 수 있다.")
    print()
    facts = _frame_facts("kb26", "r")
    zp_e = _flux_scale("kb26", "r", facts)
    sigma = facts["fwhm"] / 2.3548 / facts["scale"]
    meas = MEASURED_APERTURE["r"]
    radii = [sigma * np.sqrt(-2.0 * np.log(1.0 - p)) for p in meas["enclosed"]]
    print("  실측  " + " ".join(f"{v:.3f}" for v in meas["span"])
          + f"   1→6 배율 {meas['span'][-1] / meas['span'][0]:.2f}")
    print()
    print(f"{'z(e-)':>7}{'k':>6}{'고친뒤':>9}  예측 여섯"
          + " " * 26 + f"{'배율':>6}{'평균오차':>10}")
    print("-" * 78)
    for z, k, after in ((200, 1.30, 0.0275), (300, 1.15, 0.1013), (450, 1.15, 0.0594),
                        (600, 1.15, 0.0293), (800, 1.15, 0.0241), (100, 1.45, 0.0863),
                        (100, 1.60, 0.0587)):
        pred = [predicted_span(facts, zp_e, r, z=z, k=k) for r in radii]
        err = float(np.mean([a - b for a, b in zip(pred, meas["span"])]))
        print(f"{z:>7}{k:>6.2f}{after:>9.4f}  " + " ".join(f"{v:.3f}" for v in pred)
              + f"{pred[-1] / pred[0]:>6.2f}{err:>10.3f}")
        rows.append(dict(part=3, z=z, k=k, after_correction=after,
                         predicted=pred, ratio=pred[-1] / pred[0], mean_error=err))
    print()
    print("  **가장 잘 고치는 쌍과 가장 잘 예측하는 쌍이 다르다** — 이것이 이 식의")
    print("  꼴이 우리 자료에 딱 맞지는 않는다는 뜻이고, 보고서도 그렇게 적어 두었다.")

    OUT.write_text(json.dumps(dict(z=Z_FIT, k=K_FIT, rows=rows),
                              ensure_ascii=False, indent=2, default=str),
                   encoding="utf-8")
    print()
    print(f"썼다: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
