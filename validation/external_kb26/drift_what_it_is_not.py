"""kb26 의 밝기 치우침이 무엇이 아닌지를 다섯 번 재서 좁힌다 (2026-09-14).

앞선 검증에서 APEX 의 전처리도 측광도 아니라는 것까지는 갈랐다
(`validation/ERROR_BUDGET.md` 8 절). 그러면 남는 것은 빛이 검출기에 닿기까지의
과정이거나 검출기 자신인데, 이 둘을 더 좁히려고 **관측소가 붙여 놓은 측광
확장(`CAT`)** 을 쓴다. BANZAI 프레임은 조리개 여섯 개의 밝기와 봉우리 계수를
함께 싣고 있어서, 같은 별을 여러 조리개로 다시 읽을 수 있기 때문이다.

## 다섯 가지를 잰다

    1 조리개별 치우침   조리개를 키우면 커지나 줄면 커지나 무관한가
    2 하늘 몫 분해      1 의 값에서 조리개 넓이에 비례하는 부분만 떼어 낸다
    3 봉우리 계수       보정에 쓰인 별들이 포화·비선형 구간에 걸쳐 있나
    4 한 장씩           여섯 장을 합쳐서 생기는 것인가
    5 성장곡선          별빛이 퍼진 모양이 밝기에 따라 달라지나

**2 번이 이 스크립트의 값어치다.** 앞서 직선과 곡선을 맞춰 「하늘 항이 아니다」고
적었던 판정이 너무 거칠었다(F-336). 하늘을 더 빼는 것은 조리개 넓이에 정비례해서
치우침을 만들므로, 조리개를 바꿔 가며 재면 **하늘 몫과 조리개와 무관한 몫을 실제로
쪼갤 수 있다.** 「아니다」가 아니라 「얼마다」로 답하는 것이 맞다.

실행 (E 드라이브에 BANZAI 프레임이 있어야 한다):
    python -X utf8 validation/external_kb26/drift_what_it_is_not.py
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

BANZAI = Path("E:/APEX_validation_output/external_kb26/raw/banzai")
BANZAI_KB27 = Path("E:/APEX_validation_output/external_kb27/raw/banzai")
RESULTS = REPO / "validation/external_kb26/results/cmd_zeropoint"
OUT = REPO / "validation/external_kb26/drift_what_it_is_not.json"

MAG_EDGES = [10.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 20.0]
MIN_CELL = 8
MIN_FIT = 60
MATCH_ARCSEC = 1.5
#: `CAT` 확장이 싣고 있는 밝기 열. 앞 여섯은 반지름이 커지는 고정 조리개이고
#: 마지막 `flux` 는 Kron 조리개다(별마다 크기가 다르므로 따로 본다).
APER_COLUMNS = [f"fluxaper{i}" for i in range(1, 7)] + ["flux"]
BAND_FOLDER = {"r": "rp", "i": "ip"}


def _robust_fit(x: np.ndarray, y: np.ndarray):
    """영점 보정이 쓰는 바로 그 적합. 같은 코드 경로로 재야 비교가 성립한다."""
    return ZeropointCalibrationRunner._robust_linfit(
        _Silent(), x, y, w=np.ones(len(x)), clip_sigma=3.0, iters=5,
        slope_absmax=0.8, min_n=10)


class _Silent:
    """`_robust_linfit` 이 부르는 로그 자리만 채운다."""

    def _log(self, *a, **k) -> None:
        return None


def _median_instrumental(folder: str, column: str, root: Path = BANZAI):
    """그 밴드의 프레임들에서 별마다의 기기 등급 중앙값.

    프레임마다 검출된 별이 조금씩 다르므로 첫 프레임을 기준 목록으로 삼고
    나머지를 하늘 좌표로 잇는다. 세 장 미만에서만 보인 별은 버린다.
    """
    frames = sorted((root / folder).glob("*.fits.fz"))
    base, cols = None, []
    for path in frames:
        with fits.open(path, memmap=False) as hdul:
            cat = hdul["CAT"].data
            exptime = float(hdul["SCI"].header.get("EXPTIME") or 1.0)
        tab = pd.DataFrame({"ra": np.asarray(cat["ra"], float),
                            "dec": np.asarray(cat["dec"], float),
                            "f": np.asarray(cat[column], float),
                            "flag": np.asarray(cat["flag"], int)})
        tab = tab[(tab["flag"] == 0) & np.isfinite(tab["f"]) & (tab["f"] > 0)]
        if tab.empty:
            continue
        tab["m"] = -2.5 * np.log10(tab["f"].to_numpy(float) / exptime)
        if base is None:
            base = tab[["ra", "dec"]].reset_index(drop=True)
            cols.append(tab["m"].to_numpy(float))
            continue
        k, sep, _ = SkyCoord(base["ra"], base["dec"], unit="deg").match_to_catalog_sky(
            SkyCoord(tab["ra"].to_numpy(float), tab["dec"].to_numpy(float), unit="deg"))
        v = np.full(len(base), np.nan)
        ok = np.asarray(sep.arcsec) <= MATCH_ARCSEC
        v[ok] = tab["m"].to_numpy(float)[np.asarray(k)[ok]]
        cols.append(v)
    if base is None:
        return None
    arr = np.vstack(cols)
    with np.errstate(invalid="ignore"):
        base["m"] = np.nanmedian(arr, axis=0)
        base["n"] = np.sum(np.isfinite(arr), axis=0)
    return base[np.isfinite(base["m"]) & (base["n"] >= 3)].reset_index(drop=True)


def _calibrators():
    """영점 보정이 쓴 보정성 별 표와 적합 계수."""
    cal = pd.read_csv(RESULTS / "gaia_sdss_calibrator_by_ID.csv")
    coeff = pd.read_csv(RESULTS / "zp_fit_coefficients.csv").set_index("filter")
    try:
        good = np.asarray(gaia_quality_report(cal, cstar_nsigma=None)[0], bool)
    except (KeyError, ValueError, TypeError) as exc:
        # 품질 표가 못 만들어지면 전부 쓴다. **어느 갈래로 갔는지 적는다** (C-237).
        print(f"  [주의] Gaia 품질 표를 못 만들어 별을 안 거르고 전부 쓴다 — {exc}")
        good = np.ones(len(cal), bool)
    return cal, coeff, good


def _drift_span(resid: np.ndarray, mag: np.ndarray) -> tuple[float, float, int]:
    """등급 구간 중앙값의 폭과 어두운 쪽에서 밝은 쪽을 뺀 값."""
    meds = [float(np.median(resid[(mag >= lo) & (mag < hi)]))
            for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:])
            if int(((mag >= lo) & (mag < hi)).sum()) >= MIN_CELL]
    if len(meds) < 2:
        return float("nan"), float("nan"), len(meds)
    return float(max(meds) - min(meds)), float(meds[-1] - meds[0]), len(meds)


def _drift_for(inst_tab, cal_sub, ref, colour) -> tuple[float, float, int]:
    """그 밝기 열로 다시 맞춘 영점의 잔차가 밝기를 따라 얼마나 미끄러지나."""
    ca = SkyCoord(cal_sub["ra_deg"].to_numpy(float),
                  cal_sub["dec_deg"].to_numpy(float), unit="deg")
    k, sep, _ = ca.match_to_catalog_sky(
        SkyCoord(inst_tab["ra"].to_numpy(float),
                 inst_tab["dec"].to_numpy(float), unit="deg"))
    ok = np.asarray(sep.arcsec) <= MATCH_ARCSEC
    inst = inst_tab["m"].to_numpy(float)[np.asarray(k)[ok]]
    delta = ref[ok] - inst
    col, mag = colour[ok], ref[ok]
    g = np.isfinite(delta) & np.isfinite(col) & np.isfinite(mag)
    if int(g.sum()) < MIN_FIT:
        return float("nan"), float("nan"), int(g.sum())
    zp, ct, _, _ = _robust_fit(col[g], delta[g])
    if not np.isfinite(zp):
        return float("nan"), float("nan"), int(g.sum())
    resid = delta[g] - (zp + ct * col[g])
    resid -= np.median(resid)
    span, signed, _ = _drift_span(resid, mag[g])
    return span, signed, int(g.sum())


# ---------------------------------------------------------------------------
# 1·2. 조리개를 바꾸면 치우침이 어떻게 되나, 그중 하늘 몫은 얼마인가
# ---------------------------------------------------------------------------


def part_aperture(rows: list[dict]) -> None:
    print("=== 1. 조리개를 바꾸면 밝기 치우침이 달라지나 ===")
    print("조리개를 키울수록 커지면 하늘을 더 빼는 것이고, 줄일수록 커지면 별 중심에서")
    print("빛이 새는 것이며, 무관하면 담긴 빛의 양과 상관없다는 뜻이다.")
    print()
    cal, coeff, good = _calibrators()
    print(f"{'밴드':<5}{'조리개':<11}{'담긴 빛':>9}{'별수':>7}{'폭':>9}{'어두운−밝은':>12}")
    print("-" * 55)
    for band, folder in BAND_FOLDER.items():
        colour_col = f"color_{coeff.loc[band, 'color_col']}"
        ref_col = f"ref_{band}"
        sub = cal[good
                  & np.isfinite(pd.to_numeric(cal[ref_col], errors="coerce"))
                  & np.isfinite(pd.to_numeric(cal[colour_col], errors="coerce"))]
        ref = pd.to_numeric(sub[ref_col], errors="coerce").to_numpy(float)
        colour = pd.to_numeric(sub[colour_col], errors="coerce").to_numpy(float)
        kron = _median_instrumental(folder, "flux")
        for column in APER_COLUMNS:
            inst = kron if column == "flux" else _median_instrumental(folder, column)
            if inst is None:
                print(f"  [{band}] {column}: 프레임이 없어 건너뛴다")
                continue
            span, signed, n = _drift_for(inst, sub, ref, colour)
            if not np.isfinite(span):
                continue
            frac = _enclosed_fraction(folder, column)
            shown = "—" if column == "flux" else f"{frac:.0%}"
            print(f"{band:<5}{column:<11}{shown:>9}{n:>7}{span:>9.4f}{signed:>12.4f}")
            rows.append(dict(part=1, band=band, column=column, enclosed=frac,
                             n=n, span=span, faint_minus_bright=signed))
        print()


def _enclosed_fraction(folder: str, column: str, root: Path = BANZAI) -> float:
    """Kron 밝기를 1 로 놓았을 때 그 조리개에 담긴 비율. Kron 자신은 뜻이 없다."""
    if column == "flux":
        return float("nan")
    vals: list[np.ndarray] = []
    for path in sorted((root / folder).glob("*.fits.fz")):
        with fits.open(path, memmap=False) as hdul:
            cat = hdul["CAT"].data
        f = np.asarray(cat["flux"], float)
        a = np.asarray(cat[column], float)
        flag = np.asarray(cat["flag"], int)
        peak = np.asarray(cat["peak"], float)
        ok = (flag == 0) & np.isfinite(f) & (f > 0) & np.isfinite(a) & (peak < 50000)
        vals.append(a[ok] / f[ok])
    return float(np.nanmedian(np.concatenate(vals))) if vals else float("nan")


def part_sky_share(rows: list[dict]) -> None:
    """치우침에서 조리개 넓이에 비례하는 몫을 떼어 낸다.

    하늘을 화소마다 Δs 만큼 더 빼면 별 하나가 잃는 빛은 **Δs × 조리개 넓이** 이고,
    담긴 비율 φ 로 나눠 등급으로 옮기면 하늘 몫은 **넓이/φ 에 비례**한다. 그러니
    조리개 여섯에서 잰 치우침을 넓이/φ 에 직선으로 맞추면 **기울기가 하늘 몫이고
    절편이 조리개와 무관하게 남는 몫**이다.

    넓이는 가우스 PSF 를 가정해 담긴 비율에서 되짚는다 — φ = 1 − exp(−r²/2σ²) 이라
    넓이 ∝ −ln(1−φ) 다. 조리개끼리의 비만 쓰므로 σ 를 몰라도 된다.
    """
    print("=== 2. 그 치우침에서 하늘을 더 뺀 몫은 얼마인가 ===")
    print("하늘 몫은 조리개 넓이에 비례하므로, 넓이/담김에 직선으로 맞추면")
    print("기울기가 하늘 몫이고 절편이 조리개와 무관하게 남는 몫이다.")
    print()
    got = [r for r in rows if r["part"] == 1 and r["column"] != "flux"]
    for band in BAND_FOLDER:
        sel = [r for r in got if r["band"] == band]
        phi = np.array([r["enclosed"] for r in sel], float)
        drift = np.array([r["span"] for r in sel], float)
        ok = np.isfinite(phi) & np.isfinite(drift) & (phi > 0) & (phi < 1)
        if int(ok.sum()) < 3:
            print(f"  [{band}] 조리개가 셋도 안 남아 쪼갤 수 없다")
            continue
        area = -np.log(1.0 - phi[ok])
        x = area / phi[ok]
        slope, const = np.polyfit(x, drift[ok], 1)
        mid = x[min(2, len(x) - 1)]          # APEX 가 쓰는 것과 비슷한 중간 조리개
        total = const + slope * mid
        share = float(slope * mid / total * 100) if total else float("nan")
        print(f"== {band} 밴드")
        print(f"   직선 맞춤     치우침 = {const:.4f} + {slope:.4f} × (넓이/담김)")
        print(f"   중간 조리개   전체 {total:.4f}"
              f" = 하늘 몫 {slope*mid:.4f} + 조리개와 무관한 몫 {const:.4f}")
        print(f"   하늘이 설명하는 비율  {share:.0f} %")
        print(f"   순수 하늘 가설이면 조리개 1→6 에서 {x[-1]/x[0]:.2f} 배여야 하는데"
              f" 실제 {drift[ok][-1]/drift[ok][0]:.2f} 배")
        print()
        rows.append(dict(part=2, band=band, sky_slope=float(slope),
                         aperture_free=float(const), total_at_mid=float(total),
                         sky_share_percent=share))


# ---------------------------------------------------------------------------
# 3. 보정에 쓴 별들이 포화·비선형 구간에 걸쳐 있나
# ---------------------------------------------------------------------------


def part_peak(rows: list[dict]) -> None:
    print("=== 3. 보정성 별의 봉우리 계수가 포화 구간에 걸쳐 있나 ===")
    print("검출기의 비선형은 우물이 찰수록 커지므로, 별들이 우물 바닥에만 있으면")
    print("적어도 「흔히 말하는 포화 근처 비선형」은 원인이 될 수 없다.")
    print()
    for label, root in (("kb26 (lsc, 05-03)", BANZAI), ("kb27 (ogg, 04-15)", BANZAI_KB27)):
        files = sorted((root / "rp").glob("*.fits.fz"))
        if not files:
            print(f"  {label}: 프레임이 없다")
            continue
        with fits.open(files[0], memmap=False) as hdul:
            head, cat = hdul["SCI"].header, hdul["CAT"].data
            saturate = float(head.get("SATURATE") or 0)
            fwhm = float(head.get("L1FWHM") or 0)
            exptime = float(head.get("EXPTIME") or 1)
            peak = np.asarray(cat["peak"], float)
            flux = np.asarray(cat["flux"], float)
            flag = np.asarray(cat["flag"], int)
        ok = (flag == 0) & np.isfinite(peak) & np.isfinite(flux) & (flux > 0)
        mag, pk = -2.5 * np.log10(flux[ok] / exptime), peak[ok]
        order = np.argsort(mag)
        print(f"== {label}   SATURATE={saturate:.0f} e-  L1FWHM={fwhm:.2f}″"
              f"  EXP={exptime:.1f}s  별 {int(ok.sum())}")
        for name, sl in (("가장 밝은 20", order[:20]),
                         ("밝은 5분위", order[:max(1, len(order) // 5)]),
                         ("어두운 5분위", order[-max(1, len(order) // 5):])):
            med, top = float(np.median(pk[sl])), float(np.max(pk[sl]))
            print(f"   {name:<13} 봉우리 중앙 {med:>9.0f}  최대 {top:>9.0f}"
                  f"   (우물의 {med/saturate*100:.1f}% / {top/saturate*100:.1f}%)")
            rows.append(dict(part=3, instrument=label.split()[0], group=name,
                             peak_median=med, peak_max=top, saturate=saturate,
                             percent_of_well=float(med / saturate * 100)))
        print()


# ---------------------------------------------------------------------------
# 4. 한 장에서도 이미 있나
# ---------------------------------------------------------------------------


def part_single_frame(rows: list[dict]) -> None:
    print("=== 4. 한 장에서도 이미 있나 — 여섯 장을 합쳐서 생기는 것인가 ===")
    print()
    cal, coeff, good = _calibrators()
    band, folder = "r", "rp"
    colour_col = f"color_{coeff.loc[band, 'color_col']}"
    sub = cal[good
              & np.isfinite(pd.to_numeric(cal[f"ref_{band}"], errors="coerce"))
              & np.isfinite(pd.to_numeric(cal[colour_col], errors="coerce"))]
    ref = pd.to_numeric(sub[f"ref_{band}"], errors="coerce").to_numpy(float)
    colour = pd.to_numeric(sub[colour_col], errors="coerce").to_numpy(float)
    ca = SkyCoord(sub["ra_deg"].to_numpy(float),
                  sub["dec_deg"].to_numpy(float), unit="deg")

    print(f"{'프레임':<44}{'별수':>6}{'폭':>9}{'대기량':>8}")
    print("-" * 68)
    spans = []
    for path in sorted((BANZAI / folder).glob("*.fits.fz")):
        with fits.open(path, memmap=False) as hdul:
            cat = hdul["CAT"].data
            exptime = float(hdul["SCI"].header.get("EXPTIME") or 1.0)
            airmass = float(hdul["SCI"].header.get("AIRMASS") or np.nan)
        tab = pd.DataFrame({"ra": np.asarray(cat["ra"], float),
                            "dec": np.asarray(cat["dec"], float),
                            "f": np.asarray(cat["flux"], float),
                            "flag": np.asarray(cat["flag"], int)})
        tab = tab[(tab["flag"] == 0) & np.isfinite(tab["f"]) & (tab["f"] > 0)]
        k, sep, _ = ca.match_to_catalog_sky(
            SkyCoord(tab["ra"].to_numpy(float), tab["dec"].to_numpy(float), unit="deg"))
        ok = np.asarray(sep.arcsec) <= MATCH_ARCSEC
        inst = np.full(len(sub), np.nan)
        inst[ok] = -2.5 * np.log10(tab["f"].to_numpy(float)[np.asarray(k)[ok]] / exptime)
        delta = ref - inst
        g = np.isfinite(delta) & np.isfinite(colour) & np.isfinite(ref)
        if int(g.sum()) < MIN_FIT:
            continue
        zp, ct, _, _ = _robust_fit(colour[g], delta[g])
        resid = delta[g] - (zp + ct * colour[g])
        resid -= np.median(resid)
        span, signed, _ = _drift_span(resid, ref[g])
        spans.append(span)
        print(f"{path.name:<44}{int(g.sum()):>6}{span:>9.4f}{airmass:>8.2f}")
        rows.append(dict(part=4, frame=path.name, n=int(g.sum()), span=span,
                         airmass=airmass, faint_minus_bright=signed))
    if spans:
        print()
        print(f"한 장씩 잰 폭의 중앙값 {np.median(spans):.4f}"
              f" · 가장 작은 것 {min(spans):.4f} · 가장 큰 것 {max(spans):.4f}")
        rows.append(dict(part=4, frame="중앙값", span=float(np.median(spans)),
                         n=len(spans)))
    print()


# ---------------------------------------------------------------------------
# 5. 별빛이 퍼진 모양이 밝기에 따라 달라지나
# ---------------------------------------------------------------------------


def part_curve_of_growth(rows: list[dict]) -> None:
    print("=== 5. 별빛이 퍼진 모양이 밝기에 따라 달라지나 ===")
    print("조리개 측광은 고정된 원 안의 빛을 재므로, 밝은 별과 어두운 별의 퍼진 모양이")
    print("다르면 그것이 곧 등급 치우침이 된다. 모양이 밝기와 무관하면 여기서 갈린다.")
    print()
    for label, root in (("kb26 rp", BANZAI), ("kb27 rp", BANZAI_KB27)):
        files = sorted((root / "rp").glob("*.fits.fz"))
        if not files:
            print(f"  {label}: 프레임이 없다")
            continue
        mags, fracs, peaks = [], {i: [] for i in range(1, 7)}, []
        for path in files:
            with fits.open(path, memmap=False) as hdul:
                cat = hdul["CAT"].data
                exptime = float(hdul["SCI"].header.get("EXPTIME") or 1.0)
            flux = np.asarray(cat["flux"], float)
            peak = np.asarray(cat["peak"], float)
            flag = np.asarray(cat["flag"], int)
            ok = (flag == 0) & np.isfinite(flux) & (flux > 0) & (peak < 50000)
            mags.append(-2.5 * np.log10(flux[ok] / exptime))
            peaks.append(peak[ok])
            for i in range(1, 7):
                fracs[i].append(np.asarray(cat[f"fluxaper{i}"], float)[ok] / flux[ok])
        mag = np.concatenate(mags)
        peak = np.concatenate(peaks)
        frac = {i: np.concatenate(v) for i, v in fracs.items()}
        print(f"== {label}   (봉우리 50,000 e- 넘는 별은 뺐다, N={len(mag)})")
        print(f"{'기기 등급':<12}{'별수':>6}"
              + "".join(f"{'조리개' + str(i):>9}" for i in range(1, 7))
              + f"{'봉우리':>10}")
        print("-" * 82)
        # 기기 등급은 음수라 고정 구간이 안 맞으므로 분위로 자른다.
        edges = np.nanpercentile(mag, [0, 20, 40, 60, 80, 100])
        first = last = None
        for lo, hi in zip(edges[:-1], edges[1:]):
            k = (mag >= lo) & (mag <= hi)
            if int(k.sum()) < 20:
                continue
            vals = [float(np.nanmedian(frac[i][k])) for i in range(1, 7)]
            first = vals if first is None else first
            last = vals
            print(f"{lo:.1f}~{hi:.1f}".ljust(12) + f"{int(k.sum()):>6}"
                  + "".join(f"{v:>9.3f}" for v in vals)
                  + f"{np.nanmedian(peak[k]):>10.0f}")
            rows.append(dict(part=5, instrument=label, mag_lo=float(lo),
                             mag_hi=float(hi), n=int(k.sum()), enclosed=vals))
        if first and last:
            worst = max(abs(a - b) for a, b in zip(first, last))
            print(f"   가장 밝은 구간과 가장 어두운 구간의 담긴 비율 차이 최대 {worst:.3f}")
            rows.append(dict(part=5, instrument=label, biggest_change=float(worst)))
        print()


def main() -> int:
    rows: list[dict] = []
    part_aperture(rows)
    part_sky_share(rows)
    part_peak(rows)
    part_single_frame(rows)
    part_curve_of_growth(rows)
    OUT.write_text(json.dumps(rows, ensure_ascii=False, indent=2, default=str),
                   encoding="utf-8")
    print(f"썼다: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
