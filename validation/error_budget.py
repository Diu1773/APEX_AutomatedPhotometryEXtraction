"""보정된 등급의 오차 예산 — 성분으로 쪼개고, 합이 관측을 설명하는지 본다.

**지금까지 낸 것은 조각들이었다.** 계수가 시야마다 다르다, 절단이 그것을 흔든다,
두 기기가 0.27 어긋난다, 밝기에 따라 0.28 휜다. 각각은 재서 낸 수인데
**「그래서 이 등급의 오차는 얼마이고 어디서 왔나」를 합쳐서 낸 적이 없다.**
이 스크립트가 그것을 한다.

## 모형

    mag_cal = mag_inst + zp + ct·C          (C 는 표준 색지수)

한 별의 보정 등급이 틀릴 수 있는 경로는 넷이다.

    σ_phot     광자·읽기 잡음. 별마다 다르고 `mag_inst_err` 에 있다.
    σ_cal(C)   영점과 색항이 함께 갖는 불확실도. **색에 따라 달라진다** —
               색 지렛대의 받침점에서 멀수록 커진다.
    σ_frame    프레임마다 영점이 흔들리는 몫. 여러 장을 중앙값으로 합치므로
               √N 으로 줄어든다.
    σ_ref      Gaia→표준 변환식 자체의 산포. **같은 기준을 쓰는 두 기기를 견줄
               때는 상쇄되고, 바깥 목록과 견줄 때만 산다.**

## σ_cal 을 어떻게 재나 — 부트스트랩

색항의 불확실도를 적합 공분산에서 뽑으면 **절단이 만드는 흔들림이 안 잡힌다.**
2026-09-09 에 잰 것이 바로 그것이었다 — MuSCAT3 의 r 은 같은 자료에서 절단 하나로
계수가 −0.381 에서 +0.021 로 0.40 움직인다.

그래서 **별을 복원추출로 다시 뽑아 APEX 의 적합기로 처음부터 다시 맞춘다.** 절단도
매번 새로 돈다. 그러면 「이 자료로 이 적합기를 돌리면 계수가 얼마나 흔들리는가」가
기제를 몰라도 수치로 나온다. 300 번 돌려 색마다 16~84 % 폭을 σ_cal 로 쓴다.

## 닫는 시험

같은 별을 두 기기로 잰 차이의 산포가 **예측값과 맞는가.** 맞으면 예산이 완성된
것이고, 관측이 더 크면 **이름 없는 항이 남아 있다는 뜻**이므로 그 크기를 적는다.

실행:
    python -X utf8 validation/error_budget.py
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
sys.path.insert(0, str(REPO / "validation" / "external_muscat3"))

from apex.analysis.cmd.zeropoint_runner import (  # noqa: E402
    ZeropointCalibrationRunner,
)
from apex.utils.gaia_quality import gaia_quality_report  # noqa: E402
from apex.utils.gaia_transforms import GAIA_TO_BAND, BAND_ALIASES  # noqa: E402

N_BOOT = 300
SNR_CUT = 20.0
RNG = np.random.default_rng(20260909)

ROOTS = ("E:/APEX_validation/reprocess/*/result", "E:/observed_Analysis/*/*/result")
EXTRA = {"MuSCAT3": REPO / "validation/external_muscat3/results",
         "kb26": REPO / "validation/external_kb26/results"}
FIELDS = ("M13", "M3", "M5", "M67", "M37", "NGC6811", "NGC457")


class _Fitter:
    """APEX 의 적합기를 그대로 쓰기 위한 껍데기.

    `_robust_linfit` 은 러너의 메서드인데 `self` 에서 쓰는 것이 로그뿐이다.
    다시 구현하면 「APEX 의 오차」가 아니라 「내 재구현의 오차」를 재게 되므로
    실제 메서드를 빌려 온다.
    """

    def _log(self, *a, **k):
        pass

    fit = ZeropointCalibrationRunner._robust_linfit


FIT = _Fitter()


def field_of(path: str) -> str:
    for name in FIELDS:
        if re.search(rf"[\\/]{name}[\\/]", path):
            return name
    return "?"


def load(result_dir: Path):
    zp_dir = Path(result_dir) / "cmd_zeropoint"
    cal_p = zp_dir / "gaia_sdss_calibrator_by_ID.csv"
    co_p = zp_dir / "zp_fit_coefficients.csv"
    if not (cal_p.exists() and co_p.exists()):
        return None
    try:
        cal = pd.read_csv(cal_p)
        co = pd.read_csv(co_p).set_index("filter")
    except Exception:  # noqa: BLE001
        return None
    frame = None
    fp = zp_dir / "frame_zeropoint.csv"
    if fp.exists():
        try:
            frame = pd.read_csv(fp)
        except Exception:  # noqa: BLE001
            frame = None
    return cal, co, frame


def bootstrap_sigma(x, y, w, colours, n_boot=N_BOOT):
    """색마다의 σ_cal — 별을 다시 뽑아 매번 처음부터 맞춘 값의 16~84 % 반폭.

    반환: (색마다의 σ, ct 의 16~84 % 반폭, 성공한 재표집 수)
    """
    n = len(x)
    if n < 20:
        return np.full(len(colours), np.nan), np.nan, 0
    preds, cts = [], []
    for _ in range(n_boot):
        idx = RNG.integers(0, n, n)
        zp, ct, n_in, _ = FIT.fit(x[idx], y[idx], w=w[idx],
                                  clip_sigma=3.0, iters=5,
                                  slope_absmax=0.8, min_n=10)
        if not np.isfinite(zp) or not np.isfinite(ct):
            continue
        preds.append(zp + ct * colours)
        cts.append(ct)
    if len(preds) < 30:
        return np.full(len(colours), np.nan), np.nan, len(preds)
    P = np.asarray(preds)
    lo, hi = np.percentile(P, [16, 84], axis=0)
    ct_lo, ct_hi = np.percentile(cts, [16, 84])
    return (hi - lo) / 2.0, (ct_hi - ct_lo) / 2.0, len(preds)


def budget_for(label: str, result_dir: Path) -> list[dict]:
    got = load(result_dir)
    if got is None:
        return []
    cal, co, frame = got
    try:
        qual = np.asarray(gaia_quality_report(cal, cstar_nsigma=None)[0], bool)
    except Exception:  # noqa: BLE001
        qual = np.ones(len(cal), bool)

    rows = []
    for band in co.index:
        band = str(band)
        colour = str(co.loc[band, "color_col"])
        dc, cc, ec, sc = (f"delta_{band}", f"color_{colour}",
                          f"mag_inst_err_{band}", f"snr_{band}")
        if dc not in cal.columns or cc not in cal.columns or ec not in cal.columns:
            continue
        x = pd.to_numeric(cal[cc], errors="coerce").to_numpy(float)
        y = pd.to_numeric(cal[dc], errors="coerce").to_numpy(float)
        e = pd.to_numeric(cal[ec], errors="coerce").to_numpy(float)
        s = (pd.to_numeric(cal[sc], errors="coerce").to_numpy(float)
             if sc in cal.columns else np.full(len(cal), np.inf))
        w = 1.0 / np.maximum(e, 1e-6) ** 2
        m = (np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0)
             & np.isfinite(s) & (s >= SNR_CUT) & qual)
        if int(m.sum()) < 30:
            continue

        c_lo, c_mid, c_hi = np.percentile(x[m], [5, 50, 95])
        grid = np.array([c_lo, c_mid, c_hi])
        sig_cal, sig_ct, n_ok = bootstrap_sigma(x[m], y[m], w[m], grid)

        # σ_phot — 보정성 별의 중앙과 어두운 쪽
        sig_phot_med = float(np.median(e[m]))
        base = np.isfinite(e) & (e > 0)
        sig_phot_faint = float(np.percentile(e[base], 90))

        # σ_frame — 프레임 사이 영점 산포를 √N 으로 나눈다
        sig_frame, n_frame = np.nan, 0
        if frame is not None and "filter" in frame.columns and "zp_frame" in frame.columns:
            fz = pd.to_numeric(
                frame.loc[frame["filter"].astype(str) == band, "zp_frame"],
                errors="coerce").to_numpy(float)
            fz = fz[np.isfinite(fz)]
            if fz.size >= 2:
                n_frame = int(fz.size)
                sig_frame = float(1.4826 * np.median(np.abs(fz - np.median(fz)))
                                  / np.sqrt(fz.size))

        # σ_repeat — **별마다의 실제 재현성.** APEX 가 프레임마다 기준별들이
        # 영점 둘레로 얼마나 흩어지는지 재어 `zp_qc_summary.csv` 의
        # `median_frame_zp_scatter` 에 적어 두는데, 등급 오차에는 안 넣는다.
        # 같은 SNR 의 광자 잡음(1.0857/SNR)보다 2~10 배 크다.
        sig_repeat, snr_ref_med = np.nan, np.nan
        qp = Path(result_dir) / "cmd_zeropoint" / "zp_qc_summary.csv"
        if qp.exists():
            try:
                q = pd.read_csv(qp)
                hit = q[q["filter"].astype(str) == band]
                if len(hit):
                    sig_repeat = float(hit.iloc[0].get("median_frame_zp_scatter", np.nan))
                    snr_ref_med = float(hit.iloc[0].get("median_snr_ref", np.nan))
            except Exception:  # noqa: BLE001
                pass

        # σ_ref — 변환식 자체의 산포 (같은 기준끼리는 상쇄된다)
        entry = GAIA_TO_BAND.get(BAND_ALIASES.get(band, band))
        sig_ref = float(entry[4]) if entry else np.nan

        rows.append(dict(
            field=label, band=band, colour=colour, n_fit=int(m.sum()),
            colour_lo=float(c_lo), colour_mid=float(c_mid), colour_hi=float(c_hi),
            ct=float(co.loc[band, "ct"]), sigma_ct=float(sig_ct),
            sigma_cal_lo=float(sig_cal[0]), sigma_cal_mid=float(sig_cal[1]),
            sigma_cal_hi=float(sig_cal[2]),
            sigma_phot_med=sig_phot_med, sigma_phot_faint=sig_phot_faint,
            sigma_frame=sig_frame, n_frame=n_frame, sigma_ref=sig_ref,
            sigma_repeat=sig_repeat, snr_ref_med=snr_ref_med,
            sigma_photon_expected=(1.0857 / snr_ref_med
                                   if np.isfinite(snr_ref_med) and snr_ref_med > 0
                                   else np.nan),
            n_boot_ok=n_ok,
        ))
    return rows


# ---------------------------------------------------------------------------
# 닫는 시험 — 예산의 합이 실제로 관측된 기기 사이 차이를 설명하는가
# ---------------------------------------------------------------------------

#: 예산표는 사장님 M67 워크스페이스를 시야 이름 "M67" 로 담는다.
LABEL = {"Moravian": "M67"}

M67_TRIO = {
    "Moravian": Path("E:/APEX_validation/reprocess/M67/result"),
    "MuSCAT3": REPO / "validation/external_muscat3/results",
    "kb26": REPO / "validation/external_kb26/results",
}


def _row(d: pd.DataFrame, field: str, band: str):
    hit = d[(d["field"] == field) & (d["band"] == band)]
    return hit.iloc[0] if len(hit) else None


def closing_test(d: pd.DataFrame) -> list[dict]:
    """같은 별을 두 기기로 잰 차이의 산포 vs 예산이 예측한 값.

    **별마다 자기 값으로 예측한다.** 처음 판은 σ_phot 을 보정성 별(SNR>=20 이라
    밝다)의 중앙값으로 쓰면서 산포는 어두운 별까지 포함한 짝 전체에서 쟀다.
    그러면 예측이 낮게 나오고 「이름 없는 항」이 있는 것처럼 보인다 — 실제로는
    내가 두 표본을 어긋나게 잡은 것이다.

    별마다:

        σ_pred = sqrt( e_A² + e_B² + σ_cal_A(C)² + σ_cal_B(C)² + σ_frame_A² + σ_frame_B² )

    σ_cal 은 색에 따라 달라지므로 5·50·95 % 세 점을 선형보간해 그 별의 색에서 읽는다.
    같은 기준을 쓰므로 σ_ref 는 상쇄된다.
    """
    from compare_magnitudes import join, load_workspace, mad

    ws, raw = {}, {}
    for name, rd in M67_TRIO.items():
        try:
            ws[name] = load_workspace(Path(rd), name)
        except SystemExit:
            print(f"[{name}] 등급 표가 없어 건너뛴다")
            continue
        rp = Path(rd) / "cmd_zeropoint" / "median_by_ID_filter_wide_raw.csv"
        raw[name] = pd.read_csv(rp) if rp.exists() else None
    names = list(ws)
    out: list[dict] = []

    def sigma_cal_at(r, colour_vals):
        xs = np.array([r["colour_lo"], r["colour_mid"], r["colour_hi"]], float)
        ys = np.array([r["sigma_cal_lo"], r["sigma_cal_mid"], r["sigma_cal_hi"]], float)
        ok = np.isfinite(xs) & np.isfinite(ys)
        if ok.sum() < 2:
            return np.full(len(colour_vals), np.nan)
        return np.interp(colour_vals, xs[ok], ys[ok])

    def err_of(name, band, ids):
        """그 워크스페이스에서 별마다의 mag_inst_err. 없으면 NaN."""
        t = raw.get(name)
        c = f"mag_inst_err_{band}"
        if t is None or c not in t.columns or "ID" not in t.columns:
            return np.full(len(ids), np.nan)
        lut = dict(zip(t["ID"].astype(str), pd.to_numeric(t[c], errors="coerce")))
        return np.array([lut.get(str(i), np.nan) for i in ids], float)

    print()
    print("=== 닫는 시험 — 예산이 관측을 설명하나 (M67, 같은 별) ===")
    print("별마다 자기 오차·자기 색으로 예측한다. 같은 기준이라 σ_ref 는 상쇄된다.")
    print()
    print("{:<20}{:<3}{:>5}{:>25}{:>38}".format(
        "", "", "", "── 계통 (중앙 어긋남) ──", "──────── 산포 (별마다) ────────"))
    print("{:<20}{:<3}{:>5}{:>9}{:>9}{:>7}{:>10}{:>9}{:>9}{:>10}".format(
        "두 기기", "밴드", "별수", "예측", "관측", "비",
        "관측", "광자만", "재현성/√N", "재현성"))
    print("-" * 93)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            j = join(ws[a], ws[b])
            num = lambda c: pd.to_numeric(j[c], errors="coerce").to_numpy(float)
            for band in ("g", "r", "i", "B", "V"):
                ca, cb = f"mag_cal_{band}_A", f"mag_cal_{band}_B"
                if ca not in j.columns or cb not in j.columns:
                    continue
                ra = _row(d, LABEL.get(a, a), band)
                rb = _row(d, LABEL.get(b, b), band)
                if ra is None or rb is None:
                    continue
                diff = num(ca) - num(cb)
                # 색은 A 쪽 보정 등급에서 만든다 (적합에 쓴 색축과 같은 조합).
                col_name = str(ra["colour"])
                p1, p2 = col_name.split("_", 1)
                c1, c2 = f"mag_cal_{p1}_A", f"mag_cal_{p2}_A"
                colour = (num(c1) - num(c2)) if (c1 in j.columns and c2 in j.columns) else None
                if colour is None:
                    continue
                ea = err_of(a, band, j["ID_A"]) if "ID_A" in j.columns else np.full(len(j), np.nan)
                eb = err_of(b, band, j["ID_B"]) if "ID_B" in j.columns else np.full(len(j), np.nan)
                m = np.isfinite(diff) & np.isfinite(colour) & np.isfinite(ea) & np.isfinite(eb)
                if m.sum() < 20:
                    continue
                # **계통과 산포는 다른 것이 만든다.**
                # σ_cal 은 영점·색항의 불확실도라 그 워크스페이스의 별을 모두
                # 똑같이 민다 — 중앙값을 옮기지 별마다의 산포를 만들지 않는다.
                # 산포를 만드는 것은 별마다 다른 것, 즉 광자 잡음이다.
                sys_pred = float(np.sqrt(
                    np.nanmedian(sigma_cal_at(ra, colour[m])) ** 2
                    + np.nanmedian(sigma_cal_at(rb, colour[m])) ** 2
                    + np.nan_to_num(ra["sigma_frame"]) ** 2
                    + np.nan_to_num(rb["sigma_frame"]) ** 2))
                sys_obs = float(abs(np.median(diff[m])))
                sca_pred = float(np.nanmedian(np.sqrt(ea[m] ** 2 + eb[m] ** 2)))
                sca_obs = mad(diff[m])
                # 실제 재현성으로 예측한 두 가지 — 프레임 수로 평균되는 경우와
                # 별마다 붙박이라 안 줄어드는 경우.
                rep_a, rep_b = ra.get("sigma_repeat"), rb.get("sigma_repeat")
                na = max(int(ra.get("n_frame") or 1), 1)
                nb = max(int(rb.get("n_frame") or 1), 1)
                rep_avg = float(np.sqrt(np.nan_to_num(rep_a) ** 2 / na
                                        + np.nan_to_num(rep_b) ** 2 / nb))
                rep_fix = float(np.sqrt(np.nan_to_num(rep_a) ** 2
                                        + np.nan_to_num(rep_b) ** 2))
                print("{:<20}{:<3}{:>5}{:>9.4f}{:>9.4f}{:>7.1f}"
                      "{:>10.4f}{:>9.4f}{:>9.4f}{:>10.4f}".format(
                          f"{a} vs {b}", band, int(m.sum()),
                          sys_pred, sys_obs, sys_obs / sys_pred if sys_pred else np.nan,
                          sca_obs, sca_pred, rep_avg, rep_fix))
                out.append(dict(pair=f"{a} vs {b}", band=band, n=int(m.sum()),
                                sys_predicted=sys_pred, sys_observed=sys_obs,
                                scatter_observed=sca_obs,
                                scatter_photon_only=sca_pred,
                                scatter_repeat_averaged=rep_avg,
                                scatter_repeat_fixed=rep_fix,
                                sys_ratio=sys_obs / sys_pred if sys_pred else None,
                                scatter_ratio=sca_obs / sca_pred if sca_pred else None))
    if out:
        sr = np.array([o["sys_ratio"] for o in out if o["sys_ratio"]], float)
        cr = np.array([o["scatter_ratio"] for o in out if o["scatter_ratio"]], float)
        print()
        print(f"  계통 관측/예측 중앙 {np.median(sr):.1f} (최소 {sr.min():.1f} · 최대 {sr.max():.1f})")
        print(f"  산포 관측/예측 중앙 {np.median(cr):.1f} (최소 {cr.min():.1f} · 최대 {cr.max():.1f})")
        print("  1 에 가까우면 그 몫은 설명된 것이다. 크면 이름 없는 항이 남아 있다.")
        for key, label in (("scatter_photon_only", "광자만"),
                           ("scatter_repeat_averaged", "재현성/√N"),
                           ("scatter_repeat_fixed", "재현성(안 줄어듦)")):
            rr = np.array([o["scatter_observed"] / o[key] for o in out
                           if o.get(key)], float)
            left = np.array([np.sqrt(max(o["scatter_observed"] ** 2 - o[key] ** 2, 0.0))
                             for o in out if o.get(key)], float)
            print(f"  산포 · {label:<16} 관측/예측 중앙 {np.median(rr):>5.1f} · "
                  f"설명 안 되는 몫 중앙 {np.median(left):.4f} 등급")
    return out


def crowding_test(d: pd.DataFrame) -> list[dict]:
    """설명 안 되는 산포가 혼잡 때문인가 — 가장 가까운 이웃까지의 거리로 가른다.

    닫는 시험이 남긴 것은 **별마다 다른 0.09 등급**이다. 계통이 아니므로 영점도
    색항도 아니고, 광자 잡음보다 20 배 크다. 남는 후보 중 가장 먼저 잴 것이
    혼잡이다 — 잔차가 큰 짝이 전부 화소가 거친 kb26(0.58 초각)이 낀 짝이고,
    가장 작은 짝이 둘 다 고운 Moravian(0.393)·MuSCAT3(0.266)이다.

    **각 기기 자기 목록 안에서** 가장 가까운 이웃까지의 거리를 재고, 두 기기 중
    더 가까운 쪽을 그 별의 혼잡도로 쓴다. 혼잡이 원인이면 이웃이 멀수록 차이가
    줄어야 한다.
    """
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    from compare_magnitudes import join, load_workspace, mad

    ws = {}
    for name, rd in M67_TRIO.items():
        try:
            ws[name] = load_workspace(Path(rd), name)
        except SystemExit:
            continue

    def nn_arcsec(t):
        c = SkyCoord(pd.to_numeric(t["ra_deg"], errors="coerce").to_numpy(float),
                     pd.to_numeric(t["dec_deg"], errors="coerce").to_numpy(float),
                     unit="deg")
        _, sep, _ = c.match_to_catalog_sky(c, nthneighbor=2)
        return sep.arcsec

    for name in ws:
        ws[name] = ws[name].assign(_nn=nn_arcsec(ws[name]))

    names = list(ws)
    out: list[dict] = []
    print()
    print("=== 설명 안 되는 산포가 혼잡 때문인가 ===")
    print("가장 가까운 이웃까지의 거리로 나눈 별마다의 차이 (MAD, 등급).")
    print("두 기기 중 더 가까운 이웃 거리를 그 별의 혼잡도로 쓴다.")
    edges = [0, 3, 5, 8, 12, 1e9]
    print()
    print("{:<20}{:<3}{:>6}".format("두 기기", "밴드", "전체")
          + "".join(f"{f'{lo}~{hi:g}″':>10}" for lo, hi in zip(edges[:-1], edges[1:])))
    print("-" * 74)
    for i, a in enumerate(names):
        for b in names[i + 1:]:
            j = join(ws[a], ws[b])
            num = lambda c: pd.to_numeric(j[c], errors="coerce").to_numpy(float)
            nn = np.fmin(num("_nn_A"), num("_nn_B"))
            for band in ("g", "r", "i"):
                ca, cb = f"mag_cal_{band}_A", f"mag_cal_{band}_B"
                if ca not in j.columns or cb not in j.columns:
                    continue
                diff = num(ca) - num(cb)
                base = np.isfinite(diff) & np.isfinite(nn)
                if base.sum() < 40:
                    continue
                cells, row = [], dict(pair=f"{a} vs {b}", band=band,
                                      all_mad=mad(diff[base]), bins={})
                for lo, hi in zip(edges[:-1], edges[1:]):
                    m = base & (nn >= lo) & (nn < hi)
                    if m.sum() < 8:
                        cells.append("     -")
                        continue
                    v = mad(diff[m])
                    cells.append(f"{v:.4f}({int(m.sum())})")
                    row["bins"][f"{lo:g}~{hi:g}"] = dict(n=int(m.sum()), mad=v)
                print("{:<20}{:<3}{:>6.4f}".format(f"{a} vs {b}", band,
                                                   mad(diff[base]))
                      + "".join(f"{c:>10}" for c in cells))
                out.append(row)
    return out


def main() -> int:
    rows: list[dict] = []
    seen: set[str] = set()
    for pat in ROOTS:
        for rd in sorted(glob.glob(pat)):
            f = field_of(str(Path(rd) / "x"))
            if f == "?" or f in seen:
                continue
            got = budget_for(f, Path(rd))
            if got:
                seen.add(f)
                rows += got
    for name, rd in EXTRA.items():
        rows += budget_for(name, rd)

    d = pd.DataFrame(rows)
    if d.empty:
        print("자료를 못 찾았다.")
        return 1

    # 같은 기준을 쓰는 두 기기를 견줄 때 σ_ref 는 상쇄된다.
    d["sigma_internal_mid"] = np.sqrt(
        d["sigma_cal_mid"] ** 2 + d["sigma_phot_med"] ** 2
        + np.nan_to_num(d["sigma_frame"]) ** 2)
    d["sigma_external_mid"] = np.sqrt(
        d["sigma_internal_mid"] ** 2 + np.nan_to_num(d["sigma_ref"]) ** 2)

    print("=== 보정 등급의 오차 예산 (1σ, 등급) ===")
    print("색항 불확실도는 별을 다시 뽑아 APEX 적합기로 300 번 다시 맞춘 폭이다.")
    print("σ_cal 은 색에 따라 달라지므로 색 5 % · 50 % · 95 % 세 자리에서 낸다.")
    print()
    print("{:<9}{:<4}{:>6}{:>9}{:>9}{:>27}{:>10}{:>9}{:>9}".format(
        "시야", "밴드", "별수", "ct", "σ(ct)",
        "σ_cal (색 5%·50%·95%)", "σ_phot", "σ_frame", "σ_ref"))
    print("-" * 92)
    for _, r in d.sort_values(["field", "band"]).iterrows():
        print("{:<9}{:<4}{:>6}{:>+9.3f}{:>9.3f}"
              "{:>9.4f}{:>9.4f}{:>9.4f}{:>10.4f}{:>9.4f}{:>9.4f}".format(
                  r["field"], r["band"], r["n_fit"], r["ct"], r["sigma_ct"],
                  r["sigma_cal_lo"], r["sigma_cal_mid"], r["sigma_cal_hi"],
                  r["sigma_phot_med"], r["sigma_frame"], r["sigma_ref"]))

    print()
    print("=== 무엇이 가장 크나 — 색 중앙에서 성분의 몫 ===")
    print("{:<9}{:<4}{:>12}{:>12}{:>12}{:>14}".format(
        "시야", "밴드", "σ_cal", "σ_phot", "σ_frame", "합(기기끼리)"))
    print("-" * 63)
    for _, r in d.sort_values("sigma_internal_mid", ascending=False).iterrows():
        print("{:<9}{:<4}{:>12.4f}{:>12.4f}{:>12.4f}{:>14.4f}".format(
            r["field"], r["band"], r["sigma_cal_mid"], r["sigma_phot_med"],
            r["sigma_frame"], r["sigma_internal_mid"]))

    closing = closing_test(d)
    crowd = crowding_test(d)

    out = REPO / "validation/error_budget.json"
    out.write_text(json.dumps({"budget": d.to_dict("records"),
                               "closing_test": closing,
                               "crowding": crowd},
                              ensure_ascii=False, indent=2, default=str),
                   encoding="utf-8")
    print()
    print(f"썼다: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
