# -*- coding: utf-8 -*-
"""우주선 제거 전/후 산출물을 같은 자료로 나란히 비교한다.

같은 raw 를 두 번 처리한 결과를 견준다:
    <target>/result_nocr/   cosmetic_enable=False 로 만든 기존 산출물
    <target>/result/        cosmetic_enable=True  로 다시 만든 산출물

비교 항목
  1. 검출     — 프레임당 소스 수, 점광원 잡음(FWHM<=1.5px) 비율
  2. 마스터   — 소스 수, **1프레임 검출 수**(D-004 의 근거)
  3. 측광     — 같은 별의 등급 차, 산포, SNR
  4. PSF      — 플럭스비(psf/ap), MAD, 프레임간 재현성
  5. CMD      — 별 수, 색 분포

실행:
    .venv-deploy\\Scripts\\python -X utf8 \\
        validation/psf_crossinstrument/compare_cr_effect.py M13 [NGC6811 ...]
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPROCESS = Path(r"E:\APEX_validation\reprocess")


def _fmt(v, nd=4, pct=False):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "—"
    return f"{v:.{nd}f}{'%' if pct else ''}"


def _detect_stats(result_dir: Path) -> dict:
    """검출 캐시에서 소스 수와 점광원 잡음 비율."""
    fs = sorted((result_dir / "cache").glob("detect_*.csv"))
    if not fs:
        return {}
    n_src, spike, tot = [], 0, 0
    for f in fs:
        try:
            d = pd.read_csv(f)
        except Exception:
            continue
        n_src.append(len(d))
        if "fwhm_px" in d.columns:
            v = pd.to_numeric(d["fwhm_px"], errors="coerce").dropna()
            spike += int((v <= 1.5).sum())
            tot += len(v)
    return {
        "n_frames": len(n_src),
        "n_sources_med": float(np.median(n_src)) if n_src else np.nan,
        "spike_pct": 100.0 * spike / tot if tot else np.nan,
        "n_fwhm_measured": tot,
    }


def _master_stats(result_dir: Path) -> dict:
    """마스터 카탈로그 소스 수와 1프레임 검출 수 (D-004 근거)."""
    cands = list(result_dir.rglob("ref_catalog.tsv"))
    if not cands:
        return {}
    d = pd.read_csv(cands[0], sep="\t")
    out = {"n_master": len(d)}
    if "n_det_frames" in d.columns:
        v = pd.to_numeric(d["n_det_frames"], errors="coerce")
        out["n_1frame"] = int((v == 1).sum())
        out["pct_1frame"] = 100.0 * out["n_1frame"] / max(len(d), 1)
        out["n_allframe"] = int((v == v.max()).sum())
        out["max_frames"] = int(v.max()) if len(v) else 0
    return out


def _phot_frames(result_dir: Path) -> dict[str, pd.DataFrame]:
    out = {}
    for t in sorted((result_dir / "step7_forced_phot").glob("photometry_*.tsv")):
        try:
            out[t.name] = pd.read_csv(t, sep="\t")
        except Exception:
            continue
    return out


def _psf_stats(result_dir: Path) -> dict:
    """PSF 플럭스비와 구경 대비 산포."""
    psf_dir = result_dir / "cmd_psf"
    ap_dir = result_dir / "step7_forced_phot"
    if not psf_dir.exists():
        return {}
    ratios, meds, mads = [], [], []
    for t in sorted(psf_dir.glob("photometry_*.tsv")):
        a = ap_dir / t.name
        if not a.exists():
            continue
        try:
            p = pd.read_csv(t, sep="\t").dropna(subset=["det_uid"])
            q = pd.read_csv(a, sep="\t").dropna(subset=["det_uid"])
        except Exception:
            continue
        m = p.merge(q, on="det_uid", suffixes=("_p", "_a"))
        if "flux_psf_e" in m.columns and "flux_e" in m.columns:
            # **검출된 별만 쓴다.** 강제 측광 위치(그 프레임에서 검출되지 않은
            # 어두운 별)는 구경이 거의 0 을 재므로 비율이 10~68 로 폭주한다.
            # 섞어서 중앙값을 내면 강제 측광 비율이 높은 자료가 「PSF 이상」으로
            # 보인다 — NGC6811(강제 56%)에서 실제로 그렇게 오진했다(2026-07-30).
            sel = m
            if "forced_flag" in m.columns:
                sel = m[~m["forced_flag"].astype(bool)]
            elif "detected_flag" in m.columns:
                sel = m[m["detected_flag"].astype(bool)]
            r = pd.to_numeric(sel["flux_psf_e"], errors="coerce") / pd.to_numeric(
                sel["flux_e"], errors="coerce"
            )
            r = r[np.isfinite(r) & (r > 0)]
            if len(r):
                ratios.append(float(np.median(r)))
        ok = m[
            np.isfinite(m.get("mag_psf")) & np.isfinite(m.get("mag_inst"))
            & (m.get("snr_psf", 0) > 20) & (m.get("snr", 0) > 20)
        ]
        if len(ok) > 10:
            d = (ok["mag_psf"] - ok["mag_inst"]).to_numpy(float)
            med = float(np.median(d))
            meds.append(med)
            mads.append(float(1.4826 * np.median(np.abs(d - med))))
    return {
        "n_psf_frames": len(ratios),
        "flux_ratio_med": float(np.median(ratios)) if ratios else np.nan,
        "mad_med": float(np.median(mads)) if mads else np.nan,
        "frame_rms": float(np.std(meds)) if len(meds) > 1 else np.nan,
    }


def _cmd_stats(result_dir: Path) -> dict:
    f = result_dir / "cmd_zeropoint" / "median_by_ID_filter_wide_cmd.csv"
    if not f.exists():
        return {}
    d = pd.read_csv(f)
    cols = [c for c in d.columns if c.startswith("mag_cal_")]
    out = {"n_cmd_rows": len(d), "bands": ",".join(sorted(c[8:] for c in cols))}
    finite = np.ones(len(d), bool)
    for c in cols:
        finite &= np.isfinite(pd.to_numeric(d[c], errors="coerce"))
    out["n_cmd_complete"] = int(finite.sum())
    return out


def compare_photometry(res_a: Path, res_b: Path, tol_arcsec: float = 1.0) -> dict:
    """같은 별의 등급을 프레임별로 비교 — **하늘 좌표로 짝짓는다**.

    master_id 로 조인하면 안 된다. 재처리하면 마스터 카탈로그가 다시 만들어져
    ID 가 재할당되므로(M13: 1574 -> 1501) 같은 ID 가 다른 별을 가리킨다.
    실제로 ID 조인으로 재면 MAD 가 1.24 mag 이라는 불가능한 값이 나온다.
    """
    from scipy.spatial import cKDTree

    A, B = _phot_frames(res_a), _phot_frames(res_b)
    common = sorted(set(A) & set(B))
    if not common:
        return {}
    tol_deg = tol_arcsec / 3600.0
    dmags, n_pairs, snr_a, snr_b = [], 0, [], []
    for name in common:
        a, b = A[name], B[name]
        if not {"ra_deg", "dec_deg"} <= set(a.columns) & set(b.columns):
            continue

        def _clean(df):
            d = df[
                np.isfinite(pd.to_numeric(df["ra_deg"], errors="coerce"))
                & np.isfinite(pd.to_numeric(df["dec_deg"], errors="coerce"))
                & np.isfinite(pd.to_numeric(df["mag_inst"], errors="coerce"))
                & (pd.to_numeric(df["snr"], errors="coerce") > 20)
            ]
            return d.reset_index(drop=True)

        ca, cb = _clean(a), _clean(b)
        if len(ca) < 5 or len(cb) < 5:
            continue
        # 적경은 위도에 따라 좁아진다 — 코사인 보정 후 평면 근사로 매칭
        cosd = np.cos(np.radians(float(np.median(cb["dec_deg"]))))
        pa = np.column_stack([ca["ra_deg"] * cosd, ca["dec_deg"]])
        pb = np.column_stack([cb["ra_deg"] * cosd, cb["dec_deg"]])
        dist, idx = cKDTree(pa).query(pb, k=1)
        hit = dist < tol_deg
        if not np.any(hit):
            continue
        mb = cb.loc[hit]
        ma = ca.iloc[idx[hit]]
        dmags.append(
            pd.to_numeric(mb["mag_inst"], errors="coerce").to_numpy(float)
            - pd.to_numeric(ma["mag_inst"], errors="coerce").to_numpy(float)
        )
        snr_a.append(float(np.median(ma["snr"])))
        snr_b.append(float(np.median(mb["snr"])))
        n_pairs += int(hit.sum())
    if not dmags:
        return {}
    d = np.concatenate(dmags)
    med = float(np.median(d))
    return {
        "n_common_frames": len(common),
        "n_star_pairs": n_pairs,
        "dmag_med": med,
        "dmag_mad": float(1.4826 * np.median(np.abs(d - med))),
        "snr_med_nocr": float(np.median(snr_a)),
        "snr_med_cr": float(np.median(snr_b)),
    }


def run(target: str) -> None:
    base = REPROCESS / target
    res_nocr, res_cr = base / "result_nocr", base / "result"
    print(f"\n{'='*72}\n{target}  —  우주선 제거 전(result_nocr) vs 후(result)\n{'='*72}")
    if not res_nocr.exists():
        print(f"  기준선 없음: {res_nocr}")
        return
    if not res_cr.exists():
        print(f"  재처리 결과 없음: {res_cr}")
        return

    rows = [
        ("검출", _detect_stats(res_nocr), _detect_stats(res_cr),
         [("n_frames", "프레임 수", 0, False),
          ("n_sources_med", "프레임당 소스(중앙값)", 0, False),
          ("spike_pct", "점광원 잡음 FWHM<=1.5px", 1, True)]),
        ("마스터", _master_stats(res_nocr), _master_stats(res_cr),
         [("n_master", "마스터 소스", 0, False),
          ("n_1frame", "1프레임 검출 (D-004)", 0, False),
          ("pct_1frame", "1프레임 비율", 1, True),
          ("n_allframe", "전프레임 검출", 0, False)]),
        ("PSF", _psf_stats(res_nocr), _psf_stats(res_cr),
         [("n_psf_frames", "PSF 프레임", 0, False),
          ("flux_ratio_med", "플럭스비 psf/ap", 3, False),
          ("mad_med", "MAD (mag)", 4, False),
          ("frame_rms", "프레임간 rms (mag)", 4, False)]),
        ("CMD", _cmd_stats(res_nocr), _cmd_stats(res_cr),
         [("n_cmd_rows", "CMD 행", 0, False),
          ("n_cmd_complete", "전 밴드 유효", 0, False)]),
    ]
    for sect, a, b, fields in rows:
        if not a and not b:
            continue
        print(f"\n[{sect}]")
        print(f"  {'항목':28s} {'CR 전':>14s} {'CR 후':>14s} {'변화':>12s}")
        for key, lab, nd, pct in fields:
            va, vb = a.get(key), b.get(key)
            if va is None and vb is None:
                continue
            if isinstance(va, (int, float)) and isinstance(vb, (int, float)) \
                    and np.isfinite(va) and np.isfinite(vb):
                diff = vb - va
                ds = f"{diff:+.{nd}f}{'%p' if pct else ''}"
            else:
                ds = "—"
            print(f"  {lab:28s} {_fmt(va, nd, pct):>14s} {_fmt(vb, nd, pct):>14s} {ds:>12s}")

    ph = compare_photometry(res_nocr, res_cr)
    if ph:
        print(f"\n[측광 — 같은 별 직접 대조]")
        print(f"  공통 프레임 {ph['n_common_frames']}개 · 별-프레임 쌍 {ph['n_star_pairs']:,}개 (양쪽 SNR>20)")
        print(f"  Δmag(CR후 − CR전)  중앙값 {ph['dmag_med']:+.4f}  MAD {ph['dmag_mad']:.4f}")
        print(f"  SNR 중앙값         {ph['snr_med_nocr']:.1f} → {ph['snr_med_cr']:.1f}")


if __name__ == "__main__":
    for t in (sys.argv[1:] or ["M13"]):
        run(t)
