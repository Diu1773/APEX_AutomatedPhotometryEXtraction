# -*- coding: utf-8 -*-
"""광시야에서 단일 EPSF 가 시야 전체를 커버하는지 측정한다.

APEX 의 Step 8 은 **프레임당 EPSF 하나**를 만든다(`model_mode=per_frame`).
`psf_policy.grid_size` 는 참조별을 시야에 고르게 뽑기 위한 것이지
(`select_spatially_balanced`) 격자별 PSF 모델이 아니다. 좁은 시야에서는
PSF 공간변화가 작아 문제가 없지만, 광시야에서는 가장자리 코마·비점수차로
PSF 가 변하므로 단일 모델의 유효 범위를 실측해야 한다.

세 가지를 시야 위치의 함수로 잰다:
  1. 별 자체의 모양 변화 — Step 4 검출의 FWHM·타원율 (PSF 가 실제로 변하는가)
  2. PSF 적합 품질 — qfit/noise, reduced chi2 (단일 EPSF 가 안 맞는 곳이 있는가)
  3. PSF - 구경 등급차 — 위치에 따라 계통 오차가 생기는가

실행:
    .venv-deploy\\Scripts\\python -X utf8 \\
        validation/psf_crossinstrument/analyze_wide_epsf.py m45_wide [m67_lco]
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE = Path(r"E:\APEX_validation\psf_crossinstrument")

LABEL = {
    "m45_wide": ("M45 · LCO 0.4m coj/tfn · QHY600 full-frame 2.0°×1.3°",
                 "B/V 120 s · 0.74\"/px · 2025-01-14 · archive.lco.global"),
    "m67_lco": ("M67 · LCO 0.4m tfn · QHY600 central 30'×30'",
                "B/V 60 s · 0.74\"/px · 2024-10-29 · archive.lco.global"),
    "sinistro": ("NGC 5985 · LCO 1m elp · Sinistro 26'×26'",
                 "rp 120 s · 0.389\"/px · 2026-07-07 · archive.lco.global"),
}

NBIN = 5  # 시야를 NBIN×NBIN 격자로 나눈다


def _load(key: str) -> tuple[pd.DataFrame, dict]:
    """PSF·구경 표를 프레임별로 합치고 시야 중심으로부터의 상대 반경을 붙인다."""
    R = BASE / key / "result"
    psf_dir, ap_dir = R / "cmd_psf", R / "step7_forced_phot"
    if not psf_dir.exists():
        raise SystemExit(f"{key}: cmd_psf 없음 — Step 8 을 먼저 실행할 것")

    rows = []
    shape = None
    for t in sorted(psf_dir.glob("photometry_*.tsv")):
        p = pd.read_csv(t, sep="\t").dropna(subset=["det_uid"])
        a_path = ap_dir / t.name
        if not a_path.exists():
            continue
        a = pd.read_csv(a_path, sep="\t").dropna(subset=["det_uid"])
        m = p.merge(a, on="det_uid", suffixes=("_p", "_a"))
        m["frame"] = t.name
        rows.append(m)
    if not rows:
        raise SystemExit(f"{key}: 병합 가능한 프레임이 없음")
    df = pd.concat(rows, ignore_index=True)

    # 프레임 크기 — sci FITS 헤더에서
    from astropy.io import fits

    sci = sorted((BASE / key / "sci").glob("*.fit"))
    with fits.open(sci[0]) as h:
        shape = h[0].data.shape  # (ny, nx)
    ny, nx = shape

    x = df["x_fit_p"] if "x_fit_p" in df.columns else df["x_fit"]
    y = df["y_fit_p"] if "y_fit_p" in df.columns else df["y_fit"]
    df["_x"], df["_y"] = x.astype(float), y.astype(float)
    # 정규화 반경: 중심 0, 모서리 1
    df["_r"] = np.hypot(
        (df["_x"] - nx / 2) / (nx / 2), (df["_y"] - ny / 2) / (ny / 2)
    ) / np.sqrt(2)
    meta = {"nx": nx, "ny": ny, "n_frames": len(rows)}
    return df, meta


def _cell_stat(df: pd.DataFrame, col: str, nx: int, ny: int, nbin: int = NBIN,
               min_n: int = 8) -> np.ndarray:
    """격자 셀별 중앙값 (부족한 셀은 NaN)."""
    xi = np.clip((df["_x"] / nx * nbin).astype(int), 0, nbin - 1)
    yi = np.clip((df["_y"] / ny * nbin).astype(int), 0, nbin - 1)
    out = np.full((nbin, nbin), np.nan)
    v = pd.to_numeric(df[col], errors="coerce").to_numpy(float)
    for j in range(nbin):
        for i in range(nbin):
            sel = (xi == i) & (yi == j) & np.isfinite(v)
            if sel.sum() >= min_n:
                out[j, i] = float(np.median(v[sel]))
    return out


def run(key: str) -> dict:
    df, meta = _load(key)
    nx, ny = meta["nx"], meta["ny"]
    title, spec = LABEL.get(key, (key, ""))

    good = df[
        np.isfinite(df.get("mag_psf")) & np.isfinite(df.get("mag_inst"))
        & (df.get("snr_psf", 0) > 20) & (df.get("snr", 0) > 20)
        & (df.get("saturated_psf", 0) == 0)
    ].copy()
    good["_dmag"] = good["mag_psf"] - good["mag_inst"]
    # 프레임별 중앙값을 빼서 프레임 간 영점차를 제거 — 남는 것이 위치 의존성
    good["_dmag_rel"] = good["_dmag"] - good.groupby("frame")["_dmag"].transform("median")

    fwhm_col = next((c for c in ("fwhm_px_a", "fwhm_px", "fwhm") if c in df.columns), None)
    elong_col = next((c for c in ("elongation_a", "elongation", "elong") if c in df.columns), None)

    panels = [
        ("qfit_noise_ratio", "PSF fit quality  (qfit / expected noise)", "viridis", None),
        ("_dmag_rel", "PSF − aperture, frame median removed [mag]", "coolwarm", (-0.06, 0.06)),
    ]
    if fwhm_col:
        panels.insert(0, (fwhm_col, "star FWHM [px]  (Step 4 detection)", "magma", None))
    if elong_col:
        panels.insert(1, (elong_col, "star elongation  (Step 4 detection)", "magma", None))

    fig, axes = plt.subplots(1, len(panels) + 1, figsize=(4.1 * (len(panels) + 1), 4.3))
    for ax, (col, lab, cmap, vlim) in zip(axes, panels):
        src = good if col in good.columns else df
        grid = _cell_stat(src, col, nx, ny)
        vmin, vmax = vlim if vlim else (np.nanpercentile(grid, 5), np.nanpercentile(grid, 95))
        im = ax.imshow(grid, origin="lower", cmap=cmap, vmin=vmin, vmax=vmax,
                       extent=[0, nx, 0, ny], aspect="auto")
        ax.set_title(lab, fontsize=9)
        ax.set_xticks([]); ax.set_yticks([])
        fig.colorbar(im, ax=ax, fraction=0.046)

    # 마지막 패널: 반경 의존성 (데이터=점, 구간 중앙값=선)
    ax = axes[-1]
    ax.plot(good["_r"], good["_dmag_rel"], ".", ms=2, alpha=0.25, color="#1f77b4",
            label=f"stars SNR>20 (N={len(good)})")
    edges = np.linspace(0, good["_r"].max(), 9)
    cen, med, mad = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = good[(good["_r"] >= lo) & (good["_r"] < hi)]
        if len(s) >= 20:
            cen.append(0.5 * (lo + hi))
            m = float(np.median(s["_dmag_rel"]))
            med.append(m)
            mad.append(float(1.4826 * np.median(np.abs(s["_dmag_rel"] - m))))
    ax.plot(cen, med, "-o", color="crimson", ms=4, lw=1.4, label="median per radial bin")
    ax.fill_between(cen, np.array(med) - np.array(mad), np.array(med) + np.array(mad),
                    color="crimson", alpha=0.15, label="±MAD")
    ax.axhline(0, color="0.4", lw=0.7)
    ax.set_xlabel("normalised field radius  (0 = centre, 1 = corner)")
    ax.set_ylabel("PSF − aperture, frame median removed [mag]")
    ax.set_ylim(-0.15, 0.15)
    ax.legend(fontsize=7.5)
    ax.set_title("radial trend", fontsize=9)

    fig.suptitle(f"Single-EPSF coverage across the field — {title}", fontsize=11)
    fig.text(0.01, 0.01,
             f"{spec} | APEX Step-8 headless, model_mode=per_frame (one EPSF per frame) | "
             f"{meta['n_frames']} frames · {len(good)} star measurements SNR>20 | "
             f"grid {NBIN}×{NBIN} cells, median per cell",
             fontsize=7.5, color="0.35")
    fig.tight_layout(rect=[0, 0.045, 1, 0.93])
    out = BASE / f"fig_wide_epsf_{key}.png"
    fig.savefig(out, dpi=140)
    plt.close(fig)

    # 정량 요약: 중심부(r<0.3) vs 가장자리(r>0.7)
    inner = good[good["_r"] < 0.3]
    outer = good[good["_r"] > 0.7]
    summ = {
        "key": key, "n_frames": meta["n_frames"], "n_stars": len(good),
        "shape": f"{nx}x{ny}",
        "dmag_inner_med": float(np.median(inner["_dmag_rel"])) if len(inner) else np.nan,
        "dmag_outer_med": float(np.median(outer["_dmag_rel"])) if len(outer) else np.nan,
        "mad_inner": float(1.4826 * np.median(np.abs(inner["_dmag_rel"] - np.median(inner["_dmag_rel"])))) if len(inner) else np.nan,
        "mad_outer": float(1.4826 * np.median(np.abs(outer["_dmag_rel"] - np.median(outer["_dmag_rel"])))) if len(outer) else np.nan,
        "qfit_inner": float(np.median(inner["qfit_noise_ratio"])) if len(inner) else np.nan,
        "qfit_outer": float(np.median(outer["qfit_noise_ratio"])) if len(outer) else np.nan,
        "n_inner": len(inner), "n_outer": len(outer),
    }
    if fwhm_col:
        src = df
        summ["fwhm_inner"] = float(np.nanmedian(src.loc[src["_r"] < 0.3, fwhm_col]))
        summ["fwhm_outer"] = float(np.nanmedian(src.loc[src["_r"] > 0.7, fwhm_col]))
    if elong_col:
        src = df
        summ["elong_inner"] = float(np.nanmedian(src.loc[src["_r"] < 0.3, elong_col]))
        summ["elong_outer"] = float(np.nanmedian(src.loc[src["_r"] > 0.7, elong_col]))

    print(f"\n=== {key} ({nx}x{ny}, {meta['n_frames']} frames) ===")
    print(f"  중심 r<0.3 (N={summ['n_inner']}) vs 가장자리 r>0.7 (N={summ['n_outer']})")
    print(f"  Δmag(프레임중앙값 제거) : {summ['dmag_inner_med']:+.4f} → {summ['dmag_outer_med']:+.4f}"
          f"  (MAD {summ['mad_inner']:.4f} → {summ['mad_outer']:.4f})")
    print(f"  qfit/noise             : {summ['qfit_inner']:.3f} → {summ['qfit_outer']:.3f}")
    if "fwhm_inner" in summ:
        print(f"  별 FWHM [px]           : {summ['fwhm_inner']:.2f} → {summ['fwhm_outer']:.2f}")
    if "elong_inner" in summ:
        print(f"  별 elongation          : {summ['elong_inner']:.3f} → {summ['elong_outer']:.3f}")
    print(f"  fig -> {out}")
    return summ


if __name__ == "__main__":
    keys = sys.argv[1:] or ["m45_wide"]
    out = [run(k) for k in keys]
    if len(out) > 1:
        print("\n요약 (시야가 넓을수록 가장자리 열화가 커지는지 비교)")
        for s in out:
            print(f"  {s['key']:10s} {s['shape']:>11s}  Δmag {s['dmag_inner_med']:+.4f}→{s['dmag_outer_med']:+.4f}"
                  f"  qfit {s['qfit_inner']:.2f}→{s['qfit_outer']:.2f}")
