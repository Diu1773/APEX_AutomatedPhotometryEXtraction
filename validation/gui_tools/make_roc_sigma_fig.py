"""검출 임계 sigma 를 스캔해 completeness–purity(천문판 ROC) 곡선을 그린다.

임계를 바꿔가며 「얼마나 놓치지 않는가(completeness)」와 「잡은 것 중 진짜가
얼마인가(purity)」를 동시에 보는 표준 방법이다. 진실은 Gaia DR3 매칭으로 잡는다
— 검출이 Gaia 별과 match_radius 안에서 짝지어지면 TP, 아니면 FP 후보다.

    python make_roc_sigma_fig.py <fits> <gaia.ecsv>

주의: Gaia 보다 어두운 진짜 별은 FP 로 세어지므로 purity 는 **하한**이다.
그래도 임계에 따른 상대 변화는 그대로 읽힌다.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

import numpy as np
import sep
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from scipy.spatial import cKDTree

for _name in ("Malgun Gothic", "NanumGothic", "Gulim"):
    if any(f.name == _name for f in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = _name
        break
plt.rcParams["axes.unicode_minus"] = False

OUT = Path(__file__).absolute().parent / "fig_roc_sigma.png"
SIGMAS = (5.0, 4.0, 3.2, 2.8, 2.5, 2.2, 2.0, 1.8, 1.5, 1.2, 1.0)
MATCH_ARCSEC = 2.0


def main() -> int:
    fits_path, gaia_path = sys.argv[1], sys.argv[2]
    hdr = fits.getheader(fits_path)
    data = fits.getdata(fits_path).astype(np.float32)
    wcs = WCS(hdr).celestial
    ny, nx = data.shape

    bkg = sep.Background(data, bw=61, bh=61, fw=3, fh=3)
    sub = data - bkg.back()
    rms = float(bkg.globalrms)

    cat = Table.read(gaia_path)
    ra = np.asarray(cat["ra"], float)
    dec = np.asarray(cat["dec"], float)
    gx, gy = wcs.world_to_pixel_values(ra, dec)
    # 실제 시야 안의 Gaia 만 진실로 삼는다 (카탈로그는 radius_fudge 로 더 넓다)
    inside = (gx > 0) & (gx < nx) & (gy > 0) & (gy < ny)
    gx, gy = gx[inside], gy[inside]
    n_truth = int(inside.sum())

    # CD 행렬 WCS 는 cdelt 가 1.0 이라 그대로 쓰면 스케일이 3600″/px 로 나온다.
    scale = float(np.mean(proj_plane_pixel_scales(wcs))) * 3600.0
    tol_px = MATCH_ARCSEC / scale
    print(f"픽셀 스케일 {scale:.3f}\"/px · 매칭 {tol_px:.2f} px · "
          f"시야 내 Gaia {n_truth}개 (카탈로그 {len(ra)}개)")
    tree = cKDTree(np.c_[gx, gy])

    rows = []
    for s in SIGMAS:
        try:
            obj = sep.extract(sub, s, err=rms, minarea=3,
                              deblend_nthresh=64, deblend_cont=0.004)
        except Exception:
            continue
        dx, dy = np.asarray(obj["x"], float), np.asarray(obj["y"], float)
        d, idx = tree.query(np.c_[dx, dy], distance_upper_bound=tol_px)
        matched = np.isfinite(d)
        tp = int(len(set(idx[matched].tolist())))     # 중복 매칭은 한 번만
        n_det = int(len(dx))
        fp = n_det - int(matched.sum())
        comp = tp / n_truth if n_truth else np.nan
        pur = int(matched.sum()) / n_det if n_det else np.nan
        f1 = 2 * comp * pur / (comp + pur) if (comp + pur) > 0 else np.nan
        rows.append((s, n_det, tp, fp, comp, pur, f1))
        print(f"sigma={s:>4.1f}  검출={n_det:>6}  TP={tp:>5}  FP={fp:>6}  "
              f"완전도={comp:.3f}  순도={pur:.3f}  F1={f1:.3f}", flush=True)

    sig = [r[0] for r in rows]
    ndet = [r[1] for r in rows]
    comp = [r[4] for r in rows]
    pur = [r[5] for r in rows]
    f1 = [r[6] for r in rows]
    best = int(np.nanargmax(f1))
    # 측광에서는 가짜 검출 하나가 가짜 광도곡선·가짜 CMD 점을 만든다.
    # F1 은 완전도와 순도를 동등하게 놓으므로 순도 하한을 따로 본다.
    PUR_MIN = 0.99
    safe = max((i for i, p in enumerate(pur) if p >= PUR_MIN), default=0)

    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.6))

    ax = axes[0]
    ax.plot(pur, comp, "o-", color="#333", zorder=2)
    for s, p, c in zip(sig, pur, comp):
        ax.annotate(f"{s:g}", (p, c), textcoords="offset points", xytext=(6, -3), fontsize=8)
    ax.axvspan(PUR_MIN, 1.005, color="#2ca02c", alpha=0.12, zorder=0)
    ax.plot(pur[best], comp[best], "o", ms=14, mfc="none", mec="#d62728", mew=2,
            label=f"F1 최대 — sigma {sig[best]:g}", zorder=3)
    ax.plot(pur[safe], comp[safe], "s", ms=14, mfc="none", mec="#2ca02c", mew=2,
            label=f"순도 {PUR_MIN:.0%} 하한 — sigma {sig[safe]:g}", zorder=3)
    ax.set_xlabel("순도 purity  =  Gaia 와 짝지어진 검출 / 전체 검출")
    ax.set_ylabel("완전도 completeness  =  찾은 Gaia 별 / 시야 안 Gaia 별")
    ax.set_title("completeness–purity 곡선 (천문판 ROC)  ·  점 옆 숫자가 sigma", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc="upper left")

    ax = axes[1]
    ax.plot(sig, comp, "o-", label="완전도 (찾은 진짜 별)", color="#1f77b4")
    ax.plot(sig, pur, "o-", label="순도 (검출 중 진짜 비율)", color="#d62728")
    ax.plot(sig, f1, "o-", label="F1 (둘의 조화평균)", color="black", lw=2)
    ax.axhline(PUR_MIN, color="#2ca02c", ls=":", lw=1.2)
    ax.axvspan(min(sig), sig[safe], color="#d62728", alpha=0.08)
    ax.axvline(sig[safe], color="#2ca02c", ls="--", lw=1.4)
    ax.text(sig[safe] + 0.06, 0.55, f"순도 {PUR_MIN:.0%} 유지 한계\nsigma {sig[safe]:g}",
            color="#2ca02c", fontsize=9, ha="left")
    ax.text(1.02, 0.14, f"오염 구간\nsigma 1.2 → 검출 {ndet[-2]:,}개 중\n"
                        f"가짜 {rows[-2][3]:,}개 (33%)",
            color="#d62728", fontsize=8.5, ha="left", va="center")
    ax.invert_xaxis()
    ax.set_xlabel("검출 임계 sigma  (오른쪽이 APEX 기본값 3.2, 왼쪽으로 갈수록 낮춤)")
    ax.set_ylabel("값")
    ax.set_ylim(0, 1.05)
    ax.set_title("임계를 바꿔가며 본 세 지표", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9, loc="center left")

    fig.suptitle(
        "검출 임계 sigma 최적화 — completeness–purity(ROC) 스캔\n"
        f"실측: M13 0004-R (Moravian C3-61000, R 60 s) · 진실 = Gaia DR3 시야 내 "
        f"{n_truth}개 · 매칭 {MATCH_ARCSEC}″\n"
        "주의: Gaia 보다 어두운 진짜 별이 FP 로 세어지므로 순도는 하한이다",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    fig.savefig(OUT, dpi=150)
    print(f"[saved] {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
