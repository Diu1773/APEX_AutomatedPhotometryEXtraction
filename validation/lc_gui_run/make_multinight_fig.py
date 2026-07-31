"""다중 밤 + 다색 LC 결과 그림 — 주기도그램(별칭 구조) + 위상 접기.

    python make_multinight_fig.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

# 한글 라벨이 두부(□)로 깨지지 않게 시스템 한글 폰트를 쓴다.
for _name in ("Malgun Gothic", "NanumGothic", "Gulim", "Batang"):
    if any(f.name == _name for f in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = _name
        break
plt.rcParams["axes.unicode_minus"] = False   # 한글 폰트는 유니코드 마이너스가 없다
import numpy as np
import pandas as pd
from astropy.timeseries import LombScargle

LIT = 0.104092          # YZ Bootis 문헌 주기 (일)
LC = Path("E:/APEX_validation/reprocess/YZBoo_2n/result/lc_lightcurve/lightcurve_ID233_raw.csv")
OUT = Path(__file__).absolute().parent / "fig_multinight_period.png"

COLORS = {"g": "#2ca02c", "r": "#d62728", "i": "#ff7f0e"}


def main() -> int:
    d = pd.read_csv(LC).dropna(subset=["JD", "diff_mag"]).copy()
    # 필터마다 영점이 다르므로 각 필터를 중앙값 0 으로 맞춰서 합친다
    d["norm"] = d["diff_mag"] - d.groupby("filter")["diff_mag"].transform("median")
    t, m = d["JD"].to_numpy(), d["norm"].to_numpy()
    T = float(t.max() - t.min())

    freq = np.linspace(1 / 0.5, 1 / 0.03, 300_000)
    power = LombScargle(t, m).power(freq)
    best_f = freq[int(np.argmax(power))]
    best_P = 1.0 / best_f

    fig, axes = plt.subplots(3, 1, figsize=(9.2, 10.4))

    # (1) 원 광곡선 — 두 밤
    ax = axes[0]
    for filt, grp in d.groupby("filter"):
        ax.errorbar(grp["JD"] - t.min(), grp["norm"], yerr=grp.get("diff_err"),
                    fmt="o", ms=2.6, lw=0, elinewidth=0.5, alpha=0.75,
                    color=COLORS.get(str(filt), "#666"), label=f"{filt} (n={len(grp)})")
    # 한글 폰트에는 유니코드 마이너스(U+2212)가 없어 두부로 깨진다 — ASCII 하이픈을 쓴다
    ax.set_xlabel(f"JD - {t.min():.4f}  [일]")
    ax.set_ylabel("차등 등급 (필터별 중앙값 0)")
    ax.invert_yaxis()
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)
    ax.set_title(f"YZ Bootis — 2밤 · 364프레임 · 기저선 {T:.3f} 일 ({T/LIT:.1f} 주기)",
                 fontsize=10)
    ax.grid(alpha=0.25)

    # (2) 주기도그램 — 별칭 구조
    ax = axes[1]
    ax.plot(freq, power, lw=0.7, color="#333")
    ax.axvline(1 / LIT, color="#1f77b4", ls="--", lw=1.2,
               label=f"문헌 {1/LIT:.3f} c/d (P={LIT} d)")
    for k in (-2, -1, 1, 2):
        ax.axvline(1 / LIT + k, color="#1f77b4", ls=":", lw=0.8, alpha=0.55)
    ax.axvline(best_f, color="#d62728", ls="-", lw=1.0, alpha=0.8,
               label=f"최적 {best_f:.3f} c/d (P={best_P:.6f} d, {(best_P-LIT)/LIT*100:+.2f}%)")
    ax.set_xlim(6, 14)
    ax.set_xlabel("주파수 [c/일]   — 점선 = 1일 간격이 만드는 별칭")
    ax.set_ylabel("Lomb–Scargle power")
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)
    ax.grid(alpha=0.25)

    # (3) 위상 접기
    ax = axes[2]
    ph = ((t - t.min()) / best_P) % 1.0
    for filt, grp in d.groupby("filter"):
        sel = d["filter"] == filt
        ax.plot(np.concatenate([ph[sel], ph[sel] + 1]),
                np.concatenate([m[sel], m[sel]]), "o", ms=2.6, alpha=0.7,
                color=COLORS.get(str(filt), "#666"), label=str(filt))
    nb = 25
    edges = np.linspace(0, 1, nb + 1)
    idx = np.digitize(ph, edges) - 1
    med = np.array([np.median(m[idx == b]) if (idx == b).sum() else np.nan for b in range(nb)])
    ctr = 0.5 * (edges[:-1] + edges[1:])
    ax.plot(np.concatenate([ctr, ctr + 1]), np.concatenate([med, med]),
            "-", color="black", lw=1.6, label="위상별 중앙값")
    ax.set_xlabel(f"위상 (P = {best_P:.6f} 일)")
    ax.set_ylabel("차등 등급")
    ax.invert_yaxis()
    ax.legend(fontsize=8, ncol=4, loc="upper right", framealpha=0.9)
    ax.grid(alpha=0.25)

    fig.suptitle(
        "APEX LC 모드 — 다중 밤·다색 주기 분석\n"
        "실측: Moravian C3-61000, 0.393″/px, g·r·i 각 30 s, 2025-04-29~30 (2밤)",
        fontsize=10.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.955))
    fig.savefig(OUT, dpi=150)
    print(f"[saved] {OUT}")
    print(f"기저선 {T:.4f} 일 | 최적 P={best_P:.6f} ({(best_P-LIT)/LIT*100:+.2f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
