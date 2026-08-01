"""검출 sigma 스캔 — 실측/예측 비가 꺾이는 지점이 그 프레임의 sigma 하한.

sigma 를 낮추면 완전도가 깊어져 검출이 는다. 그 증가가 **완전도 모델로 설명되는
만큼**이면 진짜 별이고, 설명을 넘어서면 잡음이다. 두 곡선의 비를 그리면 평평한
구간과 꺾이는 지점이 눈에 보인다.

수치는 sigma_scan 실측값(M13 0004-R, Moravian C3-61000, R 60 s)이다.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

for _name in ("Malgun Gothic", "NanumGothic", "Gulim"):
    if any(f.name == _name for f in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = _name
        break
plt.rcParams["axes.unicode_minus"] = False

OUT = Path(__file__).absolute().parent / "fig_sigma_scan.png"

# sigma, S/N50, 예측, 실측  — sigma_scan.py 실행 결과 (실측)
ROWS = [
    (3.2, 4.05, 3370, 1514),
    (2.8, 3.54, 3697, 1641),
    (2.5, 3.16, 4000, 1746),
    (2.2, 2.78, 4389, 1845),
    (2.0, 2.53, 4729, 1954),
    (1.8, 2.28, 5176, 2114),
    (1.5, 1.90, 6150, 2509),
    (1.2, 1.52, 7363, 4067),
    (1.0, 1.27, 7845, 10437),
]


def main() -> int:
    sig = [r[0] for r in ROWS]
    pred = [r[2] for r in ROWS]
    obs = [r[3] for r in ROWS]
    ratio = [o / p for o, p in zip(obs, pred)]

    fig, axes = plt.subplots(2, 1, figsize=(8.6, 8.0), sharex=True)

    ax = axes[0]
    ax.plot(sig, pred, "o--", color="#1f77b4", label="예측 (완전도 모델 x Gaia 카탈로그)")
    ax.plot(sig, obs, "o-", color="#d62728", label="실측 (sep.extract)")
    ax.set_yscale("log")
    ax.set_ylabel("검출 수")
    ax.invert_xaxis()
    ax.legend(fontsize=9)
    ax.grid(alpha=0.3)
    ax.set_title("검출 sigma 를 낮추면 어디서 잡음이 섞이기 시작하나", fontsize=11)

    ax = axes[1]
    ax.plot(sig, ratio, "o-", color="black")
    ax.axhspan(0.38, 0.48, color="#2ca02c", alpha=0.14)
    ax.text(3.15, 0.43, "평평한 구간 — 늘어난 검출이\n완전도 모델로 설명된다(진짜 별)",
            fontsize=8.5, va="center", color="#2ca02c")
    ax.axvline(1.35, color="#d62728", ls="--", lw=1.2)
    ax.text(1.33, 1.05, "꺾이는 지점\n≈ 이 프레임의 sigma 하한", fontsize=8.5,
            color="#d62728", ha="right")
    ax.set_xlabel("검출 임계 sigma  (오른쪽이 기본값 3.2, 왼쪽으로 갈수록 낮춤)")
    ax.set_ylabel("실측 / 예측")
    ax.invert_xaxis()
    ax.grid(alpha=0.3)

    fig.suptitle(
        "APEX — sigma 스캔으로 검출 오염을 잡는다\n"
        "실측: M13 0004-R (Moravian C3-61000, R 60 s, sky σ 22.3 e⁻/px, FWHM 4.47 px)\n"
        "예측 = APEX 완전도 erf 모델(S/N₅₀ = 4.05 @ sigma 3.2) x Gaia FOV 7,925개",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(OUT, dpi=150)
    print(f"[saved] {OUT}")
    print("비율:", " ".join(f"{s}:{r:.2f}" for s, r in zip(sig, ratio)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
