"""Step 4 검출 QC — 음수 검출 순도 추정이 대상 종류를 가리지 않는가.

sigma_qc_scan.py 가 남긴 JSON 을 읽어 두 가지를 본다.
  왼쪽  대상 4개(구상 2·산개 2)의 순도-임계 곡선. 안전 하한은 대상마다 다르지만
        무너지는 지점은 공통이다.
  오른쪽 음수 검출로 낸 추정이 Gaia 실측과 맞는가(대각선에 붙을수록 정확).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

import numpy as np

for _name in ("Malgun Gothic", "NanumGothic", "Gulim"):
    if any(f.name == _name for f in font_manager.fontManager.ttflist):
        plt.rcParams["font.family"] = _name
        break
plt.rcParams["axes.unicode_minus"] = False

HERE = Path(__file__).absolute().parent
SRC = HERE / "sigma_qc"
OUT = HERE / "fig_sigma_qc.png"

# 대상, 표시이름, 색.  M13 은 같은 성단에서 필터만 바꾼 두 프레임을 함께 본다
# — 하한이 대상 사이뿐 아니라 한 대상 안에서도 움직인다는 걸 보이기 위해서.
TARGETS = [
    ("M13_R", "M13  구상·혼잡", "#d62728"),
    ("M13", "M13  구상·혼잡", "#d62728"),
    ("M3", "M3  구상", "#ff7f0e"),
    ("M67", "M67  산개", "#1f77b4"),
    ("NGC6811", "NGC6811  산개·성김", "#2ca02c"),
]
STYLE = {"M13": dict(ls="--", ms=4.0), "M13_R": dict(ls="-", ms=4.5)}
# 순도 하한. 0.99 는 검출이 적은 프레임(M3 B 는 298개)에서 음수 몇 개만으로
# 깨져 쓸 수 없다. 0.95 는 네 대상 모두에서 기본값 sigma 를 통과시키면서
# sigma 1.2 오염은 전부 잡는다.
PUR_MIN = 0.95


def main() -> int:
    data = {}
    for key, _, _ in TARGETS:
        path = SRC / f"{key}.json"
        if path.exists() and path.stat().st_size > 2:
            data[key] = json.loads(path.read_text(encoding="utf-8"))
    if not data:
        raise SystemExit(f"결과 JSON 이 없다: {SRC}")

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.6))

    # ── 왼쪽: 순도-임계 곡선 ──────────────────────────────────────────
    ax = axes[0]
    ax.axvspan(1.0, 1.35, color="#d62728", alpha=0.08)
    ax.text(1.33, 0.30, "무너지는 구간\n(대상 공통)", color="#d62728",
            fontsize=9, ha="right")
    ax.axhline(PUR_MIN, color="#555", ls="--", lw=1.3)
    ax.text(4.95, PUR_MIN - 0.035, f"순도 {PUR_MIN:.0%} — 제안 게이트",
            fontsize=8.5, color="#555")

    for key, label, color in TARGETS:
        d = data.get(key)
        if not d:
            continue
        sig = [r["sigma"] for r in d["rows"]]
        pur = [r["purity_est"] for r in d["rows"]]
        st = STYLE.get(key, {})
        ax.plot(sig, pur, "o", ls=st.get("ls", "-"), color=color,
                label=f"{label}  ({d['filter']})", ms=st.get("ms", 4.5))
        # 순도 게이트를 유지하는 가장 낮은 sigma
        ok = [s for s, p in zip(sig, pur) if p >= PUR_MIN]
        if ok:
            s_safe = min(ok)
            ax.plot([s_safe], [dict(zip(sig, pur))[s_safe]], "s", color=color,
                    ms=11, mfc="none", mew=2)
            ax.annotate(f"{s_safe:g}", (s_safe, dict(zip(sig, pur))[s_safe]),
                        textcoords="offset points", xytext=(2, 10),
                        fontsize=8.5, color=color, weight="bold")

    ax.invert_xaxis()
    ax.set_xlabel("검출 임계 sigma   (오른쪽이 APEX 기본값 3.2)")
    ax.set_ylabel("순도 추정  =  (양수검출 - 음수검출) / 양수검출")
    ax.set_ylim(0.2, 1.03)
    ax.set_title("안전 하한은 대상마다 다르고(□), 무너지는 지점은 같다", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8.5, loc="lower left")

    # ── 오른쪽: 추정 vs Gaia 실측 ────────────────────────────────────
    ax = axes[1]
    lim = [0.5, 2e4]
    ax.plot(lim, lim, "-", color="#888", lw=1, zorder=1)
    ax.fill_between(lim, [v * 0.5 for v in lim], [v * 2.0 for v in lim],
                    color="#888", alpha=0.12, zorder=0)
    ax.text(1.2, 3.2, "회색 띠 = 2배 이내", fontsize=8.5, color="#555", rotation=38)

    n_pt = 0
    for key, label, color in TARGETS:
        d = data.get(key)
        if not d:
            continue
        xs = [r["fp_gaia"] for r in d["rows"] if "fp_gaia" in r]
        ys = [r["n_neg"] for r in d["rows"] if "fp_gaia" in r]
        ss = [r["sigma"] for r in d["rows"] if "fp_gaia" in r]
        if not xs:
            continue
        n_pt += len(xs)
        ax.plot(xs, ys, "o" if key != "M13" else "^", color=color,
                label=f"{label.split('  ')[0]} {d['filter']}", ms=7,
                alpha=0.85, zorder=2)
        for x, y, s in zip(xs, ys, ss):
            if x >= 100:
                ax.annotate(f"{s:g}", (x, y), textcoords="offset points",
                            xytext=(7, -2), fontsize=8, color=color)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel("실측 가짜 검출 수  (Gaia 와 안 짝지어진 검출)")
    ax.set_ylabel("추정 가짜 검출 수  (음수 영상에서 검출)")
    ax.set_title(f"외부 카탈로그 없이 낸 추정이 맞는가  ·  {n_pt}점", fontsize=11)
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=9, loc="upper left")

    fig.suptitle(
        "Step 4 검출 QC — 배경 뺀 영상의 부호를 뒤집어 가짜 검출을 센다  "
        "(WCS·외부 카탈로그 없이, 검출 비용의 3.6%)\n"
        "실측: Moravian C3-61000 4788x3194 · 대상 4개 프레임 5장(60 s) · "
        "대조 진실 = Gaia DR3 시야 내 매칭 2.0\"  (M3 는 Gaia 산출물 없어 왼쪽만)",
        fontsize=9.5,
    )
    fig.tight_layout(rect=(0, 0, 1, 0.91))
    fig.savefig(OUT, dpi=150)
    print(f"[saved] {OUT}")
    for key, _, _ in TARGETS:
        d = data.get(key)
        if not d:
            continue
        rows = d["rows"]
        ok = [r["sigma"] for r in rows if r["purity_est"] >= PUR_MIN]
        print(f"  {key:<9} 순도 {PUR_MIN:.0%} 하한 sigma = "
              f"{min(ok) if ok else float('nan'):g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
