"""LC 8→11 headless against the window's saved run, on one workspace.

Both read the same Step 8 window selection (target 153, six comparisons, check
187), the same 364 frames, the same detrend mode. What the run had to show is
whether the batch chain reproduces the window's science and where it does not.
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = ["Malgun Gothic", "Segoe UI", "DejaVu Sans"]
BG, TEXT, MUTED, GRID = "#FFFFFF", "#1F2933", "#5B6573", "#D8DEE7"
WIN, HEAD, BAD = "#8A6A00", "#3A66DB", "#C73030"

REF = pathlib.Path(r"E:\APEX_validation\_yzboo_manual_ref")
NEW = pathlib.Path(r"E:\APEX_validation\reprocess\YZBoo_2n\result")
LIT = 0.104092

ref_raw = pd.read_csv(REF / "lc_lightcurve/lightcurve_ID153_raw.csv")
new_raw = pd.read_csv(NEW / "lc_lightcurve/lightcurve_ID153_raw.csv")
ref_cor = pd.read_csv(REF / "lc_detrend/lightcurve_ID153_offset.csv")
new_cor = pd.read_csv(NEW / "lc_detrend/lightcurve_ID153_offset.csv")

fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.6), dpi=110)
fig.patch.set_facecolor(BG)

# ── 1. before detrending: the ensembles differ ─────────────────────────────
ax = axes[0]
t0 = float(ref_raw["BJD_TDB"].min())
ax.scatter(ref_raw["BJD_TDB"] - t0, ref_raw["diff_mag_raw"], s=7, c=WIN,
           alpha=0.6, label="창 (합집합 11 개)", zorder=3)
ax.scatter(new_raw["BJD_TDB"] - t0, new_raw["diff_mag_raw"], s=7, c=HEAD,
           alpha=0.6, label="헤드리스 (선별 6 개)", zorder=3)
ax.invert_yaxis()
ax.set_xlabel("BJD_TDB − 2460795.09 (일)", fontsize=10)
ax.set_ylabel("diff_mag_raw (등급)", fontsize=10)
ax.set_title("추세제거 전 — 앙상블이 달라 최대 0.40 mag 어긋난다",
             fontsize=11, color=TEXT, weight="bold", loc="left", pad=9)
ax.legend(fontsize=9, frameon=False, loc="upper center", ncol=2)
ax.grid(True, color=GRID, lw=0.7, zorder=0)

# ── 2. after detrending: the difference is absorbed ────────────────────────
ax = axes[1]
ax.scatter(ref_cor["BJD_TDB"] - t0, ref_cor["diff_mag_corr"], s=7, c=WIN,
           alpha=0.6, zorder=3)
ax.scatter(new_cor["BJD_TDB"] - t0, new_cor["diff_mag_corr"], s=7, c=HEAD,
           alpha=0.6, zorder=3)
ax.invert_yaxis()
amp_r = np.percentile(ref_cor["diff_mag_corr"], 99) - np.percentile(ref_cor["diff_mag_corr"], 1)
amp_n = np.percentile(new_cor["diff_mag_corr"], 99) - np.percentile(new_cor["diff_mag_corr"], 1)
ax.set_xlabel("BJD_TDB − 2460795.09 (일)", fontsize=10)
ax.set_ylabel("diff_mag_corr (등급)", fontsize=10)
ax.set_title(f"밤·필터별 영점 보정 후 — 진폭 {amp_r:.4f} vs {amp_n:.4f} mag",
             fontsize=11, color=TEXT, weight="bold", loc="left", pad=9)
ax.text(0.02, 0.05, "6 군(밤 2 × 필터 3)에 각각 영점을 맞추면\n"
                    "앙상블 차이가 그대로 흡수된다",
        transform=ax.transAxes, fontsize=9, color=MUTED, va="bottom")
ax.grid(True, color=GRID, lw=0.7, zorder=0)

# ── 3. the periods, and the verdict that differs ───────────────────────────
ax = axes[2]
methods = ["raw_ls", "corr_ls", "raw_pdm", "corr_pdm"]
vals = [0.094468, 0.095289, 0.104295, 0.105297]      # identical in both runs
y = np.arange(len(methods))
ax.barh(y, [(v - LIT) / LIT * 100 for v in vals], color=MUTED, alpha=0.5,
        height=0.5, zorder=3)
for i, v in enumerate(vals):
    ax.text((v - LIT) / LIT * 100, i, f"  {v:.6f} d", va="center",
            fontsize=8.5, color=TEXT,
            ha="left" if v > LIT else "right")
ax.axvline(0, color=TEXT, lw=1.2, zorder=4)
ax.set_yticks(y)
ax.set_yticklabels(methods, fontsize=9)
ax.set_xlim(-12, 7.5)
# Room above the bars for the verdict box, which otherwise sits on corr_pdm.
ax.set_ylim(-0.6, 6.4)
ax.set_xlabel("문헌(0.104092 d) 대비 (%)", fontsize=10)
ax.set_title("주기 4 개는 소수 6 자리까지 같다 — 갈린 건 별칭 판정",
             fontsize=11, color=TEXT, weight="bold", loc="left", pad=9)
ax.grid(True, axis="x", color=GRID, lw=0.7, zorder=0)
ax.text(0.03, 0.97,
        "창 08-01  RESOLVED  0.104209 d\n"
        "   밤을 1 개로 봤다(night_id 전부 0, 1.10 일짜리 「하룻밤」)\n"
        "   → leave-one-night-out 시행 0 회, 공짜 통과\n\n"
        "헤드리스  AMBIGUOUS  0.104151 d\n"
        "   밤 2 개(04-29 · 04-30), 시행 2 회 중 1 회 일치 = 50 %\n"
        "   → 가장 긴 밤이 1.23 주기뿐이라 못 가린다는 판정",
        transform=ax.transAxes, fontsize=8.4, color=TEXT, va="top",
        bbox=dict(fc="#F4F6F9", ec=GRID, lw=0.8, pad=6))

note = ("실측: YZ Boo · Moravian C3-61000 · 2025-04-29/30 · 364 프레임(g 124 · r 119 · i 21)"
        " · 워크스페이스 E:\\APEX_validation\\reprocess\\YZBoo_2n\n"
        "양쪽 모두 스텝 8 창 선택(대상 ID 153, 비교성 119·166·182·199·209·226, 체크 187)에서"
        " 출발 · 창 저장본 2026-08-01 실행 · 헤드리스 2026-08-21"
        " `apex run --mode lc --steps 8-11`")
fig.text(0.5, 0.012, note, fontsize=8.3, color=MUTED, ha="center",
         linespacing=1.6)
fig.tight_layout(rect=(0, 0.085, 1, 1))
out = pathlib.Path("validation/fig_lc_headless_vs_window.png")
fig.savefig(out, dpi=110, facecolor=BG)
print(f"  {out}")
