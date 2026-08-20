"""D-012 evidence, in one page.

Left: the fitted curvature per band and field. A colour term belongs to the
filter and detector, so a band's points should sit on top of each other.
Right: what happens when one field's curvature is lent to another.
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams["font.family"] = ["Malgun Gothic", "Segoe UI", "DejaVu Sans"]
BG, TEXT, MUTED, GRID = "#FFFFFF", "#1F2933", "#5B6573", "#D8DEE7"
OK, BAD, NEUTRAL = "#247A46", "#C73030", "#3A66DB"

crit = pd.read_csv("validation/color_term_quadratic_criterion.csv")
transfer = pd.read_csv("validation/color_term_quadratic_transfer.csv")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.5, 6.2), dpi=110)
fig.patch.set_facecolor(BG)

# ── left: the coefficient per (band, colour) and field ─────────────────────
# Grouped by colour index, not just band: R is fitted against V-R on the
# Johnson fields and g-r on the SDSS ones, and those are different
# transformations. Stacking them on one tick would compare two unlike things.
def label(band, pair):
    return "{}\n({})".format(band, str(pair).replace("_", "−"))

crit["group"] = [label(b, p) for b, p in zip(crit["band"], crit["pair"])]
groups = sorted(crit["group"].unique(), key=lambda g: (g.split("\n")[0], g))
xpos = {g: i for i, g in enumerate(groups)}
markers = {"M13": "o", "M3": "s", "M5": "^", "M67": "D", "NGC6811": "v"}

for cluster, group in crit.groupby("cluster"):
    ax1.scatter([xpos[g] for g in group["group"]], group["ct2"],
                s=95, marker=markers.get(cluster, "o"), label=cluster,
                edgecolor="white", linewidth=1.1, zorder=3, alpha=0.92)

ax1.axhline(0.0, color=MUTED, lw=1.0, zorder=1)
ax1.axhspan(-0.05, 0.05, color=OK, alpha=0.10, zorder=0)
ax1.text(len(groups) - 0.45, -0.30,
         "초록 띠 = 광대역 변환 곡률의\n물리적 크기 (≈ ±0.05)",
         fontsize=8.5, color=OK, va="top", ha="right")
ax1.axhline(0.25, color=BAD, lw=1.0, ls="--", zorder=1)
ax1.axhline(-0.25, color=BAD, lw=1.0, ls="--", zorder=1)
ax1.text(-0.42, 0.256, "옛 채택 상한 |ct2| ≤ 0.25", fontsize=8.5, color=BAD, va="bottom")

ax1.set_xticks(range(len(groups)))
ax1.set_xticklabels(groups, fontsize=9)
ax1.set_xlim(-0.5, len(groups) - 0.5)
ax1.set_ylim(-0.85, 0.42)
ax1.set_xlabel("밴드 (적합에 쓴 색지수)", fontsize=10.5, color=TEXT)
ax1.set_ylabel("적합된 2차 색항 ct2", fontsize=10.5, color=TEXT)
ax1.set_title("한 대의 카메라, 다섯 시야 — 같은 밴드가 다른 곡률을 낸다",
              fontsize=11.5, color=TEXT, weight="bold", loc="left", pad=10)
ax1.grid(True, color=GRID, lw=0.7, zorder=0)
ax1.legend(fontsize=9, frameon=False, loc="lower left", ncol=2)
r_vr = crit[(crit["band"] == "R") & (crit["pair"] == "V_R")]
if not r_vr.empty:
    top = r_vr.loc[r_vr["ct2"].idxmax()]
    ax1.annotate(
        "R(V−R) 세 시야: {:+.3f} ~ {:+.3f}\n같은 색지수인데 부호가 반대".format(
            r_vr["ct2"].min(), r_vr["ct2"].max()),
        xy=(xpos[label(top["band"], top["pair"])], float(top["ct2"])),
        xytext=(0.7, -0.62), fontsize=9, color=BAD,
        arrowprops=dict(arrowstyle="->", color=BAD, lw=1.2))

# ── right: does the curvature transfer? ────────────────────────────────────
# The transfer table was already built per (band, colour); align its ticks
# with the left panel so a reader can follow one group across both.
pair_of = dict(zip(crit["band"] + "|" + crit["cluster"], crit["pair"]))
transfer["pair"] = [pair_of.get(f"{b}|{t}", "") for b, t
                    in zip(transfer["band"], transfer["target"])]
transfer["group"] = [label(b, p) for b, p
                     in zip(transfer["band"], transfer["pair"])]
transfer = transfer[transfer["group"].isin(groups)]
jitter = np.random.default_rng(0).uniform(-0.13, 0.13, len(transfer))
colors = [BAD if v > 0 else OK for v in transfer["d_mmag"]]
ax2.scatter([xpos[g] for g in transfer["group"]] + jitter,
            transfer["d_mmag"], s=85, c=colors, alpha=0.85,
            edgecolor="white", linewidth=1.0, zorder=3)
ax2.axhline(0.0, color=MUTED, lw=1.2, zorder=2)
ax2.set_xticks(range(len(groups)))
ax2.set_xticklabels(groups, fontsize=9)
ax2.set_xlim(-0.5, len(groups) - 0.5)
ax2.set_xlabel("밴드 (적합에 쓴 색지수)", fontsize=10.5, color=TEXT)
ax2.set_ylabel("남의 곡률을 빌려 썼을 때 산포 변화 (mmag)\n"
               "0 = 곡률 안 쓴 것과 같음", fontsize=10.5, color=TEXT)
ax2.set_title("빌려 쓰면 나아지는가 — B 는 그렇고 R 은 아니다",
              fontsize=11.5, color=TEXT, weight="bold", loc="left", pad=10)
ax2.grid(True, color=GRID, lw=0.7, zorder=0)
ax2.text(0.02, 0.96, "위 = 안 쓰느니만 못하다", transform=ax2.transAxes,
         fontsize=9.5, color=BAD, va="top")
ax2.text(0.98, 0.04, "아래 = 진짜 기기 곡률", transform=ax2.transAxes,
         fontsize=9.5, color=OK, va="bottom", ha="right")

note = ("실측 데이터: Moravian C3-61000 · 성단 5개(M13·M3·M5·M67·NGC 6811) · 2025 시즌 · "
        "APEX Step 10 보정별 표(gaia_sdss_calibrator_by_ID.csv), SNR≥20 + Gaia 품질컷 적용 후 337~1,705개.  "
        "왼쪽은 각 시야가 자기 자료로 적합한 ct2. 오른쪽은 그 값을 다른 시야에 적용해 잰 robust 산포 변화(24쌍).")
fig.text(0.5, 0.012, note, fontsize=8.4, color=MUTED, ha="center", wrap=True)

fig.tight_layout(rect=(0, 0.055, 1, 1))
out = pathlib.Path("validation/fig_d012_color_term.png")
out.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(out, dpi=110, facecolor=BG)
print(f"  {out}")
