"""두 기기의 색-등급도를 영점 보정 전후로 겹쳐 그린다.

**무엇을 보여 주려는 그림인가.** 서로 다른 망원경·검출기·필터로 찍은 같은 성단이
APEX 를 지나면 같은 자리에 오는가. 왼쪽은 영점을 붙이기 전(기기 등급), 오른쪽은
붙인 뒤(표준 등급)다.

**같은 별만 그린다.** 두 워크스페이스에서 `gaia_source_id` 로 이어진 별만 쓴다.
서로 다른 별을 그려 놓고 「겹친다」고 하면 그림이 거짓말을 한다.

**기기 등급은 APEX 가 낸 표에서 그대로 읽는다** (`median_by_ID_filter_wide_raw.csv`).

처음 그릴 때는 MuSCAT3 쪽 원본 표가 갤럭시북에 있고 그 기계가 꺼져 있어서 영점
모형을 거꾸로 풀어 되돌렸다(`mag_inst = mag_cal − zp − ct·색 − ct2·색²`).
2026-09-08 에 원본 표를 가져와 **되돌린 값을 실제 값으로 바꿨다.** 되돌린 값은
실제와 1~6 밀리등급 안에서 맞았지만, 되돌리는 식에 그때의 `ct2` 가 들어가 있어서
**그림이 영점 계수의 선택에 딸려 움직이는 것이 문제였다.** 이제 안 그렇다.

실행:
    python -X utf8 validation/external_muscat3/fig_cmd_overlay.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).absolute().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from apex.utils.io_utils import read_csv_int64_source_id  # noqa: E402
from compare_magnitudes import join, load_workspace  # noqa: E402

M3 = HERE / "results"
MOR = Path("E:/APEX_validation/reprocess/M67/result")

STYLE = {
    "MuSCAT3": dict(color="#1f4e79", marker="o", label="MuSCAT3  ·  LCO 2 m (ogg)"),
    "Moravian": dict(color="#a03020", marker="^", label="Moravian C3-61000  ·  사용자 망원경"),
}


def instrumental(cal: pd.DataFrame, result_dir: Path, label: str) -> pd.DataFrame:
    """APEX 가 낸 기기 등급 표를 `ID` 로 붙인다."""
    p = result_dir / "cmd_zeropoint" / "median_by_ID_filter_wide_raw.csv"
    if not p.exists():
        raise SystemExit(f"{label}: 기기 등급 표가 없다 — {p}")
    raw = read_csv_int64_source_id(p)
    cols = ["ID"] + [f"mag_inst_{b}" for b in ("g", "r", "i")
                     if f"mag_inst_{b}" in raw.columns]
    out = cal.merge(raw[cols], on="ID", how="left")
    n = int(out["mag_inst_i"].notna().sum()) if "mag_inst_i" in out else 0
    print(f"[{label}] 기기 등급이 붙은 별 {n} / {len(out)}")
    return out


def main() -> int:
    m3 = instrumental(load_workspace(M3, "MuSCAT3"), M3, "MuSCAT3")
    mv = instrumental(load_workspace(MOR, "Moravian"), MOR, "Moravian")
    j = join(m3, mv)
    num = lambda c: pd.to_numeric(j[c], errors="coerce").to_numpy(float)

    sets = {}
    for tag, suf in (("MuSCAT3", "_A"), ("Moravian", "_B")):
        sets[tag] = {
            "inst": (num(f"mag_inst_g{suf}") - num(f"mag_inst_i{suf}"),
                     num(f"mag_inst_i{suf}")),
            "cal": (num(f"mag_cal_g{suf}") - num(f"mag_cal_i{suf}"),
                    num(f"mag_cal_i{suf}")),
        }

    n = int(np.isfinite(sets["MuSCAT3"]["cal"][1]
                        + sets["Moravian"]["cal"][1]).sum())

    plt.rcParams.update({
        "font.family": ["Malgun Gothic", "DejaVu Sans"],
        "axes.unicode_minus": False, "font.size": 10,
        "axes.grid": True, "grid.alpha": 0.25, "grid.linewidth": 0.6,
        "axes.edgecolor": "#444444", "axes.linewidth": 0.9,
    })
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 6.6))

    for ax, key, title in (
        (axes[0], "inst", "영점 보정 전  —  기기 등급"),
        (axes[1], "cal", "영점 보정 후  —  표준 등급"),
    ):
        for tag, st in STYLE.items():
            x, y = sets[tag][key]
            m = np.isfinite(x) & np.isfinite(y)
            ax.scatter(x[m], y[m], s=15, alpha=0.55, linewidths=0,
                       color=st["color"], marker=st["marker"], label=st["label"])
        ax.set_xlabel("g − i" if key == "cal" else "g − i  (기기)")
        ax.set_title(title, fontsize=11.5, pad=8)
        ax.invert_yaxis()
    axes[0].set_ylabel("i  (기기 등급)")
    axes[1].set_ylabel("i  (표준 등급)")

    # **별마다의 차이**를 적는다. 두 무리의 중앙값을 빼는 것과는 다른 값이다 —
    # 분포 모양이 다르면 둘이 어긋나고(여기서 0.067 대 −0.015), 우리가 재는 것은
    # 앞쪽이다. 분석 문서의 수치도 그것이다.
    def pair_delta(key):
        (xa, ya), (xb, yb) = sets["MuSCAT3"][key], sets["Moravian"][key]
        m = np.isfinite(ya) & np.isfinite(yb) & np.isfinite(xa) & np.isfinite(xb)
        dy = (ya - yb)[m]
        return (float(np.median(dy)), float(np.median((xa - xb)[m])),
                float(1.4826 * np.median(np.abs(dy - np.median(dy)))))

    for ax, key, fc, ec in ((axes[0], "inst", "#f4f4f2", "#bbbbbb"),
                            (axes[1], "cal", "#eef3f7", "#88a8c0")):
        dmag, dcol, sd = pair_delta(key)
        fmt = ".2f" if key == "inst" else ".3f"
        ax.text(0.03, 0.03,
                "별마다의 차이 (MuSCAT3 − Moravian)\n"
                f"등급 {dmag:+{fmt}}   ·   색 {dcol:+{fmt}}   ·   흩어짐 {sd:.3f}",
                transform=ax.transAxes, fontsize=9, va="bottom",
                bbox=dict(boxstyle="round,pad=0.4", fc=fc, ec=ec))
        print(f"{key}: 등급 {dmag:+.4f} · 색 {dcol:+.4f} · MAD {sd:.4f}")

    axes[1].legend(loc="upper left", framealpha=0.92, fontsize=9.5)
    fig.suptitle(f"M67 — 같은 별 {n} 개를 두 기기로 잰 색-등급도", fontsize=13.5, y=0.975)
    fig.text(0.5, 0.040,
             "자료 : LCO MuSCAT3 (ogg 2m0a · gp·rp·ip · 2021-03-17 · 밴드당 60 장)"
             "     ·     Moravian C3-61000 (사용자 망원경 · g·r·i)"
             "     ·     둘 다 APEX 로 원본부터 처리",
             ha="center", fontsize=8.4, color="#444444")
    fig.text(0.5, 0.013,
             "별은 gaia_source_id 로 이었다     ·     양쪽 값 모두 APEX 산출물에서"
             " 그대로 읽었다 (기기 등급 median_by_ID_filter_wide_raw.csv,"
             " 표준 등급 median_by_ID_filter_wide.csv)",
             ha="center", fontsize=8.4, color="#666666")
    fig.tight_layout(rect=(0, 0.065, 1, 0.955))

    out = HERE / "figures" / "cmd_overlay_zp_before_after.png"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=150)
    print(f"같은 별 {n} 개")
    print(f"썼다: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
