"""r 밴드의 어두운 쪽 휨이 색항 때문인지 — 색항만 바꿔 놓고 확인한다.

**닫는 시험이다.** 두 기기의 r 등급이 SNR 50~100 에서 0.27 어긋나는 것을 좇아
`validation/color_term_stability.py` 가 여기까지 왔다 — MuSCAT3 의 rp 는 적합이
쓴 별(153 개)에서 ct = −0.204 가 나오는데 보정성 표 전체(241 개)에서는 +0.021 이
나온다. 스물넷 조합 중 이런 일이 벌어지는 것은 이 하나뿐이다(차이 0.304, 2 등이
0.081, 중앙이 0.0045).

**그러면 −0.204 를 +0.021 로 바꿨을 때 어긋남이 사라지는지 보면 된다.** 사라지면
색항이 원인이고, 안 사라지면 다른 것을 찾아야 한다.

## 무엇을 바꾸나

APEX 의 영점 모형은 `r_std = r_inst + zp + ct·(g−r)` 이다. 적용된 ct 와 우리가
전체로 잰 ct 의 차이만큼을 등급에서 빼면, 다시 돌리지 않고도 「그 색항을 썼다면」
을 만들 수 있다.

    고친 등급 = 보정 등급 − (적용된 ct − 전체 ct) · (g−r)

**영점은 상수라 다시 잡는다.** 색항을 바꾸면 상수도 따라 바뀌는데 여기서 보는
것은 색과 밝기에 따른 어긋남이지 전체 이동이 아니다.

**이것은 사후 보정이지 APEX 산출물이 아니다.** 「APEX 가 이만큼 낸다」가 아니라
「색항이 원인이다」의 근거로만 읽어야 한다.

실행:
    python -X utf8 validation/external_muscat3/r_band_colorterm.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).absolute().parent
sys.path.insert(0, str(HERE.parents[1]))
sys.path.insert(0, str(HERE))

from apex.analysis.cmd.zeropoint_runner import robust_weighted_polyfit  # noqa: E402
from apex.utils.gaia_quality import gaia_quality_report  # noqa: E402
from compare_magnitudes import join, load_workspace, mad  # noqa: E402

M3 = HERE / "results"
MOR = Path("E:/APEX_validation/reprocess/M67/result")
SNR_CUT = 20.0


def fit_both_ways(result_dir: Path, band: str, colour: str):
    """(적합이 쓴 별로 잰 ct, 보정성 전체로 잰 ct, 각각의 별 수)."""
    d = pd.read_csv(result_dir / "cmd_zeropoint" / "gaia_sdss_calibrator_by_ID.csv")
    x = pd.to_numeric(d[f"color_{colour}"], errors="coerce").to_numpy(float)
    y = pd.to_numeric(d[f"delta_{band}"], errors="coerce").to_numpy(float)
    e = pd.to_numeric(d[f"mag_inst_err_{band}"], errors="coerce").to_numpy(float)
    s = pd.to_numeric(d[f"snr_{band}"], errors="coerce").to_numpy(float)
    q = np.asarray(gaia_quality_report(d, cstar_nsigma=None)[0], bool)
    w = 1.0 / np.maximum(e, 1e-6) ** 2
    base = np.isfinite(x) & np.isfinite(y) & np.isfinite(e) & (e > 0)
    m_fit = base & np.isfinite(s) & (s >= SNR_CUT) & q

    def one(m):
        c, n, _ = robust_weighted_polyfit(x[m], y[m], w=w[m], degree=1,
                                          clip_sigma=3.0, iters=5, min_n=10)
        return float(c[0]), int(n)

    return one(m_fit), one(base)


def bin_table(delta, axis, edges, label, name):
    """구간별 표. **비율은 원래 값이 작을 때 쓰지 않는다.**

    0.001 이 0.020 이 되는 것을 「1178 % 나빠졌다」고 적으면 20 밀리등급짜리
    변화가 표에서 제일 큰 사건처럼 보인다. 원래 값이 0.02 를 넘을 때만 비율을
    쓰고, 아닐 때는 얼마나 움직였는지를 밀리등급으로 적는다.
    """
    print(f"\n--- {name} ---")
    print(f"{label:>12}{'N':>6}{'지금':>10}{'고친 뒤':>11}{'변화':>14}")
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = np.isfinite(axis) & (axis >= lo) & (axis < hi)
        if m.sum() < 5:
            continue
        a, b = float(np.median(delta[0][m])), float(np.median(delta[1][m]))
        if abs(a) > 0.02:
            note = f"{(1 - abs(b) / abs(a)) * 100:>+11.0f} %"
        else:
            note = f"{(abs(b) - abs(a)) * 1000:>+9.1f} mmag"
        tag = f"{lo:g}~{hi:g}" if hi < 1e8 else f"{lo:g} 이상"
        print(f"{tag:>12}{int(m.sum()):>6}{a:>+10.4f}{b:>+11.4f}{note:>16}")
        rows.append(dict(bin=tag, n=int(m.sum()), before=a, after=b))
    return rows


def main() -> int:
    (ct_fit, n_fit), (ct_all, n_all) = fit_both_ways(M3, "r", "g_r")
    applied = float(pd.read_csv(M3 / "cmd_zeropoint" / "zp_fit_coefficients.csv")
                    .set_index("filter").loc["r", "ct"])
    print(f"MuSCAT3 r 의 색항 — 산출물에 적힌 값 {applied:+.4f}")
    print(f"  적합이 쓴 별로 다시 재면 {ct_fit:+.4f} (N={n_fit})")
    print(f"  보정성 전체로 재면       {ct_all:+.4f} (N={n_all})")
    d_ct = applied - ct_all
    print(f"  → 등급에서 뺄 몫 = ({applied:+.4f}) − ({ct_all:+.4f}) = {d_ct:+.4f} × (g−r)")

    j = join(load_workspace(M3, "MuSCAT3"), load_workspace(MOR, "Moravian"))
    num = lambda c: pd.to_numeric(j[c], errors="coerce").to_numpy(float)  # noqa: E731
    col = num("mag_cal_g_A") - num("mag_cal_r_A")
    snr = num("snr_r_A")
    before = num("mag_cal_r_A") - num("mag_cal_r_B")
    after = (num("mag_cal_r_A") - d_ct * col) - num("mag_cal_r_B")

    m = np.isfinite(before) & np.isfinite(after) & np.isfinite(col) & np.isfinite(snr)
    before, after, col, snr = before[m], after[m], col[m], snr[m]
    # 색항을 바꾸면 상수도 따라 바뀐다. 전체 이동이 아니라 색·밝기에 따른
    # 어긋남을 보는 것이므로 중앙값을 원래 자리에 맞춰 둔다.
    after = after - np.median(after) + np.median(before)

    print(f"\n=== r 밴드 · 같은 별 {int(m.sum())} 개 ===")
    print(f"{'':<16}{'중앙':>10}{'MAD':>10}{'색 기울기':>12}")
    for lab, v in (("지금 그대로", before), ("색항을 전체값으로", after)):
        print(f"{lab:<16}{np.median(v):>+10.4f}{mad(v):>10.4f}"
              f"{float(np.polyfit(col, v, 1)[0]):>+12.4f}")

    snr_rows = bin_table((before, after), snr,
                         np.array([0, 20, 50, 100, 300, 1e9]),
                         "SNR", "SNR 구간별 (MuSCAT3 − Moravian)")
    col_rows = bin_table((before, after), col,
                         np.arange(0.2, 1.7, 0.3),
                         "g−r", "색 구간별 (MuSCAT3 − Moravian)")

    out = HERE / "r_band_colorterm.json"
    out.write_text(json.dumps(
        {"ct_applied": applied, "ct_fit_subset": ct_fit, "n_fit": n_fit,
         "ct_all_calibrators": ct_all, "n_all": n_all, "delta_ct": d_ct,
         "n_stars": int(m.sum()),
         "mad_before": mad(before), "mad_after": mad(after),
         "slope_before": float(np.polyfit(col, before, 1)[0]),
         "slope_after": float(np.polyfit(col, after, 1)[0]),
         "by_snr": snr_rows, "by_color": col_rows},
        ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"\n썼다: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
