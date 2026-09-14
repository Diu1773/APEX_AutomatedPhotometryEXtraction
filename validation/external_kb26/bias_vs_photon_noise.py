"""이 치우침이 광자 잡음의 몇 배인가 (2026-09-15).

## 왜 이걸 재나

11~13 절에서 kb26 의 밝기 치우침이 무엇이고 언제 나타나는지는 갈랐는데, **그래서
얼마나 중요한가**는 아직 안 냈다. 결함의 크기를 절대값(0.28 등급)으로만 말하면
쓰는 사람이 판단을 못 한다. **잣대는 광자 잡음**이다 — 그보다 훨씬 작으면 무시해도
되고, 크면 그 등급을 믿으면 안 된다.

## 두 가지로 잰다

    별 하나       계통 / 그 별의 광자 오차          한 별의 등급을 쓸 때
    별 여럿       계통 / (광자 오차 / √N)          CMD·평균·이소크론에 쓸 때

**두 번째가 진짜 문제다.** 광자 잡음은 별을 모으면 √N 으로 줄어드는데 계통은 안
줄어들기 때문이고, 그래서 별이 많아질수록 계통의 비중이 커진다. 성단 CMD 나
이소크론 맞춤처럼 수백 개를 한꺼번에 쓰는 자리에서 이 비가 결정적이다.

## 대조군

사장님의 Moravian C3-61000 을 같이 잰다. **같은 프로그램·같은 성단·같은 지표로
재는데 카메라만 다르므로**, kb26 의 비가 크게 나오면 그것이 이 결함 탓이라는 것이
바로 드러난다.

실행:
    python -X utf8 validation/external_kb26/bias_vs_photon_noise.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).absolute().parents[2]
sys.path.insert(0, str(REPO))

from apex.analysis.cmd.zeropoint_runner import ZeropointCalibrationRunner  # noqa: E402
from apex.utils.gaia_quality import gaia_quality_report  # noqa: E402
from apex.utils.constants import MAG_ERR_COEFF  # noqa: E402
from apex.utils.io_utils import read_csv_int64_source_id  # noqa: E402

OUT = REPO / "validation/external_kb26/bias_vs_photon_noise.json"

WORKSPACES = {
    "kb26 (LCO 0.4m · SBIG 6303)": REPO / "validation/external_kb26/results",
    "kb27 (같은 모델 · 보름달 밤)": REPO / "validation/external_kb27/results",
    "Moravian (사장님 · CMOS)": Path("E:/APEX_validation/reprocess/M67/result"),
}
MAG_EDGES = [11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0]
MIN_CELL = 8
#: 적합과 같은 문턱. 신호가 거의 없는 별을 넣으면 기준 등급 자체가 못 미더워진다.
SNR_CUT = 20.0


class _Silent:
    def _log(self, *a, **k) -> None:
        return None


def _robust_fit(x, y):
    return ZeropointCalibrationRunner._robust_linfit(
        _Silent(), x, y, w=np.ones(len(x)), clip_sigma=3.0, iters=5,
        slope_absmax=0.8, min_n=10)


def _table(result_dir: Path, band: str):
    """보정성 별의 잔차·기준 등급·광자 오차를 한 표로.

    광자 오차는 `mag_cal_err_phot_*` 를 쓴다 — 이 열이 **광자 통계만** 담고 있어서
    (`MAG_ERR_COEFF / SNR`), 프레임 사이 흔들림이나 계통이 섞이지 않은 순수한
    잣대가 된다.
    """
    zp_dir = result_dir / "cmd_zeropoint"
    cal_path = zp_dir / "gaia_sdss_calibrator_by_ID.csv"
    wide_path = zp_dir / "median_by_ID_filter_wide.csv"
    coeff_path = zp_dir / "zp_fit_coefficients.csv"
    if not (cal_path.exists() and wide_path.exists() and coeff_path.exists()):
        return None
    cal = read_csv_int64_source_id(cal_path)
    wide = read_csv_int64_source_id(wide_path)
    coeff = pd.read_csv(coeff_path).set_index("filter")
    if band not in coeff.index:
        return None
    colour_col = f"color_{coeff.loc[band, 'color_col']}"
    if colour_col not in cal.columns or f"ref_{band}" not in cal.columns:
        return None

    # **광자 오차는 SNR 에서 직접 만든다.** `mag_cal_err_phot_*` 열은 나중에
    # 추가된 것이라 오래 전에 돌린 워크스페이스에는 없고, 그것 때문에 대조군이
    # 통째로 빠지면 비교가 성립하지 않는다. SNR 은 모든 워크스페이스에 있고
    # APEX 가 그 열을 만들 때 쓰는 식이 `MAG_ERR_COEFF / SNR` 이므로, 같은 식을
    # 여기서 쓰면 두 워크스페이스를 같은 자로 잰다.
    tab = cal
    if f"snr_{band}" not in tab.columns:
        return None

    try:
        good = np.asarray(gaia_quality_report(tab, cstar_nsigma=None)[0], bool)
    except (KeyError, ValueError, TypeError) as exc:
        print(f"    [주의] Gaia 품질 표를 못 만들어 전부 쓴다 — {exc}")
        good = np.ones(len(tab), bool)

    ref = pd.to_numeric(tab[f"ref_{band}"], errors="coerce").to_numpy(float)
    colour = pd.to_numeric(tab[colour_col], errors="coerce").to_numpy(float)
    delta = pd.to_numeric(tab.get(f"delta_{band}"), errors="coerce").to_numpy(float)
    snr = pd.to_numeric(tab[f"snr_{band}"], errors="coerce").to_numpy(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        phot = np.where(snr > 0, MAG_ERR_COEFF / snr, np.nan)

    g = (good & np.isfinite(ref) & np.isfinite(colour) & np.isfinite(delta)
         & np.isfinite(phot) & (phot > 0))
    g &= (snr >= SNR_CUT)
    if int(g.sum()) < 60:
        return None
    zp, ct, _, _ = _robust_fit(colour[g], delta[g])
    if not np.isfinite(zp):
        return None
    resid = delta[g] - (zp + ct * colour[g])
    resid -= np.median(resid)
    return pd.DataFrame(dict(mag=ref[g], resid=resid, phot=phot[g]))


def main() -> int:
    rows: list[dict] = []
    print("=== 밝기 치우침은 광자 잡음의 몇 배인가 ===")
    print(f"보정성 별을 SNR {SNR_CUT:.0f} 이상으로 자르고, 영점을 맞춘 뒤 남은 잔차의")
    print("구간 중앙값을 그 구간의 광자 오차와 견준다.")
    print()
    for label, rd in WORKSPACES.items():
        for band in ("r", "i"):
            t = _table(Path(rd), band)
            if t is None:
                continue
            print(f"-- {label} · {band} 밴드 · 별 {len(t)} 개")
            print(f"{'등급 구간':<12}{'별수':>6}{'계통(등급)':>12}"
                  f"{'광자오차':>10}{'별 하나':>9}{'구간 평균':>11}")
            print("-" * 60)
            for lo, hi in zip(MAG_EDGES[:-1], MAG_EDGES[1:]):
                m = (t["mag"] >= lo) & (t["mag"] < hi)
                n = int(m.sum())
                if n < MIN_CELL:
                    continue
                bias = float(np.median(t["resid"][m]))
                sigma = float(np.median(t["phot"][m]))
                one = abs(bias) / sigma if sigma > 0 else float("nan")
                many = one * np.sqrt(n)
                print(f"{lo:.0f}~{hi:.0f}".ljust(12) + f"{n:>6}{bias:>12.4f}"
                      f"{sigma:>10.4f}{one:>9.1f}배{many:>10.0f}배")
                rows.append(dict(workspace=label, band=band, lo=lo, hi=hi, n=n,
                                 bias=bias, photon=sigma, ratio_one=one,
                                 ratio_mean=many))
            print()
    print("  「별 하나」는 그 별 한 개의 등급을 쓸 때, 「구간 평균」은 그 구간의 별을")
    print("  다 모아 쓸 때다. **광자 잡음은 √N 으로 줄고 계통은 안 줄기 때문에**")
    print("  성단 CMD 처럼 수백 개를 한꺼번에 쓰는 자리에서 비가 훨씬 커진다.")

    OUT.write_text(json.dumps(rows, ensure_ascii=False, indent=2, default=str),
                   encoding="utf-8")
    print()
    print(f"썼다: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
