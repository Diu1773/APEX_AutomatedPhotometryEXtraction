"""필터마다의 밝기 치우침이 산출물에 나오는가.

영점 모형은 색에 대해서만 기울기를 갖는다 — `mag_cal = mag_inst + zp + ct·C`.
**밝기에 대한 기울기는 모형에 없다.** 그래서 기기가 밝기를 따라 치우치면 그
치우침이 통째로 잔차에 남고, 잔차의 산포(`fit_scatter_rms`)에 섞여 들어가
**별마다의 흩어짐처럼 보인다.**

실제로 그런 기기가 있었다. LCO 0.4 m(kb26)의 M67 은 밝은 쪽과 어두운 쪽의 영점이
네 밴드 모두 0.28~0.38 등급 다르고, 같은 시야를 같은 기준에 견준 Moravian 과
MuSCAT3 는 0.03 안쪽이다. 관측소 자신의 파이프라인으로 보정하고 관측소 자신의
측광으로 재도 같으므로 기기 쪽 성질이다(`validation/ERROR_BUDGET.md` 7~8 절).

APEX 는 이 수치를 이미 내고 있었지만 **CMD 에 쓰는 한 밴드에 대해서만**이었다
(`gaia_cmd_drift_by_mag.csv`). 영점을 맞춘 모든 필터에 대해 내야 한다.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO))

from apex.analysis.cmd.zeropoint_runner import (  # noqa: E402
    build_zp_qc_summary,
    magnitude_drift_by_filter,
)

RNG = np.random.default_rng(20260914)


def _coeff(filt: str = "r", zp: float = -5.0, ct: float = 0.3) -> pd.DataFrame:
    return pd.DataFrame([{"filter": filt, "zp": zp, "ct": ct, "ct2": 0.0,
                          "color_col": "r_i", "N": 200, "scatter_rms": 0.02}])


def _calibrators(n: int = 240, drift: float = 0.0, noise: float = 0.01,
                 zp: float = -5.0, ct: float = 0.3) -> pd.DataFrame:
    """밝기를 따라 `drift` 만큼 미끄러지는 보정성 별 표.

    가장 밝은 별과 가장 어두운 별 사이가 정확히 `drift` 가 되도록 만든다.
    """
    ref = np.linspace(11.0, 17.0, n)
    colour = RNG.uniform(0.0, 1.2, n)
    slope = drift / (ref.max() - ref.min())
    delta = zp + ct * colour + slope * (ref - ref.min())
    delta = delta + RNG.normal(0.0, noise, n)
    return pd.DataFrame({"ref_r": ref, "color_r_i": colour, "delta_r": delta})


def test_a_flat_instrument_reports_about_zero():
    """밝기를 안 타는 기기는 0 근처가 나와야 한다."""
    got = magnitude_drift_by_filter(_calibrators(drift=0.0), _coeff())
    assert "r" in got
    assert abs(got["r"]["drift"]) < 0.01


def test_a_sliding_instrument_reports_the_size_and_the_sign():
    """미끄러지는 기기는 크기와 부호가 다 맞아야 한다.

    잰 값은 아래 5 분의 1 과 위 5 분의 1 의 **중앙값** 차이다. 아래 5 분위의
    중앙값은 전체의 10 % 지점이고 위 5 분위의 중앙값은 90 % 지점이므로, 등급에
    대해 직선으로 미끄러지면 넣은 값의 **80 %** 가 나온다.
    """
    for put in (+0.30, -0.30):
        got = magnitude_drift_by_filter(_calibrators(drift=put), _coeff())["r"]
        assert got["drift"] == pytest.approx(0.8 * put, abs=0.02)
        assert np.sign(got["drift"]) == np.sign(put)


def test_the_colour_term_is_removed_first():
    """색항이 큰 기기에서도 밝기 치우침만 잡아야 한다.

    색과 밝기가 얽히지 않게 색을 무작위로 뿌려 놓았으므로, 색항을 제대로
    빼지 않으면 색의 산포가 밝기 치우침으로 새어 든다.
    """
    got = magnitude_drift_by_filter(
        _calibrators(drift=0.0, ct=0.8), _coeff(ct=0.8))["r"]
    assert abs(got["drift"]) < 0.015


def test_too_few_stars_reports_nothing_rather_than_a_number():
    """별이 모자라면 값을 지어내지 않는다."""
    assert magnitude_drift_by_filter(_calibrators(n=30), _coeff()) == {}


def test_a_filter_without_a_reference_column_is_skipped():
    """기준 등급이 없는 필터는 건너뛴다 — kb26 의 zs 가 그렇다."""
    cal = _calibrators().drop(columns=["ref_r"])
    assert magnitude_drift_by_filter(cal, _coeff()) == {}


def test_the_summary_carries_the_drift_columns():
    """`zp_qc_summary.csv` 에 열이 실제로 실린다."""
    summary = build_zp_qc_summary(_coeff(), None, None, None,
                                  _calibrators(drift=0.30))
    assert "bright_to_faint_drift" in summary.columns
    row = summary[summary["filter"] == "r"].iloc[0]
    assert float(row["bright_to_faint_drift"]) == pytest.approx(0.24, abs=0.02)
    assert int(row["n_drift_calibrators"]) > 0
    assert float(row["drift_bright_mag"]) < float(row["drift_faint_mag"])


def test_the_summary_still_builds_without_a_calibrator_table():
    """보정성 표를 안 줘도 예전처럼 만들어진다 — 있던 열은 그대로다."""
    summary = build_zp_qc_summary(_coeff(), None, None, None)
    assert not summary.empty
    assert "global_zp" in summary.columns
    assert np.isnan(float(summary.iloc[0]["bright_to_faint_drift"]))
    assert int(summary.iloc[0]["n_drift_calibrators"]) == 0


# ---------------------------------------------------------------------------
# 표본을 적합과 같게 자르는 문턱
# ---------------------------------------------------------------------------
#
# **보정성 별 표에는 적합이 안 쓴 별까지 들어 있다.** 거기에는 신호가 거의 없는
# 별도 섞여 있고, 그 별들은 기준 등급 자체가 못 미덥다 — 어두운 쪽에서 Gaia 의
# BP 가 오염되기 때문이다. kb26 의 B 는 그 별들을 넣느냐 빼느냐로 **값의 부호가
# 뒤집혔다**: SNR 20 이상만 두면 9.9~14.5 등급에서 −0.32, 20.2 등급까지 다 넣으면
# +0.26. 그래서 적합과 같은 문턱을 걸고, 어느 등급 범위에서 잰 값인지 함께 낸다.


def _with_snr(cal: pd.DataFrame, faint_tail_snr: float = 3.0) -> pd.DataFrame:
    """어두운 쪽 4 분의 1 만 신호가 약한 표. 그 별들은 잔차도 반대로 휜다."""
    cal = cal.copy()
    n = len(cal)
    snr = np.full(n, 100.0)
    tail = cal["ref_r"] >= cal["ref_r"].quantile(0.75)
    snr[tail.to_numpy()] = faint_tail_snr
    cal["snr_r"] = snr
    cal.loc[tail, "delta_r"] = cal.loc[tail, "delta_r"] + 0.8
    return cal


def test_the_snr_gate_keeps_the_low_signal_tail_out():
    """신호 없는 어두운 꼬리가 값을 끌고 가지 못한다."""
    cal = _with_snr(_calibrators(drift=0.20))
    gated = magnitude_drift_by_filter(cal, _coeff(), snr_cut=20.0)["r"]
    assert gated["drift"] == pytest.approx(0.8 * 0.20 * 0.75, abs=0.05)


def test_without_the_gate_the_same_data_gives_a_different_answer():
    """문턱을 안 걸면 같은 자료가 다른 답을 준다 — 그래서 문턱을 함께 적는다."""
    cal = _with_snr(_calibrators(drift=0.20))
    open_ = magnitude_drift_by_filter(cal, _coeff(), snr_cut=0.0)["r"]
    gated = magnitude_drift_by_filter(cal, _coeff(), snr_cut=20.0)["r"]
    assert abs(open_["drift"] - gated["drift"]) > 0.3
    assert open_["faint_mag"] > gated["faint_mag"]


def test_the_gate_and_the_magnitude_range_are_written_down():
    """값만 내면 재현이 안 된다 — 문턱과 등급 범위를 같이 낸다."""
    summary = build_zp_qc_summary(_coeff(), None, None, None,
                                  _with_snr(_calibrators(drift=0.20)), 20.0)
    row = summary[summary["filter"] == "r"].iloc[0]
    assert float(row["drift_snr_cut"]) == 20.0
    assert float(row["drift_bright_mag"]) < float(row["drift_faint_mag"])


def test_a_table_without_an_snr_column_is_not_emptied():
    """SNR 열이 없는 기기에서 별을 전부 버리면 안 된다."""
    cal = _calibrators(drift=0.20)
    assert "snr_r" not in cal.columns
    assert "r" in magnitude_drift_by_filter(cal, _coeff(), snr_cut=20.0)


def test_a_legacy_coefficient_file_without_ct2_still_works():
    """`ct2` 가 없거나 NaN 인 옛 계수 파일에서 필터가 사라지면 안 된다.

    `float(NaN) or 0.0` 은 NaN 을 돌려준다 — NaN 이 참이기 때문이다. 그대로
    쓰면 모형이 통째로 NaN 이 되어 그 필터가 조용히 빠진다.
    """
    cal = _calibrators(drift=0.20)
    legacy = _coeff().drop(columns=["ct2"])
    assert "r" in magnitude_drift_by_filter(cal, legacy)

    nan_ct2 = _coeff()
    nan_ct2.loc[0, "ct2"] = np.nan
    got = magnitude_drift_by_filter(cal, nan_ct2)
    assert "r" in got
    assert got["r"]["drift"] == pytest.approx(0.8 * 0.20, abs=0.02)


def test_a_filter_with_no_zeropoint_is_left_out_rather_than_faked():
    """영점을 못 맞춘 필터는 값을 내지 않는다 — kb26 의 zs 가 그렇다."""
    no_zp = _coeff()
    no_zp.loc[0, "zp"] = np.nan
    assert magnitude_drift_by_filter(_calibrators(), no_zp) == {}
