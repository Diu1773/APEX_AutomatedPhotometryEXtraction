"""보고하는 등급 오차에는 광자 잡음 말고 재현성 바닥도 들어가야 한다.

Step 10 이 내던 `mag_cal_err` 는 광자·읽기 잡음뿐이었다. 별을 서로 견줄 때 그것은
실제보다 훨씬 작다 — 같은 프레임 안에서 기준별들이 영점 둘레로 흩어지는 폭이 같은
SNR 이 예측하는 광자 잡음보다 2~10 배 크고, **프레임을 쌓아도 안 줄어든다.**

2026-09-09 에 M67 을 세 기기(Moravian · MuSCAT3 · LCO kb26)로 재서 확인했다.
같은 별을 두 기기로 잰 차이의 산포를 예측하면

    광자 잡음만              관측/예측 19.7   설명 안 되는 몫 0.091 등급
    이 바닥을 √N 으로 나눠    관측/예측  3.7
    이 바닥을 그대로          관측/예측  1.3   설명 안 되는 몫 0.028 등급

**안 줄어든다는 것이 핵심이다.** 무작위 잡음이면 √N 으로 줄어야 하는데 안 주므로
별마다 붙박이인 계통이고, 그래서 나누지 않고 그대로 더한다.
근거: `validation/ERROR_BUDGET.md`.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from apex.analysis.cmd.zeropoint_runner import (
    apply_repeatability_floor,
    repeatability_floor_by_filter,
)


def _frames(**by_filter) -> pd.DataFrame:
    rows = []
    for filt, scatters in by_filter.items():
        for i, sc in enumerate(scatters):
            rows.append({"file": f"{filt}{i}.fits", "filter": filt,
                         "zp_frame": 20.0, "zp_scatter": sc, "n_ref": 100})
    return pd.DataFrame(rows)


def test_floor_is_the_median_frame_scatter_per_band():
    floors = repeatability_floor_by_filter(
        _frames(g=[0.020, 0.024, 0.028], r=[0.100, 0.110, 0.120]))
    assert floors["g"] == pytest.approx(0.024)
    assert floors["r"] == pytest.approx(0.110)


def test_no_frames_means_no_floor_rather_than_a_guess():
    assert repeatability_floor_by_filter(pd.DataFrame()) == {}
    assert repeatability_floor_by_filter(None) == {}
    # 산포가 0 이나 음수뿐이면 잴 것이 없다 — 0 을 바닥이라고 부르지 않는다.
    assert repeatability_floor_by_filter(_frames(g=[0.0, -1.0])) == {}


def test_the_floor_lands_in_the_reported_error():
    """이것이 고치기 전에 실패하는 시험이다 — 그때는 광자 값이 그대로 나왔다."""
    wide = pd.DataFrame({"mag_cal_err_g": [0.005, 0.010, np.nan],
                         "mag_cal_err_r": [0.004, 0.008, 0.012]})
    out, applied = apply_repeatability_floor(wide.copy(),
                                             {"g": 0.024, "r": 0.110})

    assert applied == {"g": pytest.approx(0.024), "r": pytest.approx(0.110)}
    assert out["mag_cal_err_g"][0] == pytest.approx(np.hypot(0.005, 0.024))
    assert out["mag_cal_err_r"][2] == pytest.approx(np.hypot(0.012, 0.110))
    # 광자만의 값은 사라지지 않는다.
    assert out["mag_cal_err_phot_g"][0] == pytest.approx(0.005)
    assert out["mag_cal_err_phot_r"][2] == pytest.approx(0.012)
    # 잴 수 없던 별은 그대로 잴 수 없다 — 바닥이 NaN 을 채우지 않는다.
    assert not np.isfinite(out["mag_cal_err_g"][2])


def test_the_floor_is_not_divided_by_the_frame_count():
    """√N 으로 나누면 M67 세 기기 시험에서 관측/예측이 1.3 이 아니라 3.7 이 된다."""
    wide = pd.DataFrame({"mag_cal_err_r": [0.004]})
    out, _ = apply_repeatability_floor(wide.copy(), {"r": 0.110})
    got = float(out["mag_cal_err_r"][0])

    assert got == pytest.approx(np.hypot(0.004, 0.110), rel=1e-6)
    for n_frames in (6, 10, 60):
        divided = np.hypot(0.004, 0.110 / np.sqrt(n_frames))
        assert got > divided, f"{n_frames} 장으로 나눈 값이 되어 버렸다"


def test_a_band_without_a_measured_floor_is_left_alone():
    wide = pd.DataFrame({"mag_cal_err_g": [0.005], "mag_cal_err_zs": [0.007]})
    out, applied = apply_repeatability_floor(wide.copy(), {"g": 0.024})

    assert set(applied) == {"g"}
    assert out["mag_cal_err_zs"][0] == pytest.approx(0.007)
