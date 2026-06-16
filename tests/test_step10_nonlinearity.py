"""Unit tests for the Step 10 bridge-star detector-linearity diagnostic.

Builds a synthetic per-measurement ``obs`` table at two exposure levels and
checks that ``_write_nonlinearity_diag`` reports a flat delta for a linear
detector, a sloped delta when a brightness-dependent offset is injected, and
skips single-exposure data. Constructs a bare worker via ``__new__`` so no Qt
event loop is needed; skipped on the no-GUI CI where PyQt5 isn't installed.
"""

import pandas as pd
import pytest

pytest.importorskip("PyQt5")

from apex.gui.workflow.cmd.step10_zeropoint_calibration import ZeropointCalibrationWorker


def _worker():
    w = ZeropointCalibrationWorker.__new__(ZeropointCalibrationWorker)
    w._log = lambda *a, **k: None
    return w


def _obs(stars_mag, exptimes, delta_fn):
    rows = []
    for sid, mtrue in enumerate(stars_mag, start=1):
        for exp in exptimes:
            for _ in range(3):  # a few measurements per (star, level)
                rows.append({
                    "ID": sid, "FILTER": "V",
                    "mag_cal": mtrue + delta_fn(mtrue, exp),
                    "exptime": float(exp), "snr": 50.0,
                })
    return pd.DataFrame(rows)


STARS = [12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0]


def test_linear_detector_flat_delta(tmp_path):
    obs = _obs(STARS, [60, 240], lambda m, e: 0.0)
    _worker()._write_nonlinearity_diag(obs, tmp_path)

    summ = pd.read_csv(tmp_path / "nonlinearity_summary.csv")
    assert len(summ) == 1
    row = summ.iloc[0]
    assert abs(row["slope_mag_per_mag"]) < 0.02
    assert not bool(row["nonlinearity_flag"])
    assert int(row["n_bridge"]) == len(STARS)
    assert (tmp_path / "nonlinearity_by_exposure.csv").exists()
    assert (tmp_path / "nonlinearity_check.png").exists()


def test_nonlinear_detector_sloped_delta(tmp_path):
    # Long (240s) exposure carries a brightness-dependent offset: a star of mag m
    # reads 0.1*(m-15) fainter. The diagnostic plots delta vs mag_l (which itself
    # includes the offset), so the recovered slope is -0.1/1.1 = -1/11.
    obs = _obs(STARS, [60, 240], lambda m, e: (0.1 * (m - 15.0) if e == 240 else 0.0))
    _worker()._write_nonlinearity_diag(obs, tmp_path)

    summ = pd.read_csv(tmp_path / "nonlinearity_summary.csv")
    row = summ.iloc[0]
    assert row["slope_mag_per_mag"] == pytest.approx(-1.0 / 11.0, abs=1e-6)
    assert bool(row["nonlinearity_flag"])


def test_low_snr_faint_bias_does_not_trigger_nonlinearity(tmp_path):
    obs = _obs(
        STARS,
        [60, 240],
        lambda m, e: (0.25 * (m - 16.5) if e == 60 and m > 16.5 else 0.0),
    )
    faint_short = (obs["exptime"] == 60.0) & (obs["mag_cal"] > 16.5)
    obs.loc[faint_short, "snr"] = 4.0

    _worker()._write_nonlinearity_diag(obs, tmp_path)

    row = pd.read_csv(tmp_path / "nonlinearity_summary.csv").iloc[0]
    assert row["snr_min"] == pytest.approx(20.0)
    assert abs(row["slope_mag_per_mag"]) < 0.02
    assert not bool(row["nonlinearity_flag"])


def test_single_exposure_skips(tmp_path):
    obs = _obs(STARS, [60], lambda m, e: 0.0)
    _worker()._write_nonlinearity_diag(obs, tmp_path)
    # Returns before writing any output for a single exposure level.
    assert not (tmp_path / "nonlinearity_summary.csv").exists()
    assert not (tmp_path / "nonlinearity_by_exposure.csv").exists()
