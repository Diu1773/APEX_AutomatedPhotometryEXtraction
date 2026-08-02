"""Phase-plot y-limits must follow the data, never the fit overlay.

A 4-harmonic fit at an unconstraining period (two nights exactly one day
apart folded at the P=1 d default) runs away by ~1e6 mag in the phase gaps.
Autoscale used to include that curve, squashing the real light curve onto a
flat line at "0" of an 800,000-magnitude axis (seen on the real YZ Boo
2-night workspace, 2026-08-02).
"""
from __future__ import annotations

import os
import types

import numpy as np
import pytest

pytest.importorskip("PyQt5")

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from PyQt5.QtWidgets import QApplication  # noqa: E402

from apex.gui.tools.variable_star import VariableStarToolWindow  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance() or QApplication([])
    yield app


def _fake_params(tmp_path):
    P = types.SimpleNamespace(
        result_dir=str(tmp_path),
        data_dir=str(tmp_path),
        target_name="TESTSTAR",
    )
    return types.SimpleNamespace(P=P, param_file=str(tmp_path / "parameters.toml"))


def _window(qapp, tmp_path, monkeypatch):
    # No workspace on disk and no worker threads: the test drives lc_data
    # directly, exactly what the loader would have produced.
    monkeypatch.setattr(
        VariableStarToolWindow, "_load_lc_from_workspace",
        lambda self, **kwargs: None,
    )
    return VariableStarToolWindow(_fake_params(tmp_path), None)


def _two_nights_one_day_apart():
    """The YZ Boo geometry that triggered the runaway: ~2 h per night,
    nights separated by almost exactly the default fold period (1 d)."""
    t0 = 2460795.0
    rng = np.random.default_rng(42)
    t = np.concatenate([
        t0 - 1.0 + np.linspace(0.0, 0.08, 60),
        t0 + np.linspace(0.0, 0.08, 64),
    ])
    mag = 14.0 + 0.3 * np.sin(2 * np.pi * t / 0.104) + rng.normal(0, 0.01, t.size)
    return t, mag


def test_runaway_fit_cannot_hijack_the_y_axis(qapp, tmp_path, monkeypatch):
    win = _window(qapp, tmp_path, monkeypatch)
    try:
        t, mag = _two_nights_one_day_apart()
        win.lc_data = {"time": t, "mag": mag, "mag_err": None,
                       "filters": None, "night_id": None, "target_id": 153}
        win.t0_edit.setValue(2460795.0)
        win.phase_p.setValue(1.0)                 # the unconstraining default
        win.phase_fit_chk.setChecked(True)
        win.phase_fit_harm.setValue(4)

        win._update_phase_plot()
        ax = win.ph_canvas.figure.axes[0]
        top, bottom = ax.get_ylim()               # inverted axis: top > bottom

        # The fit really does run away — the guard is what keeps it off-axis.
        fit_extremes = [
            float(np.nanmax(np.abs(line.get_ydata())))
            for line in ax.lines if line.get_ydata() is not None and len(line.get_ydata())
        ]
        assert fit_extremes and max(fit_extremes) > 1e3

        # ...while the axis stays pinned to the data (13.7..14.3 plus margin).
        assert bottom < mag.min() < mag.max() < top
        assert (top - bottom) < 5 * (mag.max() - mag.min() + 0.1)
        assert top > bottom                        # magnitude axis: faint down
    finally:
        win.close()


def test_axis_still_fits_the_data_without_the_fit_overlay(qapp, tmp_path, monkeypatch):
    win = _window(qapp, tmp_path, monkeypatch)
    try:
        t, mag = _two_nights_one_day_apart()
        win.lc_data = {"time": t, "mag": mag, "mag_err": None,
                       "filters": None, "night_id": None, "target_id": 153}
        win.t0_edit.setValue(2460795.0)
        win.phase_p.setValue(0.104)               # a sensible period
        win.phase_fit_chk.setChecked(False)

        win._update_phase_plot()
        ax = win.ph_canvas.figure.axes[0]
        top, bottom = ax.get_ylim()
        assert bottom < mag.min() < mag.max() < top
    finally:
        win.close()
