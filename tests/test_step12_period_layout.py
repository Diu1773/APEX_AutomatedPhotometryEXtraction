import numpy as np

import apex.analysis.light_curve.period_plot as plot_module
from apex.gui.workflow.lc.step11_period_analysis import PeriodAnalysisWindow


def test_step12_summary_uses_vertical_stack_on_narrow_canvas():
    assert PeriodAnalysisWindow._summary_uses_compact_layout(394)
    assert PeriodAnalysisWindow._summary_uses_compact_layout(899)
    assert PeriodAnalysisWindow._summary_uses_stacked_layout(394)
    assert not PeriodAnalysisWindow._summary_uses_stacked_layout(700)


def test_step12_summary_keeps_dashboard_layout_on_wide_canvas():
    assert not PeriodAnalysisWindow._summary_uses_compact_layout(900)
    assert not PeriodAnalysisWindow._summary_uses_compact_layout(1200)


def test_step12_periodogram_stacks_panels_on_narrow_canvas():
    assert PeriodAnalysisWindow._periodogram_uses_compact_layout(360)
    assert PeriodAnalysisWindow._periodogram_uses_compact_layout(619)
    assert not PeriodAnalysisWindow._periodogram_uses_compact_layout(620)


def test_step12_hides_control_column_below_desktop_breakpoint():
    assert PeriodAnalysisWindow._uses_compact_shell(900)
    assert PeriodAnalysisWindow._uses_compact_shell(1349)
    assert not PeriodAnalysisWindow._uses_compact_shell(1350)


def test_step12_reuses_check_star_periodogram(monkeypatch):
    class Spin:
        def __init__(self, value):
            self._value = value

        def value(self):
            return self._value

    window = PeriodAnalysisWindow.__new__(PeriodAnalysisWindow)
    window._check_star_plot_cache_key = ("lc.csv", 1, "V", "corr")
    window._check_star_ls_cache = {}
    window.min_period_spin = Spin(0.01)
    window.max_period_spin = Spin(2.0)
    window.samples_spin = Spin(10)
    calls = []

    def fake_compute(*args):
        calls.append(args)
        return {"best_period": 0.1}

    # `compute_ls` is called from `period_plot` since the figure moved out of
    # the window (2026-08-19). Patching the GUI module's name silently stopped
    # intercepting anything, and the test then asserted on an empty list.
    monkeypatch.setattr(plot_module, "compute_ls", fake_compute)
    check = {
        "time": np.array([1.0, 2.0]),
        "mag": np.array([0.1, -0.1]),
        "mag_err": np.array([0.01, 0.01]),
    }

    first = window._check_star_ls_result(check)
    second = window._check_star_ls_result(check)

    assert first is second
    assert len(calls) == 1
