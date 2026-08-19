"""LC step 11: the period, headless.

Nothing had to move for this one. `PeriodAnalysisWorker` is a pass-through — it
takes arrays, calls `run_period_analysis()`, emits the dict — and loading,
file choice and saving were already services. What was missing was a caller
that does not need someone to type a range into a spin box.

Measured against the analysis the window saved (YZ Boo, two nights, 364 points,
`E:/APEX_validation/reprocess/YZBoo_2n`), with the same filter grouping and the
same input file: all four methods agree to 0.0e+00, and PDM lands 0.19 % from
the literature period (0.104092 d). That needs the workspace, so it lives in
`docs/audit/LC_HEADLESS_STEPS_9_AND_11.md` rather than here; what is testable here is
the wiring and the refusals.
"""

from __future__ import annotations

import json
import logging
from types import SimpleNamespace

import pandas as pd
import pytest

from apex.analysis.light_curve.target_config import LcTarget, write_selection
from apex.pipeline.base import StepStatus
from apex.pipeline.context import RunContext
from apex.pipeline.steps.lc_period import LcPeriodStep, _methods
from apex.utils.step_paths_lc import step9_lc_dir, step11_period_dir


def _params(tmp_path, **over):
    body = dict(
        result_dir=str(tmp_path), data_dir=str(tmp_path),
        lc_target_id=-1, lc_filter="",
        lc_period_min_days=0.01, lc_period_max_days=10.0,
        lc_period_samples_per_peak=10, lc_period_methods="ls,pdm",
        lc_period_pdm_bins=10,
    )
    body.update(over)
    return SimpleNamespace(P=SimpleNamespace(**body))


def _ctx(tmp_path, params=None):
    return RunContext(mode="lc", params=params or _params(tmp_path),
                      result_dir=tmp_path, data_dir=tmp_path,
                      logger=logging.getLogger("test"))


def _curve(tmp_path, target_id=5, period=0.1, n=200):
    """A sinusoid at a known period, so the answer is checkable."""
    import numpy as np

    out = step9_lc_dir(tmp_path)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20260819)
    t = np.sort(rng.uniform(0, 3.0, n))
    mag = 12.0 + 0.30 * np.sin(2 * np.pi * t / period) + rng.normal(0, 0.004, n)
    pd.DataFrame({
        "file": [f"f{i:04d}.fit" for i in range(n)],
        "filter": ["V"] * n,
        "BJD_TDB": t + 2460000.0,
        "mag": mag,
        "mag_err": np.full(n, 0.004),
        "diff_mag_raw": mag - 12.0,
        "diff_err": np.full(n, 0.006),
    }).to_csv(out / f"lightcurve_ID{target_id}_raw.csv", index=False)


def test_it_blocks_when_no_target_was_chosen(tmp_path):
    result = LcPeriodStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.BLOCKED
    assert "Step 8" in result.message or "target_id" in result.message


def test_it_blocks_when_step_9_has_produced_no_curve(tmp_path):
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    result = LcPeriodStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.BLOCKED
    assert "Step 9" in result.message


def test_it_blocks_on_an_empty_search_window(tmp_path):
    """A period outside the range is simply not found, so an inverted range
    would return a confident wrong answer rather than nothing."""
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    _curve(tmp_path)
    params = _params(tmp_path, lc_period_min_days=5.0, lc_period_max_days=0.5)
    result = LcPeriodStep().run(_ctx(tmp_path, params))
    assert result.status == StepStatus.BLOCKED
    assert "search window" in result.message


def test_it_recovers_a_known_period(tmp_path):
    """0.1 d sinusoid in, 0.1 d out — the wiring carries real numbers."""
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    _curve(tmp_path, period=0.1)
    result = LcPeriodStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.OK, result.message

    files = sorted(step11_period_dir(tmp_path).glob("period_analysis_*.json"))
    assert files, "no analysis written"
    body = json.loads(files[0].read_text(encoding="utf-8"))
    periods = [v.get("best_period") for v in (body.get("results") or {}).values()
               if isinstance(v, dict) and v.get("best_period")]
    assert periods, f"no method reported a period: {body.get('results', {}).keys()}"
    assert any(abs(p - 0.1) < 0.002 for p in periods), (
        f"injected 0.1 d, recovered {periods}")


def test_it_runs_each_filter_separately_by_default(tmp_path):
    """A period belongs to the star, but each band measures it on its own —
    disagreement between bands is the diagnostic, so they are not merged
    unless `lightcurve.filter` says to."""
    import numpy as np

    write_selection(tmp_path, LcTarget(target_id=5), [1, 2])
    out = step9_lc_dir(tmp_path)
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(7)
    rows = []
    for flt in ("g", "r"):
        t = np.sort(rng.uniform(0, 3.0, 150))
        mag = 12.0 + 0.3 * np.sin(2 * np.pi * t / 0.1) + rng.normal(0, 0.004, 150)
        rows.append(pd.DataFrame({
            "file": [f"{flt}{i:04d}.fit" for i in range(150)], "filter": flt,
            "BJD_TDB": t + 2460000.0, "mag": mag, "mag_err": 0.004,
            "diff_mag_raw": mag - 12.0, "diff_err": 0.006,
        }))
    pd.concat(rows, ignore_index=True).to_csv(out / "lightcurve_ID5_raw.csv", index=False)

    result = LcPeriodStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.OK, result.message
    names = {p.name for p in step11_period_dir(tmp_path).glob("period_analysis_*.json")}
    assert names == {"period_analysis_g_ID5.json", "period_analysis_r_ID5.json"}, names


def test_the_configured_filter_wins(tmp_path):
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2])
    _curve(tmp_path)
    result = LcPeriodStep().run(_ctx(tmp_path, _params(tmp_path, lc_filter="all")))
    assert result.status == StepStatus.OK, result.message
    names = {p.name for p in step11_period_dir(tmp_path).glob("period_analysis_*.json")}
    assert names == {"period_analysis_all_ID5.json"}, names


def test_methods_parse_and_never_come_back_empty(tmp_path):
    """Two different empties, deliberately two different answers.

    An unset value means "I did not choose", so it gets the default pair.
    A value that is present but yields no method means something was meant and
    could not be read — that falls back to the one method that always works,
    rather than silently running the pair the user did not ask for.
    """
    assert _methods(_params(tmp_path, lc_period_methods="LS, PDM")) == ["ls", "pdm"]
    assert _methods(_params(tmp_path, lc_period_methods="ls;bls")) == ["ls", "bls"]
    assert _methods(_params(tmp_path, lc_period_methods="")) == ["ls", "pdm"]
    assert _methods(_params(tmp_path, lc_period_methods="  ,  ")) == ["ls"]
    assert _methods(_params(tmp_path, lc_period_methods=None)) == ["ls", "pdm"]


def test_the_new_settings_are_actually_read():
    """The guard that rejected six unread rows when step 8 was built."""
    from apex.config.config_audit import unread_settings

    dead = set(unread_settings(mode="lc")["dead"])
    for name in ("lc_period_min_days", "lc_period_max_days",
                 "lc_period_samples_per_peak", "lc_period_methods",
                 "lc_period_pdm_bins"):
        assert name not in dead, f"{name} 을 읽는 코드가 없다"


def test_the_lc_pipeline_now_reaches_step_11():
    """10 is deferred: detrending is the one LC stage whose calculation reads
    from widgets and writes to them, so it is a refactor rather than a move."""
    from apex.pipeline.base import DeferredStep
    from apex.pipeline.registry import get_steps

    steps = get_steps("lc")
    assert [s.index for s in steps] == list(range(1, 12))
    assert [s.key for s in steps][-4:] == [
        "lctarget", "lclightcurve", "lcdetrend", "lcperiod"]
    assert isinstance({s.index: s for s in steps}[10], DeferredStep)


def test_the_deferred_detrend_still_recognises_gui_work(tmp_path):
    """A deferred step must see a window-produced result as complete."""
    from apex.pipeline.registry import get_steps
    from apex.utils.step_paths_lc import step10_detrend_dir

    step = {s.index: s for s in get_steps("lc")}[10]
    ctx = _ctx(tmp_path)
    assert not step.is_complete(ctx)
    step10_detrend_dir(tmp_path).mkdir(parents=True, exist_ok=True)
    (step10_detrend_dir(tmp_path) / "lightcurve_ID5_offset.csv").write_text("a\n1\n")
    assert step.is_complete(ctx)


# ── the summary figure ──────────────────────────────────────────────────────
#
# The window drew it and nobody else could. It is the picture that says whether
# a period is real — a clean fold at the adopted period, or a smear that says
# the peak was an alias of the sampling window.

PLOT_SHARED = [
    "_update_summary_plot",
    "_summary_data_type",
    "_summary_period",
    "_summary_method_result",
    "_median_center_for_plot",
    "_check_star_plot_data",
    "_check_star_ls_result",
    "_two_harmonic_phase_model",
    "_phase_plot_check_time_column",
    "_summary_uses_compact_layout",
    "_summary_uses_stacked_layout",
]

# What the window supplies from widgets, and the batch path is handed. An
# accidental addition here means the batch figure silently uses a default the
# user never chose; an accidental removal means the window stops reading its
# own controls — which is exactly what happened while this was being built.
PLOT_WINDOW_ONLY = [
    "_summary_figure",
    "_summary_canvas_width",
    "_summary_redraw",
    "_search_window",
    "_show_alias_marks",
    "_load_check_star",
]


def test_the_figure_code_does_not_reach_for_qt():
    import subprocess
    import sys

    code = ("import sys;"
            "import apex.analysis.light_curve.period_plot as m;"
            "print('PyQt5' in sys.modules)")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "False"


def test_the_figure_does_not_import_the_gui_layer_at_module_scope():
    """`analysis` must not depend on `gui`. The palette is looked up lazily so
    a themed window still recolours its plot, and a batch run needs nothing."""
    import ast
    import pathlib

    src = pathlib.Path("apex/analysis/light_curve/period_plot.py").read_text(encoding="utf-8")
    for node in ast.parse(src).body:
        mod = getattr(node, "module", "") or ""
        assert not mod.startswith("apex.gui"), (
            f"module-scope import of {mod} — analysis may not depend on gui")


def test_the_window_draws_the_analysis_figure():
    pytest.importorskip("PyQt5")
    from apex.analysis.light_curve.period_plot import PeriodSummaryPlotter
    from apex.gui.workflow.lc.step11_period_analysis import PeriodAnalysisWindow

    drifted = [n for n in PLOT_SHARED
               if getattr(PeriodAnalysisWindow, n, None) is not getattr(PeriodSummaryPlotter, n)]
    assert not drifted, (
        f"the window no longer inherits {drifted} — a second copy of the "
        "summary figure has appeared")


def test_the_window_still_reads_its_own_controls():
    """The regression this test exists for: moving the figure out left the
    window inheriting the batch defaults, so its spin boxes stopped mattering."""
    pytest.importorskip("PyQt5")
    from apex.analysis.light_curve.period_plot import PeriodSummaryPlotter
    from apex.gui.workflow.lc.step11_period_analysis import PeriodAnalysisWindow

    for name in PLOT_WINDOW_ONLY:
        assert getattr(PeriodAnalysisWindow, name) is not getattr(PeriodSummaryPlotter, name), (
            f"{name} must be the window's own — otherwise the window uses the "
            "batch default instead of what the user set")


def test_the_plotter_needs_nothing_the_class_does_not_declare():
    """Same guard as the light-curve builder, for the same reason: an attribute
    the window happened to set in its constructor fails deep inside a draw."""
    import ast
    import inspect

    from apex.analysis.light_curve.period_plot import PeriodSummaryPlotter

    cls = next(n for n in ast.parse(inspect.getsource(PeriodSummaryPlotter)).body
               if isinstance(n, ast.ClassDef))
    read, written = set(), set()
    for node in ast.walk(cls):
        if (isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name)
                and node.value.id == "self"):
            (written if isinstance(node.ctx, ast.Store) else read).add(node.attr)

    provided = set(dir(PeriodSummaryPlotter)) | written | {"__dict__"}
    missing = sorted(read - provided)
    assert not missing, f"the figure reads {missing} and the class never provides them"


def test_it_draws_a_figure_for_a_real_period(tmp_path):
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    _curve(tmp_path, period=0.1)
    result = LcPeriodStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.OK, result.message

    figs = list(step11_period_dir(tmp_path).glob("period_summary_*.png"))
    assert figs, "no summary figure written"
    assert figs[0].stat().st_size > 20_000, "figure is suspiciously small"


def test_alias_resolution_runs_by_default(tmp_path):
    """Off by default it would adopt whichever peak Lomb-Scargle liked, and on
    a multi-night run that is routinely an alias of the sampling window."""
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    _curve(tmp_path, period=0.1)
    LcPeriodStep().run(_ctx(tmp_path))

    body = json.loads(sorted(step11_period_dir(tmp_path)
                             .glob("period_analysis_*.json"))[0].read_text(encoding="utf-8"))
    assert body.get("alias_analysis"), "no alias analysis in the saved result"
    assert body["alias_analysis"].get("adopted_period")


def test_alias_resolution_can_be_turned_off(tmp_path):
    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    _curve(tmp_path, period=0.1)
    params = _params(tmp_path, lc_period_resolve_aliases=False)
    LcPeriodStep().run(_ctx(tmp_path, params))

    body = json.loads(sorted(step11_period_dir(tmp_path)
                             .glob("period_analysis_*.json"))[0].read_text(encoding="utf-8"))
    assert not body.get("alias_analysis")


def test_nightly_offset_detrending_forces_the_raw_series(tmp_path):
    """Nightly-offset detrending removes the between-night baseline a
    multi-night period search needs. Using the corrected series after it moved
    YZ Boo's adopted period from 0.11 % off the literature value to 1.22 %.
    """
    import numpy as np

    from apex.pipeline.steps.lc_period import _resolve_aliases

    n = 200
    rng = np.random.default_rng(11)
    t = np.sort(rng.uniform(0, 3.0, n))
    raw = 12.0 + 0.3 * np.sin(2 * np.pi * t / 0.1) + rng.normal(0, 0.004, n)
    lc_data = {
        "time": t, "mag_raw": raw, "mag_corr": raw + 5.0, "mag_err": np.full(n, 0.004),
        "night_id": np.zeros(n, dtype=int),
        "correction_preserves_nightly_baseline": False,
    }
    out = _resolve_aliases(lc_data, {}, 0.01, 10.0, 10)
    assert out is not None
    assert out["input_series"] == "raw", (
        "the corrected series was used after a correction that does not "
        "preserve the nightly baseline")
    assert out["nightly_baseline_preserved"] is False
