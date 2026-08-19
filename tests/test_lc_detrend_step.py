"""LC step 10: detrending, and why it was not the refactor it looked like.

The step was deferred with the reason that its calculation reads its inputs
from widgets and writes its results to them. The writes were real and are hooks
now. The reads were not: the window's `_sync_state_from_controls` copied every
spin box into a plain attribute in one place, so the calculation had been
reading attributes all along and exactly four widget reads bypassed it.

Measured against the corrected curve the window saved (YZ Boo, two nights, 364
points, `E:/APEX_validation/reprocess/YZBoo_2n`): all 26 columns and all six
fit-parameter rows agree to 0.0e+00. That needs the workspace, so it lives in
`docs/audit/LC_HEADLESS_STEPS_8_TO_11.md`. What is testable here is the seam.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import pytest

from apex.analysis.light_curve.detrend_runner import DetrendRunner, HeadlessDetrendRunner
from apex.pipeline.base import StepStatus
from apex.pipeline.context import RunContext
from apex.pipeline.steps.lc_detrend import SETTINGS, LcDetrendStep, _settings

# The calculation. Named rather than discovered, so deleting one from the base
# fails here instead of quietly shrinking the check.
SHARED = [
    "fit_and_apply",
    "load_raw_data",
    "_run_sysrem",
    "_run_global_ensemble",
    "_save_comprehensive_results",
    "_save_global_results",
    "_selected_dates",
    "_update_plots",
    "_plot_global_diagnostics",
    "_apply_plot_view",
]

# What the window supplies that a batch run supplies differently. An accidental
# addition means the batch path silently uses a default nobody chose; an
# accidental removal means the window stops reading its own controls.
WINDOW_ONLY = [
    "_target_id_text",
    "_filter_selection",
    "_use_global_k2",
    "_plot_figure",
    "_plot_redraw",
    "_ensure_plot_drawn",
    "_refresh_style",
    "_tell_user",
    "log",
    "_set_busy_state",
    "_set_busy_message",
    "_sync_mode_controls_from_state",
]


def _params(tmp_path, **over):
    body = dict(result_dir=str(tmp_path), data_dir=str(tmp_path))
    body.update(over)
    return SimpleNamespace(P=SimpleNamespace(**body))


def _ctx(tmp_path, params=None):
    return RunContext(mode="lc", params=params or _params(tmp_path),
                      result_dir=tmp_path, data_dir=tmp_path,
                      logger=logging.getLogger("test"))


def test_the_detrend_is_importable_without_qt():
    """A batch install has no PyQt5, and this module must not reach for it.

    It used to: seventeen `QMessageBox` calls came across with the move. They
    are `self._tell_user(...)` now — the same trade `PsfPhotometryRunner` made
    when it stopped calling `.emit()`.
    """
    import subprocess
    import sys

    code = ("import sys;"
            "import apex.analysis.light_curve.detrend_runner as m;"
            "print('PyQt5' in sys.modules)")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "False", "importing the detrend pulled in Qt"


def test_the_detrend_does_not_import_the_gui_layer():
    import ast
    import pathlib

    src = pathlib.Path("apex/analysis/light_curve/detrend_runner.py").read_text(
        encoding="utf-8")
    for node in ast.parse(src).body:
        mod = getattr(node, "module", "") or ""
        assert not mod.startswith("apex.gui"), (
            f"module-scope import of {mod} — analysis may not depend on gui")


def test_the_window_is_the_analysis_code():
    pytest.importorskip("PyQt5")
    from apex.gui.workflow.lc.step10_detrend_merge import DetrendNightMergeWindow

    drifted = [n for n in SHARED
               if getattr(DetrendNightMergeWindow, n, None) is not getattr(DetrendRunner, n)]
    assert not drifted, (
        f"the window no longer inherits {drifted} — a second copy of the "
        "detrend has appeared, and the two can now disagree")


def test_the_window_still_reads_its_own_controls():
    """The regression this guards: moving the calculation out leaves the window
    inheriting the batch defaults, so its spin boxes stop mattering — and
    nothing fails, the figure just quietly uses 0.01-10 d."""
    pytest.importorskip("PyQt5")
    from apex.gui.workflow.lc.step10_detrend_merge import DetrendNightMergeWindow

    for name in WINDOW_ONLY:
        assert getattr(DetrendNightMergeWindow, name) is not getattr(DetrendRunner, name), (
            f"{name} must be the window's own — otherwise the window uses the "
            "batch default instead of what the user set")


def test_the_runner_needs_nothing_the_class_does_not_declare():
    import ast
    import inspect

    cls = next(n for n in ast.parse(inspect.getsource(DetrendRunner)).body
               if isinstance(n, ast.ClassDef))
    read, written = set(), set()
    for node in ast.walk(cls):
        if (isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name)
                and node.value.id == "self"):
            (written if isinstance(node.ctx, ast.Store) else read).add(node.attr)
    provided = set(dir(DetrendRunner)) | written | {"__dict__", "__class__"}
    missing = sorted(read - provided)
    assert not missing, f"the detrend reads {missing} and the class never provides them"


def test_an_unknown_setting_is_refused_rather_than_ignored():
    """A typo that silently does nothing is worse than one that stops the run."""
    with pytest.raises(TypeError, match="unknown detrend setting"):
        HeadlessDetrendRunner(None, [], target_id=1, glabal_sigma=3.0)


def test_every_config_row_maps_to_a_real_setting():
    """`SETTINGS` is the bridge between the config and the runner. A name on
    either side with nothing on the other is a setting that does nothing."""
    for attr, name in SETTINGS.items():
        assert hasattr(DetrendRunner, name), (
            f"config {attr} maps to {name}, which DetrendRunner does not have")


def test_the_batch_figure_shows_every_panel_and_says_so():
    """The window shows one panel because its canvas shares a screen with the
    controls. A file has no such constraint, and the two panels the single view
    hides — the raw curve, and Delta-mag against airmass — are the ones that
    make the figure worth keeping."""
    assert DetrendRunner._plot_view_mode == "all"


def test_it_blocks_when_no_target_was_chosen(tmp_path):
    result = LcDetrendStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.BLOCKED
    assert "Step 8" in result.message or "target_id" in result.message


def test_it_blocks_on_a_mode_it_does_not_have(tmp_path):
    from apex.analysis.light_curve.target_config import LcTarget, write_selection

    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    params = _params(tmp_path, lc_detrend_mode="polynomial")
    result = LcDetrendStep().run(_ctx(tmp_path, params))
    assert result.status == StepStatus.BLOCKED
    assert "detrend_mode" in result.message


def test_it_blocks_when_step_9_has_produced_no_curve(tmp_path):
    from apex.analysis.light_curve.target_config import LcTarget, write_selection

    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    result = LcDetrendStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.BLOCKED
    assert "Step 9" in result.message


def test_settings_come_across_from_the_config(tmp_path):
    params = _params(tmp_path, lc_detrend_mode="global", lc_detrend_clip_sigma=2.5,
                     lc_detrend_global_min_comps=7, lc_detrend_plot_view="corr")
    got = _settings(params)
    assert got["mode"] == "global"
    assert got["clip_sigma"] == 2.5
    assert got["global_min_comps"] == 7
    assert got["_plot_view_mode"] == "corr"


def test_the_new_settings_are_actually_read():
    from apex.config.config_audit import unread_settings

    dead = set(unread_settings(mode="lc")["dead"])
    for name in SETTINGS:
        assert name not in dead, f"{name} 을 읽는 코드가 없다"


def test_the_lc_pipeline_has_no_deferred_step_left():
    """LC reached 7 until 2026-08-19. It reaches 11 now, and nothing on the way
    waits for a window."""
    from apex.pipeline.base import DeferredStep
    from apex.pipeline.registry import get_steps

    steps = get_steps("lc")
    assert [s.index for s in steps] == list(range(1, 12))
    deferred = [s.key for s in steps if isinstance(s, DeferredStep)]
    assert not deferred, f"still deferred: {deferred}"
