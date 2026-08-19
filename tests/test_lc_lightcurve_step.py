"""LC step 9: the window inherits the build rather than owning a copy of it.

The light curve was 1,149 lines living on a `QMainWindow` subclass — not because
it needed a window, but because that is where it was written. Everything it took
from the window was data: the parameter model, two caches, the night
assignments. Never a widget.

Moving it to `apex.analysis.light_curve.raw_lightcurve` and having the window
inherit it means identity is structural. `LightCurveBuilderWindow` does not call
the same calculation as the batch path — it *is* it, and no edit can make the
two drift without deleting an inheritance.

Measured against the curve the window itself saved (YZ Boo, two nights, 364
frames, `E:/APEX_validation/reprocess/YZBoo_2n`): every science column matches to
0.0e+00. That comparison needs 20 GB of frames so it is not a unit test; it is
recorded in `docs/audit/LC_STEP9_HEADLESS.md`. What is testable here is the
structure that makes the comparison stay true.
"""

from __future__ import annotations

import logging
from types import SimpleNamespace

import pytest

from apex.analysis.light_curve.raw_lightcurve import (
    HeadlessLightCurveBuilder,
    RawLightCurveBuilder,
)
from apex.pipeline.base import StepStatus
from apex.pipeline.context import RunContext
from apex.pipeline.steps.lc_lightcurve import LcLightCurveStep

# Methods that carry the science. Named rather than discovered, so that deleting
# one from the base class fails here instead of silently shrinking the check.
SHARED = [
    "_build_light_curve_core",
    "_build_ensemble_series",
    "_build_check_star_series",
    "_build_star_mag_series",
    "_compute_comp_qc",
    "_compute_airmass",
    "_compute_bjd_array",
    "_get_header",
    "_get_photometry_df",
    "_load_active_photometry_index",
    "_photometry_source_for_dir",
    "_preload_photometry_cache",
    "_map_comp_source_id",
    "_get_frame_exclude_map",
    "_refresh_photometry_source_policy",
]


def _params(tmp_path, **over):
    body = dict(result_dir=str(tmp_path), data_dir=str(tmp_path))
    body.update(over)
    return SimpleNamespace(P=SimpleNamespace(**body))


def _ctx(tmp_path):
    return RunContext(mode="lc", params=_params(tmp_path), result_dir=tmp_path,
                      data_dir=tmp_path, logger=logging.getLogger("test"))


def test_the_calculation_is_importable_without_qt():
    """A batch install has no PyQt5. The module must not reach for it."""
    import subprocess
    import sys

    code = (
        "import sys;"
        "import apex.analysis.light_curve.raw_lightcurve as m;"
        "print('PyQt5' in sys.modules)"
    )
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "False", "importing the light-curve build pulled in Qt"


def test_the_window_is_the_analysis_code():
    """Not 'calls the same function' — the same object."""
    pytest.importorskip("PyQt5")
    from apex.gui.workflow.lc.step9_lightcurve_builder import LightCurveBuilderWindow

    drifted = [
        name for name in SHARED
        if getattr(LightCurveBuilderWindow, name, None) is not getattr(RawLightCurveBuilder, name)
    ]
    assert not drifted, (
        f"the window no longer inherits {drifted} — a second copy of the "
        "light-curve build has appeared, and the two can now disagree"
    )


def test_the_window_keeps_its_own_reporting():
    """The three deliberate exceptions, so an accidental fourth is visible."""
    pytest.importorskip("PyQt5")
    from apex.gui.workflow.lc.step9_lightcurve_builder import LightCurveBuilderWindow

    for name in ("log", "_preload_progress", "_preload_finished"):
        assert getattr(LightCurveBuilderWindow, name) is not getattr(RawLightCurveBuilder, name), (
            f"{name} should be the window's own — it drives Qt widgets, and the "
            "base class deliberately keeps a no-op so batch runs skip the reporting"
        )


def test_the_headless_runner_supplies_every_attribute_the_build_reads():
    """The build reads instance state the window sets across a dozen callbacks.

    Missing one shows up as `AttributeError` deep inside a frame loop — which is
    how `_frame_exclude_cache` was found, several minutes into a real run. This
    reads the class rather than exercising it, so a new attribute is caught at
    the moment it is introduced.
    """
    import ast
    import inspect

    tree = ast.parse(inspect.getsource(RawLightCurveBuilder))
    read, written = set(), set()
    for node in ast.walk(tree):
        if (isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name)
                and node.value.id == "self"):
            (written if isinstance(node.ctx, ast.Store) else read).add(node.attr)

    provided = set(dir(RawLightCurveBuilder))
    runner = HeadlessLightCurveBuilder(_params("."), [], logger=lambda _m: None)
    provided |= set(vars(runner)) | written | {"__dict__", "__class__"}

    missing = sorted(read - provided)
    assert not missing, (
        f"the build reads {missing}, and HeadlessLightCurveBuilder never sets them — "
        "a batch run will fail partway through the frame loop"
    )


def test_it_blocks_when_no_target_was_chosen(tmp_path):
    """No selection on disk means step 8 never ran. Guessing a star is worse."""
    result = LcLightCurveStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.BLOCKED
    assert "Step 8" in result.message or "target_id" in result.message


def test_it_blocks_when_the_selection_has_no_comparisons(tmp_path):
    """A differential light curve without an ensemble is not a light curve."""
    from apex.analysis.light_curve.target_config import LcTarget, write_selection

    write_selection(tmp_path, LcTarget(target_id=5, target_name="AE UMa"), [])
    result = LcLightCurveStep().run(_ctx(tmp_path))
    assert result.status == StepStatus.BLOCKED
    assert "comparison" in result.message


def test_it_is_not_complete_before_it_runs(tmp_path):
    from apex.analysis.light_curve.target_config import LcTarget, write_selection

    step = LcLightCurveStep()
    assert not step.is_complete(_ctx(tmp_path))

    write_selection(tmp_path, LcTarget(target_id=5), [1, 2, 3])
    ctx = _ctx(tmp_path)
    assert not step.is_complete(ctx)
    assert step.outputs(ctx)[0].name == "lightcurve_ID5_raw.csv"
