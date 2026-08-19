"""LC Step 10 (headless): detrending, and the figures that come with it.

This was the one LC stage still marked deferred, and the reason given was that
its calculation read from widgets and wrote to them — a refactor rather than a
move. Half of that was right. The writes are presentation, and they are hooks
now. The reads were already funnelled through a single method: the window's
`_sync_state_from_controls` copied every spin box into a plain attribute, so
the calculation had been reading attributes all along and only four widget
reads bypassed it.

What was genuinely missing is what this step supplies: somewhere for those
attributes to come from when nobody is turning a dial.

The figures come with it. They were never drawn headless because the only
figure in reach belonged to a `FigureCanvas` — not because the drawing needed
Qt, which it does not.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import List

from apex.analysis.light_curve.target_config import read_selection
from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths_lc import step9_lc_dir, step10_detrend_dir

# config attribute -> runner setting
SETTINGS = {
    "lc_detrend_mode": "mode",
    "lc_detrend_sigma_clip": "sigma_clip",
    "lc_detrend_clip_sigma": "clip_sigma",
    "lc_detrend_clip_iters": "clip_iters",
    "lc_detrend_plot_view": "_plot_view_mode",
    "lc_detrend_color_by": "color_by",
    "lc_detrend_global_min_comps": "global_min_comps",
    "lc_detrend_global_sigma": "global_sigma",
    "lc_detrend_global_iters": "global_iters",
    "lc_detrend_global_rms_pct": "global_rms_pct",
    "lc_detrend_global_rms_threshold": "global_rms_threshold",
    "lc_detrend_global_frame_sigma": "global_frame_sigma",
    "lc_detrend_global_robust": "global_robust",
    "lc_detrend_global_normalize": "global_normalize",
    "lc_detrend_global_k2": "use_global_k2",
    "lc_detrend_sysrem_iter": "sysrem_iter",
    "lc_detrend_sysrem_apply": "sysrem_apply",
}

MODES = ("offset", "color", "global", "sysrem")


def _settings(params) -> dict:
    P = getattr(params, "P", params)
    out = {}
    for attr, name in SETTINGS.items():
        value = getattr(P, attr, None)
        if value is not None:
            out[name] = value
    return out


class LcDetrendStep(PipelineStep):
    index = 10
    key = "lcdetrend"
    name = "Detrend & night merge"

    def inputs(self, ctx: RunContext) -> List[Path]:
        return [step9_lc_dir(ctx.result_dir)]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step10_detrend_dir(ctx.result_dir)]

    def is_complete(self, ctx: RunContext) -> bool:
        out = step10_detrend_dir(ctx.result_dir)
        return out.exists() and any(out.glob("lightcurve_ID*_current.csv"))

    def run(self, ctx: RunContext) -> StepResult:
        selection = read_selection(ctx.result_dir)
        if not selection:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=("no target selection — run LC Step 8, or set "
                         "lightcurve.target_id in the config"),
            )
        target_id = int(selection.get("target_id") or 0)
        if target_id <= 0:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=f"the selection names no target: {selection}",
            )

        settings = _settings(ctx.params)
        mode = str(settings.get("mode", "offset")).strip().lower()
        if mode not in MODES:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=(f"lightcurve.detrend_mode={mode!r} is not one of "
                         + ", ".join(MODES)),
            )
        settings["mode"] = mode

        from apex.analysis.light_curve.detrend_runner import HeadlessDetrendRunner

        started = time.perf_counter()
        log = getattr(ctx, "log", None)
        runner = HeadlessDetrendRunner(
            ctx.params, [ctx.result_dir], logger=log, target_id=target_id,
            project_state=getattr(ctx, "project_state", None), **settings,
        )

        if not runner.load_raw_data(silent=True):
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=(f"no raw light curve for target ID {target_id} — run "
                         f"LC Step 9 first ({step9_lc_dir(ctx.result_dir)})"),
                duration_s=time.perf_counter() - started,
            )

        try:
            # `sync_controls=False`: there are no controls to sync from, and the
            # window's version of that call would reach for spin boxes.
            runner.fit_and_apply(update_ui=False, save_outputs=True,
                                 sync_controls=False)
        except Exception as exc:                        # noqa: BLE001
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.FAILED,
                message=f"the {mode} detrend failed: {exc}",
                duration_s=time.perf_counter() - started,
            )
        elapsed = time.perf_counter() - started

        out = step10_detrend_dir(ctx.result_dir)
        written = sorted(out.glob("*")) if out.exists() else []
        corrected = getattr(runner, "corrected_df", None)
        n_corr = 0 if corrected is None else int(len(corrected))
        if n_corr == 0:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=(f"the {mode} detrend corrected no points — check that "
                         "the nights and filter in the raw curve are the ones "
                         "the settings name"),
                duration_s=elapsed,
            )

        params_df = getattr(runner, "params_df", None)
        n_fits = 0 if params_df is None else int(len(params_df))
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=(f"target ID {target_id}, {mode} — {n_corr} points, "
                     f"{n_fits} fit group(s)"),
            duration_s=elapsed,
            outputs=[str(p) for p in written],
        )
