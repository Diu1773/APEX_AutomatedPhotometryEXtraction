"""LC Step 8 (headless): resolve the target and its comparison ensemble.

The LC branch stopped at Step 7 — the pipeline registry lists steps 1-7 for LC
and nothing after. Not because the science needs a window: fourteen Qt-free
services in `apex.analysis.light_curve` do the building, detrending and period
work already. What stopped it is that nothing in the config could say which star
the light curve is of.

This step is that gate, and it blocks rather than guesses, the same way the
isochrone step does. A batch that produces nothing beats one that produces a
clean light curve of the wrong object — because a light curve of the wrong
object does not look wrong.

A gate refusing also does not look wrong, which is how this one spent its first
eight days blocked on every workspace in existence: it asked for
`master_sources.csv`, and Step 6 writes `ref_catalog.tsv`. The test suite agreed
with the mistake because its fixture invented the same name. Both now ask
`master_catalog_path()`, and a test holds the demand against what `RefBuildStep`
itself calls done.

The second thing it got wrong was quieter. `auto` mode ranked comparisons by
*catalogue order* — whatever `ref_catalog.tsv` happened to list first — because
the stability screen that ranks them by how steady they actually were lived
inside the window. On the same 364 YZ Boo frames the batch ensemble's average
came out 0.68 mag from the window's, and nothing in the output said the two runs
had used different stars. The screen is now `comparison_screening`, the window
calls the same functions, and `comparison_screening_<filter>.tsv` records every
measured star and the stage that dropped it.

The third was the gate refusing work that had already been done. The window
saves its choice as `selection_<filter>.json`, one per filter; this step read
the config and nothing else. So `apex run --mode lc` on a workspace where the
user had just picked YZ Boo in Step 8 answered that it "needs to know which star
it is of". Reading their choice is not guessing — the step now does, reports
which file it came from, and records `selected_by: window` in the output. The
config still wins wherever it speaks; the window fills its silence.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import List, Optional

import pandas as pd

from apex.analysis.light_curve.comparison_screening import (
    ScreeningResult, colors_from_catalog, screen_measurements,
    write_screening_report,
)
from apex.analysis.light_curve.photometry_source_service import (
    load_filter_photometry_timeseries,
)
from apex.analysis.light_curve.target_config import (
    LcTarget, filters_by_frame_count, missing_target_settings, read_target,
    read_window_selection, resolve_comparisons, write_selection,
)
from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.common_helpers import normalize_filter_key
from apex.utils.step_paths import forced_phot_input_dir, master_catalog_path
from apex.utils.step_paths_lc import step8_selection_dir


# Which filter stands for the workspace: the busiest one. The window has to
# answer this the same way, and it cannot import from `apex.pipeline`, so the
# rule lives in the analysis layer and both sides call it.
_available_filters = filters_by_frame_count


def _id_maps(catalog: pd.DataFrame) -> tuple[dict[int, int], dict[int, int]]:
    """`ID` ↔ `source_id`, both directions.

    The screen works in Gaia `source_id`; the selection file and every step
    after it work in the catalogue's display `ID`. Crossing that boundary by
    hand is how a target becomes a different star, so it happens once, here.
    """
    if catalog is None or catalog.empty:
        return {}, {}
    if "ID" not in catalog.columns or "source_id" not in catalog.columns:
        return {}, {}
    frame = catalog[["ID", "source_id"]].copy()
    frame["ID"] = pd.to_numeric(frame["ID"], errors="coerce")
    frame["source_id"] = pd.to_numeric(frame["source_id"], errors="coerce")
    frame = frame.dropna()
    id_to_sid = {int(r.ID): int(r.source_id) for r in frame.itertuples()}
    return id_to_sid, {sid: did for did, sid in id_to_sid.items()}


def _screen(
    ctx: RunContext, target: LcTarget, catalog: pd.DataFrame, target_sid: int
) -> tuple[Optional[ScreeningResult], str]:
    """Run the stability screen on one filter; return it and what to report.

    `lightcurve.filter` may hold the sentinel `all`, which the window writes to
    mean "these roles apply to every filter". Asking the photometry loader for a
    filter literally named `all` returns nothing, and the step then reported
    "no usable photometry" on a workspace holding 364 frames — a real answer
    hidden behind a word that is not a filter name.
    """
    available = _available_filters(ctx.result_dir)
    wanted = normalize_filter_key(target.filter_key) if target.filter_key else ""
    if wanted and wanted.lower() in ("all", "any", "*"):
        wanted = ""
    if wanted and available and wanted not in available:
        return None, (f"lightcurve.filter={wanted} has no photometry "
                      f"(present: {', '.join(available)}) — could not screen")
    if wanted:
        chosen, note = wanted, f"filter {wanted} (from lightcurve.filter)"
    elif available:
        chosen, note = available[0], f"filter {available[0]} (most frames)"
    else:
        return None, "no photometry index — could not screen"

    measurements, source_info = load_filter_photometry_timeseries(
        Path(ctx.result_dir), chosen, ctx.project_state
    )
    if measurements is None or measurements.empty:
        return None, f"no usable photometry in {chosen} — could not screen"

    P = getattr(ctx.params, "P", ctx.params)
    star_ids = pd.to_numeric(measurements["star_id"], errors="coerce").dropna()
    colors = colors_from_catalog(catalog, star_ids.astype("int64").unique())
    result = screen_measurements(
        measurements,
        int(target_sid),
        filter_key=chosen,
        source_info=source_info,
        desired_count=int(getattr(P, "comparison_auto_ensemble_max", 12)),
        pool_cap=int(getattr(P, "comparison_auto_pool_max", 30)),
        colors_by_id=colors,
        target_color=colors.get(int(target_sid), float("nan")),
    )
    return result, note


class LcTargetStep(PipelineStep):
    index = 8
    key = "lctarget"
    name = "Light-curve target"

    def inputs(self, ctx: RunContext) -> List[Path]:
        return [master_catalog_path(ctx.result_dir)]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step8_selection_dir(ctx.result_dir) / "lc_target_selection.json"]

    def is_complete(self, ctx: RunContext) -> bool:
        return self.outputs(ctx)[0].exists()

    def run(self, ctx: RunContext) -> StepResult:
        started = time.perf_counter()
        target = read_target(ctx.params)
        from_window = None

        # The window keeps its choice in `selection_<filter>.json`; this step
        # reads the config. So a user who picked a star in Step 8 and then ran
        # `apex run --mode lc` was told the run "needs to know which star it is
        # of" about a workspace where they had just said which. Reading what
        # they chose is not a guess — and it is announced, not assumed.
        if not target.is_resolved:
            from_window = read_window_selection(
                ctx.result_dir, target.filter_key,
                prefer=_available_filters(ctx.result_dir),
            )
            if from_window is not None:
                target = LcTarget(
                    target_id=int(from_window["target_id"]),
                    target_name=target.target_name,
                    comparison_ids=list(from_window["comparison_ids"]),
                    comparison_mode="manual",
                    comparison_count=target.comparison_count,
                    filter_key=str(from_window["filter"]),
                )

        missing = missing_target_settings(ctx.params) if from_window is None else []
        if missing:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=("a light curve needs to know which star it is of, and "
                         "no default is defensible: " + "; ".join(missing)
                         + " — or pick one in the Step 8 window, whose "
                           "selection_<filter>.json this step reads"),
            )

        master = master_catalog_path(ctx.result_dir)
        if not master.exists():
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=f"no master catalog at {master}",
            )

        try:
            catalog = pd.read_csv(master, sep="\t")
        except Exception as exc:                    # noqa: BLE001
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.FAILED,
                message=f"could not read the master catalog: {exc}",
            )

        known = set(pd.to_numeric(catalog.get("ID"), errors="coerce").dropna().astype(int))
        if target.target_id not in known:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=(f"lightcurve.target_id={target.target_id} is not in the "
                         f"master catalog ({len(known)} stars)"),
            )

        id_to_sid, sid_to_id = _id_maps(catalog)
        notes: list[str] = []
        outputs: list[str] = []
        comparisons: list[int] = []
        check_id: Optional[int] = None

        if from_window is not None:
            comparisons = list(from_window["comparison_ids"])
            check_id = from_window.get("check_id")
            notes.append(f"from the Step 8 window, filter {from_window['filter']}")
            notes.append(f"{Path(from_window['source_path']).name}")
        elif target.comparison_mode == "manual":
            comparisons = resolve_comparisons(target, catalog)
            notes.append("manual")
        else:
            target_sid = id_to_sid.get(int(target.target_id))
            if target_sid is None:
                notes.append("no source_id for the target — screening skipped")
            else:
                try:
                    screened, note = _screen(ctx, target, catalog, target_sid)
                except Exception as exc:            # noqa: BLE001
                    screened, note = None, f"screening failed: {exc}"
                notes.append(note)
                if screened is not None:
                    comparisons = [
                        sid_to_id[int(sid)] for sid in screened.selected_ids
                        if int(sid) in sid_to_id
                    ]
                    if screened.check_id is not None:
                        check_id = sid_to_id.get(int(screened.check_id))
                    notes.append(screened.funnel.as_text())
                    outputs.append(str(write_screening_report(ctx.result_dir, screened)))
            if not comparisons:
                # Catalogue order is not a ranking, so it is a fallback and it
                # is named as one rather than passed off as a selection.
                comparisons = resolve_comparisons(target, catalog)
                notes.append("fell back to catalogue order")

        if not comparisons:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=("no comparison stars — set lightcurve.comparison_ids, "
                         "or check that the master catalog has more than the target"),
            )

        path = write_selection(
            ctx.result_dir, target, comparisons, check_id=check_id,
            selected_by=("window" if from_window is not None
                         else "manual" if target.comparison_mode == "manual"
                         else "stability" if check_id is not None
                         else "catalog_order"),
        )
        outputs.insert(0, str(path))
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=(f"target ID {target.target_id}"
                     + (f" ({target.target_name})" if target.target_name else "")
                     + f"; {len(comparisons)} comparisons"
                     + (f" + check ID {check_id}" if check_id is not None else "")
                     + " [" + "; ".join(notes) + "]"),
            duration_s=time.perf_counter() - started,
            outputs=outputs,
        )
