"""CMD Step 12 (headless): isochrone fit.

The fitting service in ``apex.analysis.cmd.isochrone_fit_service`` has never
needed Qt — it takes a table and an ``IsochroneFitConfig`` and returns a
posterior. What kept this step deferred was not the code but the settings: the
answer is decided by which colours are fitted, how wide the age window is, and
whether there is a reddening or distance prior, and none of that was expressible
in a config file. Running on library defaults does not fail loudly, it fails
quietly — the default 0.2–6 Gyr window cannot reach a globular at all, and
without an E(B-V) prior an open cluster rails at the floor (both measured; see
``validation/psf_crossinstrument/REPORT_UB_DEGENERACY.md``).

So the settings are now config rows, and this step refuses to run until the ones
that decide the answer are written down. A batch that produces nothing is a
smaller problem than a batch that produces a confident wrong age (2026-08-17).
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import List

import pandas as pd

from apex.analysis.cmd.isochrone_config import (
    build_fit_config,
    check_bounds,
    colors_from_params as _colors,
    missing_decisive_settings,
    prior_from_params as _prior,
)
from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths_cmd import step10_zp_dir, step12_iso_dir

CMD_TABLE = "median_by_ID_filter_wide_cmd.csv"

__all__ = [
    "IsochroneStep", "build_fit_config", "check_bounds",
    "missing_decisive_settings", "_colors", "_prior",
]


class IsochroneStep(PipelineStep):
    index = 12
    key = "isochrone"
    name = "Isochrone fit"

    def inputs(self, ctx: RunContext) -> List[Path]:
        return [step10_zp_dir(ctx.result_dir) / CMD_TABLE]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step12_iso_dir(ctx.result_dir) / "isochrone_fit_summary.json"]

    def is_complete(self, ctx: RunContext) -> bool:
        return self.outputs(ctx)[0].exists()

    def run(self, ctx: RunContext) -> StepResult:
        table = step10_zp_dir(ctx.result_dir) / CMD_TABLE
        if not table.exists():
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=f"no Step 10 table at {table}",
            )

        missing = missing_decisive_settings(ctx.params)
        if missing:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=("isochrone fit needs settings that decide the answer, "
                         "and running on defaults gives a confident wrong age: "
                         + "; ".join(missing)),
            )

        inverted = check_bounds(build_fit_config(ctx.params))
        if inverted:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=("search box is inverted, so the grid would be empty: "
                         + "; ".join(inverted)),
            )

        from apex.analysis.cmd.isochrone_fit_service import fit_cluster_isochrone

        out_dir = step12_iso_dir(ctx.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        df = pd.read_csv(table)

        started = time.perf_counter()
        progress = None
        if ctx.logger is not None:
            # The service reports (fraction, message), not a bare line.
            def progress(fraction: float, message: str) -> None:
                ctx.logger.info("%3.0f%% %s", 100.0 * float(fraction), message)
        result = fit_cluster_isochrone(
            df, build_fit_config(ctx.params), make_figures=True,
            progress_cb=progress,
        )
        elapsed = time.perf_counter() - started

        summary = getattr(result, "summary", None) or {}
        import json

        # Everything needed to say what this fit was, not just what it found:
        # a posterior without its bounds, its priors and its star count cannot
        # be judged (a median sitting on a bound is the wall, not a result).
        config = build_fit_config(ctx.params)
        record = {
            "summary": summary,
            "n_stars": getattr(result, "n_stars", None),
            "member_meta": getattr(result, "member_meta", None),
            "warnings": list(getattr(result, "warnings", []) or []),
            "settings": {k: getattr(config, k) for k in vars(config)},
            "elapsed_s": elapsed,
            "wide_table": str(table),
        }
        # Make the directory again, here. It was created before the fit, but the
        # fit can run for the better part of an hour and an empty directory under
        # a temp root does not necessarily survive that — one 40-minute M13 chain
        # (32 x 6000) finished and then died with FileNotFoundError on this very
        # write, throwing away the whole posterior. Re-creating costs nothing.
        out_dir.mkdir(parents=True, exist_ok=True)

        # The service built these and the step used to drop them on the floor,
        # so a batch run produced the age and not the picture of it — and the
        # picture is what goes in the paper. Same filenames the desktop window
        # writes, so a headless run and a window run leave the same products.
        figures = []
        for attr, name, dpi in (("cmd_figure", "mcmc_cmd_isochrone.png", 160),
                                ("corner_figure", "mcmc_corner.png", 130)):
            fig = getattr(result, attr, None)
            if fig is None:
                continue
            try:
                fig.savefig(out_dir / name, dpi=dpi, bbox_inches="tight")
                figures.append(name)
            except Exception:  # noqa: BLE001 - a figure must not fail the fit
                if ctx.logger is not None:
                    ctx.logger.exception("Could not save %s", name)
            finally:
                try:
                    import matplotlib.pyplot as plt
                    plt.close(fig)
                except Exception:  # noqa: BLE001
                    pass
        record["figures"] = figures

        (out_dir / "isochrone_fit_summary.json").write_text(
            json.dumps(record, indent=1, default=str), encoding="utf-8")

        # The service calls it `convergence_ok`. Reading a key the service does
        # not write returns None, which reports every fit as a failure.
        converged = bool(summary.get("convergence_ok"))
        note = "converged" if converged else "did NOT converge — see the summary"
        made = f"; {len(figures)} figures" if figures else "; no figures"
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=f"isochrone fit {note}{made}", duration_s=elapsed,
        )
