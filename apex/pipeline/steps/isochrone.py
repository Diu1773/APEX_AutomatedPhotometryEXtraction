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

from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths_cmd import step10_zp_dir, step12_iso_dir

CMD_TABLE = "median_by_ID_filter_wide_cmd.csv"


def _colors(params) -> list[tuple[str, str]]:
    """`isochrone.colors` as pairs, e.g. "B-V,V-R" -> [("B","V"), ("V","R")]."""
    raw = str(getattr(params.P, "iso_colors", "") or "").strip()
    pairs = []
    for token in raw.replace(";", ",").split(","):
        token = token.strip()
        if not token or "-" not in token:
            continue
        left, _, right = token.partition("-")
        if left.strip() and right.strip():
            pairs.append((left.strip(), right.strip()))
    return pairs


def _prior(params, name: str) -> tuple[float, float] | None:
    """A prior written as "value,sigma"; absent means no prior."""
    raw = str(getattr(params.P, name, "") or "").strip()
    if not raw:
        return None
    parts = [p.strip() for p in raw.replace(";", ",").split(",")]
    if len(parts) != 2:
        return None
    try:
        return float(parts[0]), float(parts[1])
    except ValueError:
        return None


def build_fit_config(params, iso_file: Path | None = None):
    """Turn the config rows into an `IsochroneFitConfig`.

    Kept out of `run` so a test can check the translation without an MCMC.
    """
    from apex.analysis.cmd.isochrone_fit_service import IsochroneFitConfig

    P = params.P
    path = iso_file or (str(getattr(P, "iso_file_path", "") or "").strip() or None)
    return IsochroneFitConfig(
        colors=_colors(params),
        mag_band=str(getattr(P, "iso_mag_band", "g") or "g"),
        iso_file=str(path) if path else None,
        age_bounds=(float(P.iso_age_min), float(P.iso_age_max)),
        mh_bounds=(float(P.iso_mh_min), float(P.iso_mh_max)),
        dm_bounds=(float(P.iso_dm_min), float(P.iso_dm_max)),
        ecolor_bounds=(float(P.iso_ecolor_min), float(P.iso_ecolor_max)),
        mh_prior=_prior(params, "iso_mh_prior"),
        ecolor_prior=_prior(params, "iso_ecolor_prior"),
        dm_prior=_prior(params, "iso_dm_prior"),
        parallax_distance_prior=bool(P.iso_parallax_distance_prior),
        parallax_dm_sigma=float(P.iso_parallax_dm_sigma),
        parallax_dm_window=float(P.iso_parallax_dm_window),
        use_membership=bool(P.iso_use_membership),
        data_snr_min=float(P.iso_data_snr_min),
        fit_snr_min=float(P.iso_fit_snr_min),
        max_stars=int(P.iso_max_stars),
        n_walkers=int(P.iso_n_walkers),
        n_burn=int(P.iso_n_burn),
        n_steps=int(P.iso_n_steps),
        f_bin=float(P.iso_f_bin),
        f_field=float(P.iso_f_field),
        err_floor=float(P.iso_err_floor),
        seed=int(P.iso_seed),
    )


def missing_decisive_settings(params) -> list[str]:
    """What must be written down before a batch fit means anything.

    Only two, and both for the reason in the module docstring: the colours
    because the fit is undefined without them, and the isochrone grid because
    there is no sensible default file.
    """
    missing = []
    if not _colors(params):
        missing.append("isochrone.colors (예: \"B-V,V-R\")")
    if not str(getattr(params.P, "iso_file_path", "") or "").strip():
        missing.append("isochrone.file_path (PARSEC 격자 파일)")
    return missing


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

        (out_dir / "isochrone_fit_summary.json").write_text(
            json.dumps(summary, indent=1, default=str), encoding="utf-8")
        converged = summary.get("converged")
        note = "converged" if converged else "did NOT converge — see the summary"
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=f"isochrone fit {note}", duration_s=elapsed,
        )
