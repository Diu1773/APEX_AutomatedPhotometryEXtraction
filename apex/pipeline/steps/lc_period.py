"""LC Step 11 (headless): find the period.

This one needed almost nothing. The window's `PeriodAnalysisWorker` is a pure
pass-through — it takes arrays in, calls `run_period_analysis()`, and emits the
dict that comes back. Loading is `load_period_lightcurve_csv()`, choosing the
input file is `find_best_lightcurve_csv()`, saving is
`save_period_analysis_outputs()`. All of it already lived in the analysis layer.

So unlike step 9, nothing had to move. What was missing was a caller that does
not need someone to type a period range into a spin box, and the range is the
one setting here worth thinking about: a period outside `[period_min_days,
period_max_days]` is simply not found. The defaults (0.01–10 d) are the
window's own, and they cover δ Scuti pulsators through most eclipsing binaries.

It runs on step 9's raw curve, so it does not wait for detrending —
`find_best_lightcurve_csv` prefers a detrended curve when step 10 has produced
one and falls back to the raw one when it has not. Detrending improves the
answer; it is not a precondition for having one.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import List

from apex.analysis.light_curve.target_config import read_selection
from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths_lc import (
    find_best_lightcurve_csv,
    step9_lc_dir,
    step11_period_dir,
)


def _methods(params) -> List[str]:
    P = getattr(params, "P", params)
    raw = str(getattr(P, "lc_period_methods", "ls,pdm") or "ls,pdm")
    out = [m.strip().lower() for m in raw.replace(";", ",").split(",") if m.strip()]
    return out or ["ls"]


def _method_spread(results) -> "float | None":
    """How far apart the periodogram methods land, as a fraction of the median.

    Only meaningful when the alias analysis could not settle which peak is real.
    A large spread there is the run telling you two nights were not enough.
    """
    import numpy as np

    periods = []
    for key in ("corr_ls", "raw_ls", "corr_pdm", "raw_pdm", "corr_bls", "raw_bls"):
        found = (results or {}).get(key) or {}
        value = found.get("best_period")
        if value and float(value) > 0:
            periods.append(float(value))
    if len(periods) < 2:
        return None
    med = float(np.median(periods))
    if med <= 0:
        return None
    return float((max(periods) - min(periods)) / med)


def _write_candidates(result_dir, target_id: int, flt: str, results: dict,
                      alias_analysis: dict | None):
    """Every period this run found, in one table, ranked by strength.

    Two nights of a variable star do not determine one period — the sampling
    leaves aliases that fit as well as the truth, and which method lands on
    which is not a property of the software. Measured on two objects it went
    both ways: YZ Boo's Lomb-Scargle peak is 9.2 % from the literature value
    and its PDM peak 0.19 %; on AE UMa it is the other way round, 0.01 % against
    0.04 %. So the code has no basis to choose, and choosing anyway publishes
    one of those errors as a number.

    What it can do is lay the candidates out next to each other, with the
    aliases the resolver ranked, so the answer is picked by comparing against
    what the object is known to do. That is the table this writes.
    """
    import csv

    import numpy as np

    rows = []
    labels = {"ls": "Lomb-Scargle", "pdm": "PDM", "bls": "BLS"}
    for key, found in (results or {}).items():
        if not isinstance(found, dict):
            continue
        period = found.get("best_period")
        if not period or float(period) <= 0:
            continue
        series, _, method = key.partition("_")
        rows.append({
            "source": "periodogram",
            "method": labels.get(method, method.upper()),
            "series": "corrected" if series == "corr" else "raw",
            "period_days": float(period),
            "period_hours": float(period) * 24.0,
            "strength": found.get("best_power"),
            "rank": "",
            "note": "",
        })

    analysis = alias_analysis or {}
    status = str(analysis.get("status", "")).upper()
    for candidate in (analysis.get("candidates") or [])[:8]:
        period = candidate.get("period")
        if not period or float(period) <= 0:
            continue
        rows.append({
            "source": "alias candidate",
            "method": "alias resolver",
            "series": analysis.get("input_series", ""),
            "period_days": float(period),
            "period_hours": float(period) * 24.0,
            "strength": candidate.get("bic"),
            "rank": candidate.get("rank", ""),
            "note": ("adopted" if candidate.get("rank") == 1 and status == "RESOLVED"
                     else ""),
        })

    if not rows:
        return None
    rows.sort(key=lambda r: r["period_days"])

    out_dir = Path(step11_period_dir(result_dir))
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"period_candidates_{flt or 'all'}_ID{int(target_id)}.csv"
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "source", "method", "series", "period_days", "period_hours",
            "strength", "rank", "note"])
        writer.writeheader()
        writer.writerows(rows)

    periods = np.array([r["period_days"] for r in rows], dtype=float)
    return {
        "path": path,
        "n": len(rows),
        "min": float(periods.min()),
        "max": float(periods.max()),
        "status": status,
    }


def _resolve_aliases(lc_data, results, min_period, max_period, samples, log=None):
    """Rank the sampling-window aliases of the strongest peak.

    Same call the window's worker makes when its alias checkbox is ticked, on
    the same series it plots: the corrected magnitudes when the correction
    preserved the nightly baseline, the raw ones otherwise.
    """
    import numpy as np

    from apex.analysis.light_curve.period_alias_service import analyze_period_aliases
    from apex.analysis.light_curve.period_analysis_service import compute_ls

    # Nightly-offset detrending removes exactly the between-night baseline a
    # multi-night period search needs, so the corrected series must not be used
    # for alias work after it. `load_period_lightcurve_csv` carries the flag;
    # ignoring it moved YZ Boo's adopted period from 0.104209 d to 0.105359 d —
    # 0.11 % from the literature value to 1.22 %.
    preserves = bool(lc_data.get("correction_preserves_nightly_baseline", True))
    mag_corr = lc_data.get("mag_corr")
    use_corrected = preserves and mag_corr is not None and np.any(
        np.isfinite(np.asarray(mag_corr, float)))
    mag = mag_corr if use_corrected else lc_data["mag_raw"]

    ls = (results.get("corr_ls") if use_corrected else results.get("raw_ls"))         or results.get("raw_ls")
    if not ls or "error" in ls or "frequency" not in ls:
        ls = compute_ls(lc_data["time"], mag, lc_data.get("mag_err"), "alias-scan",
                        min_period, max_period, samples)
    if not ls or "frequency" not in ls:
        return None

    if log:
        log("Evaluating sampling-window aliases...")
    analysis = analyze_period_aliases(
        lc_data["time"], mag, lc_data.get("mag_err"), lc_data.get("night_id"),
        ls["frequency"], ls["power"], min_period, max_period, harmonics=2,
    )
    analysis["input_series"] = "corrected" if use_corrected else "raw"
    analysis["nightly_baseline_preserved"] = preserves
    return analysis


class LcPeriodStep(PipelineStep):
    index = 11
    key = "lcperiod"
    name = "Period analysis"

    def inputs(self, ctx: RunContext) -> List[Path]:
        return [step9_lc_dir(ctx.result_dir)]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step11_period_dir(ctx.result_dir)]

    def is_complete(self, ctx: RunContext) -> bool:
        out = step11_period_dir(ctx.result_dir)
        return out.exists() and any(out.glob("period_analysis_*.json"))

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

        lc_file = find_best_lightcurve_csv(ctx.result_dir, target_id)
        if lc_file is None or not Path(lc_file).exists():
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=(f"no light curve for target ID {target_id} — run LC "
                         f"Step 9 first ({step9_lc_dir(ctx.result_dir)})"),
            )

        import pandas as pd

        from apex.analysis.light_curve.period_analysis_service import run_period_analysis
        from apex.analysis.light_curve.period_io_service import (
            load_period_lightcurve_csv, save_period_analysis_outputs,
        )

        P = getattr(ctx.params, "P", ctx.params)
        min_period = float(getattr(P, "lc_period_min_days", 0.01) or 0.01)
        max_period = float(getattr(P, "lc_period_max_days", 10.0) or 10.0)
        samples = int(getattr(P, "lc_period_samples_per_peak", 10) or 10)
        pdm_bins = int(getattr(P, "lc_period_pdm_bins", 10) or 10)
        methods = _methods(ctx.params)
        resolve_aliases = bool(getattr(P, "lc_period_resolve_aliases", True))

        if not (0 < min_period < max_period):
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=(f"the search window is empty: period_min_days={min_period} "
                         f"period_max_days={max_period}"),
            )

        # One run per filter present in the curve — a period is a property of
        # the star, but each band measures it separately and disagreement
        # between bands is itself the diagnostic.
        try:
            head = pd.read_csv(lc_file)
        except Exception as exc:                        # noqa: BLE001
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.FAILED,
                message=f"could not read {Path(lc_file).name}: {exc}",
            )
        configured = str(getattr(P, "lc_filter", "") or "").strip()
        if configured:
            filters = [configured]
        elif "filter" in head.columns:
            filters = sorted({str(v) for v in head["filter"].dropna().unique() if str(v).strip()})
        else:
            filters = [""]
        if not filters:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=f"{Path(lc_file).name} names no filter to analyse",
            )

        started = time.perf_counter()
        log = getattr(ctx, "log", None)
        written: List[str] = []
        best: List[str] = []
        failed: List[str] = []
        analysed = 0                       # filters, not files — `written`
                                           # also holds figures and periodograms

        for flt in filters:
            try:
                lc_data = load_period_lightcurve_csv(Path(lc_file), flt, target_id)
                results = run_period_analysis(
                    time=lc_data["time"],
                    mag_raw=lc_data["mag_raw"],
                    mag_corr=lc_data.get("mag_corr"),
                    mag_err=lc_data.get("mag_err"),
                    min_period=min_period,
                    max_period=max_period,
                    samples_per_peak=samples,
                    methods=methods,
                    pdm_n_bins=pdm_bins,
                    progress_cb=log,
                )
                # Which peak to believe. Without this the adopted period is
                # whichever one Lomb-Scargle liked, and on a two-night run that
                # is routinely an alias of the sampling window — YZ Boo adopts
                # 0.0945 d that way and 0.1042 d with resolution, against a
                # literature 0.104092 d. The window offers it as a checkbox;
                # a batch run has nobody to tick it, so it defaults to on.
                alias_analysis = None
                if resolve_aliases:
                    alias_analysis = _resolve_aliases(
                        lc_data, results, min_period, max_period, samples, log)

                path = save_period_analysis_outputs(
                    result_dir=Path(ctx.result_dir), lc_data=lc_data, results=results,
                    min_period=min_period, max_period=max_period,
                    alias_analysis=alias_analysis,
                )
                written.append(str(path))

                # Not a pick — a list. See `_write_candidates`.
                spread_info = _write_candidates(
                    ctx.result_dir, target_id, lc_data.get("filter", flt),
                    results, alias_analysis)
                if spread_info is not None:
                    written.append(str(spread_info["path"]))

                # The same three panels the window draws — light curve,
                # periodogram, phase fold. A period without the fold is a
                # number; the fold is what says whether it is the right one.
                try:
                    from apex.analysis.light_curve.period_plot import (
                        save_period_summary_figure,
                    )
                    fig_path = step11_period_dir(ctx.result_dir) / (
                        f"period_summary_{lc_data.get('filter', flt) or 'all'}"
                        f"_ID{target_id}.png")
                    drawn = save_period_summary_figure(
                        fig_path, results=results, lc_data=lc_data,
                        alias_analysis=alias_analysis,
                        # No check star handed in: the plotter resolves it the
                        # way the window does, from the filter it actually
                        # plotted. Passing one keyed on the *requested* filter
                        # loses it whenever that is "all".
                        params=ctx.params,
                        search_window=(min_period, max_period, samples),
                    )
                    if drawn is not None:
                        written.append(str(drawn))
                except Exception as exc:              # noqa: BLE001
                    if log:
                        log(f"[figure] {flt or 'all'}: {exc}")
                # Report what the run concluded, by the same rule the figure
                # folds at: the alias candidate only when the service says it
                # resolved, the periodogram peak otherwise.
                status = str((alias_analysis or {}).get("status", "")).upper()
                adopted = (alias_analysis or {}).get("adopted_period")
                if status == "RESOLVED" and adopted:
                    best.append(f"{flt or 'all'}={float(adopted):.6f} d")
                else:
                    # Unresolved, so the fallback picks whichever method comes
                    # first — and on two nights those methods can disagree by
                    # more than the answer is worth. YZ Boo with its nights
                    # correctly labelled: LS 0.094468 d, PDM 0.104295 d, and the
                    # literature says 0.104092. Picking one silently would
                    # publish a 9 % error as a number; saying they disagree is
                    # the honest output, and it costs one line of message.
                    note = f" [alias {status.lower()}]" if status else ""
                    picked = None
                    for key in ("corr_ls", "raw_ls", "corr_pdm", "raw_pdm"):
                        found = (results or {}).get(key) or {}
                        period = found.get("best_period")
                        if period and picked is None:
                            picked = (key.split("_")[-1].upper(), float(period))
                    spread = _method_spread(results)
                    if picked:
                        line = f"{flt or 'all'} {picked[0]}={picked[1]:.6f} d{note}"
                        if spread is not None and spread > 0.02:
                            line += (f" — methods disagree by {100 * spread:.0f} %; "
                                     f"candidates in period_candidates_"
                                     f"{flt or 'all'}_ID{target_id}.csv")
                        best.append(line)
                analysed += 1
            except Exception as exc:                    # noqa: BLE001
                failed.append(f"{flt or 'all'}: {exc}")

        elapsed = time.perf_counter() - started
        if not written:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.FAILED,
                message=("period analysis produced nothing; "
                         + "; ".join(failed[:3]) if failed else "no filters analysed"),
                duration_s=elapsed,
            )

        note = f"target ID {target_id}, {analysed} filter(s)"
        if best:
            note += " — " + ", ".join(best)
        if failed:
            note += f"; {len(failed)} failed"
        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=note, duration_s=elapsed, outputs=written,
        )
