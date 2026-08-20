"""Choosing comparison stars by how they actually behaved, without a window.

`comparison_stability_service` has been Qt-free from the start: leave-one-out
residuals, stability metrics, adaptive ensemble, check-star recommendation. What
kept comparison selection inside the window was the two steps *around* it — the
per-star coverage screen that decides who is even a candidate, and the wiring
that hands the candidates to the stability search. Both lived as methods on
`TargetComparisonSelectionWindow`, reading `self.filter_catalogs`,
`self.filter_rejected_sources`, `self.params`.

None of that is presentation. Reading `self` is not the same as needing Qt, and
the difference showed up as a wrong answer: a headless run fell back to
*catalogue order* for its ensemble, and its comparison average came out 0.68 mag
from the window's on the same 364 frames. The batch run was not approximating
the window — it was doing something else entirely and saying nothing.

So the screen moves here as free functions, and the window calls them. Same
code, one copy, and the funnel counts come back as data instead of a log line.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Optional, Sequence

import numpy as np
import pandas as pd

from apex.analysis.light_curve.comparison_stability_service import (
    ComparisonSelectionConfig,
    select_adaptive_ensemble,
    select_stable_comparisons,
)

MIN_COVERAGE = 0.90
"""A candidate must be measured on at least this fraction of the frames.

Below it the leave-one-out residuals are computed on a different subset of
nights than the target's, and the resulting stability score compares two things
that were not observed together.
"""

MIN_COMPARISONS = 3
"""Fewer than three and a single misbehaving star moves the ensemble average."""


@dataclass
class ScreeningFunnel:
    """How many candidates each stage left standing.

    The endpoint alone ("3 comparisons + check 187") reads the same whether it
    came out of four thousand stars or out of five.
    """

    measured: int = 0
    coverage: int = 0
    eligible: int = 0
    pool: int = 0
    adopted: int = 0

    def as_dict(self) -> dict:
        return {
            "measured": self.measured,
            "coverage": self.coverage,
            "eligible": self.eligible,
            "pool": self.pool,
            "adopted": self.adopted,
        }

    def as_text(self) -> str:
        return (f"{self.measured} measured → {self.coverage} coverage "
                f"→ {self.eligible} eligible → {self.pool} pool "
                f"→ {self.adopted} adopted")


@dataclass
class ScreeningResult:
    """Everything one filter's screening produced, in source_id space."""

    filter_key: str
    target_id: int
    target_mag: float
    measurements: pd.DataFrame
    source_info: dict
    report: pd.DataFrame
    candidate_ids: list[int]
    metrics: pd.DataFrame
    selected_ids: list[int]
    check_id: Optional[int]
    check_metrics: dict = field(default_factory=dict)
    ensemble_trials: pd.DataFrame = field(default_factory=pd.DataFrame)
    reason: str = ""
    funnel: ScreeningFunnel = field(default_factory=ScreeningFunnel)
    stability: dict = field(default_factory=dict)
    """`select_stable_comparisons`'s own return, whole.

    Naming the two or three keys a caller happens to need makes the rest vanish
    silently: `active_ids` and `removed_ids` feed the window's stability report,
    and picking out `metrics` alone emptied both without any error.
    """


def _as_id_set(values: Optional[Iterable]) -> set[int]:
    # `if not values` raises on a numpy array ("truth value ... is ambiguous"),
    # and the caller catches it — so the screen quietly fell back to catalogue
    # order on exactly the input a pipeline step naturally hands it.
    if values is None:
        return set()
    out = set()
    for value in values:
        try:
            out.add(int(value))
        except (TypeError, ValueError):
            continue
    return out


def build_candidate_pool(
    measurements: pd.DataFrame,
    target_id: int,
    target_count: int,
    *,
    manual_rejects: Optional[Iterable] = None,
    variable_ids: Optional[Iterable] = None,
    is_variable: Optional[Callable[[int], bool]] = None,
    external_rejects: Optional[Iterable] = None,
    prefer_ids: Optional[Sequence] = None,
    pool_cap: int = 30,
    min_coverage: float = MIN_COVERAGE,
) -> tuple[list[int], pd.DataFrame, float]:
    """Who is eligible to be a comparison, and in what order to try them.

    The report keeps a row per measured star — including the rejected ones and
    why — because a run that found three candidates out of four thousand and a
    run that found three out of five need to be told apart afterwards.

    `mean_mag` holds a *median*; the name is what the on-disk stability report
    has always called it and renaming it here would break reading old reports.
    """
    total_frames = int(measurements["frame"].nunique())
    summary = (
        measurements.groupby("star_id")
        .agg(
            n=("frame", "nunique"),
            mean_mag=("mag", "median"),
            median_error=("mag_err", "median"),
        )
        .reset_index()
    )
    summary["star_id"] = pd.to_numeric(summary["star_id"], errors="coerce").astype("Int64")
    summary = summary[summary["star_id"].notna()].copy()
    summary["star_id"] = summary["star_id"].astype("int64")
    summary["coverage"] = summary["n"] / max(total_frames, 1)

    target_rows = summary[summary["star_id"] == int(target_id)]
    target_mag = (
        float(target_rows.iloc[0]["mean_mag"]) if not target_rows.empty else float("nan")
    )
    summary["d_mag"] = (
        np.abs(summary["mean_mag"] - target_mag) if np.isfinite(target_mag) else 0.0
    )

    summary["eligible"] = summary["coverage"] >= float(min_coverage)
    summary["basic_reason"] = ""
    summary.loc[summary["coverage"] < float(min_coverage), "basic_reason"] = "low_coverage"
    summary.loc[
        summary["star_id"] == int(target_id), ["eligible", "basic_reason"]
    ] = [False, "target"]

    rejects = _as_id_set(manual_rejects)
    if rejects:
        summary.loc[
            summary["star_id"].isin(rejects), ["eligible", "basic_reason"]
        ] = [False, "manual_reject"]

    # The variability screen only ever ran on stars that were still standing,
    # so a star already dropped for coverage keeps that reason rather than
    # being relabelled — the report says what actually removed it first.
    still_standing = set(summary.loc[summary["eligible"], "star_id"].astype("int64"))
    gaia_variable = _as_id_set(variable_ids) & still_standing
    if is_variable is not None:
        gaia_variable |= {sid for sid in still_standing if is_variable(int(sid))}
    if gaia_variable:
        summary.loc[
            summary["star_id"].isin(gaia_variable), ["eligible", "basic_reason"]
        ] = [False, "gaia_variable"]

    simbad_variable = _as_id_set(external_rejects)
    if simbad_variable:
        summary.loc[
            summary["star_id"].isin(simbad_variable), ["eligible", "basic_reason"]
        ] = [False, "simbad_variable"]

    eligible = summary[summary["eligible"]].copy()
    eligible = eligible.sort_values(
        ["d_mag", "median_error", "coverage", "star_id"],
        ascending=[True, True, False, True],
        na_position="last",
    )
    cap = max(3, int(pool_cap), int(target_count))
    eligible_ids = set(eligible["star_id"].astype("int64"))
    preferred = [
        int(star_id) for star_id in sorted(_as_id_set(prefer_ids))
        if int(star_id) in eligible_ids
    ]
    ranked = [int(value) for value in eligible["star_id"].tolist()]
    pool = preferred + [value for value in ranked if value not in preferred]
    return pool[:cap], summary, target_mag


def screen_measurements(
    measurements: pd.DataFrame,
    target_id: int,
    *,
    filter_key: str = "",
    source_info: Optional[dict] = None,
    desired_count: int = 12,
    pool_cap: int = 30,
    colors_by_id: Optional[dict] = None,
    color_for: Optional[Callable[[int], float]] = None,
    target_color: float = float("nan"),
    manual_rejects: Optional[Iterable] = None,
    variable_ids: Optional[Iterable] = None,
    is_variable: Optional[Callable[[int], bool]] = None,
    external_rejects: Optional[Iterable] = None,
    prefer_ids: Optional[Sequence] = None,
    should_stop: Optional[Callable[[], bool]] = None,
) -> ScreeningResult:
    """Screen, score, and pick — the whole selection for one filter.

    Raises rather than returning a half-answer: an ensemble too small to average
    is not a smaller ensemble, it is no result, and returning one lets a batch
    run write a light curve nobody can defend.
    """
    if measurements is None or measurements.empty:
        raise ValueError(f"No usable photometry time series for filter {filter_key}.")
    available = set(pd.to_numeric(measurements["star_id"], errors="coerce")
                    .dropna().astype("int64"))
    if int(target_id) not in available:
        raise ValueError(
            f"Target source {target_id} has no usable measurements in {filter_key}."
        )

    pool, report, target_mag = build_candidate_pool(
        measurements,
        int(target_id),
        int(desired_count),
        manual_rejects=manual_rejects,
        variable_ids=variable_ids,
        is_variable=is_variable,
        external_rejects=external_rejects,
        prefer_ids=prefer_ids,
        pool_cap=int(pool_cap),
    )
    if len(pool) < MIN_COMPARISONS + 1:
        raise ValueError(
            f"Only {len(pool)} candidates pass the "
            f"{int(MIN_COVERAGE * 100)}% coverage and exclusion checks in "
            f"{filter_key}; at least {MIN_COMPARISONS} comparisons plus one "
            "check star are required."
        )
    if should_stop is not None and should_stop():
        raise RuntimeError("Comparison analysis canceled.")

    # Colour is looked up *after* the pool is known, so a four-thousand-star
    # field costs thirty lookups rather than four thousand — and the pool is
    # computed once, not once here and once for the colours.
    colors = dict(colors_by_id or {})
    if color_for is not None:
        for star_id in pool:
            colors.setdefault(int(star_id), float(color_for(int(star_id))))
    config = ComparisonSelectionConfig(
        target_count=max(4, len(pool)), min_comparisons=MIN_COMPARISONS
    )
    scored = select_stable_comparisons(
        measurements,
        pool,
        target_mag=target_mag,
        target_color=target_color,
        color_by_id=colors,
        config=config,
    )
    if should_stop is not None and should_stop():
        raise RuntimeError("Comparison analysis canceled.")

    adaptive = select_adaptive_ensemble(
        measurements,
        scored.get("metrics", pd.DataFrame()),
        min_comparisons=MIN_COMPARISONS,
        max_comparisons=max(MIN_COMPARISONS, int(desired_count)),
    )
    selected = [int(value) for value in adaptive.get("selected_ids", [])]
    check_id = adaptive.get("check_id")

    funnel = ScreeningFunnel(
        measured=int(len(report)),
        coverage=int((pd.to_numeric(report["coverage"], errors="coerce")
                      >= MIN_COVERAGE).sum()),
        eligible=int(report["eligible"].astype(bool).sum()),
        pool=len(pool),
        adopted=len(selected),
    )
    return ScreeningResult(
        filter_key=str(filter_key),
        target_id=int(target_id),
        target_mag=target_mag,
        measurements=measurements,
        source_info=dict(source_info or {}),
        report=report,
        candidate_ids=list(pool),
        metrics=scored.get("metrics", pd.DataFrame()),
        selected_ids=selected,
        check_id=int(check_id) if check_id is not None else None,
        check_metrics=dict(adaptive.get("check_metrics", {})),
        ensemble_trials=adaptive.get("ensemble_trials", pd.DataFrame()),
        reason=str(adaptive.get("reason", "")),
        funnel=funnel,
        stability=dict(scored),
    )


def colors_from_catalog(
    catalog: pd.DataFrame, source_ids: Iterable
) -> dict[int, float]:
    """BP-RP per source, from whichever column the catalogue actually carries.

    Colour feeds the stability score's colour-match term, so a missing column
    must yield NaN rather than 0.0 — a fabricated 0.0 reads as "same colour as
    the target" and promotes the wrong stars.
    """
    if catalog is None or catalog.empty or "source_id" not in catalog.columns:
        return {}
    frame = catalog.copy()
    frame["source_id"] = pd.to_numeric(frame["source_id"], errors="coerce")
    frame = frame[frame["source_id"].notna()]
    frame["source_id"] = frame["source_id"].astype("int64")

    color = pd.Series(np.nan, index=frame.index, dtype="float64")
    for col in ("color_gr", "bp_rp", "BP_RP"):
        if col in frame.columns:
            color = color.fillna(pd.to_numeric(frame[col], errors="coerce"))
    if color.isna().any():
        bp = pd.Series(np.nan, index=frame.index, dtype="float64")
        rp = pd.Series(np.nan, index=frame.index, dtype="float64")
        for col in ("gaia_BP", "gaia_bp", "phot_bp_mean_mag"):
            if col in frame.columns:
                bp = bp.fillna(pd.to_numeric(frame[col], errors="coerce"))
        for col in ("gaia_RP", "gaia_rp", "phot_rp_mean_mag"):
            if col in frame.columns:
                rp = rp.fillna(pd.to_numeric(frame[col], errors="coerce"))
        color = color.fillna(bp - rp)

    wanted = _as_id_set(source_ids)
    lookup = dict(zip(frame["source_id"].tolist(), color.tolist()))
    return {sid: float(lookup[sid]) for sid in wanted if sid in lookup}


def write_screening_report(result_dir, result: ScreeningResult) -> Path:
    """Persist the per-star verdict beside the selection it produced."""
    from apex.utils.step_paths_lc import step8_selection_dir

    out_dir = Path(step8_selection_dir(result_dir))
    out_dir.mkdir(parents=True, exist_ok=True)
    safe = "".join(ch if ch.isalnum() or ch in "._-" else "_"
                   for ch in str(result.filter_key)) or "all"
    path = out_dir / f"comparison_screening_{safe}.tsv"

    report = result.report.copy()
    selected = set(result.selected_ids)
    report["role"] = ""
    report.loc[report["star_id"] == result.target_id, "role"] = "target"
    report.loc[report["star_id"].isin(selected), "role"] = "comparison"
    if result.check_id is not None:
        report.loc[report["star_id"] == int(result.check_id), "role"] = "check"
    report["in_pool"] = report["star_id"].isin(set(result.candidate_ids))

    metrics = result.metrics
    if isinstance(metrics, pd.DataFrame) and not metrics.empty:
        key = "star_id" if "star_id" in metrics.columns else None
        if key is not None:
            keep = [c for c in metrics.columns if c != key]
            report = report.merge(
                metrics[[key] + keep], on="star_id", how="left", suffixes=("", "_metric")
            )
    report.to_csv(path, sep="\t", index=False, encoding="utf-8")
    return path
