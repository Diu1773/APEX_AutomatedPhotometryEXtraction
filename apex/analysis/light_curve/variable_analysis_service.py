"""Headless advanced analysis for a Main-workflow light-curve release.

This module deliberately does not perform comparison/check-star QC or a broad
period search. Those decisions belong to the LC workflow and arrive through a
``ValidatedLightCurveBundle``. The service only refines an adopted candidate
and fits the requested single- or multi-mode model.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

import numpy as np
from scipy.optimize import minimize_scalar

from apex.analysis.light_curve.period_alias_service import fit_multimode_model
from apex.analysis.light_curve.period_io_service import load_period_lightcurve_csv
from apex.analysis.light_curve.variable_analysis_contract import (
    ReviewRequired,
    VariableAnalysisRequest,
    VariableAnalysisResult,
    compute_file_fingerprint,
)


ProgressCallback = Callable[[str, float], None]
CancelCallback = Callable[[], bool]


class AnalysisCancelled(RuntimeError):
    pass


def _report(callback: ProgressCallback | None, stage: str, fraction: float) -> None:
    if callback is not None:
        callback(str(stage), float(np.clip(fraction, 0.0, 1.0)))


def _check_cancelled(callback: CancelCallback | None) -> None:
    if callback is not None and bool(callback()):
        raise AnalysisCancelled("Advanced variable analysis was cancelled.")


def _window_frequency_offsets(alias_analysis: dict[str, Any]) -> list[float]:
    offsets: list[float] = []
    for row in alias_analysis.get("window_peaks", []):
        raw_value = row.get("freq_cd") if isinstance(row, dict) else row
        try:
            value = float(raw_value)
        except (TypeError, ValueError):
            continue
        if np.isfinite(value) and value > 0:
            offsets.append(value)
    return offsets


def fit_fixed_period_fourier(
    time: np.ndarray,
    mag: np.ndarray,
    period: float,
    harmonics: int,
    mag_err: np.ndarray | None = None,
) -> dict[str, Any]:
    """Fit a weighted cosine-series model at one fixed period."""
    t0 = np.asarray(time, dtype=float)
    y0 = np.asarray(mag, dtype=float)
    if len(t0) != len(y0):
        raise ValueError("time and mag must have the same length")
    if not np.isfinite(period) or period <= 0:
        raise ValueError("period must be positive and finite")

    n_harmonics = max(int(harmonics), 1)
    min_points = max(8, 2 * n_harmonics + 3)
    mask = np.isfinite(t0) & np.isfinite(y0)
    dy0 = None if mag_err is None else np.asarray(mag_err, dtype=float)
    if dy0 is not None:
        if len(dy0) != len(t0):
            raise ValueError("mag_err must have the same length as time")
        mask &= np.isfinite(dy0) & (dy0 > 0)

    t = t0[mask]
    y = y0[mask]
    dy = dy0[mask] if dy0 is not None else None
    if len(t) < min_points:
        raise ValueError(
            f"Not enough valid points for a {n_harmonics}-harmonic fit "
            f"({len(t)} < {min_points})."
        )

    order = np.argsort(t)
    t = t[order]
    y = y[order]
    dy = dy[order] if dy is not None else None
    time_ref = float(np.min(t))
    tau = t - time_ref
    omega = 2.0 * np.pi / float(period)
    columns = [np.ones(len(t), dtype=float)]
    for harmonic in range(1, n_harmonics + 1):
        columns.extend(
            [
                np.cos(harmonic * omega * tau),
                np.sin(harmonic * omega * tau),
            ]
        )
    design = np.column_stack(columns)

    if dy is None:
        weights = np.ones(len(t), dtype=float)
    else:
        weights = 1.0 / np.clip(dy, 1e-8, None) ** 2
        weight_scale = float(np.nanmedian(weights))
        if np.isfinite(weight_scale) and weight_scale > 0:
            weights /= weight_scale
    root_weights = np.sqrt(weights)
    coeff, _, rank, singular = np.linalg.lstsq(
        design * root_weights[:, None],
        y * root_weights,
        rcond=None,
    )
    if int(rank) < design.shape[1]:
        raise ValueError("The fixed-period Fourier fit is rank deficient.")
    singular = np.asarray(singular, dtype=float)
    condition = (
        float(singular[0] / singular[-1])
        if len(singular) >= 2 and singular[-1] > 0
        else np.inf
    )
    if not np.isfinite(condition) or condition > 1e12:
        raise ValueError(f"The fixed-period Fourier fit is ill conditioned ({condition:.3e}).")

    coeff = np.asarray(coeff, dtype=float)
    model = design @ coeff
    residual = y - model
    wrss = float(np.sum(weights * residual**2))
    rss = float(np.sum(residual**2))
    n_points = int(len(t))
    n_params = int(len(coeff))
    return {
        "time": t,
        "mag": y,
        "mag_err": dy,
        "period": float(period),
        "harmonics": n_harmonics,
        "coeff": coeff,
        "time_ref": time_ref,
        "model": model,
        "residual": residual,
        "n_points": n_points,
        "n_params": n_params,
        "rank": int(rank),
        "rmse": float(np.sqrt(rss / n_points)),
        "wrms": float(np.sqrt(wrss / n_points)),
        "rss": rss,
        "wrss": wrss,
        "bic": float(n_points * np.log(max(wrss / n_points, 1e-30)) + n_params * np.log(n_points)),
        "design_condition": condition,
        "weighted": dy is not None,
    }


def evaluate_fixed_period_fourier(
    time: np.ndarray,
    fit_result: dict[str, Any],
) -> np.ndarray:
    t = np.asarray(time, dtype=float)
    coeff = np.asarray(fit_result["coeff"], dtype=float)
    period = float(fit_result["period"])
    tau = t - float(fit_result["time_ref"])
    model = np.full(len(t), float(coeff[0]), dtype=float)
    omega = 2.0 * np.pi / period
    coefficient_index = 1
    for harmonic in range(1, int(fit_result["harmonics"]) + 1):
        model += float(coeff[coefficient_index]) * np.cos(harmonic * omega * tau)
        model += float(coeff[coefficient_index + 1]) * np.sin(harmonic * omega * tau)
        coefficient_index += 2
    return model


def fourier_shape_parameters(coeff: np.ndarray) -> dict[str, Any]:
    """Return convention-explicit shape terms for A*cos(k*omega*t + phi)."""
    values = np.asarray(coeff, dtype=float)
    if len(values) < 3 or (len(values) - 1) % 2:
        raise ValueError("Invalid Fourier coefficient vector")
    cos_coeff = values[1::2]
    sin_coeff = values[2::2]
    amplitudes = np.hypot(cos_coeff, sin_coeff)
    phases = np.arctan2(-sin_coeff, cos_coeff)

    def wrapped(value: float) -> float:
        return float((value + np.pi) % (2.0 * np.pi) - np.pi)

    return {
        "convention": "A_k*cos(k*omega*(t-time_ref)+phi_k)",
        "amplitudes": amplitudes,
        "phases": phases,
        "r21": float(amplitudes[1] / amplitudes[0])
        if len(amplitudes) > 1 and amplitudes[0] > 0
        else np.nan,
        "r31": float(amplitudes[2] / amplitudes[0])
        if len(amplitudes) > 2 and amplitudes[0] > 0
        else np.nan,
        "phi21": wrapped(float(phases[1] - 2.0 * phases[0]))
        if len(phases) > 1
        else np.nan,
        "phi31": wrapped(float(phases[2] - 3.0 * phases[0]))
        if len(phases) > 2
        else np.nan,
    }


def refine_period_local(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: np.ndarray | None,
    center_period: float,
    *,
    filter_values: np.ndarray | None = None,
    alias_frequency_offsets: list[float] | None = None,
    bootstrap_resamples: int = 300,
    harmonics: int = 4,
    random_seed: int | None = None,
    progress_callback: ProgressCallback | None = None,
    cancel_requested: CancelCallback | None = None,
) -> dict[str, Any]:
    """Refine one alias peak with a per-filter multi-harmonic objective.

    A one-harmonic Lomb-Scargle peak can shift for asymmetric pulsator light
    curves. Each filter therefore gets its own Fourier coefficients while all
    filters share the trial period. Bootstrap trials retain timestamps and
    resample residuals within each filter.
    """
    t0 = np.asarray(time, dtype=float)
    y0 = np.asarray(mag, dtype=float)
    mask = np.isfinite(t0) & np.isfinite(y0)
    dy0 = None if mag_err is None else np.asarray(mag_err, dtype=float)
    if dy0 is not None:
        if len(dy0) != len(t0):
            raise ValueError("mag_err must have the same length as time")
        mask &= np.isfinite(dy0) & (dy0 > 0)
    if filter_values is None:
        filter0 = np.full(len(t0), "all", dtype=str)
    else:
        filter0 = np.asarray(filter_values, dtype=str)
        if len(filter0) != len(t0):
            raise ValueError("filter_values must have the same length as time")
    t = t0[mask]
    y = y0[mask]
    dy = dy0[mask] if dy0 is not None else None
    filters = filter0[mask]
    if len(t) < 10:
        raise ValueError("At least 10 valid points are required for local refinement.")
    baseline = float(np.ptp(t))
    if not np.isfinite(baseline) or baseline <= 0:
        raise ValueError("The observation baseline must be positive.")
    if not np.isfinite(center_period) or center_period <= 0:
        raise ValueError("The center period must be positive and finite.")

    center_frequency = 1.0 / float(center_period)
    half_width = min(1.0 / baseline, 0.45)
    clean_alias_offsets = [
        float(offset)
        for offset in (alias_frequency_offsets or [])
        if np.isfinite(offset) and float(offset) > 0
    ]
    if clean_alias_offsets:
        half_width = min(half_width, 0.45 * min(clean_alias_offsets))
    lower = max(center_frequency - half_width, max(center_frequency * 0.25, 1e-8))
    upper = center_frequency + half_width
    frequency = np.linspace(lower, upper, 1001)
    filter_groups: list[tuple[str, np.ndarray, int]] = []
    for filter_name in dict.fromkeys(filters.tolist()):
        indices = np.flatnonzero(filters == filter_name)
        group_harmonics = min(
            max(int(harmonics), 1),
            max(1, (len(indices) - 3) // 2),
        )
        if len(indices) >= max(8, 2 * group_harmonics + 3):
            filter_groups.append((str(filter_name), indices, group_harmonics))
    if not filter_groups:
        raise ValueError("No filter has enough valid points for local refinement.")

    def objective(trial_frequency: float, values: np.ndarray) -> float:
        if not np.isfinite(trial_frequency) or trial_frequency <= 0:
            return np.inf
        total = 0.0
        for _, indices, group_harmonics in filter_groups:
            fit = fit_fixed_period_fourier(
                t[indices],
                values[indices],
                1.0 / float(trial_frequency),
                group_harmonics,
                mag_err=dy[indices] if dy is not None else None,
            )
            total += float(fit["wrss"])
        return total

    null_wrss = 0.0
    for _, indices, _ in filter_groups:
        values = y[indices]
        if dy is None:
            weights = np.ones(len(indices), dtype=float)
        else:
            weights = 1.0 / np.clip(dy[indices], 1e-8, None) ** 2
            scale = float(np.nanmedian(weights))
            if np.isfinite(scale) and scale > 0:
                weights /= scale
        center = float(np.sum(weights * values) / np.sum(weights))
        null_wrss += float(np.sum(weights * (values - center) ** 2))

    _check_cancelled(cancel_requested)
    _report(progress_callback, "Evaluating the multi-harmonic local period grid", 0.10)
    objective_grid = np.empty(len(frequency), dtype=float)
    for index, trial_frequency in enumerate(frequency):
        _check_cancelled(cancel_requested)
        objective_grid[index] = objective(float(trial_frequency), y)
        if index and index % 100 == 0:
            _report(
                progress_callback,
                f"Local period grid {index}/{len(frequency)}",
                0.10 + 0.18 * (index / len(frequency)),
            )
    power = 1.0 - objective_grid / max(null_wrss, 1e-30)
    peak_index = int(np.nanargmax(power))
    minima = np.flatnonzero(
        (power[1:-1] <= power[:-2]) & (power[1:-1] <= power[2:])
    ) + 1
    left_candidates = minima[minima < peak_index]
    right_candidates = minima[minima > peak_index]
    left_index = int(left_candidates[-1]) if len(left_candidates) else 0
    right_index = int(right_candidates[0]) if len(right_candidates) else len(frequency) - 1
    if right_index - left_index < 4:
        left_index = max(0, peak_index - 5)
        right_index = min(len(frequency) - 1, peak_index + 5)
    local_lower = float(frequency[left_index])
    local_upper = float(frequency[right_index])
    optimum = minimize_scalar(
        lambda trial_frequency: objective(float(trial_frequency), y),
        bounds=(local_lower, local_upper),
        method="bounded",
        options={"xatol": max((local_upper - local_lower) * 1e-8, 1e-10)},
    )
    refined_frequency = (
        float(optimum.x)
        if optimum.success and np.isfinite(optimum.fun)
        else float(frequency[peak_index])
    )
    refined_period = 1.0 / refined_frequency

    template_model = np.full(len(t), np.nan, dtype=float)
    residual = np.full(len(t), np.nan, dtype=float)
    for _, indices, group_harmonics in filter_groups:
        template = fit_fixed_period_fourier(
            t[indices],
            y[indices],
            refined_period,
            group_harmonics,
            mag_err=dy[indices] if dy is not None else None,
        )
        group_model = evaluate_fixed_period_fourier(t[indices], template)
        group_residual = y[indices] - group_model
        group_residual -= float(np.nanmedian(group_residual))
        template_model[indices] = group_model
        residual[indices] = group_residual
    rng = np.random.default_rng(random_seed)
    n_bootstrap = max(int(bootstrap_resamples), 0)
    bootstrap_periods: list[float] = []
    for index in range(n_bootstrap):
        _check_cancelled(cancel_requested)
        simulated_mag = template_model.copy()
        for _, indices, _ in filter_groups:
            simulated_mag[indices] += rng.choice(
                residual[indices], len(indices), replace=True
            )
        bootstrap_optimum = minimize_scalar(
            lambda trial_frequency: objective(float(trial_frequency), simulated_mag),
            bounds=(local_lower, local_upper),
            method="bounded",
            options={"xatol": max((local_upper - local_lower) * 1e-7, 1e-9)},
        )
        if bootstrap_optimum.success and np.isfinite(bootstrap_optimum.fun):
            bootstrap_periods.append(1.0 / float(bootstrap_optimum.x))
        if index == 0 or (index + 1) % max(1, n_bootstrap // 20) == 0:
            fraction = 0.30 + 0.50 * ((index + 1) / max(n_bootstrap, 1))
            _report(
                progress_callback,
                f"Local residual bootstrap {index + 1}/{n_bootstrap}",
                fraction,
            )

    bootstrap_array = np.asarray(bootstrap_periods, dtype=float)
    local_sigma = None
    if len(bootstrap_array) >= 2:
        median = float(np.nanmedian(bootstrap_array))
        mad = float(np.nanmedian(np.abs(bootstrap_array - median)))
        if np.isfinite(mad) and mad > 0:
            keep = np.abs(bootstrap_array - median) <= 5.0 * 1.4826 * mad
            values = bootstrap_array[keep] if np.count_nonzero(keep) >= 5 else bootstrap_array
        else:
            values = bootstrap_array
        local_sigma = float(np.nanstd(values, ddof=1))

    return {
        "center_period": float(center_period),
        "refined_period": float(refined_period),
        "local_period_sigma": local_sigma,
        "uncertainty_scope": "local_conditional_on_selected_alias",
        "refinement_method": "per_filter_multi_harmonic_fourier",
        "frequency_half_width_cd": float(half_width),
        "local_frequency_bounds_cd": [local_lower, local_upper],
        "fine_frequency": frequency,
        "fine_power": power,
        "bootstrap_periods": bootstrap_array,
        "bootstrap_resamples_completed": int(len(bootstrap_array)),
        "filter_harmonics": {
            filter_name: int(group_harmonics)
            for filter_name, _, group_harmonics in filter_groups
        },
    }


def _resolve_source_path(request: VariableAnalysisRequest) -> Path:
    source = Path(request.bundle.source_file)
    if source.exists():
        return source
    workspace = Path(request.bundle.workspace_dir)
    alternatives = [workspace / source, workspace / source.name]
    for alternative in alternatives:
        if alternative.exists():
            return alternative
    raise FileNotFoundError(f"Released light-curve source is missing: {source}")


def _load_released_series(request: VariableAnalysisRequest) -> dict[str, Any]:
    bundle = request.bundle
    source = _resolve_source_path(request)
    loaded = load_period_lightcurve_csv(source, bundle.analysis_filter, bundle.target_id)
    loaded_raw_col = str(loaded.get("col_raw") or "")
    loaded_corr_col = str(loaded.get("col_corr") or "")
    requested_col = str(bundle.mag_col or "")
    use_corrected = bool(
        loaded.get("mag_corr") is not None
        and (
            requested_col == loaded_corr_col
            or str(bundle.series_mode).lower() in {"corr", "corrected"}
        )
    )
    selected_column = loaded_corr_col if use_corrected else loaded_raw_col
    if requested_col and requested_col not in {loaded_raw_col, loaded_corr_col}:
        raise ValueError(
            f"Released magnitude column {requested_col!r} is not present in the loaded "
            f"series ({loaded_raw_col!r}, {loaded_corr_col!r})."
        )

    mag = loaded.get("mag_corr") if use_corrected else loaded.get("mag_raw")
    time = np.asarray(loaded.get("time", []), dtype=float)
    mag = np.asarray(mag, dtype=float)
    mag_err_raw = loaded.get("mag_err")
    mag_err = None if mag_err_raw is None else np.asarray(mag_err_raw, dtype=float)
    night_id = np.asarray(loaded.get("night_id", []), dtype=str)
    filters = np.asarray(loaded.get("filter_values", []), dtype=str)
    mask = np.isfinite(time) & np.isfinite(mag)
    if mag_err is not None:
        mask &= np.isfinite(mag_err) & (mag_err > 0)
    time = time[mask]
    mag = mag[mask]
    mag_err = mag_err[mask] if mag_err is not None else None
    night_id = night_id[mask]
    filters = filters[mask]
    if len(time) < 10:
        raise ValueError("The released target series has fewer than 10 valid points.")

    order = np.argsort(time)
    time = time[order]
    mag = mag[order]
    mag_err = mag_err[order] if mag_err is not None else None
    night_id = night_id[order]
    filters = filters[order]
    signature = dict(bundle.input_signature or {})
    mismatches: list[str] = []
    if signature.get("n_points") is not None and int(signature["n_points"]) != len(time):
        mismatches.append(f"n_points {len(time)} != {int(signature['n_points'])}")
    for key, actual in (("time_min", float(time[0])), ("time_max", float(time[-1]))):
        expected = signature.get(key)
        if expected is not None and not np.isclose(float(expected), actual, rtol=0.0, atol=1e-9):
            mismatches.append(f"{key} {actual:.12f} != {float(expected):.12f}")
    current_fingerprint = compute_file_fingerprint(source)
    fingerprint_keys = (
        ("source_sha256",)
        if signature.get("source_sha256")
        else ("source_size", "source_mtime_ns")
    )
    for key in fingerprint_keys:
        expected = signature.get(key)
        actual = current_fingerprint[key]
        if expected is not None and str(expected) != str(actual):
            mismatches.append(f"{key} changed")
    if mismatches and bundle.release_status in {"APPROVED", "OVERRIDDEN"}:
        raise ValueError(
            "The released source changed after Main-workflow validation: " + "; ".join(mismatches)
        )

    return {
        "source_file": str(source),
        "time": time,
        "mag": mag,
        "mag_err": mag_err,
        "night_id": night_id,
        "filter_values": filters,
        "selected_mag_column": selected_column,
        "signature_warnings": mismatches,
        "filter_alignment": loaded.get("filter_alignment", "none"),
    }


def _alias_review(request: VariableAnalysisRequest) -> VariableAnalysisResult | None:
    alias = dict(request.bundle.alias_analysis or {})
    status = str(alias.get("status", "UNASSESSED")).upper()
    if status == "RESOLVED" or request.adopted_period_override is not None:
        return None
    candidates = [dict(row) for row in alias.get("candidates", []) if isinstance(row, dict)]
    return VariableAnalysisResult(
        status="REVIEW_REQUIRED",
        adopted_period=float(request.bundle.adopted_period),
        diagnostics={"alias_analysis": alias},
        review=ReviewRequired(
            code="ALIAS_SELECTION_REQUIRED",
            message=(
                f"Main workflow alias status is {status}. Select one period candidate "
                "explicitly before local refinement."
            ),
            candidates=candidates,
            allowed_actions=["adopt_period_candidate", "cancel"],
        ),
    )


def _fit_single_branch(
    series: dict[str, Any],
    request: VariableAnalysisRequest,
    refined_period: float,
    progress_callback: ProgressCallback | None,
    cancel_requested: CancelCallback | None,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    filters = np.asarray(series["filter_values"], dtype=str)
    per_filter: dict[str, dict[str, Any]] = {}
    unique_filters = list(dict.fromkeys(filters.tolist()))
    for index, filter_name in enumerate(unique_filters):
        _check_cancelled(cancel_requested)
        selected = filters == filter_name
        try:
            fit = fit_fixed_period_fourier(
                series["time"][selected],
                series["mag"][selected],
                refined_period,
                request.single_harmonics,
                mag_err=(series["mag_err"][selected] if series["mag_err"] is not None else None),
            )
            fit["shape"] = fourier_shape_parameters(fit["coeff"])
            per_filter[str(filter_name)] = fit
        except ValueError as exc:
            per_filter[str(filter_name)] = {"error": str(exc), "n_points": int(np.sum(selected))}
        _report(
            progress_callback,
            f"Fitting single-mode Fourier model: {filter_name}",
            0.82 + 0.15 * ((index + 1) / max(len(unique_filters), 1)),
        )

    successful = [fit for fit in per_filter.values() if "error" not in fit]
    if not successful:
        raise ValueError("No filter has enough data for the requested single-mode model.")
    residual = np.concatenate([np.asarray(fit["residual"], dtype=float) for fit in successful])
    n_points = int(sum(int(fit["n_points"]) for fit in successful))
    model = {
        "kind": "per_filter_fourier",
        "period": float(refined_period),
        "harmonics": int(request.single_harmonics),
        "n_filters": int(len(successful)),
        "n_points": n_points,
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "bic_sum": float(sum(float(fit["bic"]) for fit in successful)),
    }
    return model, per_filter


def _fit_multi_branch(
    series: dict[str, Any],
    request: VariableAnalysisRequest,
    periods: list[float],
    progress_callback: ProgressCallback | None,
    cancel_requested: CancelCallback | None,
) -> tuple[dict[str, Any], dict[str, dict[str, Any]]]:
    filters = np.asarray(series["filter_values"], dtype=str)
    per_filter: dict[str, dict[str, Any]] = {}
    unique_filters = list(dict.fromkeys(filters.tolist()))
    for index, filter_name in enumerate(unique_filters):
        _check_cancelled(cancel_requested)
        selected = filters == filter_name
        minimum_points = max(
            10,
            2 * (1 + 2 * len(periods) * int(request.multimode_harmonics)),
        )
        if int(np.count_nonzero(selected)) < minimum_points:
            per_filter[str(filter_name)] = {
                "error": (
                    "Not enough points for a stable multi-mode fit "
                    f"({int(np.count_nonzero(selected))} < {minimum_points})."
                ),
                "n_points": int(np.count_nonzero(selected)),
            }
            _report(
                progress_callback,
                f"Skipping underconstrained multi-mode filter: {filter_name}",
                0.82 + 0.15 * ((index + 1) / max(len(unique_filters), 1)),
            )
            continue
        try:
            fit = fit_multimode_model(
                series["time"][selected],
                series["mag"][selected],
                series["mag_err"][selected] if series["mag_err"] is not None else None,
                periods=periods,
                harmonics=request.multimode_harmonics,
                night_id=series["night_id"][selected],
                include_night_offsets=request.include_night_offsets,
            )
            per_filter[str(filter_name)] = fit
        except ValueError as exc:
            per_filter[str(filter_name)] = {"error": str(exc), "n_points": int(np.sum(selected))}
        _report(
            progress_callback,
            f"Fitting multi-mode model: {filter_name}",
            0.82 + 0.15 * ((index + 1) / max(len(unique_filters), 1)),
        )

    successful = [fit for fit in per_filter.values() if "error" not in fit]
    if not successful:
        raise ValueError("No filter has enough data for the requested multi-mode model.")
    residual = np.concatenate([np.asarray(fit["residual"], dtype=float) for fit in successful])
    n_points = int(sum(int(fit["n_points"]) for fit in successful))
    model = {
        "kind": "per_filter_multimode",
        "periods": [float(period) for period in periods],
        "harmonics_per_mode": int(request.multimode_harmonics),
        "include_night_offsets": bool(request.include_night_offsets),
        "n_filters": int(len(successful)),
        "n_points": n_points,
        "rmse": float(np.sqrt(np.mean(residual**2))),
        "bic_sum": float(sum(float(fit["bic"]) for fit in successful)),
    }
    return model, per_filter


def run_variable_analysis(
    request: VariableAnalysisRequest,
    *,
    progress_callback: ProgressCallback | None = None,
    cancel_requested: CancelCallback | None = None,
) -> VariableAnalysisResult:
    """Execute one advanced-analysis request without importing Qt."""
    try:
        if not request.bundle.can_launch:
            return VariableAnalysisResult(
                status="FAILED",
                error=(
                    "Main workflow blocked advanced analysis: "
                    f"{request.bundle.release_message}"
                ),
            )
        review = _alias_review(request)
        if review is not None:
            return review

        _check_cancelled(cancel_requested)
        _report(progress_callback, "Loading the released target series", 0.03)
        series = _load_released_series(request)
        adopted_period = float(
            request.adopted_period_override
            if request.adopted_period_override is not None
            else request.bundle.adopted_period
        )
        if not np.isfinite(adopted_period) or adopted_period <= 0:
            raise ValueError("The released adopted period is invalid.")

        refinement = refine_period_local(
            series["time"],
            series["mag"],
            series["mag_err"],
            adopted_period,
            filter_values=series["filter_values"],
            alias_frequency_offsets=_window_frequency_offsets(
                request.bundle.alias_analysis
            ),
            bootstrap_resamples=request.bootstrap_resamples,
            harmonics=request.refinement_harmonics,
            random_seed=request.random_seed,
            progress_callback=progress_callback,
            cancel_requested=cancel_requested,
        )
        refined_period = float(refinement["refined_period"])
        diagnostic = dict(request.bundle.multimode_diagnostic or {})
        branch = request.analysis_branch
        if branch == "auto":
            branch = (
                "multi"
                if str(diagnostic.get("status", "")).upper() == "MULTIMODE-SUSPECT"
                else "single"
            )

        _check_cancelled(cancel_requested)
        if branch == "multi":
            candidate_period = (
                request.secondary_period_override
                if request.secondary_period_override is not None
                else diagnostic.get("candidate_period")
            )
            try:
                candidate_period = float(candidate_period)
            except (TypeError, ValueError):
                candidate_period = np.nan
            if (
                not np.isfinite(candidate_period)
                or candidate_period <= 0
                or abs(candidate_period - refined_period) / refined_period <= 1e-5
            ):
                return VariableAnalysisResult(
                    status="REVIEW_REQUIRED",
                    branch="multi",
                    adopted_period=adopted_period,
                    refined_period=refined_period,
                    local_period_sigma=refinement["local_period_sigma"],
                    refinement=refinement,
                    diagnostics={"multimode_diagnostic": diagnostic},
                    review=ReviewRequired(
                        code="SECOND_MODE_REQUIRED",
                        message="Multi-mode analysis needs an explicit independent second period.",
                        candidates=[
                            dict(row)
                            for row in diagnostic.get("candidates", [])
                            if isinstance(row, dict)
                        ],
                        allowed_actions=["select_second_mode", "use_single_mode", "cancel"],
                    ),
                )
            periods = [refined_period, candidate_period]
            model, per_filter = _fit_multi_branch(
                series,
                request,
                periods,
                progress_callback,
                cancel_requested,
            )
        else:
            model, per_filter = _fit_single_branch(
                series,
                request,
                refined_period,
                progress_callback,
                cancel_requested,
            )

        _report(progress_callback, "Advanced analysis complete", 1.0)
        baseline_days = float(np.ptp(series["time"]))
        observed_cycles = baseline_days / refined_period
        return VariableAnalysisResult(
            status="COMPLETE",
            branch=branch,
            adopted_period=adopted_period,
            refined_period=refined_period,
            local_period_sigma=refinement["local_period_sigma"],
            data_summary={
                "source_file": series["source_file"],
                "target_id": int(request.bundle.target_id),
                "analysis_filter": request.bundle.analysis_filter,
                "selected_mag_column": series["selected_mag_column"],
                "n_points": int(len(series["time"])),
                "n_nights": int(len(np.unique(series["night_id"]))),
                "filters": list(dict.fromkeys(series["filter_values"].tolist())),
                "baseline_days": baseline_days,
                "observed_cycles": float(observed_cycles),
                "filter_alignment": series["filter_alignment"],
                "signature_warnings": list(series["signature_warnings"]),
                "main_qc_release": request.bundle.release_status,
            },
            refinement=refinement,
            model=model,
            per_filter_models=per_filter,
            diagnostics={
                "alias_analysis": dict(request.bundle.alias_analysis or {}),
                "multimode_diagnostic": diagnostic,
                "local_uncertainty_excludes_alias_choice": True,
                "limited_cycle_coverage": bool(observed_cycles < 3.0),
                "main_qc_recomputed": False,
            },
        )
    except AnalysisCancelled as exc:
        return VariableAnalysisResult(status="CANCELLED", error=str(exc))
    except Exception as exc:
        return VariableAnalysisResult(status="FAILED", error=str(exc))
