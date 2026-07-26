"""Sampling-window aware period and multi-mode analysis.

The functions in this module deliberately keep alias ambiguity separate from
the local uncertainty of a periodogram peak.  Missing temporal coverage cannot
be repaired by a different periodogram, so callers receive ranked candidates
and an explicit RESOLVED / AMBIGUOUS / INSUFFICIENT status.
"""

from __future__ import annotations

from itertools import islice, product
from typing import Iterable, Optional, Sequence

import numpy as np
from astropy.timeseries import LombScargle
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks


def _finite_time(time: np.ndarray) -> np.ndarray:
    t = np.asarray(time, dtype=float)
    return np.sort(t[np.isfinite(t)])


def infer_night_ids(time: np.ndarray, gap_days: float = 0.5) -> np.ndarray:
    """Group observations into sessions when no explicit night label exists."""
    t = np.asarray(time, dtype=float)
    labels = np.full(len(t), "invalid", dtype=object)
    finite_idx = np.flatnonzero(np.isfinite(t))
    if len(finite_idx) == 0:
        return labels.astype(str)
    order = finite_idx[np.argsort(t[finite_idx])]
    session = 1
    previous = None
    threshold = max(float(gap_days), 1e-6)
    for idx in order:
        value = float(t[idx])
        if previous is not None and value - previous > threshold:
            session += 1
        labels[idx] = f"session-{session}"
        previous = value
    return labels.astype(str)


def _coerce_series(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray] = None,
    night_id: Optional[np.ndarray] = None,
) -> tuple[np.ndarray, np.ndarray, Optional[np.ndarray], np.ndarray]:
    t0 = np.asarray(time, dtype=float)
    y0 = np.asarray(mag, dtype=float)
    if len(t0) != len(y0):
        raise ValueError("time and mag must have the same length")

    mask = np.isfinite(t0) & np.isfinite(y0)
    dy0 = None if mag_err is None else np.asarray(mag_err, dtype=float)
    if dy0 is not None:
        if len(dy0) != len(t0):
            raise ValueError("mag_err must have the same length as time")
        mask &= np.isfinite(dy0) & (dy0 > 0)

    if night_id is None:
        nid0 = np.full(len(t0), "night-1", dtype=object)
    else:
        nid0 = np.asarray(night_id, dtype=object)
        if len(nid0) != len(t0):
            raise ValueError("night_id must have the same length as time")

    t = t0[mask]
    y = y0[mask]
    dy = dy0[mask] if dy0 is not None else None
    nid = nid0[mask].astype(str)
    order = np.argsort(t)
    return t[order], y[order], (dy[order] if dy is not None else None), nid[order]


def compute_spectral_window(
    time: np.ndarray,
    max_frequency: float = 4.0,
    samples_per_peak: int = 20,
    max_peaks: int = 8,
) -> dict:
    """Return the normalized sampling window and its strongest side lobes."""
    t = _finite_time(time)
    if len(t) < 3:
        return {
            "frequency": np.array([], dtype=float),
            "power": np.array([], dtype=float),
            "peaks": [],
            "baseline_days": 0.0,
            "resolution_cd": np.nan,
        }

    baseline = float(np.ptp(t))
    if baseline <= 0:
        return {
            "frequency": np.array([], dtype=float),
            "power": np.array([], dtype=float),
            "peaks": [],
            "baseline_days": baseline,
            "resolution_cd": np.nan,
        }

    resolution = 1.0 / baseline
    step = max(resolution / max(int(samples_per_peak), 1), 5e-4)
    upper = max(float(max_frequency), step)
    n_freq = max(32, min(50000, int(np.ceil(upper / step))))
    freqs = np.linspace(step, upper, n_freq, dtype=float)
    tau = t - float(np.min(t))
    power = np.empty_like(freqs)

    batch_size = max(32, min(1024, int(2_000_000 / max(len(t), 1))))
    for start in range(0, len(freqs), batch_size):
        stop = min(start + batch_size, len(freqs))
        phase = 2.0 * np.pi * freqs[start:stop, None] * tau[None, :]
        window = np.mean(np.exp(1j * phase), axis=1)
        power[start:stop] = np.abs(window) ** 2

    min_alias_freq = max(0.05, 0.25 / baseline)
    min_distance = max(1, int(round(0.02 / max(step, 1e-9))))
    peak_idx, _ = find_peaks(
        power,
        height=0.05,
        prominence=0.01,
        distance=min_distance,
    )
    peak_idx = [int(i) for i in peak_idx if freqs[int(i)] >= min_alias_freq]
    peak_idx.sort(key=lambda i: float(power[i]), reverse=True)
    peaks = [
        {"freq_cd": float(freqs[i]), "power": float(power[i])}
        for i in peak_idx[: max(int(max_peaks), 1)]
    ]
    return {
        "frequency": freqs,
        "power": power,
        "peaks": peaks,
        "baseline_days": baseline,
        "resolution_cd": resolution,
    }


def _periodogram_peaks(
    frequency: np.ndarray,
    power: np.ndarray,
    max_peaks: int = 10,
) -> list[int]:
    f = np.asarray(frequency, dtype=float)
    p = np.asarray(power, dtype=float)
    valid = np.isfinite(f) & np.isfinite(p) & (f > 0)
    if not np.any(valid):
        return []
    valid_idx = np.flatnonzero(valid)
    fv = f[valid]
    pv = p[valid]
    best_local = int(np.argmax(pv))
    idx, _ = find_peaks(pv, height=0.08 * float(pv[best_local]))
    selected = {best_local, *(int(i) for i in idx)}
    ordered = sorted(selected, key=lambda i: float(pv[i]), reverse=True)
    return [int(valid_idx[i]) for i in ordered[: max(int(max_peaks), 1)]]


def build_alias_candidates(
    time: np.ndarray,
    frequency: np.ndarray,
    power: np.ndarray,
    min_period: float,
    max_period: float,
    window: Optional[dict] = None,
    max_candidates: int = 10,
) -> list[dict]:
    """Build signal candidates from periodogram peaks and window offsets."""
    t = _finite_time(time)
    f = np.asarray(frequency, dtype=float)
    p = np.asarray(power, dtype=float)
    if len(t) < 3 or len(f) == 0 or len(f) != len(p):
        return []

    baseline = max(float(np.ptp(t)), 1e-9)
    resolution = 1.0 / baseline
    window = window or compute_spectral_window(t)
    fmin, fmax = 1.0 / float(max_period), 1.0 / float(min_period)
    direct_idx = _periodogram_peaks(f, p, max_peaks=max_candidates)
    if not direct_idx:
        return []

    grid_step = float(np.nanmedian(np.diff(np.sort(f)))) if len(f) > 1 else resolution / 10.0
    snap_radius = max(3.0 * abs(grid_step), min(0.08, 0.12 * resolution))
    dedup_tol = max(2.0 * abs(grid_step), 1e-4)
    candidates: list[dict] = []

    def snap(candidate_f: float) -> tuple[float, float]:
        near = np.flatnonzero(np.abs(f - candidate_f) <= snap_radius)
        if len(near):
            idx = int(near[np.argmax(p[near])])
        else:
            idx = int(np.argmin(np.abs(f - candidate_f)))
        return float(f[idx]), float(p[idx])

    def add(candidate_f: float, origin: str) -> None:
        if not np.isfinite(candidate_f) or not (fmin <= candidate_f <= fmax):
            return
        sf, sp = snap(float(candidate_f))
        for row in candidates:
            if abs(float(row["freq_cd"]) - sf) <= dedup_tol:
                origins = row.setdefault("origins", [])
                if origin not in origins:
                    origins.append(origin)
                row["ls_power"] = max(float(row["ls_power"]), sp)
                return
        candidates.append(
            {
                "freq_cd": sf,
                "period": 1.0 / sf,
                "ls_power": sp,
                "origins": [origin],
            }
        )

    for rank, idx in enumerate(direct_idx, start=1):
        add(float(f[idx]), f"periodogram-peak-{rank}")

    seed_f = float(f[direct_idx[0]])
    for peak in (window.get("peaks") or [])[:6]:
        offset = float(peak.get("freq_cd", np.nan))
        if not np.isfinite(offset) or offset <= 0:
            continue
        add(seed_f - offset, f"best-window-{offset:.5f}")
        add(seed_f + offset, f"best+window-{offset:.5f}")

    candidates.sort(key=lambda row: float(row["ls_power"]), reverse=True)
    best_power = max(float(candidates[0]["ls_power"]), 1e-30) if candidates else 1.0
    for rank, row in enumerate(candidates[:max_candidates], start=1):
        row["rank_by_power"] = rank
        row["relative_power"] = float(row["ls_power"]) / best_power
    return candidates[: max(int(max_candidates), 1)]


def _normalized_weights(mag_err: Optional[np.ndarray], n: int) -> np.ndarray:
    if mag_err is None:
        return np.ones(n, dtype=float)
    dy = np.asarray(mag_err, dtype=float)
    w = 1.0 / np.clip(dy, 1e-8, None) ** 2
    scale = float(np.nanmedian(w[np.isfinite(w) & (w > 0)]))
    if not np.isfinite(scale) or scale <= 0:
        return np.ones(n, dtype=float)
    return w / scale


def _solve_weighted_design(
    design: np.ndarray,
    values: np.ndarray,
    weights: np.ndarray,
) -> dict:
    """Solve a linear model and reject underconstrained numerical fits."""
    matrix = np.asarray(design, dtype=float)
    y = np.asarray(values, dtype=float)
    w = np.asarray(weights, dtype=float)
    n = int(len(y))
    k = int(matrix.shape[1])
    if n <= k:
        raise ValueError(f"Underconstrained fit: n_points={n}, n_params={k}.")

    root_w = np.sqrt(w)
    coeff, _, rank, singular = np.linalg.lstsq(
        matrix * root_w[:, None], y * root_w, rcond=None
    )
    singular = np.asarray(singular, dtype=float)
    condition = (
        float(singular[0] / singular[-1])
        if len(singular) >= 2 and singular[-1] > 0
        else np.inf
    )
    if int(rank) < k:
        raise ValueError(f"Rank-deficient fit: rank={int(rank)}, n_params={k}.")
    if not np.isfinite(condition) or condition > 1e12:
        raise ValueError(f"Ill-conditioned fit: condition={condition:.3e}.")

    coeff = np.asarray(coeff, dtype=float)
    model = matrix @ coeff
    residual = y - model
    wrss = float(np.sum(w * residual**2))
    rss = float(np.sum(residual**2))
    bic = float(n * np.log(max(wrss / n, 1e-30)) + k * np.log(n))
    return {
        "coeff": coeff,
        "model": model,
        "residual": residual,
        "rss": rss,
        "wrss": wrss,
        "rmse": float(np.sqrt(rss / n)),
        "wrms": float(np.sqrt(wrss / n)),
        "bic": bic,
        "n_points": n,
        "n_params": k,
        "rank": int(rank),
        "design_condition": condition,
    }


def _fit_frequency_model(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
    night_id: np.ndarray,
    freq_cd: float,
    harmonics: int,
    include_night_offsets: bool,
) -> dict:
    t = np.asarray(time, dtype=float)
    y = np.asarray(mag, dtype=float)
    nid = np.asarray(night_id, dtype=str)
    tau = t - float(np.min(t))
    columns = [np.ones(len(t), dtype=float)]
    labels = ["intercept"]

    unique_nights = list(dict.fromkeys(nid.tolist()))
    if include_night_offsets and len(unique_nights) > 1:
        for night in unique_nights[1:]:
            columns.append((nid == night).astype(float))
            labels.append(f"night:{night}")

    for harmonic in range(1, max(int(harmonics), 1) + 1):
        phase = 2.0 * np.pi * float(freq_cd) * harmonic * tau
        columns.extend([np.cos(phase), np.sin(phase)])
        labels.extend([f"cos:{harmonic}", f"sin:{harmonic}"])

    design = np.column_stack(columns)
    weights = _normalized_weights(mag_err, len(t))
    solved = _solve_weighted_design(design, y, weights)
    return {
        **solved,
        "labels": labels,
    }


def _refine_frequency_fit(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
    night_id: np.ndarray,
    initial_frequency: float,
    half_width: float,
    min_frequency: float,
    max_frequency: float,
    harmonics: int,
    include_night_offsets: bool,
) -> tuple[float, dict]:
    """Minimize model BIC locally without jumping to another alias family."""
    initial = float(np.clip(initial_frequency, min_frequency, max_frequency))
    best_fit = _fit_frequency_model(
        time,
        mag,
        mag_err,
        night_id,
        initial,
        harmonics,
        include_night_offsets,
    )
    lower = max(float(min_frequency), initial - float(half_width))
    upper = min(float(max_frequency), initial + float(half_width))
    if not upper > lower:
        return initial, best_fit

    def objective(frequency: float) -> float:
        return float(
            _fit_frequency_model(
                time,
                mag,
                mag_err,
                night_id,
                float(frequency),
                harmonics,
                include_night_offsets,
            )["bic"]
        )

    result = minimize_scalar(
        objective,
        bounds=(lower, upper),
        method="bounded",
        options={"xatol": max((upper - lower) * 1e-5, 1e-8)},
    )
    if result.success and np.isfinite(result.fun) and float(result.fun) < float(best_fit["bic"]):
        frequency = float(result.x)
        return frequency, _fit_frequency_model(
            time,
            mag,
            mag_err,
            night_id,
            frequency,
            harmonics,
            include_night_offsets,
        )
    return initial, best_fit


def phase_coverage_metrics(
    time: np.ndarray,
    period: float,
    n_bins: int = 20,
) -> dict:
    t = _finite_time(time)
    if len(t) == 0 or not np.isfinite(period) or period <= 0:
        return {"occupied_fraction": 0.0, "max_phase_gap": 1.0, "cycles": 0.0}
    phases = np.sort(((t - float(np.min(t))) / float(period)) % 1.0)
    wrapped = np.concatenate([phases, phases[:1] + 1.0])
    max_gap = float(np.max(np.diff(wrapped))) if len(phases) else 1.0
    occupied = np.unique(np.clip((phases * n_bins).astype(int), 0, n_bins - 1))
    return {
        "occupied_fraction": float(len(occupied) / max(int(n_bins), 1)),
        "max_phase_gap": max_gap,
        "cycles": float(np.ptp(t) / float(period)) if len(t) > 1 else 0.0,
    }


def _window_offsets(window_peaks: Optional[Sequence[dict | float]]) -> list[float]:
    offsets: list[float] = []
    for item in window_peaks or []:
        value = item.get("freq_cd") if isinstance(item, dict) else item
        try:
            freq = float(value)
        except Exception:
            continue
        if np.isfinite(freq) and freq > 0:
            offsets.append(freq)
    return offsets


def periods_are_window_aliases(
    period1: float,
    period2: float,
    window_peaks: Optional[Sequence[dict | float]],
    baseline_days: float = 1.0,
) -> bool:
    if period1 <= 0 or period2 <= 0:
        return False
    f1, f2 = 1.0 / float(period1), 1.0 / float(period2)
    delta = abs(f1 - f2)
    tol = max(0.02, min(0.12, 0.10 / max(float(baseline_days), 1.0)))
    return any(abs(delta - offset) <= tol for offset in _window_offsets(window_peaks))


def classify_frequency_relation(
    frequency: float,
    adopted_frequencies: Sequence[float],
    window_peaks: Optional[Sequence[dict | float]] = None,
    baseline_days: float = 1.0,
    max_harmonic: int = 6,
) -> tuple[str, str]:
    """Classify a frequency as duplicate, sampling alias, harmonic or new."""
    freq = float(frequency)
    previous = [float(v) for v in adopted_frequencies if np.isfinite(v) and v > 0]
    tol = max(0.01, min(0.08, 0.05 / max(float(baseline_days), 1.0)))
    for idx, prev in enumerate(previous):
        if abs(freq - prev) <= tol:
            return "duplicate", f"already adopted as M{idx + 1}"

    offsets = _window_offsets(window_peaks)
    alias_tol = max(0.02, min(0.12, 0.10 / max(float(baseline_days), 1.0)))
    for idx, prev in enumerate(previous):
        delta = abs(freq - prev)
        for offset in offsets:
            if abs(delta - offset) <= alias_tol:
                return "alias", f"window alias of M{idx + 1} ({offset:.4f} d^-1)"

    for idx, prev in enumerate(previous):
        for order in range(2, max(int(max_harmonic), 2) + 1):
            if abs(freq - order * prev) / max(abs(freq), 1e-12) < 0.01:
                return "harmonic", f"near {order}f(M{idx + 1})"
            if abs(prev - order * freq) / max(abs(prev), 1e-12) < 0.01:
                return "harmonic", f"subharmonic of M{idx + 1}"

    if len(previous) >= 2:
        for j in range(len(previous)):
            for k in range(j + 1, len(previous)):
                for sign in (-1.0, 1.0):
                    combo = previous[j] + sign * previous[k]
                    if combo > 0 and abs(freq - combo) / max(abs(freq), 1e-12) < 0.01:
                        op = "+" if sign > 0 else "-"
                        return "combination", f"near f(M{j + 1}) {op} f(M{k + 1})"
    return "new", "independent-frequency candidate"


def analyze_period_aliases(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
    night_id: Optional[np.ndarray],
    ls_frequency: np.ndarray,
    ls_power: np.ndarray,
    min_period: float,
    max_period: float,
    harmonics: int = 2,
    max_candidates: int = 10,
    n_injections: int = 8,
    random_seed: int = 12345,
    night_offset_policy: str = "diagnostic",
) -> dict:
    """Rank alias candidates and return a conservative resolution status."""
    t, y, dy, nid = _coerce_series(time, mag, mag_err, night_id)
    if len(t) < 10:
        return {
            "status": "INSUFFICIENT",
            "reason": "Not enough valid points (< 10).",
            "candidates": [],
            "window_peaks": [],
            "n_points": int(len(t)),
        }

    window = compute_spectral_window(t)
    candidates = build_alias_candidates(
        t,
        ls_frequency,
        ls_power,
        min_period,
        max_period,
        window=window,
        max_candidates=max_candidates,
    )
    if not candidates:
        return {
            "status": "INSUFFICIENT",
            "reason": "No usable period candidates.",
            "candidates": [],
            "window_peaks": window.get("peaks", []),
            "n_points": int(len(t)),
        }

    ls_frequency_arr = np.asarray(ls_frequency, dtype=float)
    finite_frequency = np.sort(ls_frequency_arr[np.isfinite(ls_frequency_arr)])
    grid_step = (
        float(np.nanmedian(np.diff(finite_frequency)))
        if len(finite_frequency) > 1
        else float(window.get("resolution_cd", 1.0)) / 10.0
    )
    refine_half_width = max(1.5 * abs(grid_step), 1e-5)
    search_fmin = 1.0 / float(max_period)
    search_fmax = 1.0 / float(min_period)
    offset_policy = str(night_offset_policy).strip().lower()
    if offset_policy not in {"none", "diagnostic", "fit"}:
        raise ValueError("night_offset_policy must be 'none', 'diagnostic', or 'fit'.")
    include_night_offsets = offset_policy == "fit" and len(np.unique(nid)) > 1

    fits: list[dict] = []
    for candidate in candidates:
        try:
            refined_frequency, fit = _refine_frequency_fit(
                t,
                y,
                dy,
                nid,
                initial_frequency=float(candidate["freq_cd"]),
                half_width=refine_half_width,
                min_frequency=search_fmin,
                max_frequency=search_fmax,
                harmonics=harmonics,
                include_night_offsets=include_night_offsets,
            )
        except (ValueError, np.linalg.LinAlgError):
            continue
        refined_period = 1.0 / refined_frequency
        coverage = phase_coverage_metrics(t, refined_period)
        row = dict(candidate)
        row.update(
            {
                "freq_cd": float(refined_frequency),
                "period": float(refined_period),
                "bic": float(fit["bic"]),
                "rmse": float(fit["rmse"]),
                "occupied_fraction": float(coverage["occupied_fraction"]),
                "max_phase_gap": float(coverage["max_phase_gap"]),
                "cycles": float(coverage["cycles"]),
                "leave_one_out_votes": 0,
                "_fit": fit,
            }
        )
        fits.append(row)

    if not fits:
        return {
            "status": "INSUFFICIENT",
            "reason": "All candidate fits were underconstrained or ill-conditioned.",
            "candidates": [],
            "window_peaks": list(window.get("peaks", [])),
            "n_points": int(len(t)),
            "night_offset_policy": offset_policy,
        }

    fits.sort(key=lambda row: float(row["bic"]))
    best_bic = float(fits[0]["bic"])
    best_period = float(fits[0]["period"])
    for rank, row in enumerate(fits, start=1):
        row["rank"] = rank
        row["delta_bic"] = float(row["bic"]) - best_bic
        if rank == 1:
            row["relation_to_best"] = "adopted"
        elif periods_are_window_aliases(
            best_period,
            float(row["period"]),
            window.get("peaks"),
            window.get("baseline_days", 1.0),
        ):
            row["relation_to_best"] = "window-alias"
        else:
            relation, _ = classify_frequency_relation(
                float(row["freq_cd"]),
                [float(fits[0]["freq_cd"])],
                window.get("peaks"),
                window.get("baseline_days", 1.0),
            )
            row["relation_to_best"] = relation

    offset_best_period = None
    offset_best_frequency = None
    offset_sensitive = False
    if offset_policy == "diagnostic" and len(np.unique(nid)) > 1:
        offset_trials: list[tuple[float, dict]] = []
        for row in fits:
            try:
                trial = _fit_frequency_model(
                    t,
                    y,
                    dy,
                    nid,
                    float(row["freq_cd"]),
                    harmonics,
                    include_night_offsets=True,
                )
            except (ValueError, np.linalg.LinAlgError):
                continue
            row["bic_with_night_offsets"] = float(trial["bic"])
            offset_trials.append((float(trial["bic"]), row))
        if offset_trials:
            _, offset_best = min(offset_trials, key=lambda item: item[0])
            offset_best_period = float(offset_best["period"])
            offset_best_frequency = float(offset_best["freq_cd"])
            offset_sensitive = int(offset_best["rank"]) != 1

    unique_nights = [night for night in np.unique(nid) if np.sum(nid == night) >= 10]
    loo_trials = 0
    if len(unique_nights) >= 2:
        for omitted in unique_nights:
            keep = nid != omitted
            if np.sum(keep) < 20:
                continue
            trial_scores = []
            for row in fits:
                try:
                    trial = _fit_frequency_model(
                        t[keep],
                        y[keep],
                        (dy[keep] if dy is not None else None),
                        nid[keep],
                        float(row["freq_cd"]),
                        harmonics,
                        include_night_offsets=include_night_offsets,
                    )
                    trial_scores.append(float(trial["bic"]))
                except (ValueError, np.linalg.LinAlgError):
                    trial_scores.append(np.inf)
            if not np.any(np.isfinite(trial_scores)):
                continue
            winner = int(np.argmin(trial_scores))
            fits[winner]["leave_one_out_votes"] += 1
            loo_trials += 1

    injection_recovery = np.nan
    if int(n_injections) > 0 and len(fits) > 1:
        rng = np.random.default_rng(int(random_seed))
        best_fit = fits[0]["_fit"]
        centered_residual = np.asarray(best_fit["residual"], dtype=float)
        centered_residual -= float(np.mean(centered_residual))
        recovered = 0
        successful_injections = 0
        for _ in range(int(n_injections)):
            simulated = np.asarray(best_fit["model"], dtype=float) + rng.choice(
                centered_residual, size=len(centered_residual), replace=True
            )
            scores = []
            for row in fits:
                try:
                    trial = _fit_frequency_model(
                        t,
                        simulated,
                        dy,
                        nid,
                        float(row["freq_cd"]),
                        harmonics,
                        include_night_offsets=include_night_offsets,
                    )
                    scores.append(float(trial["bic"]))
                except (ValueError, np.linalg.LinAlgError):
                    scores.append(np.inf)
            if not np.any(np.isfinite(scores)):
                continue
            successful_injections += 1
            recovered += int(int(np.argmin(scores)) == 0)
        if successful_injections:
            injection_recovery = float(recovered / successful_injections)

    night_spans = [float(np.ptp(t[nid == night])) for night in unique_nights]
    longest_night_span = max(night_spans, default=float(np.ptp(t)))
    longest_night_cycles = longest_night_span / best_period
    runner_delta = float(fits[1]["delta_bic"]) if len(fits) > 1 else np.inf
    loo_fraction = (
        float(fits[0]["leave_one_out_votes"] / loo_trials) if loo_trials else 1.0
    )

    global_cycles = float(fits[0]["cycles"])
    occupied_fraction = float(fits[0]["occupied_fraction"])
    max_phase_gap = float(fits[0]["max_phase_gap"])

    reasons: list[str] = []
    if len(t) < 30 or global_cycles < 1.0:
        status = "INSUFFICIENT"
        reasons.append("The full dataset covers less than one cycle or has too few points.")
    elif occupied_fraction < 0.35 or max_phase_gap > 0.55:
        status = "INSUFFICIENT"
        reasons.append(
            f"Global phase coverage is sparse (occupied={occupied_fraction:.0%}, "
            f"max gap={max_phase_gap:.2f})."
        )
    else:
        status = "RESOLVED"
        if runner_delta < 6.0:
            status = "AMBIGUOUS"
            reasons.append(f"Runner-up is close (delta BIC={runner_delta:.2f}).")
        if len(unique_nights) >= 2 and loo_fraction < 0.60:
            status = "AMBIGUOUS"
            reasons.append(f"Leave-one-night-out agreement is {loo_fraction:.0%}.")
        if np.isfinite(injection_recovery) and injection_recovery < 0.60:
            status = "AMBIGUOUS"
            reasons.append(f"Timestamp injection recovery is {injection_recovery:.0%}.")
        if longest_night_cycles < 1.5 and len(unique_nights) <= 2:
            status = "AMBIGUOUS"
            reasons.append(
                f"Only {len(unique_nights) or 1} night(s) are available and the longest "
                f"covers {longest_night_cycles:.2f} cycles."
            )
        elif (
            longest_night_cycles < 1.5
            and (occupied_fraction < 0.70 or max_phase_gap > 0.25)
        ):
            status = "AMBIGUOUS"
            reasons.append(
                "No single night spans 1.5 cycles and the combined phase coverage "
                "still has substantial gaps."
            )
        if offset_sensitive:
            status = "AMBIGUOUS"
            reasons.append(
                "The preferred alias changes when free nightly offsets are allowed."
            )
    if status == "RESOLVED":
        reasons.append("Candidate separation and subset checks are consistent.")

    serial_candidates = []
    for row in fits:
        clean = {key: value for key, value in row.items() if key != "_fit"}
        serial_candidates.append(clean)

    ls_frequency_values = np.asarray(ls_frequency, dtype=float)
    ls_power_values = np.asarray(ls_power, dtype=float)
    naive_index = int(np.nanargmax(ls_power_values))
    naive_period = float(1.0 / ls_frequency_values[naive_index])
    was_aliased = periods_are_window_aliases(
        best_period,
        naive_period,
        window.get("peaks", []),
        float(window.get("baseline_days", np.ptp(t))),
    )

    return {
        "status": status,
        "reason": " ".join(reasons),
        "adopted_period": best_period,
        "adopted_freq_cd": float(fits[0]["freq_cd"]),
        "naive_period": naive_period,
        "was_aliased": bool(was_aliased),
        "candidates": serial_candidates,
        "window_peaks": list(window.get("peaks", [])),
        "baseline_days": float(window.get("baseline_days", np.ptp(t))),
        "resolution_cd": float(window.get("resolution_cd", np.nan)),
        "n_points": int(len(t)),
        "n_nights": int(len(unique_nights) or 1),
        "longest_night_span_days": float(longest_night_span),
        "longest_night_cycles": float(longest_night_cycles),
        "global_cycles": global_cycles,
        "occupied_fraction": occupied_fraction,
        "max_phase_gap": max_phase_gap,
        "leave_one_out_trials": int(loo_trials),
        "leave_one_out_fraction": float(loo_fraction),
        "injection_recovery": (
            float(injection_recovery) if np.isfinite(injection_recovery) else None
        ),
        "harmonics": int(harmonics),
        "night_offset_policy": offset_policy,
        "offset_sensitive": bool(offset_sensitive),
        "offset_best_period": offset_best_period,
        "offset_best_frequency_cd": offset_best_frequency,
    }


def _multimode_design(
    time_rel: np.ndarray,
    periods: Sequence[float],
    harmonics: int,
) -> tuple[np.ndarray, list[dict]]:
    columns = [np.ones(len(time_rel), dtype=float)]
    terms: list[dict] = []
    for mode_index, period in enumerate(periods):
        base_freq = 1.0 / float(period)
        for harmonic in range(1, max(int(harmonics), 1) + 1):
            omega = 2.0 * np.pi * base_freq * harmonic
            for kind, values in (
                ("cos", np.cos(omega * time_rel)),
                ("sin", np.sin(omega * time_rel)),
            ):
                columns.append(values)
                terms.append(
                    {
                        "mode_index": int(mode_index),
                        "harmonic": int(harmonic),
                        "kind": kind,
                        "omega": float(omega),
                        "coefficient_index": len(columns) - 1,
                    }
                )
    return np.column_stack(columns), terms


def evaluate_multimode_result(result: dict, times: np.ndarray) -> dict:
    arr_t = np.asarray(times, dtype=float)
    tau = arr_t - float(result["time_ref"])
    coeff = np.asarray(result["coeff"], dtype=float)
    periods = [float(p) for p in result["periods"]]
    components = np.zeros((len(periods), len(arr_t)), dtype=float)
    derivatives = np.zeros_like(components)
    total = np.full(len(arr_t), float(coeff[0]), dtype=float)
    for term in result["terms"]:
        idx = int(term.get("coefficient_index", 0))
        omega = float(term["omega"])
        if term["kind"] == "cos":
            basis = np.cos(omega * tau)
            derivative = -omega * np.sin(omega * tau)
        else:
            basis = np.sin(omega * tau)
            derivative = omega * np.cos(omega * tau)
        contribution = float(coeff[idx]) * basis
        mode_idx = int(term["mode_index"])
        components[mode_idx] += contribution
        derivatives[mode_idx] += float(coeff[idx]) * derivative
        total += contribution
    return {
        "times": arr_t,
        "intercept": float(coeff[0]),
        "total": total,
        "components": components,
        "component_derivatives": derivatives,
        "total_derivative": np.sum(derivatives, axis=0),
    }


def fit_multimode_model(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
    periods: Sequence[float],
    harmonics: int = 1,
    night_id: Optional[np.ndarray] = None,
    include_night_offsets: bool = False,
) -> dict:
    """Simultaneously fit all supplied modes and their harmonics."""
    t, y, dy, nid = _coerce_series(time, mag, mag_err, night_id)
    clean_periods = [float(p) for p in periods if np.isfinite(p) and p > 0]
    if len(t) < 10:
        raise ValueError("Not enough valid data points for multi-mode fit (< 10).")
    if not clean_periods:
        raise ValueError("At least one positive period is required.")

    time_ref = float(np.min(t))
    design, terms = _multimode_design(t - time_ref, clean_periods, harmonics)
    unique_nights = list(dict.fromkeys(nid.tolist()))
    night_offset_terms: list[dict] = []
    if include_night_offsets and len(unique_nights) > 1:
        offset_columns = []
        for night in unique_nights[1:]:
            coefficient_index = design.shape[1] + len(offset_columns)
            offset_columns.append((nid == night).astype(float))
            night_offset_terms.append(
                {"night_id": str(night), "coefficient_index": coefficient_index}
            )
        design = np.column_stack([design, *offset_columns])
    weights = _normalized_weights(dy, len(t))
    solved = _solve_weighted_design(design, y, weights)
    coeff = np.asarray(solved["coeff"], dtype=float)
    evaluated = evaluate_multimode_result(
        {
            "time_ref": time_ref,
            "coeff": coeff,
            "periods": clean_periods,
            "terms": terms,
        },
        t,
    )
    night_offsets = {str(unique_nights[0]): 0.0} if unique_nights else {}
    for item in night_offset_terms:
        night_offsets[str(item["night_id"])] = float(
            coeff[int(item["coefficient_index"])]
        )
    return {
        "time": t,
        "mag": y,
        "mag_err": dy,
        "periods": clean_periods,
        "harmonics": int(harmonics),
        "coeff": coeff,
        "terms": terms,
        "night_offset_terms": night_offset_terms,
        "night_offsets": night_offsets,
        "night_offset_reference": str(unique_nights[0]) if unique_nights else None,
        "include_night_offsets": bool(night_offset_terms),
        "time_ref": time_ref,
        "intercept": float(coeff[0]),
        "model": solved["model"],
        "components": evaluated["components"],
        "component_derivatives": evaluated["component_derivatives"],
        "residual": solved["residual"],
        "baseline": float(np.ptp(t)),
        "n_points": int(solved["n_points"]),
        "n_params": int(solved["n_params"]),
        "rank": int(solved["rank"]),
        "rmse": float(solved["rmse"]),
        "wrms": float(solved["wrms"]),
        "bic": float(solved["bic"]),
        "design_condition": float(solved["design_condition"]),
    }


def _residual_frequency_relation(
    frequency: float,
    primary_frequency: float,
    window_peaks: Sequence[dict | float],
    baseline_days: float,
    max_harmonic: int,
) -> tuple[str, str]:
    resolution = 1.5 / max(float(baseline_days), 1e-9)
    alias_tol = max(0.02, min(0.12, 0.10 / max(float(baseline_days), 1.0)))
    offsets = _window_offsets(window_peaks)
    for order in range(1, max(int(max_harmonic), 1) + 1):
        harmonic_frequency = order * float(primary_frequency)
        separation = abs(float(frequency) - harmonic_frequency)
        if separation <= resolution:
            label = "primary" if order == 1 else f"{order}f_primary"
            return "unresolved", f"within 1.5/T of {label}"
        for offset in offsets:
            if abs(separation - offset) <= alias_tol:
                label = "primary" if order == 1 else f"{order}f_primary"
                return "alias", f"sampling-window alias of {label}"

    for order in range(2, max(int(max_harmonic), 2) + 1):
        subharmonic = float(primary_frequency) / order
        if abs(float(frequency) - subharmonic) <= resolution:
            return "harmonic", f"within 1.5/T of f_primary/{order}"
    return "independent", "independent residual-frequency candidate"


def _mode_fundamental_amplitude(result: dict, mode_index: int) -> float:
    coeff = np.asarray(result.get("coeff", []), dtype=float)
    values: dict[str, float] = {}
    for fallback_idx, term in enumerate(result.get("terms", []), start=1):
        if int(term.get("mode_index", -1)) != int(mode_index):
            continue
        if int(term.get("harmonic", 0)) != 1:
            continue
        coeff_idx = int(term.get("coefficient_index", fallback_idx))
        if 0 <= coeff_idx < len(coeff):
            values[str(term.get("kind", ""))] = float(coeff[coeff_idx])
    return float(np.hypot(values.get("cos", 0.0), values.get("sin", 0.0)))


def _fit_primary_with_symmetric_sidebands(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
    night_id: np.ndarray,
    primary_frequency: float,
    candidate_frequency: float,
    harmonics: int,
    include_night_offsets: bool,
) -> dict:
    """Fit a linearized amplitude/phase-modulated primary-mode hypothesis."""
    t = np.asarray(time, dtype=float)
    y = np.asarray(mag, dtype=float)
    nid = np.asarray(night_id, dtype=str)
    tau = t - float(np.min(t))
    delta = abs(float(candidate_frequency) - float(primary_frequency))
    mirror_frequency = float(primary_frequency) - np.sign(
        float(candidate_frequency) - float(primary_frequency)
    ) * delta
    if mirror_frequency <= 0:
        raise ValueError("The symmetric modulation sideband is outside positive frequencies.")

    columns = [np.ones(len(t), dtype=float)]
    for harmonic in range(1, max(int(harmonics), 1) + 1):
        phase = 2.0 * np.pi * float(primary_frequency) * harmonic * tau
        columns.extend([np.cos(phase), np.sin(phase)])
    for frequency in (float(candidate_frequency), mirror_frequency):
        phase = 2.0 * np.pi * frequency * tau
        columns.extend([np.cos(phase), np.sin(phase)])
    if include_night_offsets:
        unique_nights = list(dict.fromkeys(nid.tolist()))
        columns.extend((nid == night).astype(float) for night in unique_nights[1:])
    design = np.column_stack(columns)
    solved = _solve_weighted_design(design, y, _normalized_weights(mag_err, len(t)))
    solved["mirror_frequency_cd"] = mirror_frequency
    return solved


def _diagnose_multimode_once(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
    night_id: Optional[np.ndarray],
    alias_analysis: dict,
    min_period: float,
    max_period: float,
    harmonics: int = 2,
    samples_per_peak: int = 20,
    max_candidates: int = 8,
    include_night_offsets: bool = False,
) -> dict:
    """Test whether a single-mode model leaves a credible independent signal.

    This is deliberately a suspicion diagnostic rather than a mode counter.
    It follows the usual prewhitening workflow but requires the residual peak
    to survive sampling-window/harmonic checks, joint-model BIC comparison,
    and leave-one-night-out fits.
    """
    t, y, dy, nid = _coerce_series(time, mag, mag_err, night_id)
    primary_frequency = float(alias_analysis.get("adopted_freq_cd", np.nan))
    primary_status = str(alias_analysis.get("status", "INSUFFICIENT")).upper()
    if len(t) < 30 or not np.isfinite(primary_frequency) or primary_frequency <= 0:
        return {
            "status": "INCONCLUSIVE",
            "reason": "A usable primary frequency and at least 30 points are required.",
            "candidates": [],
        }

    include_offsets = bool(include_night_offsets) and len(np.unique(nid)) > 1
    primary_period = 1.0 / primary_frequency
    try:
        single_fit = fit_multimode_model(
            t,
            y,
            dy,
            periods=[primary_period],
            harmonics=harmonics,
            night_id=nid,
            include_night_offsets=include_offsets,
        )
    except (ValueError, np.linalg.LinAlgError) as exc:
        return {
            "status": "INCONCLUSIVE",
            "reason": f"Single-mode fit failed: {exc}",
            "candidates": [],
        }

    residual = np.asarray(single_fit["residual"], dtype=float)
    fmin = 1.0 / float(max_period)
    fmax = 1.0 / float(min_period)
    if not (np.isfinite(fmin) and np.isfinite(fmax) and 0 < fmin < fmax):
        return {
            "status": "INCONCLUSIVE",
            "reason": "The residual frequency-search range is invalid.",
            "candidates": [],
        }

    residual_ls = LombScargle(t, residual, dy) if dy is not None else LombScargle(t, residual)
    frequency, power = residual_ls.autopower(
        minimum_frequency=fmin,
        maximum_frequency=fmax,
        samples_per_peak=max(int(samples_per_peak), 5),
    )
    peak_indices = _periodogram_peaks(frequency, power, max_peaks=max_candidates)
    if not peak_indices:
        return {
            "status": "SINGLE-COMPATIBLE" if primary_status != "INSUFFICIENT" else "INCONCLUSIVE",
            "reason": "No residual periodogram peak was found in the search range.",
            "primary_period": primary_period,
            "primary_frequency_cd": primary_frequency,
            "single_bic": float(single_fit["bic"]),
            "candidates": [],
        }

    baseline = max(float(np.ptp(t)), 1e-9)
    window_peaks = list(alias_analysis.get("window_peaks", []))
    grid_step = float(np.nanmedian(np.diff(frequency))) if len(frequency) > 1 else 1.0 / baseline
    refine_half_width = max(1.5 * abs(grid_step), 1e-5)
    rows: list[dict] = []

    for peak_index in peak_indices:
        initial_frequency = float(frequency[int(peak_index)])

        def joint_bic(candidate_frequency: float) -> float:
            if not np.isfinite(candidate_frequency) or candidate_frequency <= 0:
                return np.inf
            try:
                fit = fit_multimode_model(
                    t,
                    y,
                    dy,
                    periods=[primary_period, 1.0 / float(candidate_frequency)],
                    harmonics=harmonics,
                    night_id=nid,
                    include_night_offsets=include_offsets,
                )
            except (ValueError, np.linalg.LinAlgError):
                return np.inf
            return float(fit["bic"])

        lower = max(fmin, initial_frequency - refine_half_width)
        upper = min(fmax, initial_frequency + refine_half_width)
        refined_frequency = initial_frequency
        if upper > lower:
            optimum = minimize_scalar(
                joint_bic,
                bounds=(lower, upper),
                method="bounded",
                options={"xatol": max((upper - lower) * 1e-5, 1e-8)},
            )
            if optimum.success and np.isfinite(optimum.fun):
                refined_frequency = float(optimum.x)

        try:
            joint_fit = fit_multimode_model(
                t,
                y,
                dy,
                periods=[primary_period, 1.0 / refined_frequency],
                harmonics=harmonics,
                night_id=nid,
                include_night_offsets=include_offsets,
            )
        except (ValueError, np.linalg.LinAlgError):
            continue

        residual_power = float(residual_ls.power(refined_frequency))
        try:
            residual_fap = float(
                residual_ls.false_alarm_probability(
                    residual_power,
                    method="baluev",
                    samples_per_peak=max(int(samples_per_peak), 5),
                    minimum_frequency=fmin,
                    maximum_frequency=fmax,
                )
            )
        except Exception:
            residual_fap = np.nan
        relation, note = _residual_frequency_relation(
            refined_frequency,
            primary_frequency,
            window_peaks,
            baseline,
            max_harmonic=max(4, int(harmonics) + 1),
        )
        modulation_bic = None
        modulation_mirror_frequency = None
        joint_advantage_over_modulation = None
        frequency_separation = abs(refined_frequency - primary_frequency)
        modulation_sideband_regime = (
            frequency_separation / max(primary_frequency, 1e-12) <= 0.20
        )
        if relation == "independent" and modulation_sideband_regime:
            try:
                modulation_fit = _fit_primary_with_symmetric_sidebands(
                    t,
                    y,
                    dy,
                    nid,
                    primary_frequency,
                    refined_frequency,
                    harmonics,
                    include_offsets,
                )
                modulation_bic = float(modulation_fit["bic"])
                modulation_mirror_frequency = float(
                    modulation_fit["mirror_frequency_cd"]
                )
                joint_advantage_over_modulation = float(
                    modulation_bic - float(joint_fit["bic"])
                )
            except (ValueError, np.linalg.LinAlgError):
                pass
        rows.append(
            {
                "frequency_cd": refined_frequency,
                "period": 1.0 / refined_frequency,
                "residual_power": residual_power,
                "residual_fap": residual_fap if np.isfinite(residual_fap) else None,
                "delta_bic": float(single_fit["bic"] - joint_fit["bic"]),
                "amplitude_mag": _mode_fundamental_amplitude(joint_fit, 1),
                "joint_bic": float(joint_fit["bic"]),
                "modulation_bic": modulation_bic,
                "modulation_mirror_frequency_cd": modulation_mirror_frequency,
                "joint_advantage_over_modulation_bic": joint_advantage_over_modulation,
                "modulation_sideband_regime": bool(modulation_sideband_regime),
                "relation": relation,
                "note": note,
            }
        )

    rows.sort(key=lambda row: float(row["delta_bic"]), reverse=True)
    independent = [row for row in rows if row["relation"] == "independent"]
    best = independent[0] if independent else None

    loo_trials = 0
    loo_votes = 0
    if best is not None:
        candidate_period = float(best["period"])
        usable_nights = [night for night in np.unique(nid) if np.sum(nid == night) >= 10]
        if len(usable_nights) >= 2:
            for omitted in usable_nights:
                keep = nid != omitted
                if np.sum(keep) < 30:
                    continue
                try:
                    single_subset = fit_multimode_model(
                        t[keep],
                        y[keep],
                        dy[keep] if dy is not None else None,
                        periods=[primary_period],
                        harmonics=harmonics,
                        night_id=nid[keep],
                        include_night_offsets=include_offsets,
                    )
                    joint_subset = fit_multimode_model(
                        t[keep],
                        y[keep],
                        dy[keep] if dy is not None else None,
                        periods=[primary_period, candidate_period],
                        harmonics=harmonics,
                        night_id=nid[keep],
                        include_night_offsets=include_offsets,
                    )
                except (ValueError, np.linalg.LinAlgError):
                    continue
                loo_trials += 1
                loo_votes += int(float(single_subset["bic"] - joint_subset["bic"]) >= 6.0)

    loo_fraction = float(loo_votes / loo_trials) if loo_trials else None
    primary_insufficient = primary_status == "INSUFFICIENT"
    if best is None:
        significant_nonindependent = any(
            row["residual_fap"] is not None
            and float(row["residual_fap"]) <= 0.05
            and float(row["delta_bic"]) >= 6.0
            for row in rows
        )
        if primary_insufficient or significant_nonindependent:
            status = "INCONCLUSIVE"
            reason = (
                "Residual structure is present, but it is unresolved from the primary, "
                "a harmonic, or a sampling-window alias."
            )
        else:
            status = "SINGLE-COMPATIBLE"
            reason = "No significant independent residual frequency was detected."
    else:
        fap = best.get("residual_fap")
        fap_value = float(fap) if fap is not None else np.nan
        delta_bic = float(best["delta_bic"])
        stable = loo_fraction is None or loo_fraction >= 0.60
        modulation_advantage = best.get("joint_advantage_over_modulation_bic")
        modulation_rejected = (
            modulation_advantage is None or float(modulation_advantage) >= 6.0
        )
        strong = (
            np.isfinite(fap_value)
            and fap_value <= 0.01
            and delta_bic >= 10.0
            and stable
            and modulation_rejected
        )
        weak = (np.isfinite(fap_value) and fap_value <= 0.05 and delta_bic >= 6.0) or delta_bic >= 10.0
        if strong and not primary_insufficient:
            status = "MULTIMODE-SUSPECT"
            reason = (
                f"Independent residual candidate at {best['frequency_cd']:.5f} d^-1 "
                f"improves the joint fit by delta BIC={delta_bic:.2f}."
            )
            if loo_fraction is not None:
                reason += f" Leave-one-night-out support is {loo_fraction:.0%}."
        elif weak:
            status = "INCONCLUSIVE"
            if not modulation_rejected:
                reason = (
                    "The residual sideband is also consistent with amplitude/phase "
                    "modulation of a single primary mode."
                )
            else:
                reason = (
                    "An independent residual candidate is present, but its significance, "
                    "sampling stability, or primary-period coverage is insufficient."
                )
        else:
            status = "SINGLE-COMPATIBLE"
            reason = "Independent residual candidates do not pass the significance and joint-fit thresholds."

    return {
        "status": status,
        "reason": reason,
        "primary_period": primary_period,
        "primary_frequency_cd": primary_frequency,
        "primary_alias_status": primary_status,
        "single_bic": float(single_fit["bic"]),
        "single_rmse": float(single_fit["rmse"]),
        "include_night_offsets": bool(include_offsets),
        "candidate_period": float(best["period"]) if best is not None else None,
        "candidate_frequency_cd": float(best["frequency_cd"]) if best is not None else None,
        "candidate_delta_bic": float(best["delta_bic"]) if best is not None else None,
        "candidate_fap": best.get("residual_fap") if best is not None else None,
        "leave_one_out_trials": int(loo_trials),
        "leave_one_out_fraction": loo_fraction,
        "resolution_cd": float(1.5 / baseline),
        "criteria": {
            "strong_fap_max": 0.01,
            "strong_delta_bic_min": 10.0,
            "subset_support_min": 0.60,
            "close_frequency_resolution": "1.5/T",
            "modulation_sideband_relative_max": 0.20,
            "joint_vs_modulation_delta_bic_min": 6.0,
        },
        "candidates": rows,
    }


def diagnose_multimode_suspicion(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
    night_id: Optional[np.ndarray],
    alias_analysis: dict,
    min_period: float,
    max_period: float,
    harmonics: int = 2,
    samples_per_peak: int = 20,
    max_candidates: int = 8,
    check_night_offset_sensitivity: bool = True,
) -> dict:
    """Diagnose extra modes and expose sensitivity to free nightly offsets."""
    primary = _diagnose_multimode_once(
        time,
        mag,
        mag_err,
        night_id,
        alias_analysis,
        min_period,
        max_period,
        harmonics=harmonics,
        samples_per_peak=samples_per_peak,
        max_candidates=max_candidates,
        include_night_offsets=False,
    )
    primary["offset_sensitivity"] = {"checked": False}

    nid = None if night_id is None else np.asarray(night_id)
    if (
        not check_night_offset_sensitivity
        or nid is None
        or len(nid) == 0
        or len(np.unique(nid.astype(str))) < 2
    ):
        return primary

    offset_result = _diagnose_multimode_once(
        time,
        mag,
        mag_err,
        night_id,
        alias_analysis,
        min_period,
        max_period,
        harmonics=harmonics,
        samples_per_peak=samples_per_peak,
        max_candidates=max_candidates,
        include_night_offsets=True,
    )
    status_without = str(primary.get("status", "INCONCLUSIVE"))
    status_with = str(offset_result.get("status", "INCONCLUSIVE"))
    sensitive = status_without != status_with
    primary["offset_sensitivity"] = {
        "checked": True,
        "sensitive": bool(sensitive),
        "status_without_offsets": status_without,
        "status_with_offsets": status_with,
        "candidate_frequency_without_offsets_cd": primary.get("candidate_frequency_cd"),
        "candidate_frequency_with_offsets_cd": offset_result.get("candidate_frequency_cd"),
        "reason_with_offsets": offset_result.get("reason", ""),
    }
    if sensitive:
        original_reason = str(primary.get("reason", "")).strip()
        primary["status"] = "INCONCLUSIVE"
        primary["reason"] = (
            "The mode classification changes when free nightly offsets are allowed. "
            "Treat the nightly zero point as unresolved; "
            + original_reason
        ).strip()
    return primary


def search_multimode_alias_solutions(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
    seed_periods: Sequence[float],
    harmonics: int = 1,
    window_peaks: Optional[Sequence[dict | float]] = None,
    night_id: Optional[np.ndarray] = None,
    include_night_offsets: bool = False,
    max_alias_offsets: int = 2,
    max_solutions: int = 1024,
) -> dict:
    """Compare simultaneous fits for sampling-alias combinations of modes."""
    t = _finite_time(time)
    if len(t) < 10:
        raise ValueError("Not enough data for multi-mode alias search.")
    periods = [float(p) for p in seed_periods if np.isfinite(p) and p > 0]
    if not periods:
        raise ValueError("At least one seed period is required.")

    if window_peaks is None:
        window_peaks = compute_spectral_window(t).get("peaks", [])
    offsets = _window_offsets(window_peaks)[: max(int(max_alias_offsets), 0)]
    alternatives: list[list[float]] = []
    for period in periods:
        base_freq = 1.0 / period
        mode_freqs = [base_freq]
        for offset in offsets:
            for sign in (-1.0, 1.0):
                candidate = base_freq + sign * offset
                if candidate > 0 and all(abs(candidate - old) > 1e-5 for old in mode_freqs):
                    mode_freqs.append(candidate)
        alternatives.append(mode_freqs)

    summaries: list[dict] = []
    best_fit: Optional[dict] = None
    best_bic = np.inf
    index_combinations = product(*(range(len(values)) for values in alternatives))
    ordered_indices = sorted(
        index_combinations,
        key=lambda combo: (
            sum((idx + 1) // 2 for idx in combo),
            sum(idx != 0 for idx in combo),
            combo,
        ),
    )
    total_combinations = int(len(ordered_indices))
    solution_limit = max(int(max_solutions), 1)
    search_complete = total_combinations <= solution_limit
    combinations: Iterable[tuple[float, ...]] = (
        tuple(alternatives[mode_idx][choice] for mode_idx, choice in enumerate(combo))
        for combo in ordered_indices
    )
    for combo in islice(combinations, solution_limit):
        if len(combo) != len({round(float(freq), 6) for freq in combo}):
            continue
        trial_periods = [1.0 / float(freq) for freq in combo]
        try:
            fit = fit_multimode_model(
                time,
                mag,
                mag_err,
                periods=trial_periods,
                harmonics=harmonics,
                night_id=night_id,
                include_night_offsets=bool(include_night_offsets),
            )
        except (ValueError, np.linalg.LinAlgError):
            continue
        summary = {
            "periods": [float(v) for v in trial_periods],
            "frequencies_cd": [float(v) for v in combo],
            "bic": float(fit["bic"]),
            "rmse": float(fit["rmse"]),
        }
        summaries.append(summary)
        if float(fit["bic"]) < best_bic:
            best_bic = float(fit["bic"])
            best_fit = fit

    if best_fit is None or not summaries:
        raise ValueError("No valid multi-mode alias solution could be fitted.")
    summaries.sort(key=lambda row: float(row["bic"]))
    for rank, row in enumerate(summaries, start=1):
        row["rank"] = rank
        row["delta_bic"] = float(row["bic"] - best_bic)
    runner_delta = float(summaries[1]["delta_bic"]) if len(summaries) > 1 else np.inf
    status = "AMBIGUOUS" if runner_delta < 6.0 else "RESOLVED"
    status_reasons = []
    if runner_delta < 6.0:
        status_reasons.append(f"Runner-up delta BIC={runner_delta:.2f}.")
    if not search_complete:
        status = "AMBIGUOUS"
        status_reasons.append(
            f"Alias search tested only {solution_limit} of {total_combinations} combinations."
        )

    t_valid, _, _, nid = _coerce_series(time, mag, mag_err, night_id)
    night_spans = [float(np.ptp(t_valid[nid == night])) for night in np.unique(nid)]
    longest_night_span = max(night_spans, default=float(np.ptp(t_valid)))
    n_nights = int(len(np.unique(nid)))
    baseline = float(np.ptp(t_valid))
    best_periods = [float(p) for p in best_fit["periods"]]
    phase_coverages = [phase_coverage_metrics(t_valid, period) for period in best_periods]
    sparse_modes = [
        idx
        for idx, coverage in enumerate(phase_coverages)
        if float(coverage["occupied_fraction"]) < 0.35
        or float(coverage["max_phase_gap"]) > 0.55
    ]
    if best_periods and baseline < max(best_periods):
        status = "INSUFFICIENT"
        status_reasons.append("The full dataset does not cover the longest fitted mode once.")
    elif sparse_modes:
        status = "INSUFFICIENT"
        labels = ", ".join(f"M{idx + 1}" for idx in sparse_modes)
        status_reasons.append(f"Combined phase coverage is sparse for {labels}.")
    elif (
        best_periods
        and longest_night_span < 1.5 * max(best_periods)
        and n_nights <= 2
    ):
        status = "AMBIGUOUS"
        status_reasons.append(
            "Only one or two nights are available and no night covers 1.5 cycles "
            "of the longest fitted mode."
        )

    best_freqs = [1.0 / period for period in best_periods]
    beat_periods = []
    for idx in range(len(best_freqs)):
        for jdx in range(idx + 1, len(best_freqs)):
            delta = abs(best_freqs[idx] - best_freqs[jdx])
            if delta > 0:
                beat_periods.append(1.0 / delta)
    strongest_window = max(
        (float(item.get("power", 0.0)) for item in (window_peaks or []) if isinstance(item, dict)),
        default=0.0,
    )
    beat_not_constrained = bool(
        beat_periods
        and strongest_window >= 0.5
        and (
            baseline < max(beat_periods)
            or (n_nights <= 2 and longest_night_span < max(beat_periods))
        )
    )
    if beat_not_constrained:
        status = "AMBIGUOUS" if status != "INSUFFICIENT" else status
        status_reasons.append(
            "The available nights do not constrain the beat cycle while window aliases are strong."
        )
    best_fit["alias_solutions"] = summaries[:10]
    best_fit["alias_status"] = status
    best_fit["alias_status_reason"] = " ".join(status_reasons)
    best_fit["alias_runner_delta_bic"] = runner_delta
    best_fit["longest_night_span_days"] = float(longest_night_span)
    best_fit["n_nights"] = n_nights
    best_fit["phase_coverages"] = phase_coverages
    best_fit["beat_periods_days"] = [float(value) for value in beat_periods]
    best_fit["window_peaks"] = [
        dict(item) if isinstance(item, dict) else {"freq_cd": float(item)}
        for item in (window_peaks or [])
    ]
    best_fit["alias_search_complete"] = bool(search_complete)
    best_fit["alias_total_combinations"] = total_combinations
    best_fit["alias_considered_combinations"] = int(
        min(solution_limit, total_combinations)
    )
    best_fit["alias_fitted_combinations"] = int(len(summaries))
    return best_fit
