"""Reusable period-analysis service functions."""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Callable, List, Optional, Tuple
import warnings

import numpy as np
from astropy.timeseries import LombScargle
from scipy.signal import find_peaks

from apex.utils.constants import get_parallel_workers


def _pdm_theta_vectorized(
    dt: np.ndarray,
    y: np.ndarray,
    y2: np.ndarray,
    trial_periods: np.ndarray,
    n_bins: int,
    var_total: float,
    batch_size: int = 2000,
) -> np.ndarray:
    """Vectorised Stellingwerf (1978) PDM theta over all trial periods.

    Replaces the pure-Python per-period loop with batched numpy bincount
    operations.  Processes ``batch_size`` periods at once to bound memory use.
    Typical speedup: 20-100× over the scalar loop.
    """
    n_trials = len(trial_periods)
    n_stars = len(dt)
    theta = np.ones(n_trials, dtype=float)

    # Auto-shrink batch_size so each batch uses ≤ 50 MB.
    # Three arrays of shape (B, n_stars) are live per batch:
    # phases (float64=8), bins (int32≈8 after broadcast), global_idx (int64=8).
    bytes_per_period = max(n_stars * 24, 1)
    batch_size = max(100, min(batch_size, (50 * 1024 * 1024) // bytes_per_period))

    for b_start in range(0, n_trials, batch_size):
        b_end = min(b_start + batch_size, n_trials)
        B = b_end - b_start

        # phases: (B, n_stars) — fraction of period elapsed
        phases = (dt[None, :] / trial_periods[b_start:b_end, None]) % 1.0
        bins = np.clip((phases * n_bins).astype(np.int32), 0, n_bins - 1)

        # Encode as flat (period_index * n_bins + bin_index) for bincount
        global_idx = (np.arange(B, dtype=np.int64)[:, None] * n_bins + bins).ravel()

        counts_flat = np.bincount(global_idx, minlength=B * n_bins)
        sums_flat   = np.bincount(global_idx, weights=np.tile(y,  B), minlength=B * n_bins)
        sumsq_flat  = np.bincount(global_idx, weights=np.tile(y2, B), minlength=B * n_bins)

        counts = counts_flat.reshape(B, n_bins)
        sums   = sums_flat.reshape(B, n_bins)
        sumsq  = sumsq_flat.reshape(B, n_bins)

        good = counts >= 2
        c_safe = np.where(good, counts, 1)
        ss  = np.where(good, sumsq - sums ** 2 / c_safe, 0.0)
        dof = np.where(good, counts - 1, 0)

        ss_sum  = ss.sum(axis=1)
        dof_sum = dof.sum(axis=1)
        valid = dof_sum > 0
        theta[b_start:b_end] = np.where(valid, ss_sum / dof_sum / var_total, 1.0)

    return theta


def _find_top_periods(
    power: np.ndarray,
    periods: np.ndarray,
    best_idx: int,
    n_top: int = 5,
) -> tuple[list[float], list[float]]:
    """Return (top_periods, top_powers) from a power spectrum.

    Finds peaks above 10% of max power; falls back to best_idx when none found.
    """
    peak_indices, _ = find_peaks(power, height=0.1 * power[best_idx])
    if len(peak_indices) == 0:
        peak_indices = [best_idx]
    sorted_peaks = sorted(peak_indices, key=lambda i: power[i], reverse=True)[:n_top]
    return [float(periods[i]) for i in sorted_peaks], [float(power[i]) for i in sorted_peaks]


def filter_valid(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
) -> Tuple[np.ndarray, np.ndarray, Optional[np.ndarray]]:
    mask = np.isfinite(time) & np.isfinite(mag)
    if mag_err is not None:
        mask &= np.isfinite(mag_err)
    t = time[mask]
    y = mag[mask]
    dy = mag_err[mask] if mag_err is not None else None
    return t, y, dy


def compute_ls(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
    label: str,
    min_period: float,
    max_period: float,
    samples_per_peak: int,
) -> dict:
    t, y, dy = filter_valid(time, mag, mag_err)

    if len(t) < 10:
        return {
            "label": label, "method": "ls",
            "error": "Not enough data points (< 10)",
            "best_period": np.nan, "best_power": np.nan,
        }

    if dy is not None and np.any(dy > 0):
        ls = LombScargle(t, y, dy)
    else:
        ls = LombScargle(t, y)

    frequency, power = ls.autopower(
        minimum_frequency=1.0 / max_period,
        maximum_frequency=1.0 / min_period,
        samples_per_peak=samples_per_peak,
    )

    best_idx = np.argmax(power)
    best_freq = float(frequency[best_idx])
    best_period = 1.0 / best_freq
    best_power = float(power[best_idx])

    try:
        with np.errstate(invalid="ignore", divide="ignore", over="ignore"):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                fap = float(ls.false_alarm_probability(best_power))
    except Exception:
        fap = np.nan

    top_periods, top_powers = _find_top_periods(
        power, 1.0 / np.asarray(frequency, float), best_idx
    )

    return {
        "label": label, "method": "ls",
        "frequency": np.array(frequency, dtype=float),
        "power": np.array(power, dtype=float),
        "best_period": best_period, "best_power": best_power, "fap": fap,
        "top_periods": top_periods, "top_powers": top_powers,
        "n_points": len(t), "time": t, "mag": y, "mag_err": dy,
    }


def compute_pdm(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
    label: str,
    min_period: float,
    max_period: float,
    samples_per_peak: int,
    pdm_n_bins: int,
) -> dict:
    t, y, dy = filter_valid(time, mag, mag_err)

    if len(t) < 10:
        return {
            "label": label, "method": "pdm",
            "error": "Not enough data points (< 10)",
            "best_period": np.nan, "best_power": np.nan,
        }

    baseline = t.max() - t.min()
    n_trials_raw = int(samples_per_peak * baseline / min_period)
    n_trials = min(50000, n_trials_raw)
    if n_trials_raw > 50000:
        warnings.warn(
            f"PDM [{label}]: requested {n_trials_raw} trial periods capped to 50000. "
            f"Consider increasing min_period (current={min_period:.4f} d) or "
            f"reducing baseline ({baseline:.1f} d) to improve period resolution.",
            UserWarning,
            stacklevel=3,
        )
    # Linear frequency spacing gives uniform resolution across all periods.
    # Period-linear spacing over-samples long periods and under-samples short ones.
    trial_freqs = np.linspace(1.0 / max_period, 1.0 / min_period, n_trials)
    trial_periods = 1.0 / trial_freqs[::-1]

    # Stellingwerf (1978) defines σ² as the sample variance (ddof=1).
    var_total = float(np.var(y, ddof=1))
    if var_total == 0:
        return {
            "label": label, "method": "pdm",
            "error": "Zero variance in data",
            "best_period": np.nan, "best_power": np.nan,
        }

    n_bins = pdm_n_bins
    dt = t - t.min()
    y2 = y * y
    theta = _pdm_theta_vectorized(dt, y, y2, trial_periods, n_bins, var_total)

    power = 1.0 - theta

    best_idx = np.argmax(power)
    best_period = float(trial_periods[best_idx])
    best_power = float(power[best_idx])
    best_theta = float(theta[best_idx])

    top_periods, top_powers = _find_top_periods(power, trial_periods, best_idx)

    return {
        "label": label, "method": "pdm",
        "trial_periods": trial_periods,
        "theta": theta, "power": power,
        "best_period": best_period, "best_power": best_power,
        "best_theta": best_theta, "fap": np.nan,
        "top_periods": top_periods, "top_powers": top_powers,
        "n_points": len(t), "time": t, "mag": y, "mag_err": dy,
    }


def compute_bls(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
    label: str,
    min_period: float,
    max_period: float,
) -> dict:
    from astropy.timeseries import BoxLeastSquares

    t, y, dy = filter_valid(time, mag, mag_err)

    if len(t) < 10:
        return {
            "label": label, "method": "bls",
            "error": "Not enough data points (< 10)",
            "best_period": np.nan, "best_power": np.nan,
        }

    if dy is not None and np.any(dy > 0):
        bls = BoxLeastSquares(t, y, dy=dy)
    else:
        bls = BoxLeastSquares(t, y)

    baseline = t.max() - t.min()
    # max_dur: up to 25 % of the longest trial period or 25 % of baseline.
    # Removed the old min_period * 0.5 cap which prevented detecting long
    # transits (e.g. 1-2 day duration around a 5-10 day period).
    max_dur = min(max_period * 0.25, baseline * 0.25)
    min_dur = max(min_period * 0.01, max_dur * 0.05)
    if min_dur >= max_dur:
        min_dur = max_dur * 0.1
    durations = np.linspace(min_dur, max_dur, 10)

    try:
        result = bls.autopower(
            durations,
            minimum_period=min_period,
            maximum_period=max_period,
        )
    except Exception as e:
        return {
            "label": label, "method": "bls",
            "error": f"BLS failed: {e}",
            "best_period": np.nan, "best_power": np.nan,
        }

    power = result.power
    periods = result.period

    best_idx = np.argmax(power)
    best_period = float(periods[best_idx])
    best_power = float(power[best_idx])

    top_periods, top_powers = _find_top_periods(
        np.asarray(power, float), np.asarray(periods, float), best_idx
    )

    return {
        "label": label, "method": "bls",
        "trial_periods": np.array(periods, dtype=float),
        "power": np.array(power, dtype=float),
        "best_period": best_period, "best_power": best_power,
        "fap": np.nan,
        "top_periods": top_periods, "top_powers": top_powers,
        "n_points": len(t), "time": t, "mag": y, "mag_err": dy,
    }


def bootstrap_fap(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
    best_power: float,
    min_period: float,
    max_period: float,
    samples_per_peak: int = 10,
    n_bootstrap: int = 1000,
    seed: Optional[int] = None,
    progress_cb: Optional[Callable[[int, int], None]] = None,
) -> float:
    """Bootstrap false alarm probability for Lomb-Scargle.

    Permutes the magnitude values *n_bootstrap* times (keeping the time
    sampling fixed) and computes the LS peak power for each permutation.
    The FAP is the fraction of permutations whose peak power ≥ *best_power*.

    This is more reliable than the analytical FAP (which assumes Gaussian
    noise and specific sampling patterns) especially for short, uneven,
    or correlated time series.

    Parameters
    ----------
    time, mag, mag_err : observed time series (NaN filtered internally).
    best_power : the LS peak power to test (from ``compute_ls``).
    min_period, max_period : period search range (days).
    samples_per_peak : LS frequency grid density.
    n_bootstrap : number of Monte-Carlo permutations (default 1000).
    seed : random seed for reproducibility.
    progress_cb : optional callback(current, total) for UI progress.

    Returns
    -------
    float : bootstrap FAP in [0, 1].  Returns NaN if fewer than 10 points.
    """
    t, y, dy = filter_valid(time, mag, mag_err)
    if len(t) < 10:
        return np.nan

    rng = np.random.default_rng(seed)
    freq_kw = dict(
        minimum_frequency=1.0 / max_period,
        maximum_frequency=1.0 / min_period,
        samples_per_peak=samples_per_peak,
    )

    # Pre-compute the frequency grid once from the *observed* data so all
    # permutations use the identical grid — avoids recomputing it 1000 times.
    _ls_ref = LombScargle(t, y, dy) if dy is not None else LombScargle(t, y)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        _frequency_grid, _ = _ls_ref.autopower(**freq_kw)

    # Pre-generate permutations sequentially so a fixed seed gives reproducible
    # results regardless of worker count or completion order.
    perms = [rng.permutation(y) for _ in range(n_bootstrap)]

    def _peak_power(y_perm: np.ndarray) -> float:
        ls = LombScargle(t, y_perm, dy) if dy is not None else LombScargle(t, y_perm)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # ls.power() uses the pre-computed grid — no grid recomputation per perm
            power_perm = ls.power(_frequency_grid)
        return float(np.max(power_perm))

    # LombScargle._unscaled_periodogram is a C extension that releases the GIL,
    # so threads scale near-linearly with cores for this workload.
    max_workers = max(1, min(get_parallel_workers(), n_bootstrap))
    n_exceed = 0
    n_done = 0
    if max_workers == 1:
        for y_perm in perms:
            if _peak_power(y_perm) >= best_power:
                n_exceed += 1
            n_done += 1
            if progress_cb is not None:
                progress_cb(n_done, n_bootstrap)
    else:
        with ThreadPoolExecutor(max_workers=max_workers) as ex:
            futures = [ex.submit(_peak_power, p) for p in perms]
            for fut in as_completed(futures):
                if fut.result() >= best_power:
                    n_exceed += 1
                n_done += 1
                if progress_cb is not None:
                    progress_cb(n_done, n_bootstrap)

    return float(n_exceed / n_bootstrap)


def run_period_analysis(
    time: np.ndarray,
    mag_raw: np.ndarray,
    mag_corr: Optional[np.ndarray],
    mag_err: Optional[np.ndarray],
    min_period: float,
    max_period: float,
    samples_per_peak: int = 10,
    methods: Optional[List[str]] = None,
    pdm_n_bins: int = 10,
    progress_cb: Optional[Callable[[str], None]] = None,
) -> dict:
    results = {}
    for method in (methods or ["ls"]):
        if progress_cb is not None:
            progress_cb(f"Computing {method.upper()} for raw magnitudes...")
        results[f"raw_{method}"] = _compute_method(
            time, mag_raw, mag_err, "raw", method,
            min_period, max_period, samples_per_peak, pdm_n_bins,
        )
        if mag_corr is not None and np.any(np.isfinite(mag_corr)):
            if progress_cb is not None:
                progress_cb(f"Computing {method.upper()} for corrected magnitudes...")
            results[f"corr_{method}"] = _compute_method(
                time, mag_corr, mag_err, "corr", method,
                min_period, max_period, samples_per_peak, pdm_n_bins,
            )
    return results


def _compute_method(
    time: np.ndarray,
    mag: np.ndarray,
    mag_err: Optional[np.ndarray],
    label: str,
    method: str,
    min_period: float,
    max_period: float,
    samples_per_peak: int,
    pdm_n_bins: int,
) -> dict:
    if method == "ls":
        return compute_ls(time, mag, mag_err, label, min_period, max_period, samples_per_peak)
    if method == "pdm":
        return compute_pdm(time, mag, mag_err, label, min_period, max_period, samples_per_peak, pdm_n_bins)
    if method == "bls":
        return compute_bls(time, mag, mag_err, label, min_period, max_period)
    return {"label": label, "method": method, "error": f"Unknown method: {method}"}
