"""Photon-transfer curve (PTC) measurement of gain and read noise.

Janesick's photon-transfer method (Janesick 2007) measured from the detector's
own calibration frames, with no dependence on FITS ``EGAIN``/``RDNOISE``
keywords.  On the Moravian C3-61000 the header ``EGAIN`` is the sensor's nominal
register value rather than the realised conversion gain, so the measured value
is the only usable one (see ``apex.benchmark.detector_characterize``).

The measurement:

* **Read noise** comes from the difference of two bias frames.  The difference
  of two independent reads has twice the single-frame variance, so
  ``RN = sqrt(var(bias_a - bias_b) / 2)``.
* **Gain** comes from the slope of ``var((flat_a - flat_b) / 2)`` against
  signal.  Differencing a same-level flat pair cancels the fixed pattern
  (PRNU, vignetting) and leaves shot noise, whose variance is
  ``S / g + RN²`` in ADU.  Fitting that line gives ``g = 1 / slope``.

All statistics are sigma-clipped so cosmic rays, hot pixels and dust motes do
not enter the variance.

Binned frames: the stored pixel of an ``n x n`` average-binned frame collects
``n²`` photosites, so its effective gain is ``n²`` times the per-photosite
value.  Both are reported — photometry runs on the stored pixels, so
``gain_eff`` is the number that belongs in the error model.

Qt-free: the GUI tool in ``apex.gui.tools.detector_ptc`` and the paper
benchmark in ``apex.benchmark.detector_characterize`` both call in here.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Sequence

import numpy as np

__all__ = [
    "PTCResult", "PTCPoint", "FrameInfo",
    "characterize_detector", "measure_ptc", "fit_ptc_points",
    "scan_calibration_frames", "measure_levels", "build_flat_pairs",
    "pairs_from_frames", "pairs_by_level", "read_box",
    "robust_center", "robust_diff_variance",
]

_CLIP_K = 4.0
_CLIP_ITERS = 4
_MIN_CLIP_SAMPLES = 64
_DEFAULT_BOX_RADIUS = 300
# Minimum (max-min)/max signal span before the slope is called weakly
# constrained.  Below this the formal error is optimistic.
_MIN_LEVER_ARM = 0.35


# ---------------------------------------------------------------------------
# Robust statistics
# ---------------------------------------------------------------------------

def _sigma_clip(values: np.ndarray, k: float = _CLIP_K) -> np.ndarray:
    """Iteratively drop points more than ``k`` robust sigma from the median."""
    v = values[np.isfinite(values)]
    if v.size == 0:
        return v
    for _ in range(_CLIP_ITERS):
        med = np.median(v)
        scale = np.median(np.abs(v - med)) * 1.4826 or 1.0
        keep = np.abs(v - med) <= k * scale
        if keep.sum() < _MIN_CLIP_SAMPLES:
            break
        v = v[keep]
    return v


def robust_center(values: np.ndarray, k: float = _CLIP_K) -> float:
    """Sigma-clipped median level."""
    v = _sigma_clip(np.asarray(values, dtype=np.float64).ravel(), k)
    return float(np.median(v)) if v.size else float("nan")


def robust_diff_variance(a: np.ndarray, b: np.ndarray, k: float = _CLIP_K) -> float:
    """Sigma-clipped variance of ``a - b``.

    Resistant to cosmic rays, hot pixels and any source that appears in only
    one of the two frames.
    """
    d = (np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)).ravel()
    d = _sigma_clip(d, k)
    return float(np.var(d)) if d.size else float("nan")


# ---------------------------------------------------------------------------
# Results
# ---------------------------------------------------------------------------

@dataclass
class PTCPoint:
    """One flat pair reduced to a (signal, variance) point."""

    signal_adu: float
    variance_adu2: float          # var of (a-b)/2, i.e. per-frame shot variance
    filter_name: str = ""
    source: tuple[str, str] = ("", "")

    def pair_gain(self, read_noise_adu: float) -> float:
        """Gain implied by this point alone, given the read noise."""
        excess = self.variance_adu2 - read_noise_adu ** 2
        if not np.isfinite(excess) or excess <= 0:
            return float("nan")
        return float(self.signal_adu / excess)


@dataclass
class PTCResult:
    ok: bool
    gain_eff: float                       # e-/ADU for the stored (binned) pixel
    gain_eff_err: float                   # 1-sigma from the fit
    read_noise_eff: float                 # e- for the stored pixel
    read_noise_adu: float
    gain_pixel: float                     # e-/ADU per photosite
    read_noise_pixel: float
    binning: int
    binning_mode: str
    fit_slope: float
    fit_slope_err: float
    fit_intercept: float
    n_pairs: int
    signal_min: float
    signal_max: float
    r_squared: float
    max_residual_frac: float              # worst |resid| / variance — linearity
    lever_arm: float = float("nan")       # (max-min)/max signal span
    points: list[PTCPoint] = field(default_factory=list)
    header_egain: float | None = None
    header_ratio: float | None = None
    message: str = ""
    log: list[str] = field(default_factory=list)

    def summary(self) -> str:
        lines = [
            f"pairs           : {self.n_pairs}",
            f"signal range    : {self.signal_min:,.0f} – {self.signal_max:,.0f} ADU",
            f"gain (stored px): {self.gain_eff:.4f} +- {self.gain_eff_err:.4f} e-/ADU",
            f"read noise      : {self.read_noise_eff:.3f} e-  ({self.read_noise_adu:.3f} ADU)",
            f"per photosite   : gain {self.gain_pixel:.4f} e-/ADU, RN {self.read_noise_pixel:.3f} e-",
            f"fit R^2         : {self.r_squared:.5f}",
            f"max |residual|  : {100 * self.max_residual_frac:.2f} % of variance",
            f"lever arm       : {100 * self.lever_arm:.0f} % of the top signal",
        ]
        if self.header_ratio is not None:
            lines.append(
                f"header EGAIN    : {self.header_egain:.5g}  "
                f"(measured / header = {self.header_ratio:.2f}x)"
            )
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Measurement
# ---------------------------------------------------------------------------

def _crop(arr: np.ndarray, box: tuple[int, int, int, int] | None) -> np.ndarray:
    if box is None:
        return arr
    y0, y1, x0, x1 = box
    return arr[y0:y1, x0:x1]


def _default_box(shape: tuple[int, ...], radius: int = _DEFAULT_BOX_RADIUS):
    h, w = shape[:2]
    cy, cx = h // 2, w // 2
    r = min(radius, h // 2 - 1, w // 2 - 1)
    return (cy - r, cy + r, cx - r, cx + r)


def read_box(path, radius: int = _DEFAULT_BOX_RADIUS, hdu_index: int = 0) -> np.ndarray:
    """Read only a centred square from a FITS image.

    Uses ``HDU.section``, which pulls just the requested slice off disk instead
    of materialising the whole frame.  On a 4788x3194 flat a 600x600 box is
    about 1/40 of the pixels, so a sweep over hundreds of pairs stops being
    I/O bound.

    Camera frames are unsigned integers carrying BZERO/BSCALE, and astropy
    refuses to memory-map those through its scaling path.  Scaling is therefore
    disabled on read and applied here, which keeps the section (and its small
    read) available.
    """
    from astropy.io import fits  # local import keeps the module import light

    with fits.open(path, memmap=True, do_not_scale_image_data=True) as hdul:
        hdu = hdul[hdu_index]
        bzero = float(hdu.header.get("BZERO", 0.0) or 0.0)
        bscale = float(hdu.header.get("BSCALE", 1.0) or 1.0)
        h, w = hdu.shape[:2]
        cy, cx = h // 2, w // 2
        r = min(radius, h // 2 - 1, w // 2 - 1)
        raw = hdu.section[cy - r:cy + r, cx - r:cx + r]
        return np.asarray(raw, dtype=np.float64) * bscale + bzero


def measure_ptc(
    pairs: Sequence[tuple],
    *,
    bias_pairs: Sequence[tuple] | None = None,
    bias_level: float | None = None,
    binning: int = 1,
    binning_mode: str = "average",
    box: tuple[int, int, int, int] | None = None,
    box_radius: int = _DEFAULT_BOX_RADIUS,
    signal_limits: tuple[float, float] | None = None,
    fix_intercept: bool = False,
    header_egain: float | None = None,
    loader: Callable[[object], np.ndarray | None] | None = None,
    filters: Sequence[str] | None = None,
    log_fn: Callable[[str], None] | None = None,
) -> PTCResult:
    """Measure gain and read noise from same-level frame pairs.

    Parameters
    ----------
    pairs
        Same-level frame pairs.  Each element is ``(a, b)`` of ndarrays or
        anything ``loader`` accepts.  Flats at many levels give the lever arm
        the slope needs.
    bias_pairs
        Bias frame pairs.  Read noise is the mean over pairs; more pairs give a
        tighter value.  Without them the read noise falls back to the fit
        intercept, which is far noisier.
    bias_level
        Pedestal in ADU removed from the flat level to form the signal axis.
        Defaults to the measured bias-pair level, else 0.
    binning, binning_mode
        ``n`` in ``n x n`` binning and how the camera combines photosites.
    box, box_radius
        Region used for statistics.  Defaults to a centred square of
        ``2 * box_radius`` pixels, which avoids vignetted corners.
    signal_limits
        ``(low, high)`` ADU bounds; pairs outside are dropped.  Use the high
        bound to stay below the onset of saturation, where the variance rolls
        over and biases the slope.
    fix_intercept
        Pin the intercept to the bias-measured ``RN_adu²`` and fit the slope
        alone.  The free two-parameter fit needs points near zero signal to
        constrain the intercept; without them the intercept runs away (it can
        even go negative, which is unphysical) and drags the slope with it.
        Since the read noise is already measured independently from bias pairs,
        pinning it is the better-determined choice whenever the flats sit far
        from the origin.  Requires ``bias_pairs``.
    filters
        Optional per-pair filter names, recorded on each point.
    """

    log: list[str] = []

    def _log(message: str) -> None:
        log.append(message)
        if log_fn is not None:
            log_fn(message)

    def _read(x) -> np.ndarray | None:
        if isinstance(x, np.ndarray):
            return np.asarray(x, dtype=np.float64)
        if loader is None:
            raise ValueError("a loader is required for non-array inputs")
        arr = loader(x)
        return None if arr is None else np.asarray(arr, dtype=np.float64)

    def _label(x) -> str:
        return Path(str(x)).name if not isinstance(x, np.ndarray) else "<array>"

    # ── region ────────────────────────────────────────────────────────────
    if box is None:
        for pr in pairs:
            sample = _read(pr[0])
            if sample is not None and sample.ndim == 2:
                box = _default_box(sample.shape, box_radius)
                break

    # ── read noise from bias pairs ────────────────────────────────────────
    rn_adu = float("nan")
    if bias_pairs:
        rn_values: list[float] = []
        levels: list[float] = []
        for ba_path, bb_path in bias_pairs:
            ba, bb = _read(ba_path), _read(bb_path)
            if ba is None or bb is None or ba.shape != bb.shape:
                continue
            ba, bb = _crop(ba, box), _crop(bb, box)
            var = robust_diff_variance(ba, bb)
            if np.isfinite(var) and var > 0:
                rn_values.append(math.sqrt(var / 2.0))
                levels.append(0.5 * (robust_center(ba) + robust_center(bb)))
        if rn_values:
            rn_adu = float(np.median(rn_values))
            spread = float(np.std(rn_values)) if len(rn_values) > 1 else 0.0
            if bias_level is None and levels:
                bias_level = float(np.median(levels))
            _log(
                f"[PTC] bias: {len(rn_values)} pairs, RN = {rn_adu:.4f} +- {spread:.4f} ADU, "
                f"level = {bias_level:.1f} ADU"
            )
    if bias_level is None:
        bias_level = 0.0

    # ── PTC points ────────────────────────────────────────────────────────
    points: list[PTCPoint] = []
    for i, (a_path, b_path) in enumerate(pairs):
        a, b = _read(a_path), _read(b_path)
        if a is None or b is None or a.shape != b.shape:
            continue
        a, b = _crop(a, box), _crop(b, box)
        signal = 0.5 * (robust_center(a) + robust_center(b)) - bias_level
        # var of the difference is twice the per-frame variance
        variance = robust_diff_variance(a, b) / 2.0
        if not (np.isfinite(signal) and np.isfinite(variance) and variance > 0):
            continue
        if signal_limits is not None and not (signal_limits[0] <= signal <= signal_limits[1]):
            continue
        points.append(PTCPoint(
            signal_adu=signal,
            variance_adu2=variance,
            filter_name=(filters[i] if filters is not None and i < len(filters) else ""),
            source=(_label(a_path), _label(b_path)),
        ))

    n = len(points)
    if n < 3:
        return PTCResult(
            ok=False, gain_eff=float("nan"), gain_eff_err=float("nan"),
            read_noise_eff=float("nan"), read_noise_adu=rn_adu,
            gain_pixel=float("nan"), read_noise_pixel=float("nan"),
            binning=binning, binning_mode=binning_mode,
            fit_slope=float("nan"), fit_slope_err=float("nan"),
            fit_intercept=float("nan"), n_pairs=n,
            signal_min=float("nan"), signal_max=float("nan"),
            r_squared=float("nan"), max_residual_frac=float("nan"),
            points=points, header_egain=header_egain,
            message="need at least 3 usable pairs spanning a signal range",
            log=log,
        )

    signal = np.array([p.signal_adu for p in points])
    variance = np.array([p.variance_adu2 for p in points])

    result = fit_ptc_points(
        points, read_noise_adu=rn_adu, binning=binning, binning_mode=binning_mode,
        fix_intercept=fix_intercept, header_egain=header_egain,
    )
    result.log = log + result.log
    for line in result.log[len(log):]:
        if log_fn is not None:
            log_fn(line)
    return result


def fit_ptc_points(
    points: Sequence[PTCPoint],
    *,
    read_noise_adu: float = float("nan"),
    binning: int = 1,
    binning_mode: str = "average",
    fix_intercept: bool = False,
    header_egain: float | None = None,
) -> PTCResult:
    """Fit gain from already-measured ``(signal, variance)`` points.

    Split out from :func:`measure_ptc` so a caller can read each frame once,
    cache the points, and then explore fit variants — signal cuts, pinned vs
    free intercept, per-epoch subsets — without touching the disk again.
    """
    log: list[str] = []
    n = len(points)
    signal = np.array([p.signal_adu for p in points], dtype=np.float64)
    variance = np.array([p.variance_adu2 for p in points], dtype=np.float64)

    if n < 3:
        return PTCResult(
            ok=False, gain_eff=float("nan"), gain_eff_err=float("nan"),
            read_noise_eff=float("nan"), read_noise_adu=read_noise_adu,
            gain_pixel=float("nan"), read_noise_pixel=float("nan"),
            binning=binning, binning_mode=binning_mode,
            fit_slope=float("nan"), fit_slope_err=float("nan"),
            fit_intercept=float("nan"), n_pairs=n,
            signal_min=float("nan"), signal_max=float("nan"),
            r_squared=float("nan"), max_residual_frac=float("nan"),
            points=list(points), header_egain=header_egain,
            message="need at least 3 usable pairs spanning a signal range", log=log,
        )

    rn_adu = read_noise_adu
    if fix_intercept:
        if not np.isfinite(rn_adu):
            raise ValueError("fix_intercept needs a measured read noise")
        intercept = rn_adu ** 2
        # single-parameter least squares through the pinned intercept
        offset = variance - intercept
        slope = float(signal @ offset / (signal @ signal))
        residual = variance - (signal * slope + intercept)
        dof = max(n - 1, 1)
        resid_var = float(residual @ residual) / dof
        slope_err = float(math.sqrt(resid_var / (signal @ signal)))
        log.append(f"[PTC] intercept pinned to RN^2 = {intercept:.3f} ADU^2")
    else:
        design = np.vstack([signal, np.ones_like(signal)]).T
        (slope, intercept), *_ = np.linalg.lstsq(design, variance, rcond=None)
        residual = variance - design @ np.array([slope, intercept])
        dof = max(n - 2, 1)
        resid_var = float(residual @ residual) / dof
        covariance = resid_var * np.linalg.inv(design.T @ design)
        slope_err = float(math.sqrt(max(covariance[0, 0], 0.0)))

    ss_tot = float(np.sum((variance - variance.mean()) ** 2))
    r_squared = 1.0 - float(residual @ residual) / ss_tot if ss_tot > 0 else float("nan")
    max_resid_frac = float(np.max(np.abs(residual) / np.maximum(variance, 1e-12)))

    gain_eff = 1.0 / slope if slope > 0 else float("nan")
    # g = 1/m  =>  sigma_g = sigma_m / m^2
    gain_eff_err = slope_err / (slope ** 2) if slope > 0 else float("nan")

    if not np.isfinite(rn_adu) and intercept > 0:
        rn_adu = math.sqrt(intercept)
        log.append("[PTC] read noise taken from the fit intercept (no bias pairs given)")
    read_noise_eff = gain_eff * rn_adu if np.isfinite(gain_eff) else float("nan")

    nb = max(int(binning), 1)
    gain_pixel = gain_eff / (nb * nb) if binning_mode == "average" else gain_eff
    read_noise_pixel = read_noise_eff / nb

    ratio = None
    if header_egain and header_egain > 0 and np.isfinite(gain_eff):
        ratio = float(gain_eff / header_egain)

    log.append(
        f"[PTC] {n} pairs, signal {signal.min():,.0f}-{signal.max():,.0f} ADU, "
        f"R^2 = {r_squared:.5f}"
    )
    log.append(
        f"[PTC] gain = {gain_eff:.4f} +- {gain_eff_err:.4f} e-/ADU, "
        f"RN = {read_noise_eff:.3f} e-"
    )

    # The quoted error is the precision of the fit, not its accuracy.  A short
    # signal range gives a short lever arm, so any systematic in the variance
    # estimate tilts the slope while the formal error stays small and
    # reassuring.  Warn rather than let that pass silently.
    span = float(signal.max() - signal.min())
    lever = span / float(signal.max()) if signal.max() > 0 else 0.0
    warning = ""
    if lever < _MIN_LEVER_ARM:
        warning = (
            f"signal range spans only {100 * lever:.0f}% of the top level "
            f"({span:,.0f} ADU) — the slope is weakly constrained and the quoted "
            f"error understates the true uncertainty; add flats at more exposure levels"
        )
        log.append(f"[PTC] WARNING: {warning}")

    return PTCResult(
        ok=bool(np.isfinite(gain_eff)),
        gain_eff=gain_eff, gain_eff_err=gain_eff_err,
        read_noise_eff=read_noise_eff, read_noise_adu=rn_adu,
        gain_pixel=gain_pixel, read_noise_pixel=read_noise_pixel,
        binning=nb, binning_mode=binning_mode,
        fit_slope=float(slope), fit_slope_err=slope_err,
        fit_intercept=float(intercept), n_pairs=n,
        signal_min=float(signal.min()), signal_max=float(signal.max()),
        r_squared=r_squared, max_residual_frac=max_resid_frac,
        lever_arm=lever,
        points=list(points), header_egain=header_egain, header_ratio=ratio,
        message=warning or "ok", log=log,
    )


# ---------------------------------------------------------------------------
# High-level entry point: calibration frames in, detector constants out
# ---------------------------------------------------------------------------

@dataclass
class FrameInfo:
    """A calibration frame with the header fields the pairing needs."""

    path: str
    image_type: str = ""
    filter_name: str = ""
    exptime: float = 0.0
    binning: int = 1
    header_egain: float | None = None
    level_adu: float = float("nan")     # filled in lazily by measure_levels


def scan_calibration_frames(paths: Iterable[str]) -> list[FrameInfo]:
    """Read headers only — cheap enough to run over a whole calibration folder."""
    from astropy.io import fits

    out: list[FrameInfo] = []
    for p in paths:
        try:
            header = fits.getheader(str(p))
        except Exception:  # noqa: BLE001 — an unreadable frame is simply skipped
            continue
        try:
            egain = float(header.get("EGAIN"))
        except (TypeError, ValueError):
            egain = None
        try:
            binning = int(header.get("XBINNING", 1) or 1)
        except (TypeError, ValueError):
            binning = 1
        out.append(FrameInfo(
            path=str(p),
            image_type=str(header.get("IMAGETYP", "")).strip(),
            filter_name=str(header.get("FILTER", "")).strip(),
            exptime=float(header.get("EXPTIME", header.get("EXPOSURE", 0)) or 0),
            binning=binning,
            header_egain=egain,
        ))
    return out


def measure_levels(
    frames: Sequence[FrameInfo],
    *,
    box_radius: int = _DEFAULT_BOX_RADIUS,
    progress: Callable[[int, int], None] | None = None,
) -> list[FrameInfo]:
    """Fill in each frame's central level.  Reads only the central box."""
    for i, frame in enumerate(frames, 1):
        try:
            frame.level_adu = robust_center(read_box(frame.path, box_radius))
        except Exception:  # noqa: BLE001
            frame.level_adu = float("nan")
        if progress is not None:
            progress(i, len(frames))
    return list(frames)


def build_flat_pairs(
    frames: Sequence[FrameInfo],
    *,
    bias_level: float = 0.0,
    tolerance: float = 0.02,
    signal_floor: float = 0.0,
    signal_ceiling: float | None = None,
) -> list[tuple[FrameInfo, FrameInfo]]:
    """Pair flats of the same filter whose levels agree to ``tolerance``.

    Only frames whose bias-subtracted signal sits inside
    ``[signal_floor, signal_ceiling]`` take part, so callers can stay clear of
    both the noisy low end and the onset of saturation.
    """
    by_filter: dict[str, list[FrameInfo]] = {}
    for frame in frames:
        signal = frame.level_adu - bias_level
        if not np.isfinite(signal) or signal < signal_floor:
            continue
        if signal_ceiling is not None and signal > signal_ceiling:
            continue
        by_filter.setdefault(frame.filter_name, []).append(frame)

    pairs: list[tuple[FrameInfo, FrameInfo]] = []
    for group in by_filter.values():
        group.sort(key=lambda f: f.level_adu)
        i = 0
        while i < len(group) - 1:
            a, b = group[i], group[i + 1]
            if abs(b.level_adu - a.level_adu) <= tolerance * max(abs(a.level_adu), 1.0):
                pairs.append((a, b))
                i += 2
            else:
                i += 1
    return pairs


def characterize_detector(
    bias_paths: Sequence[str],
    flat_paths: Sequence[str],
    *,
    binning: int | None = None,
    binning_mode: str = "average",
    box_radius: int = _DEFAULT_BOX_RADIUS,
    signal_floor: float = 20000.0,
    signal_ceiling: float | None = None,
    tolerance: float = 0.02,
    max_bias_pairs: int = 24,
    fix_intercept: bool = True,
    progress: Callable[[str, int, int], None] | None = None,
    log_fn: Callable[[str], None] | None = None,
) -> PTCResult:
    """Measure gain and read noise from a set of bias and flat frames.

    This is the entry point a GUI or CLI should call: hand it the calibration
    frames the user already has and it does the classification, pairing and
    fitting.

    ``fix_intercept`` defaults to True because flats useful for flat-fielding
    sit far from zero signal, where a free intercept is unconstrained.  The
    read noise measured from the bias pairs pins it instead.
    """
    def _progress(stage: str, i: int, n: int) -> None:
        if progress is not None:
            progress(stage, i, n)

    bias_frames = scan_calibration_frames(bias_paths)
    flat_frames = scan_calibration_frames(flat_paths)
    if len(bias_frames) < 2:
        raise ValueError("at least two bias frames are required")
    if len(flat_frames) < 6:
        raise ValueError("at least six flat frames are required")

    if binning is None:
        binning = bias_frames[0].binning or 1
    header_egain = next(
        (f.header_egain for f in flat_frames if f.header_egain), None)

    def _log(message: str) -> None:
        if log_fn is not None:
            log_fn(message)

    # ── bias: pedestal and read noise ─────────────────────────────────────
    bias_pairs = pairs_from_frames([f.path for f in bias_frames])[:max_bias_pairs]
    rn_values: list[float] = []
    levels: list[float] = []
    _progress("bias", 0, len(bias_pairs))
    for i, (pa, pb) in enumerate(bias_pairs, 1):
        a, b = read_box(pa, box_radius), read_box(pb, box_radius)
        var = robust_diff_variance(a, b)
        if np.isfinite(var) and var > 0:
            rn_values.append(math.sqrt(var / 2.0))
            levels.append(0.5 * (robust_center(a) + robust_center(b)))
        _progress("bias", i, len(bias_pairs))
    if not rn_values:
        raise ValueError("no usable bias pair — every difference had zero variance")
    rn_adu = float(np.median(rn_values))
    bias_level = float(np.median(levels))
    _log(f"[PTC] bias: {len(rn_values)} pairs, RN = {rn_adu:.4f} +- "
         f"{np.std(rn_values):.4f} ADU, level = {bias_level:.1f} ADU")

    # ── flats: levels, then same-level pairs ──────────────────────────────
    _progress("flat-levels", 0, len(flat_frames))
    measure_levels(flat_frames, box_radius=box_radius,
                   progress=lambda i, n: _progress("flat-levels", i, n))

    pairs = build_flat_pairs(
        flat_frames, bias_level=bias_level, tolerance=tolerance,
        signal_floor=signal_floor, signal_ceiling=signal_ceiling,
    )
    if len(pairs) < 3:
        raise ValueError(
            f"only {len(pairs)} flat pairs found above {signal_floor:,.0f} ADU; "
            "lower the signal floor or supply flats at more exposure levels"
        )
    _log(f"[PTC] {len(pairs)} flat pairs across "
         f"{len({a.filter_name for a, _ in pairs})} filters")

    # ── one point per pair ────────────────────────────────────────────────
    points: list[PTCPoint] = []
    _progress("ptc", 0, len(pairs))
    for i, (fa, fb) in enumerate(pairs, 1):
        a, b = read_box(fa.path, box_radius), read_box(fb.path, box_radius)
        signal = 0.5 * (robust_center(a) + robust_center(b)) - bias_level
        variance = robust_diff_variance(a, b) / 2.0
        if np.isfinite(signal) and np.isfinite(variance) and variance > 0:
            points.append(PTCPoint(
                signal_adu=signal, variance_adu2=variance,
                filter_name=fa.filter_name,
                source=(Path(fa.path).name, Path(fb.path).name),
            ))
        _progress("ptc", i, len(pairs))

    dropped = len(pairs) - len(points)
    if dropped:
        _log(f"[PTC] {dropped} pairs dropped (zero or non-finite difference variance)")

    result = fit_ptc_points(
        points, read_noise_adu=rn_adu, binning=binning,
        binning_mode=binning_mode, fix_intercept=fix_intercept,
        header_egain=header_egain,
    )
    for line in result.log:
        _log(line)
    return result


def pairs_from_frames(frames: Sequence) -> list[tuple]:
    """Pair consecutive frames as (0,1), (2,3), ….

    Consecutive exposures are closest in level and time, so their difference is
    the cleanest noise estimate.
    """
    return [(frames[i], frames[i + 1]) for i in range(0, len(frames) - 1, 2)]


def pairs_by_level(
    frames: Iterable[tuple[object, float]],
    *,
    tolerance: float = 0.02,
) -> list[tuple]:
    """Pair frames whose levels agree to within ``tolerance`` (fractional).

    ``frames`` is an iterable of ``(handle, level)``.  Frames are sorted by
    level and greedily paired, so a set of flats taken at many exposure times
    yields the widest possible signal range.
    """
    ordered = sorted(frames, key=lambda item: item[1])
    out: list[tuple] = []
    i = 0
    while i < len(ordered) - 1:
        (ha, la), (hb, lb) = ordered[i], ordered[i + 1]
        if abs(lb - la) <= tolerance * max(abs(la), 1.0):
            out.append((ha, hb))
            i += 2
        else:
            i += 1
    return out
