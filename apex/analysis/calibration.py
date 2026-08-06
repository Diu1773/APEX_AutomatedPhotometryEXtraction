"""Detector calibration (bias / dark / flat) — the science core of Step 0.

Qt-free, pure numpy/astropy (analysis layer → utils only).  Ports the
*science-legitimate* subset of the AstralImage/AIPPI individual-frame
preprocessing (``core/light_process.py::process_single_light``) plus a simple
master-combine (median / sigma-clipped mean).  Deliberately EXCLUDES the
imaging-only stages that are unsafe for photometry:

  * debayer (CFA demosaic) — interpolation destroys per-pixel photometric
    independence; APEX targets mono science CCDs, so it never applies anyway;
  * light-frame integration / registration / local normalisation ("LN") —
    APEX does forced photometry per frame and never coadds science frames;
  * linear-defect (row/column) correction — off by default (can flatten real
    sky gradients); not implemented here.

Calibration model (per science frame)::

    calibrated = ((raw - bias) - k * dark) / flat  + pedestal

where ``k`` is the dark scale (exposure ratio by default; an optional
noise-minimising fit refines it), pixels with ``flat < flat_min`` are marked
NaN (dead), and ``pedestal`` is an additive constant recorded in the header.
Because photometry subtracts a local sky annulus, a uniform pedestal cancels
and does not bias the measured magnitudes.
"""

from __future__ import annotations

from dataclasses import dataclass, fields as dataclass_fields
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np
from astropy.io import fits

from apex.utils import fast_stats
from apex.utils.constants import MAD_TO_SIGMA, EXPTIME_HEADER_KEYS
from apex.analysis.overscan import (
    correct_overscan,
    correct_overscan_from_header,
    update_header_for_overscan,
)

PathLike = Any

# Stack dtype for master combination — float32 halves peak RAM vs float64 and
# is ample precision for DN-scale calibration frames.
_STACK_DTYPE = np.float32
_OUTPUT_DTYPE = np.float32


@dataclass(frozen=True)
class CalibrationOptions:
    """Detector-calibration parameters (primitive values only, no Qt/config).

    The GUI / pipeline maps a ``CalibrationConfig`` onto this dataclass so the
    analysis layer stays free of config-model imports.
    """

    # --- master combine ---
    combine_method: str = "median"        # "median" | "mean" | "sigmaclip_mean"
    sigma_low: float = 3.0
    sigma_high: float = 3.0
    maxiters: int = 5

    # --- dark ---
    dark_scale: bool = True               # scale dark by exposure-time ratio
    dark_optimize: bool = False           # optional noise-min k-fit (off: imaging heuristic)
    # How far a dark's sensor temperature may sit from the light's before the
    # match is called into question. Dark current roughly doubles every ~6 C, so
    # 1 C is a few percent — tolerable for most work, but observers chasing
    # faint signal want it tighter, hence a configurable tolerance rather than
    # the old fixed 1 C rounding bucket. Beyond the tolerance the match is
    # reported (delta_T is logged and written to calibration.json) and, with
    # strict_temp on, refused outright.
    temp_match_tol_c: float = 1.0
    strict_temp: bool = False

    # --- flat ---
    flat_min: float = 0.01                # flat values below this -> dead pixel (NaN)

    # --- pedestal ---
    pedestal_mode: str = "adaptive"       # "none" | "adaptive" | "fixed"
    pedestal_value: float = 100.0         # DN, used by "fixed"
    pedestal_max: float = 5000.0          # DN cap for "adaptive"

    # --- overscan (off by default) ---
    overscan_enable: bool = False
    overscan_edge: str = "right"          # "left" | "right" | "top" | "bottom"
    overscan_width: int = 32
    overscan_trim: bool = True

    # --- cosmetic: cosmic-ray + hot pixel (L.A.Cosmic / astroscrappy) ---
    # On by default: standard practice in reduction pipelines (BANZAI removes
    # ~47% of the single-pixel spikes in its own e91 products; van Dokkum 2001).
    # Leaving it off let point-like noise dominate the ePSF reference stars on
    # CMOS detectors, where 16-40% of detections are 1-pixel spikes — the ePSF
    # then came out 2.75x too narrow and PSF fluxes fell to 32% of aperture
    # (measured 2026-07-29 on M67/QHY600; a CCD frame with 0% spikes was fine).
    # Measured cost on real stars: 94% are untouched in the core and shift only
    # because their sky annulus got cleaner (median +0.07%, more go up than
    # down); the 6% whose core is masked lose a median 0.44%.
    cosmetic_enable: bool = True
    cr_sigclip: float = 4.5               # L.A.Cosmic detection sigma
    cr_objlim: float = 5.0                # contrast limit (protects real sources)
    hot_sigma: float = 6.0                # hot-pixel threshold on the master dark
    gain: float = 1.0                     # e-/ADU fallback (header EGAIN preferred)
    readnoise: float = 6.5                # e- fallback (header RDNOISE preferred)

    # Skip frames whose header shows they were ALREADY bias/dark/flat-corrected
    # (CALBIAS/CALFLAT/CALSOFT keys or "CAL BIAS/FLAT" HISTORY, from a prior APEX
    # run or an external pipeline like MaxIm/AstralImage). Re-calibrating them
    # double-subtracts bias/dark and double-divides by flat, producing artifacts.
    skip_precalibrated: bool = True

    @classmethod
    def field_names(cls) -> Tuple[str, ...]:
        return tuple(f.name for f in dataclass_fields(cls))

    @classmethod
    def from_mapping(cls, mapping: Optional[Dict[str, Any]]) -> "CalibrationOptions":
        """Build options from a plain dict, ignoring unknown keys.

        The bridge between the TOML ``[calibration]`` table (or the GUI's
        settings dict) and the analysis layer, which must not import the config
        models.  Values that cannot be coerced fall back to the default rather
        than raising, so one bad line in a hand-edited TOML cannot block a run.
        """
        if not mapping:
            return cls()
        defaults = cls()
        kwargs: Dict[str, Any] = {}
        for field in dataclass_fields(cls):
            if field.name not in mapping:
                continue
            value = mapping[field.name]
            if value is None:
                continue
            current = getattr(defaults, field.name)
            try:
                if isinstance(current, bool):
                    if isinstance(value, str):
                        value = value.strip().lower() in ("1", "true", "yes", "on")
                    else:
                        value = bool(value)
                elif isinstance(current, int):
                    value = int(float(value))
                elif isinstance(current, float):
                    value = float(value)
                else:
                    value = str(value)
            except (TypeError, ValueError):
                continue
            kwargs[field.name] = value
        return cls(**kwargs)

    def to_mapping(self) -> Dict[str, Any]:
        return {f.name: getattr(self, f.name) for f in dataclass_fields(self)}


# ---------------------------------------------------------------------------
# FITS I/O + header helpers
# ---------------------------------------------------------------------------

def _load_fits(path: PathLike) -> Tuple[Optional[np.ndarray], Any]:
    """Return ``(data float64, header)`` from the first HDU with image data."""
    with fits.open(path, memmap=False) as hdul:
        for hdu in hdul:
            if getattr(hdu, "data", None) is not None:
                return np.asarray(hdu.data, dtype=np.float64), hdu.header.copy()
    return None, None


def _header_exptime(header, default: float = 1.0) -> float:
    if header is not None:
        for key in EXPTIME_HEADER_KEYS:
            if key in header:
                try:
                    return float(header[key])
                except (TypeError, ValueError):
                    continue
    return float(default)


def _header_float(header, keys, default: float) -> float:
    if header is not None:
        for key in keys:
            if key in header:
                try:
                    return float(header[key])
                except (TypeError, ValueError):
                    continue
    return float(default)


def _apply_overscan(data: np.ndarray, header, opts: CalibrationOptions) -> np.ndarray:
    """Apply overscan (header BIASSEC/DATASEC first, then manual edge/width)."""
    if not opts.overscan_enable:
        return data
    corrected, info = correct_overscan_from_header(data, header, trim=opts.overscan_trim)
    if not info.get("applied"):
        corrected, info = correct_overscan(
            data, edge=opts.overscan_edge, width=opts.overscan_width,
            trim=opts.overscan_trim,
        )
    if info.get("applied"):
        update_header_for_overscan(header, info)
        return corrected
    return data


def load_frame(path: PathLike, opts: CalibrationOptions) -> Tuple[np.ndarray, Any]:
    """Load a raw frame and apply overscan (shared by masters + science).

    Applying overscan at load keeps master frames and science frames on the
    same (optionally trimmed) pixel grid, so subtraction/division line up.
    """
    data, header = _load_fits(path)
    if data is None:
        raise ValueError(f"Unsupported or unreadable image: {path}")
    data = _apply_overscan(data, header, opts)
    return data, header


# ---------------------------------------------------------------------------
# Combine (stacking + simple sigma-clip rejection)
# ---------------------------------------------------------------------------

def _sigma_clip_mask(stack: np.ndarray, sigma_low: float, sigma_high: float,
                     maxiters: int) -> np.ndarray:
    """Per-pixel sigma-clip mask (median centre, robust MAD spread).

    Adapted from AstralImage ``core/rejection.py::reject_sigma_clip``.  Returns
    a boolean array (True = rejected) matching ``stack``.
    """
    mask = ~np.isfinite(stack)
    for _ in range(max(1, int(maxiters))):
        work = np.where(mask, np.nan, stack)
        med = fast_stats.nanmedian(work, axis=0)[np.newaxis]
        mad = fast_stats.nanmedian(np.abs(work - med), axis=0)[np.newaxis] * MAD_TO_SIGMA
        mad = np.where((~np.isfinite(mad)) | (mad == 0), 1e-10, mad)
        residual = stack - med
        finite = np.isfinite(residual)
        new_mask = mask | (
            finite & ((residual > sigma_high * mad) | (residual < -sigma_low * mad))
        )
        if np.array_equal(new_mask, mask):
            break
        mask = new_mask
    return mask


# Peak working set allowed for one combine chunk.  sigmaclip_mean holds the
# stack, a boolean mask and a clipped copy at once, so the true peak is a few
# times this; 384 MB keeps a 30-frame full-frame combine near 1 GB instead of
# the ~6 GB that stacking everything at once needed.
_COMBINE_CHUNK_BYTES = 384 * 1024 * 1024


def _combine_stack(stack: np.ndarray, method: str, sigma_low: float,
                   sigma_high: float, maxiters: int) -> np.ndarray:
    """Collapse an already-stacked (N, …) block along axis 0."""
    if method == "mean":
        return np.asarray(fast_stats.nanmean(stack, axis=0), dtype=np.float64)
    if method == "sigmaclip_mean":
        mask = _sigma_clip_mask(stack, sigma_low, sigma_high, maxiters)
        work = np.where(mask, np.nan, stack)
        return np.asarray(fast_stats.nanmean(work, axis=0), dtype=np.float64)
    return np.asarray(fast_stats.nanmedian(stack, axis=0), dtype=np.float64)


def combine_frames(arrays: Sequence[np.ndarray], method: str = "median",
                   sigma_low: float = 3.0, sigma_high: float = 3.0,
                   maxiters: int = 5) -> np.ndarray:
    """Combine a list of 2-D frames into one master (float64).

    ``method``: "median" (robust default), "mean", or "sigmaclip_mean"
    (per-pixel sigma clip, then NaN-mean of the survivors).  NaNs in the inputs
    (e.g. dead pixels) are treated as missing.

    Combining is done in row bands rather than on one big stack.  Every method
    here is per-pixel, so banding is exact — it only bounds the working set,
    which otherwise scales with frame count and overflowed on a full-frame
    30-flat combine.
    """
    if not arrays:
        raise ValueError("combine_frames: empty array list")
    method = str(method).lower()

    first = np.asarray(arrays[0], dtype=_STACK_DTYPE)
    if len(arrays) == 1:
        return first.astype(np.float64)
    if first.ndim != 2:
        stack = np.stack([np.asarray(a, dtype=_STACK_DTYPE) for a in arrays], axis=0)
        return _combine_stack(stack, method, sigma_low, sigma_high, maxiters)

    height, width = first.shape
    row_bytes = len(arrays) * width * np.dtype(_STACK_DTYPE).itemsize
    rows_per_chunk = max(1, int(_COMBINE_CHUNK_BYTES // max(row_bytes, 1)))

    out = np.empty((height, width), dtype=np.float64)
    for y0 in range(0, height, rows_per_chunk):
        y1 = min(y0 + rows_per_chunk, height)
        band = np.stack(
            [np.asarray(a[y0:y1], dtype=_STACK_DTYPE) for a in arrays], axis=0)
        out[y0:y1] = _combine_stack(band, method, sigma_low, sigma_high, maxiters)
    return out


# ---------------------------------------------------------------------------
# Pre-built master detection
# ---------------------------------------------------------------------------
# Mirrors AstralImage/AIPPI: a frame whose IMAGETYP reads "MASTER <kind>" is a
# pre-built master and is used verbatim instead of being re-stacked — so a user
# can drop a single ready-made master bias/dark/flat into the scan and have it
# applied directly (no dedicated UI field, exactly like AIPPI's auto-detection).

_IMAGETYP_KEYS = ("IMAGETYP", "IMGTYPE", "FRAMETYP", "OBSTYPE")


def _is_master_header(header, kind: str) -> bool:
    """True when the header marks this frame as a pre-built master of ``kind``."""
    if header is None:
        return False
    for key in _IMAGETYP_KEYS:
        if key in header:
            it = str(header[key]).upper()
            if "MASTER" in it and kind.upper() in it:
                return True
    return False


def _find_master(paths: Sequence[PathLike], kind: str, opts: CalibrationOptions):
    """If any path is a pre-built master of ``kind``, load it (overscan applied
    for grid parity) and return ``(data, header, path)``; else ``None``."""
    for p in paths:
        try:
            data, header = _load_fits(p)
        except Exception:
            continue
        if data is not None and _is_master_header(header, kind):
            return _apply_overscan(data, header, opts), header, str(p)
    return None


# ---------------------------------------------------------------------------
# Master builders
# ---------------------------------------------------------------------------

def build_master_bias(paths: Sequence[PathLike],
                      opts: CalibrationOptions) -> Tuple[np.ndarray, Dict]:
    """Combine bias frames into a master bias (or use a pre-built master)."""
    m = _find_master(paths, "BIAS", opts)
    if m is not None:
        data, _hdr, src = m
        return data, {"type": "bias", "n_frames": 1, "master_input": True,
                      "source": Path(src).name,
                      "median": float(fast_stats.finite_nanmedian(data, 0.0))}
    arrays = [load_frame(p, opts)[0] for p in paths]
    master = combine_frames(arrays, opts.combine_method,
                            opts.sigma_low, opts.sigma_high, opts.maxiters)
    prov = {
        "type": "bias",
        "n_frames": len(paths),
        "method": opts.combine_method,
        "median": float(fast_stats.finite_nanmedian(master, 0.0)),
    }
    return master, prov


def build_master_dark(paths: Sequence[PathLike], opts: CalibrationOptions,
                      master_bias: Optional[np.ndarray] = None
                      ) -> Tuple[np.ndarray, float, Dict]:
    """Combine dark frames (bias-subtracted) into a master dark.

    Returns ``(master_dark, dark_exptime, provenance)`` where ``dark_exptime``
    is the reference (median) exposure used to scale to science frames.
    """
    m = _find_master(paths, "DARK", opts)
    if m is not None:
        data, header, src = m           # a master dark is already bias-subtracted
        dark_exp = _header_exptime(header)
        return data, dark_exp, {
            "type": "dark", "n_frames": 1, "master_input": True,
            "source": Path(src).name, "exptime": dark_exp,
            "bias_subtracted": True,
            "median": float(fast_stats.finite_nanmedian(data, 0.0))}
    arrays: List[np.ndarray] = []
    exps: List[float] = []
    for p in paths:
        data, header = load_frame(p, opts)
        if master_bias is not None:
            data = data - master_bias
        arrays.append(data)
        exps.append(_header_exptime(header))
    master = combine_frames(arrays, opts.combine_method,
                            opts.sigma_low, opts.sigma_high, opts.maxiters)
    dark_exp = float(np.median(exps)) if exps else 1.0
    prov = {
        "type": "dark",
        "n_frames": len(paths),
        "method": opts.combine_method,
        "exptime": dark_exp,
        "bias_subtracted": master_bias is not None,
        "median": float(fast_stats.finite_nanmedian(master, 0.0)),
    }
    return master, dark_exp, prov


def build_master_flat(paths: Sequence[PathLike], opts: CalibrationOptions,
                      master_bias: Optional[np.ndarray] = None,
                      master_dark: Optional[np.ndarray] = None,
                      dark_exp: float = 1.0) -> Tuple[np.ndarray, Dict]:
    """Combine flat frames into a master flat normalised to unit median.

    Each flat is bias/dark-subtracted, then divided by its own median so
    illumination/exposure differences cancel before combining; the master is
    re-normalised to a median of 1.0.
    """
    m = _find_master(paths, "FLAT", opts)
    if m is not None:
        data, _hdr, src = m             # a master flat is already reduced
        mm = float(fast_stats.finite_nanmedian(data, 1.0))
        if mm and np.isfinite(mm) and abs(mm - 1.0) > 0.05:
            data = data / mm            # defensively renormalise to unit median
        return data, {"type": "flat", "n_frames": 1, "master_input": True,
                      "source": Path(src).name, "bias_subtracted": True,
                      "dark_subtracted": True,
                      "median": float(fast_stats.finite_nanmedian(data, 0.0))}
    arrays: List[np.ndarray] = []
    for p in paths:
        data, header = load_frame(p, opts)
        if master_bias is not None:
            data = data - master_bias
        if master_dark is not None:
            exp = _header_exptime(header)
            ratio = (exp / dark_exp) if (opts.dark_scale and dark_exp > 0) else 1.0
            data = data - master_dark * ratio
        med = float(fast_stats.finite_nanmedian(data, 0.0))
        if med and np.isfinite(med):
            data = data / med
        arrays.append(data)
    master = combine_frames(arrays, opts.combine_method,
                            opts.sigma_low, opts.sigma_high, opts.maxiters)
    mm = float(fast_stats.finite_nanmedian(master, 1.0))
    if mm and np.isfinite(mm):
        master = master / mm
    prov = {
        "type": "flat",
        "n_frames": len(paths),
        "method": opts.combine_method,
        "bias_subtracted": master_bias is not None,
        "dark_subtracted": master_dark is not None,
        "median": float(fast_stats.finite_nanmedian(master, 0.0)),
    }
    return master, prov


# ---------------------------------------------------------------------------
# Per-frame apply
# ---------------------------------------------------------------------------

def _dark_scale(light_bias_sub: np.ndarray, master_dark: np.ndarray,
                exposure_ratio: float, optimize: bool) -> float:
    """Dark scale ``k``.  Exposure ratio unless ``optimize`` refines it.

    The optional fit minimises calibrated background noise via
    ``k = Σ(L·D) / Σ(D²)`` on median-centred pixels, clamped around the
    exposure ratio (matched to AstralImage's bounded PI-style optimisation).
    """
    if not optimize:
        return float(exposure_ratio)
    finite = np.isfinite(light_bias_sub) & np.isfinite(master_dark)
    if not np.any(finite):
        return float(exposure_ratio)
    L = light_bias_sub[finite].astype(np.float64)
    D = master_dark[finite].astype(np.float64)
    Lc = L - np.median(L)
    Dc = D - np.median(D)
    denom = float(np.dot(Dc, Dc))
    if denom <= 0:
        return float(exposure_ratio)
    k = float(np.dot(Lc, Dc)) / denom
    k_low = max(0.1, exposure_ratio * 0.5)
    k_high = min(3.0, exposure_ratio * 1.5)
    if k_low > k_high:
        k_low, k_high = k_high, k_low
    k = float(np.clip(k, k_low, k_high))
    return k if np.isfinite(k) else float(exposure_ratio)


def _apply_pedestal(data: np.ndarray, mode: str, value: float,
                    pmax: float) -> Tuple[np.ndarray, float]:
    """Add an output pedestal (additive constant) after flat calibration."""
    mode = str(mode).lower()
    if mode == "none":
        return data, 0.0
    if mode == "fixed":
        ped = float(value)
        data += ped
        return data, ped
    # adaptive: lift only if the low tail went negative
    low = float(np.nanpercentile(data, 0.01)) if np.any(np.isfinite(data)) else 0.0
    ped = 0.0
    if low < 0:
        ped = float(min(np.ceil(abs(low)) + 100.0, pmax))
        data += ped
    return data, ped


def frame_stats(data: np.ndarray) -> Dict[str, float]:
    """Simple per-frame QC statistics (the 'simple statistics' of Step 0)."""
    arr = np.asarray(data, dtype=np.float64)
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return {"finite_pct": 0.0, "min": 0.0, "max": 0.0, "mean": 0.0,
                "median": 0.0, "mad_sigma": 0.0, "neg": 0, "neg_pct": 0.0}
    med = float(np.median(finite))
    mad = float(np.median(np.abs(finite - med)) * MAD_TO_SIGMA)
    neg = int(np.sum(finite < 0))
    return {
        "finite_pct": float(finite.size / max(1, arr.size) * 100.0),
        "min": float(np.min(finite)),
        "max": float(np.max(finite)),
        "mean": float(np.mean(finite)),
        "median": med,
        "mad_sigma": mad,
        "neg": neg,
        "neg_pct": float(neg / max(1, arr.size) * 100.0),
    }


def frame_is_calibrated(header) -> bool:
    """True if the header records prior bias/dark/flat calibration — from a
    previous APEX run or an external pipeline (MaxIm DL, AstralImage). Detected
    from the CALSOFT/CALBIAS/CALFLAT/CALDARK keywords (APEX writes booleans,
    external tools write the master-frame filename) or "CAL BIAS/FLAT" HISTORY
    lines. Used to skip re-calibration and avoid double bias/dark/flat."""
    if header is None:
        return False
    for key in ("CALSOFT", "CALBIAS", "CALFLAT", "CALDARK"):
        if key in header:
            val = header[key]
            if isinstance(val, bool):
                if val:
                    return True
            elif str(val).strip().lower() not in ("", "false", "0", "none"):
                return True
    try:
        hist = " ".join(str(h) for h in header.get("HISTORY", [])).upper()
        if "CAL BIAS" in hist or "CAL FLAT" in hist or "APEX CALIBRATION" in hist:
            return True
    except Exception:
        pass
    return False


def calibrate_light(data: np.ndarray, header, opts: CalibrationOptions,
                    master_bias: Optional[np.ndarray] = None,
                    master_dark: Optional[np.ndarray] = None,
                    dark_exp: float = 1.0,
                    master_flat: Optional[np.ndarray] = None
                    ) -> Tuple[np.ndarray, Any, Dict]:
    """Apply detector calibration to one (overscan-corrected) science frame.

    ``data`` is assumed already overscan-corrected (use :func:`load_frame`).
    Returns ``(calibrated float32, header, qc)``.  No debayer, no LDC.
    """
    work = np.array(data, dtype=np.float64, copy=True)
    prov: Dict[str, Any] = {}

    # Guard: never re-calibrate an already-calibrated frame. A prior APEX run or
    # an external pipeline (MaxIm/AstralImage) leaves CALBIAS/CALFLAT/CALSOFT in
    # the header; re-applying bias/dark/flat would double-subtract and double-
    # divide, leaving artifacts (e.g. i-band fringe residuals that SEP then
    # detects as thousands of spurious sources). Pass such frames through as-is.
    if opts.skip_precalibrated and frame_is_calibrated(header):
        qc = frame_stats(work)
        qc.update({"already_calibrated": True, "skipped_recalibration": True})
        if header is not None:
            header.add_history("APEX Step 0: skipped (already calibrated)")
        return work.astype(_OUTPUT_DTYPE), header, qc

    # A. bias
    if master_bias is not None:
        work -= master_bias
        prov["bias"] = True

    # B. dark (exposure-ratio scale by default; optional noise-min refine)
    if master_dark is not None:
        exp = _header_exptime(header)
        ratio = (exp / dark_exp) if (opts.dark_scale and dark_exp > 0) else 1.0
        k = _dark_scale(work, master_dark, ratio, opts.dark_optimize)
        work -= master_dark * k
        prov["dark_scale"] = float(k)
        prov["dark_exptime"] = float(dark_exp)

    # C. flat (dead pixels -> NaN)
    if master_flat is not None:
        bad = master_flat < opts.flat_min
        safe = np.where(bad, 1.0, master_flat)
        work /= safe
        work[bad] = np.nan
        prov["flat"] = True
        prov["flat_bad_pct"] = float(np.count_nonzero(bad) / max(1, bad.size) * 100.0)

    # C2. cosmetic — cosmic-ray (L.A.Cosmic/astroscrappy, star-protected) + hot
    # pixels (dark-derived mask). Optional; standard CCD reduction, off by default.
    if opts.cosmetic_enable:
        from apex.analysis.cosmetic import clean_frame, hot_pixel_mask, HAS_ASTROSCRAPPY
        if HAS_ASTROSCRAPPY:
            # Use the measured/configured gain, NOT the FITS EGAIN keyword: on
            # this camera (Moravian C3-61000) the header EGAIN is the nominal
            # max-gain value and is ~14x off the actual delivered gain, so it
            # must not drive the CR noise model. RDNOISE is usually absent, so it
            # falls back to the configured value.
            gain = float(opts.gain)
            rdn = _header_float(header, ("RDNOISE", "READNOIS", "RDNOIS"), opts.readnoise)
            sat = _header_float(header, ("SATURATE", "DATAMAX"), 65535.0)
            hmask = hot_pixel_mask(master_dark, opts.hot_sigma)
            work, _cmask, ncorr = clean_frame(
                work, gain=gain, readnoise=rdn, satlevel=sat,
                sigclip=opts.cr_sigclip, objlim=opts.cr_objlim, hot_mask=hmask)
            prov["cosmetic_pixels"] = ncorr

    # D. pedestal (additive constant, recorded — cancels under sky subtraction)
    work, pedestal = _apply_pedestal(work, opts.pedestal_mode,
                                     opts.pedestal_value, opts.pedestal_max)
    prov["pedestal"] = float(pedestal)

    # E. header provenance
    if header is not None:
        header["PEDESTAL"] = int(round(pedestal))
        header["CALBIAS"] = bool(master_bias is not None)
        header["CALDARK"] = bool(master_dark is not None)
        header["CALFLAT"] = bool(master_flat is not None)
        header["CALSOFT"] = ("APEX", "Detector calibration performed by APEX Step 0")
        header.add_history(
            "APEX calibration: ((raw-bias)-k*dark)/flat + pedestal"
        )

    qc = frame_stats(work)
    qc.update(prov)
    return work.astype(_OUTPUT_DTYPE), header, qc


def calibrate_light_file(path: PathLike, opts: CalibrationOptions,
                         master_bias: Optional[np.ndarray] = None,
                         master_dark: Optional[np.ndarray] = None,
                         dark_exp: float = 1.0,
                         master_flat: Optional[np.ndarray] = None
                         ) -> Tuple[np.ndarray, Any, Dict]:
    """Convenience: load (with overscan) then :func:`calibrate_light`."""
    data, header = load_frame(path, opts)
    return calibrate_light(data, header, opts, master_bias=master_bias,
                           master_dark=master_dark, dark_exp=dark_exp,
                           master_flat=master_flat)
