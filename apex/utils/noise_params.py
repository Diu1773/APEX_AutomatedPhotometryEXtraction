"""Resolve effective detector gain/read-noise for a FITS frame.

The pipeline stores measured gain/read-noise in parameters.toml. FITS files can
carry per-frame metadata such as EGAIN, RDNOISE, and X/YBINNING, but those
header values are not always calibrated to the detector's measured PTC values.
This module keeps the interpretation in one place so Step 7/8 do not each
invent slightly different rules.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping


@dataclass(frozen=True)
class EffectiveNoiseParams:
    """Effective noise parameters for the stored image pixels."""

    gain_e_per_adu: float
    rdnoise_e: float
    bin_x: int
    bin_y: int
    bin_factor: int
    gain_source: str
    rdnoise_source: str


_GAIN_KEYS = ("EGAIN", "EPADU", "GAIN_EADU", "GAIN_E_PER_ADU")
_RDNOISE_KEYS = ("RDNOISE", "READNOIS", "READNOISE", "ENOISE")
_XBIN_KEYS = ("XBINNING", "XBIN", "BINX", "CCDXBIN", "HBIN")
_YBIN_KEYS = ("YBINNING", "YBIN", "BINY", "CCDYBIN", "VBIN")


def _header_get(header: Mapping[str, Any] | None, key: str) -> Any:
    if header is None:
        return None
    try:
        return header.get(key)
    except Exception:
        return None


def _positive_float(value: Any) -> float | None:
    try:
        out = float(str(value).strip())
    except Exception:
        return None
    return out if math.isfinite(out) and out > 0 else None


def _positive_int(value: Any) -> int | None:
    try:
        out = int(round(float(str(value).strip())))
    except Exception:
        return None
    return out if out > 0 else None


def _first_header_float(header: Mapping[str, Any] | None, keys: tuple[str, ...]) -> tuple[float | None, str | None]:
    for key in keys:
        val = _positive_float(_header_get(header, key))
        if val is not None:
            return val, key
    return None, None


def _parse_binning_pair(value: Any) -> tuple[int, int] | None:
    if value is None:
        return None
    text = str(value).strip()
    if not text:
        return None
    nums = re.findall(r"\d+(?:\.\d+)?", text)
    if len(nums) >= 2:
        bx = _positive_int(nums[0])
        by = _positive_int(nums[1])
        if bx and by:
            return bx, by
    if len(nums) == 1:
        b = _positive_int(nums[0])
        if b:
            return b, b
    return None


def infer_binning(
    header: Mapping[str, Any] | None,
    default_binning: Any = 1,
) -> tuple[int, int]:
    """Return image binning as (xbin, ybin)."""

    bx = None
    for key in _XBIN_KEYS:
        bx = _positive_int(_header_get(header, key))
        if bx is not None:
            break

    by = None
    for key in _YBIN_KEYS:
        by = _positive_int(_header_get(header, key))
        if by is not None:
            break

    if bx is not None and by is not None:
        return bx, by

    for key in ("CCDSUM", "BINNING", "BIN"):
        parsed = _parse_binning_pair(_header_get(header, key))
        if parsed is not None:
            return parsed

    fallback = _positive_int(default_binning) or 1
    return fallback, fallback


def resolve_effective_noise_params(
    params: Any,
    header: Mapping[str, Any] | None = None,
) -> EffectiveNoiseParams:
    """Resolve gain/RDNOISE for the stored image pixels.

    Rules:
    - Manual measured values are authoritative by default when present.
    - If a manual gain/read-noise value is missing, FITS EGAIN/RDNOISE-like
      keywords are used as fallback.
    - FITS EGAIN/RDNOISE-like keywords override manual values only when
      ``noise_use_fits_header`` is explicitly enabled.
    - Manual values are scaled by binning only when
      ``noise_scale_by_binning`` is explicitly enabled and
      ``noise_reference_binning`` describes the binning of the measured values.
      Otherwise manual values are treated as already effective for the stored
      FITS pixels.
    """

    default_binning = getattr(params, "binning_default", getattr(params, "binning", 1))
    bin_x, bin_y = infer_binning(header, default_binning=default_binning)
    bin_factor = max(int(bin_x) * int(bin_y), 1)

    ref_binning = _positive_int(getattr(params, "noise_reference_binning", None))
    if ref_binning is None:
        ref_binning = _positive_int(default_binning) or max(int(round(math.sqrt(bin_factor))), 1)
    ref_factor = max(ref_binning * ref_binning, 1)

    use_header = bool(getattr(params, "noise_use_fits_header", False))
    scale_manual = bool(getattr(params, "noise_scale_by_binning", False))

    manual_gain = _positive_float(getattr(params, "gain_e_per_adu", None))
    manual_rdnoise = _positive_float(getattr(params, "rdnoise_e", None))

    gain_scale = (float(bin_factor) / float(ref_factor)) if scale_manual else 1.0
    rd_scale = math.sqrt(float(bin_factor) / float(ref_factor)) if scale_manual else 1.0

    header_gain, gain_key = _first_header_float(header, _GAIN_KEYS)
    if use_header and header_gain is not None:
        gain = header_gain
        gain_source = f"header:{gain_key}"
    elif manual_gain is not None:
        gain = manual_gain * gain_scale
        gain_source = "manual" if abs(gain_scale - 1.0) < 1e-12 else f"manual*bin({bin_factor}/{ref_factor})"
    elif header_gain is not None:
        gain = header_gain
        gain_source = f"header:{gain_key}"
    else:
        gain = 1.0
        gain_source = "fallback:1.0"

    header_rdnoise, rdnoise_key = _first_header_float(header, _RDNOISE_KEYS)
    if use_header and header_rdnoise is not None:
        rdnoise = header_rdnoise
        rdnoise_source = f"header:{rdnoise_key}"
    elif manual_rdnoise is not None:
        rdnoise = manual_rdnoise * rd_scale
        rdnoise_source = "manual" if abs(rd_scale - 1.0) < 1e-12 else f"manual*sqrtbin({bin_factor}/{ref_factor})"
    elif header_rdnoise is not None:
        rdnoise = header_rdnoise
        rdnoise_source = f"header:{rdnoise_key}"
    else:
        rdnoise = 0.0
        rdnoise_source = "fallback:0.0"

    return EffectiveNoiseParams(
        gain_e_per_adu=float(gain),
        rdnoise_e=float(rdnoise),
        bin_x=int(bin_x),
        bin_y=int(bin_y),
        bin_factor=int(bin_factor),
        gain_source=gain_source,
        rdnoise_source=rdnoise_source,
    )
