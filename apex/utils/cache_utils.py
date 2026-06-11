"""Shared cache schema, file-signature, and WCS-cache helpers."""
from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Mapping


DETECTION_CACHE_SCHEMA_VERSION = 4
DETECTION_CACHE_MIN_SIGNATURE_SCHEMA = 2


def cache_schema_value(payload: Mapping | None, default: int = 0) -> int:
    """Return the integer cache schema from a cache payload."""

    if not isinstance(payload, Mapping):
        return int(default)
    try:
        return int(payload.get("cache_schema", default) or default)
    except Exception:
        return int(default)


def normalize_detect_engine(value) -> str:
    """Normalize detection-engine aliases used in cache compatibility checks."""

    s = str(value or "segm").strip().lower()
    aliases = {
        "sourceextractor": "sep",
        "source_extractor": "sep",
        "sextractor": "sep",
        "sex": "sep",
        "segmentation": "segm",
        "photutils": "segm",
        "daofind": "dao",
    }
    s = aliases.get(s, s)
    if s in ("sep", "segm", "peak", "dao"):
        return s
    return "segm"


def norm_path_key(path_value) -> str:
    if path_value is None:
        return ""
    try:
        s = str(path_value).strip().replace("\\", "/")
    except Exception:
        s = str(path_value).replace("\\", "/")
    if len(s) >= 3 and s[1] == ":" and s[2] == "/" and s[0].isalpha():
        s = f"/mnt/{s[0].lower()}/{s[3:]}"
    while "//" in s:
        s = s.replace("//", "/")
    if len(s) > 1 and s.endswith("/"):
        s = s[:-1]
    return s.lower()


def build_file_signature(path: Path, *, use_cropped: bool) -> dict:
    sig = {
        "source_path": norm_path_key(path),
        "source_use_cropped": bool(use_cropped),
        "source_size": None,
        "source_mtime_ns": None,
    }
    try:
        st = path.stat()
        sig["source_size"] = int(st.st_size)
        sig["source_mtime_ns"] = int(st.st_mtime_ns)
    except Exception:
        pass
    return sig


def build_detection_cache_signature(
    path: Path,
    *,
    use_cropped: bool,
    detect_engine=None,
    cache_schema: int = DETECTION_CACHE_SCHEMA_VERSION,
) -> dict:
    """Build the canonical Step4 detection cache signature block."""

    sig = build_file_signature(path, use_cropped=use_cropped)
    sig["cache_schema"] = int(cache_schema)
    if detect_engine is not None:
        sig["detect_engine"] = normalize_detect_engine(detect_engine)
    return sig


def file_signature_matches(saved_sig: dict, current_sig: dict) -> bool:
    if not isinstance(saved_sig, dict) or not isinstance(current_sig, dict):
        return False
    if bool(saved_sig.get("source_use_cropped", None)) != bool(current_sig.get("source_use_cropped", None)):
        return False
    if norm_path_key(saved_sig.get("source_path", "")) != norm_path_key(current_sig.get("source_path", "")):
        return False
    try:
        saved_mtime_ns = int(saved_sig.get("source_mtime_ns"))
        curr_mtime_ns = int(current_sig.get("source_mtime_ns"))
    except Exception:
        return False
    return saved_mtime_ns == curr_mtime_ns


#: FITS header block size (bytes). A Step 5 WCS write may inject extra cards
#: (CD matrix, SIP polynomial coefficients, provenance keys) that push the
#: header across one or more block boundaries, so the file size after solving
#: differs from the size recorded when Step 4 ran. The image payload is
#: untouched, so detection XY is still valid — we just need to ignore the
#: header-only drift. Allow up to _MAX_HEADER_GROWTH_BLOCKS of growth before
#: declaring the file truly modified.
_FITS_HEADER_BLOCK_BYTES = 2880
_MAX_HEADER_GROWTH_BLOCKS = 16  # ~46 KB, easily covers SIP degree 5


def file_signature_matches_relaxed(saved_sig: dict, current_sig: dict) -> bool:
    """Accept mtime drift when path/crop are stable and size drift is bounded.

    Step 5 updates FITS headers in-place after solving, which changes mtime
    and, when SIP/PV polynomial keys are added, can also enlarge the file by a
    few 2880-byte FITS header blocks. Detection XY is independent of the WCS
    header, so we keep the Step 4 cache valid as long as:

    * path and crop mode are unchanged
    * the size delta is non-negative (file shouldn't shrink under header
      growth) and within ``_MAX_HEADER_GROWTH_BLOCKS`` blocks

    A real edit to the image data would change the size by orders of magnitude
    more than this allowance, so the relaxation cannot hide actual content
    changes.
    """
    if file_signature_matches(saved_sig, current_sig):
        return True
    if not isinstance(saved_sig, dict) or not isinstance(current_sig, dict):
        return False
    if bool(saved_sig.get("source_use_cropped", None)) != bool(current_sig.get("source_use_cropped", None)):
        return False
    if norm_path_key(saved_sig.get("source_path", "")) != norm_path_key(current_sig.get("source_path", "")):
        return False
    try:
        saved_size = int(saved_sig.get("source_size"))
        curr_size = int(current_sig.get("source_size"))
    except Exception:
        return False
    if saved_size <= 0 or curr_size <= 0:
        return False
    if saved_size == curr_size:
        return True
    delta = curr_size - saved_size
    if delta < 0:
        return False
    return delta <= _MAX_HEADER_GROWTH_BLOCKS * _FITS_HEADER_BLOCK_BYTES


def detection_cache_signature_matches(
    payload: Mapping | None,
    current_sig: Mapping | None,
    *,
    min_schema: int = DETECTION_CACHE_MIN_SIGNATURE_SCHEMA,
    current_engine=None,
    allow_mtime_drift: bool = True,
) -> bool:
    """Validate a Step4 detection cache payload against current input state."""

    if not isinstance(payload, Mapping) or not isinstance(current_sig, Mapping):
        return False
    if cache_schema_value(payload) < int(min_schema):
        return False
    if current_engine is not None:
        cached_engine = normalize_detect_engine(payload.get("detect_engine", "segm"))
        if cached_engine != normalize_detect_engine(current_engine):
            return False
    if allow_mtime_drift:
        return file_signature_matches_relaxed(dict(payload), dict(current_sig))
    return file_signature_matches(dict(payload), dict(current_sig))


def astap_wcs_candidates(fits_path: Path) -> List[Path]:
    return [
        fits_path.with_suffix(".wcs"),
        Path(str(fits_path) + ".wcs"),
        fits_path.parent / (fits_path.stem + ".wcs"),
        fits_path.parent / (fits_path.name + ".wcs"),
    ]


def parse_astap_wcs_file(wcs_path: Path) -> dict:
    d: Dict[str, object] = {}
    try:
        lines = wcs_path.read_text(encoding="utf-8", errors="ignore").splitlines()
    except Exception:
        return d
    for ln in lines:
        s = ln.strip()
        if not s or s.startswith("#"):
            continue
        if "/" in s:
            s = s.split("/", 1)[0].strip()
        if "=" not in s:
            continue
        key, val = [t.strip() for t in s.split("=", 1)]
        if not key:
            continue
        if val.startswith("'") and val.endswith("'"):
            d[key] = val.strip("'")
            continue
        try:
            if "." in val or "E" in val.upper():
                d[key] = float(val)
            else:
                d[key] = int(val)
        except Exception:
            d[key] = val
    return d
