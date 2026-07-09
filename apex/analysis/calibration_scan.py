"""Auto-scan and classify raw frames for detector calibration (Step 0).

Qt-free (analysis layer).  Walks a root folder, reads FITS headers, classifies
each frame as bias / dark / flat / light by ``IMAGETYP`` (with a filename
fallback), and buckets frames by observing night, exposure, temperature and
filter — mirroring the AstralImage/AIPPI auto-detection so the GUI can show a
"scan folder → tree of nights/types → run" flow instead of manual folder picks.

Grouping keys (matched to AIPPI):
  * bias  — global pool (session-agnostic)
  * dark  — (exposure, temperature bucket[±1 °C])
  * flat  — (filter)
  * light — (filter); matched to a dark by exp+temp and a flat by filter
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from astropy.io import fits

from apex.utils.constants import (
    FITS_EXTENSIONS, FILTER_HEADER_KEYS, EXPTIME_HEADER_KEYS,
)

FRAME_TYPES = ("bias", "dark", "flat", "light")

_TEMP_KEYS = ("CCD-TEMP", "CCD_TEMP", "CCDTEMP", "SET-TEMP", "SETTEMP",
              "SENSORTEMP", "CAMTEMP", "TEMP")
_IMAGETYP_KEYS = ("IMAGETYP", "IMGTYPE", "FRAMETYP", "OBSTYPE")
_DATE_KEYS = ("DATE-OBS", "DATE", "DATEOBS")
_PATH_DATE_RE = re.compile(r"(20\d{6})")


# ---------------------------------------------------------------------------
# classification helpers
# ---------------------------------------------------------------------------

def classify_type(imagetyp: Optional[str], filename: str) -> Optional[str]:
    """Classify a frame as bias/dark/flat/light from IMAGETYP, else filename."""
    t = str(imagetyp or "").upper()
    if "BIAS" in t or "ZERO" in t:
        return "bias"
    if "DARK" in t:
        return "dark"
    if "FLAT" in t:
        return "flat"
    if "LIGHT" in t or "SCIENCE" in t or "OBJECT" in t:
        return "light"
    name = os.path.basename(str(filename)).lower()
    for key in ("bias", "dark", "flat", "light"):
        if key in name:
            return key
    return None


def _header_first(header, keys, default=None):
    for k in keys:
        if header is not None and k in header:
            val = header[k]
            if val not in (None, ""):
                return val
    return default


def night_from_path(path: str) -> str:
    """YYYYMMDD from the path (folders named ..._20260515), else ''."""
    m = None
    for part in reversed(Path(path).parts):
        m = _PATH_DATE_RE.search(part)
        if m:
            return m.group(1)
    return ""


def night_from_dateobs(dateobs: Optional[str]) -> str:
    """Observing night as YYYYMMDD from DATE-OBS, split at local noon.

    Evening and the following pre-dawn hours share one night: subtract 12 h and
    take the calendar date, so 20:00 on day N and 03:00 on day N+1 both bucket
    to night N.
    """
    if not dateobs:
        return ""
    s = str(dateobs).strip().replace("Z", "")
    for fmt in ("%Y-%m-%dT%H:%M:%S.%f", "%Y-%m-%dT%H:%M:%S",
                "%Y-%m-%dT%H:%M", "%Y-%m-%d"):
        try:
            dt = datetime.strptime(s.split(".")[0] if fmt.endswith("S") else s, fmt)
            return (dt - timedelta(hours=12)).strftime("%Y%m%d")
        except ValueError:
            continue
    return ""


def temp_bucket(temp_c: Optional[float]) -> Optional[int]:
    """Round temperature to the nearest integer °C (dark grouping)."""
    if temp_c is None:
        return None
    try:
        return int(round(float(temp_c)))
    except (TypeError, ValueError):
        return None


def _to_float(val, default=None):
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


# ---------------------------------------------------------------------------
# scan
# ---------------------------------------------------------------------------

@dataclass
class FrameInfo:
    path: str
    ftype: str                       # bias | dark | flat | light
    exp: float = 0.0
    temp: Optional[float] = None
    filt: str = ""
    night: str = ""                  # YYYYMMDD or "" (global)
    name: str = field(default="")

    @property
    def temp_bucket(self) -> Optional[int]:
        return temp_bucket(self.temp)


def read_frame_info(path: str) -> Optional[FrameInfo]:
    """Read one FITS header and return its classified :class:`FrameInfo`."""
    try:
        header = fits.getheader(path)
    except Exception:
        return None
    ftype = classify_type(_header_first(header, _IMAGETYP_KEYS), path)
    if ftype is None:
        return None
    exp = _to_float(_header_first(header, EXPTIME_HEADER_KEYS), 0.0) or 0.0
    temp = _to_float(_header_first(header, _TEMP_KEYS), None)
    filt = str(_header_first(header, FILTER_HEADER_KEYS, "") or "").strip()
    night = night_from_path(path) or night_from_dateobs(_header_first(header, _DATE_KEYS))
    return FrameInfo(path=str(path), ftype=ftype, exp=exp, temp=temp,
                     filt=filt, night=night, name=os.path.basename(path))


def find_fits(root: str) -> List[str]:
    """Recursively collect FITS files under ``root`` (sorted)."""
    out: List[str] = []
    for dirpath, _dirs, files in os.walk(root):
        for f in files:
            if os.path.splitext(f)[1] in FITS_EXTENSIONS:
                out.append(os.path.join(dirpath, f))
    return sorted(out)


def scan_folder(root: str,
                progress: Optional[Callable[[int, int, str], None]] = None,
                stop: Optional[Callable[[], bool]] = None) -> List[FrameInfo]:
    """Scan ``root`` and return classified :class:`FrameInfo` for every frame."""
    paths = find_fits(root)
    total = len(paths)
    frames: List[FrameInfo] = []
    for i, p in enumerate(paths):
        if stop is not None and stop():
            break
        if progress is not None:
            progress(i, total, os.path.basename(p))
        info = read_frame_info(p)
        if info is not None:
            frames.append(info)
    if progress is not None:
        progress(total, total, "done")
    return frames


# ---------------------------------------------------------------------------
# grouping for display + processing
# ---------------------------------------------------------------------------

def nights(frames: List[FrameInfo]) -> List[str]:
    """Sorted list of distinct nights that contain LIGHT frames, newest first."""
    ns = {f.night for f in frames if f.ftype == "light"}
    return sorted(ns, reverse=True)


def group_for_night(frames: List[FrameInfo], night: str) -> Dict:
    """Frames usable for calibrating ``night``.

    Bias is a global pool (any night). Darks/flats prefer the same night but
    fall back to the global pool if the night has none. Returns a dict with
    ``bias`` (list), ``dark`` ({(exp,temp): list}), ``flat`` ({filter: list}),
    ``light`` (list for this night).
    """
    def _pool(ftype):
        same = [f for f in frames if f.ftype == ftype and f.night == night]
        return same if same else [f for f in frames if f.ftype == ftype]

    bias = _pool("bias")
    darks: Dict[Tuple[float, Optional[int]], List[FrameInfo]] = {}
    for f in _pool("dark"):
        darks.setdefault((round(f.exp, 3), f.temp_bucket), []).append(f)
    flats: Dict[str, List[FrameInfo]] = {}
    for f in _pool("flat"):
        flats.setdefault(f.filt, []).append(f)
    lights = [f for f in frames if f.ftype == "light" and f.night == night]
    return {"bias": bias, "dark": darks, "flat": flats, "light": lights}


def match_dark(darks: Dict[Tuple[float, Optional[int]], List[FrameInfo]],
               exp: float, temp: Optional[float]) -> Optional[Tuple]:
    """Pick the dark group key best matching a light's exposure + temperature.

    Prefers exact exposure and nearest temperature; falls back to nearest
    exposure. Returns the group key or None.
    """
    if not darks:
        return None
    tb = temp_bucket(temp)
    exact = [k for k in darks if abs(k[0] - exp) < 1e-3]
    if exact:
        if tb is None:
            return exact[0]
        return min(exact, key=lambda k: abs((k[1] if k[1] is not None else tb) - tb))
    # nearest exposure, then nearest temp
    return min(darks, key=lambda k: (abs(k[0] - exp),
                                     abs((k[1] if k[1] is not None else (tb or 0)) - (tb or 0))))


def match_flat(flats: Dict[str, List[FrameInfo]], filt: str) -> Optional[str]:
    """Pick the flat whose filter matches the light's filter.

    Flats MUST match by filter — a V flat cannot correct a B light — so there is
    no "use the only flat" fallback; an unmatched filter returns None and the
    caller skips flat-fielding for that light (with a warning).
    """
    if not flats:
        return None
    if filt in flats:
        return filt
    for k in flats:                       # case-insensitive
        if k.lower() == str(filt).lower():
            return k
    return None
