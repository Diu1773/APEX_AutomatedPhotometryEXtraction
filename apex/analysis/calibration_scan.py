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
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from astropy.io import fits

from apex.utils.constants import (
    FITS_EXTENSIONS, FILTER_HEADER_KEYS, EXPTIME_HEADER_KEYS,
)
from apex.utils.night_utils import (
    DATE_HEADER_KEYS as _DATE_KEYS,
    LON_HEADER_KEYS as _LON_KEYS,
    has_local_reference,
    observing_night_detail,
    parse_lon_east as _parse_lon_east,
)

FRAME_TYPES = ("bias", "dark", "flat", "light")

_TEMP_KEYS = ("CCD-TEMP", "CCD_TEMP", "CCDTEMP", "SET-TEMP", "SETTEMP",
              "SENSORTEMP", "CAMTEMP", "TEMP")
_IMAGETYP_KEYS = ("IMAGETYP", "IMGTYPE", "FRAMETYP", "OBSTYPE")
_PATH_DATE_RE = re.compile(r"(20\d{6})")

NIGHT_METHOD_PATH = "path"


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


def night_from_dateobs(dateobs: Optional[str],
                       lon_east_deg: Optional[float] = None,
                       tz_offset_hours: Optional[float] = None) -> str:
    """Observing night as YYYYMMDD from DATE-OBS, split at local noon.

    Thin wrapper over :mod:`apex.utils.night_utils` — the single definition of
    an observing night shared with LC Step 1.  Evening and the following
    pre-dawn hours share one night because the cut falls at local noon (while
    the Sun is up), never at midnight.
    """
    return observing_night_detail(dateobs, lon_east_deg, tz_offset_hours)[0]


def resolve_night(path: str, header,
                  tz_offset_hours: Optional[float] = None) -> Tuple[str, str, str]:
    """Decide one frame's observing night.

    Returns ``(night, method, conflicting_path_date)``.  Priority:

    1. ``DATE-OBS`` + site longitude (header)      → local solar noon split
    2. ``DATE-OBS`` + configured tz offset         → local civil noon split
    3. a ``YYYYMMDD`` date in the file path
    4. ``DATE-OBS`` − 12 h (Greenwich rule, last resort)

    Steps 1-2 outrank the path date because capture software often stamps the
    *next* day onto frames taken after midnight, which would tear one night in
    two.  But they are only reachable when a real local reference exists — with
    neither longitude nor tz offset the noon split degenerates to the Greenwich
    rule, which mis-buckets East-Asian evenings, so the path date is preferred
    over it.  ``conflicting_path_date`` is set when 1-2 succeeded and the path
    disagrees, so the caller can warn.
    """
    lon = _parse_lon_east(_header_first(header, _LON_KEYS))
    dateobs = _header_first(header, _DATE_KEYS)
    path_night = night_from_path(path)

    if has_local_reference(lon, tz_offset_hours):
        night, method = observing_night_detail(dateobs, lon, tz_offset_hours)
        if night:
            conflict = path_night if (path_night and path_night != night) else ""
            return night, method, conflict

    if path_night:
        return path_night, NIGHT_METHOD_PATH, ""

    night, method = observing_night_detail(dateobs, lon, tz_offset_hours)
    return night, method, ""


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
    is_master: bool = False          # pre-built master (IMAGETYP "MASTER …")
    night_method: str = ""           # solar | civil | path | utc
    night_conflict: str = ""         # path date that disagrees with 1-2, if any

    @property
    def temp_bucket(self) -> Optional[int]:
        return temp_bucket(self.temp)


def read_frame_info(path: str,
                    tz_offset_hours: Optional[float] = None) -> Optional[FrameInfo]:
    """Read one FITS header and return its classified :class:`FrameInfo`.

    ``tz_offset_hours`` (``params.P.site_tz_offset_hours``) is the fallback
    local reference for the observing-night split when the header carries no
    site longitude — see :func:`resolve_night`.
    """
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
    night, night_method, night_conflict = resolve_night(path, header, tz_offset_hours)
    is_master = "MASTER" in str(_header_first(header, _IMAGETYP_KEYS) or "").upper()
    return FrameInfo(path=str(path), ftype=ftype, exp=exp, temp=temp,
                     filt=filt, night=night, name=os.path.basename(path),
                     is_master=is_master, night_method=night_method,
                     night_conflict=night_conflict)


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
                stop: Optional[Callable[[], bool]] = None,
                tz_offset_hours: Optional[float] = None,
                warn: Optional[Callable[[str], None]] = None) -> List[FrameInfo]:
    """Scan ``root`` and return classified :class:`FrameInfo` for every frame.

    ``warn`` receives one line per scan-level anomaly (currently: frames whose
    path date disagrees with the night derived from DATE-OBS, which is normal
    for post-midnight frames but worth showing once).
    """
    paths = find_fits(root)
    total = len(paths)
    frames: List[FrameInfo] = []
    for i, p in enumerate(paths):
        if stop is not None and stop():
            break
        if progress is not None:
            progress(i, total, os.path.basename(p))
        info = read_frame_info(p, tz_offset_hours)
        if info is not None:
            frames.append(info)
    if progress is not None:
        progress(total, total, "done")
    if warn is not None:
        conflicts = [f for f in frames if f.night_conflict]
        if conflicts:
            sample = conflicts[0]
            warn(f"[night] {len(conflicts)} frame(s): path date differs from "
                 f"DATE-OBS night (e.g. {sample.name}: path {sample.night_conflict} "
                 f"→ night {sample.night}). Using DATE-OBS ({sample.night_method}).")
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
