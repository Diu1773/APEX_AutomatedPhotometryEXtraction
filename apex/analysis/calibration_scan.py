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

import json
import os
import re
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from astropy.io import fits

from apex.utils.common_helpers import normalize_filter_key
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
# Frames whose type could not be read. They used to be dropped from the scan
# without a word, so a folder of 300 files could show 280 and the user had no
# way to tell, let alone fix it.
TYPE_UNKNOWN = "unknown"
OVERRIDES_FILENAME = "classification_overrides.json"

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


def _to_int(val) -> Optional[int]:
    try:
        return int(round(float(val)))
    except (TypeError, ValueError):
        return None


def _header_shape(header) -> Optional[Tuple[int, int]]:
    """Detector geometry as ``(NAXIS2, NAXIS1)`` — numpy order, no data read."""
    ny = _to_int(_header_first(header, ("NAXIS2",)))
    nx = _to_int(_header_first(header, ("NAXIS1",)))
    return (ny, nx) if (ny and nx) else None


def _header_binning(header) -> Optional[Tuple[int, int]]:
    bx = _to_int(_header_first(header, ("XBINNING", "BINX", "CCDXBIN", "BINNING")))
    by = _to_int(_header_first(header, ("YBINNING", "BINY", "CCDYBIN")))
    if bx and not by:
        by = bx                       # square binning written once
    return (bx, by) if (bx and by) else None


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
    shape: Optional[Tuple[int, int]] = None    # (NAXIS2, NAXIS1)
    binning: Optional[Tuple[int, int]] = None  # (XBINNING, YBINNING)

    @property
    def temp_bucket(self) -> Optional[int]:
        return temp_bucket(self.temp)

    @property
    def geometry_label(self) -> str:
        if self.shape:
            base = f"{self.shape[1]}×{self.shape[0]}"
        elif self.binning:
            base = ""
        else:
            return ""
        if self.binning:
            binned = f"bin{self.binning[0]}×{self.binning[1]}"
            return f"{base} {binned}".strip()
        return base


def compatible_geometry(a, b) -> bool:
    """Can these two frames be combined pixel-for-pixel?

    False only when both sides declare a geometry and it differs — a 2×2-binned
    dark cannot be subtracted from a 1×1 light (numpy would raise mid-run, or
    silently broadcast), so the mismatch has to be caught while matching.
    Unknown geometry never blocks a match.
    """
    shape_a, shape_b = getattr(a, "shape", None), getattr(b, "shape", None)
    if shape_a and shape_b and tuple(shape_a) != tuple(shape_b):
        return False
    bin_a, bin_b = getattr(a, "binning", None), getattr(b, "binning", None)
    if bin_a and bin_b and tuple(bin_a) != tuple(bin_b):
        return False
    return True


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
    # An unreadable IMAGETYP yields TYPE_UNKNOWN rather than dropping the frame:
    # the GUI shows those in their own bucket so they can be reclassified by
    # hand instead of vanishing from the scan.
    ftype = classify_type(_header_first(header, _IMAGETYP_KEYS), path) or TYPE_UNKNOWN
    exp = _to_float(_header_first(header, EXPTIME_HEADER_KEYS), 0.0) or 0.0
    temp = _to_float(_header_first(header, _TEMP_KEYS), None)
    filt = str(_header_first(header, FILTER_HEADER_KEYS, "") or "").strip()
    night, night_method, night_conflict = resolve_night(path, header, tz_offset_hours)
    is_master = "MASTER" in str(_header_first(header, _IMAGETYP_KEYS) or "").upper()
    return FrameInfo(path=str(path), ftype=ftype, exp=exp, temp=temp,
                     filt=filt, night=night, name=os.path.basename(path),
                     is_master=is_master, night_method=night_method,
                     night_conflict=night_conflict,
                     shape=_header_shape(header), binning=_header_binning(header))


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
        unknown = [f for f in frames if f.ftype == TYPE_UNKNOWN]
        if unknown:
            warn(f"[type] {len(unknown)} frame(s) could not be classified "
                 f"(no usable IMAGETYP, no keyword in the filename) — "
                 f"e.g. {unknown[0].name}. Set their type in the tree to use them.")
    return frames


# ---------------------------------------------------------------------------
# manual reclassification
# ---------------------------------------------------------------------------

def load_overrides(path) -> Dict[str, Dict]:
    """Read the saved manual reclassifications, keyed by frame path."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except Exception:
        return {}
    return data.get("overrides", {}) if isinstance(data, dict) else {}


def save_overrides(path, overrides: Dict[str, Dict]) -> None:
    """Persist manual reclassifications so a re-scan keeps them."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"overrides": overrides}, indent=2,
                               ensure_ascii=False), encoding="utf-8")


def apply_overrides(frames: List[FrameInfo],
                    overrides: Dict[str, Dict]) -> List[FrameInfo]:
    """Return ``frames`` with any saved manual type/filter/night applied."""
    if not overrides:
        return frames
    out: List[FrameInfo] = []
    for frame in frames:
        override = overrides.get(frame.path)
        if not override:
            out.append(frame)
            continue
        ftype = str(override.get("ftype") or frame.ftype)
        filt = override.get("filt")
        night = override.get("night")
        out.append(replace(
            frame,
            ftype=ftype if ftype in FRAME_TYPES or ftype == TYPE_UNKNOWN else frame.ftype,
            filt=frame.filt if filt is None else str(filt),
            night=frame.night if night is None else str(night),
        ))
    return out


# ---------------------------------------------------------------------------
# grouping for display + processing
# ---------------------------------------------------------------------------

def nights(frames: List[FrameInfo]) -> List[str]:
    """Sorted list of distinct nights that contain LIGHT frames, newest first."""
    ns = {f.night for f in frames if f.ftype == "light"}
    return sorted(ns, reverse=True)


def _by_exp_temp(frames: List[FrameInfo]) -> Dict[Tuple[float, Optional[int]], List[FrameInfo]]:
    groups: Dict[Tuple[float, Optional[int]], List[FrameInfo]] = {}
    for f in frames:
        groups.setdefault((round(f.exp, 3), f.temp_bucket), []).append(f)
    return groups


def group_for_night(frames: List[FrameInfo], night: str) -> Dict:
    """Frames usable for calibrating ``night``.

    Bias is a global pool (any night). Darks/flats prefer the same night but
    fall back to the global pool if the night has none. Returns a dict with
    ``bias`` (list), ``dark`` ({(exp,temp): list}), ``flat`` ({filter: list}),
    ``light`` (list for this night).

    ``dark_library`` holds the darks that are *not* from this night, kept
    separate so a same-night dark of the wrong exposure or temperature cannot
    shut out a better one from a shared dark library: preferring the night is
    right, but only among matches that are otherwise equally good.  Callers
    pass it to :func:`match_dark_detail` as ``fallback``.  The pools stay
    separate rather than merged so a master is never stacked from frames taken
    months apart.
    """
    def _pool(ftype):
        same = [f for f in frames if f.ftype == ftype and f.night == night]
        return same if same else [f for f in frames if f.ftype == ftype]

    bias = _pool("bias")
    darks = _by_exp_temp(_pool("dark"))
    library = _by_exp_temp([f for f in frames
                            if f.ftype == "dark" and f.night != night])
    flats: Dict[str, List[FrameInfo]] = {}
    for f in _pool("flat"):
        flats.setdefault(f.filt, []).append(f)
    lights = [f for f in frames if f.ftype == "light" and f.night == night]
    return {"bias": bias, "dark": darks, "dark_library": library,
            "flat": flats, "light": lights}


@dataclass(frozen=True)
class DarkMatch:
    """The dark group chosen for one light, with how far off the match is."""

    key: Tuple[float, Optional[int]]
    delta_temp_c: Optional[float]     # |T_dark - T_light|, None if either unknown
    delta_exp_s: float                # |exp_dark - exp_light|
    within_temp_tol: bool             # False -> caller warns (or refuses)
    frames: Tuple[FrameInfo, ...] = ()   # the darks themselves
    source: str = "night"             # "night" | "library"

    @property
    def exp(self) -> float:
        return self.key[0]

    @property
    def temp_bucket(self) -> Optional[int]:
        return self.key[1]

    @property
    def night(self) -> str:
        return self.frames[0].night if self.frames else ""

    def rank(self) -> Tuple:
        """Sort key for "which of two matches is better" (smaller is better)."""
        return (
            0 if self.delta_exp_s <= 1e-3 else 1,     # exposure-exact first
            0 if self.within_temp_tol else 1,         # then temperature
            self.delta_exp_s,
            self.delta_temp_c if self.delta_temp_c is not None else 0.0,
        )


def group_temperature(frames: List[FrameInfo],
                      key: Tuple[float, Optional[int]]) -> Optional[float]:
    """Mean *actual* sensor temperature of a dark group.

    The grouping key only carries the 1 °C display bucket, which is far coarser
    than observers who care about dark current can accept, so matching reads the
    frames' real temperatures and falls back to the bucket only when they are
    missing.
    """
    temps: List[float] = []
    for f in frames or ():
        value = getattr(f, "temp", None)
        if value is None:
            continue
        try:
            temps.append(float(value))
        except (TypeError, ValueError):
            continue
    if temps:
        return sum(temps) / len(temps)
    return float(key[1]) if key[1] is not None else None


def _geometry_filtered(groups: Dict, light: Optional[FrameInfo]) -> Dict:
    """Drop calibration groups that cannot be combined with ``light``.

    Returns the original mapping when the light's geometry is unknown or when
    nothing survives (so the caller reports "no dark/flat" rather than picking
    a frame that would blow up on subtraction).
    """
    if light is None:
        return groups
    kept = {k: v for k, v in groups.items()
            if any(compatible_geometry(light, f) for f in (v or ()))}
    return kept


def _match_one_pool(darks, exp, temp, tol_c, light, source) -> Optional[DarkMatch]:
    if not darks:
        return None
    darks = _geometry_filtered(darks, light)
    if not darks:
        return None

    def _delta_temp(key) -> Optional[float]:
        group_temp = group_temperature(darks.get(key), key)
        if group_temp is None or temp is None:
            return None
        try:
            return abs(group_temp - float(temp))
        except (TypeError, ValueError):
            return None

    exact = [k for k in darks if abs(k[0] - exp) < 1e-3]
    pool = exact or list(darks)
    best = min(pool, key=lambda k: (
        abs(k[0] - exp),
        float("inf") if _delta_temp(k) is None else _delta_temp(k),
    ))
    delta_temp = _delta_temp(best)
    return DarkMatch(
        key=best,
        delta_temp_c=delta_temp,
        delta_exp_s=abs(best[0] - exp),
        within_temp_tol=(delta_temp is None or delta_temp <= float(tol_c)),
        frames=tuple(darks.get(best) or ()),
        source=source,
    )


def match_dark_detail(darks: Dict[Tuple[float, Optional[int]], List[FrameInfo]],
                      exp: float, temp: Optional[float],
                      tol_c: float = 1.0,
                      light: Optional[FrameInfo] = None,
                      fallback: Optional[Dict] = None) -> Optional[DarkMatch]:
    """Pick the dark group best matching a light's exposure + temperature.

    Prefers an exact exposure, then the nearest *actual* temperature. The
    returned :class:`DarkMatch` reports the residual mismatch so the caller can
    show it, log it, or refuse it (``strict_temp``) — silently accepting an
    arbitrarily large ΔT was the old behaviour.

    When ``light`` is given, darks of an incompatible geometry (different frame
    size or binning) are excluded; ``None`` is returned if that leaves nothing.

    ``fallback`` (``group_for_night()["dark_library"]``) is consulted when it
    offers a strictly better match than the night's own darks — a same-night
    dark of the wrong exposure or temperature must not shut out a shared dark
    library that has the right one. Ties go to the night's own darks.
    """
    primary = _match_one_pool(darks, exp, temp, tol_c, light, "night")
    spare = _match_one_pool(fallback, exp, temp, tol_c, light, "library")
    if primary is None:
        return spare
    if spare is not None and spare.rank() < primary.rank():
        return spare
    return primary


def match_dark(darks: Dict[Tuple[float, Optional[int]], List[FrameInfo]],
               exp: float, temp: Optional[float]) -> Optional[Tuple]:
    """Group key of the best dark match (see :func:`match_dark_detail`)."""
    match = match_dark_detail(darks, exp, temp)
    return match.key if match is not None else None


def match_flat(flats: Dict[str, List[FrameInfo]], filt: str,
               light: Optional[FrameInfo] = None) -> Optional[str]:
    """Pick the flat whose filter matches the light's filter.

    Flats MUST match by filter — a V flat cannot correct a B light — so there is
    no "use the only flat" fallback; an unmatched filter returns None and the
    caller skips flat-fielding for that light (with a warning).

    Comparison goes through :func:`normalize_filter_key`, the same canonical
    form the rest of the pipeline uses, so a header writing ``HA`` for the
    flats and ``Ha`` for the lights still matches instead of reporting
    "NO FLAT".  ``light`` additionally excludes flats of an incompatible
    geometry (different frame size or binning).
    """
    if not flats:
        return None
    flats = _geometry_filtered(flats, light)
    if not flats:
        return None
    if filt in flats:
        return filt
    want = normalize_filter_key(filt)
    if want:
        for k in flats:
            if normalize_filter_key(k) == want:
                return k
    for k in flats:                       # last resort: plain case-insensitive
        if str(k).lower() == str(filt).lower():
            return k
    return None
