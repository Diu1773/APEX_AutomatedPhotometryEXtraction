"""CCD/CMOS overscan correction at the raw calibration boundary.

Vendored (near-verbatim) from the AstralImage/AIPPI engine
(``core/overscan.py``) as part of the APEX detector-calibration Step 0.  Pure
numpy, no GUI/engine coupling — belongs in the analysis layer.  Only the
science-legitimate row/column overscan subtraction is kept; nothing here is
imaging-specific.
"""

from __future__ import annotations

import re

import numpy as np


VALID_EDGES = ("left", "right", "top", "bottom")


_SECTION_RE = re.compile(
    r"^\s*\[\s*([+-]?\d+)\s*:\s*([+-]?\d+)\s*,\s*"
    r"([+-]?\d+)\s*:\s*([+-]?\d+)\s*\]\s*$"
)


def parse_fits_section(value, shape=None):
    """Parse a FITS image section like ``[x1:x2,y1:y2]``.

    FITS sections are 1-based and inclusive.  Returns half-open Python
    ``(y0, y1, x0, x1)`` coordinates, or ``None`` if unparseable/out of range.
    """
    if value is None:
        return None
    match = _SECTION_RE.match(str(value).strip())
    if not match:
        return None
    x1, x2, y1, y2 = (int(group) for group in match.groups())
    x0 = min(x1, x2) - 1
    xh = max(x1, x2)
    y0 = min(y1, y2) - 1
    yh = max(y1, y2)
    if x0 < 0 or y0 < 0 or xh <= x0 or yh <= y0:
        return None
    if shape is not None:
        h, w = shape
        if xh > w or yh > h:
            return None
    return y0, yh, x0, xh


def _header_value(header, key: str):
    if header is None:
        return None
    try:
        return header.get(key)
    except AttributeError:
        try:
            return header[key]
        except Exception:
            return None


def _infer_edge(bsec, dsec) -> str:
    by0, by1, bx0, bx1 = bsec
    dy0, dy1, dx0, dx1 = dsec
    if bx0 >= dx1:
        return "right"
    if bx1 <= dx0:
        return "left"
    if by0 >= dy1:
        return "bottom"
    if by1 <= dy0:
        return "top"
    return "section"


def correct_overscan(data: np.ndarray, *, edge: str = "right",
                     width: int = 32, trim: bool = True):
    """Subtract row/column overscan bias and optionally trim the region.

    Left/right regions estimate one robust bias level per sensor row; top/bottom
    regions estimate one per column.  Removes the electronic pedestal and slow
    readout banding before master construction.

    Returns ``(corrected, info)``.  Non-2D or invalid configs return a copy with
    ``info['applied'] == False``.
    """
    arr = np.asarray(data)
    info = {
        "applied": False,
        "edge": str(edge).lower(),
        "width": int(width),
        "trimmed": False,
        "level": 0.0,
        "noise": 0.0,
    }
    if arr.ndim != 2:
        return arr.copy(), info

    edge = str(edge).strip().lower()
    width = int(width)
    h, w = arr.shape
    limit = w if edge in ("left", "right") else h
    if edge not in VALID_EDGES or width < 2 or width >= limit:
        return arr.copy(), info

    work = arr.astype(np.float64, copy=True)
    if edge == "left":
        region = work[:, :width]
        levels = np.nanmedian(region, axis=1)
        work -= levels[:, np.newaxis]
        if trim:
            work = work[:, width:]
    elif edge == "right":
        region = work[:, -width:]
        levels = np.nanmedian(region, axis=1)
        work -= levels[:, np.newaxis]
        if trim:
            work = work[:, :-width]
    elif edge == "top":
        region = work[:width, :]
        levels = np.nanmedian(region, axis=0)
        work -= levels[np.newaxis, :]
        if trim:
            work = work[width:, :]
    else:
        region = work[-width:, :]
        levels = np.nanmedian(region, axis=0)
        work -= levels[np.newaxis, :]
        if trim:
            work = work[:-width, :]

    finite = region[np.isfinite(region)]
    median = float(np.nanmedian(finite)) if finite.size else 0.0
    noise = (
        float(np.nanmedian(np.abs(finite - median)) * 1.4826)
        if finite.size else 0.0
    )
    info.update({
        "applied": True,
        "edge": edge,
        "width": width,
        "trimmed": bool(trim),
        "level": float(np.nanmedian(levels)),
        "noise": noise,
    })
    return work.astype(arr.dtype, copy=False), info


def correct_overscan_from_header(data: np.ndarray, header, *, trim: bool = True):
    """Apply FITS ``BIASSEC``/``DATASEC`` overscan correction when available.

    Supports single-amplifier vertical (left/right) and horizontal (top/bottom)
    overscan strips.  More complex multi-amplifier sections deliberately return
    ``applied=False`` so the caller can fall back to explicit edge/width.
    """
    arr = np.asarray(data)
    biassec_raw = _header_value(header, "BIASSEC")
    datasec_raw = _header_value(header, "DATASEC") or _header_value(header, "TRIMSEC")
    info = {
        "applied": False,
        "source": "header",
        "edge": "section",
        "width": 0,
        "trimmed": False,
        "level": 0.0,
        "noise": 0.0,
        "biassec": str(biassec_raw or ""),
        "datasec": str(datasec_raw or ""),
    }
    if arr.ndim != 2 or biassec_raw in (None, "") or datasec_raw in (None, ""):
        return arr.copy(), info

    bsec = parse_fits_section(biassec_raw, arr.shape)
    dsec = parse_fits_section(datasec_raw, arr.shape)
    if bsec is None or dsec is None:
        return arr.copy(), info

    by0, by1, bx0, bx1 = bsec
    dy0, dy1, dx0, dx1 = dsec
    work = arr.astype(np.float64, copy=False)
    data_region = work[dy0:dy1, dx0:dx1]
    bias_region = work[by0:by1, bx0:bx1]
    if data_region.size == 0 or bias_region.size == 0:
        return arr.copy(), info

    vertical = (by0 <= dy0 and by1 >= dy1 and (bx1 - bx0) < (dx1 - dx0))
    horizontal = (bx0 <= dx0 and bx1 >= dx1 and (by1 - by0) < (dy1 - dy0))
    if vertical:
        bias_rows = bias_region[(dy0 - by0):(dy1 - by0), :]
        if bias_rows.shape[0] != data_region.shape[0]:
            return arr.copy(), info
        levels = np.nanmedian(bias_rows, axis=1)
        corrected_data = data_region - levels[:, np.newaxis]
        width = bx1 - bx0
    elif horizontal:
        bias_cols = bias_region[:, (dx0 - bx0):(dx1 - bx0)]
        if bias_cols.shape[1] != data_region.shape[1]:
            return arr.copy(), info
        levels = np.nanmedian(bias_cols, axis=0)
        corrected_data = data_region - levels[np.newaxis, :]
        width = by1 - by0
    else:
        return arr.copy(), info

    if trim:
        out = corrected_data.copy()
        trim_origin = (dx0, dy0)
    else:
        out = work.astype(np.float64, copy=True)
        out[dy0:dy1, dx0:dx1] = corrected_data
        trim_origin = (0, 0)

    finite = bias_region[np.isfinite(bias_region)]
    median = float(np.nanmedian(finite)) if finite.size else 0.0
    noise = (
        float(np.nanmedian(np.abs(finite - median)) * 1.4826)
        if finite.size else 0.0
    )
    info.update({
        "applied": True,
        "source": "header",
        "edge": _infer_edge(bsec, dsec),
        "width": int(width),
        "trimmed": bool(trim),
        "level": float(np.nanmedian(levels)),
        "noise": noise,
        "biassec": str(biassec_raw),
        "datasec": str(datasec_raw),
        "trim_origin": trim_origin,
    })
    return out.astype(arr.dtype, copy=False), info


def update_header_for_overscan(header, info: dict) -> None:
    """Record overscan provenance and adjust WCS CRPIX for leading-edge trims."""
    if header is None or not info.get("applied"):
        return
    edge = info["edge"]
    width = int(info["width"])
    if info.get("trimmed") and info.get("source") == "header":
        x0, y0 = info.get("trim_origin", (0, 0))
        if x0 and "CRPIX1" in header:
            header["CRPIX1"] = float(header["CRPIX1"]) - int(x0)
        if y0 and "CRPIX2" in header:
            header["CRPIX2"] = float(header["CRPIX2"]) - int(y0)
    else:
        if info.get("trimmed") and edge == "left" and "CRPIX1" in header:
            header["CRPIX1"] = float(header["CRPIX1"]) - width
        if info.get("trimmed") and edge == "top" and "CRPIX2" in header:
            header["CRPIX2"] = float(header["CRPIX2"]) - width
    header["OVRSCAN"] = (True, "Overscan correction applied")
    header["OVRSSRC"] = (str(info.get("source", "manual")), "Overscan source")
    header["OVRSEDGE"] = (edge, "Overscan region edge")
    header["OVRSWID"] = (width, "Overscan region width in pixels")
    header["OVRSTRIM"] = (bool(info.get("trimmed")), "Overscan region trimmed")
    if info.get("biassec"):
        header["OVRSBIAS"] = (str(info.get("biassec")), "BIASSEC used")
    if info.get("datasec"):
        header["OVRSDATA"] = (str(info.get("datasec")), "DATASEC/TRIMSEC used")
    header.add_history(
        f"OVERSCAN source={str(info.get('source', 'manual'))} "
        f"edge={edge} width={width} trim={bool(info.get('trimmed'))} "
        f"level={float(info.get('level', 0.0)):.4f}"
    )
