"""Shared helpers for adopting target coordinates from scanned FITS headers."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np


def select_header_target(
    headers_df,
    eligible_filenames: Iterable[str],
    preferred_filenames: Iterable[str] = (),
) -> dict | None:
    """Return the first valid header target, preferring explicitly selected rows."""
    if headers_df is None or getattr(headers_df, "empty", True):
        return None
    if "Filename" not in headers_df.columns:
        return None

    eligible = [str(name) for name in eligible_filenames if str(name)]
    eligible_set = set(eligible)
    preferred = [
        str(name)
        for name in preferred_filenames
        if str(name) and str(name) in eligible_set
    ]
    ordered = list(dict.fromkeys(preferred + eligible))
    if not ordered:
        return None

    rows_by_name = {}
    for _, row in headers_df.iterrows():
        filename = str(row.get("Filename", ""))
        if filename and filename not in rows_by_name:
            rows_by_name[filename] = row

    for filename in ordered:
        row = rows_by_name.get(filename)
        if row is None:
            continue
        try:
            ra_deg = float(row.get("RA_DEG", np.nan))
            dec_deg = float(row.get("DEC_DEG", np.nan))
        except (TypeError, ValueError):
            continue
        if not (np.isfinite(ra_deg) and np.isfinite(dec_deg)):
            continue
        if not (0.0 <= ra_deg < 360.0 and -90.0 <= dec_deg <= 90.0):
            continue

        object_name = str(row.get("OBJECT", "") or "").strip()
        if object_name.lower() in {"nan", "none", "n/a", "unknown"}:
            object_name = ""
        return {
            "filename": filename,
            "object_name": object_name or Path(filename).stem,
            "ra_deg": ra_deg,
            "dec_deg": dec_deg,
        }
    return None
