"""Community-standard photometry output exporters.

Converts an APEX light-curve :class:`pandas.DataFrame` into the on-disk formats
accepted by AAVSO, ExoClock, and ExoFOP/TFOP.

APEX light-curve column conventions
-----------------------------------
The canonical raw light-curve table produced by the LC pipeline (Step 9,
:mod:`apex.analysis.light_curve.lightcurve_output_service` and
``lc/step9_lightcurve_builder.py``) has these columns::

    file, filter, date, night_id, JD, BJD_TDB, rel_time_hr,
    mag, mag_err, comp_avg, comp_err, diff_mag_raw, diff_err, airmass

The detrended table (Step 10) additionally carries ``diff_mag_corr`` and
``correction_mode`` (and an optional ``ID`` column). The exporters below accept
either table and pick the best available column (corrected differential
magnitude preferred over raw, BJD_TDB preferred over JD, etc.).

LAYER RULE: Qt-free. No ``apex.gui`` / ``PyQt5`` imports.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

# AAVSO's sentinel for "no data" in the Extended File Format.
NA = "na"

PathLike = "Path | str"


# ── shared helpers ────────────────────────────────────────────────────────────

def _as_path(out_path) -> Path:
    if out_path is None:
        raise ValueError("out_path is required")
    return Path(out_path)


def _require_dataframe(lc_table) -> pd.DataFrame:
    if not isinstance(lc_table, pd.DataFrame):
        raise TypeError(f"lc_table must be a pandas DataFrame, got {type(lc_table)!r}")
    if lc_table.empty:
        raise ValueError("lc_table is empty; nothing to export")
    return lc_table


def _first_column(df: pd.DataFrame, candidates: Sequence[str]) -> Optional[str]:
    """Return the first column name in *candidates* present in *df*."""
    for col in candidates:
        if col in df.columns:
            return col
    return None


def _pick_time_column(df: pd.DataFrame, prefer_bjd: bool) -> str:
    if prefer_bjd:
        order = ("BJD_TDB", "BJD", "bjd", "JD", "jd", "HJD", "hjd")
    else:
        order = ("JD", "jd", "BJD_TDB", "BJD", "bjd", "HJD", "hjd")
    col = _first_column(df, order)
    if col is None:
        raise ValueError(
            "No time column found (looked for BJD_TDB/BJD/JD/HJD). "
            f"Available columns: {list(df.columns)}"
        )
    return col


def _time_system_label(time_col: str) -> str:
    upper = time_col.upper()
    if upper.startswith("BJD"):
        return "BJD"
    if upper.startswith("HJD"):
        return "HJD"
    return "JD"


def _pick_mag_column(df: pd.DataFrame) -> str:
    """Prefer corrected differential mag, then raw differential, then instrumental."""
    col = _first_column(
        df,
        (
            "diff_mag_corr",
            "diff_mag",
            "mag_corr",
            "calibrated_mag",
            "diff_mag_raw",
            "mag",
            "inst_mag",
        ),
    )
    if col is None:
        raise ValueError(
            "No magnitude column found (looked for diff_mag_corr/diff_mag_raw/mag). "
            f"Available columns: {list(df.columns)}"
        )
    return col


def _pick_err_column(df: pd.DataFrame) -> Optional[str]:
    return _first_column(
        df,
        ("diff_err", "diff_err_corr", "mag_err", "err", "sigma", "diff_mag_err", "comp_err"),
    )


def _meta_get(meta: Optional[Mapping], *keys: str, default=None):
    if not meta:
        return default
    for key in keys:
        if key in meta and meta[key] not in (None, ""):
            return meta[key]
    return default


def _fmt_float(value, fmt: str, fallback: str = NA) -> str:
    try:
        v = float(value)
    except (TypeError, ValueError):
        return fallback
    if not np.isfinite(v):
        return fallback
    return format(v, fmt)


def _join_ids(ids: Optional[Iterable]) -> Optional[str]:
    if ids is None:
        return None
    if isinstance(ids, (str, bytes)):
        s = ids.decode() if isinstance(ids, bytes) else ids
        s = s.strip()
        return s or None
    parts = [str(x).strip() for x in ids if str(x).strip()]
    return ",".join(parts) if parts else None


def _software_label(meta: Optional[Mapping]) -> str:
    explicit = _meta_get(meta, "software")
    if explicit:
        return str(explicit)
    try:
        from apex import __version__

        return f"APEX {__version__}"
    except Exception:  # pragma: no cover - apex always importable here
        return "APEX"


# ── AAVSO Extended File Format ────────────────────────────────────────────────

# Canonical AAVSO Extended File Format data columns (order is REQUIRED).
# Ref: AAVSO "Extended File Format" spec, https://www.aavso.org/aavso-extended-file-format
AAVSO_COLUMNS = [
    "NAME",
    "DATE",
    "MAG",
    "MERR",
    "FILT",
    "TRANS",
    "MTYPE",
    "CNAME",
    "CMAG",
    "KNAME",
    "KMAG",
    "AMASS",
    "GROUP",
    "CHART",
    "NOTES",
]


def export_aavso_extended(lc_table, meta: Optional[Mapping], out_path) -> Path:
    """Write an AAVSO Extended File Format file (AAVSO Visual/CCD observations).

    The AAVSO Extended File Format is a comma-delimited text file with a small
    header block of ``#KEY=value`` lines followed by one data row per
    observation. APEX writes the six standard header keys (``#TYPE=EXTENDED``,
    ``#OBSCODE``, ``#SOFTWARE``, ``#DELIM``, ``#DATE``, ``#OBSTYPE``) and the 15
    required data columns: NAME, DATE, MAG, MERR, FILT, TRANS, MTYPE, CNAME,
    CMAG, KNAME, KMAG, AMASS, GROUP, CHART, NOTES. Unknown fields use ``na``
    per the spec.

    Parameters
    ----------
    lc_table : pandas.DataFrame
        APEX light-curve table (raw or detrended). Must contain a time column
        (``BJD_TDB``/``JD``/``HJD``) and a magnitude column
        (``diff_mag_corr``/``diff_mag_raw``/``mag``).
    meta : Mapping
        Submission metadata. Recognized keys:
        ``obscode`` (AAVSO observer code; required by AAVSO),
        ``target`` / ``name`` (star designation; required),
        ``filter`` (AAVSO band code; falls back to the table's ``filter`` column),
        ``software`` (defaults to ``"APEX {version}"``),
        ``obstype`` (default ``"CCD"``),
        ``mtype`` (magnitude type, default ``"STD"``),
        ``trans`` (transformed, ``"YES"``/``"NO"``; default ``"NO"``),
        ``cname``/``cmag`` (comparison star name/mag),
        ``kname``/``kmag`` (check star name/mag),
        ``chart`` (AAVSO chart/sequence id),
        ``group`` (observation grouping id),
        ``notes`` (free text appended to each row).
    out_path : Path | str
        Destination file (created/overwritten).

    Returns
    -------
    Path
        The written file path.

    Raises
    ------
    ValueError / TypeError
        On empty/invalid input or missing required columns.
    """
    df = _require_dataframe(lc_table)
    out_path = _as_path(out_path)

    target = _meta_get(meta, "target", "name", "star", "object")
    if not target:
        raise ValueError("meta['target'] (star designation) is required for AAVSO export")
    obscode = _meta_get(meta, "obscode", "obs_code", default="XXX")

    time_col = _pick_time_column(df, prefer_bjd=False)  # AAVSO accepts JD or BJD
    mag_col = _pick_mag_column(df)
    err_col = _pick_err_column(df)
    filt_default = _meta_get(meta, "filter", "filt", "band")
    has_filter_col = "filter" in df.columns

    obstype = str(_meta_get(meta, "obstype", default="CCD"))
    mtype = str(_meta_get(meta, "mtype", default="STD"))
    trans = str(_meta_get(meta, "trans", "transformed", default="NO"))
    cname = _meta_get(meta, "cname", "comp_name", default=None)
    cmag = _meta_get(meta, "cmag", "comp_mag", default=None)
    kname = _meta_get(meta, "kname", "check_name", default=None)
    kmag = _meta_get(meta, "kmag", "check_mag", default=None)
    chart = _meta_get(meta, "chart", "sequence", default=None)
    group_default = _meta_get(meta, "group", "grouping", default=None)
    notes_default = _meta_get(meta, "notes", "comment", default=None)
    date_system = _time_system_label(time_col)

    header_lines = [
        "#TYPE=EXTENDED",
        f"#OBSCODE={obscode}",
        f"#SOFTWARE={_software_label(meta)}",
        "#DELIM=,",
        f"#DATE={date_system}",
        f"#OBSTYPE={obstype}",
    ]

    rows: list[list[str]] = []
    for _, row in df.iterrows():
        t = row.get(time_col)
        m = row.get(mag_col)
        try:
            if not (np.isfinite(float(t)) and np.isfinite(float(m))):
                continue
        except (TypeError, ValueError):
            continue

        if has_filter_col and pd.notna(row.get("filter")):
            filt = str(row["filter"]).strip()
        else:
            filt = str(filt_default) if filt_default else NA
        filt = filt or NA

        rows.append(
            [
                str(target),
                _fmt_float(t, ".5f"),
                _fmt_float(m, ".4f"),
                _fmt_float(row.get(err_col) if err_col else None, ".4f"),
                filt,
                trans,
                mtype,
                str(cname) if cname else NA,
                _fmt_float(cmag, ".4f") if cmag is not None else NA,
                str(kname) if kname else NA,
                _fmt_float(kmag, ".4f") if kmag is not None else NA,
                _fmt_float(row.get("airmass"), ".4f") if "airmass" in df.columns else NA,
                str(group_default) if group_default else NA,
                str(chart) if chart else NA,
                str(notes_default) if notes_default else NA,
            ]
        )

    if not rows:
        raise ValueError("No finite (time, magnitude) rows to export after filtering")

    lines = list(header_lines)
    lines.append(",".join(AAVSO_COLUMNS))
    lines.extend(",".join(r) for r in rows)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path


# ── ExoClock ──────────────────────────────────────────────────────────────────

EXOCLOCK_COLUMNS = ["BJD_TDB", "flux", "flux_err"]


def export_exoclock(lc_table, meta: Optional[Mapping], out_path) -> Path:
    """Write an ExoClock-style transit photometry file.

    ExoClock (https://www.exoclock.eu) ingests time-series transit photometry as
    a comma-delimited table of ``BJD_TDB, flux, flux_err`` (normalized,
    out-of-transit flux ≈ 1), with a small comment header describing the target,
    observer, telescope, and filter. APEX converts its differential magnitudes
    to relative flux via ``flux = 10**(-0.4 * (mag - mag_ref))`` where ``mag_ref``
    is the median magnitude (so the baseline normalizes to ~1.0), and propagates
    the magnitude error to flux error via ``flux_err = ln(10)/2.5 * flux * mag_err``.

    Parameters
    ----------
    lc_table : pandas.DataFrame
        APEX light-curve table. Requires a BJD/JD time column and a magnitude
        column. ExoClock requires BJD_TDB; if only JD is present the JD values
        are written under the BJD_TDB column and a warning note is added to the
        header (the caller is responsible for a proper barycentric correction).
    meta : Mapping
        Recognized keys: ``target``/``name`` (required), ``filter``,
        ``observer``, ``telescope``, ``obscode``, ``software``,
        ``exposure``/``exptime``.
    out_path : Path | str
        Destination file.

    Returns
    -------
    Path
    """
    df = _require_dataframe(lc_table)
    out_path = _as_path(out_path)

    target = _meta_get(meta, "target", "name", "star", "object")
    if not target:
        raise ValueError("meta['target'] is required for ExoClock export")

    time_col = _pick_time_column(df, prefer_bjd=True)
    mag_col = _pick_mag_column(df)
    err_col = _pick_err_column(df)

    times = pd.to_numeric(df[time_col], errors="coerce").to_numpy(float)
    mags = pd.to_numeric(df[mag_col], errors="coerce").to_numpy(float)
    if err_col is not None:
        merr = pd.to_numeric(df[err_col], errors="coerce").to_numpy(float)
    else:
        merr = np.full(mags.shape, np.nan)

    finite = np.isfinite(times) & np.isfinite(mags)
    if not finite.any():
        raise ValueError("No finite (time, magnitude) rows to export after filtering")
    times, mags, merr = times[finite], mags[finite], merr[finite]

    # Sort by time for a deterministic, monotonic transit series.
    order = np.argsort(times, kind="stable")
    times, mags, merr = times[order], mags[order], merr[order]

    mag_ref = float(np.nanmedian(mags))
    flux = np.power(10.0, -0.4 * (mags - mag_ref))
    ln10_over_2p5 = np.log(10.0) / 2.5
    flux_err = ln10_over_2p5 * flux * merr  # NaN where merr is NaN

    header = [
        f"# target={target}",
        f"# software={_software_label(meta)}",
        f"# time_system={_time_system_label(time_col)}",
    ]
    for key, label in (
        ("filter", "filter"),
        ("observer", "observer"),
        ("telescope", "telescope"),
        ("obscode", "obscode"),
        ("exposure", "exposure"),
        ("exptime", "exposure"),
    ):
        val = _meta_get(meta, key)
        if val is not None and f"# {label}=" not in "".join(header):
            header.append(f"# {label}={val}")
    if _time_system_label(time_col) != "BJD":
        header.append(
            "# WARNING: time column is not BJD_TDB; values copied verbatim - "
            "apply a barycentric correction before submission"
        )

    lines = list(header)
    lines.append(",".join(EXOCLOCK_COLUMNS))
    for t, fl, fe in zip(times, flux, flux_err):
        lines.append(
            ",".join(
                [
                    _fmt_float(t, ".6f"),
                    _fmt_float(fl, ".6f"),
                    _fmt_float(fe, ".6f"),
                ]
            )
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path


# ── ExoFOP / TFOP ─────────────────────────────────────────────────────────────

EXOFOP_COLUMNS = ["BJD_TDB", "rel_flux", "rel_flux_err", "airmass"]


def export_exofop(lc_table, meta: Optional[Mapping], out_path) -> Path:
    """Write an ExoFOP/TFOP-SG1 seeing-limited photometry table.

    Spec note
    ---------
    The APEX interop spec (``.wf_specs/interop_citability.md``) sketches ExoFOP
    as an ISO-8601-timestamped magnitude table but does NOT pin down the
    authoritative TFOP measurements format, which in practice is the
    AstroImageJ "Measurements" table (a tab/space-delimited file keyed on
    ``BJD_TDB`` with ``rel_flux_T1`` / ``rel_flux_err_T1`` columns). To stay
    useful and deterministic without over-claiming conformance, APEX emits the
    minimal, machine-readable subset that TFOP SG1 tooling actually consumes:
    ``BJD_TDB, rel_flux, rel_flux_err, airmass`` with a descriptive comment
    header. Treat this as a TFOP-compatible relative-flux table rather than a
    byte-for-byte AIJ Measurements file.

    Parameters
    ----------
    lc_table : pandas.DataFrame
        APEX light-curve table (needs a BJD/JD time column and a magnitude
        column; ``airmass`` is used when present).
    meta : Mapping
        Recognized keys: ``target``/``name`` (required), ``filter``,
        ``observer``, ``telescope``, ``obscode``, ``aperture_radius``,
        ``software``.
    out_path : Path | str

    Returns
    -------
    Path
    """
    df = _require_dataframe(lc_table)
    out_path = _as_path(out_path)

    target = _meta_get(meta, "target", "name", "star", "object")
    if not target:
        raise ValueError("meta['target'] is required for ExoFOP export")

    time_col = _pick_time_column(df, prefer_bjd=True)
    mag_col = _pick_mag_column(df)
    err_col = _pick_err_column(df)

    times = pd.to_numeric(df[time_col], errors="coerce").to_numpy(float)
    mags = pd.to_numeric(df[mag_col], errors="coerce").to_numpy(float)
    if err_col is not None:
        merr = pd.to_numeric(df[err_col], errors="coerce").to_numpy(float)
    else:
        merr = np.full(mags.shape, np.nan)
    if "airmass" in df.columns:
        amass = pd.to_numeric(df["airmass"], errors="coerce").to_numpy(float)
    else:
        amass = np.full(mags.shape, np.nan)

    finite = np.isfinite(times) & np.isfinite(mags)
    if not finite.any():
        raise ValueError("No finite (time, magnitude) rows to export after filtering")
    times, mags, merr, amass = times[finite], mags[finite], merr[finite], amass[finite]

    order = np.argsort(times, kind="stable")
    times, mags, merr, amass = times[order], mags[order], merr[order], amass[order]

    mag_ref = float(np.nanmedian(mags))
    rel_flux = np.power(10.0, -0.4 * (mags - mag_ref))
    ln10_over_2p5 = np.log(10.0) / 2.5
    rel_flux_err = ln10_over_2p5 * rel_flux * merr

    header = [
        f"# target={target}",
        f"# software={_software_label(meta)}",
        f"# time_system={_time_system_label(time_col)}",
        "# format=TFOP-compatible relative-flux table (see export_exofop docstring)",
    ]
    for key, label in (
        ("filter", "filter"),
        ("observer", "observer"),
        ("telescope", "telescope"),
        ("obscode", "obscode"),
        ("aperture_radius", "aperture_radius"),
    ):
        val = _meta_get(meta, key)
        if val is not None:
            header.append(f"# {label}={val}")

    lines = list(header)
    lines.append(",".join(EXOFOP_COLUMNS))
    for t, fl, fe, am in zip(times, rel_flux, rel_flux_err, amass):
        lines.append(
            ",".join(
                [
                    _fmt_float(t, ".6f"),
                    _fmt_float(fl, ".6f"),
                    _fmt_float(fe, ".6f"),
                    _fmt_float(am, ".4f"),
                ]
            )
        )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path
