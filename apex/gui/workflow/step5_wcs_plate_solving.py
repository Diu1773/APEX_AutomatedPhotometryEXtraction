"""
Step 5: WCS Plate Solving (ASTAP)
WCS solving window and workers.
"""

from __future__ import annotations

import json
import sys
import time
import subprocess
import threading
import tempfile
import warnings
import shlex
import shutil
from pathlib import Path, PureWindowsPath
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord
import astropy.units as u

from scipy.spatial import cKDTree as KDTree

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QGroupBox, QMessageBox,
    QTextEdit, QComboBox, QDialog, QFormLayout, QLineEdit, QDialogButtonBox,
    QProgressBar, QCheckBox, QSpinBox, QDoubleSpinBox, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView, QWidget, QTabWidget,
    QScrollArea, QFileDialog
)

from PyQt5.QtCore import Qt, QThread, pyqtSignal

from .step_window_base import StepWindowBase
from .run_control import RunControlBar, format_duration, progress_status_text
from .log_panel import WorkflowLogWindow, WorkerStatusPanel, append_timestamped_log, show_raised
from .ui_helpers import (
    add_parameter_reset_button,
    create_cache_checkbox,
    create_collapsible_section,
    create_parameter_button,
    configure_parameter_dialog,
    set_table_row_background,
    status_row_background,
)
from apex.core.cache_manager import StepCacheManager
from apex.utils.step_paths import (
    step2_cropped_dir,
    crop_is_active,
    crop_rect_path,
    step4_dir,
    step5_wcs_dir,
    step1_dir,
)
from apex.utils.constants import get_parallel_workers, MAD_TO_SIGMA
from apex.utils.gaia_catalog_service import (
    GaiaCatalogService,
    gaia_runtime_available,
)
from apex.utils.qc_utils import filter_files_by_qc
from apex.utils.cache_utils import (
    norm_path_key,
    build_file_signature,
    detection_cache_signature_matches,
    file_signature_matches_relaxed,
    astap_wcs_candidates,
    parse_astap_wcs_file,
)


# Windows defaults subprocess text pipes to the active ANSI code page (often cp949),
# which crashes on UTF-8 output from WSL/solver tools. Decode explicitly.
# CREATE_NO_WINDOW prevents console flicker in frozen (PyInstaller) apps.
_SUBPROCESS_TEXT_KWARGS: dict = {
    "text": True,
    "encoding": "utf-8",
    "errors": "replace",
}
if sys.platform == "win32":
    _SUBPROCESS_TEXT_KWARGS["creationflags"] = subprocess.CREATE_NO_WINDOW
    _SUBPROCESS_TEXT_KWARGS["stdin"] = subprocess.DEVNULL


def _tail_text(value: str | None, limit: int = 800, max_lines: int = 8) -> str:
    if value is None:
        return ""
    s = str(value).replace("\r\n", "\n").replace("\r", "\n").strip()
    if not s:
        return ""
    lines = [ln.strip() for ln in s.split("\n") if ln.strip()]
    if max_lines > 0 and len(lines) > max_lines:
        lines = lines[-max_lines:]
    one_line = " | ".join(lines)
    if len(one_line) > limit:
        one_line = "..." + one_line[-limit:]
    return one_line


def _exc_brief(exc: Exception, limit: int = 260) -> str:
    return _tail_text(f"{type(exc).__name__}: {exc}", limit=limit, max_lines=4)


def _is_explicit_stopped_message(value: str | None) -> bool:
    s = str(value or "").strip().lower()
    return s == "stopped" or s.startswith("stopped |")


def _coord_sep_deg(a: SkyCoord | None, b: SkyCoord | None) -> float:
    if a is None or b is None:
        return float("nan")
    try:
        return float(a.separation(b).deg)
    except Exception:
        return float("nan")


def _format_coord_hint(coord: SkyCoord | None) -> str:
    if coord is None:
        return "blind"
    try:
        return f"{coord.ra.deg:.6f},{coord.dec.deg:.6f}"
    except Exception:
        return "invalid"


def _header_pointing_coord(hdr) -> SkyCoord | None:
    """Return a SkyCoord from a frame header's mount-pointing keys, or None.

    Reuses the shared ``astro_utils._parse_ra_dec_from_header`` key-priority
    parser (OBJCTRA/OBJRA/RA/RA_OBJ, sexagesimal or decimal) so the internal
    solver reads pointing the same way the rest of APEX does.
    """
    try:
        from apex.utils.astro_utils import _parse_ra_dec_from_header
        rd = _parse_ra_dec_from_header(hdr)
        if rd is None:
            return None
        ra_deg, dec_deg = rd
        if not (np.isfinite(ra_deg) and np.isfinite(dec_deg)):
            return None
        return SkyCoord(float(ra_deg) * u.deg, float(dec_deg) * u.deg)
    except Exception:
        return None


def _strip_outer_quotes(value: str) -> str:
    value = str(value or "").strip()
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {"'", '"'}:
        return value[1:-1]
    return value


def _split_command(value: str) -> list[str]:
    value = str(value or "").strip()
    if not value:
        return []
    try:
        return shlex.split(value)
    except Exception:
        return [value]


def _local_executable_available(exe: str) -> tuple[bool, str]:
    exe = _strip_outer_quotes(exe)
    if not exe:
        return False, "empty command"
    path = Path(exe)
    if path.exists():
        return True, f"found at {path}"
    found = shutil.which(exe)
    if found:
        return True, f"found on PATH: {found}"
    return False, f"not found: {exe}"


# ── APEX-managed WCS provenance/quality keys ─────────────────────────────────
# These are written by each solver to identify itself and surface quality stats.
# We clear them at the start of every per-frame solve so a re-run with a
# different solver cannot leave stale metadata behind (e.g. ASTAP RMS lingering
# in a frame re-solved by the Internal Python engine).
_APEX_WCS_META_KEYS: tuple[str, ...] = (
    "WCSSRC",   # solver tag: APEX_INTERNAL / ASTAP / ASTNET_WSL / ASTNET_REFINED
    "WCS_OK",   # bool — last solve success
    "WCSPIXI",  # input pixel scale (arcsec/px)
    "WCSPIXF",  # fitted pixel scale
    "WCSREFN",  # refine note
    "WCSRMD",   # residual median (arcsec)
    "WCSRMAX",  # residual max (arcsec)
    "WCSROT",   # rotation (deg, E of N)
    "WCSCRA",   # center RA (deg)
    "WCSCDEC",  # center Dec (deg)
    "WCSMOD",   # TAN / SIPn
    "WCSSIP",   # SIP order
    "WCSNST",   # number of matched stars
)


def _load_simbad_target_coord(result_dir) -> "SkyCoord | None":
    """Read Step 1's targets_simbad.tsv and return the first row's SkyCoord.

    Step 1 writes the SIMBAD-resolved target to
    ``<result_dir>/step1_file_selection/targets_simbad.tsv``. We trust this
    over a FITS-header OBJCTRA/OBJCTDEC because mount logs are sometimes
    wrong (e.g. a previous slew target gets stamped into the headers).
    Returns None on any failure so callers can fall back.
    """
    try:
        path = step1_dir(Path(result_dir)) / "targets_simbad.tsv"
        if not path.exists():
            return None
        df = pd.read_csv(path, sep="\t")
        if df.empty:
            return None
        cols = {c.lower(): c for c in df.columns}
        ra_col = cols.get("ra_deg")
        dec_col = cols.get("dec_deg")
        if ra_col is None or dec_col is None:
            return None
        ra = float(df[ra_col].iloc[0])
        dec = float(df[dec_col].iloc[0])
        if not (np.isfinite(ra) and np.isfinite(dec)):
            return None
        return SkyCoord(ra * u.deg, dec * u.deg)
    except Exception:
        return None


def _reset_apex_wcs_meta(hdr) -> None:
    """Remove APEX-managed WCS provenance/quality keys from *hdr* in place.

    Raw WCS keywords (CRVAL/CRPIX/CTYPE/CD/PC/CDELT/PV/...) are intentionally
    not cleared — the calling solver will overwrite them with its own values.
    """
    for key in _APEX_WCS_META_KEYS:
        try:
            if key in hdr:
                del hdr[key]
        except Exception:
            pass


def _shared_wcs_center_coords(w: WCS, nx: int, ny: int) -> tuple[float, float]:
    try:
        if not w.has_celestial:
            return (float("nan"), float("nan"))
        sky = w.pixel_to_world(nx / 2.0, ny / 2.0)
        return (float(sky.ra.deg), float(sky.dec.deg))
    except Exception:
        return (float("nan"), float("nan"))


def _shared_empty_wcs_qc_metrics(n_detect: int = 0) -> dict:
    return {
        "n_detect": int(max(0, n_detect)),
        "n_catalog_in_fov": 0,
        "n_match": 0,
        "n_inlier": 0,
        "match_rate": np.nan,
        "match_rate_cat": np.nan,
        "match_rate_eff": np.nan,
        "match_radius_arcsec": np.nan,
        "match_radius_px": np.nan,
        "dx_med_px": np.nan,
        "dy_med_px": np.nan,
        "resid_med_px": np.nan,
        "resid_mad_px": np.nan,
        "resid_peak_px": np.nan,
        "resid_p99_px": np.nan,
        "rms_px": np.nan,
        "inlier_rate": np.nan,
        "resid_vs_radius_slope": np.nan,
        "edge_resid_ratio": np.nan,
        "center_offset_arcsec": np.nan,
        "pix_scale_input_arcsec": np.nan,
        "pix_scale_fit_arcsec": np.nan,
        "scale_delta_pct": np.nan,
        "gaia_available": False,
    }


def _shared_compute_wcs_qc_metrics(
    params,
    *,
    w: WCS | None,
    det_xy: np.ndarray,
    nx: int,
    ny: int,
    gaia_ra_deg: np.ndarray,
    gaia_dec_deg: np.ndarray,
    pix_input_arcsec: float,
    pix_fit_arcsec: float,
    center_coord: SkyCoord | None,
) -> dict:
    det_xy = np.asarray(det_xy, dtype=float)
    if det_xy.ndim != 2 or det_xy.shape[1] != 2:
        det_xy = np.zeros((0, 2), dtype=float)
    out = _shared_empty_wcs_qc_metrics(n_detect=len(det_xy))
    gaia_ra_deg = np.asarray(gaia_ra_deg, dtype=float)
    gaia_dec_deg = np.asarray(gaia_dec_deg, dtype=float)
    out["gaia_available"] = bool(len(gaia_ra_deg) > 0 and len(gaia_dec_deg) > 0)
    if np.isfinite(pix_input_arcsec):
        out["pix_scale_input_arcsec"] = float(pix_input_arcsec)
    if np.isfinite(pix_fit_arcsec):
        out["pix_scale_fit_arcsec"] = float(pix_fit_arcsec)
    if np.isfinite(pix_input_arcsec) and pix_input_arcsec > 0 and np.isfinite(pix_fit_arcsec):
        out["scale_delta_pct"] = float((pix_fit_arcsec - pix_input_arcsec) / pix_input_arcsec * 100.0)

    if w is not None and w.has_celestial and center_coord is not None:
        try:
            c_ra, c_dec = _shared_wcs_center_coords(w, nx, ny)
            if np.isfinite(c_ra) and np.isfinite(c_dec):
                c_sky = SkyCoord(c_ra * u.deg, c_dec * u.deg, frame="icrs")
                out["center_offset_arcsec"] = float(c_sky.separation(center_coord).arcsec)
        except Exception:
            pass

    if w is None or (not w.has_celestial) or len(det_xy) == 0:
        return out
    if gaia_ra_deg.size == 0 or gaia_dec_deg.size == 0:
        return out

    try:
        xg, yg = w.celestial.all_world2pix(gaia_ra_deg, gaia_dec_deg, 0)
        xg = np.asarray(xg, float)
        yg = np.asarray(yg, float)
    except Exception:
        return out

    ok_g = (
        np.isfinite(xg)
        & np.isfinite(yg)
        & (xg >= 0.0)
        & (xg < float(nx))
        & (yg >= 0.0)
        & (yg < float(ny))
    )
    if not np.any(ok_g):
        return out

    gaia_xy = np.column_stack((xg[ok_g], yg[ok_g]))
    out["n_catalog_in_fov"] = int(len(gaia_xy))

    pix_use = pix_fit_arcsec if np.isfinite(pix_fit_arcsec) and pix_fit_arcsec > 0 else pix_input_arcsec
    match_r_arcsec = float(getattr(params.P, "wcs_qc_match_radius_arcsec", 2.0))
    if not np.isfinite(match_r_arcsec) or match_r_arcsec <= 0:
        match_r_arcsec = 2.0
    if np.isfinite(pix_use) and pix_use > 0:
        match_r_px = float(match_r_arcsec / pix_use)
    else:
        match_r_px = float(getattr(params.P, "wcs_qc_match_radius_px", 2.5))
    match_r_px = float(np.clip(match_r_px, 1.0, 25.0))
    out["match_radius_arcsec"] = float(match_r_arcsec)
    out["match_radius_px"] = float(match_r_px)

    tree = KDTree(gaia_xy)
    d, j = tree.query(det_xy, k=1)
    d = np.asarray(d, float)
    j = np.asarray(j, int)
    ok = np.isfinite(d) & (d <= match_r_px) & (j >= 0) & (j < len(gaia_xy))
    if not np.any(ok):
        return out

    det_candidates = np.where(ok)[0]
    order = np.argsort(d[det_candidates])
    used_gaia = set()
    keep_det = []
    keep_gaia = []
    for ord_idx in order:
        det_i = int(det_candidates[ord_idx])
        gaia_i = int(j[det_i])
        if gaia_i in used_gaia:
            continue
        used_gaia.add(gaia_i)
        keep_det.append(det_i)
        keep_gaia.append(gaia_i)
    if not keep_det:
        return out

    det_keep = np.asarray(keep_det, dtype=int)
    gaia_keep = np.asarray(keep_gaia, dtype=int)
    dx = det_xy[det_keep, 0] - gaia_xy[gaia_keep, 0]
    dy = det_xy[det_keep, 1] - gaia_xy[gaia_keep, 1]
    r = np.hypot(dx, dy)
    finite_r = np.isfinite(r)
    if not np.any(finite_r):
        return out
    if not np.all(finite_r):
        dx = dx[finite_r]
        dy = dy[finite_r]
        r = r[finite_r]
        gaia_keep = gaia_keep[finite_r]

    n_match = int(len(r))
    out["n_match"] = n_match
    out["match_rate"] = float(n_match / max(int(len(det_xy)), 1))
    out["match_rate_cat"] = float(n_match / max(int(len(gaia_xy)), 1))
    out["match_rate_eff"] = float(max(out["match_rate"], out["match_rate_cat"]))
    if n_match == 0:
        return out

    out["dx_med_px"] = float(np.nanmedian(dx)) if len(dx) else np.nan
    out["dy_med_px"] = float(np.nanmedian(dy)) if len(dy) else np.nan
    resid_med = float(np.nanmedian(r))
    resid_mad = float(MAD_TO_SIGMA * np.nanmedian(np.abs(r - resid_med)))
    resid_p99 = float(np.nanpercentile(r, 99))
    out["resid_med_px"] = resid_med
    out["resid_mad_px"] = resid_mad
    out["resid_p99_px"] = resid_p99
    out["resid_peak_px"] = resid_p99

    clip_sigma = float(getattr(params.P, "wcs_qc_clip_sigma", 3.0))
    if not np.isfinite(clip_sigma) or clip_sigma <= 0:
        clip_sigma = 3.0
    if np.isfinite(resid_mad) and resid_mad > 0:
        inlier = np.abs(r - resid_med) <= clip_sigma * resid_mad
    else:
        resid_std = float(np.nanstd(r))
        if np.isfinite(resid_std) and resid_std > 0:
            inlier = np.abs(r - float(np.nanmean(r))) <= clip_sigma * resid_std
        else:
            inlier = np.ones(len(r), dtype=bool)
    n_inlier = int(np.sum(inlier))
    r_in = r[inlier] if n_inlier > 0 else r
    out["n_inlier"] = n_inlier
    out["inlier_rate"] = float(n_inlier / max(n_match, 1))
    out["rms_px"] = float(np.sqrt(np.nanmean(r_in ** 2))) if len(r_in) else np.nan

    if len(det_keep) >= 8:
        cx = float(nx) / 2.0
        cy = float(ny) / 2.0
        rr = np.hypot(gaia_xy[gaia_keep, 0] - cx, gaia_xy[gaia_keep, 1] - cy)
        max_rr = max(float(np.hypot(max(cx, 1.0), max(cy, 1.0))), 1.0)
        rho = rr / max_rr
        if np.isfinite(np.nanstd(rho)) and float(np.nanstd(rho)) > 1e-6:
            try:
                out["resid_vs_radius_slope"] = float(np.polyfit(rho, r, 1)[0])
            except Exception:
                out["resid_vs_radius_slope"] = np.nan
        core = r[rho <= 0.4]
        edge = r[rho >= 0.8]
        if len(core) >= 3 and len(edge) >= 3:
            core_med = float(np.nanmedian(core))
            if np.isfinite(core_med) and core_med > 1e-9:
                out["edge_resid_ratio"] = float(np.nanmedian(edge) / core_med)

    return out


def _shared_evaluate_wcs_qc_pass(params, metrics: dict, *, wcs_ok: bool) -> tuple[bool, list[str]]:
    def _num(key: str) -> float:
        try:
            return float(metrics.get(key, np.nan))
        except Exception:
            return np.nan

    reasons: list[str] = []
    gaia_available = bool(metrics.get("gaia_available", True))

    require_wcs_ok = bool(getattr(params.P, "wcs_qc_require_wcs_ok", True))
    if require_wcs_ok and not wcs_ok:
        reasons.append("wcs_fail")

    n_detect = int(metrics.get("n_detect", 0) or 0)
    n_match = int(metrics.get("n_match", 0) or 0)
    if n_detect <= 0:
        reasons.append("no_detect_data")

    if not gaia_available:
        reasons.append("gaia_unavailable")
    else:
        min_match_n = int(getattr(params.P, "wcs_qc_min_match_n", 20))
        if min_match_n > 0 and n_match < min_match_n:
            reasons.append("low_match_n")

        min_match_rate = float(getattr(params.P, "wcs_qc_min_match_rate", 0.20))
        mrate_det = _num("match_rate")
        mrate_cat = _num("match_rate_cat")
        mrate_eff = _num("match_rate_eff")
        if not np.isfinite(mrate_eff):
            finite_rates = [v for v in [mrate_det, mrate_cat] if np.isfinite(v)]
            if finite_rates:
                mrate_eff = float(np.nanmax(finite_rates))
        if np.isfinite(min_match_rate) and min_match_rate > 0:
            if (not np.isfinite(mrate_eff)) or (mrate_eff < min_match_rate):
                reasons.append("low_match_rate")

    max_rms_px = float(getattr(params.P, "wcs_qc_max_rms_px", 2.5))
    if n_match > 0 and np.isfinite(max_rms_px) and max_rms_px > 0:
        rms_px = _num("rms_px")
        if (not np.isfinite(rms_px)) or (rms_px > max_rms_px):
            reasons.append("high_rms")

    max_p99_px = float(getattr(params.P, "wcs_qc_max_p99_px", 5.0))
    if n_match > 0 and np.isfinite(max_p99_px) and max_p99_px > 0:
        p99_px = _num("resid_p99_px")
        if (not np.isfinite(p99_px)) or (p99_px > max_p99_px):
            reasons.append("high_p99")

    min_inlier_rate = float(getattr(params.P, "wcs_qc_min_inlier_rate", 0.50))
    if n_match > 0 and np.isfinite(min_inlier_rate) and min_inlier_rate > 0:
        inlier_rate = _num("inlier_rate")
        if (not np.isfinite(inlier_rate)) or (inlier_rate < min_inlier_rate):
            reasons.append("low_inlier")

    max_edge_ratio = float(getattr(params.P, "wcs_qc_max_edge_ratio", 0.0))
    edge_ratio = _num("edge_resid_ratio")
    if np.isfinite(max_edge_ratio) and max_edge_ratio > 0:
        if np.isfinite(edge_ratio) and edge_ratio > max_edge_ratio:
            reasons.append("edge_resid")

    max_center_off = float(getattr(params.P, "wcs_qc_max_center_offset_arcsec", 0.0))
    center_off = _num("center_offset_arcsec")
    if np.isfinite(max_center_off) and max_center_off > 0:
        if (not np.isfinite(center_off)) or (center_off > max_center_off):
            reasons.append("center_offset")

    return len(reasons) == 0, reasons


def _check_astap_available(params) -> tuple[bool, str]:
    exe = str(getattr(params.P, "astap_exe", "astap_cli.exe") or "astap_cli.exe").strip()
    ok, detail = _local_executable_available(exe)
    if ok:
        db = str(getattr(params.P, "astap_database", "") or "").strip()
        suffix = f"; selected ASTAP DB={db}" if db else ""
        return True, f"ASTAP executable {detail}{suffix}"
    return (
        False,
        f"ASTAP executable {detail}. Install ASTAP and set Step 5 > ASTAP Parameters > ASTAP CLI Path.",
    )


def _classify_astap_failure(rc: int, stdout: str | None, stderr: str | None) -> str:
    text = f"{stdout or ''}\n{stderr or ''}".lower()
    if rc == -996:
        return "astap_not_found"
    if rc == -999 or "timeout" in text:
        return "astap_timeout"
    if rc == -997 or _is_explicit_stopped_message(stderr):
        return "astap_stopped"
    if any(token in text for token in (
        "star database",
        "database not found",
        "no database",
        "could not find database",
        "error reading database",
        "d50",
        "d80",
    )):
        return "astap_database_missing_or_mismatch"
    if any(token in text for token in (
        "no solution",
        "not solved",
        "solution not found",
        "solving failed",
        "not enough stars",
    )):
        return "astap_no_solution"
    if any(token in text for token in (
        "fov",
        "field of view",
        "pixel scale",
        "search radius",
    )):
        return "astap_fov_or_scale"
    return f"astap_rc_{rc}"


def _check_astnet_available(params) -> tuple[bool, str]:
    cmd_text = str(getattr(params.P, "astnet_local_command", "solve-field") or "solve-field").strip()
    cmd_base = _split_command(cmd_text)
    if not cmd_base:
        return False, "astrometry.net command is empty; set Step 5 > Astrometry.net Parameters > solve-field Command."

    use_wsl = bool(getattr(params.P, "astnet_local_use_wsl", True))
    if use_wsl:
        wsl_ok, wsl_detail = _local_executable_available("wsl")
        if not wsl_ok:
            return False, f"WSL executable {wsl_detail}. Install WSL/Ubuntu or disable Use WSL."

        inner_cmd = cmd_base[1:] if cmd_base[0].lower() == "wsl" else cmd_base
        if not inner_cmd:
            return False, "WSL command has no solve-field executable."
        inner_exe = inner_cmd[0]
        probe = f"command -v {shlex.quote(inner_exe)} >/dev/null 2>&1 || test -x {shlex.quote(inner_exe)}"
        try:
            cp = subprocess.run(
                ["wsl", "bash", "-lc", probe],
                capture_output=True,
                timeout=8.0,
                **_SUBPROCESS_TEXT_KWARGS,
            )
            if cp.returncode == 0:
                return True, f"WSL available ({wsl_detail}); {inner_exe} found inside WSL"
            err = _tail_text((cp.stderr or "") + "\n" + (cp.stdout or ""), limit=500, max_lines=4)
            extra = f" Details: {err}" if err else ""
            return (
                False,
                f"WSL is available, but {inner_exe} was not found inside WSL. "
                f"Install astrometry.net/index files in WSL or update solve-field Command.{extra}",
            )
        except Exception as exc:
            return False, f"Could not probe WSL astrometry.net command: {_exc_brief(exc)}"

    exe = cmd_base[0]
    ok, detail = _local_executable_available(exe)
    if ok:
        return True, f"astrometry.net command {detail}"
    return (
        False,
        f"astrometry.net command {detail}. Install solve-field or update Step 5 > Astrometry.net Parameters.",
    )


def _check_gaia_runtime_available() -> tuple[bool, str]:
    return gaia_runtime_available()


def _astnet_center_candidates(
    header_coord: SkyCoord | None,
    target_coord: SkyCoord | None,
) -> list[tuple[str, SkyCoord | None]]:
    attempts: list[tuple[str, SkyCoord | None]] = []

    def _add(label: str, coord: SkyCoord | None) -> None:
        for _, existing in attempts:
            if coord is None and existing is None:
                return
            sep_deg = _coord_sep_deg(coord, existing)
            if coord is not None and existing is not None and np.isfinite(sep_deg) and sep_deg < 1e-6:
                return
        attempts.append((label, coord))

    if header_coord is not None:
        _add("header", header_coord)

    sep_deg = _coord_sep_deg(header_coord, target_coord)
    if target_coord is not None and (
        header_coord is None or (not np.isfinite(sep_deg)) or sep_deg > 0.25
    ):
        _add("target", target_coord)

    _add("blind", None)
    return attempts


def _wsl_path_exists_probe(wsl_path: str) -> bool:
    try:
        cp = subprocess.run(
            ["wsl", "bash", "-lc", f"test -e {shlex.quote(wsl_path)}"],
            capture_output=True,
            timeout=8.0,
            **_SUBPROCESS_TEXT_KWARGS,
        )
        return cp.returncode == 0
    except Exception:
        return False


def _wsl_ensure_writable_dir_probe(wsl_path: str) -> bool:
    try:
        cp = subprocess.run(
            ["wsl", "bash", "-lc", f"mkdir -p {shlex.quote(wsl_path)} && test -w {shlex.quote(wsl_path)}"],
            capture_output=True,
            timeout=8.0,
            **_SUBPROCESS_TEXT_KWARGS,
        )
        return cp.returncode == 0
    except Exception:
        return False


def _preferred_astnet_solution_path(new_path: Path, wcs_path: Path) -> Path | None:
    # Prefer the lightweight header-only WCS artifact over the full solved FITS copy.
    if wcs_path.exists():
        return wcs_path
    if new_path.exists():
        return new_path
    return None


def _astnet_solution_artifacts_ready(new_path: Path, solved_path: Path, wcs_path: Path) -> bool:
    # A header-only .wcs sidecar is sufficient; otherwise require the .new + .solved pair.
    if wcs_path.exists():
        return True
    return new_path.exists() and solved_path.exists()


def _cleanup_redundant_astnet_new(new_path: Path, wcs_path: Path) -> None:
    if new_path.exists() and wcs_path.exists():
        try:
            new_path.unlink()
        except Exception:
            pass


def _solution_header_shape(solution_hdr: fits.Header, fallback_hdr: fits.Header | None = None) -> tuple[float, float]:
    nx = int(solution_hdr.get("NAXIS1", 0) or 0)
    ny = int(solution_hdr.get("NAXIS2", 0) or 0)
    if (nx <= 0 or ny <= 0) and fallback_hdr is not None:
        nx = int(fallback_hdr.get("NAXIS1", 0) or 0)
        ny = int(fallback_hdr.get("NAXIS2", 0) or 0)
    return float(nx), float(ny)


class WcsWorker(QThread):
    """Worker thread for ASTAP WCS solving"""
    progress = pyqtSignal(int, int, str)
    file_done = pyqtSignal(str, dict)
    worker_status = pyqtSignal(int, str, str, int)
    log_message = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str, str)

    def __init__(self, file_list, params, data_dir, result_dir, cache_dir,
                 use_cropped=False, target_coord=None):
        super().__init__()
        self.file_list = file_list
        self.params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(result_dir)
        self.cache_dir = Path(cache_dir)
        self.use_cropped = use_cropped
        self.target_coord = target_coord
        self._stop_requested = False

    def stop(self):
        self._stop_requested = True

    def _resolve_exe(self, exe: str) -> str:
        p = Path(exe)
        if p.exists():
            return str(p)
        return exe

    def _inject_wcs_into_header(self, hdr: fits.Header, wcs_dict: dict):
        for k, v in wcs_dict.items():
            try:
                hdr[k] = v
            except Exception:
                pass

    def _copy_wcs_keywords(self, src_hdr: fits.Header, dst_hdr: fits.Header):
        prefixes = (
            "CRVAL", "CRPIX", "CTYPE", "CUNIT", "CDELT",
            "CD1_", "CD2_", "PC1_", "PC2_", "CROTA",
            "PV", "LONPOLE", "LATPOLE", "RADESYS", "EQUINOX", "WCSAXES",
        )
        for key in src_hdr.keys():
            if key.startswith(prefixes):
                try:
                    dst_hdr[key] = src_hdr[key]
                except Exception:
                    pass

    def _win_to_wsl_path(self, path: Path) -> str:
        try:
            wp = PureWindowsPath(str(path))
            if wp.drive:
                drive = wp.drive.rstrip(":").lower()
                parts = "/".join(wp.parts[1:])
                return f"/mnt/{drive}/{parts}"
        except Exception:
            pass
        return str(path).replace("\\", "/")

    def _wsl_path_exists(self, wsl_path: str) -> bool:
        try:
            cp = subprocess.run(
                ["wsl", "bash", "-lc", f"test -e {shlex.quote(wsl_path)}"],
                capture_output=True,
                timeout=8.0,
                **_SUBPROCESS_TEXT_KWARGS,
            )
            return cp.returncode == 0
        except Exception:
            return False

    def _wsl_ensure_writable_dir(self, wsl_path: str) -> bool:
        try:
            cp = subprocess.run(
                ["wsl", "bash", "-lc", f"mkdir -p {shlex.quote(wsl_path)} && test -w {shlex.quote(wsl_path)}"],
                capture_output=True,
                timeout=8.0,
                **_SUBPROCESS_TEXT_KWARGS,
            )
            return cp.returncode == 0
        except Exception:
            return False

    def _run_solve_field(
        self,
        fits_path: Path,
        center_coord: SkyCoord | None,
        scale_low: float,
        scale_high: float,
        radius_deg: float,
        downsample: int,
        timeout_s: float,
        outdir: Path,
        use_wsl: bool,
        stage_in_outdir: bool = True,
        use_cache: bool = False,
        max_objs: int | None = None,
        cpulimit_s: float | None = None,
    ):
        outdir.mkdir(parents=True, exist_ok=True)
        stem = fits_path.stem
        new_path = outdir / f"{stem}.new"
        solved_path = outdir / f"{stem}.solved"
        wcs_path_out = outdir / f"{stem}.wcs"
        sig_path = outdir / f"{stem}.input.json"
        run_outdir = outdir
        run_new_path = new_path
        run_solved_path = solved_path
        run_wcs_path_out = wcs_path_out
        stage_root = None
        stage_note = ""
        source_sig = build_file_signature(fits_path, use_cropped=bool(self.use_cropped))
        cached_solution_path = _preferred_astnet_solution_path(new_path, wcs_path_out)
        cache_ready = _astnet_solution_artifacts_ready(new_path, solved_path, wcs_path_out)
        if use_cache and cache_ready and cached_solution_path is not None:
            cache_ok = False
            try:
                if sig_path.exists():
                    saved_sig = json.loads(sig_path.read_text(encoding="utf-8"))
                    cache_ok = file_signature_matches_relaxed(saved_sig, source_sig)
            except Exception:
                cache_ok = False
            if cache_ok:
                return True, 0.0, "cache_hit", "", [], cached_solution_path
        for p in outdir.glob(f"{stem}.*"):
            try:
                p.unlink()
            except Exception:
                pass
        staged_path = fits_path
        if stage_in_outdir and (not use_wsl):
            try:
                staged_path = outdir / fits_path.name
                if staged_path != fits_path:
                    shutil.copy2(fits_path, staged_path)
            except Exception:
                staged_path = fits_path
        cmd_str = str(getattr(self.params.P, "astnet_local_command", "solve-field"))
        cmd_base = shlex.split(cmd_str) if cmd_str.strip() else ["solve-field"]
        if use_wsl and cmd_base and cmd_base[0].lower() != "wsl":
            cmd = ["wsl"] + cmd_base
        else:
            cmd = cmd_base

        outdir_arg = self._win_to_wsl_path(run_outdir) if use_wsl else str(run_outdir)
        fits_arg = self._win_to_wsl_path(staged_path) if use_wsl else str(staged_path)
        if not staged_path.exists():
            return False, 0.0, "", f"input_missing:{staged_path}", cmd, None

        if use_wsl:
            input_visible = _wsl_path_exists_probe(fits_arg)
            outdir_writable = _wsl_ensure_writable_dir_probe(outdir_arg)
            if not input_visible or not outdir_writable:
                try:
                    stage_root = (
                        Path(tempfile.gettempdir())
                        / "apex_astnet_wsl"
                        / f"{stem}_{int(time.time() * 1000)}_{threading.get_ident()}"
                    )
                    stage_root.mkdir(parents=True, exist_ok=True)
                    run_outdir = stage_root
                    run_new_path = run_outdir / f"{stem}.new"
                    run_solved_path = run_outdir / f"{stem}.solved"
                    run_wcs_path_out = run_outdir / f"{stem}.wcs"
                    staged_path = run_outdir / fits_path.name
                    shutil.copy2(fits_path, staged_path)
                    outdir_arg = self._win_to_wsl_path(run_outdir)
                    fits_arg = self._win_to_wsl_path(staged_path)
                    stage_note = (
                        f"wsl_stage_fallback:input_visible={input_visible},"
                        f"outdir_writable={outdir_writable},"
                        f"src={fits_path},stage={run_outdir}"
                    )
                except Exception as e:
                    return False, 0.0, "", f"wsl_stage_failed:{e}", cmd, None

        cmd += [
            "--dir", outdir_arg,
            "--scale-units", "arcsecperpix",
            "--scale-low", f"{scale_low:.5f}",
            "--scale-high", f"{scale_high:.5f}",
            "--downsample", str(int(downsample)),
            "--no-verify",
            "--no-plots",
            "--overwrite",
            fits_arg,
        ]
        if max_objs is not None and int(max_objs) > 0:
            cmd += ["--objs", str(int(max_objs))]
        if cpulimit_s is not None and float(cpulimit_s) > 0:
            cmd += ["--cpulimit", f"{float(cpulimit_s):.1f}"]

        if center_coord is not None:
            cmd += [
                "--ra", f"{center_coord.ra.deg:.6f}",
                "--dec", f"{center_coord.dec.deg:.6f}",
                "--radius", f"{radius_deg:.3f}",
            ]

        try:
            start = time.time()
            cp = subprocess.run(
                cmd,
                capture_output=True,
                timeout=timeout_s,
                **_SUBPROCESS_TEXT_KWARGS,
            )
            dt = time.time() - start
            # Some solve-field builds return non-zero even when solution artifacts are created.
            # Prefer artifact existence over process return code.
            solution_path = _preferred_astnet_solution_path(run_new_path, run_wcs_path_out)
            ok = bool(
                solution_path is not None
                and _astnet_solution_artifacts_ready(run_new_path, run_solved_path, run_wcs_path_out)
            )
            if ok and run_outdir != outdir:
                try:
                    if run_wcs_path_out.exists():
                        shutil.copy2(run_wcs_path_out, wcs_path_out)
                    elif run_new_path.exists():
                        shutil.copy2(run_new_path, new_path)
                    if run_solved_path.exists():
                        shutil.copy2(run_solved_path, solved_path)
                except Exception as e:
                    ok = False
                    cp_stderr = cp.stderr or ""
                    cp.stderr = _tail_text(f"stage_copy_failed:{e} | {cp_stderr}", limit=2000, max_lines=12)
            if ok:
                try:
                    sig_path.write_text(json.dumps(source_sig, indent=2), encoding="utf-8")
                except Exception:
                    pass
                _cleanup_redundant_astnet_new(new_path, wcs_path_out)
            if staged_path != fits_path:
                try:
                    staged_path.unlink()
                except Exception:
                    pass
            solution_path = _preferred_astnet_solution_path(new_path, wcs_path_out)
            return ok, dt, cp.stdout, cp.stderr, cmd, solution_path
        except subprocess.TimeoutExpired as e:
            if staged_path != fits_path:
                try:
                    staged_path.unlink()
                except Exception:
                    pass
            out_s = e.stdout or ""
            err_s = e.stderr or ""
            err_msg = "timeout"
            err_tail = _tail_text(err_s, limit=1000, max_lines=10)
            if err_tail:
                err_msg = f"timeout | {err_tail}"
            return False, timeout_s, out_s, err_msg, cmd, None
        except Exception as e:
            if staged_path != fits_path:
                try:
                    staged_path.unlink()
                except Exception:
                    pass
            return False, 0.0, "", str(e), cmd, None

    def _try_ingest_wcs(self, fits_path: Path, hdr: fits.Header) -> bool:
        for wp in astap_wcs_candidates(fits_path):
            if not wp.exists():
                continue
            wcsd = parse_astap_wcs_file(wp)
            if not wcsd:
                continue
            self._inject_wcs_into_header(hdr, wcsd)
            try:
                w = WCS(hdr, relax=True)
                if w.has_celestial:
                    return True
            except Exception:
                pass
        # Fallback: ASTAP이 직접 FITS에 WCS를 쓴 경우 (헤더 재로드)
        try:
            with fits.open(fits_path, memmap=False) as hdul_check:
                hdr_check = hdul_check[0].header
                w_check = WCS(hdr_check, relax=True)
                if w_check.has_celestial:
                    # WCS 키워드 복사
                    for key in hdr_check.keys():
                        if key.startswith(('CRVAL', 'CRPIX', 'CDELT', 'CD1_', 'CD2_',
                                          'CTYPE', 'CUNIT', 'CROTA', 'PC1_', 'PC2_')):
                            hdr[key] = hdr_check[key]
                    return True
        except Exception:
            pass
        return False

    def _pixscale_from_wcs(self, w: WCS) -> float:
        try:
            sc = proj_plane_pixel_scales(w.celestial) * 3600.0
            return float(np.mean(sc))
        except Exception:
            return float("nan")

    def _wcs_rotation_deg(self, w: WCS) -> float:
        """Extract rotation angle from WCS CD matrix (degrees, E of N)."""
        try:
            if not w.has_celestial:
                return float("nan")
            # Get CD matrix or compute from PC+CDELT
            if hasattr(w.wcs, 'cd') and w.wcs.cd is not None:
                cd = w.wcs.cd
            elif hasattr(w.wcs, 'pc') and w.wcs.pc is not None:
                pc = w.wcs.pc
                cdelt = w.wcs.cdelt
                cd = pc * cdelt[:, np.newaxis]
            else:
                return float("nan")
            # Rotation from CD matrix: theta = atan2(-CD1_2, CD2_2)
            rot_rad = np.arctan2(-cd[0, 1], cd[1, 1])
            return float(np.degrees(rot_rad))
        except Exception:
            return float("nan")

    def _wcs_center_coords(self, w: WCS, nx: int, ny: int) -> tuple:
        """Get center RA/Dec from WCS."""
        try:
            if not w.has_celestial:
                return (float("nan"), float("nan"))
            cx, cy = nx / 2.0, ny / 2.0
            sky = w.pixel_to_world(cx, cy)
            return (float(sky.ra.deg), float(sky.dec.deg))
        except Exception:
            return (float("nan"), float("nan"))

    def _wcs_sip_order(self, hdr) -> int:
        """Get SIP distortion polynomial order (0 if none)."""
        try:
            # Check for SIP keywords
            a_order = hdr.get("A_ORDER", 0)
            b_order = hdr.get("B_ORDER", 0)
            return max(int(a_order), int(b_order))
        except Exception:
            return 0

    def _resolve_source_fits_path(self, fname: str):
        if self.use_cropped:
            cand = step2_cropped_dir(self.result_dir) / fname
            if cand.exists():
                return cand
        try:
            orig = self.params.get_file_path(fname)
            orig = Path(orig)
            if orig.exists():
                return orig
        except Exception:
            pass
        return None

    def _compatible_detect_signature(self, fname: str):
        src = self._resolve_source_fits_path(fname)
        if src is None or not src.exists():
            return None
        return build_file_signature(src, use_cropped=bool(self.use_cropped))

    def _detect_meta_matches(self, payload: dict, sig_now: dict, meta_path: Path) -> bool:
        return detection_cache_signature_matches(
            payload,
            sig_now,
            min_schema=2,
            allow_mtime_drift=True,
        )

    def _schema1_detect_cache_allowed(self, marker_path: Path) -> bool:
        try:
            marker_mtime = int(marker_path.stat().st_mtime_ns)
        except Exception:
            return False
        if self.use_cropped:
            rect_path = crop_rect_path(self.result_dir)
            if rect_path.exists():
                try:
                    rect_mtime = int(rect_path.stat().st_mtime_ns)
                    if marker_mtime < rect_mtime:
                        return False
                except Exception:
                    return False
        return True

    def _refresh_detect_cache_signature(self, fits_path: Path, fname: str) -> None:
        """After writing WCS headers into a FITS file, update size+mtime in detect JSONs."""
        try:
            st = fits_path.stat()
            new_size = int(st.st_size)
            new_mtime = int(st.st_mtime_ns)
        except Exception:
            return
        candidates = [
            self.cache_dir / f"detect_{fname}.json",
            step4_dir(self.result_dir) / f"detect_{fname}.json",
        ]
        for p in candidates:
            if not p.exists():
                continue
            try:
                payload = json.loads(p.read_text(encoding="utf-8"))
                payload["source_size"] = new_size
                payload["source_mtime_ns"] = new_mtime
                p.write_text(json.dumps(payload, indent=2), encoding="utf-8")
            except Exception:
                pass

    def _load_fwhm_for_frame(self, fname: str):
        sig_now = self._compatible_detect_signature(fname)
        if sig_now is None:
            return float(getattr(self.params.P, "fwhm_seed_px", 6.0)), np.nan
        candidates = [
            self.cache_dir / f"detect_{fname}.json",
            step4_dir(self.result_dir) / f"detect_{fname}.json",
        ]
        candidates = [p for p in candidates if p.exists()]
        candidates.sort(key=lambda p: p.stat().st_mtime_ns, reverse=True)
        for meta_json in candidates:
            try:
                meta = json.loads(meta_json.read_text(encoding="utf-8"))
            except Exception:
                continue
            if not self._detect_meta_matches(meta, sig_now, meta_json):
                continue
            try:
                fpx = float(meta.get("fwhm_med_rad_px", meta.get("fwhm_med_px", np.nan)))
                farc = float(meta.get("fwhm_med_rad_arcsec", meta.get("fwhm_med_arc", np.nan)))
                return fpx, farc
            except Exception:
                continue
        # Backward compatibility: schema<2 detect cache can be used
        # when it is newer than the current crop selection marker.
        for meta_json in candidates:
            try:
                meta = json.loads(meta_json.read_text(encoding="utf-8"))
            except Exception:
                continue
            try:
                schema = int(meta.get("cache_schema", 0) or 0)
            except Exception:
                schema = 0
            if schema >= 2:
                continue
            if not self._schema1_detect_cache_allowed(meta_json):
                continue
            try:
                fpx = float(meta.get("fwhm_med_rad_px", meta.get("fwhm_med_px", np.nan)))
                farc = float(meta.get("fwhm_med_rad_arcsec", meta.get("fwhm_med_arc", np.nan)))
                return fpx, farc
            except Exception:
                continue
        return float(getattr(self.params.P, "fwhm_seed_px", 6.0)), np.nan

    def _load_detect_xy(self, fname: str) -> np.ndarray:
        sig_now = self._compatible_detect_signature(fname)
        if sig_now is None:
            return np.zeros((0, 2), float)
        candidates = [
            self.cache_dir / f"detect_{fname}.csv",
            step4_dir(self.result_dir) / f"detect_{fname}.csv",
        ]
        for path in candidates:
            if not path.exists():
                continue
            meta_path = path.with_suffix(".json")
            if not meta_path.exists():
                continue
            try:
                payload = json.loads(meta_path.read_text(encoding="utf-8"))
            except Exception:
                continue
            if not self._detect_meta_matches(payload, sig_now, meta_path):
                continue
            try:
                df = pd.read_csv(path)
            except Exception:
                continue
            if not {"x", "y"} <= set(df.columns):
                continue
            xy = df[["x", "y"]].to_numpy(float)
            xy = xy[np.isfinite(xy).all(axis=1)]
            return xy
        # Backward compatibility for schema<2 detect cache.
        for csv_path in candidates:
            if not csv_path.exists():
                continue
            meta_path = csv_path.with_suffix(".json")
            schema = 0
            if meta_path.exists():
                try:
                    payload = json.loads(meta_path.read_text(encoding="utf-8"))
                    schema = int(payload.get("cache_schema", 0) or 0)
                except Exception:
                    continue
            if schema >= 2:
                continue
            marker_path = meta_path if meta_path.exists() else csv_path
            if not self._schema1_detect_cache_allowed(marker_path):
                continue
            try:
                df = pd.read_csv(csv_path)
            except Exception:
                continue
            if not {"x", "y"} <= set(df.columns):
                continue
            xy = df[["x", "y"]].to_numpy(float)
            xy = xy[np.isfinite(xy).all(axis=1)]
            return xy
        return np.zeros((0, 2), float)

    def _empty_wcs_qc_metrics(self, n_detect: int = 0) -> dict:
        return {
            "n_detect": int(max(0, n_detect)),
            "n_catalog_in_fov": 0,
            "n_match": 0,
            "n_inlier": 0,
            "match_rate": np.nan,
            "match_rate_cat": np.nan,
            "match_rate_eff": np.nan,
            "match_radius_arcsec": np.nan,
            "match_radius_px": np.nan,
            "dx_med_px": np.nan,
            "dy_med_px": np.nan,
            "resid_med_px": np.nan,
            "resid_mad_px": np.nan,
            "resid_peak_px": np.nan,
            "resid_p99_px": np.nan,
            "rms_px": np.nan,
            "inlier_rate": np.nan,
            "resid_vs_radius_slope": np.nan,
            "edge_resid_ratio": np.nan,
            "center_offset_arcsec": np.nan,
            "pix_scale_input_arcsec": np.nan,
            "pix_scale_fit_arcsec": np.nan,
            "scale_delta_pct": np.nan,
            "gaia_available": False,
        }

    def _compute_wcs_qc_metrics(
        self,
        *,
        w: WCS | None,
        det_xy: np.ndarray,
        nx: int,
        ny: int,
        gaia_ra_deg: np.ndarray,
        gaia_dec_deg: np.ndarray,
        pix_input_arcsec: float,
        pix_fit_arcsec: float,
        center_coord: SkyCoord | None,
    ) -> dict:
        out = self._empty_wcs_qc_metrics(n_detect=len(det_xy))
        out["gaia_available"] = bool(len(gaia_ra_deg) > 0 and len(gaia_dec_deg) > 0)
        if np.isfinite(pix_input_arcsec):
            out["pix_scale_input_arcsec"] = float(pix_input_arcsec)
        if np.isfinite(pix_fit_arcsec):
            out["pix_scale_fit_arcsec"] = float(pix_fit_arcsec)
        if np.isfinite(pix_input_arcsec) and pix_input_arcsec > 0 and np.isfinite(pix_fit_arcsec):
            out["scale_delta_pct"] = float((pix_fit_arcsec - pix_input_arcsec) / pix_input_arcsec * 100.0)

        if w is not None and w.has_celestial and center_coord is not None:
            try:
                c_ra, c_dec = self._wcs_center_coords(w, nx, ny)
                if np.isfinite(c_ra) and np.isfinite(c_dec):
                    c_sky = SkyCoord(c_ra * u.deg, c_dec * u.deg, frame="icrs")
                    out["center_offset_arcsec"] = float(c_sky.separation(center_coord).arcsec)
            except Exception:
                pass

        if w is None or (not w.has_celestial):
            return out
        if len(det_xy) == 0:
            return out
        if gaia_ra_deg.size == 0 or gaia_dec_deg.size == 0:
            return out

        try:
            xg, yg = w.celestial.all_world2pix(gaia_ra_deg, gaia_dec_deg, 0)
            xg = np.asarray(xg, float)
            yg = np.asarray(yg, float)
        except Exception:
            return out

        ok_g = (
            np.isfinite(xg)
            & np.isfinite(yg)
            & (xg >= 0.0)
            & (xg < float(nx))
            & (yg >= 0.0)
            & (yg < float(ny))
        )
        if not np.any(ok_g):
            return out

        gaia_xy = np.column_stack((xg[ok_g], yg[ok_g]))
        out["n_catalog_in_fov"] = int(len(gaia_xy))

        pix_use = pix_fit_arcsec if np.isfinite(pix_fit_arcsec) and pix_fit_arcsec > 0 else pix_input_arcsec
        match_r_arcsec = float(getattr(self.params.P, "wcs_qc_match_radius_arcsec", 2.0))
        if not np.isfinite(match_r_arcsec) or match_r_arcsec <= 0:
            match_r_arcsec = 2.0
        if np.isfinite(pix_use) and pix_use > 0:
            match_r_px = float(match_r_arcsec / pix_use)
        else:
            match_r_px = float(getattr(self.params.P, "wcs_qc_match_radius_px", 2.5))
        match_r_px = float(np.clip(match_r_px, 1.0, 25.0))
        out["match_radius_arcsec"] = float(match_r_arcsec)
        out["match_radius_px"] = float(match_r_px)

        tree = KDTree(gaia_xy)
        d, j = tree.query(det_xy, k=1)
        d = np.asarray(d, float)
        j = np.asarray(j, int)
        ok = np.isfinite(d) & (d <= match_r_px) & (j >= 0) & (j < len(gaia_xy))
        if not np.any(ok):
            return out

        det_candidates = np.where(ok)[0]
        order = np.argsort(d[det_candidates])
        used_gaia = set()
        keep_det = []
        keep_gaia = []
        for ord_idx in order:
            det_i = int(det_candidates[ord_idx])
            gaia_i = int(j[det_i])
            if gaia_i in used_gaia:
                continue
            used_gaia.add(gaia_i)
            keep_det.append(det_i)
            keep_gaia.append(gaia_i)
        if not keep_det:
            return out

        det_keep = np.asarray(keep_det, dtype=int)
        gaia_keep = np.asarray(keep_gaia, dtype=int)
        dx = det_xy[det_keep, 0] - gaia_xy[gaia_keep, 0]
        dy = det_xy[det_keep, 1] - gaia_xy[gaia_keep, 1]
        r = np.hypot(dx, dy)
        finite_r = np.isfinite(r)
        if not np.any(finite_r):
            return out
        if not np.all(finite_r):
            dx = dx[finite_r]
            dy = dy[finite_r]
            r = r[finite_r]
            gaia_keep = gaia_keep[finite_r]

        n_match = int(len(r))
        out["n_match"] = n_match
        out["match_rate"] = float(n_match / max(int(len(det_xy)), 1))
        out["match_rate_cat"] = float(n_match / max(int(len(gaia_xy)), 1))
        out["match_rate_eff"] = float(max(out["match_rate"], out["match_rate_cat"]))
        if n_match == 0:
            return out

        dx_med = float(np.nanmedian(dx)) if len(dx) else np.nan
        dy_med = float(np.nanmedian(dy)) if len(dy) else np.nan
        out["dx_med_px"] = dx_med
        out["dy_med_px"] = dy_med

        resid_med = float(np.nanmedian(r))
        resid_mad = float(MAD_TO_SIGMA * np.nanmedian(np.abs(r - resid_med)))
        resid_p99 = float(np.nanpercentile(r, 99))
        out["resid_med_px"] = resid_med
        out["resid_mad_px"] = resid_mad
        out["resid_p99_px"] = resid_p99
        out["resid_peak_px"] = resid_p99

        clip_sigma = float(getattr(self.params.P, "wcs_qc_clip_sigma", 3.0))
        if not np.isfinite(clip_sigma) or clip_sigma <= 0:
            clip_sigma = 3.0
        if np.isfinite(resid_mad) and resid_mad > 0:
            inlier = np.abs(r - resid_med) <= clip_sigma * resid_mad
        else:
            resid_std = float(np.nanstd(r))
            if np.isfinite(resid_std) and resid_std > 0:
                inlier = np.abs(r - float(np.nanmean(r))) <= clip_sigma * resid_std
            else:
                inlier = np.ones(len(r), dtype=bool)
        n_inlier = int(np.sum(inlier))
        r_in = r[inlier] if n_inlier > 0 else r
        out["n_inlier"] = n_inlier
        out["inlier_rate"] = float(n_inlier / max(n_match, 1))
        out["rms_px"] = float(np.sqrt(np.nanmean(r_in ** 2))) if len(r_in) else np.nan

        if len(det_keep) >= 8:
            cx = float(nx) / 2.0
            cy = float(ny) / 2.0
            rr = np.hypot(gaia_xy[gaia_keep, 0] - cx, gaia_xy[gaia_keep, 1] - cy)
            max_rr = max(float(np.hypot(max(cx, 1.0), max(cy, 1.0))), 1.0)
            rho = rr / max_rr
            if np.isfinite(np.nanstd(rho)) and float(np.nanstd(rho)) > 1e-6:
                try:
                    out["resid_vs_radius_slope"] = float(np.polyfit(rho, r, 1)[0])
                except Exception:
                    out["resid_vs_radius_slope"] = np.nan
            core = r[rho <= 0.4]
            edge = r[rho >= 0.8]
            if len(core) >= 3 and len(edge) >= 3:
                core_med = float(np.nanmedian(core))
                if np.isfinite(core_med) and core_med > 1e-9:
                    out["edge_resid_ratio"] = float(np.nanmedian(edge) / core_med)

        return out

    def _evaluate_wcs_qc_pass(self, metrics: dict, *, wcs_ok: bool) -> tuple[bool, list[str]]:
        def _num(key: str) -> float:
            try:
                return float(metrics.get(key, np.nan))
            except Exception:
                return np.nan

        reasons: list[str] = []
        gaia_available = bool(metrics.get("gaia_available", True))

        require_wcs_ok = bool(getattr(self.params.P, "wcs_qc_require_wcs_ok", True))
        if require_wcs_ok and not wcs_ok:
            reasons.append("wcs_fail")

        n_detect = int(metrics.get("n_detect", 0) or 0)
        n_match = int(metrics.get("n_match", 0) or 0)
        if n_detect <= 0:
            reasons.append("no_detect_data")

        if not gaia_available:
            reasons.append("gaia_unavailable")
        else:
            min_match_n = int(getattr(self.params.P, "wcs_qc_min_match_n", 20))
            if min_match_n > 0 and n_match < min_match_n:
                reasons.append("low_match_n")

            min_match_rate = float(getattr(self.params.P, "wcs_qc_min_match_rate", 0.20))
            mrate_det = _num("match_rate")
            mrate_cat = _num("match_rate_cat")
            mrate_eff = _num("match_rate_eff")
            if not np.isfinite(mrate_eff):
                if np.isfinite(mrate_det) or np.isfinite(mrate_cat):
                    mrate_eff = float(np.nanmax([v for v in [mrate_det, mrate_cat] if np.isfinite(v)] or [np.nan]))
            if np.isfinite(min_match_rate) and min_match_rate > 0:
                if (not np.isfinite(mrate_eff)) or (mrate_eff < min_match_rate):
                    reasons.append("low_match_rate")

        max_rms_px = float(getattr(self.params.P, "wcs_qc_max_rms_px", 2.5))
        if n_match > 0 and np.isfinite(max_rms_px) and max_rms_px > 0:
            rms_px = _num("rms_px")
            if (not np.isfinite(rms_px)) or (rms_px > max_rms_px):
                reasons.append("high_rms")

        max_p99_px = float(getattr(self.params.P, "wcs_qc_max_p99_px", 5.0))
        if n_match > 0 and np.isfinite(max_p99_px) and max_p99_px > 0:
            p99_px = _num("resid_p99_px")
            if (not np.isfinite(p99_px)) or (p99_px > max_p99_px):
                reasons.append("high_p99")

        min_inlier_rate = float(getattr(self.params.P, "wcs_qc_min_inlier_rate", 0.50))
        if n_match > 0 and np.isfinite(min_inlier_rate) and min_inlier_rate > 0:
            inlier_rate = _num("inlier_rate")
            if (not np.isfinite(inlier_rate)) or (inlier_rate < min_inlier_rate):
                reasons.append("low_inlier")

        max_edge_ratio = float(getattr(self.params.P, "wcs_qc_max_edge_ratio", 0.0))
        edge_ratio = _num("edge_resid_ratio")
        if np.isfinite(max_edge_ratio) and max_edge_ratio > 0:
            if np.isfinite(edge_ratio) and edge_ratio > max_edge_ratio:
                reasons.append("edge_resid")

        max_center_off = float(getattr(self.params.P, "wcs_qc_max_center_offset_arcsec", 0.0))
        center_off = _num("center_offset_arcsec")
        if np.isfinite(max_center_off) and max_center_off > 0:
            if (not np.isfinite(center_off)) or (center_off > max_center_off):
                reasons.append("center_offset")

        return len(reasons) == 0, reasons

    def _load_or_query_gaia(self, center: SkyCoord, radius_deg: float):
        service = GaiaCatalogService(
            self.params,
            self.result_dir,
            log_fn=self.log_message.emit,
            stop_fn=lambda: self._stop_requested,
        )
        return service.load_or_query(center, radius_deg)

    def _refine_crpix_by_match(self, w: WCS, hdr: fits.Header, det_xy: np.ndarray,
                               gaia_df: pd.DataFrame, fwhm_px: float, max_match: int):
        if w is None or (not w.has_celestial):
            return False, "no_wcs", np.nan, np.nan, 0
        if det_xy.size == 0:
            return False, "no_det", np.nan, np.nan, 0
        if gaia_df is None or len(gaia_df) == 0:
            return False, "gaia_unavailable", np.nan, np.nan, 0

        try:
            ra = gaia_df["ra"].to_numpy(float)
            dec = gaia_df["dec"].to_numpy(float)
        except Exception:
            return False, "gaia_cols_missing", np.nan, np.nan, 0

        nx = int(hdr.get("NAXIS1", 0))
        ny = int(hdr.get("NAXIS2", 0))
        if nx <= 0 or ny <= 0:
            return False, "bad_shape", np.nan, np.nan, 0

        try:
            xg, yg = w.celestial.world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
            xg = np.asarray(xg, float)
            yg = np.asarray(yg, float)
            okb = np.isfinite(xg) & np.isfinite(yg) & (xg >= 0) & (xg < nx) & (yg >= 0) & (yg < ny)
            if okb.sum() == 0:
                return False, "gaia_outside", np.nan, np.nan, 0
            gaia_xy = np.vstack([xg[okb], yg[okb]]).T
        except Exception as e:
            return False, f"world2pix_fail:{e}", np.nan, np.nan, 0

        r_match = max(3.0, float(fwhm_px) * float(getattr(self.params.P, "wcs_refine_match_r_fwhm", 1.6)))
        tree = KDTree(gaia_xy)
        d, j = tree.query(det_xy, k=1)
        m = np.isfinite(d) & (d <= r_match)
        min_match = int(getattr(self.params.P, "wcs_refine_min_match", 50))
        if m.sum() < min_match:
            return False, f"match_too_small:{m.sum()}", np.nan, np.nan, int(m.sum())

        det_m = det_xy[m]
        gaia_m = gaia_xy[j[m]]

        if det_m.shape[0] > int(max_match):
            order = np.argsort(d[m])[:int(max_match)]
            det_m, gaia_m = det_m[order], gaia_m[order]

        dx = det_m[:, 0] - gaia_m[:, 0]
        dy = det_m[:, 1] - gaia_m[:, 1]
        dx_med = float(np.median(dx))
        dy_med = float(np.median(dy))

        if "CRPIX1" in hdr and "CRPIX2" in hdr:
            hdr["CRPIX1"] = float(hdr["CRPIX1"]) + dx_med
            hdr["CRPIX2"] = float(hdr["CRPIX2"]) + dy_med
        else:
            return False, "no_crpix", np.nan, np.nan, det_m.shape[0]

        w2 = WCS(hdr, relax=True)
        pix_fit = self._pixscale_from_wcs(w2)
        if not np.isfinite(pix_fit):
            pix_fit = float(getattr(self.params.P, "pixel_scale_arcsec", np.nan))

        resid_arc = np.hypot(dx - dx_med, dy - dy_med) * float(pix_fit)
        resid_med = float(np.median(resid_arc)) if resid_arc.size else np.nan
        resid_max = float(np.max(resid_arc)) if resid_arc.size else np.nan
        return True, f"m1={det_m.shape[0]}", resid_med, resid_max, int(det_m.shape[0])

    def _run_astap(self, fits_path: Path, fov_deg: float, radius_deg: float, timeout_s: float):
        exe = self._resolve_exe(str(getattr(self.params.P, "astap_exe", "astap_cli.exe")))
        db = str(getattr(self.params.P, "astap_database", "") or "").strip()
        z = int(getattr(self.params.P, "astap_downsample_z", 2))
        s = int(getattr(self.params.P, "astap_max_stars_s", 500))
        cmd = [
            exe,
            "-f", str(fits_path),
            "-fov", f"{fov_deg:.6f}",
            "-r", f"{radius_deg:.3f}",
            "-z", str(z),
            "-s", str(s),
        ]
        if db:
            cmd += ["-D", db]
        try:
            start = time.time()
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                **_SUBPROCESS_TEXT_KWARGS,
            )
            stdout_s = ""
            stderr_s = ""
            while True:
                if self._stop_requested:
                    try:
                        proc.terminate()
                        proc.wait(timeout=2.0)
                    except Exception:
                        try:
                            proc.kill()
                        except Exception:
                            pass
                    try:
                        out_s, err_s = proc.communicate(timeout=0.5)
                        stdout_s = out_s or ""
                        stderr_s = err_s or ""
                    except Exception:
                        pass
                    err_msg = "stopped"
                    err_tail = _tail_text(stderr_s, limit=1000, max_lines=10)
                    if err_tail:
                        err_msg = f"stopped | {err_tail}"
                    return False, -997, time.time() - start, stdout_s, err_msg, cmd

                elapsed = time.time() - start
                if elapsed >= timeout_s:
                    try:
                        proc.terminate()
                        proc.wait(timeout=2.0)
                    except Exception:
                        try:
                            proc.kill()
                        except Exception:
                            pass
                    try:
                        out_s, err_s = proc.communicate(timeout=0.5)
                        stdout_s = out_s or ""
                        stderr_s = err_s or ""
                    except Exception:
                        pass
                    err_msg = "timeout"
                    err_tail = _tail_text(stderr_s, limit=1000, max_lines=10)
                    if err_tail:
                        err_msg = f"timeout | {err_tail}"
                    return False, -999, timeout_s, stdout_s, err_msg, cmd

                try:
                    out_s, err_s = proc.communicate(timeout=0.2)
                    stdout_s = out_s or ""
                    stderr_s = err_s or ""
                    break
                except subprocess.TimeoutExpired:
                    continue

            dt = time.time() - start
            rc = int(proc.returncode if proc.returncode is not None else -998)
            ok = (rc == 0)
            return ok, rc, dt, stdout_s, stderr_s, cmd
        except Exception as e:
            return False, -998, 0.0, "", str(e), cmd

    def run(self):
        try:
            warnings.filterwarnings("ignore", message="Keyword name*HIERARCH*")
            warnings.filterwarnings("ignore", message="The WCS transformation has more axes*")

            results = []
            files = list(self.file_list)
            total = len(files)
            pix_arc = float(getattr(self.params.P, "pixel_scale_arcsec", np.nan))
            if not np.isfinite(pix_arc) or pix_arc <= 0:
                raise RuntimeError("pixel_scale_arcsec is not set; run instrument setup first.")

            if not files:
                raise RuntimeError("No files to process (QC filter removed all files).")

            astap_timeout = float(getattr(self.params.P, "astap_timeout_s", 120.0))
            astap_radius = float(getattr(self.params.P, "astap_search_radius_deg", 8.0))
            astap_db = str(getattr(self.params.P, "astap_database", "D80") or "").strip()
            astap_fov_fudge = float(getattr(self.params.P, "astap_fov_fudge", 1.0))
            astnet_local_enable = bool(getattr(self.params.P, "astnet_local_enable", False))
            astnet_use_wsl = bool(getattr(self.params.P, "astnet_local_use_wsl", True))
            astnet_timeout_s = float(getattr(self.params.P, "astnet_local_timeout_s", 300.0))
            astnet_downsample = int(getattr(self.params.P, "astnet_local_downsample", 2))
            astnet_scale_low = float(getattr(self.params.P, "astnet_local_scale_low", 0.0))
            astnet_scale_high = float(getattr(self.params.P, "astnet_local_scale_high", 0.0))
            astnet_radius_deg = float(getattr(self.params.P, "astnet_local_radius_deg", 8.0))
            astnet_keep_outputs = bool(getattr(self.params.P, "astnet_local_keep_outputs", True))
            astnet_use_cache = bool(getattr(self.params.P, "astnet_local_use_cache", True))
            astnet_max_objs = int(getattr(self.params.P, "astnet_local_max_objs", 2000))
            astnet_cpulimit_s = float(getattr(self.params.P, "astnet_local_cpulimit_s", 30.0))
            astnet_blind_retry = bool(getattr(self.params.P, "astnet_blind_retry_on_fail", True))
            astnet_blind_cpulimit_s = float(getattr(self.params.P, "astnet_blind_cpulimit_s", 120.0))
            astap_available, astap_preflight = _check_astap_available(self.params)
            astnet_available = False
            astnet_preflight = (
                "local astrometry.net fallback deferred until ASTAP fails"
                if astnet_local_enable
                else "local astrometry.net fallback disabled"
            )
            astnet_preflight_checked = False
            astnet_preflight_lock = threading.Lock()
            if astnet_local_enable and not astap_available:
                astnet_available, astnet_preflight = _check_astnet_available(self.params)
                astnet_preflight_checked = True
            meta_dir = self.cache_dir / "wcs_solve"
            meta_dir.mkdir(parents=True, exist_ok=True)
            log_path = self.cache_dir / "wcs_solve.log"

            def L(msg, *, gui: bool = False):
                ts = time.strftime("%Y-%m-%d %H:%M:%S")
                line = f"{ts} {msg}"
                try:
                    with open(log_path, "a", encoding="utf-8") as fh:
                        fh.write(line + "\n")
                except Exception:
                    pass
                if gui:
                    self.log_message.emit(msg)

            def log_cmd_failure(tag, fname, reason, cmd=None, stdout=None, stderr=None):
                # Surface the one-line reason in the GUI panel too, otherwise a
                # silent ASTAP failure looks like "nothing happened, then astnet
                # results appeared". Verbose cmd / stdout / stderr tails stay
                # file-only to keep the panel readable.
                L(f"{fname}: {tag} fail reason={reason}", gui=True)
                if cmd:
                    L(f"{fname}: {tag} cmd={' '.join(str(c) for c in cmd)}")
                out_tail = _tail_text(stdout, limit=1600, max_lines=12)
                err_tail = _tail_text(stderr, limit=1600, max_lines=12)
                if out_tail:
                    L(f"{fname}: {tag} stdout_tail={out_tail}")
                if err_tail:
                    L(f"{fname}: {tag} stderr_tail={err_tail}")

            def ensure_astnet_available() -> bool:
                nonlocal astnet_available, astnet_preflight, astnet_preflight_checked
                if not astnet_local_enable:
                    return False
                if astnet_preflight_checked:
                    return bool(astnet_available)
                with astnet_preflight_lock:
                    if not astnet_preflight_checked:
                        astnet_available, astnet_preflight = _check_astnet_available(self.params)
                        astnet_preflight_checked = True
                        L(f"[WCS][Fallback] astrometry.net: {astnet_preflight}", gui=True)
                return bool(astnet_available)

            L("=" * 60)
            L(f"[WCS] start files={len(files)} use_cropped={self.use_cropped} cache_dir={self.cache_dir}")
            L(f"[WCS][Preflight] ASTAP: {astap_preflight}")
            L(f"[WCS][Preflight] astrometry.net: {astnet_preflight}", gui=(not astap_available))
            L(f"[WCS] astap_timeout_s={astap_timeout} astap_radius_deg={astap_radius} astap_db={astap_db or 'default'} astap_fov_fudge={astap_fov_fudge}")
            L(
                f"[WCS] astnet_local_enable={astnet_local_enable} use_wsl={astnet_use_wsl} "
                f"timeout_s={astnet_timeout_s} downsample={astnet_downsample} "
                f"scale=[{astnet_scale_low:.5f},{astnet_scale_high:.5f}] radius_deg={astnet_radius_deg} "
                f"blind_retry={astnet_blind_retry} blind_cpulimit_s={astnet_blind_cpulimit_s}"
            )

            # Determine Gaia center.  PRIORITY:
            #   1. self.target_coord (project_state, set by Step 1)
            #   2. targets_simbad.tsv (Step 1's SIMBAD resolution sidecar)
            #   3. FITS header OBJCTRA/OBJCTDEC (mount log — least trusted)
            #
            # Mount logs occasionally carry the previous slew target's
            # coords (we have seen M3 stamped into M13 frames).  SIMBAD-
            # resolved values always describe the science target, so they
            # are the authoritative fallback.
            header_coord = None
            try:
                sample = files[0]
                if self.use_cropped:
                    fits_path = step2_cropped_dir(self.result_dir) / sample
                else:
                    fits_path = self.params.get_file_path(sample)
                with fits.open(fits_path, memmap=False) as hdul:
                    hdr = hdul[0].header
                    ra0 = hdr.get("OBJCTRA", None)
                    dec0 = hdr.get("OBJCTDEC", None)
                    if ra0 is not None and dec0 is not None:
                        header_coord = SkyCoord(str(ra0), str(dec0), unit=(u.hourangle, u.deg))
            except Exception:
                header_coord = None

            simbad_coord = _load_simbad_target_coord(self.result_dir)

            # Choose center, layered fallback with cross-checks against header
            center_coord = None
            center_source = "missing"

            if self.target_coord is not None:
                center_coord = self.target_coord
                center_source = "project_state"
            elif simbad_coord is not None:
                center_coord = simbad_coord
                center_source = "simbad_tsv"
            elif header_coord is not None:
                center_coord = header_coord
                center_source = "fits_header"

            # If a trusted target (project_state/simbad) and header disagree
            # by more than 5 deg, warn loudly — keep the trusted target.
            trusted_coord = self.target_coord or simbad_coord
            if trusted_coord is not None and header_coord is not None:
                sep_deg = float(header_coord.separation(trusted_coord).deg)
                if sep_deg > 5.0:
                    L(
                        f"[WCS] WARNING: FITS header (OBJCTRA/DEC) differs by "
                        f"{sep_deg:.2f}deg from {center_source}. Using "
                        f"{center_source} center; header is likely a stale "
                        f"mount log.",
                        gui=True,
                    )

            if center_coord is None:
                raise RuntimeError(
                    "Target coordinate not set. Resolve target in Step 1 "
                    "(writes targets_simbad.tsv) or set parameters.toml "
                    "target.ra_deg / target.dec_deg."
                )
            L(f"[WCS] Gaia center source = {center_source} ({center_coord.ra.deg:.6f}, {center_coord.dec.deg:.6f})")

            # Gaia query/cache
            if self._stop_requested:
                self.finished.emit({"stopped": True, "total": 0, "ok": 0, "wcs_qc_pass": 0})
                return
            gaia_fudge = float(getattr(self.params.P, "gaia_radius_fudge", 1.35))
            sample = files[0]
            if self.use_cropped:
                sample_path = step2_cropped_dir(self.result_dir) / sample
            else:
                sample_path = self.params.get_file_path(sample)
            with fits.open(sample_path, memmap=False, ignore_missing_simple=True) as hdul:
                data0 = hdul[0].data
                if data0 is None:
                    raise RuntimeError("First frame data is None")
                ny0, nx0 = data0.shape
            fov_w = (nx0 * pix_arc) / 3600.0
            fov_h = (ny0 * pix_arc) / 3600.0
            diag_deg = float(np.hypot(fov_w, fov_h))
            gaia_r = float(0.5 * diag_deg * gaia_fudge)
            L(
                f"[WCS] sample={sample} shape={nx0}x{ny0} pix_arcsec={pix_arc:.5f} "
                f"fov_w={fov_w:.5f}deg fov_h={fov_h:.5f}deg diag={diag_deg:.5f}deg "
                f"gaia_r={gaia_r:.5f}deg"
            )
            self.progress.emit(0, total, "Querying Gaia catalog...")
            self.worker_status.emit(0, "Gaia catalog", "Querying", 5)
            L(
                f"[Gaia] querying catalog before WCS refinement: "
                f"center=({center_coord.ra.deg:.6f},{center_coord.dec.deg:.6f}) "
                f"r={gaia_r:.4f}deg",
                gui=True,
            )
            gaia_df, gaia_src = self._load_or_query_gaia(center_coord, gaia_r)
            if self._stop_requested:
                self.finished.emit({"stopped": True, "total": 0, "ok": 0, "wcs_qc_pass": 0})
                return
            L(f"[Gaia] center=({center_coord.ra.deg:.6f},{center_coord.dec.deg:.6f}) r={gaia_r:.4f}deg source={gaia_src} N={len(gaia_df)}")
            if len(gaia_df) == 0:
                self.worker_status.emit(0, "Gaia catalog", "Unavailable", 100)
                L(
                    f"[Gaia][WARN] catalog unavailable ({gaia_src}); "
                    "WCS refine/residual match stats will be skipped.",
                    gui=True,
                )
            else:
                self.worker_status.emit(0, "Gaia catalog", "Loaded", 100)
                L(f"[Gaia] catalog loaded: {len(gaia_df)} sources ({gaia_src})", gui=True)
            gaia_ra_vals = np.array([], dtype=float)
            gaia_dec_vals = np.array([], dtype=float)
            if isinstance(gaia_df, pd.DataFrame) and (not gaia_df.empty) and {"ra", "dec"} <= set(gaia_df.columns):
                gaia_ra_vals = pd.to_numeric(gaia_df["ra"], errors="coerce").to_numpy(float)
                gaia_dec_vals = pd.to_numeric(gaia_df["dec"], errors="coerce").to_numpy(float)
                ok_rd = np.isfinite(gaia_ra_vals) & np.isfinite(gaia_dec_vals)
                gaia_ra_vals = gaia_ra_vals[ok_rd]
                gaia_dec_vals = gaia_dec_vals[ok_rd]

            def solve_one(filename):
                if self._stop_requested:
                    return filename, None

                if self.use_cropped:
                    fits_path = step2_cropped_dir(self.result_dir) / filename
                else:
                    fits_path = self.params.get_file_path(filename)

                status = "fail"
                fail_reason = ""
                pix_fit = np.nan
                wcs_ok = False
                refine_note = ""
                resid_med = np.nan
                resid_max = np.nan
                match_n = 0
                rc = -998
                dt_astap = 0.0
                ok_astap = False
                astap_cmd = []
                astap_stdout = ""
                astap_stderr = ""
                astnet_reason = ""
                det_xy = self._load_detect_xy(filename)
                qc_metrics = self._empty_wcs_qc_metrics(n_detect=len(det_xy))
                qc_metrics["pix_scale_input_arcsec"] = float(pix_arc) if np.isfinite(pix_arc) else np.nan

                with fits.open(fits_path, memmap=False, ignore_missing_simple=True) as hdul:
                    hdr = hdul[0].header
                    data = hdul[0].data
                    if data is None:
                        return filename, {
                            "fname": filename,
                            "ok": False,
                            "status": "data_none",
                            "wcs_ok": False,
                            "pix_fit": None,
                            "elapsed": 0.0,
                            "fail_reason": "data_none",
                            "refine": "",
                            "solver": "astap",
                            "wcs_qc_pass": False,
                            "wcs_qc_reason": "data_none,wcs_fail",
                            **qc_metrics,
                        }
                    ny, nx = data.shape

                if self._stop_requested:
                    return filename, None
                fov_w_deg = (nx * pix_arc) / 3600.0 * astap_fov_fudge
                fov_h_deg = (ny * pix_arc) / 3600.0 * astap_fov_fudge
                fov_deg = float(max(fov_w_deg, fov_h_deg))

                if astap_available:
                    ok_astap, rc, dt_astap, astap_stdout, astap_stderr, astap_cmd = self._run_astap(
                        fits_path, fov_deg=fov_deg, radius_deg=astap_radius, timeout_s=astap_timeout
                    )
                else:
                    ok_astap = False
                    rc = -996
                    dt_astap = 0.0
                    astap_stdout = ""
                    astap_stderr = astap_preflight
                    astap_cmd = [str(getattr(self.params.P, "astap_exe", "astap_cli.exe") or "astap_cli.exe")]
                if self._stop_requested:
                    return filename, None
                astap_cmd_str = " ".join(str(c) for c in astap_cmd)
                if not ok_astap:
                    fail_reason = _classify_astap_failure(rc, astap_stdout, astap_stderr)
                    log_cmd_failure(
                        "ASTAP",
                        filename,
                        f"{fail_reason}, dt={dt_astap:.1f}s",
                        cmd=astap_cmd,
                        stdout=astap_stdout,
                        stderr=astap_stderr,
                    )
                    if not astnet_local_enable or not ensure_astnet_available():
                        if astnet_local_enable and not astnet_available:
                            fail_reason = f"{fail_reason}; astnet_not_found"
                        qc_metrics["pix_scale_fit_arcsec"] = np.nan
                        status_text = "astap_not_found" if rc == -996 else f"astap_fail rc={rc}"
                        return filename, {
                            "fname": filename,
                            "ok": False,
                            "status": status_text,
                            "wcs_ok": False,
                            "pix_fit": pix_fit,
                            "elapsed": float(dt_astap),
                            "fail_reason": fail_reason,
                            "refine": refine_note,
                            "solver": "astap",
                            "astap_cmd": astap_cmd_str,
                            "astap_stdout": _tail_text(astap_stdout, limit=2000, max_lines=12),
                            "astap_stderr": _tail_text(astap_stderr, limit=2000, max_lines=12),
                            "wcs_qc_pass": False,
                            "wcs_qc_reason": "astap_fail,wcs_fail",
                            **qc_metrics,
                        }

                astnet_ok = False
                astnet_dt = np.nan
                astnet_stdout = ""
                astnet_stderr = ""
                astnet_cmd = []
                astnet_cmd_str = ""
                astnet_solution_path = None
                astnet_hint_source = ""
                astnet_center_hint = ""
                astnet_attempts = ""
                solver = "astap"
                used_elapsed = float(dt_astap)

                # WCS를 FITS 파일에 저장 (writeto 사용 - Windows 호환성)
                with fits.open(fits_path, memmap=False, ignore_missing_simple=True) as hdul:
                    hdr = hdul[0].header
                    data = hdul[0].data
                    w_final = None
                    wcs_ok = False
                    if ok_astap:
                        try:
                            w0 = WCS(hdr, relax=True)
                            wcs_ok = w0.has_celestial
                        except Exception:
                            wcs_ok = False

                        if not wcs_ok:
                            wcs_ok = self._try_ingest_wcs(fits_path, hdr)

                    if (not ok_astap or not wcs_ok) and astnet_local_enable and ensure_astnet_available():
                        scale_low = astnet_scale_low
                        scale_high = astnet_scale_high
                        if scale_low <= 0 or scale_high <= 0:
                            scale_low = float(pix_arc) * 0.85
                            scale_high = float(pix_arc) * 1.15

                        outdir = meta_dir / "astnet_local"
                        frame_header_coord = None
                        try:
                            ra0 = hdr.get("OBJCTRA", None)
                            dec0 = hdr.get("OBJCTDEC", None)
                            if ra0 is not None and dec0 is not None:
                                frame_header_coord = SkyCoord(str(ra0), str(dec0), unit=(u.hourangle, u.deg))
                        except Exception:
                            frame_header_coord = None

                        target_attempt_coord = self.target_coord if self.target_coord is not None else center_coord
                        attempt_specs = _astnet_center_candidates(
                            header_coord=frame_header_coord,
                            target_coord=target_attempt_coord,
                        )
                        has_hint_attempt = any(coord is not None for _, coord in attempt_specs)
                        if not astnet_blind_retry and has_hint_attempt:
                            attempt_specs = [(label, coord) for label, coord in attempt_specs if coord is not None]
                        if not attempt_specs:
                            attempt_specs = [("blind", None)]
                        attempt_notes: list[str] = []
                        for attempt_idx, (astnet_hint_source, attempt_coord) in enumerate(attempt_specs, start=1):
                            astnet_center_hint = _format_coord_hint(attempt_coord)
                            is_blind_retry = attempt_coord is None and has_hint_attempt
                            attempt_cpulimit_s = astnet_blind_cpulimit_s if is_blind_retry else astnet_cpulimit_s
                            L(
                                f"[ASTNET_WSL] {filename} attempt={attempt_idx}/{len(attempt_specs)} "
                                f"hint={astnet_hint_source} center={astnet_center_hint} "
                                f"cpulimit={attempt_cpulimit_s:.0f}s"
                            )
                            astnet_ok, astnet_dt, astnet_stdout, astnet_stderr, astnet_cmd, astnet_solution_path = (
                                self._run_solve_field(
                                    fits_path,
                                    center_coord=attempt_coord,
                                    scale_low=scale_low,
                                    scale_high=scale_high,
                                    radius_deg=astnet_radius_deg,
                                    downsample=astnet_downsample,
                                    timeout_s=astnet_timeout_s,
                                    outdir=outdir,
                                    use_wsl=astnet_use_wsl,
                                    use_cache=astnet_use_cache,
                                    max_objs=astnet_max_objs,
                                    cpulimit_s=attempt_cpulimit_s,
                                )
                            )
                            attempt_notes.append(
                                f"{astnet_hint_source}@{astnet_center_hint}:{'ok' if astnet_ok else 'fail'}:{float(astnet_dt):.1f}s:cpu={attempt_cpulimit_s:.0f}"
                            )
                            if astnet_ok:
                                break
                        astnet_attempts = " | ".join(attempt_notes)
                        astnet_cmd_str = " ".join(str(c) for c in astnet_cmd)
                        if astnet_ok:
                            L(
                                f"[ASTNET_WSL] {filename} success dt={astnet_dt:.1f}s "
                                f"hint={astnet_hint_source} center={astnet_center_hint}"
                            )
                            if astnet_stderr:
                                L(f"{filename}: ASTNET_WSL stderr_tail={_tail_text(astnet_stderr, limit=600, max_lines=6)}")
                        else:
                            if "timeout" in str(astnet_stderr).lower():
                                astnet_reason = "astnet_timeout"
                            elif self._stop_requested or _is_explicit_stopped_message(astnet_stderr):
                                astnet_reason = "astnet_stopped"
                            elif (
                                "wsl_path_unavailable" in str(astnet_stderr).lower()
                                or "wsl_not_found" in str(astnet_stderr).lower()
                                or "wsl_input_missing" in str(astnet_stderr).lower()
                                or "input_missing:" in str(astnet_stderr).lower()
                            ):
                                astnet_reason = "astnet_wsl_path_unavailable"
                            elif (
                                "cannot open `" in str(astnet_stderr).lower()
                                and "no such file or directory" in str(astnet_stderr).lower()
                            ):
                                astnet_reason = "astnet_missing_file_or_index"
                            elif "cache_hit" in str(astnet_stdout):
                                astnet_reason = "astnet_cache_miss_solved_marker_missing"
                            else:
                                astnet_reason = "astnet_fail_no_solution"
                            if not fail_reason:
                                fail_reason = astnet_reason
                            log_cmd_failure(
                                f"ASTNET_WSL[{astnet_hint_source or 'unknown'}]",
                                filename,
                                f"{astnet_reason}, dt={astnet_dt:.1f}s",
                                cmd=astnet_cmd,
                                stdout=astnet_stdout,
                                stderr=astnet_stderr,
                            )
                            if astnet_attempts:
                                L(f"{filename}: ASTNET_WSL attempts={astnet_attempts}")

                        if astnet_ok and astnet_solution_path is not None and astnet_solution_path.exists():
                            try:
                                with fits.open(astnet_solution_path, memmap=False) as hdul_new:
                                    new_hdr = hdul_new[0].header
                                    self._copy_wcs_keywords(new_hdr, hdr)
                                w_new = WCS(hdr, relax=True)
                                wcs_ok = w_new.has_celestial
                            except Exception:
                                wcs_ok = False

                            if wcs_ok:
                                hdr["WCSSRC"] = ("ASTNET_WSL", "WCS source")
                                solver = "astnet_wsl"
                                used_elapsed = float(astnet_dt)

                        if not astnet_keep_outputs:
                            for p in outdir.glob(f"{fits_path.stem}.*"):
                                try:
                                    p.unlink()
                                except Exception:
                                    pass

                    if wcs_ok:
                        w = WCS(hdr, relax=True)
                        pix_fit = self._pixscale_from_wcs(w)

                        refine_enable = bool(getattr(self.params.P, "wcs_refine_enable", True))
                        if not refine_enable:
                            refine_note = "disabled"
                        elif gaia_df is None or len(gaia_df) == 0:
                            refine_note = f"gaia_unavailable:{gaia_src}"
                        elif len(det_xy) == 0:
                            refine_note = "no_det"
                        else:
                            fwhm_px, _ = self._load_fwhm_for_frame(filename)
                            ok_ref, note, rmed, rmax, nmatch = self._refine_crpix_by_match(
                                w, hdr, det_xy, gaia_df,
                                fwhm_px=float(fwhm_px),
                                max_match=int(getattr(self.params.P, "wcs_refine_max_match", 600))
                            )
                            refine_note = note
                            match_n = int(nmatch)
                            if ok_ref:
                                w2 = WCS(hdr, relax=True)
                                pix_fit = self._pixscale_from_wcs(w2)
                                resid_med = rmed
                                resid_max = rmax

                        # Use final WCS (refined if available)
                        w_final = WCS(hdr, relax=True)
                        qc_metrics = self._compute_wcs_qc_metrics(
                            w=w_final,
                            det_xy=det_xy,
                            nx=int(hdr.get("NAXIS1", nx)),
                            ny=int(hdr.get("NAXIS2", ny)),
                            gaia_ra_deg=gaia_ra_vals,
                            gaia_dec_deg=gaia_dec_vals,
                            pix_input_arcsec=float(pix_arc),
                            pix_fit_arcsec=float(pix_fit) if np.isfinite(pix_fit) else np.nan,
                            center_coord=center_coord,
                        )

                        # Wipe APEX-managed meta from any previous solver before
                        # writing this run's tags (raw WCS keys are overwritten
                        # naturally by the solver).
                        _reset_apex_wcs_meta(hdr)
                        if solver == "astnet_wsl":
                            hdr["WCSSRC"] = ("ASTNET_WSL", "WCS source (apex)")
                        else:
                            hdr["WCSSRC"] = ("ASTAP", "WCS source (apex)")
                        hdr["WCS_OK"] = (True, "WCS solve success")
                        hdr["WCSPIXI"] = (float(pix_arc), "pixscale input (arcsec/pix)")
                        if np.isfinite(pix_fit):
                            hdr["WCSPIXF"] = (float(pix_fit), "pixscale fit (arcsec/pix)")
                        if isinstance(match_n, (int, float)) and int(match_n) > 0:
                            hdr["WCSNST"] = (int(match_n), "WCS matched stars")
                        if refine_note:
                            hdr["WCSREFN"] = (str(refine_note)[:68], "refine note")
                        if np.isfinite(resid_med):
                            hdr["WCSRMD"] = (float(resid_med), "ref resid med (arcsec)")
                        if np.isfinite(resid_max):
                            hdr["WCSRMAX"] = (float(resid_max), "ref resid max (arcsec)")

                        # Additional WCS quality stats
                        wcs_rot_deg = self._wcs_rotation_deg(w_final)
                        nx = int(hdr.get("NAXIS1", 0))
                        ny = int(hdr.get("NAXIS2", 0))
                        center_ra, center_dec = self._wcs_center_coords(w_final, nx, ny)
                        sip_order = self._wcs_sip_order(hdr)

                        if np.isfinite(wcs_rot_deg):
                            hdr["WCSROT"] = (float(wcs_rot_deg), "WCS rotation (deg, E of N)")
                        if np.isfinite(center_ra):
                            hdr["WCSCRA"] = (float(center_ra), "WCS center RA (deg)")
                        if np.isfinite(center_dec):
                            hdr["WCSCDEC"] = (float(center_dec), "WCS center Dec (deg)")
                        if sip_order > 0:
                            hdr["WCSSIP"] = (int(sip_order), "SIP distortion order")

                        status = "ok_astnet_wsl" if solver == "astnet_wsl" else "ok"
                    else:
                        hdr["WCS_OK"] = (False, "WCS solve failed")
                        if rc == -996:
                            status = "astap_not_found"
                        elif not ok_astap:
                            status = f"astap_fail rc={rc}"
                        else:
                            status = "wcs_missing"
                        if not fail_reason:
                            if not ok_astap:
                                fail_reason = _classify_astap_failure(rc, astap_stdout, astap_stderr)
                            elif astnet_local_enable and astnet_reason:
                                fail_reason = astnet_reason
                            elif astnet_local_enable:
                                fail_reason = "astnet_ran_but_wcs_missing"
                            else:
                                fail_reason = "wcs_keywords_missing_after_astap"
                        L(
                            f"{filename}: final_wcs_fail status={status} fail_reason={fail_reason} "
                            f"solver={solver} astap_ok={ok_astap} astnet_ok={astnet_ok}"
                        )
                        # Set defaults for failed WCS
                        wcs_rot_deg = np.nan
                        center_ra = np.nan
                        center_dec = np.nan
                        sip_order = 0

                    qc_pass, qc_reasons = self._evaluate_wcs_qc_pass(qc_metrics, wcs_ok=bool(wcs_ok))
                    qc_reason = ",".join(qc_reasons)

                if self._stop_requested:
                    return filename, None
                # writeto로 확실하게 저장 (Windows 호환)
                fits.writeto(fits_path, data, hdr, overwrite=True)
                self._refresh_detect_cache_signature(fits_path, filename)
                if self._stop_requested:
                    return filename, None

                meta = {
                    "fname": filename,
                    "ok": bool(wcs_ok),
                    "wcs_ok": bool(wcs_ok),
                    "status": status,
                    "fail_reason": fail_reason,
                    "pix_fit": float(pix_fit) if np.isfinite(pix_fit) else None,
                    "elapsed": float(used_elapsed),
                    "refine": refine_note,
                    "resid_med": float(resid_med) if np.isfinite(resid_med) else None,
                    "resid_max": float(resid_max) if np.isfinite(resid_max) else None,
                    "match_n": int(match_n),
                    "gaia_source": str(gaia_src),
                    "solver": solver,
                    # New WCS quality stats
                    "wcs_rot_deg": float(wcs_rot_deg) if np.isfinite(wcs_rot_deg) else None,
                    "center_ra_deg": float(center_ra) if np.isfinite(center_ra) else None,
                    "center_dec_deg": float(center_dec) if np.isfinite(center_dec) else None,
                    "sip_order": int(sip_order) if sip_order > 0 else None,
                    # Solver details
                    "astap_ok": bool(ok_astap),
                    "astap_rc": int(rc),
                    "astap_elapsed": float(dt_astap),
                    "astap_cmd": astap_cmd_str,
                    "astap_stdout": _tail_text(astap_stdout, limit=2000, max_lines=12),
                    "astap_stderr": _tail_text(astap_stderr, limit=2000, max_lines=12),
                    "astnet_wsl_ok": bool(astnet_ok),
                    "astnet_wsl_elapsed": float(astnet_dt) if np.isfinite(astnet_dt) else None,
                    "astnet_wsl_fail_reason": astnet_reason,
                    "astnet_wsl_hint_source": astnet_hint_source,
                    "astnet_wsl_center_hint": astnet_center_hint,
                    "astnet_wsl_attempts": astnet_attempts,
                    "astnet_wsl_cmd": astnet_cmd_str,
                    "astnet_wsl_stdout": _tail_text(astnet_stdout, limit=2000, max_lines=12),
                    "astnet_wsl_stderr": _tail_text(astnet_stderr, limit=2000, max_lines=12),
                    "wcs_qc_pass": bool(qc_pass),
                    "wcs_qc_reason": qc_reason,
                }
                meta.update(qc_metrics)
                rms_px_val = qc_metrics.get("rms_px", np.nan)
                try:
                    rms_px_val = float(rms_px_val)
                except Exception:
                    rms_px_val = np.nan
                rms_px_str = f"{rms_px_val:.3f}" if np.isfinite(rms_px_val) else "-"
                L(
                    f"{filename}: {status} pix_fit={pix_fit:.4f} dt={used_elapsed:.1f}s "
                    f"refine={refine_note or '-'} resid_med={resid_med if np.isfinite(resid_med) else '-'} "
                    f"match_n={match_n} wcs_qc={'PASS' if qc_pass else 'FAIL'} "
                    f"n_det={int(qc_metrics.get('n_detect', 0) or 0)} "
                    f"n_match={int(qc_metrics.get('n_match', 0) or 0)} "
                    f"rms_px={rms_px_str} "
                    f"reason={qc_reason or '-'} "
                    f"fail_reason={fail_reason or '-'}"
                )
                if self._stop_requested:
                    return filename, None
                (meta_dir / f"wcs_{filename}.json").write_text(
                    json.dumps(meta, indent=2), encoding="utf-8"
                )
                return filename, meta

            completed = 0
            _default_workers = get_parallel_workers(self.params)
            max_workers = int(getattr(self.params.P, "wcs_max_workers", _default_workers))
            from queue import Queue as _Queue
            _slot_q: _Queue = _Queue()
            for _i in range(max(1, max_workers)):
                _slot_q.put(_i)

            def _slotted_solve(filename: str):
                slot = _slot_q.get()
                try:
                    self.worker_status.emit(slot, filename, "Solving", 5)
                    out = solve_one(filename)
                    self.worker_status.emit(slot, filename, "Done", 100)
                    return out
                finally:
                    _slot_q.put(slot)

            # Explicit ownership allows queued ASTAP solves to be cancelled.
            ex = ThreadPoolExecutor(max_workers=max(1, max_workers))
            try:
                futures = {ex.submit(_slotted_solve, f): f for f in files}
                for fut in as_completed(futures):
                    if self._stop_requested:
                        for f_cancel in futures:
                            f_cancel.cancel()
                        break
                    fname = futures[fut]
                    completed += 1
                    try:
                        fn, res = fut.result()
                        if res is not None:
                            results.append(res)
                            self.file_done.emit(fn, res)
                        else:
                            self.error.emit(fname, "stopped")
                    except Exception as e:
                        L(f"{fname}: worker_exception={_exc_brief(e)}")
                        self.error.emit(fname, str(e))
                    self.progress.emit(completed, len(files), fname)
            finally:
                ex.shutdown(wait=True, cancel_futures=True)

            # Save summary CSV
            try:
                df = pd.DataFrame(results)
                step5_out = step5_wcs_dir(self.result_dir)
                step5_out.mkdir(parents=True, exist_ok=True)
                df.to_csv(step5_out / "wcs_solve_summary.csv", index=False)
                qc_cols = [
                    "fname",
                    "status",
                    "fail_reason",
                    "ok",
                    "wcs_ok",
                    "solver",
                    "elapsed",
                    "n_detect",
                    "n_catalog_in_fov",
                    "n_match",
                    "n_inlier",
                    "match_rate",
                    "match_rate_cat",
                    "match_rate_eff",
                    "match_radius_arcsec",
                    "match_radius_px",
                    "dx_med_px",
                    "dy_med_px",
                    "resid_med_px",
                    "resid_mad_px",
                    "resid_p99_px",
                    "resid_peak_px",
                    "rms_px",
                    "inlier_rate",
                    "resid_vs_radius_slope",
                    "edge_resid_ratio",
                    "pix_scale_input_arcsec",
                    "pix_scale_fit_arcsec",
                    "scale_delta_pct",
                    "wcs_rot_deg",
                    "center_ra_deg",
                    "center_dec_deg",
                    "center_offset_arcsec",
                    "sip_order",
                    "wcs_qc_pass",
                    "wcs_qc_reason",
                ]
                qc_df = df[[c for c in qc_cols if c in df.columns]].copy()
                if "fname" in qc_df.columns and "file" not in qc_df.columns:
                    qc_df = qc_df.rename(columns={"fname": "file"})
                qc_df.to_csv(step5_out / "frame_wcs_qc.csv", index=False)
            except Exception:
                pass

            n_qc_pass = sum(1 for r in results if bool(r.get("wcs_qc_pass", False)))
            n_qc_not_eval = sum(
                1
                for r in results
                if "gaia_unavailable" in str(r.get("wcs_qc_reason", ""))
            )
            summary = {
                "total": len(results),
                "ok": sum(1 for r in results if r.get("ok")),
                "wcs_qc_pass": int(n_qc_pass),
                "wcs_qc_not_evaluated": int(n_qc_not_eval),
                "stopped": bool(self._stop_requested),
            }
            self.finished.emit(summary)
        except Exception as e:
            self.error.emit("WORKER", str(e))
            self.finished.emit({})


class AstrometryNetWorker(QThread):
    """Worker thread for local astrometry.net (solve-field) WCS solving"""
    progress = pyqtSignal(int, int, str)
    file_done = pyqtSignal(str, dict)
    refine_done = pyqtSignal(str, dict)  # 파일명, refine 결과
    log_message = pyqtSignal(str)
    worker_status = pyqtSignal(int, str, str, int)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str, str)

    def __init__(self, file_list, params, data_dir, result_dir, cache_dir,
                 use_cropped=False, target_coord=None):
        super().__init__()
        self.file_list = file_list
        self.params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(result_dir)
        self.cache_dir = Path(cache_dir)
        self.use_cropped = use_cropped
        self.target_coord = target_coord
        self._stop_requested = False

    def stop(self):
        self._stop_requested = True

    def _pixscale_from_wcs(self, w: WCS):
        try:
            sc = proj_plane_pixel_scales(w.celestial) * 3600.0
            return float(np.mean(sc))
        except Exception:
            return float("nan")

    def _load_or_query_gaia(self, center: SkyCoord, radius_deg: float):
        service = GaiaCatalogService(
            self.params,
            self.result_dir,
            log_fn=self.log_message.emit,
            stop_fn=lambda: self._stop_requested,
        )
        return service.load_or_query(center, radius_deg)

    def _load_fwhm_for_frame(self, fname: str):
        meta_json = self.cache_dir / f"detect_{fname}.json"
        if not meta_json.exists():
            fallback = step4_dir(self.result_dir) / f"detect_{fname}.json"
            if fallback.exists():
                meta_json = fallback
        if meta_json.exists():
            try:
                meta = json.loads(meta_json.read_text(encoding="utf-8"))
                fpx = float(
                    meta.get(
                        "fwhm_med_rad_px",
                        meta.get("fwhm_med_px", meta.get("fwhm_px", np.nan)),
                    )
                )
                farc = float(
                    meta.get(
                        "fwhm_med_rad_arcsec",
                        meta.get("fwhm_med_arc", meta.get("fwhm_arcsec", np.nan)),
                    )
                )
                return fpx, farc
            except Exception:
                pass
        return float(getattr(self.params.P, "fwhm_seed_px", 6.0)), np.nan

    def _load_detect_xy(self, fname: str):
        csv_path = self.cache_dir / f"detect_{fname}.csv"
        if not csv_path.exists():
            fallback = step4_dir(self.result_dir) / f"detect_{fname}.csv"
            if fallback.exists():
                csv_path = fallback
            else:
                return np.empty((0, 2)), None
        try:
            df = pd.read_csv(csv_path)
            xy = df[["x", "y"]].values
            return xy, df
        except Exception:
            return np.empty((0, 2)), None

    def _wcs_rotation_deg(self, w: WCS) -> float:
        try:
            if not w.has_celestial:
                return float("nan")
            if hasattr(w.wcs, "cd") and w.wcs.cd is not None:
                cd = w.wcs.cd
            elif hasattr(w.wcs, "pc") and w.wcs.pc is not None:
                pc = w.wcs.pc
                cdelt = w.wcs.cdelt
                cd = pc * cdelt[:, np.newaxis]
            else:
                return float("nan")
            rot_rad = np.arctan2(-cd[0, 1], cd[1, 1])
            return float(np.degrees(rot_rad))
        except Exception:
            return float("nan")

    def _wcs_center_coords(self, w: WCS, nx: int, ny: int) -> tuple:
        try:
            if not w.has_celestial:
                return (float("nan"), float("nan"))
            cx, cy = nx / 2.0, ny / 2.0
            sky = w.pixel_to_world(cx, cy)
            return (float(sky.ra.deg), float(sky.dec.deg))
        except Exception:
            return (float("nan"), float("nan"))

    def _empty_wcs_qc_metrics(self, n_detect: int = 0) -> dict:
        return {
            "n_detect": int(max(0, n_detect)),
            "n_catalog_in_fov": 0,
            "n_match": 0,
            "n_inlier": 0,
            "match_rate": np.nan,
            "match_rate_cat": np.nan,
            "match_rate_eff": np.nan,
            "match_radius_arcsec": np.nan,
            "match_radius_px": np.nan,
            "dx_med_px": np.nan,
            "dy_med_px": np.nan,
            "resid_med_px": np.nan,
            "resid_mad_px": np.nan,
            "resid_peak_px": np.nan,
            "resid_p99_px": np.nan,
            "rms_px": np.nan,
            "inlier_rate": np.nan,
            "resid_vs_radius_slope": np.nan,
            "edge_resid_ratio": np.nan,
            "center_offset_arcsec": np.nan,
            "pix_scale_input_arcsec": np.nan,
            "pix_scale_fit_arcsec": np.nan,
            "scale_delta_pct": np.nan,
            "gaia_available": False,
        }

    def _compute_wcs_qc_metrics(
        self,
        *,
        w: WCS | None,
        det_xy: np.ndarray,
        nx: int,
        ny: int,
        gaia_ra_deg: np.ndarray,
        gaia_dec_deg: np.ndarray,
        pix_input_arcsec: float,
        pix_fit_arcsec: float,
        center_coord: SkyCoord | None,
    ) -> dict:
        out = self._empty_wcs_qc_metrics(n_detect=len(det_xy))
        out["gaia_available"] = bool(len(gaia_ra_deg) > 0 and len(gaia_dec_deg) > 0)
        if np.isfinite(pix_input_arcsec):
            out["pix_scale_input_arcsec"] = float(pix_input_arcsec)
        if np.isfinite(pix_fit_arcsec):
            out["pix_scale_fit_arcsec"] = float(pix_fit_arcsec)
        if np.isfinite(pix_input_arcsec) and pix_input_arcsec > 0 and np.isfinite(pix_fit_arcsec):
            out["scale_delta_pct"] = float((pix_fit_arcsec - pix_input_arcsec) / pix_input_arcsec * 100.0)

        if w is not None and w.has_celestial and center_coord is not None:
            try:
                c_ra, c_dec = self._wcs_center_coords(w, nx, ny)
                if np.isfinite(c_ra) and np.isfinite(c_dec):
                    c_sky = SkyCoord(c_ra * u.deg, c_dec * u.deg, frame="icrs")
                    out["center_offset_arcsec"] = float(c_sky.separation(center_coord).arcsec)
            except Exception:
                pass

        if w is None or (not w.has_celestial):
            return out
        if len(det_xy) == 0:
            return out
        if gaia_ra_deg.size == 0 or gaia_dec_deg.size == 0:
            return out

        try:
            xg, yg = w.celestial.all_world2pix(gaia_ra_deg, gaia_dec_deg, 0)
            xg = np.asarray(xg, float)
            yg = np.asarray(yg, float)
        except Exception:
            return out

        ok_g = (
            np.isfinite(xg)
            & np.isfinite(yg)
            & (xg >= 0.0)
            & (xg < float(nx))
            & (yg >= 0.0)
            & (yg < float(ny))
        )
        if not np.any(ok_g):
            return out

        gaia_xy = np.column_stack((xg[ok_g], yg[ok_g]))
        out["n_catalog_in_fov"] = int(len(gaia_xy))

        pix_use = pix_fit_arcsec if np.isfinite(pix_fit_arcsec) and pix_fit_arcsec > 0 else pix_input_arcsec
        match_r_arcsec = float(getattr(self.params.P, "wcs_qc_match_radius_arcsec", 2.0))
        if not np.isfinite(match_r_arcsec) or match_r_arcsec <= 0:
            match_r_arcsec = 2.0
        if np.isfinite(pix_use) and pix_use > 0:
            match_r_px = float(match_r_arcsec / pix_use)
        else:
            match_r_px = float(getattr(self.params.P, "wcs_qc_match_radius_px", 2.5))
        match_r_px = float(np.clip(match_r_px, 1.0, 25.0))
        out["match_radius_arcsec"] = float(match_r_arcsec)
        out["match_radius_px"] = float(match_r_px)

        tree = KDTree(gaia_xy)
        d, j = tree.query(det_xy, k=1)
        d = np.asarray(d, float)
        j = np.asarray(j, int)
        ok = np.isfinite(d) & (d <= match_r_px) & (j >= 0) & (j < len(gaia_xy))
        if not np.any(ok):
            return out

        det_candidates = np.where(ok)[0]
        order = np.argsort(d[det_candidates])
        used_gaia = set()
        keep_det = []
        keep_gaia = []
        for ord_idx in order:
            det_i = int(det_candidates[ord_idx])
            gaia_i = int(j[det_i])
            if gaia_i in used_gaia:
                continue
            used_gaia.add(gaia_i)
            keep_det.append(det_i)
            keep_gaia.append(gaia_i)
        if not keep_det:
            return out

        det_keep = np.asarray(keep_det, dtype=int)
        gaia_keep = np.asarray(keep_gaia, dtype=int)
        dx = det_xy[det_keep, 0] - gaia_xy[gaia_keep, 0]
        dy = det_xy[det_keep, 1] - gaia_xy[gaia_keep, 1]
        r = np.hypot(dx, dy)
        finite_r = np.isfinite(r)
        if not np.any(finite_r):
            return out
        if not np.all(finite_r):
            dx = dx[finite_r]
            dy = dy[finite_r]
            r = r[finite_r]
            gaia_keep = gaia_keep[finite_r]

        n_match = int(len(r))
        out["n_match"] = n_match
        out["match_rate"] = float(n_match / max(int(len(det_xy)), 1))
        out["match_rate_cat"] = float(n_match / max(int(len(gaia_xy)), 1))
        out["match_rate_eff"] = float(max(out["match_rate"], out["match_rate_cat"]))
        if n_match == 0:
            return out

        out["dx_med_px"] = float(np.nanmedian(dx)) if len(dx) else np.nan
        out["dy_med_px"] = float(np.nanmedian(dy)) if len(dy) else np.nan
        resid_med = float(np.nanmedian(r))
        resid_mad = float(MAD_TO_SIGMA * np.nanmedian(np.abs(r - resid_med)))
        out["resid_med_px"] = resid_med
        out["resid_mad_px"] = resid_mad
        out["resid_p99_px"] = float(np.nanpercentile(r, 99))
        out["resid_peak_px"] = out["resid_p99_px"]

        clip_sigma = float(getattr(self.params.P, "wcs_qc_clip_sigma", 3.0))
        if not np.isfinite(clip_sigma) or clip_sigma <= 0:
            clip_sigma = 3.0
        if np.isfinite(resid_mad) and resid_mad > 0:
            inlier = np.abs(r - resid_med) <= clip_sigma * resid_mad
        else:
            rstd = float(np.nanstd(r))
            inlier = np.abs(r - float(np.nanmean(r))) <= clip_sigma * rstd if np.isfinite(rstd) and rstd > 0 else np.ones(len(r), dtype=bool)
        n_inlier = int(np.sum(inlier))
        r_in = r[inlier] if n_inlier > 0 else r
        out["n_inlier"] = n_inlier
        out["inlier_rate"] = float(n_inlier / max(n_match, 1))
        out["rms_px"] = float(np.sqrt(np.nanmean(r_in ** 2))) if len(r_in) else np.nan

        if len(det_keep) >= 8:
            cx = float(nx) / 2.0
            cy = float(ny) / 2.0
            rr = np.hypot(gaia_xy[gaia_keep, 0] - cx, gaia_xy[gaia_keep, 1] - cy)
            max_rr = max(float(np.hypot(max(cx, 1.0), max(cy, 1.0))), 1.0)
            rho = rr / max_rr
            if np.isfinite(np.nanstd(rho)) and float(np.nanstd(rho)) > 1e-6:
                try:
                    out["resid_vs_radius_slope"] = float(np.polyfit(rho, r, 1)[0])
                except Exception:
                    out["resid_vs_radius_slope"] = np.nan
            core = r[rho <= 0.4]
            edge = r[rho >= 0.8]
            if len(core) >= 3 and len(edge) >= 3:
                core_med = float(np.nanmedian(core))
                if np.isfinite(core_med) and core_med > 1e-9:
                    out["edge_resid_ratio"] = float(np.nanmedian(edge) / core_med)

        return out

    def _evaluate_wcs_qc_pass(self, metrics: dict, *, wcs_ok: bool) -> tuple[bool, list[str]]:
        def _num(key: str) -> float:
            try:
                return float(metrics.get(key, np.nan))
            except Exception:
                return np.nan

        reasons: list[str] = []
        if bool(getattr(self.params.P, "wcs_qc_require_wcs_ok", True)) and not wcs_ok:
            reasons.append("wcs_fail")
        gaia_available = bool(metrics.get("gaia_available", True))
        if not gaia_available:
            reasons.append("gaia_unavailable")
            return False, reasons
        n_match = int(metrics.get("n_match", 0) or 0)
        if int(getattr(self.params.P, "wcs_qc_min_match_n", 20)) > 0 and n_match < int(getattr(self.params.P, "wcs_qc_min_match_n", 20)):
            reasons.append("low_match_n")
        mrate = _num("match_rate")
        min_rate = float(getattr(self.params.P, "wcs_qc_min_match_rate", 0.20))
        if np.isfinite(min_rate) and min_rate > 0 and ((not np.isfinite(mrate)) or (mrate < min_rate)):
            reasons.append("low_match_rate")
        rms_px = _num("rms_px")
        max_rms = float(getattr(self.params.P, "wcs_qc_max_rms_px", 2.5))
        if np.isfinite(max_rms) and max_rms > 0 and ((not np.isfinite(rms_px)) or (rms_px > max_rms)):
            reasons.append("high_rms")
        p99_px = _num("resid_p99_px")
        max_p99 = float(getattr(self.params.P, "wcs_qc_max_p99_px", 5.0))
        if np.isfinite(max_p99) and max_p99 > 0 and ((not np.isfinite(p99_px)) or (p99_px > max_p99)):
            reasons.append("high_p99")
        inlier = _num("inlier_rate")
        min_inlier = float(getattr(self.params.P, "wcs_qc_min_inlier_rate", 0.50))
        if np.isfinite(min_inlier) and min_inlier > 0 and ((not np.isfinite(inlier)) or (inlier < min_inlier)):
            reasons.append("low_inlier")
        edge_ratio = _num("edge_resid_ratio")
        max_edge = float(getattr(self.params.P, "wcs_qc_max_edge_ratio", 0.0))
        if np.isfinite(max_edge) and max_edge > 0 and np.isfinite(edge_ratio) and edge_ratio > max_edge:
            reasons.append("edge_resid")
        center_off = _num("center_offset_arcsec")
        max_center = float(getattr(self.params.P, "wcs_qc_max_center_offset_arcsec", 0.0))
        if np.isfinite(max_center) and max_center > 0 and ((not np.isfinite(center_off)) or (center_off > max_center)):
            reasons.append("center_offset")
        return len(reasons) == 0, reasons

    def _refine_crpix_by_match(self, w: WCS, hdr: fits.Header, det_xy: np.ndarray,
                               gaia_df: pd.DataFrame, fwhm_px: float, max_match: int):
        if w is None or (not w.has_celestial):
            return False, "no_wcs", np.nan, np.nan, 0
        if det_xy.size == 0:
            return False, "no_det", np.nan, np.nan, 0
        if gaia_df is None or len(gaia_df) == 0:
            return False, "gaia_unavailable", np.nan, np.nan, 0

        try:
            ra = gaia_df["ra"].to_numpy(float)
            dec = gaia_df["dec"].to_numpy(float)
        except Exception:
            return False, "gaia_cols_missing", np.nan, np.nan, 0

        nx = int(hdr.get("NAXIS1", 0))
        ny = int(hdr.get("NAXIS2", 0))
        if nx <= 0 or ny <= 0:
            return False, "bad_shape", np.nan, np.nan, 0

        try:
            xg, yg = w.celestial.world_to_pixel(SkyCoord(ra * u.deg, dec * u.deg))
            xg = np.asarray(xg, float)
            yg = np.asarray(yg, float)
            okb = np.isfinite(xg) & np.isfinite(yg) & (xg >= 0) & (xg < nx) & (yg >= 0) & (yg < ny)
            if okb.sum() == 0:
                return False, "gaia_outside", np.nan, np.nan, 0
            gaia_xy = np.vstack([xg[okb], yg[okb]]).T
        except Exception as e:
            return False, f"world2pix_fail:{e}", np.nan, np.nan, 0

        r_match = max(3.0, float(fwhm_px) * float(getattr(self.params.P, "wcs_refine_match_r_fwhm", 1.6)))
        tree = KDTree(gaia_xy)
        d, j = tree.query(det_xy, k=1)
        m = np.isfinite(d) & (d <= r_match)
        min_match = int(getattr(self.params.P, "wcs_refine_min_match", 50))
        if m.sum() < min_match:
            return False, f"match_too_small:{m.sum()}", np.nan, np.nan, int(m.sum())

        det_m = det_xy[m]
        gaia_m = gaia_xy[j[m]]

        if det_m.shape[0] > int(max_match):
            order = np.argsort(d[m])[:int(max_match)]
            det_m, gaia_m = det_m[order], gaia_m[order]

        dx = det_m[:, 0] - gaia_m[:, 0]
        dy = det_m[:, 1] - gaia_m[:, 1]
        dx_med = float(np.median(dx))
        dy_med = float(np.median(dy))

        if "CRPIX1" in hdr and "CRPIX2" in hdr:
            hdr["CRPIX1"] = float(hdr["CRPIX1"]) + dx_med
            hdr["CRPIX2"] = float(hdr["CRPIX2"]) + dy_med
        else:
            return False, "no_crpix", np.nan, np.nan, det_m.shape[0]

        w2 = WCS(hdr, relax=True)
        pix_fit = self._pixscale_from_wcs(w2)
        if not np.isfinite(pix_fit):
            pix_fit = float(getattr(self.params.P, "pixel_scale_arcsec", np.nan))

        resid_arc = np.hypot(dx - dx_med, dy - dy_med) * float(pix_fit)
        resid_med = float(np.median(resid_arc)) if resid_arc.size else np.nan
        resid_max = float(np.max(resid_arc)) if resid_arc.size else np.nan
        return True, f"m1={det_m.shape[0]}", resid_med, resid_max, int(det_m.shape[0])

    def _win_to_wsl_path(self, path: Path) -> str:
        try:
            wp = PureWindowsPath(str(path))
            if wp.drive:
                drive = wp.drive.rstrip(":").lower()
                parts = "/".join(wp.parts[1:])
                return f"/mnt/{drive}/{parts}"
        except Exception:
            pass
        return str(path).replace("\\", "/")

    def _run_solve_field(
        self,
        fits_path: Path,
        center_coord: SkyCoord | None,
        scale_low: float,
        scale_high: float,
        radius_deg: float,
        downsample: int,
        timeout_s: float,
        outdir: Path,
        use_wsl: bool,
        stage_in_outdir: bool = True,
        use_cache: bool = False,
        max_objs: int | None = None,
        cpulimit_s: float | None = None,
    ):
        outdir.mkdir(parents=True, exist_ok=True)
        stem = fits_path.stem
        new_path = outdir / f"{stem}.new"
        solved_path = outdir / f"{stem}.solved"
        wcs_path_out = outdir / f"{stem}.wcs"
        sig_path = outdir / f"{stem}.input.json"
        run_outdir = outdir
        run_new_path = new_path
        run_solved_path = solved_path
        run_wcs_path_out = wcs_path_out
        stage_root = None
        stage_note = ""
        source_sig = build_file_signature(fits_path, use_cropped=bool(self.use_cropped))
        cached_solution_path = _preferred_astnet_solution_path(new_path, wcs_path_out)
        cache_ready = _astnet_solution_artifacts_ready(new_path, solved_path, wcs_path_out)
        if use_cache and cache_ready and cached_solution_path is not None:
            cache_ok = False
            try:
                if sig_path.exists():
                    saved_sig = json.loads(sig_path.read_text(encoding="utf-8"))
                    cache_ok = file_signature_matches_relaxed(saved_sig, source_sig)
            except Exception:
                cache_ok = False
            if cache_ok:
                return True, 0.0, "cache_hit", "", [], cached_solution_path
        for p in outdir.glob(f"{stem}.*"):
            try:
                p.unlink()
            except Exception:
                pass
        staged_path = fits_path
        if stage_in_outdir and (not use_wsl):
            try:
                staged_path = outdir / fits_path.name
                if staged_path != fits_path:
                    shutil.copy2(fits_path, staged_path)
            except Exception:
                staged_path = fits_path
        cmd_str = str(getattr(self.params.P, "astnet_local_command", "solve-field"))
        cmd_base = shlex.split(cmd_str) if cmd_str.strip() else ["solve-field"]
        if use_wsl and cmd_base and cmd_base[0].lower() != "wsl":
            cmd = ["wsl"] + cmd_base
        else:
            cmd = cmd_base

        outdir_arg = self._win_to_wsl_path(run_outdir) if use_wsl else str(run_outdir)
        fits_arg = self._win_to_wsl_path(staged_path) if use_wsl else str(staged_path)
        if not staged_path.exists():
            return False, 0.0, "", f"input_missing:{staged_path}", cmd, None

        if use_wsl:
            input_visible = _wsl_path_exists_probe(fits_arg)
            outdir_writable = _wsl_ensure_writable_dir_probe(outdir_arg)
            if not input_visible or not outdir_writable:
                try:
                    stage_root = (
                        Path(tempfile.gettempdir())
                        / "apex_astnet_wsl"
                        / f"{stem}_{int(time.time() * 1000)}_{threading.get_ident()}"
                    )
                    stage_root.mkdir(parents=True, exist_ok=True)
                    run_outdir = stage_root
                    run_new_path = run_outdir / f"{stem}.new"
                    run_solved_path = run_outdir / f"{stem}.solved"
                    run_wcs_path_out = run_outdir / f"{stem}.wcs"
                    staged_path = run_outdir / fits_path.name
                    shutil.copy2(fits_path, staged_path)
                    outdir_arg = self._win_to_wsl_path(run_outdir)
                    fits_arg = self._win_to_wsl_path(staged_path)
                    stage_note = (
                        f"wsl_stage_fallback:input_visible={input_visible},"
                        f"outdir_writable={outdir_writable},"
                        f"src={fits_path},stage={run_outdir}"
                    )
                except Exception as e:
                    return False, 0.0, "", f"wsl_stage_failed:{e}", cmd, None

        cmd += [
            "--dir", outdir_arg,
            "--scale-units", "arcsecperpix",
            "--scale-low", f"{scale_low:.5f}",
            "--scale-high", f"{scale_high:.5f}",
            "--downsample", str(int(downsample)),
            "--no-verify",
            "--no-plots",
            "--overwrite",
            fits_arg,
        ]
        if max_objs is not None and int(max_objs) > 0:
            cmd += ["--objs", str(int(max_objs))]
        if cpulimit_s is not None and float(cpulimit_s) > 0:
            cmd += ["--cpulimit", f"{float(cpulimit_s):.1f}"]

        if center_coord is not None:
            cmd += [
                "--ra", f"{center_coord.ra.deg:.6f}",
                "--dec", f"{center_coord.dec.deg:.6f}",
                "--radius", f"{radius_deg:.3f}",
            ]

        try:
            start = time.time()
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                **_SUBPROCESS_TEXT_KWARGS,
            )
            stdout_s = ""
            stderr_s = ""
            while True:
                if self._stop_requested:
                    try:
                        proc.terminate()
                        proc.wait(timeout=2.0)
                    except Exception:
                        try:
                            proc.kill()
                        except Exception:
                            pass
                    err_msg = "stopped"
                    try:
                        out_s, err_s = proc.communicate(timeout=2.0)
                        stdout_s = out_s or ""
                        stderr_s = err_s or ""
                    except Exception:
                        pass
                    err_tail = _tail_text(stderr_s, limit=1000, max_lines=10)
                    if err_tail:
                        err_msg = f"stopped | {err_tail}"
                    if staged_path != fits_path:
                        try:
                            staged_path.unlink()
                        except Exception:
                            pass
                    if stage_root is not None:
                        try:
                            shutil.rmtree(stage_root, ignore_errors=True)
                        except Exception:
                            pass
                    return False, time.time() - start, stdout_s, err_msg, cmd, None

                elapsed = time.time() - start
                if elapsed >= timeout_s:
                    try:
                        proc.terminate()
                        proc.wait(timeout=2.0)
                    except Exception:
                        try:
                            proc.kill()
                        except Exception:
                            pass
                    try:
                        out_s, err_s = proc.communicate(timeout=0.5)
                        stdout_s = out_s or ""
                        stderr_s = err_s or ""
                    except Exception:
                        pass
                    err_msg = "timeout"
                    err_tail = _tail_text(stderr_s, limit=1000, max_lines=10)
                    if err_tail:
                        err_msg = f"timeout | {err_tail}"
                    if staged_path != fits_path:
                        try:
                            staged_path.unlink()
                        except Exception:
                            pass
                    if stage_root is not None:
                        try:
                            shutil.rmtree(stage_root, ignore_errors=True)
                        except Exception:
                            pass
                    return False, timeout_s, stdout_s, err_msg, cmd, None

                try:
                    out_s, err_s = proc.communicate(timeout=0.5)
                    stdout_s = out_s or ""
                    stderr_s = err_s or ""
                    break
                except subprocess.TimeoutExpired:
                    continue

            dt = time.time() - start
            rc = int(proc.returncode if proc.returncode is not None else -998)
            # Some solve-field builds return non-zero or omit the .solved marker
            # even when usable WCS artifacts were produced.
            solution_path = _preferred_astnet_solution_path(run_new_path, run_wcs_path_out)
            ok = bool(
                solution_path is not None
                and _astnet_solution_artifacts_ready(run_new_path, run_solved_path, run_wcs_path_out)
            )
            if ok and run_outdir != outdir:
                try:
                    if run_wcs_path_out.exists():
                        shutil.copy2(run_wcs_path_out, wcs_path_out)
                    elif run_new_path.exists():
                        shutil.copy2(run_new_path, new_path)
                    if run_solved_path.exists():
                        shutil.copy2(run_solved_path, solved_path)
                except Exception as e:
                    ok = False
                    stderr_s = _tail_text(f"stage_copy_failed:{e} | {stderr_s}", limit=2000, max_lines=12)
            if ok:
                try:
                    sig_path.write_text(json.dumps(source_sig, indent=2), encoding="utf-8")
                except Exception:
                    pass
                _cleanup_redundant_astnet_new(new_path, wcs_path_out)
            if stage_note:
                stdout_s = f"{stage_note}\n{stdout_s}".strip()
            if staged_path != fits_path:
                try:
                    staged_path.unlink()
                except Exception:
                    pass
            if stage_root is not None:
                try:
                    shutil.rmtree(stage_root, ignore_errors=True)
                except Exception:
                    pass
            solution_path = _preferred_astnet_solution_path(new_path, wcs_path_out)
            return ok, dt, stdout_s, stderr_s, cmd, solution_path
        except Exception as e:
            if staged_path != fits_path:
                try:
                    staged_path.unlink()
                except Exception:
                    pass
            if stage_root is not None:
                try:
                    shutil.rmtree(stage_root, ignore_errors=True)
                except Exception:
                    pass
            return False, 0.0, "", str(e), cmd, None

    def _refresh_detect_cache_signature(self, fits_path: Path, fname: str) -> None:
        """After writing WCS headers into a FITS file, update size+mtime in detect JSONs.

        Step 5 writes WCS keywords into the original FITS (changing file size/mtime).
        Without this, Step 6 treats all detection caches as invalid (signature mismatch).
        """
        try:
            st = fits_path.stat()
            new_size = int(st.st_size)
            new_mtime = int(st.st_mtime_ns)
        except Exception:
            return
        candidates = [
            self.cache_dir / f"detect_{fname}.json",
            step4_dir(self.result_dir) / f"detect_{fname}.json",
        ]
        for p in candidates:
            if not p.exists():
                continue
            try:
                payload = json.loads(p.read_text(encoding="utf-8"))
                payload["source_size"] = new_size
                payload["source_mtime_ns"] = new_mtime
                p.write_text(json.dumps(payload, indent=2), encoding="utf-8")
            except Exception:
                pass

    def _safe_header_update(self, fits_path: Path, new_hdr: fits.Header) -> None:
        for i in range(5):
            try:
                with fits.open(fits_path, mode="update", memmap=False) as hdul:
                    hdr = hdul[0].header
                    # Wipe APEX-managed meta from a previous solver before
                    # stamping our own. Raw WCS keys are overwritten below.
                    _reset_apex_wcs_meta(hdr)
                    for key in new_hdr.keys():
                        if key.startswith(("CRVAL", "CRPIX", "CTYPE", "CUNIT", "CDELT",
                                           "CD1_", "CD2_", "PC1_", "PC2_", "CROTA",
                                           "PV", "LONPOLE", "LATPOLE", "RADESYS", "EQUINOX", "WCSAXES")):
                            try:
                                hdr[key] = new_hdr[key]
                            except Exception:
                                pass
                    hdr["WCS_OK"] = (True, "WCS solve by astrometry.net (local)")
                    hdr["WCSSRC"] = ("ASTNET_WSL", "WCS source (apex)")
                return
            except Exception:
                time.sleep(0.4 * (i + 1))
        raise RuntimeError("failed to update FITS header (locked)")

    def run(self):
        results = []
        total = len(self.file_list)

        # --- 파라미터 로드 ---
        pix_arc = float(getattr(self.params.P, "pixel_scale_arcsec", np.nan))
        if not np.isfinite(pix_arc) or pix_arc <= 0:
            self.error.emit("WORKER", "pixel_scale_arcsec is not set; run instrument setup first.")
            self.finished.emit({})
            return

        use_wsl = bool(getattr(self.params.P, "astnet_local_use_wsl", True))
        timeout_s = float(getattr(self.params.P, "astnet_local_timeout_s", 300.0))

        downsample = int(getattr(self.params.P, "astnet_local_downsample", 2))
        max_objs = int(getattr(self.params.P, "astnet_local_max_objs", 2000))
        scale_low = float(getattr(self.params.P, "astnet_local_scale_low", 0.0))
        scale_high = float(getattr(self.params.P, "astnet_local_scale_high", 0.0))
        radius_deg = float(getattr(self.params.P, "astnet_local_radius_deg", 8.0))
        keep_outputs = bool(getattr(self.params.P, "astnet_local_keep_outputs", True))
        use_cache = bool(getattr(self.params.P, "astnet_local_use_cache", True))
        cpulimit_s = float(getattr(self.params.P, "astnet_local_cpulimit_s", 30.0))
        blind_retry = bool(getattr(self.params.P, "astnet_blind_retry_on_fail", True))
        blind_cpulimit_s = float(getattr(self.params.P, "astnet_blind_cpulimit_s", 120.0))
        max_workers = get_parallel_workers(self.params)

        if scale_low <= 0 or scale_high <= 0:
            scale_low = float(pix_arc) * 0.85
            scale_high = float(pix_arc) * 1.15

        outdir = self.cache_dir / "wcs_solve" / "astnet_local"
        log_path = self.cache_dir / "astnet_solve.log"

        def LOG(msg, emit=False):
            ts = time.strftime("%Y-%m-%d %H:%M:%S")
            line = f"{ts} {msg}"
            try:
                with open(log_path, "a", encoding="utf-8") as fh:
                    fh.write(line + "\n")
            except Exception:
                pass
            if emit:
                self.log_message.emit(msg)

        LOG("=" * 60)
        LOG(
            f"[ASTNET] start files={len(self.file_list)} use_cropped={self.use_cropped} "
            f"workers={max_workers} outdir={outdir}"
        )
        LOG(
            f"[ASTNET] scale=[{scale_low:.5f},{scale_high:.5f}] downsample={downsample} "
            f"max_objs={max_objs} radius_deg={radius_deg} timeout_s={timeout_s} "
            f"cpulimit_s={cpulimit_s} use_wsl={use_wsl} use_cache={use_cache} "
            f"blind_retry={blind_retry} blind_cpulimit_s={blind_cpulimit_s}"
        )

        self.log_message.emit(f"Starting parallel plate solving with {max_workers} workers...")
        self.log_message.emit(f"  Scale: {scale_low:.4f} - {scale_high:.4f} arcsec/px")
        self.log_message.emit(f"  Radius: {radius_deg:.1f} deg")
        self.log_message.emit(f"  Downsample: {downsample}, Max objs: {max_objs}")
        self.log_message.emit(
            f"  Blind retry: {'on' if blind_retry else 'off'}"
            + (f", blind CPU limit: {blind_cpulimit_s:.0f}s" if blind_retry else "")
        )
        self.log_message.emit(f"  Debug log: {log_path}")

        # 내부 함수: 단일 파일 처리 로직 (스레드에서 실행됨)
        def process_single_file(filename):
            if self._stop_requested:
                LOG(f"{filename}: stop_requested_before_start")
                return filename, {"ok": False, "status": "stopped", "fail_reason": "stopped"}

            if self.use_cropped:
                fits_path = step2_cropped_dir(self.result_dir) / filename
            else:
                fits_path = self.params.get_file_path(filename)

            if not fits_path.exists():
                LOG(f"{filename}: file_not_found path={fits_path}")
                return filename, {"ok": False, "status": "file_not_found", "fail_reason": "file_not_found"}

            # 헤더에서 중심 좌표 읽기
            center_coord = None
            try:
                # memmap=False로 읽어서 파일 핸들 즉시 반환 유도
                with fits.open(fits_path, memmap=False) as hdul:
                    hdr0 = hdul[0].header
                    ra0 = hdr0.get("OBJCTRA", None)
                    dec0 = hdr0.get("OBJCTDEC", None)
                    if ra0 is not None and dec0 is not None:
                        center_coord = SkyCoord(str(ra0), str(dec0), unit=(u.hourangle, u.deg))
            except Exception:
                center_coord = None

            attempt_specs = _astnet_center_candidates(header_coord=center_coord, target_coord=self.target_coord)
            has_hint_attempt = any(coord is not None for _, coord in attempt_specs)
            if not blind_retry and has_hint_attempt:
                attempt_specs = [(label, coord) for label, coord in attempt_specs if coord is not None]
            if not attempt_specs:
                attempt_specs = [("blind", None)]
            attempt_notes: list[str] = []
            ok = False
            dt = 0.0
            out_s = ""
            err_s = ""
            cmd = []
            solution_path = None
            attempt_label = ""
            attempt_coord_hint = ""
            for idx, (attempt_label, attempt_coord) in enumerate(attempt_specs, start=1):
                attempt_coord_hint = _format_coord_hint(attempt_coord)
                is_blind_retry = attempt_coord is None and has_hint_attempt
                attempt_cpulimit_s = blind_cpulimit_s if is_blind_retry else cpulimit_s
                LOG(
                    f"{filename}: attempt {idx}/{len(attempt_specs)} "
                    f"hint={attempt_label} center={attempt_coord_hint} "
                    f"cpulimit={attempt_cpulimit_s:.0f}s"
                )
                ok, dt, out_s, err_s, cmd, solution_path = self._run_solve_field(
                    fits_path, attempt_coord, scale_low, scale_high, radius_deg,
                    downsample, timeout_s, outdir, use_wsl, True, use_cache, max_objs, attempt_cpulimit_s
                )
                attempt_notes.append(
                    f"{attempt_label}@{attempt_coord_hint}:{'ok' if ok else 'fail'}:{float(dt):.1f}s:cpu={attempt_cpulimit_s:.0f}"
                )
                if ok:
                    break

            cmd_str = " ".join(str(c) for c in cmd) if cmd else ""
            stdout_tail = _tail_text(out_s, limit=2000, max_lines=12)
            stderr_tail = _tail_text(err_s, limit=2000, max_lines=12)
            attempt_summary = " | ".join(attempt_notes)
            fail_reason = ""
            if not ok:
                err_l = str(err_s).lower()
                out_l = str(out_s).lower()
                if "timeout" in err_l:
                    fail_reason = "astnet_timeout"
                elif self._stop_requested or _is_explicit_stopped_message(err_s):
                    fail_reason = "astnet_stopped"
                elif (
                    "wsl_path_unavailable" in err_l
                    or "wsl_not_found" in err_l
                    or "wsl_input_missing" in err_l
                    or "input_missing:" in err_l
                    or "wsl_stage_failed:" in err_l
                ):
                    fail_reason = "astnet_wsl_path_unavailable"
                elif (
                    "backend.cfg" in err_l
                    or "index-" in err_l
                    or "index file" in err_l
                    or "/usr/share/astrometry" in err_l
                ):
                    fail_reason = "astnet_index_missing"
                elif "cannot open `" in err_l and "no such file or directory" in err_l:
                    if filename.lower() in err_l:
                        fail_reason = "astnet_wsl_input_unavailable"
                    else:
                        fail_reason = "astnet_missing_file_or_index"
                elif "cache_hit" in str(out_s):
                    fail_reason = "astnet_cache_hit_without_solution_marker"
                else:
                    fail_reason = "astnet_fail_no_solution"

            result = {
                "fname": filename,
                "file": filename,
                "ok": False,
                "wcs_ok": False,
                "status": "fail",
                "fail_reason": fail_reason,
                "ra": 0.0, "dec": 0.0, "pixscale": 0.0,
                "elapsed_s": float(dt),
                "solver": "astnet_wsl",
                "astnet_wsl_hint_source": attempt_label,
                "astnet_wsl_center_hint": attempt_coord_hint,
                "astnet_wsl_attempts": attempt_summary,
                "astnet_wsl_cmd": cmd_str,
                "astnet_wsl_stdout": stdout_tail,
                "astnet_wsl_stderr": stderr_tail,
            }

            if ok and solution_path is not None and solution_path.exists():
                try:
                    with fits.open(solution_path, memmap=False) as hdul_new:
                        new_hdr = hdul_new[0].header
                        w = WCS(new_hdr, relax=True)
                        if w.has_celestial:
                            self._safe_header_update(fits_path, new_hdr)
                            self._refresh_detect_cache_signature(fits_path, filename)

                            pix_fit = float(np.mean(proj_plane_pixel_scales(w.celestial) * 3600.0))
                            with fits.open(fits_path, memmap=False) as hdul_src:
                                src_hdr = hdul_src[0].header
                            nx, ny = _solution_header_shape(new_hdr, src_hdr)
                            cx = nx / 2.0
                            cy = ny / 2.0
                            ra_dec = w.pixel_to_world(cx, cy)

                            result = {
                                "fname": filename,
                                "file": filename,
                                "ok": True,
                                "wcs_ok": True,
                                "status": "solved",
                                "ra": float(ra_dec.ra.deg),
                                "dec": float(ra_dec.dec.deg),
                                "pixscale": pix_fit,
                                "elapsed_s": float(dt),
                                "solver": "astnet_wsl",
                                "fail_reason": "",
                                "astnet_wsl_hint_source": attempt_label,
                                "astnet_wsl_center_hint": attempt_coord_hint,
                                "astnet_wsl_attempts": attempt_summary,
                                "astnet_wsl_cmd": cmd_str,
                                "astnet_wsl_stdout": stdout_tail,
                                "astnet_wsl_stderr": stderr_tail,
                                "wcs_header": dict(new_hdr),
                                "fits_path": str(fits_path),
                            }
                            LOG(
                                f"{filename}: solved dt={dt:.1f}s RA={result['ra']:.6f} "
                                f"Dec={result['dec']:.6f} pix={pix_fit:.5f} "
                                f"hint={attempt_label} center={attempt_coord_hint}"
                            )
                            if stderr_tail:
                                LOG(f"{filename}: solver_stderr_tail={stderr_tail}")
                        else:
                            fail_reason = "wcs_header_not_celestial"
                            result["status"] = "wcs_not_celestial"
                            result["fail_reason"] = fail_reason
                            LOG(f"{filename}: fail reason={fail_reason} dt={dt:.1f}s")
                            if cmd_str:
                                LOG(f"{filename}: cmd={cmd_str}")
                            if stdout_tail:
                                LOG(f"{filename}: stdout_tail={stdout_tail}")
                            if stderr_tail:
                                LOG(f"{filename}: stderr_tail={stderr_tail}")
                except Exception as e:
                    fail_reason = f"header_update_error:{_exc_brief(e, limit=160)}"
                    result = {
                        "fname": filename,
                        "file": filename,
                        "ok": False,
                        "wcs_ok": False,
                        "status": f"error: {e}",
                        "fail_reason": fail_reason,
                        "elapsed_s": float(dt),
                        "solver": "astnet_wsl",
                        "astnet_wsl_hint_source": attempt_label,
                        "astnet_wsl_center_hint": attempt_coord_hint,
                        "astnet_wsl_attempts": attempt_summary,
                        "astnet_wsl_cmd": cmd_str,
                        "astnet_wsl_stdout": stdout_tail,
                        "astnet_wsl_stderr": stderr_tail,
                    }
                    LOG(f"{filename}: {fail_reason}")
                    if cmd_str:
                        LOG(f"{filename}: cmd={cmd_str}")
                    if stdout_tail:
                        LOG(f"{filename}: stdout_tail={stdout_tail}")
                    if stderr_tail:
                        LOG(f"{filename}: stderr_tail={stderr_tail}")
            else:
                if not fail_reason:
                    fail_reason = "astnet_no_solution_file"
                result["fail_reason"] = fail_reason
                LOG(f"{filename}: fail reason={fail_reason} dt={dt:.1f}s")
                if cmd_str:
                    LOG(f"{filename}: cmd={cmd_str}")
                if stdout_tail:
                    LOG(f"{filename}: stdout_tail={stdout_tail}")
                if stderr_tail:
                    LOG(f"{filename}: stderr_tail={stderr_tail}")
                if attempt_summary:
                    LOG(f"{filename}: attempts={attempt_summary}")

            # 임시 파일 정리
            if not keep_outputs:
                for p in outdir.glob(f"{fits_path.stem}.*"):
                    try: p.unlink()
                    except Exception: pass

            return filename, result

        file_results = {}  # filename -> result 매핑
        from queue import Queue as _Queue
        _slot_q: _Queue = _Queue()
        for _i in range(max(1, max_workers)):
            _slot_q.put(_i)

        def _slotted_process(filename: str):
            slot = _slot_q.get()
            try:
                self.worker_status.emit(slot, filename, "Astnet", 5)
                out = process_single_file(filename)
                self.worker_status.emit(slot, filename, "Done", 100)
                return out
            finally:
                _slot_q.put(slot)

        # Explicit ownership allows queued solves to be cancelled on Stop.
        executor = ThreadPoolExecutor(max_workers=max_workers)
        try:
            # 모든 작업을 큐에 등록
            future_to_file = {executor.submit(_slotted_process, f): f for f in self.file_list}

            completed_count = 0
            for future in as_completed(future_to_file):
                if self._stop_requested:
                    for f_cancel in future_to_file:
                        f_cancel.cancel()
                    break

                fname = future_to_file[future]
                try:
                    filename, res = future.result()
                    res["filename"] = filename  # 결과에 파일명 저장
                    results.append(res)
                    file_results[filename] = res

                    # UI 업데이트 시그널
                    self.file_done.emit(filename, res)

                    # 로그 출력 (성공/실패 여부에 따라)
                    if res.get("ok"):
                        ra_val = res.get('ra', 0)
                        dec_val = res.get('dec', 0)
                        self.log_message.emit(f"[OK] {filename} (RA={ra_val:.4f}, Dec={dec_val:.4f})")
                    else:
                        status_txt = str(res.get("status", "fail"))
                        fail_txt = str(res.get("fail_reason", "")).strip()
                        if fail_txt:
                            self.log_message.emit(f"[FAIL] {filename}: {status_txt} ({fail_txt})")
                        else:
                            self.log_message.emit(f"[FAIL] {filename}: {status_txt}")
                        try:
                            elapsed_val = float(res.get("elapsed_s", np.nan))
                        except Exception:
                            elapsed_val = np.nan
                        LOG(
                            f"{filename}: status={status_txt} fail_reason={fail_txt or '-'} "
                            f"elapsed={elapsed_val:.1f}s"
                        )

                except Exception as e:
                    LOG(f"{fname}: worker_exception={_exc_brief(e)}")
                    self.error.emit(fname, str(e))

                completed_count += 1
                self.progress.emit(completed_count, total, f"Solved {completed_count}/{total}")
        finally:
            executor.shutdown(wait=True, cancel_futures=True)

        # --- Gaia 쿼리 및 WCS Refine ---
        if self._stop_requested:
            self.finished.emit({
                "total": len(results),
                "ok": sum(1 for r in results if r.get("ok")),
                "wcs_qc_pass": sum(1 for r in results if bool(r.get("wcs_qc_pass", False))),
                "stopped": True,
            })
            return

        # 성공한 프레임에서 중심 좌표 얻기
        solved_frames = [r for r in results if r.get("ok")]
        center_coord = None
        if solved_frames:
            first_ok = solved_frames[0]
            ra_center = float(first_ok.get("ra", 0))
            dec_center = float(first_ok.get("dec", 0))
            if np.isfinite(ra_center) and np.isfinite(dec_center):
                center_coord = SkyCoord(ra=ra_center * u.deg, dec=dec_center * u.deg)

        gaia_df = None
        if center_coord is not None:
            self.log_message.emit("[Gaia] Querying Gaia catalog for ID matching...")
            self.progress.emit(total, total, "Querying Gaia catalog...")

            # FOV 계산
            sample_fname = self.file_list[0]
            if self.use_cropped:
                sample_path = step2_cropped_dir(self.result_dir) / sample_fname
            else:
                sample_path = self.params.get_file_path(sample_fname)
            try:
                with fits.open(sample_path, memmap=False) as hdul:
                    ny, nx = hdul[0].data.shape
                fov_w = (nx * pix_arc) / 3600.0
                fov_h = (ny * pix_arc) / 3600.0
                diag_deg = float(np.hypot(fov_w, fov_h))
                gaia_fudge = float(getattr(self.params.P, "gaia_radius_fudge", 1.35))
                gaia_r = float(0.5 * diag_deg * gaia_fudge)

                gaia_df, gaia_src = self._load_or_query_gaia(center_coord, gaia_r)
                self.log_message.emit(f"[Gaia] center=({center_coord.ra.deg:.6f},{center_coord.dec.deg:.6f}) r={gaia_r:.4f}deg source={gaia_src} N={len(gaia_df)}")
            except Exception as e:
                msg = _exc_brief(e, limit=240)
                self.log_message.emit(f"[Gaia] Query error: {msg}")
                LOG(f"[Gaia] Query error: {msg}")
                gaia_df = pd.DataFrame()

        # WCS Refine 수행
        if gaia_df is not None and len(gaia_df) > 0:
            self.log_message.emit("[Refine] Starting WCS refinement with Gaia matching...")
            refine_max_match = int(getattr(self.params.P, "wcs_refine_max_match", 300))

            for i, res in enumerate(solved_frames):
                if self._stop_requested:
                    break

                filename = res.get("filename", "")
                fits_path_str = res.get("fits_path", "")
                if not fits_path_str:
                    continue

                fits_path = Path(fits_path_str)
                if not fits_path.exists():
                    continue

                try:
                    # Detection 로드
                    det_xy, _ = self._load_detect_xy(filename)
                    fwhm_px, _ = self._load_fwhm_for_frame(filename)

                    with fits.open(fits_path, mode="update", memmap=False) as hdul:
                        hdr = hdul[0].header
                        w = WCS(hdr, relax=True)

                        ok_refine, refine_note, resid_med, resid_max, match_n = self._refine_crpix_by_match(
                            w, hdr, det_xy, gaia_df, fwhm_px, refine_max_match
                        )

                        # 결과 업데이트
                        res["refine"] = refine_note
                        res["resid_med"] = resid_med
                        res["resid_max"] = resid_max
                        res["match_n"] = match_n

                        if ok_refine:
                            hdr["WCSSRC"] = ("ASTNET_REFINED", "WCS source (refined with Gaia)")
                            self.log_message.emit(f"  [Refine] {filename}: {refine_note}, resid_med={resid_med:.2f}\"")
                        else:
                            self.log_message.emit(f"  [Refine] {filename}: skip - {refine_note}")

                    # Refine 결과 시그널 emit
                    self.refine_done.emit(filename, {
                        "refine": res.get("refine", "-"),
                        "resid_med": res.get("resid_med", np.nan),
                        "resid_max": res.get("resid_max", np.nan),
                        "match_n": res.get("match_n", 0),
                    })

                except Exception as e:
                    self.log_message.emit(f"  [Refine] {filename}: error - {e}")
                    res["refine"] = f"error:{e}"
                    res["resid_med"] = np.nan
                    res["resid_max"] = np.nan
                    res["match_n"] = 0
                    self.refine_done.emit(filename, {
                        "refine": res["refine"],
                        "resid_med": res["resid_med"],
                        "resid_max": res["resid_max"],
                        "match_n": res["match_n"],
                    })

                self.progress.emit(total + i + 1, total + len(solved_frames), f"Refining {i+1}/{len(solved_frames)}")
        else:
            self.log_message.emit("[Refine] Skipped - no Gaia data available")

        gaia_ra_vals = np.array([], dtype=float)
        gaia_dec_vals = np.array([], dtype=float)
        if isinstance(gaia_df, pd.DataFrame) and (not gaia_df.empty) and {"ra", "dec"} <= set(gaia_df.columns):
            gaia_ra_vals = pd.to_numeric(gaia_df["ra"], errors="coerce").to_numpy(float)
            gaia_dec_vals = pd.to_numeric(gaia_df["dec"], errors="coerce").to_numpy(float)
            ok_rd = np.isfinite(gaia_ra_vals) & np.isfinite(gaia_dec_vals)
            gaia_ra_vals = gaia_ra_vals[ok_rd]
            gaia_dec_vals = gaia_dec_vals[ok_rd]

        for res in results:
            filename = str(res.get("filename") or res.get("fname") or res.get("file") or "").strip()
            if not filename:
                continue

            det_xy, _ = self._load_detect_xy(filename)
            qc_metrics = self._empty_wcs_qc_metrics(n_detect=len(det_xy))
            qc_metrics["pix_scale_input_arcsec"] = float(pix_arc) if np.isfinite(pix_arc) else np.nan

            fits_path = None
            fp = str(res.get("fits_path", "")).strip()
            if fp:
                fits_path = Path(fp)
            if fits_path is None or not fits_path.exists():
                if self.use_cropped:
                    cp = step2_cropped_dir(self.result_dir) / filename
                    if cp.exists():
                        fits_path = cp
                if (fits_path is None or not fits_path.exists()):
                    try:
                        op = Path(self.params.get_file_path(filename))
                        if op.exists():
                            fits_path = op
                    except Exception:
                        fits_path = None

            w = None
            wcs_ok = False
            nx = 0
            ny = 0
            wcs_rot_deg = np.nan
            center_ra = np.nan
            center_dec = np.nan
            pix_fit = float(res.get("pixscale", np.nan))
            if fits_path is not None and fits_path.exists():
                try:
                    with fits.open(fits_path, memmap=False) as hdul:
                        hdr = hdul[0].header
                        nx = int(hdr.get("NAXIS1", 0))
                        ny = int(hdr.get("NAXIS2", 0))
                        w = WCS(hdr, relax=True)
                    wcs_ok = bool(w is not None and w.has_celestial)
                    if wcs_ok:
                        pix_fit = self._pixscale_from_wcs(w)
                        wcs_rot_deg = self._wcs_rotation_deg(w)
                        center_ra, center_dec = self._wcs_center_coords(w, nx, ny)
                except Exception:
                    w = None
                    wcs_ok = False

            qc_metrics = self._compute_wcs_qc_metrics(
                w=w,
                det_xy=det_xy,
                nx=nx,
                ny=ny,
                gaia_ra_deg=gaia_ra_vals,
                gaia_dec_deg=gaia_dec_vals,
                pix_input_arcsec=float(pix_arc),
                pix_fit_arcsec=float(pix_fit) if np.isfinite(pix_fit) else np.nan,
                center_coord=center_coord,
            )
            qc_pass, qc_reasons = self._evaluate_wcs_qc_pass(qc_metrics, wcs_ok=bool(wcs_ok))

            elapsed_s = float(res.get("elapsed_s", np.nan))
            elapsed = elapsed_s if np.isfinite(elapsed_s) else np.nan
            res.update({
                "fname": filename,
                "file": filename,
                "wcs_ok": bool(wcs_ok),
                "pix_fit": float(pix_fit) if np.isfinite(pix_fit) else None,
                "elapsed": float(elapsed) if np.isfinite(elapsed) else None,
                "wcs_rot_deg": float(wcs_rot_deg) if np.isfinite(wcs_rot_deg) else None,
                "center_ra_deg": float(center_ra) if np.isfinite(center_ra) else None,
                "center_dec_deg": float(center_dec) if np.isfinite(center_dec) else None,
                "wcs_qc_pass": bool(qc_pass),
                "wcs_qc_reason": ",".join(qc_reasons),
            })
            res.update(qc_metrics)

        try:
            step5_out = step5_wcs_dir(self.result_dir)
            step5_out.mkdir(parents=True, exist_ok=True)
            df = pd.DataFrame(results)
            df.to_csv(step5_out / "wcs_solve_summary.csv", index=False)
            qc_cols = [
                "file",
                "status",
                "fail_reason",
                "ok",
                "wcs_ok",
                "solver",
                "elapsed",
                "elapsed_s",
                "n_detect",
                "n_catalog_in_fov",
                "n_match",
                "n_inlier",
                "match_rate",
                "match_rate_cat",
                "match_rate_eff",
                "match_radius_arcsec",
                "match_radius_px",
                "dx_med_px",
                "dy_med_px",
                "resid_med_px",
                "resid_mad_px",
                "resid_p99_px",
                "resid_peak_px",
                "rms_px",
                "inlier_rate",
                "resid_vs_radius_slope",
                "edge_resid_ratio",
                "pix_scale_input_arcsec",
                "pixscale",
                "pix_fit",
                "scale_delta_pct",
                "wcs_rot_deg",
                "center_ra_deg",
                "center_dec_deg",
                "center_offset_arcsec",
                "wcs_qc_pass",
                "wcs_qc_reason",
            ]
            qc_df = df[[c for c in qc_cols if c in df.columns]].copy()
            qc_df.to_csv(step5_out / "frame_wcs_qc.csv", index=False)
        except Exception as e:
            self.log_message.emit(f"[WCS-QC] Failed to write QC CSV: {e}")

        # --- 마무리 ---
        n_qc_pass = sum(1 for r in results if bool(r.get("wcs_qc_pass", False)))
        n_qc_not_eval = sum(
            1
            for r in results
            if "gaia_unavailable" in str(r.get("wcs_qc_reason", ""))
        )
        summary = {
            "total": len(results),
            "ok": sum(1 for r in results if r.get("ok")),
            "wcs_qc_pass": int(n_qc_pass),
            "wcs_qc_not_evaluated": int(n_qc_not_eval),
            "stopped": bool(self._stop_requested),
        }
        self.finished.emit(summary)


class InternalWcsWorker(QThread):
    """WCS solver using the built-in Python engine (no external executables)."""

    progress      = pyqtSignal(int, int, str)
    log_message   = pyqtSignal(str)
    frame_done    = pyqtSignal(str, dict)   # (fname, result_info)
    finished      = pyqtSignal(dict)
    error         = pyqtSignal(str)
    worker_status = pyqtSignal(int, str, str, int)  # (slot, fname, status, pct)

    def __init__(self, file_list, params, data_dir, result_dir,
                 use_cropped=False,
                 min_matches=8,
                 rms_max_arcsec=2.0,
                 advanced_params=None):
        super().__init__()
        self.file_list         = list(file_list)
        self.params            = params
        self.data_dir          = Path(data_dir)
        self.result_dir        = Path(result_dir)
        self.use_cropped       = use_cropped
        self.min_matches       = int(min_matches)
        self.rms_max_arcsec    = float(rms_max_arcsec)
        self.advanced_params   = dict(advanced_params) if advanced_params else {}
        self._stop             = False

    def stop(self):
        self._stop = True

    def _log(self, msg: str) -> None:
        self.log_message.emit(str(msg))

    def _load_or_fetch_gaia(
        self,
        step5_out: Path,
        approx_ra: float,
        approx_dec: float,
        P,
        radius_deg: float | None = None,
    ):
        """Load Gaia catalog from cache, or query it if not present.

        Returns ``(gaia_ra, gaia_dec, gaia_g, gaia_pmra, gaia_pmdec)`` numpy
        arrays, or all-None on failure. PM columns are NaN-replaced with 0.0
        so the solver can apply ``ra + pmra·dt/cos(dec)`` without branching on
        missing data per star — Gaia DR3 has ≳20 % of faint sources without
        PM, and treating those as zero-motion is the safe default at the
        catalogue-fit residual level.
        """
        try:
            service = GaiaCatalogService(
                self.params,
                self.result_dir,
                log_fn=self._log,
                stop_fn=lambda: self._stop,
            )
            center = SkyCoord(float(approx_ra) * u.deg, float(approx_dec) * u.deg)
            radius_use = float(radius_deg) if radius_deg is not None else float("nan")
            if not (0.05 < radius_use < 5.0):
                fov_w = float(getattr(P, "fov_w_deg", 0.0))
                fov_h = float(getattr(P, "fov_h_deg", 0.0))
                radius_use = float(np.hypot(fov_w, fov_h)) / 2.0
            if not (0.05 < radius_use < 5.0):
                radius_use = 0.4
            df, source = service.load_or_query(center, radius_use)
            self._log(f"[Internal WCS] Gaia catalog source={source} N={len(df)}")
            if df is None or df.empty:
                return None, None, None, None, None
            col_ra = "ra" if "ra" in df.columns else "ra_deg" if "ra_deg" in df.columns else None
            col_dec = "dec" if "dec" in df.columns else "dec_deg" if "dec_deg" in df.columns else None
            if col_ra is None or col_dec is None:
                self._log("[Internal WCS] Gaia catalog has no RA/Dec columns.")
                return None, None, None, None, None
            ra = pd.to_numeric(df[col_ra], errors="coerce").to_numpy(float)
            dec = pd.to_numeric(df[col_dec], errors="coerce").to_numpy(float)
            if "phot_g_mean_mag" in df.columns:
                g = pd.to_numeric(df["phot_g_mean_mag"], errors="coerce").to_numpy(float)
            else:
                g = None

            def _col_zero_filled(name: str):
                if name not in df.columns:
                    return None
                v = pd.to_numeric(df[name], errors="coerce").to_numpy(float)
                v = np.where(np.isfinite(v), v, 0.0)
                return v

            pmra = _col_zero_filled("pmra")
            pmdec = _col_zero_filled("pmdec")

            ok = np.isfinite(ra) & np.isfinite(dec)
            if g is not None and len(g) == len(ok):
                g = g[ok]
            if pmra is not None and len(pmra) == len(ok):
                pmra = pmra[ok]
            if pmdec is not None and len(pmdec) == len(ok):
                pmdec = pmdec[ok]
            return ra[ok], dec[ok], g, pmra, pmdec
        except Exception as exc:
            self._log(f"[Internal WCS] Shared Gaia service failed: {exc}")
            return None, None, None, None, None

    def _fits_path(self, fname: str) -> Path:
        if self.use_cropped:
            p = step2_cropped_dir(self.result_dir) / fname
            if p.exists():
                return p
        return self.data_dir / fname

    def _detect_sources_from_fits(self, fits_path: Path, max_sources: int = 2500):
        """Small internal fallback detector for the WCS solver.

        Step 5 normally consumes Step 4 detections, but a plate solver should not
        become unusable just because the detection cache was not produced or was
        too sparse.  This intentionally stays simple: robust background,
        local-maxima, and brightest peaks.
        """
        try:
            from scipy.ndimage import maximum_filter

            data = fits.getdata(fits_path)
            data = np.asarray(data, dtype=float)
            if data.ndim > 2:
                data = np.squeeze(data)
            if data.ndim != 2:
                return None, None

            finite = np.isfinite(data)
            if not finite.any():
                return None, None
            vals = data[finite]
            med = float(np.nanmedian(vals))
            mad = float(MAD_TO_SIGMA * np.nanmedian(np.abs(vals - med)))
            if not np.isfinite(mad) or mad <= 0:
                mad = float(np.nanstd(vals))
            if not np.isfinite(mad) or mad <= 0:
                return None, None

            threshold = med + 5.0 * mad
            local_max = maximum_filter(data, size=5, mode="nearest")
            mask = finite & (data == local_max) & (data > threshold)
            mask[:8, :] = False
            mask[-8:, :] = False
            mask[:, :8] = False
            mask[:, -8:] = False
            yy, xx = np.nonzero(mask)
            if len(xx) == 0:
                return None, None

            flux = data[yy, xx] - med
            order = np.argsort(np.where(np.isfinite(flux), flux, -np.inf))[::-1]
            order = order[:max(1, int(max_sources))]
            xy = np.column_stack([xx[order].astype(float), yy[order].astype(float)])
            flux = flux[order].astype(float)
            return xy, flux
        except Exception as exc:
            self._log(f"[Internal WCS] fallback detector failed: {exc}")
            return None, None

    def run(self):
        """Solve WCS for every queued frame.

        The entire body is wrapped in try/finally so that the worker
        ALWAYS emits ``finished`` even when something silently explodes
        (e.g. a missing PyInstaller hidden import for wcs_solver).
        Without this guarantee the UI's Run button stays disabled and
        the user has no way to recover short of closing the window.
        """
        import traceback as _tb_internal

        crash_reason: str | None = None
        n_ok = 0
        total = len(self.file_list)
        try:
            self._run_body()
            return
        except ImportError as exc:
            crash_reason = (
                f"Internal solver dependency missing: {exc}.  "
                "If this is an installed build, the wcs_solver module or one of "
                "its scientific deps (scipy/astropy) was not bundled."
            )
            self._log(f"[Internal WCS] {crash_reason}")
            self._log(_tb_internal.format_exc())
        except Exception as exc:
            crash_reason = f"Internal solver crashed: {exc}"
            self._log(f"[Internal WCS] {crash_reason}")
            self._log(_tb_internal.format_exc())
        finally:
            if crash_reason is not None:
                try:
                    self.error.emit(crash_reason)
                except Exception:
                    pass
                try:
                    self.finished.emit({
                        "total": total,
                        "ok": n_ok,
                        "stopped": self._stop,
                        "crashed": True,
                    })
                except Exception:
                    pass

    def _run_body(self):
        from apex.analysis.astrometry import solve as solve_wcs_internal
        from apex.utils.io_utils import read_ecsv_int64_source_id

        P = self.params.P
        result_dir = self.result_dir
        step5_out  = step5_wcs_dir(result_dir)
        step5_out.mkdir(parents=True, exist_ok=True)

        # ── Pixel scale ───────────────────────────────────────────────────
        approx_scale = float(getattr(P, "pixel_scale_arcsec", 0.0))
        if not (approx_scale > 0):
            self.error.emit(
                "pixel_scale_arcsec not set in parameters.\n"
                "Configure telescope/camera in parameters.toml first."
            )
            return

        # ── Target coords ─────────────────────────────────────────────────
        approx_ra  = float(getattr(P, "target_ra_deg",  float("nan")))
        approx_dec = float(getattr(P, "target_dec_deg", float("nan")))

        # Fallback: load from Step 1 SIMBAD-resolved targets file
        if not (np.isfinite(approx_ra) and np.isfinite(approx_dec)):
            from apex.utils.step_paths import step1_dir
            simbad_path = step1_dir(result_dir) / "targets_simbad.tsv"
            if simbad_path.exists():
                try:
                    sdf = pd.read_csv(simbad_path, sep="\t")
                    if {"ra_deg", "dec_deg"} <= set(sdf.columns) and not sdf.empty:
                        approx_ra  = float(sdf["ra_deg"].iloc[0])
                        approx_dec = float(sdf["dec_deg"].iloc[0])
                        name = str(sdf.get("name", pd.Series([""])).iloc[0])
                        self._log(
                            f"[Internal WCS] Target loaded from Step 1: "
                            f"{name} ({approx_ra:.4f}, {approx_dec:.4f})"
                        )
                except Exception as exc:
                    self._log(f"[Internal WCS] Step 1 SIMBAD read error: {exc}")

        if not (np.isfinite(approx_ra) and np.isfinite(approx_dec)):
            self.error.emit(
                "Target RA/Dec not found.\n\n"
                "Resolve the target in Step 1 first (writes step1_file_selection/targets_simbad.tsv)."
            )
            return

        # ── Load Gaia catalog (download if not cached) ───────────────────
        gaia_radius_deg = None
        try:
            for probe_name in self.file_list:
                probe_path = self._fits_path(probe_name)
                if not probe_path.exists():
                    continue
                with fits.open(probe_path, memmap=False) as hdul:
                    ph = hdul[0].header
                    pny = int(ph.get("NAXIS2", 0))
                    pnx = int(ph.get("NAXIS1", 0))
                if pnx > 0 and pny > 0:
                    diag_deg = float(np.hypot(pnx, pny) * approx_scale / 3600.0)
                    gaia_fudge = float(getattr(P, "gaia_radius_fudge", 1.35))
                    # Use the SAME radius formula as the ASTAP worker
                    # (0.5·diag·fudge).  The old ×1.25 made the Internal
                    # tab download ~56 % more stars AND mismatch the ASTAP
                    # cache, forcing a fresh Gaia query every time the user
                    # switched tabs.  Shared radius = shared gaia_fov.ecsv.
                    gaia_radius_deg = max(0.4, 0.5 * diag_deg * max(gaia_fudge, 1.0))
                    self._log(f"[Internal WCS] Gaia query radius={gaia_radius_deg:.4f}deg")
                    break
        except Exception as exc:
            self._log(f"[Internal WCS] Gaia radius estimate failed: {exc}")

        # Surface "querying Gaia" feedback immediately, like the ASTAP
        # worker does — otherwise the Worker Status panel sits blank for
        # the whole (blocking) download and the step looks frozen.
        try:
            self.progress.emit(0, len(self.file_list), "Querying Gaia catalog...")
            self.worker_status.emit(0, "Gaia catalog", "Querying", 5)
        except Exception:
            pass
        gaia_ra, gaia_dec, gaia_g, gaia_pmra, gaia_pmdec = self._load_or_fetch_gaia(
            step5_out, approx_ra, approx_dec, P, radius_deg=gaia_radius_deg,
        )
        try:
            self.worker_status.emit(0, "Gaia catalog", "Ready", 100)
        except Exception:
            pass
        if gaia_ra is None:
            self.error.emit(
                "Gaia catalog unavailable and download failed.\n\n"
                "Check internet connection or run the ASTAP tab once to cache "
                "the Gaia catalogue."
            )
            return

        # ── Per-frame solve (parallel across frames) ─────────────────────
        # Matches the ASTAP / astnet worker pattern: extract per-frame work
        # into a method, submit to ThreadPoolExecutor sized by
        # get_parallel_workers(self.params). The heavy CPU inside the solver
        # (KDTree builds, RANSAC, lstsq) is numpy/scipy and releases the GIL,
        # so threads actually overlap. Each frame has its own FITS path so
        # the in-place header write is race-free.
        total   = len(self.file_list)
        n_ok    = 0
        n_qc_pass = 0
        n_qc_not_evaluated = 0
        results = {}

        max_workers = get_parallel_workers(self.params)
        max_workers = max(1, int(max_workers))
        self._log(f"[Internal WCS] Solving {total} frames with {max_workers} worker(s)")

        # Target SkyCoord for the per-frame hint chain (header -> target ->
        # blind). The single Gaia cone above stays target-centred (mosaic is
        # out of scope), so for non-mosaic data every candidate resolves
        # against the same catalogue; the chain just rescues frames whose
        # header pointing is stale/missing or whose target is slightly off.
        target_coord = None
        if np.isfinite(approx_ra) and np.isfinite(approx_dec):
            try:
                target_coord = SkyCoord(float(approx_ra) * u.deg, float(approx_dec) * u.deg)
            except Exception:
                target_coord = None

        bundle = dict(
            gaia_ra=gaia_ra, gaia_dec=gaia_dec, gaia_g=gaia_g,
            gaia_pmra=gaia_pmra, gaia_pmdec=gaia_pmdec,
            approx_ra=approx_ra, approx_dec=approx_dec,
            approx_scale=approx_scale,
            result_dir=result_dir,
            solve_fn=solve_wcs_internal,
            target_coord=target_coord,
        )

        # Per-thread slot index for the Worker Status panel, same pattern as the
        # ASTAP / astnet workers: a bounded queue hands each running thread a
        # stable slot number so the panel shows one row per worker.
        from queue import Queue as _Queue
        _slot_q: _Queue = _Queue()
        for _i in range(max_workers):
            _slot_q.put(_i)

        def _slotted_solve(filename: str):
            slot = _slot_q.get()
            try:
                self.worker_status.emit(slot, filename, "Solving", 5)
                out = self._solve_one_internal_frame(filename, bundle)
                status = "Done" if out.get("ok") else "Failed"
                self.worker_status.emit(slot, filename, status, 100)
                return out
            finally:
                _slot_q.put(slot)

        done_count = 0
        # Keep ownership explicit so queued work can be cancelled on Stop.
        ex = ThreadPoolExecutor(max_workers=max_workers)
        try:
            fut_to_name = {
                ex.submit(_slotted_solve, f): f
                for f in self.file_list
            }
            for fut in as_completed(fut_to_name):
                if self._stop:
                    # Running frames are joined below before finished emits.
                    for f_cancel in fut_to_name:
                        f_cancel.cancel()
                    break
                fname = fut_to_name[fut]
                try:
                    info = fut.result()
                except Exception as exc:
                    info = {"ok": False, "reason": f"unhandled: {exc}"}
                    self._log(f"[Internal WCS] Unhandled error {fname}: {exc}")
                results[fname] = info
                if info.get("ok"):
                    n_ok += 1
                    if not bool(info.get("gaia_available", True)):
                        n_qc_not_evaluated += 1
                    elif bool(info.get("wcs_qc_pass", False)):
                        n_qc_pass += 1
                done_count += 1
                self.progress.emit(done_count, total, fname)
                self.frame_done.emit(fname, info)
        finally:
            ex.shutdown(wait=True, cancel_futures=True)

        self.finished.emit({
            "total": total,
            "ok": n_ok,
            "wcs_qc_pass": n_qc_pass,
            "wcs_qc_not_evaluated": n_qc_not_evaluated,
            "stopped": self._stop,
        })

    def _solve_one_internal_frame(self, fname: str, bundle: dict) -> dict:
        """Solve WCS for one frame. Thread-safe — only touches the per-frame
        FITS file and uses signal-based logging back to the main thread.
        """
        gaia_ra      = bundle["gaia_ra"]
        gaia_dec     = bundle["gaia_dec"]
        gaia_g       = bundle["gaia_g"]
        gaia_pmra    = bundle["gaia_pmra"]
        gaia_pmdec   = bundle["gaia_pmdec"]
        approx_ra    = bundle["approx_ra"]
        approx_dec   = bundle["approx_dec"]
        approx_scale = bundle["approx_scale"]
        result_dir   = bundle["result_dir"]
        solve_wcs_internal = bundle["solve_fn"]
        target_coord = bundle.get("target_coord")

        _t_frame_start = time.time()

        fits_path = self._fits_path(fname)
        if not fits_path.exists():
            self._log(f"[Internal WCS] Missing: {fname}")
            return {"ok": False, "reason": "fits_not_found",
                    "elapsed_s": time.time() - _t_frame_start}

        # Load FITS shape + observation epoch (J-year) for Gaia PM update,
        # plus this frame's own pointing for the hint chain.
        obs_epoch_jyear: float | None = None
        header_coord = None
        try:
            with fits.open(fits_path, memmap=False) as hdul:
                hdr = hdul[0].header.copy()
                _date_obs = hdr.get("DATE-OBS")
                if _date_obs:
                    try:
                        from astropy.time import Time as _AT
                        obs_epoch_jyear = float(_AT(_date_obs, scale="utc").jyear)
                    except Exception:
                        obs_epoch_jyear = None
                ny = int(hdr.get("NAXIS2", 0))
                nx = int(hdr.get("NAXIS1", 0))
                header_coord = _header_pointing_coord(hdr)
        except Exception as exc:
            self._log(f"[Internal WCS] FITS read error {fname}: {exc}")
            return {"ok": False, "reason": "fits_error",
                    "elapsed_s": time.time() - _t_frame_start}

        if ny == 0 or nx == 0:
            return {"ok": False, "reason": "bad_shape",
                    "elapsed_s": time.time() - _t_frame_start}

        # Load Step 4 detections and build an ORDERED LADDER of source sets to
        # try: isolated "anchor" stars first (cleanest WCS — best RMS on normal
        # fields), then the raw brightest detections as a safety net. The
        # benchmark showed that in dense globular cores the isolation filter can
        # decimate the core and starve quad matching, while the raw brightest
        # set still solves; the ladder gives the clean path when it works and
        # the robust path when it doesn't.
        det_candidates = [
            self.params.P.cache_dir / f"detect_{fname}.csv",
            step4_dir(result_dir) / f"detect_{fname}.csv",
        ]
        src_sets: list[tuple[str, np.ndarray, np.ndarray | None]] = []
        det_xy_for_qc = np.zeros((0, 2), dtype=float)
        ANCHOR_MIN = 12  # below this an anchor-only set is too sparse to bother
        for dp in det_candidates:
            if not dp.exists():
                continue
            try:
                ddf = pd.read_csv(dp)
                if not ({"x", "y"} <= set(ddf.columns)):
                    continue
                xy_all = ddf[["x", "y"]].to_numpy(float)
                valid_all = np.isfinite(xy_all).all(axis=1)
                det_xy_for_qc = xy_all[valid_all]
                flux_col = next(
                    (c for c in ("dao_flux", "flux", "peak_adu") if c in ddf.columns),
                    None,
                )
                flux_all = (
                    pd.to_numeric(ddf[flux_col], errors="coerce").to_numpy(float)
                    if flux_col else np.ones(len(ddf))
                )
                nn_all = (
                    pd.to_numeric(ddf["nearest_neighbor_fwhm"], errors="coerce").to_numpy(float)
                    if "nearest_neighbor_fwhm" in ddf.columns else None
                )

                # Anchor set: Step 4 anchor_candidate, else a relaxed
                # nearest-neighbour isolation cut. Flux is isolation-weighted so
                # the most isolated stars dominate quad construction.
                anchor_mask = None
                if "anchor_candidate" in ddf.columns:
                    anchor_mask = (
                        ddf["anchor_candidate"].astype(str).str.strip().str.lower()
                        .isin({"1", "true", "t", "yes", "y"}).to_numpy(bool)
                    )
                elif nn_all is not None:
                    anchor_mask = nn_all >= 3.0
                if anchor_mask is not None:
                    am = anchor_mask & valid_all
                    if int(am.sum()) >= ANCHOR_MIN:
                        a_xy = xy_all[am]
                        a_fl = flux_all[am]
                        if nn_all is not None:
                            nn = nn_all[am]
                            nn = np.where(np.isfinite(nn) & (nn > 0), nn, 1.0)
                            a_fl = a_fl * np.clip(nn / np.nanmedian(nn), 0.3, 3.0)
                        src_sets.append(("anchor", a_xy, a_fl))

                # Raw set: all valid detections, plain flux. Always available as
                # the fallback — keeps the crowded core that isolation removes.
                r_xy = xy_all[valid_all]
                r_fl = flux_all[valid_all] if flux_col else None
                if len(r_xy) >= 8:
                    # Avoid a duplicate pass when anchors == all detections.
                    if not (src_sets and len(src_sets[0][1]) == len(r_xy)):
                        src_sets.append(("raw", r_xy, r_fl))
                break
            except Exception:
                pass

        if not src_sets:
            src_xy_fb, src_flux_fb = self._detect_sources_from_fits(fits_path)
            if src_xy_fb is not None and len(src_xy_fb) >= 8:
                src_sets.append(("fitsdetect", src_xy_fb, src_flux_fb))
                det_xy_for_qc = src_xy_fb
                self._log(
                    f"[Internal WCS] Step 4 detections missing/sparse; "
                    f"fallback detector found {len(src_xy_fb)} peaks for {fname}"
                )
            else:
                self._log(f"[Internal WCS] No usable detections for {fname}")
                return {"ok": False, "reason": "no_detections",
                        "elapsed_s": time.time() - _t_frame_start}

        # Initial binding for the closures below; the source ladder rebinds
        # these per attempt.
        src_xy = src_sets[0][1]
        src_flux = src_sets[0][2]
        if len(det_xy_for_qc) == 0:
            det_xy_for_qc = np.asarray(src_xy, dtype=float)

        # Solve. Crowded fields are sensitive to how many Gaia stars enter quad
        # construction: too many bright off-frame/cluster-core stars can swamp
        # the correct pattern. Try an AstralImage-style catalogue-size ladder
        # before falling back to the wider local-blind search.
        solve_params = dict(self.advanced_params)
        local_blind_fallback = bool(solve_params.pop("local_blind_fallback", True))
        solve_params["local_blind"] = bool(solve_params.get("local_blind", False))
        try:
            def _cat_ladder_values(base_params: dict) -> list[int]:
                configured = int(base_params.get("n_brightest_cat", 900) or 900)
                src_cap = int(base_params.get("n_brightest_src", len(src_xy)) or len(src_xy))
                src_count = max(int(self.min_matches), min(src_cap, len(src_xy)))
                cat_limit = max(int(self.min_matches), len(gaia_ra))
                raw_values = [
                    src_count,
                    int(round(src_count * 1.5)),
                    configured,
                    int(round(src_count * 2.5)),
                    int(round(src_count * 4.0)),
                ]
                out: list[int] = []
                for value in raw_values:
                    value = max(int(self.min_matches), min(int(value), cat_limit))
                    if value not in out:
                        out.append(value)
                return out

            def _acceptable_result(res) -> bool:
                if not getattr(res, "converged", False):
                    return False
                rms = float(getattr(res, "rms_arcsec", np.nan))
                return (not np.isfinite(rms)) or rms <= float(self.rms_max_arcsec)

            def _run_solver_ladder(base_params: dict, *, local_blind: bool,
                                   cra: float, cdec: float):
                ladder = _cat_ladder_values(base_params)
                self._log(
                    f"[Internal WCS] Cat ladder for {fname}"
                    f"{' (local-blind)' if local_blind else ''}: {ladder}"
                )
                best = None
                for n_cat in ladder:
                    attempt_params = dict(base_params)
                    attempt_params["local_blind"] = bool(local_blind)
                    attempt_params["n_brightest_cat"] = int(n_cat)
                    result_i = solve_wcs_internal(
                        src_xy,
                        src_flux,
                        gaia_ra, gaia_dec, gaia_g,
                        cra, cdec, approx_scale,
                        (ny, nx),
                        min_matches=self.min_matches,
                        log_fn=self._log,
                        gaia_pmra=gaia_pmra,
                        gaia_pmdec=gaia_pmdec,
                        obs_epoch_jyear=obs_epoch_jyear,
                        **attempt_params,
                    )
                    if _acceptable_result(result_i):
                        return result_i
                    if getattr(result_i, "converged", False):
                        if best is None:
                            best = result_i
                        else:
                            best_rms = float(getattr(best, "rms_arcsec", np.inf))
                            this_rms = float(getattr(result_i, "rms_arcsec", np.inf))
                            if np.isfinite(this_rms) and this_rms < best_rms:
                                best = result_i
                    elif best is None:
                        best = result_i
                return best

            # ── Hint chain: header pointing → Step 1 target → blind ───────
            # Each candidate only changes the seed/projection CENTER fed to the
            # solver (the solver is translation-invariant, so a wrong center
            # just selects the wrong catalogue subset). The first candidate
            # that produces an acceptable solution wins, and its label is
            # recorded as provenance. Runs against the CURRENT src_xy/src_flux
            # binding so the outer source ladder can re-run it with a different
            # detection set.
            def _solve_one_source():
                candidates = _astnet_center_candidates(header_coord, target_coord)
                if header_coord is not None and target_coord is not None:
                    sep = _coord_sep_deg(header_coord, target_coord)
                    if np.isfinite(sep) and sep > 5.0:
                        self._log(
                            f"[Internal WCS] {fname}: header pointing is {sep:.2f}° "
                            f"from Step 1 target — possible stale mount header."
                        )

                def _candidate_center(coord):
                    if coord is not None:
                        try:
                            return float(coord.ra.deg), float(coord.dec.deg)
                        except Exception:
                            return None
                    if np.isfinite(approx_ra) and np.isfinite(approx_dec):
                        return float(approx_ra), float(approx_dec)
                    if header_coord is not None:
                        try:
                            return float(header_coord.ra.deg), float(header_coord.dec.deg)
                        except Exception:
                            return None
                    return None

                res_l = None
                hint_l = "none"
                best_l = None
                best_l_hint = "none"
                for label, coord in candidates:
                    center = _candidate_center(coord)
                    if center is None:
                        continue
                    cra, cdec = center
                    use_blind = (label == "blind") or bool(solve_params.get("local_blind", False))
                    cand_res = _run_solver_ladder(
                        solve_params, local_blind=use_blind, cra=cra, cdec=cdec,
                    )
                    if cand_res is None:
                        continue
                    if _acceptable_result(cand_res):
                        return cand_res, label
                    if best_l is None or (
                        getattr(cand_res, "converged", False)
                        and not getattr(best_l, "converged", False)
                    ):
                        best_l = cand_res
                        best_l_hint = label

                if (local_blind_fallback
                        and not bool(solve_params.get("local_blind", False))):
                    center = _candidate_center(None)
                    if center is not None:
                        self._log(f"[Internal WCS] Local-blind retry for {fname}")
                        retry_result = _run_solver_ladder(
                            solve_params, local_blind=True, cra=center[0], cdec=center[1],
                        )
                        if retry_result is not None and _acceptable_result(retry_result):
                            return retry_result, "blind"
                        if retry_result is not None and best_l is None:
                            best_l = retry_result
                            best_l_hint = "blind"
                return best_l, best_l_hint

            # ── Source ladder: anchor (clean) → raw (robust) ──────────────
            result = None
            hint_source = "none"
            src_set_used = "none"
            overall_best = None
            overall_best_hint = "none"
            overall_best_src = "none"
            for _si, (_src_label, _sxy, _sfl) in enumerate(src_sets):
                src_xy = _sxy
                src_flux = _sfl
                if len(src_xy) < self.min_matches:
                    continue
                if _si > 0:
                    self._log(
                        f"[Internal WCS] {fname}: '{src_sets[0][0]}' set unsolved — "
                        f"retrying with '{_src_label}' source set ({len(src_xy)} stars)"
                    )
                r_i, h_i = _solve_one_source()
                if r_i is not None and _acceptable_result(r_i):
                    result, hint_source, src_set_used = r_i, h_i, _src_label
                    break
                if r_i is not None and (
                    overall_best is None
                    or (getattr(r_i, "converged", False)
                        and not getattr(overall_best, "converged", False))
                ):
                    overall_best, overall_best_hint, overall_best_src = r_i, h_i, _src_label

            if result is None:
                result, hint_source, src_set_used = overall_best, overall_best_hint, overall_best_src
            if result is None:
                # Synthesize a non-converged result so downstream handling is
                # uniform (no candidate produced anything at all).
                from apex.analysis.astrometry import WCSSolveResult as _WSR
                result = _WSR(
                    converged=False, wcs=None, n_matches=0, rms_arcsec=float("nan"),
                    rotation_deg=float("nan"), scale_arcsec_per_px=float("nan"),
                    flip_x=False, flip_y=False,
                )
        except Exception as exc:
            self._log(f"[Internal WCS] Solver error {fname}: {exc}")
            return {"ok": False, "reason": f"solver_exception: {exc}",
                    "elapsed_s": time.time() - _t_frame_start}

        if not result.converged:
            reason_line = ""
            try:
                for entry in reversed(result.log or []):
                    s = str(entry).strip()
                    low = s.lower()
                    if any(tok in low for tok in (
                        "ransac", "match_too_small", "no candidates",
                        "no quad", "code matches", "rejected",
                        "fit failed", "too few",
                    )):
                        reason_line = s
                        break
                if not reason_line and result.log:
                    reason_line = str(result.log[-1]).strip()
            except Exception:
                reason_line = ""
            reason_short = reason_line[:200] if reason_line else "no_solution"
            self._log(f"[Internal WCS] No solution: {fname}  ({reason_short})")
            return {"ok": False, "reason": reason_short, "hint_source": hint_source,
                    "elapsed_s": time.time() - _t_frame_start}

        if np.isfinite(result.rms_arcsec) and result.rms_arcsec > self.rms_max_arcsec:
            self._log(
                f"[Internal WCS] RMS={result.rms_arcsec:.3f}\" > limit={self.rms_max_arcsec}\" → rejected: {fname}"
            )
            return {"ok": False, "reason": f"rms_exceeded_{result.rms_arcsec:.2f}",
                    "hint_source": hint_source,
                    "elapsed_s": time.time() - _t_frame_start}

        center_ra_deg = float("nan")
        center_dec_deg = float("nan")
        qc_metrics = _shared_empty_wcs_qc_metrics(n_detect=len(det_xy_for_qc))
        qc_pass = False
        qc_reason = "not_evaluated"
        try:
            wcs_hdr = result.wcs.to_header(relax=True)
            try:
                cx = (nx - 1) / 2.0
                cy = (ny - 1) / 2.0
                sky_c = result.wcs.pixel_to_world(cx, cy)
                center_ra_deg = float(sky_c.ra.deg)
                center_dec_deg = float(sky_c.dec.deg)
            except Exception:
                pass
            try:
                target_coord = SkyCoord(approx_ra * u.deg, approx_dec * u.deg, frame="icrs")
                qc_metrics = _shared_compute_wcs_qc_metrics(
                    self.params,
                    w=result.wcs,
                    det_xy=det_xy_for_qc,
                    nx=nx,
                    ny=ny,
                    gaia_ra_deg=gaia_ra,
                    gaia_dec_deg=gaia_dec,
                    pix_input_arcsec=float(approx_scale),
                    pix_fit_arcsec=float(result.scale_arcsec_per_px),
                    center_coord=target_coord,
                )
                qc_pass, qc_reasons = _shared_evaluate_wcs_qc_pass(
                    self.params, qc_metrics, wcs_ok=True,
                )
                qc_reason = ",".join(qc_reasons)
            except Exception as exc:
                qc_reason = f"qc_error:{exc}"
            with fits.open(fits_path, mode="update", memmap=False) as hdul:
                hdr_out = hdul[0].header
                _reset_apex_wcs_meta(hdr_out)
                for key, val in wcs_hdr.items():
                    try:
                        hdr_out[key] = val
                    except Exception:
                        pass
                hdr_out["WCSSRC"] = ("APEX_INTERNAL", "WCS source (apex)")
                hdr_out["WCS_OK"] = (True, "WCS solve success (APEX Internal)")
                hdr_out["WCSPIXI"] = (float(approx_scale), "pixscale input (arcsec/pix)")
                if np.isfinite(result.scale_arcsec_per_px):
                    hdr_out["WCSPIXF"] = (float(result.scale_arcsec_per_px), "pixscale fit (arcsec/pix)")
                if np.isfinite(result.rms_arcsec):
                    hdr_out["WCSRMD"] = (float(result.rms_arcsec), "ref resid med (arcsec)")
                if np.isfinite(result.rotation_deg):
                    hdr_out["WCSROT"] = (float(result.rotation_deg), "WCS rotation (deg, E of N)")
                if np.isfinite(center_ra_deg):
                    hdr_out["WCSCRA"] = (float(center_ra_deg), "WCS center RA (deg)")
                if np.isfinite(center_dec_deg):
                    hdr_out["WCSCDEC"] = (float(center_dec_deg), "WCS center Dec (deg)")
                if result.n_matches > 0:
                    hdr_out["WCSNST"] = (int(result.n_matches), "WCS matched stars")
                if getattr(result, "sip_order", 0):
                    hdr_out["WCSSIP"] = (int(result.sip_order), "SIP distortion order")
                hdr_out["WCSMOD"] = (str(getattr(result, "model", "TAN"))[:16], "WCS model")
                hdul.flush()
        except Exception as exc:
            self._log(f"[Internal WCS] Header write error {fname}: {exc}")
            return {"ok": False, "reason": f"header_write: {exc}",
                    "elapsed_s": time.time() - _t_frame_start}

        self._log(
            f"[Internal WCS] ✓ {fname}  hint={hint_source}  "
            f"matches={result.n_matches}  RMS={result.rms_arcsec:.3f}\"  "
            f"rot={result.rotation_deg:.1f}°  "
            f"model={getattr(result, 'model', 'TAN')}  "
            f"wcs_qc={'PASS' if qc_pass else 'FAIL'}"
        )
        info = {
            "ok": True,
            "wcs_ok": True,
            "status": "ok",
            "hint_source": hint_source,
            "src_set": src_set_used,
            "n_matches": result.n_matches,
            "rms_arcsec": result.rms_arcsec,
            "resid_med": result.rms_arcsec,
            "rotation_deg": result.rotation_deg,
            "flip_x": result.flip_x,
            "flip_y": result.flip_y,
            "model": getattr(result, "model", "TAN"),
            "sip_order": int(getattr(result, "sip_order", 0) or 0),
            "scale_arcsec_per_px": result.scale_arcsec_per_px,
            "pix_fit": result.scale_arcsec_per_px,
            "refine": f"m={int(result.n_matches)}",
            "match_n": int(qc_metrics.get("n_match", 0) or 0),
            "center_ra_deg": center_ra_deg,
            "center_dec_deg": center_dec_deg,
            "solver": "apex_internal",
            "wcssrc": "APEX_INTERNAL",
            "wcs_qc_pass": bool(qc_pass),
            "wcs_qc_reason": qc_reason,
            "elapsed_s": time.time() - _t_frame_start,
        }
        info.update(qc_metrics)
        self._write_internal_sidecar(fname, info, approx_scale)
        return info

    def _write_internal_sidecar(self, fname: str, info: dict, pix_input: float) -> None:
        """Write per-frame WCS metadata sidecar to cache_dir/wcs_solve/.

        Mirrors the ASTAP/AstNet path so downstream steps (Step 6 QC,
        diagnostics) can read solver provenance & residuals uniformly.
        """
        try:
            meta_dir = Path(self.params.P.cache_dir) / "wcs_solve"
            meta_dir.mkdir(parents=True, exist_ok=True)
            payload = {
                "fname": fname,
                "ok": bool(info.get("ok", False)),
                "wcs_ok": bool(info.get("ok", False)),
                "solver": "apex_internal",
                "wcssrc": info.get("wcssrc", "APEX_INTERNAL"),
                "model": str(info.get("model", "TAN") or "TAN"),
                "sip_order": int(info.get("sip_order", 0) or 0),
                "pix_input_arcsec": float(pix_input) if np.isfinite(pix_input) else None,
                "pix_fit_arcsec": (float(info["scale_arcsec_per_px"])
                                    if np.isfinite(info.get("scale_arcsec_per_px", np.nan)) else None),
                "rms_arcsec": (float(info["rms_arcsec"])
                                if np.isfinite(info.get("rms_arcsec", np.nan)) else None),
                "rotation_deg": (float(info["rotation_deg"])
                                  if np.isfinite(info.get("rotation_deg", np.nan)) else None),
                "center_ra_deg": (float(info["center_ra_deg"])
                                   if np.isfinite(info.get("center_ra_deg", np.nan)) else None),
                "center_dec_deg": (float(info["center_dec_deg"])
                                    if np.isfinite(info.get("center_dec_deg", np.nan)) else None),
                "n_matches": int(info.get("n_matches", 0) or 0),
                "flip_x": bool(info.get("flip_x", False)),
                "flip_y": bool(info.get("flip_y", False)),
                "wcs_qc_pass": bool(info.get("wcs_qc_pass", False)),
                "wcs_qc_reason": str(info.get("wcs_qc_reason", "") or ""),
            }
            for key in (
                "n_detect", "n_catalog_in_fov", "n_match", "n_inlier",
                "match_rate", "match_rate_cat", "match_rate_eff",
                "match_radius_arcsec", "match_radius_px",
                "resid_med_px", "resid_mad_px", "resid_p99_px",
                "rms_px", "inlier_rate", "center_offset_arcsec",
                "scale_delta_pct", "gaia_available",
            ):
                if key in info:
                    value = info.get(key)
                    try:
                        if isinstance(value, (np.integer,)):
                            value = int(value)
                        elif isinstance(value, (np.floating, float)):
                            value = float(value) if np.isfinite(value) else None
                    except Exception:
                        pass
                    payload[key] = value
            (meta_dir / f"wcs_{fname}.json").write_text(
                json.dumps(payload, indent=2), encoding="utf-8"
            )
        except Exception as exc:
            self._log(f"[Internal WCS] sidecar write skipped for {fname}: {exc}")


class WcsPlateSolvingWindow(StepWindowBase):
    """Step 5: WCS Plate Solving"""

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.worker = None
        self.results = {}
        self.stop_requested = False
        self.log_window = None
        self._wcs_cache_mgr: StepCacheManager | None = None

        self.file_list = []
        self.use_cropped = False

        super().__init__(
            step_index=4,
            step_name="WCS Plate Solving",
            params=params,
            project_state=project_state,
            main_window=main_window
        )

        self.setup_step_ui()
        self.restore_state()

    def setup_step_ui(self):
        # ── Persistent "which solver ran each file?" banner ─────────────────
        # The three tabs each only render their own solver's runs, so after a
        # project re-open it is otherwise invisible which engine produced the
        # current WCS for each frame. This banner reads self.results (loaded
        # from wcs_solve_summary.csv) and prints a per-solver tally above the
        # tabs so the answer is one glance away, regardless of which tab is
        # active.
        self.solver_breakdown_label = QLabel("")
        self.solver_breakdown_label.setWordWrap(True)
        self.solver_breakdown_label.setStyleSheet(
            "QLabel { background-color: #ECEFF1; padding: 6px 10px;"
            " border-radius: 4px; color: #37474F; }"
        )
        self.content_layout.addWidget(self.solver_breakdown_label)

        # Create tab widget
        self.tab_widget = QTabWidget()
        self.content_layout.addWidget(self.tab_widget)

        # Internal (Python) Tab — first, no external tools required
        self.internal_tab = QWidget()
        self.setup_internal_tab()
        self.tab_widget.addTab(self.internal_tab, "Internal (Python)")

        # ASTAP Tab
        self.astap_tab = QWidget()
        self.setup_astap_tab()
        self.tab_widget.addTab(self.astap_tab, "ASTAP (Local)")

        # Astrometry.net Tab
        self.astrometrynet_tab = QWidget()
        self.setup_astrometrynet_tab()
        self.tab_widget.addTab(self.astrometrynet_tab, "Astrometry.net (Local)")

        self.setup_log_window()
        self.populate_file_list()

    def setup_internal_tab(self):
        """Internal (Python) WCS solver — no external executables needed."""
        layout = QVBoxLayout(self.internal_tab)

        info = QLabel(
            "Solve WCS using the built-in astnet-style quad-hash matcher and "
            "Gaia refinement (numpy/scipy/astropy only — no external binary). "
            "Needs Step 1 target RA/Dec and the Instrument pixel scale."
        )
        info.setWordWrap(True)
        info.setStyleSheet("QLabel { background-color: #E8F5E9; padding: 10px; border-radius: 5px; }")
        layout.addWidget(info)

        # All Internal-solver parameters live in one dict edited by the
        # parameter dialog — same pattern as ASTAP/astnet so the three tabs
        # behave consistently.
        self._internal_params = {
            "min_matches": 8,
            "rms_max_arcsec": 2.0,
            "n_brightest_src": 350,
            "n_brightest_cat": 900,
            "quad_k_neighbor": 10,
            "quad_neighbor_pool_factor": 3,
            "quad_code_tol": 0.020,
            "quad_scale_ratio_tol": 0.25,
            "quad_max_per_side": 3000,
            "ransac_inlier_radius_px": 4.0,
            "ransac_max_trials": 4000,
            "ransac_keep_candidates": 8,
            "allow_reflection": True,
            "local_blind_fallback": True,
            "local_blind_radius_factor": 2.5,
            "sip_degree": 3,
            "sip_min_pairs": 30,
            "sip_holdout_fraction": 0.25,
            "sip_min_improvement": 0.10,
        }

        control_layout = QHBoxLayout()
        btn_int_params = create_parameter_button("Internal Parameters")
        btn_int_params.clicked.connect(self._open_internal_parameters_dialog)
        control_layout.addWidget(btn_int_params)

        self.run_bar_internal = RunControlBar(
            "Run Internal Solver", "Log",
            run_cb=self.run_wcs_internal_solver,
            stop_cb=self.stop_wcs_internal_solver,
            log_cb=self.show_log_window,
        )
        control_layout.addWidget(self.run_bar_internal)
        layout.addLayout(control_layout)

        progress_layout = QHBoxLayout()
        self.internal_progress = QProgressBar()
        self.internal_progress.setMinimum(0)
        self.internal_progress.setMaximum(100)
        self.internal_progress.setValue(0)
        progress_layout.addWidget(self.internal_progress)

        self.internal_status = QLabel("Ready")
        self.internal_status.setMinimumWidth(350)
        progress_layout.addWidget(self.internal_status)
        layout.addLayout(progress_layout)

        results_group = QGroupBox("Internal Solver Results")
        results_layout = QVBoxLayout(results_group)

        self.internal_results_table = QTableWidget()
        self.internal_results_table.setColumnCount(9)
        self.internal_results_table.setHorizontalHeaderLabels([
            "File", "Status", "Hint", "RA", "Dec", "PixScale", "Refine", "Resid(\")", "Elapsed (s)"
        ])
        self.internal_results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.internal_results_table.horizontalHeader().setStretchLastSection(True)
        self.internal_results_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.internal_results_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        results_layout.addWidget(self.internal_results_table)

        layout.addWidget(results_group)

    # ── Shared run front-matter (all three solver tabs) ──────────────────
    def _qc_filter_and_log(self, file_list=None):
        """Apply the Step 4 QC frame filter and log the outcome.

        Shared by the ASTAP / astrometry.net / Internal tabs so the QC line is
        computed and worded identically everywhere. Returns the filtered file
        list, or ``None`` when nothing remains (the caller should abort — a
        warning dialog has already been shown).
        """
        src = list(self.file_list if file_list is None else file_list)
        require_qc = bool(getattr(self.params.P, "wcs_require_qc_pass", True))
        filtered, qc_info = filter_files_by_qc(
            Path(self.params.P.result_dir), src, require_qc=require_qc,
        )
        if require_qc:
            if qc_info.get("applied"):
                self.log(f"[WCS][QC] Frame QC filter: {qc_info['kept']}/{qc_info['total']} kept.")
            elif qc_info.get("path") is None:
                self.log("[WCS][QC] frame_quality.csv not found; using all frames.")
            else:
                self.log(f"[WCS][QC] frame_quality.csv ignored ({qc_info['reason']}); using all frames.")
        if not filtered:
            QMessageBox.warning(self, "Warning", "No frames remain after Step 4 QC filtering.")
            return None
        return filtered

    def _log_preflight(self, *, astap=False, astnet=False, gaia=False):
        """Log a uniform ``[WCS][Preflight]`` header and return availability.

        Each tab requests only the backends it actually uses, so no slow probe
        runs unnecessarily (e.g. the Internal tab never spins up the WSL
        astrometry.net check). The wording/format is shared so the preflight
        header reads the same regardless of which tab launched the run.

        Returns a dict like ``{"astap": (ok, detail), ...}`` for requested
        backends only.
        """
        out: dict[str, tuple[bool, str]] = {}
        if astap:
            ok, detail = _check_astap_available(self.params)
            self.log(f"[WCS][Preflight] ASTAP: {detail}")
            out["astap"] = (ok, detail)
        if astnet:
            ok, detail = _check_astnet_available(self.params)
            self.log(f"[WCS][Preflight] astrometry.net: {detail}")
            out["astnet"] = (ok, detail)
        if gaia:
            ok, detail = _check_gaia_runtime_available()
            self.log(f"[WCS][Preflight] Gaia: {detail}")
            out["gaia"] = (ok, detail)
        return out

    def run_wcs_internal_solver(self):
        if not self.file_list:
            QMessageBox.warning(self, "Warning", "No files to process")
            return
        if hasattr(self, "_internal_worker") and self._internal_worker and \
                self._internal_worker.isRunning():
            return

        self.log_text.clear()
        file_list = self._qc_filter_and_log()
        if file_list is None:
            return

        # Internal can still solve from a cached gaia_fov.ecsv even when the
        # runtime probe is unavailable, so this is informational — no early
        # return.
        self._log_preflight(gaia=True)

        self.internal_results_table.setRowCount(0)
        # Drop any prior run's rows so a fresh Internal run starts from a
        # clean slate. ASTAP / astnet handlers reset implicitly through their
        # own paths; we do it explicitly here so cache writers below see only
        # the current solver's output and never accidentally pick up stale
        # ASTAP/astnet rows that share filenames.
        self.results = {}
        self._internal_run_results = {}

        p = dict(self._internal_params)
        min_matches = int(p.pop("min_matches"))
        rms_max = float(p.pop("rms_max_arcsec"))

        use_cropped = bool(getattr(self, "use_cropped", False))
        self._internal_worker = InternalWcsWorker(
            file_list,
            self.params,
            self.params.P.data_dir,
            self.params.P.result_dir,
            use_cropped=use_cropped,
            min_matches=min_matches,
            rms_max_arcsec=rms_max,
            advanced_params=p,
        )
        self._internal_worker.progress.connect(self._on_internal_progress)
        self._internal_worker.log_message.connect(self.log)
        self._internal_worker.frame_done.connect(self._on_internal_frame_done)
        self._internal_worker.finished.connect(self._on_internal_finished)
        self._internal_worker.error.connect(self._on_internal_error)

        # Wire the shared Worker Status panel so the Internal tab shows one row
        # per parallel worker, just like the ASTAP / astnet tabs.
        self.setup_log_window()
        if hasattr(self, "_worker_panel") and self._worker_panel is not None:
            self._worker_panel.clear()
            self._internal_worker.worker_status.connect(self._worker_panel.update_worker)

        self.run_bar_internal.set_running(True)
        self.internal_progress.setValue(0)
        self.internal_progress.setMaximum(len(file_list))
        self.internal_status.setText("Solving…")
        self._internal_worker.start()

    def _on_internal_progress(self, i: int, total: int, fname: str):
        self.internal_progress.setMaximum(max(total, 1))
        self.internal_progress.setValue(int(i))
        self.internal_status.setText(f"{i}/{total}  {fname}")

    def stop_wcs_internal_solver(self):
        if hasattr(self, "_internal_worker") and self._internal_worker:
            self._internal_worker.stop()

    def _open_internal_parameters_dialog(self):
        """Edit Internal Solver parameters in a single dialog.

        Same styling pattern as ASTAP / Astrometry.net parameter dialogs
        (configure_parameter_dialog + grouped form).
        """
        from PyQt5.QtWidgets import (
            QDialog, QFormLayout, QVBoxLayout, QSpinBox, QDoubleSpinBox,
            QCheckBox, QDialogButtonBox, QGroupBox,
        )

        d = QDialog(self)
        configure_parameter_dialog(d, "Internal Parameters", 560, 780)
        root = QVBoxLayout(d)
        p = self._internal_params

        basic = QGroupBox("Basic")
        bf = QFormLayout(basic)
        sp_min = QSpinBox(); sp_min.setRange(4, 100); sp_min.setValue(int(p["min_matches"]))
        sp_min.setToolTip("Reject solutions with fewer than this many matched stars")
        bf.addRow("min_matches:", sp_min)
        sp_rms = QDoubleSpinBox(); sp_rms.setRange(0.1, 10.0); sp_rms.setDecimals(2)
        sp_rms.setSingleStep(0.1); sp_rms.setSuffix(" \"")
        sp_rms.setValue(float(p["rms_max_arcsec"]))
        sp_rms.setToolTip("Reject solutions with residual RMS above this")
        bf.addRow("rms_max:", sp_rms)
        root.addWidget(basic)

        quad = QGroupBox("Quad matching")
        qf = QFormLayout(quad)
        sp_src = QSpinBox(); sp_src.setRange(20, 2000); sp_src.setValue(int(p["n_brightest_src"]))
        sp_src.setToolTip("Brightest detected sources fed into quad construction")
        qf.addRow("n_brightest_src:", sp_src)
        sp_cat = QSpinBox(); sp_cat.setRange(20, 5000); sp_cat.setValue(int(p["n_brightest_cat"]))
        sp_cat.setToolTip("Brightest in-field Gaia stars fed into quad construction. "
                          "Crowded fields may need ≥500 so that halo stars enter the pool.")
        qf.addRow("n_brightest_cat:", sp_cat)
        sp_k = QSpinBox(); sp_k.setRange(4, 30); sp_k.setValue(int(p["quad_k_neighbor"]))
        sp_k.setToolTip("Neighbour sample size per star for quad construction")
        qf.addRow("quad_k_neighbor:", sp_k)
        sp_pool = QSpinBox(); sp_pool.setRange(1, 8); sp_pool.setValue(int(p["quad_neighbor_pool_factor"]))
        sp_pool.setToolTip("Sample neighbours from a wider pool; higher values add long baselines in crowded fields")
        qf.addRow("quad_neighbor_pool_factor:", sp_pool)
        sp_tol = QDoubleSpinBox(); sp_tol.setRange(0.001, 0.20); sp_tol.setDecimals(4)
        sp_tol.setSingleStep(0.005); sp_tol.setValue(float(p["quad_code_tol"]))
        sp_tol.setToolTip("L2 tolerance in 4-D code space for two quads to be a match")
        qf.addRow("quad_code_tol:", sp_tol)
        sp_sr = QDoubleSpinBox(); sp_sr.setRange(0.02, 1.0); sp_sr.setDecimals(3)
        sp_sr.setSingleStep(0.05); sp_sr.setValue(float(p["quad_scale_ratio_tol"]))
        sp_sr.setToolTip("Allowed |log(src_side / cat_side)| between matched quads")
        qf.addRow("quad_scale_ratio_tol:", sp_sr)
        sp_mq = QSpinBox(); sp_mq.setRange(100, 10000); sp_mq.setValue(int(p["quad_max_per_side"]))
        sp_mq.setSingleStep(100)
        sp_mq.setToolTip("Cap on quads kept per side (larger = slower but more thorough)")
        qf.addRow("quad_max_per_side:", sp_mq)
        root.addWidget(quad)

        ransac = QGroupBox("RANSAC verification")
        rf = QFormLayout(ransac)
        sp_ir = QDoubleSpinBox(); sp_ir.setRange(0.5, 20.0); sp_ir.setDecimals(2)
        sp_ir.setSingleStep(0.5); sp_ir.setValue(float(p["ransac_inlier_radius_px"]))
        sp_ir.setSuffix(" px"); sp_ir.setToolTip("Inlier radius for the candidate similarity")
        rf.addRow("ransac_inlier_radius_px:", sp_ir)
        sp_mt = QSpinBox(); sp_mt.setRange(50, 20000); sp_mt.setValue(int(p["ransac_max_trials"]))
        sp_mt.setSingleStep(100)
        sp_mt.setToolTip("Maximum quad-pair candidates to test (runtime cap)")
        rf.addRow("ransac_max_trials:", sp_mt)
        sp_keep = QSpinBox(); sp_keep.setRange(1, 20); sp_keep.setValue(int(p["ransac_keep_candidates"]))
        sp_keep.setToolTip("Distinct RANSAC candidates carried into final WCS fitting")
        rf.addRow("ransac_keep_candidates:", sp_keep)
        cb_refl = QCheckBox("Allow mirrored transform (image/sky parity flip)")
        cb_refl.setChecked(bool(p["allow_reflection"]))
        rf.addRow(cb_refl)
        root.addWidget(ransac)

        fallback = QGroupBox("Fallback")
        lf = QFormLayout(fallback)
        cb_lb = QCheckBox("Retry local-blind when normal solve fails")
        cb_lb.setChecked(bool(p.get("local_blind_fallback", True)))
        cb_lb.setToolTip(
            "Retry once with a wider Gaia tangent-plane window. This is not "
            "all-sky blind solving; Step 1 target RA/Dec and pixel scale are "
            "still required."
        )
        lf.addRow(cb_lb)
        sp_lb = QDoubleSpinBox(); sp_lb.setRange(1.0, 8.0); sp_lb.setDecimals(2)
        sp_lb.setSingleStep(0.25); sp_lb.setValue(float(p.get("local_blind_radius_factor", 2.5)))
        sp_lb.setToolTip("Gaia window radius in half-image-diagonal units for local-blind retry")
        lf.addRow("local_blind_radius_factor:", sp_lb)
        root.addWidget(fallback)

        refine = QGroupBox("Refine")
        ff = QFormLayout(refine)
        sp_sip = QSpinBox(); sp_sip.setRange(0, 5); sp_sip.setValue(int(p["sip_degree"]))
        sp_sip.setToolTip("SIP polynomial degree for distortion correction. "
                          "0 = pure TAN. 2–3 typical. Needs ≥30 matched pairs to engage.")
        ff.addRow("sip_degree:", sp_sip)
        sp_sip_min = QSpinBox(); sp_sip_min.setRange(12, 300)
        sp_sip_min.setValue(int(p.get("sip_min_pairs", 30)))
        sp_sip_min.setToolTip("Minimum matched pairs required before SIP validation is attempted")
        ff.addRow("sip_min_pairs:", sp_sip_min)
        sp_sip_hold = QDoubleSpinBox(); sp_sip_hold.setRange(0.05, 0.50)
        sp_sip_hold.setDecimals(2); sp_sip_hold.setSingleStep(0.05)
        sp_sip_hold.setValue(float(p.get("sip_holdout_fraction", 0.25)))
        sp_sip_hold.setToolTip("Fraction of matched pairs held out for SIP-vs-TAN validation")
        ff.addRow("sip_holdout_fraction:", sp_sip_hold)
        sp_sip_imp = QDoubleSpinBox(); sp_sip_imp.setRange(0.00, 0.50)
        sp_sip_imp.setDecimals(2); sp_sip_imp.setSingleStep(0.05)
        sp_sip_imp.setValue(float(p.get("sip_min_improvement", 0.10)))
        sp_sip_imp.setToolTip("Minimum holdout RMS improvement required to adopt SIP")
        ff.addRow("sip_min_improvement:", sp_sip_imp)
        root.addWidget(refine)

        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        add_parameter_reset_button(
            bb,
            [
                (sp_min, 8),
                (sp_rms, 2.0),
                (sp_src, 350),
                (sp_cat, 900),
                (sp_k, 10),
                (sp_pool, 3),
                (sp_tol, 0.020),
                (sp_sr, 0.25),
                (sp_mq, 3000),
                (sp_ir, 4.0),
                (sp_mt, 4000),
                (sp_keep, 8),
                (cb_refl, True),
                (cb_lb, True),
                (sp_lb, 2.5),
                (sp_sip, 3),
                (sp_sip_min, 30),
                (sp_sip_hold, 0.25),
                (sp_sip_imp, 0.10),
            ],
        )
        bb.accepted.connect(d.accept)
        bb.rejected.connect(d.reject)
        root.addWidget(bb)

        if d.exec_() != QDialog.Accepted:
            return

        self._internal_params.update({
            "min_matches": int(sp_min.value()),
            "rms_max_arcsec": float(sp_rms.value()),
            "n_brightest_src": int(sp_src.value()),
            "n_brightest_cat": int(sp_cat.value()),
            "quad_k_neighbor": int(sp_k.value()),
            "quad_neighbor_pool_factor": int(sp_pool.value()),
            "quad_code_tol": float(sp_tol.value()),
            "quad_scale_ratio_tol": float(sp_sr.value()),
            "quad_max_per_side": int(sp_mq.value()),
            "ransac_inlier_radius_px": float(sp_ir.value()),
            "ransac_max_trials": int(sp_mt.value()),
            "ransac_keep_candidates": int(sp_keep.value()),
            "allow_reflection": bool(cb_refl.isChecked()),
            "local_blind_fallback": bool(cb_lb.isChecked()),
            "local_blind_radius_factor": float(sp_lb.value()),
            "sip_degree": int(sp_sip.value()),
            "sip_min_pairs": int(sp_sip_min.value()),
            "sip_holdout_fraction": float(sp_sip_hold.value()),
            "sip_min_improvement": float(sp_sip_imp.value()),
        })

    def _on_internal_frame_done(self, fname: str, info: dict):
        """Add one row to the results table as each frame completes.

        Column layout mirrors the Astrometry.net tab so the three solver tabs
        present comparable data side-by-side. Also records the result into
        ``self.results`` so that ``validate_step()`` recognises the step as
        complete and the Next Step button enables — the ASTAP and astnet
        worker callbacks already do the same.
        """
        from PyQt5.QtWidgets import QTableWidgetItem

        if not hasattr(self, "_internal_run_results"):
            self._internal_run_results = {}
        self._internal_run_results[fname] = info
        if bool(info.get("ok", False)):
            self.results[fname] = info
        else:
            self.results.pop(fname, None)
        # Cache writes (parity with ASTAP/astnet on_file_done):
        #   * StepCacheManager manifest per successful frame
        #   * Navigation refresh so the Next Step button enables as soon as
        #     at least one frame is solved.
        # Both are wrapped in broad try/except — a cache hiccup must never
        # stop the per-frame loop from displaying results.
        try:
            if bool(info.get("ok", False)):
                self._write_wcs_manifest(fname, info)
        except Exception:
            pass
        try:
            self.update_navigation_buttons()
        except Exception:
            pass
        try:
            self._refresh_solver_breakdown_label()
        except Exception:
            pass

        t = self.internal_results_table
        row = t.rowCount()
        t.insertRow(row)
        ok = bool(info.get("ok", False))

        def _fmt_num(key, fmt):
            v = info.get(key)
            if v is None or not np.isfinite(float(v)):
                return ""
            return fmt.format(float(v))

        if ok:
            ra_s = _fmt_num("center_ra_deg", "{:.4f}")
            dec_s = _fmt_num("center_dec_deg", "{:.4f}")
            pix_s = _fmt_num("scale_arcsec_per_px", "{:.4f}\"/px")
            model_s = str(info.get("model", "") or "").strip()
            refine_s = f"{model_s} m={int(info.get('n_matches', 0) or 0)}".strip()
            resid_s = _fmt_num("rms_arcsec", "{:.3f}")
            status_s = "✓ QC" if self._boolish(info.get("wcs_qc_pass")) else "✓ QC?"
        else:
            ra_s = dec_s = pix_s = refine_s = resid_s = ""
            status_s = "✗ " + str(info.get("reason", "fail"))

        elapsed_s = _fmt_num("elapsed_s", "{:.2f}")
        hint_s = str(info.get("hint_source", "") or "")

        values = [fname, status_s, hint_s, ra_s, dec_s, pix_s, refine_s, resid_s, elapsed_s]
        for col, val in enumerate(values):
            t.setItem(row, col, QTableWidgetItem(val))
        set_table_row_background(t, row, self._wcs_result_row_background(info))
        t.scrollToBottom()

    def _on_internal_finished(self, summary: dict):
        self.run_bar_internal.set_running(False)
        ok = int(summary.get("ok", 0) or 0)
        total = int(summary.get("total", 0) or 0)
        stopped = bool(summary.get("stopped", False))
        msg = f"Internal solver: {ok}/{total} solved"
        if stopped:
            msg += " (stopped)"
        self.internal_status.setText(msg)
        qc_not_eval = int(summary.get("wcs_qc_not_evaluated", 0) or 0)
        qc_text = (
            f"WCS-QC not evaluated: {qc_not_eval} (Gaia unavailable)"
            if qc_not_eval
            else f"WCS-QC pass: {summary.get('wcs_qc_pass', 0)}"
        )
        self.log(f"[Internal WCS] {msg} | {qc_text}")

        # Final cache layer: write the solver-agnostic summary CSV. Without
        # this, _restore_success_results_from_summary() on the next project
        # open finds no rows for Internal-solved frames, the Next Step button
        # stays disabled, and downstream QC tools that read this CSV miss the
        # frames entirely. Wrapped so a write failure is logged but does not
        # break the user-visible "N/N solved" completion.
        try:
            self._write_internal_summary_csv()
        except Exception as exc:
            self.log(f"[Internal WCS] summary CSV write failed: {exc}")

    def _refresh_solver_breakdown_label(self) -> None:
        """Re-tally ``self.results`` by WCS solver and update the banner.

        Source of truth for the per-frame tag is, in order:
          1. ``info["wcssrc"]``  (the FITS-header keyword written by the
             worker — survives across project re-opens because it's also in
             the summary CSV)
          2. ``info["solver"]`` (legacy / Internal-native key)

        Frames whose result row doesn't carry either key are bucketed under
        ``"unknown"`` rather than silently dropped, so a missing tag is
        visible instead of pretending everything was solved.
        """
        label = getattr(self, "solver_breakdown_label", None)
        if label is None:
            return
        results = getattr(self, "results", {}) or {}
        if not results:
            label.setText("No WCS results yet — run a solver tab below.")
            return
        counts: dict[str, int] = {}
        for info in results.values():
            if not isinstance(info, dict):
                continue
            tag = info.get("wcssrc") or info.get("solver") or "unknown"
            tag = str(tag).strip() or "unknown"
            counts[tag] = counts.get(tag, 0) + 1
        total = sum(counts.values())
        parts = [f"<b>{tag}</b>: {n}" for tag, n in sorted(counts.items())]
        label.setText(
            f"Solved frames ({total}) — " + ",&nbsp; ".join(parts)
        )

    @staticmethod
    def _finite_result_float(value, default=np.nan) -> float:
        try:
            out = float(value)
        except Exception:
            return float(default)
        return out if np.isfinite(out) else float(default)

    def _append_restored_internal_row(self, fname: str, info: dict) -> None:
        if not hasattr(self, "internal_results_table"):
            return
        t = self.internal_results_table
        row = t.rowCount()
        t.insertRow(row)
        ok = self._is_successful_wcs_result(info)

        def _fmt(keys, fmt, default=""):
            for key in keys:
                val = self._finite_result_float(info.get(key))
                if np.isfinite(val):
                    return fmt.format(val)
            return default

        if ok:
            ra_s = _fmt(("center_ra_deg", "ra"), "{:.4f}")
            dec_s = _fmt(("center_dec_deg", "dec"), "{:.4f}")
            pix_s = _fmt(("scale_arcsec_per_px", "pix_fit", "pixscale"), "{:.4f}\"/px")
            model_s = str(info.get("model", "") or "").strip()
            refine_s = str(info.get("refine", "") or "").strip()
            if not refine_s and info.get("n_matches") is not None:
                try:
                    refine_s = f"{model_s} m={int(info.get('n_matches', 0) or 0)}".strip()
                except Exception:
                    refine_s = model_s
            resid_s = _fmt(("rms_arcsec", "resid_med", "resid_med_px"), "{:.3f}")
            status_s = str(info.get("status", "") or "ok")
        else:
            ra_s = dec_s = pix_s = refine_s = resid_s = ""
            status_s = str(info.get("status", "") or info.get("reason", "") or "fail")

        elapsed_s = _fmt(("elapsed_s", "elapsed"), "{:.2f}")
        hint_s = str(info.get("hint_source", "") or "")
        values = [fname, status_s, hint_s, ra_s, dec_s, pix_s, refine_s, resid_s, elapsed_s]
        for col, val in enumerate(values):
            t.setItem(row, col, QTableWidgetItem(val))
        set_table_row_background(t, row, self._wcs_result_row_background(info))

    def _append_restored_astnet_row(self, fname: str, info: dict) -> None:
        if not hasattr(self, "astrometrynet_results_table"):
            return
        t = self.astrometrynet_results_table
        row = t.rowCount()
        t.insertRow(row)
        pixscale = self._finite_result_float(info.get("pixscale", info.get("pix_fit")))
        resid_med = self._finite_result_float(info.get("resid_med", info.get("resid_med_px")))
        elapsed = self._finite_result_float(info.get("elapsed_s", info.get("elapsed")))
        values = [
            fname,
            str(info.get("status", "") or "solved"),
            f"{pixscale:.4f}" if np.isfinite(pixscale) and pixscale > 0 else "-",
            str(info.get("refine", "") or "-"),
            f"{resid_med:.2f}" if np.isfinite(resid_med) else "-",
            f"{elapsed:.1f}" if np.isfinite(elapsed) and elapsed > 0 else "-",
        ]
        for col, val in enumerate(values):
            t.setItem(row, col, QTableWidgetItem(val))
        set_table_row_background(t, row, self._wcs_result_row_background(info))

    def _append_restored_astap_row(self, fname: str, info: dict) -> None:
        if not hasattr(self, "results_table"):
            return
        row = self.results_table.rowCount()
        self.results_table.insertRow(row)
        self.results_table.setItem(row, 0, QTableWidgetItem(fname))
        self.results_table.setItem(row, 1, QTableWidgetItem(str(info.get("status", "") or "ok")))
        pix_fit = self._finite_result_float(info.get("pix_fit", info.get("pixscale")))
        self.results_table.setItem(
            row, 2,
            QTableWidgetItem(f"{pix_fit:.4f}" if np.isfinite(pix_fit) and pix_fit > 0 else "-"),
        )
        self.results_table.setItem(row, 3, QTableWidgetItem(str(info.get("refine", "") or "-")))
        resid_med = self._finite_result_float(info.get("resid_med"))
        if np.isfinite(resid_med):
            resid_str = f"{resid_med:.3f}\""
        else:
            resid_px = self._finite_result_float(info.get("resid_med_px"))
            resid_str = f"{resid_px:.3f}px" if np.isfinite(resid_px) else "-"
        self.results_table.setItem(row, 4, QTableWidgetItem(resid_str))
        elapsed = self._finite_result_float(info.get("elapsed", info.get("elapsed_s")), 0.0)
        self.results_table.setItem(row, 5, QTableWidgetItem(f"{elapsed:.1f}" if np.isfinite(elapsed) else "0.0"))
        set_table_row_background(self.results_table, row, self._wcs_result_row_background(info))

    def _populate_restored_wcs_tables(self, restored: dict[str, dict]) -> dict[str, int]:
        counts = {"internal": 0, "astnet": 0, "astap": 0}
        for table_name in ("internal_results_table", "results_table", "astrometrynet_results_table"):
            table = getattr(self, table_name, None)
            if table is not None:
                table.setRowCount(0)

        for fname, info in restored.items():
            tag = str(info.get("wcssrc") or info.get("solver") or "").strip().lower()
            if "internal" in tag or "apex_internal" in tag:
                self._append_restored_internal_row(fname, info)
                counts["internal"] += 1
            elif "astnet" in tag or "astrometry" in tag or "solve-field" in tag:
                self._append_restored_astnet_row(fname, info)
                counts["astnet"] += 1
            else:
                self._append_restored_astap_row(fname, info)
                counts["astap"] += 1
        return counts

    def _write_internal_summary_csv(self):
        """Persist ``self.results`` to ``step5/wcs_solve_summary.csv``.

        Schema-aligned with the ASTAP and astrometry.net writers so the same
        downstream loaders (project re-open, Step 6 QC, diagnostics) work
        unchanged. The Internal solver's native result dict uses slightly
        different field names than ASTAP's, so we map the common columns:

        ============================ ==========================================
        Internal key                  Aliased to (for cross-solver consumers)
        ============================ ==========================================
        ``rms_arcsec``                ``resid_med``
        ``scale_arcsec_per_px``       ``pix_fit``
        ``elapsed_s``                 ``elapsed``
        ``n_matches``                 ``refine`` (formatted as ``m={N}``)
        ============================ ==========================================

        We never drop the original field — the alias is added alongside so
        Internal-specific consumers and the cross-solver loader both see what
        they expect. ``file`` and ``fname`` are both written for the same
        reason (``_restore_success_results_from_summary`` accepts either).
        """
        source_results = getattr(self, "_internal_run_results", None) or self.results
        if not source_results:
            return
        rows: list[dict] = []
        for fname, info in source_results.items():
            row = dict(info) if isinstance(info, dict) else {}
            row.setdefault("file", str(fname))
            row.setdefault("fname", str(fname))
            ok = bool(row.get("ok", False))
            row["ok"] = ok
            row.setdefault("status", "ok" if ok else "fail")
            row.setdefault("solver", "apex_internal")
            row.setdefault("wcssrc", "APEX_INTERNAL")
            row.setdefault("hint_source", "")
            if "rms_arcsec" in row and "resid_med" not in row:
                row["resid_med"] = row["rms_arcsec"]
            if "scale_arcsec_per_px" in row and "pix_fit" not in row:
                row["pix_fit"] = row["scale_arcsec_per_px"]
            if "elapsed_s" in row and "elapsed" not in row:
                row["elapsed"] = row["elapsed_s"]
            if "n_matches" in row and "refine" not in row:
                try:
                    row["refine"] = f"m={int(row.get('n_matches', 0) or 0)}"
                except Exception:
                    row["refine"] = ""
            rows.append(row)
        step5_out = step5_wcs_dir(self.params.P.result_dir)
        step5_out.mkdir(parents=True, exist_ok=True)
        df = pd.DataFrame(rows)
        out_path = step5_out / "wcs_solve_summary.csv"
        df.to_csv(out_path, index=False)
        self.log(f"[Internal WCS] Summary CSV: {len(df)} rows → {out_path.name}")

    def _on_internal_error(self, msg: str):
        self.run_bar_internal.set_running(False)
        self.log(f"[Internal WCS] ERROR: {msg}")
        QMessageBox.warning(self, "Internal WCS Solver Error", msg)

    def setup_astap_tab(self):
        """Setup ASTAP tab UI"""
        layout = QVBoxLayout(self.astap_tab)

        info = QLabel(
            "Solve WCS for all frames using ASTAP (local). "
            "ASTAP and its D80/D50 star database are installed separately; "
            "see Help > WCS Solver Installation Help."
        )
        info.setWordWrap(True)
        info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; border-radius: 5px; }")
        layout.addWidget(info)

        control_layout = QHBoxLayout()
        btn_params = create_parameter_button("ASTAP Parameters")
        btn_params.clicked.connect(self.open_parameters_dialog)
        control_layout.addWidget(btn_params)

        self.run_bar_astap = RunControlBar(
            "Run ASTAP", "Log",
            run_cb=self.run_wcs,
            stop_cb=self.stop_wcs,
            log_cb=self.show_log_window,
        )
        control_layout.addWidget(self.run_bar_astap)
        self.btn_run = self.run_bar_astap.btn_run
        self.btn_stop = self.run_bar_astap.btn_stop

        layout.addLayout(control_layout)

        progress_layout = QHBoxLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        progress_layout.addWidget(self.progress_bar)

        self.progress_label = QLabel("Ready")
        self.progress_label.setMinimumWidth(350)
        progress_layout.addWidget(self.progress_label)
        layout.addLayout(progress_layout)

        results_group = QGroupBox("WCS Results")
        results_layout = QVBoxLayout(results_group)

        self.results_table = QTableWidget()
        self.results_table.setColumnCount(6)
        self.results_table.setHorizontalHeaderLabels([
            "File", "Status", "PixScale Fit", "Refine", "Resid Med", "Elapsed (s)"
        ])
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.results_table.horizontalHeader().setStretchLastSection(True)
        self.results_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.results_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        results_layout.addWidget(self.results_table)

        layout.addWidget(results_group)

    def setup_astrometrynet_tab(self):
        """Setup Astrometry.net tab UI"""
        layout = QVBoxLayout(self.astrometrynet_tab)

        info_text = (
            "Solve WCS for all frames using local astrometry.net (solve-field). "
            "This optional fallback needs a local solve-field installation and "
            "matching astrometry.net index files; see Help > WCS Solver Installation Help."
        )
        info_style = "QLabel { background-color: #E8F5E9; padding: 10px; border-radius: 5px; }"
        info = QLabel(info_text)
        info.setWordWrap(True)
        info.setStyleSheet(info_style)
        layout.addWidget(info)

        control_layout = QHBoxLayout()

        btn_astnet_params = create_parameter_button("Astrometry.net Parameters")
        btn_astnet_params.clicked.connect(self.open_astrometrynet_parameters_dialog)
        control_layout.addWidget(btn_astnet_params)

        self.run_bar_astnet = RunControlBar(
            "Solve All Frames", "Log",
            run_cb=self.run_astrometrynet_solve,
            stop_cb=self.stop_astrometrynet_solve,
            log_cb=self.show_log_window,
        )
        control_layout.addWidget(self.run_bar_astnet)
        self.btn_solve_astrometrynet = self.run_bar_astnet.btn_run
        self.btn_stop_astrometrynet = self.run_bar_astnet.btn_stop

        layout.addLayout(control_layout)

        progress_layout = QHBoxLayout()
        self.astrometrynet_progress = QProgressBar()
        self.astrometrynet_progress.setMinimum(0)
        self.astrometrynet_progress.setMaximum(100)
        self.astrometrynet_progress.setValue(0)
        progress_layout.addWidget(self.astrometrynet_progress)

        self.astrometrynet_status = QLabel("Ready")
        self.astrometrynet_status.setMinimumWidth(350)
        progress_layout.addWidget(self.astrometrynet_status)
        layout.addLayout(progress_layout)

        results_group = QGroupBox("Astrometry.net Results")
        results_layout = QVBoxLayout(results_group)

        self.astrometrynet_results_table = QTableWidget()
        self.astrometrynet_results_table.setColumnCount(6)
        self.astrometrynet_results_table.setHorizontalHeaderLabels([
            "File", "Status", "PixScale", "Refine", "Resid(\")", "Elapsed (s)"
        ])
        self.astrometrynet_results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.astrometrynet_results_table.horizontalHeader().setStretchLastSection(True)
        self.astrometrynet_results_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.astrometrynet_results_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        results_layout.addWidget(self.astrometrynet_results_table)

        layout.addWidget(results_group)

    def run_astrometrynet_solve(self):
        if not self.file_list:
            QMessageBox.warning(self, "Warning", "No frames found to solve")
            return

        self.log_text.clear()
        file_list = self._qc_filter_and_log()
        if file_list is None:
            return

        pf = self._log_preflight(astnet=True, gaia=True)
        astnet_ok, astnet_detail = pf["astnet"]
        gaia_ok, gaia_detail = pf["gaia"]
        if not astnet_ok:
            self.astrometrynet_status.setText("solve-field not found or not configured")
            QMessageBox.warning(
                self,
                "Astrometry.net Not Found",
                f"{astnet_detail}\n\n"
                "Install local astrometry.net with matching index files, or update Step 5 > "
                "Astrometry.net Parameters. See Help > WCS Solver Installation Help.",
            )
            return

        if not gaia_ok:
            QMessageBox.warning(
                self,
                "Gaia Runtime Not Available",
                f"{gaia_detail}\n\n"
                "Local astrometry.net can still solve WCS, but Gaia attach, "
                "refine, residual medians, and Step 6 Gaia stats will be absent "
                "unless a compatible gaia_fov.ecsv cache already exists.",
            )

        target_coord = None
        ra = getattr(self.params.P, "target_ra_deg", None)
        dec = getattr(self.params.P, "target_dec_deg", None)
        if ra is not None and dec is not None:
            try:
                target_coord = SkyCoord(float(ra) * u.deg, float(dec) * u.deg)
            except Exception:
                target_coord = None

        self.astrometrynet_results_table.setRowCount(0)
        self.results = {}

        # Start solving in thread
        self.astrometrynet_worker = AstrometryNetWorker(
            file_list,
            self.params,
            self.params.P.data_dir,
            self.params.P.result_dir,
            self.params.P.cache_dir,
            self.use_cropped,
            target_coord=target_coord,
        )
        self.astrometrynet_worker.progress.connect(self.on_astrometrynet_progress)
        self.astrometrynet_worker.file_done.connect(self.on_astrometrynet_file_done)
        self.astrometrynet_worker.refine_done.connect(self.on_astrometrynet_refine_done)
        self.astrometrynet_worker.log_message.connect(self.log)
        self.astrometrynet_worker.finished.connect(self.on_astrometrynet_finished)
        self.astrometrynet_worker.error.connect(self.on_astrometrynet_error)
        self.setup_log_window()
        if hasattr(self, "_worker_panel") and self._worker_panel is not None:
            self._worker_panel.clear()
            self.astrometrynet_worker.worker_status.connect(self._worker_panel.update_worker)

        self.run_bar_astnet.set_running(True)
        self.astrometrynet_progress.setValue(0)
        self._astnet_start_time = time.monotonic()
        self.astrometrynet_status.setText(
            progress_status_text(0, len(file_list), self._astnet_start_time, message="Starting local astrometry.net...")
        )
        self.log("=" * 50)
        self.log("Starting local astrometry.net (solve-field) plate solving...")
        self.log(f"Frames: {len(file_list)}")
        self.astrometrynet_worker.start()
        self.show_log_window()

    def stop_astrometrynet_solve(self):
        if self.astrometrynet_worker and self.astrometrynet_worker.isRunning():
            self.run_bar_astnet.set_stopping()
            self.astrometrynet_status.setText("Stopping...")
            self.log("Astrometry.net stop requested...")
            self.astrometrynet_worker.stop()

    def on_astrometrynet_progress(self, current, total, status):
        pct = int(100 * current / max(1, total))
        self.astrometrynet_progress.setValue(pct)
        self.astrometrynet_status.setText(
            progress_status_text(current, total, getattr(self, "_astnet_start_time", None), message=status)
        )

    def on_astrometrynet_file_done(self, filename, result):
        row = self.astrometrynet_results_table.rowCount()
        self.astrometrynet_results_table.insertRow(row)
        self.astrometrynet_results_table.setItem(row, 0, QTableWidgetItem(filename))
        self.astrometrynet_results_table.setItem(row, 1, QTableWidgetItem(result.get("status", "")))
        pixscale = float(result.get("pixscale", 0.0))
        refine = result.get("refine", "-")
        resid_med = result.get("resid_med", np.nan)
        elapsed = float(result.get("elapsed_s", 0.0))
        self.astrometrynet_results_table.setItem(row, 2, QTableWidgetItem(f"{pixscale:.4f}" if np.isfinite(pixscale) and pixscale > 0 else "-"))
        self.astrometrynet_results_table.setItem(row, 3, QTableWidgetItem(str(refine) if refine else "-"))
        self.astrometrynet_results_table.setItem(row, 4, QTableWidgetItem(f"{resid_med:.2f}" if np.isfinite(resid_med) else "-"))
        self.astrometrynet_results_table.setItem(row, 5, QTableWidgetItem(f"{elapsed:.1f}" if np.isfinite(elapsed) and elapsed > 0 else "-"))
        set_table_row_background(
            self.astrometrynet_results_table,
            row,
            self._wcs_result_row_background(result),
        )

        if result.get("ok"):
            self.results[filename] = result
            self.log(f"Astrometry.net solved: {filename} (RA={result.get('ra', 0):.4f}, Dec={result.get('dec', 0):.4f})")
        try:
            self._refresh_solver_breakdown_label()
        except Exception:
            pass

    def on_astrometrynet_error(self, filename, error):
        self.log(f"Astrometry.net ERROR {filename}: {error}")

    def on_astrometrynet_refine_done(self, filename, result):
        """Refine 결과로 테이블 업데이트"""
        # 파일명으로 해당 행 찾기
        for row in range(self.astrometrynet_results_table.rowCount()):
            item = self.astrometrynet_results_table.item(row, 0)
            if item and item.text() == filename:
                refine = result.get("refine", "-")
                resid_med = result.get("resid_med", np.nan)
                self.astrometrynet_results_table.setItem(row, 3, QTableWidgetItem(str(refine) if refine else "-"))
                self.astrometrynet_results_table.setItem(row, 4, QTableWidgetItem(f"{resid_med:.2f}" if np.isfinite(resid_med) else "-"))
                # results에도 업데이트
                if filename in self.results:
                    self.results[filename].update(result)
                    row_result = self.results[filename]
                else:
                    row_result = result
                set_table_row_background(
                    self.astrometrynet_results_table,
                    row,
                    self._wcs_result_row_background(row_result),
                )
                break

    def on_astrometrynet_finished(self, summary):
        self.run_bar_astnet.set_running(False)
        stopped = bool(summary.get("stopped")) if isinstance(summary, dict) else False
        n_ok = summary.get("ok", 0)
        n_qc = summary.get("wcs_qc_pass", 0)
        n_qc_not_eval = int(summary.get("wcs_qc_not_evaluated", 0) or 0)
        if not stopped:
            self.astrometrynet_progress.setValue(100)
            self.astrometrynet_status.setText(f"Done: {n_ok}/{summary.get('total', 0)} solved")
        else:
            self.astrometrynet_status.setText(f"Stopped: {n_ok}/{summary.get('total', 0)} solved")
        if n_ok > 0:
            qc_text = (
                f"WCS-QC not evaluated: {n_qc_not_eval} (Gaia unavailable)"
                if n_qc_not_eval
                else f"WCS-QC pass: {n_qc}"
            )
            self.log(f"Astrometry.net: {n_ok} frames solved successfully | {qc_text}")
        if stopped:
            self.log("Astrometry.net solve stopped by user")
        self.save_state()
        self.update_navigation_buttons()

    def setup_log_window(self):
        if self.log_window is not None:
            return

        worker_group = QGroupBox("Workers")
        worker_group.setMinimumWidth(430)
        wg_layout = QVBoxLayout(worker_group)
        wg_layout.setContentsMargins(5, 5, 5, 5)
        self._worker_panel = WorkerStatusPanel(worker_group)
        wg_layout.addWidget(self._worker_panel)

        self.log_window = WorkflowLogWindow(
            self, "WCS Log & Workers", width=900, height=500,
            side_widget=worker_group,
        )
        self.log_text = self.log_window.log_text

    def show_log_window(self):
        if self.log_window is None:
            self.setup_log_window()
        show_raised(self.log_window)

    def log(self, message: str):
        append_timestamped_log(self.log_text, message)

    def populate_file_list(self):
        crop_active = crop_is_active(self.params.P.result_dir)
        cropped_dir = step2_cropped_dir(self.params.P.result_dir)
        excluded = getattr(self.file_manager, "excluded_files", set()) if self.file_manager else set()

        if crop_active and cropped_dir.exists() and list(cropped_dir.glob("*.fit*")):
            files = sorted([f.name for f in cropped_dir.glob("*.fit*")
                            if f.name not in excluded])
            self.use_cropped = True
        else:
            files = self.file_manager.get_file_list() if self.file_manager else []
            self.use_cropped = False

        self.file_list = list(files)

    def open_parameters_dialog(self):
        dialog = QDialog(self)
        configure_parameter_dialog(dialog, "WCS Parameters", 560, 720)

        layout = QVBoxLayout(dialog)
        scroll = QScrollArea(dialog)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        content = QWidget()
        body = QVBoxLayout(content)
        body.setContentsMargins(4, 4, 4, 4)
        body.setSpacing(8)
        scroll.setWidget(content)
        layout.addWidget(scroll, 1)

        _info = QLabel(
            "Adjust WCS plate-solving parameters. ASTAP and the selected D80/D50 "
            "star database must be installed outside APEX. See Help > WCS Solver "
            "Installation Help."
        )
        _info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; margin-bottom: 10px; }")
        _info.setWordWrap(True)
        body.addWidget(_info)

        def _add_group(title: str, *, expanded: bool = False) -> QFormLayout:
            group, container = create_collapsible_section(title, initial_expanded=expanded)
            form = QFormLayout(container)
            form.setLabelAlignment(Qt.AlignRight)
            form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
            body.addWidget(group)
            return form

        astap_form = _add_group("ASTAP Solver", expanded=True)
        refine_form = _add_group("WCS Refinement", expanded=True)
        gaia_form = _add_group("Gaia / Hybrid ID")

        self.param_astap_exe = QLineEdit(str(getattr(self.params.P, "astap_exe", "astap_cli.exe")))
        astap_exe_row = QWidget()
        astap_exe_layout = QHBoxLayout(astap_exe_row)
        astap_exe_layout.setContentsMargins(0, 0, 0, 0)
        astap_exe_layout.setSpacing(6)
        astap_exe_layout.addWidget(self.param_astap_exe, 1)
        btn_browse_astap = QPushButton("Browse...")
        btn_browse_astap.clicked.connect(self.browse_astap_cli)
        astap_exe_layout.addWidget(btn_browse_astap)
        astap_form.addRow("ASTAP CLI Path:", astap_exe_row)

        self.param_timeout = QDoubleSpinBox()
        self.param_timeout.setRange(10, 1000)
        self.param_timeout.setValue(float(getattr(self.params.P, "astap_timeout_s", 120.0)))
        astap_form.addRow("Timeout (s):", self.param_timeout)

        self.param_radius = QDoubleSpinBox()
        self.param_radius.setRange(0.5, 30.0)
        self.param_radius.setValue(float(getattr(self.params.P, "astap_search_radius_deg", 8.0)))
        astap_form.addRow("Search Radius (deg):", self.param_radius)

        self.param_astap_db = QComboBox()
        self.param_astap_db.addItems(["D80", "D50"])
        current_db = str(getattr(self.params.P, "astap_database", "D80"))
        idx = self.param_astap_db.findText(current_db)
        if idx >= 0:
            self.param_astap_db.setCurrentIndex(idx)
        astap_form.addRow("ASTAP Star DB:", self.param_astap_db)

        self.param_annotate_variables = QCheckBox("Enable")
        self.param_annotate_variables.setChecked(bool(getattr(self.params.P, "astap_annotate_variables", False)))
        self.param_annotate_variables.setToolTip("ASTAP 변광성 데이터베이스로 변광성 주석 표시 (별도 설치 필요)")
        astap_form.addRow("Annotate Variable Stars:", self.param_annotate_variables)

        self.param_fov_fudge = QDoubleSpinBox()
        self.param_fov_fudge.setRange(0.5, 2.0)
        self.param_fov_fudge.setSingleStep(0.05)
        self.param_fov_fudge.setValue(float(getattr(self.params.P, "astap_fov_fudge", 1.0)))
        astap_form.addRow("FOV Fudge:", self.param_fov_fudge)

        self.param_downsample = QSpinBox()
        self.param_downsample.setRange(1, 8)
        self.param_downsample.setValue(int(getattr(self.params.P, "astap_downsample_z", 2)))
        astap_form.addRow("Downsample Z:", self.param_downsample)

        self.param_max_stars = QSpinBox()
        self.param_max_stars.setRange(50, 5000)
        self.param_max_stars.setValue(int(getattr(self.params.P, "astap_max_stars_s", 500)))
        astap_form.addRow("Max Stars (S):", self.param_max_stars)

        self.param_max_workers = QSpinBox()
        self.param_max_workers.setRange(1, 16)
        self.param_max_workers.setValue(int(getattr(self.params.P, "wcs_max_workers", 1)))
        astap_form.addRow("Max Workers:", self.param_max_workers)

        self.param_require_qc = QCheckBox("Enable")
        self.param_require_qc.setChecked(bool(getattr(self.params.P, "wcs_require_qc_pass", True)))
        refine_form.addRow("QC Pass Only:", self.param_require_qc)

        self.param_refine_enable = QCheckBox("Enable")
        self.param_refine_enable.setChecked(bool(getattr(self.params.P, "wcs_refine_enable", True)))
        refine_form.addRow("Refine CRPIX:", self.param_refine_enable)

        self.param_refine_max_match = QSpinBox()
        self.param_refine_max_match.setRange(50, 5000)
        self.param_refine_max_match.setValue(int(getattr(self.params.P, "wcs_refine_max_match", 600)))
        refine_form.addRow("Refine Max Match:", self.param_refine_max_match)

        self.param_refine_match_r = QDoubleSpinBox()
        self.param_refine_match_r.setRange(0.5, 5.0)
        self.param_refine_match_r.setSingleStep(0.1)
        self.param_refine_match_r.setValue(float(getattr(self.params.P, "wcs_refine_match_r_fwhm", 1.6)))
        refine_form.addRow("Refine Match R (×FWHM):", self.param_refine_match_r)

        self.param_refine_min_match = QSpinBox()
        self.param_refine_min_match.setRange(5, 500)
        self.param_refine_min_match.setValue(int(getattr(self.params.P, "wcs_refine_min_match", 50)))
        refine_form.addRow("Refine Min Match:", self.param_refine_min_match)

        self.param_gaia_fudge = QDoubleSpinBox()
        self.param_gaia_fudge.setRange(0.5, 3.0)
        self.param_gaia_fudge.setSingleStep(0.05)
        self.param_gaia_fudge.setValue(float(getattr(self.params.P, "gaia_radius_fudge", 1.35)))
        gaia_form.addRow("Gaia Radius Fudge:", self.param_gaia_fudge)

        self.param_gaia_mag_max = QDoubleSpinBox()
        self.param_gaia_mag_max.setRange(10.0, 25.0)
        self.param_gaia_mag_max.setSingleStep(0.5)
        self.param_gaia_mag_max.setValue(float(getattr(self.params.P, "gaia_mag_max", 18.0)))
        gaia_form.addRow("Gaia Mag Max:", self.param_gaia_mag_max)

        self.param_gaia_wcs_mag_max = QDoubleSpinBox()
        self.param_gaia_wcs_mag_max.setRange(10.0, 21.0)
        self.param_gaia_wcs_mag_max.setSingleStep(0.5)
        self.param_gaia_wcs_mag_max.setValue(float(getattr(self.params.P, "gaia_wcs_mag_max", 18.0)))
        self.param_gaia_wcs_mag_max.setToolTip(
            "Server-side Gaia G magnitude cap used by Step 5 WCS refinement/QC. "
            "Keep this near 18-20 to avoid TAP timeouts."
        )
        gaia_form.addRow("Gaia WCS Mag Max:", self.param_gaia_wcs_mag_max)

        self.param_ref_gaia_match_tol = QDoubleSpinBox()
        self.param_ref_gaia_match_tol.setRange(0.1, 30.0)
        self.param_ref_gaia_match_tol.setDecimals(2)
        self.param_ref_gaia_match_tol.setSingleStep(0.1)
        self.param_ref_gaia_match_tol.setValue(float(getattr(self.params.P, "ref_wcs_match_radius_arcsec", 2.0)))
        gaia_form.addRow("Gaia Match Tol (Ref, arcsec):", self.param_ref_gaia_match_tol)

        self.param_gaia_g_limit = QDoubleSpinBox()
        self.param_gaia_g_limit.setRange(10.0, 25.0)
        self.param_gaia_g_limit.setDecimals(2)
        self.param_gaia_g_limit.setSingleStep(0.5)
        self.param_gaia_g_limit.setValue(
            float(
                getattr(
                    self.params.P,
                    "idmatch_gaia_g_limit",
                    getattr(self.params.P, "gaia_mag_max", 18.0),
                )
            )
        )
        gaia_form.addRow("Gaia G limit (Hybrid ID):", self.param_gaia_g_limit)

        self.param_gaia_retry = QSpinBox()
        self.param_gaia_retry.setRange(0, 10)
        self.param_gaia_retry.setValue(int(getattr(self.params.P, "gaia_retry", 2)))
        gaia_form.addRow("Gaia Retry:", self.param_gaia_retry)

        self.param_gaia_timeout = QDoubleSpinBox()
        self.param_gaia_timeout.setRange(5.0, 300.0)
        self.param_gaia_timeout.setSingleStep(5.0)
        self.param_gaia_timeout.setValue(float(getattr(self.params.P, "gaia_timeout_s", 30.0)))
        self.param_gaia_timeout.setSuffix(" s")
        gaia_form.addRow("Gaia Timeout:", self.param_gaia_timeout)

        self.param_gaia_backoff = QDoubleSpinBox()
        self.param_gaia_backoff.setRange(0.0, 30.0)
        self.param_gaia_backoff.setSingleStep(1.0)
        self.param_gaia_backoff.setValue(float(getattr(self.params.P, "gaia_backoff_s", 6.0)))
        gaia_form.addRow("Gaia Backoff (s):", self.param_gaia_backoff)

        self.param_gaia_allow_no_cache = QCheckBox("Allow query when cache missing")
        self.param_gaia_allow_no_cache.setChecked(bool(getattr(self.params.P, "gaia_allow_no_cache", True)))
        self.param_gaia_allow_no_cache.setToolTip(
            "This is not a reuse-cache toggle. It controls whether Step 5 may query Gaia online "
            "when no local Gaia cache exists."
        )
        gaia_form.addRow("Gaia cache miss:", self.param_gaia_allow_no_cache)

        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        add_parameter_reset_button(
            buttons,
            [
                (self.param_astap_exe, "C:/Program Files/astap/astap_cli.exe"),
                (self.param_timeout, 120.0),
                (self.param_radius, 8.0),
                (self.param_astap_db, "D50"),
                (self.param_annotate_variables, False),
                (self.param_fov_fudge, 1.0),
                (self.param_downsample, 2),
                (self.param_max_stars, 500),
                (self.param_max_workers, 4),
                (self.param_require_qc, True),
                (self.param_refine_enable, True),
                (self.param_refine_max_match, 600),
                (self.param_refine_match_r, 2.0),
                (self.param_refine_min_match, 25),
                (self.param_gaia_fudge, 1.35),
                (self.param_gaia_mag_max, 25.0),
                (self.param_gaia_wcs_mag_max, 18.0),
                (self.param_ref_gaia_match_tol, 2.0),
                (self.param_gaia_g_limit, 25.0),
                (self.param_gaia_retry, 2),
                (self.param_gaia_timeout, 30.0),
                (self.param_gaia_backoff, 6.0),
                (self.param_gaia_allow_no_cache, True),
            ],
        )
        buttons.accepted.connect(lambda: self.save_parameters(dialog))
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        dialog.exec_()

    def open_astrometrynet_parameters_dialog(self):
        dialog = QDialog(self)
        configure_parameter_dialog(dialog, "Astrometry.net Parameters", 540, 620)

        layout = QVBoxLayout(dialog)
        scroll = QScrollArea(dialog)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        content = QWidget()
        body = QVBoxLayout(content)
        body.setContentsMargins(4, 4, 4, 4)
        body.setSpacing(8)
        scroll.setWidget(content)
        layout.addWidget(scroll, 1)

        _info = QLabel(
            "Configure Astrometry.net local solver as a fallback for ASTAP failures. "
            "Local solve-field and matching index files must be installed outside APEX. "
            "See Help > WCS Solver Installation Help."
        )
        _info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; margin-bottom: 10px; }")
        _info.setWordWrap(True)
        body.addWidget(_info)

        def _add_group(title: str, *, expanded: bool = False) -> QFormLayout:
            group, container = create_collapsible_section(title, initial_expanded=expanded)
            form = QFormLayout(container)
            form.setLabelAlignment(Qt.AlignRight)
            form.setFieldGrowthPolicy(QFormLayout.AllNonFixedFieldsGrow)
            body.addWidget(group)
            return form

        form = _add_group("Local Solve", expanded=True)
        retry_form = _add_group("Blind Retry")

        self.param_astnet_enable = QCheckBox("Enable")
        self.param_astnet_enable.setChecked(bool(getattr(self.params.P, "astnet_local_enable", False)))
        form.addRow("Enable Local Solve:", self.param_astnet_enable)

        self.param_astnet_use_wsl = QCheckBox("Use WSL")
        self.param_astnet_use_wsl.setChecked(bool(getattr(self.params.P, "astnet_local_use_wsl", True)))
        form.addRow("Use WSL:", self.param_astnet_use_wsl)

        self.param_astnet_command = QLineEdit(str(getattr(self.params.P, "astnet_local_command", "solve-field")))
        form.addRow("solve-field Command:", self.param_astnet_command)

        self.param_astnet_timeout = QDoubleSpinBox()
        self.param_astnet_timeout.setRange(30, 3600)
        self.param_astnet_timeout.setValue(float(getattr(self.params.P, "astnet_local_timeout_s", 300.0)))
        form.addRow("Timeout (s):", self.param_astnet_timeout)

        self.param_astnet_downsample = QSpinBox()
        self.param_astnet_downsample.setRange(1, 8)
        self.param_astnet_downsample.setValue(int(getattr(self.params.P, "astnet_local_downsample", 2)))
        form.addRow("Downsample:", self.param_astnet_downsample)

        self.param_astnet_scale_low = QDoubleSpinBox()
        self.param_astnet_scale_low.setRange(0.0, 10.0)
        self.param_astnet_scale_low.setDecimals(5)
        self.param_astnet_scale_low.setValue(float(getattr(self.params.P, "astnet_local_scale_low", 0.0)))
        form.addRow("Scale Low (arcsec/pix):", self.param_astnet_scale_low)

        self.param_astnet_scale_high = QDoubleSpinBox()
        self.param_astnet_scale_high.setRange(0.0, 10.0)
        self.param_astnet_scale_high.setDecimals(5)
        self.param_astnet_scale_high.setValue(float(getattr(self.params.P, "astnet_local_scale_high", 0.0)))
        form.addRow("Scale High (arcsec/pix):", self.param_astnet_scale_high)

        self.param_astnet_radius = QDoubleSpinBox()
        self.param_astnet_radius.setRange(0.1, 30.0)
        self.param_astnet_radius.setValue(float(getattr(self.params.P, "astnet_local_radius_deg", 8.0)))
        form.addRow("Radius (deg):", self.param_astnet_radius)

        self.param_astnet_keep_outputs = QCheckBox("Keep .wcs/debug files")
        self.param_astnet_keep_outputs.setChecked(bool(getattr(self.params.P, "astnet_local_keep_outputs", True)))
        self.param_astnet_keep_outputs.setToolTip(
            "Keep .wcs and debug sidecars; redundant .new FITS copies are removed automatically."
        )
        form.addRow("Keep Debug Outputs:", self.param_astnet_keep_outputs)

        self.param_astnet_use_cache = create_cache_checkbox(
            "Use Cache",
            bool(getattr(self.params.P, "astnet_local_use_cache", True)),
            "Reuse compatible local astrometry.net sidecar outputs instead of running solve-field again.",
        )
        form.addRow("Use Cached Outputs:", self.param_astnet_use_cache)

        self.param_astnet_max_objs = QSpinBox()
        self.param_astnet_max_objs.setRange(100, 20000)
        self.param_astnet_max_objs.setValue(int(getattr(self.params.P, "astnet_local_max_objs", 2000)))
        form.addRow("Max Objects:", self.param_astnet_max_objs)

        self.param_astnet_cpulimit = QDoubleSpinBox()
        self.param_astnet_cpulimit.setRange(1, 300)
        self.param_astnet_cpulimit.setValue(float(getattr(self.params.P, "astnet_local_cpulimit_s", 30.0)))
        form.addRow("CPU Limit (s):", self.param_astnet_cpulimit)

        self.param_astnet_blind_retry = QCheckBox("Retry blind when hint-based solve fails")
        self.param_astnet_blind_retry.setChecked(bool(getattr(self.params.P, "astnet_blind_retry_on_fail", True)))
        retry_form.addRow("Blind Retry:", self.param_astnet_blind_retry)

        self.param_astnet_blind_cpulimit = QDoubleSpinBox()
        self.param_astnet_blind_cpulimit.setRange(10, 600)
        self.param_astnet_blind_cpulimit.setValue(float(getattr(self.params.P, "astnet_blind_cpulimit_s", 120.0)))
        self.param_astnet_blind_cpulimit.setSuffix(" s")
        retry_form.addRow("Blind CPU Limit:", self.param_astnet_blind_cpulimit)

        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        add_parameter_reset_button(
            buttons,
            [
                (self.param_astnet_enable, True),
                (self.param_astnet_use_wsl, True),
                (self.param_astnet_command, "solve-field"),
                (self.param_astnet_timeout, 300.0),
                (self.param_astnet_downsample, 2),
                (self.param_astnet_scale_low, 0.3),
                (self.param_astnet_scale_high, 0.5),
                (self.param_astnet_radius, 30.0),
                (self.param_astnet_keep_outputs, True),
                (self.param_astnet_use_cache, True),
                (self.param_astnet_max_objs, 2000),
                (self.param_astnet_cpulimit, 30.0),
                (self.param_astnet_blind_retry, True),
                (self.param_astnet_blind_cpulimit, 120.0),
            ],
        )
        buttons.accepted.connect(lambda: self.save_astrometrynet_parameters(dialog))
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        dialog.exec_()

    def browse_astap_cli(self):
        current = _strip_outer_quotes(self.param_astap_exe.text())
        start_dir = Path(current).parent if current and Path(current).parent.exists() else Path.home()
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select ASTAP CLI executable",
            str(start_dir),
            "ASTAP CLI (astap_cli.exe);;Executables (*.exe);;All Files (*)",
        )
        if path:
            self.param_astap_exe.setText(path)

    def save_parameters(self, dialog):
        self.params.P.astap_exe = self.param_astap_exe.text().strip()
        self.params.P.astap_timeout_s = self.param_timeout.value()
        self.params.P.astap_search_radius_deg = self.param_radius.value()
        self.params.P.astap_database = self.param_astap_db.currentText()
        self.params.P.astap_annotate_variables = self.param_annotate_variables.isChecked()
        self.params.P.astap_fov_fudge = self.param_fov_fudge.value()
        self.params.P.astap_downsample_z = self.param_downsample.value()
        self.params.P.astap_max_stars_s = self.param_max_stars.value()
        self.params.P.wcs_max_workers = self.param_max_workers.value()
        self.params.P.wcs_require_qc_pass = self.param_require_qc.isChecked()
        self.params.P.wcs_refine_enable = self.param_refine_enable.isChecked()
        self.params.P.wcs_refine_max_match = self.param_refine_max_match.value()
        self.params.P.wcs_refine_match_r_fwhm = self.param_refine_match_r.value()
        self.params.P.wcs_refine_min_match = self.param_refine_min_match.value()
        self.params.P.gaia_radius_fudge = self.param_gaia_fudge.value()
        self.params.P.gaia_mag_max = self.param_gaia_mag_max.value()
        self.params.P.gaia_wcs_mag_max = self.param_gaia_wcs_mag_max.value()
        self.params.P.ref_wcs_match_radius_arcsec = self.param_ref_gaia_match_tol.value()
        self.params.P.idmatch_gaia_g_limit = self.param_gaia_g_limit.value()
        self.params.P.gaia_retry = self.param_gaia_retry.value()
        self.params.P.gaia_timeout_s = self.param_gaia_timeout.value()
        self.params.P.gaia_backoff_s = self.param_gaia_backoff.value()
        self.params.P.gaia_allow_no_cache = self.param_gaia_allow_no_cache.isChecked()
        self.persist_params()
        self.save_state()
        QMessageBox.information(dialog, "Success", "Parameters saved!")
        dialog.accept()

    def save_astrometrynet_parameters(self, dialog):
        self.params.P.astnet_local_enable = self.param_astnet_enable.isChecked()
        self.params.P.astnet_local_use_wsl = self.param_astnet_use_wsl.isChecked()
        self.params.P.astnet_local_command = self.param_astnet_command.text().strip()
        self.params.P.astnet_local_timeout_s = self.param_astnet_timeout.value()
        self.params.P.astnet_local_downsample = self.param_astnet_downsample.value()
        self.params.P.astnet_local_scale_low = self.param_astnet_scale_low.value()
        self.params.P.astnet_local_scale_high = self.param_astnet_scale_high.value()
        self.params.P.astnet_local_radius_deg = self.param_astnet_radius.value()
        self.params.P.astnet_local_keep_outputs = self.param_astnet_keep_outputs.isChecked()
        self.params.P.astnet_local_use_cache = self.param_astnet_use_cache.isChecked()
        self.params.P.astnet_local_max_objs = self.param_astnet_max_objs.value()
        self.params.P.astnet_local_cpulimit_s = self.param_astnet_cpulimit.value()
        self.params.P.astnet_blind_retry_on_fail = self.param_astnet_blind_retry.isChecked()
        self.params.P.astnet_blind_cpulimit_s = self.param_astnet_blind_cpulimit.value()
        self.persist_params()
        self.save_state()
        QMessageBox.information(dialog, "Success", "Astrometry.net parameters saved!")
        dialog.accept()

    def run_wcs(self):
        if not self.file_list:
            QMessageBox.warning(self, "Warning", "No files to process")
            return
        if self.worker and self.worker.isRunning():
            return

        self.results = {}
        self.results_table.setRowCount(0)
        self.log_text.clear()
        file_list = self._qc_filter_and_log()
        if file_list is None:
            return
        self.stop_requested = False

        astap_ok, astap_detail = self._log_preflight(astap=True)["astap"]
        astnet_enabled = bool(getattr(self.params.P, "astnet_local_enable", False))
        # astrometry.net is a fallback for ASTAP: only probe it (slow WSL spin-up)
        # when ASTAP is unavailable, otherwise log a deferred note instead.
        if astnet_enabled and not astap_ok:
            astnet_ok, astnet_detail = self._log_preflight(astnet=True)["astnet"]
        else:
            astnet_ok = False
            astnet_detail = (
                "local astrometry.net fallback deferred until ASTAP fails"
                if astnet_enabled
                else "local astrometry.net fallback disabled"
            )
            if not astap_ok:
                self.log(f"[WCS][Preflight] astrometry.net: {astnet_detail}")

        gaia_ok, gaia_detail = self._log_preflight(gaia=True)["gaia"]
        if not gaia_ok:
            QMessageBox.warning(
                self,
                "Gaia Runtime Not Available",
                f"{gaia_detail}\n\n"
                "ASTAP can still solve WCS, but Gaia attach, refine, residual "
                "medians, and Step 6 Gaia stats will be absent unless a compatible "
                "gaia_fov.ecsv cache already exists.",
            )

        if not astap_ok and not astnet_ok:
            self.progress_label.setText("WCS solver not found or not configured")
            QMessageBox.warning(
                self,
                "WCS Solver Not Found",
                f"{astap_detail}\n\n{astnet_detail}\n\n"
                "Install ASTAP with a D80/D50 star database, or enable/install local "
                "astrometry.net with matching index files. See Help > WCS Solver "
                "Installation Help.",
            )
            return

        if not astap_ok and astnet_ok:
            self.progress_label.setText("ASTAP not found; using local astrometry.net fallback")
            QMessageBox.information(
                self,
                "ASTAP Not Found",
                f"{astap_detail}\n\n"
                "Local astrometry.net is enabled and available, so APEX will continue "
                "with solve-field fallback. See Help > WCS Solver Installation Help "
                "for ASTAP/D80/D50 setup.",
            )

        # Get target from params.P (loaded from TOML)
        target_coord = None
        ra = getattr(self.params.P, "target_ra_deg", None)
        dec = getattr(self.params.P, "target_dec_deg", None)
        if ra is not None and dec is not None:
            try:
                target_coord = SkyCoord(float(ra) * u.deg, float(dec) * u.deg)
            except Exception:
                target_coord = None

        self.worker = WcsWorker(
            file_list,
            self.params,
            self.params.P.data_dir,
            self.params.P.result_dir,
            self.params.P.cache_dir,
            self.use_cropped,
            target_coord=target_coord
        )
        self.worker.progress.connect(self.on_progress)
        self.worker.file_done.connect(self.on_file_done)
        self.worker.log_message.connect(self.log)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)
        self.setup_log_window()
        if hasattr(self, "_worker_panel") and self._worker_panel is not None:
            self._worker_panel.clear()
            self.worker.worker_status.connect(self._worker_panel.update_worker)

        self.run_bar_astap.set_running(True)
        self.progress_bar.setValue(0)
        self.progress_bar.setMaximum(len(file_list))
        self._wcs_start_time = time.monotonic()
        self.progress_label.setText(
            progress_status_text(0, len(file_list), self._wcs_start_time, message="Starting...")
        )
        self.worker.start()
        self.show_log_window()

    def stop_wcs(self):
        if self.worker and self.worker.isRunning():
            self.stop_requested = True
            self.run_bar_astap.set_stopping()
            self.progress_label.setText("Stopping...")
            self.log("Stop requested...")
            self.worker.stop()

    def on_progress(self, current, total, filename):
        self.progress_bar.setValue(current)
        self.progress_label.setText(
            progress_status_text(current, total, getattr(self, "_wcs_start_time", None), message=filename)
        )

    @staticmethod
    def _boolish(value) -> bool:
        if isinstance(value, (bool, np.bool_)):
            return bool(value)
        if value is None:
            return False
        if isinstance(value, (int, float, np.integer, np.floating)):
            return bool(value) and np.isfinite(value)
        return str(value).strip().lower() in {"1", "true", "t", "yes", "y", "ok", "solved"}

    @classmethod
    def _is_successful_wcs_result(cls, result: dict) -> bool:
        if not isinstance(result, dict):
            return False
        if cls._boolish(result.get("ok")):
            return True
        status = str(result.get("status", "")).strip().lower()
        return status in {"ok", "ok_astnet_wsl", "solved"}

    @classmethod
    def _wcs_result_row_background(cls, result: dict) -> str:
        ok = cls._is_successful_wcs_result(result)
        warning = False
        if ok and isinstance(result, dict):
            if "wcs_qc_pass" in result and not cls._boolish(result.get("wcs_qc_pass")):
                warning = True
            refine = str(result.get("refine", "") or "").strip().lower()
            if refine.startswith("match_too_small") or refine.startswith("error"):
                warning = True
        return status_row_background(ok, warning=warning)

    def _restore_success_results_from_summary(self):
        summary_path = step5_wcs_dir(self.params.P.result_dir) / "wcs_solve_summary.csv"
        if not summary_path.exists():
            return
        try:
            df = pd.read_csv(summary_path)
        except Exception:
            return
        if df.empty:
            return

        file_col = "file" if "file" in df.columns else "fname" if "fname" in df.columns else None
        if not file_col:
            return

        current_files = {Path(str(f)).name for f in getattr(self, "file_list", []) or []}
        restored = {}
        skipped_foreign = 0
        for row in df.to_dict("records"):
            if not self._is_successful_wcs_result(row):
                continue
            filename = row.get(file_col)
            if filename is None or pd.isna(filename):
                continue
            fname = Path(str(filename)).name
            if current_files and fname not in current_files:
                skipped_foreign += 1
                continue
            row["file"] = fname
            row.setdefault("fname", fname)
            restored[fname] = row

        if restored:
            self.results = restored
            counts = self._populate_restored_wcs_tables(restored)
            total = len(restored)
            parts = []
            if counts.get("internal"):
                parts.append(f"Internal {counts['internal']}")
                if hasattr(self, "internal_status"):
                    self.internal_status.setText(f"Loaded previous results: {counts['internal']} solved")
            if counts.get("astap"):
                parts.append(f"ASTAP {counts['astap']}")
            if counts.get("astnet"):
                parts.append(f"Astrometry.net {counts['astnet']}")
                if hasattr(self, "astrometrynet_status"):
                    self.astrometrynet_status.setText(f"Loaded previous results: {counts['astnet']} solved")
            if hasattr(self, "progress_label"):
                detail = ", ".join(parts) if parts else f"{total} solved"
                self.progress_label.setText(f"Loaded previous WCS summary: {detail}")
            try:
                self.log(f"[WCS][CACHE] Loaded previous WCS summary: {total} solved frame(s).")
            except Exception:
                pass
        elif skipped_foreign:
            try:
                self.log(f"[WCS][CACHE] Ignored {skipped_foreign} WCS row(s) not in the current file list.")
            except Exception:
                pass
        try:
            self._refresh_solver_breakdown_label()
        except Exception:
            pass

    def _get_wcs_cache_mgr(self) -> StepCacheManager:
        if self._wcs_cache_mgr is None:
            self._wcs_cache_mgr = StepCacheManager(
                self.params.P.cache_dir, "wcs_plate_solving", cache_schema_version=1
            )
        return self._wcs_cache_mgr

    def _write_wcs_manifest(self, filename: str, result: dict) -> None:
        try:
            from apex.utils.step_paths import step5_wcs_dir
            fits_path = (
                step2_cropped_dir(self.params.P.result_dir) / filename
                if self.use_cropped
                else Path(self.params.P.data_dir) / filename
            )
            wcs_out = step5_wcs_dir(self.params.P.result_dir)
            stem = Path(filename).stem
            mgr = self._get_wcs_cache_mgr()
            manifest = mgr.build_manifest(
                input_paths=[fits_path] if fits_path.exists() else [],
                payload_paths={"wcs_summary": wcs_out / "wcs_solve_summary.csv"},
                extra={"status": result.get("status", ""), "filename": filename},
            )
            mgr.write_manifest(filename, manifest)
        except Exception:
            pass

    @staticmethod
    def _wcs_failure_detail(result: dict, *, limit: int = 260) -> str:
        reason = str(result.get("fail_reason", "") or "").strip()
        for key in (
            "astap_stderr",
            "astap_stdout",
            "astnet_wsl_stderr",
            "astnet_wsl_stdout",
        ):
            tail = _tail_text(result.get(key), limit=limit, max_lines=2)
            if tail:
                return f"{reason} | {tail}" if reason else tail
        return reason

    @staticmethod
    def _wcs_failure_tooltip(result: dict) -> str:
        parts: list[str] = []
        reason = str(result.get("fail_reason", "") or "").strip()
        if reason:
            parts.append(f"reason: {reason}")
        for label, key in (
            ("ASTAP stderr", "astap_stderr"),
            ("ASTAP stdout", "astap_stdout"),
            ("astrometry.net stderr", "astnet_wsl_stderr"),
            ("astrometry.net stdout", "astnet_wsl_stdout"),
        ):
            tail = _tail_text(result.get(key), limit=900, max_lines=5)
            if tail:
                parts.append(f"{label}: {tail}")
        return "\n".join(parts)

    def on_file_done(self, filename, result):
        if self._is_successful_wcs_result(result):
            self.results[filename] = result
        else:
            self.results.pop(filename, None)
        try:
            self._refresh_solver_breakdown_label()
        except Exception:
            pass
        row = self.results_table.rowCount()
        self.results_table.insertRow(row)
        self.results_table.setItem(row, 0, QTableWidgetItem(filename))
        status_item = QTableWidgetItem(str(result.get("status", "")))
        if not self._is_successful_wcs_result(result):
            tooltip = self._wcs_failure_tooltip(result)
            if tooltip:
                status_item.setToolTip(tooltip)
        self.results_table.setItem(row, 1, status_item)
        pix_fit = result.get("pix_fit")
        pix_str = f"{pix_fit:.4f}" if isinstance(pix_fit, float) and np.isfinite(pix_fit) else "-"
        self.results_table.setItem(row, 2, QTableWidgetItem(pix_str))
        refine = result.get("refine", "")
        self.results_table.setItem(row, 3, QTableWidgetItem(str(refine)))
        resid_med = result.get("resid_med")
        if isinstance(resid_med, float) and np.isfinite(resid_med):
            resid_str = f"{resid_med:.3f}\""
        else:
            resid_px = result.get("resid_med_px")
            resid_str = f"{resid_px:.3f}px" if isinstance(resid_px, float) and np.isfinite(resid_px) else "-"
        self.results_table.setItem(row, 4, QTableWidgetItem(resid_str))
        elapsed = result.get("elapsed", 0.0)
        self.results_table.setItem(row, 5, QTableWidgetItem(f"{elapsed:.1f}"))
        set_table_row_background(self.results_table, row, self._wcs_result_row_background(result))
        self.results_table.scrollToBottom()
        pix_fit_log = result.get("pix_fit")
        pix_log = f"{pix_fit_log:.4f}" if isinstance(pix_fit_log, float) and np.isfinite(pix_fit_log) else "-"
        resid_log = result.get("resid_med")
        if isinstance(resid_log, float) and np.isfinite(resid_log):
            resid_str = f"{resid_log:.3f}\""
        else:
            resid_px = result.get("resid_med_px")
            resid_str = f"{resid_px:.3f}px" if isinstance(resid_px, float) and np.isfinite(resid_px) else "-"
        if self._is_successful_wcs_result(result):
            self.log(f"{filename}: {result.get('status', '')} pix={pix_log} refine={refine or '-'} resid_med={resid_str}")
            self._write_wcs_manifest(filename, result)
        else:
            detail = self._wcs_failure_detail(result)
            detail_txt = f" reason={detail}" if detail else ""
            self.log(
                f"{filename}: {result.get('status', '')} pix={pix_log} "
                f"refine={refine or '-'} resid_med={resid_str}{detail_txt}"
            )

    def on_error(self, filename, error):
        self.log(f"ERROR {filename}: {error}")

    def on_finished(self, summary):
        self.run_bar_astap.set_running(False)
        self.stop_requested = False
        stopped = bool(summary.get("stopped")) if isinstance(summary, dict) else False
        elapsed_txt = ""
        if hasattr(self, "_wcs_start_time"):
            elapsed_txt = f" | elapsed {format_duration(time.monotonic() - self._wcs_start_time)}"
        self.progress_label.setText(("Stopped" if stopped else "Done") + elapsed_txt)
        if summary:
            qc_not_eval = int(summary.get("wcs_qc_not_evaluated", 0) or 0)
            qc_text = (
                f"WCS-QC not evaluated: {qc_not_eval} (Gaia unavailable)"
                if qc_not_eval
                else f"WCS-QC pass: {summary.get('wcs_qc_pass', 0)}"
            )
            self.log(f"WCS done: {summary.get('ok', 0)}/{summary.get('total', 0)} OK | {qc_text}")
        self.save_state()
        self.update_navigation_buttons()

    def validate_step(self) -> bool:
        return len(self.results) > 0

    def save_state(self):
        state_data = {
            "wcs_complete": len(self.results) > 0,
            "n_files": len(self.results),
            "use_cropped": self.use_cropped,
            "astap_exe": getattr(self.params.P, "astap_exe", "astap_cli.exe"),
            "astap_timeout_s": getattr(self.params.P, "astap_timeout_s", 120.0),
            "astap_search_radius_deg": getattr(self.params.P, "astap_search_radius_deg", 8.0),
            "astap_database": getattr(self.params.P, "astap_database", "D80"),
            "astap_annotate_variables": getattr(self.params.P, "astap_annotate_variables", False),
            "astap_fov_fudge": getattr(self.params.P, "astap_fov_fudge", 1.0),
            "astap_downsample_z": getattr(self.params.P, "astap_downsample_z", 2),
            "astap_max_stars_s": getattr(self.params.P, "astap_max_stars_s", 500),
            "wcs_max_workers": getattr(self.params.P, "wcs_max_workers", 1),
            "wcs_require_qc_pass": getattr(self.params.P, "wcs_require_qc_pass", True),
            "wcs_refine_enable": getattr(self.params.P, "wcs_refine_enable", True),
            "wcs_refine_max_match": getattr(self.params.P, "wcs_refine_max_match", 600),
            "wcs_refine_match_r_fwhm": getattr(self.params.P, "wcs_refine_match_r_fwhm", 1.6),
            "wcs_refine_min_match": getattr(self.params.P, "wcs_refine_min_match", 50),
            "gaia_radius_fudge": getattr(self.params.P, "gaia_radius_fudge", 1.35),
            "gaia_mag_max": getattr(self.params.P, "gaia_mag_max", 18.0),
            "gaia_wcs_mag_max": getattr(self.params.P, "gaia_wcs_mag_max", 18.0),
            "gaia_retry": getattr(self.params.P, "gaia_retry", 2),
            "gaia_timeout_s": getattr(self.params.P, "gaia_timeout_s", 30.0),
            "gaia_backoff_s": getattr(self.params.P, "gaia_backoff_s", 6.0),
            "gaia_allow_no_cache": getattr(self.params.P, "gaia_allow_no_cache", True),
            "internal_params": dict(getattr(self, "_internal_params", {}) or {}),
        }
        self.project_state.store_step_data("wcs_plate_solve", state_data)

    def restore_state(self):
        state_data = self.project_state.get_step_data("wcs_plate_solve")
        if state_data:
            for key, val in state_data.items():
                if hasattr(self.params.P, key):
                    setattr(self.params.P, key, val)
            internal_params = state_data.get("internal_params")
            if isinstance(internal_params, dict) and hasattr(self, "_internal_params"):
                for key, val in internal_params.items():
                    if key in self._internal_params:
                        self._internal_params[key] = val
        self._restore_success_results_from_summary()
        self.update_navigation_buttons()

    def closeEvent(self, event):
        workers = [
            (getattr(self, "worker", None), self.stop_wcs),
            (getattr(self, "astrometrynet_worker", None), self.stop_astrometrynet_solve),
            (getattr(self, "_internal_worker", None), self.stop_wcs_internal_solver),
        ]
        for worker, stop_fn in workers:
            if worker is not None and worker.isRunning():
                stop_fn()

        for worker, _ in workers:
            if worker is not None and worker.isRunning() and not worker.wait(10000):
                QMessageBox.warning(
                    self,
                    "Background Task Running",
                    "A WCS solver is still stopping. Please wait and close again.",
                )
                event.ignore()
                return
        super().closeEvent(event)
