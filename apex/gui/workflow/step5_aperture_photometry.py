"""
Step 5 – Aperture Photometry  (detection-based, det_uid keyed)

Detection-first philosophy: photometry is performed at positions
from Step 4 detect_*.csv without requiring any WCS or ID matching.

Reads:
    step4_detection/detect_{fname}.csv   (det_uid, x, y)
    step4_detection/detect_{fname}.json  (FWHM metadata)
    <FITS image>

Apcorr (growth curve) writes to step5_aperture/:
    aperture_by_frame.csv   – optimal aperture per frame
    apcorr_summary.csv      – per-frame Apcorr summary
    growth_curve.csv        – full growth curve data

Photometry writes to step5_aperture/:
    photometry_{fname}.tsv  – det_uid keyed photometry results
    photometry_index.csv    – per-frame summary
"""
from __future__ import annotations

import json
import re
import traceback
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

import numpy as np
import pandas as pd
from astropy.io import fits

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QGroupBox, QMessageBox,
    QDialog, QFormLayout, QDialogButtonBox, QProgressBar,
    QCheckBox, QSpinBox, QDoubleSpinBox, QLineEdit, QWidget,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView, QTabWidget,
    QSplitter, QListWidget,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from apex.gui.workflow.step_window_base import StepWindowBase
from apex.gui.workflow.log_panel import WorkflowLogWindow, append_timestamped_log, show_raised
from apex.gui.workflow.step5_aperture_worker import ApertureWorker as Step5ApertureWorker  # noqa: F401
from apex.gui.workflow.param_dialog import ParamSpec, build_param_form, read_param_form
from apex.utils.step_paths import (
    step2_cropped_dir, step4_dir, step5_aperture_dir, crop_is_active,
)
from apex.utils.constants import get_parallel_workers
from apex.utils.photometry_utils import (
    circle_mask as _circle_mask,
    refine_local_centroid as _refine_local_centroid,
    phot_one_star as _phot_one_target,
)
from apex.utils.astro_utils import normalize_filter_name


# ── Scalar helpers ────────────────────────────────────────────────────────────

def _as_bool(x, default=False):
    if isinstance(x, bool):
        return x
    if x is None:
        return default
    if isinstance(x, (int, np.integer)):
        return bool(x)
    s = str(x).strip().lower()
    if s in ("1", "true", "t", "yes", "y", "on"):
        return True
    if s in ("0", "false", "f", "no", "n", "off"):
        return False
    return default


def _safe_float(x, default=np.nan):
    try:
        return float(x)
    except Exception:
        return default


def _to_float(val, default):
    try:
        if val is None:
            return float(default)
        out = float(val)
        return out if np.isfinite(out) else float(default)
    except Exception:
        return float(default)


def _to_int(val, default):
    try:
        if val is None:
            return int(default)
        return int(float(val))
    except Exception:
        return int(default)


# ── FITS helpers ───────────────────────────────────────────────────────────────

def _get_filter_lower(fits_path: Path):
    try:
        h = fits.getheader(fits_path)
        f = h.get("FILTER", None)
        if f is None:
            return "unknown"
        return normalize_filter_name(f)
    except Exception:
        return "unknown"


def _get_exptime(fits_path: Path, default=1.0):
    try:
        h = fits.getheader(fits_path)
        for k in ("EXPTIME", "EXPOSURE", "ITIME", "ELAPTIME"):
            if k in h:
                v = float(h[k])
                if np.isfinite(v) and v > 0:
                    return v
    except Exception:
        pass
    return float(default)


# ── Detect CSV / JSON helpers ──────────────────────────────────────────────────

def _load_detect_positions(fname: str, cache_dir: Path, result_dir: Path):
    """Load det_uid, x, y from detect_{fname}.csv.  Returns DataFrame or None."""
    candidates = [
        cache_dir / f"detect_{fname}.csv",
        step4_dir(result_dir) / f"detect_{fname}.csv",
        result_dir / f"detect_{fname}.csv",
    ]
    for p in candidates:
        if p.exists() and p.stat().st_size > 0:
            try:
                df = pd.read_csv(p)
                x_col = "x" if "x" in df.columns else ("xcenter" if "xcenter" in df.columns else None)
                y_col = "y" if "y" in df.columns else ("ycenter" if "ycenter" in df.columns else None)
                if x_col is None or y_col is None:
                    continue
                out = pd.DataFrame({"x": pd.to_numeric(df[x_col], errors="coerce"),
                                    "y": pd.to_numeric(df[y_col], errors="coerce")})
                if "det_uid" in df.columns:
                    out["det_uid"] = pd.to_numeric(df["det_uid"], errors="coerce").fillna(
                        pd.Series(range(len(df)))).astype(int)
                else:
                    out["det_uid"] = range(len(df))
                out = out.dropna(subset=["x", "y"])
                return out[["det_uid", "x", "y"]]
            except Exception:
                continue
    return None


def _load_fwhm_from_meta(fname: str, cache_dir: Path, result_dir: Path,
                          params_fwhm_guess=6.0) -> float:
    candidates = [
        cache_dir / f"detect_{fname}.json",
        step4_dir(result_dir) / f"detect_{fname}.json",
    ]
    for p in sorted([c for c in candidates if c.exists()],
                    key=lambda q: q.stat().st_mtime_ns, reverse=True):
        try:
            meta = json.loads(p.read_text(encoding="utf-8"))
        except Exception:
            continue
        for k in ("fwhm_med_rad_px", "fwhm_med_px", "fwhm_px", "fwhm_med"):
            v = meta.get(k, None)
            if v is not None:
                try:
                    v = float(v)
                    if np.isfinite(v) and v > 0:
                        return v
                except Exception:
                    continue
    return float(params_fwhm_guess)



# ── Step5PhotWorker – actual aperture photometry ──────────────────────────────

class Step5PhotWorker(QThread):
    """Per-frame aperture photometry at detected positions (det_uid keyed).

    Output TSV columns:
        det_uid, x_init, y_init, xcenter, ycenter, FILTER,
        r_ap_px, r_in_px, r_out_px,
        flux_net_adu, flux_e, mag_ap, mag_ap_err, snr, sky_med, sky_std,
        n_sky, peak_adu, flags_ap, exptime
    """
    progress = pyqtSignal(int, int, str)
    frame_done = pyqtSignal(str, dict)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str, str)
    log = pyqtSignal(str)

    # flags_ap bitmask
    FLAG_SAT = 1
    FLAG_NONLINEAR = 2
    FLAG_RECENTER_SHIFT = 4
    FLAG_LOW_SNR = 8

    def __init__(self, file_list, params, data_dir, result_dir, cache_dir, use_cropped=False):
        super().__init__()
        self.file_list = list(file_list)
        self.params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(result_dir)
        self.cache_dir = Path(cache_dir)
        self.use_cropped = use_cropped
        self.max_workers = get_parallel_workers(params)
        self._stop_requested = False
        self._write_lock = Lock()

    def stop(self):
        self._stop_requested = True

    def _log(self, msg):
        self.log.emit(msg)

    def _resolve_fits_path(self, fname: str) -> Path | None:
        if self.use_cropped and crop_is_active(self.result_dir):
            cdir = step2_cropped_dir(self.result_dir)
            cpath = cdir / fname
            if cpath.exists():
                return cpath
        fpath = self.data_dir / fname
        return fpath if fpath.exists() else None

    def run(self):  # noqa: C901
        try:
            P = self.params.P
            output_dir = step5_aperture_dir(self.result_dir)
            output_dir.mkdir(parents=True, exist_ok=True)

            GAIN = _to_float(getattr(P, "gain_e_per_adu", 1.0), 1.0)
            ZP = _to_float(getattr(P, "zp_initial", 25.0), 25.0)
            rn_e = _to_float(getattr(P, "rdnoise_e", 7.5), 7.5)
            ann_sigma = _to_float(getattr(P, "annulus_sigma_clip", 3.0), 3.0)
            ann_maxiter = _to_int(getattr(P, "fitsky_max_iter", 5), 5)
            sat_adu = _to_float(getattr(P, "saturation_adu", 60000.0), 60000.0)
            datamax_adu = _to_float(getattr(P, "datamax_adu", np.nan), np.nan)
            recenter = bool(getattr(P, "recenter_aperture", True))
            max_recenter_shift = _to_float(getattr(P, "max_recenter_shift", 2.0), 2.0)
            cbox_scale = _to_float(getattr(P, "center_cbox_scale", 1.5), 1.5)
            min_snr = _to_float(getattr(P, "min_snr_for_mag", 3.0), 3.0)
            sky_sigma_mode = str(getattr(P, "sky_sigma_mode", "local")).strip().lower()
            sky_sigma_incl_rn = bool(getattr(P, "sky_sigma_includes_rn", True))
            min_n_sky = _to_int(getattr(P, "sky_sigma_min_n_sky", 50), 50)
            force = bool(getattr(P, "force_rephot", False))
            fwhm_guess = _to_float(getattr(P, "fwhm_pix_guess", 6.0), 6.0)

            # Load aperture table
            ap_path = output_dir / "aperture_by_frame.csv"
            if not ap_path.exists():
                raise RuntimeError(
                    f"aperture_by_frame.csv not found at {ap_path}. "
                    "Run Apcorr first (Apcorr QC tab → Run Apcorr)."
                )
            df_ap = pd.read_csv(ap_path)
            df_ap["file"] = df_ap["file"].astype(str).map(lambda s: Path(s).name)

            # Load apcorr summary (may not exist if apcorr was disabled)
            apcorr_map: dict[str, float] = {}  # fname -> apcorr value (negative mag offset)
            apcorr_sum_path = output_dir / "apcorr_summary.csv"
            if apcorr_sum_path.exists():
                try:
                    df_apc = pd.read_csv(apcorr_sum_path)
                    df_apc["file"] = df_apc["file"].astype(str).map(lambda s: Path(s).name)
                    for _, row in df_apc.iterrows():
                        if bool(row.get("apply", False)) and np.isfinite(_to_float(row.get("apcorr"), np.nan)):
                            apcorr_map[str(row["file"])] = float(row["apcorr"])
                    self._log(f"[APCORR] Loaded correction for {len(apcorr_map)}/{len(df_apc)} frames")
                except Exception as e:
                    self._log(f"[APCORR] Failed to load apcorr_summary.csv: {e}")

            frames = list(self.file_list)
            total = len(frames)
            index_rows = []
            counters = {"processed": 0, "no_detect": 0, "no_aperture": 0, "stopped": 0}
            completed_count = [0]

            def process_frame(fname):
                if self._stop_requested:
                    return fname, None, "stopped"

                img_path = self._resolve_fits_path(fname)
                if img_path is None or not img_path.exists():
                    return fname, None, "no_fits"

                det_df = _load_detect_positions(fname, self.cache_dir, self.result_dir)
                if det_df is None or len(det_df) == 0:
                    return fname, None, "no_detect"

                ap_row = df_ap[df_ap["file"] == fname]
                frame_apcorr = apcorr_map.get(fname, np.nan)
                if ap_row.empty:
                    # Fallback: compute from FWHM
                    fwhm_used = _load_fwhm_from_meta(fname, self.cache_dir, self.result_dir, fwhm_guess)
                    ap_scale = _to_float(getattr(P, "phot_aperture_scale", 1.0), 1.0)
                    ann_in_scale = _to_float(getattr(P, "fitsky_annulus_scale", 4.0), 4.0)
                    ann_out_scale = _to_float(getattr(P, "fitsky_dannulus_scale", 2.0), 2.0)
                    ann_gap = _to_float(getattr(P, "annulus_min_gap_px", 6.0), 6.0)
                    r_ap = max(ap_scale * fwhm_used, 4.0)
                    r_in = max(ann_in_scale * fwhm_used, r_ap + ann_gap)
                    r_out = r_in + ann_out_scale * fwhm_used
                else:
                    r_ap = float(ap_row["r_ap"].values[0])
                    r_in = float(ap_row["r_in"].values[0])
                    r_out = float(ap_row["r_out"].values[0])
                    fwhm_used = _load_fwhm_from_meta(fname, self.cache_dir, self.result_dir, fwhm_guess)

                try:
                    img = fits.getdata(img_path).astype(float)
                except Exception as e:
                    return fname, None, f"fits_read_error:{e}"

                this_filter = _get_filter_lower(img_path)
                exptime = _get_exptime(img_path, default=1.0)
                h, w = img.shape
                cut_r = int(np.ceil(r_out + 2))

                phot_rows = []
                for _, row_d in det_df.iterrows():
                    det_uid = int(row_d["det_uid"])
                    x_init = float(row_d["x"])
                    y_init = float(row_d["y"])

                    # Centroid refinement
                    xi, yi = int(round(x_init)), int(round(y_init))
                    x0 = max(0, xi - cut_r)
                    x1 = min(w, xi + cut_r + 1)
                    y0 = max(0, yi - cut_r)
                    y1 = min(h, yi + cut_r + 1)
                    if (x1 - x0) < 5 or (y1 - y0) < 5:
                        continue

                    img_cut = img[y0:y1, x0:x1]
                    lxc = x_init - x0
                    lyc = y_init - y0
                    xcenter = x_init
                    ycenter = y_init
                    flags = 0

                    if recenter:
                        xc_new, yc_new = _refine_local_centroid(img_cut, lxc, lyc, fwhm_used, cbox_scale)
                        shift = np.sqrt((xc_new - lxc) ** 2 + (yc_new - lyc) ** 2)
                        if shift <= max_recenter_shift:
                            xcenter = x0 + xc_new
                            ycenter = y0 + yc_new
                            lxc, lyc = xc_new, yc_new
                        else:
                            flags |= self.FLAG_RECENTER_SHIFT

                    try:
                        (flux_e, sigma_e, snr, ap_sum, bkg_med, bkg_std,
                         ap_area, n_sky, is_sat, is_nonlinear) = _phot_one_target(
                            img_cut, lxc, lyc, r_ap, r_in, r_out,
                            sigma_clip_val=ann_sigma, maxiters=ann_maxiter,
                            gain=GAIN, rn_param_e=rn_e,
                            sky_sigma_mode=sky_sigma_mode,
                            sky_sigma_includes_rn=sky_sigma_incl_rn,
                            min_n_sky_for_local=min_n_sky,
                            sat_adu=sat_adu, datamax_adu=datamax_adu,
                        )
                    except Exception:
                        continue

                    if is_sat:
                        flags |= self.FLAG_SAT
                    if is_nonlinear:
                        flags |= self.FLAG_NONLINEAR

                    flux_net_adu = flux_e / GAIN if GAIN > 0 else np.nan

                    if np.isfinite(snr) and snr >= min_snr and np.isfinite(flux_e) and flux_e > 0:
                        mag_ap = ZP - 2.5 * np.log10(flux_e / exptime)
                        mag_ap_err = 2.5 / np.log(10) * sigma_e / flux_e if flux_e > 0 else np.nan
                    else:
                        mag_ap = np.nan
                        mag_ap_err = np.nan
                        flags |= self.FLAG_LOW_SNR

                    # Apply aperture correction (flux ratio: flux_corrected = flux * apcorr)
                    apcorr_val = frame_apcorr
                    if (np.isfinite(flux_e) and flux_e > 0
                            and np.isfinite(apcorr_val) and apcorr_val > 0):
                        flux_corr = flux_e * apcorr_val
                        mag_ap_corr = ZP - 2.5 * np.log10(flux_corr / exptime)
                    else:
                        mag_ap_corr = np.nan

                    peak_adu_val = float(np.nanmax(img_cut[_circle_mask(img_cut.shape, lxc, lyc, r_ap)])) \
                        if ap_area > 0 else np.nan

                    phot_rows.append({
                        "det_uid": det_uid,
                        "x_init": round(x_init, 4),
                        "y_init": round(y_init, 4),
                        "xcenter": round(xcenter, 4),
                        "ycenter": round(ycenter, 4),
                        "FILTER": this_filter,
                        "r_ap_px": round(r_ap, 3),
                        "r_in_px": round(r_in, 3),
                        "r_out_px": round(r_out, 3),
                        "flux_net_adu": round(flux_net_adu, 4) if np.isfinite(flux_net_adu) else np.nan,
                        "flux_e": round(flux_e, 4) if np.isfinite(flux_e) else np.nan,
                        "mag_ap": round(mag_ap, 6) if np.isfinite(mag_ap) else np.nan,
                        "mag_ap_err": round(mag_ap_err, 6) if np.isfinite(mag_ap_err) else np.nan,
                        "mag_ap_corr": round(mag_ap_corr, 6) if np.isfinite(mag_ap_corr) else np.nan,
                        "apcorr": round(apcorr_val, 6) if np.isfinite(apcorr_val) else np.nan,
                        "snr": round(snr, 3) if np.isfinite(snr) else np.nan,
                        "sky_med": round(bkg_med, 4) if np.isfinite(bkg_med) else np.nan,
                        "sky_std": round(bkg_std, 4) if np.isfinite(bkg_std) else np.nan,
                        "n_sky": int(n_sky),
                        "peak_adu": round(peak_adu_val, 2) if np.isfinite(peak_adu_val) else np.nan,
                        "flags_ap": int(flags),
                        "exptime": round(exptime, 4),
                    })

                if not phot_rows:
                    return fname, None, "no_phot"

                df_out = pd.DataFrame(phot_rows)
                out_tsv = output_dir / f"photometry_{fname}.tsv"
                with self._write_lock:
                    df_out.to_csv(out_tsv, sep="\t", index=False, encoding="utf-8-sig")

                n_rows = len(phot_rows)
                n_good = int(df_out["mag_ap"].notna().sum())
                n_fail = n_rows - n_good
                return fname, {
                    "file": fname, "filter": this_filter,
                    "n": n_rows, "n_goodmag": n_good, "n_fail": n_fail,
                    "r_ap": round(r_ap, 2),
                }, "processed"

            with ThreadPoolExecutor(max_workers=max(1, self.max_workers)) as ex:
                future_map = {ex.submit(process_frame, f): f for f in frames}
                for future in as_completed(future_map):
                    fname_done = future_map[future]
                    try:
                        fname_r, result_d, status = future.result()
                    except Exception as e:
                        self._log(f"Error {fname_done}: {e}")
                        continue

                    completed_count[0] += 1
                    self.progress.emit(completed_count[0], total, fname_r)

                    if status == "processed" and result_d:
                        counters["processed"] += 1
                        index_rows.append(result_d)
                        self.frame_done.emit(fname_r, result_d)
                        self._log(
                            f"[{completed_count[0]}/{total}] OK {fname_r} | "
                            f"f={result_d['filter']} n={result_d['n']} "
                            f"good={result_d['n_goodmag']} fail={result_d['n_fail']}"
                        )
                    elif status == "no_detect":
                        counters["no_detect"] += 1
                        self._log(f"[{completed_count[0]}/{total}] SKIP {fname_r} | reason=no detect csv")
                    elif status == "no_aperture":
                        counters["no_aperture"] += 1
                        self._log(f"[{completed_count[0]}/{total}] SKIP {fname_r} | reason=no aperture row")
                    elif status == "stopped":
                        counters["stopped"] += 1
                    else:
                        self._log(f"[{completed_count[0]}/{total}] SKIP {fname_r} | reason={status}")

            if index_rows:
                idx_path = output_dir / "photometry_index.csv"
                pd.DataFrame(index_rows).to_csv(idx_path, index=False)

            summary = {
                "frames": len(frames),
                "processed": counters["processed"],
                "no_detect": counters["no_detect"],
                "stopped": counters["stopped"],
            }
            self._log(
                f"Done | processed={counters['processed']} | "
                f"no_detect={counters['no_detect']} | stopped={counters['stopped']}"
            )
            self.finished.emit(summary)
        except Exception as e:
            err = f"{type(e).__name__}: {e}\n{traceback.format_exc()}"
            self.error.emit("WORKER", err)
            self.finished.emit({})


# ── Step5 parameter dialog specs ─────────────────────────────────────────────

_STEP5_MAIN_SPECS: tuple[ParamSpec, ...] = (
    ParamSpec("── Photometry ─────────────────────────", kind="sep"),
    ParamSpec("Zero Point", "zp_initial", "float", 10.0, 40.0, default=25.0),
    ParamSpec("Min SNR", "min_snr_for_mag", "float", 0.5, 20.0, default=3.0),
    ParamSpec("Recenter", "recenter_aperture", "bool", default=True),
    ParamSpec("── Aperture/Annulus scales ─────────────", kind="sep"),
    ParamSpec("Aperture (×FWHM)", "phot_aperture_scale", "float", 0.5, 5.0, 0.1, default=1.0),
    ParamSpec("Annulus inner (×FWHM)", "fitsky_annulus_scale", "float", 1.0, 10.0, 0.5, default=4.0),
    ParamSpec("Annulus width (×FWHM)", "fitsky_dannulus_scale", "float", 0.5, 10.0, 0.5, default=2.0),
    ParamSpec("Min annulus gap (px)", "annulus_min_gap_px", "float", 0.0, 50.0, 0.5, default=6.0),
    ParamSpec("Apcorr apply", "apcorr_apply", "bool", default=True),
    ParamSpec("Force", "force_rephot", "bool", default=False),
)

_STEP5_APCORR_SPECS: tuple[ParamSpec, ...] = (
    ParamSpec(
        "Max rel_scatter", "apcorr_scatter_max", "float", 0.01, 0.50, 0.01, 3,
        tooltip="Max allowed rel_scatter (1.4826×MAD/apcorr). 희소 시야는 0.10–0.20으로 올릴 것.",
        default=0.05,
    ),
    ParamSpec(
        "Large ref scale (×FWHM)", "apcorr_large_scale", "float", 1.5, 8.0, 0.5,
        tooltip="Large ref aperture (×FWHM). 낮출수록 PSF wing noise 감소.",
        write_also=("apcorr_large_ref_scale",), default=5.0,
    ),
    ParamSpec(
        "Min stars", "apcorr_use_min_n", "int", 3, 100,
        tooltip="apcorr 계산에 필요한 최소 기준성 수.",
        default=20,
    ),
)


# ── AperturePhotometryWindow (Step 5 GUI) ─────────────────────────────────────

class AperturePhotometryWindow(StepWindowBase):
    """Step 5 – Aperture Photometry (detection-based, det_uid keyed)."""

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.worker = None
        self.apcorr_worker = None
        self._pending_phot_after_apcorr = False
        self.file_list = []
        self.use_cropped = False
        self.log_window = None
        self._apcorr_summary_df = pd.DataFrame()
        self._growth_curve_df = pd.DataFrame()
        self._filter_cache = {}

        super().__init__(
            step_index=4,
            step_name="Aperture Photometry",
            params=params,
            project_state=project_state,
            main_window=main_window,
        )
        self.btn_complete.hide()
        self.setup_step_ui()
        self.restore_state()

    def setup_step_ui(self):
        info = QLabel(
            "Aperture photometry at detected positions (det_uid keyed). "
            "Run Apcorr first to find optimal aperture per frame, then Run Photometry."
        )
        info.setWordWrap(True)
        info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; border-radius: 5px; }")
        self.content_layout.addWidget(info)

        prog_layout = QHBoxLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setValue(0)
        prog_layout.addWidget(self.progress_bar)
        self.progress_label = QLabel("Ready")
        self.progress_label.setMinimumWidth(350)
        prog_layout.addWidget(self.progress_label)
        self.content_layout.addLayout(prog_layout)

        self.main_tabs = QTabWidget()
        self.content_layout.addWidget(self.main_tabs)

        # ── Tab 0: Apcorr QC ──────────────────────────────────────────────────
        apcorr_tab = QWidget()
        apcorr_layout = QVBoxLayout(apcorr_tab)

        apcorr_ctrl = QHBoxLayout()
        self.btn_apcorr_params = QPushButton("Parameters")
        self.btn_apcorr_params.setStyleSheet(
            "QPushButton { background-color: #795548; color: white; font-weight: bold; padding: 8px 15px; }"
        )
        self.btn_apcorr_params.clicked.connect(self.open_parameters_dialog)
        apcorr_ctrl.addWidget(self.btn_apcorr_params)
        apcorr_ctrl.addStretch()
        self.btn_run_apcorr = QPushButton("Run Apcorr")
        self.btn_run_apcorr.setStyleSheet(
            "QPushButton { background-color: #2E7D32; color: white; font-weight: bold; padding: 8px 20px; }"
        )
        self.btn_run_apcorr.clicked.connect(self.run_apcorr)
        apcorr_ctrl.addWidget(self.btn_run_apcorr)
        self.btn_stop_apcorr = QPushButton("Stop")
        self.btn_stop_apcorr.setStyleSheet(
            "QPushButton { background-color: #c62828; color: white; font-weight: bold; padding: 8px 15px; }"
        )
        self.btn_stop_apcorr.clicked.connect(self.stop_apcorr)
        self.btn_stop_apcorr.setEnabled(False)
        apcorr_ctrl.addWidget(self.btn_stop_apcorr)
        self.btn_log_apcorr = QPushButton("Log")
        self.btn_log_apcorr.setStyleSheet(
            "QPushButton { background-color: #607D8B; color: white; font-weight: bold; padding: 8px 15px; }"
        )
        self.btn_log_apcorr.clicked.connect(self.show_log_window)
        apcorr_ctrl.addWidget(self.btn_log_apcorr)
        apcorr_layout.addLayout(apcorr_ctrl)

        apcorr_note = QLabel(
            "Growth curve analysis: optimal aperture = radius with minimum magnitude error.\n"
            "Apcorr = correction from optimal aperture to large reference aperture."
        )
        apcorr_note.setStyleSheet("QLabel { background-color: #FFF8E1; padding: 8px; border-radius: 4px; }")
        apcorr_layout.addWidget(apcorr_note)

        apcorr_status_row = QHBoxLayout()
        self.apcorr_status = QLabel("growth curve: not loaded")
        apcorr_status_row.addWidget(self.apcorr_status)
        apcorr_status_row.addStretch()
        self.btn_refresh_apcorr = QPushButton("Refresh")
        self.btn_refresh_apcorr.setStyleSheet(
            "QPushButton { background-color: #607D8B; color: white; font-weight: bold; padding: 6px 12px; }"
        )
        self.btn_refresh_apcorr.clicked.connect(self.update_apcorr_tables)
        apcorr_status_row.addWidget(self.btn_refresh_apcorr)
        apcorr_layout.addLayout(apcorr_status_row)

        self.apcorr_splitter = QSplitter(Qt.Horizontal)

        plot_panel = QWidget()
        plot_layout = QVBoxLayout(plot_panel)
        self.apcorr_fig = Figure(figsize=(6.4, 5.6))
        self.apcorr_canvas = FigureCanvas(self.apcorr_fig)
        self.apcorr_ax_mag = self.apcorr_fig.add_subplot(211)
        self.apcorr_ax_err = self.apcorr_fig.add_subplot(212)
        self.apcorr_fig.subplots_adjust(hspace=0.35)
        plot_layout.addWidget(self.apcorr_canvas)
        self.apcorr_splitter.addWidget(plot_panel)

        side_panel = QWidget()
        side_layout = QVBoxLayout(side_panel)
        side_layout.addWidget(QLabel("Frames"))
        self.apcorr_file_list = QListWidget()
        self.apcorr_file_list.currentTextChanged.connect(self._on_apcorr_file_changed)
        side_layout.addWidget(self.apcorr_file_list, stretch=1)

        side_layout.addWidget(QLabel("Growth Curve (selected frame)"))
        self.apcorr_file_cand_table = QTableWidget()
        self.apcorr_file_cand_table.setColumnCount(7)
        self.apcorr_file_cand_table.setHorizontalHeaderLabels([
            "scale", "r (px)", "med_mag", "med_mag_err", "med_SNR", "n_stars", "selected",
        ])
        for col in range(7):
            mode = QHeaderView.Stretch if col in (0, 1) else QHeaderView.ResizeToContents
            self.apcorr_file_cand_table.horizontalHeader().setSectionResizeMode(col, mode)
        side_layout.addWidget(self.apcorr_file_cand_table, stretch=2)
        self.apcorr_splitter.addWidget(side_panel)
        self.apcorr_splitter.setStretchFactor(0, 3)
        self.apcorr_splitter.setStretchFactor(1, 2)
        apcorr_layout.addWidget(self.apcorr_splitter, stretch=1)

        apcorr_sum_group = QGroupBox("Growth Curve Summary")
        apcorr_sum_layout = QVBoxLayout(apcorr_sum_group)
        self.apcorr_table = QTableWidget()
        self.apcorr_table.setColumnCount(8)
        self.apcorr_table.setHorizontalHeaderLabels([
            "Frame", "FWHM(px)", "Opt r(px)", "Opt scale",
            "Apcorr", "SNR_opt", "mag_err_opt", "apply",
        ])
        self.apcorr_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        for col in range(1, 8):
            self.apcorr_table.horizontalHeader().setSectionResizeMode(col, QHeaderView.ResizeToContents)
        apcorr_sum_layout.addWidget(self.apcorr_table)
        apcorr_layout.addWidget(apcorr_sum_group)
        self.main_tabs.addTab(apcorr_tab, "Apcorr QC")

        # ── Tab 1: Photometry ─────────────────────────────────────────────────
        phot_tab = QWidget()
        phot_layout = QVBoxLayout(phot_tab)

        phot_ctrl = QHBoxLayout()
        self.btn_phot_params = QPushButton("Parameters")
        self.btn_phot_params.setStyleSheet(
            "QPushButton { background-color: #9C27B0; color: white; font-weight: bold; padding: 8px 15px; }"
        )
        self.btn_phot_params.clicked.connect(self.open_parameters_dialog)
        phot_ctrl.addWidget(self.btn_phot_params)
        phot_ctrl.addStretch()
        self.btn_run_phot = QPushButton("Run Photometry")
        self.btn_run_phot.setStyleSheet(
            "QPushButton { background-color: #4CAF50; color: white; font-weight: bold; padding: 8px 20px; }"
        )
        self.btn_run_phot.clicked.connect(self.run_photometry)
        phot_ctrl.addWidget(self.btn_run_phot)
        self.btn_stop_phot = QPushButton("Stop")
        self.btn_stop_phot.setStyleSheet(
            "QPushButton { background-color: #f44336; color: white; font-weight: bold; padding: 8px 15px; }"
        )
        self.btn_stop_phot.clicked.connect(self.stop_photometry)
        self.btn_stop_phot.setEnabled(False)
        phot_ctrl.addWidget(self.btn_stop_phot)
        self.btn_log_phot = QPushButton("Log")
        self.btn_log_phot.setStyleSheet(
            "QPushButton { background-color: #607D8B; color: white; font-weight: bold; padding: 8px 15px; }"
        )
        self.btn_log_phot.clicked.connect(self.show_log_window)
        phot_ctrl.addWidget(self.btn_log_phot)
        phot_layout.addLayout(phot_ctrl)

        frame_group = QGroupBox("Per-Frame Summary")
        frame_layout = QVBoxLayout(frame_group)
        self.frame_table = QTableWidget()
        self.frame_table.setColumnCount(6)
        self.frame_table.setHorizontalHeaderLabels(
            ["Frame", "Filter", "N_det", "N_goodmag", "N_fail", "r_ap (px)"]
        )
        self.frame_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.Stretch)
        for c in range(1, 6):
            self.frame_table.horizontalHeader().setSectionResizeMode(c, QHeaderView.ResizeToContents)
        self.frame_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        frame_layout.addWidget(self.frame_table)
        phot_layout.addWidget(frame_group)
        self.main_tabs.addTab(phot_tab, "Photometry")
        self.main_tabs.setCurrentIndex(0)

        # ── Log window ────────────────────────────────────────────────────────
        self.log_window = WorkflowLogWindow(
            self,
            "Step5 Aperture Photometry Log",
            width=800,
            height=400,
        )
        self.log_text = self.log_window.log_text

        self.populate_file_list()
        self.update_frame_table()
        self.update_apcorr_tables()

    # ── File list ─────────────────────────────────────────────────────────────

    def populate_file_list(self):
        crop_active = crop_is_active(self.params.P.result_dir)
        cropped_dir = step2_cropped_dir(self.params.P.result_dir)
        if crop_active and cropped_dir.exists() and list(cropped_dir.glob("*.fit*")):
            files = sorted([f.name for f in cropped_dir.glob("*.fit*")])
            self.use_cropped = True
        else:
            if not self.file_manager.filenames:
                try:
                    self.file_manager.scan_files()
                except Exception:
                    pass
            files = self.file_manager.filenames
            self.use_cropped = False
        self.file_list = list(files)

    # ── Actions ───────────────────────────────────────────────────────────────

    def run_apcorr(self):
        if not self.file_list:
            QMessageBox.warning(self, "Warning", "No files to process")
            return
        if self.worker and self.worker.isRunning():
            QMessageBox.warning(self, "Busy", "Photometry is running. Stop it first.")
            return
        if self.apcorr_worker and self.apcorr_worker.isRunning():
            return

        self.log_text.clear()
        self.log(f"Start Apcorr | {len(self.file_list)} frames")

        self.apcorr_worker = Step5ApertureWorker(
            self.file_list, self.params,
            self.params.P.data_dir, self.params.P.result_dir,
            self.params.P.cache_dir, self.use_cropped,
        )
        self.apcorr_worker.progress.connect(self.on_apcorr_progress)
        self.apcorr_worker.finished.connect(self.on_apcorr_finished)
        self.apcorr_worker.error.connect(self.on_apcorr_error)
        self.apcorr_worker.log.connect(self.log)

        self.btn_run_apcorr.setEnabled(False)
        self.btn_stop_apcorr.setEnabled(True)
        self.btn_run_phot.setEnabled(False)
        self.progress_bar.setMaximum(len(self.file_list))
        self.progress_bar.setValue(0)
        self.progress_label.setText(f"0/{len(self.file_list)} | Apcorr...")
        self.apcorr_worker.start()
        self.show_log_window()

    def stop_apcorr(self):
        if self.apcorr_worker and self.apcorr_worker.isRunning():
            self.apcorr_worker.stop()

    def run_photometry(self):
        if not self.file_list:
            QMessageBox.warning(self, "Warning", "No files to process")
            return
        if self.apcorr_worker and self.apcorr_worker.isRunning():
            QMessageBox.warning(self, "Busy", "Apcorr is running. Stop it first.")
            return
        if self.worker and self.worker.isRunning():
            return

        # Check prereq
        ap_path = step5_aperture_dir(self.params.P.result_dir) / "aperture_by_frame.csv"
        if not ap_path.exists():
            self.log("[AUTO] aperture_by_frame.csv missing – running Apcorr first.")
            self._pending_phot_after_apcorr = True
            self.run_apcorr()
            return

        files_to_run = list(self.file_list)
        apcorr_sum_path = step5_aperture_dir(self.params.P.result_dir) / "apcorr_summary.csv"
        if apcorr_sum_path.exists():
            try:
                df_apc = pd.read_csv(apcorr_sum_path)
                if (not df_apc.empty) and {"file", "apply"} <= set(df_apc.columns):
                    ok_vals = df_apc["apply"].astype(str).str.strip().str.lower().isin(
                        {"true", "1", "yes", "y", "on"}
                    )
                    ok_set = set(df_apc.loc[ok_vals, "file"].astype(str).map(lambda s: Path(str(s)).name))
                    before_n = len(files_to_run)
                    files_to_run = [f for f in files_to_run if f in ok_set]
                    self.log(f"[APCORR] apply=True frame filter: {len(files_to_run)}/{before_n} kept")
            except Exception as e:
                self.log(f"[APCORR] filter skipped ({e})")
        if not files_to_run:
            QMessageBox.warning(self, "Warning", "No frames with apcorr apply=True")
            return

        self._pending_phot_after_apcorr = False
        self.log_text.clear()
        self.log(f"Start photometry | {len(files_to_run)} frames | use_cropped={self.use_cropped}")

        self.worker = Step5PhotWorker(
            files_to_run, self.params,
            self.params.P.data_dir, self.params.P.result_dir,
            self.params.P.cache_dir, self.use_cropped,
        )
        self.worker.progress.connect(self.on_progress)
        self.worker.frame_done.connect(self.on_frame_done)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)
        self.worker.log.connect(self.log)

        self.btn_run_phot.setEnabled(False)
        self.btn_stop_phot.setEnabled(True)
        self.btn_run_apcorr.setEnabled(False)
        self.progress_bar.setMaximum(len(files_to_run))
        self.progress_bar.setValue(0)
        self.progress_label.setText(f"0/{len(files_to_run)} | Starting...")
        self.frame_table.setRowCount(0)
        self.worker.start()
        self.show_log_window()

    def stop_photometry(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()

    # ── Callbacks ─────────────────────────────────────────────────────────────

    def on_apcorr_progress(self, current, total, filename):
        self.progress_bar.setValue(current)
        self.progress_label.setText(f"{current}/{total} | Apcorr {filename}")

    def on_apcorr_finished(self, summary):
        self.btn_run_apcorr.setEnabled(True)
        self.btn_stop_apcorr.setEnabled(False)
        self.btn_run_phot.setEnabled(True)
        self.progress_label.setText("Apcorr done")
        self.log(f"Apcorr done: {summary}")
        self._cleanup_apcorr_worker()
        self.update_apcorr_tables()
        self.save_state()
        self.update_navigation_buttons()

        if self._pending_phot_after_apcorr:
            ap_path = step5_aperture_dir(self.params.P.result_dir) / "aperture_by_frame.csv"
            if ap_path.exists():
                self._pending_phot_after_apcorr = False
                self.run_photometry()
            else:
                self._pending_phot_after_apcorr = False
                QMessageBox.warning(self, "Apcorr Failed",
                                    "Apcorr finished but aperture_by_frame.csv not found.")

    def on_apcorr_error(self, src, err):
        self.log(f"APCORR ERROR {src}: {err}")

    def on_progress(self, current, total, filename):
        self.progress_bar.setValue(current)
        self.progress_label.setText(f"{current}/{total} | {filename}")

    def on_frame_done(self, filename, result):
        if not hasattr(self, "frame_table"):
            return
        r = self.frame_table.rowCount()
        self.frame_table.insertRow(r)
        self.frame_table.setItem(r, 0, QTableWidgetItem(str(result.get("file", filename))))
        self.frame_table.setItem(r, 1, QTableWidgetItem(str(result.get("filter", ""))))
        self.frame_table.setItem(r, 2, QTableWidgetItem(str(result.get("n", 0))))
        self.frame_table.setItem(r, 3, QTableWidgetItem(str(result.get("n_goodmag", 0))))
        self.frame_table.setItem(r, 4, QTableWidgetItem(str(result.get("n_fail", 0))))
        self.frame_table.setItem(r, 5, QTableWidgetItem(str(result.get("r_ap", ""))))
        self.frame_table.scrollToBottom()

    def on_finished(self, summary):
        self.btn_run_phot.setEnabled(True)
        self.btn_stop_phot.setEnabled(False)
        self.btn_run_apcorr.setEnabled(True)
        self.progress_label.setText("Done")
        self.log(f"Photometry done: {summary}")
        self._cleanup_worker()
        self.save_state()
        self.update_navigation_buttons()

    def on_error(self, src, err):
        self.log(f"ERROR {src}: {err}")

    # ── Thread cleanup ────────────────────────────────────────────────────────

    def _cleanup_worker(self, timeout_ms=5000):
        if not self.worker:
            return
        if self.worker.isRunning():
            self.worker.stop()
            self.worker.quit()
            self.worker.wait(timeout_ms)
        try:
            self.worker.deleteLater()
        except Exception:
            pass
        self.worker = None

    def _cleanup_apcorr_worker(self, timeout_ms=5000):
        if not self.apcorr_worker:
            return
        if self.apcorr_worker.isRunning():
            self.apcorr_worker.stop()
            self.apcorr_worker.quit()
            self.apcorr_worker.wait(timeout_ms)
        try:
            self.apcorr_worker.deleteLater()
        except Exception:
            pass
        self.apcorr_worker = None

    # ── Log ───────────────────────────────────────────────────────────────────

    def log(self, message: str):
        append_timestamped_log(self.log_text, message)

    def show_log_window(self):
        show_raised(self.log_window)

    # ── Table refresh ─────────────────────────────────────────────────────────

    def update_frame_table(self):
        idx_path = step5_aperture_dir(self.params.P.result_dir) / "photometry_index.csv"
        if not idx_path.exists() or not hasattr(self, "frame_table"):
            return
        try:
            idx = pd.read_csv(idx_path)
        except Exception:
            return
        self.frame_table.setRowCount(0)
        for _, row in idx.iterrows():
            r = self.frame_table.rowCount()
            self.frame_table.insertRow(r)
            self.frame_table.setItem(r, 0, QTableWidgetItem(str(row.get("file", ""))))
            self.frame_table.setItem(r, 1, QTableWidgetItem(str(row.get("filter", ""))))
            self.frame_table.setItem(r, 2, QTableWidgetItem(str(int(_safe_float(row.get("n", 0), 0)))))
            self.frame_table.setItem(r, 3, QTableWidgetItem(str(int(_safe_float(row.get("n_goodmag", 0), 0)))))
            self.frame_table.setItem(r, 4, QTableWidgetItem(str(int(_safe_float(row.get("n_fail", 0), 0)))))
            self.frame_table.setItem(r, 5, QTableWidgetItem(str(row.get("r_ap", ""))))

    def update_apcorr_tables(self):
        if not hasattr(self, "apcorr_table"):
            return
        base = step5_aperture_dir(self.params.P.result_dir)
        sum_path = base / "apcorr_summary.csv"
        gc_path = base / "growth_curve.csv"

        df_sum = pd.read_csv(sum_path) if sum_path.exists() else pd.DataFrame()
        df_gc = pd.read_csv(gc_path) if gc_path.exists() else pd.DataFrame()
        self._apcorr_summary_df = df_sum.copy()
        self._growth_curve_df = df_gc.copy()

        self.apcorr_table.setRowCount(0)
        for _, row in df_sum.iterrows():
            r = self.apcorr_table.rowCount()
            self.apcorr_table.insertRow(r)
            vals = [
                str(row.get("file", "")),
                f"{row.get('fwhm_used', np.nan):.2f}" if np.isfinite(_safe_float(row.get("fwhm_used"), np.nan)) else "",
                f"{row.get('r_optimal', np.nan):.2f}" if np.isfinite(_safe_float(row.get("r_optimal"), np.nan)) else "",
                f"{row.get('optimal_scale', np.nan):.2f}" if np.isfinite(_safe_float(row.get("optimal_scale"), np.nan)) else "",
                f"{row.get('apcorr', np.nan):.4f}" if np.isfinite(_safe_float(row.get("apcorr"), np.nan)) else "",
                f"{row.get('snr_optimal', np.nan):.1f}" if np.isfinite(_safe_float(row.get("snr_optimal"), np.nan)) else "",
                f"{row.get('mag_err_optimal', np.nan):.4f}" if np.isfinite(_safe_float(row.get("mag_err_optimal"), np.nan)) else "",
                str(row.get("apply", False)),
            ]
            for c, text in enumerate(vals):
                self.apcorr_table.setItem(r, c, QTableWidgetItem(text))

        if hasattr(self, "apcorr_status"):
            if not df_sum.empty:
                self.apcorr_status.setText(f"growth curve: {len(df_sum)} frames")
            else:
                self.apcorr_status.setText("growth curve: not found")

        if hasattr(self, "apcorr_file_list"):
            files = []
            if not df_gc.empty and "file" in df_gc.columns:
                files = sorted(pd.unique(df_gc["file"].astype(str)).tolist())
            elif not df_sum.empty and "file" in df_sum.columns:
                files = sorted(pd.unique(df_sum["file"].astype(str)).tolist())
            self.apcorr_file_list.blockSignals(True)
            self.apcorr_file_list.clear()
            if files:
                self.apcorr_file_list.addItems(files)
            self.apcorr_file_list.blockSignals(False)
            if files:
                self.apcorr_file_list.setCurrentRow(0)
                self._on_apcorr_file_changed(files[0])
            else:
                self._clear_apcorr_plot("No growth curve data")

    def _clear_apcorr_plot(self, message: str):
        for ax in (self.apcorr_ax_mag, self.apcorr_ax_err):
            ax.clear()
            ax.text(0.5, 0.5, message, ha="center", va="center",
                    transform=ax.transAxes, fontsize=9, color="#888")
            ax.set_axis_off()
        self.apcorr_canvas.draw_idle()

    def _on_apcorr_file_changed(self, filename: str):
        if not filename or not hasattr(self, "apcorr_file_cand_table"):
            return
        df_gc = getattr(self, "_growth_curve_df", None)
        if df_gc is None or df_gc.empty or "file" not in df_gc.columns:
            self.apcorr_file_cand_table.setRowCount(0)
            self._clear_apcorr_plot("No growth curve data")
            return

        sub = df_gc[df_gc["file"].astype(str) == str(filename)].copy()
        if sub.empty:
            self.apcorr_file_cand_table.setRowCount(0)
            self._clear_apcorr_plot("No data for this frame")
            return

        sub = sub.sort_values("r_px", kind="stable")
        self.apcorr_file_cand_table.setRowCount(0)
        for _, row in sub.iterrows():
            r = self.apcorr_file_cand_table.rowCount()
            self.apcorr_file_cand_table.insertRow(r)
            vals = [
                f"{_safe_float(row.get('scale'), np.nan):.2f}",
                f"{_safe_float(row.get('r_px'), np.nan):.2f}",
                f"{_safe_float(row.get('median_mag'), np.nan):.4f}",
                f"{_safe_float(row.get('median_mag_err'), np.nan):.4f}",
                f"{_safe_float(row.get('median_snr'), np.nan):.1f}",
                str(int(_safe_float(row.get("n_stars", 0), 0))),
                str(row.get("selected", False)),
            ]
            for c, text in enumerate(vals):
                self.apcorr_file_cand_table.setItem(r, c, QTableWidgetItem(text))

        # Growth curve plot
        ax_mag = self.apcorr_ax_mag
        ax_err = self.apcorr_ax_err
        ax_mag.clear()
        ax_err.clear()

        r_px = pd.to_numeric(sub["r_px"], errors="coerce").to_numpy(float)
        med_mag = pd.to_numeric(sub.get("median_mag", np.nan), errors="coerce").to_numpy(float)
        med_err = pd.to_numeric(sub.get("median_mag_err", np.nan), errors="coerce").to_numpy(float)
        sel_mask = sub.get("selected", False)
        if hasattr(sel_mask, "map"):
            sel_mask = sel_mask.map(lambda v: _as_bool(v, False)).to_numpy(bool)
        else:
            sel_mask = np.zeros(len(sub), dtype=bool)

        finite_mag = np.isfinite(r_px) & np.isfinite(med_mag)
        finite_err = np.isfinite(r_px) & np.isfinite(med_err)

        if not np.any(finite_mag) and not np.any(finite_err):
            self._clear_apcorr_plot("No finite growth curve points")
            return

        # Get summary info
        df_sum = getattr(self, "_apcorr_summary_df", pd.DataFrame())
        r_optimal = np.nan
        fwhm_used = np.nan
        apcorr_val = np.nan
        apply_on = False
        if not df_sum.empty and "file" in df_sum.columns:
            row_s_df = df_sum[df_sum["file"].astype(str) == str(filename)]
            if not row_s_df.empty:
                row_s = row_s_df.iloc[0]
                r_optimal = float(_safe_float(row_s.get("r_optimal"), np.nan))
                fwhm_used = float(_safe_float(row_s.get("fwhm_used"), np.nan))
                apcorr_val = float(_safe_float(row_s.get("apcorr"), np.nan))
                apply_on = _as_bool(row_s.get("apply", False), False)

        color = "#1565C0"
        if np.any(finite_mag):
            ax_mag.plot(r_px[finite_mag], med_mag[finite_mag], "-o", color=color,
                        linewidth=1.5, markersize=5, markeredgecolor="white", markeredgewidth=0.5)
        if np.isfinite(r_optimal):
            ax_mag.axvline(r_optimal, color="#E53935", linewidth=1.5, linestyle="--",
                           alpha=0.8, label=f"opt r={r_optimal:.1f}px")
        if np.isfinite(fwhm_used) and fwhm_used > 0:
            ax_mag.axvline(fwhm_used, color="#6D4C41", linewidth=1.2, linestyle=":",
                           alpha=0.8, label=f"FWHM={fwhm_used:.2f}px")
        ax_mag.invert_yaxis()
        ax_mag.set_ylabel("Inst Magnitude", fontsize=9)
        title = filename
        if np.isfinite(apcorr_val):
            title += f" | apcorr={apcorr_val:.4f}"
        title += f" | apply={'ON' if apply_on else 'OFF'}"
        ax_mag.set_title(title, fontsize=9)
        ax_mag.grid(True, alpha=0.25)
        ax_mag.legend(fontsize=8, frameon=False)

        if np.any(finite_err):
            ax_err.plot(r_px[finite_err], med_err[finite_err], "-s", color="#E53935",
                        linewidth=1.8, markersize=5, markeredgecolor="white", markeredgewidth=0.5)
        if np.isfinite(r_optimal):
            ax_err.axvline(r_optimal, color="#E53935", linewidth=1.5, linestyle="--", alpha=0.8)
        if np.isfinite(fwhm_used) and fwhm_used > 0:
            ax_err.axvline(fwhm_used, color="#6D4C41", linewidth=1.2, linestyle=":", alpha=0.8)
        ax_err.set_xlabel("Aperture radius (px)", fontsize=9)
        ax_err.set_ylabel("Median mag_err", fontsize=9)
        ax_err.set_title("Error vs Aperture (U-shape)", fontsize=9)
        ax_err.grid(True, alpha=0.25)

        self.apcorr_fig.tight_layout()
        self.apcorr_canvas.draw_idle()

    # ── Parameters dialog ─────────────────────────────────────────────────────

    def open_parameters_dialog(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Step5 Aperture Photometry Parameters")
        dialog.resize(460, 480)
        layout = QVBoxLayout(dialog)
        form = QFormLayout()

        form.addRow(QLabel("── Photometry ─────────────────────────"))

        self.param_zp = QDoubleSpinBox()
        self.param_zp.setRange(10.0, 40.0)
        self.param_zp.setValue(float(getattr(self.params.P, "zp_initial", 25.0)))
        form.addRow("Zero Point:", self.param_zp)

        self.param_snr = QDoubleSpinBox()
        self.param_snr.setRange(0.5, 20.0)
        self.param_snr.setValue(float(getattr(self.params.P, "min_snr_for_mag", 3.0)))
        form.addRow("Min SNR:", self.param_snr)

        self.param_recenter = QCheckBox("Enable")
        self.param_recenter.setChecked(bool(getattr(self.params.P, "recenter_aperture", True)))
        form.addRow("Recenter:", self.param_recenter)

        form.addRow(QLabel("── Aperture/Annulus scales ─────────────"))

        self.param_ap_scale = QDoubleSpinBox()
        self.param_ap_scale.setRange(0.5, 5.0)
        self.param_ap_scale.setSingleStep(0.1)
        self.param_ap_scale.setValue(float(getattr(self.params.P, "phot_aperture_scale", 1.0)))
        form.addRow("Aperture (×FWHM):", self.param_ap_scale)

        self.param_ann_in = QDoubleSpinBox()
        self.param_ann_in.setRange(1.0, 10.0)
        self.param_ann_in.setSingleStep(0.5)
        self.param_ann_in.setValue(float(getattr(self.params.P, "fitsky_annulus_scale", 4.0)))
        form.addRow("Annulus inner (×FWHM):", self.param_ann_in)

        self.param_ann_out = QDoubleSpinBox()
        self.param_ann_out.setRange(0.5, 10.0)
        self.param_ann_out.setSingleStep(0.5)
        self.param_ann_out.setValue(float(getattr(self.params.P, "fitsky_dannulus_scale", 2.0)))
        form.addRow("Annulus width (×FWHM):", self.param_ann_out)

        self.param_ann_gap = QDoubleSpinBox()
        self.param_ann_gap.setRange(0.0, 50.0)
        self.param_ann_gap.setSingleStep(0.5)
        self.param_ann_gap.setValue(float(getattr(self.params.P, "annulus_min_gap_px", 6.0)))
        form.addRow("Min annulus gap (px):", self.param_ann_gap)

        self.param_apcorr_apply = QCheckBox("Enable")
        self.param_apcorr_apply.setChecked(bool(getattr(self.params.P, "apcorr_apply", True)))
        form.addRow("Apcorr apply:", self.param_apcorr_apply)

        self.param_force = QCheckBox("Force re-phot")
        self.param_force.setChecked(bool(getattr(self.params.P, "force_rephot", False)))
        form.addRow("Force:", self.param_force)

        layout.addLayout(form)

        apcorr_group = QGroupBox("Aperture Correction (Apcorr)")
        apcorr_form = QFormLayout(apcorr_group)

        self.param_apcorr_scatter = QDoubleSpinBox()
        self.param_apcorr_scatter.setRange(0.01, 0.50)
        self.param_apcorr_scatter.setSingleStep(0.01)
        self.param_apcorr_scatter.setDecimals(3)
        self.param_apcorr_scatter.setValue(float(getattr(self.params.P, "apcorr_scatter_max", 0.05)))
        self.param_apcorr_scatter.setToolTip("Max allowed rel_scatter (1.4826×MAD/apcorr). 희소 시야는 0.10–0.20으로 올릴 것.")
        apcorr_form.addRow("Max rel_scatter:", self.param_apcorr_scatter)

        self.param_apcorr_scale = QDoubleSpinBox()
        self.param_apcorr_scale.setRange(1.5, 8.0)
        self.param_apcorr_scale.setSingleStep(0.5)
        self.param_apcorr_scale.setValue(float(
            getattr(self.params.P, "apcorr_large_scale", getattr(self.params.P, "apcorr_large_ref_scale", 5.0))
        ))
        self.param_apcorr_scale.setToolTip("Large ref aperture (×FWHM). 낮출수록 PSF wing noise 감소.")
        apcorr_form.addRow("Large ref scale (×FWHM):", self.param_apcorr_scale)

        self.param_apcorr_min_n = QSpinBox()
        self.param_apcorr_min_n.setRange(3, 100)
        self.param_apcorr_min_n.setValue(int(getattr(self.params.P, "apcorr_use_min_n", 20)))
        self.param_apcorr_min_n.setToolTip("apcorr 계산에 필요한 최소 기준성 수.")
        apcorr_form.addRow("Min stars:", self.param_apcorr_min_n)

        layout.addWidget(apcorr_group)

        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        buttons.accepted.connect(lambda: self._save_params(dialog))
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)
        dialog.exec_()

    def _save_params(self, dialog):
        self.params.P.zp_initial = self.param_zp.value()
        self.params.P.min_snr_for_mag = self.param_snr.value()
        self.params.P.recenter_aperture = self.param_recenter.isChecked()
        self.params.P.phot_aperture_scale = self.param_ap_scale.value()
        self.params.P.fitsky_annulus_scale = self.param_ann_in.value()
        self.params.P.fitsky_dannulus_scale = self.param_ann_out.value()
        self.params.P.annulus_min_gap_px = self.param_ann_gap.value()
        self.params.P.apcorr_apply = self.param_apcorr_apply.isChecked()
        self.params.P.force_rephot = self.param_force.isChecked()
        self.params.P.apcorr_scatter_max = self.param_apcorr_scatter.value()
        self.params.P.apcorr_large_scale = self.param_apcorr_scale.value()
        self.params.P.apcorr_large_ref_scale = self.param_apcorr_scale.value()
        self.params.P.apcorr_use_min_n = self.param_apcorr_min_n.value()
        self.save_state()
        saved = self.persist_params()
        msg = "Parameters saved." if saved else "Parameters updated, but TOML save failed."
        QMessageBox.information(dialog, "Saved", msg)
        dialog.accept()

    # ── Validation / State ────────────────────────────────────────────────────

    def update_navigation_buttons(self):
        self._sync_completion_from_outputs()
        super().update_navigation_buttons()
        self.btn_complete.hide()

    def _csv_has_rows(self, path: Path, required_columns=()) -> bool:
        if not path.exists():
            return False
        try:
            df = pd.read_csv(path)
        except Exception:
            return False
        if df.empty:
            return False
        return set(required_columns) <= set(df.columns)

    def _sync_completion_from_outputs(self):
        is_complete = self.validate_step()
        was_complete = self.project_state.is_step_completed(self.step_index)
        self.completed = is_complete
        if is_complete == was_complete:
            return

        if is_complete:
            self.project_state.mark_step_completed(self.step_index)
            if hasattr(self.main_window, "update_step_buttons"):
                self.main_window.update_step_buttons()
            if hasattr(self.main_window, "append_log"):
                self.main_window.append_log(
                    f"✓ Step {self.step_index + 1} completed: {self.step_name}"
                )
            return

        self.project_state.mark_step_incomplete(self.step_index)
        if hasattr(self.main_window, "update_step_buttons"):
            self.main_window.update_step_buttons()

    def validate_step(self) -> bool:
        base = step5_aperture_dir(self.params.P.result_dir)
        apcorr_sum_path = base / "apcorr_summary.csv"
        idx_path = base / "photometry_index.csv"
        return (
            self._csv_has_rows(apcorr_sum_path, required_columns=("file",))
            and self._csv_has_rows(idx_path, required_columns=("file",))
        )

    def save_state(self):
        state = {
            "zp_initial": getattr(self.params.P, "zp_initial", 25.0),
            "min_snr_for_mag": getattr(self.params.P, "min_snr_for_mag", 3.0),
            "recenter_aperture": getattr(self.params.P, "recenter_aperture", True),
            "phot_aperture_scale": getattr(self.params.P, "phot_aperture_scale", 1.0),
            "fitsky_annulus_scale": getattr(self.params.P, "fitsky_annulus_scale", 4.0),
            "fitsky_dannulus_scale": getattr(self.params.P, "fitsky_dannulus_scale", 2.0),
            "annulus_min_gap_px": getattr(self.params.P, "annulus_min_gap_px", 6.0),
            "apcorr_apply": getattr(self.params.P, "apcorr_apply", True),
            "force_rephot": getattr(self.params.P, "force_rephot", False),
        }
        self.project_state.store_step_data("aperture_photometry", state)

    def restore_state(self):
        state = self.project_state.get_step_data("aperture_photometry")
        if state:
            for k, v in state.items():
                if hasattr(self.params.P, k):
                    setattr(self.params.P, k, v)

    def closeEvent(self, event):
        if self.worker and self.worker.isRunning():
            self.stop_photometry()
        if self.apcorr_worker and self.apcorr_worker.isRunning():
            self.stop_apcorr()
        self._cleanup_worker(timeout_ms=10000)
        self._cleanup_apcorr_worker(timeout_ms=10000)
        super().closeEvent(event)
