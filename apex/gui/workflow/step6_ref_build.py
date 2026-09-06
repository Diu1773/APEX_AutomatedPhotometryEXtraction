"""
Step 6: Master Catalog Build (WCS-based reference catalog)

- Select reference frame using detection-based quality metrics
- Build a fixed master star list from the reference frame detections
- Write master catalogs for downstream steps
"""

from __future__ import annotations

import json
import hashlib
import re
import time
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd

from astropy.wcs import WCS

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QGroupBox, QMessageBox,
    QTextEdit, QFormLayout, QProgressBar, QDoubleSpinBox, QSpinBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView, QWidget,
    QDialog, QDialogButtonBox, QTabWidget, QCheckBox, QComboBox
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from .step_window_base import StepWindowBase
from .run_control import RunControlBar, format_duration, progress_status_text
from .param_dialog import ParamSpec, run_param_dialog, specs_from_map
from .log_panel import WorkflowLogWindow, WorkerStatusPanel, append_timestamped_log, show_raised
from .ui_helpers import (
    create_output_reuse_checkbox,
    create_parameter_button,
    set_table_row_background,
    status_row_background,
)
from apex.utils.step_paths import (
    step5_wcs_dir,
    step6_refbuild_dir,
    step4_dir,
)
from apex.utils.common_helpers import normalize_filter_key, safe_float as _safe_float
from apex.gui.layout_rules import scroll_wrap
from apex.utils.qc_utils import filter_files_by_qc, should_use_frame_quality_qc


_FILTER_RE = re.compile(r"[-_]([ugrizbvUGRIZBV])[-_.]", re.IGNORECASE)
_DATE_RE = re.compile(r"(20\d{6})")
_REF_SIGNATURE_FILE = "ref_build_signature.json"
_REF_SIGNATURE_VERSION = 2
_REF_SIGNATURE_PARAMS = (
    "wcs_require_qc_pass",
    "global_ref_filter",
    "ref_select_sat_pct",
    "ref_select_elong_pct",
    "ref_cat_max_sources",
    "ref_cat_min_sources",
    "ref_cat_max_elong",
    "ref_cat_max_abs_round",
    "ref_cat_sharp_min",
    "ref_cat_sharp_max",
    "ref_cat_min_peak_adu",
    "ref_wcs_match_radius_arcsec",
    "ref_wcs_min_match_rate",
    "ref_wcs_min_match_n",
    "ref_wcs_max_sep_med_arcsec",
    "ref_wcs_max_sep_p90_arcsec",
    "ref_wcs_max_dup_rate",
    "ref_per_date",
    "ref_master_union",
    "ref_union_min_frames",
    "idmatch_gaia_g_limit",
    "gaia_mag_max",
)


def _get_filter_from_filename(filename: str) -> Optional[str]:
    match = _FILTER_RE.search(str(filename))
    return normalize_filter_key(match.group(1)) if match else None


def _parse_date_key(value: str, params) -> Optional[str]:
    mode = str(getattr(params.P, "night_parse_mode", "regex") or "regex").strip().lower()
    if mode == "split":
        delim = str(getattr(params.P, "night_parse_split_delim", "_"))
        parts = value.split(delim) if delim else [value]
        idx = int(getattr(params.P, "night_parse_split_index", -1))
        if idx < 0:
            idx = len(parts) + idx
        if idx < 0 or idx >= len(parts):
            return None
        return parts[idx]
    if mode == "last_digits":
        n_digits = max(1, int(getattr(params.P, "night_parse_last_digits", 8)))
        m = re.search(rf"(\\d{{{n_digits}}})$", value)
        return m.group(1) if m else None
    try:
        pattern = str(getattr(params.P, "night_parse_regex", r".*_(\d{8})"))
        m = re.search(pattern, value)
    except re.error:
        return None
    if not m:
        return None
    if m.groupdict().get("date"):
        return m.group("date")
    if m.groups():
        return m.group(1)
    return m.group(0)


def _extract_date_key(filename: str, params=None) -> str:
    if params is None or not hasattr(params, "P"):
        match = _DATE_RE.search(str(filename))
        return match.group(1) if match else "unknown_date"
    date_key = None
    try:
        data_dir = Path(getattr(params.P, "data_dir", "."))
        file_path = Path(params.get_file_path(filename))
        if file_path.parent != data_dir:
            date_key = _parse_date_key(file_path.parent.name, params)
        if not date_key:
            date_key = _parse_date_key(file_path.name, params)
    except Exception:
        date_key = None
    if not date_key:
        date_key = _parse_date_key(str(filename), params)
    return date_key or "unknown_date"


# _safe_float imported from utils.common_helpers


class RefBuildWorker(QThread):
    progress = pyqtSignal(int, int, str)
    log = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str, str)

    def __init__(
        self,
        params,
        data_dir: Path,
        result_dir: Path,
        cache_dir: Path,
        file_list: List[str],
        ref_filter: str,
        sat_drop_pct: float,
        elong_drop_pct: float,
        ref_cat_max_sources: int,
        ref_cat_min_sources: int,
        ref_cat_max_elong: float,
        ref_cat_max_abs_round: float,
        ref_cat_sharp_min: float,
        ref_cat_sharp_max: float,
        ref_cat_min_peak_adu: float,
        wcs_match_radius_arcsec: float,
        wcs_min_match_rate: float,
        wcs_min_match_n: int,
        wcs_max_sep_med_arcsec: float,
        wcs_max_sep_p90_arcsec: float,
        wcs_max_dup_rate: float,
        ref_per_date: bool,
        gaia_mag_limit: float = 18.0,
        ref_master_union: bool = True,
        ref_union_min_frames: int = 1,
    ):
        super().__init__()
        self.params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(result_dir)
        self.cache_dir = Path(cache_dir)
        self.file_list = list(file_list)
        self.ref_filter = ref_filter
        self.sat_drop_pct = float(sat_drop_pct)
        self.elong_drop_pct = float(elong_drop_pct)
        self.ref_cat_max_sources = int(ref_cat_max_sources)
        self.ref_cat_min_sources = int(ref_cat_min_sources)
        self.ref_cat_max_elong = float(ref_cat_max_elong)
        self.ref_cat_max_abs_round = float(ref_cat_max_abs_round)
        self.ref_cat_sharp_min = float(ref_cat_sharp_min)
        self.ref_cat_sharp_max = float(ref_cat_sharp_max)
        self.ref_cat_min_peak_adu = float(ref_cat_min_peak_adu)
        self.wcs_match_radius_arcsec = float(wcs_match_radius_arcsec)
        self.wcs_min_match_rate = float(wcs_min_match_rate)
        self.wcs_min_match_n = int(wcs_min_match_n)
        self.wcs_max_sep_med_arcsec = float(wcs_max_sep_med_arcsec)
        self.wcs_max_sep_p90_arcsec = float(wcs_max_sep_p90_arcsec)
        self.wcs_max_dup_rate = float(wcs_max_dup_rate)
        self.ref_per_date = bool(ref_per_date)
        self.gaia_mag_limit = float(gaia_mag_limit)
        self.ref_master_union = bool(ref_master_union)
        self.ref_union_min_frames = max(1, int(ref_union_min_frames))
        self._stop_requested = False

    def stop(self):
        self._stop_requested = True

    def _log(self, msg: str):
        self.log.emit(msg)

    def run(self):
        """Run the master catalog build.

        Thin wrapper that delegates the compute to the Qt-free
        ``apex.analysis.refbuild.run_refbuild``, re-emitting its callbacks
        through the existing Qt signals. Behavior is identical to the prior
        inline ``_run_impl`` implementation.
        """
        from apex.analysis.refbuild import run_refbuild

        try:
            summary = run_refbuild(
                params=self.params,
                data_dir=self.data_dir,
                result_dir=self.result_dir,
                cache_dir=self.cache_dir,
                file_list=self.file_list,
                ref_filter=self.ref_filter,
                sat_drop_pct=self.sat_drop_pct,
                elong_drop_pct=self.elong_drop_pct,
                ref_cat_max_sources=self.ref_cat_max_sources,
                ref_cat_min_sources=self.ref_cat_min_sources,
                ref_cat_max_elong=self.ref_cat_max_elong,
                ref_cat_max_abs_round=self.ref_cat_max_abs_round,
                ref_cat_sharp_min=self.ref_cat_sharp_min,
                ref_cat_sharp_max=self.ref_cat_sharp_max,
                ref_cat_min_peak_adu=self.ref_cat_min_peak_adu,
                wcs_match_radius_arcsec=self.wcs_match_radius_arcsec,
                wcs_min_match_rate=self.wcs_min_match_rate,
                wcs_min_match_n=self.wcs_min_match_n,
                wcs_max_sep_med_arcsec=self.wcs_max_sep_med_arcsec,
                wcs_max_sep_p90_arcsec=self.wcs_max_sep_p90_arcsec,
                wcs_max_dup_rate=self.wcs_max_dup_rate,
                ref_per_date=self.ref_per_date,
                gaia_mag_limit=self.gaia_mag_limit,
                ref_master_union=self.ref_master_union,
                ref_union_min_frames=self.ref_union_min_frames,
                progress_cb=self.progress.emit,
                log_cb=self._log,
                error_cb=self.error.emit,
                should_stop=lambda: self._stop_requested,
            )
        except Exception as e:
            import traceback
            self._log(f"[ERROR] {e}\n{traceback.format_exc()}")
            self.error.emit("WORKER", str(e))
            self.finished.emit({})
            return
        self.finished.emit(summary if summary else {})


# The window names the rows it shows, in order. Everything else about each
# setting — where it lives in the config file, its type, its default, its label
# and range — is the map row, so a widget can only exist for a setting the
# loader builds and `save_toml` persists. Before 2026-08-16 this list repeated
# the type and default, and thirty settings across the app had a widget and no
# row: editing one showed "Parameters saved." and wrote nothing.
_STEP6_ATTRS: tuple[str, ...] = (
    "ref_select_sat_pct",
    "ref_select_elong_pct",
    "ref_per_date",
    "ref_master_union",
    "ref_union_min_frames",
    "ref_cat_max_sources",
    "ref_cat_min_sources",
    "ref_cat_max_elong",
    "ref_cat_max_abs_round",
    "ref_cat_sharp_min",
    "ref_cat_sharp_max",
    "ref_cat_min_peak_adu",
    "ref_wcs_match_radius_arcsec",
    "ref_wcs_min_match_rate",
    "ref_wcs_min_match_n",
    "ref_wcs_max_sep_med_arcsec",
    "ref_wcs_max_sep_p90_arcsec",
    "ref_wcs_max_dup_rate",
    "idmatch_gaia_g_limit",
)

_STEP6_SPECS: tuple[ParamSpec, ...] = specs_from_map(_STEP6_ATTRS)


class RefBuildWindow(StepWindowBase):
    """Step 6: Master Catalog Build (WCS-based)."""

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.worker = None
        self.log_window = None
        self.results = {}
        self._current_ref_signature: dict | None = None

        super().__init__(
            step_index=5,
            step_name="Master Catalog Build",
            params=params,
            project_state=project_state,
            main_window=main_window,
        )

        self.setup_step_ui()
        self.restore_state()

    def setup_step_ui(self):
        info = QLabel(
            "Build a fixed master catalog using WCS-solved frames.\n"
            "Selection prefers good WCS match stats, then saturation/elongation/FWHM."
        )
        info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; border-radius: 5px; }")
        info.setWordWrap(True)
        self.content_layout.addWidget(info)

        status_group = QGroupBox("WCS/Detection Status")
        status_layout = QVBoxLayout(status_group)
        self.status_label = QLabel("Checking...")
        status_layout.addWidget(self.status_label)
        self.content_layout.addWidget(status_group)

        control_layout = QHBoxLayout()

        btn_params = create_parameter_button("Master Catalog Parameters")
        btn_params.clicked.connect(self.open_parameters_dialog)
        control_layout.addWidget(btn_params)

        self.chk_use_existing_output = create_output_reuse_checkbox(
            not bool(getattr(self.params.P, "force_master_build", False)),
            "When enabled, Step 6 loads the existing reference-build summary and master catalogs "
            "if they are complete. Disable to rebuild the master catalog.",
        )
        control_layout.addWidget(self.chk_use_existing_output)

        self.run_bar = RunControlBar(
            "Run Master Catalog Build", "Show Log",
            run_cb=self.run_ref_build,
            stop_cb=self.stop_ref_build,
            log_cb=self.show_log_window,
        )
        control_layout.addWidget(self.run_bar)
        self.btn_run = self.run_bar.btn_run
        self.btn_stop = self.run_bar.btn_stop

        self.content_layout.addLayout(control_layout)

        progress_layout = QHBoxLayout()
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setValue(0)
        progress_layout.addWidget(self.progress_bar)

        self.progress_label = QLabel("Ready")
        self.progress_label.setMinimumWidth(300)
        progress_layout.addWidget(self.progress_label)
        self.content_layout.addLayout(progress_layout)

        self.tabs = QTabWidget()

        summary_tab = QWidget()
        summary_layout = QVBoxLayout(summary_tab)
        self.results_table = QTableWidget()
        self.results_table.setColumnCount(6)
        self.results_table.setHorizontalHeaderLabels([
            "Date", "Filter", "Ref Frame", "Sources", "FWHM (px)", "Sat Count"
        ])
        self.results_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.results_table.horizontalHeader().setStretchLastSection(True)
        self.results_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.results_table.setMinimumHeight(120)
        summary_layout.addWidget(self.results_table)
        self.tabs.addTab(summary_tab, "Summary")

        stats_tab = QWidget()
        stats_layout = QVBoxLayout(stats_tab)
        self.stats_table = QTableWidget()
        self.stats_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.stats_table.setMinimumHeight(220)
        stats_layout.addWidget(self.stats_table)
        self.tabs.addTab(stats_tab, "Stats")

        plot_tab = QWidget()
        plot_layout = QVBoxLayout(plot_tab)
        plot_controls = QHBoxLayout()
        plot_controls.addWidget(QLabel("Date:"))
        self.plot_date_combo = QComboBox()
        self.plot_date_combo.addItem("All")
        self.plot_date_combo.currentIndexChanged.connect(self._on_plot_date_changed)
        plot_controls.addWidget(self.plot_date_combo)
        plot_controls.addStretch()
        plot_layout.addLayout(plot_controls)
        self.plot_canvas = FigureCanvas(Figure(figsize=(8, 4)))
        self.plot_canvas.setMinimumHeight(260)
        plot_layout.addWidget(self.plot_canvas, 1)
        # Tallest page (318 px): scroll it so it does not set the window's
        # minimum height — see layout_rules.scroll_wrap.
        self.tabs.addTab(scroll_wrap(plot_tab), "Plot")

        self.content_layout.addWidget(self.tabs)

        _worker_group = QGroupBox("Workers")
        _worker_group.setMinimumWidth(300)
        _wg_layout = QVBoxLayout(_worker_group)
        _wg_layout.setContentsMargins(5, 5, 5, 5)
        self.worker_panel = WorkerStatusPanel(_worker_group)
        _wg_layout.addWidget(self.worker_panel)

        self.log_window = WorkflowLogWindow(
            self, "Master Catalog Build Log", width=900, height=500,
            side_widget=_worker_group,
        )
        self.log_text = self.log_window.log_text

        self.check_detection_status()

    def check_detection_status(self):
        cache_dir = Path(self.params.P.cache_dir)
        metas = list(cache_dir.glob("detect_*.json"))
        if not metas:
            step4_out = step4_dir(self.params.P.result_dir)
            metas = list(step4_out.glob("detect_*.json"))
        if not metas:
            self.status_label.setText("No detection cache found. Run Source Detection first.")
            self.status_label.setStyleSheet("color: red;")
            return
        self.status_label.setText(f"Detection cache found: {len(metas)} frames")
        self.status_label.setStyleSheet("color: green;")

    def open_parameters_dialog(self):
        run_param_dialog(
            self, "Master Catalog Build Parameters", _STEP6_SPECS,
            on_save=self.persist_params,
            overrides={"idmatch_gaia_g_limit": getattr(self.params.P, "idmatch_gaia_g_limit",
                        getattr(self.params.P, "gaia_mag_max", 18.0))},
            resize=(460, 560),
            info_text="Adjust reference catalog and WCS quality filter parameters. Changes apply to the next build run.",
        )

    @staticmethod
    def _signature_value(value):
        if isinstance(value, np.generic):
            value = value.item()
        if isinstance(value, Path):
            return str(value)
        if isinstance(value, float):
            return float(value) if np.isfinite(value) else None
        if isinstance(value, (bool, int, str)) or value is None:
            return value
        if isinstance(value, (list, tuple, set)):
            return [RefBuildWindow._signature_value(v) for v in value]
        if isinstance(value, dict):
            return {
                str(k): RefBuildWindow._signature_value(v)
                for k, v in sorted(value.items(), key=lambda item: str(item[0]))
            }
        try:
            if pd.isna(value):
                return None
        except Exception:
            pass
        return str(value)

    @staticmethod
    def _file_signature(path: Path | None) -> dict | None:
        if path is None:
            return None
        try:
            p = Path(path)
            if not p.exists():
                return None
            st = p.stat()
            try:
                path_text = str(p.resolve())
            except Exception:
                path_text = str(p)
            return {
                "path": path_text,
                "size": int(st.st_size),
                "mtime_ns": int(st.st_mtime_ns),
            }
        except Exception:
            return None

    @staticmethod
    def _first_existing(candidates: list[Path]) -> Path | None:
        for path in candidates:
            try:
                if path.exists() and path.stat().st_size > 0:
                    return path
            except Exception:
                continue
        return None

    def _build_ref_output_signature(self, files: list[str]) -> dict:
        result_dir = Path(self.params.P.result_dir)
        cache_dir = Path(self.params.P.cache_dir)
        s4_dir = step4_dir(result_dir)
        s5_dir = step5_wcs_dir(result_dir)

        frame_inputs = []
        for fname in files:
            detect_json = self._first_existing([
                cache_dir / f"detect_{fname}.json",
                s4_dir / f"detect_{fname}.json",
            ])
            detect_csv = self._first_existing([
                cache_dir / f"detect_{fname}.csv",
                s4_dir / f"detect_{fname}.csv",
            ])
            frame_inputs.append({
                "file": str(fname),
                "detect_json": self._file_signature(detect_json),
                "detect_csv": self._file_signature(detect_csv),
            })

        payload = {
            "signature_version": _REF_SIGNATURE_VERSION,
            "step": "step6_ref_build",
            "frames": [str(f) for f in files],
            "params": {
                k: self._signature_value(getattr(self.params.P, k, None))
                for k in _REF_SIGNATURE_PARAMS
            },
            "inputs": {
                "frame_quality": self._file_signature(s4_dir / "frame_quality.csv"),
                "wcs_summary": self._file_signature(s5_dir / "wcs_solve_summary.csv"),
                "gaia_fov": self._file_signature(s5_dir / "gaia_fov.ecsv"),
                "frames": frame_inputs,
            },
        }
        encoded = json.dumps(payload, sort_keys=True, ensure_ascii=False, allow_nan=False)
        payload["signature_hash"] = hashlib.sha1(encoded.encode("utf-8")).hexdigest()
        return payload

    def _stored_ref_signature(self) -> dict | None:
        sig_path = step6_refbuild_dir(self.params.P.result_dir) / _REF_SIGNATURE_FILE
        if not sig_path.exists():
            return None
        try:
            data = json.loads(sig_path.read_text(encoding="utf-8"))
            return data if isinstance(data, dict) else None
        except Exception:
            return None

    def _write_ref_signature(self, signature: dict) -> None:
        out_dir = step6_refbuild_dir(self.params.P.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / _REF_SIGNATURE_FILE).write_text(
            json.dumps(signature, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False),
            encoding="utf-8",
        )

    def _ref_signature_matches(self, signature: dict) -> tuple[bool, str]:
        stored = self._stored_ref_signature()
        if not stored:
            return False, "missing signature"
        if stored.get("signature_version") != _REF_SIGNATURE_VERSION:
            return False, "signature version mismatch"
        if stored.get("signature_hash") != signature.get("signature_hash"):
            return False, "signature hash mismatch"
        return True, "ok"

    def _current_ref_cache_status(self) -> tuple[bool, str, dict | None]:
        try:
            files = list(self.file_manager.get_file_list())
        except Exception as exc:
            return False, f"file list unavailable: {exc}", None
        if not files:
            return False, "no current frames", None
        use_qc = should_use_frame_quality_qc(
            Path(self.params.P.result_dir),
            self.params.P,
            "wcs_require_qc_pass",
            default=True,
        )
        files, _ = filter_files_by_qc(
            Path(self.params.P.result_dir),
            files,
            require_qc=use_qc,
        )
        if not files:
            return False, "no current frames after QC", None
        signature = self._build_ref_output_signature(list(files))
        matches, reason = self._ref_signature_matches(signature)
        if not matches:
            return False, reason, None
        summary = self._load_existing_summary()
        if not summary or int(summary.get("n_sources", 0) or 0) <= 0:
            return False, "reference output missing or empty", None
        return True, "ok", summary

    def _sync_cache_params(self) -> None:
        use_existing = bool(
            getattr(self, "chk_use_existing_output", None)
            and self.chk_use_existing_output.isChecked()
        )
        self.params.P.force_master_build = not use_existing
        if hasattr(self, "persist_params"):
            self.persist_params()

    def run_ref_build(self):
        if self.worker and self.worker.isRunning():
            return
        self._sync_cache_params()
        files = self.file_manager.get_file_list() if self.file_manager else []
        if not files:
            QMessageBox.warning(self, "Warning", "No frames found")
            return

        use_qc = should_use_frame_quality_qc(
            Path(self.params.P.result_dir),
            self.params.P,
            "wcs_require_qc_pass",
            default=True,
        )
        files, qc_info = filter_files_by_qc(Path(self.params.P.result_dir), files, require_qc=use_qc)
        if use_qc:
            if qc_info.get("applied"):
                self.log(f"[REF][QC] Frame QC filter: {qc_info['kept']}/{qc_info['total']} kept.")
            elif qc_info.get("path") is None:
                self.log("[REF][QC] frame_quality.csv not found; using all frames.")
            else:
                self.log(f"[REF][QC] frame_quality.csv ignored ({qc_info['reason']}); using all frames.")
        if not files:
            QMessageBox.warning(self, "Warning", "No frames after QC filter.")
            return

        ref_filter = normalize_filter_key(getattr(self.params.P, "global_ref_filter", ""))
        self.params.P.ref_select_sat_pct = float(getattr(self.params.P, "ref_select_sat_pct", 20.0))
        self.params.P.ref_select_elong_pct = float(getattr(self.params.P, "ref_select_elong_pct", 20.0))
        self.params.P.ref_cat_max_sources = int(getattr(self.params.P, "ref_cat_max_sources", 0))
        self.params.P.ref_cat_min_sources = int(getattr(self.params.P, "ref_cat_min_sources", 50))
        self.params.P.ref_cat_max_elong = float(getattr(self.params.P, "ref_cat_max_elong", 1.5))
        self.params.P.ref_cat_max_abs_round = float(getattr(self.params.P, "ref_cat_max_abs_round", 0.4))
        self.params.P.ref_cat_sharp_min = float(getattr(self.params.P, "ref_cat_sharp_min", 0.2))
        self.params.P.ref_cat_sharp_max = float(getattr(self.params.P, "ref_cat_sharp_max", 1.0))
        self.params.P.ref_cat_min_peak_adu = float(getattr(self.params.P, "ref_cat_min_peak_adu", 0.0))
        self.params.P.ref_wcs_match_radius_arcsec = float(getattr(self.params.P, "ref_wcs_match_radius_arcsec", 2.0))
        self.params.P.ref_wcs_min_match_rate = float(getattr(self.params.P, "ref_wcs_min_match_rate", 0.2))
        self.params.P.ref_wcs_min_match_n = int(getattr(self.params.P, "ref_wcs_min_match_n", 50))
        self.params.P.ref_wcs_max_sep_med_arcsec = float(getattr(self.params.P, "ref_wcs_max_sep_med_arcsec", 1.5))
        self.params.P.ref_wcs_max_sep_p90_arcsec = float(getattr(self.params.P, "ref_wcs_max_sep_p90_arcsec", 2.5))
        self.params.P.ref_wcs_max_dup_rate = float(getattr(self.params.P, "ref_wcs_max_dup_rate", 0.1))
        self.params.P.ref_per_date = bool(getattr(self.params.P, "ref_per_date", True))
        self.params.P.ref_master_union = bool(getattr(self.params.P, "ref_master_union", True))
        self.params.P.ref_union_min_frames = int(getattr(self.params.P, "ref_union_min_frames", 1))

        signature = self._build_ref_output_signature(files)
        self._current_ref_signature = signature
        if not bool(getattr(self.params.P, "force_master_build", False)):
            sig_ok, sig_reason = self._ref_signature_matches(signature)
            if sig_ok:
                cached_summary = self._load_existing_summary()
                if cached_summary:
                    self.log("[REF][CACHE] Existing Step 6 reference build matches current inputs; using cached output.")
                    self.results = cached_summary
                    self._update_results_table(cached_summary)
                    self._update_stats_table()
                    self._update_plot_tab(cached_summary)
                    self.save_state(cached_summary)
                    self.progress_label.setText("Cached Step 6 output loaded")
                    self.update_navigation_buttons()
                    self._current_ref_signature = None
                    return
            self.log(f"[REF][CACHE] Existing Step 6 output not reusable ({sig_reason}); rebuilding.")

        self.worker = RefBuildWorker(
            params=self.params,
            data_dir=self.params.P.data_dir,
            result_dir=self.params.P.result_dir,
            cache_dir=self.params.P.cache_dir,
            file_list=files,
            ref_filter=ref_filter,
            sat_drop_pct=self.params.P.ref_select_sat_pct,
            elong_drop_pct=self.params.P.ref_select_elong_pct,
            ref_cat_max_sources=self.params.P.ref_cat_max_sources,
            ref_cat_min_sources=self.params.P.ref_cat_min_sources,
            ref_cat_max_elong=self.params.P.ref_cat_max_elong,
            ref_cat_max_abs_round=self.params.P.ref_cat_max_abs_round,
            ref_cat_sharp_min=self.params.P.ref_cat_sharp_min,
            ref_cat_sharp_max=self.params.P.ref_cat_sharp_max,
            ref_cat_min_peak_adu=self.params.P.ref_cat_min_peak_adu,
            wcs_match_radius_arcsec=self.params.P.ref_wcs_match_radius_arcsec,
            wcs_min_match_rate=self.params.P.ref_wcs_min_match_rate,
            wcs_min_match_n=self.params.P.ref_wcs_min_match_n,
            wcs_max_sep_med_arcsec=self.params.P.ref_wcs_max_sep_med_arcsec,
            wcs_max_sep_p90_arcsec=self.params.P.ref_wcs_max_sep_p90_arcsec,
            wcs_max_dup_rate=self.params.P.ref_wcs_max_dup_rate,
            ref_per_date=self.params.P.ref_per_date,
            ref_master_union=self.params.P.ref_master_union,
            ref_union_min_frames=self.params.P.ref_union_min_frames,
            gaia_mag_limit=float(
                getattr(
                    self.params.P,
                    "idmatch_gaia_g_limit",
                    getattr(self.params.P, "gaia_mag_max", 18.0),
                )
            ),
        )
        self.worker.progress.connect(self.on_progress)
        self.worker.log.connect(self.log)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)

        self.run_bar.set_running(True)
        self.progress_bar.setValue(0)
        self.progress_bar.setMaximum(len(files))
        self._ref_start_time = time.monotonic()
        self.progress_label.setText(
            progress_status_text(0, len(files), self._ref_start_time, message="Starting...")
        )
        self.worker.start()
        self.show_log_window()

    def stop_ref_build(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()

    def on_progress(self, current, total, filename):
        self.progress_bar.setValue(current)
        if hasattr(self, "worker_panel"):
            pct = int(100 * current / max(1, total))
            self.worker_panel.update_worker(0, filename, f"{current}/{total}", pct)
        self.progress_label.setText(
            progress_status_text(current, total, getattr(self, "_ref_start_time", None), message=filename)
        )

    def on_finished(self, summary: dict):
        self.run_bar.set_running(False)
        elapsed_txt = ""
        if hasattr(self, "_ref_start_time"):
            elapsed_txt = f" | elapsed {format_duration(time.monotonic() - self._ref_start_time)}"
        self.progress_label.setText(f"Done{elapsed_txt}")
        if summary:
            self.results = summary
            self._update_results_table(summary)
            self._update_stats_table()
            self._update_plot_tab(summary)
            self.save_state(summary)
            if self._current_ref_signature:
                try:
                    self._write_ref_signature(self._current_ref_signature)
                    self.log("[REF][CACHE] Output signature saved for future reuse.")
                except Exception as exc:
                    self.log(f"[REF][CACHE] Signature write failed: {exc}")
        self._current_ref_signature = None
        self.update_navigation_buttons()

    def on_error(self, filename, error):
        self.log(f"ERROR {filename}: {error}")

    def _update_results_table(self, summary: dict):
        self.results_table.setRowCount(0)
        ref_frame = summary.get("ref_frame", "")
        ref_filter = summary.get("ref_filter", "")
        n_sources = summary.get("n_sources", 0)
        ref_frames_by_date = summary.get("ref_frames_by_date", {}) or {}

        metrics_path = step6_refbuild_dir(self.params.P.result_dir) / "ref_frame_stats.csv"
        metrics_df = None
        if metrics_path.exists():
            try:
                metrics_df = pd.read_csv(metrics_path)
            except Exception:
                metrics_df = None

        def _row_stats(fname: str):
            fwhm = "-"
            sat = "-"
            if metrics_df is not None and fname:
                row = metrics_df[metrics_df["file"] == fname]
                if not row.empty:
                    fwhm = f"{float(row.iloc[0].get('fwhm_px', np.nan)):.2f}" if np.isfinite(row.iloc[0].get('fwhm_px', np.nan)) else "-"
                    sat = str(int(row.iloc[0].get('sat_star_count', 0) or 0))
            return fwhm, sat

        def _has_sources(value) -> bool:
            try:
                return int(float(value)) > 0
            except (TypeError, ValueError):
                return False

        if ref_frames_by_date:
            for date_key, fname in sorted(ref_frames_by_date.items()):
                fwhm, sat = _row_stats(fname)
                row = self.results_table.rowCount()
                self.results_table.insertRow(row)
                self.results_table.setItem(row, 0, QTableWidgetItem(str(date_key)))
                self.results_table.setItem(row, 1, QTableWidgetItem(str(ref_filter)))
                self.results_table.setItem(row, 2, QTableWidgetItem(str(fname)))
                self.results_table.setItem(row, 3, QTableWidgetItem(str(n_sources)))
                self.results_table.setItem(row, 4, QTableWidgetItem(str(fwhm)))
                self.results_table.setItem(row, 5, QTableWidgetItem(str(sat)))
                set_table_row_background(
                    self.results_table,
                    row,
                    status_row_background(bool(fname) and _has_sources(n_sources)),
                )
        else:
            fwhm, sat = _row_stats(ref_frame)
            row = self.results_table.rowCount()
            self.results_table.insertRow(row)
            self.results_table.setItem(row, 0, QTableWidgetItem(_extract_date_key(ref_frame, self.params)))
            self.results_table.setItem(row, 1, QTableWidgetItem(str(ref_filter)))
            self.results_table.setItem(row, 2, QTableWidgetItem(str(ref_frame)))
            self.results_table.setItem(row, 3, QTableWidgetItem(str(n_sources)))
            self.results_table.setItem(row, 4, QTableWidgetItem(str(fwhm)))
            self.results_table.setItem(row, 5, QTableWidgetItem(str(sat)))
            set_table_row_background(
                self.results_table,
                row,
                status_row_background(bool(ref_frame) and _has_sources(n_sources)),
            )

    def _update_stats_table(self):
        stats_path = step6_refbuild_dir(self.params.P.result_dir) / "ref_frame_stats.csv"
        if not stats_path.exists():
            self.stats_table.setRowCount(0)
            self.stats_table.setColumnCount(0)
            return
        try:
            df = pd.read_csv(stats_path)
        except Exception:
            self.stats_table.setRowCount(0)
            self.stats_table.setColumnCount(0)
            return

        preferred_cols = [
            "date_key",
            "file",
            "filter",
            "wcs_ok",
            "match_rate_eff",
            "match_rate",
            "match_rate_cat",
            "n_match",
            "sep_med_arcsec",
            "sep_p90_arcsec",
            "dup_rate",
            "wcs_resid_med",
            "wcs_resid_med_px",
            "wcs_rms_px",
            "wcs_center_offset_arcsec",
            "wcs_match_n",
            "wcs_qc_pass",
            "wcs_qc_reason",
            "gaia_source",
            "fwhm_px",
            "sat_star_count",
            "n_sources",
            "selected",
            "selected_date",
        ]
        cols = [c for c in preferred_cols if c in df.columns]
        if not cols:
            cols = list(df.columns)
        df = df[cols].copy()

        self.stats_table.setColumnCount(len(cols))
        self.stats_table.setRowCount(len(df))
        self.stats_table.setHorizontalHeaderLabels(cols)
        self.stats_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.stats_table.horizontalHeader().setStretchLastSection(True)

        for r, (_, row) in enumerate(df.iterrows()):
            for c, col in enumerate(cols):
                val = row.get(col, "")
                if isinstance(val, float):
                    if "rate" in col or "sep_" in col:
                        text = f"{val:.3f}" if np.isfinite(val) else ""
                    else:
                        text = f"{val:.2f}" if np.isfinite(val) else ""
                else:
                    text = str(val)
                self.stats_table.setItem(r, c, QTableWidgetItem(text))

    def _update_plot_tab(self, summary: Optional[dict] = None) -> None:
        """The reference-build overview, drawn by the code the batch also uses.

        This was 202 lines here, which is why a headless run wrote
        `ref_frame_stats.csv` and no picture of it. The window keeps the date
        selector; the drawing moved next to the calculation.
        """
        if not hasattr(self, "plot_canvas") or self.plot_canvas is None:
            return

        from apex.analysis.refbuild_qc import draw_refbuild_overview

        self._refresh_plot_date_combo()
        date_key = (self.plot_date_combo.currentText()
                    if getattr(self, "plot_date_combo", None) is not None else "All")
        draw_refbuild_overview(self.plot_canvas.figure, self.params,
                               date_key=date_key, summary=summary)
        self.plot_canvas.draw_idle()

    def _refresh_plot_date_combo(self) -> None:
        """Keep the date selector in step with what the stats table holds."""
        combo = getattr(self, "plot_date_combo", None)
        if combo is None:
            return
        from apex.utils.step_paths import step6_refbuild_dir

        stats = Path(step6_refbuild_dir(self.params.P.result_dir)) / "ref_frame_stats.csv"
        if not stats.exists():
            return
        try:
            df = pd.read_csv(stats)
            if "date_key" not in df.columns:
                return
            dates = sorted(df["date_key"].fillna("unknown_date").astype(str).unique().tolist())
            items = ["All"] + dates
            prev = combo.currentText() if combo.count() else "All"
            combo.blockSignals(True)
            combo.clear()
            combo.addItems(items)
            combo.setCurrentText(prev) if prev in items else combo.setCurrentIndex(0)
            combo.blockSignals(False)
        except Exception:  # noqa: BLE001 - a stale selector must not stop the plot
            pass

    def _on_plot_date_changed(self, index: int) -> None:
        if index < 0:
            return
        summary = self.results if isinstance(self.results, dict) else None
        self._update_plot_tab(summary)

    def _load_existing_summary(self) -> Optional[dict]:
        out_dir = step6_refbuild_dir(self.params.P.result_dir)
        meta_path = out_dir / "ref_build_meta.json"
        if not meta_path.exists():
            return None
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
        except Exception as exc:
            self.log(f"[REF][CACHE] Existing ref_build_meta.json ignored: {exc}")
            return None

        ref_filter = str(meta.get("ref_filter", "") or "")
        catalog_candidates = [out_dir / "master_catalog.tsv"]
        if ref_filter:
            catalog_candidates.append(out_dir / f"ref_catalog_{ref_filter}.tsv")
        catalog_candidates.append(out_dir / "ref_catalog.tsv")
        catalog_path = next((p for p in catalog_candidates if p.exists()), None)
        if catalog_path is None:
            self.log("[REF][CACHE] ref_build_meta.json exists but no reference catalog was found; rebuilding.")
            return None

        n_sources = int(meta.get("n_ref_used") or 0)
        if n_sources <= 0:
            try:
                n_sources = len(pd.read_csv(catalog_path, sep="\t"))
            except Exception:
                n_sources = 0

        return {
            "ref_frame": meta.get("ref_frame", ""),
            "ref_filter": ref_filter,
            "n_sources": n_sources,
            "filters": meta.get("filters", []),
            "ref_per_date": bool(meta.get("ref_per_date", False)),
            "ref_frames_by_date": meta.get("ref_frames_by_date", {}) or {},
        }

    def validate_step(self) -> bool:
        valid, _, _ = self._current_ref_cache_status()
        return valid

    def save_state(self, summary: Optional[dict] = None):
        if summary is None:
            summary = self.results if isinstance(self.results, dict) else {}
        summary = summary or {}

        ref_frame = summary.get("ref_frame")
        ref_filter = summary.get("ref_filter")
        n_sources = summary.get("n_sources", 0)

        if not ref_frame or not ref_filter:
            meta_path = step6_refbuild_dir(self.params.P.result_dir) / "ref_build_meta.json"
            if meta_path.exists():
                try:
                    meta = json.loads(meta_path.read_text(encoding="utf-8"))
                    ref_frame = ref_frame or meta.get("ref_frame")
                    ref_filter = ref_filter or meta.get("ref_filter")
                except Exception:
                    pass

        if (not n_sources) and ref_filter:
            ref_path = step6_refbuild_dir(self.params.P.result_dir) / f"ref_catalog_{ref_filter}.tsv"
            if ref_path.exists():
                try:
                    n_sources = len(pd.read_csv(ref_path, sep="\t"))
                except Exception:
                    n_sources = 0

        state_data = {
            "ref_frame": ref_frame,
            "ref_filter": ref_filter,
            "n_sources": n_sources,
        }
        self.project_state.store_step_data("ref_build", state_data)
        if ref_frame:
            self.file_manager.ref_filename = ref_frame

    def restore_state(self):
        state = self.project_state.get_step_data("ref_build")
        valid, reason, summary = self._current_ref_cache_status()
        if valid and summary:
            self.results = summary
            ref_frame = summary.get("ref_frame")
            if ref_frame:
                self.file_manager.ref_filename = ref_frame
            self._update_results_table(summary)
            self._update_stats_table()
            self._update_plot_tab(summary)
            n_sources = summary.get("n_sources", 0)
            self.progress_label.setText(f"Loaded previous reference build ({n_sources} sources)")
            self.update_navigation_buttons()
            try:
                self.log("[REF][CACHE] Loaded previous Step 6 reference build from disk.")
            except Exception:
                pass
            return

        if state or self._load_existing_summary():
            try:
                self.log(f"[REF][CACHE] Previous Step 6 output not restored ({reason}).")
            except Exception:
                pass

    def log(self, msg: str):
        append_timestamped_log(self.log_text, msg)

    def show_log_window(self):
        show_raised(self.log_window)
