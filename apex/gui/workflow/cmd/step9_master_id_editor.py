"""
Step 9: Master ID Editor
"""

from __future__ import annotations

import time
import json
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
import astropy.units as u
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QGroupBox, QMessageBox,
    QTextEdit, QDialog, QFormLayout, QDialogButtonBox, QProgressBar,
    QCheckBox, QSpinBox, QDoubleSpinBox, QLineEdit, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView, QWidget, QComboBox,
    QColorDialog, QFrame
)
from PyQt5.QtGui import QColor
from PyQt5.QtCore import Qt, QPoint

from apex.gui.widgets.fits_viewer import FITSViewerWidget, OverlayMarker

from apex.gui.workflow.step_window_base import StepWindowBase
from apex.gui.workflow.ui_helpers import create_parameter_button, build_scroll_param_dialog
from apex.utils.step_paths import (
    crop_is_active,
    step2_cropped_dir,
    step7_forced_phot_dir,
    step5_wcs_dir,
    step6_refbuild_dir,
)
from apex.utils.step_paths_cmd import step8_psf_dir, step9_selection_dir
from apex.utils.io_utils import (
    parse_int64_scalar,
    parse_int64_series,
    read_csv_int64_source_id,
    read_ecsv_int64_source_id,
)


class MasterIdEditorWindow(StepWindowBase):
    """Step 9: Master ID Editor"""

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.file_list = []
        self.use_cropped = False
        self.current_filename = None
        self.image_data = None
        self.header = None
        self._file_filter_map = {}
        self.phot_df = None
        self.master_ids = set()
        self._master_ids_persisted = False
        self.selected_source_id = None
        self.selected_phot_row = None
        self.selected_phot_file = None
        self.last_click_xy = None
        self.gaia_df = None  # Gaia catalog cache
        self._gaia_gmag_map = {}
        self.membership_by_source = {}
        self.membership_by_gaia = {}
        self.membership_by_id = {}
        self._membership_loaded = False
        self._membership_log_once = False
        self.master_gmag_map = {}
        # Fixed ID map (source_id -> display ID); IDs must stay stable across sessions.
        self.internal_id_map = {}
        self.source_id_from_internal = {}
        self._global_id_map = {}
        self._auto_master_dirty = False

        # GPU-accelerated FITS viewer
        self._viewer: FITSViewerWidget | None = None

        self.hover_xy = None  # Track mouse hover position for G key

        # PSF iter2 source tracking (det_uid < 0)
        self.psf_iter2_ids: set = set()

        # CMD ROI (circle in image pixel coords, saved to cmd_roi.json)
        self._roi_circle: dict | None = None   # {ra_deg, dec_deg, radius_arcsec}
        self._roi_preview: tuple | None = None # (cx, cy, r_px) live preview while dragging
        self._roi_mode = False                 # draw-mode toggle
        self._roi_drag_start: tuple | None = None

        # Overlay color customization (simplified: 3 categories + selected)
        self._overlay_colors: dict = {
            "gaia":     "#00FF00",   # step4 detected + Gaia + in master
            "local":    "#FFC107",   # step4 detected + non-Gaia + in master
            "forced":   "#00BFFF",   # forced phot only (step4 미검출) + in master
            "removed":  "#9E9E9E",   # not in master
            "selected": "#FF0000",
        }
        self._color_window: "QWidget | None" = None
        self._color_btns: dict = {}
        self._overlay_visible: dict = {k: True for k in (
            "gaia", "local", "forced", "removed", "selected"
        )}
        self._show_members: bool = False  # toggle for member center dot

        # Frame data caches: filename -> data  (LRU, max _FITS_CACHE_SIZE entries)
        self._fits_cache: dict = {}        # filename -> (image_data, header)
        self._fits_cache_order: list = []  # LRU order
        self._FITS_CACHE_SIZE = max(3, int(getattr(params.P, "step8_fits_cache_size", 8)))
        self._phot_cache: dict = {}     # filename -> phot_df (all frames, small)
        self._phot_arr_cache: dict = {} # filename -> (x, y, sids) numpy arrays
        self._prefetch_lock = threading.Lock()
        self._prefetch_pending: set[str] = set()
        self._prefetch_executor = ThreadPoolExecutor(max_workers=1)

        super().__init__(
            step_index=8,
            step_name="Master ID Editor",
            params=params,
            project_state=project_state,
            main_window=main_window
        )

        self.setup_step_ui()
        self.restore_state()
        self.setFocusPolicy(Qt.StrongFocus)

    def setup_step_ui(self):
        info = QLabel(
            "Inspect forced-photometry results from step7. Click a star to see info.\n"
            "Shortcuts: G=Radial profile (at cursor), [ / ]=Prev/Next frame\n"
            "Green=Gaia in master, Yellow=Local, Cyan=forced-only. Members ● = Gaia member dot."
        )
        info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; border-radius: 5px; }")
        self.content_layout.addWidget(info)

        control_layout = QHBoxLayout()
        btn_params = create_parameter_button("Editor Parameters")
        btn_params.clicked.connect(self.open_parameters_dialog)
        control_layout.addWidget(btn_params)

        btn_log = QPushButton("Log")
        btn_log.setStyleSheet("QPushButton { background-color: #607D8B; color: white; font-weight: bold; padding: 8px 15px; }")
        btn_log.clicked.connect(self.show_log_window)
        control_layout.addWidget(btn_log)

        self.btn_colors = QPushButton("🎨 Colors")
        self.btn_colors.setCheckable(True)
        self.btn_colors.setStyleSheet(
            "QPushButton { background-color: #455A64; color: white; font-weight: bold; padding: 8px 12px; }"
            "QPushButton:checked { background-color: #78909C; }"
        )
        self.btn_colors.toggled.connect(self._toggle_color_panel)
        control_layout.addWidget(self.btn_colors)

        self.chk_members = QCheckBox("Members ●")
        self.chk_members.setToolTip("Show center dot on Gaia member stars (Pmem ≥ threshold)")
        self.chk_members.toggled.connect(self._on_members_toggled)
        control_layout.addWidget(self.chk_members)

        self.content_layout.addLayout(control_layout)

        # Selected info
        select_layout = QHBoxLayout()
        self.selected_label = QLabel("Selected source_id: (none)")
        select_layout.addWidget(self.selected_label)
        select_layout.addStretch()
        self.content_layout.addLayout(select_layout)

        # Viewer + table
        main_layout = QHBoxLayout()

        viewer_group = QGroupBox("Preview")
        viewer_layout = QVBoxLayout(viewer_group)

        file_layout = QHBoxLayout()
        file_layout.addWidget(QLabel("File:"))
        self.file_combo = QComboBox()
        self.file_combo.currentIndexChanged.connect(self.on_file_changed)
        file_layout.addWidget(self.file_combo)
        viewer_layout.addLayout(file_layout)


        # CMD ROI controls
        roi_layout = QHBoxLayout()
        self.btn_set_roi = QPushButton("Set CMD ROI")
        self.btn_set_roi.setCheckable(True)
        self.btn_set_roi.setToolTip("Click to enter ROI draw mode, then click+drag on the image to define a circle.\nThis ROI filters sources shown in the CMD (ZP/photometry are unaffected).")
        self.btn_set_roi.setStyleSheet(
            "QPushButton { background-color: #37474F; color: white; font-weight: bold; padding: 4px 10px; }"
            "QPushButton:checked { background-color: #00BCD4; color: black; }"
        )
        self.btn_set_roi.toggled.connect(self._on_set_roi_toggled)
        roi_layout.addWidget(self.btn_set_roi)

        self.btn_clear_roi = QPushButton("Clear ROI")
        self.btn_clear_roi.setStyleSheet("QPushButton { background-color: #546E7A; color: white; padding: 4px 10px; }")
        self.btn_clear_roi.clicked.connect(self._on_clear_roi)
        roi_layout.addWidget(self.btn_clear_roi)

        self.roi_info_label = QLabel("No ROI set")
        self.roi_info_label.setStyleSheet("QLabel { color: #90A4AE; font-size: 9pt; }")
        roi_layout.addWidget(self.roi_info_label)
        roi_layout.addStretch()
        viewer_layout.addLayout(roi_layout)

        self._viewer = FITSViewerWidget(self)
        self._viewer.gl.mouse_pressed.connect(self._on_viewer_press)
        self._viewer.gl.mouse_released.connect(self._on_viewer_release)
        self._viewer.gl.mouse_moved.connect(self._on_viewer_move)
        viewer_layout.addWidget(self._viewer)
        main_layout.addWidget(viewer_group, stretch=2)

        table_group = QGroupBox("Master IDs")
        table_layout = QVBoxLayout(table_group)

        self.master_table = QTableWidget()
        self.master_table.setColumnCount(4)
        self.master_table.setHorizontalHeaderLabels(["ID", "source_id", "G mag", "Pmem"])
        self.master_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.master_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.master_table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeToContents)
        self.master_table.horizontalHeader().setSectionResizeMode(3, QHeaderView.ResizeToContents)
        self.master_table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.master_table.itemSelectionChanged.connect(self.on_table_selection_changed)
        table_layout.addWidget(self.master_table)

        main_layout.addWidget(table_group, stretch=1)

        self.content_layout.addLayout(main_layout)

        # Floating color-legend window (built now, shown on demand)
        self._color_window = self._build_color_window()

        self.log_window = QWidget(self, Qt.Window)
        self.log_window.setWindowTitle("Master ID Log")
        self.log_window.resize(700, 350)
        log_layout = QVBoxLayout(self.log_window)
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setStyleSheet("QTextEdit { font-family: monospace; font-size: 9pt; }")
        log_layout.addWidget(self.log_text)

        self.populate_file_list()
        self.load_master_ids()
        self.load_psf_new_sources()
        self.load_gaia_catalog()
        self._load_roi()

    def log(self, message: str):
        timestamp = time.strftime("%H:%M:%S")
        self.log_text.append(f"[{timestamp}] {message}")

    def populate_file_list(self):
        self._file_filter_map = {}
        cropped_dir = step2_cropped_dir(self.params.P.result_dir)
        crop_active = crop_is_active(self.params.P.result_dir)
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

        base_count = len(files)
        # Hide unsolved frames: keep only rows with wcs_ok=True when available.
        forced_dir = step7_forced_phot_dir(self.params.P.result_dir)
        stats_path = next(
            (
                p for p in (
                    forced_dir / "frame_stats.csv",
                    forced_dir / "step8_frame_stats.csv",  # legacy
                )
                if p.exists()
            ),
            forced_dir / "frame_stats.csv",
        )
        if stats_path.exists():
            try:
                s7 = pd.read_csv(stats_path)
            except Exception:
                s7 = pd.DataFrame()
            if (not s7.empty) and ("file" in s7.columns):
                if "filter" in s7.columns:
                    for fname, filt in zip(s7["file"].astype(str), s7["filter"].astype(str)):
                        self._file_filter_map[fname] = self._normalize_filter_key(filt)
                if "wcs_ok" in s7.columns:
                    wcs_ok = s7["wcs_ok"].astype(str).str.strip().str.lower().isin({"true", "1", "yes", "y", "on"})
                else:
                    wcs_ok = pd.Series(True, index=s7.index, dtype=bool)
                keep = set(s7.loc[wcs_ok, "file"].astype(str).tolist())
                files = [f for f in files if f in keep]
                self.log(f"Step 9 frame filter (wcs_ok): {len(files)}/{base_count} kept")

        # Filter: keep frames that have forced phot output.
        forced_dir = step7_forced_phot_dir(self.params.P.result_dir)
        if forced_dir.exists():
            before_idm = len(files)
            missing_forced = []
            kept_files = []
            for f in files:
                if self._resolve_forced_phot_path(f) is not None:
                    kept_files.append(f)
                else:
                    missing_forced.append(f)
            files = kept_files
            self.log(f"Step 9 frame filter (forced phot tsv): {len(files)}/{before_idm} kept")
            if missing_forced:
                sample = ", ".join(missing_forced[:5])
                suffix = "..." if len(missing_forced) > 5 else ""
                self.log(f"Step 9 skipped frames without forced TSV: {sample}{suffix}")

        self.file_list = list(files)
        self.file_combo.clear()
        self.file_combo.addItems(self.file_list)
        self._load_filter_map_from_index()
        with self._prefetch_lock:
            self._fits_cache.clear()
            self._fits_cache_order.clear()
            self._phot_arr_cache.clear()
            self._prefetch_pending.clear()

    @staticmethod
    def _normalize_filter_key(value: str | None) -> str:
        if value is None:
            return ""
        return str(value).strip().upper()

    def _extract_filter_from_header(self, header) -> str:
        if header is None:
            return ""
        for key in ("FILTER", "FILTER1", "FILTER2", "FILTNAM"):
            val = header.get(key)
            if val:
                return self._normalize_filter_key(val)
        return ""

    @staticmethod
    def _infer_filter_from_filename(fname: str) -> str:
        base = Path(fname).name
        for ext in (".fits", ".fit", ".fts", ".fz", ".gz"):
            if base.lower().endswith(ext):
                base = base[: -len(ext)]
        parts = [p for p in base.replace(".", "_").replace("-", "_").split("_") if p]
        for token in reversed(parts):
            cand = token.lower()
            if 1 <= len(cand) <= 3 and cand.isalpha():
                return cand.upper()
        return ""

    def _load_filter_map_from_index(self):
        idx_path = next(
            (
                p for p in (
                    step7_forced_phot_dir(self.params.P.result_dir) / "photometry_index.csv",
                    step8_psf_dir(self.params.P.result_dir) / "photometry_index.csv",
                    self.params.P.result_dir / "photometry_index.csv",
                )
                if p.exists()
            ),
            None,
        )
        if idx_path is None:
            return
        try:
            df = pd.read_csv(idx_path)
        except Exception:
            return
        if "file" not in df.columns or "filter" not in df.columns:
            return
        for fname, filt in zip(df["file"].astype(str), df["filter"].astype(str)):
            fkey = self._normalize_filter_key(filt)
            if not fkey:
                fkey = self._infer_filter_from_filename(fname)
            self._file_filter_map[fname] = fkey
            if fname.lower().startswith("crop_"):
                self._file_filter_map[fname[5:]] = fkey
            if fname.lower().startswith("cropped_"):
                self._file_filter_map[fname[8:]] = fkey

    def _get_filter_for_file(self, fname: str) -> str:
        fkey = self._file_filter_map.get(fname, "")
        if fkey:
            return fkey
        if fname == self.current_filename and self.header is not None:
            fkey = self._extract_filter_from_header(self.header)
        if not fkey:
            fkey = self._infer_filter_from_filename(fname)
        self._file_filter_map[fname] = fkey
        return fkey

    def _forced_phot_candidates(self, filename: str) -> list[Path]:
        forced_dir = step7_forced_phot_dir(self.params.P.result_dir)
        names = [str(filename)]
        low = str(filename).lower()
        if low.startswith("crop_"):
            names.append(str(filename)[5:])
        if low.startswith("cropped_"):
            names.append(str(filename)[8:])
        names.append(f"crop_{filename}")
        names.append(f"cropped_{filename}")

        out = []
        seen = set()
        for name in names:
            if not name or name in seen:
                continue
            seen.add(name)
            out.append(forced_dir / f"photometry_{name}.tsv")
        return out

    def _resolve_forced_phot_path(self, filename: str) -> Path | None:
        for path in self._forced_phot_candidates(filename):
            if path.exists():
                return path
        return None

    def load_master_ids(self):
        master_path = step9_selection_dir(self.params.P.result_dir) / "master_star_ids.csv"
        self.internal_id_map = {}
        self.source_id_from_internal = {}
        self._master_ids_persisted = False
        self._load_global_id_map()
        if master_path.exists():
            self._master_ids_persisted = True
            try:
                df = pd.read_csv(master_path)
                if "source_id" in df.columns:
                    sid_vals = parse_int64_series(df["source_id"]).dropna().astype("int64")
                    self.master_ids = set(sid_vals.tolist())
                    if "g_mag" in df.columns:
                        gmag_series = pd.to_numeric(df["g_mag"], errors="coerce")
                        sid_series = parse_int64_series(df["source_id"])
                        valid_sid = sid_series.notna()
                        sid_clean = sid_series.loc[valid_sid].astype("int64")
                        gmag_clean = gmag_series.loc[valid_sid]
                        self.master_gmag_map = dict(zip(sid_clean.tolist(), gmag_clean.tolist()))
                    # Prefer saved fixed ID mapping when available.
                    if "ID" in df.columns:
                        id_vals = pd.to_numeric(df["ID"], errors="coerce")
                    else:
                        id_vals = pd.Series([np.nan] * len(df))
                    sid_series = parse_int64_series(df["source_id"])
                    for sid_v, id_v in zip(sid_series.tolist(), id_vals):
                        if not (pd.notna(sid_v) and np.isfinite(id_v)):
                            continue
                        sid_i = int(sid_v)
                        id_i = int(id_v)
                        if sid_i in self.internal_id_map:
                            continue
                        if id_i in self.source_id_from_internal and self.source_id_from_internal[id_i] != sid_i:
                            continue
                        self.internal_id_map[sid_i] = id_i
                        self.source_id_from_internal[id_i] = sid_i
                    for sid_i in sorted(self.master_ids):
                        self._ensure_stable_id(sid_i)
                    self._load_refbuild_gmag_map()
                    self.log(f"Loaded {len(self.master_ids)} master IDs from {master_path.name}")
                    self.update_master_table()
            except Exception as e:
                self.log(f"Error loading master IDs: {e}")

    def load_psf_new_sources(self):
        """Step 8 PSF iter2에서 발견된 새 소스(det_uid < 0)를 마스터 목록에 추가.

        PSF 측광 iter2에서 잔차 이미지로부터 새로 검출된 별들은 음수 det_uid를
        가진다. 이 소스들은 RefBuild ref catalog에 없으므로 별도로 master_id를
        부여하여 마스터 목록에 추가한다.
        """
        psf_dir = step8_psf_dir(self.params.P.result_dir)
        if not psf_dir.exists():
            return

        rows = []
        for tsv in sorted(psf_dir.glob("photometry_*.tsv")):
            try:
                df = pd.read_csv(tsv, sep="\t")
                if "det_uid" not in df.columns or "iter_found" not in df.columns:
                    continue
                new = df[(pd.to_numeric(df["det_uid"], errors="coerce") < 0) &
                         (pd.to_numeric(df["iter_found"], errors="coerce") > 1)].copy()
                if not new.empty:
                    fname = tsv.name.replace("photometry_", "").replace(".tsv", "")
                    new["_source_file"] = fname
                    rows.append(new)
            except Exception:
                continue

        if not rows:
            return

        all_new = pd.concat(rows, ignore_index=True)
        # per det_uid 대표값 (첫 번째 등장 프레임 기준)
        unique_new = all_new.groupby("det_uid", as_index=False).first()
        n_new = len(unique_new)
        if n_new == 0:
            return

        # 기존 master_id 최대값 다음부터 부여 (Gaia source_id와 충돌 없이 음수 psf_uid 사용)
        # PSF 새 소스는 source_id = det_uid (음수) 로 내부 식별
        added = 0
        for _, row in unique_new.iterrows():
            det_uid = int(row["det_uid"])
            if det_uid not in self.internal_id_map:
                self._ensure_stable_id(det_uid)
                self.master_ids.add(det_uid)
                added += 1
            self.psf_iter2_ids.add(det_uid)

        if added > 0:
            self.log(f"PSF new sources: {added} iter2 소스 master_id 부여 완료 (det_uid<0)")
            self.update_master_table()

    def _load_refbuild_gmag_map(self):
        """Load G magnitudes from the APEX RefBuild catalog."""
        candidates = [
            step6_refbuild_dir(self.params.P.result_dir) / "ref_catalog.tsv",
            self.params.P.result_dir / "ref_catalog.tsv",
        ]
        for path in candidates:
            if not path.exists():
                continue
            try:
                df = pd.read_csv(path, sep="\t")
            except Exception:
                continue
            if "source_id" not in df.columns:
                continue
            g_col = None
            for cand in ("phot_g_mean_mag", "gaia_G", "g_mag"):
                if cand in df.columns:
                    g_col = cand
                    break
            if g_col is None:
                continue
            sid = parse_int64_series(df["source_id"])
            g = pd.to_numeric(df[g_col], errors="coerce")
            valid = sid.notna() & g.notna()
            if not valid.any():
                continue
            sid_i = sid.loc[valid].astype("int64")
            g_v = g.loc[valid].astype(float)
            loaded = 0
            for s, gv in zip(sid_i.tolist(), g_v.tolist()):
                if not np.isfinite(gv):
                    continue
                self.master_gmag_map[int(s)] = float(gv)
                loaded += 1
            if loaded > 0:
                self.log(f"Gmag fallback loaded from {path.name}: {loaded} rows")
                break

    def _load_global_id_map(self):
        """Load source_id -> ID map generated by RefBuild or this editor."""
        self._global_id_map = {}
        candidates = [
            step9_selection_dir(self.params.P.result_dir) / "sourceid_to_ID.csv",
            step6_refbuild_dir(self.params.P.result_dir) / "sourceid_to_ID.csv",
            self.params.P.result_dir / "sourceid_to_ID.csv",
        ]
        for path in candidates:
            if not path.exists():
                continue
            try:
                df = read_csv_int64_source_id(path)
            except Exception:
                continue
            if not {"source_id", "ID"} <= set(df.columns):
                continue
            sid_vals = parse_int64_series(df["source_id"])
            id_vals = pd.to_numeric(df["ID"], errors="coerce")
            for sid_v, id_v in zip(sid_vals, id_vals):
                if pd.notna(sid_v) and np.isfinite(id_v):
                    self._global_id_map[int(sid_v)] = int(id_v)
            if self._global_id_map:
                break

    def _next_available_id(self) -> int:
        used = {int(v) for v in self.source_id_from_internal.keys() if int(v) > 0}
        used.update(int(v) for v in self._global_id_map.values() if int(v) > 0)
        if not used:
            return 1
        return max(used) + 1

    def _ensure_stable_id(self, source_id: int) -> int:
        sid = int(source_id)
        if sid in self.internal_id_map:
            return int(self.internal_id_map[sid])

        # Prefer global ID mapping from upstream ID match when conflict-free.
        gid = self._global_id_map.get(sid)
        if gid is not None:
            gid = int(gid)
            owner = self.source_id_from_internal.get(gid)
            if owner is None or int(owner) == sid:
                self.internal_id_map[sid] = gid
                self.source_id_from_internal[gid] = sid
                return gid

        # Allocate new stable ID (never per-frame re-number).
        new_id = self._next_available_id()
        while new_id in self.source_id_from_internal:
            new_id += 1
        self.internal_id_map[sid] = int(new_id)
        self.source_id_from_internal[int(new_id)] = sid
        return int(new_id)

    def load_gaia_catalog(self):
        """Load Gaia catalog for source info lookup"""
        candidates = [
            step5_wcs_dir(self.params.P.result_dir) / "gaia_derived.csv",
            self.params.P.result_dir / "gaia_derived.csv",
            step5_wcs_dir(self.params.P.result_dir) / "gaia_fov.ecsv",
            self.params.P.result_dir / "gaia_fov.ecsv",
        ]
        for gaia_path in candidates:
            if not gaia_path.exists():
                continue
            try:
                if gaia_path.suffix.lower() == ".ecsv":
                    self.gaia_df = read_ecsv_int64_source_id(gaia_path)
                else:
                    self.gaia_df = pd.read_csv(gaia_path)
                if "source_id" in self.gaia_df.columns:
                    sid_series = parse_int64_series(self.gaia_df["source_id"])
                    valid = sid_series.notna()
                    dropped = int((~valid).sum())
                    if dropped > 0:
                        self.gaia_df = self.gaia_df.loc[valid].copy()
                        sid_series = sid_series.loc[valid]
                        self.log(f"Gaia catalog: dropped {dropped} rows with invalid source_id")
                    self.gaia_df["source_id"] = sid_series.astype("int64")
                self._rebuild_gaia_gmag_map()
                self.log(f"Gaia catalog loaded: {len(self.gaia_df)} sources ({gaia_path.name})")
                if self.master_ids:
                    self.update_master_table()
                    self.update_overlay()
                return
            except Exception as e:
                self.log(f"Failed to load Gaia catalog ({gaia_path.name}): {e}")
                self.gaia_df = None
                self._gaia_gmag_map = {}

    def _rebuild_gaia_gmag_map(self):
        self._gaia_gmag_map = {}
        if self.gaia_df is None or len(self.gaia_df) == 0:
            return
        if "source_id" not in self.gaia_df.columns:
            return
        g_col = None
        for cand in ("phot_g_mean_mag", "gaia_G", "g_mag"):
            if cand in self.gaia_df.columns:
                g_col = cand
                break
        if g_col is None:
            return
        sid = parse_int64_series(self.gaia_df["source_id"])
        g = pd.to_numeric(self.gaia_df[g_col], errors="coerce")
        valid = sid.notna() & g.notna()
        if not valid.any():
            return
        sid_i = sid.loc[valid].astype("int64").to_numpy(np.int64, copy=False)
        g_v = g.loc[valid].astype(float).to_numpy(float, copy=False)
        for s, gv in zip(sid_i.tolist(), g_v.tolist()):
            if np.isfinite(gv):
                self._gaia_gmag_map[int(s)] = float(gv)

    @staticmethod
    def _pick_membership_col(cols) -> str | None:
        for c in ("gaia_pmem", "pmem_gaia", "membership_prob_gaia", "membership_prob", "pmem"):
            if c in cols:
                return c
        return None

    @staticmethod
    def _logpdf_gauss(x: np.ndarray, mu: np.ndarray, cov: np.ndarray) -> np.ndarray:
        d = x.shape[1]
        cov_r = np.asarray(cov, float) + np.eye(d) * 1e-6
        sign, logdet = np.linalg.slogdet(cov_r)
        if sign <= 0:
            return np.full(x.shape[0], -np.inf, dtype=float)
        try:
            inv = np.linalg.inv(cov_r)
        except Exception:
            return np.full(x.shape[0], -np.inf, dtype=float)
        diff = x - mu[None, :]
        q = np.einsum("ni,ij,nj->n", diff, inv, diff)
        return -0.5 * (d * np.log(2.0 * np.pi) + logdet + q)

    def _fit_two_component_gmm(self, x_fit: np.ndarray):
        n, d = x_fit.shape
        if n < max(30, d * 8):
            return None

        center = np.nanmedian(x_fit, axis=0)
        mad = np.nanmedian(np.abs(x_fit - center), axis=0)
        mad = np.where(np.isfinite(mad) & (mad > 1e-6), mad, 1.0)
        z = (x_fit - center[None, :]) / mad[None, :]
        d2 = np.sum(z * z, axis=1)
        q40 = float(np.nanquantile(d2, 0.40))
        m0 = d2 <= q40
        if m0.sum() < max(12, d * 3) or m0.sum() > (n - max(12, d * 3)):
            order = np.argsort(d2)
            m0 = np.zeros(n, dtype=bool)
            m0[order[: max(n // 2, 1)]] = True
        m1 = ~m0

        def _cov(arr: np.ndarray) -> np.ndarray:
            if arr.shape[0] < (d + 1):
                return np.eye(d, dtype=float)
            c = np.cov(arr, rowvar=False)
            if np.ndim(c) == 0:
                c = np.eye(d, dtype=float) * float(c)
            return np.asarray(c, float) + np.eye(d, dtype=float) * 1e-4

        pi = np.array([max(m0.mean(), 1e-3), max(m1.mean(), 1e-3)], dtype=float)
        pi /= pi.sum()
        mu = np.vstack([
            np.nanmean(x_fit[m0], axis=0),
            np.nanmean(x_fit[m1], axis=0),
        ])
        cov = np.stack([_cov(x_fit[m0]), _cov(x_fit[m1])], axis=0)
        last_ll = -np.inf

        for _ in range(80):
            lp0 = np.log(max(pi[0], 1e-9)) + self._logpdf_gauss(x_fit, mu[0], cov[0])
            lp1 = np.log(max(pi[1], 1e-9)) + self._logpdf_gauss(x_fit, mu[1], cov[1])
            m = np.maximum(lp0, lp1)
            e0 = np.exp(lp0 - m)
            e1 = np.exp(lp1 - m)
            den = e0 + e1 + 1e-12
            r0 = e0 / den
            r1 = e1 / den
            nk = np.array([r0.sum(), r1.sum()], dtype=float)
            if np.any(nk < (d + 2)):
                break
            pi = nk / float(n)
            mu[0] = (r0[:, None] * x_fit).sum(axis=0) / nk[0]
            mu[1] = (r1[:, None] * x_fit).sum(axis=0) / nk[1]
            for k, rk in enumerate((r0, r1)):
                diff = x_fit - mu[k][None, :]
                cov[k] = (diff.T * rk).dot(diff) / max(nk[k], 1.0)
                cov[k] += np.eye(d, dtype=float) * 1e-4
            ll = float(np.sum(m + np.log(den)))
            if np.isfinite(last_ll):
                if abs(ll - last_ll) < 1e-4 * max(1.0, abs(last_ll)):
                    break
            last_ll = ll

        det0 = abs(float(np.linalg.det(cov[0])))
        det1 = abs(float(np.linalg.det(cov[1])))
        cluster_idx = 0 if det0 <= det1 else 1
        return {
            "pi": pi,
            "mu": mu,
            "cov": cov,
            "cluster_idx": int(cluster_idx),
        }

    def _compute_membership_from_master(self) -> bool:
        # Prefer gaia_fov.ecsv (written by step5 WCS) — always has pmra/pmdec/parallax
        # keyed by Gaia source_id.  master_catalog.tsv only carries photometry (G/BP/RP),
        # so it never has the astrometric columns needed for GMM membership.
        result_dir = self.params.P.result_dir
        req = ("pmra", "pmdec", "parallax")

        df = None
        df_is_gaia_keyed = False  # True when source_id column = Gaia source_id

        gaia_fov_path = step5_wcs_dir(result_dir) / "gaia_fov.ecsv"
        if gaia_fov_path.exists():
            try:
                from astropy.table import Table as _ATable
                tab = _ATable.read(str(gaia_fov_path), format="ascii.ecsv")
                cols_lower = [c.lower() for c in tab.colnames]
                if cols_lower != list(tab.colnames):
                    tab.rename_columns(tab.colnames, cols_lower)
                gaia_df = tab.to_pandas()
                if all(c in gaia_df.columns for c in req):
                    # source_id here is the Gaia source_id
                    df = gaia_df.rename(columns={"source_id": "gaia_source_id"})
                    df_is_gaia_keyed = True
            except Exception:
                pass

        # Fallback: master_catalog.tsv (only works if astrometric cols were written)
        if df is None:
            master_candidates = [
                step6_refbuild_dir(result_dir) / "master_catalog.tsv",
                result_dir / "master_catalog.tsv",
            ]
            master_path = next((p for p in master_candidates if p.exists()), None)
            if master_path is not None:
                try:
                    master_df = pd.read_csv(master_path, sep="\t")
                    if not master_df.empty and all(c in master_df.columns for c in req):
                        df = master_df
                except Exception:
                    pass

        if df is None or df.empty:
            return False

        pmra = pd.to_numeric(df["pmra"], errors="coerce").to_numpy(float)
        pmdec = pd.to_numeric(df["pmdec"], errors="coerce").to_numpy(float)
        plx = pd.to_numeric(df["parallax"], errors="coerce").to_numpy(float)
        finite = np.isfinite(pmra) & np.isfinite(pmdec) & np.isfinite(plx)
        if int(finite.sum()) < 30:
            return False

        fit_mask = finite.copy()
        if "ruwe" in df.columns:
            ruwe = pd.to_numeric(df["ruwe"], errors="coerce").to_numpy(float)
            fit_mask &= (~np.isfinite(ruwe)) | (ruwe <= 2.0)
        if "visibility_periods_used" in df.columns:
            vpu = pd.to_numeric(df["visibility_periods_used"], errors="coerce").to_numpy(float)
            fit_mask &= (~np.isfinite(vpu)) | (vpu >= 8.0)
        if int(fit_mask.sum()) < 25:
            fit_mask = finite.copy()

        x_fit = np.column_stack([pmra[fit_mask], pmdec[fit_mask], plx[fit_mask]])
        model = self._fit_two_component_gmm(x_fit)
        if model is None:
            return False

        x_all = np.column_stack([pmra[finite], pmdec[finite], plx[finite]])
        pi = np.asarray(model["pi"], float)
        mu = np.asarray(model["mu"], float)
        cov = np.asarray(model["cov"], float)
        k_cluster = int(model["cluster_idx"])

        lp0 = np.log(max(pi[0], 1e-9)) + self._logpdf_gauss(x_all, mu[0], cov[0])
        lp1 = np.log(max(pi[1], 1e-9)) + self._logpdf_gauss(x_all, mu[1], cov[1])
        m = np.maximum(lp0, lp1)
        e0 = np.exp(lp0 - m)
        e1 = np.exp(lp1 - m)
        den = e0 + e1 + 1e-12
        p_cluster = (e0 if k_cluster == 0 else e1) / den

        pmem_arr = np.full(len(df), np.nan, dtype=float)
        pmem_arr[finite] = np.clip(p_cluster, 0.0, 1.0)

        gaia_ids = parse_int64_series(df["gaia_source_id"]).tolist() if "gaia_source_id" in df.columns else [pd.NA] * len(df)
        src_ids  = parse_int64_series(df["source_id"]).tolist()     if "source_id"      in df.columns else [pd.NA] * len(df)
        int_ids  = pd.to_numeric(df["ID"], errors="coerce").tolist() if "ID"             in df.columns else [np.nan] * len(df)

        # When we loaded from gaia_fov (gaia_keyed), source_id col was renamed to
        # gaia_source_id, so src_ids is empty — use gaia_ids for membership_by_source
        # as well (hybrid mode assigns Gaia source_id as source_id for matched stars).
        if df_is_gaia_keyed:
            src_ids = gaia_ids

        # Load master catalog to resolve gaia_source_id → source_id / ID mapping
        if df_is_gaia_keyed:
            master_candidates = [
                step6_refbuild_dir(result_dir) / "master_catalog.tsv",
                result_dir / "master_catalog.tsv",
            ]
            master_path = next((p for p in master_candidates if p.exists()), None)
            if master_path is not None:
                try:
                    mdf = pd.read_csv(master_path, sep="\t")
                    if "gaia_source_id" in mdf.columns and "source_id" in mdf.columns:
                        g2s = dict(zip(
                            parse_int64_series(mdf["gaia_source_id"]).tolist(),
                            parse_int64_series(mdf["source_id"]).tolist(),
                        ))
                        src_ids = [g2s.get(g, g) for g in gaia_ids]
                    if "gaia_source_id" in mdf.columns and "ID" in mdf.columns:
                        g2id = dict(zip(
                            parse_int64_series(mdf["gaia_source_id"]).tolist(),
                            pd.to_numeric(mdf["ID"], errors="coerce").tolist(),
                        ))
                        int_ids = [g2id.get(g, np.nan) for g in gaia_ids]
                except Exception:
                    pass

        n_src = n_gaia = n_id = 0
        for s, g, i, p in zip(src_ids, gaia_ids, int_ids, pmem_arr.tolist()):
            if not np.isfinite(p):
                continue
            pp = float(np.clip(p, 0.0, 1.0))
            if pd.notna(s):
                self.membership_by_source[int(s)] = pp
                n_src += 1
            if pd.notna(g):
                self.membership_by_gaia[int(g)] = pp
                n_gaia += 1
            if np.isfinite(i):
                self.membership_by_id[int(i)] = pp
                n_id += 1

        src_label = "gaia_fov.ecsv" if df_is_gaia_keyed else "master_catalog.tsv"
        self.log(
            f"Membership map computed from {src_label} astrometry "
            f"(source={n_src}, gaia={n_gaia}, id={n_id}, fit={int(fit_mask.sum())})"
        )
        return (n_src + n_gaia + n_id) > 0

    def _load_membership_map(self):
        self.membership_by_source = {}
        self.membership_by_gaia = {}
        self.membership_by_id = {}
        self._membership_loaded = True

        result_dir = self.params.P.result_dir
        candidates = [
            step5_wcs_dir(result_dir) / "gaia_derived.csv",
            result_dir / "gaia_derived.csv",
            result_dir / "cmd_with_gaia_membership.csv",
            step9_selection_dir(result_dir) / "cmd_with_gaia_membership.csv",
            result_dir / "cmd_with_membership.csv",
            step9_selection_dir(result_dir) / "cmd_with_membership.csv",
            result_dir / "median_by_ID_filter_wide_cmd.csv",
            step9_selection_dir(result_dir) / "median_by_ID_filter_wide_cmd.csv",
            step6_refbuild_dir(result_dir) / "master_catalog.tsv",
            result_dir / "master_catalog.tsv",
        ]

        for path in candidates:
            if not path.exists():
                continue
            try:
                if path.suffix.lower() == ".tsv":
                    df = pd.read_csv(path, sep="\t")
                else:
                    df = pd.read_csv(path)
            except Exception:
                continue
            if df.empty:
                continue

            p_col = self._pick_membership_col(df.columns)
            if p_col is None:
                continue
            prob = pd.to_numeric(df[p_col], errors="coerce").to_numpy(float)

            n_src = 0
            n_gaia = 0
            n_id = 0
            if "source_id" in df.columns:
                src = parse_int64_series(df["source_id"])
                for sid, p in zip(src.tolist(), prob.tolist()):
                    if pd.isna(sid) or (not np.isfinite(p)):
                        continue
                    self.membership_by_source[int(sid)] = float(np.clip(p, 0.0, 1.0))
                    n_src += 1

            if "gaia_source_id" in df.columns:
                gs = parse_int64_series(df["gaia_source_id"])
                for sid, p in zip(gs.tolist(), prob.tolist()):
                    if pd.isna(sid) or (not np.isfinite(p)):
                        continue
                    self.membership_by_gaia[int(sid)] = float(np.clip(p, 0.0, 1.0))
                    n_gaia += 1

            if "ID" in df.columns:
                ids = pd.to_numeric(df["ID"], errors="coerce")
                for id_v, p in zip(ids.tolist(), prob.tolist()):
                    if (not np.isfinite(id_v)) or (not np.isfinite(p)):
                        continue
                    self.membership_by_id[int(id_v)] = float(np.clip(p, 0.0, 1.0))
                    n_id += 1

            if (n_src + n_gaia + n_id) > 0:
                self.log(
                    f"Membership map loaded: source={n_src}, gaia={n_gaia}, id={n_id} "
                    f"from {path.name} ({p_col})"
                )
                return

        # Compute membership from RefBuild astrometric columns when no saved map exists.
        if self._compute_membership_from_master():
            return

        if not self._membership_log_once:
            self.log("Membership map not found (gaia_pmem column/file missing).")
            self._membership_log_once = True

    def _ensure_membership_map(self, force: bool = False):
        if (not force) and (not bool(getattr(self.params.P, "step8_membership_overlay_enable", False))):
            return
        if not self._membership_loaded:
            self._load_membership_map()

    def _membership_prob_for_sid(self, sid: int) -> float:
        s = int(sid)
        p = self.membership_by_source.get(s, None)
        if p is not None:
            return float(p)
        if s > 0:
            p = self.membership_by_gaia.get(s, None)
            if p is not None:
                return float(p)
        id_v = self.internal_id_map.get(s, self._global_id_map.get(s, None))
        if id_v is not None:
            p = self.membership_by_id.get(int(id_v), None)
            if p is not None:
                return float(p)
        return np.nan

    def get_gaia_info(self, source_id: int) -> str:
        """Get Gaia info string for a source_id"""
        if self.gaia_df is None or len(self.gaia_df) == 0:
            return ""
        try:
            row = self.gaia_df[self.gaia_df["source_id"] == source_id]
            if len(row) == 0:
                return ""
            row = row.iloc[0]
            parts = []
            # Magnitudes
            g = row.get("phot_g_mean_mag", np.nan)
            bp = row.get("phot_bp_mean_mag", np.nan)
            rp = row.get("phot_rp_mean_mag", np.nan)
            if np.isfinite(g):
                parts.append(f"G={g:.3f}")
            if np.isfinite(bp):
                parts.append(f"BP={bp:.3f}")
            if np.isfinite(rp):
                parts.append(f"RP={rp:.3f}")
            # Color
            if np.isfinite(bp) and np.isfinite(rp):
                parts.append(f"BP-RP={bp-rp:.3f}")
            # Coordinates
            ra = row.get("ra", np.nan)
            dec = row.get("dec", np.nan)
            if np.isfinite(ra) and np.isfinite(dec):
                parts.append(f"RA={ra:.5f}")
                parts.append(f"Dec={dec:.5f}")
            # Proper motion
            pmra = row.get("pmra", np.nan)
            pmdec = row.get("pmdec", np.nan)
            if np.isfinite(pmra) and np.isfinite(pmdec):
                parts.append(f"PM=({pmra:.2f},{pmdec:.2f})mas/yr")
            return " | ".join(parts)
        except Exception:
            return ""

    def update_master_table(self):
        """Update master table using stable ID (not frame-dependent numbering)."""
        self._ensure_membership_map(force=True)
        rows = []
        for sid in sorted(self.master_ids):
            sid_i = int(sid)
            fixed_id = self._ensure_stable_id(sid_i)
            g_mag = self._get_gmag_for_source(sid_i)
            pmem = self._membership_prob_for_sid(sid_i)
            rows.append((int(fixed_id), sid_i, g_mag, pmem))

        rows.sort(key=lambda x: int(x[0]))

        self.master_table.blockSignals(True)
        self.master_table.setUpdatesEnabled(False)
        try:
            self.master_table.setRowCount(len(rows))
            for i, (fixed_id, sid, g_mag, pmem) in enumerate(rows):
                self.master_table.setItem(i, 0, QTableWidgetItem(str(fixed_id)))
                self.master_table.setItem(i, 1, QTableWidgetItem(str(sid)))
                if np.isfinite(g_mag):
                    g_str = f"{g_mag:.3f}"
                elif int(sid) < 0:
                    g_str = "local"
                else:
                    g_str = "-"
                self.master_table.setItem(i, 2, QTableWidgetItem(g_str))
                p_str = f"{pmem:.3f}" if np.isfinite(pmem) else "-"
                self.master_table.setItem(i, 3, QTableWidgetItem(p_str))
        finally:
            self.master_table.setUpdatesEnabled(True)
            self.master_table.blockSignals(False)

        n_total = len(rows)
        n_gmag = sum(1 for _, _, g, _ in rows if np.isfinite(g))
        n_pmem = sum(1 for _, _, _, p in rows if np.isfinite(p))
        if hasattr(self, "log_text"):
            self.log(f"Master IDs: {n_total} | Gmag: {n_gmag} | Pmem: {n_pmem} | Stable IDs")

        self.update_navigation_buttons()

    def on_table_selection_changed(self):
        """Handle table selection change"""
        rows = self.master_table.selectionModel().selectedRows()
        if rows:
            row_idx = rows[0].row()
            sid_item = self.master_table.item(row_idx, 1)
            if sid_item:
                try:
                    self.selected_source_id = int(sid_item.text())
                    self.selected_phot_row = None
                    self.selected_phot_file = None
                    in_master = "✓ IN MASTER"
                    internal_id = self.internal_id_map.get(self.selected_source_id, "?")
                    self.selected_label.setText(
                        f"Selected: ID {internal_id} | source_id: {self.selected_source_id} ({in_master})"
                    )
                    self.update_overlay()
                except ValueError:
                    pass

    def select_source_in_table(self, source_id: int):
        """Select a source_id in the master table"""
        for row in range(self.master_table.rowCount()):
            sid_item = self.master_table.item(row, 1)
            if sid_item and int(sid_item.text()) == source_id:
                self.master_table.blockSignals(True)
                self.master_table.selectRow(row)
                self.master_table.scrollToItem(sid_item)
                self.master_table.blockSignals(False)
                break

    def on_file_changed(self, index):
        if index < 0 or index >= len(self.file_list):
            return
        # Fast switch: avoid full redraw if shape is unchanged.
        self.load_and_display(quick_switch=True)

    def _get_fits_path(self, filename):
        if self.use_cropped:
            return step2_cropped_dir(self.params.P.result_dir) / filename
        return self.params.P.data_dir / filename

    def _load_fits_cached(self, filename):
        """Return (image_data, header) from cache or disk. LRU eviction."""
        with self._prefetch_lock:
            if filename in self._fits_cache:
                # Move to end (most recently used)
                try:
                    self._fits_cache_order.remove(filename)
                except ValueError:
                    pass
                self._fits_cache_order.append(filename)
                return self._fits_cache[filename]
        file_path = self._get_fits_path(filename)
        with fits.open(file_path, memmap=False) as hdul:
            # Use float32 to reduce memory bandwidth on frame switching.
            data_raw = hdul[0].data
            if data_raw is None:
                raise ValueError(f"FITS image data is empty: {file_path.name}")
            data = np.asarray(data_raw, dtype=np.float32)
            header = hdul[0].header.copy()
        with self._prefetch_lock:
            # In case prefetch completed first, return existing cache entry.
            if filename in self._fits_cache:
                try:
                    self._fits_cache_order.remove(filename)
                except ValueError:
                    pass
                self._fits_cache_order.append(filename)
                return self._fits_cache[filename]
            # Evict oldest if over limit
            while len(self._fits_cache_order) >= self._FITS_CACHE_SIZE:
                oldest = self._fits_cache_order.pop(0)
                self._fits_cache.pop(oldest, None)
            self._fits_cache[filename] = (data, header)
            self._fits_cache_order.append(filename)
            return data, header

    def _prefetch_fits_worker(self, filename: str):
        try:
            file_path = self._get_fits_path(filename)
            if not file_path.exists():
                return
            with fits.open(file_path, memmap=False) as hdul:
                data_raw = hdul[0].data
                if data_raw is None:
                    return
                data = np.asarray(data_raw, dtype=np.float32)
                header = hdul[0].header.copy()
            with self._prefetch_lock:
                if filename not in self._fits_cache:
                    while len(self._fits_cache_order) >= self._FITS_CACHE_SIZE:
                        oldest = self._fits_cache_order.pop(0)
                        self._fits_cache.pop(oldest, None)
                    self._fits_cache[filename] = (data, header)
                    self._fits_cache_order.append(filename)
        except Exception:
            pass
        finally:
            with self._prefetch_lock:
                self._prefetch_pending.discard(filename)

    def _schedule_prefetch_neighbors(self):
        if not self.file_list:
            return
        idx = self.file_combo.currentIndex()
        if idx < 0:
            return
        candidates = []
        if idx + 1 < len(self.file_list):
            candidates.append(self.file_list[idx + 1])
        if idx - 1 >= 0:
            candidates.append(self.file_list[idx - 1])
        for fname in candidates:
            with self._prefetch_lock:
                if fname in self._fits_cache or fname in self._prefetch_pending:
                    continue
                self._prefetch_pending.add(fname)
            try:
                self._prefetch_executor.submit(self._prefetch_fits_worker, fname)
            except Exception:
                with self._prefetch_lock:
                    self._prefetch_pending.discard(fname)

    def load_and_display(self, quick_switch=False):
        filename = self.file_combo.currentText()
        if not filename:
            return
        try:
            new_data, new_header = self._load_fits_cached(filename)
            shape_changed = (self.image_data is None or new_data.shape != self.image_data.shape)
            self.image_data = new_data
            self.header = new_header
            self.current_filename = filename
            self.load_photometry_for_file(filename)
            self._update_viewer_image(keep_view=quick_switch and not shape_changed)
            self.update_overlay()
            if self._auto_master_dirty:
                self.save_master_ids(log_action="auto_add")
                self._auto_master_dirty = False
            self._schedule_prefetch_neighbors()
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load: {str(e)}")

    def closeEvent(self, event):
        try:
            self._prefetch_executor.shutdown(wait=False, cancel_futures=True)
        except Exception:
            pass
        super().closeEvent(event)

    @staticmethod
    def _truthy_series(series: pd.Series, default: bool = False) -> pd.Series:
        if series is None:
            return pd.Series([], dtype=bool)
        if pd.api.types.is_bool_dtype(series):
            return series.fillna(default).astype(bool)
        num = pd.to_numeric(series, errors="coerce")
        if pd.api.types.is_numeric_dtype(series) or num.notna().any():
            return num.fillna(1 if default else 0).astype(float) != 0.0
        text = series.astype(str).str.strip().str.lower()
        true_values = {"1", "true", "t", "yes", "y", "on", "forced", "forced_only"}
        false_values = {"0", "false", "f", "no", "n", "off", "detected", "matched", ""}
        out = pd.Series(default, index=series.index, dtype=bool)
        out.loc[text.isin(true_values)] = True
        out.loc[text.isin(false_values)] = False
        return out

    @classmethod
    def _forced_mask_from_photometry(cls, df: pd.DataFrame) -> np.ndarray:
        forced = np.zeros(len(df), dtype=bool)
        if "forced_flag" in df.columns:
            forced |= cls._truthy_series(df["forced_flag"], default=False).to_numpy(bool, copy=False)
        if "detected_flag" in df.columns:
            detected = cls._truthy_series(df["detected_flag"], default=True).to_numpy(bool, copy=False)
            forced |= ~detected
        if "det_uid" in df.columns:
            det_uid = pd.to_numeric(df["det_uid"], errors="coerce").to_numpy(float)
            forced |= np.isfinite(det_uid) & (det_uid < 0)
        for col in ("centering_quality", "source_type", "source_role", "match_status", "phot_status"):
            if col not in df.columns:
                continue
            text = df[col].astype(str).str.strip().str.lower()
            forced |= text.str.contains("forced", na=False).to_numpy(bool, copy=False)
        return forced

    def load_photometry_for_file(self, filename):
        # Check cache first
        if filename in self._phot_cache:
            self.phot_df = self._phot_cache[filename]
            arr = self._phot_arr_cache.get(filename, None)
            if arr is not None:
                self._auto_add_detections_to_master(self.phot_df, arr[2])
            else:
                self._auto_add_detections_to_master(self.phot_df)
            return
        phot_path = self._resolve_forced_phot_path(filename)
        if phot_path is not None and phot_path.exists():
            try:
                sep = "\t" if phot_path.suffix.lower() == ".tsv" else ","
                df = read_csv_int64_source_id(phot_path, sep=sep)
                if {"x_fit", "y_fit"} <= set(df.columns) and not {"x", "y"} <= set(df.columns):
                    df = df.rename(columns={"x_fit": "x", "y_fit": "y"})
                if {"x", "y", "source_id"} <= set(df.columns):
                    extra_cols = [
                        c for c in (
                            "forced_flag", "detected_flag", "det_uid",
                            "centering_quality", "source_type", "source_role",
                            "match_status", "phot_status",
                        )
                        if c in df.columns
                    ]
                    clean = df[["x", "y", "source_id"] + extra_cols].copy()
                    clean["x"] = pd.to_numeric(clean["x"], errors="coerce")
                    clean["y"] = pd.to_numeric(clean["y"], errors="coerce")
                    valid = (
                        np.isfinite(clean["x"].to_numpy(float))
                        & np.isfinite(clean["y"].to_numpy(float))
                    )
                    dropped = int((~valid).sum())
                    if dropped > 0:
                        self.log(f"[{filename}] photometry: dropped {dropped} invalid rows (x/y)")
                    clean = clean.loc[valid].copy()
                    unmatched = int(clean["source_id"].apply(
                        lambda s: pd.isna(s)
                    ).sum())
                    if unmatched > 0:
                        self.log(f"[{filename}] photometry: unmatched rows kept={unmatched}")
                    clean["source_id"] = parse_int64_series(clean["source_id"]).fillna(0).astype("int64")
                    result = clean.reset_index(drop=True)
                    x_arr = result["x"].to_numpy(float, copy=False)
                    y_arr = result["y"].to_numpy(float, copy=False)
                    sid_arr = result["source_id"].to_numpy(np.int64, copy=False)
                    # forced=True means forced-photometry-only, regardless of master inclusion.
                    forced_arr = self._forced_mask_from_photometry(result)
                    self._phot_cache[filename] = result
                    self._phot_arr_cache[filename] = (x_arr, y_arr, sid_arr, forced_arr)
                    self.phot_df = result
                    self._auto_add_detections_to_master(self.phot_df, sid_arr)
                    return
            except Exception as exc:
                self.log(f"[{filename}] failed to load forced photometry {phot_path.name}: {exc}")
        else:
            tried = ", ".join(p.name for p in self._forced_phot_candidates(filename)[:3])
            self.log(f"[{filename}] forced photometry TSV not found (tried: {tried})")
        empty = pd.DataFrame(columns=["x", "y", "source_id"])
        self._phot_cache[filename] = empty
        self._phot_arr_cache[filename] = (
            np.array([], dtype=float),
            np.array([], dtype=float),
            np.array([], dtype=np.int64),
            np.array([], dtype=bool),
        )
        self.phot_df = empty

    def _auto_add_detections_to_master(self, df: pd.DataFrame, sid_arr: np.ndarray | None = None):
        if self._master_ids_persisted:
            return
        if sid_arr is None:
            try:
                sid_arr = parse_int64_series(df["source_id"]).dropna().astype("int64").to_numpy(np.int64, copy=False)
            except Exception:
                return
        if sid_arr.size == 0:
            return
        sids = set(np.unique(sid_arr[sid_arr != 0]).tolist())
        new_ids = sids - self.master_ids
        if not new_ids:
            return
        self.master_ids |= new_ids
        self._auto_master_dirty = True

    def _update_viewer_image(self, keep_view=False):
        if self._viewer is None or self.image_data is None:
            return
        self._viewer.set_data(self.image_data)
        self._viewer.auto_stf()
        if not keep_view:
            self._viewer.fit_in_view()

    def update_overlay(self):
        if self._viewer is None:
            return

        markers: list[OverlayMarker] = []
        c = self._overlay_colors
        vis = self._overlay_visible

        arr = self._phot_arr_cache.get(self.current_filename or "", None)
        if arr is None and self.phot_df is not None and not self.phot_df.empty:
            x_f = pd.to_numeric(self.phot_df["x"], errors="coerce").to_numpy(float)
            y_f = pd.to_numeric(self.phot_df["y"], errors="coerce").to_numpy(float)
            sid_f = parse_int64_series(self.phot_df["source_id"])
            valid = np.isfinite(x_f) & np.isfinite(y_f) & sid_f.notna().to_numpy(bool)
            if np.any(valid):
                x_f = x_f[valid]
                y_f = y_f[valid]
                sids_f = sid_f.loc[valid].astype("int64").to_numpy(np.int64, copy=False)
                forced_f = self._forced_mask_from_photometry(self.phot_df).astype(bool, copy=False)[valid]
                self._phot_arr_cache[self.current_filename or ""] = (x_f, y_f, sids_f, forced_f)
                arr = (x_f, y_f, sids_f, forced_f)

        if arr is not None:
            if len(arr) == 4:
                x, y, sids, forced = arr
            else:
                x, y, sids = arr
                forced = np.zeros(len(sids), dtype=bool)
            if sids.size > 0:
                if self.master_ids:
                    master_vals = np.fromiter((int(v) for v in self.master_ids), dtype=np.int64, count=len(self.master_ids))
                    in_master = np.isin(sids, master_vals)
                else:
                    in_master = np.zeros_like(sids, dtype=bool)
                is_matched   = sids != 0
                is_gaia_sid  = sids > 0
                is_local_sid = sids < 0
                is_forced    = is_matched & forced                     # step4 미검출 forced
                is_removed   = is_matched & ~in_master & ~forced
                is_gaia      = is_matched & in_master & is_gaia_sid & ~forced
                is_local     = is_matched & in_master & is_local_sid & ~forced

                is_member = np.zeros(len(sids), dtype=bool)
                if self._show_members:
                    self._ensure_membership_map()
                    thr = float(getattr(self.params.P, "step8_membership_threshold", 0.5))
                    idx_g = np.where(is_gaia)[0]
                    if idx_g.size > 0:
                        pmem = np.fromiter(
                            (self._membership_prob_for_sid(int(s)) for s in sids[idx_g]),
                            dtype=float, count=idx_g.size,
                        )
                        is_member[idx_g] = np.isfinite(pmem) & (pmem >= thr)

                sel_id  = self.selected_source_id
                sel_row = self.selected_phot_row
                sel_file = self.selected_phot_file
                sel_vis = vis.get("selected", True)
                sel_col = QColor(c["selected"])

                def _add_markers(mask, key, radius):
                    if not vis.get(key, True):
                        return
                    cat_col = QColor(c[key])
                    for idx in np.flatnonzero(mask):
                        sid = int(sids[idx])
                        row_selected = (
                            sel_vis
                            and sel_row is not None
                            and sel_file == self.current_filename
                            and int(idx) == int(sel_row)
                        )
                        sid_selected = (
                            sel_vis
                            and sel_id is not None
                            and sid == int(sel_id)
                            and not row_selected
                        )
                        is_sel = row_selected or sid_selected
                        col = sel_col if is_sel else cat_col
                        markers.append(OverlayMarker(
                            col=float(x[idx]), row=float(y[idx]),
                            radius=radius, color=col,
                            member=bool(is_member[idx]),
                        ))

                _add_markers(is_removed, "removed", 5.0)
                _add_markers(is_forced,  "forced",  5.5)
                _add_markers(is_local,     "local",   5.5)
                _add_markers(is_gaia,      "gaia",    6.5)

        # ROI circle (saved, projected to current frame pixel coords)
        pix = self._sky_roi_to_pixels()
        if pix is not None:
            cx, cy, r_px = pix
            markers.append(OverlayMarker(col=cx, row=cy, radius=r_px, color=QColor(0, 229, 255, 200)))

        # ROI drag preview
        if self._roi_preview is not None:
            cx, cy, r_px = self._roi_preview
            markers.append(OverlayMarker(col=cx, row=cy, radius=max(r_px, 1.0), color=QColor(0, 229, 255, 130)))

        self._viewer.set_overlay_markers(markers)

    def on_click_xy(self, x: float, y: float):
        self.setFocus()
        self.last_click_xy = (x, y)

        # Try to find nearby detected source
        search_r = float(getattr(self.params.P, "search_radius_px", 7.0))
        found_detected = False

        if self.phot_df is not None and not self.phot_df.empty:
            dx = self.phot_df["x"].to_numpy(float) - x
            dy = self.phot_df["y"].to_numpy(float) - y
            dist2 = dx * dx + dy * dy
            if dist2.size > 0:
                i = int(np.argmin(dist2))
                if dist2[i] <= search_r * search_r:
                    found_detected = True
                    sid = int(self.phot_df.iloc[i]["source_id"])
                    if sid == 0:
                        self.selected_source_id = None
                        self.selected_phot_row = int(i)
                        self.selected_phot_file = self.current_filename
                        frame = self.current_filename or "?"
                        px_x = float(self.phot_df.iloc[i]["x"])
                        px_y = float(self.phot_df.iloc[i]["y"])
                        self.selected_label.setText(
                            f"Selected: unmatched detection at ({px_x:.1f}, {px_y:.1f})"
                        )
                        self.log(f"[{frame}] Unmatched detection selected | px=({px_x:.1f}, {px_y:.1f})")
                        self.update_overlay()
                        return
                    self.selected_source_id = sid
                    self.selected_phot_row = int(i)
                    self.selected_phot_file = self.current_filename
                    in_master = "✓ IN MASTER" if self.selected_source_id in self.master_ids else "✗ not in master"
                    internal_id = self.internal_id_map.get(
                        self.selected_source_id,
                        self._global_id_map.get(self.selected_source_id, "-")
                    )
                    self.selected_label.setText(
                        f"Selected: ID {internal_id} | source_id: {self.selected_source_id} ({in_master})"
                    )

                    # Log with frame name
                    gaia_info = self.get_gaia_info(self.selected_source_id)
                    row = self.phot_df.iloc[i]
                    px_x, px_y = row["x"], row["y"]
                    frame = self.current_filename or "?"
                    self.log(f"[{frame}] Selected ID {internal_id}: {self.selected_source_id} | px=({px_x:.1f}, {px_y:.1f}) | {in_master}")
                    if gaia_info:
                        self.log(f"  Gaia: {gaia_info}")

                    # Auto-select in table if in master
                    if self.selected_source_id in self.master_ids:
                        self.select_source_in_table(self.selected_source_id)

                    self.update_overlay()

        if not found_detected:
            # No detected source nearby
            self.selected_source_id = None
            self.selected_phot_row = None
            self.selected_phot_file = None
            self.selected_label.setText(f"No detection at ({x:.1f}, {y:.1f}) - click on a circled star")
            self.update_overlay()

    def _on_viewer_press(self, x: float, y: float, btn: int):
        if btn == int(Qt.LeftButton):
            if self._roi_mode:
                self._roi_drag_start = (x, y)
                self._roi_preview = (x, y, 0.0)
                self.update_overlay()
            else:
                self.on_click_xy(x, y)

    def _on_viewer_release(self, x: float, y: float, btn: int):
        if btn == int(Qt.LeftButton) and self._roi_mode and self._roi_drag_start is not None:
            x0, y0 = self._roi_drag_start
            r_px = float(np.hypot(x - x0, y - y0))
            self._roi_drag_start = None
            self._roi_preview = None
            if r_px >= 1.0:
                roi = self._pixels_to_sky_roi(x0, y0, r_px)
                if roi is not None:
                    self._roi_circle = roi
                    self._save_roi()
                    self._update_roi_info_label()
                else:
                    self.log("ROI: no WCS in current frame header, cannot convert to sky coords")
            self.btn_set_roi.blockSignals(True)
            self.btn_set_roi.setChecked(False)
            self.btn_set_roi.setText("Set CMD ROI")
            self.btn_set_roi.blockSignals(False)
            self._roi_mode = False
            self.update_overlay()

    def _on_viewer_move(self, x: float, y: float, val: float):
        self.hover_xy = (x, y)
        if self._roi_mode and self._roi_drag_start is not None:
            x0, y0 = self._roi_drag_start
            r_px = float(np.hypot(x - x0, y - y0))
            self._roi_preview = (x0, y0, r_px)
            self.update_overlay()

    def save_master_ids(self, log_action: str = None):
        output_dir = step9_selection_dir(self.params.P.result_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        master_path = output_dir / "master_star_ids.csv"
        sid2id_path = output_dir / "sourceid_to_ID.csv"
        backup_path = output_dir / "master_star_ids.orig.csv"
        if master_path.exists() and (not backup_path.exists()):
            try:
                backup_path.write_text(master_path.read_text(encoding="utf-8"), encoding="utf-8")
            except Exception:
                pass

        # Build stable-ID rows and persist mapping.
        rows = []
        for sid in sorted(int(s) for s in self.master_ids):
            fixed_id = self._ensure_stable_id(sid)
            g_mag = self._get_gmag_for_source(sid)
            rows.append({
                "ID": int(fixed_id),
                "source_id": sid,
                "g_mag": g_mag if np.isfinite(g_mag) else None
            })
        df = pd.DataFrame(rows)
        if len(df):
            df = df.sort_values(["ID", "source_id"]).reset_index(drop=True)
        df.to_csv(master_path, index=False)
        if len(df):
            df[["source_id", "ID"]].to_csv(sid2id_path, index=False)
        self._master_ids_persisted = True

        # Only log save summary if no specific action (e.g., on load or undo)
        if log_action is None:
            n_sources = len(self.master_ids)
            self.log(f"Saved {n_sources} sources to {master_path.name}")

        self.save_state()
        self.update_master_table()
        self.update_overlay()
        self.update_navigation_buttons()

    def _get_gmag_for_source(self, source_id: int) -> float:
        """Get G magnitude for a source_id from Gaia catalog"""
        sid = int(source_id)
        val = self._gaia_gmag_map.get(sid, np.nan)
        if np.isfinite(val):
            return float(val)
        return float(self.master_gmag_map.get(int(source_id), np.nan))

    def _get_wcs(self) -> WCS | None:
        """Get WCS from header or .wcs file"""
        # First try header
        if self.header is not None:
            try:
                w = WCS(self.header, relax=True)
                if w.has_celestial:
                    return w
            except Exception:
                pass

        # Try .wcs file
        if self.current_filename:
            if self.use_cropped:
                fits_path = step2_cropped_dir(self.params.P.result_dir) / self.current_filename
            else:
                fits_path = self.params.P.data_dir / self.current_filename
            wcs_path = fits_path.with_suffix(".wcs")
            if wcs_path.exists():
                try:
                    from astropy.io import fits as afits
                    with afits.open(wcs_path) as whdul:
                        w = WCS(whdul[0].header, relax=True)
                        if w.has_celestial:
                            return w
                except Exception:
                    pass

        return None


    def open_parameters_dialog(self):
        dialog, layout, buttons = build_scroll_param_dialog(
            self, "Editor Parameters",
            info_text="Adjust star editor interaction parameters. Changes apply immediately.",
            size=(460, 420),
        )

        form = QFormLayout()

        self.param_search = QDoubleSpinBox()
        self.param_search.setRange(1.0, 50.0)
        self.param_search.setValue(float(getattr(self.params.P, "search_radius_px", 7.0)))
        form.addRow("Search Radius (px):", self.param_search)

        self.param_box = QSpinBox()
        self.param_box.setRange(10, 2000)
        self.param_box.setValue(int(getattr(self.params.P, "bulk_drop_box_px", 200)))
        form.addRow("Remove Box Size (px):", self.param_box)

        self.param_gaia_sep = QDoubleSpinBox()
        self.param_gaia_sep.setRange(0.1, 10.0)
        self.param_gaia_sep.setDecimals(1)
        self.param_gaia_sep.setValue(float(getattr(self.params.P, "gaia_add_max_sep_arcsec", 2.0)))
        form.addRow("Gaia Add Max Sep (\"):", self.param_gaia_sep)

        self.param_mem_overlay = QCheckBox("Enable membership color overlay")
        self.param_mem_overlay.setChecked(bool(getattr(self.params.P, "step8_membership_overlay_enable", False)))
        form.addRow("Membership Overlay:", self.param_mem_overlay)

        self.param_mem_thr = QDoubleSpinBox()
        self.param_mem_thr.setRange(0.0, 1.0)
        self.param_mem_thr.setDecimals(2)
        self.param_mem_thr.setSingleStep(0.05)
        self.param_mem_thr.setValue(float(getattr(self.params.P, "step8_membership_threshold", 0.5)))
        form.addRow("Membership P threshold:", self.param_mem_thr)

        layout.addLayout(form)
        layout.addStretch(1)
        buttons.accepted.connect(lambda: self.save_parameters(dialog))
        buttons.rejected.connect(dialog.reject)
        dialog.exec_()

    def save_parameters(self, dialog):
        self.params.P.search_radius_px = self.param_search.value()
        self.params.P.bulk_drop_box_px = self.param_box.value()
        self.params.P.gaia_add_max_sep_arcsec = self.param_gaia_sep.value()
        self.params.P.step8_membership_overlay_enable = self.param_mem_overlay.isChecked()
        self.params.P.step8_membership_threshold = self.param_mem_thr.value()
        self._membership_loaded = False
        self._ensure_membership_map(force=True)
        self.update_overlay()
        self.save_state()
        QMessageBox.information(dialog, "Success", "Parameters saved!")
        dialog.accept()

    def validate_step(self) -> bool:
        return (step9_selection_dir(self.params.P.result_dir) / "master_star_ids.csv").exists()

    # ── Color legend floating window ─────────────────────────────────────────

    def _build_color_window(self) -> QWidget:
        """Build the overlay color/visibility window (floating, Qt.Tool)."""
        win = QWidget(None, Qt.Tool | Qt.WindowTitleHint | Qt.WindowCloseButtonHint)
        win.setWindowTitle("Overlay Colors")
        win.setFixedWidth(200)
        win.closeEvent = lambda e: (self.btn_colors.setChecked(False), e.accept())

        layout = QVBoxLayout(win)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(4)

        hint = QLabel("Left-click: toggle on/off\nRight-click: change color")
        hint.setStyleSheet("color: #888; font-size: 8pt;")
        layout.addWidget(hint)

        entries = [
            ("gaia",     "In master · Gaia"),
            ("local",    "In master · Local"),
            ("forced",   "In master · Forced"),
            ("removed",  "Removed by D-key (press A to restore)"),
            ("selected", "Selected"),
        ]
        for key, label_text in entries:
            row_w = QWidget()
            row = QHBoxLayout(row_w)
            row.setContentsMargins(0, 2, 0, 2)
            row.setSpacing(8)

            btn = QPushButton()
            btn.setFixedSize(24, 24)
            btn.setCheckable(True)
            btn.setChecked(True)
            btn.setToolTip("Left-click: toggle  |  Right-click: change color")
            btn.setContextMenuPolicy(Qt.CustomContextMenu)
            btn.customContextMenuRequested.connect(lambda _, k=key: self._on_pick_color(k))
            btn.toggled.connect(lambda on, k=key: self._toggle_layer_visible(k, on))
            self._color_btns[key] = btn
            self._refresh_color_btn(key)

            lbl = QLabel(label_text)
            row.addWidget(btn)
            row.addWidget(lbl)
            row.addStretch()
            layout.addWidget(row_w)

        layout.addStretch()
        return win

    def _refresh_color_btn(self, key: str):
        """Update a color button's appearance based on current color + visibility."""
        btn = self._color_btns.get(key)
        if btn is None:
            return
        visible = self._overlay_visible.get(key, True)
        color = self._overlay_colors.get(key, "#ffffff")
        if visible:
            btn.setStyleSheet(
                f"QPushButton {{ background-color: {color}; border: 2px solid #aaa; border-radius: 3px; }}"
                f"QPushButton:hover {{ border: 2px solid #fff; }}"
            )
        else:
            btn.setStyleSheet(
                "QPushButton { background-color: #333; border: 2px solid #555; border-radius: 3px; "
                "color: #666; }"
            )

    def _toggle_layer_visible(self, key: str, on: bool):
        """Toggle a specific overlay layer on/off."""
        self._overlay_visible[key] = on
        self._refresh_color_btn(key)
        self.update_overlay()

    def _on_pick_color(self, key: str):
        cur = QColor(self._overlay_colors.get(key, "#ffffff"))
        col = QColorDialog.getColor(cur, self, f"Pick color — {key}")
        if not col.isValid():
            return
        self._overlay_colors[key] = col.name()
        self._refresh_color_btn(key)
        self.update_overlay()

    def _on_members_toggled(self, checked: bool):
        self._show_members = bool(checked)
        self.update_overlay()

    def _toggle_color_panel(self, checked: bool):
        if self._color_window is None:
            return
        if checked:
            # Position to the left of this widget
            top_left = self.mapToGlobal(QPoint(0, 0))
            win_w = self._color_window.sizeHint().width() or 200
            x = max(0, top_left.x() - win_w - 8)
            y = top_left.y()
            self._color_window.move(x, y)
            self._color_window.show()
            self._color_window.raise_()
        else:
            self._color_window.hide()

    @staticmethod
    def _apply_color_btn_style(btn: QPushButton, color: str):
        """Apply the current swatch color to a button."""
        btn.setStyleSheet(
            f"QPushButton {{ background-color: {color}; border: 2px solid #aaa; border-radius: 3px; }}"
            f"QPushButton:hover {{ border: 2px solid #fff; }}"
        )

    # ────────────────────────────────────────────────────────────────────────

    def save_state(self):
        state_data = {
            "search_radius_px": getattr(self.params.P, "search_radius_px", 7.0),
            "bulk_drop_box_px": getattr(self.params.P, "bulk_drop_box_px", 200),
            "gaia_add_max_sep_arcsec": getattr(self.params.P, "gaia_add_max_sep_arcsec", 2.0),
            "step8_membership_overlay_enable": getattr(self.params.P, "step8_membership_overlay_enable", False),
            "step8_membership_threshold": getattr(self.params.P, "step8_membership_threshold", 0.5),
        }
        self.project_state.store_step_data("master_id_editor", state_data)

    def restore_state(self):
        state_data = self.project_state.get_step_data("master_id_editor")
        if state_data:
            for key, val in state_data.items():
                if hasattr(self.params.P, key):
                    setattr(self.params.P, key, val)

    # ------------------------------------------------------------------
    # CMD ROI helpers
    # ------------------------------------------------------------------
    def _roi_path(self) -> Path:
        return step9_selection_dir(self.params.P.result_dir) / "cmd_roi.json"

    def _load_roi(self):
        p = self._roi_path()
        try:
            if p.exists():
                self._roi_circle = json.loads(p.read_text())
            else:
                self._roi_circle = None
        except Exception:
            self._roi_circle = None
        self._update_roi_info_label()

    def _save_roi(self):
        p = self._roi_path()
        try:
            p.parent.mkdir(parents=True, exist_ok=True)
            if self._roi_circle:
                p.write_text(json.dumps(self._roi_circle))
            elif p.exists():
                p.unlink()
        except Exception:
            pass

    def _pixels_to_sky_roi(self, cx_px: float, cy_px: float, r_px: float) -> dict | None:
        """Convert pixel-space ROI to sky ROI {ra_deg, dec_deg, radius_arcsec} using current frame WCS."""
        if self.header is None:
            return None
        try:
            wcs = WCS(self.header, naxis=2)
            ra, dec = wcs.all_pix2world([[cx_px, cy_px]], 0)[0]
            # pixel scale: sqrt(|det(CD matrix)|) in deg/px → arcsec/px
            ps_deg = float(np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix))))
            ps_arcsec = ps_deg * 3600.0
            if ps_arcsec <= 0:
                ps_arcsec = float(getattr(self.params.P, "pixel_scale_arcsec", 1.0))
            return {
                "ra_deg": float(ra),
                "dec_deg": float(dec),
                "radius_arcsec": r_px * ps_arcsec,
            }
        except Exception:
            return None

    def _sky_roi_to_pixels(self) -> tuple[float, float, float] | None:
        """Convert stored sky ROI to pixel coords for the current frame. Returns (cx, cy, r_px) or None."""
        if self._roi_circle is None or self.header is None:
            return None
        try:
            wcs = WCS(self.header, naxis=2)
            ra = self._roi_circle["ra_deg"]
            dec = self._roi_circle["dec_deg"]
            r_arcsec = self._roi_circle["radius_arcsec"]
            cx_px, cy_px = wcs.all_world2pix([[ra, dec]], 0)[0]
            ps_deg = float(np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix))))
            ps_arcsec = ps_deg * 3600.0
            if ps_arcsec <= 0:
                ps_arcsec = float(getattr(self.params.P, "pixel_scale_arcsec", 1.0))
            r_px = r_arcsec / ps_arcsec
            return float(cx_px), float(cy_px), float(r_px)
        except Exception:
            return None

    def _update_roi_info_label(self):
        if self._roi_circle:
            ra = self._roi_circle["ra_deg"]
            dec = self._roi_circle["dec_deg"]
            r = self._roi_circle["radius_arcsec"]
            self.roi_info_label.setText(f"ROI: RA={ra:.4f} Dec={dec:.4f}  r={r:.0f}\"")
            self.roi_info_label.setStyleSheet("QLabel { color: #00E5FF; font-size: 9pt; }")
        else:
            self.roi_info_label.setText("No ROI set")
            self.roi_info_label.setStyleSheet("QLabel { color: #90A4AE; font-size: 9pt; }")

    def _redraw_roi_patch(self):
        """Refresh the ROI overlay circle for the current frame."""
        self.update_overlay()

    def _on_set_roi_toggled(self, checked: bool):
        self._roi_mode = checked
        if checked:
            self.btn_set_roi.setText("Cancel (click+drag to draw)")
            self._roi_drag_start = None
        else:
            self.btn_set_roi.setText("Set CMD ROI")
            self._roi_drag_start = None
            self._roi_preview = None
            self.update_overlay()

    def _on_clear_roi(self):
        self._roi_circle = None
        self._save_roi()
        self._update_roi_info_label()
        self._redraw_roi_patch()

    # ------------------------------------------------------------------

    def show_log_window(self):
        self.log_window.show()
        self.log_window.raise_()
        self.log_window.activateWindow()

    def step_frame(self, delta: int):
        if not self.file_list:
            return
        idx = self.file_combo.currentIndex()
        if idx < 0:
            return
        current_fname = self.file_combo.itemText(idx)
        current_filter = self._get_filter_for_file(current_fname)
        if not current_filter:
            current_filter = ""

        filter_indices = []
        for i, fname in enumerate(self.file_list):
            if self._get_filter_for_file(fname) == current_filter:
                filter_indices.append(i)

        if len(filter_indices) <= 1:
            new_idx = (idx + delta) % len(self.file_list)
        else:
            try:
                pos = filter_indices.index(idx)
            except ValueError:
                pos = 0
            new_idx = filter_indices[(pos + delta) % len(filter_indices)]

        if new_idx != idx:
            # Disconnect signal temporarily to call load_and_display with quick_switch flag
            self.file_combo.blockSignals(True)
            self.file_combo.setCurrentIndex(new_idx)
            self.file_combo.blockSignals(False)
            self.load_and_display(quick_switch=True)

    def keyPressEvent(self, event):
        if event.key() == Qt.Key_G:
            self.show_radial_profile()
            return
        if event.key() == Qt.Key_BracketLeft or event.key() == Qt.Key_Comma:
            self.step_frame(-1)
            return
        if event.key() == Qt.Key_BracketRight or event.key() == Qt.Key_Period:
            self.step_frame(1)
            return
        super().keyPressEvent(event)

    def show_radial_profile(self):
        """Show radial profile at mouse hover position (G key) - updates if already open"""
        if self.image_data is None:
            self.log("No image loaded for radial profile")
            return

        # Use hover position (preferred) or last click position
        if self.hover_xy is not None:
            x, y = self.hover_xy
        elif self.last_click_xy is not None:
            x, y = self.last_click_xy
        else:
            self.log("Move mouse over image first for radial profile")
            return
        frame = self.current_filename or "?"

        try:
            from astropy.stats import sigma_clipped_stats
            from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

            xc, yc = float(x), float(y)

            # Calculate radial profile
            rmax = 50
            ry0 = max(0, int(yc - rmax))
            ry1 = min(self.image_data.shape[0], int(yc + rmax))
            rx0 = max(0, int(xc - rmax))
            rx1 = min(self.image_data.shape[1], int(xc + rmax))

            region = self.image_data[ry0:ry1, rx0:rx1]
            yy, xx = np.mgrid[ry0:ry1, rx0:rx1]
            rr = np.sqrt((xx - xc)**2 + (yy - yc)**2)

            # Background subtraction
            _, reg_median, _ = sigma_clipped_stats(region, sigma=3.0)
            region_sub = region - reg_median

            # Radial bins
            dr = 0.5
            edges = np.arange(0, rmax, dr)
            centers = 0.5 * (edges[:-1] + edges[1:])
            profile = np.full_like(centers, np.nan)

            for i in range(len(centers)):
                mask = (rr >= edges[i]) & (rr < edges[i+1])
                if np.any(mask):
                    vals = region_sub[mask]
                    vals = vals[np.isfinite(vals)]
                    if vals.size > 0:
                        profile[i] = np.mean(vals)

            # Create dialog if not exists
            if not hasattr(self, 'radial_dialog') or self.radial_dialog is None or not self.radial_dialog.isVisible():
                self.radial_dialog = QDialog(self)
                self.radial_dialog.setWindowTitle("Radial Profile (G key to update)")
                self.radial_dialog.resize(600, 400)

                layout = QVBoxLayout(self.radial_dialog)
                self.prof_fig = Figure(figsize=(6, 4))
                self.prof_canvas = FigureCanvas(self.prof_fig)
                self.prof_ax = self.prof_fig.add_subplot(111)
                layout.addWidget(NavigationToolbar(self.prof_canvas, self.radial_dialog))
                layout.addWidget(self.prof_canvas)

            # Plot (always update)
            self.prof_ax.clear()
            self.prof_ax.plot(centers, profile, 'o-', color='steelblue', markersize=3, linewidth=1)
            self.prof_ax.axhline(0, color='red', linestyle='--', linewidth=1, alpha=0.5)
            self.prof_ax.set_xlabel('Radius (pixels)')
            self.prof_ax.set_ylabel('Pixel Value - Background (ADU)')
            self.prof_ax.set_title(f'{frame} | Position ({int(xc)}, {int(yc)})')
            self.prof_ax.grid(True, alpha=0.3)

            # Estimate FWHM
            fwhm_note = ""
            peak = np.nanmax(profile) if np.isfinite(profile).any() else 0
            if peak > 0:
                half = 0.5 * peak
                idx = np.where((profile[:-1] >= half) & (profile[1:] < half))[0]
                if len(idx) > 0:
                    i = idx[0]
                    x1_p, y1_p = centers[i], profile[i]
                    x2_p, y2_p = centers[i+1], profile[i+1]
                    if y1_p != y2_p:
                        r_half = x1_p + (half - y1_p) * (x2_p - x1_p) / (y2_p - y1_p)
                        fwhm_px = 2.0 * r_half
                        pixscale = float(getattr(self.params.P, "pixel_scale_arcsec", 0.4))
                        fwhm_arcsec = fwhm_px * pixscale
                        self.prof_ax.axvline(r_half, color='orange', linestyle='--', linewidth=2)
                        fwhm_note = f'FWHM: {fwhm_arcsec:.2f}" ({fwhm_px:.2f} px)'
                        self.prof_ax.text(0.02, 0.95, fwhm_note,
                                         transform=self.prof_ax.transAxes, ha='left', va='top',
                                         fontsize=10, bbox=dict(facecolor='white', alpha=0.8))
                        self.log(f"[{frame}] Radial profile at ({int(xc)},{int(yc)}): FWHM={fwhm_arcsec:.2f}\" ({fwhm_px:.2f}px)")

            self.prof_fig.tight_layout()
            self.prof_canvas.draw()
            self.radial_dialog.show()
            self.radial_dialog.raise_()

        except Exception as e:
            self.log(f"[{frame}] Radial profile error: {e}")
