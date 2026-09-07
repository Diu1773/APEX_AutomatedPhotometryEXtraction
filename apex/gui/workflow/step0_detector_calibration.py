"""Step 0: Detector Calibration (optional, off-chain).

GUI front-end for :mod:`apex.analysis.calibration`.  The user points at one or
more raw folders; the window auto-scans them, reads FITS headers on a background
thread, classifies each frame as bias/dark/flat/light and groups them by night /
exposure / temperature / filter (:mod:`apex.analysis.calibration_scan`), and
shows the result as a tree.  On Run it builds master frames per group, matches
each light to the right dark (exposure+temperature) and flat (filter), and
writes calibrated FITS + masters + a QC summary under ``step0_calibration/``.

Off-chain: never sets ``current_step`` / touches ``completed_steps``.  On
success records ``mark_calibration("done")`` and offers the calibrated folder as
the File-Selection input; Skip records ``"skipped"``.  Heavy work runs on
QThreads (workers emit signals; the main thread touches widgets).
"""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Dict, List, Optional

from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtWidgets import (
    QCheckBox, QComboBox, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
    QFormLayout, QGroupBox, QHBoxLayout, QInputDialog, QLabel, QLineEdit,
    QMenu, QMessageBox, QProgressBar, QPushButton, QScrollArea, QSpinBox,
    QTableWidget, QTableWidgetItem, QTreeWidget, QTreeWidgetItem, QVBoxLayout,
    QWidget,
)

from apex.analysis import calibration_scan as scan
from apex.analysis.calibration import CalibrationOptions
from apex.analysis.calibration_run import run_calibration
from apex.analysis.calibration_scan import FrameInfo
from apex.config.calibration_section import calibration_toml_sections
from apex.gui.layout_rules import FittedDialog
from apex.gui.theme import ICON, Tokens, style_button
from apex.gui.tools.tool_window_base import ToolWindowBase
from apex.utils.astro_utils import (
    filter_aliases_in_effect,
    normalize_filter_name,
    register_filter_aliases,
)
from apex.utils.param_file import update_param_file
from apex.utils.step_paths import step0_calibration_dir


# ---------------------------------------------------------------------------
# Workers
# ---------------------------------------------------------------------------

class _ScanWorker(QThread):
    progress = pyqtSignal(int, int, str)
    logline = pyqtSignal(str)
    done = pyqtSignal(list)                       # List[FrameInfo]

    def __init__(self, roots: List[str], tz_offset_hours: Optional[float] = None):
        super().__init__()
        self.roots = roots
        # Fallback local reference for the observing-night split when the FITS
        # header carries no site longitude (apex.utils.night_utils).
        self.tz_offset_hours = tz_offset_hours

    def run(self):
        frames: List[FrameInfo] = []
        for root in self.roots:
            frames.extend(scan.scan_folder(
                root,
                progress=lambda i, t, m: self.progress.emit(i, t, m),
                stop=self.isInterruptionRequested,
                tz_offset_hours=self.tz_offset_hours,
                warn=self.logline.emit,
            ))
        self.done.emit(frames)


class _CalibrationWorker(QThread):
    progress = pyqtSignal(int, int, str)
    logline = pyqtSignal(str)
    finished_ok = pyqtSignal(dict)
    failed = pyqtSignal(str)

    def __init__(self, frames: List[FrameInfo], night: str, out_dir: Path,
                 opts: CalibrationOptions):
        super().__init__()
        self.frames = frames
        self.night = night
        self.out_dir = Path(out_dir)
        self.opts = opts

    def run(self):
        try:
            self._run()
        except Exception as exc:                  # pragma: no cover - surfaced to UI
            self.failed.emit(f"{type(exc).__name__}: {exc}")

    def _run(self):
        # All the work lives in the Qt-free core, which the headless pipeline
        # step and the reprocess scripts call too — one code path, so a batch
        # run cannot drift from what this window shows.
        summary = run_calibration(
            self.frames, self.night, self.out_dir, self.opts,
            progress=lambda done, total, msg: self.progress.emit(done, total, msg),
            log=self.logline.emit,
            stop=self.isInterruptionRequested,
        )
        self.finished_ok.emit(summary)


# ---------------------------------------------------------------------------
# Window
# ---------------------------------------------------------------------------

class DetectorCalibrationWindow(ToolWindowBase):
    """Optional off-chain Step 0 window (detector calibration)."""

    _TYPE_ORDER = ("bias", "dark", "flat", "light", scan.TYPE_UNKNOWN)

    def __init__(self, params, project_state, parent=None):
        super().__init__("Step 0: Detector Calibration", params=params,
                         project_state=project_state, parent=parent,
                         min_size=(560, 640))
        self._frames: List[FrameInfo] = []
        self._roots: List[str] = []
        self._scan_worker: Optional[_ScanWorker] = None
        self._worker: Optional[_CalibrationWorker] = None
        # Manual type/filter/night fixes keyed by frame path, persisted so a
        # re-scan of the same folder keeps them.
        self._overrides: Dict[str, Dict] = {}
        # All calibration parameters live in one mutable holder, edited via the
        # Parameters dialog (AIPPI-style single settings panel) and persisted
        # back to the TOML's [calibration] table so they survive a restart.
        P = getattr(params, "P", None)
        self._settings: Dict = asdict(
            CalibrationOptions.from_mapping(getattr(P, "calibration", None)))
        # Seed cosmetic's noise model from the CONFIGURED (measured) gain/read
        # noise, never the FITS EGAIN header (unreliable on this camera).
        if P is not None:
            try:
                self._settings["gain"] = float(getattr(P, "gain_e_per_adu", 1.0) or 1.0)
                self._settings["readnoise"] = float(
                    getattr(P, "rdnoise_e", getattr(P, "readnoise", 6.5)) or 6.5)
            except (TypeError, ValueError):
                pass
        self.add_header_action(f"{ICON['params']} Parameters", self._open_parameters,
                               tooltip="Calibration parameters")
        self.add_log_dock(popup=True)
        self.content_layout.addWidget(self._build_body())

    # -- UI ----------------------------------------------------------------

    def _build_body(self) -> QWidget:
        body = QWidget()
        v = QVBoxLayout(body)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(Tokens.S3)

        intro = QLabel(
            "Scan a raw folder — APEX reads the FITS headers, sorts frames into "
            "bias/dark/flat/light by night, and calibrates the science frames "
            "(mono CCD). Optional — skip if frames are already calibrated."
        )
        intro.setWordWrap(True)
        intro.setStyleSheet(f"QLabel {{ color: {Tokens.TEXT_MUTED}; }}")
        v.addWidget(intro)

        # scan controls
        scan_row = QHBoxLayout()
        scan_row.setSpacing(Tokens.GAP)
        self.btn_scan = QPushButton("Scan folder…")
        style_button(self.btn_scan, "primary", height=Tokens.H_BUTTON)
        self.btn_scan.clicked.connect(self._on_scan)
        self.btn_add = QPushButton("Add folder…")
        style_button(self.btn_add, height=Tokens.H_BUTTON)
        self.btn_add.clicked.connect(lambda: self._on_scan(add=True))
        self.scan_label = QLabel("No folder scanned yet.")
        self.scan_label.setStyleSheet(f"QLabel {{ color: {Tokens.TEXT_MUTED}; }}")
        scan_row.addWidget(self.btn_scan)
        scan_row.addWidget(self.btn_add)
        scan_row.addWidget(self.scan_label, 1)
        v.addLayout(scan_row)

        # classified tree
        self.tree = QTreeWidget()
        self.tree.setColumnCount(3)
        self.tree.setHeaderLabels(["Detected frames", "Count", "Calibration"])
        self.tree.setRootIsDecorated(True)
        self.tree.header().setStretchLastSection(True)
        self.tree.setColumnWidth(0, 260)
        self.tree.setColumnWidth(1, 60)
        self.tree.setContextMenuPolicy(Qt.CustomContextMenu)
        self.tree.customContextMenuRequested.connect(self._on_tree_menu)
        self.tree.setToolTip("Right-click a row to correct its frame type, "
                             "filter or night when the FITS header is wrong.")
        v.addWidget(self.tree, 1)

        # 필터 이름 — 관측소마다 FILTER 헤더를 자기 방식으로 적는다. 스캔이
        # 찾아낸 값을 그대로 보여 주고 그 자리에서 고치게 한다. 설정 파일을
        # 직접 열게 하는 것은 자료를 앞에 둔 사람에게 할 일이 아니다.
        self.filter_box = QGroupBox("Filter names in these headers")
        fb = QVBoxLayout(self.filter_box)
        fb.setContentsMargins(Tokens.S3, Tokens.S3, Tokens.S3, Tokens.S3)
        fb.setSpacing(Tokens.GAP)
        hint = QLabel(
            "Left is what the FITS header says, right is the key APEX uses. "
            "Edit the right column when this observatory spells a filter its "
            "own way. Johnson is upper case (B V R I), SDSS lower (g r i z).")
        hint.setWordWrap(True)
        fb.addWidget(hint)
        self.filter_table = QTableWidget(0, 3)
        self.filter_table.setHorizontalHeaderLabels(
            ["FILTER header", "Frames", "APEX key"])
        self.filter_table.verticalHeader().setVisible(False)
        self.filter_table.setEditTriggers(
            QTableWidget.DoubleClicked | QTableWidget.SelectedClicked)
        self.filter_table.horizontalHeader().setStretchLastSection(True)
        self.filter_table.setColumnWidth(0, 200)
        self.filter_table.setColumnWidth(1, 70)
        self.filter_table.itemChanged.connect(self._on_filter_edited)
        fb.addWidget(self.filter_table)
        self.filter_box.setVisible(False)
        v.addWidget(self.filter_box)

        # process-night row (all tunables live in the Parameters dialog)
        proc_row = QHBoxLayout()
        proc_row.setSpacing(Tokens.GAP)
        proc_row.addWidget(QLabel("Process night"))
        self.cmb_night = QComboBox()
        self.cmb_night.setMinimumHeight(Tokens.H_BUTTON)
        self.cmb_night.currentIndexChanged.connect(lambda _i: self._rebuild_tree())
        proc_row.addWidget(self.cmb_night, 1)
        v.addLayout(proc_row)

        # output
        out_box = QGroupBox("Output folder")
        ov = QHBoxLayout(out_box)
        ov.setContentsMargins(Tokens.S3, Tokens.S3, Tokens.S3, Tokens.S3)
        self.out_edit = QLineEdit(str(self._default_out_dir()))
        self.out_edit.setMinimumHeight(Tokens.H_BUTTON)
        ov.addWidget(self.out_edit, 1)
        v.addWidget(out_box)

        # progress
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimumHeight(20)
        self.progress_bar.setTextVisible(True)
        v.addWidget(self.progress_bar)

        # actions
        actions = QHBoxLayout()
        self.btn_run = QPushButton("Run Calibration")
        style_button(self.btn_run, "primary", height=Tokens.H_ACTION)
        self.btn_run.setEnabled(False)
        self.btn_run.clicked.connect(self._on_run)
        self.btn_skip = QPushButton("Skip")
        style_button(self.btn_skip, height=Tokens.H_ACTION)
        self.btn_skip.clicked.connect(self._on_skip)
        self.btn_close = QPushButton("Close")
        style_button(self.btn_close, height=Tokens.H_ACTION)
        self.btn_close.clicked.connect(self.close)
        actions.addWidget(self.btn_skip)
        actions.addStretch(1)
        actions.addWidget(self.btn_close)
        actions.addWidget(self.btn_run)
        v.addLayout(actions)
        return body

    # -- helpers -----------------------------------------------------------

    def _site_tz_offset(self) -> Optional[float]:
        """Configured site tz offset, the civil fallback for the night split.

        Used only when a frame's header has no site longitude; 0.0 (the config
        default, indistinguishable from a genuine Greenwich site) is treated as
        unset by :mod:`apex.utils.night_utils`.
        """
        P = getattr(self.params, "P", None)
        try:
            return float(getattr(P, "site_tz_offset_hours", 0.0) or 0.0)
        except (TypeError, ValueError):
            return None

    def _default_out_dir(self) -> Path:
        result_dir = getattr(getattr(self.params, "P", None), "result_dir", "") or ""
        base = Path(result_dir) if result_dir else Path.cwd()
        return step0_calibration_dir(base)

    def _current_options(self) -> CalibrationOptions:
        return CalibrationOptions(**self._settings)

    def _persist_settings(self) -> None:
        """Write the settings back to the TOML's [calibration] table.

        Without this the whole panel resets every time the window closes; the
        section was also inert until the loader learned to read it.
        """
        P = getattr(self.params, "P", None)
        if P is not None:
            try:
                P.calibration = dict(self._settings)
            except Exception:
                pass
        # gain/readnoise come from [instrument]; writing them here too would
        # create a second source of truth for the same numbers.
        payload = {k: v for k, v in self._settings.items()
                   if k not in ("gain", "readnoise")}
        saved = update_param_file(
            getattr(self.params, "param_file", None),
            calibration_toml_sections(payload),
            params=self.params,
        )
        self.log("Parameters updated and saved to the parameter file."
                 if saved else
                 "Parameters updated (this session only — no writable parameter file).")

    # -- parameters dialog (single settings panel) -------------------------

    def _open_parameters(self):
        s = self._settings
        dlg = FittedDialog(self)
        dlg.setWindowTitle("Detector Calibration — Parameters")
        outer = QVBoxLayout(dlg)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.NoFrame)
        inner = QWidget()
        col = QVBoxLayout(inner)
        col.setSpacing(Tokens.S3)

        def _dspin(lo, hi, step, dec, val):
            w = QDoubleSpinBox(); w.setRange(lo, hi); w.setSingleStep(step)
            w.setDecimals(dec); w.setValue(float(val)); w.setMinimumHeight(Tokens.H_BUTTON)
            return w

        def _ispin(lo, hi, val):
            w = QSpinBox(); w.setRange(lo, hi); w.setValue(int(val))
            w.setMinimumHeight(Tokens.H_BUTTON); return w

        def _combo(items, val):
            w = QComboBox(); w.addItems(items); w.setCurrentText(str(val))
            w.setMinimumHeight(Tokens.H_BUTTON); return w

        def _lbl(text):                       # fixed-width label keeps the form
            lab = QLabel(text); lab.setFixedWidth(150); return lab  # column tight

        def _group(title):
            g = QGroupBox(title)
            f = QFormLayout(g)
            f.setContentsMargins(Tokens.S3, Tokens.S3, Tokens.S3, Tokens.S3)
            f.setSpacing(Tokens.GAP)
            col.addWidget(g)
            return f

        # Master combine
        f = _group("Master combine")
        w_method = _combo(["median", "mean", "sigmaclip_mean"], s["combine_method"])
        w_slow = _dspin(0.5, 10.0, 0.5, 1, s["sigma_low"])
        w_shigh = _dspin(0.5, 10.0, 0.5, 1, s["sigma_high"])
        w_iter = _ispin(1, 10, s["maxiters"])
        f.addRow(_lbl("Method"), w_method)
        f.addRow(_lbl("Sigma low (clip)"), w_slow)
        f.addRow(_lbl("Sigma high (clip)"), w_shigh)
        f.addRow(_lbl("Max iterations"), w_iter)

        # Dark
        f = _group("Dark")
        w_dscale = QCheckBox("Scale master dark by exposure ratio"); w_dscale.setChecked(s["dark_scale"])
        w_dopt = QCheckBox("Optimise scale (noise-min k-fit)"); w_dopt.setChecked(s["dark_optimize"])
        w_dtol = _dspin(0.1, 10.0, 0.1, 1, s["temp_match_tol_c"])
        w_dtol.setToolTip(
            "How far the dark's sensor temperature may sit from the light's.\n"
            "Dark current roughly doubles every ~6 °C, so tighten this if you\n"
            "are chasing faint signal. Exceeding it is always logged."
        )
        w_dstrict = QCheckBox("Refuse a dark outside the tolerance (no dark subtraction)")
        w_dstrict.setChecked(s["strict_temp"])
        f.addRow(w_dscale); f.addRow(w_dopt)
        f.addRow(_lbl("Temp tolerance (°C)"), w_dtol)
        f.addRow(w_dstrict)

        # Flat
        f = _group("Flat")
        w_fmin = _dspin(0.0, 0.5, 0.01, 3, s["flat_min"])
        f.addRow(_lbl("Dead-pixel threshold"), w_fmin)

        # Pedestal
        f = _group("Pedestal")
        w_pmode = _combo(["adaptive", "none", "fixed"], s["pedestal_mode"])
        w_pval = _dspin(0.0, 10000.0, 10.0, 1, s["pedestal_value"])
        w_pmax = _dspin(0.0, 50000.0, 100.0, 1, s["pedestal_max"])
        f.addRow(_lbl("Mode"), w_pmode)
        f.addRow(_lbl("Fixed value (DN)"), w_pval)
        f.addRow(_lbl("Adaptive max (DN)"), w_pmax)

        # Overscan
        f = _group("Overscan")
        w_oen = QCheckBox("Subtract overscan"); w_oen.setChecked(s["overscan_enable"])
        w_oedge = _combo(["left", "right", "top", "bottom"], s["overscan_edge"])
        w_owid = _ispin(2, 512, s["overscan_width"])
        w_otrim = QCheckBox("Trim overscan region"); w_otrim.setChecked(s["overscan_trim"])
        f.addRow(w_oen)
        f.addRow(_lbl("Edge"), w_oedge)
        f.addRow(_lbl("Width (px)"), w_owid)
        f.addRow(w_otrim)

        # Cosmetic (cosmic-ray + hot pixel)
        f = _group("Cosmetic (cosmic-ray + hot pixel)")
        w_cc = QCheckBox("Remove cosmic rays + hot pixels (L.A.Cosmic, star-protected)")
        w_cc.setChecked(s["cosmetic_enable"])
        w_crs = _dspin(1.0, 15.0, 0.5, 1, s["cr_sigclip"])
        w_cro = _dspin(1.0, 20.0, 0.5, 1, s["cr_objlim"])
        w_hot = _dspin(2.0, 20.0, 0.5, 1, s["hot_sigma"])
        f.addRow(w_cc)
        f.addRow(_lbl("CR sigma (sigclip)"), w_crs)
        f.addRow(_lbl("CR object limit"), w_cro)
        f.addRow(_lbl("Hot-pixel sigma"), w_hot)
        col.addStretch(1)

        scroll.setWidget(inner)
        outer.addWidget(scroll, 1)
        buttons = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Cancel)
        buttons.accepted.connect(dlg.accept)
        buttons.rejected.connect(dlg.reject)
        outer.addWidget(buttons)

        if dlg.exec_() == FittedDialog.Accepted:
            s.update({
                "combine_method": w_method.currentText(),
                "sigma_low": w_slow.value(), "sigma_high": w_shigh.value(),
                "maxiters": w_iter.value(),
                "dark_scale": w_dscale.isChecked(), "dark_optimize": w_dopt.isChecked(),
                "temp_match_tol_c": w_dtol.value(), "strict_temp": w_dstrict.isChecked(),
                "flat_min": w_fmin.value(),
                "pedestal_mode": w_pmode.currentText(),
                "pedestal_value": w_pval.value(), "pedestal_max": w_pmax.value(),
                "overscan_enable": w_oen.isChecked(), "overscan_edge": w_oedge.currentText(),
                "overscan_width": w_owid.value(), "overscan_trim": w_otrim.isChecked(),
                "cosmetic_enable": w_cc.isChecked(), "cr_sigclip": w_crs.value(),
                "cr_objlim": w_cro.value(), "hot_sigma": w_hot.value(),
            })
            self._persist_settings()

    # -- scan flow ---------------------------------------------------------

    def _on_scan(self, add: bool = False):
        if self._scan_worker is not None and self._scan_worker.isRunning():
            return
        folder = QFileDialog.getExistingDirectory(self, "Select a raw folder to scan")
        if not folder:
            return
        if not add:
            self._frames = []
            self._roots = []
        self._roots.append(folder)
        self.btn_scan.setEnabled(False)
        self.btn_add.setEnabled(False)
        self.scan_label.setText(f"Scanning {folder}…")
        self.log(f"Scanning {folder}")
        self._scan_worker = _ScanWorker([folder], self._site_tz_offset())
        self._scan_worker.progress.connect(self._on_scan_progress)
        self._scan_worker.logline.connect(self.log)
        self._scan_worker.done.connect(self._on_scan_done)
        self._scan_worker.start()

    def _on_scan_progress(self, done: int, total: int, message: str):
        if total > 0:
            self.progress_bar.setMaximum(total)
            self.progress_bar.setValue(done)
        self.progress_bar.setFormat(f"Reading headers… {message} (%p%)")

    def _on_scan_done(self, frames: List[FrameInfo]):
        self.btn_scan.setEnabled(True)
        self.btn_add.setEnabled(True)
        self._frames.extend(frames)
        self._load_overrides()
        self._apply_overrides()
        self.progress_bar.reset()
        self._refresh_scan_label()
        self.log(f"Scanned {len(frames)} new frames ({len(self._frames)} total)")
        self._populate_nights()
        self._rebuild_filter_table()
        self._rebuild_tree()

    # -- filter-name mapping ------------------------------------------------

    def _rebuild_filter_table(self):
        """스캔이 찾은 FILTER 값을 그대로 보여 준다 (많은 것부터)."""
        counts: dict[str, int] = {}
        for f in self._frames:
            # FrameInfo.filt 은 헤더에 적힌 값 그대로다(정규화는 뒤 단계의 일).
            raw = str(f.filt or "").strip()
            if raw:
                counts[raw] = counts.get(raw, 0) + 1
        rows = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0].lower()))

        self.filter_table.blockSignals(True)
        self.filter_table.setSortingEnabled(False)
        try:
            self.filter_table.setRowCount(len(rows))
            for i, (raw, n) in enumerate(rows):
                left = QTableWidgetItem(raw)
                left.setFlags(left.flags() & ~Qt.ItemIsEditable)
                count = QTableWidgetItem(str(n))
                count.setFlags(count.flags() & ~Qt.ItemIsEditable)
                count.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                key = QTableWidgetItem(normalize_filter_name(raw))
                key.setData(Qt.UserRole, raw)
                self.filter_table.setItem(i, 0, left)
                self.filter_table.setItem(i, 1, count)
                self.filter_table.setItem(i, 2, key)
        finally:
            self.filter_table.blockSignals(False)
        # 줄 수에 맞춰 키우되 여섯 줄에서 멈춘다 — 필터가 많은 자료에서 이 표가
        # 창을 다 먹으면 안 된다. 행 높이는 defaultSectionSize 로 잡는다.
        # sizeHintForRow 는 배치 전에 부르면 작은 값을 돌려줘서 표가 잘린다.
        vh = self.filter_table.verticalHeader()
        row_h = vh.defaultSectionSize() or 24
        wanted = (self.filter_table.horizontalHeader().sizeHint().height()
                  + row_h * min(max(len(rows), 1), 6)
                  + 2 * self.filter_table.frameWidth())
        self.filter_table.setMinimumHeight(wanted)
        self.filter_table.setMaximumHeight(wanted)
        self.filter_box.setVisible(bool(rows))

    def _on_filter_edited(self, item):
        if item is None or item.column() != 2:
            return
        raw = str(item.data(Qt.UserRole) or "").strip()
        new_key = str(item.text()).strip()
        if not raw:
            return
        if not new_key:
            # 빈 칸으로 지우면 그 별칭을 없애고 기본 판정으로 되돌린다.
            aliases = {k: v for k, v in self._filter_aliases().items()
                       if k.lower() != raw.lower()}
        else:
            aliases = dict(self._filter_aliases())
            aliases[raw] = new_key
        self._persist_filter_aliases(aliases)
        self.filter_table.blockSignals(True)
        item.setText(normalize_filter_name(raw))
        self.filter_table.blockSignals(False)
        self.log(f"Filter '{raw}' now reads as '{normalize_filter_name(raw)}'")
        self._rebuild_tree()

    def _filter_aliases(self) -> dict:
        P = getattr(self.params, "P", None)
        return dict(getattr(P, "filter_aliases", None) or {})

    def _persist_filter_aliases(self, aliases: dict) -> None:
        # 지운 별칭은 「빠진 키」가 아니라 None 으로 보내야 파일에서도 없어진다.
        # 새 표만 쓰면 파일에는 옛 줄이 그대로 남아, 다음에 열 때 되살아난다.
        previous = self._filter_aliases()
        payload = dict(aliases)
        for key in previous:
            if key not in payload:
                payload[key] = None

        register_filter_aliases(aliases)
        P = getattr(self.params, "P", None)
        if P is not None:
            try:
                P.filter_aliases = dict(aliases)
            except Exception:
                pass
        saved = update_param_file(
            getattr(self.params, "param_file", None),
            {("filters", "aliases"): payload},
            params=self.params,
        )
        if not saved:
            self.log("Filter names updated (this session only — "
                     "no writable parameter file).")

    def _refresh_scan_label(self):
        by_type = {t: sum(1 for f in self._frames if f.ftype == t)
                   for t in self._TYPE_ORDER}
        parts = [f"{t} {by_type[t]}" for t in self._TYPE_ORDER
                 if t != scan.TYPE_UNKNOWN or by_type[t]]
        self.scan_label.setText(f"{len(self._frames)} frames  ·  " + "  ".join(parts))

    # -- manual reclassification -------------------------------------------

    def _overrides_path(self) -> Path:
        return Path(self.out_edit.text().strip() or self._default_out_dir()) \
            / scan.OVERRIDES_FILENAME

    def _load_overrides(self):
        saved = scan.load_overrides(self._overrides_path())
        if saved:
            self._overrides.update(saved)
            self.log(f"Loaded {len(saved)} saved classification override(s)")

    def _apply_overrides(self):
        self._frames = scan.apply_overrides(self._frames, self._overrides)

    def _set_override(self, paths: List[str], **changes):
        for path in paths:
            entry = dict(self._overrides.get(path, {}))
            entry.update({k: v for k, v in changes.items() if v is not None})
            self._overrides[path] = entry
        scan.save_overrides(self._overrides_path(), self._overrides)
        self._apply_overrides()
        self._refresh_scan_label()
        self._populate_nights()
        self._rebuild_tree()

    def _clear_override(self, paths: List[str]):
        removed = sum(1 for p in paths if self._overrides.pop(p, None) is not None)
        if not removed:
            return
        scan.save_overrides(self._overrides_path(), self._overrides)
        # The stored FrameInfo already carries the override, so re-read the
        # headers of just those frames to get their original classification.
        tz = self._site_tz_offset()
        by_path = {p: scan.read_frame_info(p, tz) for p in paths}
        self._frames = [by_path.get(f.path) or f for f in self._frames]
        self._apply_overrides()
        self.log(f"Cleared {removed} override(s)")
        self._refresh_scan_label()
        self._populate_nights()
        self._rebuild_tree()

    def _on_tree_menu(self, pos):
        item = self.tree.itemAt(pos)
        if item is None:
            return
        paths = item.data(0, Qt.UserRole) or []
        if not paths:
            return
        n = len(paths)
        menu = QMenu(self)
        type_menu = menu.addMenu(f"Set frame type  ({n} frame(s))")
        for ftype in scan.FRAME_TYPES:
            action = type_menu.addAction(ftype.capitalize())
            action.triggered.connect(
                lambda _checked, t=ftype: self._set_override(paths, ftype=t))
        act_filter = menu.addAction("Set filter…")
        act_night = menu.addAction("Set night (YYYYMMDD)…")
        menu.addSeparator()
        act_clear = menu.addAction("Clear manual override")
        act_clear.setEnabled(any(p in self._overrides for p in paths))

        chosen = menu.exec_(self.tree.viewport().mapToGlobal(pos))
        if chosen is act_filter:
            current = next((f.filt for f in self._frames if f.path == paths[0]), "")
            value, ok = QInputDialog.getText(self, "Set filter", "Filter name:",
                                             text=current)
            if ok:
                self._set_override(paths, filt=value.strip())
        elif chosen is act_night:
            current = next((f.night for f in self._frames if f.path == paths[0]), "")
            value, ok = QInputDialog.getText(
                self, "Set night", "Observing night (YYYYMMDD):", text=current)
            if ok:
                self._set_override(paths, night=value.strip())
        elif chosen is act_clear:
            self._clear_override(paths)

    def _populate_nights(self):
        self.cmb_night.blockSignals(True)
        self.cmb_night.clear()
        ns = scan.nights(self._frames)
        if ns:
            # Default: process every night at once (dump folders → all done).
            self.cmb_night.addItem(f"All nights ({len(ns)})", "__all__")
            for n in ns:
                self.cmb_night.addItem(f"Night {n}" if n else "Undated", n)
        self.cmb_night.blockSignals(False)
        self.btn_run.setEnabled(bool(ns))

    def _selected_night(self) -> Optional[str]:
        if self.cmb_night.count() == 0:
            return None
        return self.cmb_night.currentData()

    def _rebuild_tree(self):
        self.tree.clear()
        night = self._selected_night()
        expand_all = night == "__all__"
        all_nights = sorted({f.night for f in self._frames}, reverse=True)
        for n in all_nights:
            g = scan.group_for_night(self._frames, n)
            n_lights = sum(1 for f in self._frames if f.ftype == "light" and f.night == n)
            title = (f"Night {n}" if n else "Undated / global pool")
            node = QTreeWidgetItem([f"☾ {title}", "", f"{n_lights} lights" if n_lights else ""])
            self.tree.addTopLevelItem(node)
            expand_unknown = False
            for t in self._TYPE_ORDER:
                items = [f for f in self._frames if f.ftype == t and f.night == n]
                if not items:
                    continue
                if t == scan.TYPE_UNKNOWN:
                    label = "Unclassified"
                    hint = "right-click → Set type"
                    expand_unknown = True
                else:
                    label = t.capitalize()
                    n_master = sum(1 for f in items if getattr(f, "is_master", False))
                    hint = "pre-built master → used directly" if n_master else ""
                type_node = QTreeWidgetItem([label, str(len(items)), hint])
                self._tag_paths(type_node, items)
                node.addChild(type_node)
                # Always pass the resolved pool so each night's light rows show
                # which dark/flat (incl. shared/global) will be used.
                self._add_subgroups(type_node, t, items, g)
                if t == scan.TYPE_UNKNOWN:
                    type_node.setExpanded(True)
            node.setExpanded(expand_all or expand_unknown or n == night)

    @staticmethod
    def _tag_paths(item: QTreeWidgetItem, frames) -> None:
        """Attach the frames a row stands for, so the context menu can act."""
        item.setData(0, Qt.UserRole, [f.path for f in frames])

    def _add_subgroups(self, parent, ftype, items, group):
        if ftype == scan.TYPE_UNKNOWN:
            # One row per file: these need individual attention, and the user
            # has to see the names to know what they are.
            for f in sorted(items, key=lambda x: x.name):
                child = QTreeWidgetItem([f.name, "1", f.geometry_label])
                self._tag_paths(child, [f])
                parent.addChild(child)
            return
        if ftype == "dark":
            keys: Dict = {}
            geoms: Dict = {}
            for f in items:
                key = (round(f.exp, 3), f.temp_bucket)
                keys[key] = keys.get(key, 0) + 1
                geoms.setdefault(key, set()).add(f.geometry_label)
            for (exp, tb), c in sorted(keys.items()):
                labels = sorted(x for x in geoms[(exp, tb)] if x)
                geom = " / ".join(labels) if len(labels) > 1 else (labels[0] if labels else "")
                parent.addChild(QTreeWidgetItem([f"{exp:g}s · {tb}°C", str(c), geom]))
        elif ftype in ("flat", "light"):
            keys = {}
            for f in items:
                keys.setdefault(f.filt or "(none)", 0)
                keys[f.filt or "(none)"] += 1
            for filt, c in sorted(keys.items()):
                detail = ""
                if ftype == "light" and group is not None:
                    detail = self._match_detail(group, filt, items)
                parent.addChild(QTreeWidgetItem([f"filter {filt}", str(c), detail]))
        # bias: no subgroups

    def _match_detail(self, group, filt, items) -> str:
        sample = next((f for f in items if (f.filt or "(none)") == filt), None)
        if sample is None:
            return ""
        tol = float(self._settings.get("temp_match_tol_c", 1.0))
        library = group.get("dark_library")
        dm = scan.match_dark_detail(group["dark"], sample.exp, sample.temp, tol,
                                    light=sample, fallback=library)
        ff = scan.match_flat(group["flat"], sample.filt, light=sample)
        if dm is None:
            # Distinguish "no dark at all" from "none of this geometry".
            blocked = bool(group["dark"]) and scan.match_dark_detail(
                group["dark"], sample.exp, sample.temp, tol,
                fallback=library) is not None
            d = "no dark (geometry mismatch)" if blocked else "no dark"
        else:
            d = f"dark {dm.exp:g}s/{dm.temp_bucket}°C"
            if dm.source == "library" and dm.night:
                d += f" [{dm.night}]"        # not this night's own dark
            # Show the residual mismatch up front — the old preview implied an
            # exact match no matter how far off the nearest dark actually was.
            if dm.delta_temp_c is not None and not dm.within_temp_tol:
                d += f" [ΔT={dm.delta_temp_c:.2f}°C > {tol:g}]"
            if dm.delta_exp_s > 1e-3:
                d += f" [Δexp={dm.delta_exp_s:g}s]"
        fl = f"flat {ff}" if ff else "NO FLAT"
        return f"→ {d}, {fl}"

    # -- run flow ----------------------------------------------------------

    def _on_run(self):
        if self._worker is not None and self._worker.isRunning():
            return
        night = self._selected_night()
        if night is None:
            QMessageBox.warning(self, "Detector Calibration",
                                "Scan a folder containing light frames first.")
            return
        if not any(f.ftype in ("bias", "dark", "flat") for f in self._frames):
            QMessageBox.warning(self, "Detector Calibration",
                                "No bias / dark / flat frames found. Add a folder "
                                "with calibration frames (or a pre-built master).")
            return
        out_dir = Path(self.out_edit.text().strip() or self._default_out_dir())
        if not self._confirm_rerun(out_dir):
            return

        self.btn_run.setEnabled(False)
        self.btn_skip.setEnabled(False)
        self.btn_scan.setEnabled(False)
        scope = "all nights" if night == "__all__" else f"night {night}"
        self.log(f"Starting calibration ({scope})…")
        opts = self._current_options()
        self._worker = _CalibrationWorker(self._frames, night, out_dir, opts)
        self._worker.progress.connect(self._on_progress)
        self._worker.logline.connect(self.log)
        self._worker.finished_ok.connect(self._on_finished)
        self._worker.failed.connect(self._on_failed)
        self._worker.start()

    def _confirm_rerun(self, out_dir: Path) -> bool:
        """Warn before overwriting a previous run in the same output folder.

        Calibration is the slowest step in the chain (tens of minutes for a
        night), so silently redoing it — or silently overwriting a good run —
        is worth one question.
        """
        summary_path = out_dir / "calibration.json"
        if not summary_path.exists():
            return True
        try:
            prev = json.loads(summary_path.read_text(encoding="utf-8"))
            n_prev = int(prev.get("n_calibrated", 0))
            nights = ", ".join(sorted(prev.get("nights", {}))) or "?"
        except Exception:
            n_prev, nights = 0, "?"
        resp = QMessageBox.question(
            self, "Detector Calibration",
            f"{out_dir} already holds a calibration run "
            f"({n_prev} frame(s), night(s): {nights}).\n\n"
            "Run again and overwrite it?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
        )
        return resp == QMessageBox.Yes

    def _on_progress(self, done: int, total: int, message: str):
        if total > 0:
            self.progress_bar.setMaximum(total)
            self.progress_bar.setValue(done)
        self.progress_bar.setFormat(f"{message} (%p%)" if total else message)

    def _on_failed(self, message: str):
        self._reset_buttons()
        self.log(f"[ERROR] {message}")
        QMessageBox.critical(self, "Detector Calibration", message)

    def _on_finished(self, summary: Dict):
        self._reset_buttons()
        n = summary.get("n_calibrated", 0)
        miss = summary.get("n_missing_flat", 0)
        nn = summary.get("n_nights", 1)
        temp_miss = summary.get("n_temp_mismatch", 0)
        refused = summary.get("n_dark_refused", 0)
        nights = summary.get("nights", {})
        # Single night → point File Selection at that night's folder; multi-night
        # → the calibrated root (per-night subfolders inside).
        if len(nights) == 1:
            calibrated_dir = next(iter(nights.values())).get("calibrated_dir", "")
        else:
            calibrated_dir = summary.get("calibrated_root", "")
        self.log(f"Done — {n} frames across {nn} night(s) calibrated"
                 + (f" ({miss} without a matching flat)" if miss else "")
                 + (f" ({temp_miss} with a dark outside the temperature tolerance"
                    + (f", {refused} refused" if refused else "") + ")"
                    if temp_miss else "")
                 + f" → {calibrated_dir}")
        if self.project_state is not None:
            self.project_state.mark_calibration("done")
            self._refresh_main_window()
        resp = QMessageBox.question(
            self, "Detector Calibration",
            f"Calibrated {n} frames"
            + (f" ({miss} had no matching flat)" if miss else "")
            + (f"\n{temp_miss} frame(s) had no dark within the temperature "
               f"tolerance" + (f"; {refused} were calibrated without a dark."
                               if refused else ".") if temp_miss else "")
            + "\n\nUse the calibrated folder as the File-Selection input directory?",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.Yes,
        )
        if resp == QMessageBox.Yes and calibrated_dir:
            self._set_data_dir(calibrated_dir)

    def _on_skip(self):
        if self.project_state is not None:
            self.project_state.mark_calibration("skipped")
            self._refresh_main_window()
        self.log("Calibration skipped.")
        self.close()

    def _reset_buttons(self):
        self.btn_run.setEnabled(self.cmb_night.count() > 0)
        self.btn_skip.setEnabled(True)
        self.btn_scan.setEnabled(True)

    def _set_data_dir(self, calibrated_dir: str):
        P = getattr(self.params, "P", None)
        if P is None:
            return
        try:
            setattr(P, "data_dir", calibrated_dir)
            # Calibrated frames are written as ``pp_<original>``, so a prefix
            # left over from the raw files matches nothing and Step 1 comes up
            # with an empty file list. Keep it only if it still fits.
            prefix = str(getattr(P, "filename_prefix", "") or "")
            if prefix and not prefix.startswith("pp_"):
                setattr(P, "filename_prefix", "")
                self.log(f"Filename prefix '{prefix}' cleared "
                         f"(calibrated frames are named pp_*)")
            self.log(f"File-Selection input set to {calibrated_dir}")
        except Exception:
            pass

    def _refresh_main_window(self):
        mw = self.parent()
        if mw is not None and hasattr(mw, "update_step_buttons"):
            try:
                mw.update_step_buttons()
            except Exception:
                pass

    def closeEvent(self, event):
        for w in (self._scan_worker, self._worker):
            if w is not None and w.isRunning():
                w.requestInterruption()
                w.wait(3000)
        super().closeEvent(event)
