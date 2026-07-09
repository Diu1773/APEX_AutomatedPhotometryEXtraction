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

import numpy as np
from astropy.io import fits
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtWidgets import (
    QCheckBox, QComboBox, QDialogButtonBox, QDoubleSpinBox, QFileDialog,
    QFormLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit, QMessageBox,
    QProgressBar, QPushButton, QScrollArea, QSpinBox, QTreeWidget,
    QTreeWidgetItem, QVBoxLayout, QWidget,
)

from apex.analysis import calibration as cal
from apex.analysis import calibration_scan as scan
from apex.analysis.calibration import CalibrationOptions
from apex.analysis.calibration_scan import FrameInfo
from apex.gui.layout_rules import FittedDialog
from apex.gui.theme import ICON, Tokens, style_button
from apex.gui.tools.tool_window_base import ToolWindowBase
from apex.utils.step_paths import step0_calibration_dir


# ---------------------------------------------------------------------------
# Workers
# ---------------------------------------------------------------------------

class _ScanWorker(QThread):
    progress = pyqtSignal(int, int, str)
    done = pyqtSignal(list)                       # List[FrameInfo]

    def __init__(self, roots: List[str]):
        super().__init__()
        self.roots = roots

    def run(self):
        frames: List[FrameInfo] = []
        for root in self.roots:
            frames.extend(scan.scan_folder(
                root,
                progress=lambda i, t, m: self.progress.emit(i, t, m),
                stop=self.isInterruptionRequested,
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
        opts = self.opts
        masters_dir = self.out_dir / "masters"
        masters_dir.mkdir(parents=True, exist_ok=True)

        nights = (scan.nights(self.frames) if self.night == "__all__"
                  else [self.night])
        nights = [n for n in nights
                  if scan.group_for_night(self.frames, n)["light"]]
        if not nights:
            raise ValueError("No light frames found to calibrate")

        # Master cache keyed by the exact source files, so a master built from a
        # shared bias/dark pool (or a single pre-built master frame) is built
        # ONCE and reused across every night; per-night flats build per night.
        cache: Dict = {}
        written: set = set()
        summary: Dict = {"nights": {}, "masters": [], "options": opts.__dict__.copy()}

        def _write_master(fname, arr, prov, label):
            if fname not in written:
                fits.writeto(masters_dir / fname, arr.astype(np.float32), overwrite=True)
                written.add(fname)
                summary["masters"].append({"file": fname, **prov})
                self.logline.emit(f"{label}: {prov['n_frames']} frame(s) → {fname}")

        def _bias(paths):
            if not paths:
                return None
            key = ("bias", frozenset(paths))
            if key not in cache:
                mb, prov = cal.build_master_bias(paths, opts)
                cache[key] = mb
                _write_master("master_bias.fits", mb, prov, "master bias")
            return cache[key]

        def _dark(paths, dk, mbias, night_tag):
            if not paths:
                return None, 1.0
            key = ("dark", frozenset(paths))
            if key not in cache:
                md, dexp, prov = cal.build_master_dark(paths, opts, master_bias=mbias)
                cache[key] = (md, dexp)
                exp, tb = dk
                tag = night_tag or "global"   # keep per-night darks distinct on disk
                _write_master(f"master_dark_{exp:g}s_{tb}C_{tag}.fits", md, prov,
                              f"master dark {exp:g}s/{tb}°C [{tag}]")
            return cache[key]

        def _flat(fis, filt, mbias):
            if not fis:
                return None
            paths = [f.path for f in fis]
            key = ("flat", frozenset(paths))
            if key not in cache:
                mf, prov = cal.build_master_flat(paths, opts, master_bias=mbias)
                cache[key] = mf
                tag = fis[0].night or "global"   # the flat frames' own night
                _write_master(f"master_flat_{filt}_{tag}.fits", mf, prov,
                              f"master flat {filt} [{tag}]")
            return cache[key]

        total = sum(len(scan.group_for_night(self.frames, n)["light"]) for n in nights)
        done = 0
        no_flat = 0
        self.logline.emit(f"Calibrating {total} light frames across {len(nights)} night(s)…")

        for night in nights:
            g = scan.group_for_night(self.frames, night)
            calibrated_dir = self.out_dir / "calibrated" / (night or "undated")
            calibrated_dir.mkdir(parents=True, exist_ok=True)
            mbias = _bias([f.path for f in g["bias"]])
            night_recs: List[Dict] = []
            for l in g["light"]:
                if self.isInterruptionRequested():
                    self.logline.emit("Stopped by user.")
                    break
                self.progress.emit(done, total, f"[{night}] {l.name}")
                dk = scan.match_dark(g["dark"], l.exp, l.temp)
                md, dexp = (_dark([f.path for f in g["dark"][dk]], dk, mbias,
                                  g["dark"][dk][0].night) if dk else (None, 1.0))
                ff = scan.match_flat(g["flat"], l.filt)
                mf = _flat(g["flat"][ff], ff, mbias) if ff else None
                if mf is None:
                    no_flat += 1
                    self.logline.emit(f"  [warn] {l.name}: no flat for '{l.filt}' — flat skipped")
                data, header, qc = cal.calibrate_light_file(
                    l.path, opts, master_bias=mbias, master_dark=md,
                    dark_exp=dexp, master_flat=mf)
                fits.writeto(calibrated_dir / f"pp_{l.name}", data, header=header, overwrite=True)
                done += 1
                night_recs.append({
                    "input": l.name, "output": f"pp_{l.name}", "filter": l.filt,
                    "dark": f"{dk[0]:g}s/{dk[1]}C" if dk else None, "flat": ff,
                    "median": qc.get("median"), "neg_pct": qc.get("neg_pct"),
                })
            summary["nights"][night] = {
                "calibrated_dir": str(calibrated_dir),
                "n_calibrated": len(night_recs), "frames": night_recs,
            }
        self.progress.emit(total, total, "Done")

        summary["n_calibrated"] = done
        summary["n_missing_flat"] = no_flat
        summary["n_nights"] = len(nights)
        summary["calibrated_root"] = str(self.out_dir / "calibrated")
        (self.out_dir / "calibration.json").write_text(
            json.dumps(summary, indent=2, default=float), encoding="utf-8")
        self.finished_ok.emit(summary)


# ---------------------------------------------------------------------------
# Window
# ---------------------------------------------------------------------------

class DetectorCalibrationWindow(ToolWindowBase):
    """Optional off-chain Step 0 window (detector calibration)."""

    _TYPE_ORDER = ("bias", "dark", "flat", "light")

    def __init__(self, params, project_state, parent=None):
        super().__init__("Step 0: Detector Calibration", params=params,
                         project_state=project_state, parent=parent,
                         min_size=(560, 640))
        self._frames: List[FrameInfo] = []
        self._roots: List[str] = []
        self._scan_worker: Optional[_ScanWorker] = None
        self._worker: Optional[_CalibrationWorker] = None
        # All calibration parameters live in one mutable holder, edited via the
        # Parameters dialog (AIPPI-style single settings panel).
        self._settings: Dict = asdict(CalibrationOptions())
        # Seed cosmetic's noise model from the CONFIGURED (measured) gain/read
        # noise, never the FITS EGAIN header (unreliable on this camera).
        P = getattr(params, "P", None)
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
        v.addWidget(self.tree, 1)

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

    def _default_out_dir(self) -> Path:
        result_dir = getattr(getattr(self.params, "P", None), "result_dir", "") or ""
        base = Path(result_dir) if result_dir else Path.cwd()
        return step0_calibration_dir(base)

    def _current_options(self) -> CalibrationOptions:
        return CalibrationOptions(**self._settings)

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
        f.addRow(w_dscale); f.addRow(w_dopt)

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
                "flat_min": w_fmin.value(),
                "pedestal_mode": w_pmode.currentText(),
                "pedestal_value": w_pval.value(), "pedestal_max": w_pmax.value(),
                "overscan_enable": w_oen.isChecked(), "overscan_edge": w_oedge.currentText(),
                "overscan_width": w_owid.value(), "overscan_trim": w_otrim.isChecked(),
                "cosmetic_enable": w_cc.isChecked(), "cr_sigclip": w_crs.value(),
                "cr_objlim": w_cro.value(), "hot_sigma": w_hot.value(),
            })
            self.log("Parameters updated.")

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
        self._scan_worker = _ScanWorker([folder])
        self._scan_worker.progress.connect(self._on_scan_progress)
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
        self.progress_bar.reset()
        by_type = {t: sum(1 for f in self._frames if f.ftype == t) for t in self._TYPE_ORDER}
        self.scan_label.setText(
            f"{len(self._frames)} frames  ·  "
            + "  ".join(f"{t} {by_type[t]}" for t in self._TYPE_ORDER)
        )
        self.log(f"Scanned {len(frames)} new frames ({len(self._frames)} total)")
        self._populate_nights()
        self._rebuild_tree()

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
            for t in self._TYPE_ORDER:
                items = [f for f in self._frames if f.ftype == t and f.night == n]
                if not items:
                    continue
                type_node = QTreeWidgetItem([t.capitalize(), str(len(items)), ""])
                node.addChild(type_node)
                # Always pass the resolved pool so each night's light rows show
                # which dark/flat (incl. shared/global) will be used.
                self._add_subgroups(type_node, t, items, g)
            node.setExpanded(expand_all or n == night)

    def _add_subgroups(self, parent, ftype, items, group):
        if ftype == "dark":
            keys: Dict = {}
            for f in items:
                keys.setdefault((round(f.exp, 3), f.temp_bucket), 0)
                keys[(round(f.exp, 3), f.temp_bucket)] += 1
            for (exp, tb), c in sorted(keys.items()):
                parent.addChild(QTreeWidgetItem([f"{exp:g}s · {tb}°C", str(c), ""]))
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
        dk = scan.match_dark(group["dark"], sample.exp, sample.temp)
        ff = scan.match_flat(group["flat"], sample.filt)
        d = f"dark {dk[0]:g}s/{dk[1]}°C" if dk else "no dark"
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
        nights = summary.get("nights", {})
        # Single night → point File Selection at that night's folder; multi-night
        # → the calibrated root (per-night subfolders inside).
        if len(nights) == 1:
            calibrated_dir = next(iter(nights.values())).get("calibrated_dir", "")
        else:
            calibrated_dir = summary.get("calibrated_root", "")
        self.log(f"Done — {n} frames across {nn} night(s) calibrated"
                 + (f" ({miss} without a matching flat)" if miss else "")
                 + f" → {calibrated_dir}")
        if self.project_state is not None:
            self.project_state.mark_calibration("done")
            self._refresh_main_window()
        resp = QMessageBox.question(
            self, "Detector Calibration",
            f"Calibrated {n} frames"
            + (f" ({miss} had no matching flat)" if miss else "")
            + ".\n\nUse the calibrated folder as the File-Selection input directory?",
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
        if P is not None:
            try:
                setattr(P, "data_dir", calibrated_dir)
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
