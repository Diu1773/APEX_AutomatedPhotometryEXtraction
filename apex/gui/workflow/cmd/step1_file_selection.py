"""
Step 1: File Selection Window
Popup window with Previous/Next navigation
"""

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QTableWidget, QTableWidgetItem, QGroupBox, QLineEdit,
    QFileDialog, QMessageBox, QHeaderView, QWidget
)
from PyQt5.QtCore import Qt, QUrl
from PyQt5.QtGui import QDesktopServices
from pathlib import Path
from urllib.parse import quote_plus

from apex.gui.workflow.step_window_base import StepWindowBase
from apex.analysis.header_target import select_header_target
from apex.gui.workflow.target_resolver import TargetResolveWorker, target_failure_message


class FileSelectionWindow(StepWindowBase):
    """Step 1: File Selection and Header Reading"""

    def __init__(self, params, file_manager, project_state, main_window):
        """
        Initialize file selection window

        Args:
            params: Parameters object
            file_manager: FileManager instance
            project_state: ProjectState object
            main_window: Main window reference
        """
        self.file_manager = file_manager
        self.simbad_worker = None

        # Initialize base class
        super().__init__(
            step_index=0,
            step_name="File Selection",
            params=params,
            project_state=project_state,
            main_window=main_window
        )

        # Setup step-specific UI
        self.setup_step_ui()

        # Restore state AFTER UI is created
        self.restore_state()

    def setup_step_ui(self):
        """Setup step-specific UI components"""

        # === Directory Selection ===
        dir_group = QGroupBox("Data Directory")
        dir_layout = QHBoxLayout(dir_group)

        dir_layout.addWidget(QLabel("Directory:"))

        self.dir_edit = QLineEdit(str(self.params.P.data_dir))
        self.dir_edit.setReadOnly(True)
        dir_layout.addWidget(self.dir_edit)

        btn_browse = QPushButton("Browse...")
        btn_browse.clicked.connect(self.browse_directory)
        dir_layout.addWidget(btn_browse)

        dir_layout.addStretch()
        self.file_count_label = QLabel("Files: 0")
        dir_layout.addWidget(self.file_count_label)

        self.content_layout.addWidget(dir_group)

        # === SIMBAD Target ===
        target_group = QGroupBox("Target Coordinates")
        target_layout = QVBoxLayout(target_group)

        target_row = QHBoxLayout()
        target_row.addWidget(QLabel("Target Name:"))
        self.target_edit = QLineEdit()
        self.target_edit.setPlaceholderText("e.g., M31 or WASP-12")
        self.target_edit.setMaximumWidth(200)
        target_row.addWidget(self.target_edit)

        self.btn_resolve_simbad = QPushButton("Resolve SIMBAD")
        self.btn_resolve_simbad.clicked.connect(self.resolve_target)
        target_row.addWidget(self.btn_resolve_simbad)

        btn_open_simbad = QPushButton("Open SIMBAD")
        btn_open_simbad.setToolTip("Open the target name in the SIMBAD web page.")
        btn_open_simbad.clicked.connect(self.open_simbad_page)
        target_row.addWidget(btn_open_simbad)

        btn_manual = QPushButton("Manual RA/Dec")
        btn_manual.setToolTip("Enter target coordinates manually when SIMBAD cannot resolve the object.")
        btn_manual.clicked.connect(self._toggle_manual_radec)
        target_row.addWidget(btn_manual)

        btn_header_radec = QPushButton("Use Header RA/Dec")
        btn_header_radec.setToolTip(
            "선택한 사용 프레임의 FITS header RA/Dec를 target 좌표로 사용합니다. "
            "선택 행이 없으면 첫 번째 사용 프레임을 사용합니다."
        )
        btn_header_radec.clicked.connect(self._use_header_radec)
        target_row.addWidget(btn_header_radec)

        self.target_result = QLabel("(not resolved)")
        self.target_result.setStyleSheet("QLabel { font-weight: bold; color: #4CAF50; }")
        target_row.addWidget(self.target_result)

        target_row.addStretch()
        target_layout.addLayout(target_row)

        self._manual_radec_widget = QWidget()
        manual_row = QHBoxLayout(self._manual_radec_widget)
        manual_row.setContentsMargins(0, 0, 0, 0)
        manual_row.addWidget(QLabel("RA:"))
        self.manual_ra_edit = QLineEdit()
        self.manual_ra_edit.setPlaceholderText("HH:MM:SS or deg")
        self.manual_ra_edit.setMaximumWidth(140)
        manual_row.addWidget(self.manual_ra_edit)
        manual_row.addWidget(QLabel("Dec:"))
        self.manual_dec_edit = QLineEdit()
        self.manual_dec_edit.setPlaceholderText("±DD:MM:SS or deg")
        self.manual_dec_edit.setMaximumWidth(140)
        manual_row.addWidget(self.manual_dec_edit)
        btn_apply_manual = QPushButton("Apply")
        btn_apply_manual.clicked.connect(self._apply_manual_radec)
        manual_row.addWidget(btn_apply_manual)
        manual_row.addStretch()
        self._manual_radec_widget.setVisible(False)
        target_layout.addWidget(self._manual_radec_widget)
        self.content_layout.addWidget(target_group)

        # === Header Table ===
        table_group = QGroupBox("FITS Headers")
        table_layout = QVBoxLayout(table_group)

        # Select / deselect buttons
        btn_row = QHBoxLayout()
        btn_select_all = QPushButton("전체 선택")
        btn_select_all.clicked.connect(self._select_all_files)
        btn_row.addWidget(btn_select_all)
        btn_deselect_all = QPushButton("선택 해제")
        btn_deselect_all.clicked.connect(self._deselect_all_files)
        btn_row.addWidget(btn_deselect_all)
        btn_row.addStretch()
        table_layout.addLayout(btn_row)

        self.header_table = QTableWidget()
        self.header_table.setColumnCount(10)
        self.header_table.setHorizontalHeaderLabels([
            "Use", "Filename", "DATE-OBS", "FILTER", "EXPTIME", "AIRMASS",
            "OBJECT", "RA_DEG", "DEC_DEG", "IMAGETYP"
        ])
        self.header_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.header_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.header_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.header_table.itemChanged.connect(self._on_file_use_changed)
        self._file_table_loading = False

        table_layout.addWidget(self.header_table)

        self.content_layout.addWidget(table_group)

        # Don't auto-load files - user must click "Rescan Files" button
        # (Files will be loaded from restore_state() if previously scanned)

    def validate_step(self) -> bool:
        """Validate if step can be completed"""
        has_files = len(self.file_manager.filenames) > 0
        has_target = self._has_target_coordinates()

        return has_files and has_target

    def _has_target_coordinates(self) -> bool:
        try:
            ra = float(getattr(self.params.P, "target_ra_deg", float("nan")))
            dec = float(getattr(self.params.P, "target_dec_deg", float("nan")))
        except (TypeError, ValueError):
            return False
        return (0.0 <= ra < 360.0) and (-90.0 <= dec <= 90.0)

    @staticmethod
    def _format_header_cell(row, key: str) -> str:
        try:
            value = row.get(key, "")
        except Exception:
            value = ""
        if value is None:
            return ""
        try:
            if value != value:
                return ""
        except Exception:
            pass
        if key in {"RA_DEG", "DEC_DEG"}:
            try:
                return f"{float(value):.5f}"
            except (TypeError, ValueError):
                return str(value)
        return str(value)

    def open_simbad_page(self) -> None:
        name = self.target_edit.text().strip()
        url = "https://simbad.cds.unistra.fr/simbad/sim-id"
        if name:
            url = f"{url}?Ident={quote_plus(name)}"
        QDesktopServices.openUrl(QUrl(url))

    def _clear_target_coordinates(self) -> None:
        self.params.P.target_name = ""
        self.params.P.target_ra_deg = None
        self.params.P.target_dec_deg = None
        try:
            inst = self.main_window.instrument
            inst.targets_resolved = []
            inst.primary_target = None
            inst.primary_coord = None
        except Exception:
            pass
        self.target_edit.clear()
        self.target_result.setText("(not resolved)")

    def _clear_resolved_target_state(self, name: str | None = None, *, persist: bool = False) -> None:
        target_name = str(name if name is not None else self.target_edit.text()).strip()
        self.params.P.target_name = target_name
        self.params.P.target_ra_deg = None
        self.params.P.target_dec_deg = None
        try:
            inst = self.main_window.instrument
            inst.targets_resolved = []
            inst.primary_target = None
            inst.primary_coord = None
        except Exception:
            pass
        if persist and hasattr(self.params, 'save_toml'):
            self.params.save_toml()
        try:
            self._invalidate_workflow_from_here()
        except Exception:
            pass
        self.update_navigation_buttons()

    def _invalidate_workflow_from_here(self) -> None:
        for idx in list(self.project_state.state.get("completed_steps", [])):
            if int(idx) >= self.step_index:
                self.project_state.mark_step_incomplete(int(idx))
        if hasattr(self.main_window, "update_step_buttons"):
            self.main_window.update_step_buttons()

    def browse_directory(self):
        """Browse for data directory"""
        dir_path = QFileDialog.getExistingDirectory(
            self, "Select Data Directory",
            str(self.params.P.data_dir)
        )
        if dir_path:
            data_path = Path(dir_path)
            self.params.P.data_dir = data_path
            self.params.P.filename_prefix = ""
            self._clear_target_coordinates()
            self._invalidate_workflow_from_here()
            self.dir_edit.setText(dir_path)

            # Update result_dir to match new data_dir
            self.params.P.result_dir = data_path / "result"
            self.params.P.result_dir.mkdir(parents=True, exist_ok=True)
            self.params.P.cache_dir = self.params.P.result_dir / "cache"
            self.params.P.cache_dir.mkdir(parents=True, exist_ok=True)

            # Save to TOML
            if hasattr(self.params, 'save_toml'):
                self.params.save_toml()

            # Clear file manager state and auto-scan new directory
            self.file_manager.filenames = []
            self.file_manager.df_headers = None
            self.file_manager.ref_filename = None
            self.load_files()

    def resolve_target(self):
        """Resolve target coordinates via SIMBAD"""
        name = self.target_edit.text().strip()
        if not name:
            QMessageBox.warning(self, "Target Missing", "Please enter a target name.")
            return

        if self.simbad_worker is not None and self.simbad_worker.isRunning():
            QMessageBox.information(self, "SIMBAD", "SIMBAD lookup is already running.")
            return

        self.target_result.setText("(resolving...)")
        self.btn_resolve_simbad.setEnabled(False)
        self._clear_resolved_target_state(name, persist=True)
        self.simbad_worker = TargetResolveWorker(self.main_window.instrument, name, self)
        self.simbad_worker.log_message.connect(print)
        self.simbad_worker.result_ready.connect(self._on_simbad_result)
        self.simbad_worker.error_message.connect(self._on_simbad_error)
        self.simbad_worker.finished.connect(self._on_simbad_finished)
        self.simbad_worker.start()

    def _on_simbad_result(self, result: dict) -> None:
        if not result.get("ok"):
            name = str(result.get("name") or self.target_edit.text().strip())
            self._clear_resolved_target_state(name, persist=True)
            self.target_result.setText(f"(not resolved: {name})")
            QMessageBox.warning(self, "SIMBAD", target_failure_message(result))
            return

        name = str(result.get("name") or self.target_edit.text().strip())
        ra_deg = float(result["ra_deg"])
        dec_deg = float(result["dec_deg"])
        self.target_result.setText(f"{result['ra_hms']}, {result['dec_dms']}")
        self.params.P.target_name = name
        self.params.P.target_ra_deg = ra_deg
        self.params.P.target_dec_deg = dec_deg
        if hasattr(self.params, 'save_toml'):
            self.params.save_toml()
        self.update_navigation_buttons()

    def _on_simbad_error(self, message: str) -> None:
        name = self.target_edit.text().strip()
        self._clear_resolved_target_state(name, persist=True)
        self.target_result.setText(f"(not resolved: {name})" if name else "(not resolved)")
        QMessageBox.warning(self, "SIMBAD Error", message)
        self.update_navigation_buttons()

    def _on_simbad_finished(self) -> None:
        self.btn_resolve_simbad.setEnabled(True)
        if self.simbad_worker is not None:
            self.simbad_worker.deleteLater()
            self.simbad_worker = None

    def _toggle_manual_radec(self) -> None:
        self._manual_radec_widget.setVisible(not self._manual_radec_widget.isVisible())

    def _apply_manual_radec(self) -> None:
        from astropy.coordinates import SkyCoord
        import astropy.units as u

        ra_text = self.manual_ra_edit.text().strip()
        dec_text = self.manual_dec_edit.text().strip()
        name = self.target_edit.text().strip() or "manual"
        if not ra_text or not dec_text:
            QMessageBox.warning(self, "Manual Coordinates", "Enter both RA and Dec.")
            return
        try:
            try:
                coord = SkyCoord(ra_text, dec_text, unit=(u.hourangle, u.deg))
            except Exception:
                coord = SkyCoord(float(ra_text), float(dec_text), unit="deg")
        except Exception as exc:
            QMessageBox.warning(
                self,
                "Manual Coordinates",
                "Failed to parse coordinates.\n\n"
                "RA examples: 16:41:41.63 or 250.4234\n"
                "Dec examples: +36:27:40 or 36.461\n\n"
                f"{exc}",
            )
            return

        ra_deg = float(coord.ra.deg)
        dec_deg = float(coord.dec.deg)
        ra_hms = coord.ra.to_string(unit="hour", sep=":", precision=2)
        dec_dms = coord.dec.to_string(unit="deg", sep=":", precision=1, alwayssign=True)

        try:
            inst = self.main_window.instrument
            inst.targets_resolved = [dict(
                name=name,
                ra_deg=ra_deg,
                dec_deg=dec_deg,
                ra_str=ra_hms,
                dec_str=dec_dms,
                vmag=float("nan"),
                otype="",
                simbad_query="manual",
            )]
            inst.primary_target = name
            inst.primary_coord = coord
        except Exception:
            pass

        self.target_result.setText(f"{ra_hms}, {dec_dms}")
        self.params.P.target_name = name
        self.params.P.target_ra_deg = ra_deg
        self.params.P.target_dec_deg = dec_deg
        if hasattr(self.params, 'save_toml'):
            self.params.save_toml()
        self._manual_radec_widget.setVisible(False)
        self.update_navigation_buttons()

    def _use_header_radec(self) -> None:
        headers_df = getattr(self.file_manager, "df_headers", None)
        excluded = set(getattr(self.file_manager, "excluded_files", set()) or set())
        eligible = [
            str(name)
            for name in getattr(self.file_manager, "filenames", [])
            if str(name) not in excluded
        ]
        preferred = []
        for row in sorted({index.row() for index in self.header_table.selectedIndexes()}):
            item = self.header_table.item(row, 1)
            if item is not None:
                preferred.append(item.text().strip())

        header_target = select_header_target(headers_df, eligible, preferred)
        if header_target is None:
            QMessageBox.warning(
                self,
                "Header RA/Dec",
                "사용할 수 있는 FITS header RA/Dec를 찾지 못했습니다.\n\n"
                "먼저 FITS 파일을 불러오고 RA_DEG/DEC_DEG 열을 확인하세요.",
            )
            return

        from astropy.coordinates import SkyCoord
        import astropy.units as u

        name = self.target_edit.text().strip() or header_target["object_name"]
        coord = SkyCoord(
            header_target["ra_deg"],
            header_target["dec_deg"],
            unit="deg",
            frame="icrs",
        )
        ra_hms = coord.ra.to_string(unit="hour", sep=":", precision=2)
        dec_dms = coord.dec.to_string(
            unit="deg",
            sep=":",
            precision=1,
            alwayssign=True,
        )

        try:
            inst = self.main_window.instrument
            inst.targets_resolved = [dict(
                name=name,
                ra_deg=float(coord.ra.deg),
                dec_deg=float(coord.dec.deg),
                ra_str=ra_hms,
                dec_str=dec_dms,
                vmag=float("nan"),
                otype="",
                simbad_query=f"fits_header:{header_target['filename']}",
            )]
            inst.primary_target = name
            inst.primary_coord = coord
        except Exception:
            pass

        self.target_edit.setText(name)
        self.target_result.setText(f"{ra_hms}, {dec_dms} (header)")
        self.target_result.setToolTip(f"Coordinates from {header_target['filename']}")
        self.params.P.target_name = name
        self.params.P.target_ra_deg = float(coord.ra.deg)
        self.params.P.target_dec_deg = float(coord.dec.deg)
        if hasattr(self.params, "save_toml"):
            self.params.save_toml()
        self.update_navigation_buttons()

    def load_files(self):
        """Load files and headers"""
        try:
            filenames = self.file_manager.scan_files()
            self.file_count_label.setText(f"Files: {len(filenames)}")

            df_headers = self.file_manager.read_headers()

            excluded = getattr(self.file_manager, "excluded_files", set())

            self._file_table_loading = True
            self.header_table.blockSignals(True)
            self.header_table.setRowCount(len(df_headers))

            for i, row in df_headers.iterrows():
                fname = str(row["Filename"])

                use_item = QTableWidgetItem()
                use_item.setFlags(use_item.flags() | Qt.ItemIsUserCheckable)
                use_item.setCheckState(Qt.Unchecked if fname in excluded else Qt.Checked)
                use_item.setData(Qt.UserRole, fname)
                self.header_table.setItem(i, 0, use_item)

                labels = [
                    "Filename", "DATE-OBS", "FILTER", "EXPTIME", "AIRMASS",
                    "OBJECT", "RA_DEG", "DEC_DEG", "IMAGETYP",
                ]
                for col, label in enumerate(labels, start=1):
                    self.header_table.setItem(
                        i,
                        col,
                        QTableWidgetItem(self._format_header_cell(row, label)),
                    )

            self.header_table.blockSignals(False)
            self._file_table_loading = False

            self.update_navigation_buttons()

        except Exception as e:
            self._file_table_loading = False
            self.header_table.blockSignals(False)
            QMessageBox.critical(self, "Error", f"Failed to load files:\n{str(e)}")

    def _on_file_use_changed(self, item: QTableWidgetItem):
        if self._file_table_loading or item.column() != 0:
            return
        self._sync_excluded_files()

    def _sync_excluded_files(self):
        excluded = set()
        for row in range(self.header_table.rowCount()):
            use_item = self.header_table.item(row, 0)
            if use_item and use_item.checkState() == Qt.Unchecked:
                fname = use_item.data(Qt.UserRole)
                if fname:
                    excluded.add(fname)
        self.file_manager.excluded_files = excluded

    def _select_all_files(self):
        self._file_table_loading = True
        self.header_table.blockSignals(True)
        for row in range(self.header_table.rowCount()):
            item = self.header_table.item(row, 0)
            if item:
                item.setCheckState(Qt.Checked)
        self.header_table.blockSignals(False)
        self._file_table_loading = False
        self._sync_excluded_files()

    def _deselect_all_files(self):
        self._file_table_loading = True
        self.header_table.blockSignals(True)
        for row in range(self.header_table.rowCount()):
            item = self.header_table.item(row, 0)
            if item:
                item.setCheckState(Qt.Unchecked)
        self.header_table.blockSignals(False)
        self._file_table_loading = False
        self._sync_excluded_files()

    def save_state(self):
        """Save step state to project"""
        self._sync_excluded_files()
        state_data = {
            "data_dir": str(self.params.P.data_dir),
            "result_dir": str(self.params.P.result_dir),
            "file_count": len(self.file_manager.filenames),
            "excluded_files": sorted(self.file_manager.excluded_files),
        }

        self.project_state.store_step_data("file_selection", state_data)

    def restore_state(self):
        """Restore step state from project"""
        state_data = self.project_state.get_step_data("file_selection")

        if state_data:
            if "data_dir" in state_data:
                self.params.P.data_dir = Path(state_data["data_dir"])
                self.dir_edit.setText(str(state_data["data_dir"]))

            # Restore exclusions before loading so checkboxes reflect saved state
            if "excluded_files" in state_data:
                self.file_manager.excluded_files = set(state_data["excluded_files"])

            try:
                self.load_files()
            except:
                pass

        # Load target from params.P (which comes from TOML)
        name = getattr(self.params.P, "target_name", None)
        ra = getattr(self.params.P, "target_ra_deg", None)
        dec = getattr(self.params.P, "target_dec_deg", None)
        if name:
            self.target_edit.setText(str(name))
        if ra is not None and dec is not None:
            try:
                self.target_result.setText(f"{float(ra):.6f}, {float(dec):.6f}")
            except (ValueError, TypeError):
                pass
        self.update_navigation_buttons()
