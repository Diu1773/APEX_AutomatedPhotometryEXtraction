"""
Step 1: File Selection Window
Popup window with Previous/Next navigation
"""

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QTableWidget, QTableWidgetItem, QGroupBox, QLineEdit,
    QFileDialog, QMessageBox, QHeaderView
)
from PyQt5.QtCore import Qt
from pathlib import Path

from apex.gui.workflow.step_window_base import StepWindowBase


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
        target_group = QGroupBox("SIMBAD Target")
        target_layout = QHBoxLayout(target_group)

        target_layout.addWidget(QLabel("Target Name:"))
        self.target_edit = QLineEdit()
        self.target_edit.setPlaceholderText("e.g., M31 or WASP-12")
        self.target_edit.setMaximumWidth(200)
        target_layout.addWidget(self.target_edit)

        btn_resolve = QPushButton("Resolve SIMBAD")
        btn_resolve.clicked.connect(self.resolve_target)
        target_layout.addWidget(btn_resolve)

        self.target_result = QLabel("(not resolved)")
        self.target_result.setStyleSheet("QLabel { font-weight: bold; color: #4CAF50; }")
        target_layout.addWidget(self.target_result)

        target_layout.addStretch()
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
        self.header_table.setColumnCount(7)
        self.header_table.setHorizontalHeaderLabels([
            "Use", "Filename", "DATE-OBS", "FILTER", "EXPTIME", "AIRMASS", "IMAGETYP"
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
        # Only check if files are loaded
        # Reference frame selection will be done after detection/ID matching in later steps
        has_files = len(self.file_manager.filenames) > 0

        return has_files

    def browse_directory(self):
        """Browse for data directory"""
        dir_path = QFileDialog.getExistingDirectory(
            self, "Select Data Directory",
            str(self.params.P.data_dir)
        )
        if dir_path:
            data_path = Path(dir_path)
            self.params.P.data_dir = data_path
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

        try:
            inst = self.main_window.instrument
            inst.targets_resolved = []
            inst.primary_target = None
            inst.primary_coord = None
            inst.resolve_targets([name])
            if inst.primary_coord is None:
                self.target_result.setText("(not resolved)")
                QMessageBox.warning(self, "SIMBAD", f"Target not found: {name}")
                return

            ra_deg = float(inst.primary_coord.ra.deg)
            dec_deg = float(inst.primary_coord.dec.deg)
            ra_hms = inst.primary_coord.ra.to_string(unit="hour", sep=":", precision=2)
            dec_dms = inst.primary_coord.dec.to_string(unit="deg", sep=":", precision=1, alwayssign=True)
            self.target_result.setText(f"{ra_hms}, {dec_dms}")

            # Save to params.P and TOML (primary source of truth)
            self.params.P.target_name = name
            self.params.P.target_ra_deg = ra_deg
            self.params.P.target_dec_deg = dec_deg
            if hasattr(self.params, 'save_toml'):
                self.params.save_toml()

        except Exception as e:
            QMessageBox.warning(self, "SIMBAD Error", str(e))

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

                self.header_table.setItem(i, 1, QTableWidgetItem(fname))
                self.header_table.setItem(i, 2, QTableWidgetItem(str(row["DATE-OBS"])))
                self.header_table.setItem(i, 3, QTableWidgetItem(str(row["FILTER"])))
                self.header_table.setItem(i, 4, QTableWidgetItem(str(row["EXPTIME"])))
                self.header_table.setItem(i, 5, QTableWidgetItem(str(row["AIRMASS"])))
                self.header_table.setItem(i, 6, QTableWidgetItem(str(row["IMAGETYP"])))

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
