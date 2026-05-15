"""
Shared Step 1 file intake window.

This module owns the mode-independent workflow: choose input files, scan FITS
headers, resolve target coordinates, choose usable frames, and persist the common
Step 1 state. CMD and LC wrappers extend only the parts that truly differ.
"""

from __future__ import annotations

import json
from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QAbstractItemView,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QToolButton,
    QVBoxLayout,
)

from apex.gui.workflow.step_window_base import StepWindowBase
from apex.utils.step_paths import step1_dir

try:  # Python 3.11+
    import tomllib  # type: ignore
except Exception:  # Python 3.10
    import tomli as tomllib  # type: ignore

try:
    import tomli_w  # type: ignore
except Exception:
    tomli_w = None


class CommonFileSelectionWindow(StepWindowBase):
    """Mode-independent Step 1 window for FITS file intake."""

    directory_label = "Directory:"
    header_labels = ["Filename", "DATE-OBS", "FILTER", "EXPTIME", "AIRMASS", "IMAGETYP"]

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.root_dir = Path(params.P.data_dir)
        self.last_data_dir = Path(params.P.data_dir)
        self._all_filenames: list[str] = []
        self._all_file_path_map: dict[str, Path] = {}
        self._all_headers_df = None
        self._selected_frame_filenames: set[str] | None = None
        self._header_table_loading = False
        self.init_mode_state()

        super().__init__(
            step_index=0,
            step_name="File Selection",
            params=params,
            project_state=project_state,
            main_window=main_window,
        )

        self.setup_step_ui()
        self.restore_state()

    # ------------------------------------------------------------------
    # Hooks for mode-specific behavior
    # ------------------------------------------------------------------

    def init_mode_state(self) -> None:
        """Initialize subclass-owned state before Qt widgets are built."""

    def setup_directory_extension(self, dir_layout: QVBoxLayout) -> None:
        """Add mode-specific directory controls below the root directory row."""

    def configure_header_table(self) -> None:
        self.header_table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeToContents)
        self.header_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.Stretch)
        self.header_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.header_table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.header_table.setEditTriggers(QTableWidget.NoEditTriggers)

    def prepare_scan_context(self) -> None:
        """Prepare file_manager state immediately before scan_files()."""
        if hasattr(self.file_manager, "clear_multi_night_dirs"):
            self.file_manager.clear_multi_night_dirs()

    def after_headers_loaded(self, df_headers) -> None:
        self.populate_header_table(df_headers)

    def extra_state(self) -> dict:
        return {"multi_night": False, "night_dirs": []}

    def restore_extra_state_before_load(self, state_data: dict) -> None:
        self.prepare_scan_context()

    def after_common_state_saved(self) -> None:
        """Persist mode-specific sidecar files after project_state is updated."""

    def after_target_resolved(self) -> None:
        """Mode-specific target side effects such as LC run manifest updates."""

    def after_frame_selection_applied(self, selected_filenames: set[str]) -> None:
        """Mode-specific side effects after frame Use checkboxes change."""

    def reset_mode_state_for_new_directory(self) -> None:
        if hasattr(self.file_manager, "clear_multi_night_dirs"):
            self.file_manager.clear_multi_night_dirs()

    def _update_result_workspace(self) -> None:
        data_dir = Path(self.params.P.data_dir)
        self.params.P.result_dir = data_dir / "result"
        self.params.P.cache_dir = self.params.P.result_dir / "cache"

    def _make_help_button(self, title: str, text: str) -> QToolButton:
        btn = QToolButton()
        btn.setText("?")
        btn.setAutoRaise(True)
        btn.setFixedSize(22, 22)
        btn.setToolTip("도움말 열기")
        btn.clicked.connect(lambda: QMessageBox.information(self, title, text))
        return btn

    # ------------------------------------------------------------------
    # UI setup
    # ------------------------------------------------------------------

    def setup_step_ui(self) -> None:
        dir_group = QGroupBox("Data Directory")
        dir_layout = QVBoxLayout(dir_group)

        root_row = QHBoxLayout()
        root_row.addWidget(QLabel(self.directory_label))
        self.dir_edit = QLineEdit(str(self.root_dir))
        self.dir_edit.setReadOnly(True)
        root_row.addWidget(self.dir_edit)

        btn_browse = QPushButton("Browse...")
        btn_browse.clicked.connect(self.browse_directory)
        root_row.addWidget(btn_browse)
        root_row.addWidget(self._make_help_button(
            "Data Directory",
            "FITS 파일이 있는 입력 폴더를 선택합니다.\n"
            "선택 후 Rescan Files를 눌러 파일 목록과 FITS 헤더를 갱신합니다.",
        ))
        dir_layout.addLayout(root_row)

        self.setup_directory_extension(dir_layout)
        self.content_layout.addWidget(dir_group)

        filter_group = QGroupBox("File Filter")
        filter_layout = QHBoxLayout(filter_group)
        filter_layout.addWidget(QLabel("Filename Prefix:"))

        self.prefix_edit = QLineEdit(self.params.P.filename_prefix)
        self.prefix_edit.setMaximumWidth(150)
        filter_layout.addWidget(self.prefix_edit)
        filter_layout.addWidget(self._make_help_button(
            "Filename Prefix",
            "읽을 FITS 파일 이름의 시작 문자열입니다.\n"
            "예: pp_0001.fits, pp_0002.fits만 읽으려면 pp_ 를 입력합니다.\n"
            "폴더 안의 모든 FITS를 읽으려면 비워 둡니다.",
        ))

        btn_rescan = QPushButton("Rescan Files")
        btn_rescan.setToolTip("현재 폴더와 Prefix 설정으로 FITS 파일과 헤더를 다시 읽습니다.")
        btn_rescan.clicked.connect(self.rescan_files)
        filter_layout.addWidget(btn_rescan)

        filter_layout.addStretch()
        self.file_count_label = QLabel("Files: 0")
        filter_layout.addWidget(self.file_count_label)
        self.content_layout.addWidget(filter_group)

        target_group = QGroupBox("SIMBAD Target")
        target_layout = QHBoxLayout(target_group)
        target_layout.addWidget(QLabel("Target Name:"))

        self.target_edit = QLineEdit()
        self.target_edit.setPlaceholderText("e.g., M31 or WASP-12")
        self.target_edit.setMaximumWidth(200)
        target_layout.addWidget(self.target_edit)

        btn_resolve = QPushButton("Resolve SIMBAD")
        btn_resolve.setToolTip("타겟 이름을 SIMBAD에서 찾아 RA/Dec 좌표를 저장합니다.")
        btn_resolve.clicked.connect(self.resolve_target)
        target_layout.addWidget(btn_resolve)
        target_layout.addWidget(self._make_help_button(
            "SIMBAD Target",
            "타겟 이름을 입력하고 Resolve SIMBAD를 누르면 Step 9 등에서 사용할 좌표를 저장합니다.\n"
            "인터넷 또는 SIMBAD 응답 문제가 있으면 나중 단계에서 수동 선택으로 진행할 수 있습니다.",
        ))

        self.target_result = QLabel("(not resolved)")
        self.target_result.setStyleSheet("QLabel { font-weight: bold; color: #4CAF50; }")
        target_layout.addWidget(self.target_result)
        target_layout.addStretch()
        self.content_layout.addWidget(target_group)

        table_group = QGroupBox("FITS Headers")
        table_layout = QVBoxLayout(table_group)
        frame_select_row = QHBoxLayout()
        frame_select_row.addWidget(QLabel("Frame Use:"))

        btn_use_all = QPushButton("전체 사용")
        btn_use_all.setToolTip("스캔된 모든 프레임을 이후 단계 입력으로 사용합니다.")
        btn_use_all.clicked.connect(self.select_all_frames)
        frame_select_row.addWidget(btn_use_all)

        btn_use_rows = QPushButton("선택 행만 사용")
        btn_use_rows.setToolTip("표에서 선택한 행만 Use 체크하고 나머지는 제외합니다. Ctrl/Shift로 여러 행을 선택할 수 있습니다.")
        btn_use_rows.clicked.connect(self.use_selected_table_rows_only)
        frame_select_row.addWidget(btn_use_rows)

        btn_clear_frames = QPushButton("전체 해제")
        btn_clear_frames.setToolTip("모든 프레임의 Use 체크를 끕니다.")
        btn_clear_frames.clicked.connect(self.clear_all_frames)
        frame_select_row.addWidget(btn_clear_frames)
        frame_select_row.addWidget(self._make_help_button(
            "Frame Use",
            "Use가 켜진 프레임만 이후 단계에서 처리합니다.\n"
            "특정 프레임을 빼고 싶으면 체크를 끄세요.\n"
            "딱 몇 장만 쓰려면 표에서 행을 Ctrl/Shift로 고른 뒤 선택 행만 사용을 누르세요.",
        ))
        frame_select_row.addStretch()
        table_layout.addLayout(frame_select_row)

        self.header_table = QTableWidget()
        self.header_table.setColumnCount(len(self.header_labels) + 1)
        self.header_table.setHorizontalHeaderLabels(["Use"] + self.header_labels)
        self.configure_header_table()
        self.header_table.itemChanged.connect(self._on_header_item_changed)
        table_layout.addWidget(self.header_table)
        self.content_layout.addWidget(table_group)

    # ------------------------------------------------------------------
    # Directory / scanning
    # ------------------------------------------------------------------

    def validate_step(self) -> bool:
        return len(self.file_manager.filenames) > 0

    def browse_directory(self) -> None:
        dir_path = QFileDialog.getExistingDirectory(
            self,
            "Select Data Directory",
            str(self.root_dir),
        )
        if not dir_path:
            return

        self.root_dir = Path(dir_path)
        self.reset_mode_state_for_new_directory()
        self._set_data_dir(self.root_dir)
        self.dir_edit.setText(dir_path)
        self._clear_loaded_file_state()

    def _clear_loaded_file_state(self) -> None:
        self.file_manager.filenames = []
        self.file_manager.df_headers = None
        self.file_manager.path_map = {}
        self.file_manager.ref_filename = None
        self._all_filenames = []
        self._all_file_path_map = {}
        self._all_headers_df = None
        self._selected_frame_filenames = None
        self.file_count_label.setText("Files: 0")
        self.header_table.setRowCount(0)

    def rescan_files(self) -> None:
        self.params.P.filename_prefix = self.prefix_edit.text()
        self._persist_param_file(io_updates={"filename_prefix": self.params.P.filename_prefix})

        try:
            self.prepare_scan_context()
            self._update_result_workspace()
            self._ensure_result_dirs()
            self._persist_param_file(
                io_updates={
                    "data_dir": str(self.params.P.data_dir),
                    "result_dir": str(self.params.P.result_dir),
                    "cache_dir": str(self.params.P.cache_dir.name),
                }
            )
            self.load_files()
            self.update_navigation_buttons()
        except Exception as e:
            QMessageBox.warning(
                self,
                "Scan Error",
                f"Failed to scan files:\n{str(e)}",
            )

    def load_files(self) -> None:
        try:
            filenames = self.file_manager.scan_files()
            df_headers = self.file_manager.read_headers()
            self._all_filenames = list(filenames)
            self._all_file_path_map = dict(getattr(self.file_manager, "path_map", {}) or {})
            self._all_headers_df = df_headers.copy()
            self.after_headers_loaded(df_headers)
            self._apply_frame_selection_from_table(write_sidecars=True)
            self.update_navigation_buttons()
        except Exception as e:
            QMessageBox.critical(
                self,
                "Error",
                f"Failed to load files:\n{str(e)}",
            )

    def populate_header_table(self, df_headers) -> None:
        self._header_table_loading = True
        self.header_table.blockSignals(True)
        self.header_table.setRowCount(len(df_headers))
        for i, row in df_headers.iterrows():
            filename = str(row["Filename"])
            self.header_table.setItem(i, 0, self._make_frame_use_item(filename))
            self.header_table.setItem(i, 1, QTableWidgetItem(filename))
            self.header_table.setItem(i, 2, QTableWidgetItem(str(row["DATE-OBS"])))
            self.header_table.setItem(i, 3, QTableWidgetItem(str(row["FILTER"])))
            self.header_table.setItem(i, 4, QTableWidgetItem(str(row["EXPTIME"])))
            self.header_table.setItem(i, 5, QTableWidgetItem(str(row["AIRMASS"])))
            self.header_table.setItem(i, 6, QTableWidgetItem(str(row["IMAGETYP"])))
        self.header_table.blockSignals(False)
        self._header_table_loading = False

    def _make_frame_use_item(self, filename: str) -> QTableWidgetItem:
        item = QTableWidgetItem()
        item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
        item.setData(Qt.UserRole, filename)
        selected = self._selected_frame_filenames
        item.setCheckState(Qt.Checked if selected is None or filename in selected else Qt.Unchecked)
        item.setToolTip("체크된 프레임만 이후 단계에서 처리합니다.")
        return item

    def _selected_filenames_from_table(self) -> set[str]:
        selected = set()
        for row in range(self.header_table.rowCount()):
            item = self.header_table.item(row, 0)
            if item is None or item.checkState() != Qt.Checked:
                continue
            filename = item.data(Qt.UserRole)
            if filename:
                selected.add(str(filename))
        return selected

    def _ordered_selected_filenames(self, selected: set[str]) -> list[str]:
        if self._all_filenames:
            return [fn for fn in self._all_filenames if fn in selected]
        return sorted(selected)

    def _filtered_headers(self, selected_order: list[str]):
        if self._all_headers_df is None:
            return None
        df = self._all_headers_df
        if "Filename" not in df.columns:
            return df.copy()
        order = {fn: i for i, fn in enumerate(selected_order)}
        out = df[df["Filename"].astype(str).isin(set(selected_order))].copy()
        if not out.empty:
            out["_apex_order"] = out["Filename"].astype(str).map(order)
            out = out.sort_values("_apex_order").drop(columns=["_apex_order"])
        return out

    def _apply_frame_selection_from_table(self, *, write_sidecars: bool = False) -> None:
        self._apply_frame_selection(
            self._selected_filenames_from_table(),
            update_table=False,
            write_sidecars=write_sidecars,
        )

    def _apply_frame_selection(
        self,
        selected: set[str],
        *,
        update_table: bool = True,
        write_sidecars: bool = True,
    ) -> None:
        selected_order = self._ordered_selected_filenames(selected)
        selected_set = set(selected_order)
        self._selected_frame_filenames = selected_set

        self.file_manager.filenames = selected_order
        self.file_manager.path_map = {
            fn: self._all_file_path_map[fn]
            for fn in selected_order
            if fn in self._all_file_path_map
        }
        self.params.P.file_path_map = {
            fn: str(path) for fn, path in self.file_manager.path_map.items()
        }
        filtered_headers = self._filtered_headers(selected_order)
        self.file_manager.df_headers = filtered_headers

        if update_table:
            self._set_frame_checks(selected_set)

        self.after_frame_selection_applied(selected_set)
        if write_sidecars:
            self._write_common_step1_sidecars()
        self._update_file_count_label()
        self.update_navigation_buttons()

    def _set_frame_checks(self, selected: set[str]) -> None:
        self._header_table_loading = True
        self.header_table.blockSignals(True)
        for row in range(self.header_table.rowCount()):
            item = self.header_table.item(row, 0)
            if item is None:
                continue
            filename = item.data(Qt.UserRole)
            item.setCheckState(Qt.Checked if str(filename) in selected else Qt.Unchecked)
        self.header_table.blockSignals(False)
        self._header_table_loading = False

    def _on_header_item_changed(self, item: QTableWidgetItem) -> None:
        if self._header_table_loading or item.column() != 0:
            return
        self._apply_frame_selection_from_table(write_sidecars=True)

    def _update_file_count_label(self) -> None:
        selected = len(self.file_manager.filenames)
        total = len(self._all_filenames) if self._all_filenames else selected
        if total and selected != total:
            self.file_count_label.setText(f"Files: {selected}/{total} selected")
        else:
            self.file_count_label.setText(f"Files: {selected}")

    def select_all_frames(self) -> None:
        self._apply_frame_selection(set(self._all_filenames), update_table=True, write_sidecars=True)

    def clear_all_frames(self) -> None:
        self._apply_frame_selection(set(), update_table=True, write_sidecars=True)

    def use_selected_table_rows_only(self) -> None:
        rows = sorted({idx.row() for idx in self.header_table.selectedIndexes()})
        if not rows:
            QMessageBox.information(self, "No Rows Selected", "사용할 프레임 행을 먼저 선택하세요.")
            return
        selected = set()
        for row in rows:
            item = self.header_table.item(row, 0)
            if item is not None and item.data(Qt.UserRole):
                selected.add(str(item.data(Qt.UserRole)))
        self._apply_frame_selection(selected, update_table=True, write_sidecars=True)

    def _set_data_dir(self, path: Path) -> None:
        self.params.P.data_dir = Path(path)
        self._update_result_workspace()
        self._ensure_result_dirs()

        if hasattr(self, "file_count_label"):
            self.file_count_label.setText("Files: 0")
        self.last_data_dir = self.params.P.data_dir

        self._persist_param_file(
            io_updates={
                "data_dir": str(self.params.P.data_dir),
                "result_dir": str(self.params.P.result_dir),
                "cache_dir": str(self.params.P.cache_dir.name),
            }
        )

    def _ensure_result_dirs(self) -> None:
        try:
            self.params.P.result_dir.mkdir(parents=True, exist_ok=True)
            self.params.P.cache_dir.mkdir(parents=True, exist_ok=True)
        except OSError:
            pass

    def _persist_param_file(self, io_updates=None, target_updates=None) -> None:
        if tomli_w is None:
            if hasattr(self.params, "save_toml"):
                try:
                    self.params.save_toml()
                except Exception:
                    pass
            return

        param_path = Path(getattr(self.params, "param_file", "parameters.toml"))
        if not param_path.exists():
            return

        try:
            with param_path.open("rb") as f:
                data = tomllib.load(f)
        except Exception:
            return

        io_updates = io_updates or {}
        target_updates = target_updates or {}

        if io_updates:
            io_block = data.get("io", {}) if isinstance(data.get("io", {}), dict) else {}
            io_block.update(io_updates)
            data["io"] = io_block

        if target_updates:
            target_block = data.get("target", {}) if isinstance(data.get("target", {}), dict) else {}
            target_block.update(target_updates)
            data["target"] = target_block

        try:
            with param_path.open("wb") as f:
                tomli_w.dump(data, f)
        except Exception:
            return

    # ------------------------------------------------------------------
    # Target lookup
    # ------------------------------------------------------------------

    def resolve_target(self) -> None:
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
                tried = getattr(inst, "last_target_attempts", [])
                errors = getattr(inst, "last_target_errors", [])
                tried_text = f"\nTried: {', '.join(tried[:8])}" if tried else ""
                error_text = f"\nLast error: {errors[-1]}" if errors else ""
                QMessageBox.warning(self, "SIMBAD", f"Target not resolved: {name}{tried_text}{error_text}")
                return

            ra_deg = float(inst.primary_coord.ra.deg)
            dec_deg = float(inst.primary_coord.dec.deg)
            ra_hms = inst.primary_coord.ra.to_string(unit="hour", sep=":", precision=2)
            dec_dms = inst.primary_coord.dec.to_string(
                unit="deg",
                sep=":",
                precision=1,
                alwayssign=True,
            )
            self.target_result.setText(f"{ra_hms}, {dec_dms}")

            self.params.P.target_name = name
            self.params.P.target_ra_deg = ra_deg
            self.params.P.target_dec_deg = dec_deg
            self._persist_param_file(
                target_updates={
                    "name": name,
                    "ra_deg": ra_deg,
                    "dec_deg": dec_deg,
                }
            )
            self.after_target_resolved()

        except Exception as e:
            QMessageBox.warning(self, "SIMBAD Error", str(e))

    # ------------------------------------------------------------------
    # State persistence
    # ------------------------------------------------------------------

    def save_state(self) -> None:
        state_data = {
            "data_dir": str(self.params.P.data_dir),
            "result_dir": str(self.params.P.result_dir),
            "root_dir": str(self.root_dir),
            "filename_prefix": self.params.P.filename_prefix,
            "file_count": len(self.file_manager.filenames),
            "all_frame_count": len(self._all_filenames),
            "selected_frame_filenames": list(self.file_manager.filenames),
            "file_path_map": {
                k: str(v) for k, v in getattr(self.file_manager, "path_map", {}).items()
            },
        }
        state_data.update(self.extra_state())

        self.project_state.store_step_data("file_selection", state_data)
        self._write_common_step1_sidecars()
        self.after_common_state_saved()

    def _write_common_step1_sidecars(self) -> None:
        try:
            s1_dir = step1_dir(self.params.P.result_dir)
            s1_dir.mkdir(parents=True, exist_ok=True)
            path_map_data = {
                k: str(v) for k, v in getattr(self.file_manager, "path_map", {}).items() if v
            }
            (s1_dir / "target_list.txt").write_text(
                "\n".join(self.file_manager.filenames),
                encoding="utf-8",
            )
            if self.file_manager.df_headers is not None:
                self.file_manager.df_headers.to_csv(
                    s1_dir / "headers.csv",
                    index=False,
                    encoding="utf-8",
                )
            (s1_dir / "file_path_map.json").write_text(
                json.dumps(path_map_data, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
        except Exception:
            pass

    def restore_state(self) -> None:
        state_data = self.project_state.get_step_data("file_selection")

        if state_data:
            saved_result_dir = state_data.get("result_dir")
            if "data_dir" in state_data:
                self.params.P.data_dir = Path(state_data["data_dir"])
                self.root_dir = Path(state_data.get("root_dir") or state_data["data_dir"])
                self.dir_edit.setText(str(self.root_dir))
                if saved_result_dir:
                    self.params.P.result_dir = Path(saved_result_dir)
                    self.params.P.cache_dir = self.params.P.result_dir / "cache"
                else:
                    self._update_result_workspace()

            if "filename_prefix" in state_data:
                self.params.P.filename_prefix = state_data["filename_prefix"]
                self.prefix_edit.setText(state_data["filename_prefix"])

            selected_frames = state_data.get("selected_frame_filenames")
            if isinstance(selected_frames, list) and selected_frames:
                self._selected_frame_filenames = {str(fn) for fn in selected_frames if fn}

            file_path_map = state_data.get("file_path_map")
            if isinstance(file_path_map, dict) and file_path_map:
                self.params.P.file_path_map = {
                    str(k): str(v) for k, v in file_path_map.items() if v
                }
                self.file_manager.path_map = {
                    str(k): Path(v) for k, v in file_path_map.items() if v
                }

            self.restore_extra_state_before_load(state_data)
            if not saved_result_dir and "data_dir" in state_data:
                self._update_result_workspace()
            self._ensure_result_dirs()

            try:
                self.load_files()
            except Exception:
                pass

        self._restore_target_from_params()

    def _restore_target_from_params(self) -> None:
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
