"""LC-only multi-night setup for the shared Step 1 file intake window."""

from __future__ import annotations

import json
import os
from collections import defaultdict
from pathlib import Path

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import (
    QCheckBox,
    QDoubleSpinBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QAbstractItemView,
    QLabel,
    QListView,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QTreeView,
)

from apex.utils.run_workspace import build_result_workspace_dir, write_run_manifest
from apex.utils.step_paths_lc import step1_dir


LC_STEP1_GUIDE = (
    "단일 night: Browse...로 FITS 폴더를 고르고 Filename Prefix를 확인한 뒤 Rescan Files를 누릅니다.\n"
    "여러 night: 입력 폴더 추가...로 사용할 night 폴더만 고르거나, Root 선택 후 하위폴더 포함을 켜고 Rescan Files를 누릅니다.\n"
    "Night gap은 JD 시간 간격으로 night를 나누는 기준입니다. 잘못 나뉘면 값을 조정하고 다시 Rescan Files를 누르세요.\n"
    "스캔 후 FITS Headers와 Night 표를 확인하고 사용할 프레임만 Use로 남긴 뒤 Next Step으로 이동합니다."
)


def _classify_nights_by_jd_gap(records: list, gap_hours: float = 8.0) -> list:
    """Assign 1-based night_id values using sorted JD gaps."""
    gap_days = gap_hours / 24.0
    valid = [(i, r) for i, r in enumerate(records) if r.get("jd") is not None]
    valid.sort(key=lambda x: x[1]["jd"])

    night_id = 1
    prev_jd = None
    id_map = {}
    for i, r in valid:
        if prev_jd is not None and (r["jd"] - prev_jd) > gap_days:
            night_id += 1
        id_map[i] = night_id
        prev_jd = r["jd"]

    return [{**r, "night_id": id_map.get(i, 1)} for i, r in enumerate(records)]


def _build_night_summary(records_with_night: list, tz_offset_hours: float = 0.0) -> list:
    """Build display rows for the LC night-selection table."""
    nights = defaultdict(lambda: {
        "jd_min": None,
        "jd_max": None,
        "files": [],
        "filters": defaultdict(int),
    })

    for r in records_with_night:
        night_id = r["night_id"]
        jd = r.get("jd")
        if jd is not None:
            if nights[night_id]["jd_min"] is None or jd < nights[night_id]["jd_min"]:
                nights[night_id]["jd_min"] = jd
            if nights[night_id]["jd_max"] is None or jd > nights[night_id]["jd_max"]:
                nights[night_id]["jd_max"] = jd
        nights[night_id]["files"].append(r["filename"])
        filt = str(r.get("filter", "") or "").strip()
        if filt:
            nights[night_id]["filters"][filt] += 1

    tz_days = float(tz_offset_hours) / 24.0

    def _jd_to_timestr(jd):
        if jd is None:
            return "?"
        try:
            from astropy.time import Time

            t = Time(jd + tz_days, format="jd", scale="utc")
            return t.datetime.strftime("%Y-%m-%d %H:%M")
        except Exception:
            return "?"

    result = []
    for night_id in sorted(nights):
        night = nights[night_id]
        first = _jd_to_timestr(night["jd_min"])
        last = _jd_to_timestr(night["jd_max"])
        date_label = first.split(" ")[0] if " " in first else first
        time_range = f"{first.split(' ')[-1]} ~ {last.split(' ')[-1]}" if first != "?" else "?"
        filter_str = "  ".join(f"{k}:{v}" for k, v in sorted(night["filters"].items()))
        result.append({
            "night_id": night_id,
            "label": f"N{night_id}",
            "date": date_label,
            "time_range": time_range,
            "file_count": len(night["files"]),
            "filter_str": filter_str or "-",
            "filenames": night["files"],
        })
    return result


class LightCurveNightSetupMixin:
    """Adds manual folders, subfolder scanning, and night selection to Step 1."""

    directory_label = "Root:"
    header_labels = [
        "Filename", "Night", "DATE-OBS", "FILTER", "EXPTIME", "AIRMASS",
        "OBJECT", "RA_DEG", "DEC_DEG", "IMAGETYP",
    ]

    def init_mode_state(self) -> None:
        self._manual_input_dirs: list[Path] = []
        self._night_records: list = []
        self._all_night_assignments: dict[str, int] = {}
        self._excluded_nights: set[int] = set()
        self._night_table_loading = False

    # ------------------------------------------------------------------
    # UI extension
    # ------------------------------------------------------------------

    def setup_directory_extension(self, dir_layout) -> None:
        guide_row = QHBoxLayout()
        guide_title = QLabel("LC Step 1 Guide")
        guide_title.setStyleSheet("QLabel { font-weight: bold; color: #263238; }")
        guide_row.addWidget(guide_title)
        guide_row.addWidget(self._make_help_button("LC Step 1 Guide", LC_STEP1_GUIDE))
        guide_row.addStretch()
        dir_layout.addLayout(guide_row)

        guide_label = QLabel(
            "단일: Browse... -> Rescan Files -> 사용할 프레임 Use 확인.  "
            "여러 night: 입력 폴더 추가... 또는 하위폴더 포함 -> Rescan Files."
        )
        guide_label.setWordWrap(True)
        guide_label.setStyleSheet(
            "QLabel { background:#E3F2FD; color:#263238; padding:6px; border-radius:4px; }"
        )
        dir_layout.addWidget(guide_label)

        pick_row = QHBoxLayout()

        btn_add_inputs = QPushButton("입력 폴더 추가...")
        btn_add_inputs.setToolTip("여러 night 폴더를 직접 선택합니다.")
        btn_add_inputs.clicked.connect(self.add_input_directories)
        pick_row.addWidget(btn_add_inputs)

        btn_remove_input = QPushButton("선택 제거")
        btn_remove_input.setToolTip("입력 폴더 목록에서 선택한 폴더를 제거합니다.")
        btn_remove_input.clicked.connect(self.remove_selected_input_directory)
        pick_row.addWidget(btn_remove_input)

        btn_clear_inputs = QPushButton("입력 초기화")
        btn_clear_inputs.setToolTip("직접 선택한 입력 폴더를 모두 지우고 Root 단일 폴더 사용으로 되돌립니다.")
        btn_clear_inputs.clicked.connect(self.clear_input_directories)
        pick_row.addWidget(btn_clear_inputs)
        pick_row.addWidget(self._make_help_button(
            "입력 폴더",
            "단일 night는 Browse...로 선택한 Root 폴더를 그대로 씁니다.\n"
            "여러 night는 입력 폴더 추가...로 각 night 폴더를 직접 고르는 방식이 가장 안전합니다.\n"
            "직접 선택한 입력 폴더가 있으면 하위폴더 포함 옵션보다 우선합니다.",
        ))
        pick_row.addStretch()
        dir_layout.addLayout(pick_row)

        self.input_dir_list = QListWidget()
        self.input_dir_list.setMaximumHeight(90)
        self.input_dir_list.setSelectionMode(QListWidget.SingleSelection)
        self.input_dir_list.setToolTip("직접 선택한 night 입력 폴더 목록입니다.")
        dir_layout.addWidget(self.input_dir_list)

        self.input_dir_info = QLabel("입력 폴더: 루트 폴더 단일 사용")
        self.input_dir_info.setStyleSheet("QLabel { color: #555; }")
        dir_layout.addWidget(self.input_dir_info)

        options_row = QHBoxLayout()
        self.include_subfolders_check = QCheckBox("하위폴더 포함")
        self.include_subfolders_check.setToolTip("Root 바로 아래 1단계 하위폴더 중 FITS가 있는 폴더를 night 입력으로 사용합니다.")
        options_row.addWidget(self.include_subfolders_check)
        options_row.addWidget(self._make_help_button(
            "하위폴더 포함",
            "Root 바로 아래 폴더만 검사합니다. 예: Root/2025-01-01, Root/2025-01-02.\n"
            "더 깊은 재귀 탐색은 하지 않습니다.\n"
            "폴더를 정확히 통제하려면 입력 폴더 추가...를 사용하세요.",
        ))
        options_row.addWidget(QLabel("Night gap:"))

        self.night_gap_spinbox = QDoubleSpinBox()
        self.night_gap_spinbox.setRange(1.0, 24.0)
        self.night_gap_spinbox.setSingleStep(0.5)
        self.night_gap_spinbox.setDecimals(1)
        self.night_gap_spinbox.setValue(float(getattr(self.params.P, "night_gap_hours", 8.0)))
        self.night_gap_spinbox.setMaximumWidth(70)
        self.night_gap_spinbox.setToolTip("연속 프레임 JD 간격이 이 시간보다 크면 다음 night로 분리합니다.")
        options_row.addWidget(self.night_gap_spinbox)
        options_row.addWidget(QLabel("h"))
        options_row.addWidget(self._make_help_button(
            "Night gap",
            "기본값 8시간은 관측일 사이의 긴 공백을 기준으로 night를 나눕니다.\n"
            "한 밤이 둘로 갈라지면 값을 키우고, 서로 다른 밤이 합쳐지면 값을 줄인 뒤 Rescan Files를 다시 누르세요.",
        ))
        options_row.addStretch()
        dir_layout.addLayout(options_row)

        night_btn_row = QHBoxLayout()
        btn_select_all_nights = QPushButton("전체 선택")
        btn_select_all_nights.setToolTip("Night 표의 Use 체크를 모두 켭니다.")
        btn_select_all_nights.clicked.connect(self._select_all_displayed_nights)
        night_btn_row.addWidget(btn_select_all_nights)

        btn_clear_nights = QPushButton("선택 해제")
        btn_clear_nights.setToolTip("Night 표의 Use 체크를 모두 끄고 해당 night 프레임을 제외합니다.")
        btn_clear_nights.clicked.connect(self._clear_all_displayed_nights)
        night_btn_row.addWidget(btn_clear_nights)
        night_btn_row.addWidget(self._make_help_button(
            "Night 표",
            "Rescan Files 후 JD 기준으로 분류된 night 요약입니다.\n"
            "관측일, 시간범위, 파일수, 필터가 예상과 맞는지 확인하세요.\n"
            "Use를 끄면 해당 night의 프레임들이 FITS Headers 표에서도 제외됩니다.",
        ))
        night_btn_row.addStretch()
        dir_layout.addLayout(night_btn_row)

        self.night_table = QTableWidget()
        self.night_table.setColumnCount(6)
        self.night_table.setHorizontalHeaderLabels(
            ["Use", "Night", "관측일", "시간범위", "파일수", "필터"]
        )
        self.night_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.night_table.horizontalHeader().setStretchLastSection(True)
        self.night_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.night_table.itemChanged.connect(self._on_night_item_changed)
        self.night_table.setMaximumHeight(160)
        self.night_table.setToolTip("Rescan Files 후 분류된 night 요약입니다. 예상과 다르면 Night gap 또는 입력 폴더를 조정하세요.")
        use_header = self.night_table.horizontalHeaderItem(0)
        if use_header is not None:
            use_header.setToolTip("체크된 night의 프레임만 이후 단계에서 처리합니다.")
        dir_layout.addWidget(self.night_table)

    def configure_header_table(self) -> None:
        self.header_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.header_table.horizontalHeader().setStretchLastSection(True)
        self.header_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.header_table.setSelectionMode(QAbstractItemView.ExtendedSelection)
        self.header_table.setEditTriggers(QTableWidget.NoEditTriggers)

    # ------------------------------------------------------------------
    # Directory and scan extension
    # ------------------------------------------------------------------

    def _pick_multiple_directories(self, start_dir: Path) -> list[Path]:
        dialog = QFileDialog(self, "입력 폴더 선택", str(start_dir))
        dialog.setFileMode(QFileDialog.Directory)
        dialog.setOption(QFileDialog.ShowDirsOnly, True)
        dialog.setOption(QFileDialog.DontUseNativeDialog, True)
        for view in dialog.findChildren((QListView, QTreeView)):
            view.setSelectionMode(view.ExtendedSelection)
        if dialog.exec_() != QFileDialog.Accepted:
            return []
        return [Path(p) for p in dialog.selectedFiles() if p]

    def _sync_input_dir_widgets(self) -> None:
        self.input_dir_list.clear()
        for path in self._manual_input_dirs:
            self.input_dir_list.addItem(QListWidgetItem(str(path)))
        if self._manual_input_dirs:
            self.input_dir_info.setText(
                f"입력 폴더: {len(self._manual_input_dirs)}개 직접 선택 "
                f"(하위폴더 포함보다 우선)"
            )
        else:
            self.input_dir_info.setText("입력 폴더: 루트 폴더 단일 사용")

    def _apply_manual_input_dirs(self) -> None:
        if self._manual_input_dirs:
            common_root = Path(os.path.commonpath([str(p) for p in self._manual_input_dirs]))
            self.root_dir = common_root
            self.dir_edit.setText(str(self.root_dir))
            self.params.P.data_dir = self.root_dir
            self.file_manager.set_multi_night_dirs(self.root_dir, self._manual_input_dirs)
        else:
            self.file_manager.clear_multi_night_dirs()

        self._sync_input_dir_widgets()
        self._update_result_workspace()
        self._ensure_result_dirs()
        self._persist_param_file(
            io_updates={
                "data_dir": str(self.params.P.data_dir),
                "result_dir": str(self.params.P.result_dir),
                "cache_dir": str(self.params.P.cache_dir.name),
            }
        )

        try:
            self.save_state()
            self.update_navigation_buttons()
        except Exception:
            pass

    def add_input_directories(self) -> None:
        picked = self._pick_multiple_directories(self.root_dir)
        if not picked:
            return

        existing = {str(p.resolve()) for p in self._manual_input_dirs}
        added = False
        for path in picked:
            try:
                key = str(path.resolve())
            except Exception:
                key = str(path)
            if key in existing:
                continue
            self._manual_input_dirs.append(Path(path))
            existing.add(key)
            added = True

        if added:
            self._manual_input_dirs = sorted(self._manual_input_dirs, key=lambda p: str(p))
            self._apply_manual_input_dirs()

    def remove_selected_input_directory(self) -> None:
        row = self.input_dir_list.currentRow()
        if 0 <= row < len(self._manual_input_dirs):
            self._manual_input_dirs.pop(row)
            self._apply_manual_input_dirs()

    def clear_input_directories(self) -> None:
        if not self._manual_input_dirs:
            return

        self._manual_input_dirs = []
        self._clear_loaded_file_state()
        self.file_manager.night_assignments = {}
        self.file_manager.excluded_nights = set()
        self._night_records = []
        self._all_night_assignments = {}
        self._excluded_nights = set()
        self.night_table.setRowCount(0)

        self._apply_manual_input_dirs()
        self.rescan_files()

    def reset_mode_state_for_new_directory(self) -> None:
        self._manual_input_dirs = []
        self._night_records = []
        self._all_night_assignments = {}
        self._excluded_nights = set()
        self.file_manager.night_assignments = {}
        self.file_manager.excluded_nights = set()
        self.file_manager.clear_multi_night_dirs()
        self._sync_input_dir_widgets()
        self.night_table.setRowCount(0)

    def prepare_scan_context(self) -> None:
        gap = self.night_gap_spinbox.value()
        self.params.P.night_gap_hours = gap
        self._persist_param_file(io_updates={"night_gap_hours": gap})

        if self._manual_input_dirs:
            self.file_manager.set_multi_night_dirs(self.root_dir, self._manual_input_dirs)
        elif self.include_subfolders_check.isChecked():
            self._setup_subdirectory_scan()
        else:
            self.file_manager.clear_multi_night_dirs()

    def _setup_subdirectory_scan(self) -> None:
        root = Path(self.root_dir)
        if not root.exists():
            return

        prefix = self.prefix_edit.text().strip().lower()
        suffixes = (".fit", ".fits", ".fit.fz", ".fits.fz")

        sub_dirs = []
        for sub in sorted(p for p in root.iterdir() if p.is_dir()):
            has_fits = any(
                f.is_file()
                and f.name.lower().startswith(prefix)
                and f.name.lower().endswith(suffixes)
                for f in sub.iterdir()
            )
            if has_fits:
                sub_dirs.append(sub)

        if sub_dirs:
            self.file_manager.set_multi_night_dirs(root, sub_dirs)
        else:
            self.file_manager.clear_multi_night_dirs()

    def _current_input_dirs(self) -> list[Path]:
        return [Path(p) for p in getattr(self.file_manager, "selected_dirs", []) if p]

    def _update_result_workspace(self) -> None:
        root_dir = Path(self.root_dir if getattr(self, "root_dir", None) else self.params.P.data_dir)
        self.params.P.result_dir = build_result_workspace_dir(root_dir, self._current_input_dirs())
        self.params.P.cache_dir = self.params.P.result_dir / "cache"

    # ------------------------------------------------------------------
    # Night classification
    # ------------------------------------------------------------------

    def after_headers_loaded(self, df_headers) -> None:
        self._classify_and_display(df_headers)

    def _classify_and_display(self, df_headers) -> None:
        import numpy as np

        records = []
        for _, row in df_headers.iterrows():
            filename = str(row["Filename"])
            date_obs = str(row.get("DATE-OBS", "") or "")
            filt = str(row.get("FILTER", "") or "")

            jd = None
            if "JD" in df_headers.columns:
                try:
                    jd_raw = float(row["JD"])
                    if np.isfinite(jd_raw):
                        jd = jd_raw
                except (TypeError, ValueError):
                    pass

            if jd is None and date_obs and date_obs not in ("N/A", "nan", ""):
                try:
                    from astropy.time import Time

                    jd = float(Time(date_obs.strip(), format="isot", scale="utc").jd)
                except Exception:
                    jd = None

            records.append({
                "filename": filename,
                "jd": jd,
                "date_obs": date_obs,
                "filter": filt,
            })

        records_with_night = _classify_nights_by_jd_gap(
            records,
            self.night_gap_spinbox.value(),
        )
        self._night_records = records_with_night
        self._all_night_assignments = {
            r["filename"]: r["night_id"] for r in records_with_night
        }
        self.file_manager.night_assignments = dict(self._all_night_assignments)

        tz = float(getattr(self.params.P, "site_tz_offset_hours", 0.0))
        self._populate_night_table(_build_night_summary(records_with_night, tz_offset_hours=tz))
        self.populate_header_table(df_headers)

    def _populate_night_table(self, summaries: list) -> None:
        self._night_table_loading = True
        self.night_table.blockSignals(True)
        self.night_table.setRowCount(0)

        for summary in summaries:
            night_id = summary["night_id"]
            row = self.night_table.rowCount()
            self.night_table.insertRow(row)

            use_item = QTableWidgetItem()
            use_item.setFlags(use_item.flags() | Qt.ItemIsUserCheckable)
            checked = Qt.Unchecked if night_id in self._excluded_nights else Qt.Checked
            use_item.setCheckState(checked)
            use_item.setData(Qt.UserRole, night_id)

            self.night_table.setItem(row, 0, use_item)
            self.night_table.setItem(row, 1, QTableWidgetItem(summary["label"]))
            self.night_table.setItem(row, 2, QTableWidgetItem(summary["date"]))
            self.night_table.setItem(row, 3, QTableWidgetItem(summary["time_range"]))
            self.night_table.setItem(row, 4, QTableWidgetItem(str(summary["file_count"])))
            self.night_table.setItem(row, 5, QTableWidgetItem(summary["filter_str"]))

        self.night_table.blockSignals(False)
        self._night_table_loading = False

    def populate_header_table(self, df_headers) -> None:
        night_map = {r["filename"]: r["night_id"] for r in self._night_records}
        self._header_table_loading = True
        self.header_table.blockSignals(True)
        self.header_table.setRowCount(len(df_headers))
        for i, row in df_headers.iterrows():
            filename = str(row["Filename"])
            night_id = night_map.get(filename, "-")
            use_item = self._make_frame_use_item(filename)
            if night_id != "-" and int(night_id) in self._excluded_nights:
                use_item.setCheckState(Qt.Unchecked)
            self.header_table.setItem(i, 0, use_item)
            values = {
                "Filename": filename,
                "Night": f"N{night_id}",
            }
            for label in self.header_labels:
                values.setdefault(label, self._format_header_cell(row, label))
            for col, label in enumerate(self.header_labels, start=1):
                self.header_table.setItem(i, col, QTableWidgetItem(values[label]))
        self.header_table.blockSignals(False)
        self._header_table_loading = False

    def after_frame_selection_applied(self, selected_filenames: set[str]) -> None:
        source = self._all_night_assignments or self.file_manager.night_assignments
        self.file_manager.night_assignments = {
            fn: night_id for fn, night_id in source.items() if fn in selected_filenames
        }
        self._sync_night_table_checks(selected_filenames)

    def _filenames_for_night(self, night_id: int) -> set[str]:
        return {
            fn for fn, assigned_night in self._all_night_assignments.items()
            if int(assigned_night) == int(night_id)
        }

    def _sync_night_table_checks(self, selected_filenames: set[str]) -> None:
        if not hasattr(self, "night_table") or not self._all_night_assignments:
            return

        excluded = set()
        self._night_table_loading = True
        self.night_table.blockSignals(True)
        for row in range(self.night_table.rowCount()):
            item = self.night_table.item(row, 0)
            if item is None:
                continue
            night_id = int(item.data(Qt.UserRole))
            night_filenames = self._filenames_for_night(night_id)
            selected_in_night = night_filenames & selected_filenames
            if not selected_in_night:
                item.setCheckState(Qt.Unchecked)
                excluded.add(night_id)
            elif selected_in_night == night_filenames:
                item.setCheckState(Qt.Checked)
            else:
                item.setCheckState(Qt.PartiallyChecked)
        self.night_table.blockSignals(False)
        self._night_table_loading = False
        self._excluded_nights = excluded
        self.file_manager.excluded_nights = excluded

    def _on_night_item_changed(self, item: QTableWidgetItem) -> None:
        if self._night_table_loading or item.column() != 0:
            return
        night_id = item.data(Qt.UserRole)
        if night_id is None:
            return
        selected = self._selected_filenames_from_table()
        night_filenames = self._filenames_for_night(int(night_id))
        if item.checkState() == Qt.Unchecked:
            selected.difference_update(night_filenames)
        else:
            selected.update(night_filenames)
        self._apply_frame_selection(selected, update_table=True, write_sidecars=True)

    def _select_all_displayed_nights(self) -> None:
        self._apply_frame_selection(
            set(self._all_night_assignments.keys()),
            update_table=True,
            write_sidecars=True,
        )

    def _clear_all_displayed_nights(self) -> None:
        self._apply_frame_selection(set(), update_table=True, write_sidecars=True)

    # ------------------------------------------------------------------
    # State and LC sidecars
    # ------------------------------------------------------------------

    def extra_state(self) -> dict:
        return {
            "root_dir": str(self.root_dir),
            "include_subfolders": bool(self.include_subfolders_check.isChecked()),
            "multi_night": bool(self.file_manager.selected_dirs),
            "night_dirs": [str(p) for p in self.file_manager.selected_dirs],
            "manual_input_dirs": [str(p) for p in self._manual_input_dirs],
            "night_gap_hours": self.night_gap_spinbox.value(),
            "excluded_nights": sorted(self._excluded_nights),
            "night_assignments": {
                k: v for k, v in self.file_manager.night_assignments.items()
            },
        }

    def restore_extra_state_before_load(self, state_data: dict) -> None:
        if "include_subfolders" in state_data:
            self.include_subfolders_check.setChecked(bool(state_data["include_subfolders"]))

        if "night_gap_hours" in state_data:
            self.night_gap_spinbox.setValue(float(state_data["night_gap_hours"]))

        if "excluded_nights" in state_data:
            self._excluded_nights = set(int(n) for n in state_data["excluded_nights"])
            self.file_manager.excluded_nights = self._excluded_nights

        if "night_assignments" in state_data:
            self.file_manager.night_assignments = {
                k: int(v) for k, v in state_data["night_assignments"].items()
            }

        self._manual_input_dirs = [
            Path(p) for p in state_data.get("manual_input_dirs", []) if p
        ]
        self._sync_input_dir_widgets()

        root_path = Path(
            state_data.get("root_dir")
            or state_data.get("data_dir")
            or self.params.P.data_dir
        )
        self.root_dir = root_path
        self.dir_edit.setText(str(self.root_dir))

        night_dirs = [Path(p) for p in state_data.get("night_dirs", []) if p]
        if self._manual_input_dirs:
            self.file_manager.set_multi_night_dirs(root_path, self._manual_input_dirs)
        elif bool(state_data.get("multi_night")) and night_dirs:
            self.file_manager.set_multi_night_dirs(root_path, night_dirs)
        else:
            self.file_manager.clear_multi_night_dirs()

    def after_common_state_saved(self) -> None:
        try:
            s1_dir = step1_dir(self.params.P.result_dir)
            s1_dir.mkdir(parents=True, exist_ok=True)
            night_data = {
                "night_assignments": {
                    k: v for k, v in self.file_manager.night_assignments.items()
                },
                "excluded_nights": sorted(self._excluded_nights),
            }
            (s1_dir / "night_assignments.json").write_text(
                json.dumps(night_data, indent=2, ensure_ascii=False),
                encoding="utf-8",
            )
            self._write_run_manifest()
        except Exception:
            pass

    def _write_run_manifest(self) -> None:
        write_run_manifest(
            self.params.P.result_dir,
            run_type="result",
            root_dir=Path(self.root_dir if getattr(self, "root_dir", None) else self.params.P.data_dir),
            input_dirs=self._current_input_dirs(),
            target_name=getattr(self.params.P, "target_name", None),
        )

    def after_target_resolved(self) -> None:
        try:
            self._write_run_manifest()
        except Exception:
            pass
