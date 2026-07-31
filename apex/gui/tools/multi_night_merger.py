"""Tools-native multi-night merger workflow."""

from __future__ import annotations

import copy
from pathlib import Path

import numpy as np
import pandas as pd

from PyQt5.QtWidgets import (
    QApplication,
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QPushButton,
    QLabel,
    QGroupBox,
    QListWidget,
    QListWidgetItem,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QTextEdit,
    QFileDialog,
    QMessageBox,
    QLineEdit,
    QComboBox,
    QAbstractItemView,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QColor

from apex.analysis.merge.id_match import reconcile_workspace_catalogs
from apex.analysis.merge.reference_store import STORAGE_FULL, STORAGE_REFERENCE
from apex.analysis.merge.workspace_build import materialize_merged_workspace
from apex.analysis.merge.workspace_scan import (
    default_merged_output_dir as _default_output_dir,
    folder_tag as _folder_tag,
    load_master_catalogs_by_filter as _load_master_catalogs_by_filter,
    load_selection_payloads as _load_selection_payloads,
    read_step7_index as _read_step7_index,
    scan_merge_input_workspace,
    workspace_scan_signature,
)
from apex.core.project_state import ProjectState
from apex.gui.layout_rules import AutoFitMixin
from apex.gui.theme import Tokens, style_button
from apex.utils.io_utils import (
    coerce_int64_source_id,
)
from apex.utils.run_workspace import (
    build_merged_workspace_dir,
)
from apex.utils.step_paths_lc import (
    step8_selection_dir,
)


def compute_selection_defaults(
    merged_catalogs: dict[str, pd.DataFrame],
    selection_payloads: "list[dict[str, dict]] | dict[str, dict]",
) -> tuple[dict[str, "int | None"], dict[str, "set[int]"], dict[str, "int | None"]]:
    """Derive target/comparison/check selections per filter from the input
    workspaces' selection payloads, keeping only source_ids that survive into
    the merged catalog. Pure — safe to run inside the merge worker thread.

    ``selection_payloads`` is the per-folder payload list in priority order
    (base first). The base folder wins where it has a selection; a target or
    check it lacks is taken from the next folder that has one, and comparison
    stars are unioned — previously only the base folder was consulted, so a
    target picked on the second night was silently dropped.
    """
    if isinstance(selection_payloads, dict):        # single-folder (legacy) form
        selection_payloads = [selection_payloads]

    sel_t: dict[str, int | None] = {}
    sel_c: dict[str, set[int]] = {}
    sel_k: dict[str, int | None] = {}
    for flt, df in merged_catalogs.items():
        if "source_id" in df.columns:
            available = set(
                coerce_int64_source_id(df["source_id"]).dropna().astype("int64").tolist()
            )
        else:
            available = set()

        def _keep(sid) -> "int | None":
            if sid is None:
                return None
            try:
                sid = int(sid)
            except (TypeError, ValueError):
                return None
            return sid if sid in available else None

        target: int | None = None
        check: int | None = None
        comp_sids: set[int] = set()
        for payload_by_filter in selection_payloads:
            payload = payload_by_filter.get(flt, {})
            if target is None:
                target = _keep(payload.get("target_source_id"))
            if check is None:
                check = _keep(payload.get("check_source_id"))
            for sid in payload.get("comparison_source_ids", []) or ():
                kept = _keep(sid)
                if kept is not None:
                    comp_sids.add(kept)
        sel_t[flt] = target
        sel_c[flt] = comp_sids
        sel_k[flt] = check
    return sel_t, sel_c, sel_k


class _MergeWorker(QThread):
    """Runs the heavy reconcile + materialize off the GUI thread.

    Emits ``progress(str)`` log lines, ``finished_ok(dict)`` with the reconcile
    result + selection defaults + materialize build_info, or ``failed(str)``.
    """

    progress = pyqtSignal(str)
    finished_ok = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(
        self,
        *,
        folders,
        catalogs_by_folder,
        folder_tags,
        match_radius,
        base_selection_by_filter,
        out_dir,
        cached_reconcile,
        storage_mode=STORAGE_FULL,
        parent=None,
    ):
        super().__init__(parent)
        self._folders = folders
        self._catalogs_by_folder = catalogs_by_folder
        self._folder_tags = folder_tags
        self._match_radius = match_radius
        self._base_selection = base_selection_by_filter
        self._out_dir = out_dir
        self._cached_reconcile = cached_reconcile
        self._storage_mode = storage_mode

    def run(self):  # QThread entry point
        try:
            cached = self._cached_reconcile
            if cached is None:
                cached = reconcile_workspace_catalogs(
                    self._folders,
                    self._catalogs_by_folder,
                    self._folder_tags,
                    self._match_radius,
                    logger=self.progress.emit,
                )
            merged_catalogs = cached["canonical_by_filter"]
            sel_t, sel_c, sel_k = compute_selection_defaults(
                merged_catalogs, self._base_selection
            )
            build_info = materialize_merged_workspace(
                out_dir=self._out_dir,
                folders=self._folders,
                folder_tags=self._folder_tags,
                local_id_maps=cached["local_id_maps"],
                merged_catalogs=merged_catalogs,
                selection_target_by_filter=sel_t,
                selection_comp_by_filter=sel_c,
                selection_check_by_filter=sel_k,
                match_records=cached["match_records"],
                storage_mode=self._storage_mode,
            )
            self.finished_ok.emit(
                {
                    "cached": cached,
                    "sel_t": sel_t,
                    "sel_c": sel_c,
                    "sel_k": sel_k,
                    "build_info": build_info,
                    "out_dir": self._out_dir,
                }
            )
        except Exception as exc:  # surface any failure to the GUI thread
            import traceback

            self.failed.emit(f"{exc}\n\n{traceback.format_exc()}")


class _MergedParamsProxy:
    """Read-only-ish params wrapper for merged runtime windows."""

    def __init__(self, base_params, merged_result_dir: Path):
        self._base = base_params
        self.P = copy.deepcopy(base_params.P)
        self.P.result_dir = Path(merged_result_dir)
        self.P.cache_dir = Path(merged_result_dir) / "cache"
        if not hasattr(self.P, "file_path_map"):
            self.P.file_path_map = {}

    def __getattr__(self, name):
        return getattr(self._base, name)

    def get_file_path(self, filename: str) -> Path:
        """Resolve merged-runtime FITS paths against the merged path map."""
        path_map = getattr(self.P, "file_path_map", None)
        if isinstance(path_map, dict):
            mapped = path_map.get(filename)
            if mapped:
                return Path(mapped)
        return Path(self.P.data_dir) / filename

    def save_toml(self):
        return False


class _MergedFileManagerProxy:
    def __init__(self, filenames: list[str], night_assignments: dict[str, int], path_map: dict[str, str] | None = None):
        self.filenames = list(filenames)
        self.night_assignments = dict(night_assignments)
        self.path_map = dict(path_map or {})

    def get_file_path(self, filename: str) -> Path | None:
        p = self.path_map.get(filename)
        return Path(p) if p else None


class MultiNightMergerWindow(AutoFitMixin, QMainWindow):
    """Merger workflow for previously processed result folders."""

    STEP_TITLES = [
        "Step 1  폴더 선택",
        "Step 2  ID 매칭",
        "Step 3  선택",
        "Step 4  Light Curve",
        "Step 5  Detrend",
        "Step 6  Period",
    ]

    # For StepWindowBase child reuse.
    step_names = [
        "File Selection",
        "Crop",
        "Sky Preview",
        "Source Detection",
        "WCS",
        "Master Catalog Build",
        "Forced Phot",
        "PSF Photometry",
        "Selection",
        "Light Curve",
        "Detrend",
        "Period Analysis",
    ]

    def __init__(self, params, project_state, main_window):
        super().__init__()
        self.params = params
        self.project_state = project_state
        self.main_window = main_window

        self.folders: list[Path] = [Path(params.P.result_dir)]
        self.folder_tags: dict[str, str] = {}
        self.folder_scan_rows: list[dict] = []
        self._workspace_scan_cache: dict[str, tuple[tuple, dict]] = {}
        self._catalog_cache: dict[str, tuple[tuple, dict[str, pd.DataFrame]]] = {}
        self._selection_payload_cache: dict[str, tuple[tuple, dict[str, dict]]] = {}
        self._id_match_cache: dict[tuple, dict] = {}

        self.match_summary_rows: list[dict] = []
        self.match_records: list[dict] = []
        self.merged_catalogs: dict[str, pd.DataFrame] = {}
        self.local_id_maps: dict[str, dict[str, dict[int, dict[str, int]]]] = {}
        self.base_selection_by_filter: dict[str, dict] = {}

        self.selection_target_by_filter: dict[str, int | None] = {}
        self.selection_comp_by_filter: dict[str, set[int]] = {}
        self.selection_check_by_filter: dict[str, int | None] = {}

        self.merged_result_dir: Path | None = None
        self.merged_runtime_params = None
        self.merged_runtime_project_state: ProjectState | None = None
        self.merged_runtime_file_manager = None
        self.step9_embedded_window = None
        self.step10_embedded_window = None
        self.step11_embedded_window = None
        self.step12_embedded_window = None

        self._merge_worker: "_MergeWorker | None" = None
        self._merge_busy = False
        self._pending_id_sig: "tuple | None" = None

        self.setWindowTitle("Multi-Night Merger Workflow")
        self.resize(1200, 820)
        self._setup_ui()
        # AutoFitMixin fits to content and clamps to the monitor on first show.

    # ───────────────────────── UI ─────────────────────────

    def _setup_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        header = QHBoxLayout()
        header.setSpacing(Tokens.GAP)
        btn_back = QPushButton("← 메인으로")
        style_button(btn_back, height=Tokens.H_COMPACT)
        # sizeHint, not a fixed 110 px: at that width the arrow was clipped.
        back_width = btn_back.sizeHint().width()
        btn_back.setFixedWidth(back_width)
        btn_back.clicked.connect(self._go_back)
        header.addWidget(btn_back)

        title = QLabel("Multi-Night Merger Workflow")
        title.setProperty("role", "title")
        title.setAlignment(Qt.AlignCenter)
        header.addWidget(title, 1)
        header.addSpacing(back_width)     # keeps the title optically centred
        root.addLayout(header)

        step_bar = QHBoxLayout()
        step_bar.setSpacing(Tokens.GAP // 2)
        self._step_btns: list[QPushButton] = []
        for i, text in enumerate(self.STEP_TITLES):
            btn = QPushButton(text)
            btn.setCheckable(True)
            style_button(btn, height=Tokens.H_BUTTON)
            btn.clicked.connect(lambda checked, idx=i: self._go_to_step(idx))
            step_bar.addWidget(btn)
            self._step_btns.append(btn)
        root.addLayout(step_bar)

        self._pages: list[QWidget] = [
            self._make_step1(),
            self._make_step2(),
            self._make_step3(),
            self._make_step4(),
            self._make_step5(),
            self._make_step6(),
        ]

        self._page_container = QWidget()
        self._page_layout = QVBoxLayout(self._page_container)
        self._page_layout.setContentsMargins(0, 0, 0, 0)
        for page in self._pages:
            self._page_layout.addWidget(page)
            page.hide()
        root.addWidget(self._page_container, 1)

        nav = QHBoxLayout()
        nav.setSpacing(Tokens.GAP)
        self.btn_prev = QPushButton("◀ 이전")
        style_button(self.btn_prev, height=Tokens.H_ACTION)
        self.btn_prev.clicked.connect(self._prev_step)
        self.btn_next = QPushButton("다음 ▶")
        style_button(self.btn_next, "primary", height=Tokens.H_ACTION)
        self.btn_next.clicked.connect(self._next_step)
        nav.addWidget(self.btn_prev)
        nav.addStretch()
        nav.addWidget(self.btn_next)
        root.addLayout(nav)

        self._current_step = 0
        self._refresh_output_dir_default(force=True)
        self._refresh_folder_list()
        self._go_to_step(0)

    def _make_info_label(self, text: str) -> QLabel:
        label = QLabel(text)
        label.setWordWrap(True)
        label.setProperty("role", "info")     # themed chip, follows dark mode
        return label

    def _make_status_label(self, text: str) -> QLabel:
        label = QLabel(text)
        label.setProperty("role", "caption")
        return label

    def _make_step1(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)

        layout.addWidget(self._make_info_label(
            "RESULT_* 또는 MERGED_* workspace 폴더들을 선택합니다.\n"
            "각 폴더에 Step 7(forced photometry)과 Step 8(마스터 카탈로그·선택)이 "
            "있으면 됩니다. Step 9 광곡선은 merged workspace에서 다시 만들므로 "
            "밤마다 미리 돌릴 필요가 없습니다.\n"
            "이후 MERGED_<target>_<start>_<end> workspace를 새로 생성합니다."
        ))

        grp = QGroupBox("입력 result 폴더")
        grp_layout = QVBoxLayout(grp)
        self.folder_list = QListWidget()
        self.folder_list.setMinimumHeight(180)
        grp_layout.addWidget(self.folder_list)

        btn_row = QHBoxLayout()
        btn_add = QPushButton("폴더 추가")
        style_button(btn_add, height=Tokens.H_BUTTON)
        btn_add.clicked.connect(self._on_add_folder)
        btn_remove = QPushButton("선택 제거")
        style_button(btn_remove, height=Tokens.H_BUTTON)
        btn_remove.clicked.connect(self._on_remove_folder)
        btn_up = QPushButton("위로")
        style_button(btn_up, height=Tokens.H_BUTTON)
        btn_up.clicked.connect(lambda: self._move_selected_folder(-1))
        btn_down = QPushButton("아래로")
        style_button(btn_down, height=Tokens.H_BUTTON)
        btn_down.clicked.connect(lambda: self._move_selected_folder(+1))
        btn_row.addWidget(btn_add)
        btn_row.addWidget(btn_remove)
        btn_row.addWidget(btn_up)
        btn_row.addWidget(btn_down)
        btn_row.addStretch()
        grp_layout.addLayout(btn_row)
        layout.addWidget(grp)

        out_grp = QGroupBox("출력 merged result 폴더")
        out_outer = QVBoxLayout(out_grp)
        out_layout = QHBoxLayout()
        self.output_dir_edit = QLineEdit()
        btn_browse_out = QPushButton("찾기...")
        style_button(btn_browse_out, height=Tokens.H_BUTTON)
        btn_browse_out.clicked.connect(self._browse_output_dir)
        out_layout.addWidget(self.output_dir_edit, 1)
        out_layout.addWidget(btn_browse_out)
        out_outer.addLayout(out_layout)

        storage_row = QHBoxLayout()
        storage_row.setSpacing(Tokens.GAP)
        storage_row.addWidget(QLabel("측광 저장 방식:"))
        self.storage_mode_combo = QComboBox()
        self.storage_mode_combo.addItem("전체 복사 (독립 실행 가능)", STORAGE_FULL)
        self.storage_mode_combo.addItem("원본 참조 (디스크 절약)", STORAGE_REFERENCE)
        self.storage_mode_combo.setMinimumHeight(Tokens.H_BUTTON)
        self.storage_mode_combo.setToolTip(
            "전체 복사: 프레임별 측광 TSV를 merged workspace 안으로 복사합니다. "
            "입력 폴더를 옮기거나 지워도 동작합니다.\n"
            "원본 참조: 복사 대신 어느 입력 workspace에서 왔는지만 기록하고 "
            "읽을 때 원본을 봅니다. 디스크를 아끼지만 "
            "입력 폴더를 옮기면 깨집니다."
        )
        self.storage_mode_combo.currentIndexChanged.connect(
            lambda _i: self._invalidate_merger_state())
        storage_row.addWidget(self.storage_mode_combo, 1)
        out_outer.addLayout(storage_row)
        layout.addWidget(out_grp)

        info_grp = QGroupBox("폴더 스캔")
        info_layout = QVBoxLayout(info_grp)
        self.folder_info_table = QTableWidget(0, 10)
        self.folder_info_table.setHorizontalHeaderLabels(["폴더", "Label", "Type", "Start", "End", "Step 7", "Step 8", "Step 9", "필터", "상태"])
        self.folder_info_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.folder_info_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.folder_info_table.setMinimumHeight(180)
        info_layout.addWidget(self.folder_info_table)
        btn_scan = QPushButton("폴더 스캔")
        style_button(btn_scan, height=Tokens.H_BUTTON)
        btn_scan.clicked.connect(self._scan_folders)
        info_layout.addWidget(btn_scan)
        layout.addWidget(info_grp)
        layout.addStretch()
        return page

    def _make_step2(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)

        layout.addWidget(self._make_info_label(
            "base 폴더의 selection / master catalog를 기준으로 canonical merged ID를 만듭니다.\n"
            "매칭 우선순위는 Gaia source_id → 기존 canonical source_id → 좌표 근접 매칭입니다."
        ))

        top = QHBoxLayout()
        top.addWidget(QLabel("Position match radius (arcsec):"))
        self.match_radius_combo = QComboBox()
        for value in ("1.0", "1.5", "2.0", "3.0", "5.0"):
            self.match_radius_combo.addItem(value)
        self.match_radius_combo.setCurrentText("2.0")
        top.addWidget(self.match_radius_combo)
        top.addStretch()
        self.btn_match = QPushButton("ID 매칭 실행")
        style_button(self.btn_match, "primary", height=Tokens.H_BUTTON)
        self.btn_match.clicked.connect(self._run_id_match)
        top.addWidget(self.btn_match)
        layout.addLayout(top)

        self.match_status_label = self._make_status_label("매칭 결과: 아직 실행 안 됨")
        layout.addWidget(self.match_status_label)

        self.match_table = QTableWidget(0, 7)
        self.match_table.setHorizontalHeaderLabels(
            ["폴더", "필터", "Exact SID", "Positional", "New", "총 rows", "상태"]
        )
        self.match_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.match_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.match_table.setMinimumHeight(240)
        layout.addWidget(self.match_table)

        self.match_log = QTextEdit()
        self.match_log.setReadOnly(True)
        self.match_log.setMinimumHeight(160)
        self.match_log.setObjectName("Log")     # themed mono log surface
        layout.addWidget(self.match_log)
        return page

    def _make_embedded_step_page(self, description: str, status_attr: str, host_attr: str) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.addWidget(self._make_info_label(description))
        status_label = self._make_status_label("Merged workspace 없음")
        layout.addWidget(status_label)
        setattr(self, status_attr, status_label)

        host = QWidget()
        host_layout = QVBoxLayout(host)
        host_layout.setContentsMargins(0, 0, 0, 0)
        host_layout.setSpacing(0)
        layout.addWidget(host, 1)
        setattr(self, host_attr, host_layout)
        return page

    def _make_step3(self) -> QWidget:
        return self._make_embedded_step_page(
            "Merged workspace의 Step 8 Selection UI를 그대로 보여줍니다.\n"
            "여기서 target / comparison / check를 바로 수정하면 merged workspace에 즉시 저장됩니다.",
            "selection_status_label",
            "step9_host_layout",
        )

    def _make_step4(self) -> QWidget:
        return self._make_embedded_step_page(
            "Merged workspace의 Step 9 Light Curve Builder UI를 그대로 보여줍니다.",
            "step10_status_label",
            "step10_host_layout",
        )

    def _make_step5(self) -> QWidget:
        return self._make_embedded_step_page(
            "Merged workspace의 Step 10 Detrend UI를 그대로 보여줍니다.",
            "step11_status_label",
            "step11_host_layout",
        )

    def _make_step6(self) -> QWidget:
        return self._make_embedded_step_page(
            "Merged workspace의 Step 11 Period Analysis UI를 그대로 보여줍니다.",
            "step12_status_label",
            "step12_host_layout",
        )

    # ───────────────────────── folder scan ─────────────────────────

    def _refresh_folder_list(self):
        self.folder_list.clear()
        for i, p in enumerate(self.folders):
            label = f"[BASE] {p}" if i == 0 else str(p)
            item = QListWidgetItem(label)
            if i == 0:
                item.setForeground(QColor(Tokens.ACCENT_TEXT))
            self.folder_list.addItem(item)

    def _invalidate_merger_state(self):
        self.folder_scan_rows = []
        self.match_summary_rows = []
        self.match_records = []
        self.merged_catalogs = {}
        self.local_id_maps = {}
        self.base_selection_by_filter = {}
        self.selection_target_by_filter = {}
        self.selection_comp_by_filter = {}
        self.selection_check_by_filter = {}
        self.merged_result_dir = None
        self.merged_runtime_params = None
        self.merged_runtime_project_state = None
        self.merged_runtime_file_manager = None
        self._clear_embedded_step_windows()
        if hasattr(self, "folder_info_table"):
            self.folder_info_table.setRowCount(0)
        if hasattr(self, "match_table"):
            self.match_table.setRowCount(0)
        if hasattr(self, "match_status_label"):
            self.match_status_label.setText("매칭 결과: 아직 실행 안 됨")
        if hasattr(self, "match_log"):
            self.match_log.clear()
        self._refresh_step_host_status_labels()

    def _detach_embedded_window(self, attr_name: str, host_layout_attr: str):
        window = getattr(self, attr_name, None)
        if window is None:
            return
        host_layout = getattr(self, host_layout_attr, None)
        if host_layout is not None:
            host_layout.removeWidget(window)
        window.setParent(None)
        window.deleteLater()
        setattr(self, attr_name, None)

    def _clear_embedded_step_windows(self):
        self._detach_embedded_window("step9_embedded_window", "step9_host_layout")
        self._detach_embedded_window("step10_embedded_window", "step10_host_layout")
        self._detach_embedded_window("step11_embedded_window", "step11_host_layout")
        self._detach_embedded_window("step12_embedded_window", "step12_host_layout")

    def _prepare_embedded_step_window(self, window):
        window.setParent(self._page_container)
        window.setWindowFlags(Qt.Widget)
        if hasattr(window, "title_label"):
            window.title_label.hide()
        for btn_name in ("btn_previous", "btn_complete", "btn_next"):
            btn = getattr(window, btn_name, None)
            if btn is not None:
                btn.hide()
        window.show()
        return window

    def _ensure_step9_embedded_window(self):
        if self.step9_embedded_window is not None:
            return self.step9_embedded_window
        if self.merged_result_dir is None:
            return None
        from apex.gui.workflow.lc.step8_target_selection import TargetComparisonSelectionWindow

        window = TargetComparisonSelectionWindow(
            self.merged_runtime_params,
            self.merged_runtime_file_manager,
            self.merged_runtime_project_state,
            self,
        )
        window = self._prepare_embedded_step_window(window)
        self.step9_host_layout.addWidget(window)
        self.step9_embedded_window = window
        return window

    def _ensure_step10_embedded_window(self):
        if self.step10_embedded_window is not None:
            return self.step10_embedded_window
        if self.merged_result_dir is None:
            return None
        from apex.gui.workflow.lc.step9_lightcurve_builder import LightCurveBuilderWindow

        window = LightCurveBuilderWindow(
            self.merged_runtime_params,
            self.merged_runtime_file_manager,
            self.merged_runtime_project_state,
            self,
        )
        window = self._prepare_embedded_step_window(window)
        self.step10_host_layout.addWidget(window)
        self.step10_embedded_window = window
        return window

    def _ensure_step11_embedded_window(self):
        if self.step11_embedded_window is not None:
            return self.step11_embedded_window
        if self.merged_result_dir is None:
            return None
        from apex.gui.workflow.lc.step10_detrend_merge import DetrendNightMergeWindow

        window = DetrendNightMergeWindow(
            self.merged_runtime_params,
            self.merged_runtime_file_manager,
            self.merged_runtime_project_state,
            self,
        )
        window = self._prepare_embedded_step_window(window)
        self.step11_host_layout.addWidget(window)
        self.step11_embedded_window = window
        return window

    def _ensure_step12_embedded_window(self):
        if self.step12_embedded_window is not None:
            return self.step12_embedded_window
        if self.merged_result_dir is None:
            return None
        from apex.gui.workflow.lc.step11_period_analysis import PeriodAnalysisWindow

        window = PeriodAnalysisWindow(
            self.merged_runtime_params,
            self.merged_runtime_file_manager,
            self.merged_runtime_project_state,
            self,
        )
        window = self._prepare_embedded_step_window(window)
        self.step12_host_layout.addWidget(window)
        self.step12_embedded_window = window
        return window

    def _refresh_embedded_step_windows_for_current_workspace(self):
        self._clear_embedded_step_windows()
        if self.merged_result_dir is None:
            return
        if getattr(self, "_current_step", 0) >= 2:
            self._sync_embedded_step_page(self._current_step)

    def _sync_embedded_step_page(self, idx: int):
        if self.merged_result_dir is None:
            return
        if idx == 2:
            window = self._ensure_step9_embedded_window()
            if window is not None:
                window.check_step7_status()
        elif idx == 3:
            window = self._ensure_step10_embedded_window()
            if window is not None:
                window._auto_load_ids()
        elif idx == 4:
            window = self._ensure_step11_embedded_window()
            if window is not None:
                window.restore_state()
                window._auto_load_ids()
                window.load_raw_data(silent=True)
        elif idx == 5:
            window = self._ensure_step12_embedded_window()
            if window is not None:
                window._auto_load_target_id()
                window._scan_available_data()
                window._load_light_curve(silent=True)
        self._refresh_step_host_status_labels()

    def _refresh_step_host_status_labels(self):
        if self.merged_result_dir is None:
            self.selection_status_label.setText("Merged workspace 없음")
            self.step10_status_label.setText("Merged workspace 없음")
            self.step11_status_label.setText("Merged workspace 없음")
            self.step12_status_label.setText("Merged workspace 없음")
            return

        merged_dir = Path(self.merged_result_dir)
        status_text = f"Workspace: {merged_dir}"
        self.selection_status_label.setText(status_text)
        self.step10_status_label.setText(status_text)
        self.step11_status_label.setText(status_text)
        self.step12_status_label.setText(status_text)

    def _step9_signature(self, result_dir: Path) -> tuple:
        s9 = step8_selection_dir(result_dir)
        return (
            str(result_dir.resolve()),
            max((p.stat().st_mtime for p in s9.glob("master_catalog_*.tsv")), default=None),
            max((p.stat().st_mtime for p in s9.glob("selection_*.json")), default=None),
        )

    def _get_cached_selection_payloads(self, result_dir: Path) -> dict[str, dict]:
        cache_key = str(result_dir.resolve())
        signature = self._step9_signature(result_dir)
        cached = self._selection_payload_cache.get(cache_key)
        if cached and cached[0] == signature:
            return cached[1]
        payloads = _load_selection_payloads(result_dir)
        self._selection_payload_cache[cache_key] = (signature, payloads)
        return payloads

    def _get_cached_master_catalogs(self, result_dir: Path) -> dict[str, pd.DataFrame]:
        cache_key = str(result_dir.resolve())
        signature = self._step9_signature(result_dir)
        cached = self._catalog_cache.get(cache_key)
        if cached and cached[0] == signature:
            return cached[1]
        catalogs = _load_master_catalogs_by_filter(result_dir)
        self._catalog_cache[cache_key] = (signature, catalogs)
        return catalogs

    def _id_match_signature(self) -> tuple:
        return (
            tuple(str(folder.resolve()) for folder in self.folders),
            tuple(self._step9_signature(folder) for folder in self.folders),
            float(self.match_radius_combo.currentText()),
        )

    def _refresh_output_dir_default(self, force: bool = False):
        if not self.folders:
            self.output_dir_edit.clear()
            return
        new_default = build_merged_workspace_dir(self.folders)
        current_text = self.output_dir_edit.text().strip()
        if force or not current_text:
            self.output_dir_edit.setText(str(new_default))
            return
        try:
            current_path = Path(current_text)
        except Exception:
            current_path = None
        old_default = _default_output_dir(self.folders[0])
        if current_path is not None and current_path == old_default:
            self.output_dir_edit.setText(str(new_default))

    def _on_add_folder(self):
        start_dir = self.folders[0].parent if self.folders else Path(self.params.P.result_dir).parent
        folder = QFileDialog.getExistingDirectory(self, "result 폴더 선택", str(start_dir))
        if not folder:
            return
        p = Path(folder)
        if any(existing.resolve() == p.resolve() for existing in self.folders):
            return
        self.folders.append(p)
        self._invalidate_merger_state()
        self._refresh_folder_list()
        self._refresh_output_dir_default(force=True)

    def _on_remove_folder(self):
        row = self.folder_list.currentRow()
        if row < 0:
            return
        self.folders.pop(row)
        self._invalidate_merger_state()
        self._refresh_folder_list()
        self._refresh_output_dir_default(force=True)

    def _move_selected_folder(self, delta: int):
        row = self.folder_list.currentRow()
        if row < 0:
            return
        new_row = row + delta
        if new_row < 0 or new_row >= len(self.folders):
            return
        self.folders[row], self.folders[new_row] = self.folders[new_row], self.folders[row]
        self._invalidate_merger_state()
        self._refresh_folder_list()
        self.folder_list.setCurrentRow(new_row)
        self._refresh_output_dir_default(force=True)

    def _browse_output_dir(self):
        path = QFileDialog.getExistingDirectory(
            self,
            "MERGED workspace 폴더 선택",
            str(Path(self.output_dir_edit.text()).parent if self.output_dir_edit.text().strip() else build_merged_workspace_dir(self.folders).parent),
        )
        if path:
            current_name = build_merged_workspace_dir(self.folders).name if self.folders else "MERGED_workspace"
            self.output_dir_edit.setText(str(Path(path) / current_name))

    def _scan_folders(self):
        self.folder_scan_rows = []
        self.folder_info_table.setRowCount(0)
        for folder in self.folders:
            cache_key = str(folder.resolve())
            signature = workspace_scan_signature(folder)
            cached = self._workspace_scan_cache.get(cache_key)
            if cached and cached[0] == signature:
                row_info = dict(cached[1])
            else:
                row_info = scan_merge_input_workspace(folder)
                self._workspace_scan_cache[cache_key] = (signature, dict(row_info))
            self.folder_scan_rows.append(row_info)

            row = self.folder_info_table.rowCount()
            self.folder_info_table.insertRow(row)
            self.folder_info_table.setItem(row, 0, QTableWidgetItem(folder.name))
            self.folder_info_table.setItem(row, 1, QTableWidgetItem(str(row_info["label"])))
            self.folder_info_table.setItem(row, 2, QTableWidgetItem(str(row_info["run_type"])))
            self.folder_info_table.setItem(row, 3, QTableWidgetItem(str(row_info["date_start"])))
            self.folder_info_table.setItem(row, 4, QTableWidgetItem(str(row_info["date_end"])))
            for col_idx, key in enumerate(("has_step7", "has_step8", "has_step9"), start=5):
                ok = bool(row_info[key])
                item = QTableWidgetItem("OK" if ok else "없음")
                item.setForeground(QColor(Tokens.OK if ok else Tokens.ERROR))
                self.folder_info_table.setItem(row, col_idx, item)
            filters = list(row_info.get("filters") or [])
            merge_ready = bool(row_info.get("merge_ready"))
            self.folder_info_table.setItem(row, 8, QTableWidgetItem(", ".join(filters) if filters else "—"))
            # Step 9 is optional — the merged workspace rebuilds the light
            # curves — so its absence is a note, not a blocker.
            if merge_ready and not row_info.get("has_step9"):
                status_item = QTableWidgetItem("사용 가능 (Step 9 없음)")
                status_item.setForeground(QColor(Tokens.WARN))
                status_item.setToolTip(
                    "Step 9 광곡선은 merge 에 필요하지 않습니다. "
                    "merged workspace 에서 다시 만듭니다."
                )
            else:
                status_item = QTableWidgetItem("사용 가능" if merge_ready else "입력 부족")
                status_item.setForeground(QColor(Tokens.OK if merge_ready else Tokens.ERROR))
            self.folder_info_table.setItem(row, 9, status_item)

    def _validate_selected_workspaces(self) -> tuple[bool, str]:
        if len(self.folders) < 2:
            return False, "Merge하려면 최소 2개의 RESULT/MERGED workspace를 선택하세요."
        if not self.folder_scan_rows:
            self._scan_folders()
        invalid_rows = [row for row in self.folder_scan_rows if not row.get("merge_ready")]
        if invalid_rows:
            missing = ", ".join(
                f"{Path(row['folder']).name}("
                + "+".join(
                    step for step, key in (("Step 7", "has_step7"), ("Step 8", "has_step8"))
                    if not row.get(key)
                ) + " 없음)"
                for row in invalid_rows
            )
            return False, ("Step 7(forced photometry)과 Step 8(마스터 카탈로그·선택)이 "
                           f"있는 workspace만 머저할 수 있습니다.\n부족: {missing}")
        labels = []
        seen = set()
        for row in self.folder_scan_rows:
            label = str(row.get("label") or "").strip()
            key = label.lower()
            if label and key not in seen:
                labels.append(label)
                seen.add(key)
        if len(labels) > 1:
            return False, "서로 다른 target/label의 RESULT를 한 번에 머저할 수 없습니다."
        return True, ""

    # ───────────────────────── Step 2: ID match ─────────────────────────

    def _run_id_match(self):
        if self._merge_busy:
            return
        self.match_log.clear()
        self.match_table.setRowCount(0)
        self.match_summary_rows = []
        self.match_records = []
        self.merged_catalogs = {}
        self.local_id_maps = {}
        self.base_selection_by_filter = {}

        if not self.folder_scan_rows:
            self._scan_folders()
        ok, msg = self._validate_selected_workspaces()
        if not ok:
            QMessageBox.warning(self, "ID Match", msg)
            return

        output_dir_text = self.output_dir_edit.text().strip()
        if not output_dir_text:
            QMessageBox.warning(self, "Merged Workspace", "출력 폴더를 지정하세요.")
            return
        out_dir = Path(output_dir_text)

        # Base folder first, then the rest: the base decides the target/check,
        # but a selection only present on a later night is no longer lost.
        selection_payloads = [self._get_cached_selection_payloads(folder)
                              for folder in self.folders]
        self.base_selection_by_filter = selection_payloads[0]
        catalogs_by_folder = {str(folder): self._get_cached_master_catalogs(folder) for folder in self.folders}
        self.folder_tags = {str(folder): _folder_tag(i, folder) for i, folder in enumerate(self.folders)}

        all_filters = sorted({
            flt for folder in self.folders
            for flt in catalogs_by_folder.get(str(folder), {}).keys()
        })
        if not all_filters:
            self.match_status_label.setText("매칭 실패: master_catalog 없음")
            self.match_log.append("[ERR] 어떤 폴더에서도 master_catalog_*.tsv 를 찾지 못했습니다.")
            return

        id_sig = self._id_match_signature()
        cached = self._id_match_cache.get(id_sig)
        if cached is not None:
            self.match_log.append("[MATCH] Using cached ID match result")

        # Heavy reconcile + materialize run on a worker thread so the UI stays
        # responsive; GUI wiring happens in _on_merge_finished on the main thread.
        self._pending_id_sig = id_sig
        worker = _MergeWorker(
            folders=list(self.folders),
            catalogs_by_folder=catalogs_by_folder,
            folder_tags=dict(self.folder_tags),
            match_radius=float(self.match_radius_combo.currentText()),
            base_selection_by_filter=selection_payloads,
            out_dir=out_dir,
            cached_reconcile=cached,
            storage_mode=self.storage_mode_combo.currentData(),
            parent=self,
        )
        worker.progress.connect(self.match_log.append)
        worker.finished_ok.connect(self._on_merge_finished)
        worker.failed.connect(self._on_merge_failed)
        self._merge_worker = worker
        self._set_merge_busy(True)
        self.match_status_label.setText("머징 중… (백그라운드 처리)")
        worker.start()

    def _set_merge_busy(self, busy: bool):
        self._merge_busy = busy
        self.btn_match.setEnabled(not busy)

    def _on_merge_finished(self, result: dict):
        cached = result["cached"]
        if self._pending_id_sig is not None:
            self._id_match_cache[self._pending_id_sig] = cached
        self.merged_catalogs = cached["canonical_by_filter"]
        self.local_id_maps = cached["local_id_maps"]
        self.match_summary_rows = cached["match_summary_rows"]
        self.match_records = cached["match_records"]
        self.selection_target_by_filter = result["sel_t"]
        self.selection_comp_by_filter = result["sel_c"]
        self.selection_check_by_filter = result["sel_k"]
        self._update_match_table()
        n_rows = sum(len(df) for df in self.merged_catalogs.values())
        self.match_status_label.setText(
            f"매칭 완료: filters={len(self.merged_catalogs)} canonical rows={n_rows} mapping rows={len(self.match_records)}"
        )

        self.merged_result_dir = result["out_dir"]
        self._build_merged_runtime_context(
            result["build_info"]["night_assignments"],
            result["build_info"]["path_map"],
        )
        self._refresh_embedded_step_windows_for_current_workspace()
        self._refresh_step_host_status_labels()
        self._set_merge_busy(False)
        self._merge_worker = None
        self._go_to_step(2)

    def _on_merge_failed(self, message: str):
        self._set_merge_busy(False)
        self._merge_worker = None
        self.match_status_label.setText("머징 실패")
        QMessageBox.warning(self, "Merged Workspace", message)

    def _update_match_table(self):
        self.match_table.setRowCount(0)
        for row_data in self.match_summary_rows:
            row = self.match_table.rowCount()
            self.match_table.insertRow(row)
            self.match_table.setItem(row, 0, QTableWidgetItem(str(row_data["folder"])))
            self.match_table.setItem(row, 1, QTableWidgetItem(str(row_data["filter"])))
            self.match_table.setItem(row, 2, QTableWidgetItem(str(row_data["exact"])))
            self.match_table.setItem(row, 3, QTableWidgetItem(str(row_data["pos"])))
            self.match_table.setItem(row, 4, QTableWidgetItem(str(row_data["new"])))
            self.match_table.setItem(row, 5, QTableWidgetItem(str(row_data["total"])))
            status_item = QTableWidgetItem(str(row_data["status"]))
            ok = row_data["status"] in ("OK", "base")
            status_item.setForeground(QColor(Tokens.OK if ok else Tokens.ERROR))
            self.match_table.setItem(row, 6, status_item)

    def _build_merged_runtime_context(self, night_assignments: dict[str, int], path_map: dict[str, str]):
        if self.merged_result_dir is None:
            return
        out_dir = Path(self.merged_result_dir)
        self.merged_runtime_params = _MergedParamsProxy(self.params, out_dir)
        self.merged_runtime_params.P.file_path_map = dict(path_map)
        idx = _read_step7_index(out_dir)
        filenames = idx["file"].astype(str).tolist() if not idx.empty and "file" in idx.columns else []
        self.merged_runtime_file_manager = _MergedFileManagerProxy(filenames, night_assignments, path_map)
        self.merged_runtime_project_state = ProjectState(out_dir)
        self.merged_runtime_project_state.state["project_name"] = out_dir.name
        self.merged_runtime_project_state.state["completed_steps"] = sorted(set(range(9)))
        self.merged_runtime_project_state.state["current_step"] = 8
        self.merged_runtime_project_state.store_step_data(
            "psf_photometry",
            {"skip_psf": True, "runtime_reason": "merged aperture workspace"},
        )
        self.merged_runtime_project_state.store_step_data(
            "file_selection",
            {
                "data_dir": str(out_dir.parent),
                "filename_prefix": getattr(self.params.P, "filename_prefix", ""),
            },
        )
        self.merged_runtime_project_state.save()

    def on_step_completed(self, step_index: int):
        self._refresh_step_host_status_labels()
        if 8 <= step_index <= 11:
            self._sync_embedded_step_page(step_index - 6)

    def open_step(self, step_index: int):
        step_map = {
            8: 2,
            9: 3,
            10: 4,
            11: 5,
        }
        idx = step_map.get(step_index)
        if idx is None:
            if step_index < 8:
                idx = 2
            else:
                return
        if idx >= 2 and self.merged_result_dir is None:
            QMessageBox.warning(self, "Merged Workflow", "Merged workspace를 먼저 생성하세요.")
            return
        self.show()
        self.raise_()
        self.activateWindow()
        self._go_to_step(idx)

    # ───────────────────────── navigation ─────────────────────────

    def _go_to_step(self, idx: int):
        for i, (page, btn) in enumerate(zip(self._pages, self._step_btns)):
            page.setVisible(i == idx)
            # The checked state alone marks the current step; the theme paints
            # it, so no per-button stylesheet is needed (and dark mode works).
            btn.setChecked(i == idx)
        self._current_step = idx
        self.btn_prev.setEnabled(idx > 0)
        self.btn_next.setEnabled(idx < len(self._pages) - 1)
        self._refresh_step_host_status_labels()
        if idx >= 2:
            self._sync_embedded_step_page(idx)

    def _prev_step(self):
        if self._current_step > 0:
            self._go_to_step(self._current_step - 1)

    def _next_step(self):
        if self._current_step == 0:
            ok, msg = self._validate_selected_workspaces()
            if not ok:
                QMessageBox.warning(self, "Merger", msg)
                return
        if self._current_step == 0 and not self.folders:
            QMessageBox.warning(self, "Merger", "폴더를 먼저 선택하세요.")
            return
        if self._current_step == 1 and not self.merged_catalogs:
            QMessageBox.warning(self, "Merger", "Step 2 ID 매칭을 먼저 실행하세요.")
            return
        if self._current_step == 2 and self.merged_result_dir is None:
            QMessageBox.warning(self, "Merger", "Merged workspace를 먼저 생성하세요.")
            return
        if self._current_step < len(self._pages) - 1:
            self._go_to_step(self._current_step + 1)

    # ───────────────────────── back ─────────────────────────

    def _go_back(self):
        self.hide()
        self.main_window.show()
        self.main_window.raise_()
        self.main_window.activateWindow()

    def closeEvent(self, event):
        self._go_back()
        event.ignore()
