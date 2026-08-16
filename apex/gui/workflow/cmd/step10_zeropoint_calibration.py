"""
Step 10: Zeropoint & Standardization
"""

from __future__ import annotations

import hashlib
import json
import time
import traceback
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.spatial import cKDTree

from apex.utils.constants import MAG_ERR_COEFF, MAD_TO_SIGMA

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.colors import LinearSegmentedColormap, Normalize
import matplotlib as mpl

from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QGroupBox, QMessageBox,
    QTextEdit, QDialog, QFormLayout, QDialogButtonBox, QDoubleSpinBox,
    QSpinBox, QCheckBox, QComboBox, QWidget, QTabWidget, QFileDialog, QLineEdit
)


from apex.analysis.cmd.zeropoint_runner import (  # noqa: F401  (re-exported)
    ZeropointCalibrationRunner,
    build_cmd_qc_summary,
    build_gaia_cmd_comparison,
    build_gaia_cmd_drift_table,
    build_gaia_cmd_snr_sweep,
    build_zp_qc_summary,
    draw_cmd_qc_overview,
    draw_gaia_cmd_comparison,
    draw_gaia_cmd_snr_sweep,
    draw_zp_qc_overview,
    export_cmd_qc_products,
    export_gaia_cmd_comparison_products,
    export_zp_qc_products,
    resolve_cmd_photometry_input,
    robust_weighted_polyfit,
    select_cmd_qc_axes,
    solve_standard_colors,
)
from apex.gui.workflow.step_window_base import StepWindowBase
from apex.analysis.light_curve.photometry_source_service import (
    resolve_lightcurve_photometry_source,
)
from apex.gui.theme import mono_note_style


def _set_label_role(label, prop: str, value) -> None:
    """Swap a theme role property at runtime and repolish so the QSS re-runs.

    Roles are static in QSS; a label whose meaning changes (ROI set vs unset)
    must clear the other property or both selectors stay live.
    """
    for name in ("role", "status", "banner"):
        label.setProperty(name, value if name == prop else None)
    style = label.style()
    style.unpolish(label)
    style.polish(label)


from apex.gui.workflow.run_control import RunControlBar
from apex.gui.workflow.log_panel import WorkflowLogWindow, WorkerStatusPanel, append_timestamped_log, show_raised
from apex.gui.workflow.ui_helpers import (
    add_parameter_reset_button,
    build_scroll_param_dialog,
    create_collapsible_section,
    create_parameter_button,
    install_parameter_wheel_guard,
)
from apex.utils.astro_utils import normalize_filter_name
from apex.utils.common_helpers import format_cmd_title, photometric_system_label
from apex.utils.cmd_gaia_enrichment import (
    load_master_table as _load_master_table,
    merge_gaia_columns_from_catalog as _merge_gaia_columns_from_catalog,
)
from apex.utils.step_paths import (
    step2_cropped_dir,
    crop_is_active,
    step7_forced_phot_dir,
    step5_wcs_dir,
    tool_extinction_dir,
)
from apex.utils.step_paths_cmd import step8_psf_dir, step9_selection_dir, step10_zp_dir
from apex.utils.io_utils import parse_int64_series, read_ecsv_int64_source_id
from apex.utils.gaia_quality import (
    gaia_corrected_excess_factor,
    gaia_cstar_sigma,
    gaia_quality_mask,
    gaia_quality_report,
)
from apex.utils.qc_utils import filter_frame_df_by_qc, should_use_frame_quality_qc
from apex.utils.photometry_provenance import (
    build_photometry_provenance,
    collapse_provenance_values,
    format_photometry_provenance,
    summarize_photometry_table,
)


from apex.utils.gaia_transforms import (
    GAIA_TO_BAND       as _GAIA_TO_BAND,
    FILTER_COLOR_PREF  as _FILTER_COLOR_PREF,
    BAND_ALIASES       as _BAND_ALIASES,
    build_color_pairs  as _build_color_pairs,
    teff_from_color    as _teff_from_color,
    TEFF_COLOR_ANCHORS as _TEFF_COLOR_ANCHORS,
    filter_bands_from_columns as _filter_bands_from_columns,
)

_ZP_SIGNATURE_FILE = "zeropoint_signature.json"
_ZP_SIGNATURE_VERSION = 3
_ZP_SIGNATURE_PARAMS = (
    "match_tol_px",
    "min_master_gaia_matches",
    "cmd_snr_calib_min",
    "frame_zp_min_n",
    "cmd_apply_extinction",
    "cmd_extinction_mode",
    "zp_clip_sigma",
    "zp_fit_iters",
    "zp_slope_absmax",
    "gaia_snr_calib_min",
    "gaia_gi_min",
    "gaia_gi_max",
    "gaia_zp_slope_absmax",
    "gaia_color_slope_absmax",
    "min_snr_for_mag",
    "phot_ref_apcorr_min_keep",
    "phot_ref_require_apcorr_candidate",
    "phot_use_qc_pass_only",
    "ref_frame",
    "site_lat_deg",
    "site_lon_deg",
    "site_alt_m",
    "site_tz_offset_hours",
)




def resolve_cmd_photometry_provenance(result_dir: Path | str) -> dict[str, str]:
    """Describe the last CMD product, or the input that Step 10 would use."""
    root = Path(result_dir)
    output_dir = step10_zp_dir(root)
    for name in (
        "median_by_ID_filter_wide_cmd.csv",
        "median_by_ID_filter_wide.csv",
        "median_by_ID_filter_wide_raw.csv",
    ):
        path = output_dir / name
        if not path.exists() or path.stat().st_size <= 0:
            continue
        try:
            info = summarize_photometry_table(pd.read_csv(path, nrows=500))
            if info["source"] != "unknown":
                return info
        except Exception:
            continue

    source_info = resolve_lightcurve_photometry_source(root)
    index_path = Path(source_info["index_path"])
    if index_path.exists() and index_path.stat().st_size > 0:
        return build_photometry_provenance(
            source_info.get("source"),
            source_info.get("mag_column"),
            source_info.get("mag_error_column"),
        )

    index_path = next(
        (
            path
            for path in _cmd_photometry_index_candidates(root)
            if path.exists() and path.stat().st_size > 0
        ),
        None,
    )
    return build_photometry_provenance("aperture")
















# Minimum sigma-clipped calibrators before the quadratic color term is
# attempted; below this the linear fit is kept (a poorly constrained curvature
# does more harm than the ±0.02-0.03 mag it corrects on rich fields).
























































# The calculation lives in `apex.analysis.cmd.zeropoint_runner` so a script can
# run Step 10 without PyQt5. This adds the thread and the real signals; because
# `pyqtSignal` declarations are class attributes, `SignalHost` finds them and
# does not replace them, so the window and the headless runner emit through the
# very same code (2026-08-16).
class ZeropointCalibrationWorker(QThread, ZeropointCalibrationRunner):
    progress = pyqtSignal(int, int, str)
    log = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    # `QThread` comes first in the MRO and has its own empty `run`, which would
    # win and quietly do nothing. Bind the calculation explicitly.
    run = ZeropointCalibrationRunner.run

    def __init__(self, params, data_dir: Path, result_dir: Path, cache_dir: Path):
        QThread.__init__(self)
        ZeropointCalibrationRunner.__init__(self, params, data_dir, result_dir, cache_dir)
        # This is the whole adapter: the calculation announces on its channels,
        # and here each one becomes a Qt signal, so it crosses to the GUI thread
        # exactly as it did when the worker was a QThread subclass.
        self.on_progress.subscribe(self.progress.emit)
        self.on_log.subscribe(self.log.emit)
        self.on_finished.subscribe(self.finished.emit)
        self.on_error.subscribe(self.error.emit)


class CmdViewerWindow(QWidget):
    """Interactive CMD viewer (Qt)."""

    def __init__(self, df: pd.DataFrame, result_dir: Path, parent=None, embedded: bool = False, params=None):
        super().__init__(parent)
        self.df = self._with_calibrated_aliases(df)
        self.photometry_provenance = summarize_photometry_table(self.df)
        self.result_dir = Path(result_dir)
        self.params = params

        self.setWindowTitle("CMD Viewer")
        if embedded:
            self.setWindowFlags(Qt.Widget)
            # Embedded in Step 11, this is a *preferred* size, not a floor:
            # a 600 px minimum exceeded the room Step 11 had, so the host
            # scroll showed only ~150 px of the plot at a time. Keep it low
            # enough to fit a laptop screen and let the canvas expand.
            self.setMinimumSize(760, 420)
        else:
            self.setWindowFlag(Qt.Window, True)
            self.resize(1200, 900)
            self.setMinimumSize(1000, 720)

        # View mode is selected after available magnitude products are detected.
        self.view_mode = 0

        self.inst_bands = _filter_bands_from_columns(self.df.columns, "mag_inst_")
        self.std_bands  = _filter_bands_from_columns(self.df.columns, "mag_std_")
        std_value_cols = [f"mag_std_{band}" for band in self.std_bands if f"mag_std_{band}" in self.df.columns]
        self.has_std = bool(std_value_cols) and np.isfinite(self.df[std_value_cols].to_numpy(float)).any()

        self.all_bands = sorted(set(self.inst_bands) | set(self.std_bands))

        # X axis: adjacent color pairs only (e.g. B-V, V-R — standard CMD indices)
        # Y axis: scalar magnitudes only (CMD viewer convention)
        axis_bands = self.inst_bands or self.std_bands
        x_allowed = _build_color_pairs(axis_bands, adjacent_only=True)
        self.x_allowed         = x_allowed
        self.y_allowed_scalars = axis_bands  # already wavelength-sorted
        self.y_allowed_colors  = []

        self.x_pairs       = x_allowed
        self.y_scalar_opts = axis_bands
        self.y_color_pairs = []

        self.snr_cols = [c for c in self.df.columns if c.startswith("snr_")]
        self.has_snr = len(self.snr_cols) > 0

        self.has_gaia_inst = (
            {"gaia_G_inst", "gaia_BP_RP_inst"}.issubset(df.columns)
            and np.isfinite(df["gaia_G_inst"].to_numpy(float)).any()
            and np.isfinite(df["gaia_BP_RP_inst"].to_numpy(float)).any()
        )
        self.has_gaia_syn = (
            {"gaia_G_syn", "gaia_BP_RP_syn"}.issubset(df.columns)
            and np.isfinite(df["gaia_G_syn"].to_numpy(float)).any()
            and np.isfinite(df["gaia_BP_RP_syn"].to_numpy(float)).any()
        )
        # Gaia CMD is a diagnostic-only view.  Keep Gaia columns available for
        # membership and click details, but do not add a third CMD panel.
        self.gaia_mode = None

        self.teff_vmin = 2400.0
        self.teff_vmax = 40000.0
        self.ob_norm = Normalize(vmin=self.teff_vmin, vmax=self.teff_vmax, clip=True)

        anchors = [
            (2400, "#E53935"),
            (3200, "#FF6A3D"),
            (4500, "#FFB84D"),
            (5800, "#FFE36A"),
            (6500, "#FFF6C7"),
            (8000, "#FFFFFF"),
            (10000, "#FFFFFF"),
            (20000, "#2D5BFF"),
            (40000, "#7A3CFF"),
        ]
        anchors = sorted(anchors, key=lambda x: x[0])
        pos = [(t - self.teff_vmin) / (self.teff_vmax - self.teff_vmin) for t, _ in anchors]
        pos[0] = 0.0
        pos[-1] = 1.0

        self.ob_cmap = LinearSegmentedColormap.from_list(
            "obafgkm_like",
            list(zip(pos, [c for _, c in anchors])),
            N=256
        )
        self.ob_cmap.set_bad("#777777")

        self.color_anchors = _TEFF_COLOR_ANCHORS

        # Determine available views
        self.available_views = []
        if self.inst_bands:
            self.available_views.append("inst")
        if self.has_std:
            self.available_views.append("std")
        if self.gaia_mode is not None:
            self.available_views.append("gaia")
        if len(self.available_views) > 1:
            self.available_views.append("all")
        if not self.available_views:
            self.available_views = ["inst"]
        self.view_mode = self.available_views.index("std") if "std" in self.available_views else 0

        self._plot_cache = {}
        self.last_pick_info = None
        self.pick_log = []
        self._membership_prob = None
        self._membership_source = "none"
        self._membership_note = ""
        self._membership_ready = False
        self._parallax_range_initialized = False
        self._roi_data: dict | None = None
        self._build_ui()
        self._update_view_label()
        self._load_roi()
        self._initialize_parallax_range(force=True)
        self._build_figure()
        self._redraw()
        self.setFocusPolicy(Qt.StrongFocus)

    @staticmethod
    def _with_calibrated_aliases(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        for band in _filter_bands_from_columns(out.columns, "mag_cal_"):
            cal_col = f"mag_cal_{band}"
            std_col = f"mag_std_{band}"
            if std_col not in out.columns and cal_col in out.columns:
                out[std_col] = out[cal_col]
            cal_err_col = f"mag_cal_err_{band}"
            std_err_col = f"mag_std_err_{band}"
            if std_err_col not in out.columns and cal_err_col in out.columns:
                out[std_err_col] = out[cal_err_col]
        return out

    def _update_view_label(self) -> None:
        view_name = self.available_views[self.view_mode] if self.available_views else "inst"
        view_labels = {"inst": "Instrumental", "std": "Calibrated", "gaia": "Gaia", "all": "All CMDs"}
        if hasattr(self, "view_label"):
            self.view_label.setText(f"View: {view_labels.get(view_name, view_name)}")

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # Top controls split across two rows so labels don't get clipped
        # in narrow windows.  Row 1 = selection axes & filters; Row 2 =
        # export buttons + current view label.
        controls = QHBoxLayout()
        controls.addWidget(QLabel("X(color):"))
        self.x_combo = QComboBox()
        self.x_combo.addItems([f"{a}-{b}" for (a, b) in self.x_pairs] or ["(none)"])
        controls.addWidget(self.x_combo)

        controls.addWidget(QLabel("Y:"))
        y_opts = self.y_scalar_opts + [f"{a}-{b}" for (a, b) in self.y_color_pairs]
        self.y_combo = QComboBox()
        self.y_combo.addItems(y_opts or ["(none)"])
        controls.addWidget(self.y_combo)

        controls.addWidget(QLabel("SNR >="))
        self.snr_spin = QSpinBox()
        self.snr_spin.setRange(0, 100)
        self.snr_spin.setValue(20)
        controls.addWidget(self.snr_spin)

        self.invert_y = QCheckBox("Invert Y")
        self.invert_y.setChecked(True)
        controls.addWidget(self.invert_y)

        # Optional extra ZP nudge for the Instrumental view. Instrumental
        # magnitudes already carry the IRAF Z=25 convention baked in at Step 7
        # (mag_inst = 25 - 2.5*log10(flux_e/exptime)), so they read in the
        # usual positive range without any shift — this control defaults to 0
        # and is only for manual fine-tuning. Colors (X = a-b) are unaffected
        # because a constant ZP cancels in a difference.
        self.manual_zp_check = QCheckBox("Manual ZP")
        self.manual_zp_check.setToolTip(
            "Add an extra constant zeropoint to Instrumental magnitudes for display only.\n"
            "mag_inst already includes the IRAF Z=25 convention, so leave at 0 normally.\n"
            "Colors are unchanged."
        )
        controls.addWidget(self.manual_zp_check)
        self.manual_zp_spin = QDoubleSpinBox()
        self.manual_zp_spin.setRange(0.0, 50.0)
        self.manual_zp_spin.setDecimals(3)
        self.manual_zp_spin.setSingleStep(0.1)
        self.manual_zp_spin.setValue(0.0)
        self.manual_zp_spin.setToolTip("Extra Instrumental-view zeropoint added to Y (display only).")
        controls.addWidget(self.manual_zp_spin)

        controls.addSpacing(8)
        controls.addWidget(QLabel("Membership:"))
        self.member_mode_combo = QComboBox()
        self.member_mode_combo.addItems([
            "Off",
            "Loose (P>=0.30)",
            "Normal (P>=0.50)",
            "Strict (P>=0.80)",
        ])
        mode_raw = "off"
        if self.params is not None and hasattr(self.params, "P"):
            mode_raw = str(getattr(self.params.P, "cmd_membership_mode", "off")).strip().lower()
        mode_to_idx = {"off": 0, "loose": 1, "normal": 2, "strict": 3}
        self.member_mode_combo.setCurrentIndex(mode_to_idx.get(mode_raw, 0))
        controls.addWidget(self.member_mode_combo)

        self.member_compare = QCheckBox("Compare")
        cmp_default = True
        if self.params is not None and hasattr(self.params, "P"):
            cmp_default = bool(getattr(self.params.P, "cmd_membership_compare", True))
        self.member_compare.setChecked(cmp_default)
        controls.addWidget(self.member_compare)
        controls.addStretch()
        layout.addLayout(controls)

        controls_row2 = QHBoxLayout()
        self.btn_save_membership = QPushButton("Save Pmem CSV")
        controls_row2.addWidget(self.btn_save_membership)

        self.btn_reset_filters = QPushButton("Reset View Filters")
        self.btn_reset_filters.setToolTip("Restore CMD viewer filters to the project defaults.")
        controls_row2.addWidget(self.btn_reset_filters)

        self.save_btn = QPushButton("Save PNG")
        controls_row2.addWidget(self.save_btn)

        self.xerr_check = QCheckBox("X err")
        self.xerr_check.setToolTip("Show color-index error bars for foreground CMD points.")
        self.xerr_check.setChecked(False)
        controls_row2.addWidget(self.xerr_check)

        self.yerr_check = QCheckBox("Y err")
        self.yerr_check.setToolTip("Show magnitude error bars for foreground CMD points.")
        self.yerr_check.setChecked(False)
        controls_row2.addWidget(self.yerr_check)

        controls_row2.addStretch()
        self.photometry_source_label = QLabel(
            format_photometry_provenance(self.photometry_provenance)
        )
        self.photometry_source_label.setProperty("role", "caption")
        _f = self.photometry_source_label.font(); _f.setBold(True)
        self.photometry_source_label.setFont(_f)
        controls_row2.addWidget(self.photometry_source_label)

        self.view_label = QLabel("View: Instrumental")
        self.view_label.setProperty("role", "info")
        controls_row2.addWidget(self.view_label)
        layout.addLayout(controls_row2)

        self.info_text = QTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setFixedHeight(90)
        self.info_text.setObjectName("Log")     # themed mono surface
        layout.addWidget(self.info_text)

        # 10×6 instead of 12×6 — the wider aspect made the CMD look stretched
        # on 1400px+ windows.
        self.figure = Figure(figsize=(10, 6), dpi=100)
        self.figure.subplots_adjust(bottom=0.14)
        self.canvas = FigureCanvas(self.figure)
        # Canvas is the only stretch=1 widget below; let it absorb ALL
        # spare vertical space.  (A maxHeight here caused leftover space
        # to be distributed as ugly gaps between the control rows in the
        # standalone CMD+Isochrone tool.)
        self.canvas.setMinimumSize(640, 300)

        # Prev/Next View are mounted next to the matplotlib navigation
        # toolbar (above the canvas) so they cannot visually collide with
        # the X-axis label below.
        self.btn_prev_view = QPushButton("Prev View")
        self.btn_next_view = QPushButton("Next View")
        self.btn_prev_view.setFixedHeight(28)
        self.btn_next_view.setFixedHeight(28)
        self.toolbar = NavigationToolbar(self.canvas, self)
        toolbar_row = QHBoxLayout()
        toolbar_row.setContentsMargins(0, 0, 0, 0)
        toolbar_row.addWidget(self.toolbar)
        toolbar_row.addStretch()
        toolbar_row.addWidget(self.btn_prev_view)
        toolbar_row.addWidget(self.btn_next_view)
        layout.addLayout(toolbar_row)
        layout.addWidget(self.canvas, stretch=1)

        self.x_combo.currentTextChanged.connect(self._redraw)
        self.y_combo.currentTextChanged.connect(self._redraw)
        self.snr_spin.valueChanged.connect(self._redraw)
        self.invert_y.stateChanged.connect(self._redraw)
        self.manual_zp_check.stateChanged.connect(self._redraw)
        self.manual_zp_spin.valueChanged.connect(self._redraw)
        self.xerr_check.stateChanged.connect(self._redraw)
        self.yerr_check.stateChanged.connect(self._redraw)
        self.member_mode_combo.currentIndexChanged.connect(self._on_membership_ui_changed)
        self.member_compare.stateChanged.connect(self._on_membership_ui_changed)
        self.btn_save_membership.clicked.connect(self._save_membership_csv)
        self.btn_reset_filters.clicked.connect(self._reset_view_filters)
        self.save_btn.clicked.connect(self._save_png)
        self.btn_prev_view.clicked.connect(lambda: self._switch_view(-1))
        self.btn_next_view.clicked.connect(lambda: self._switch_view(1))
        self.canvas.mpl_connect("button_press_event", self._on_plot_click)

        controls2 = QHBoxLayout()
        self.plx_check = QCheckBox("Parallax filter")
        self.plx_check.setChecked(False)
        controls2.addWidget(self.plx_check)
        controls2.addWidget(QLabel("min:"))
        self.plx_min_spin = QDoubleSpinBox()
        self.plx_min_spin.setRange(-5.0, 20.0)
        self.plx_min_spin.setDecimals(3)
        self.plx_min_spin.setSingleStep(0.05)
        self.plx_min_spin.setValue(-0.5)
        self.plx_min_spin.setSuffix(" mas")
        controls2.addWidget(self.plx_min_spin)
        controls2.addWidget(QLabel("max:"))
        self.plx_max_spin = QDoubleSpinBox()
        self.plx_max_spin.setRange(-5.0, 20.0)
        self.plx_max_spin.setDecimals(3)
        self.plx_max_spin.setSingleStep(0.05)
        self.plx_max_spin.setValue(0.5)
        self.plx_max_spin.setSuffix(" mas")
        controls2.addWidget(self.plx_max_spin)
        controls2.addSpacing(16)
        self.roi_check = QCheckBox("ROI filter")
        self.roi_check.setChecked(False)
        self.roi_check.setToolTip("Filter CMD sources by the spatial ROI circle set in Step 9.\nDoes not affect ZP calibration.")
        controls2.addWidget(self.roi_check)
        self.roi_info_label = QLabel("(no ROI)")
        self.roi_info_label.setProperty("role", "caption")
        controls2.addWidget(self.roi_info_label)
        self.btn_reload_roi = QPushButton("Reload")
        # sizeHint, not a hard 56 px: with the theme's button padding that
        # width clipped the label to "elo:".
        self.btn_reload_roi.setFixedWidth(self.btn_reload_roi.sizeHint().width())
        self.btn_reload_roi.setToolTip("Re-read cmd_roi.json from Step 9 output directory")
        controls2.addWidget(self.btn_reload_roi)
        controls2.addSpacing(16)

        # Gaia astrometric/photometric quality, as *display* filters. These are
        # the same two cuts step10 can apply when choosing zero-point
        # calibrators, offered here so the CMD can be inspected with and
        # without them — in a globular the C* cut alone removes about half the
        # Gaia-matched sources, and that is worth seeing rather than assuming.
        # Filtering here never touches the calibration; it only changes what
        # is drawn.
        self.ruwe_check = QCheckBox("RUWE ≤")
        self.ruwe_check.setChecked(False)
        self.ruwe_check.setToolTip(
            "Hide sources whose Gaia astrometric fit is poor (RUWE above the\n"
            "threshold) — typically unresolved binaries and blends.\n"
            "Display only; does not affect ZP calibration.")
        controls2.addWidget(self.ruwe_check)
        self.ruwe_spin = QDoubleSpinBox()
        self.ruwe_spin.setRange(1.0, 10.0)
        self.ruwe_spin.setDecimals(2)
        self.ruwe_spin.setSingleStep(0.1)
        self.ruwe_spin.setValue(1.4)
        self.ruwe_spin.setFixedWidth(self.ruwe_spin.sizeHint().width())
        controls2.addWidget(self.ruwe_spin)

        self.cstar_check = QCheckBox("C* ≤")
        self.cstar_check.setChecked(False)
        self.cstar_check.setToolTip(
            "Hide sources whose BP/RP flux excess is inconsistent with G by\n"
            "more than N sigma (Riello+2021) — BP/RP window contamination in\n"
            "crowded fields. Needs phot_bp_rp_excess_factor in the CMD table.\n"
            "Display only; does not affect ZP calibration.")
        controls2.addWidget(self.cstar_check)
        self.cstar_spin = QDoubleSpinBox()
        self.cstar_spin.setRange(1.0, 10.0)
        self.cstar_spin.setDecimals(1)
        self.cstar_spin.setSingleStep(0.5)
        self.cstar_spin.setValue(3.0)
        self.cstar_spin.setSuffix(" σ")
        self.cstar_spin.setFixedWidth(self.cstar_spin.sizeHint().width())
        controls2.addWidget(self.cstar_spin)

        self.quality_info_label = QLabel("")
        self.quality_info_label.setProperty("role", "caption")
        controls2.addWidget(self.quality_info_label)

        controls2.addStretch()
        layout.addLayout(controls2)

        self.plx_check.stateChanged.connect(self._on_plx_filter_changed)
        self.plx_min_spin.valueChanged.connect(self._redraw)
        self.plx_max_spin.valueChanged.connect(self._redraw)
        self.roi_check.stateChanged.connect(self._redraw)
        self.ruwe_check.stateChanged.connect(self._on_quality_filter_changed)
        self.ruwe_spin.valueChanged.connect(self._redraw)
        self.cstar_check.stateChanged.connect(self._on_quality_filter_changed)
        self.cstar_spin.valueChanged.connect(self._redraw)
        self.btn_reload_roi.clicked.connect(self._on_reload_roi)
        install_parameter_wheel_guard(self)

    def _reset_view_filters(self):
        for widget in (
            self.snr_spin,
            self.invert_y,
            self.member_mode_combo,
            self.member_compare,
            self.xerr_check,
            self.yerr_check,
            self.plx_check,
            self.plx_min_spin,
            self.plx_max_spin,
            self.roi_check,
        ):
            widget.blockSignals(True)

        self.snr_spin.setValue(20)
        self.invert_y.setChecked(True)
        self.member_mode_combo.setCurrentIndex(2)  # Normal (P>=0.50)
        self.member_compare.setChecked(True)
        self.xerr_check.setChecked(False)
        self.yerr_check.setChecked(False)
        self.plx_check.setChecked(False)
        self.plx_min_spin.setValue(-0.5)
        self.plx_max_spin.setValue(0.5)
        self._initialize_parallax_range(force=True)
        self.roi_check.setChecked(False)

        for widget in (
            self.snr_spin,
            self.invert_y,
            self.member_mode_combo,
            self.member_compare,
            self.xerr_check,
            self.yerr_check,
            self.plx_check,
            self.plx_min_spin,
            self.plx_max_spin,
            self.roi_check,
        ):
            widget.blockSignals(False)

        if self.params is not None and hasattr(self.params, "P"):
            self.params.P.cmd_membership_mode = "normal"
            self.params.P.cmd_membership_compare = True
            if hasattr(self.params, "save_toml"):
                try:
                    self.params.save_toml()
                except Exception:
                    pass
        self._redraw()

    def _load_roi(self):
        """Load cmd_roi.json from step8 output directory and update UI."""
        roi_path = step9_selection_dir(self.result_dir) / "cmd_roi.json"
        try:
            if roi_path.exists():
                self._roi_data = json.loads(roi_path.read_text())
            else:
                self._roi_data = None
        except Exception:
            self._roi_data = None
        if hasattr(self, "roi_check"):
            self.roi_check.setEnabled(self._roi_data is not None)
        if hasattr(self, "roi_info_label"):
            if self._roi_data:
                ra = self._roi_data.get("ra_deg", 0.0)
                dec = self._roi_data.get("dec_deg", 0.0)
                r = self._roi_data.get("radius_arcsec", 0.0)
                self.roi_info_label.setText(f"RA={ra:.4f} Dec={dec:.4f}  r={r:.0f}\"")
                _set_label_role(self.roi_info_label, "status", "ok")
            else:
                self.roi_info_label.setText("(no ROI)")
                _set_label_role(self.roi_info_label, "role", "caption")

    def _on_reload_roi(self):
        self._load_roi()
        self._redraw()

    def _parallax_values(self):
        if "parallax" not in self.df.columns:
            self._ensure_membership_columns_from_master()
        if "parallax" not in self.df.columns:
            return None
        return pd.to_numeric(self.df["parallax"], errors="coerce").to_numpy(float)

    def _set_parallax_range(self, plx_min: float, plx_max: float):
        if not np.isfinite(plx_min) or not np.isfinite(plx_max):
            return
        if plx_min > plx_max:
            plx_min, plx_max = plx_max, plx_min
        lo = max(float(self.plx_min_spin.minimum()), float(plx_min))
        hi = min(float(self.plx_max_spin.maximum()), float(plx_max))
        if lo >= hi:
            pad = max(0.05, abs(float(plx_min)) * 0.05)
            lo = max(float(self.plx_min_spin.minimum()), float(plx_min) - pad)
            hi = min(float(self.plx_max_spin.maximum()), float(plx_max) + pad)
        for widget in (self.plx_min_spin, self.plx_max_spin):
            widget.blockSignals(True)
        self.plx_min_spin.setValue(lo)
        self.plx_max_spin.setValue(hi)
        for widget in (self.plx_min_spin, self.plx_max_spin):
            widget.blockSignals(False)

    def _auto_set_parallax_range(self, preferred_mask=None) -> bool:
        plx = self._parallax_values()
        if plx is None:
            return False
        finite = np.isfinite(plx)
        base = finite.copy()
        if preferred_mask is not None and len(preferred_mask) == len(plx):
            preferred = np.asarray(preferred_mask, bool) & finite
            if int(preferred.sum()) >= 10:
                base = preferred
        if int(base.sum()) == 0:
            return False

        vals = plx[base]
        center = float(np.nanmedian(vals))
        mad = float(np.nanmedian(np.abs(vals - center)))
        robust_sigma = MAD_TO_SIGMA * mad if np.isfinite(mad) and mad > 0 else 0.0
        half_width = max(0.5, 4.0 * robust_sigma)
        half_width = min(5.0, half_width)
        self._set_parallax_range(center - half_width, center + half_width)
        return True

    def _initialize_parallax_range(self, force: bool = False) -> bool:
        if self._parallax_range_initialized and not force:
            return False
        if self._auto_set_parallax_range():
            self._parallax_range_initialized = True
            return True
        return False

    def _roi_mask(self):
        """Returns boolean mask selecting sources inside the CMD ROI circle (sky coords), or None if disabled."""
        if not self.roi_check.isChecked() or self._roi_data is None:
            return None
        roi_ra = float(self._roi_data["ra_deg"])
        roi_dec = float(self._roi_data["dec_deg"])
        roi_r_arcsec = float(self._roi_data["radius_arcsec"])
        # Prefer RA/Dec angular distance (correct across frames)
        if "ra_deg" in self.df.columns and "dec_deg" in self.df.columns:
            ra = pd.to_numeric(self.df["ra_deg"], errors="coerce").to_numpy(float)
            dec = pd.to_numeric(self.df["dec_deg"], errors="coerce").to_numpy(float)
            valid = np.isfinite(ra) & np.isfinite(dec)
            # Small-field approximation (accurate to ~0.01% within 1 deg)
            cos_dec = np.cos(np.radians(roi_dec))
            d_ra = (ra - roi_ra) * cos_dec * 3600.0   # arcsec
            d_dec = (dec - roi_dec) * 3600.0           # arcsec
            return valid & (d_ra ** 2 + d_dec ** 2 <= roi_r_arcsec ** 2)
        return None

    def _on_plx_filter_changed(self):
        if self.plx_check.isChecked():
            plx = self._parallax_values()
            if plx is None:
                self.plx_check.blockSignals(True)
                self.plx_check.setChecked(False)
                self.plx_check.blockSignals(False)
                from PyQt5.QtWidgets import QMessageBox
                QMessageBox.warning(self, "Parallax Unavailable",
                    "parallax column not found in CMD data or master_catalog.\n"
                    "Rerun Step 6 (Master Catalog Build) and Step 10 (ZP Calibration).")
                return
            plx_min = float(self.plx_min_spin.value())
            plx_max = float(self.plx_max_spin.value())
            finite = np.isfinite(plx)
            selected = finite & (plx >= plx_min) & (plx <= plx_max)
            if int(finite.sum()) > 0 and int(selected.sum()) < 5:
                member_mask, _, _ = self._current_membership_mask()
                self._auto_set_parallax_range(member_mask)
        self._redraw()

    def _parallax_mask(self):
        """Returns boolean mask for parallax range filter, or None if disabled/unavailable."""
        if not self.plx_check.isChecked():
            return None
        plx = self._parallax_values()
        if plx is None:
            return None
        plx_min = float(self.plx_min_spin.value())
        plx_max = float(self.plx_max_spin.value())
        mask = np.isfinite(plx) & (plx >= plx_min) & (plx <= plx_max)
        n_finite = int(np.isfinite(plx).sum())
        if n_finite == 0:
            return None
        return mask

    def _quality_column(self, name: str):
        """A Gaia quality column from the CMD table, as floats, or None."""
        df = getattr(self, "df", None)
        if df is None or name not in df.columns:
            return None
        vals = pd.to_numeric(df[name], errors="coerce").to_numpy(dtype=float)
        return vals if np.isfinite(vals).any() else None

    def _quality_mask(self):
        """Combined RUWE / C* display filter, or None when neither is on.

        Missing values are KEPT, matching `gaia_quality_mask`: a star with no
        RUWE has not failed the cut, it was never measured. Rejecting those
        would quietly drop every non-Gaia source from the plot.
        """
        mask = None
        n_before = None
        parts: list[str] = []

        if self.ruwe_check.isChecked():
            ruwe = self._quality_column("ruwe")
            if ruwe is not None:
                thr = float(self.ruwe_spin.value())
                keep = ~(np.isfinite(ruwe) & (ruwe > thr))
                n_before = keep.size
                mask = keep if mask is None else (mask & keep)
                parts.append(f"RUWE −{int((~keep).sum())}")
            else:
                parts.append("RUWE n/a")

        if self.cstar_check.isChecked():
            excess = self._quality_column("phot_bp_rp_excess_factor")
            bp_rp = self._quality_column("gaia_BP_RP")
            if bp_rp is None:
                bp = self._quality_column("gaia_BP")
                rp = self._quality_column("gaia_RP")
                bp_rp = (bp - rp) if (bp is not None and rp is not None) else None
            gmag = self._quality_column("gaia_G")
            if excess is not None and bp_rp is not None and gmag is not None:
                cstar = gaia_corrected_excess_factor(bp_rp, excess)
                sigma = gaia_cstar_sigma(gmag)
                nsig = float(self.cstar_spin.value())
                bad = (np.isfinite(cstar) & np.isfinite(sigma)
                       & (np.abs(cstar) > nsig * sigma))
                keep = ~bad
                n_before = keep.size
                mask = keep if mask is None else (mask & keep)
                parts.append(f"C* −{int(bad.sum())}")
            else:
                parts.append("C* n/a (no excess factor)")

        label = getattr(self, "quality_info_label", None)
        if label is not None:
            if mask is not None and n_before:
                label.setText(f"({' · '.join(parts)} → {int(mask.sum())}/{n_before})")
            else:
                label.setText(f"({' · '.join(parts)})" if parts else "")
        return mask

    def _on_quality_filter_changed(self):
        """Warn once when a cut is asked for but its column is not there."""
        if self.cstar_check.isChecked() and \
                self._quality_column("phot_bp_rp_excess_factor") is None:
            self.cstar_check.blockSignals(True)
            self.cstar_check.setChecked(False)
            self.cstar_check.blockSignals(False)
            from PyQt5.QtWidgets import QMessageBox
            QMessageBox.information(
                self, "C* Unavailable",
                "phot_bp_rp_excess_factor is not in this CMD table.\n\n"
                "Catalogues fetched before 2026-08-11 did not request it. "
                "Re-run Step 5 (Gaia query) and Step 6 (Master Catalog Build) "
                "to make the C* filter available.")
        self._redraw()

    def _membership_mode_key(self) -> str:
        idx = int(self.member_mode_combo.currentIndex())
        if idx == 1:
            return "loose"
        if idx == 2:
            return "normal"
        if idx == 3:
            return "strict"
        return "off"

    def _membership_threshold(self, mode_key: str) -> float:
        if mode_key == "loose":
            return 0.30
        if mode_key == "strict":
            return 0.80
        return 0.50

    def _on_membership_ui_changed(self):
        if self.params is not None and hasattr(self.params, "P"):
            self.params.P.cmd_membership_mode = self._membership_mode_key()
            self.params.P.cmd_membership_compare = bool(self.member_compare.isChecked())
            if hasattr(self.params, "save_toml"):
                try:
                    self.params.save_toml()
                except Exception:
                    pass
        self._redraw()

    def _pick_existing_membership_col(self):
        cands = [
            "gaia_pmem",
            "pmem_gaia",
            "membership_prob_gaia",
            "membership_prob",
            "pmem",
        ]
        for c in cands:
            if c not in self.df.columns:
                continue
            v = pd.to_numeric(self.df[c], errors="coerce").to_numpy(float)
            if np.isfinite(v).sum() >= 10:
                return c, v
        return None, None

    def _merge_columns_from_gaia_derived(self, needed_cols):
        self.df = _merge_gaia_columns_from_catalog(self.df, self.result_dir, needed_cols)

    def _ensure_membership_columns_from_master(self):
        needed = [
            "source_id", "gaia_source_id",
            "pmra", "pmdec", "parallax",
            "pmra_error", "pmdec_error", "parallax_error",
            "ruwe", "visibility_periods_used",
            "gaia_pmem", "pmem_gaia", "membership_prob_gaia",
        ]
        missing = [c for c in needed if c not in self.df.columns]
        if not missing:
            return self._merge_columns_from_gaia_derived(needed)

        try:
            master, _, _ = _load_master_table(self.result_dir)
        except Exception:
            return
        if master.empty:
            return

        key = None
        if "ID" in self.df.columns and "ID" in master.columns:
            key = "ID"
        elif "source_id" in self.df.columns and "source_id" in master.columns:
            key = "source_id"
            self.df["source_id"] = parse_int64_series(self.df["source_id"]).astype("Int64")
            master["source_id"] = parse_int64_series(master["source_id"]).astype("Int64")
        elif "gaia_source_id" in self.df.columns and "gaia_source_id" in master.columns:
            key = "gaia_source_id"
            self.df["gaia_source_id"] = parse_int64_series(self.df["gaia_source_id"]).astype("Int64")
            master["gaia_source_id"] = parse_int64_series(master["gaia_source_id"]).astype("Int64")
        if key is None:
            return

        add_cols = [c for c in missing if c in master.columns and c != key]
        if add_cols:
            use = master[[key] + add_cols].copy()
            use = use.drop_duplicates(subset=[key], keep="first")
            try:
                self.df = self.df.merge(use, on=key, how="left")
            except Exception:
                return

        self._merge_columns_from_gaia_derived(needed)

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

    def _compute_gaia_membership_prob(self):
        req = ("pmra", "pmdec", "parallax")
        if not all(c in self.df.columns for c in req):
            return None, "pm/parallax columns missing"

        pmra = pd.to_numeric(self.df["pmra"], errors="coerce").to_numpy(float)
        pmdec = pd.to_numeric(self.df["pmdec"], errors="coerce").to_numpy(float)
        plx = pd.to_numeric(self.df["parallax"], errors="coerce").to_numpy(float)
        finite = np.isfinite(pmra) & np.isfinite(pmdec) & np.isfinite(plx)
        if int(finite.sum()) < 30:
            return None, "too few stars with finite pm/parallax"

        fit_mask = finite.copy()
        if "ruwe" in self.df.columns:
            ruwe = pd.to_numeric(self.df["ruwe"], errors="coerce").to_numpy(float)
            fit_mask &= (~np.isfinite(ruwe)) | (ruwe <= 2.0)
        if "visibility_periods_used" in self.df.columns:
            vpu = pd.to_numeric(self.df["visibility_periods_used"], errors="coerce").to_numpy(float)
            fit_mask &= (~np.isfinite(vpu)) | (vpu >= 8.0)
        if int(fit_mask.sum()) < 25:
            fit_mask = finite.copy()

        x_fit = np.column_stack([pmra[fit_mask], pmdec[fit_mask], plx[fit_mask]])
        model = self._fit_two_component_gmm(x_fit)
        if model is None:
            return None, "GMM fit failed"

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
        r0 = e0 / den
        r1 = e1 / den
        p_cluster = r0 if k_cluster == 0 else r1

        pmem = np.full(len(self.df), np.nan, dtype=float)
        pmem[finite] = p_cluster
        note = f"gaia_gmm_3d(valid={int(finite.sum())}, fit={int(fit_mask.sum())})"
        return pmem, note

    def _ensure_membership_prob(self):
        if self._membership_ready:
            return self._membership_prob
        self._membership_ready = True

        c_exist, v_exist = self._pick_existing_membership_col()
        if c_exist is not None:
            self._membership_prob = np.clip(np.asarray(v_exist, float), 0.0, 1.0)
            self._membership_source = c_exist
            self._membership_note = f"existing:{c_exist}"
            return self._membership_prob

        self._ensure_membership_columns_from_master()
        c_exist, v_exist = self._pick_existing_membership_col()
        if c_exist is not None:
            self._membership_prob = np.clip(np.asarray(v_exist, float), 0.0, 1.0)
            self._membership_source = c_exist
            self._membership_note = f"existing:{c_exist}"
            return self._membership_prob

        p_auto, note = self._compute_gaia_membership_prob()
        if p_auto is not None:
            self._membership_prob = np.clip(np.asarray(p_auto, float), 0.0, 1.0)
            self.df["gaia_pmem"] = self._membership_prob
            self._membership_source = "gaia_gmm_3d"
            self._membership_note = note
            return self._membership_prob

        self._membership_prob = None
        self._membership_source = "none"
        self._membership_note = note
        return None

    def _current_membership_mask(self):
        mode = self._membership_mode_key()
        if mode == "off":
            return None, mode, np.nan
        prob = self._ensure_membership_prob()
        if prob is None:
            return None, mode, self._membership_threshold(mode)
        thr = self._membership_threshold(mode)
        mask = np.isfinite(prob) & (prob >= thr)
        return mask, mode, thr

    def _save_membership_csv(self):
        prob = self._ensure_membership_prob()
        if prob is None:
            self.info_text.setPlainText(
                "Membership not available.\n"
                f"Reason: {self._membership_note}"
            )
            return
        out = self.result_dir / "cmd_with_gaia_membership.csv"
        save_cols = ["ID", "source_id", "gaia_source_id", "pmra", "pmdec", "parallax"]
        keep = [c for c in save_cols if c in self.df.columns]
        df_out = self.df[keep].copy() if keep else pd.DataFrame(index=self.df.index)
        df_out["gaia_pmem"] = prob
        df_out.to_csv(out, index=False, na_rep="NaN")
        self.info_text.setPlainText(f"Saved: {out}")

    def _build_figure(self):
        self.figure.clear()
        view_name = self.available_views[self.view_mode]

        # Single view modes
        # Padded right margin + cax pulled inward keeps the colorbar
        # label "Teff (K) + OBAFGKM-like color" fully on-canvas even
        # with long tick labels like "35000 K (O)".  top=0.88 leaves
        # breathing room for the plot title.
        if view_name == "inst":
            self.ax_inst = self.figure.add_subplot(1, 1, 1)
            self.ax_std = None
            self.ax_gaia = None
            self.figure.subplots_adjust(left=0.10, right=0.80, bottom=0.18, top=0.88)
            self.cax = self.figure.add_axes([0.84, 0.18, 0.018, 0.70])
        elif view_name == "std":
            self.ax_inst = None
            self.ax_std = self.figure.add_subplot(1, 1, 1)
            self.ax_gaia = None
            self.figure.subplots_adjust(left=0.10, right=0.80, bottom=0.18, top=0.88)
            self.cax = self.figure.add_axes([0.84, 0.18, 0.018, 0.70])
        elif view_name == "gaia":
            self.ax_inst = None
            self.ax_std = None
            self.ax_gaia = self.figure.add_subplot(1, 1, 1)
            self.figure.subplots_adjust(left=0.10, right=0.80, bottom=0.18, top=0.88)
            self.cax = self.figure.add_axes([0.84, 0.18, 0.018, 0.70])
        elif view_name == "all":
            # Show all available CMDs
            if self.has_std and self.gaia_mode is not None:
                self.ax_inst = self.figure.add_subplot(1, 3, 1)
                self.ax_std = self.figure.add_subplot(1, 3, 2)
                self.ax_gaia = self.figure.add_subplot(1, 3, 3)
                self.figure.subplots_adjust(left=0.055, right=0.85, bottom=0.16, top=0.85, wspace=0.30)
                self.cax = self.figure.add_axes([0.88, 0.16, 0.015, 0.68])
            elif self.has_std:
                self.ax_inst = self.figure.add_subplot(1, 2, 1)
                self.ax_std = self.figure.add_subplot(1, 2, 2)
                self.ax_gaia = None
                self.figure.subplots_adjust(left=0.07, right=0.85, bottom=0.16, top=0.85, wspace=0.30)
                self.cax = self.figure.add_axes([0.88, 0.16, 0.015, 0.68])
            elif self.gaia_mode is not None:
                self.ax_inst = self.figure.add_subplot(1, 2, 1)
                self.ax_std = None
                self.ax_gaia = self.figure.add_subplot(1, 2, 2)
                self.figure.subplots_adjust(left=0.07, right=0.85, bottom=0.16, top=0.85, wspace=0.30)
                self.cax = self.figure.add_axes([0.88, 0.16, 0.015, 0.68])
            else:
                self.ax_inst = self.figure.add_subplot(1, 1, 1)
                self.ax_std = None
                self.ax_gaia = None
                self.figure.subplots_adjust(left=0.13, right=0.80, bottom=0.16, top=0.85)
                self.cax = self.figure.add_axes([0.84, 0.18, 0.018, 0.66])
        else:
            self.ax_inst = self.figure.add_subplot(1, 1, 1)
            self.ax_std = None
            self.ax_gaia = None
            self.figure.subplots_adjust(left=0.13, right=0.80, bottom=0.14, top=0.85)
            self.cax = self.figure.add_axes([0.84, 0.18, 0.018, 0.66])

        self.figure.patch.set_facecolor("black")
        for ax in [self.ax_inst, self.ax_std, self.ax_gaia]:
            if ax is None:
                continue
            ax.set_facecolor("black")
            for sp in ax.spines.values():
                sp.set_color("white")
            ax.tick_params(colors="white")
            ax.xaxis.label.set_color("white")
            ax.yaxis.label.set_color("white")
            ax.title.set_color("white")

        self.cax.set_facecolor("black")
        sm = mpl.cm.ScalarMappable(norm=self.ob_norm, cmap=self.ob_cmap)
        sm.set_array([])
        cbar = self.figure.colorbar(sm, cax=self.cax)
        cbar.set_label("Teff (K) + OBAFGKM-like color", fontsize=9, color="white", labelpad=14)
        ticks = [35000, 20000, 10000, 7500, 6000, 4500, 3200]
        labels = ["35000 K (O)", "20000 K (B)", "10000 K (A)", " 7500 K (F)", " 6000 K (G)", " 4500 K (K)", " 3200 K (M)"]
        cbar.set_ticks(ticks)
        cbar.set_ticklabels(labels)
        cbar.ax.tick_params(colors="white")
        for sp in cbar.ax.spines.values():
            sp.set_color("white")

    def _safe_float(self, series):
        return pd.to_numeric(series, errors="coerce").to_numpy(float)

    def _teff_from_color_index(self, color_x: np.ndarray, mode: str):
        if color_x.size == 0:
            return np.full_like(color_x, np.nan, dtype=float)
        parts = mode.split("-", 1)
        if len(parts) == 2:
            return _teff_from_color(color_x, parts[0], parts[1], self.teff_vmin, self.teff_vmax)
        return np.full_like(color_x, np.nan, dtype=float)

    def _get_y_mode(self, yval):
        if yval in self.y_scalar_opts:
            return ("scalar", yval)
        if isinstance(yval, str) and "-" in yval:
            a, b = yval.split("-", 1)
            if (a, b) in self.y_color_pairs:
                return ("color", (a, b))
        return (None, None)

    def _compute_arrays_and_mask(self, system: str, x_pair, y_choice, snr_cut: float, membership_mask=None):
        a, b = x_pair
        col_ax = f"mag_{system}_{a}"
        col_bx = f"mag_{system}_{b}"
        if (col_ax not in self.df.columns) or (col_bx not in self.df.columns):
            return np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])

        Ax = self._safe_float(self.df[col_ax])
        Bx = self._safe_float(self.df[col_bx])
        xcolor = Ax - Bx
        x = xcolor

        y_mode, y_param = self._get_y_mode(y_choice)
        involved = set([a, b])

        if y_mode == "scalar":
            by = y_param
            col_y = f"mag_{system}_{by}"
            if col_y not in self.df.columns:
                return np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])
            y = self._safe_float(self.df[col_y])
            # Display-only manual ZP: optional extra shift of the Instrumental
            # magnitude axis. mag_inst already carries the IRAF Z=25 convention
            # (baked in at Step 7), so this defaults to 0; it only adds a manual
            # nudge on top. Only the magnitude (scalar) axis of the Instrumental
            # view is shifted; colors and the Std view are untouched (a constant
            # ZP cancels in any a-b color).
            if system == "inst" and getattr(self, "manual_zp_check", None) is not None \
                    and self.manual_zp_check.isChecked():
                y = y + float(self.manual_zp_spin.value())
            involved.add(by)
        elif y_mode == "color":
            ya, yb = y_param
            col_ya = f"mag_{system}_{ya}"
            col_yb = f"mag_{system}_{yb}"
            if (col_ya not in self.df.columns) or (col_yb not in self.df.columns):
                return np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])
            Ya = self._safe_float(self.df[col_ya])
            Yb = self._safe_float(self.df[col_yb])
            y = Ya - Yb
            involved.update([ya, yb])
        else:
            return np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])

        mask = np.isfinite(x) & np.isfinite(y)

        if snr_cut > 0 and self.has_snr:
            for band in involved:
                sc = f"snr_{band}"
                if sc in self.df.columns:
                    sv = self._safe_float(self.df[sc])
                    # Only exclude stars with known (finite) SNR below threshold;
                    # NaN SNR means unmeasured → keep (do not reject unknowns)
                    mask &= ~(np.isfinite(sv) & (sv < snr_cut))

        if membership_mask is not None and len(membership_mask) == len(mask):
            mask &= np.asarray(membership_mask, bool)

        return x[mask], y[mask], mask, xcolor[mask]

    def _mag_error_array(self, system: str, band: str) -> np.ndarray:
        candidates = []
        if system == "std":
            candidates.append(f"mag_std_err_{band}")
        candidates.append(f"mag_inst_err_{band}")
        for col in candidates:
            if col in self.df.columns:
                arr = self._safe_float(self.df[col])
                return np.where(np.isfinite(arr) & (arr >= 0), arr, np.nan)
        return np.full(len(self.df), np.nan, dtype=float)

    def _quadrature_error(self, *arrays: np.ndarray) -> np.ndarray:
        if not arrays:
            return np.full(len(self.df), np.nan, dtype=float)
        stack = np.vstack([np.asarray(a, dtype=float) for a in arrays])
        finite = np.isfinite(stack).all(axis=0)
        out = np.full(stack.shape[1], np.nan, dtype=float)
        out[finite] = np.sqrt(np.sum(stack[:, finite] ** 2, axis=0))
        return out

    def _compute_cmd_error_arrays(self, system: str, x_pair, y_choice, mask: np.ndarray):
        if mask is None or len(mask) != len(self.df):
            return np.array([]), np.array([])
        a, b = x_pair
        xerr_full = self._quadrature_error(
            self._mag_error_array(system, a),
            self._mag_error_array(system, b),
        )

        y_mode, y_param = self._get_y_mode(y_choice)
        if y_mode == "scalar":
            yerr_full = self._mag_error_array(system, y_param)
        elif y_mode == "color":
            ya, yb = y_param
            yerr_full = self._quadrature_error(
                self._mag_error_array(system, ya),
                self._mag_error_array(system, yb),
            )
        else:
            yerr_full = np.full(len(self.df), np.nan, dtype=float)

        return xerr_full[mask], yerr_full[mask]

    def _plot_cmd_errorbars(self, ax, x, y, mask: np.ndarray, system: str, x_pair, y_choice):
        show_x = getattr(self, "xerr_check", None) is not None and self.xerr_check.isChecked()
        show_y = getattr(self, "yerr_check", None) is not None and self.yerr_check.isChecked()
        if not (show_x or show_y) or len(x) == 0:
            return

        xerr, yerr = self._compute_cmd_error_arrays(system, x_pair, y_choice, mask)
        if len(xerr) != len(x) or len(yerr) != len(y):
            return

        finite = np.isfinite(x) & np.isfinite(y)
        if show_x:
            finite &= np.isfinite(xerr) & (xerr > 0)
        if show_y:
            finite &= np.isfinite(yerr) & (yerr > 0)
        idx = np.flatnonzero(finite)
        if idx.size == 0:
            return

        max_bars = 400
        if idx.size > max_bars:
            idx = idx[np.linspace(0, idx.size - 1, max_bars).astype(int)]

        ax.errorbar(
            np.asarray(x)[idx],
            np.asarray(y)[idx],
            xerr=np.asarray(xerr)[idx] if show_x else None,
            yerr=np.asarray(yerr)[idx] if show_y else None,
            fmt="none",
            ecolor="#DDE7F0",
            elinewidth=0.55,
            alpha=0.35,
            capsize=0,
            zorder=1,
        )

    def _compute_gaia_arrays_and_mask(self, snr_cut: float, membership_mask=None):
        if self.gaia_mode is None:
            return np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])

        if self.gaia_mode == "inst":
            G = self._safe_float(self.df["gaia_G_inst"])
            C = self._safe_float(self.df["gaia_BP_RP_inst"])
        else:
            G = self._safe_float(self.df["gaia_G_syn"])
            C = self._safe_float(self.df["gaia_BP_RP_syn"])

        mask = np.isfinite(G) & np.isfinite(C)

        if snr_cut > 0 and self.has_snr:
            for band in ("g", "r", "i"):
                sc = f"snr_{band}"
                if sc in self.df.columns:
                    sv = self._safe_float(self.df[sc])
                    mask &= ~(np.isfinite(sv) & (sv < snr_cut))

        if membership_mask is not None and len(membership_mask) == len(mask):
            mask &= np.asarray(membership_mask, bool)

        return C[mask], G[mask], mask, C[mask]

    def _style_axis(self, ax):
        ax.set_facecolor("black")
        for sp in ax.spines.values():
            sp.set_color("white")
        ax.tick_params(colors="white")
        ax.xaxis.label.set_color("white")
        ax.yaxis.label.set_color("white")
        ax.title.set_color("white")

    def _apply_y_orientation(self):
        axes = [self.ax_inst, self.ax_std, self.ax_gaia]
        for ax in axes:
            if ax is None:
                continue
            ymin, ymax = ax.get_ylim()
            if self.invert_y.isChecked():
                if ymin < ymax:
                    ax.set_ylim(ymax, ymin)
            else:
                if ymin > ymax:
                    ax.set_ylim(ymax, ymin)

    def _redraw(self):
        ax_primary = self.ax_inst or self.ax_std or self.ax_gaia
        self._plot_cache = {}
        if not self.x_pairs:
            if ax_primary is not None:
                ax_primary.clear()
                ax_primary.set_title("No available X color index", fontsize=10, color="white")
            self.canvas.draw_idle()
            return

        x_text = self.x_combo.currentText()
        if x_text not in [f"{a}-{b}" for (a, b) in self.x_pairs]:
            self.canvas.draw_idle()
            return
        a, b = x_text.split("-", 1)
        x_pair = (a, b)

        yval = self.y_combo.currentText()
        snr_cut = float(self.snr_spin.value())
        member_mask, member_mode, member_thr = self._current_membership_mask()
        member_active = (member_mode != "off") and (member_mask is not None)
        member_compare = bool(self.member_compare.isChecked()) and member_active

        plx_mask = self._parallax_mask()
        plx_active = plx_mask is not None

        roi_mask = self._roi_mask()
        roi_active = roi_mask is not None

        # background (grey dots): spatial/parallax pre-filter (everything inside ROI or parallax range)
        if plx_active and roi_active:
            bg_mask = plx_mask & roi_mask
        elif plx_active:
            bg_mask = plx_mask
        elif roi_active:
            bg_mask = roi_mask
        else:
            bg_mask = None

        # Gaia quality (RUWE / C*) filters the *drawn* sample, background
        # included: the point of enabling them is to see the CMD without the
        # contaminated stars, so leaving them in the grey layer would defeat it.
        qual_mask = self._quality_mask()
        if qual_mask is not None:
            bg_mask = qual_mask if bg_mask is None else (bg_mask & qual_mask)

        # foreground mask: member & spatial filters combined
        spatial_mask = bg_mask  # reuse combined spatial pre-filter
        spatial_active = spatial_mask is not None
        if member_active and spatial_active:
            effective_mask = np.asarray(member_mask, bool) & spatial_mask
        elif member_active:
            effective_mask = member_mask
        elif spatial_active:
            effective_mask = spatial_mask
        else:
            effective_mask = None

        if self.ax_inst is not None:
            self.ax_inst.clear()
        if self.ax_std is not None:
            self.ax_std.clear()
        if self.ax_gaia is not None:
            self.ax_gaia.clear()

        if self.ax_inst is not None:
            x_i_all, y_i_all, mask_i_all, xcol_i_all = self._compute_arrays_and_mask("inst", x_pair, yval, snr_cut, bg_mask)
            if member_active or plx_active or roi_active:
                x_i, y_i, mask_i, xcol_i = self._compute_arrays_and_mask("inst", x_pair, yval, snr_cut, effective_mask)
            else:
                x_i, y_i, mask_i, xcol_i = x_i_all, y_i_all, mask_i_all, xcol_i_all
            teff_i = self._teff_from_color_index(xcol_i, f"{a}-{b}")
        else:
            x_i_all, y_i_all = np.array([]), np.array([])
            x_i, y_i, mask_i, teff_i = np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])

        if self.has_std and self.ax_std is not None:
            x_s_all, y_s_all, mask_s_all, xcol_s_all = self._compute_arrays_and_mask("std", x_pair, yval, snr_cut, bg_mask)
            if member_active or plx_active or roi_active:
                x_s, y_s, mask_s, xcol_s = self._compute_arrays_and_mask("std", x_pair, yval, snr_cut, effective_mask)
            else:
                x_s, y_s, mask_s, xcol_s = x_s_all, y_s_all, mask_s_all, xcol_s_all
            teff_s = self._teff_from_color_index(xcol_s, f"{a}-{b}")
        else:
            x_s_all, y_s_all = np.array([]), np.array([])
            x_s, y_s, mask_s, teff_s = np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])

        if self.gaia_mode is not None and self.ax_gaia is not None:
            x_g_all, y_g_all, mask_g_all, xcol_g_all = self._compute_gaia_arrays_and_mask(snr_cut, bg_mask)
            if member_active or plx_active or roi_active:
                x_g, y_g, mask_g, xcol_g = self._compute_gaia_arrays_and_mask(snr_cut, effective_mask)
            else:
                x_g, y_g, mask_g, xcol_g = x_g_all, y_g_all, mask_g_all, xcol_g_all
            teff_g = self._teff_from_color_index(xcol_g, "BP-RP")
        else:
            x_g_all, y_g_all = np.array([]), np.array([])
            x_g, y_g, mask_g, teff_g = np.array([]), np.array([]), np.zeros(len(self.df), bool), np.array([])

        color_title = f"{a}-{b}"

        def _count_text(n, n_all, active):
            return f"N={n}/{n_all}" if active else f"N={n}"

        if self.ax_inst is not None:
            self._style_axis(self.ax_inst)
            if member_compare and len(x_i_all) > 0:
                self.ax_inst.scatter(x_i_all, y_i_all, s=10, alpha=0.22, linewidths=0, rasterized=True, c="#9E9E9E")
            if len(x_i) > 0:
                self._plot_cmd_errorbars(self.ax_inst, x_i, y_i, mask_i, "inst", x_pair, yval)
                self.ax_inst.scatter(x_i, y_i, s=12, alpha=0.92, linewidths=0, rasterized=True, c=teff_i, cmap=self.ob_cmap, norm=self.ob_norm)
                self.ax_inst.set_title(
                    format_cmd_title(
                        self.params,
                        yval,
                        color_title,
                        system_label="Instrumental",
                        count_text=_count_text(len(x_i), len(x_i_all), member_active or plx_active or roi_active),
                        result_dir=self.result_dir,
                    ),
                    fontsize=11,
                    color="white",
                )
                self._plot_cache[self.ax_inst] = {
                    "system": "inst",
                    "x": x_i,
                    "y": y_i,
                    "df_index": np.where(mask_i)[0],
                }
            else:
                self.ax_inst.set_title(
                    format_cmd_title(
                        self.params,
                        yval,
                        color_title,
                        system_label="Instrumental",
                        count_text=_count_text(0, len(x_i_all), member_active or plx_active or roi_active),
                        result_dir=self.result_dir,
                    ),
                    fontsize=11,
                    color="white",
                )

        if self.ax_std is not None:
            std_label = photometric_system_label(a, b, yval)
            self._style_axis(self.ax_std)
            if member_compare and len(x_s_all) > 0:
                self.ax_std.scatter(x_s_all, y_s_all, s=10, alpha=0.22, linewidths=0, rasterized=True, c="#9E9E9E")
            if len(x_s) > 0:
                self._plot_cmd_errorbars(self.ax_std, x_s, y_s, mask_s, "std", x_pair, yval)
                self.ax_std.scatter(x_s, y_s, s=12, alpha=0.92, linewidths=0, rasterized=True, c=teff_s, cmap=self.ob_cmap, norm=self.ob_norm)
                self.ax_std.set_title(
                    format_cmd_title(
                        self.params,
                        yval,
                        color_title,
                        system_label=std_label,
                        count_text=_count_text(len(x_s), len(x_s_all), member_active or plx_active or roi_active),
                        result_dir=self.result_dir,
                    ),
                    fontsize=11,
                    color="white",
                )
                self._plot_cache[self.ax_std] = {
                    "system": "std",
                    "x": x_s,
                    "y": y_s,
                    "df_index": np.where(mask_s)[0],
                }
            else:
                self.ax_std.set_title(
                    format_cmd_title(
                        self.params,
                        yval,
                        color_title,
                        system_label=std_label,
                        count_text=_count_text(0, len(x_s_all), member_active or plx_active or roi_active),
                        result_dir=self.result_dir,
                    ),
                    fontsize=11,
                    color="white",
                )

        if self.ax_gaia is not None:
            self._style_axis(self.ax_gaia)
            if member_compare and len(x_g_all) > 0:
                self.ax_gaia.scatter(x_g_all, y_g_all, s=10, alpha=0.22, linewidths=0, rasterized=True, c="#9E9E9E")
            if len(x_g) > 0:
                self.ax_gaia.scatter(x_g, y_g, s=12, alpha=0.92, linewidths=0, rasterized=True, c=teff_g, cmap=self.ob_cmap, norm=self.ob_norm)
                gaia_label = "Gaia instrumental" if self.gaia_mode == "inst" else "Gaia synthetic"
                self.ax_gaia.set_title(
                    format_cmd_title(
                        self.params,
                        "G",
                        "BP-RP",
                        system_label=gaia_label,
                        count_text=_count_text(len(x_g), len(x_g_all), member_active or plx_active or roi_active),
                        result_dir=self.result_dir,
                    ),
                    fontsize=11,
                    color="white",
                )
                self._plot_cache[self.ax_gaia] = {
                    "system": "gaia",
                    "x": x_g,
                    "y": y_g,
                    "df_index": np.where(mask_g)[0],
                }
            else:
                gaia_label = "Gaia instrumental" if self.gaia_mode == "inst" else "Gaia synthetic"
                self.ax_gaia.set_title(
                    format_cmd_title(
                        self.params,
                        "G",
                        "BP-RP",
                        system_label=gaia_label,
                        count_text=_count_text(0, len(x_g_all), member_active or plx_active or roi_active),
                        result_dir=self.result_dir,
                    ),
                    fontsize=11,
                    color="white",
                )

        x_label = f"{a}-{b} (mag)"
        y_label = f"{yval} (mag)" if yval else ""
        if self.ax_inst is not None:
            self.ax_inst.set_xlabel(x_label)
            self.ax_inst.set_ylabel(y_label)
        if self.ax_std is not None:
            self.ax_std.set_xlabel(x_label)
            self.ax_std.set_ylabel(y_label)
        if self.ax_gaia is not None:
            self.ax_gaia.set_xlabel("BP-RP (mag)")
            self.ax_gaia.set_ylabel("G (mag)")

        self._apply_y_orientation()

        def _rng(arr):
            arr = np.asarray(arr, float)
            arr = arr[np.isfinite(arr)]
            if arr.size == 0:
                return "n/a"
            return f"{arr.min():.0f}-{arr.max():.0f} K"

        lines = [
            f"X={a}-{b}, Y={yval}, SNR>={snr_cut:.0f}",
        ]
        if member_mode == "off":
            lines.append("Membership: Off")
        elif member_mask is None:
            lines.append(f"Membership: {member_mode} (P>={member_thr:.2f}) unavailable | {self._membership_note}")
        else:
            n_valid = int(np.isfinite(self._membership_prob).sum()) if self._membership_prob is not None else 0
            n_mem = int(np.asarray(member_mask, bool).sum())
            lines.append(
                f"Membership: {member_mode} (P>={member_thr:.2f}) | source={self._membership_source} | "
                f"selected={n_mem}/{n_valid} | compare={member_compare}"
            )
        if self.plx_check.isChecked():
            plx = self._parallax_values()
            if plx is None:
                lines.append("Parallax: unavailable")
            else:
                finite = np.isfinite(plx)
                n_sel = int(plx_mask.sum()) if plx_mask is not None else 0
                lines.append(
                    f"Parallax: {self.plx_min_spin.value():.3f}..{self.plx_max_spin.value():.3f} mas | "
                    f"selected={n_sel}/{int(finite.sum())}"
                )
        if self.inst_bands:
            lines.append(f"[Inst] N={len(x_i)}{'/' + str(len(x_i_all)) if member_active else ''} | Teff range: {_rng(teff_i)}")
        if self.has_std:
            lines.append(f"[Cal]  N={len(x_s)}{'/' + str(len(x_s_all)) if member_active else ''} | Teff range: {_rng(teff_s)}")
        if self.gaia_mode is not None:
            lines.append(f"[Gaia:{self.gaia_mode}] N={len(x_g)}{'/' + str(len(x_g_all)) if member_active else ''} | Teff range: {_rng(teff_g)}")
        if not self.has_snr:
            lines.append("(snr_* columns missing: SNR cut disabled)")
        if self.last_pick_info:
            lines.append(self.last_pick_info)
        if self.pick_log:
            lines.append("Pick log (latest 5):")
            lines.extend(self.pick_log[-5:])
        self.info_text.setPlainText("\n".join(lines))

        self.canvas.draw_idle()

    def _fmt_val(self, v, nd=3):
        try:
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                return "NaN"
            return f"{float(v):.{nd}f}"
        except Exception:
            return str(v)

    def _on_plot_click(self, event):
        if event.inaxes is None or event.xdata is None or event.ydata is None:
            return
        if not getattr(event, "dblclick", False):
            return
        if event.inaxes not in self._plot_cache:
            return
        cache = self._plot_cache[event.inaxes]
        x = cache["x"]
        y = cache["y"]
        if x.size == 0:
            return
        xy_disp = event.inaxes.transData.transform(np.column_stack([x, y]))
        click = np.array([event.x, event.y])
        d2 = np.sum((xy_disp - click) ** 2, axis=1)
        idx = int(np.argmin(d2))
        if d2[idx] > (12.0 ** 2):
            return
        df_idx = int(cache["df_index"][idx])
        row = self.df.iloc[df_idx]

        parts = [f"Pick[{cache['system']}] ID={row.get('ID', 'n/a')}", f"source_id={row.get('source_id', 'n/a')}"]
        for band in ("g", "r", "i"):
            c_inst = f"mag_inst_{band}"
            if c_inst in row:
                parts.append(f"{band}_inst={self._fmt_val(row.get(c_inst))}")
        for band in ("g", "r", "i"):
            c_cal = f"mag_cal_{band}"
            if c_cal in row:
                parts.append(f"{band}_cal={self._fmt_val(row.get(c_cal))}")
        for band in ("g", "r", "i"):
            c_std = f"mag_std_{band}"
            if c_std in row:
                parts.append(f"{band}_std={self._fmt_val(row.get(c_std))}")
        if "gaia_G" in row:
            parts.append(f"gaia_G={self._fmt_val(row.get('gaia_G'))}")
        if "gaia_BP" in row:
            parts.append(f"gaia_BP={self._fmt_val(row.get('gaia_BP'))}")
        if "gaia_RP" in row:
            parts.append(f"gaia_RP={self._fmt_val(row.get('gaia_RP'))}")
        for c in ("pmra", "pmdec", "parallax"):
            if c in row:
                parts.append(f"{c}={self._fmt_val(row.get(c))}")
        if self._membership_prob is not None and 0 <= df_idx < len(self._membership_prob):
            pm = self._membership_prob[df_idx]
            if np.isfinite(pm):
                parts.append(f"Pmem={float(pm):.3f}")
        msg = " | ".join(parts)
        self.last_pick_info = msg
        self.pick_log.append(msg)
        self._redraw()

    def _save_png(self):
        if not self.x_pairs:
            self.info_text.setPlainText("No available X color index")
            return
        a, b = self.x_combo.currentText().split("-", 1)
        yv = self.y_combo.currentText().replace(" ", "")

        if self.has_std and self.gaia_mode is not None:
            mode = f"inst_std_gaia{self.gaia_mode}"
        elif self.has_std:
            mode = "inst_std"
        elif self.gaia_mode is not None:
            mode = f"inst_gaia{self.gaia_mode}"
        else:
            mode = "inst_only"

        mem_tag = ""
        mkey = self._membership_mode_key()
        if mkey != "off":
            mem_tag = f"_mem{mkey}"
            if self.member_compare.isChecked():
                mem_tag += "_cmp"

        out = self.result_dir / f"cmd_{mode}_{a}-{b}_vs_{yv}_snr{int(self.snr_spin.value())}{mem_tag}_OBcolor_dark.png"
        self.figure.savefig(out, dpi=170, bbox_inches="tight", facecolor=self.figure.get_facecolor(), edgecolor="none")
        self.info_text.setPlainText(f"Saved: {out}")

    def keyPressEvent(self, event):
        super().keyPressEvent(event)

    def _switch_view(self, delta: int):
        """Switch between views: inst, std, gaia, all"""
        self.view_mode = (self.view_mode + delta) % len(self.available_views)
        self._update_view_label()
        self._build_figure()
        self._redraw()


class ZPFitPlotWidget(QWidget):
    """ZP linear fit (delta vs color) + per-frame ZP timeline plots."""

    FILT_COLORS = {"g": "#1976D2", "r": "#D32F2F", "i": "#388E3C",
                   "u": "#7B1FA2", "z": "#E65100", "y": "#00695C",
                   "ha": "#AD1457", "r_spec": "#B71C1C"}

    def __init__(self, result_dir: Path, parent=None):
        super().__init__(parent)
        self.result_dir = Path(result_dir)
        self._cal_df = None
        self._frame_df = None
        self._coeff_df = None
        self._artist_map = {}   # id(artist) -> list of row dicts
        self._setup_ui()

    def _filter_color(self, filt: str) -> str:
        """Return a consistent color for any filter name."""
        if filt in self.FILT_COLORS:
            return self.FILT_COLORS[filt]
        import hashlib
        h = int(hashlib.md5(filt.encode()).hexdigest()[:6], 16)
        # Keep it reasonably saturated (not too dark/light)
        r = 80 + (h >> 16 & 0xFF) % 140
        g = 80 + (h >> 8 & 0xFF) % 140
        b = 80 + (h & 0xFF) % 140
        return f"#{r:02X}{g:02X}{b:02X}"

    def _detected_filters(self) -> list[str]:
        """Collect all filter names present in loaded data."""
        filts: set[str] = set()
        if self._cal_df is not None:
            filts.update(c[len("delta_"):] for c in self._cal_df.columns if c.startswith("delta_"))
        if self._frame_df is not None and "filter" in self._frame_df.columns:
            filts.update(str(f) for f in self._frame_df["filter"].dropna().unique())
        if self._coeff_df is not None and "filter" in self._coeff_df.columns:
            filts.update(str(f) for f in self._coeff_df["filter"].dropna().unique())
        return sorted(filts)

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)

        ctrl = QHBoxLayout()
        ctrl.addWidget(QLabel("Filter:"))
        self._filt_combo = QComboBox()
        self._filt_combo.addItem("All")
        self._filt_combo.currentTextChanged.connect(self._redraw)
        ctrl.addWidget(self._filt_combo)
        ctrl.addSpacing(12)
        ctrl.addWidget(QLabel("Date:"))
        self._date_combo = QComboBox()
        self._date_combo.addItem("All")
        self._date_combo.currentTextChanged.connect(self._redraw)
        ctrl.addWidget(self._date_combo)
        ctrl.addSpacing(12)
        btn_reload = QPushButton("Reload")
        btn_reload.clicked.connect(lambda: self.reload())
        ctrl.addWidget(btn_reload)
        btn_save = QPushButton("Save PNG")
        btn_save.clicked.connect(self._save_png)
        ctrl.addWidget(btn_save)
        ctrl.addStretch()
        layout.addLayout(ctrl)

        self._fig = Figure(figsize=(12, 8))
        gs = self._fig.add_gridspec(2, 2, hspace=0.35, wspace=0.30)
        self._ax_fit = self._fig.add_subplot(gs[0, 0])
        self._ax_hist = self._fig.add_subplot(gs[0, 1])
        self._ax_frame = self._fig.add_subplot(gs[1, :])
        self._canvas = FigureCanvas(self._fig)
        self._canvas.mpl_connect("pick_event", self._on_pick)
        self._toolbar = NavigationToolbar(self._canvas, self)
        layout.addWidget(self._toolbar)
        layout.addWidget(self._canvas, 1)

        self._info_label = QLabel("Click a data point to see star info.")
        self._info_label.setStyleSheet(mono_note_style())
        layout.addWidget(self._info_label)

    def reload(self, result_dir: Path = None):
        if result_dir is not None:
            self.result_dir = Path(result_dir)
        out_dir = step10_zp_dir(self.result_dir)
        self._cal_df = None
        self._frame_df = None
        self._coeff_df = None

        cal_path = out_dir / "gaia_sdss_calibrator_by_ID.csv"
        if cal_path.exists() and cal_path.stat().st_size > 0:
            try:
                self._cal_df = pd.read_csv(cal_path)
            except Exception:
                self._cal_df = None

        frame_path = out_dir / "frame_zeropoint.csv"
        if frame_path.exists() and frame_path.stat().st_size > 0:
            try:
                self._frame_df = pd.read_csv(frame_path)
            except Exception:
                self._frame_df = None

        coeff_path = out_dir / "zp_fit_coefficients.csv"
        if coeff_path.exists() and coeff_path.stat().st_size > 0:
            try:
                self._coeff_df = pd.read_csv(coeff_path)
            except Exception:
                self._coeff_df = None

        self._update_filter_combo()
        self._update_date_combo()
        self._redraw()

    def _update_filter_combo(self) -> None:
        """Repopulate filter combo from currently loaded data."""
        prev = self._filt_combo.currentText()
        self._filt_combo.blockSignals(True)
        self._filt_combo.clear()
        self._filt_combo.addItem("All")
        for f in self._detected_filters():
            self._filt_combo.addItem(f)
        idx = self._filt_combo.findText(prev)
        self._filt_combo.setCurrentIndex(max(idx, 0))
        self._filt_combo.blockSignals(False)

    def _update_date_combo(self):
        import re
        prev = self._date_combo.currentText()
        self._date_combo.blockSignals(True)
        self._date_combo.clear()
        self._date_combo.addItem("All")
        if self._frame_df is not None and "file" in self._frame_df.columns:
            def _extract_date(fname):
                m = re.search(r"\d{4}-\d{2}-\d{2}", str(fname))
                return m.group(0) if m else None
            dates = sorted(set(d for d in self._frame_df["file"].apply(_extract_date) if d))
            for d in dates:
                self._date_combo.addItem(d)
        idx = self._date_combo.findText(prev)
        self._date_combo.setCurrentIndex(max(idx, 0))
        self._date_combo.blockSignals(False)

    def _redraw(self):
        self._ax_fit.cla()
        self._ax_hist.cla()
        self._ax_frame.cla()
        self._artist_map.clear()
        self._draw_fit_plot()
        self._draw_zp_hist()
        self._draw_frame_zp()
        self._fig.tight_layout(pad=2.0)
        self._canvas.draw_idle()

    def _draw_fit_plot(self):
        ax = self._ax_fit
        filt_sel = self._filt_combo.currentText()

        if self._cal_df is None:
            ax.set_title("ZP Linear Fit — run ZP calibration first")
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes,
                    fontsize=14, color="gray")
            return

        cal = self._cal_df
        filts = self._detected_filters() if filt_sel == "All" else [filt_sel]

        color_labels = set()
        for filt in filts:
            fc = self._filter_color(filt)

            # Determine color index column from coefficients if available, else guess
            color_col = None
            if self._coeff_df is not None:
                row_c = self._coeff_df[self._coeff_df["filter"] == filt]
                if len(row_c) and "color_col" in row_c.columns:
                    cc = str(row_c["color_col"].iloc[0])
                    if cc != "none":
                        # coeff_df stores bare name ("B_V"); cal column has "color_" prefix
                        color_col = cc if cc.startswith("color_") else f"color_{cc}"
            if color_col is None or color_col not in cal.columns:
                # Fallback: scan cal for any color_* column
                for cand in cal.columns:
                    if cand.startswith("color_"):
                        color_col = cand
                        break
            if color_col is None or color_col not in cal.columns:
                continue

            ref_col = f"ref_{filt}"
            inst_col = f"mag_inst_{filt}"
            err_col = f"mag_inst_err_{filt}"
            delta_col = f"delta_{filt}"

            if delta_col in cal.columns:
                delta = cal[delta_col].to_numpy(float)
            elif ref_col in cal.columns and inst_col in cal.columns:
                delta = cal[ref_col].to_numpy(float) - cal[inst_col].to_numpy(float)
            else:
                continue

            color_x = cal[color_col].to_numpy(float)
            mask = np.isfinite(delta) & np.isfinite(color_x)
            if mask.sum() == 0:
                continue

            x_plot = color_x[mask]
            y_plot = delta[mask]
            stars_idx = np.where(mask)[0]

            err_arr = None
            if err_col in cal.columns:
                err_raw = cal[err_col].to_numpy(float)[mask]
                err_arr = np.where(np.isfinite(err_raw), err_raw, np.nan)

            # Sigma-clip outlier mask re-computed from saved fit line
            inlier = np.ones(len(x_plot), dtype=bool)
            zp_val = ct_val = None
            ct2_val = 0.0
            fit_label = f"{filt} (N={mask.sum()})"
            if self._coeff_df is not None:
                row = self._coeff_df[self._coeff_df["filter"] == filt]
                if len(row):
                    zp_val = float(row["zp"].iloc[0])
                    ct_val = float(row["ct"].iloc[0])
                    if "ct2" in row.columns and np.isfinite(row["ct2"].iloc[0]):
                        ct2_val = float(row["ct2"].iloc[0])
                    N_val = int(row["N"].iloc[0]) if "N" in row.columns else mask.sum()
                    sc_val = float(row["scatter_rms"].iloc[0]) if "scatter_rms" in row.columns else np.nan
                    resid = y_plot - (zp_val + ct_val * x_plot + ct2_val * x_plot**2)
                    med_r = np.nanmedian(resid)
                    mad_r = np.nanmedian(np.abs(resid - med_r)) + 1e-12
                    sig_r = MAD_TO_SIGMA * mad_r
                    inlier = np.abs(resid - med_r) <= 3.0 * sig_r
                    sc_str = f"σ={sc_val:.4f}" if np.isfinite(sc_val) else ""
                    ct2_str = f" CT2={ct2_val:+.3f}" if ct2_val else ""
                    fit_label = f"{filt}: ZP={zp_val:.3f} CT={ct_val:+.3f}{ct2_str} {sc_str} (N={N_val})"

            # Outliers first (behind inliers)
            if (~inlier).any():
                ax.scatter(x_plot[~inlier], y_plot[~inlier],
                           marker="x", c="gray", s=25, alpha=0.35, linewidths=1.0, zorder=2)

            sc = ax.scatter(
                x_plot[inlier], y_plot[inlier],
                c=fc, s=12, alpha=0.60,
                label=fit_label,
                picker=True, pickradius=6, zorder=3,
            )

            if err_arr is not None and np.isfinite(err_arr[inlier]).any():
                ax.errorbar(
                    x_plot[inlier], y_plot[inlier], yerr=err_arr[inlier],
                    fmt="none", ecolor=fc, elinewidth=0.7, capsize=2,
                    alpha=0.35, zorder=2,
                )

            row_list = []
            inlier_idx = np.where(inlier)[0]
            for i, si_local in enumerate(inlier_idx):
                si = stars_idx[si_local]
                row_list.append({
                    "filt": filt,
                    "ID": cal["ID"].iloc[si] if "ID" in cal.columns else "?",
                    "color": float(x_plot[si_local]),
                    "delta": float(y_plot[si_local]),
                    "ref_mag": float(cal[ref_col].iloc[si]) if ref_col in cal.columns else np.nan,
                    "inst_mag": float(cal[inst_col].iloc[si]) if inst_col in cal.columns else np.nan,
                    "err": float(err_arr[si_local]) if err_arr is not None and np.isfinite(err_arr[si_local]) else np.nan,
                })
            self._artist_map[id(sc)] = row_list

            if zp_val is not None and ct_val is not None:
                x_fit = np.linspace(np.nanmin(x_plot) - 0.05, np.nanmax(x_plot) + 0.05, 200)
                y_fit = zp_val + ct_val * x_fit + ct2_val * x_fit**2
                ax.plot(x_fit, y_fit, "-", color=fc, linewidth=2.0, zorder=4)

            color_labels.add(color_col)

        # X-axis label based on color columns used
        if len(color_labels) == 1:
            cc = next(iter(color_labels))
            color_label = cc.replace("color_", "").replace("_", " - ") + " (inst)"
        else:
            color_label = "color index (inst)"
        ax.set_xlabel(color_label)
        ax.set_ylabel("gaia_ref − inst (mag)")
        ax.set_title("ZP Linear Fit: gaia_ref − inst vs color  [× = sigma-clipped outliers]")
        if ax.collections or ax.get_lines():
            ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    def _draw_zp_hist(self):
        """ZP distribution histogram per filter (논문 AutoPHOT Fig.8 방식)."""
        from scipy.stats import gaussian_kde
        ax = self._ax_hist
        filt_sel = self._filt_combo.currentText()

        if self._frame_df is None:
            ax.set_title("ZP Distribution")
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes, fontsize=12, color="gray")
            return

        df = self._frame_df.copy()
        date_sel = self._date_combo.currentText()
        if date_sel != "All" and "file" in df.columns:
            import re
            df = df[df["file"].apply(lambda f: bool(re.search(re.escape(date_sel), str(f))))].copy()

        filts = (sorted(df["filter"].unique()) if "filter" in df.columns else []) if filt_sel == "All" else [filt_sel]
        has_any = False
        for filt in filts:
            sub = df[df["filter"] == filt] if "filter" in df.columns else pd.DataFrame()
            if sub.empty:
                continue
            zp_vals = pd.to_numeric(sub["zp_frame"], errors="coerce").dropna().to_numpy(float)
            zp_vals = zp_vals[np.isfinite(zp_vals)]
            if len(zp_vals) < 3:
                continue
            has_any = True
            fc = self._filter_color(filt)
            med = float(np.median(zp_vals))
            mad = float(np.median(np.abs(zp_vals - med)))
            sigma = MAD_TO_SIGMA * mad
            n_bins = max(8, min(30, len(zp_vals) // 2))
            ax.hist(zp_vals, bins=n_bins, color=fc, alpha=0.5,
                    label=f"{filt}: μ={med:.3f} σ={sigma:.3f} N={len(zp_vals)}")
            # KDE curve
            if len(zp_vals) >= 5:
                try:
                    kde = gaussian_kde(zp_vals, bw_method="scott")
                    xg = np.linspace(zp_vals.min() - 3 * sigma, zp_vals.max() + 3 * sigma, 200)
                    yk = kde(xg) * len(zp_vals) * (zp_vals.max() - zp_vals.min()) / n_bins
                    ax.plot(xg, yk, color=fc, linewidth=1.5)
                except Exception:
                    pass
            ax.axvline(med, color=fc, linestyle="--", linewidth=1.2, alpha=0.9)

        ax.set_xlabel("ZP (mag)")
        ax.set_ylabel("Count")
        ax.set_title("ZP Distribution")
        if has_any:
            ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    def _draw_frame_zp(self):
        import re
        ax = self._ax_frame
        filt_sel = self._filt_combo.currentText()
        date_sel = self._date_combo.currentText()

        if self._frame_df is None:
            ax.set_title("Per-Frame ZP — run ZP calibration first")
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes,
                    fontsize=14, color="gray")
            return

        df = self._frame_df.copy()
        if date_sel != "All" and "file" in df.columns:
            df = df[df["file"].apply(lambda f: bool(re.search(re.escape(date_sel), str(f))))].copy()

        if "filter" not in df.columns:
            ax.set_title("Per-Frame ZP — no filter column")
            return

        filts = (sorted(df["filter"].unique())) if filt_sel == "All" else [filt_sel]

        # Build shared x-tick labels from the first filter's file list
        _label_sub = None
        for filt in filts:
            sub0 = df[df["filter"] == filt].reset_index(drop=True)
            if not sub0.empty:
                _label_sub = sub0
                break

        for filt in filts:
            sub = df[df["filter"] == filt].reset_index(drop=True)
            if sub.empty:
                continue
            fc = self._filter_color(filt)
            x = np.arange(len(sub))
            y = sub["zp_frame"].to_numpy(float)
            yerr = sub["zp_scatter"].to_numpy(float) if "zp_scatter" in sub.columns else None

            # Point size ∝ n_ref
            has_nref = "n_ref" in sub.columns
            s_size = np.clip(sub["n_ref"].to_numpy(float) * 3.0, 12, 90) if has_nref else np.full(len(sub), 25)

            # Connection line
            ax.plot(x, y, "-", color=fc, alpha=0.45, linewidth=1.2, zorder=2)
            # Error bars (scatter)
            if yerr is not None and np.isfinite(yerr).any():
                ax.errorbar(x, y, yerr=yerr, fmt="none", ecolor=fc,
                           elinewidth=0.7, capsize=2, alpha=0.5, zorder=2)
            # Points with variable size
            n_label = int(sub["n_ref"].mean()) if has_nref else len(sub)
            ax.scatter(x, y, s=s_size, c=fc, alpha=0.85, zorder=3,
                      label=f"{filt}  N_frame={len(sub)}  <n_ref>≈{n_label}",
                      picker=True, pickradius=6)

            # Median ± σ band
            yf = y[np.isfinite(y)]
            if len(yf) >= 3:
                med_zp = float(np.median(yf))
                mad_zp = float(np.median(np.abs(yf - med_zp)))
                sigma_zp = MAD_TO_SIGMA * mad_zp
                ax.axhline(med_zp, color=fc, linestyle="--", linewidth=0.9, alpha=0.7)
                ax.axhspan(med_zp - sigma_zp, med_zp + sigma_zp,
                           color=fc, alpha=0.07, zorder=1)

        # Date-based x-tick labels from the first filter
        if _label_sub is not None and "file" in _label_sub.columns:
            n = len(_label_sub)
            tick_every = max(1, n // 18)
            tick_x = list(range(0, n, tick_every))
            tick_labels = []
            for i in tick_x:
                fname = str(_label_sub["file"].iloc[i])
                m = re.search(r"\d{4}-\d{2}-\d{2}", fname)
                if m:
                    tick_labels.append(m.group(0))
                else:
                    m2 = re.search(r"\d{8}", fname)
                    tick_labels.append(m2.group(0) if m2 else str(i))
            ax.set_xticks(tick_x)
            ax.set_xticklabels(tick_labels, rotation=35, ha="right", fontsize=7)
            xlabel = "Frame  (date from filename;  point size ∝ n_ref)"
        else:
            xlabel = "Frame index  (point size ∝ n_ref)"

        ax.set_xlabel(xlabel)
        ax.set_ylabel("ZP (mag)")
        ax.set_title("Per-Frame Zeropoint")
        if ax.get_lines() or ax.collections:
            ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    def _on_pick(self, event):
        artist = event.artist
        aid = id(artist)
        if aid not in self._artist_map:
            return
        rows = self._artist_map[aid]
        ind = event.ind
        if len(ind) == 0:
            return
        info = rows[ind[0]]
        parts = [
            f"Filter={info['filt']}",
            f"ID={info['ID']}",
            f"color={info['color']:.3f}",
            f"delta={info['delta']:.3f}",
            f"ref_mag={info['ref_mag']:.3f}",
            f"inst_mag={info['inst_mag']:.3f}",
        ]
        if np.isfinite(info["err"]):
            parts.append(f"err={info['err']:.4f}")
        self._info_label.setText(" | ".join(parts))

    def _save_png(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Plot", str(self.result_dir), "PNG Images (*.png)"
        )
        if path:
            self._fig.savefig(path, dpi=150, bbox_inches="tight")


class _AnchorDiscoveryWorker(QThread):
    """Background VizieR catalog discovery for the standard-anchor dialog.

    Network-bound (keyword search + per-candidate probe, ~10-60 s) so it must
    never run on the GUI thread; results come back via signals."""

    found = pyqtSignal(list)
    failed = pyqtSignal(str)

    def __init__(self, target_name, ra_deg, dec_deg, bands, parent=None):
        super().__init__(parent)
        self._args = (str(target_name), float(ra_deg), float(dec_deg), list(bands))

    def run(self):
        try:
            from apex.analysis.cmd.standard_anchor import discover_standard_catalogs

            name, ra, dec, bands = self._args
            self.found.emit(discover_standard_catalogs(name, ra, dec, bands))
        except Exception as exc:
            self.failed.emit(f"{type(exc).__name__}: {exc}")


class ZeropointCalibrationWindow(StepWindowBase):
    """Step 10: Zeropoint & Standardization"""

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.worker = None
        self.viewer = None
        self._current_zp_signature = None
        self._zp_cache_validation_result = None

        super().__init__(
            step_index=9,
            step_name="Zeropoint Calibration",
            params=params,
            project_state=project_state,
            main_window=main_window
        )

        self.setup_step_ui()
        self.restore_state()

    def setup_step_ui(self):
        # ── Controls ─────────────────────────────────────────────────────────
        info = QLabel("Build per-frame ZP calibration and standardized catalogs.")
        info.setProperty("role", "info")
        self.content_layout.addWidget(info)

        self.photometry_source_label = QLabel()
        self.photometry_source_label.setProperty("role", "caption")
        _f = self.photometry_source_label.font(); _f.setBold(True)
        self.photometry_source_label.setFont(_f)
        self.content_layout.addWidget(self.photometry_source_label)
        self._refresh_photometry_source_label()

        control_layout = QHBoxLayout()
        btn_params = create_parameter_button("Calibration Parameters")
        btn_params.clicked.connect(self.open_parameters_dialog)
        control_layout.addWidget(btn_params)

        self.run_bar = RunControlBar(
            "Run ZP Calibration", "Open Log",
            run_cb=self.run_analysis,
            stop_cb=self.stop_analysis,
            log_cb=self.show_log_window,
        )
        control_layout.addWidget(self.run_bar)
        self.btn_run = self.run_bar.btn_run
        self.btn_stop = self.run_bar.btn_stop
        self.content_layout.addLayout(control_layout)

        self.progress_label = QLabel("Idle")
        self.progress_label.setStyleSheet("QLabel { padding: 4px; }")
        self.content_layout.addWidget(self.progress_label)

        # ── ZP Fit Plot (takes remaining space) ──────────────────────────────
        self.fit_tab = ZPFitPlotWidget(self.params.P.result_dir)
        self.content_layout.addWidget(self.fit_tab, 1)

        # ── Log window with Workers panel ────────────────────────────────────
        _worker_group = QGroupBox("Workers")
        _worker_group.setMinimumWidth(300)
        _wg_layout = QVBoxLayout(_worker_group)
        _wg_layout.setContentsMargins(5, 5, 5, 5)
        self.worker_panel = WorkerStatusPanel(_worker_group)
        _wg_layout.addWidget(self.worker_panel)

        self.log_window = WorkflowLogWindow(
            self, "Calibration Log", width=850, height=420,
            side_widget=_worker_group,
        )
        self.log_text = self.log_window.log_text

    def _refresh_photometry_source_label(self, info: dict | None = None) -> None:
        provenance = info or resolve_cmd_photometry_input(
            self.params.P.result_dir,
            self.project_state,
        )
        self.photometry_source_label.setText(
            format_photometry_provenance(provenance)
        )

    def log(self, message: str):
        append_timestamped_log(self.log_text, message)

    def show_log_window(self):
        show_raised(self.log_window)

    def open_parameters_dialog(self):
        dialog, layout, buttons = build_scroll_param_dialog(
            self, "Calibration Parameters",
            info_text="Adjust zero-point calibration parameters. Changes apply to the next calibration run.",
            size=(500, 620),
        )

        match_group, match_container = create_collapsible_section("Matching", initial_expanded=True)
        match_form = QFormLayout(match_container)
        match_form.setContentsMargins(0, 0, 0, 0)

        self.param_pix = QDoubleSpinBox()
        self.param_pix.setRange(0.0, 50.0)
        self.param_pix.setValue(float(getattr(self.params.P, "pixel_scale_arcsec", 0.0) or 0.0))
        self.param_pix.setEnabled(False)
        match_form.addRow("Pixel scale (arcsec):", self.param_pix)

        self.param_match = QDoubleSpinBox()
        self.param_match.setRange(0.1, 20.0)
        self.param_match.setValue(float(getattr(self.params.P, "match_tol_px", 5.0)))
        match_form.addRow("Match tol (px):", self.param_match)

        self.param_min_match = QSpinBox()
        self.param_min_match.setRange(3, 1000)
        self.param_min_match.setValue(int(getattr(self.params.P, "min_master_gaia_matches", 10)))
        match_form.addRow("Min Gaia matches:", self.param_min_match)

        layout.addWidget(match_group)

        zp_group, zp_container = create_collapsible_section("Zero-point Fit", initial_expanded=True)
        zp_form = QFormLayout(zp_container)
        zp_form.setContentsMargins(0, 0, 0, 0)

        self.param_cmd_snr = QDoubleSpinBox()
        self.param_cmd_snr.setRange(0.0, 200.0)
        self.param_cmd_snr.setValue(float(getattr(self.params.P, "cmd_snr_calib_min", 20.0)))
        zp_form.addRow("CMD calib SNR min:", self.param_cmd_snr)

        self.param_frame_min = QSpinBox()
        self.param_frame_min.setRange(1, 1000)
        self.param_frame_min.setValue(int(getattr(self.params.P, "frame_zp_min_n", 5)))
        zp_form.addRow("Frame ZP min refs:", self.param_frame_min)

        self.param_apply_ext = QCheckBox("Enable")
        self.param_apply_ext.setChecked(bool(getattr(self.params.P, "cmd_apply_extinction", False)))
        zp_form.addRow("Apply extinction (k·X):", self.param_apply_ext)

        self.param_ext_mode = QComboBox()
        self.param_ext_mode.addItems(["absorb", "two_step"])
        self.param_ext_mode.setCurrentText(str(getattr(self.params.P, "cmd_extinction_mode", "absorb")))
        zp_form.addRow("Extinction mode:", self.param_ext_mode)

        self.param_clip = QDoubleSpinBox()
        self.param_clip.setRange(0.5, 10.0)
        self.param_clip.setValue(float(getattr(self.params.P, "zp_clip_sigma", 3.0)))
        zp_form.addRow("ZP clip sigma:", self.param_clip)

        self.param_iters = QSpinBox()
        self.param_iters.setRange(1, 20)
        self.param_iters.setValue(int(getattr(self.params.P, "zp_fit_iters", 5)))
        zp_form.addRow("ZP fit iters:", self.param_iters)

        self.param_slope = QDoubleSpinBox()
        self.param_slope.setRange(0.1, 5.0)
        self.param_slope.setValue(float(getattr(self.params.P, "zp_slope_absmax", 1.0)))
        zp_form.addRow("ZP slope abs max:", self.param_slope)

        layout.addWidget(zp_group)

        gaia_group, gaia_container = create_collapsible_section("Calibration Star Selection")
        gaia_form = QFormLayout(gaia_container)
        gaia_form.setContentsMargins(0, 0, 0, 0)

        self.param_gaia_snr = QDoubleSpinBox()
        self.param_gaia_snr.setRange(0.0, 200.0)
        self.param_gaia_snr.setValue(float(getattr(self.params.P, "gaia_snr_calib_min", 20.0)))
        self.param_gaia_snr.setToolTip("SNR threshold for calibration reference stars (applied in both global fit and per-frame ZP)")
        gaia_form.addRow("Calib star SNR min:", self.param_gaia_snr)

        self.param_gi_min = QDoubleSpinBox()
        self.param_gi_min.setRange(-2.0, 5.0)
        self.param_gi_min.setDecimals(2)
        self.param_gi_min.setValue(float(getattr(self.params.P, "gaia_gi_min", -0.5)))
        self.param_gi_min.setToolTip("Global lower bound for Gaia BP-RP color. Stars outside this range are excluded from all Jordi fits.")
        gaia_form.addRow("BP-RP min:", self.param_gi_min)

        self.param_gi_max = QDoubleSpinBox()
        self.param_gi_max.setRange(0.5, 6.0)
        self.param_gi_max.setDecimals(2)
        self.param_gi_max.setValue(float(getattr(self.params.P, "gaia_gi_max", 3.5)))
        self.param_gi_max.setToolTip("Global upper bound for Gaia BP-RP color. Stars outside this range are excluded from all Jordi fits.")
        gaia_form.addRow("BP-RP max:", self.param_gi_max)

        layout.addWidget(gaia_group)

        anchor_group, anchor_container = create_collapsible_section(
            "External Standard Anchor")
        anchor_form = QFormLayout(anchor_container)
        anchor_form.setContentsMargins(0, 0, 0, 0)

        anchor_info = QLabel(
            "Gaia 변환 참조가 부정확한 밴드(특히 Johnson U, σ≈0.20)를 "
            "VizieR 표준성 카탈로그(Landolt 계열, Gaia 독립)로 재앵커합니다.\n"
            "M67 실측: U 영점이 -0.13 mag 틀어져 이소크론 [M/H]가 "
            "-0.83으로 rail — 재앵커 후 +0.06(문헌 일치). 모든 밴드가 한 "
            "표준계에 앉아야 U-B가 축퇴를 풉니다."
        )
        anchor_info.setWordWrap(True)
        anchor_form.addRow(anchor_info)

        self.param_anchor_enable = QCheckBox("Enable")
        self.param_anchor_enable.setChecked(
            bool(getattr(self.params.P, "std_anchor_enable", False)))
        self.param_anchor_enable.setToolTip(
            "켜면 Step 10이 wide CMD 테이블 저장 직전에 mag_std_* 를 외부 "
            "표준성 오프셋만큼 이동합니다 (계측 측광은 불변). 오프셋 QC는 "
            "cmd_zeropoint/standard_anchor_offsets.csv 에 남습니다.")
        anchor_form.addRow("Anchor mag_std to standards:", self.param_anchor_enable)

        self.param_anchor_catalog = QLineEdit()
        self.param_anchor_catalog.setText(
            str(getattr(self.params.P, "std_anchor_catalog", "") or ""))
        self.param_anchor_catalog.setPlaceholderText(
            "VizieR catalog id — e.g. J/AJ/106/181 (M67 UBVRI, Montgomery+1993)")
        self.param_anchor_catalog.setToolTip(
            "대상 시야를 덮는 표준성 측광 카탈로그의 VizieR ID. 직접 "
            "<band>mag 컬럼 또는 Vmag+색(B-V, U-B, V-R, V-I)이 있어야 하며, "
            "Gaia에서 유도된 카탈로그는 쓰면 안 됩니다(독립성 상실).")
        anchor_form.addRow("VizieR catalog:", self.param_anchor_catalog)

        self.param_anchor_find = QPushButton("Find catalogs for this field (VizieR)")
        self.param_anchor_find.setToolTip(
            "대상 이름·좌표로 VizieR를 검색해 이 시야를 덮는 표준 측광 "
            "카탈로그 후보를 찾습니다 (네트워크, 10~60초). 후보를 고르면 "
            "위 칸이 채워집니다.")
        self.param_anchor_find.clicked.connect(self._discover_anchor_catalogs)
        anchor_form.addRow("", self.param_anchor_find)

        self.param_anchor_candidates = QComboBox()
        self.param_anchor_candidates.setVisible(False)
        self.param_anchor_candidates.activated.connect(
            self._pick_anchor_candidate)
        anchor_form.addRow("", self.param_anchor_candidates)

        self.param_anchor_radius = QDoubleSpinBox()
        self.param_anchor_radius.setRange(0.1, 10.0)
        self.param_anchor_radius.setDecimals(1)
        self.param_anchor_radius.setValue(
            float(getattr(self.params.P, "std_anchor_match_radius", 1.5)))
        anchor_form.addRow("Match radius (arcsec):", self.param_anchor_radius)

        self.param_anchor_min_stars = QSpinBox()
        self.param_anchor_min_stars.setRange(5, 1000)
        self.param_anchor_min_stars.setValue(
            int(getattr(self.params.P, "std_anchor_min_stars", 20)))
        anchor_form.addRow("Min matched stars:", self.param_anchor_min_stars)

        layout.addWidget(anchor_group)
        layout.addStretch(1)
        add_parameter_reset_button(
            buttons,
            [
                (self.param_match, 1.0),
                (self.param_min_match, 10),
                (self.param_cmd_snr, 50.0),
                (self.param_frame_min, 5),
                (self.param_apply_ext, False),
                (self.param_ext_mode, "absorb"),
                (self.param_clip, 3.0),
                (self.param_iters, 5),
                (self.param_slope, 1.0),
                (self.param_gaia_snr, 20.0),
                (self.param_gi_min, -0.5),
                (self.param_gi_max, 4.5),
                (self.param_anchor_enable, False),
                (self.param_anchor_radius, 1.5),
                (self.param_anchor_min_stars, 20),
            ],
        )
        buttons.accepted.connect(lambda: self.save_parameters(dialog))
        buttons.rejected.connect(dialog.reject)
        dialog.exec_()

    def _discover_anchor_catalogs(self):
        P = self.params.P
        ra = getattr(P, "target_ra_deg", None) or getattr(P, "ra_deg", None)
        dec = getattr(P, "target_dec_deg", None) or getattr(P, "dec_deg", None)
        if ra is None or dec is None:
            QMessageBox.warning(self, "Standard anchor",
                                "대상 좌표(target.ra_deg/dec_deg)가 설정에 없어 "
                                "탐색할 수 없습니다.")
            return
        name = str(getattr(P, "target_name", "") or "").strip()
        if not name:
            QMessageBox.warning(self, "Standard anchor",
                                "대상 이름(target.name)이 설정에 없어 키워드 "
                                "탐색을 할 수 없습니다.")
            return
        self.param_anchor_find.setEnabled(False)
        self.param_anchor_find.setText("Searching VizieR…")
        bands = ["U", "B", "V", "R", "I", "g", "r", "i"]
        self._anchor_discovery = _AnchorDiscoveryWorker(name, ra, dec, bands, self)
        self._anchor_discovery.found.connect(self._anchor_candidates_found)
        self._anchor_discovery.failed.connect(self._anchor_discovery_failed)
        self._anchor_discovery.start()

    def _anchor_discovery_reset_button(self):
        self.param_anchor_find.setEnabled(True)
        self.param_anchor_find.setText("Find catalogs for this field (VizieR)")

    def _anchor_candidates_found(self, candidates):
        self._anchor_discovery_reset_button()
        combo = self.param_anchor_candidates
        combo.clear()
        if not candidates:
            combo.setVisible(False)
            QMessageBox.information(
                self, "Standard anchor",
                "이 시야를 덮는 표준 측광 카탈로그를 VizieR에서 찾지 "
                "못했습니다. Gaia 참조만으로 진행하거나, 표준장 전이 "
                "(같은 밤 표준장 관측)나 분광 [M/H] prior를 고려하세요.")
            return
        self._anchor_candidate_ids = [c.catalog_id for c in candidates]
        for c in candidates:
            field = "시야 내" if c.in_field else "시야 밖?"
            combo.addItem(f"{c.catalog_id} — {'/'.join(c.bands)} ({field}) "
                          f"{c.description[:40]}")
        combo.setVisible(True)
        combo.showPopup()

    def _pick_anchor_candidate(self, index):
        ids = getattr(self, "_anchor_candidate_ids", [])
        if 0 <= index < len(ids):
            self.param_anchor_catalog.setText(ids[index])
            self.param_anchor_enable.setChecked(True)

    def _anchor_discovery_failed(self, message):
        self._anchor_discovery_reset_button()
        QMessageBox.warning(self, "Standard anchor",
                            f"VizieR 탐색 실패: {message}")

    def save_parameters(self, dialog):
        self.params.P.match_tol_px = self.param_match.value()
        self.params.P.min_master_gaia_matches = self.param_min_match.value()
        self.params.P.cmd_snr_calib_min = self.param_cmd_snr.value()
        self.params.P.frame_zp_min_n = self.param_frame_min.value()
        self.params.P.cmd_apply_extinction = self.param_apply_ext.isChecked()
        self.params.P.cmd_extinction_mode = self.param_ext_mode.currentText().strip()
        self.params.P.zp_clip_sigma = self.param_clip.value()
        self.params.P.zp_fit_iters = self.param_iters.value()
        self.params.P.zp_slope_absmax = self.param_slope.value()
        self.params.P.gaia_snr_calib_min = self.param_gaia_snr.value()
        self.params.P.gaia_gi_min = self.param_gi_min.value()
        self.params.P.gaia_gi_max = self.param_gi_max.value()
        self.params.P.std_anchor_enable = self.param_anchor_enable.isChecked()
        self.params.P.std_anchor_catalog = self.param_anchor_catalog.text().strip()
        self.params.P.std_anchor_match_radius = self.param_anchor_radius.value()
        self.params.P.std_anchor_min_stars = self.param_anchor_min_stars.value()
        self._zp_cache_validation_result = None
        self.save_state()
        saved = self.persist_params()
        msg = "Parameters saved to TOML." if saved else "Parameters saved (TOML save failed)."
        QMessageBox.information(dialog, "Success", msg)
        dialog.accept()

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
            return [ZeropointCalibrationWindow._signature_value(v) for v in value]
        if isinstance(value, dict):
            return {
                str(k): ZeropointCalibrationWindow._signature_value(v)
                for k, v in sorted(value.items(), key=lambda item: str(item[0]))
            }
        return str(value)

    @staticmethod
    def _file_signature(path: Path | None) -> dict | None:
        if path is None:
            return None
        try:
            path = Path(path)
            if not path.is_file():
                return None
            stat = path.stat()
            return {
                "path": str(path.resolve()),
                "size": int(stat.st_size),
                "mtime_ns": int(stat.st_mtime_ns),
            }
        except OSError:
            return None

    def _build_zp_output_signature(self) -> dict:
        result_dir = Path(self.params.P.result_dir)
        source_info = resolve_cmd_photometry_input(
            result_dir,
            self.project_state,
        )
        active_phot_dir = Path(source_info["directory"])
        upstream_paths: list[Path] = [
            step5_wcs_dir(result_dir) / "wcs_solve_summary.csv",
            Path(source_info["index_path"]),
        ]
        if source_info["source"] == "psf":
            upstream_paths.extend([
                step7_forced_phot_dir(result_dir) / "photometry_index.csv",
                active_phot_dir / "psf_output_signature.json",
            ])
            active_patterns = ("photometry_*.tsv",)
        else:
            active_patterns = ("photometry_*.tsv", "apcorr_summary.csv")
        for directory, patterns in (
            (active_phot_dir, active_patterns),
            (step9_selection_dir(result_dir), ("*.csv", "*.tsv", "*.json")),
            (tool_extinction_dir(result_dir), ("*.csv", "*.json")),
        ):
            if directory.exists():
                for pattern in patterns:
                    upstream_paths.extend(sorted(directory.glob(pattern)))

        frame_paths: list[Path] = []
        try:
            for filename in self.file_manager.get_file_list():
                frame_paths.append(Path(self.file_manager.get_file_path(filename)))
        except Exception:
            frame_paths = []

        def _unique_signatures(paths: list[Path]) -> list[dict]:
            signatures: list[dict] = []
            seen: set[str] = set()
            for path in paths:
                signature = self._file_signature(path)
                if not signature:
                    continue
                key = signature["path"]
                if key in seen:
                    continue
                seen.add(key)
                signatures.append(signature)
            return sorted(signatures, key=lambda item: item["path"])

        payload = {
            "signature_version": _ZP_SIGNATURE_VERSION,
            "step": "cmd_step10_zeropoint_calibration",
            "photometry_input": {
                "source": source_info["source"],
                "mag_column": source_info["mag_column"],
                "mag_error_column": source_info["mag_error_column"],
            },
            "params": {
                key: self._signature_value(getattr(self.params.P, key, None))
                for key in _ZP_SIGNATURE_PARAMS
            },
            "inputs": {
                "upstream": _unique_signatures(upstream_paths),
                "frames": _unique_signatures(frame_paths),
            },
        }
        encoded = json.dumps(payload, sort_keys=True, ensure_ascii=False, allow_nan=False)
        payload["signature_hash"] = hashlib.sha1(encoded.encode("utf-8")).hexdigest()
        return payload

    def _stored_zp_signature(self) -> dict | None:
        path = step10_zp_dir(self.params.P.result_dir) / _ZP_SIGNATURE_FILE
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
            return data if isinstance(data, dict) else None
        except Exception:
            return None

    def _write_zp_signature(self, signature: dict) -> None:
        out_dir = step10_zp_dir(self.params.P.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        path = out_dir / _ZP_SIGNATURE_FILE
        temp_path = path.with_name(path.name + ".tmp")
        temp_path.write_text(
            json.dumps(signature, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False),
            encoding="utf-8",
        )
        temp_path.replace(path)

    def _remove_zp_signature(self) -> None:
        path = step10_zp_dir(self.params.P.result_dir) / _ZP_SIGNATURE_FILE
        try:
            path.unlink()
        except FileNotFoundError:
            pass
        except OSError:
            pass

    def _zp_cache_status(self) -> tuple[bool, str, dict | None]:
        if self._zp_cache_validation_result is not None:
            return self._zp_cache_validation_result
        stored = self._stored_zp_signature()
        if not stored:
            result = (False, "missing signature", None)
            self._zp_cache_validation_result = result
            return result
        if stored.get("signature_version") != _ZP_SIGNATURE_VERSION:
            result = (False, "signature version mismatch", None)
            self._zp_cache_validation_result = result
            return result
        current = self._build_zp_output_signature()
        if stored.get("signature_hash") != current.get("signature_hash"):
            result = (False, "signature hash mismatch", None)
            self._zp_cache_validation_result = result
            return result
        summary = self._existing_output_summary()
        if not summary:
            result = (False, "calibration output missing or empty", None)
            self._zp_cache_validation_result = result
            return result
        result = (True, "ok", summary)
        self._zp_cache_validation_result = result
        return result

    def run_analysis(self):
        if self.worker and self.worker.isRunning():
            return
        self._refresh_photometry_source_label()
        self.log_text.clear()
        self.progress_label.setText("Starting...")
        self._zp_cache_validation_result = None
        self._current_zp_signature = self._build_zp_output_signature()
        self._remove_zp_signature()

        self.worker = ZeropointCalibrationWorker(
            self.params,
            self.params.P.data_dir,
            self.params.P.result_dir,
            self.params.P.cache_dir,
        )
        self.worker.progress.connect(self.on_progress)
        self.worker.log.connect(self.log)
        self.worker.finished.connect(self.on_finished)
        self.worker.error.connect(self.on_error)

        self.run_bar.set_running(True)
        self.worker.start()
        self.show_log_window()

    def stop_analysis(self):
        if self.worker and self.worker.isRunning():
            self.worker.stop()
            self.log("Stop requested")

    def on_progress(self, current, total, filename):
        self.progress_label.setText(f"{current}/{total} | {filename}")
        if hasattr(self, "worker_panel"):
            pct = int(100 * current / max(1, total))
            self.worker_panel.update_worker(0, filename, f"{current}/{total}", pct)

    def on_finished(self, summary):
        self.run_bar.set_running(False)
        if summary.get("stopped"):
            self.progress_label.setText("Stopped")
            self.log("Analysis stopped")
        else:
            self.progress_label.setText("Done")
            self.log("ZP calibration complete")
            self._refresh_photometry_source_label(summary)
            if self._current_zp_signature and self._existing_output_summary():
                try:
                    self._write_zp_signature(self._current_zp_signature)
                    self.log("[ZP][CACHE] Output signature saved for future reuse.")
                except Exception as exc:
                    self.log(f"[ZP][CACHE] Signature write failed: {exc}")
            self._zp_cache_validation_result = None
            self.save_state()
            self.update_navigation_buttons()
            self.fit_tab.reload()
        self._current_zp_signature = None
        self._cleanup_worker()

    def on_error(self, message):
        self.run_bar.set_running(False)
        self.progress_label.setText("Error")
        self.log(f"ERROR: {message}")
        self._current_zp_signature = None
        self._cleanup_worker()

    def _cleanup_worker(self, timeout_ms=5000):
        if not self.worker:
            return True

        worker = self.worker
        if worker.isRunning():
            try:
                worker.stop()
            except Exception:
                pass
            worker.quit()
            if not worker.wait(int(timeout_ms)):
                self.log("Calibration worker is still running; close is deferred.")
                return False

        try:
            worker.deleteLater()
        except Exception:
            pass
        self.worker = None
        return True

    def closeEvent(self, event):
        if self.worker and self.worker.isRunning():
            self.stop_analysis()
        if not self._cleanup_worker(timeout_ms=10000):
            QMessageBox.warning(
                self,
                "Background Task Running",
                "Calibration worker is still stopping. Please wait a few seconds and close again.",
            )
            event.ignore()
            return
        super().closeEvent(event)

    def _existing_output_summary(self) -> dict | None:
        out_dir = step10_zp_dir(self.params.P.result_dir)
        wide_cmd = out_dir / "median_by_ID_filter_wide_cmd.csv"
        wide = out_dir / "median_by_ID_filter_wide.csv"
        if not wide_cmd.exists() and not wide.exists():
            legacy_cmd = self.params.P.result_dir / "median_by_ID_filter_wide_cmd.csv"
            legacy = self.params.P.result_dir / "median_by_ID_filter_wide.csv"
            wide_cmd = legacy_cmd
            wide = legacy
        main_path = wide_cmd if wide_cmd.exists() else wide if wide.exists() else None
        if main_path is None:
            return None

        n_sources = 0
        try:
            n_sources = len(pd.read_csv(main_path, nrows=1000000))
        except Exception:
            return None
        if n_sources <= 0:
            return None

        coeff_path = out_dir / "zp_fit_coefficients.csv"
        frame_path = out_dir / "frame_zeropoint.csv"
        cal_path = out_dir / "gaia_sdss_calibrator_by_ID.csv"
        n_coeff = 0
        n_frames = 0
        n_cal = 0
        try:
            if coeff_path.exists():
                n_coeff = len(pd.read_csv(coeff_path))
        except Exception:
            n_coeff = 0
        try:
            if frame_path.exists():
                n_frames = len(pd.read_csv(frame_path))
        except Exception:
            n_frames = 0
        try:
            if cal_path.exists():
                n_cal = len(pd.read_csv(cal_path))
        except Exception:
            n_cal = 0

        return {
            "main_path": str(main_path),
            "n_sources": int(n_sources),
            "n_coeff": int(n_coeff),
            "n_frames": int(n_frames),
            "n_calibrators": int(n_cal),
        }

    def _try_load_existing_results(self) -> bool:
        valid, reason, summary = self._zp_cache_status()
        if not valid or not summary:
            if self._existing_output_summary():
                try:
                    self.log(f"[ZP][CACHE] Previous output not restored ({reason}).")
                except Exception:
                    pass
            return False
        try:
            self.fit_tab.reload(self.params.P.result_dir)
        except Exception:
            pass
        try:
            export_zp_qc_products(step10_zp_dir(self.params.P.result_dir), self.log)
            export_cmd_qc_products(step10_zp_dir(self.params.P.result_dir), self.log)
            export_gaia_cmd_comparison_products(step10_zp_dir(self.params.P.result_dir), self.log)
        except Exception:
            pass
        parts = [f"{summary.get('n_sources', 0)} sources"]
        if summary.get("n_frames", 0):
            parts.append(f"{summary['n_frames']} frame ZPs")
        if summary.get("n_coeff", 0):
            parts.append(f"{summary['n_coeff']} fit coeffs")
        if summary.get("n_calibrators", 0):
            parts.append(f"{summary['n_calibrators']} calibrators")
        self.progress_label.setText("Loaded previous ZP calibration (" + ", ".join(parts) + ")")
        try:
            self.log("[ZP][CACHE] Loaded previous Step 10 ZP calibration from disk.")
        except Exception:
            pass
        self.update_navigation_buttons()
        return True

    def validate_step(self) -> bool:
        valid, _, _ = self._zp_cache_status()
        return valid

    def save_state(self):
        state_data = {
            "match_tol_px": getattr(self.params.P, "match_tol_px", 5.0),
            "min_master_gaia_matches": getattr(self.params.P, "min_master_gaia_matches", 10),
            "cmd_snr_calib_min": getattr(self.params.P, "cmd_snr_calib_min", 20.0),
            "frame_zp_min_n": getattr(self.params.P, "frame_zp_min_n", 5),
            "cmd_apply_extinction": getattr(self.params.P, "cmd_apply_extinction", False),
            "cmd_extinction_mode": getattr(self.params.P, "cmd_extinction_mode", "absorb"),
            "zp_clip_sigma": getattr(self.params.P, "zp_clip_sigma", 3.0),
            "zp_fit_iters": getattr(self.params.P, "zp_fit_iters", 5),
            "zp_slope_absmax": getattr(self.params.P, "zp_slope_absmax", 1.0),
            "gaia_snr_calib_min": getattr(self.params.P, "gaia_snr_calib_min", 20.0),
            "gaia_gi_min": getattr(self.params.P, "gaia_gi_min", -0.5),
            "gaia_gi_max": getattr(self.params.P, "gaia_gi_max", 4.5),
            "gaia_zp_slope_absmax": getattr(self.params.P, "gaia_zp_slope_absmax", 1.0),
            "gaia_color_slope_absmax": getattr(self.params.P, "gaia_color_slope_absmax", 2.0),
        }
        self.project_state.store_step_data("zeropoint_calibration", state_data)

    def restore_state(self):
        state_data = self.project_state.get_step_data("zeropoint_calibration")
        if not state_data:
            state_data = self.project_state.get_step_data("cmd_analysis")
        if state_data:
            for key, val in state_data.items():
                if key == "pixel_scale_arcsec":
                    continue
                if hasattr(self.params.P, key):
                    setattr(self.params.P, key, val)
        self._try_load_existing_results()
