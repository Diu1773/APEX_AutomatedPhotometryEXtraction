"""
Step 12: Isochrone Model

Extended with automatic isochrone fitting:
- AutoFit mode: Global search + local refinement for initial parameter estimation
"""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import matplotlib as mpl
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from matplotlib.colors import LinearSegmentedColormap, Normalize
from scipy.spatial import cKDTree

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel, QGroupBox, QMessageBox,
    QTextEdit, QFormLayout, QDoubleSpinBox, QComboBox,
    QCheckBox, QLineEdit, QWidget, QFileDialog, QProgressBar,
    QApplication, QSlider, QSpinBox,
    QTabWidget, QSizePolicy, QScrollArea, QFrame
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal, QTimer

from apex.gui.workflow.step_window_base import StepWindowBase
from apex.gui.layout_rules import scroll_wrap
from apex.gui.theme import Tokens, refresh, style_button
from apex.utils.step_paths_cmd import step9_selection_dir, step10_zp_dir, step12_iso_dir
from apex.utils.common_helpers import format_cmd_title, target_display_name
from apex.utils.cmd_gaia_enrichment import merge_gaia_columns_from_catalog
from apex.utils.gaia_transforms import teff_from_color
from apex.utils.parsec_columns import read_parsec_header
from apex.utils.photometry_provenance import (
    format_photometry_provenance,
    summarize_photometry_table,
)

# Extinction coefficients R = A/E(B-V) (Cardelli+1989 / Fitzpatrick+1999)
EXTINCTION_R: dict[str, float] = {
    # Johnson-Cousins
    "U": 4.902, "B": 4.035, "V": 3.116, "R": 2.634, "I": 1.903,
    # SDSS (u = metallicity-sensitive blue band; A_u/E(B-V), SF2011)
    "u": 4.239, "g": 3.303, "r": 2.285, "i": 1.698, "z": 1.263,
}
# Backwards compat alias
SDSS_R = EXTINCTION_R

# Default color pairs/bands (overridden dynamically from loaded data at runtime)
_DEFAULT_COLOR_PAIRS = [("g", "r"), ("g", "i"), ("r", "i"),
                        ("B", "V"), ("V", "R"), ("V", "I"), ("B", "R")]
_DEFAULT_MAG_BANDS   = ["g", "r", "i", "V", "B", "R", "I"]


def _bands_from_df(df: "pd.DataFrame") -> tuple[list[tuple[str, str]], list[str]]:
    """Return (adjacent color_pairs, mag_bands) from calibrated columns when available."""
    from apex.utils.gaia_transforms import build_color_pairs, filter_bands_from_columns
    bands = (
        filter_bands_from_columns(df.columns, "mag_cal_")
        or filter_bands_from_columns(df.columns, "mag_std_")
        or filter_bands_from_columns(df.columns, "mag_inst_")
    )
    if not bands:
        return _DEFAULT_COLOR_PAIRS[:3], _DEFAULT_MAG_BANDS[:3]
    return build_color_pairs(bands, adjacent_only=True), bands


def _preferred_mag_col(columns, band: str) -> tuple[str | None, str]:
    for prefix, label in (
        ("mag_cal_", "Calibrated"),
        ("mag_std_", "Standardized"),
        ("mag_inst_", "Instrumental"),
    ):
        col = f"{prefix}{band}"
        if col in columns:
            return col, label
    return None, "Missing"


def _preferred_err_col(columns, band: str) -> str | None:
    for prefix in ("mag_cal_err_", "mag_std_err_", "mag_inst_err_"):
        col = f"{prefix}{band}"
        if col in columns:
            return col
    return None


def _extinction_magnitude_shift(
    color_excess: float,
    band_1: str,
    band_2: str,
    magnitude_band: str,
) -> float:
    """Convert E(band_1-band_2) to extinction in the magnitude band."""
    r1 = EXTINCTION_R.get(band_1)
    r2 = EXTINCTION_R.get(band_2)
    rmag = EXTINCTION_R.get(magnitude_band)
    if r1 is None or r2 is None or rmag is None:
        return 0.0
    denominator = float(r1) - float(r2)
    if abs(denominator) < 1e-12:
        return 0.0
    return float(rmag) * float(color_excess) / denominator


def _segmented_isochrone_line(
    color: np.ndarray,
    magnitude: np.ndarray,
    *,
    min_jump: float = 0.5,
    jump_factor: float = 8.0,
    local_window: int = 8,
) -> tuple[np.ndarray, np.ndarray]:
    """Insert NaN breaks where adjacent model rows are discontinuous.

    Isochrone tables can append remnant or separately computed evolutionary
    blocks after the visible stellar sequence. Connecting those rows with a
    normal line creates long diagonals across the CMD. The original points
    remain untouched for fitting; only the display line is segmented.
    """
    x = np.asarray(color, dtype=float)
    y = np.asarray(magnitude, dtype=float)
    if x.shape != y.shape:
        raise ValueError("Isochrone color and magnitude arrays must match")
    if x.ndim != 1 or len(x) < 2:
        return x.copy(), y.copy()

    steps = np.hypot(np.diff(x), np.diff(y))
    break_after = ~(
        np.isfinite(x[:-1])
        & np.isfinite(y[:-1])
        & np.isfinite(x[1:])
        & np.isfinite(y[1:])
    )

    for index, step in enumerate(steps):
        if break_after[index] or not np.isfinite(step) or step <= min_jump:
            continue
        start = max(0, index - local_window)
        stop = min(len(steps), index + local_window + 1)
        nearby = np.concatenate((steps[start:index], steps[index + 1:stop]))
        nearby = nearby[np.isfinite(nearby) & (nearby > 0)]
        local_step = float(np.median(nearby)) if len(nearby) else 0.0
        if local_step == 0.0 or step > jump_factor * local_step:
            break_after[index] = True

    if not break_after.any():
        return x.copy(), y.copy()

    extra = int(break_after.sum())
    line_x = np.empty(len(x) + extra, dtype=float)
    line_y = np.empty(len(y) + extra, dtype=float)
    output_index = 0
    for index in range(len(x)):
        line_x[output_index] = x[index]
        line_y[output_index] = y[index]
        output_index += 1
        if index < len(break_after) and break_after[index]:
            line_x[output_index] = np.nan
            line_y[output_index] = np.nan
            output_index += 1
    return line_x, line_y


def _legacy_b_calibration_reason(input_dir: str | Path, bands) -> str | None:
    """Return a rerun message when Johnson B still uses legacy calibration."""
    if "B" not in set(bands):
        return None
    coeff_path = Path(input_dir) / "zp_fit_coefficients.csv"
    if not coeff_path.exists():
        return None
    try:
        coeff = pd.read_csv(coeff_path)
    except Exception:
        return None
    if "filter" not in coeff.columns:
        return None
    b_rows = coeff[coeff["filter"].astype(str) == "B"]
    if b_rows.empty:
        return None
    if "ref_source" not in b_rows.columns:
        return (
            "The Johnson B calibration predates the corrected Gaia-to-B "
            "transformation. Rerun Step 10 before using a B-based CMD."
        )
    sources = b_rows["ref_source"].fillna("").astype(str)
    if not sources.str.startswith("Pancino+2022").all():
        return (
            "Johnson B was calibrated with a legacy Gaia transformation. "
            "Rerun Step 10 before using a B-based CMD."
        )
    return None


class IsoCacheWorker(QThread):
    """Build the ``.apex.npy`` sidecar for an isochrone table in the
    background so opening Step 12 doesn't freeze the GUI for the
    ~30 s it takes to parse a 400 MB text file.
    """

    progress_msg = pyqtSignal(str)
    finished_ok = pyqtSignal(str)         # iso_file path
    failed = pyqtSignal(str, str)         # iso_file path, error message

    def __init__(self, iso_file: Path):
        super().__init__()
        self.iso_file = Path(iso_file)
        self._cancel_requested = False

    def cancel(self) -> None:
        self._cancel_requested = True

    def run(self):  # type: ignore[override]
        temp_path = None
        try:
            cache_path = self.iso_file.with_suffix(self.iso_file.suffix + ".apex.npy")
            self.progress_msg.emit(
                f"Parsing isochrone (one-time, ~30 s): {self.iso_file.name}"
            )
            arr = np.loadtxt(self.iso_file, comments="#")
            if self._cancel_requested:
                return
            if arr.ndim == 1:
                arr = arr.reshape(1, -1)
            arr = arr[~np.isnan(arr).any(axis=1)]
            self.progress_msg.emit(f"Writing binary cache: {cache_path.name}")
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            temp_path = cache_path.with_name(
                f"{cache_path.name}.{id(self)}.tmp"
            )
            with temp_path.open("wb") as handle:
                np.save(handle, arr, allow_pickle=False)
            if self._cancel_requested:
                return
            temp_path.replace(cache_path)
        except Exception as exc:
            if not self._cancel_requested:
                self.failed.emit(str(self.iso_file), str(exc))
            return
        finally:
            if temp_path is not None:
                try:
                    temp_path.unlink()
                except FileNotFoundError:
                    pass
                except OSError:
                    pass
        if not self._cancel_requested:
            self.finished_ok.emit(str(self.iso_file))


class McmcFitWorker(QThread):
    """Background worker for the automatic MCMC isochrone fit (Bayesian).

    Runs the full generative-mixture MCMC via the shared, Qt-free service
    :func:`apex.analysis.cmd.isochrone_fit_service.fit_cluster_isochrone` and
    returns posterior summary + paper-quality CMD/corner figures.
    """

    progress = pyqtSignal(float, str)
    finished = pyqtSignal(object)  # IsochroneFitOutput or Exception

    def __init__(self, df, config):
        super().__init__()
        self._df = df
        self._config = config

    def run(self):  # type: ignore[override]
        try:
            from apex.analysis.cmd.isochrone_fit_service import fit_cluster_isochrone
            # make_figures=False: figures are built on the main (GUI) thread in the
            # finished-slot, because Matplotlib/pyplot is not thread-safe.
            out = fit_cluster_isochrone(
                self._df, self._config, make_figures=False,
                progress_cb=lambda f, m: self.progress.emit(float(f), str(m)),
            )
            self.finished.emit(out)
        except Exception as e:  # pragma: no cover - surfaced in the GUI
            self.finished.emit(e)


class FloatSlider(QWidget):
    """QSlider wrapper that emits float values.

    matplotlib Slider objects float around the figure's pixel coords and
    collide with subplots once the canvas resizes; a Qt QSlider mounted in
    its own row never overlaps the plot.
    """

    valueChanged = pyqtSignal(float)

    def __init__(self, label: str, vmin: float, vmax: float, step: float,
                 init_val: float, value_fmt: str = "{:.3f}",
                 snap_values: Optional[np.ndarray] = None, parent=None):
        super().__init__(parent)
        self._scale = max(1.0, 1.0 / max(step, 1e-9))
        self._init = float(init_val)
        self._value_fmt = value_fmt
        # When provided, slider output snaps to the nearest grid value
        # (used for log_age / [Fe/H] so grid lookup hits an actual row
        # in the isochrone table — otherwise get_iso_points() returns
        # empty arrays and the viewer shows "No isochrone points").
        self._snap = None
        if snap_values is not None and len(snap_values) > 0:
            self._snap = np.sort(np.asarray(snap_values, dtype=float))

        layout = QHBoxLayout(self)
        layout.setContentsMargins(2, 0, 2, 0)
        layout.setSpacing(8)

        self.label = QLabel(label)
        self.label.setMinimumWidth(80)
        self.label.setStyleSheet("QLabel { font-weight: bold; }")
        layout.addWidget(self.label)

        self._slider = QSlider(Qt.Horizontal)
        self._slider.setMinimum(int(round(vmin * self._scale)))
        self._slider.setMaximum(int(round(vmax * self._scale)))
        self._slider.setValue(int(round(init_val * self._scale)))
        self._slider.setSingleStep(1)
        self._slider.setPageStep(max(1, int(round((vmax - vmin) * self._scale / 20))))
        self._slider.valueChanged.connect(self._on_changed)
        layout.addWidget(self._slider, stretch=1)

        self.value_label = QLabel()
        self.value_label.setMinimumWidth(70)
        self.value_label.setStyleSheet("QLabel { font-family: monospace; }")
        layout.addWidget(self.value_label)
        self._update_label()

    def _on_changed(self, _: int) -> None:
        self._update_label()
        self.valueChanged.emit(self.value())

    def value(self) -> float:
        raw = float(self._slider.value()) / self._scale
        if self._snap is not None:
            idx = int(np.argmin(np.abs(self._snap - raw)))
            return float(self._snap[idx])
        return raw

    # Compatibility shim for older matplotlib-Slider call sites.
    @property
    def val(self) -> float:
        return self.value()

    def setValue(self, val: float) -> None:
        self._slider.blockSignals(True)
        self._slider.setValue(int(round(float(val) * self._scale)))
        self._slider.blockSignals(False)
        self._update_label()

    def reset(self) -> None:
        self.setValue(self._init)
        self.valueChanged.emit(self.value())

    def _update_label(self) -> None:
        self.value_label.setText(self._value_fmt.format(self.value()))


class IsochroneViewerWindow(QWidget):
    """Interactive isochrone viewer using Qt sliders.

    Sliders live in their own row below the canvas so they cannot overlap
    the CMD/residual subplots when the window is resized.
    """

    def __init__(self, df: pd.DataFrame, iso_raw: np.ndarray, params,
                 parent=None, embedded=False,
                 band_mag: str = "g", band_color: tuple = ("g", "r"),
                 iso_col_1: int = 29, iso_col_2: int = 30,
                 iso_col_mag: Optional[int] = None):
        super().__init__(parent)
        self.df = df
        self.iso_raw = iso_raw
        self.params = params
        self.embedded = bool(embedded)
        self.band_mag = band_mag
        self.band_color = band_color
        self.iso_col_1 = iso_col_1
        self.iso_col_2 = iso_col_2
        self.iso_col_mag = iso_col_1 if iso_col_mag is None else iso_col_mag
        self.color_label = f"{band_color[0]}-{band_color[1]}"
        self.extinction_label = f"E({self.color_label})"

        if not self.embedded:
            self.setWindowTitle("Isochrone Viewer")
            self.resize(1200, 900)
            self.setMinimumSize(900, 700)
        else:
            self.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
            # Force a minimum so step12's CMD-Viewer tab actually allocates
            # enough vertical room for toolbar + canvas + slider row.
            # toolbar 34 + canvas 280 min + sliders 200 + margins ≈ 540.
            self.setMinimumHeight(560)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(5, 5, 5, 5)
        # Slightly more square aspect (was 11×8 → too wide on 1400×900 windows).
        fig_size = (14, 10) if not self.embedded else (9, 7)
        self.figure = Figure(figsize=fig_size)
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        if not self.embedded:
            self.canvas.setMinimumSize(800, 600)
        # Step 12's outer page (Isochrone Source / Band Selection / Source
        # Filters / Tabs / Open Log / nav buttons) leaves the embedded
        # viewer roughly 400–500 px of vertical room.  Use a tight canvas
        # range so the slider row below never overlaps the plot.
        self.canvas.setMinimumHeight(280)
        self.canvas.setMaximumHeight(520)
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas, stretch=1)

        self._build_plot()

    def current_slider_values(self):
        """Return current manual slider state as a fit-initial guess dict."""
        if not all(hasattr(self, attr) for attr in ("s_age", "s_mh", "s_vshift", "s_hshift")):
            return None
        return {
            "log_age": float(self.s_age.val),
            "metallicity": float(self.s_mh.val),
            "distance_mod": float(self.s_vshift.val),
            "extinction_gr": float(self.s_hshift.val),
        }

    def _build_plot(self):
        mpl.rcParams['axes.unicode_minus'] = False
        self.figure.clear()

        teff_vmin = 2400.0
        teff_vmax = 40000.0
        ob_norm = Normalize(vmin=teff_vmin, vmax=teff_vmax, clip=True)

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
        pos = [(t - teff_vmin) / (teff_vmax - teff_vmin) for t, _ in anchors]
        pos[0] = 0.0
        pos[-1] = 1.0

        ob_cmap = LinearSegmentedColormap.from_list(
            "obafgkm_like",
            list(zip(pos, [c for _, c in anchors])),
            N=256
        )
        ob_cmap.set_bad("#777777")

        available_ages = np.unique(self.iso_raw[:, 2])
        available_mhs = np.unique(self.iso_raw[:, 1])

        b1, b2 = self.band_color
        bm = self.band_mag
        col1, color_system1 = _preferred_mag_col(self.df.columns, b1)
        col2, color_system2 = _preferred_mag_col(self.df.columns, b2)
        colm, mag_system = _preferred_mag_col(self.df.columns, bm)

        if col1 is None or col2 is None:
            obs_band1 = np.asarray([], dtype=float)
            obs_band2 = np.asarray([], dtype=float)
            color_system = "Missing"
        else:
            obs_band1 = self.df[col1].to_numpy(float)
            obs_band2 = self.df[col2].to_numpy(float)
            color_system = color_system1 if color_system1 == color_system2 else "Mixed"
        if colm is not None:
            obs_mag = self.df[colm].to_numpy(float)
        else:
            obs_mag = obs_band1.copy()
            mag_system = color_system
        obs_system = color_system if color_system == mag_system else "Mixed"

        obs_color = obs_band1 - obs_band2
        mask = np.isfinite(obs_mag) & np.isfinite(obs_color)
        obs_mag, obs_color = obs_mag[mask], obs_color[mask]
        obs_pts = np.c_[obs_color, obs_mag]
        obs_teff = teff_from_color(obs_color, b1, b2, teff_vmin, teff_vmax)

        # Wider gap between subplots (was wspace=0.2/hspace=0.3) so the
        # right histogram doesn't kiss the CMD's right tick labels and
        # the residual subplot has breathing room below the CMD title.
        gs = self.figure.add_gridspec(
            2, 2, width_ratios=[2.5, 1], height_ratios=[3, 1],
            hspace=0.45, wspace=0.40,
        )
        ax_cmd = self.figure.add_subplot(gs[0, 0])
        ax_hist = self.figure.add_subplot(gs[0, 1])
        ax_res = self.figure.add_subplot(gs[1, 0])

        # Two-line CMD title needs more headroom (top=0.84); left=0.07
        # accommodates a long "M67 - Standardized CMD g vs g-r (N=…)"
        # without clipping the first character.
        self.figure.subplots_adjust(left=0.08, right=0.82, bottom=0.13, top=0.84)

        self.figure.patch.set_facecolor("black")
        for ax in (ax_cmd, ax_hist, ax_res):
            ax.set_facecolor("black")
            for sp in ax.spines.values():
                sp.set_color("white")
            ax.tick_params(colors="white")
            ax.xaxis.label.set_color("white")
            ax.yaxis.label.set_color("white")
            ax.title.set_color("white")

        sc_obs = ax_cmd.scatter(obs_color, obs_mag, s=3, alpha=0.85, linewidths=0, c=obs_teff, cmap=ob_cmap, norm=ob_norm, label=f"Observed ({obs_system})")
        iso_line, = ax_cmd.plot(
            [],
            [],
            color="red",
            linewidth=1.4,
            alpha=0.95,
            label="Isochrone",
            zorder=6,
            solid_capstyle="round",
            solid_joinstyle="round",
        )
        empty_offsets = np.empty((0, 2), dtype=float)

        ax_cmd.invert_yaxis()
        ax_cmd.set_xlabel(f"{obs_system} ({self.color_label})")
        ax_cmd.set_ylabel(f"{obs_system} {bm}")
        ax_cmd.grid(True, linestyle=":", alpha=0.35)
        ax_cmd.legend(loc="upper right")

        res_scat = ax_res.scatter([], [], s=3, alpha=0.75, linewidths=0, color="cyan")
        ax_res.axhline(0, color="white", lw=1, ls="--", alpha=0.6)
        ax_res.set_xlabel(f"{obs_system} {bm}")
        ax_res.set_ylabel("Residual (NN dist in CMD)")

        sm = mpl.cm.ScalarMappable(norm=ob_norm, cmap=ob_cmap)
        sm.set_array([])
        # Manual cax keeps the 'Teff (K) + OBAFGKM-like color' label fully
        # visible — the auto-positioned colorbar was getting shoved off
        # the right edge by the long Teff (e.g. "35000 K (O)") tick labels.
        cax = self.figure.add_axes([0.86, 0.13, 0.018, 0.75])
        cbar = self.figure.colorbar(sm, cax=cax)
        cbar.set_label("Teff (K) + OBAFGKM-like color", color="white", labelpad=14)
        cbar.ax.tick_params(colors="white")
        for sp in cbar.ax.spines.values():
            sp.set_color("white")

        c1 = self.iso_col_1
        c2 = self.iso_col_2
        cm = self.iso_col_mag

        def get_iso_points(age, mh, h_shift, v_shift):
            m = (self.iso_raw[:, 2] == age) & (self.iso_raw[:, 1] == mh)
            filtered = self.iso_raw[m]
            if len(filtered) == 0:
                return np.array([]), np.array([])
            extinction_mag = _extinction_magnitude_shift(
                h_shift,
                b1,
                b2,
                bm,
            )
            mag_model = filtered[:, cm] + v_shift + extinction_mag
            color_model = (filtered[:, c1] - filtered[:, c2]) + h_shift
            return color_model, mag_model

        def style_axis_dark(ax):
            ax.set_facecolor("black")
            for sp in ax.spines.values():
                sp.set_color("white")
            ax.tick_params(colors="white")
            ax.xaxis.label.set_color("white")
            ax.yaxis.label.set_color("white")
            ax.title.set_color("white")

        age_init = float(getattr(self.params.P, "iso_age_init", 9.7))
        mh_init = float(getattr(self.params.P, "iso_mh_init", -0.1))
        if len(available_ages) > 0:
            age_init = float(available_ages[np.argmin(np.abs(available_ages - age_init))])
        if len(available_mhs) > 0:
            mh_init = float(available_mhs[np.argmin(np.abs(available_mhs - mh_init))])

        # Qt sliders mounted in their own row below the canvas.  Each
        # slider snaps to the same grid the old matplotlib widget used.
        age_min_v = float(available_ages.min()) if len(available_ages) else 6.0
        age_max_v = float(available_ages.max()) if len(available_ages) else 10.5
        mh_min_v = float(available_mhs.min()) if len(available_mhs) else -2.2
        mh_max_v = float(available_mhs.max()) if len(available_mhs) else 1.0
        # Step sizes inferred from grid spacing (fall back to sensible defaults).
        age_step = float(np.median(np.diff(available_ages))) if len(available_ages) > 1 else 0.05
        mh_step = float(np.median(np.diff(available_mhs))) if len(available_mhs) > 1 else 0.1
        if not np.isfinite(age_step) or age_step <= 0:
            age_step = 0.05
        if not np.isfinite(mh_step) or mh_step <= 0:
            mh_step = 0.1

        # Pass the discrete isochrone age/[M/H] grid to the sliders so each
        # tick lands on a row that exists in the isochrone table; otherwise
        # interpolation between rows yields empty point lists and the
        # viewer shows "No isochrone points".
        s_age = FloatSlider("log Age", age_min_v, age_max_v, age_step,
                             age_init, value_fmt="{:.2f}",
                             snap_values=available_ages)
        s_mh = FloatSlider("[Fe/H]", mh_min_v, mh_max_v, mh_step,
                            mh_init, value_fmt="{:.2f}",
                            snap_values=available_mhs)
        s_hshift = FloatSlider(self.extinction_label, -0.1, 0.8, 0.0001,
                                float(getattr(self.params.P, "iso_eg_r_init", 0.0033)),
                                value_fmt="{:.4f}")
        s_vshift = FloatSlider("Dist. Mod", 5.0, 20.0, 0.01,
                                float(getattr(self.params.P, "iso_dm_init", 9.46)),
                                value_fmt="{:.2f}")

        sliders_box = QGroupBox()
        sliders_box.setFlat(True)
        sliders_box.setStyleSheet("QGroupBox { border: 0; }")
        sliders_layout = QVBoxLayout(sliders_box)
        sliders_layout.setContentsMargins(4, 4, 4, 4)
        sliders_layout.setSpacing(4)
        sliders_layout.addWidget(s_age)
        sliders_layout.addWidget(s_mh)
        sliders_layout.addWidget(s_hshift)
        sliders_layout.addWidget(s_vshift)

        reset_row = QHBoxLayout()
        reset_row.addStretch()
        self.reset_button = QPushButton("Reset")
        self.reset_button.setFixedWidth(80)
        style_button(self.reset_button)
        reset_row.addWidget(self.reset_button)
        sliders_layout.addLayout(reset_row)
        # _build_plot runs after __init__ set up `self.layout()` (QVBoxLayout
        # containing the canvas).  Add the slider row to that same layout
        # so it sits below the canvas.
        self.layout().addWidget(sliders_box)

        self.s_age = s_age
        self.s_mh = s_mh
        self.s_hshift = s_hshift
        self.s_vshift = s_vshift

        def update(_=None):
            age, mh = s_age.value(), s_mh.value()
            h_s, v_s = s_hshift.value(), s_vshift.value()

            # Keep current manual values in sync for subsequent auto-fit runs.
            self.params.P.iso_age_init = float(age)
            self.params.P.iso_mh_init = float(mh)
            self.params.P.iso_eg_r_init = float(h_s)
            self.params.P.iso_dm_init = float(v_s)

            new_gr, new_g = get_iso_points(age, mh, h_s, v_s)

            if len(new_gr) > 0:
                # Isochrone rows are already ordered along the model sequence.
                # Preserve that order, but break appended/discontinuous model
                # blocks so they do not draw artificial diagonals across CMD.
                line_color, line_mag = _segmented_isochrone_line(new_gr, new_g)
                iso_line.set_data(line_color, line_mag)
            else:
                iso_line.set_data([], [])

            if len(new_gr) > 0 and len(obs_pts) > 0:
                iso_pts = np.c_[new_gr, new_g]
                tree = cKDTree(iso_pts)
                dist, _ = tree.query(obs_pts)

                res_scat.set_offsets(np.c_[obs_mag, dist])
                ax_res.set_xlim(ax_cmd.get_ylim())
                ax_res.set_ylim(0, np.percentile(dist, 95))

                ax_hist.clear()
                style_axis_dark(ax_hist)
                hi = np.percentile(dist, 98)
                ax_hist.hist(dist, bins=30, range=(0, hi), color="deepskyblue", edgecolor="white", alpha=0.75)
                ax_hist.set_title(f"Mean Res: {np.mean(dist):.4f}", color="white", pad=10)
            else:
                res_scat.set_offsets(empty_offsets)
                ax_hist.clear()
                style_axis_dark(ax_hist)
                ax_hist.set_title("No isochrone points", color="white", pad=10)

            base_title = format_cmd_title(
                self.params,
                bm,
                self.color_label,
                system_label=obs_system,
                count_text=f"N={len(obs_mag)}",
            )
            # Split into two lines so the long parameter readout
            # doesn't push the title off the left edge of the canvas.
            ax_cmd.set_title(
                f"{base_title}\n"
                f"log(Age)={age:.2f}, [M/H]={mh:.2f}, DM={v_s:.2f}, "
                f"{self.extinction_label}={h_s:.4f}",
                color="white",
                pad=10,
                fontsize=10,
            )
            self.canvas.draw_idle()

        def reset(_=None):
            s_age.reset()
            s_mh.reset()
            s_hshift.reset()
            s_vshift.reset()
            update()

        s_age.valueChanged.connect(update)
        s_mh.valueChanged.connect(update)
        s_hshift.valueChanged.connect(update)
        s_vshift.valueChanged.connect(update)
        self.reset_button.clicked.connect(reset)

        update()


class IsochroneModelWindow(StepWindowBase):
    """Step 12: Isochrone Model"""

    def __init__(self, params, file_manager, project_state, main_window):
        self.file_manager = file_manager
        self.viewer = None
        self.iso_path_edit = None
        self.iso_status_label = None
        self._iso_cache_worker = None
        self._pending_iso_cache_path = None
        # Data cache – avoids re-reading files on every slider change
        self._cache_key = None          # (iso_path, csv_path)
        self._cached_df = None          # pd.DataFrame
        self._cached_iso_raw = None     # np.ndarray
        self._cached_iso_file = None    # Path
        self._iso_band_columns: dict[str, int] = {}
        self._iso_photometric_system = "Unknown"
        self._iso_header_key = None
        self._pending_band_color_text = None
        self._pending_band_mag_text = None

        super().__init__(
            step_index=11,
            step_name="Isochrone Model",
            params=params,
            project_state=project_state,
            main_window=main_window
        )

        self.setup_step_ui()
        self.restore_state()

    def setup_step_ui(self):
        info = QLabel("Load isochrone data, explore with sliders, or run automatic fitting.")
        info.setProperty("role", "info")
        self.content_layout.addWidget(info)

        self.photometry_source_label = QLabel()
        self.photometry_source_label.setProperty("role", "caption")
        _f = self.photometry_source_label.font(); _f.setBold(True)
        self.photometry_source_label.setFont(_f)
        self.photometry_source_label.setToolTip(
            "Original measurement source used to build the calibrated CMD magnitudes."
        )
        self.content_layout.addWidget(self.photometry_source_label)
        self._refresh_photometry_source_label()

        # === File Selection ===
        file_group = QGroupBox("Isochrone Source")
        file_layout = QVBoxLayout(file_group)
        file_row = QHBoxLayout()
        self.iso_path_edit = QLineEdit()
        self.iso_path_edit.setPlaceholderText("Select isochrone file or folder")
        file_row.addWidget(self.iso_path_edit)
        btn_browse = QPushButton("Browse")
        btn_browse.clicked.connect(self.browse_iso_file)
        file_row.addWidget(btn_browse)
        btn_folder = QPushButton("Folder")
        btn_folder.clicked.connect(self.browse_iso_folder)
        file_row.addWidget(btn_folder)
        file_layout.addLayout(file_row)
        self.iso_status_label = QLabel("Single file mode")
        self.iso_status_label.setProperty("role", "caption")
        file_layout.addWidget(self.iso_status_label)
        self.content_layout.addWidget(file_group)

        # === Filter Selection ===
        filter_group = QGroupBox("Band Selection")
        filter_layout = QHBoxLayout(filter_group)
        filter_layout.addWidget(QLabel("Color (X):"))
        self.color_combo = QComboBox()
        self.color_combo.addItems([f"{a}-{b}" for a, b in _DEFAULT_COLOR_PAIRS[:3]])
        self.color_combo.setCurrentIndex(0)
        self.color_combo.currentIndexChanged.connect(self._on_band_changed)
        filter_layout.addWidget(self.color_combo)
        filter_layout.addWidget(QLabel("Mag (Y):"))
        self.mag_combo = QComboBox()
        self.mag_combo.addItems(_DEFAULT_MAG_BANDS[:3])
        self.mag_combo.setCurrentIndex(0)
        self.mag_combo.currentIndexChanged.connect(self._on_band_changed)
        filter_layout.addWidget(self.mag_combo)
        filter_layout.addStretch()
        self.content_layout.addWidget(filter_group)

        # === Source Filters ===
        sf_group = QGroupBox("Source Filters")
        sf_layout = QHBoxLayout(sf_group)

        # Parallax filter
        self.plx_check = QCheckBox("Parallax filter")
        self.plx_check.setChecked(False)
        self.plx_check.setToolTip("Filter sources by Gaia parallax range.")
        sf_layout.addWidget(self.plx_check)

        self.plx_min_spin = QDoubleSpinBox()
        self.plx_min_spin.setRange(-5.0, 20.0)
        self.plx_min_spin.setDecimals(3)
        self.plx_min_spin.setSingleStep(0.05)
        self.plx_min_spin.setValue(-0.5)
        self.plx_min_spin.setSuffix(" mas")
        sf_layout.addWidget(self.plx_min_spin)

        sf_layout.addWidget(QLabel("–"))

        self.plx_max_spin = QDoubleSpinBox()
        self.plx_max_spin.setRange(-5.0, 20.0)
        self.plx_max_spin.setDecimals(3)
        self.plx_max_spin.setSingleStep(0.05)
        self.plx_max_spin.setValue(0.5)
        self.plx_max_spin.setSuffix(" mas")
        sf_layout.addWidget(self.plx_max_spin)

        sf_layout.addSpacing(20)

        # ROI filter
        self.roi_check = QCheckBox("ROI filter")
        self.roi_check.setChecked(False)
        self.roi_check.setEnabled(False)
        self.roi_check.setToolTip(
            "Filter sources by the spatial ROI circle saved in Step 9.\n"
            "Enable Step 9's ROI tool first."
        )
        sf_layout.addWidget(self.roi_check)

        self.roi_label = QLabel("(no ROI)")
        self.roi_label.setProperty("status", "idle")
        _f = self.roi_label.font(); _f.setItalic(True); self.roi_label.setFont(_f)
        sf_layout.addWidget(self.roi_label)

        sf_layout.addSpacing(20)

        # SNR display filter
        self.snr_display_check = QCheckBox("SNR >=")
        self.snr_display_check.setChecked(True)
        self.snr_display_check.setToolTip("Filter CMD display sources by minimum SNR.")
        sf_layout.addWidget(self.snr_display_check)

        self.snr_display_spin = QDoubleSpinBox()
        self.snr_display_spin.setRange(0.0, 200.0)
        self.snr_display_spin.setDecimals(1)
        self.snr_display_spin.setSingleStep(1.0)
        self.snr_display_spin.setValue(20.0)
        sf_layout.addWidget(self.snr_display_spin)

        sf_layout.addStretch()
        self.content_layout.addWidget(sf_group)

        # Internal ROI state
        self._roi_data: dict | None = None
        self._parallax_auto_range_pending = False
        self._parallax_unavailable_logged = False
        self._parallax_range_initialized = False
        self._load_roi_data()

        # Connect filter signals
        self.plx_check.stateChanged.connect(self._on_parallax_filter_changed)
        self.plx_min_spin.valueChanged.connect(self._on_source_filter_changed)
        self.plx_max_spin.valueChanged.connect(self._on_source_filter_changed)
        self.roi_check.stateChanged.connect(self._on_source_filter_changed)
        self.snr_display_check.stateChanged.connect(self._on_source_filter_changed)
        self.snr_display_spin.valueChanged.connect(self._on_source_filter_changed)

        # === Tabs ===
        self.tabs = QTabWidget()
        self.content_layout.addWidget(self.tabs, stretch=1)

        # --- Tab: Auto-fit (MCMC) — Bayesian, paper-quality figures ---
        self._build_mcmc_autofit_tab()

        # --- Tab 2: CMD Viewer (default tab) ---
        # Wrap in QScrollArea — the embedded viewer needs ~560px of
        # vertical room for toolbar + canvas + slider row; on a 900px
        # Step 12 window the tab area is only ~500px, so a scrollbar
        # keeps the slider row reachable without forcing the user to
        # enlarge the window.
        manual_tab = QWidget()
        manual_outer = QVBoxLayout(manual_tab)
        manual_outer.setContentsMargins(0, 0, 0, 0)
        manual_scroll = QScrollArea()
        manual_scroll.setWidgetResizable(True)
        manual_scroll.setFrameShape(QFrame.NoFrame)
        manual_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAsNeeded)
        manual_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarAsNeeded)

        manual_inner = QWidget()
        manual_layout = QVBoxLayout(manual_inner)
        manual_layout.setContentsMargins(6, 6, 6, 6)

        self.viewer_container = QWidget()
        self.viewer_container.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.viewer_layout = QVBoxLayout(self.viewer_container)
        self.viewer_layout.setContentsMargins(0, 0, 0, 0)
        self.viewer_layout.setSpacing(0)
        manual_layout.addWidget(self.viewer_container, stretch=1)
        self.viewer_placeholder = None
        self._show_viewer_placeholder(
            "Select an isochrone file to render CMD + sliders.\n"
            "Auto-fit (MCMC) results and sliders share this viewer."
        )

        manual_scroll.setWidget(manual_inner)
        manual_outer.addWidget(manual_scroll)
        # Workflow order: fit first (Auto-fit tab, index 0, default), then
        # inspect the result on the CMD Viewer.
        self.cmd_viewer_tab_index = self.tabs.addTab(manual_tab, "CMD Viewer")
        self.tabs.setCurrentIndex(self.mcmc_tab_index)
        # Lazy load: defer the (slow) CMD Viewer first render until the
        # user actually opens that tab.  See _on_tab_changed below.
        self._cmd_viewer_pending = False
        self.tabs.currentChanged.connect(self._on_tab_changed)

        # --- Log Window ---
        log_row = QHBoxLayout()
        btn_log = QPushButton("Open Log")
        style_button(btn_log, "ghost")
        btn_log.clicked.connect(self.show_log_window)
        log_row.addWidget(btn_log)
        log_row.addStretch()
        self.content_layout.addLayout(log_row)

        self.log_window = QWidget(self, Qt.Window)
        self.log_window.setWindowTitle("Isochrone Log")
        self.log_window.resize(700, 350)
        log_layout = QVBoxLayout(self.log_window)
        self.log_text = QTextEdit()
        self.log_text.setReadOnly(True)
        self.log_text.setObjectName("Log")
        log_layout.addWidget(self.log_text)

        # Internal state
        self.cmd_df: Optional[pd.DataFrame] = None

    def log(self, message: str):
        timestamp = time.strftime("%H:%M:%S")
        self.log_text.append(f"[{timestamp}] {message}")

    # =========================================================================
    # Source Filters
    # =========================================================================

    def _set_roi_state(self, status: str, *, italic: bool) -> None:
        """Colour the ROI readout by state, italic while nothing is loaded."""
        self.roi_label.setProperty("status", status)
        font = self.roi_label.font()
        font.setItalic(italic)
        self.roi_label.setFont(font)
        refresh(self.roi_label)

    def _load_roi_data(self):
        """Try to load ROI circle from step10's saved JSON."""
        try:
            roi_path = step9_selection_dir(self.params.P.result_dir) / "cmd_roi.json"
            if roi_path.exists():
                self._roi_data = json.loads(roi_path.read_text())
                ra  = self._roi_data.get("ra_deg", 0.0)
                dec = self._roi_data.get("dec_deg", 0.0)
                r   = self._roi_data.get("radius_arcsec", 0.0)
                self.roi_label.setText(f"RA={ra:.4f}° Dec={dec:.4f}° r={r:.1f}\"")
                self._set_roi_state("ok", italic=False)
                self.roi_check.setEnabled(True)
            else:
                self._roi_data = None
                self.roi_label.setText("(no ROI saved)")
                self._set_roi_state("idle", italic=True)
                self.roi_check.setEnabled(False)
        except Exception:
            self._roi_data = None
            self.roi_check.setEnabled(False)

    def _enrich_cmd_gaia_columns(self, df: pd.DataFrame) -> pd.DataFrame:
        if df is None or df.empty:
            return df
        needed = [
            "source_id", "gaia_source_id",
            "pmra", "pmdec", "parallax",
            "pmra_error", "pmdec_error", "parallax_error",
            "ruwe", "visibility_periods_used",
            "gaia_pmem", "pmem_gaia", "membership_prob_gaia",
        ]
        before_plx = 0
        if "parallax" in df.columns:
            before_plx = int(pd.to_numeric(df["parallax"], errors="coerce").notna().sum())
        try:
            out = merge_gaia_columns_from_catalog(df, Path(self.params.P.result_dir), needed)
        except Exception as exc:
            self.log(f"[Source filters] Gaia enrichment skipped: {exc}")
            return df
        after_plx = 0
        if "parallax" in out.columns:
            after_plx = int(pd.to_numeric(out["parallax"], errors="coerce").notna().sum())
        if after_plx > before_plx:
            self.log(f"[Source filters] Gaia astrometry attached: parallax={after_plx}/{len(out)}")
            self._parallax_unavailable_logged = False
        return out

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

    def _auto_set_parallax_range(self, plx: np.ndarray) -> bool:
        finite = np.isfinite(plx)
        if int(finite.sum()) == 0:
            return False
        vals = plx[finite]
        center = float(np.nanmedian(vals))
        mad = float(np.nanmedian(np.abs(vals - center)))
        robust_sigma = 1.4826 * mad if np.isfinite(mad) and mad > 0 else 0.0
        half_width = max(0.5, 4.0 * robust_sigma)
        half_width = min(5.0, half_width)
        self._set_parallax_range(center - half_width, center + half_width)
        return True

    def _initialize_parallax_range(self, df: pd.DataFrame, force: bool = False) -> bool:
        if self._parallax_range_initialized and not force:
            return False
        if df is None or df.empty or "parallax" not in df.columns:
            return False
        plx = pd.to_numeric(df["parallax"], errors="coerce").to_numpy(float)
        if self._auto_set_parallax_range(plx):
            self._parallax_range_initialized = True
            return True
        return False

    def _on_parallax_filter_changed(self):
        self._parallax_auto_range_pending = self.plx_check.isChecked()
        self._on_source_filter_changed()

    def _apply_source_filters(self, df: pd.DataFrame) -> pd.DataFrame:
        """Return a filtered copy of df based on parallax and ROI settings."""
        if df is None:
            return df
        mask = np.ones(len(df), dtype=bool)

        # Parallax filter
        if self.plx_check.isChecked():
            if "parallax" not in df.columns:
                self.plx_check.blockSignals(True)
                self.plx_check.setChecked(False)
                self.plx_check.blockSignals(False)
                if not self._parallax_unavailable_logged:
                    self.log(
                        "[Source filters] Parallax filter disabled: parallax column is unavailable. "
                        "Rerun Step 10 or ensure Step 5 Gaia cache exists."
                    )
                    self._parallax_unavailable_logged = True
                self._parallax_auto_range_pending = False
            else:
                plx = pd.to_numeric(df["parallax"], errors="coerce").to_numpy(float)
                plx_min = float(self.plx_min_spin.value())
                plx_max = float(self.plx_max_spin.value())
                plx_mask = np.isfinite(plx) & (plx >= plx_min) & (plx <= plx_max)
                default_range = abs(plx_min + 0.5) < 1e-6 and abs(plx_max - 0.5) < 1e-6
                if int(np.isfinite(plx).sum()) > 0 and int(plx_mask.sum()) < 5 and (
                    self._parallax_auto_range_pending or default_range
                ):
                    if self._auto_set_parallax_range(plx):
                        plx_min = float(self.plx_min_spin.value())
                        plx_max = float(self.plx_max_spin.value())
                        plx_mask = np.isfinite(plx) & (plx >= plx_min) & (plx <= plx_max)
                        self.log(
                            "[Source filters] Parallax range auto-centered: "
                            f"{plx_min:.3f}..{plx_max:.3f} mas ({int(plx_mask.sum())} stars)"
                        )
                self._parallax_auto_range_pending = False
                mask &= plx_mask

        # ROI filter
        if self.roi_check.isChecked() and self._roi_data is not None:
            roi_ra       = float(self._roi_data["ra_deg"])
            roi_dec      = float(self._roi_data["dec_deg"])
            roi_r_arcsec = float(self._roi_data["radius_arcsec"])
            if "ra_deg" in df.columns and "dec_deg" in df.columns:
                ra  = pd.to_numeric(df["ra_deg"],  errors="coerce").to_numpy(float)
                dec = pd.to_numeric(df["dec_deg"], errors="coerce").to_numpy(float)
                cos_dec = np.cos(np.radians(roi_dec))
                d_ra    = (ra  - roi_ra)  * cos_dec * 3600.0
                d_dec   = (dec - roi_dec) * 3600.0
                dist    = np.sqrt(d_ra**2 + d_dec**2)
                mask &= np.isfinite(dist) & (dist <= roi_r_arcsec)

        # SNR display filter (applies to all available snr_* columns)
        if self.snr_display_check.isChecked():
            snr_cut = float(self.snr_display_spin.value())
            if snr_cut > 0:
                for sc in [c for c in df.columns if c.startswith("snr_")]:
                    sv = pd.to_numeric(df[sc], errors="coerce").to_numpy(float)
                    mask &= ~(np.isfinite(sv) & (sv < snr_cut))

        if not mask.all():
            df = df[mask].reset_index(drop=True)
        return df

    def _on_source_filter_changed(self):
        """Called when any source filter widget changes."""
        self.refresh_cmd_viewer(show_error=False)

    def _get_iso_path(self) -> str:
        iso_path = ""
        if self.iso_path_edit is not None:
            iso_path = self.iso_path_edit.text().strip()
        if not iso_path:
            iso_path = str(getattr(self.params.P, "iso_file_path", ""))
        return iso_path

    # ===================================================================
    # Auto-fit (MCMC) tab — Bayesian fit + paper-quality figures
    # ===================================================================
    def _build_mcmc_autofit_tab(self):
        """Build the 'Auto-fit (MCMC)' tab (separate from the manual fit tab)."""
        tab = QWidget()
        outer = QVBoxLayout(tab)
        outer.setContentsMargins(6, 6, 6, 6)

        intro = QLabel(
            "<b>Automatic Bayesian isochrone fit</b> (generative-mixture MCMC) → "
            "paper-quality CMD + corner figures.<br>"
            "<b>Filters / colours:</b> the <i>colours</i> row shows a checkbox per "
            "colour found in your data (e.g. BVR → <tt>B-V</tt>, <tt>V-R</tt>) — tick "
            "the ones to fit, untick for a single colour. All ticked = use every "
            "filter.<br>"
            "<b>What each colour does:</b> ONLY a blue/UV colour (<tt>u-g</tt>, "
            "<tt>U-B</tt>, Strömgren m1) constrains <b>[M/H]</b> — red colours "
            "(r-i, i-z) mainly add temperature/age. So <tt>gri</tt>/<tt>BVR</tt> alone "
            "≈ an <b>age</b> meter; for metallicity you need a blue band or an [M/H] prior.<br>"
            "<b>Depth caveat (age):</b> a deep, dense faint main sequence dilutes the "
            "age constraint — which lives in the sparse turn-off — so broadband ages can "
            "read ~10–20% young for old clusters. A membership-clean, turn-off-sampled CMD "
            "fits best (validated on M67/NGC 6811, docs §10.8–10.9).<br>"
            "<b>Priors (how 4D is made reliable):</b> distance ← Gaia parallax "
            "(auto ✓); [M/H] ← spectroscopy; E(B−V) ← dust map. One-stop sources: "
            "Cantat-Gaudin+2020 / Dias+2021 / WEBDA cluster catalogues. Set σ to the "
            "real measurement precision (too loose can't break the degeneracy). "
            "Details: docs §10."
        )
        intro.setWordWrap(True)
        intro.setTextFormat(Qt.RichText)
        intro.setProperty("role", "caption")
        outer.addWidget(intro)

        controls = QHBoxLayout()

        # -- data / sampler options --
        opt_group = QGroupBox("Fit options")
        opt_form = QFormLayout(opt_group)
        self.mcmc_membership_chk = QCheckBox("Gaia PM+parallax membership filter")
        self.mcmc_membership_chk.setChecked(True)
        opt_form.addRow(self.mcmc_membership_chk)
        # Colour selection: one checkbox per colour available in the data
        # (auto-detected). Explicit — tick the colours to fit, untick for a single
        # colour. Populated when the tab is opened (and via the 'detect' button).
        self.mcmc_color_checks: dict = {}
        colors_row = QHBoxLayout()
        colors_row.setContentsMargins(0, 0, 0, 0)
        self._mcmc_colors_holder = QWidget()
        self.mcmc_colors_layout = QHBoxLayout(self._mcmc_colors_holder)
        self.mcmc_colors_layout.setContentsMargins(0, 0, 0, 0)
        self.mcmc_colors_layout.addWidget(QLabel("(open tab to detect)"))
        btn_detect_colors = QPushButton("detect")
        btn_detect_colors.setMaximumWidth(60)
        btn_detect_colors.setToolTip("Re-scan the data for available colours and rebuild the checkboxes.")
        btn_detect_colors.clicked.connect(self._mcmc_refresh_color_checks)
        colors_row.addWidget(self._mcmc_colors_holder, 1)
        colors_row.addWidget(btn_detect_colors)
        opt_form.addRow("colours:", self._wrap(colors_row))
        self.mcmc_parallax_chk = QCheckBox("Gaia parallax distance prior (auto)")
        self.mcmc_parallax_chk.setChecked(True)
        self.mcmc_parallax_chk.setToolTip(
            "Derive a Gaussian (m-M)0 prior from the membership clump's median Gaia "
            "parallax. This external distance breaks the gri degeneracy and is what "
            "recovers the correct age + distance (see docs sec 10.6).")
        opt_form.addRow(self.mcmc_parallax_chk)
        self.mcmc_walkers_spin = QSpinBox(); self.mcmc_walkers_spin.setRange(8, 256)
        self.mcmc_walkers_spin.setValue(32)
        opt_form.addRow("walkers:", self.mcmc_walkers_spin)
        self.mcmc_steps_spin = QSpinBox(); self.mcmc_steps_spin.setRange(200, 20000)
        self.mcmc_steps_spin.setSingleStep(200); self.mcmc_steps_spin.setValue(2000)
        opt_form.addRow("steps:", self.mcmc_steps_spin)
        self.mcmc_burn_spin = QSpinBox(); self.mcmc_burn_spin.setRange(50, 5000)
        self.mcmc_burn_spin.setSingleStep(100); self.mcmc_burn_spin.setValue(600)
        opt_form.addRow("burn-in:", self.mcmc_burn_spin)
        self.mcmc_maxstars_spin = QSpinBox(); self.mcmc_maxstars_spin.setRange(50, 5000)
        self.mcmc_maxstars_spin.setSingleStep(50); self.mcmc_maxstars_spin.setValue(500)
        opt_form.addRow("max stars:", self.mcmc_maxstars_spin)
        # Age search range (Gyr) — self-contained so this tab does NOT depend on
        # the old 'Auto Fit' tab. Set it to bracket the cluster's expected age
        # (too narrow/low rails the fit; e.g. NGC 6811 ≈ 1 Gyr).
        age_row = QHBoxLayout(); age_row.setContentsMargins(0, 0, 0, 0)
        self.mcmc_age_min = QDoubleSpinBox(); self.mcmc_age_min.setRange(0.001, 15.0)
        self.mcmc_age_min.setDecimals(2); self.mcmc_age_min.setSingleStep(0.1); self.mcmc_age_min.setValue(0.20)
        self.mcmc_age_max = QDoubleSpinBox(); self.mcmc_age_max.setRange(0.002, 15.0)
        self.mcmc_age_max.setDecimals(2); self.mcmc_age_max.setSingleStep(0.5); self.mcmc_age_max.setValue(6.0)
        age_row.addWidget(self.mcmc_age_min); age_row.addWidget(QLabel("–")); age_row.addWidget(self.mcmc_age_max)
        opt_form.addRow("age range (Gyr):", self._wrap(age_row))
        controls.addWidget(opt_group)

        # -- priors (Gaussian; tiny sigma == fixed) --
        prior_group = QGroupBox("External priors (Gaussian: mean ± σ; small σ ≈ fixed)")
        prior_form = QFormLayout(prior_group)

        # Shown only when the data has no blue/UV band: without u/U the
        # likelihood cannot constrain [M/H] (documented ~0.2 dex metal-poor
        # floor), so the spectroscopic prior is not optional garnish — it IS
        # the [M/H] measurement. Populated in _mcmc_refresh_color_checks.
        self.mcmc_blue_warning = QLabel(
            "⚠ 데이터에 청색 밴드(u/U)가 없습니다 — 이 색만으로 [M/H]는 "
            "~0.2 dex 금속결핍으로 치우칩니다. 아래 [M/H] prior에 분광 "
            "문헌값을 넣는 것이 공식 경로입니다.")
        self.mcmc_blue_warning.setWordWrap(True)
        self.mcmc_blue_warning.setProperty("banner", "warn")
        self.mcmc_blue_warning.setVisible(False)
        prior_form.addRow(self.mcmc_blue_warning)

        self.mcmc_mh_prior_chk = QCheckBox("[M/H] prior (spectroscopy)")
        self.mcmc_mh_prior_chk.setToolTip(
            "분광 문헌값을 입력하세요 — 예: APOGEE/GALAH의 성단 [Fe/H].\n"
            "U/u 청색 밴드가 없는 자료에서 [M/H]를 정하는 표준 경로입니다\n"
            "(gri·BVR 색만으로는 ~0.2 dex 금속결핍으로 치우침, docs §10.9).\n"
            "σ는 그 측정의 실제 정밀도(보통 0.03~0.10 dex).")
        mh_row = QHBoxLayout()
        self.mcmc_mh_mean = QDoubleSpinBox(); self.mcmc_mh_mean.setRange(-2.5, 1.0)
        self.mcmc_mh_mean.setDecimals(3); self.mcmc_mh_mean.setSingleStep(0.05)
        self.mcmc_mh_mean.setToolTip("문헌 [Fe/H]/[M/H] 값 [dex]")
        self.mcmc_mh_sigma = QDoubleSpinBox(); self.mcmc_mh_sigma.setRange(0.001, 1.0)
        self.mcmc_mh_sigma.setDecimals(3); self.mcmc_mh_sigma.setSingleStep(0.01)
        self.mcmc_mh_sigma.setValue(0.10)
        self.mcmc_mh_sigma.setToolTip("측정 정밀도 [dex] — 너무 크면 축퇴를 못 깹니다")
        mh_row.addWidget(QLabel("μ")); mh_row.addWidget(self.mcmc_mh_mean)
        mh_row.addWidget(QLabel("σ")); mh_row.addWidget(self.mcmc_mh_sigma)
        prior_form.addRow(self.mcmc_mh_prior_chk, self._wrap(mh_row))

        self.mcmc_ebv_prior_chk = QCheckBox("E(B−V) prior (dust map)")
        self.mcmc_ebv_prior_chk.setToolTip(
            "먼지지도/문헌 적색화 — 예: SFD98·Green+2019(Bayestar), 성단\n"
            "카탈로그(Cantat-Gaudin+2020, Dias+2021). 젊은 산개성단의 단색\n"
            "적합은 이 prior 없이는 나이·E가 축퇴로 rail 합니다(NGC 6811 실측).")
        ebv_row = QHBoxLayout()
        self.mcmc_ebv_mean = QDoubleSpinBox(); self.mcmc_ebv_mean.setRange(0.0, 3.0)
        self.mcmc_ebv_mean.setDecimals(3); self.mcmc_ebv_mean.setSingleStep(0.01)
        self.mcmc_ebv_mean.setToolTip("문헌 E(B−V) [mag]")
        self.mcmc_ebv_sigma = QDoubleSpinBox(); self.mcmc_ebv_sigma.setRange(0.001, 1.0)
        self.mcmc_ebv_sigma.setDecimals(3); self.mcmc_ebv_sigma.setSingleStep(0.01)
        self.mcmc_ebv_sigma.setValue(0.02)
        self.mcmc_ebv_sigma.setToolTip("불확실성 [mag] — 차등적색화 성단이면 크게")
        ebv_row.addWidget(QLabel("μ")); ebv_row.addWidget(self.mcmc_ebv_mean)
        ebv_row.addWidget(QLabel("σ")); ebv_row.addWidget(self.mcmc_ebv_sigma)
        prior_form.addRow(self.mcmc_ebv_prior_chk, self._wrap(ebv_row))

        # visual affordance: μ/σ fields live only while their prior is on
        def _bind_prior_row(chk, widgets):
            def _sync(_state=None):
                for w in widgets:
                    w.setEnabled(chk.isChecked())
            chk.stateChanged.connect(_sync)
            _sync()
        _bind_prior_row(self.mcmc_mh_prior_chk,
                        (self.mcmc_mh_mean, self.mcmc_mh_sigma))
        _bind_prior_row(self.mcmc_ebv_prior_chk,
                        (self.mcmc_ebv_mean, self.mcmc_ebv_sigma))

        hint = QLabel("Gaia parallax → (m−M)₀ is set from the membership clump when available.")
        hint.setProperty("role", "caption")
        hint.setWordWrap(True)
        prior_form.addRow(hint)
        controls.addWidget(prior_group)
        outer.addLayout(controls)

        # -- run row --
        run_row = QHBoxLayout()
        self.btn_run_mcmc = QPushButton("Run MCMC Auto-Fit")
        style_button(self.btn_run_mcmc, "primary")
        self.btn_run_mcmc.clicked.connect(self._run_mcmc_autofit)
        run_row.addWidget(self.btn_run_mcmc)
        self.btn_save_mcmc = QPushButton("Save Figures")
        self.btn_save_mcmc.setEnabled(False)
        self.btn_save_mcmc.clicked.connect(self._save_mcmc_figures)
        run_row.addWidget(self.btn_save_mcmc)
        self.mcmc_progress = QProgressBar(); self.mcmc_progress.setRange(0, 100)
        run_row.addWidget(self.mcmc_progress, stretch=1)
        outer.addLayout(run_row)

        self.mcmc_status = QLabel("")
        _f = self.mcmc_status.font(); _f.setBold(True); self.mcmc_status.setFont(_f)
        outer.addWidget(self.mcmc_status)

        # -- results: summary + figures in a scroll area --
        self.mcmc_summary = QTextEdit(); self.mcmc_summary.setReadOnly(True)
        self.mcmc_summary.setMaximumHeight(120)
        outer.addWidget(self.mcmc_summary)

        fig_scroll = QScrollArea(); fig_scroll.setWidgetResizable(True)
        fig_holder = QWidget()
        self.mcmc_fig_layout = QVBoxLayout(fig_holder)
        fig_scroll.setWidget(fig_holder)
        outer.addWidget(fig_scroll, stretch=1)

        self._mcmc_worker = None
        self._mcmc_last_output = None
        saved = getattr(self, "_mcmc_saved_state", None)
        if saved:
            self._restore_mcmc_state(saved)
        # This page stacks fit options, the priors form and the run/progress
        # rows: 844 px against the CMD Viewer's 105. A QTabWidget takes the
        # *maximum* minimum over its pages, so unwrapped it demanded a
        # 1558 px-tall window — nearly twice a laptop screen — and Qt squeezed
        # the page to a 191 px sliver with its headings and Run button
        # overlapping. Wrapped, the page scrolls and the window fits.
        self.mcmc_tab_index = self.tabs.addTab(
            scroll_wrap(tab, horizontal=True), "Auto-fit (MCMC)")
        # Auto-detect available colours when the user opens this tab.
        self.tabs.currentChanged.connect(self._on_mcmc_tab_changed)

    def _on_mcmc_tab_changed(self, index):
        if index == getattr(self, "mcmc_tab_index", -1) and not self.mcmc_color_checks:
            try:
                self._mcmc_refresh_color_checks()
            except Exception as exc:  # pragma: no cover - best-effort UI refresh
                self.log(f"[MCMC] colour detect failed: {exc}")

    @staticmethod
    def _wrap(layout):
        """Wrap a layout in a QWidget (for QFormLayout rows)."""
        w = QWidget(); w.setLayout(layout); return w

    # Wavelength order for auto-detecting an adjacent colour chain. Limited to the
    # bands the bundled PARSEC grids carry (SDSS ugriz + Johnson UBVRI); more bands
    # (y, JHK, …) need a PARSEC grid with those columns + a PARSEC_ISO_COL mapping.
    _MCMC_BAND_ORDER = ["u", "g", "r", "i", "z", "U", "B", "V", "R", "I"]

    def _mcmc_available_bands(self, df) -> list:
        return [b for b in self._MCMC_BAND_ORDER if f"mag_std_{b}" in df.columns]

    def _mcmc_refresh_color_checks(self):
        """(Re)build one checkbox per adjacent colour available in the data.

        All ticked by default (use every colour). Untick to drop a colour; leave
        one ticked for a single-colour fit. Called when the tab is opened and via
        the 'detect' button.
        """
        # clear existing widgets
        while self.mcmc_colors_layout.count():
            it = self.mcmc_colors_layout.takeAt(0)
            w = it.widget()
            if w is not None:
                w.setParent(None)
        self.mcmc_color_checks = {}
        df, _iso_raw, _iso_file = self._load_cmd_and_iso_data(show_error=False)
        if df is None:
            self.mcmc_colors_layout.addWidget(QLabel("load CMD data first"))
            return
        bands = self._mcmc_available_bands(df)
        if len(bands) < 2:
            self.mcmc_colors_layout.addWidget(QLabel("need ≥2 bands"))
            return
        for i in range(len(bands) - 1):
            b1, b2 = bands[i], bands[i + 1]
            chk = QCheckBox(f"{b1}-{b2}")
            saved_colors = getattr(self, "_mcmc_saved_colors", None)
            chk.setChecked(True if not saved_colors
                           else f"{b1}-{b2}" in saved_colors)
            if b1 in ("u", "U"):
                chk.setToolTip("Blue/UV colour — this is the one that constrains [M/H].")
            self.mcmc_color_checks[(b1, b2)] = chk
            self.mcmc_colors_layout.addWidget(chk)
        self.mcmc_colors_layout.addStretch()

        # No blue band -> surface the spectroscopic-prior guidance up front
        # (the run-time log warning is too late to change the plan).
        has_blue = any(b in ("u", "U") for b in bands)
        self.mcmc_blue_warning.setVisible(not has_blue)
        self.mcmc_mh_prior_chk.setStyleSheet(
            "" if has_blue
            else f"QCheckBox {{ color: {Tokens.WARN}; font-weight: bold; }}")

    def _mcmc_parse_colors(self, bc, df) -> list:
        """Return the (b1, b2) colours ticked in the checkboxes.

        Falls back to the Band-Selection colour if the checkboxes are not built
        yet or none are ticked.
        """
        sel = [c for c, chk in self.mcmc_color_checks.items() if chk.isChecked()]
        if sel:
            return sel
        bands = self._mcmc_available_bands(df)
        if len(bands) >= 2:
            return [(bands[i], bands[i + 1]) for i in range(len(bands) - 1)]
        return [tuple(bc["band_color"])]

    def _run_mcmc_autofit(self):
        if self._mcmc_worker is not None and self._mcmc_worker.isRunning():
            return
        self.save_state()   # persist priors/settings the moment they are used
        bc = self._get_band_config()
        if not self._require_iso_config(bc):
            return
        df, iso_raw, iso_file = self._load_cmd_and_iso_data(show_error=True)
        if df is None:
            return
        try:
            from apex.analysis.cmd.isochrone_fit_service import IsochroneFitConfig
        except Exception as exc:
            QMessageBox.critical(self, "Error", f"Fit service unavailable: {exc}")
            return

        mag_band = bc["band_mag"]
        try:
            colors = self._mcmc_parse_colors(bc, df)
        except Exception as exc:
            QMessageBox.warning(self, "Colors", str(exc))
            return
        if not colors:
            QMessageBox.warning(self, "Colors", "No valid colours specified.")
            return
        has_blue = any(b1 in ("u", "U") for b1, _ in colors)
        if not has_blue:
            self.log("[MCMC] note: no blue/UV colour (u-g, U-B) -> [M/H] is "
                     "degenerate; rely on the [M/H] prior (see tab help).")
        b1, b2 = colors[0]
        r_color = EXTINCTION_R.get(b1, 3.303) - EXTINCTION_R.get(b2, 2.285)
        # Self-contained search range: age from this tab's Gyr controls; [M/H]/dm/E
        # left wide here because their priors (below) tighten them. No dependency
        # on the old 'Auto Fit' tab.
        a_min = max(self.mcmc_age_min.value(), 1e-3)
        a_max = max(self.mcmc_age_max.value(), a_min * 1.01)
        age_bounds = (9.0 + float(np.log10(a_min)), 9.0 + float(np.log10(a_max)))

        mh_prior = None
        if self.mcmc_mh_prior_chk.isChecked():
            mh_prior = (self.mcmc_mh_mean.value(), max(self.mcmc_mh_sigma.value(), 1e-3))
        ecolor_prior = None
        if self.mcmc_ebv_prior_chk.isChecked():
            # convert E(B-V) prior into E(color)=E(B-V)*(R_b1-R_b2)
            ecolor_prior = (self.mcmc_ebv_mean.value() * r_color,
                            max(self.mcmc_ebv_sigma.value() * abs(r_color), 1e-3))

        cfg = IsochroneFitConfig(
            colors=colors, mag_band=mag_band, iso_file=Path(iso_file),
            age_bounds=age_bounds,
            mh_bounds=(-1.0, 0.5), dm_bounds=(5.0, 18.0), ecolor_bounds=(0.0, 1.0),
            mh_prior=mh_prior, ecolor_prior=ecolor_prior,
            use_membership=self.mcmc_membership_chk.isChecked(),
            parallax_distance_prior=self.mcmc_parallax_chk.isChecked(),
            max_stars=self.mcmc_maxstars_spin.value(),
            n_walkers=self.mcmc_walkers_spin.value(),
            n_steps=self.mcmc_steps_spin.value(),
            n_burn=self.mcmc_burn_spin.value(),
        )

        # Record the settings actually used, to show alongside the results.
        self._mcmc_run_info = {
            "colors": "+".join(f"{a}-{b}" for a, b in colors),
            "age": (a_min, a_max),
            "mh_prior": (self.mcmc_mh_mean.value(), self.mcmc_mh_sigma.value())
            if self.mcmc_mh_prior_chk.isChecked() else None,
            "ebv_prior": (self.mcmc_ebv_mean.value(), self.mcmc_ebv_sigma.value())
            if self.mcmc_ebv_prior_chk.isChecked() else None,
            "parallax": self.mcmc_parallax_chk.isChecked(),
            "steps": self.mcmc_steps_spin.value(),
        }

        self.btn_run_mcmc.setEnabled(False)
        self.btn_save_mcmc.setEnabled(False)
        self.mcmc_progress.setValue(0)
        self.mcmc_status.setText("Starting…")
        self.log(f"[MCMC] auto-fit colors={colors} mag={mag_band} age={a_min}-{a_max}Gyr "
                 f"mh_prior={self._mcmc_run_info['mh_prior']} "
                 f"ebv_prior={self._mcmc_run_info['ebv_prior']} membership={cfg.use_membership}")
        self._mcmc_worker = McmcFitWorker(df, cfg)
        self._mcmc_worker.progress.connect(self._on_mcmc_progress)
        self._mcmc_worker.finished.connect(self._on_mcmc_autofit_done)
        self._mcmc_worker.start()

    def _on_mcmc_progress(self, frac, msg):
        self.mcmc_progress.setValue(int(max(0.0, min(1.0, frac)) * 100))
        self.mcmc_status.setText(str(msg))

    def _on_mcmc_autofit_done(self, out):
        self.btn_run_mcmc.setEnabled(True)
        if isinstance(out, Exception):
            self.mcmc_status.setText(f"Failed: {out}")
            self.log(f"[MCMC] failed: {out}")
            QMessageBox.critical(self, "MCMC fit failed", str(out))
            return
        self._mcmc_last_output = out
        self.mcmc_progress.setValue(100)
        s = out.summary

        def fmt(key, scale=1.0, name=None):
            p = s.get(key)
            if not p:
                return ""
            med, lo, hi = p[1] * scale, p[0] * scale, p[2] * scale
            return f"{name or key} = {med:.3f} (+{hi - med:.3f} / −{med - lo:.3f})\n"

        # Settings actually used (so the priors / age range are visible).
        ri = getattr(self, "_mcmc_run_info", {}) or {}
        txt = ""
        if ri:
            mhp = ri.get("mh_prior"); ebvp = ri.get("ebv_prior"); ag = ri.get("age", (None, None))
            dist = "Gaia π" if ri.get("parallax") else "free"
            if dist == "Gaia π" and out.member_meta.get("dm_from_parallax") is not None:
                dist += f" (m−M={out.member_meta['dm_from_parallax']})"
            txt += (f"SETTINGS — colours: {ri.get('colors')} | age range: "
                    f"{ag[0]:g}–{ag[1]:g} Gyr | steps: {ri.get('steps')}\n"
                    f"  priors → [M/H]: {('%.3f±%.3f' % mhp) if mhp else 'off'}; "
                    f"E(B−V): {('%.3f±%.3f' % ebvp) if ebvp else 'off'}; "
                    f"distance: {dist}\n"
                    "──────────────────────────────\n")
        txt += fmt("age_gyr", name="Age [Gyr]")
        txt += fmt("metallicity", name="[M/H]")
        txt += fmt("distance_mod", name="(m−M)₀")
        if s.get("e_bv"):
            txt += fmt("e_bv", name="E(B−V)")
        txt += f"acceptance = {s['acceptance_fraction']:.3f}   converged = {s['convergence_ok']}\n"
        mm = out.member_meta
        if mm.get("applied"):
            txt += f"members = {mm.get('n_members')} (PM {mm.get('pm_center')}, π {mm.get('parallax_center')})\n"
        txt += f"fit stars = {out.n_stars}\n"
        for w in out.warnings:
            txt += f"⚠ {w}\n"
        self.mcmc_summary.setPlainText(txt)
        self.mcmc_status.setText("Done" + (" (not converged — see warnings)" if not s["convergence_ok"] else ""))

        # Build the figures here, on the main thread (pyplot is not thread-safe).
        from apex.analysis.cmd import isochrone_plots as P
        import numpy as _np
        self._mcmc_cmd_fig = None
        self._mcmc_corner_fig = None
        try:
            if out.obs_color is not None and out.iso_color is not None:
                try:
                    _target = target_display_name(self.params)
                except Exception:
                    _target = "Cluster"
                _cl = out.color_label or "color"
                _ml = out.mag_label or "mag"
                cmd_title = f"{_target} — {_ml} vs {_cl} CMD"
                self._mcmc_cmd_fig = P.cmd_figure(
                    out.obs_color, out.obs_mag, out.iso_color, out.iso_mag,
                    color_label=_cl, mag_label=_ml,
                    title=cmd_title, annotations=out.annotations or None,
                )
            r = out.result
            if getattr(r, "flat_chain", None) is not None and getattr(r, "labels", None):
                self._mcmc_corner_fig = P.corner_figure(r.flat_chain, list(r.labels))
        except Exception as exc:
            self.log(f"[MCMC] figure render failed: {exc}")

        while self.mcmc_fig_layout.count():
            item = self.mcmc_fig_layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.setParent(None)
        for fig in (self._mcmc_cmd_fig, self._mcmc_corner_fig):
            if fig is None:
                continue
            cv = FigureCanvas(fig)
            # Give each canvas its full natural pixel size so the scroll area shows
            # the whole figure (instead of squishing/clipping it to the viewport).
            w_in, h_in = fig.get_size_inches()
            dpi = fig.get_dpi()
            cv.setMinimumSize(int(w_in * dpi * 0.95), int(h_in * dpi * 0.95))
            self.mcmc_fig_layout.addWidget(cv)
        self.btn_save_mcmc.setEnabled(self._mcmc_cmd_fig is not None)
        self.log(f"[MCMC] done: age={s.get('age_gyr',[None,'?'])[1]} [M/H]={s.get('metallicity',[None,'?'])[1]} acc={s['acceptance_fraction']:.2f}")

    def _save_mcmc_figures(self):
        cmd_fig = getattr(self, "_mcmc_cmd_fig", None)
        corner_fig = getattr(self, "_mcmc_corner_fig", None)
        if cmd_fig is None and corner_fig is None:
            return
        try:
            outdir = Path(step12_iso_dir(self.params.P.result_dir))
        except Exception:
            outdir = Path(self.params.P.result_dir)
        outdir.mkdir(parents=True, exist_ok=True)
        saved = []
        # savefig directly (do NOT close — the figures are live in the canvases).
        if cmd_fig is not None:
            p = outdir / "mcmc_cmd_isochrone.png"
            cmd_fig.savefig(p, dpi=160, bbox_inches="tight"); saved.append(p)
        if corner_fig is not None:
            p = outdir / "mcmc_corner.png"
            corner_fig.savefig(p, dpi=130, bbox_inches="tight"); saved.append(p)
        QMessageBox.information(self, "Saved", "Saved:\n" + "\n".join(str(p) for p in saved))
        self.log(f"[MCMC] figures saved to {outdir}")

    def _get_band_config(self):
        """Return current band selection. Iso/R values are ``None`` for unsupported filters.

        Use :meth:`_require_iso_config` to validate before fitting / iso-dependent work.
        """
        color_text = self.color_combo.currentText()
        parts = color_text.split("-")
        b1 = parts[0].strip() if parts else ""
        b2 = parts[1].strip() if len(parts) > 1 else ""
        bm = self.mag_combo.currentText().strip()
        iso_columns = self._iso_band_columns
        return {
            "band_color":  (b1, b2),
            "band_mag":    bm,
            "iso_col_1":   iso_columns.get(b1),
            "iso_col_2":   iso_columns.get(b2),
            "iso_col_mag": iso_columns.get(bm),
            "R_band1":     EXTINCTION_R.get(b1),
            "R_band2":     EXTINCTION_R.get(b2),
            "R_mag":       EXTINCTION_R.get(bm),
            "iso_system":  self._iso_photometric_system,
            "iso_bands":   tuple(iso_columns),
        }

    def _require_iso_config(self, bc: dict) -> bool:
        """Validate band config for iso fitting. Show user dialog and return False on miss."""
        b1, b2 = bc["band_color"]
        bm     = bc["band_mag"]
        missing_iso = [b for b, k in ((b1,"iso_col_1"), (b2,"iso_col_2"), (bm,"iso_col_mag"))
                       if bc[k] is None]
        missing_R   = [b for b, k in ((b1,"R_band1"),   (b2,"R_band2"),   (bm,"R_mag"))
                       if bc[k] is None]
        if not missing_iso and not missing_R:
            return True
        parts = []
        if missing_iso:
            available = ", ".join(bc.get("iso_bands", ())) or "none detected"
            parts.append(
                "Selected observed bands are absent from this isochrone file: "
                f"{', '.join(missing_iso)}\n"
                f"Isochrone system: {bc.get('iso_system', 'Unknown')}\n"
                f"Available isochrone bands: {available}"
            )
        if missing_R:
            parts.append(f"Missing extinction coefficients: {', '.join(missing_R)}")
        parts.append(
            "Choose an isochrone file generated for the same photometric "
            "system as the observations."
        )
        QMessageBox.warning(self, "Isochrone Filter Mismatch", "\n\n".join(parts))
        return False

    def _default_iso_dir(self) -> Path:
        preferred = Path.cwd() / "isochrone"
        if preferred.exists():
            return preferred
        return Path.cwd()

    def _set_iso_status(self, message: str):
        if self.iso_status_label is not None:
            self.iso_status_label.setText(message)

    def _list_iso_files(self, iso_dir: Path) -> list[Path]:
        files = sorted(
            p for p in iso_dir.iterdir()
            if p.is_file()
            and p.suffix.lower() in {".dat", ".txt"}
            and not p.name.startswith(".apex_")
        )
        chunk_files = [p for p in files if "_chunk_" in p.stem]
        if chunk_files:
            return chunk_files
        return files

    def _folder_cache_paths(self, iso_dir: Path) -> tuple[Path, Path]:
        cache_dir = iso_dir / ".apex_cache"
        return cache_dir / "combined_isochrones.dat", cache_dir / "combined_isochrones_meta.json"

    def _compute_iso_signature(self, paths: list[Path]) -> str:
        digest = hashlib.sha1()
        for path in paths:
            stat = path.stat()
            digest.update(path.name.encode("utf-8"))
            digest.update(str(stat.st_size).encode("ascii"))
            digest.update(str(stat.st_mtime_ns).encode("ascii"))
        return digest.hexdigest()

    def _ensure_folder_iso_cache(self, iso_dir: Path, files: list[Path], signature: str) -> Path:
        merged_path, meta_path = self._folder_cache_paths(iso_dir)
        merged_path.parent.mkdir(parents=True, exist_ok=True)

        if merged_path.exists() and meta_path.exists():
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                if meta.get("signature") == signature:
                    self.log(f"Using cached folder isochrone: {merged_path}")
                    self._set_iso_status(f"Folder mode: {len(files)} files cached")
                    return merged_path
            except Exception:
                pass

        self.log(f"Merging {len(files)} isochrone files from {iso_dir}")
        self._set_iso_status(f"Folder mode: merging {len(files)} files...")
        app = QApplication.instance()
        if app is not None:
            app.setOverrideCursor(Qt.WaitCursor)
            app.processEvents()

        tmp_path = merged_path.with_suffix(".tmp")
        try:
            with tmp_path.open("w", encoding="utf-8", newline="\n") as out_handle:
                out_handle.write(f"# APEX merged {len(files)} isochrone files from {iso_dir}\n")
                for idx, src_path in enumerate(files, start=1):
                    self.log(f"[iso merge] {idx}/{len(files)} {src_path.name}")
                    with src_path.open("r", encoding="utf-8", errors="replace") as in_handle:
                        for line in in_handle:
                            if idx > 1 and line.lstrip().startswith("#"):
                                continue
                            out_handle.write(line)
                    if app is not None and (idx == len(files) or idx % 5 == 0):
                        app.processEvents()

            tmp_path.replace(merged_path)
            meta_path.write_text(
                json.dumps(
                    {
                        "signature": signature,
                        "file_count": len(files),
                        "files": [p.name for p in files],
                    },
                    indent=2,
                ),
                encoding="utf-8",
            )
        finally:
            if tmp_path.exists():
                tmp_path.unlink(missing_ok=True)
            if app is not None:
                app.restoreOverrideCursor()

        self.log(f"Merged folder isochrone cache ready: {merged_path}")
        self._set_iso_status(f"Folder mode: {len(files)} files cached")
        return merged_path

    def _resolve_iso_source(self, show_error=True):
        iso_path = self._get_iso_path()
        if not iso_path:
            if show_error:
                QMessageBox.warning(self, "Missing File", "Select an isochrone file or folder first")
            return None, None, None

        source_path = Path(iso_path)
        if not source_path.exists():
            if show_error:
                QMessageBox.warning(self, "Missing File", f"Isochrone source not found: {source_path}")
            return None, None, None

        if source_path.is_dir():
            files = self._list_iso_files(source_path)
            if not files:
                if show_error:
                    QMessageBox.warning(self, "Missing File", f"No .dat/.txt isochrone files found in: {source_path}")
                self._set_iso_status("Folder mode: no .dat/.txt files found")
                return None, None, None
            signature = self._compute_iso_signature(files)
            merged_path = self._ensure_folder_iso_cache(source_path, files, signature)
            cache_token = f"dir::{source_path.resolve()}::{signature}"
            return merged_path, source_path, cache_token

        self._set_iso_status("Single file mode")
        stat = source_path.stat()
        cache_token = f"file::{source_path.resolve()}::{stat.st_size}::{stat.st_mtime_ns}"
        return source_path, source_path, cache_token

    def _on_band_changed(self):
        """Refresh CMD viewer when band selection changes."""
        if self._get_iso_path():
            self.refresh_cmd_viewer(show_error=False)

    def _load_cmd_and_iso_data(self, show_error=True):
        iso_file, source_path, iso_cache_token = self._resolve_iso_source(show_error=show_error)
        if iso_file is None or source_path is None or iso_cache_token is None:
            return None, None, None

        try:
            header_info = read_parsec_header(iso_file)
        except Exception as exc:
            if show_error:
                QMessageBox.critical(
                    self,
                    "Isochrone Header Error",
                    f"Failed to read isochrone column metadata: {exc}",
                )
            return None, None, None
        if not header_info.band_columns:
            if show_error:
                QMessageBox.warning(
                    self,
                    "Isochrone Header Error",
                    "No supported photometric columns were found in the "
                    "isochrone header. Use a file that includes its column "
                    "name comment.",
                )
            return None, None, None

        self._iso_band_columns = dict(header_info.band_columns)
        self._iso_photometric_system = header_info.photometric_system
        header_key = (
            str(iso_file),
            header_info.photometric_system,
            tuple(header_info.band_columns.items()),
        )
        if header_key != self._iso_header_key:
            self._iso_header_key = header_key
            bands_text = ", ".join(header_info.available_bands)
            self.log(
                "Isochrone header detected: "
                f"{header_info.photometric_system} | bands={bands_text}"
            )
            self._set_iso_status(
                f"{header_info.photometric_system} | bands: {bands_text}"
            )

        input_dir = step10_zp_dir(self.params.P.result_dir)
        if not input_dir.exists():
            input_dir = self.params.P.result_dir
        selected = self._get_band_config()
        selected_bands = (*selected["band_color"], selected["band_mag"])
        legacy_reason = _legacy_b_calibration_reason(input_dir, selected_bands)
        if legacy_reason:
            self.log(f"[CMD viewer] {legacy_reason}")
            if show_error:
                QMessageBox.warning(
                    self,
                    "Step 10 Rerun Required",
                    legacy_reason,
                )
            return None, None, None
        wide_path = input_dir / "median_by_ID_filter_wide_cmd.csv"
        if not wide_path.exists():
            wide_path = input_dir / "median_by_ID_filter_wide.csv"
        if not wide_path.exists():
            if show_error:
                QMessageBox.warning(self, "Missing Data", "CMD wide CSV not found")
            return None, None, None

        # ---- cache check ----
        wide_stat = wide_path.stat()
        cache_key = (iso_cache_token, str(wide_path), wide_stat.st_size, wide_stat.st_mtime_ns)
        if (
            self._cache_key == cache_key
            and self._cached_df is not None
            and self._cached_iso_raw is not None
        ):
            return self._cached_df, self._cached_iso_raw, self._cached_iso_file

        try:
            df = pd.read_csv(wide_path, dtype={"source_id": str, "gaia_source_id": str})
        except Exception as e:
            if show_error:
                QMessageBox.critical(self, "Error", f"Failed to load CMD data: {e}")
            return None, None, None
        self._refresh_photometry_source_label(df)
        df = self._enrich_cmd_gaia_columns(df)
        self._initialize_parallax_range(df)

        try:
            app = QApplication.instance()
            if app is not None:
                app.setOverrideCursor(Qt.WaitCursor)
                app.processEvents()
            iso_raw = self._load_iso_raw_fast(Path(iso_file))
            if iso_raw is None or iso_raw.size == 0:
                if show_error:
                    QMessageBox.warning(self, "Data Error", "Isochrone file is empty")
                return None, None, None
            max_column = max(self._iso_band_columns.values())
            if max_column >= iso_raw.shape[1]:
                raise ValueError(
                    "Isochrone header/data column mismatch: "
                    f"header requires column {max_column}, data has "
                    f"{iso_raw.shape[1]} columns"
                )
        except Exception as e:
            if show_error:
                QMessageBox.critical(self, "Error", f"Failed to parse isochrone file: {e}")
            return None, None, None
        finally:
            app = QApplication.instance()
            if app is not None:
                app.restoreOverrideCursor()

        # ---- store in cache ----
        self._cache_key = cache_key
        self._cached_df = df
        self._cached_iso_raw = iso_raw
        self._cached_iso_file = iso_file

        return df, iso_raw, iso_file

    def _refresh_photometry_source_label(
        self,
        df: Optional[pd.DataFrame] = None,
    ) -> None:
        if df is None:
            input_dir = step10_zp_dir(self.params.P.result_dir)
            candidates = (
                input_dir / "median_by_ID_filter_wide_cmd.csv",
                input_dir / "median_by_ID_filter_wide.csv",
            )
            wide_path = next((path for path in candidates if path.exists()), None)
            if wide_path is not None:
                try:
                    df = pd.read_csv(wide_path, nrows=500)
                except Exception:
                    df = None
        provenance = summarize_photometry_table(df) if df is not None else None
        self.photometry_source_label.setText(
            format_photometry_provenance(provenance)
        )

    def _load_iso_raw_fast(self, iso_file: Path) -> Optional[np.ndarray]:
        """Load an isochrone table quickly, with an on-disk .npy cache.

        Large text grids can contain millions of rows; running
        ``np.genfromtxt`` on every Step 12 open burned ~30 s and froze the
        window.  This helper:

          * Uses ``np.loadtxt`` (≈3× faster than genfromtxt on the same data).
          * Saves a binary ``.apex.npy`` sidecar next to the source file —
            subsequent loads are <100 ms regardless of file size, and the
            cache is auto-invalidated when the source's mtime/size changes.

        Falls back to the slow text parse if the cache is unreadable or
        the sidecar cannot be written (read-only directory).
        """
        cache_path = iso_file.with_suffix(iso_file.suffix + ".apex.npy")
        try:
            src_stat = iso_file.stat()
        except OSError:
            return None

        if cache_path.exists():
            try:
                cache_stat = cache_path.stat()
                if (
                    cache_stat.st_mtime_ns >= src_stat.st_mtime_ns
                    and cache_stat.st_size > 0
                ):
                    arr = np.load(cache_path, allow_pickle=False)
                    if arr.ndim == 2 and arr.size > 0:
                        return arr
            except Exception:
                pass  # corrupt cache → fall through to rebuild

        arr = np.loadtxt(iso_file, comments="#")
        if arr.ndim == 1:
            arr = arr.reshape(1, -1)
        arr = arr[~np.isnan(arr).any(axis=1)]

        try:
            cache_path.parent.mkdir(parents=True, exist_ok=True)
            np.save(cache_path, arr, allow_pickle=False)
        except Exception:
            pass  # cache write failure is non-fatal

        return arr

    def _invalidate_cache(self):
        """Force reload on next access (call after file change)."""
        self._cache_key = None
        self._cached_df = None
        self._cached_iso_raw = None
        self._cached_iso_file = None
        self._iso_band_columns = {}
        self._iso_photometric_system = "Unknown"
        self._iso_header_key = None

    def _show_viewer_placeholder(self, message: str):
        self._clear_viewer_widget()
        placeholder = QLabel(message)
        placeholder.setAlignment(Qt.AlignCenter)
        placeholder.setProperty("role", "caption")
        placeholder.setStyleSheet(
            f"QLabel {{ border: 1px dashed {Tokens.BORDER}; padding: {Tokens.MARGIN}px; }}"
        )
        self.viewer_layout.addWidget(placeholder, stretch=1)
        self.viewer_placeholder = placeholder

    def _clear_viewer_widget(self):
        while self.viewer_layout.count():
            item = self.viewer_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
                widget.deleteLater()
        self.viewer = None

    def _repopulate_band_combos(self, df: "pd.DataFrame") -> None:
        """Repopulate color/mag combos from calibrated columns in df."""
        pairs, bands = _bands_from_df(df)
        pending_color = getattr(self, "_pending_band_color_text", None)
        pending_mag = getattr(self, "_pending_band_mag_text", None)
        prev_color = pending_color or self.color_combo.currentText()
        prev_mag   = pending_mag or self.mag_combo.currentText()

        self.color_combo.blockSignals(True)
        self.mag_combo.blockSignals(True)
        self.color_combo.clear()
        self.mag_combo.clear()
        self.color_combo.addItems([f"{a}-{b}" for a, b in pairs])
        self.mag_combo.addItems(bands)

        idx_c = self.color_combo.findText(prev_color)
        idx_m = self.mag_combo.findText(prev_mag)
        self.color_combo.setCurrentIndex(max(idx_c, 0))
        self.mag_combo.setCurrentIndex(max(idx_m, 0))
        self.color_combo.blockSignals(False)
        self.mag_combo.blockSignals(False)
        self._pending_band_color_text = None
        self._pending_band_mag_text = None

    def refresh_cmd_viewer(self, show_error=True) -> bool:
        df, iso_raw, iso_file = self._load_cmd_and_iso_data(show_error=show_error)
        if df is None or iso_raw is None:
            self._show_viewer_placeholder(
                "CMD viewer not ready.\nSelect an isochrone file and ensure Step 10 output exists."
            )
            return False

        self._repopulate_band_combos(df)
        df = self._apply_source_filters(df)
        self._clear_viewer_widget()
        bc = self._get_band_config()
        missing_iso = [
            b for b, k in (
                (bc["band_color"][0], "iso_col_1"),
                (bc["band_color"][1], "iso_col_2"),
                (bc["band_mag"], "iso_col_mag"),
            )
            if bc[k] is None
        ]
        if missing_iso:
            available = ", ".join(self._iso_band_columns) or "none"
            message = (
                "Isochrone/filter mismatch.\n"
                f"Observed selection: {bc['band_color'][0]}-"
                f"{bc['band_color'][1]}, {bc['band_mag']}\n"
                f"Isochrone: {self._iso_photometric_system} "
                f"({available})\n\n"
                "Select an isochrone file generated for the observation filters."
            )
            self.log(f"[CMD viewer] {message.replace(chr(10), ' ')}")
            self._show_viewer_placeholder(message)
            if show_error:
                self._require_iso_config(bc)
            return False
        viewer = IsochroneViewerWindow(
            df, iso_raw, self.params, self.viewer_container, embedded=True,
            band_mag=bc["band_mag"], band_color=bc["band_color"],
            iso_col_1=bc["iso_col_1"], iso_col_2=bc["iso_col_2"],
            iso_col_mag=bc["iso_col_mag"],
        )
        self.viewer_layout.addWidget(viewer, stretch=1)
        self.viewer = viewer
        n_total = len(self._cached_df) if self._cached_df is not None else len(df)
        self.log(
            f"CMD viewer updated: {iso_file.name}  "
            f"bands={bc['band_color'][0]}-{bc['band_color'][1]}, mag={bc['band_mag']}  "
            f"N={len(df)}/{n_total}"
        )
        return True

    def show_log_window(self):
        self.log_window.show()
        self.log_window.raise_()
        self.log_window.activateWindow()

    def browse_iso_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Isochrone File",
            str(self._default_iso_dir()),
            "Data Files (*.dat *.txt);;All Files (*.*)",
        )
        if path:
            self.iso_path_edit.setText(path)
            self.params.P.iso_file_path = path
            self._set_iso_status("Single file mode")
            self._invalidate_cache()
            self.save_state()
            self.persist_params()
            self.refresh_cmd_viewer(show_error=True)
            self.update_navigation_buttons()

    def browse_iso_folder(self):
        path = QFileDialog.getExistingDirectory(
            self,
            "Select Isochrone Folder",
            str(self._default_iso_dir()),
        )
        if path:
            self.iso_path_edit.setText(path)
            self.params.P.iso_file_path = path
            self._set_iso_status("Folder mode")
            self._invalidate_cache()
            self.save_state()
            self.persist_params()
            self.refresh_cmd_viewer(show_error=True)
            self.update_navigation_buttons()

    def open_viewer(self):
        if self.refresh_cmd_viewer(show_error=True):
            self.tabs.setCurrentIndex(self.cmd_viewer_tab_index)

    def validate_step(self) -> bool:
        iso_path = ""
        if getattr(self, "iso_path_edit", None) is not None:
            iso_path = self.iso_path_edit.text().strip()
        if not iso_path:
            iso_path = str(getattr(self.params.P, "iso_file_path", ""))
        if not iso_path:
            return False
        if not Path(iso_path).exists():
            return False
        input_dir = step10_zp_dir(self.params.P.result_dir)
        if not input_dir.exists():
            input_dir = self.params.P.result_dir
        return (input_dir / "median_by_ID_filter_wide_cmd.csv").exists() or (input_dir / "median_by_ID_filter_wide.csv").exists()

    def save_state(self):
        state_data = {
            "iso_file_path": self.iso_path_edit.text().strip() or str(getattr(self.params.P, "iso_file_path", "")),
            "iso_age_init": getattr(self.params.P, "iso_age_init", 9.7),
            "iso_mh_init": getattr(self.params.P, "iso_mh_init", -0.1),
            "iso_eg_r_init": getattr(self.params.P, "iso_eg_r_init", 0.0033),
            "iso_dm_init": getattr(self.params.P, "iso_dm_init", 9.46),
            "band_color_idx": self.color_combo.currentIndex(),
            "band_mag_idx": self.mag_combo.currentIndex(),
            "band_color_text": self.color_combo.currentText(),
            "band_mag_text": self.mag_combo.currentText(),
            # Auto-fit (MCMC) settings — priors survive a restart (a re-typed
            # prior is a re-typed literature lookup; persist it).
            # "workspace" guards against the shared primary state file
            # (apex/.state/<mode> is one file for ALL workspaces): restoring
            # cluster A's [M/H]/E priors into cluster B would silently
            # contaminate a science fit, so restore only for the same
            # result_dir.
            "mcmc": {
                "workspace": str(getattr(self.params.P, "result_dir", "")),
                "membership": self.mcmc_membership_chk.isChecked(),
                "parallax_prior": self.mcmc_parallax_chk.isChecked(),
                "mh_prior_on": self.mcmc_mh_prior_chk.isChecked(),
                "mh_mean": self.mcmc_mh_mean.value(),
                "mh_sigma": self.mcmc_mh_sigma.value(),
                "ebv_prior_on": self.mcmc_ebv_prior_chk.isChecked(),
                "ebv_mean": self.mcmc_ebv_mean.value(),
                "ebv_sigma": self.mcmc_ebv_sigma.value(),
                "age_min": self.mcmc_age_min.value(),
                "age_max": self.mcmc_age_max.value(),
                "walkers": self.mcmc_walkers_spin.value(),
                "steps": self.mcmc_steps_spin.value(),
                "burn": self.mcmc_burn_spin.value(),
                "max_stars": self.mcmc_maxstars_spin.value(),
                "colors_checked": [f"{a}-{b}" for (a, b), chk
                                   in self.mcmc_color_checks.items()
                                   if chk.isChecked()],
            },
        }
        self.project_state.store_step_data("isochrone_model", state_data)

    def _restore_mcmc_state(self, m: dict):
        """Apply saved Auto-fit (MCMC) settings back onto the widgets."""
        try:
            self.mcmc_membership_chk.setChecked(bool(m.get("membership", True)))
            self.mcmc_parallax_chk.setChecked(bool(m.get("parallax_prior", True)))
            self.mcmc_mh_prior_chk.setChecked(bool(m.get("mh_prior_on", False)))
            self.mcmc_mh_mean.setValue(float(m.get("mh_mean", 0.0)))
            self.mcmc_mh_sigma.setValue(float(m.get("mh_sigma", 0.10)))
            self.mcmc_ebv_prior_chk.setChecked(bool(m.get("ebv_prior_on", False)))
            self.mcmc_ebv_mean.setValue(float(m.get("ebv_mean", 0.0)))
            self.mcmc_ebv_sigma.setValue(float(m.get("ebv_sigma", 0.02)))
            self.mcmc_age_min.setValue(float(m.get("age_min", 0.20)))
            self.mcmc_age_max.setValue(float(m.get("age_max", 6.0)))
            self.mcmc_walkers_spin.setValue(int(m.get("walkers", 32)))
            self.mcmc_steps_spin.setValue(int(m.get("steps", 2000)))
            self.mcmc_burn_spin.setValue(int(m.get("burn", 600)))
            self.mcmc_maxstars_spin.setValue(int(m.get("max_stars", 500)))
            # colour checkboxes are rebuilt on tab open; remember the wish list
            self._mcmc_saved_colors = list(m.get("colors_checked", []))
        except Exception as exc:
            self.log(f"[MCMC] saved settings restore skipped: {exc}")

    def restore_state(self):
        state_data = self.project_state.get_step_data("isochrone_model")
        if state_data:
            for key, val in state_data.items():
                if hasattr(self.params.P, key):
                    setattr(self.params.P, key, val)
            if state_data.get("iso_file_path"):
                if self.iso_path_edit is not None:
                    self.iso_path_edit.setText(state_data["iso_file_path"])
            self._pending_band_color_text = state_data.get("band_color_text")
            self._pending_band_mag_text = state_data.get("band_mag_text")
            if "band_color_idx" in state_data:
                idx = int(state_data["band_color_idx"])
                if 0 <= idx < self.color_combo.count():
                    self.color_combo.setCurrentIndex(idx)
            if "band_mag_idx" in state_data:
                idx = int(state_data["band_mag_idx"])
                if 0 <= idx < self.mag_combo.count():
                    self.mag_combo.setCurrentIndex(idx)
            # restore_state runs AFTER setup_step_ui() (see __init__), so the
            # MCMC widgets already exist here — apply directly. The stash +
            # build-site apply stays as a guard for any future order change.
            # Workspace guard: the primary state file is shared by every
            # workspace of this mode, so drop saved priors that belong to a
            # different result_dir (cross-cluster prior contamination).
            cur_ws = str(getattr(self.params.P, "result_dir", "") or "")

            def _ws_of(m: dict) -> str:
                return str((m or {}).get("workspace", "") or "")

            saved_mcmc = state_data.get("mcmc") or {}
            if saved_mcmc and _ws_of(saved_mcmc) and _ws_of(saved_mcmc) != cur_ws:
                self.log(f"[MCMC] primary saved settings belong to another "
                         f"workspace ({_ws_of(saved_mcmc)}) — ignored")
                saved_mcmc = {}
            if not saved_mcmc and cur_ws:
                # The shared primary holds only the LAST-closed workspace's
                # settings; this workspace's own copy lives in its mirror
                # (result_dir/project_state.json, written on every save).
                # The mirror's location proves ownership, so a legacy entry
                # without a workspace tag is accepted from here.
                try:
                    mirror = Path(cur_ws) / "project_state.json"
                    if mirror.exists():
                        md = json.loads(mirror.read_text(encoding="utf-8"))
                        mm = ((md.get("step_data") or {})
                              .get("isochrone_model") or {}).get("mcmc") or {}
                        if mm and _ws_of(mm) in ("", cur_ws):
                            saved_mcmc = mm
                            self.log("[MCMC] settings restored from the "
                                     "workspace mirror")
                except Exception as exc:
                    self.log(f"[MCMC] mirror restore skipped: {exc}")
            self._mcmc_saved_state = saved_mcmc
            if self._mcmc_saved_state and hasattr(self, "mcmc_mh_prior_chk"):
                self._restore_mcmc_state(self._mcmc_saved_state)
        if self.iso_path_edit is not None and not self.iso_path_edit.text().strip():
            iso_path = str(getattr(self.params.P, "iso_file_path", "") or "")
            if iso_path:
                self.iso_path_edit.setText(iso_path)
        current_iso = self._get_iso_path()
        if current_iso:
            self._set_iso_status("Folder mode" if Path(current_iso).is_dir() else "Single file mode")
        if self._get_iso_path():
            # The Color-Color plot is light enough to render up front.
            # The full CMD viewer (isochrone load + matplotlib scatter) can
            # take a few seconds on cold cache, so we mark it pending and
            # let _on_tab_changed render it the first time the user clicks
            # the CMD Viewer tab.  Without this the Step 12 window
            # appears frozen during open.
            if self.tabs.currentIndex() == getattr(self, "cmd_viewer_tab_index", -1):
                self._cmd_viewer_pending = False
            else:
                self._cmd_viewer_pending = True
            self._maybe_start_iso_cache_build()
        self.update_navigation_buttons()

    def _maybe_start_iso_cache_build(self) -> None:
        """Render plots immediately if the binary ``.apex.npy`` cache is
        ready, otherwise spawn an ``IsoCacheWorker`` to build it in the
        background while the window stays responsive.

        A large text grid can take tens of seconds to
        parse on cold cache; freezing the UI during that period is what
        the user reported as "step 12가 켜지는데 오래 걸림".
        """
        iso_path_str = self._get_iso_path() or ""
        if not iso_path_str:
            return
        iso_path = Path(iso_path_str)
        if not iso_path.exists():
            return
        if iso_path.is_dir():
            QTimer.singleShot(0, self._render_pending_plots)
            return

        cache_path = iso_path.with_suffix(iso_path.suffix + ".apex.npy")
        cache_fresh = False
        try:
            if cache_path.exists():
                if cache_path.stat().st_mtime_ns >= iso_path.stat().st_mtime_ns:
                    cache_fresh = True
        except OSError:
            cache_fresh = False

        if cache_fresh:
            QTimer.singleShot(0, self._render_pending_plots)
            return

        worker = self._iso_cache_worker
        if worker is not None and worker.isRunning():
            if worker.iso_file == iso_path:
                return
            worker.cancel()
            self._pending_iso_cache_path = iso_path
            self._set_iso_status("Waiting for previous cache build to stop.")
            return

        self._set_iso_status(
            f"Building isochrone cache in background — first open takes "
            f"~30 s ({iso_path.name})."
        )
        self.log(f"[iso] background cache build started: {iso_path.name}")
        self._iso_cache_worker = IsoCacheWorker(iso_path)
        self._iso_cache_worker.progress_msg.connect(self.log)
        self._iso_cache_worker.finished_ok.connect(self._on_iso_cache_ready)
        self._iso_cache_worker.failed.connect(self._on_iso_cache_failed)
        self._iso_cache_worker.finished.connect(self._on_iso_cache_worker_finished)
        self._iso_cache_worker.start()

    def _render_pending_plots(self) -> None:
        if not getattr(self, "_cmd_viewer_pending", False):
            try:
                self.refresh_cmd_viewer(show_error=False)
            except Exception as exc:
                self.log(f"[iso] viewer error: {exc}")

    def _on_iso_cache_ready(self, iso_path_str: str) -> None:
        if Path(iso_path_str) != Path(self._get_iso_path() or ""):
            return
        self.log(f"[iso] cache ready: {Path(iso_path_str).name}")
        self._set_iso_status("Single file mode (cache ready)")
        QTimer.singleShot(0, self._render_pending_plots)

    def _on_iso_cache_failed(self, iso_path_str: str, err: str) -> None:
        if Path(iso_path_str) != Path(self._get_iso_path() or ""):
            return
        self.log(f"[iso] background cache build failed: {err}")
        # Fall back to synchronous load (slow but at least works).
        QTimer.singleShot(0, self._render_pending_plots)

    def _on_iso_cache_worker_finished(self) -> None:
        worker = self._iso_cache_worker
        if worker is not None:
            worker.deleteLater()
        self._iso_cache_worker = None
        pending = self._pending_iso_cache_path
        self._pending_iso_cache_path = None
        if pending is not None and Path(self._get_iso_path() or "") == pending:
            QTimer.singleShot(0, self._maybe_start_iso_cache_build)

    def closeEvent(self, event):
        worker = self._iso_cache_worker
        if worker is not None and worker.isRunning():
            worker.cancel()
            if not worker.wait(10000):
                QMessageBox.warning(
                    self,
                    "Background Task Running",
                    "Isochrone cache generation is still stopping. Please wait and close again.",
                )
                event.ignore()
                return
        super().closeEvent(event)

    def _on_tab_changed(self, index: int) -> None:
        if (
            index == getattr(self, "cmd_viewer_tab_index", -1)
            and getattr(self, "_cmd_viewer_pending", False)
        ):
            self._cmd_viewer_pending = False
            # singleShot so the tab UI paints before the heavy refresh.
            QTimer.singleShot(0, lambda: self.refresh_cmd_viewer(show_error=False))
