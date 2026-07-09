"""
Step 3: Sky Preview & QC Window
Interactive image viewer with click-to-measure stats
Sky preview workflow with imexamine-like keyboard shortcuts.
"""

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QGroupBox, QMessageBox, QTextEdit, QComboBox, QDialog,
    QFormLayout, QLineEdit, QDialogButtonBox, QSplitter, QApplication,
    QSlider, QSpinBox, QDoubleSpinBox
)
from PyQt5.QtCore import Qt, pyqtSignal
from pathlib import Path
import json
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.visualization import ZScaleInterval, ImageNormalize
from apex.utils.constants import MAG_ERR_COEFF
from photutils.aperture import CircularAperture, CircularAnnulus, ApertureStats
from photutils.detection import DAOStarFinder
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
from PyQt5.QtGui import QColor
from apex.gui.layout_rules import FittedDialog, tame_canvas
from apex.gui.widgets.fits_viewer import FITSViewerWidget, OverlayMarker

from .step_window_base import StepWindowBase
from .ui_helpers import (
    add_parameter_reset_button,
    build_scroll_param_dialog,
    create_parameter_button,
)
from apex.utils.step_paths import (
    step1_dir,
    step2_cropped_dir,
    step7_forced_phot_dir,
    crop_is_active,
)
from apex.utils.common_helpers import normalize_filter_key as _canonical_filter_key
from apex.utils.noise_params import resolve_effective_noise_params


_MIN_R_AP_PX = 4.0
_MIN_R_IN_PX = 12.0
_MIN_R_OUT_PX = 20.0


class SkyPreviewWindow(StepWindowBase):
    """
    Step 3: Sky Preview & QC
    Click on stars to measure FWHM, SNR, magnitude
    """

    def __init__(self, params, file_manager, project_state, main_window):
        """Initialize sky preview window"""
        self.file_manager = file_manager

        # Image data
        self.current_filename = None
        self.image_data = None
        self.header = None
        self.current_file_index = 0

        # Image viewer (FITSViewerWidget — GPU-accelerated)
        self._viewer: FITSViewerWidget | None = None

        # Cursor position in image coords (updated by mouse_moved signal)
        self.cursor_x = None
        self.cursor_y = None

        # Photometry parameters from parameters.toml (hud5x section)
        hud5 = getattr(params.P, "_hud5", {})
        self.aperture_scale = float(hud5.get("5x.aperture_scale", "") or 1.0)
        self.annulus_in_scale = float(hud5.get("5x.annulus_in_scale", "") or 4.0)
        self.annulus_out_scale = float(hud5.get("5x.annulus_out_scale", "") or 2.0)  # width
        # Hidden safety floors only. User-facing control stays scale-based.
        self.min_r_ap_px = _MIN_R_AP_PX
        self.min_r_in_px = _MIN_R_IN_PX
        self.min_r_out_px = _MIN_R_OUT_PX
        self.sigma_clip_value = float(hud5.get("5x.sigma_clip", "") or 3.0)
        # FWHM seed in arcsec (converted to px using pixel_scale)
        self.fwhm_seed_arcsec = float(getattr(params.P, "fwhm_guess_arcsec", None) or 2.5)

        # File list
        self.file_list = []
        self.use_cropped = False
        self.cropped_dir = None

        # Persistent analysis windows
        self.histogram_dialog = None
        self.radial_profile_dialog = None
        self.last_radial_fwhm_px = None
        self.last_radial_center = None
        self.hist_fig = None
        self.hist_canvas = None
        self.hist_ax = None
        self.hist_stats_label = None
        self.prof_fig = None
        self.prof_canvas = None
        self.prof_ax = None
        self.last_measurement = None
        self._file_filter_map = {}
        self._file_frame_key_map = {}
        self._frame_key_map = {}
        self._frame_keys_by_filter = {}
        self._filter_order = []
        # FITS data cache for fast frame switching
        self._fits_cache: dict = {}
        self._fits_cache_order: list = []
        self._FITS_CACHE_SIZE = max(3, int(getattr(params.P, "step3_fits_cache_size", 6)))

        # Stretch plot window (2D Plot)
        self.stretch_plot_dialog = None
        self.stretch_plot_canvas = None
        self.stretch_plot_ax = None
        self.stretch_plot_fig = None
        self.stretch_plot_info_label = None
        self._stretch_vmin = None
        self._stretch_vmax = None
        self._stretch_data_range = None
        self._stretch_dragging = False
        self._stretch_drag_target = None  # "min" or "max"
        self._stretch_marker_min_line = None
        self._stretch_marker_max_line = None

        # Initialize base class
        super().__init__(
            step_index=2,
            step_name="Sky Preview & QC",
            params=params,
            project_state=project_state,
            main_window=main_window
        )

        # Setup step-specific UI
        self.setup_step_ui()

        # Restore state AFTER UI is created
        self.restore_state()

    def _sigma_clip_sigma(self) -> float:
        try:
            value = float(self.sigma_clip_value)
        except Exception:
            return 3.0
        return value if np.isfinite(value) and value > 0 else 3.0

    def setup_step_ui(self):
        """Setup step-specific UI components"""

        # === Info Label ===
        info_label = QLabel(
            "Keyboard: [m] measure at cursor | [h] histogram | [g] radial profile | [.] next filter | [[/]] prev/next frame\n"
            "Mouse: Wheel to zoom | Right-click drag to pan"
        )
        info_label.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; border-radius: 5px; }")
        self.content_layout.addWidget(info_label)

        # === File Selector ===
        file_layout = QHBoxLayout()
        file_layout.addWidget(QLabel("Image:"))

        self.file_combo = QComboBox()
        self.file_combo.currentIndexChanged.connect(self.on_file_changed)
        file_layout.addWidget(self.file_combo)

        btn_params = create_parameter_button("Photometry Parameters")
        btn_params.clicked.connect(self.open_parameters_dialog)
        file_layout.addWidget(btn_params)

        file_layout.addStretch()
        self.content_layout.addLayout(file_layout)

        # === Image Viewer ===
        viewer_group = QGroupBox("Image Viewer")
        viewer_layout = QVBoxLayout(viewer_group)

        top_bar = QHBoxLayout()
        btn_2d_plot = QPushButton("2D Plot")
        btn_2d_plot.setStyleSheet("QPushButton { background-color: #FF9800; color: white; font-weight: bold; }")
        btn_2d_plot.clicked.connect(self.open_stretch_plot)
        top_bar.addWidget(btn_2d_plot)
        top_bar.addStretch()
        viewer_layout.addLayout(top_bar)

        self._viewer = FITSViewerWidget(self)
        self._viewer.mouse_moved.connect(self._on_viewer_hover)
        self._viewer.mouse_pressed.connect(self._on_viewer_click)
        viewer_layout.addWidget(self._viewer)

        # === Create splitter for viewer and stats side-by-side ===
        main_splitter = QSplitter(Qt.Horizontal)
        main_splitter.addWidget(viewer_group)

        # === Stats Display (right side) ===
        stats_group = QGroupBox("Measurement Statistics")
        stats_layout = QVBoxLayout(stats_group)

        self.stats_text = QTextEdit()
        self.stats_text.setReadOnly(True)
        self.stats_text.setStyleSheet("QTextEdit { font-family: monospace; font-size: 10pt; }")
        self.stats_text.setPlainText("Press 'm' at cursor to measure\n'h' for histogram | 'g' for radial profile")

        self.stats_text.setMinimumWidth(360)
        stats_layout.addWidget(self.stats_text)
        main_splitter.addWidget(stats_group)
        main_splitter.setStretchFactor(0, 3)
        main_splitter.setStretchFactor(1, 2)
        main_splitter.setChildrenCollapsible(False)

        self.content_layout.addWidget(main_splitter, 1)

        # Populate file list
        self.populate_file_list()

    def populate_file_list(self):
        """Populate file combo box with cropped files"""
        # Check for cropped files
        crop_active = crop_is_active(self.params.P.result_dir)
        cropped_dir = step2_cropped_dir(self.params.P.result_dir)

        if crop_active and cropped_dir.exists() and list(cropped_dir.glob("*.fit*")):
            # Use cropped files
            files = sorted([f.name for f in cropped_dir.glob("*.fit*")])
            self.use_cropped = True
            self.cropped_dir = cropped_dir
        else:
            # Use original files
            if not self.file_manager.filenames:
                try:
                    self.file_manager.scan_files()
                except Exception:
                    pass
            files = self.file_manager.filenames
            self.use_cropped = False
            self.cropped_dir = None

        self.file_list = list(files)
        self.file_combo.clear()
        self.file_combo.addItems(files)
        self._fits_cache.clear()
        self._fits_cache_order.clear()
        self._file_filter_map = {}
        self._file_frame_key_map = {}
        self._frame_key_map = {}
        self._frame_keys_by_filter = {}
        self._filter_order = []
        self._load_filter_map_from_index()

    def load_selected_image(self, keep_view=False):
        """Load selected image"""
        filename = self.file_combo.currentText()
        if not filename:
            return

        try:
            new_data, new_header = self._load_fits_cached(filename)
            first_load = self.image_data is None or new_data.shape != self.image_data.shape
            self.image_data = new_data
            self.header = new_header

            self.current_filename = filename
            self.current_file_index = self.file_combo.currentIndex()
            self._file_filter_map.setdefault(filename, self._extract_filter_from_header(self.header))
            self._file_frame_key_map.setdefault(
                filename,
                self._extract_frame_key(filename, self._file_filter_map.get(filename, ""))
            )

            self.reset_stretch_plot_values()
            self._update_viewer_image(keep_view=keep_view and not first_load)
            self._viewer.setFocus()

            if not (keep_view and self.last_measurement is not None):
                self.stats_text.setPlainText(
                    f"Loaded: {filename}\n"
                    f"Size: {self.image_data.shape}\n"
                    f"Press 'm' at cursor to measure | 'h' for histogram | 'g' for radial profile"
                )

            if keep_view:
                self._restore_measurement_overlay()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load image:\n{str(e)}")

    def _get_fits_path(self, filename: str) -> Path:
        if self.use_cropped:
            cropped_dir = self.cropped_dir if self.cropped_dir is not None else step2_cropped_dir(self.params.P.result_dir)
            return cropped_dir / filename
        try:
            return Path(self.params.get_file_path(filename))
        except Exception:
            return self.params.P.data_dir / filename

    def _load_fits_cached(self, filename: str):
        if filename in self._fits_cache:
            try:
                self._fits_cache_order.remove(filename)
            except ValueError:
                pass
            self._fits_cache_order.append(filename)
            return self._fits_cache[filename]

        file_path = self._get_fits_path(filename)
        with fits.open(file_path, memmap=False) as hdul:
            data_raw = hdul[0].data
            if data_raw is None:
                raise ValueError(f"FITS image data is empty: {file_path.name}")
            data = np.asarray(data_raw, dtype=np.float32)
            header = hdul[0].header.copy()

        while len(self._fits_cache_order) >= self._FITS_CACHE_SIZE:
            oldest = self._fits_cache_order.pop(0)
            self._fits_cache.pop(oldest, None)
        self._fits_cache[filename] = (data, header)
        self._fits_cache_order.append(filename)
        return data, header

    def _update_viewer_image(self, keep_view: bool = False):
        """Upload raw data and refresh auto-STF. keep_view preserves zoom/pan only."""
        if self.image_data is None or self._viewer is None:
            return
        self._viewer.set_data(self.image_data)
        self._viewer.auto_stf()
        if not keep_view:
            self._viewer.fit_in_view()

    def _on_viewer_hover(self, x: float, y: float, val: float):
        self.cursor_x = x
        self.cursor_y = y

    def _on_viewer_click(self, x: float, y: float, btn: int):
        from PyQt5.QtCore import Qt as _Qt
        if btn == int(_Qt.MiddleButton):
            self.measure_star(int(x), int(y))

    def on_file_changed(self, index):
        """Handle file selection change"""
        if index < 0 or index >= len(self.file_list):
            return
        filename = self.file_combo.itemText(index)
        if filename == self.current_filename and self.image_data is not None:
            return
        self.load_selected_image(keep_view=False)

    def measure_star(self, x, y):
        """Measure star properties at (x, y)"""
        if self.image_data is None:
            return

        try:

            # Get local region for FWHM estimation
            size = 25
            y0 = max(0, y - size)
            y1 = min(self.image_data.shape[0], y + size)
            x0 = max(0, x - size)
            x1 = min(self.image_data.shape[1], x + size)

            cutout = self.image_data[y0:y1, x0:x1]
            cutout_min = float(np.min(cutout)) if cutout.size else float("nan")
            cutout_max = float(np.max(cutout)) if cutout.size else float("nan")

            # Estimate local background
            mean, median, std = sigma_clipped_stats(cutout, sigma=self._sigma_clip_sigma())

            # Convert FWHM seed from arcsec to pixels
            pixscale = self.params.P.pixel_scale_arcsec
            fwhm_seed_px = self.fwhm_seed_arcsec / pixscale if np.isfinite(pixscale) and pixscale > 0 else 6.0

            # Estimate FWHM using radial profile method (more accurate per-star)
            fwhm = fwhm_seed_px
            fwhm_method = "seed"
            try:
                # Find centroid using DAOStarFinder
                import warnings
                from photutils.utils.exceptions import NoDetectionsWarning

                xc, yc = float(x), float(y)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", NoDetectionsWarning)
                    finder = DAOStarFinder(threshold=median + 3*std, fwhm=fwhm_seed_px)
                    sources = finder(cutout - median)

                if sources is not None and len(sources) > 0:
                    # Find closest source to click
                    dx = sources['xcentroid'] - (x - x0)
                    dy = sources['ycentroid'] - (y - y0)
                    dist = np.sqrt(dx**2 + dy**2)
                    idx = np.argmin(dist)
                    if dist[idx] < 10:
                        xc = float(x0 + sources['xcentroid'][idx])
                        yc = float(y0 + sources['ycentroid'][idx])

                # Calculate FWHM from radial profile
                rmax = 30
                ry0 = max(0, int(yc - rmax))
                ry1 = min(self.image_data.shape[0], int(yc + rmax))
                rx0 = max(0, int(xc - rmax))
                rx1 = min(self.image_data.shape[1], int(xc + rmax))

                region = self.image_data[ry0:ry1, rx0:rx1]
                yy, xx = np.mgrid[ry0:ry1, rx0:rx1]
                rr = np.sqrt((xx - xc)**2 + (yy - yc)**2)

                # Background subtraction
                _, reg_median, _ = sigma_clipped_stats(region, sigma=self._sigma_clip_sigma())
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

                # Estimate FWHM from profile
                peak = np.nanmax(profile) if np.isfinite(profile).any() else 0
                if peak > 0:
                    half = 0.5 * peak
                    idx_half = np.where((profile[:-1] >= half) & (profile[1:] < half))[0]
                    if len(idx_half) > 0:
                        i = idx_half[0]
                        x1_p, y1_p = centers[i], profile[i]
                        x2_p, y2_p = centers[i+1], profile[i+1]
                        if y1_p != y2_p:
                            r_half = x1_p + (half - y1_p) * (x2_p - x1_p) / (y2_p - y1_p)
                            fwhm = 2.0 * r_half
                            fwhm_method = "radial"
            except:
                pass

            # Define aperture and annulus using FWHM seed × scale (fixed size, not per-star)
            # Measured FWHM is for display only, aperture size uses seed value
            r_ap = max(self.min_r_ap_px, fwhm_seed_px * self.aperture_scale)
            r_in = max(self.min_r_in_px, fwhm_seed_px * self.annulus_in_scale)
            r_out = max(self.min_r_out_px, r_in + fwhm_seed_px * self.annulus_out_scale)
            self.last_measurement = {
                "x": float(x),
                "y": float(y),
                "r_ap": float(r_ap),
                "r_in": float(r_in),
                "r_out": float(r_out),
            }

            aperture = CircularAperture((x, y), r=r_ap)
            annulus = CircularAnnulus((x, y), r_in=r_in, r_out=r_out)

            # Calculate stats
            ap_stats = ApertureStats(self.image_data, aperture)
            an_stats = ApertureStats(self.image_data, annulus, sigma_clip=None)

            # Sky background - handle both scalar and array returns
            sky_median = float(an_stats.median[0]) if hasattr(an_stats.median, '__len__') else float(an_stats.median)
            sky_std = float(an_stats.std[0]) if hasattr(an_stats.std, '__len__') else float(an_stats.std)

            # Source flux - handle both scalar and array returns
            flux_total = float(ap_stats.sum[0]) if hasattr(ap_stats.sum, '__len__') else float(ap_stats.sum)
            flux_sky = sky_median * float(ap_stats.sum_aper_area.value)
            flux_source = flux_total - flux_sky

            # SNR calculation
            noise = resolve_effective_noise_params(self.params.P, self.header)
            gain = noise.gain_e_per_adu
            rdnoise = noise.rdnoise_e

            flux_e = flux_source * gain
            sky_e = flux_sky * gain
            n_pix = ap_stats.sum_aper_area.value

            snr = flux_e / np.sqrt(flux_e + sky_e + n_pix * rdnoise**2)

            # Magnitude
            exptime = self.header.get('EXPTIME', 1.0)
            zp = self.params.P.zp_initial
            mag = -2.5 * np.log10(flux_source / exptime) + zp if flux_source > 0 else float("nan")
            mag_err = MAG_ERR_COEFF / snr if snr > 0 else float("nan")  # 2.5/ln(10)
            mag_str = f"{mag:.2f}" if np.isfinite(mag) else "NaN"
            mag_err_str = f"{mag_err:.3f}" if np.isfinite(mag_err) else "NaN"

            # FWHM in arcsec
            pixscale = self.params.P.pixel_scale_arcsec
            fwhm_arcsec = fwhm * pixscale if np.isfinite(pixscale) else np.nan

            # Display results with FWHM in arcsec first
            stats_text = f"""
Measurement at ({x}, {y})
─────────────────────────
File: {self.current_filename}
EXPTIME: {exptime:.1f} s
Local Pixel Min/Max: {cutout_min:.2f} / {cutout_max:.2f} ADU

Sky Background:
  Median: {sky_median:.2f} ADU
  StdDev: {sky_std:.2f} ADU

Aperture (r={r_ap:.1f} px):
  Total: {flux_total:.0f} ADU
  Sky:   {flux_sky:.0f} ADU
  Source: {flux_source:.0f} ADU

FWHM: {fwhm_arcsec:.2f}" ({fwhm:.2f} px; {fwhm_method})
SNR:  {snr:.1f}
Mag:  {mag_str} ± {mag_err_str}
"""
            self.stats_text.setPlainText(stats_text)

            # Draw aperture/annulus overlay markers on GL viewer
            if self._viewer is not None:
                self._viewer.set_overlay_markers([
                    OverlayMarker(col=float(x), row=float(y), radius=r_ap,
                                  color=QColor(255, 60,  60, 220)),
                    OverlayMarker(col=float(x), row=float(y), radius=r_in,
                                  color=QColor(0,  200, 220, 160)),
                    OverlayMarker(col=float(x), row=float(y), radius=r_out,
                                  color=QColor(0,  200, 220, 100)),
                ])

        except Exception as e:
            self.stats_text.setPlainText(f"Measurement failed:\n{str(e)}")

    def _restore_measurement_overlay(self):
        if self.last_measurement is None or self.image_data is None:
            return
        x = int(self.last_measurement["x"])
        y = int(self.last_measurement["y"])
        if x < 0 or y < 0 or y >= self.image_data.shape[0] or x >= self.image_data.shape[1]:
            return
        self.measure_star(x, y)

    # ========== Keyboard Event Handling ==========

    def keyPressEvent(self, event):
        """Handle Qt keyboard events (overridden from QWidget)"""
        key = event.key()
        text = event.text().lower()

        if text == 'm' and self.cursor_x is not None and self.cursor_y is not None:
            self.measure_star_at_cursor()
        elif text == 'h' and self.cursor_x is not None and self.cursor_y is not None:
            self.show_histogram()
        elif text == 'g' and self.cursor_x is not None and self.cursor_y is not None:
            self.show_radial_profile()
        elif key == Qt.Key_Period or text == '.':  # . key
            self.cycle_filter()
        elif key == Qt.Key_BracketLeft:  # '['
            self.navigate_frame(-1)
        elif key == Qt.Key_BracketRight:  # ']'
            self.navigate_frame(1)
        else:
            # Call parent implementation for other keys
            super().keyPressEvent(event)

    def measure_star_at_cursor(self):
        """Measure star properties at current cursor position (m key)"""
        if self.image_data is None or self.cursor_x is None or self.cursor_y is None:
            return

        x = int(self.cursor_x)
        y = int(self.cursor_y)

        # Use the existing measure_star method
        self.measure_star(x, y)

    def navigate_frame(self, direction):
        """Navigate to previous/next frame within the SAME filter"""
        if not self.file_list:
            return

        # Ensure filter map is built
        self._build_frame_key_map()

        # Get current filter from header if possible
        current_filter = self._file_filter_map.get(self.current_filename)
        if current_filter is None and self.header is not None:
            current_filter = self._extract_filter_from_header(self.header)
            if current_filter:
                self._file_filter_map[self.current_filename] = current_filter

        if not current_filter:
            current_filter = ""

        # Build list of file indices for current filter
        filter_indices = []
        for idx, fname in enumerate(self.file_list):
            fkey = self._file_filter_map.get(fname, "")
            if fkey == current_filter:
                filter_indices.append(idx)

        # Debug: if filter_indices is empty or has only 1 item, try filename-based filter
        if len(filter_indices) <= 1 and current_filter:
            # Try matching by filename pattern (e.g., files ending with _G, _B, _I, etc.)
            filter_indices = []
            filter_suffix = f"_{current_filter.lower()}"
            for idx, fname in enumerate(self.file_list):
                base = fname.lower()
                for ext in ('.fits', '.fit', '.fts'):
                    if base.endswith(ext):
                        base = base[:-len(ext)]
                        break
                if base.endswith(filter_suffix) or f"_{current_filter.lower()}." in fname.lower():
                    filter_indices.append(idx)

        if not filter_indices or len(filter_indices) <= 1:
            # Fallback: cycle through all files
            self.current_file_index = (self.current_file_index + direction) % len(self.file_list)
        else:
            # Find current position within filter's files
            try:
                pos = filter_indices.index(self.current_file_index)
            except ValueError:
                pos = 0
            # Move within filter only
            pos = (pos + direction) % len(filter_indices)
            self.current_file_index = filter_indices[pos]

        self.file_combo.setCurrentIndex(self.current_file_index)
        self.load_selected_image(keep_view=True)

    def cycle_filter(self):
        """Cycle to next filter, keeping the same frame position (index within filter)"""
        if not self.file_list or self.image_data is None:
            return
        self._build_frame_key_map()
        current_filter = self._file_filter_map.get(self.current_filename, "")
        filters = [f for f in self._filter_order if f]
        if len(filters) <= 1:
            self.stats_text.append("\n[Filter Cycle] Only one filter found")
            return

        # Find next filter
        try:
            current_filter_idx = filters.index(current_filter)
            next_filter = filters[(current_filter_idx + 1) % len(filters)]
        except ValueError:
            next_filter = filters[0]

        # Build file indices for current and next filter
        current_filter_indices = []
        next_filter_indices = []
        for idx, fname in enumerate(self.file_list):
            fkey = self._file_filter_map.get(fname, "")
            if fkey == current_filter:
                current_filter_indices.append(idx)
            if fkey == next_filter:
                next_filter_indices.append(idx)

        if not next_filter_indices:
            return

        # Find current position within current filter's files
        try:
            pos_in_filter = current_filter_indices.index(self.current_file_index)
        except ValueError:
            pos_in_filter = 0

        # Go to same position in next filter (or last if shorter)
        target_pos = min(pos_in_filter, len(next_filter_indices) - 1)
        self.current_file_index = next_filter_indices[target_pos]

        self.file_combo.setCurrentIndex(self.current_file_index)
        self.load_selected_image(keep_view=True)
        self.stats_text.append(f"\n[Filter Cycle] Switched to {next_filter}")

    @staticmethod
    def _normalize_filter_key(value: str | None) -> str:
        return _canonical_filter_key(value)

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
            cand = _canonical_filter_key(token)
            if cand and 1 <= len(cand) <= 3 and str(cand).isalpha():
                return cand
        return ""

    def _extract_frame_key(self, fname: str, filter_key: str) -> str:
        name = Path(fname).name
        base = name
        for ext in (".fits", ".fit", ".fts", ".fz", ".gz"):
            if base.lower().endswith(ext):
                base = base[: -len(ext)]
        base_lower = base.lower()
        fkey = self._normalize_filter_key(filter_key) or self._infer_filter_from_filename(fname)
        if fkey:
            filt = fkey.lower()
            for sep in ("-", "_", "."):
                suffix = f"{sep}{filt}"
                if base_lower.endswith(suffix):
                    return base[: -len(suffix)]
        return base

    def _load_filter_map_from_index(self):
        # 1) headers.csv에서 먼저 로드 (Step 1에서 생성됨 - FITS 안 읽어도 됨)
        headers_path = step1_dir(self.params.P.result_dir) / "headers.csv"
        if not headers_path.exists():
            headers_path = self.params.P.result_dir / "headers.csv"
        if headers_path.exists():
            try:
                hdf = pd.read_csv(headers_path)
                if "Filename" in hdf.columns and "FILTER" in hdf.columns:
                    for fname, filt in zip(hdf["Filename"].astype(str), hdf["FILTER"].astype(str)):
                        fkey = self._normalize_filter_key(filt)
                        if fkey:
                            self._file_filter_map[fname] = fkey
                            if fname.lower().startswith("crop_"):
                                self._file_filter_map[fname[5:]] = fkey
                            if fname.lower().startswith("cropped_"):
                                self._file_filter_map[fname[8:]] = fkey
            except Exception:
                pass

        # 2) photometry_index.csv에서 보충 (Step 7 forced phot 완료 후)
        _result = self.params.P.result_dir
        idx_path = step7_forced_phot_dir(_result) / "photometry_index.csv"
        if not idx_path.exists():
            return
        try:
            df = pd.read_csv(idx_path)
        except Exception:
            return
        if "file" not in df.columns or "filter" not in df.columns:
            return
        for fname, filt in zip(df["file"].astype(str), df["filter"].astype(str)):
            if fname in self._file_filter_map:
                continue  # 이미 있으면 스킵
            fkey = self._normalize_filter_key(filt)
            if not fkey:
                fkey = self._infer_filter_from_filename(fname)
            self._file_filter_map[fname] = fkey
            if fname.lower().startswith("crop_"):
                self._file_filter_map[fname[5:]] = fkey
            if fname.lower().startswith("cropped_"):
                self._file_filter_map[fname[8:]] = fkey

    def _build_frame_key_map(self):
        if self._frame_key_map:
            return
        for idx, fname in enumerate(self.file_list):
            fkey = self._file_filter_map.get(fname)
            if not fkey:
                try:
                    if self.use_cropped and self.cropped_dir is not None:
                        fpath = self.cropped_dir / fname
                    else:
                        fpath = self.params.get_file_path(fname)
                    hdr = fits.getheader(fpath)
                    fkey = self._extract_filter_from_header(hdr)
                except Exception:
                    fkey = ""
                if not fkey:
                    fkey = self._infer_filter_from_filename(fname)
                self._file_filter_map[fname] = fkey
            frame_key = self._extract_frame_key(fname, fkey)
            self._file_frame_key_map[fname] = frame_key
            self._frame_key_map[(frame_key, fkey)] = idx
            self._frame_keys_by_filter.setdefault(fkey, []).append(frame_key)
        self._filter_order = [f for f in sorted(self._frame_keys_by_filter.keys()) if f]
    @staticmethod
    def _normalize_filter_key(value: str | None) -> str:
        if value is None:
            return ""
        return str(value).strip().upper()

    def show_histogram(self):
        """Show histogram of cursor region (h key) - persistent window"""
        if self.image_data is None or self.cursor_x is None or self.cursor_y is None:
            return

        x = int(self.cursor_x)
        y = int(self.cursor_y)

        # Get region around cursor
        size = 50
        y0 = max(0, y - size)
        y1 = min(self.image_data.shape[0], y + size)
        x0 = max(0, x - size)
        x1 = min(self.image_data.shape[1], x + size)

        region = self.image_data[y0:y1, x0:x1]
        region_finite = region[np.isfinite(region)]

        if region_finite.size == 0:
            return

        # Create or update histogram window
        if self.histogram_dialog is None or not self.histogram_dialog.isVisible():
            # Create new dialog
            self.histogram_dialog = FittedDialog(self)
            self.histogram_dialog.setWindowTitle("Histogram (Press 'h' to update)")
            self.histogram_dialog.resize(700, 500)
            self.histogram_dialog.setWindowFlags(Qt.Window)  # Make it resizable

            layout = QVBoxLayout(self.histogram_dialog)

            # Create matplotlib figure
            self.hist_fig = Figure(figsize=(7, 5))
            self.hist_canvas = FigureCanvas(self.hist_fig)
            self.hist_ax = self.hist_fig.add_subplot(111)

            layout.addWidget(NavigationToolbar(self.hist_canvas, self.histogram_dialog))
            layout.addWidget(tame_canvas(self.hist_canvas), 1)

            # Stats text
            self.hist_stats_label = QLabel()
            self.hist_stats_label.setStyleSheet("QLabel { font-family: monospace; padding: 10px; background-color: #f0f0f0; }")
            layout.addWidget(self.hist_stats_label)

            self.histogram_dialog.show()
            self._move_dialog_to_side(self.histogram_dialog)

        # Update plot
        self.hist_ax.clear()
        mean, median, std = sigma_clipped_stats(region_finite, sigma=self._sigma_clip_sigma())

        self.hist_ax.hist(region_finite.flatten(), bins=50, alpha=0.7, color='steelblue', edgecolor='black')
        self.hist_ax.axvline(median, color='red', linestyle='--', linewidth=2, label=f'Median: {median:.2f}')
        self.hist_ax.axvline(mean, color='orange', linestyle='--', linewidth=2, label=f'Mean: {mean:.2f}')

        self.hist_ax.set_xlabel('Pixel Value (ADU)')
        self.hist_ax.set_ylabel('Count')
        self.hist_ax.set_title(f'Histogram at ({x}, {y}) | {x0}:{x1}, {y0}:{y1} ({region_finite.size} pixels)')
        self.hist_ax.legend()
        self.hist_ax.grid(True, alpha=0.3)

        self.hist_fig.tight_layout()
        self.hist_canvas.draw()

        # Update stats
        self.hist_stats_label.setText(
            f"Mean: {mean:.2f} ADU  |  Median: {median:.2f} ADU  |  StdDev: {std:.2f} ADU  |  "
            f"Min: {np.min(region_finite):.2f} ADU  |  Max: {np.max(region_finite):.2f} ADU"
        )

    def show_radial_profile(self):
        """Show radial profile plot (g key)"""
        if self.image_data is None or self.cursor_x is None or self.cursor_y is None:
            return

        x = int(self.cursor_x)
        y = int(self.cursor_y)

        try:
            # Recenter on star
            size = 25
            y0 = max(0, y - size)
            y1 = min(self.image_data.shape[0], y + size)
            x0 = max(0, x - size)
            x1 = min(self.image_data.shape[1], x + size)

            cutout = self.image_data[y0:y1, x0:x1]
            mean, median, std = sigma_clipped_stats(cutout, sigma=self._sigma_clip_sigma())

            # Convert FWHM seed from arcsec to pixels
            pixscale = self.params.P.pixel_scale_arcsec
            fwhm_seed_px = self.fwhm_seed_arcsec / pixscale if np.isfinite(pixscale) and pixscale > 0 else 6.0

            # Try to find center
            xc, yc = float(x), float(y)
            try:
                import warnings
                from photutils.utils.exceptions import NoDetectionsWarning

                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", NoDetectionsWarning)
                    finder = DAOStarFinder(threshold=median + 3*std, fwhm=fwhm_seed_px)
                    sources = finder(cutout - median)

                if sources is not None and len(sources) > 0:
                    dx = sources['xcentroid'] - (x - x0)
                    dy = sources['ycentroid'] - (y - y0)
                    dist = np.sqrt(dx**2 + dy**2)
                    idx = np.argmin(dist)
                    if dist[idx] < 10:
                        xc = float(x0 + sources['xcentroid'][idx])
                        yc = float(y0 + sources['ycentroid'][idx])
            except:
                pass

            # Calculate radial profile
            rmax = 50
            y0 = max(0, int(yc - rmax))
            y1 = min(self.image_data.shape[0], int(yc + rmax))
            x0 = max(0, int(xc - rmax))
            x1 = min(self.image_data.shape[1], int(xc + rmax))

            self.last_radial_center = (float(xc), float(yc))

            region = self.image_data[y0:y1, x0:x1]
            yy, xx = np.mgrid[y0:y1, x0:x1]
            rr = np.sqrt((xx - xc)**2 + (yy - yc)**2)

            # Background subtraction
            mean, median, std = sigma_clipped_stats(region, sigma=self._sigma_clip_sigma())
            region_sub = region - median

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

            # Create or reuse radial profile window
            if self.radial_profile_dialog is None or not self.radial_profile_dialog.isVisible():
                self.radial_profile_dialog = FittedDialog(self)
                self.radial_profile_dialog.setWindowTitle("Radial Profile (Press 'g' to update)")
                self.radial_profile_dialog.resize(700, 500)
                self.radial_profile_dialog.setWindowFlags(Qt.Window)  # Make it resizable

                layout = QVBoxLayout(self.radial_profile_dialog)

                self.prof_fig = Figure(figsize=(7, 5))
                self.prof_canvas = FigureCanvas(self.prof_fig)
                self.prof_ax = self.prof_fig.add_subplot(111)

                layout.addWidget(NavigationToolbar(self.prof_canvas, self.radial_profile_dialog))
                layout.addWidget(tame_canvas(self.prof_canvas), 1)

                close_btn = QPushButton("Close")
                close_btn.clicked.connect(self.radial_profile_dialog.close)
                layout.addWidget(close_btn)

                self.radial_profile_dialog.show()
                self._move_dialog_to_side(self.radial_profile_dialog)

            # Plot profile
            self.prof_ax.clear()
            self.prof_ax.plot(centers, profile, 'o-', color='steelblue', markersize=4, linewidth=1.5)
            self.prof_ax.axhline(0, color='red', linestyle='--', linewidth=1, alpha=0.5)
            self.prof_ax.set_xlabel('Radius (pixels)')
            self.prof_ax.set_ylabel('Pixel Value - Background (ADU)')
            self.prof_ax.set_title(f'Radial Profile (centered at {int(xc)}, {int(yc)})')
            self.prof_ax.grid(True, alpha=0.3)

            # Estimate FWHM from profile
            self.last_radial_fwhm_px = None
            fwhm_note = "FWHM: N/A"
            pixscale = self.params.P.pixel_scale_arcsec
            peak = np.nanmax(profile) if np.isfinite(profile).any() else 0
            if peak > 0:
                half = 0.5 * peak
                idx = np.where((profile[:-1] >= half) & (profile[1:] < half))[0]
                if len(idx) > 0:
                    i = idx[0]
                    x1, y1 = centers[i], profile[i]
                    x2, y2 = centers[i+1], profile[i+1]
                    if y1 != y2:
                        r_half = x1 + (half - y1) * (x2 - x1) / (y2 - y1)
                        fwhm_est = 2.0 * r_half
                        self.last_radial_fwhm_px = float(fwhm_est)
                        fwhm_arcsec = fwhm_est * pixscale if np.isfinite(pixscale) else np.nan
                        if np.isfinite(fwhm_arcsec):
                            fwhm_note = f'FWHM: {fwhm_arcsec:.2f}" ({fwhm_est:.2f} px; radial)'
                        else:
                            fwhm_note = f"FWHM: {fwhm_est:.2f} px (radial)"
                        self.prof_ax.axvline(r_half, color='orange', linestyle='--', linewidth=2, label=f'FWHM ~{fwhm_arcsec:.2f}" ({fwhm_est:.2f} px)')
                        self.prof_ax.axhline(half, color='green', linestyle=':', linewidth=1, alpha=0.7)
                        self.prof_ax.legend()

            self.prof_fig.tight_layout()
            self.prof_ax.text(0.02, 0.95, fwhm_note, transform=self.prof_ax.transAxes,
                              ha='left', va='top', fontsize=10,
                              bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))
            self.prof_canvas.draw()

        except Exception as e:
            self.stats_text.append(f"\n[Radial Profile] Error: {str(e)}")

    def _move_dialog_to_side(self, dialog):
        """Place tool dialog to the right of the main window."""
        if dialog is None:
            return
        main_geo = self.frameGeometry()
        screen = QApplication.primaryScreen()
        if screen is None:
            return
        avail = screen.availableGeometry()
        x = min(main_geo.topRight().x() + 10, avail.right() - dialog.width())
        y = min(main_geo.top(), avail.bottom() - dialog.height())
        dialog.move(x, y)

    # ========== Parameter Dialog ==========

    def open_parameters_dialog(self):
        dialog, layout, buttons = build_scroll_param_dialog(
            self, "Photometry Parameters",
            info_text="Adjust aperture photometry parameters. Changes apply immediately to measurements.",
            size=(440, 400),
        )

        form = QFormLayout()

        spin_fwhm = QDoubleSpinBox()
        spin_fwhm.setRange(0.1, 20.0)
        spin_fwhm.setSingleStep(0.1)
        spin_fwhm.setDecimals(2)
        spin_fwhm.setValue(self.fwhm_seed_arcsec)
        form.addRow("FWHM Seed (arcsec):", spin_fwhm)

        spin_ap = QDoubleSpinBox()
        spin_ap.setRange(0.5, 10.0)
        spin_ap.setSingleStep(0.1)
        spin_ap.setDecimals(2)
        spin_ap.setValue(self.aperture_scale)
        form.addRow("Aperture Scale (× FWHM):", spin_ap)

        spin_ann_in = QDoubleSpinBox()
        spin_ann_in.setRange(1.0, 20.0)
        spin_ann_in.setSingleStep(0.1)
        spin_ann_in.setDecimals(2)
        spin_ann_in.setValue(self.annulus_in_scale)
        form.addRow("Annulus Inner Scale (× FWHM):", spin_ann_in)

        spin_ann_out = QDoubleSpinBox()
        spin_ann_out.setRange(0.5, 20.0)
        spin_ann_out.setSingleStep(0.1)
        spin_ann_out.setDecimals(2)
        spin_ann_out.setValue(self.annulus_out_scale)
        form.addRow("Annulus Outer Width (× FWHM):", spin_ann_out)

        spin_sigma = QDoubleSpinBox()
        spin_sigma.setRange(1.0, 10.0)
        spin_sigma.setSingleStep(0.1)
        spin_sigma.setDecimals(1)
        spin_sigma.setValue(self.sigma_clip_value)
        form.addRow("Sigma Clipping (σ):", spin_sigma)

        spin_zp = QDoubleSpinBox()
        spin_zp.setRange(0.0, 50.0)
        spin_zp.setSingleStep(0.1)
        spin_zp.setDecimals(2)
        spin_zp.setValue(float(getattr(self.params.P, "zp_initial", 25.0)))
        spin_zp.setToolTip("Display-only zero point for Step 3 preview magnitude.")
        form.addRow("Preview ZP:", spin_zp)

        layout.addLayout(form)
        layout.addStretch(1)

        add_parameter_reset_button(
            buttons,
            [
                (spin_fwhm, 2.5),
                (spin_ap, 1.0),
                (spin_ann_in, 4.0),
                (spin_ann_out, 2.0),
                (spin_sigma, 3.0),
                (spin_zp, 25.0),
            ],
        )

        def _save():
            self.fwhm_seed_arcsec = spin_fwhm.value()
            self.aperture_scale = spin_ap.value()
            self.annulus_in_scale = spin_ann_in.value()
            self.annulus_out_scale = spin_ann_out.value()
            self.min_r_ap_px = _MIN_R_AP_PX
            self.min_r_in_px = _MIN_R_IN_PX
            self.min_r_out_px = _MIN_R_OUT_PX
            self.sigma_clip_value = spin_sigma.value()
            self.params.P.fwhm_guess_arcsec = self.fwhm_seed_arcsec
            self.params.P.phot_aperture_scale = self.aperture_scale
            self.params.P.fitsky_annulus_scale = self.annulus_in_scale
            self.params.P.fitsky_dannulus_scale = self.annulus_out_scale
            self.params.P.annulus_sigma_clip = self.sigma_clip_value
            self.params.P.zp_initial = spin_zp.value()
            self.persist_params()
            self.save_state()
            QMessageBox.information(dialog, "Saved", "Parameters updated and saved!")
            dialog.accept()

        buttons.accepted.connect(_save)
        buttons.rejected.connect(dialog.reject)
        dialog.exec_()

    def open_stretch_plot(self):
        """Open stretch plot window showing histogram with draggable min/max markers"""
        if self.image_data is None:
            QMessageBox.warning(self, "Warning", "Load an image first")
            return

        # If dialog exists and is visible, just raise it
        if self.stretch_plot_dialog is not None and self.stretch_plot_dialog.isVisible():
            self.stretch_plot_dialog.raise_()
            self.stretch_plot_dialog.activateWindow()
            self._update_stretch_plot()
            return

        # Create dialog
        self.stretch_plot_dialog = FittedDialog(self)
        self.stretch_plot_dialog.setWindowTitle("2D Plot - Stretch Control")
        self.stretch_plot_dialog.resize(500, 250)

        layout = QVBoxLayout(self.stretch_plot_dialog)

        # Info label
        self.stretch_plot_info_label = QLabel("Drag min/max markers to adjust stretch")
        self.stretch_plot_info_label.setStyleSheet("QLabel { padding: 5px; background-color: #E3F2FD; border-radius: 3px; }")
        layout.addWidget(self.stretch_plot_info_label)

        # Matplotlib figure for histogram
        self.stretch_plot_fig = Figure(figsize=(6, 2.5))
        self.stretch_plot_canvas = FigureCanvas(self.stretch_plot_fig)
        self.stretch_plot_ax = self.stretch_plot_fig.add_subplot(111)
        self.stretch_plot_fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15, top=0.9)

        # Connect mouse events
        self.stretch_plot_canvas.mpl_connect('button_press_event', self._on_stretch_plot_press)
        self.stretch_plot_canvas.mpl_connect('motion_notify_event', self._on_stretch_plot_motion)
        self.stretch_plot_canvas.mpl_connect('button_release_event', self._on_stretch_plot_release)

        layout.addWidget(tame_canvas(self.stretch_plot_canvas, min_h=140), 1)

        # Hint label
        hint_label = QLabel("Click and drag < > markers to adjust min/max | Changes apply in real-time")
        hint_label.setStyleSheet("QLabel { color: #666; font-size: 10px; }")
        layout.addWidget(hint_label)

        self.stretch_plot_dialog.show()
        self._update_stretch_plot()

    def _update_stretch_plot(self):
        """Update the stretch plot histogram and markers"""
        if self.stretch_plot_ax is None or self.image_data is None:
            return

        ax = self.stretch_plot_ax
        ax.clear()

        # Get image data statistics
        data = self.image_data.copy()
        finite_mask = np.isfinite(data)
        if not finite_mask.any():
            return

        flat = data[finite_mask].flatten()

        # Use percentile range for histogram display
        p_low, p_high = np.percentile(flat, [1, 99])
        display_data = flat[(flat >= p_low) & (flat <= p_high)]

        if len(display_data) == 0:
            display_data = flat

        # Store data range for marker calculations
        self._stretch_data_range = (float(p_low), float(p_high))

        # Initial vmin/vmax: use current viewer STF params if available, else compute
        if self._stretch_vmin is None or self._stretch_vmax is None:
            if self._viewer is not None:
                s, h, _ = self._viewer.get_stf_params()
                vmin, vmax = s, h
            else:
                _, median_val, std_val = sigma_clipped_stats(
                    flat, sigma=self._sigma_clip_sigma(), maxiters=5
                )
                vmin = max(np.min(flat), median_val - 2.8 * std_val)
                vmax = min(np.max(flat), np.percentile(flat, 99.9))
            if vmax <= vmin:
                vmin = np.min(flat)
                vmax = np.max(flat)
            self._stretch_vmin = float(vmin)
            self._stretch_vmax = float(vmax)

        # Plot histogram
        ax.hist(display_data, bins=128, color='#3a6ea5', edgecolor='none', alpha=0.7)
        ax.set_xlim(p_low, p_high)

        # Draw min/max markers
        vmin = self._stretch_vmin
        vmax = self._stretch_vmax

        # Clamp markers to display range for visualization
        vmin_display = max(p_low, min(p_high, vmin))
        vmax_display = max(p_low, min(p_high, vmax))

        self._stretch_marker_min_line = ax.axvline(vmin_display, color='#FF5722', linewidth=2, linestyle='-', label=f"Min: {vmin:.1f}")
        self._stretch_marker_max_line = ax.axvline(vmax_display, color='#4CAF50', linewidth=2, linestyle='-', label=f"Max: {vmax:.1f}")

        # Add marker text
        y_max = ax.get_ylim()[1]
        ax.text(vmin_display, y_max * 0.95, '◀', color='#FF5722', fontsize=14, ha='center', va='top', fontweight='bold')
        ax.text(vmax_display, y_max * 0.95, '▶', color='#4CAF50', fontsize=14, ha='center', va='top', fontweight='bold')

        ax.set_xlabel('Pixel Value')
        ax.set_ylabel('Count')
        ax.set_title('Image Histogram')
        ax.legend(loc='upper right', fontsize=8)

        # Update info label
        if self.stretch_plot_info_label:
            mode = self._viewer._mode_combo.currentText() if self._viewer else "STF"
            self.stretch_plot_info_label.setText(
                f"Stretch: {mode} | Min: {vmin:.2f} | Max: {vmax:.2f}"
            )

        self.stretch_plot_canvas.draw_idle()

    def _on_stretch_plot_press(self, event):
        """Handle mouse press on stretch plot"""
        if event.inaxes != self.stretch_plot_ax or event.xdata is None:
            return

        if self._stretch_vmin is None or self._stretch_vmax is None:
            return

        x = event.xdata
        dist_to_min = abs(x - self._stretch_vmin)
        dist_to_max = abs(x - self._stretch_vmax)

        # Select closer marker
        if dist_to_min < dist_to_max:
            self._stretch_drag_target = "min"
        else:
            self._stretch_drag_target = "max"

        self._stretch_dragging = True

    def _on_stretch_plot_motion(self, event):
        """Handle mouse motion on stretch plot (dragging)"""
        if not self._stretch_dragging or event.xdata is None:
            return

        x = event.xdata

        if self._stretch_drag_target == "min":
            # Don't let min exceed max
            new_val = min(x, self._stretch_vmax - 1)
            self._stretch_vmin = new_val
        else:
            # Don't let max go below min
            new_val = max(x, self._stretch_vmin + 1)
            self._stretch_vmax = new_val

        # Update plot and image in real-time
        self._update_stretch_plot()
        self._apply_custom_stretch()

    def _on_stretch_plot_release(self, event):
        """Handle mouse release on stretch plot"""
        self._stretch_dragging = False
        self._stretch_drag_target = None

    def _apply_custom_stretch(self):
        """Apply custom vmin/vmax from stretch_plot drag markers."""
        if self._viewer is None or self._stretch_vmin is None or self._stretch_vmax is None:
            return
        vmin = self._stretch_vmin
        vmax = self._stretch_vmax if self._stretch_vmax > self._stretch_vmin else self._stretch_vmin + 1
        self._viewer.set_shadow_highlight(vmin, vmax)

    def reset_stretch_plot_values(self):
        """Reset stretch plot values when changing image or stretch mode"""
        self._stretch_vmin = None
        self._stretch_vmax = None
        if self.stretch_plot_dialog and self.stretch_plot_dialog.isVisible():
            self._update_stretch_plot()

    def validate_step(self) -> bool:
        """Validate if step can be completed"""
        # Step 3 is optional QC, always allow completion
        return True

    def save_state(self):
        """Save step state including photometry parameters"""
        state_data = {
            "current_file": self.current_filename,
            "use_cropped": self.use_cropped,
            # Photometry parameters (scale-based)
            "fwhm_seed_arcsec": self.fwhm_seed_arcsec,
            "aperture_scale": self.aperture_scale,
            "annulus_in_scale": self.annulus_in_scale,
            "annulus_out_scale": self.annulus_out_scale,
            "sigma_clip_value": self.sigma_clip_value,
        }
        self.project_state.store_step_data("sky_preview", state_data)

    def restore_state(self):
        """Restore step state including photometry parameters"""
        state_data = self.project_state.get_step_data("sky_preview")
        if state_data:
            # Restore photometry parameters
            if "fwhm_seed_arcsec" in state_data:
                self.fwhm_seed_arcsec = float(state_data["fwhm_seed_arcsec"])
            if "aperture_scale" in state_data:
                self.aperture_scale = float(state_data["aperture_scale"])
            if "annulus_in_scale" in state_data:
                self.annulus_in_scale = float(state_data["annulus_in_scale"])
            if "annulus_out_scale" in state_data:
                self.annulus_out_scale = float(state_data["annulus_out_scale"])
            self.min_r_ap_px = _MIN_R_AP_PX
            self.min_r_in_px = _MIN_R_IN_PX
            self.min_r_out_px = _MIN_R_OUT_PX
            if "sigma_clip_value" in state_data:
                self.sigma_clip_value = float(state_data["sigma_clip_value"])
