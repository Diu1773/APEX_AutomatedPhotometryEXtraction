"""
Step 2: Image Crop Window
Allows user to select crop region on reference image and apply crop to all files
"""

from PyQt5.QtWidgets import (
    QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QGroupBox, QMessageBox, QProgressDialog, QApplication
)
from PyQt5.QtCore import Qt
from pathlib import Path
import json
import shutil
import numpy as np
from astropy.io import fits
from concurrent.futures import ThreadPoolExecutor, as_completed

from .step_window_base import StepWindowBase
from apex.gui.widgets.fits_viewer import FITSViewerWidget
from apex.utils.constants import get_parallel_workers
from apex.utils.step_paths import crop_rect_path, step2_dir, step2_cropped_dir


class CropSelectorWindow(StepWindowBase):
    """Step 2: Image Crop Selection"""

    def __init__(self, params, file_manager, project_state, main_window):
        """
        Initialize crop selector window

        Args:
            params: Parameters object
            file_manager: FileManager instance
            project_state: ProjectState object
            main_window: Main window reference
        """
        self.file_manager = file_manager

        # Crop rectangle coordinates
        self.crop_x0 = None
        self.crop_y0 = None
        self.crop_x1 = None
        self.crop_y1 = None

        # Viewer + state
        self._viewer: FITSViewerWidget | None = None
        self.original_image_data = None  # Original uncropped image
        self.displayed_image_data = None  # Currently displayed (may be cropped)
        self.ref_filename = None
        self.is_cropped = False  # Track if images have been cropped
        self.crop_skipped = False

        # Initialize base class
        super().__init__(
            step_index=1,
            step_name="Image Crop",
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

        # === Info Label ===
        info_label = QLabel(
            "1. Left-click drag on the image to draw a crop rectangle\n"
            "   (Right-click drag to pan, mouse wheel to zoom)\n"
            "2. Click 'Apply Crop' to save cropped images to 'result/step2_crop/cropped/' folder"
        )
        info_label.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 10px; border-radius: 5px; }")
        self.content_layout.addWidget(info_label)

        # === Image Viewer ===
        viewer_group = QGroupBox("Reference Image Viewer")
        viewer_layout = QVBoxLayout(viewer_group)

        self._viewer = FITSViewerWidget(self)
        self._viewer.set_selection_mode(True)
        self._viewer.selection_finished.connect(self._on_selection_finished)
        viewer_layout.addWidget(self._viewer)

        # Controls below image
        image_controls_layout = QHBoxLayout()

        self.btn_crop = QPushButton("Apply Crop")
        self.btn_crop.setEnabled(False)
        self.btn_crop.setMinimumHeight(35)
        self.btn_crop.setStyleSheet("""
            QPushButton {
                background-color: #FF5722;
                color: white;
                font-weight: bold;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #E64A19;
            }
            QPushButton:disabled {
                background-color: #CCCCCC;
                color: #666666;
            }
        """)
        self.btn_crop.clicked.connect(self.apply_crop)
        image_controls_layout.addWidget(self.btn_crop)

        self.btn_reset = QPushButton("Reset to Original")
        self.btn_reset.setMinimumHeight(35)
        self.btn_reset.clicked.connect(self.reset_to_original)
        image_controls_layout.addWidget(self.btn_reset)

        self.btn_skip = QPushButton("Skip Crop")
        self.btn_skip.setMinimumHeight(35)
        self.btn_skip.setStyleSheet("""
            QPushButton {
                background-color: #607D8B;
                color: white;
                font-weight: bold;
                border-radius: 5px;
            }
            QPushButton:hover {
                background-color: #546E7A;
            }
        """)
        self.btn_skip.clicked.connect(self.skip_crop)
        image_controls_layout.addWidget(self.btn_skip)

        image_controls_layout.addStretch()

        viewer_layout.addLayout(image_controls_layout)
        self.content_layout.addWidget(viewer_group)

        # === Crop Info ===
        crop_info_group = QGroupBox("Status")
        crop_info_layout = QVBoxLayout(crop_info_group)

        self.crop_info_label = QLabel("No image loaded")
        self.crop_info_label.setStyleSheet("QLabel { font-weight: bold; color: blue; }")
        crop_info_layout.addWidget(self.crop_info_label)

        self.content_layout.addWidget(crop_info_group)

    def load_reference_image(self):
        """Load reference image from file manager"""
        try:
            # Check if files are available from Step 1
            if not self.file_manager.filenames:
                # Try to scan files
                try:
                    self.file_manager.scan_files()
                except:
                    QMessageBox.warning(
                        self, "No Files",
                        "Run Step 1 first to scan files."
                    )
                    return

            # Use first available file
            if self.ref_filename and self.ref_filename in self.file_manager.filenames:
                pass  # keep previously used file for this step
            else:
                self.ref_filename = self.file_manager.filenames[0]

            file_path = self.params.get_file_path(self.ref_filename)

            # Load FITS data from original
            data = None
            with fits.open(file_path) as hdul:
                for hdu in hdul:
                    if hdu.data is not None:
                        data = hdu.data
                        break
            if data is None:
                raise RuntimeError("No image data found in FITS file")
            self.original_image_data = data.astype(float).copy()

            self.displayed_image_data = self.original_image_data.copy()
            self.is_cropped = False

            # Display image
            self.display_image()

            # Update status
            self.crop_info_label.setText(f"Loaded: {self.ref_filename}\nReady to draw crop region")

            self.update_navigation_buttons()

        except Exception as e:
            QMessageBox.critical(
                self, "Error",
                f"Failed to load reference image:\n{str(e)}"
            )

    def display_image(self):
        """Push current displayed_image_data to FITSViewerWidget."""
        if self.displayed_image_data is None or self._viewer is None:
            return
        data = np.asarray(self.displayed_image_data, dtype=np.float32)
        self._viewer.set_data(data)
        self._viewer.auto_stf()
        self._viewer.fit_in_view()

        # Restore saved selection rectangle (only on uncropped image)
        if self.crop_x0 is not None and not self.is_cropped:
            self._viewer.set_selection_rect(
                self.crop_x0, self.crop_y0, self.crop_x1, self.crop_y1
            )
            self._viewer.set_selection_mode(True)
        else:
            self._viewer.clear_selection_rect()
            # Disable selection drawing once cropped
            self._viewer.set_selection_mode(not self.is_cropped)

    def _on_selection_finished(self, x0: float, y0: float, x1: float, y1: float):
        """Handle rectangle finalize from FITSViewerWidget."""
        if self.is_cropped or self.original_image_data is None:
            return

        x0i, y0i = int(round(x0)), int(round(y0))
        x1i, y1i = int(round(x1)), int(round(y1))

        # Clamp to image bounds
        h, w = self.original_image_data.shape[:2]
        x0i = max(0, min(w, x0i))
        x1i = max(0, min(w, x1i))
        y0i = max(0, min(h, y0i))
        y1i = max(0, min(h, y1i))

        if (x1i - x0i) < 50 or (y1i - y0i) < 50:
            QMessageBox.warning(self, "Too Small", "Crop region must be at least 50×50 pixels")
            self.crop_x0 = None
            self.crop_y0 = None
            self.crop_x1 = None
            self.crop_y1 = None
            self._viewer.clear_selection_rect()
            return

        self.crop_x0, self.crop_y0 = x0i, y0i
        self.crop_x1, self.crop_y1 = x1i, y1i

        width  = self.crop_x1 - self.crop_x0
        height = self.crop_y1 - self.crop_y0
        self.crop_info_label.setText(
            f"Crop Region Selected: ({self.crop_x0}, {self.crop_y0}) → ({self.crop_x1}, {self.crop_y1})\n"
            f"Size: {width} × {height} pixels\n"
            f"Click 'Apply Crop' to crop all images"
        )
        self.btn_crop.setEnabled(True)
        self._viewer.set_selection_rect(self.crop_x0, self.crop_y0, self.crop_x1, self.crop_y1)

    def apply_crop(self):
        """Apply crop to all FITS files"""
        if not self.validate_crop_region():
            return
        if not self.confirm_recrop_warning():
            return

        # Confirm with user
        reply = QMessageBox.question(
            self, "Apply Crop",
            f"This will create cropped files in 'result/step2_crop/cropped' folder.\n"
            f"Original files will not be modified.\n\n"
            f"Crop {len(self.file_manager.filenames)} images\n"
            f"Region: ({self.crop_x0}, {self.crop_y0}) → ({self.crop_x1}, {self.crop_y1})\n\n"
            f"Continue?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )

        if reply != QMessageBox.Yes:
            return
        self.crop_skipped = False

        # Create progress dialog
        progress = QProgressDialog("Cropping images...", "Cancel", 0, len(self.file_manager.filenames), self)
        progress.setWindowTitle("Processing")
        progress.setWindowModality(Qt.WindowModal)

        try:
            # Create cropped directory only
            cropped_dir = step2_cropped_dir(self.params.P.result_dir)
            cropped_dir.mkdir(parents=True, exist_ok=True)

            files = list(self.file_manager.filenames)

            def crop_one(filename: str) -> str:
                original_path = self.params.get_file_path(filename)
                cropped_path = cropped_dir / filename

                if self._is_crop_cache_valid(original_path, cropped_path):
                    return filename  # reuse existing cropped file

                with fits.open(original_path, memmap=False) as hdul:
                    data = None
                    header = None
                    for hdu in hdul:
                        if hdu.data is not None:
                            data = hdu.data
                            header = hdu.header
                            break
                    if data is None:
                        raise RuntimeError(f"No image data found in FITS file: {original_path.name}")
                    header = header.copy() if header is not None else hdul[0].header.copy()

                cropped_data = data[self.crop_y0:self.crop_y1, self.crop_x0:self.crop_x1]

                header['CROPPED'] = (True, 'Image has been cropped')
                header['CROP_X0'] = (self.crop_x0, 'Crop region start X')
                header['CROP_Y0'] = (self.crop_y0, 'Crop region start Y')
                header['CROP_X1'] = (self.crop_x1, 'Crop region end X')
                header['CROP_Y1'] = (self.crop_y1, 'Crop region end Y')

                fits.PrimaryHDU(data=cropped_data, header=header).writeto(
                    cropped_path, overwrite=True
                )
                return filename

            max_workers = get_parallel_workers(self.params)
            max_workers = min(max_workers, len(files)) if files else 1

            errors = []
            if max_workers <= 1 or len(files) <= 1:
                for i, filename in enumerate(files, 1):
                    if progress.wasCanceled():
                        break
                    progress.setValue(i)
                    progress.setLabelText(f"Cropping {filename}...")
                    QApplication.processEvents()
                    try:
                        crop_one(filename)
                    except Exception as e:
                        errors.append((filename, str(e)))
            else:
                with ThreadPoolExecutor(max_workers=max_workers) as executor:
                    futures = {executor.submit(crop_one, fn): fn for fn in files}
                    for i, future in enumerate(as_completed(futures), 1):
                        if progress.wasCanceled():
                            for f in futures:
                                f.cancel()
                            break
                        filename = futures[future]
                        try:
                            future.result()
                        except Exception as e:
                            errors.append((filename, str(e)))
                        progress.setValue(i)
                        progress.setLabelText(f"Cropping {filename}...")
                        QApplication.processEvents()

            progress.setValue(len(files))

            if errors:
                self.is_cropped = False
                first_err = errors[0][1]
                QMessageBox.warning(
                    self, "Crop Warning",
                    f"Finished with {len(errors)} failures.\nFirst error: {first_err}"
                )
                self.update_navigation_buttons()
                return

            # Update displayed image to show cropped result
            self.displayed_image_data = self.original_image_data[self.crop_y0:self.crop_y1, self.crop_x0:self.crop_x1]

            # Redisplay cropped image (without rectangle)
            self.display_image()

            self.is_cropped = self._cropped_outputs_ready(require_marker=False)
            if not self.is_cropped:
                QMessageBox.warning(
                    self, "Crop Warning",
                    "Crop outputs are incomplete. Re-apply crop or skip crop before proceeding."
                )
                self.update_navigation_buttons()
                return

            # Save crop info only after all cropped outputs are confirmed.
            self.save_crop_info()

            # Update status
            width = self.crop_x1 - self.crop_x0
            height = self.crop_y1 - self.crop_y0
            self.crop_info_label.setText(
                f"Crop Applied!\n"
                f"Cropped files saved to: result/step2_crop/cropped/\n"
                f"Size: {width} × {height} pixels"
            )

            # Disable crop button
            self.btn_crop.setEnabled(False)

            self.update_navigation_buttons()
            self.go_next()

        except Exception as e:
            QMessageBox.critical(
                self, "Crop Error",
                f"Failed to crop images:\n{str(e)}"
            )

    def _is_crop_cache_valid(self, original_path: Path, cropped_path: Path) -> bool:
        """Return True if existing cropped file can be reused for current crop rectangle."""
        if not cropped_path.exists():
            return False
        try:
            src_stat = original_path.stat()
            dst_stat = cropped_path.stat()
            if int(dst_stat.st_mtime_ns) < int(src_stat.st_mtime_ns):
                return False

            header = fits.getheader(cropped_path, ext=0)
            if int(header.get("CROP_X0", -1)) != int(self.crop_x0):
                return False
            if int(header.get("CROP_Y0", -1)) != int(self.crop_y0):
                return False
            if int(header.get("CROP_X1", -1)) != int(self.crop_x1):
                return False
            if int(header.get("CROP_Y1", -1)) != int(self.crop_y1):
                return False

            expected_w = int(self.crop_x1 - self.crop_x0)
            expected_h = int(self.crop_y1 - self.crop_y0)
            if int(header.get("NAXIS1", -1)) != expected_w:
                return False
            if int(header.get("NAXIS2", -1)) != expected_h:
                return False
            return True
        except Exception:
            return False

    def has_completed_step4_or_later(self) -> bool:
        """Return True if Step 4+ has already been completed."""
        completed_steps = self.project_state.state.get("completed_steps", [])
        return any(step_index >= 3 for step_index in completed_steps)

    def confirm_recrop_warning(self) -> bool:
        """Warn that re-cropping after Step 4 invalidates downstream coordinates."""
        if not self.has_completed_step4_or_later():
            return True

        reply = QMessageBox.question(
            self, "Re-crop Warning",
            "Step 4 (Source Detection) or later is already completed.\n"
            "Applying crop again can shift pixel coordinates and invalidate downstream results.\n\n"
            "After re-cropping, rerun Step 4 and all later steps.\n\n"
            "Continue anyway?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )
        return reply == QMessageBox.Yes

    def reset_to_original(self):
        """Reset selection and reload original image"""
        if self.is_cropped:
            # Just reload original
            try:
                file_path = self.params.get_file_path(self.ref_filename)
                hdul = fits.open(file_path)
                data = None
                for hdu in hdul:
                    if hdu.data is not None:
                        data = hdu.data
                        break
                if data is None:
                    raise RuntimeError("No image data found in FITS file")
                self.original_image_data = data.astype(float).copy()
                hdul.close()

                self.displayed_image_data = self.original_image_data.copy()
                self.is_cropped = False

                # Clear crop coordinates
                self.crop_x0 = None
                self.crop_y0 = None
                self.crop_x1 = None
                self.crop_y1 = None

                # Redisplay
                self.display_image()

                self.crop_info_label.setText(f"Loaded: {self.ref_filename}\nReady to draw crop region")
                self.btn_crop.setEnabled(False)
                self.update_navigation_buttons()
            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to reload:\n{str(e)}")
            return

        # Just reset selection if not cropped
        # Reset crop coordinates
        self.crop_x0 = None
        self.crop_y0 = None
        self.crop_x1 = None
        self.crop_y1 = None

        # Disable crop button
        self.btn_crop.setEnabled(False)

        # Reset display to original
        if self.original_image_data is not None:
            self.displayed_image_data = self.original_image_data.copy()
            self.display_image()

        self.crop_info_label.setText(f"Loaded: {self.ref_filename}\nReady to draw crop region")
        self.update_navigation_buttons()

    def skip_crop(self):
        """Mark crop step as skipped (use original images)."""
        reply = QMessageBox.question(
            self, "Skip Crop",
            "Skip cropping and use original images?\n"
            "Any existing cropped outputs will be ignored.",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.No
        )
        if reply != QMessageBox.Yes:
            return

        self.crop_skipped = True
        self.is_cropped = False
        self.crop_x0 = None
        self.crop_y0 = None
        self.crop_x1 = None
        self.crop_y1 = None
        if self._viewer is not None:
            self._viewer.clear_selection_rect()

        if self.original_image_data is not None:
            self.displayed_image_data = self.original_image_data.copy()
            self.display_image()

        self.btn_crop.setEnabled(False)
        self.save_skip_info()
        self.crop_info_label.setText("Crop skipped. Original images will be used.")
        self.update_navigation_buttons()
        self.go_next()

    def validate_crop_region(self) -> bool:
        """Validate if crop region is valid"""
        if self.crop_x0 is None or self.crop_y0 is None or self.crop_x1 is None or self.crop_y1 is None:
            QMessageBox.warning(
                self, "No Crop Region",
                "Please draw a crop region first by clicking and dragging on the image."
            )
            return False
        return True

    def validate_step(self) -> bool:
        """Validate if step can be completed"""
        return self.crop_skipped or self._cropped_outputs_ready()

    def _cropped_outputs_ready(self, require_marker: bool = True) -> bool:
        """Return True only when the saved crop outputs are actually usable."""
        marker = crop_rect_path(self.params.P.result_dir)
        cropped_dir = step2_cropped_dir(self.params.P.result_dir)
        if require_marker and not marker.exists():
            return False
        if not cropped_dir.exists():
            return False

        files = list(self.file_manager.filenames) if self.file_manager else []
        if files:
            return all((cropped_dir / fname).exists() for fname in files)
        return any(cropped_dir.glob("*.fit*"))

    def save_crop_info(self):
        """Save crop info to JSON"""
        crop_data = {
            "x0": int(self.crop_x0),
            "y0": int(self.crop_y0),
            "x1": int(self.crop_x1),
            "y1": int(self.crop_y1),
            "ref_file": self.ref_filename or "",
            "skipped": False
        }

        # Save to result directory
        step2_out = step2_dir(self.params.P.result_dir)
        step2_out.mkdir(parents=True, exist_ok=True)
        crop_rect_path = step2_out / "crop_rect.json"
        with open(crop_rect_path, 'w') as f:
            json.dump(crop_data, f, indent=2)

        # Also save to project state
        self.project_state.store_step_data("crop", crop_data)

        print(f"Crop region saved: {crop_rect_path}")

    def save_skip_info(self):
        """Persist skip state and disable crop marker."""
        crop_data = {
            "skipped": True,
            "ref_file": self.ref_filename or ""
        }
        self.project_state.store_step_data("crop", crop_data)
        crop_rect_path = step2_dir(self.params.P.result_dir) / "crop_rect.json"
        try:
            crop_rect_path.unlink()
        except FileNotFoundError:
            pass
        except Exception:
            pass

    def save_state(self):
        """Save step state to project"""
        if self._cropped_outputs_ready():
            self.is_cropped = True
            self.crop_skipped = False
            self.save_crop_info()
        elif self.crop_skipped:
            self.save_skip_info()

    def restore_state(self):
        """Restore step state from project"""
        # Try to restore crop state
        state_data = self.project_state.get_step_data("crop")

        if state_data:
            if state_data.get("skipped"):
                self.crop_skipped = True
                self.is_cropped = False
                if state_data.get("ref_file"):
                    self.ref_filename = state_data.get("ref_file")
                try:
                    self.load_reference_image()
                except Exception:
                    pass
                self.crop_info_label.setText("Crop skipped. Original images will be used.")
                self.update_navigation_buttons()
                return

            self.crop_x0 = state_data.get("x0")
            self.crop_y0 = state_data.get("y0")
            self.crop_x1 = state_data.get("x1")
            self.crop_y1 = state_data.get("y1")

            # Use crop's ref_file if available
            if state_data.get("ref_file"):
                self.ref_filename = state_data.get("ref_file")

            self.crop_skipped = False
            self.is_cropped = self._cropped_outputs_ready()

            # Try to load image
            try:
                self.load_reference_image()

                # Update info label
                if self.is_cropped and self.crop_x0 is not None:
                    width = self.crop_x1 - self.crop_x0
                    height = self.crop_y1 - self.crop_y0
                    self.crop_info_label.setText(
                        f"Crop Previously Applied\n"
                        f"Cropped files in: result/step2_crop/cropped/\n"
                        f"Size: {width} × {height} pixels"
                    )
                elif self.crop_x0 is not None:
                    self.crop_info_label.setText(
                        "Saved crop region found, but cropped outputs are missing.\n"
                        "Re-apply crop or skip crop to continue."
                    )
            except:
                pass

            self.update_navigation_buttons()
        else:
            # No saved state — auto-load first frame
            try:
                self.load_reference_image()
            except Exception:
                pass
