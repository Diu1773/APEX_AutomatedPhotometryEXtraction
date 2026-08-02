"""
IRAF Photometry Tool Window
Comprehensive IRAF/DAOPHOT photometry interface with all parameters.
Uses WSL subprocess to run PyRAF on Windows.
"""

from __future__ import annotations

import os
import sys
import subprocess
import tempfile
import json
import time
import shutil
try:
    import tomllib
except ImportError:
    import tomli as tomllib  # type: ignore
try:
    import tomli_w  # type: ignore
except Exception:
    tomli_w = None
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from apex.utils.constants import MAD_TO_SIGMA

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QColor, QKeySequence
from PyQt5.QtWidgets import (
    QWidget, QMainWindow, QVBoxLayout, QHBoxLayout, QPushButton, QLabel,
    QTextEdit, QGroupBox, QFormLayout, QLineEdit, QDoubleSpinBox, QSpinBox,
    QFileDialog, QMessageBox, QProgressBar, QTabWidget, QComboBox,
    QCheckBox, QSplitter, QTableWidget, QTableWidgetItem, QScrollArea,
    QFrame, QSlider, QShortcut, QDialog, QHeaderView, QGridLayout
)

from apex.gui.layout_rules import AutoFitMixin, FittedDialog, clamp_to_screen, tame_canvas
from apex.gui.widgets.fits_viewer import FITSViewerWidget, OverlayMarker
from apex.gui.workflow.run_control import RunControlBar
from apex.gui.theme import Tokens, mono_note_style
from apex.gui.workflow.ui_helpers import create_cache_action_button, create_parameter_button
from apex.utils.common_helpers import normalize_filter_key
from apex.utils.step_paths import step4_dir, step7_forced_phot_dir
from apex.utils.step_paths_lc import step2_cropped_dir


_SUBPROCESS_TEXT_KWARGS: dict = {
    "text": True,
    "encoding": "utf-8",
    "errors": "replace",
}
if sys.platform == "win32":
    _SUBPROCESS_TEXT_KWARGS["creationflags"] = subprocess.CREATE_NO_WINDOW


def windows_to_wsl_path(win_path: str) -> str:
    """Convert Windows path to WSL path."""
    path = str(win_path).replace("\\", "/")
    if len(path) >= 2 and path[1] == ":":
        drive = path[0].lower()
        path = f"/mnt/{drive}{path[2:]}"
    return path


def wsl_to_windows_path(wsl_path: str) -> str:
    """Convert WSL path to Windows path."""
    if wsl_path.startswith("/mnt/") and len(wsl_path) > 6:
        drive = wsl_path[5].upper()
        return f"{drive}:{wsl_path[6:]}".replace("/", "\\")
    return wsl_path


def _format_subprocess_output(text: str | None, max_chars: int = 4000) -> str:
    out = str(text or "").strip()
    if len(out) > max_chars:
        out = out[-max_chars:]
    return out


def _iraf_runtime_cmd(is_windows: bool) -> list[str]:
    return ["wsl", "python3"] if is_windows else ["python3"]


def _unbuffered_python_cmd(runtime_cmd: list[str]) -> list[str]:
    return list(runtime_cmd) + ["-u"]


def _parse_probe_statuses(output: str) -> dict[str, tuple[str, str]]:
    statuses: dict[str, tuple[str, str]] = {}
    for raw in str(output or "").splitlines():
        line = raw.strip()
        if line.startswith("[CHECK-OK] "):
            state = "ok"
            payload = line[len("[CHECK-OK] "):]
        elif line.startswith("[CHECK-FAIL] "):
            state = "fail"
            payload = line[len("[CHECK-FAIL] "):]
        elif line.startswith("[CHECK-SKIP] "):
            state = "warn"
            payload = line[len("[CHECK-SKIP] "):]
        else:
            continue
        key, _, message = payload.partition("|")
        statuses[key.strip()] = (state, message.strip() or state.upper())
    return statuses


def _run_iraf_environment_probe(
    runtime_cmd: list[str],
    data_dir: str,
    output_dir: str,
    log_callback=None,
    timeout: int = 30,
) -> dict:
    code = f"""
import os
import sys

failed = False
data_dir = {data_dir!r}
output_dir = {output_dir!r}
print(f"[CHECK-OK] runtime|Python: {{sys.executable}}")
if not os.path.isdir(data_dir):
    print(f"[CHECK-FAIL] data|Data directory is not visible: {{data_dir}}")
    failed = True
else:
    print(f"[CHECK-OK] data|Data directory visible")
try:
    os.makedirs(output_dir, exist_ok=True)
    print(f"[CHECK-OK] output|Output directory writable")
except Exception as exc:
    print(f"[CHECK-FAIL] output|{{type(exc).__name__}}: {{exc}}")
    failed = True
try:
    from pyraf import iraf
    print("[CHECK-OK] pyraf|PyRAF import OK")
except Exception as exc:
    iraf = None
    print(f"[CHECK-FAIL] pyraf|{{type(exc).__name__}}: {{exc}}")
    failed = True
if iraf is None:
    print("[CHECK-SKIP] daophot|Skipped because PyRAF import failed")
else:
    try:
        iraf.noao(); iraf.digiphot(); iraf.daophot()
        print("[CHECK-OK] daophot|IRAF DAOPHOT ready")
    except Exception as exc:
        print(f"[CHECK-FAIL] daophot|{{type(exc).__name__}}: {{exc}}")
        failed = True
sys.exit(1 if failed else 0)
"""
    proc = subprocess.run(
        runtime_cmd + ["-c", code],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=timeout,
        **_SUBPROCESS_TEXT_KWARGS,
    )
    out = _format_subprocess_output(proc.stdout)
    if out and log_callback:
        for line in out.splitlines():
            log_callback(line)
    return {
        "ok": proc.returncode == 0,
        "returncode": proc.returncode,
        "output": out,
        "statuses": _parse_probe_statuses(out),
        "command": runtime_cmd + ["-c", "<probe>"],
    }


# ============================================================================
# IRAF Parameter Dataclass-like storage
# ============================================================================
class IRAFParameters:
    """Storage for all IRAF DAOPHOT parameters."""

    def __init__(self):
        # === DATAPARS ===
        self.scale = 1.0           # Image scale in units/pixel
        self.emission = True       # Emission features (True) or absorption (False)
        self.datamax = 60000.0     # Maximum good data value
        self.noise = "poisson"     # Noise model: poisson, constant, file
        self.readnoise = 2.5       # CCD readout noise in electrons
        self.epadu = 0.689         # Gain in electrons per ADU
        self.exposure = "EXPTIME"  # Exposure time header keyword
        self.itime = 1.0           # Integration time (if not in header)

        # === Filter-specific FWHM (seeing in arcsec) ===
        self.seeing_g = 2.5
        self.seeing_r = 2.5
        self.seeing_i = 2.5
        self.seeing_default = 2.5

        # === Filter-specific sigma (sky background std) ===
        self.sigma_g = 50.0
        self.sigma_r = 50.0
        self.sigma_i = 50.0
        self.sigma_default = 50.0

        # === FINDPARS ===
        # Filter-specific thresholds (g, r, i)
        self.threshold_g = 4.0
        self.threshold_r = 4.5
        self.threshold_i = 5.0
        self.threshold_default = 5.0  # For unknown filters
        self.nsigma = 1.5          # Width of convolution kernel in sigma
        self.ratio = 1.0           # Ratio of minor to major axis of Gaussian
        self.theta = 0.0           # Position angle of major axis
        # Filter-specific sharplo
        self.sharplo_g = 0.2
        self.sharplo_r = 0.2
        self.sharplo_i = 0.4
        self.sharplo_default = 0.2
        self.sharphi = 1.0         # Upper bound on sharpness
        self.roundlo = -1.0        # Lower bound on roundness
        self.roundhi = 1.0         # Upper bound on roundness
        # Filter-specific datamin
        self.datamin_g = -100.0
        self.datamin_r = -100.0
        self.datamin_i = 0.0
        self.datamin_default = -100.0

        # === CENTERPARS ===
        self.calgorithm = "centroid"  # Centering algorithm
        self.cbox_mult = 2.0       # Centering box = FWHM * cbox_mult
        self.cthreshold = 0.0      # Centering threshold in sigma
        self.minsnratio = 1.0      # Minimum signal-to-noise ratio
        self.cmaxiter = 10         # Maximum iterations for centering
        self.maxshift = 1.0        # Maximum shift in scale units
        self.clean = False         # Symmetry clean before centering
        self.rclean = 1.0          # Cleaning radius in scale units
        self.rclip = 2.0           # Clipping radius in scale units
        self.kclean = 3.0          # K-sigma rejection criterion

        # === FITSKYPARS (FWHM multipliers) ===
        self.salgorithm = "mode"   # Sky algorithm
        self.annulus_mult = 4.0    # Inner radius = FWHM * annulus_mult
        self.dannulus_mult = 2.0   # Width = FWHM * dannulus_mult
        self.skyvalue = 0.0        # User sky value (for constant algorithm)
        self.smaxiter = 10         # Maximum iterations for sky fitting
        self.sloclip = 0.0         # Lower clipping factor (sigma)
        self.shiclip = 0.0         # Upper clipping factor (sigma)
        self.snreject = 50         # Maximum number of rejection iterations
        self.sloreject = 3.0       # Lower K-sigma rejection limit
        self.shireject = 3.0       # Upper K-sigma rejection limit
        self.khist = 3.0           # Half-width of histogram in sigma
        self.binsize = 0.1         # Binsize of histogram in sigma
        self.smooth = False        # Smooth histogram before fitting
        self.rgrow = 0.0           # Region growing radius in scale units

        # === PHOTPARS (FWHM multiplier) ===
        self.aperture_mult = 1.0   # Aperture = FWHM * aperture_mult
        self.zmag = 25.0           # Zero point of magnitude scale
        self.mkapert = False       # Make aperture plots

        # === Convenience settings ===
        self.pix_scale = 0.392     # Pixel scale in arcsec/pixel
        self.sigma_ref = 50.0      # Reference sigma for threshold scaling

    def to_dict(self):
        """Convert to dictionary for script generation."""
        return {k: v for k, v in self.__dict__.items()}

    def from_dict(self, d: dict):
        """Load parameters from dictionary."""
        for k, v in d.items():
            if hasattr(self, k):
                setattr(self, k, v)


# ============================================================================
# PyRAF Script Template
# ============================================================================
PYRAF_SCRIPT_TEMPLATE = '''#!/usr/bin/env python3
"""Auto-generated PyRAF photometry script with filter-specific parameters"""

import os, glob, sys, time, shutil, atexit, signal
import numpy as np

try:
    sys.stdout.reconfigure(line_buffering=True)
except Exception:
    pass

def log(message):
    print(message, flush=True)

log("[INIT] Importing PyRAF...")
from pyraf import iraf
log("[INIT] Loading IRAF packages...")
iraf.noao(); iraf.digiphot(); iraf.daophot()

# =====================
# Paths
# =====================
DATA_DIR = "{data_dir}"
OUTDIR = "{output_dir}"
FILE_PATTERN = "{file_pattern}"
SKIP_EXISTING = {skip_existing}
USE_WSL_SCRATCH = {use_wsl_scratch}
SCRATCH_ROOT = "{scratch_root}"
SCRATCH_RUN_DIR_REQUESTED = "{scratch_run_dir}"

# =====================
# Filter config (from TOML)
# =====================
FILTER_PARAMS = {filter_params}
FILTER_ALIASES = {filter_aliases}
BUILTIN_FILTER_ALIASES = {{
    "b": "B", "v": "V",
    "rc": "R", "r_c": "R", "r_cousins": "R", "rj": "R", "r_j": "R",
    "ic": "I", "i_c": "I", "i_cousins": "I", "ij": "I", "i_j": "I",
    "bj": "B", "b_j": "B", "b_john": "B",
    "vj": "V", "v_j": "V", "v_john": "V",
    "sdss-g": "g", "sdss_g": "g", "g'": "g", "gprime": "g",
    "sdss-r": "r", "sdss_r": "r", "r'": "r", "rprime": "r",
    "sdss-i": "i", "sdss_i": "i", "i'": "i", "iprime": "i",
}}
PARAM_DEFAULTS = {param_defaults}
APEX_FRAME_PARAMS = {apex_frame_params}
USE_APEX_FRAME_PARAMS = {use_apex_frame_params}

# =====================
# Pixel scale
# =====================
PIX_SCALE = {pix_scale}  # arcsec/pixel

# =====================
# Filter-specific SEEING (arcsec)
# =====================
SEEING = {{"g": {seeing_g}, "r": {seeing_r}, "i": {seeing_i}, "default": {seeing_default}}}

# =====================
# Filter-specific SIGMA (sky background std)
# =====================
SIGMA = {{"g": {sigma_g}, "r": {sigma_r}, "i": {sigma_i}, "default": {sigma_default}}}
SIGMA_REF = {sigma_ref}

# =====================
# DATAPARS (base values)
# =====================
DATAPARS = {{
    "scale": {scale},
    "emission": {emission},
    "datamax": {datamax},
    "noise": "{noise}",
    "readnoise": {readnoise},
    "epadu": {epadu},
    "exposure": "{exposure}",
    "itime": {itime},
}}

# =====================
# FINDPARS - Filter-specific
# =====================
THRESHOLD = {{"g": {threshold_g}, "r": {threshold_r}, "i": {threshold_i}, "default": {threshold_default}}}
SHARPLO = {{"g": {sharplo_g}, "r": {sharplo_r}, "i": {sharplo_i}, "default": {sharplo_default}}}
DATAMIN = {{"g": {datamin_g}, "r": {datamin_r}, "i": {datamin_i}, "default": {datamin_default}}}

FINDPARS_BASE = {{
    "nsigma": {nsigma},
    "ratio": {ratio},
    "theta": {theta},
    "sharphi": {sharphi},
    "roundlo": {roundlo},
    "roundhi": {roundhi},
}}

# =====================
# CENTERPARS
# =====================
CENTERPARS = {{
    "calgorithm": "{calgorithm}",
    "cbox_mult": {cbox_mult},
    "cthreshold": {cthreshold},
    "minsnratio": {minsnratio},
    "cmaxiter": {cmaxiter},
    "maxshift": {maxshift},
    "clean": {clean},
    "rclean": {rclean},
    "rclip": {rclip},
    "kclean": {kclean},
}}

# =====================
# FITSKYPARS (FWHM multipliers)
# =====================
FITSKYPARS = {{
    "salgorithm": "{salgorithm}",
    "annulus_mult": {annulus_mult},
    "dannulus_mult": {dannulus_mult},
    "skyvalue": {skyvalue},
    "smaxiter": {smaxiter},
    "sloclip": {sloclip},
    "shiclip": {shiclip},
    "snreject": {snreject},
    "sloreject": {sloreject},
    "shireject": {shireject},
    "khist": {khist},
    "binsize": {binsize},
    "smooth": {smooth},
    "rgrow": {rgrow},
}}

# =====================
# PHOTPARS (FWHM multiplier)
# =====================
PHOTPARS = {{
    "aperture_mult": {aperture_mult},
    "zmag": {zmag},
    "mkapert": {mkapert},
}}

# =====================
# Helper functions
# =====================
def safe_rm(p):
    try: os.remove(p)
    except FileNotFoundError: pass

def get_header(im, key):
    try:
        out = iraf.hselect(im, key, "yes", Stdout=1)
        v = out[0].strip() if out else ""
        if v in ["", "INDEF", "indef"]:
            return None
        return v
    except:
        return None

def normalize_filter(val):
    raw = str(val).strip()
    key = raw.lower()
    if key in FILTER_ALIASES:
        return str(FILTER_ALIASES[key]).strip()
    if key in BUILTIN_FILTER_ALIASES:
        return str(BUILTIN_FILTER_ALIASES[key]).strip()
    return raw

def filename_filter_tokens(im):
    stem = os.path.basename(im).lower()
    for ext in [".fits", ".fit", ".fts", ".fz", ".gz"]:
        if stem.endswith(ext):
            stem = stem[:-len(ext)]
    tokens = []
    buf = []
    for ch in stem:
        if ch.isalnum():
            buf.append(ch)
        elif buf:
            tokens.append("".join(buf))
            buf = []
    if buf:
        tokens.append("".join(buf))
    return tokens

def guess_filter(im):
    for k in ["FILTER", "FILTER1", "FILTER2", "FILTNAM"]:
        v = get_header(im, k)
        if v:
            return normalize_filter(v)
    low = os.path.basename(im).lower()
    tokens = filename_filter_tokens(im)
    for alias, canon in FILTER_ALIASES.items():
        alias_key = str(alias).strip().lower()
        if alias_key and (alias_key in tokens or (len(alias_key) > 1 and alias_key in low)):
            return str(canon).strip()
    for key in FILTER_PARAMS.keys():
        key_low = str(key).strip().lower()
        if key != "default" and key_low and (key_low in tokens or (len(key_low) > 1 and key_low in low)):
            return key
    if "b" in tokens: return "B"
    if "v" in tokens: return "V"
    if "g" in tokens: return "g"
    if "r" in tokens: return "r"
    if "i" in tokens: return "i"
    return "unknown"

def get_param(name, band, default=None):
    b = normalize_filter(band)
    if b in FILTER_PARAMS and name in FILTER_PARAMS[b]:
        return FILTER_PARAMS[b][name]
    if "default" in FILTER_PARAMS and name in FILTER_PARAMS["default"]:
        return FILTER_PARAMS["default"][name]
    if name in PARAM_DEFAULTS:
        return PARAM_DEFAULTS[name]
    return default

def estimate_sigma(image):
    try:
        out = iraf.imstat(image, fields="stddev", format="no",
                         nclip=5, lsigma=3., usigma=3., Stdout=1)
        return float(out[0].strip())
    except:
        return None

def clean_frame_key(name):
    stem = os.path.basename(str(name))
    low = stem.lower()
    for ext in [".fits.gz", ".fit.gz", ".fits", ".fit", ".fts", ".fz", ".gz"]:
        if low.endswith(ext):
            return stem[:-len(ext)]
    return os.path.splitext(stem)[0]

def frame_params_for(im):
    if not USE_APEX_FRAME_PARAMS or not APEX_FRAME_PARAMS:
        return {{}}
    fname = os.path.basename(im)
    keys = [fname, os.path.splitext(fname)[0], clean_frame_key(fname)]
    for key in keys:
        if key in APEX_FRAME_PARAMS:
            return APEX_FRAME_PARAMS[key]
    return {{}}

def finite_positive(value):
    try:
        out = float(value)
        if np.isfinite(out) and out > 0:
            return out
    except:
        pass
    return None

def fwhm_pix_from_header_or_default(im, band):
    # Try arcsec keys from header
    for k in ["FWHMARC", "SEEING", "FWHM_AS"]:
        v = get_header(im, k)
        if v:
            try:
                val = float(v)
                if val > 0:
                    return val / PIX_SCALE
            except:
                pass
    # Try pixel keys from header
    for k in ["FWHMPSF", "FWHM_PIX"]:
        v = get_header(im, k)
        if v:
            try:
                val = float(v)
                if val > 0:
                    return val
            except:
                pass
    # Use filter-specific seeing
    seeing_arcsec = get_param("seeing", band, SEEING.get(band, SEEING["default"]))
    return seeing_arcsec / PIX_SCALE

def format_bytes(nbytes):
    try:
        value = float(nbytes)
    except:
        value = 0.0
    units = ["B", "KB", "MB", "GB", "TB"]
    for unit in units:
        if value < 1024.0 or unit == units[-1]:
            return f"{{value:.1f}} {{unit}}" if unit != "B" else f"{{int(value)}} B"
        value /= 1024.0

SCRATCH_RUN_DIR = None

def cleanup_scratch():
    try:
        if SCRATCH_RUN_DIR and os.path.isdir(SCRATCH_RUN_DIR):
            shutil.rmtree(SCRATCH_RUN_DIR, ignore_errors=True)
            log(f"[SCRATCH] Cleaned scratch run dir: {{SCRATCH_RUN_DIR}}")
    except Exception as exc:
        log(f"[SCRATCH] Cleanup warning: {{exc}}")

atexit.register(cleanup_scratch)

def handle_stop_signal(signum, _frame):
    log(f"[STOP] Received signal {{signum}}; cleaning scratch and exiting.")
    cleanup_scratch()
    sys.exit(128 + int(signum))

try:
    signal.signal(signal.SIGTERM, handle_stop_signal)
    signal.signal(signal.SIGINT, handle_stop_signal)
except Exception:
    pass

# =====================
# Unlearn all tasks
# =====================
for t in ["daofind", "phot", "datapars", "findpars", "centerpars", "fitskypars", "photpars"]:
    try: iraf.unlearn(t)
    except: pass

# =====================
# Apply base DATAPARS
# =====================
iraf.datapars.scale = DATAPARS["scale"]
iraf.datapars.emission = "yes" if DATAPARS["emission"] else "no"
iraf.datapars.datamax = DATAPARS["datamax"]
iraf.datapars.noise = DATAPARS["noise"]
iraf.datapars.readnoise = DATAPARS["readnoise"]
iraf.datapars.epadu = DATAPARS["epadu"]
iraf.datapars.exposure = DATAPARS["exposure"]
iraf.datapars.itime = DATAPARS["itime"]

# =====================
# Apply base FINDPARS
# =====================
iraf.findpars.nsigma = FINDPARS_BASE["nsigma"]
iraf.findpars.ratio = FINDPARS_BASE["ratio"]
iraf.findpars.theta = FINDPARS_BASE["theta"]
iraf.findpars.sharphi = FINDPARS_BASE["sharphi"]
iraf.findpars.roundlo = FINDPARS_BASE["roundlo"]
iraf.findpars.roundhi = FINDPARS_BASE["roundhi"]

# =====================
# Apply CENTERPARS (base)
# =====================
iraf.centerpars.calgorithm = CENTERPARS["calgorithm"]
iraf.centerpars.cthreshold = CENTERPARS["cthreshold"]
iraf.centerpars.minsnratio = CENTERPARS["minsnratio"]
iraf.centerpars.cmaxiter = CENTERPARS["cmaxiter"]
iraf.centerpars.maxshift = CENTERPARS["maxshift"]
iraf.centerpars.clean = "yes" if CENTERPARS["clean"] else "no"
iraf.centerpars.rclean = CENTERPARS["rclean"]
iraf.centerpars.rclip = CENTERPARS["rclip"]
iraf.centerpars.kclean = CENTERPARS["kclean"]

# =====================
# Apply FITSKYPARS (base)
# =====================
iraf.fitskypars.salgorithm = FITSKYPARS["salgorithm"]
iraf.fitskypars.skyvalue = FITSKYPARS["skyvalue"]
iraf.fitskypars.smaxiter = FITSKYPARS["smaxiter"]
iraf.fitskypars.sloclip = FITSKYPARS["sloclip"]
iraf.fitskypars.shiclip = FITSKYPARS["shiclip"]
iraf.fitskypars.snreject = FITSKYPARS["snreject"]
iraf.fitskypars.sloreject = FITSKYPARS["sloreject"]
iraf.fitskypars.shireject = FITSKYPARS["shireject"]
iraf.fitskypars.khist = FITSKYPARS["khist"]
iraf.fitskypars.binsize = FITSKYPARS["binsize"]
iraf.fitskypars.smooth = "yes" if FITSKYPARS["smooth"] else "no"
iraf.fitskypars.rgrow = FITSKYPARS["rgrow"]

# =====================
# Apply PHOTPARS (base)
# =====================
iraf.photpars.zmag = PHOTPARS["zmag"]
iraf.photpars.mkapert = "yes" if PHOTPARS["mkapert"] else "no"

log("[INIT] PyRAF initialized")
log(f"[PARAM] pix_scale={{PIX_SCALE:.3f}} arcsec/pix")
log(f"[PARAM] epadu={{DATAPARS['epadu']}}, rdnoise={{DATAPARS['readnoise']}}")
log(f"[PARAM] FWHM multipliers: ap={{PHOTPARS['aperture_mult']}}x, ann={{FITSKYPARS['annulus_mult']}}x, dann={{FITSKYPARS['dannulus_mult']}}x")
log(f"[PARAM] Filter-specific seeing: g={{SEEING['g']}}, r={{SEEING['r']}}, i={{SEEING['i']}} arcsec")
log(f"[PARAM] Filter-specific sigma: g={{SIGMA['g']}}, r={{SIGMA['r']}}, i={{SIGMA['i']}}")

FINAL_OUTDIR = OUTDIR
WORK_DATA_DIR = DATA_DIR
WORK_OUTDIR = OUTDIR
SCRATCH_SOURCE_BY_NAME = {{}}

os.makedirs(FINAL_OUTDIR, exist_ok=True)
os.chdir(DATA_DIR)
source_imgs = sorted(glob.glob(FILE_PATTERN))

if not source_imgs:
    log(f"[ERROR] No images found matching {{FILE_PATTERN}}")
    sys.exit(1)

pending_imgs = []
skipped = 0
for src in source_imgs:
    base = os.path.splitext(os.path.basename(src))[0]
    final_txt = os.path.join(FINAL_OUTDIR, f"{{base}}.txt")
    if SKIP_EXISTING and os.path.exists(final_txt):
        skipped += 1
        continue
    pending_imgs.append(src)

if not pending_imgs:
    log(f"[DONE] Processed: 0, Skipped: {{skipped}}, Output: {{FINAL_OUTDIR}}")
    sys.exit(0)

imgs = pending_imgs
if USE_WSL_SCRATCH:
    scratch_base = os.path.abspath(SCRATCH_ROOT or "/tmp/apex_iraf")
    os.makedirs(scratch_base, exist_ok=True)

    input_bytes = 0
    for src in pending_imgs:
        try:
            input_bytes += os.path.getsize(src)
        except:
            pass
    required_bytes = int(input_bytes * 2.0 + 512 * 1024 * 1024)
    usage = shutil.disk_usage(scratch_base)
    log(
        f"[SCRATCH] Enabled root={{scratch_base}} | input={{format_bytes(input_bytes)}} "
        f"| free={{format_bytes(usage.free)}} | required~={{format_bytes(required_bytes)}}"
    )
    if usage.free < required_bytes:
        log("[ERROR] Not enough free space in WSL scratch. Disable scratch or free WSL disk space.")
        sys.exit(2)

    if SCRATCH_RUN_DIR_REQUESTED:
        SCRATCH_RUN_DIR = os.path.abspath(SCRATCH_RUN_DIR_REQUESTED)
    else:
        SCRATCH_RUN_DIR = os.path.join(scratch_base, f"apex_iraf_{{int(time.time())}}_{{os.getpid()}}")
    WORK_DATA_DIR = os.path.join(SCRATCH_RUN_DIR, "data")
    WORK_OUTDIR = os.path.join(SCRATCH_RUN_DIR, "out")
    os.makedirs(WORK_DATA_DIR, exist_ok=True)
    os.makedirs(WORK_OUTDIR, exist_ok=True)
    log(f"[SCRATCH] Run dir: {{SCRATCH_RUN_DIR}}")

    SCRATCH_SOURCE_BY_NAME = {{
        os.path.basename(src): os.path.abspath(src)
        for src in pending_imgs
    }}
    imgs = list(SCRATCH_SOURCE_BY_NAME.keys())
    log(f"[SCRATCH] Copy mode: per-frame copy for {{len(imgs)}} pending FITS file(s)")
    log(f"[SCRATCH] IRAF work data: {{WORK_DATA_DIR}}")
    log(f"[SCRATCH] IRAF work output: {{WORK_OUTDIR}}")
    log(f"[SCRATCH] Final output: {{FINAL_OUTDIR}}")
    os.chdir(WORK_DATA_DIR)
else:
    os.chdir(WORK_DATA_DIR)

log(f"[INFO] Found {{len(source_imgs)}} images; pending {{len(imgs)}}; skipped existing {{skipped}}")

processed = 0

for idx, im in enumerate(imgs):
    base = os.path.splitext(os.path.basename(im))[0]
    final_coo = os.path.join(FINAL_OUTDIR, f"{{base}}.coo")
    final_mag = os.path.join(FINAL_OUTDIR, f"{{base}}.mag")
    final_txt = os.path.join(FINAL_OUTDIR, f"{{base}}.txt")

    # Skip if already processed
    if SKIP_EXISTING and os.path.exists(final_txt):
        log(f"[SKIP] {{idx+1}}/{{len(imgs)}} {{im}} (already exists)")
        skipped += 1
        continue

    copyin_seconds = 0.0
    work_image_path = None
    if USE_WSL_SCRATCH:
        src_path = SCRATCH_SOURCE_BY_NAME.get(os.path.basename(im))
        work_image_path = os.path.join(WORK_DATA_DIR, os.path.basename(im))
        if not src_path:
            log(f"  -> ERROR: scratch source missing for {{im}}")
            continue
        try:
            t_copyin = time.perf_counter()
            shutil.copy2(src_path, work_image_path)
            copyin_seconds = time.perf_counter() - t_copyin
            log(f"[SCRATCH] Copy-in {{idx+1}}/{{len(imgs)}} {{os.path.basename(im)}}: {{copyin_seconds:.1f}}s")
            im = os.path.basename(work_image_path)
        except Exception as exc:
            log(f"  -> ERROR: scratch copy-in failed for {{im}}: {{exc}}")
            continue

    band = guess_filter(im)
    frame_params = frame_params_for(im)

    # Get filter-specific FWHM
    apex_fwhm = finite_positive(frame_params.get("fwhm_px"))
    if apex_fwhm is not None:
        fwhm_pix = apex_fwhm
        fwhm_source = "apex-step4"
    else:
        fwhm_pix = fwhm_pix_from_header_or_default(im, band)
        fwhm_source = "header/default"
    iraf.datapars.fwhmpsf = float(fwhm_pix)

    # Get filter-specific sigma (or auto-estimate)
    apex_sigma = finite_positive(frame_params.get("sigma"))
    if apex_sigma is not None:
        sig = apex_sigma
        sigma_source = "apex-step4"
    elif {auto_sigma}:
        sig = estimate_sigma(im)
        if sig is None or sig <= 0:
            sig = get_param("sigma", band, SIGMA.get(band, SIGMA["default"]))
            sigma_source = "filter-default"
        else:
            sigma_source = "imstat"
    else:
        sig = get_param("sigma", band, SIGMA.get(band, SIGMA["default"]))
        sigma_source = "filter-default"
    iraf.datapars.sigma = float(sig)

    # Filter-specific threshold with sigma scaling
    base_thr = get_param("threshold", band, THRESHOLD.get(band, THRESHOLD["default"]))
    sigma_ref = float(get_param("sigma_ref", band, SIGMA_REF))
    thr = base_thr * np.clip(sig / sigma_ref, 0.8, 1.6)
    thr = float(np.clip(thr, 3.5, 15.0))
    iraf.findpars.threshold = thr

    # Filter-specific sharplo
    sharplo = get_param("sharplo", band, SHARPLO.get(band, SHARPLO["default"]))
    iraf.findpars.sharplo = float(sharplo)

    # Filter-specific datamin
    datamin = get_param("datamin", band, DATAMIN.get(band, DATAMIN["default"]))
    iraf.datapars.datamin = float(datamin)

    # Per-filter overrides
    iraf.datapars.scale = float(get_param("scale", band, DATAPARS["scale"]))
    iraf.datapars.emission = "yes" if get_param("emission", band, DATAPARS["emission"]) else "no"
    iraf.datapars.datamax = float(get_param("datamax", band, DATAPARS["datamax"]))
    iraf.datapars.noise = str(get_param("noise", band, DATAPARS["noise"]))
    iraf.datapars.readnoise = float(get_param("readnoise", band, DATAPARS["readnoise"]))
    iraf.datapars.epadu = float(get_param("epadu", band, DATAPARS["epadu"]))
    iraf.datapars.exposure = str(get_param("exposure", band, DATAPARS["exposure"]))
    iraf.datapars.itime = float(get_param("itime", band, DATAPARS["itime"]))

    iraf.findpars.nsigma = float(get_param("nsigma", band, FINDPARS_BASE["nsigma"]))
    iraf.findpars.ratio = float(get_param("ratio", band, FINDPARS_BASE["ratio"]))
    iraf.findpars.theta = float(get_param("theta", band, FINDPARS_BASE["theta"]))
    iraf.findpars.sharphi = float(get_param("sharphi", band, FINDPARS_BASE["sharphi"]))
    iraf.findpars.roundlo = float(get_param("roundlo", band, FINDPARS_BASE["roundlo"]))
    iraf.findpars.roundhi = float(get_param("roundhi", band, FINDPARS_BASE["roundhi"]))

    iraf.centerpars.calgorithm = str(get_param("calgorithm", band, CENTERPARS["calgorithm"]))
    iraf.centerpars.cthreshold = float(get_param("cthreshold", band, CENTERPARS["cthreshold"]))
    iraf.centerpars.minsnratio = float(get_param("minsnratio", band, CENTERPARS["minsnratio"]))
    iraf.centerpars.cmaxiter = int(get_param("cmaxiter", band, CENTERPARS["cmaxiter"]))
    iraf.centerpars.maxshift = float(get_param("maxshift", band, CENTERPARS["maxshift"]))
    iraf.centerpars.clean = "yes" if get_param("clean", band, CENTERPARS["clean"]) else "no"
    iraf.centerpars.rclean = float(get_param("rclean", band, CENTERPARS["rclean"]))
    iraf.centerpars.rclip = float(get_param("rclip", band, CENTERPARS["rclip"]))
    iraf.centerpars.kclean = float(get_param("kclean", band, CENTERPARS["kclean"]))

    iraf.fitskypars.salgorithm = str(get_param("salgorithm", band, FITSKYPARS["salgorithm"]))
    iraf.fitskypars.skyvalue = float(get_param("skyvalue", band, FITSKYPARS["skyvalue"]))
    iraf.fitskypars.smaxiter = int(get_param("smaxiter", band, FITSKYPARS["smaxiter"]))
    iraf.fitskypars.sloclip = float(get_param("sloclip", band, FITSKYPARS["sloclip"]))
    iraf.fitskypars.shiclip = float(get_param("shiclip", band, FITSKYPARS["shiclip"]))
    iraf.fitskypars.snreject = int(get_param("snreject", band, FITSKYPARS["snreject"]))
    iraf.fitskypars.sloreject = float(get_param("sloreject", band, FITSKYPARS["sloreject"]))
    iraf.fitskypars.shireject = float(get_param("shireject", band, FITSKYPARS["shireject"]))
    iraf.fitskypars.khist = float(get_param("khist", band, FITSKYPARS["khist"]))
    iraf.fitskypars.binsize = float(get_param("binsize", band, FITSKYPARS["binsize"]))
    iraf.fitskypars.smooth = "yes" if get_param("smooth", band, FITSKYPARS["smooth"]) else "no"
    iraf.fitskypars.rgrow = float(get_param("rgrow", band, FITSKYPARS["rgrow"]))

    # FWHM-based parameters
    cbox = fwhm_pix * float(get_param("cbox_mult", band, CENTERPARS["cbox_mult"]))
    iraf.centerpars.cbox = float(cbox)

    aperture = fwhm_pix * float(get_param("aperture_mult", band, PHOTPARS["aperture_mult"]))
    iraf.photpars.apertures = f"{{aperture:.2f}}"

    annulus = fwhm_pix * float(get_param("annulus_mult", band, FITSKYPARS["annulus_mult"]))
    dannulus = fwhm_pix * float(get_param("dannulus_mult", band, FITSKYPARS["dannulus_mult"]))
    iraf.fitskypars.annulus = float(annulus)
    iraf.fitskypars.dannulus = float(dannulus)

    iraf.photpars.zmag = float(get_param("zmag", band, PHOTPARS["zmag"]))
    iraf.photpars.mkapert = "yes" if get_param("mkapert", band, PHOTPARS["mkapert"]) else "no"

    coo = os.path.join(WORK_OUTDIR, f"{{base}}.coo")
    mag = os.path.join(WORK_OUTDIR, f"{{base}}.mag")
    txt = os.path.join(WORK_OUTDIR, f"{{base}}.txt")
    safe_rm(coo); safe_rm(mag); safe_rm(txt)
    safe_rm(final_coo); safe_rm(final_mag); safe_rm(final_txt)

    log(f"[PROGRESS] {{idx+1}}/{{len(imgs)}} {{im}} band={{band}}")
    log(f"  fwhm={{fwhm_pix:.2f}}px({{fwhm_source}}) sigma={{sig:.1f}}({{sigma_source}}) thr={{thr:.2f}} sharplo={{sharplo}} datamin={{datamin}}")
    log(f"  ap={{aperture:.2f}}px ann={{annulus:.2f}}px dann={{dannulus:.2f}}px cbox={{cbox:.2f}}px")

    try:
        t0 = time.perf_counter()
        iraf.daofind(im, output=coo, verify="no", interactive="no", verbose="no")
        t_daofind = time.perf_counter()
        iraf.phot(im, coords=coo, output=mag, verify="no", interactive="no", verbose="no", Stdout=os.devnull)
        t_phot = time.perf_counter()
        iraf.txdump(mag, fields="ID,XCENTER,YCENTER,MAG,MERR,MSKY,STDEV,NSKY",
                   expr="yes", headers="no", Stdout=txt)
        t_txdump = time.perf_counter()

        t_copyout = t_txdump
        if os.path.abspath(WORK_OUTDIR) != os.path.abspath(FINAL_OUTDIR):
            for work_path, final_path in ((coo, final_coo), (mag, final_mag), (txt, final_txt)):
                if os.path.exists(work_path):
                    shutil.copy2(work_path, final_path)
            t_copyout = time.perf_counter()

        n_stars = 0
        try:
            with open(txt, 'r') as f:
                n_stars = sum(1 for _ in f)
        except:
            pass
        try:
            apex_n_sources = float(frame_params.get("apex_n_sources", 0) or 0)
        except:
            apex_n_sources = 0.0
        if apex_n_sources > 0 and n_stars > max(2000, apex_n_sources * 1.5):
            log(
                f"  [QC] IRAF detected many sources: {{n_stars}} vs APEX Step4 {{int(apex_n_sources)}}. "
                f"Consider raising FINDPARS threshold (current {{thr:.2f}}) for faster comparison runs."
            )
        elif n_stars > 2500:
            log(
                f"  [QC] IRAF detected many sources: {{n_stars}}. "
                f"Consider raising FINDPARS threshold (current {{thr:.2f}}) for faster comparison runs."
            )
        copy_part = f" copyout={{t_copyout-t_txdump:.1f}}s" if os.path.abspath(WORK_OUTDIR) != os.path.abspath(FINAL_OUTDIR) else ""
        copyin_part = f"copyin={{copyin_seconds:.1f}}s " if USE_WSL_SCRATCH else ""
        log(f"  -> {{n_stars}} stars detected ({{copyin_part}}daofind={{t_daofind-t0:.1f}}s phot={{t_phot-t_daofind:.1f}}s txdump={{t_txdump-t_phot:.1f}}s{{copy_part}})")
        processed += 1
    except Exception as e:
        log(f"  -> ERROR: {{e}}")
    if work_image_path:
        safe_rm(work_image_path)

log(f"[DONE] Processed: {{processed}}, Skipped: {{skipped}}, Output: {{OUTDIR}}")
'''


# ============================================================================
# Worker Thread
# ============================================================================
class IRAFPhotometryWorker(QThread):
    """Worker thread for running IRAF photometry via WSL PyRAF."""

    progress = pyqtSignal(int, int, str)
    log = pyqtSignal(str)
    finished = pyqtSignal(dict)
    error = pyqtSignal(str)

    def __init__(self, data_dir: Path, output_dir: Path, file_pattern: str,
                 params: IRAFParameters, auto_sigma: bool = True, skip_existing: bool = True,
                 filter_params: dict | None = None, filter_aliases: dict | None = None,
                 param_defaults: dict | None = None, apex_frame_params: dict | None = None,
                 use_apex_frame_params: bool = True, use_wsl_scratch: bool = False,
                 scratch_root: str = "/tmp/apex_iraf"):
        super().__init__()
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)
        self.file_pattern = file_pattern
        self.params = params
        self.auto_sigma = auto_sigma
        self.skip_existing = skip_existing
        self.filter_params = filter_params or {}
        self.filter_aliases = filter_aliases or {}
        self.param_defaults = param_defaults or {}
        self.apex_frame_params = apex_frame_params or {}
        self.use_apex_frame_params = bool(use_apex_frame_params)
        self.use_wsl_scratch = bool(use_wsl_scratch)
        self.scratch_root = str(scratch_root or "/tmp/apex_iraf")
        self._stop_requested = False
        self._process = None
        self._script_path = None
        self._wsl_script_path = None
        self._scratch_run_dir = None

    def stop(self):
        self._stop_requested = True
        if sys.platform == "win32" and self._wsl_script_path:
            try:
                subprocess.run(
                    ["wsl", "pkill", "-TERM", "-f", self._wsl_script_path],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    timeout=5,
                )
            except Exception:
                pass
        if self._process:
            self._process.terminate()

    def _log(self, msg: str):
        self.log.emit(msg)

    def _runtime_cmd(self, is_windows: bool) -> list[str]:
        return _iraf_runtime_cmd(is_windows)

    def _cleanup_scratch_run_dir(self, is_windows: bool):
        scratch_dir = str(self._scratch_run_dir or "").replace("\\", "/").rstrip("/")
        scratch_root = str(self.scratch_root or "/tmp/apex_iraf").replace("\\", "/").rstrip("/")
        if not self.use_wsl_scratch or not scratch_dir:
            return
        name = scratch_dir.rsplit("/", 1)[-1]
        if not scratch_dir.startswith(scratch_root + "/") or not name.startswith("apex_iraf_"):
            self._log(f"[SCRATCH] Cleanup skipped for unexpected path: {scratch_dir}")
            return
        try:
            if is_windows:
                subprocess.run(
                    ["wsl", "rm", "-rf", "--", scratch_dir],
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                    timeout=30,
                )
            else:
                shutil.rmtree(scratch_dir, ignore_errors=True)
            self._log(f"[SCRATCH] Worker cleanup removed: {scratch_dir}")
        except Exception as exc:
            self._log(f"[SCRATCH] Worker cleanup failed: {exc}")

    def _check_runtime_environment(self, runtime_cmd: list[str], data_dir: str, output_dir: str):
        try:
            result = _run_iraf_environment_probe(
                runtime_cmd, data_dir, output_dir, log_callback=self._log, timeout=30
            )
        except subprocess.TimeoutExpired as exc:
            out = _format_subprocess_output(exc.stdout)
            details = f"\n\nLast output:\n{out}" if out else ""
            raise RuntimeError(f"IRAF environment check timed out after 30s.{details}") from exc

        if not result["ok"]:
            raise RuntimeError(
                "IRAF/PyRAF environment check failed before photometry started.\n"
                "Install/configure PyRAF and IRAF DAOPHOT in the runtime used by this tool "
                "(WSL python3 on Windows, local python3 otherwise).\n\n"
                f"Command: {' '.join(result['command'])}\n"
                f"Exit code: {result['returncode']}\n\n{result['output']}"
            )

    @staticmethod
    def _looks_like_iraf_source_row(line: str) -> bool:
        parts = line.split()
        if len(parts) < 5:
            return False

        image_name = parts[0].split("[", 1)[0].lower()
        if not image_name.endswith((".fit", ".fits", ".fts", ".fit.gz", ".fits.gz", ".fts.gz")):
            return False

        try:
            float(parts[1])
            float(parts[2])
            float(parts[3])
        except ValueError:
            return False
        return True

    def run(self):
        try:
            self._script_path = None
            is_windows = sys.platform == "win32"

            if is_windows:
                data_dir_str = windows_to_wsl_path(str(self.data_dir))
                output_dir_str = windows_to_wsl_path(str(self.output_dir))
            else:
                data_dir_str = str(self.data_dir)
                output_dir_str = str(self.output_dir)

            runtime_cmd = self._runtime_cmd(is_windows)
            self._log("[CHECK] Verifying PyRAF/IRAF runtime...")
            self._check_runtime_environment(runtime_cmd, data_dir_str, output_dir_str)

            if self.use_wsl_scratch:
                root = self.scratch_root.replace("\\", "/").rstrip("/") or "/tmp/apex_iraf"
                self._scratch_run_dir = f"{root}/apex_iraf_{int(time.time())}_{os.getpid()}_{id(self) & 0xffff:x}"
            else:
                self._scratch_run_dir = None

            p = self.params
            script_content = PYRAF_SCRIPT_TEMPLATE.format(
                data_dir=data_dir_str,
                output_dir=output_dir_str,
                file_pattern=self.file_pattern,
                skip_existing="True" if self.skip_existing else "False",
                use_wsl_scratch="True" if self.use_wsl_scratch else "False",
                scratch_root=self.scratch_root.replace("\\", "/"),
                scratch_run_dir=(self._scratch_run_dir or ""),
                filter_params=repr(self.filter_params),
                filter_aliases=repr(self.filter_aliases),
                param_defaults=repr(self.param_defaults),
                apex_frame_params=repr(self.apex_frame_params),
                use_apex_frame_params="True" if self.use_apex_frame_params else "False",
                # Pixel scale
                pix_scale=p.pix_scale,
                # Filter-specific seeing
                seeing_g=p.seeing_g,
                seeing_r=p.seeing_r,
                seeing_i=p.seeing_i,
                seeing_default=p.seeing_default,
                # Filter-specific sigma
                sigma_g=p.sigma_g,
                sigma_r=p.sigma_r,
                sigma_i=p.sigma_i,
                sigma_default=p.sigma_default,
                sigma_ref=p.sigma_ref,
                # DATAPARS
                scale=p.scale,
                emission="True" if p.emission else "False",
                datamax=p.datamax,
                noise=p.noise,
                readnoise=p.readnoise,
                epadu=p.epadu,
                exposure=p.exposure,
                itime=p.itime,
                # Filter-specific FINDPARS
                threshold_g=p.threshold_g,
                threshold_r=p.threshold_r,
                threshold_i=p.threshold_i,
                threshold_default=p.threshold_default,
                sharplo_g=p.sharplo_g,
                sharplo_r=p.sharplo_r,
                sharplo_i=p.sharplo_i,
                sharplo_default=p.sharplo_default,
                datamin_g=p.datamin_g,
                datamin_r=p.datamin_r,
                datamin_i=p.datamin_i,
                datamin_default=p.datamin_default,
                # FINDPARS base
                nsigma=p.nsigma,
                ratio=p.ratio,
                theta=p.theta,
                sharphi=p.sharphi,
                roundlo=p.roundlo,
                roundhi=p.roundhi,
                # CENTERPARS
                calgorithm=p.calgorithm,
                cbox_mult=p.cbox_mult,
                cthreshold=p.cthreshold,
                minsnratio=p.minsnratio,
                cmaxiter=p.cmaxiter,
                maxshift=p.maxshift,
                clean="True" if p.clean else "False",
                rclean=p.rclean,
                rclip=p.rclip,
                kclean=p.kclean,
                # FITSKYPARS (FWHM multipliers)
                salgorithm=p.salgorithm,
                annulus_mult=p.annulus_mult,
                dannulus_mult=p.dannulus_mult,
                skyvalue=p.skyvalue,
                smaxiter=p.smaxiter,
                sloclip=p.sloclip,
                shiclip=p.shiclip,
                snreject=p.snreject,
                sloreject=p.sloreject,
                shireject=p.shireject,
                khist=p.khist,
                binsize=p.binsize,
                smooth="True" if p.smooth else "False",
                rgrow=p.rgrow,
                # PHOTPARS (FWHM multiplier)
                aperture_mult=p.aperture_mult,
                zmag=p.zmag,
                mkapert="True" if p.mkapert else "False",
                # Additional
                auto_sigma="True" if self.auto_sigma else "False",
            )

            # Save script
            if is_windows:
                script_dir = Path(self.output_dir)
                script_dir.mkdir(parents=True, exist_ok=True)
                script_path = script_dir / "_pyraf_photometry.py"
                script_path.write_text(script_content, encoding="utf-8")
                self._script_path = script_path
                wsl_script_path = windows_to_wsl_path(str(script_path))
                self._wsl_script_path = wsl_script_path
                self._log(f"Script: {script_path}")
                cmd = _unbuffered_python_cmd(runtime_cmd) + [wsl_script_path]
            else:
                with tempfile.NamedTemporaryFile(mode='w', suffix='.py', delete=False) as f:
                    f.write(script_content)
                script_path = f.name
                self._script_path = Path(script_path)
                self._wsl_script_path = script_path
                self._log(f"Script: {script_path}")
                cmd = _unbuffered_python_cmd(runtime_cmd) + [script_path]

            self._log("Starting PyRAF via WSL (unbuffered)..." if is_windows else "Starting PyRAF (unbuffered)...")

            self._process = subprocess.Popen(
                cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                bufsize=1, **_SUBPROCESS_TEXT_KWARGS,
            )

            total_images = 0
            results = []
            suppressed_source_rows = 0
            current_image_index = 0
            current_image_name = ""
            output_tail: list[str] = []

            for line in iter(self._process.stdout.readline, ''):
                if self._stop_requested:
                    self._process.terminate()
                    break

                line = line.strip()
                if not line:
                    continue

                if self._looks_like_iraf_source_row(line):
                    suppressed_source_rows += 1
                    if suppressed_source_rows == 1:
                        self._log("[IRAF] Suppressing per-source phot rows; detailed rows are in phot output files.")
                    continue

                self._log(line)
                output_tail.append(line)
                if len(output_tail) > 80:
                    output_tail = output_tail[-80:]

                if line.startswith("[INFO] Found"):
                    try:
                        total_images = int(line.split()[2])
                        self.progress.emit(0, total_images, "Starting")
                    except:
                        pass
                elif line.startswith("[SCRATCH] Run dir:"):
                    self._scratch_run_dir = line.split(":", 1)[1].strip()
                elif line.startswith("[PROGRESS]"):
                    try:
                        parts = line.split(maxsplit=3)
                        frame_count = parts[1].split("/")
                        current_image_index = int(frame_count[0])
                        if len(frame_count) > 1:
                            total_images = int(frame_count[1])
                        current_image_name = parts[2] if len(parts) > 2 else ""
                        completed = max(current_image_index - 1, 0)
                        self.progress.emit(
                            completed,
                            total_images,
                            f"Running {current_image_index}/{total_images}: {current_image_name}",
                        )
                    except:
                        pass
                elif line.startswith("[SKIP]"):
                    try:
                        parts = line.split(maxsplit=3)
                        frame_count = parts[1].split("/")
                        current = int(frame_count[0])
                        if len(frame_count) > 1:
                            total_images = int(frame_count[1])
                        skipped_name = parts[2] if len(parts) > 2 else ""
                        self.progress.emit(current, total_images, f"Skipped {skipped_name}")
                    except:
                        pass
                elif "stars detected" in line:
                    try:
                        n_stars = int(line.split("->")[1].split()[0])
                        results.append({"n_stars": n_stars})
                        self.progress.emit(
                            current_image_index,
                            total_images,
                            f"Completed {current_image_name} ({n_stars} stars)",
                        )
                    except:
                        pass
                elif line.startswith("-> ERROR"):
                    if current_image_index:
                        self.progress.emit(
                            current_image_index,
                            total_images,
                            f"Error on {current_image_name}",
                        )

            if self._stop_requested:
                try:
                    self._process.wait(timeout=10)
                except subprocess.TimeoutExpired:
                    self._log("[STOP] PyRAF did not exit after terminate; killing process.")
                    self._process.kill()
                    self._process.wait()
            else:
                self._process.wait()

            if suppressed_source_rows:
                self._log(f"[IRAF] Suppressed {suppressed_source_rows} per-source phot rows.")

            if self._stop_requested:
                self._log(f"[STOP] IRAF photometry stopped by user (exit code {self._process.returncode}).")
                self._cleanup_scratch_run_dir(is_windows)
                if self._script_path:
                    try:
                        self._script_path.unlink()
                        self._log(f"[CLEANUP] Removed script: {self._script_path}")
                    except Exception as e:
                        self._log(f"[CLEANUP] Failed to remove script: {e}")
                self.finished.emit({
                    "results": results,
                    "output_dir": str(self.output_dir),
                    "stopped": True,
                    "returncode": self._process.returncode,
                })
            elif self._process.returncode == 0:
                self._log("\n[SUCCESS] PyRAF photometry completed!")
                self._cleanup_scratch_run_dir(is_windows)
                if self._script_path:
                    try:
                        self._script_path.unlink()
                        self._log(f"[CLEANUP] Removed script: {self._script_path}")
                    except Exception as e:
                        self._log(f"[CLEANUP] Failed to remove script: {e}")
                self.finished.emit({"results": results, "output_dir": str(self.output_dir)})
            else:
                self._cleanup_scratch_run_dir(is_windows)
                tail = "\n".join(output_tail[-12:])
                details = f"\n\nLast output:\n{tail}" if tail else ""
                self.error.emit(f"PyRAF exited with code {self._process.returncode}{details}")

        except FileNotFoundError as e:
            missing = str(getattr(e, "filename", "") or "")
            if missing == "wsl" or sys.platform == "win32":
                self.error.emit("WSL not found. Install WSL, then install PyRAF/IRAF inside WSL.")
            else:
                self.error.emit(f"Python runtime not found: {missing or 'python3'}")
        except Exception as e:
            self.error.emit(f"Error: {e}")


class IRAFEnvironmentCheckWorker(QThread):
    """Worker for IRAF dependency and path status checks."""

    status = pyqtSignal(str, str, str)
    log = pyqtSignal(str)
    finished = pyqtSignal(dict)

    def __init__(self, data_dir: Path, output_dir: Path, file_pattern: str):
        super().__init__()
        self.data_dir = Path(data_dir)
        self.output_dir = Path(output_dir)
        self.file_pattern = str(file_pattern or "*.fit*").strip() or "*.fit*"

    def _emit_status(self, key: str, state: str, message: str):
        self.status.emit(key, state, message)

    def run(self):
        statuses: dict[str, tuple[str, str]] = {}

        def set_status(key: str, state: str, message: str):
            statuses[key] = (state, message)
            self._emit_status(key, state, message)

        is_windows = sys.platform == "win32"
        runtime_cmd = _iraf_runtime_cmd(is_windows)

        try:
            if not self.data_dir.exists():
                set_status("images", "fail", "Data directory missing")
            elif not self.data_dir.is_dir():
                set_status("images", "fail", "Data path is not a directory")
            else:
                try:
                    image_count = sum(1 for _ in self.data_dir.glob(self.file_pattern))
                    if image_count:
                        set_status("images", "ok", f"{image_count} file(s)")
                    else:
                        set_status("images", "fail", f"No match: {self.file_pattern}")
                except Exception as exc:
                    set_status("images", "fail", f"Pattern error: {exc}")

            data_dir_str = windows_to_wsl_path(str(self.data_dir)) if is_windows else str(self.data_dir)
            output_dir_str = windows_to_wsl_path(str(self.output_dir)) if is_windows else str(self.output_dir)
            result = _run_iraf_environment_probe(
                runtime_cmd, data_dir_str, output_dir_str, log_callback=self.log.emit, timeout=30
            )
            for key, (state, message) in result["statuses"].items():
                set_status(key, state, message)
            if "runtime" not in statuses:
                set_status("runtime", "ok" if result["ok"] else "fail", "Runtime checked")
            overall_ok = bool(result["ok"]) and not any(state == "fail" for state, _ in statuses.values())
            self.finished.emit({"ok": overall_ok, "statuses": statuses, "output": result["output"]})
        except FileNotFoundError as exc:
            missing = str(getattr(exc, "filename", "") or "python3")
            label = "WSL not found" if is_windows else f"{missing} not found"
            set_status("runtime", "fail", label)
            set_status("pyraf", "warn", "Not checked")
            set_status("daophot", "warn", "Not checked")
            self.finished.emit({"ok": False, "statuses": statuses, "output": str(exc)})
        except subprocess.TimeoutExpired as exc:
            out = _format_subprocess_output(exc.stdout)
            set_status("runtime", "fail", "Check timed out")
            self.finished.emit({"ok": False, "statuses": statuses, "output": out})
        except Exception as exc:
            set_status("runtime", "fail", str(exc))
            self.finished.emit({"ok": False, "statuses": statuses, "output": str(exc)})


# ============================================================================
# Comparison Functions (from iraf_comparison_window.py)
# ============================================================================
def _read_iraf_txt(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r"\s+",
                     names=["ID", "x", "y", "mag", "merr", "msky", "stdev", "nsky"],
                     engine="python")
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if "x" in df.columns:
        df["x"] = df["x"] - 1.0
    if "y" in df.columns:
        df["y"] = df["y"] - 1.0
    return df


def _read_iraf_coo(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        names=["x", "y", "mag", "sharp", "sround", "ground", "ID"],
        engine="python",
    )
    for col in ["ID", "x", "y", "mag"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df[["ID", "x", "y", "mag"]].copy()
    df["x"] = df["x"] - 1.0
    df["y"] = df["y"] - 1.0
    return df


def _read_apex_tsv(path: Path) -> pd.DataFrame:
    try:
        return pd.read_csv(path, sep="\t")
    except:
        return pd.read_csv(path)


def _normalize_frame_key(stem: str) -> str:
    key = str(stem)
    if key.startswith("photometry_"):
        key = key[len("photometry_"):]
    if key.startswith("detect_"):
        key = key[len("detect_"):]
    if key.endswith("_photometry"):
        key = key[:-len("_photometry")]
    if key.startswith("Crop_"):
        key = key[len("Crop_"):]
    low = key.lower()
    for ext in (".fits.gz", ".fit.gz", ".fits", ".fit", ".fts", ".fz", ".gz"):
        if low.endswith(ext):
            key = key[:-len(ext)]
            break
    return key


def _iter_apex_photometry_files(apex_dir: Path):
    patterns = ("*_photometry.tsv", "photometry_*.tsv")
    seen: set[Path] = set()
    for pattern in patterns:
        for path in apex_dir.rglob(pattern):
            if not path.is_file() or path in seen:
                continue
            seen.add(path)
            yield path


BASE_IRAF_SHIFT = -1.0
AUTO_SHIFT_THRESHOLD = 0.6


def _auto_axis_shift(delta_med: float) -> float:
    if not np.isfinite(delta_med):
        return 0.0
    if abs(delta_med - 1.0) <= AUTO_SHIFT_THRESHOLD:
        return 1.0
    if abs(delta_med + 1.0) <= AUTO_SHIFT_THRESHOLD:
        return -1.0
    return 0.0


def _pick_first(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    return None


# ============================================================================
# Main Window
# ============================================================================
class IRAFPhotometryWindow(AutoFitMixin, QMainWindow):
    """Comprehensive IRAF Photometry tool with parameters and comparison."""

    def __init__(self, params, data_dir: Path, result_dir: Path, project_state=None, parent=None):
        super().__init__(parent)
        self.app_params = params
        self.data_dir = Path(data_dir)
        self.result_dir = Path(result_dir)
        self.project_state = project_state
        self.iraf_params = IRAFParameters()
        self.worker = None
        self.env_worker = None
        self.run_log_window = None
        self.run_log_text = None
        self.env_status_labels = {}
        # Only the integration essentials — is the runtime/PyRAF/DAOPHOT stack
        # actually wired up. Path/image checks still run in the worker but the
        # paths are already visible in their own fields, so they're not shown
        # as status cards (keeps the window uncluttered).
        self.env_status_defs = (
            ("runtime", "WSL/Python" if sys.platform == "win32" else "Python"),
            ("pyraf", "PyRAF"),
            ("daophot", "IRAF/DAOPHOT"),
        )

        # Comparison data
        self.frame_rows = []
        self.frame_matches = {}
        self.matched_all = None

        # Overlay data (comparison)
        self.overlay_file_list = []
        self.overlay_keys = []
        self.overlay_image_map = {}
        self.overlay_key_to_index = {}
        self.overlay_last_image_dir = None
        self.overlay_current_index = 0
        self.overlay_image_data = None
        self.overlay_header = None
        self.overlay_apex_map = {}
        self.overlay_iraf_map = {}
        self.overlay_filter_cache = {}
        self._overlay_normalized_cache = None
        self.overlay_viewer: FITSViewerWidget | None = None
        self.overlay_xlim_original = None
        self.overlay_ylim_original = None
        self.overlay_panning = False
        self.overlay_pan_start = None
        self._overlay_imshow_obj = None

        # Stretch plot window (2D Plot)
        self.overlay_stretch_plot_dialog = None
        self.overlay_stretch_plot_canvas = None
        self.overlay_stretch_plot_ax = None
        self.overlay_stretch_plot_fig = None
        self.overlay_stretch_plot_info_label = None
        self._overlay_stretch_vmin = None
        self._overlay_stretch_vmax = None
        self._overlay_stretch_data_range = None
        self._overlay_stretch_dragging = False
        self._overlay_stretch_drag_target = None
        self._overlay_stretch_marker_min_line = None
        self._overlay_stretch_marker_max_line = None

        self.setWindowTitle("IRAF/DAOPHOT Photometry Tool")
        # Minimum clamped to the monitor; AutoFitMixin fits to content on show.
        self.setMinimumSize(*clamp_to_screen(1200, 900, self))
        self._init_project_path_widgets()
        self.setup_ui()

    def _default_project_image_dir(self) -> Path:
        cropped = step2_cropped_dir(self.result_dir)
        if cropped.exists():
            return cropped
        legacy_cropped = step2_cropped_dir(self.data_dir / "result")
        if legacy_cropped.exists():
            return legacy_cropped
        return self.data_dir

    def _init_project_path_widgets(self):
        default_images = self._default_project_image_dir()
        default_iraf_dir = self.result_dir / "iraf_phot"

        self.data_edit = QLineEdit(str(default_images))
        self.out_edit = QLineEdit(str(default_iraf_dir))
        self.pattern_edit = QLineEdit("*.fit*")

        self.auto_sigma_check = QCheckBox("Auto-estimate sigma per image")
        self.auto_sigma_check.setChecked(True)
        self.skip_existing_check = QCheckBox("Skip already processed files")
        self.skip_existing_check.setChecked(True)
        self.use_apex_step4_check = QCheckBox("Use APEX Step 4 frame FWHM/sky when available")
        self.use_apex_step4_check.setChecked(True)
        self.use_wsl_scratch_check = QCheckBox("Use WSL scratch copy for IRAF run")
        self.use_wsl_scratch_check.setChecked(sys.platform == "win32")
        self.use_wsl_scratch_check.setToolTip(
            "Copy FITS files to WSL /tmp before IRAF, then copy results back. "
            "This is faster for external drives and /mnt paths but needs WSL disk space."
        )

        self.cmp_apex_edit = QLineEdit(str(self.result_dir))
        self.cmp_iraf_edit = QLineEdit(str(default_iraf_dir))
        self.cmp_image_edit = QLineEdit(str(default_images))
        self.cmp_tol = QDoubleSpinBox()
        self.cmp_tol.setRange(0.1, 20.0)
        self.cmp_tol.setValue(1.5)
        self.cmp_tol.valueChanged.connect(lambda _value: self._refresh_path_summaries())

        for edit in (self.data_edit, self.out_edit, self.pattern_edit):
            edit.editingFinished.connect(self._on_photometry_path_settings_changed)
        self.use_wsl_scratch_check.stateChanged.connect(lambda _value: self._refresh_path_summaries())
        for edit in (self.cmp_apex_edit, self.cmp_iraf_edit):
            edit.editingFinished.connect(self._on_comparison_path_settings_changed)
        self.cmp_image_edit.editingFinished.connect(self._on_comparison_image_path_changed)

    def _path_row(self, line_edit: QLineEdit) -> QWidget:
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(line_edit)
        btn = QPushButton("Browse")
        btn.clicked.connect(lambda: self._browse_dir(line_edit))
        layout.addWidget(btn)
        return row

    def _reset_project_paths(self):
        default_images = self._default_project_image_dir()
        default_iraf_dir = self.result_dir / "iraf_phot"
        self.data_edit.setText(str(default_images))
        self.out_edit.setText(str(default_iraf_dir))
        self.pattern_edit.setText("*.fit*")
        self.use_wsl_scratch_check.setChecked(sys.platform == "win32")
        self.cmp_apex_edit.setText(str(self.result_dir))
        self.cmp_iraf_edit.setText(str(default_iraf_dir))
        self.cmp_image_edit.setText(str(default_images))
        self._on_photometry_path_settings_changed()
        self._on_comparison_image_path_changed()

    def _on_photometry_path_settings_changed(self):
        if hasattr(self, "env_status_labels"):
            self._mark_environment_pending()
        self._refresh_path_summaries()

    def _on_comparison_path_settings_changed(self):
        self._refresh_path_summaries()
        if hasattr(self, "overlay_file_combo"):
            self._overlay_refresh_maps()

    def _on_comparison_image_path_changed(self):
        self._refresh_path_summaries()
        if hasattr(self, "overlay_file_combo"):
            self._overlay_reload()

    def _refresh_path_summaries(self):
        if hasattr(self, "run_path_summary"):
            self.run_path_summary.setText(
                "Data Directory: "
                f"{self.data_edit.text()}\n"
                "Output Directory: "
                f"{self.out_edit.text()}\n"
                f"File Pattern: {self.pattern_edit.text().strip() or '*.fit*'}\n"
                "WSL Scratch: "
                f"{'enabled' if getattr(self, 'use_wsl_scratch_check', None) and self.use_wsl_scratch_check.isChecked() else 'disabled'}"
            )
        if hasattr(self, "cmp_path_summary"):
            self.cmp_path_summary.setText(
                f"APEX: {self.cmp_apex_edit.text()}\n"
                f"IRAF: {self.cmp_iraf_edit.text()}\n"
                f"Images: {self.cmp_image_edit.text()} | Tol: {self.cmp_tol.value():.2f} px | Scratch: temp"
            )

    @staticmethod
    def _format_bytes(nbytes: int | float | None) -> str:
        try:
            value = float(nbytes or 0)
        except Exception:
            value = 0.0
        units = ("B", "KB", "MB", "GB", "TB")
        for unit in units:
            if value < 1024.0 or unit == units[-1]:
                return f"{int(value)} B" if unit == "B" else f"{value:.1f} {unit}"
            value /= 1024.0
        return f"{value:.1f} TB"

    @staticmethod
    def _scratch_required_bytes(input_bytes: int) -> int:
        return int(max(0, input_bytes) * 2.0 + 512 * 1024 * 1024)

    @staticmethod
    def _sum_file_sizes(paths: list[Path]) -> int:
        total = 0
        for path in paths:
            try:
                total += path.stat().st_size
            except OSError:
                pass
        return total

    def _wsl_free_bytes(self, wsl_path: str = "/tmp") -> int | None:
        if sys.platform != "win32":
            return None
        try:
            result = subprocess.run(
                ["wsl", "df", "-Pk", wsl_path],
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                timeout=10,
                **_SUBPROCESS_TEXT_KWARGS,
            )
        except Exception as exc:
            self._scratch_last_check = {"error": str(exc)}
            return None
        output = _format_subprocess_output(result.stdout, max_chars=2000)
        self._scratch_last_check = {"df_output": output}
        if result.returncode != 0:
            return None
        lines = [line for line in output.splitlines() if line.strip()]
        if len(lines) < 2:
            return None
        parts = lines[-1].split()
        if len(parts) < 4:
            return None
        try:
            return int(parts[3]) * 1024
        except ValueError:
            return None

    def _confirm_wsl_scratch_run(self, image_count: int, input_bytes: int) -> bool:
        self._scratch_last_check = {}
        required = self._scratch_required_bytes(input_bytes)
        free = self._wsl_free_bytes("/tmp")
        self._scratch_last_check.update({
            "enabled": True,
            "image_count": image_count,
            "input_bytes": input_bytes,
            "required_bytes": required,
            "free_bytes": free,
            "scratch_root": "/tmp/apex_iraf",
        })

        if free is not None and free < required:
            QMessageBox.critical(
                self,
                "WSL Scratch Space",
                "Not enough free space in WSL /tmp for scratch photometry.\n\n"
                f"Input FITS: {self._format_bytes(input_bytes)} ({image_count} files)\n"
                f"Required: {self._format_bytes(required)}\n"
                f"Free: {self._format_bytes(free)}\n\n"
                "Disable WSL scratch or free WSL disk space.",
            )
            return False

        if free is None:
            choice = QMessageBox.warning(
                self,
                "WSL Scratch Space",
                "Could not verify free space in WSL /tmp.\n\n"
                f"Input FITS: {self._format_bytes(input_bytes)} ({image_count} files)\n"
                f"Estimated required: {self._format_bytes(required)}\n"
                "Scratch path: /tmp/apex_iraf\n\n"
                "Continue anyway?",
                QMessageBox.Yes | QMessageBox.No,
                QMessageBox.No,
            )
            return choice == QMessageBox.Yes

        choice = QMessageBox.warning(
            self,
            "WSL Scratch Copy",
            "IRAF will run from WSL scratch, then copy results back to the project output folder.\n\n"
            f"Input FITS: {self._format_bytes(input_bytes)} ({image_count} files)\n"
            f"Estimated required: {self._format_bytes(required)}\n"
            f"WSL /tmp free: {self._format_bytes(free)}\n"
            "Scratch path: /tmp/apex_iraf\n\n"
            "Continue?",
            QMessageBox.Yes | QMessageBox.No,
            QMessageBox.Yes,
        )
        return choice == QMessageBox.Yes

    def show_parameters_tab(self):
        if hasattr(self, "tabs"):
            self.tabs.setCurrentIndex(1)
        if hasattr(self, "param_tabs"):
            self.param_tabs.setCurrentIndex(0)

    def setup_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        layout = QVBoxLayout(central)

        # Main tab widget
        self.tabs = QTabWidget()

        # Tab 1: Run Photometry
        self.tabs.addTab(self._create_run_tab(), "Run Photometry")

        # Tab 2: IRAF Parameters
        self.tabs.addTab(self._create_params_tab(), "IRAF Parameters")

        # Tab 3: Comparison
        self.tabs.addTab(self._create_comparison_tab(), "Comparison")

        layout.addWidget(self.tabs)

    def closeEvent(self, event):
        self._auto_save_params()
        super().closeEvent(event)

    def _create_environment_status_group(self):
        group = QGroupBox("Environment Status")
        layout = QVBoxLayout(group)

        grid = QGridLayout()
        grid.setSpacing(8)
        self.env_status_labels = {}
        for idx, (key, title) in enumerate(self.env_status_defs):
            label = QLabel()
            label.setAlignment(Qt.AlignCenter)
            label.setWordWrap(True)
            label.setMinimumHeight(58)
            self.env_status_labels[key] = label
            grid.addWidget(label, idx // 3, idx % 3)
            self._set_environment_status(key, "pending", "Not checked")
        layout.addLayout(grid)

        btn_row = QHBoxLayout()
        self.check_env_btn = QPushButton("Check Environment")
        self.check_env_btn.clicked.connect(self.check_iraf_environment)
        btn_row.addWidget(self.check_env_btn)
        self.env_summary_label = QLabel("Run this before photometry on a new PC or after changing paths.")
        self.env_summary_label.setProperty("role", "caption")
        btn_row.addWidget(self.env_summary_label, stretch=1)
        layout.addLayout(btn_row)
        return group

    def _set_environment_status(self, key: str, state: str, message: str):
        label = self.env_status_labels.get(key)
        if label is None:
            return
        title = dict(self.env_status_defs).get(key, key)
        state_key = str(state or "pending").lower()
        # State colours from the live theme tokens — the fixed light pastels
        # glared on every dark preset.
        palette = {
            "ok": (Tokens.OK_SOFT, Tokens.OK, "OK"),
            "fail": (Tokens.ERROR_SOFT, Tokens.ERROR, "FAIL"),
            "warn": (Tokens.WARN_SOFT, Tokens.WARN, "CHECK"),
            "running": (Tokens.ACCENT_SOFT, Tokens.ACCENT_TEXT, "CHECKING"),
            "pending": (Tokens.SURFACE_ALT, Tokens.TEXT_MUTED, "NOT CHECKED"),
        }
        bg, fg, status_text = palette.get(state_key, palette["pending"])
        label.setText(f"{title}\n{status_text}: {message}")
        label.setStyleSheet(
            f"QLabel {{ background: {bg}; color: {fg}; border: 1px solid {fg}; "
            "border-radius: 6px; padding: 8px; font-weight: 600; }}"
        )

    def _mark_environment_pending(self):
        for key, _ in self.env_status_defs:
            self._set_environment_status(key, "pending", "Not checked")
        if hasattr(self, "env_summary_label"):
            self.env_summary_label.setText("Status is stale. Run Check Environment.")

    def check_iraf_environment(self):
        data_dir = Path(self.data_edit.text())
        output_dir = Path(self.out_edit.text())
        pattern = self.pattern_edit.text().strip() or "*.fit*"

        for key, _ in self.env_status_defs:
            self._set_environment_status(key, "running", "Checking")
        self.env_summary_label.setText("Checking IRAF runtime and paths...")
        self.check_env_btn.setEnabled(False)
        self.run_btn.setEnabled(False)

        self.env_worker = IRAFEnvironmentCheckWorker(data_dir, output_dir, pattern)
        self.env_worker.status.connect(self._on_environment_status)
        self.env_worker.log.connect(self._on_log)
        self.env_worker.finished.connect(self._on_environment_finished)
        self.env_worker.start()

    def _on_environment_status(self, key: str, state: str, message: str):
        self._set_environment_status(key, state, message)

    def _on_environment_finished(self, result: dict):
        ok = bool(result.get("ok"))
        self.check_env_btn.setEnabled(True)
        self.run_btn.setEnabled(True)
        self.env_summary_label.setText(
            "Environment ready." if ok else "Environment check failed. See red status boxes and log."
        )

    # ========================================================================
    # Tab 1: Run Photometry
    # ========================================================================
    def _create_run_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        # Info
        info = QLabel(
            "Run IRAF DAOPHOT photometry via PyRAF (WSL python3 on Windows, local python3 otherwise).\n"
            "Configure parameters in the 'IRAF Parameters' tab before running."
        )
        info.setProperty("role", "caption")
        layout.addWidget(info)

        controls = QHBoxLayout()
        btn_params = create_parameter_button("IRAF Parameters")
        btn_params.clicked.connect(self.show_parameters_tab)
        controls.addWidget(btn_params)
        controls.addStretch()
        layout.addLayout(controls)

        paths_group = QGroupBox("Project Paths")
        paths_layout = QVBoxLayout(paths_group)
        self.run_path_summary = QLabel()
        self.run_path_summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.run_path_summary.setStyleSheet(mono_note_style())
        paths_layout.addWidget(self.run_path_summary)
        layout.addWidget(paths_group)

        layout.addWidget(self._create_environment_status_group())

        # Buttons
        self.run_bar = RunControlBar(
            "Run IRAF Photometry",
            "Log",
            run_cb=self.run_photometry,
            stop_cb=self.stop_photometry,
            log_cb=self.show_run_log_window,
        )
        layout.addWidget(self.run_bar)
        self.run_btn = self.run_bar.btn_run
        self.stop_btn = self.run_bar.btn_stop

        # Progress
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        layout.addWidget(self.progress_bar)

        self.progress_label = QLabel("")
        layout.addWidget(self.progress_label)

        # Log
        log_group = QGroupBox("Execution Log")
        log_layout = QVBoxLayout()
        self.log_text = QTextEdit()
        self.log_text.setObjectName("Log")     # themed mono surface
        self.log_text.setReadOnly(True)
        self.log_text.setMinimumHeight(300)
        log_layout.addWidget(self.log_text)
        log_group.setLayout(log_layout)
        layout.addWidget(log_group)

        self._refresh_path_summaries()
        return tab

    # ========================================================================
    # Tab 2: IRAF Parameters
    # ========================================================================
    def _create_params_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)

        # Sub-tabs for each parameter group
        param_tabs = QTabWidget()
        self.param_tabs = param_tabs

        param_tabs.addTab(self._create_project_paths_panel(), "Project Paths")
        param_tabs.addTab(self._create_datapars_panel(), "DATAPARS")
        param_tabs.addTab(self._create_findpars_panel(), "FINDPARS")
        param_tabs.addTab(self._create_centerpars_panel(), "CENTERPARS")
        param_tabs.addTab(self._create_fitskypars_panel(), "FITSKYPARS")
        param_tabs.addTab(self._create_photpars_panel(), "PHOTPARS")

        layout.addWidget(param_tabs)

        btn_layout = QHBoxLayout()

        note = QLabel("Parameters auto-save on run/close.")
        note.setProperty("role", "caption")
        btn_layout.addWidget(note)
        btn_layout.addStretch()

        defaults_btn = QPushButton("Reset to Defaults")
        defaults_btn.clicked.connect(self._load_defaults)
        btn_layout.addWidget(defaults_btn)

        layout.addLayout(btn_layout)

        # Auto-load saved parameters if file exists
        self._auto_load_params()

        return tab

    def _create_project_paths_panel(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        panel = QWidget()
        layout = QVBoxLayout(panel)

        project_group = QGroupBox("Main Project")
        project_layout = QFormLayout(project_group)
        project_result = QLabel(str(self.result_dir))
        project_result.setTextInteractionFlags(Qt.TextSelectableByMouse)
        project_data = QLabel(str(self.data_dir))
        project_data.setTextInteractionFlags(Qt.TextSelectableByMouse)
        project_layout.addRow("APEX Result Dir:", project_result)
        project_layout.addRow("Raw Data Dir:", project_data)
        btn_defaults = create_cache_action_button("Use Project Defaults")
        btn_defaults.clicked.connect(self._reset_project_paths)
        project_layout.addRow("", btn_defaults)
        btn_sync = create_parameter_button("Sync From APEX Step 4/7")
        btn_sync.clicked.connect(self.sync_iraf_params_from_apex)
        project_layout.addRow("", btn_sync)
        self.apex_sync_status = QLabel("Step 4/7 sync has not been run in this session.")
        self.apex_sync_status.setWordWrap(True)
        self.apex_sync_status.setProperty("role", "caption")
        project_layout.addRow("APEX Sync:", self.apex_sync_status)
        layout.addWidget(project_group)

        run_group = QGroupBox("IRAF Photometry")
        run_layout = QFormLayout(run_group)
        run_layout.addRow("Data Directory:", self._path_row(self.data_edit))
        run_layout.addRow("Output Directory:", self._path_row(self.out_edit))
        run_layout.addRow("File Pattern:", self.pattern_edit)
        run_layout.addRow("", self.auto_sigma_check)
        run_layout.addRow("", self.skip_existing_check)
        run_layout.addRow("", self.use_apex_step4_check)
        run_layout.addRow("", self.use_wsl_scratch_check)
        layout.addWidget(run_group)

        compare_group = QGroupBox("Comparison")
        compare_layout = QFormLayout(compare_group)
        compare_layout.addRow("APEX Result Dir:", self._path_row(self.cmp_apex_edit))
        compare_layout.addRow("IRAF Result Dir:", self._path_row(self.cmp_iraf_edit))
        compare_layout.addRow("Image Dir:", self._path_row(self.cmp_image_edit))
        compare_layout.addRow("Match Tolerance (px):", self.cmp_tol)
        layout.addWidget(compare_group)

        layout.addStretch()
        scroll.setWidget(panel)
        return scroll

    @staticmethod
    def _first_column(df: pd.DataFrame, candidates: list[str]) -> str | None:
        names = {str(col).lower(): str(col) for col in df.columns}
        for name in candidates:
            col = names.get(str(name).lower())
            if col is not None:
                return col
        return None

    @staticmethod
    def _finite_median(values) -> float:
        arr = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(float)
        arr = arr[np.isfinite(arr)]
        return float(np.nanmedian(arr)) if arr.size else float("nan")

    @staticmethod
    def _frame_key_variants(filename) -> set[str]:
        raw = Path(str(filename)).name
        variants = {raw, Path(raw).stem}
        low = raw.lower()
        for ext in (".fits.gz", ".fit.gz", ".fits", ".fit", ".fts", ".fz", ".gz"):
            if low.endswith(ext):
                variants.add(raw[:-len(ext)])
                break
        return {v for v in variants if v}

    def _load_apex_step4_frame_table(self) -> tuple[pd.DataFrame, str]:
        s4_dir = step4_dir(self.result_dir)
        frame_quality = s4_dir / "frame_quality.csv"
        if frame_quality.exists():
            try:
                return pd.read_csv(frame_quality), str(frame_quality)
            except Exception:
                pass

        rows = []
        if s4_dir.exists():
            for path in sorted(s4_dir.glob("detect_*.json")):
                if path.name.startswith("detect_peak_"):
                    continue
                try:
                    payload = json.loads(path.read_text(encoding="utf-8"))
                except Exception:
                    continue
                fname = path.name[len("detect_"):-len(".json")]
                rows.append({
                    "file": fname,
                    "filter": payload.get("filter", ""),
                    "fwhm_px": payload.get("fwhm_px", np.nan),
                    "fwhm_arcsec": payload.get("fwhm_arcsec", np.nan),
                    "sky_sigma_med_adu": payload.get("bkg_rms", np.nan),
                    "sigma_used": payload.get("sigma_used", np.nan),
                    "n_sources": payload.get("n_sources", np.nan),
                })
        if rows:
            return pd.DataFrame(rows), str(s4_dir / "detect_*.json")
        return pd.DataFrame(), str(frame_quality)

    def _load_apex_frame_params(self) -> dict:
        df, _source = self._load_apex_step4_frame_table()
        if df.empty:
            return {}
        file_col = self._first_column(df, ["file", "filename", "File", "Filename"])
        fwhm_col = self._first_column(df, ["fwhm_med", "fwhm_px", "fwhm_med_px"])
        sigma_col = self._first_column(df, ["sky_sigma_med_adu", "sky_sigma", "bkg_rms"])
        n_sources_col = self._first_column(df, ["n_sources", "sources", "source_count"])
        if file_col is None or (fwhm_col is None and sigma_col is None and n_sources_col is None):
            return {}

        frame_params = {}
        for _, row in df.iterrows():
            entry = {}
            if fwhm_col is not None:
                fwhm = pd.to_numeric(pd.Series([row.get(fwhm_col)]), errors="coerce").iloc[0]
                if np.isfinite(float(fwhm)) and float(fwhm) > 0:
                    entry["fwhm_px"] = float(fwhm)
            if sigma_col is not None:
                sigma = pd.to_numeric(pd.Series([row.get(sigma_col)]), errors="coerce").iloc[0]
                if np.isfinite(float(sigma)) and float(sigma) > 0:
                    entry["sigma"] = float(sigma)
            if n_sources_col is not None:
                n_sources = pd.to_numeric(pd.Series([row.get(n_sources_col)]), errors="coerce").iloc[0]
                if np.isfinite(float(n_sources)) and float(n_sources) > 0:
                    entry["apex_n_sources"] = int(n_sources)
            if not entry:
                continue
            for key in self._frame_key_variants(row.get(file_col, "")):
                frame_params[key] = dict(entry)
        return frame_params

    def sync_iraf_params_from_apex(self):
        self._apply_params_silent()
        df, source = self._load_apex_step4_frame_table()
        if df.empty:
            message = f"No Step 4 frame quality data found at {source}"
            if hasattr(self, "apex_sync_status"):
                self.apex_sync_status.setText(message)
            QMessageBox.information(self, "APEX Sync", message)
            return

        p = self.iraf_params
        updates = []

        pix_scale = getattr(getattr(self.app_params, "P", object()), "pixel_scale_arcsec", np.nan)
        try:
            pix_scale = float(pix_scale)
        except Exception:
            pix_scale = np.nan
        if np.isfinite(pix_scale) and pix_scale > 0:
            p.pix_scale = pix_scale
            updates.append(f"pix_scale={pix_scale:.3f}\"/px")

        fwhm_col = self._first_column(df, ["fwhm_med", "fwhm_px", "fwhm_med_px"])
        sigma_col = self._first_column(df, ["sky_sigma_med_adu", "sky_sigma", "bkg_rms"])
        filter_col = self._first_column(df, ["filter", "FILTER"])

        if fwhm_col is not None:
            fwhm_med_px = self._finite_median(df[fwhm_col])
            if np.isfinite(fwhm_med_px) and fwhm_med_px > 0:
                p.seeing_default = fwhm_med_px * p.pix_scale
                updates.append(f"default seeing={p.seeing_default:.2f}\" ({fwhm_med_px:.2f}px)")

        if sigma_col is not None:
            sigma_med = self._finite_median(df[sigma_col])
            if np.isfinite(sigma_med) and sigma_med > 0:
                p.sigma_default = sigma_med
                p.sigma_ref = sigma_med
                updates.append(f"default sigma={sigma_med:.1f}")

        if filter_col is not None:
            work = df.copy()
            work["_iraf_filter_key"] = work[filter_col].map(normalize_filter_key)
            for filt_key, group in work.groupby("_iraf_filter_key"):
                if filt_key not in {"g", "r", "i"}:
                    continue
                if fwhm_col is not None:
                    fwhm_px = self._finite_median(group[fwhm_col])
                    if np.isfinite(fwhm_px) and fwhm_px > 0:
                        setattr(p, f"seeing_{filt_key}", fwhm_px * p.pix_scale)
                if sigma_col is not None:
                    sigma_val = self._finite_median(group[sigma_col])
                    if np.isfinite(sigma_val) and sigma_val > 0:
                        setattr(p, f"sigma_{filt_key}", sigma_val)

        P = getattr(self.app_params, "P", object())
        detect_sigma = getattr(P, "detect_sigma", None)
        try:
            detect_sigma = float(detect_sigma)
        except Exception:
            detect_sigma = np.nan
        if np.isfinite(detect_sigma) and detect_sigma > 0:
            p.threshold_default = float(detect_sigma)
        sigma_by_filter = getattr(P, "detect_sigma_by_filter", {}) or {}
        if isinstance(sigma_by_filter, dict):
            for filt, value in sigma_by_filter.items():
                key = normalize_filter_key(filt)
                if key not in {"g", "r", "i"}:
                    continue
                try:
                    val = float(value)
                except Exception:
                    continue
                if np.isfinite(val) and val > 0:
                    setattr(p, f"threshold_{key}", val)

        for attr, target in (
            ("forced_r_ap_scale", "aperture_mult"),
            ("fitsky_annulus_scale", "annulus_mult"),
            ("fitsky_dannulus_scale", "dannulus_mult"),
            ("center_cbox_scale", "cbox_mult"),
        ):
            try:
                val = float(getattr(P, attr))
            except Exception:
                continue
            if np.isfinite(val) and val > 0:
                setattr(p, target, val)
        updates.append(
            f"geometry ap={p.aperture_mult:.2f}x ann={p.annulus_mult:.2f}x "
            f"dann={p.dannulus_mult:.2f}x cbox={p.cbox_mult:.2f}x"
        )

        step7_path = step7_forced_phot_dir(self.result_dir) / "apcorr_summary.csv"
        if step7_path.exists():
            updates.append("Step 7 apcorr summary found for comparison context")

        self._update_ui_from_params()
        self.use_apex_step4_check.setChecked(True)
        self._auto_save_params()

        frame_count = len(self._load_apex_frame_params())
        message = f"Synced from {source}. Frame overrides available: {frame_count}. " + "; ".join(updates)
        if hasattr(self, "apex_sync_status"):
            self.apex_sync_status.setText(message)
        QMessageBox.information(self, "APEX Sync", message)

    def _create_datapars_panel(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        panel = QWidget()
        layout = QFormLayout(panel)

        # --- Pixel Scale ---
        layout.addRow(QLabel("--- Pixel Scale ---"))
        self.dp_pix_scale = QDoubleSpinBox()
        self.dp_pix_scale.setRange(0.01, 10.0)
        self.dp_pix_scale.setDecimals(3)
        self.dp_pix_scale.setValue(self.iraf_params.pix_scale)
        layout.addRow("pix_scale (arcsec/pix):", self.dp_pix_scale)

        # --- Filter-specific SEEING (FWHM in arcsec) ---
        layout.addRow(QLabel(""))
        layout.addRow(QLabel("--- Filter-specific Seeing (arcsec) ---"))

        self.dp_seeing_g = QDoubleSpinBox()
        self.dp_seeing_g.setRange(0.1, 20.0)
        self.dp_seeing_g.setDecimals(2)
        self.dp_seeing_g.setValue(self.iraf_params.seeing_g)
        layout.addRow("seeing (g band):", self.dp_seeing_g)

        self.dp_seeing_r = QDoubleSpinBox()
        self.dp_seeing_r.setRange(0.1, 20.0)
        self.dp_seeing_r.setDecimals(2)
        self.dp_seeing_r.setValue(self.iraf_params.seeing_r)
        layout.addRow("seeing (r band):", self.dp_seeing_r)

        self.dp_seeing_i = QDoubleSpinBox()
        self.dp_seeing_i.setRange(0.1, 20.0)
        self.dp_seeing_i.setDecimals(2)
        self.dp_seeing_i.setValue(self.iraf_params.seeing_i)
        layout.addRow("seeing (i band):", self.dp_seeing_i)

        self.dp_seeing_default = QDoubleSpinBox()
        self.dp_seeing_default.setRange(0.1, 20.0)
        self.dp_seeing_default.setDecimals(2)
        self.dp_seeing_default.setValue(self.iraf_params.seeing_default)
        layout.addRow("seeing (default):", self.dp_seeing_default)

        # --- Filter-specific SIGMA (sky background std) ---
        layout.addRow(QLabel(""))
        layout.addRow(QLabel("--- Filter-specific Sigma (sky std) ---"))

        self.dp_sigma_g = QDoubleSpinBox()
        self.dp_sigma_g.setRange(0.1, 10000.0)
        self.dp_sigma_g.setDecimals(1)
        self.dp_sigma_g.setValue(self.iraf_params.sigma_g)
        layout.addRow("sigma (g band):", self.dp_sigma_g)

        self.dp_sigma_r = QDoubleSpinBox()
        self.dp_sigma_r.setRange(0.1, 10000.0)
        self.dp_sigma_r.setDecimals(1)
        self.dp_sigma_r.setValue(self.iraf_params.sigma_r)
        layout.addRow("sigma (r band):", self.dp_sigma_r)

        self.dp_sigma_i = QDoubleSpinBox()
        self.dp_sigma_i.setRange(0.1, 10000.0)
        self.dp_sigma_i.setDecimals(1)
        self.dp_sigma_i.setValue(self.iraf_params.sigma_i)
        layout.addRow("sigma (i band):", self.dp_sigma_i)

        self.dp_sigma_default = QDoubleSpinBox()
        self.dp_sigma_default.setRange(0.1, 10000.0)
        self.dp_sigma_default.setDecimals(1)
        self.dp_sigma_default.setValue(self.iraf_params.sigma_default)
        layout.addRow("sigma (default):", self.dp_sigma_default)

        # Sigma reference (for threshold scaling)
        self.dp_sigma_ref = QDoubleSpinBox()
        self.dp_sigma_ref.setRange(1.0, 1000.0)
        self.dp_sigma_ref.setDecimals(1)
        self.dp_sigma_ref.setValue(self.iraf_params.sigma_ref)
        layout.addRow("sigma_ref (thr scaling):", self.dp_sigma_ref)

        # --- Other DATAPARS ---
        layout.addRow(QLabel(""))
        layout.addRow(QLabel("--- Other DATAPARS ---"))

        # Scale
        self.dp_scale = QDoubleSpinBox()
        self.dp_scale.setRange(0.01, 100.0)
        self.dp_scale.setDecimals(3)
        self.dp_scale.setValue(self.iraf_params.scale)
        layout.addRow("scale (units/pixel):", self.dp_scale)

        # Emission
        self.dp_emission = QCheckBox()
        self.dp_emission.setChecked(self.iraf_params.emission)
        layout.addRow("emission:", self.dp_emission)

        # Datamax
        self.dp_datamax = QDoubleSpinBox()
        self.dp_datamax.setRange(0, 1000000)
        self.dp_datamax.setDecimals(0)
        self.dp_datamax.setValue(self.iraf_params.datamax)
        layout.addRow("datamax:", self.dp_datamax)

        # Noise model
        self.dp_noise = QComboBox()
        self.dp_noise.addItems(["poisson", "constant", "file"])
        self.dp_noise.setCurrentText(self.iraf_params.noise)
        layout.addRow("noise model:", self.dp_noise)

        # Readnoise
        self.dp_readnoise = QDoubleSpinBox()
        self.dp_readnoise.setRange(0.0, 1000.0)
        self.dp_readnoise.setDecimals(2)
        self.dp_readnoise.setValue(self.iraf_params.readnoise)
        layout.addRow("readnoise (e-):", self.dp_readnoise)

        # Gain (epadu)
        self.dp_epadu = QDoubleSpinBox()
        self.dp_epadu.setRange(0.001, 1000.0)
        self.dp_epadu.setDecimals(3)
        self.dp_epadu.setValue(self.iraf_params.epadu)
        layout.addRow("epadu (e-/ADU):", self.dp_epadu)

        # Exposure keyword
        self.dp_exposure = QLineEdit(self.iraf_params.exposure)
        layout.addRow("exposure keyword:", self.dp_exposure)

        # Integration time
        self.dp_itime = QDoubleSpinBox()
        self.dp_itime.setRange(0.001, 100000.0)
        self.dp_itime.setDecimals(2)
        self.dp_itime.setValue(self.iraf_params.itime)
        layout.addRow("itime (default):", self.dp_itime)

        scroll.setWidget(panel)
        return scroll

    def _create_findpars_panel(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        panel = QWidget()
        layout = QFormLayout(panel)

        # --- Filter-specific Threshold ---
        layout.addRow(QLabel("--- Filter-specific Threshold ---"))
        layout.addRow(QLabel("(threshold = base × clip(sigma/sigma_ref, 0.8, 1.6))"))

        self.fp_threshold_g = QDoubleSpinBox()
        self.fp_threshold_g.setRange(0.1, 100.0)
        self.fp_threshold_g.setDecimals(2)
        self.fp_threshold_g.setValue(self.iraf_params.threshold_g)
        layout.addRow("threshold (g band):", self.fp_threshold_g)

        self.fp_threshold_r = QDoubleSpinBox()
        self.fp_threshold_r.setRange(0.1, 100.0)
        self.fp_threshold_r.setDecimals(2)
        self.fp_threshold_r.setValue(self.iraf_params.threshold_r)
        layout.addRow("threshold (r band):", self.fp_threshold_r)

        self.fp_threshold_i = QDoubleSpinBox()
        self.fp_threshold_i.setRange(0.1, 100.0)
        self.fp_threshold_i.setDecimals(2)
        self.fp_threshold_i.setValue(self.iraf_params.threshold_i)
        layout.addRow("threshold (i band):", self.fp_threshold_i)

        self.fp_threshold_default = QDoubleSpinBox()
        self.fp_threshold_default.setRange(0.1, 100.0)
        self.fp_threshold_default.setDecimals(2)
        self.fp_threshold_default.setValue(self.iraf_params.threshold_default)
        layout.addRow("threshold (default):", self.fp_threshold_default)

        # --- Filter-specific Sharplo ---
        layout.addRow(QLabel(""))
        layout.addRow(QLabel("--- Filter-specific Sharplo ---"))

        self.fp_sharplo_g = QDoubleSpinBox()
        self.fp_sharplo_g.setRange(-10.0, 10.0)
        self.fp_sharplo_g.setDecimals(2)
        self.fp_sharplo_g.setValue(self.iraf_params.sharplo_g)
        layout.addRow("sharplo (g band):", self.fp_sharplo_g)

        self.fp_sharplo_r = QDoubleSpinBox()
        self.fp_sharplo_r.setRange(-10.0, 10.0)
        self.fp_sharplo_r.setDecimals(2)
        self.fp_sharplo_r.setValue(self.iraf_params.sharplo_r)
        layout.addRow("sharplo (r band):", self.fp_sharplo_r)

        self.fp_sharplo_i = QDoubleSpinBox()
        self.fp_sharplo_i.setRange(-10.0, 10.0)
        self.fp_sharplo_i.setDecimals(2)
        self.fp_sharplo_i.setValue(self.iraf_params.sharplo_i)
        layout.addRow("sharplo (i band):", self.fp_sharplo_i)

        self.fp_sharplo_default = QDoubleSpinBox()
        self.fp_sharplo_default.setRange(-10.0, 10.0)
        self.fp_sharplo_default.setDecimals(2)
        self.fp_sharplo_default.setValue(self.iraf_params.sharplo_default)
        layout.addRow("sharplo (default):", self.fp_sharplo_default)

        # --- Filter-specific Datamin ---
        layout.addRow(QLabel(""))
        layout.addRow(QLabel("--- Filter-specific Datamin ---"))

        self.fp_datamin_g = QDoubleSpinBox()
        self.fp_datamin_g.setRange(-100000, 100000)
        self.fp_datamin_g.setDecimals(1)
        self.fp_datamin_g.setValue(self.iraf_params.datamin_g)
        layout.addRow("datamin (g band):", self.fp_datamin_g)

        self.fp_datamin_r = QDoubleSpinBox()
        self.fp_datamin_r.setRange(-100000, 100000)
        self.fp_datamin_r.setDecimals(1)
        self.fp_datamin_r.setValue(self.iraf_params.datamin_r)
        layout.addRow("datamin (r band):", self.fp_datamin_r)

        self.fp_datamin_i = QDoubleSpinBox()
        self.fp_datamin_i.setRange(-100000, 100000)
        self.fp_datamin_i.setDecimals(1)
        self.fp_datamin_i.setValue(self.iraf_params.datamin_i)
        layout.addRow("datamin (i band):", self.fp_datamin_i)

        self.fp_datamin_default = QDoubleSpinBox()
        self.fp_datamin_default.setRange(-100000, 100000)
        self.fp_datamin_default.setDecimals(1)
        self.fp_datamin_default.setValue(self.iraf_params.datamin_default)
        layout.addRow("datamin (default):", self.fp_datamin_default)

        # --- Other FINDPARS ---
        layout.addRow(QLabel(""))
        layout.addRow(QLabel("--- Other FINDPARS ---"))

        # Nsigma
        self.fp_nsigma = QDoubleSpinBox()
        self.fp_nsigma.setRange(0.1, 10.0)
        self.fp_nsigma.setDecimals(2)
        self.fp_nsigma.setValue(self.iraf_params.nsigma)
        layout.addRow("nsigma:", self.fp_nsigma)

        # Ratio
        self.fp_ratio = QDoubleSpinBox()
        self.fp_ratio.setRange(0.0, 1.0)
        self.fp_ratio.setDecimals(2)
        self.fp_ratio.setValue(self.iraf_params.ratio)
        layout.addRow("ratio (minor/major):", self.fp_ratio)

        # Theta
        self.fp_theta = QDoubleSpinBox()
        self.fp_theta.setRange(-180.0, 180.0)
        self.fp_theta.setDecimals(1)
        self.fp_theta.setValue(self.iraf_params.theta)
        layout.addRow("theta (degrees):", self.fp_theta)

        # Sharphi
        self.fp_sharphi = QDoubleSpinBox()
        self.fp_sharphi.setRange(-10.0, 10.0)
        self.fp_sharphi.setDecimals(2)
        self.fp_sharphi.setValue(self.iraf_params.sharphi)
        layout.addRow("sharphi:", self.fp_sharphi)

        # Roundlo
        self.fp_roundlo = QDoubleSpinBox()
        self.fp_roundlo.setRange(-10.0, 10.0)
        self.fp_roundlo.setDecimals(2)
        self.fp_roundlo.setValue(self.iraf_params.roundlo)
        layout.addRow("roundlo:", self.fp_roundlo)

        # Roundhi
        self.fp_roundhi = QDoubleSpinBox()
        self.fp_roundhi.setRange(-10.0, 10.0)
        self.fp_roundhi.setDecimals(2)
        self.fp_roundhi.setValue(self.iraf_params.roundhi)
        layout.addRow("roundhi:", self.fp_roundhi)

        scroll.setWidget(panel)
        return scroll

    def _create_centerpars_panel(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        panel = QWidget()
        layout = QFormLayout(panel)

        # Centering algorithm
        self.cp_calgorithm = QComboBox()
        self.cp_calgorithm.addItems(["none", "centroid", "gauss", "ofilter"])
        self.cp_calgorithm.setCurrentText(self.iraf_params.calgorithm)
        layout.addRow("calgorithm:", self.cp_calgorithm)

        # Cbox multiplier (cbox = FWHM * cbox_mult)
        layout.addRow(QLabel("--- FWHM Multiplier ---"))
        self.cp_cbox_mult = QDoubleSpinBox()
        self.cp_cbox_mult.setRange(0.5, 10.0)
        self.cp_cbox_mult.setDecimals(1)
        self.cp_cbox_mult.setValue(self.iraf_params.cbox_mult)
        layout.addRow("cbox_mult (cbox = FWHM ×):", self.cp_cbox_mult)
        layout.addRow(QLabel(""))

        # Cthreshold
        self.cp_cthreshold = QDoubleSpinBox()
        self.cp_cthreshold.setRange(0.0, 100.0)
        self.cp_cthreshold.setDecimals(2)
        self.cp_cthreshold.setValue(self.iraf_params.cthreshold)
        layout.addRow("cthreshold:", self.cp_cthreshold)

        # Min SNR
        self.cp_minsnratio = QDoubleSpinBox()
        self.cp_minsnratio.setRange(0.0, 100.0)
        self.cp_minsnratio.setDecimals(2)
        self.cp_minsnratio.setValue(self.iraf_params.minsnratio)
        layout.addRow("minsnratio:", self.cp_minsnratio)

        # Max iterations
        self.cp_cmaxiter = QSpinBox()
        self.cp_cmaxiter.setRange(1, 100)
        self.cp_cmaxiter.setValue(self.iraf_params.cmaxiter)
        layout.addRow("cmaxiter:", self.cp_cmaxiter)

        # Max shift
        self.cp_maxshift = QDoubleSpinBox()
        self.cp_maxshift.setRange(0.0, 100.0)
        self.cp_maxshift.setDecimals(2)
        self.cp_maxshift.setValue(self.iraf_params.maxshift)
        layout.addRow("maxshift:", self.cp_maxshift)

        # Clean
        self.cp_clean = QCheckBox()
        self.cp_clean.setChecked(self.iraf_params.clean)
        layout.addRow("clean:", self.cp_clean)

        # Rclean
        self.cp_rclean = QDoubleSpinBox()
        self.cp_rclean.setRange(0.0, 100.0)
        self.cp_rclean.setDecimals(2)
        self.cp_rclean.setValue(self.iraf_params.rclean)
        layout.addRow("rclean:", self.cp_rclean)

        # Rclip
        self.cp_rclip = QDoubleSpinBox()
        self.cp_rclip.setRange(0.0, 100.0)
        self.cp_rclip.setDecimals(2)
        self.cp_rclip.setValue(self.iraf_params.rclip)
        layout.addRow("rclip:", self.cp_rclip)

        # Kclean
        self.cp_kclean = QDoubleSpinBox()
        self.cp_kclean.setRange(0.0, 100.0)
        self.cp_kclean.setDecimals(2)
        self.cp_kclean.setValue(self.iraf_params.kclean)
        layout.addRow("kclean:", self.cp_kclean)

        scroll.setWidget(panel)
        return scroll

    def _create_fitskypars_panel(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        panel = QWidget()
        layout = QFormLayout(panel)

        # Sky algorithm
        self.sp_salgorithm = QComboBox()
        self.sp_salgorithm.addItems([
            "constant", "file", "mean", "median", "mode",
            "centroid", "gauss", "ofilter", "crosscor", "histplot"
        ])
        self.sp_salgorithm.setCurrentText(self.iraf_params.salgorithm)
        layout.addRow("salgorithm:", self.sp_salgorithm)

        # --- FWHM Multipliers ---
        layout.addRow(QLabel("--- FWHM Multipliers ---"))

        # Annulus multiplier (annulus = FWHM * annulus_mult)
        self.sp_annulus_mult = QDoubleSpinBox()
        self.sp_annulus_mult.setRange(1.0, 20.0)
        self.sp_annulus_mult.setDecimals(1)
        self.sp_annulus_mult.setValue(self.iraf_params.annulus_mult)
        layout.addRow("annulus_mult (ann = FWHM ×):", self.sp_annulus_mult)

        # Dannulus multiplier (dannulus = FWHM * dannulus_mult)
        self.sp_dannulus_mult = QDoubleSpinBox()
        self.sp_dannulus_mult.setRange(0.5, 10.0)
        self.sp_dannulus_mult.setDecimals(1)
        self.sp_dannulus_mult.setValue(self.iraf_params.dannulus_mult)
        layout.addRow("dannulus_mult (width = FWHM ×):", self.sp_dannulus_mult)

        layout.addRow(QLabel(""))

        # Skyvalue
        self.sp_skyvalue = QDoubleSpinBox()
        self.sp_skyvalue.setRange(-100000, 100000)
        self.sp_skyvalue.setDecimals(2)
        self.sp_skyvalue.setValue(self.iraf_params.skyvalue)
        layout.addRow("skyvalue:", self.sp_skyvalue)

        # Smaxiter
        self.sp_smaxiter = QSpinBox()
        self.sp_smaxiter.setRange(1, 100)
        self.sp_smaxiter.setValue(self.iraf_params.smaxiter)
        layout.addRow("smaxiter:", self.sp_smaxiter)

        # Sloclip
        self.sp_sloclip = QDoubleSpinBox()
        self.sp_sloclip.setRange(0.0, 100.0)
        self.sp_sloclip.setDecimals(2)
        self.sp_sloclip.setValue(self.iraf_params.sloclip)
        layout.addRow("sloclip:", self.sp_sloclip)

        # Shiclip
        self.sp_shiclip = QDoubleSpinBox()
        self.sp_shiclip.setRange(0.0, 100.0)
        self.sp_shiclip.setDecimals(2)
        self.sp_shiclip.setValue(self.iraf_params.shiclip)
        layout.addRow("shiclip:", self.sp_shiclip)

        # Snreject
        self.sp_snreject = QSpinBox()
        self.sp_snreject.setRange(0, 1000)
        self.sp_snreject.setValue(self.iraf_params.snreject)
        layout.addRow("snreject:", self.sp_snreject)

        # Sloreject
        self.sp_sloreject = QDoubleSpinBox()
        self.sp_sloreject.setRange(0.0, 100.0)
        self.sp_sloreject.setDecimals(2)
        self.sp_sloreject.setValue(self.iraf_params.sloreject)
        layout.addRow("sloreject (sigma):", self.sp_sloreject)

        # Shireject
        self.sp_shireject = QDoubleSpinBox()
        self.sp_shireject.setRange(0.0, 100.0)
        self.sp_shireject.setDecimals(2)
        self.sp_shireject.setValue(self.iraf_params.shireject)
        layout.addRow("shireject (sigma):", self.sp_shireject)

        # Khist
        self.sp_khist = QDoubleSpinBox()
        self.sp_khist.setRange(0.0, 100.0)
        self.sp_khist.setDecimals(2)
        self.sp_khist.setValue(self.iraf_params.khist)
        layout.addRow("khist:", self.sp_khist)

        # Binsize
        self.sp_binsize = QDoubleSpinBox()
        self.sp_binsize.setRange(0.001, 10.0)
        self.sp_binsize.setDecimals(3)
        self.sp_binsize.setValue(self.iraf_params.binsize)
        layout.addRow("binsize:", self.sp_binsize)

        # Smooth
        self.sp_smooth = QCheckBox()
        self.sp_smooth.setChecked(self.iraf_params.smooth)
        layout.addRow("smooth:", self.sp_smooth)

        # Rgrow
        self.sp_rgrow = QDoubleSpinBox()
        self.sp_rgrow.setRange(0.0, 100.0)
        self.sp_rgrow.setDecimals(2)
        self.sp_rgrow.setValue(self.iraf_params.rgrow)
        layout.addRow("rgrow:", self.sp_rgrow)

        scroll.setWidget(panel)
        return scroll

    def _create_photpars_panel(self):
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        panel = QWidget()
        layout = QFormLayout(panel)

        # --- FWHM Multiplier ---
        layout.addRow(QLabel("--- FWHM Multiplier ---"))

        # Aperture multiplier (aperture = FWHM * aperture_mult)
        self.pp_aperture_mult = QDoubleSpinBox()
        self.pp_aperture_mult.setRange(0.5, 10.0)
        self.pp_aperture_mult.setDecimals(1)
        self.pp_aperture_mult.setValue(self.iraf_params.aperture_mult)
        layout.addRow("aperture_mult (ap = FWHM ×):", self.pp_aperture_mult)

        layout.addRow(QLabel(""))
        layout.addRow(QLabel("--- Other PHOTPARS ---"))

        # Zmag
        self.pp_zmag = QDoubleSpinBox()
        self.pp_zmag.setRange(0.0, 50.0)
        self.pp_zmag.setDecimals(2)
        self.pp_zmag.setValue(self.iraf_params.zmag)
        layout.addRow("zmag (zero point):", self.pp_zmag)

        # Mkapert
        self.pp_mkapert = QCheckBox()
        self.pp_mkapert.setChecked(self.iraf_params.mkapert)
        layout.addRow("mkapert:", self.pp_mkapert)

        # Info note
        layout.addRow(QLabel(""))
        layout.addRow(QLabel("--- Summary ---"))
        layout.addRow(QLabel("FWHM is auto-calculated per image from header"))
        layout.addRow(QLabel("(SEEING, FWHMARC) or defaults to seeing/pix_scale"))
        layout.addRow(QLabel(""))
        layout.addRow(QLabel("aperture = FWHM × aperture_mult"))
        layout.addRow(QLabel("annulus = FWHM × annulus_mult (in FITSKYPARS)"))
        layout.addRow(QLabel("dannulus = FWHM × dannulus_mult (in FITSKYPARS)"))
        layout.addRow(QLabel("cbox = FWHM × cbox_mult (in CENTERPARS)"))

        scroll.setWidget(panel)
        return scroll

    def _get_config_path(self) -> Path:
        """Get path for IRAF config file (separate from project state)."""
        return self.result_dir / "iraf_config.json"

    def _load_iraf_params_from_toml(self) -> bool:
        toml_path = Path("parameters.toml")
        if not toml_path.exists():
            return False
        try:
            data = toml_path.read_text(encoding="utf-8")
            cfg = tomllib.loads(data)
        except Exception:
            return False

        tools = cfg.get("tools", {}) if isinstance(cfg, dict) else {}
        iraf_cfg = tools.get("iraf", {}) if isinstance(tools, dict) else {}
        if not iraf_cfg:
            iraf_cfg = cfg.get("iraf", {})
        params = iraf_cfg.get("params", {})
        if isinstance(params, dict) and params:
            try:
                self.iraf_params.from_dict(params)
                self._update_ui_from_params()
                return True
            except Exception:
                return False
        return False

    def _save_iraf_params_to_toml(self, notify: bool = False):
        if tomli_w is None:
            if notify:
                QMessageBox.warning(self, "Error", "tomli_w is required to write parameters.toml")
            return
        toml_path = Path("parameters.toml")
        try:
            if toml_path.exists():
                cfg = tomllib.loads(toml_path.read_text(encoding="utf-8"))
            else:
                cfg = {}
            tools = cfg.get("tools", {})
            if not isinstance(tools, dict):
                tools = {}
            iraf_cfg = tools.get("iraf", {})
            if not isinstance(iraf_cfg, dict):
                iraf_cfg = {}
            iraf_cfg["params"] = self.iraf_params.to_dict()
            tools["iraf"] = iraf_cfg
            cfg["tools"] = tools
            toml_path.write_text(tomli_w.dumps(cfg), encoding="utf-8")
            if notify:
                QMessageBox.information(self, "Saved", f"IRAF parameters saved to:\n{toml_path}")
        except Exception as e:
            if notify:
                QMessageBox.warning(self, "Error", f"Failed to save: {e}")

    def _save_params_to_file(self):
        """Save current parameters to separate config file."""
        self._apply_params_silent()  # Apply UI values to params object first
        self._save_iraf_params_to_toml(notify=True)

    def _auto_save_params(self):
        """Auto-save current parameters without UI prompts."""
        self._apply_params_silent()
        self._save_iraf_params_to_toml(notify=False)

    def _load_params_from_file(self):
        """Load parameters from config file."""
        if self._load_iraf_params_from_toml():
            QMessageBox.information(self, "Loaded", "IRAF parameters loaded from parameters.toml")
            return
        config_path = self._get_config_path()
        if not config_path.exists():
            QMessageBox.information(self, "Not Found", "No IRAF config found")
            return
        try:
            with open(config_path, 'r', encoding='utf-8') as f:
                data = json.load(f)
            self.iraf_params.from_dict(data)
            self._update_ui_from_params()
            QMessageBox.information(self, "Loaded", f"IRAF parameters loaded from:\n{config_path}")
        except Exception as e:
            QMessageBox.warning(self, "Error", f"Failed to load: {e}")

    def _auto_load_params(self):
        """Auto-load parameters from config file if exists."""
        if self._load_iraf_params_from_toml():
            return
        config_path = self._get_config_path()
        if config_path.exists():
            try:
                with open(config_path, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                self.iraf_params.from_dict(data)
                self._update_ui_from_params()
            except Exception:
                pass  # Silently ignore errors on auto-load

    def _apply_params_silent(self):
        """Apply UI values to params object without showing message."""
        p = self.iraf_params
        # DATAPARS - Pixel scale
        p.pix_scale = self.dp_pix_scale.value()
        # DATAPARS - Filter-specific seeing
        p.seeing_g = self.dp_seeing_g.value()
        p.seeing_r = self.dp_seeing_r.value()
        p.seeing_i = self.dp_seeing_i.value()
        p.seeing_default = self.dp_seeing_default.value()
        # DATAPARS - Filter-specific sigma
        p.sigma_g = self.dp_sigma_g.value()
        p.sigma_r = self.dp_sigma_r.value()
        p.sigma_i = self.dp_sigma_i.value()
        p.sigma_default = self.dp_sigma_default.value()
        p.sigma_ref = self.dp_sigma_ref.value()
        # DATAPARS - Other
        p.scale = self.dp_scale.value()
        p.emission = self.dp_emission.isChecked()
        p.datamax = self.dp_datamax.value()
        p.noise = self.dp_noise.currentText()
        p.readnoise = self.dp_readnoise.value()
        p.epadu = self.dp_epadu.value()
        p.exposure = self.dp_exposure.text()
        p.itime = self.dp_itime.value()
        # FINDPARS - Filter-specific
        p.threshold_g = self.fp_threshold_g.value()
        p.threshold_r = self.fp_threshold_r.value()
        p.threshold_i = self.fp_threshold_i.value()
        p.threshold_default = self.fp_threshold_default.value()
        p.sharplo_g = self.fp_sharplo_g.value()
        p.sharplo_r = self.fp_sharplo_r.value()
        p.sharplo_i = self.fp_sharplo_i.value()
        p.sharplo_default = self.fp_sharplo_default.value()
        p.datamin_g = self.fp_datamin_g.value()
        p.datamin_r = self.fp_datamin_r.value()
        p.datamin_i = self.fp_datamin_i.value()
        p.datamin_default = self.fp_datamin_default.value()
        # FINDPARS - Other
        p.nsigma = self.fp_nsigma.value()
        p.ratio = self.fp_ratio.value()
        p.theta = self.fp_theta.value()
        p.sharphi = self.fp_sharphi.value()
        p.roundlo = self.fp_roundlo.value()
        p.roundhi = self.fp_roundhi.value()
        # CENTERPARS
        p.calgorithm = self.cp_calgorithm.currentText()
        p.cbox_mult = self.cp_cbox_mult.value()
        p.cthreshold = self.cp_cthreshold.value()
        p.minsnratio = self.cp_minsnratio.value()
        p.cmaxiter = self.cp_cmaxiter.value()
        p.maxshift = self.cp_maxshift.value()
        p.clean = self.cp_clean.isChecked()
        p.rclean = self.cp_rclean.value()
        p.rclip = self.cp_rclip.value()
        p.kclean = self.cp_kclean.value()
        # FITSKYPARS
        p.salgorithm = self.sp_salgorithm.currentText()
        p.annulus_mult = self.sp_annulus_mult.value()
        p.dannulus_mult = self.sp_dannulus_mult.value()
        p.skyvalue = self.sp_skyvalue.value()
        p.smaxiter = self.sp_smaxiter.value()
        p.sloclip = self.sp_sloclip.value()
        p.shiclip = self.sp_shiclip.value()
        p.snreject = self.sp_snreject.value()
        p.sloreject = self.sp_sloreject.value()
        p.shireject = self.sp_shireject.value()
        p.khist = self.sp_khist.value()
        p.binsize = self.sp_binsize.value()
        p.smooth = self.sp_smooth.isChecked()
        p.rgrow = self.sp_rgrow.value()
        # PHOTPARS
        p.aperture_mult = self.pp_aperture_mult.value()
        p.zmag = self.pp_zmag.value()
        p.mkapert = self.pp_mkapert.isChecked()

    def _load_defaults(self):
        self.iraf_params = IRAFParameters()
        self._update_ui_from_params()
        QMessageBox.information(self, "Defaults Loaded", "Default IRAF parameters restored.")

    def _update_ui_from_params(self):
        p = self.iraf_params
        # DATAPARS - Pixel scale
        self.dp_pix_scale.setValue(p.pix_scale)
        # DATAPARS - Filter-specific seeing
        self.dp_seeing_g.setValue(p.seeing_g)
        self.dp_seeing_r.setValue(p.seeing_r)
        self.dp_seeing_i.setValue(p.seeing_i)
        self.dp_seeing_default.setValue(p.seeing_default)
        # DATAPARS - Filter-specific sigma
        self.dp_sigma_g.setValue(p.sigma_g)
        self.dp_sigma_r.setValue(p.sigma_r)
        self.dp_sigma_i.setValue(p.sigma_i)
        self.dp_sigma_default.setValue(p.sigma_default)
        self.dp_sigma_ref.setValue(p.sigma_ref)
        # DATAPARS - Other
        self.dp_scale.setValue(p.scale)
        self.dp_emission.setChecked(p.emission)
        self.dp_datamax.setValue(p.datamax)
        self.dp_noise.setCurrentText(p.noise)
        self.dp_readnoise.setValue(p.readnoise)
        self.dp_epadu.setValue(p.epadu)
        self.dp_exposure.setText(p.exposure)
        self.dp_itime.setValue(p.itime)
        # FINDPARS - Filter-specific
        self.fp_threshold_g.setValue(p.threshold_g)
        self.fp_threshold_r.setValue(p.threshold_r)
        self.fp_threshold_i.setValue(p.threshold_i)
        self.fp_threshold_default.setValue(p.threshold_default)
        self.fp_sharplo_g.setValue(p.sharplo_g)
        self.fp_sharplo_r.setValue(p.sharplo_r)
        self.fp_sharplo_i.setValue(p.sharplo_i)
        self.fp_sharplo_default.setValue(p.sharplo_default)
        self.fp_datamin_g.setValue(p.datamin_g)
        self.fp_datamin_r.setValue(p.datamin_r)
        self.fp_datamin_i.setValue(p.datamin_i)
        self.fp_datamin_default.setValue(p.datamin_default)
        # FINDPARS - Other
        self.fp_nsigma.setValue(p.nsigma)
        self.fp_ratio.setValue(p.ratio)
        self.fp_theta.setValue(p.theta)
        self.fp_sharphi.setValue(p.sharphi)
        self.fp_roundlo.setValue(p.roundlo)
        self.fp_roundhi.setValue(p.roundhi)
        # CENTERPARS
        self.cp_calgorithm.setCurrentText(p.calgorithm)
        self.cp_cbox_mult.setValue(p.cbox_mult)
        self.cp_cthreshold.setValue(p.cthreshold)
        self.cp_minsnratio.setValue(p.minsnratio)
        self.cp_cmaxiter.setValue(p.cmaxiter)
        self.cp_maxshift.setValue(p.maxshift)
        self.cp_clean.setChecked(p.clean)
        self.cp_rclean.setValue(p.rclean)
        self.cp_rclip.setValue(p.rclip)
        self.cp_kclean.setValue(p.kclean)
        # FITSKYPARS
        self.sp_salgorithm.setCurrentText(p.salgorithm)
        self.sp_annulus_mult.setValue(p.annulus_mult)
        self.sp_dannulus_mult.setValue(p.dannulus_mult)
        self.sp_skyvalue.setValue(p.skyvalue)
        self.sp_smaxiter.setValue(p.smaxiter)
        self.sp_sloclip.setValue(p.sloclip)
        self.sp_shiclip.setValue(p.shiclip)
        self.sp_snreject.setValue(p.snreject)
        self.sp_sloreject.setValue(p.sloreject)
        self.sp_shireject.setValue(p.shireject)
        self.sp_khist.setValue(p.khist)
        self.sp_binsize.setValue(p.binsize)
        self.sp_smooth.setChecked(p.smooth)
        self.sp_rgrow.setValue(p.rgrow)
        # PHOTPARS
        self.pp_aperture_mult.setValue(p.aperture_mult)
        self.pp_zmag.setValue(p.zmag)
        self.pp_mkapert.setChecked(p.mkapert)

    def _apply_params(self):
        p = self.iraf_params
        # DATAPARS - Pixel scale
        p.pix_scale = self.dp_pix_scale.value()
        # DATAPARS - Filter-specific seeing
        p.seeing_g = self.dp_seeing_g.value()
        p.seeing_r = self.dp_seeing_r.value()
        p.seeing_i = self.dp_seeing_i.value()
        p.seeing_default = self.dp_seeing_default.value()
        # DATAPARS - Filter-specific sigma
        p.sigma_g = self.dp_sigma_g.value()
        p.sigma_r = self.dp_sigma_r.value()
        p.sigma_i = self.dp_sigma_i.value()
        p.sigma_default = self.dp_sigma_default.value()
        p.sigma_ref = self.dp_sigma_ref.value()
        # DATAPARS - Other
        p.scale = self.dp_scale.value()
        p.emission = self.dp_emission.isChecked()
        p.datamax = self.dp_datamax.value()
        p.noise = self.dp_noise.currentText()
        p.readnoise = self.dp_readnoise.value()
        p.epadu = self.dp_epadu.value()
        p.exposure = self.dp_exposure.text()
        p.itime = self.dp_itime.value()
        # FINDPARS - Filter-specific
        p.threshold_g = self.fp_threshold_g.value()
        p.threshold_r = self.fp_threshold_r.value()
        p.threshold_i = self.fp_threshold_i.value()
        p.threshold_default = self.fp_threshold_default.value()
        p.sharplo_g = self.fp_sharplo_g.value()
        p.sharplo_r = self.fp_sharplo_r.value()
        p.sharplo_i = self.fp_sharplo_i.value()
        p.sharplo_default = self.fp_sharplo_default.value()
        p.datamin_g = self.fp_datamin_g.value()
        p.datamin_r = self.fp_datamin_r.value()
        p.datamin_i = self.fp_datamin_i.value()
        p.datamin_default = self.fp_datamin_default.value()
        # FINDPARS - Other
        p.nsigma = self.fp_nsigma.value()
        p.ratio = self.fp_ratio.value()
        p.theta = self.fp_theta.value()
        p.sharphi = self.fp_sharphi.value()
        p.roundlo = self.fp_roundlo.value()
        p.roundhi = self.fp_roundhi.value()
        # CENTERPARS
        p.calgorithm = self.cp_calgorithm.currentText()
        p.cbox_mult = self.cp_cbox_mult.value()
        p.cthreshold = self.cp_cthreshold.value()
        p.minsnratio = self.cp_minsnratio.value()
        p.cmaxiter = self.cp_cmaxiter.value()
        p.maxshift = self.cp_maxshift.value()
        p.clean = self.cp_clean.isChecked()
        p.rclean = self.cp_rclean.value()
        p.rclip = self.cp_rclip.value()
        p.kclean = self.cp_kclean.value()
        # FITSKYPARS
        p.salgorithm = self.sp_salgorithm.currentText()
        p.annulus_mult = self.sp_annulus_mult.value()
        p.dannulus_mult = self.sp_dannulus_mult.value()
        p.skyvalue = self.sp_skyvalue.value()
        p.smaxiter = self.sp_smaxiter.value()
        p.sloclip = self.sp_sloclip.value()
        p.shiclip = self.sp_shiclip.value()
        p.snreject = self.sp_snreject.value()
        p.sloreject = self.sp_sloreject.value()
        p.shireject = self.sp_shireject.value()
        p.khist = self.sp_khist.value()
        p.binsize = self.sp_binsize.value()
        p.smooth = self.sp_smooth.isChecked()
        p.rgrow = self.sp_rgrow.value()
        # PHOTPARS
        p.aperture_mult = self.pp_aperture_mult.value()
        p.zmag = self.pp_zmag.value()
        p.mkapert = self.pp_mkapert.isChecked()

        self._auto_save_params()
        QMessageBox.information(self, "Applied", "Parameters applied to current session.")

    # ========================================================================
    # Tab 3: Comparison
    # ========================================================================
    def _create_comparison_tab(self):
        tab = QWidget()
        layout = QVBoxLayout(tab)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(6)

        btn_layout = QHBoxLayout()
        btn_params = create_parameter_button("IRAF Parameters")
        btn_params.clicked.connect(self.show_parameters_tab)
        btn_layout.addWidget(btn_params)

        self.cmp_run_btn = QPushButton("Run Comparison")
        self.cmp_run_btn.clicked.connect(self.run_comparison)
        btn_layout.addWidget(self.cmp_run_btn)

        self.cmp_export_btn = QPushButton("Export CSV")
        self.cmp_export_btn.clicked.connect(self.export_comparison)
        btn_layout.addWidget(self.cmp_export_btn)

        self.cmp_log_btn = create_cache_action_button("Log")
        self.cmp_log_btn.clicked.connect(self.show_comparison_log)
        btn_layout.addWidget(self.cmp_log_btn)

        btn_layout.addStretch()
        layout.addLayout(btn_layout)

        paths_group = QGroupBox("Project Paths")
        paths_group.setMaximumHeight(82)
        paths_layout = QVBoxLayout(paths_group)
        paths_layout.setContentsMargins(8, 4, 8, 4)
        self.cmp_path_summary = QLabel()
        self.cmp_path_summary.setTextInteractionFlags(Qt.TextSelectableByMouse)
        self.cmp_path_summary.setStyleSheet(mono_note_style())
        paths_layout.addWidget(self.cmp_path_summary)
        layout.addWidget(paths_group)
        self._refresh_path_summaries()

        # Summary
        self.cmp_summary = QLabel("No comparison run yet.")
        _f = self.cmp_summary.font(); _f.setBold(True)
        self.cmp_summary.setFont(_f)
        self.cmp_summary.setMaximumHeight(30)
        layout.addWidget(self.cmp_summary)

        # Tabs: results + log
        tabs = QTabWidget()
        tabs.setMinimumHeight(430)
        self.cmp_tabs = tabs

        result_tab = QWidget()
        result_layout = QVBoxLayout(result_tab)
        result_layout.setContentsMargins(0, 0, 0, 0)

        # Splitter: table + plot
        splitter = QSplitter(Qt.Horizontal)

        # Table
        self.cmp_table = QTableWidget(0, 13)
        self.cmp_table.setHorizontalHeaderLabels([
            "Frame", "Matched", "dmag_med", "resid_std", "dx_med", "dy_med",
            "dist_med", "dist_p95", "frac<=tol", "shift_x", "shift_y",
            "N_iraf", "N_apex"
        ])
        self.cmp_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        self.cmp_table.horizontalHeader().setStretchLastSection(True)
        self.cmp_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.cmp_table.setSelectionMode(QTableWidget.SingleSelection)
        self.cmp_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.cmp_table.itemSelectionChanged.connect(self._plot_comparison)
        self.cmp_table.setMinimumWidth(360)
        self.cmp_table.setMaximumWidth(540)
        splitter.addWidget(self.cmp_table)

        # Plot
        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)
        plot_layout.setContentsMargins(2, 2, 2, 2)
        plot_layout.setSpacing(2)
        self.cmp_fig = Figure(figsize=(9, 6.5), tight_layout=True)
        self.cmp_canvas = FigureCanvas(self.cmp_fig)
        self.cmp_canvas.setMinimumSize(620, 390)
        self.cmp_toolbar = NavigationToolbar(self.cmp_canvas, self)
        plot_layout.addWidget(self.cmp_toolbar)
        plot_layout.addWidget(self.cmp_canvas, stretch=1)
        splitter.addWidget(plot_widget)

        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 4)
        splitter.setSizes([420, 850])
        splitter.setChildrenCollapsible(False)

        result_layout.addWidget(splitter)
        tabs.addTab(result_tab, "Results")

        log_tab = QWidget()
        log_layout = QVBoxLayout(log_tab)
        self.cmp_log_text = QTextEdit()
        self.cmp_log_text.setReadOnly(True)
        self.cmp_log_text.setObjectName("Log")
        log_layout.addWidget(self.cmp_log_text)
        tabs.addTab(log_tab, "Log")

        overlay_widget = self._create_overlay_widget()
        tabs.addTab(overlay_widget, "Overlay")

        layout.addWidget(tabs, stretch=1)

        return tab

    def _create_overlay_widget(self):
        widget = QWidget()
        layout = QVBoxLayout(widget)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        select_group = QGroupBox("Frame Selection")
        select_layout = QHBoxLayout(select_group)
        select_layout.setContentsMargins(8, 4, 8, 4)
        select_layout.setSpacing(6)
        select_layout.addWidget(QLabel("Index:"))
        self.overlay_index_spin = QSpinBox()
        self.overlay_index_spin.setRange(0, 0)
        self.overlay_index_spin.valueChanged.connect(self._overlay_on_index_changed)
        select_layout.addWidget(self.overlay_index_spin)

        select_layout.addWidget(QLabel("File:"))
        self.overlay_file_combo = QComboBox()
        self.overlay_file_combo.currentIndexChanged.connect(self._overlay_on_file_changed)
        select_layout.addWidget(self.overlay_file_combo, stretch=1)

        btn_reload = QPushButton("Reload")
        btn_reload.clicked.connect(self._overlay_reload)
        select_layout.addWidget(btn_reload)
        select_group.setMaximumHeight(58)
        layout.addWidget(select_group)

        control_group = QGroupBox("Display Controls")
        control_layout = QHBoxLayout(control_group)
        control_layout.setContentsMargins(8, 4, 8, 4)
        control_layout.setSpacing(8)
        btn_reset_zoom = QPushButton("Reset Zoom")
        btn_reset_zoom.clicked.connect(self._overlay_reset_zoom)
        control_layout.addWidget(btn_reset_zoom)

        self.overlay_show_apex = QCheckBox("APEX TSV")
        self.overlay_show_apex.setChecked(True)
        self.overlay_show_apex.stateChanged.connect(self._overlay_redisplay)
        control_layout.addWidget(self.overlay_show_apex)

        self.overlay_show_iraf = QCheckBox("IRAF daofind")
        self.overlay_show_iraf.setChecked(True)
        self.overlay_show_iraf.stateChanged.connect(self._overlay_redisplay)
        control_layout.addWidget(self.overlay_show_iraf)

        control_layout.addStretch()
        control_group.setMaximumHeight(58)
        layout.addWidget(control_group)

        self.overlay_viewer = FITSViewerWidget(self)
        self.overlay_viewer.setMinimumHeight(360)
        layout.addWidget(self.overlay_viewer, stretch=1)

        self.overlay_status = QLabel("No frame loaded.")
        self.overlay_status.setProperty("role", "caption")
        self.overlay_status.setMaximumHeight(28)
        layout.addWidget(self.overlay_status)

        sc_prev = QShortcut(QKeySequence(Qt.Key_BracketLeft), widget)
        sc_prev.setContext(Qt.WidgetWithChildrenShortcut)
        sc_prev.activated.connect(lambda: self._overlay_navigate(-1))

        sc_next = QShortcut(QKeySequence(Qt.Key_BracketRight), widget)
        sc_next.setContext(Qt.WidgetWithChildrenShortcut)
        sc_next.activated.connect(lambda: self._overlay_navigate(1))

        sc_filter = QShortcut(QKeySequence("."), widget)
        sc_filter.setContext(Qt.WidgetWithChildrenShortcut)
        sc_filter.activated.connect(self._overlay_cycle_filter)

        self._overlay_reload()
        return widget

    def _build_iraf_param_defaults(self) -> dict:
        p = self.iraf_params
        return {
            "scale": p.scale,
            "emission": p.emission,
            "datamax": p.datamax,
            "noise": p.noise,
            "readnoise": p.readnoise,
            "epadu": p.epadu,
            "exposure": p.exposure,
            "itime": p.itime,
            "seeing_g": p.seeing_g,
            "seeing_r": p.seeing_r,
            "seeing_i": p.seeing_i,
            "seeing_default": p.seeing_default,
            "sigma_g": p.sigma_g,
            "sigma_r": p.sigma_r,
            "sigma_i": p.sigma_i,
            "sigma_default": p.sigma_default,
            "threshold_g": p.threshold_g,
            "threshold_r": p.threshold_r,
            "threshold_i": p.threshold_i,
            "threshold_default": p.threshold_default,
            "nsigma": p.nsigma,
            "ratio": p.ratio,
            "theta": p.theta,
            "sharplo_g": p.sharplo_g,
            "sharplo_r": p.sharplo_r,
            "sharplo_i": p.sharplo_i,
            "sharplo_default": p.sharplo_default,
            "sharphi": p.sharphi,
            "roundlo": p.roundlo,
            "roundhi": p.roundhi,
            "datamin_g": p.datamin_g,
            "datamin_r": p.datamin_r,
            "datamin_i": p.datamin_i,
            "datamin_default": p.datamin_default,
            "calgorithm": p.calgorithm,
            "cbox_mult": p.cbox_mult,
            "cthreshold": p.cthreshold,
            "minsnratio": p.minsnratio,
            "cmaxiter": p.cmaxiter,
            "maxshift": p.maxshift,
            "clean": p.clean,
            "rclean": p.rclean,
            "rclip": p.rclip,
            "kclean": p.kclean,
            "salgorithm": p.salgorithm,
            "annulus_mult": p.annulus_mult,
            "dannulus_mult": p.dannulus_mult,
            "skyvalue": p.skyvalue,
            "smaxiter": p.smaxiter,
            "sloclip": p.sloclip,
            "shiclip": p.shiclip,
            "snreject": p.snreject,
            "sloreject": p.sloreject,
            "shireject": p.shireject,
            "khist": p.khist,
            "binsize": p.binsize,
            "smooth": p.smooth,
            "rgrow": p.rgrow,
            "aperture_mult": p.aperture_mult,
            "zmag": p.zmag,
            "mkapert": p.mkapert,
            "pix_scale": p.pix_scale,
            "sigma_ref": p.sigma_ref,
        }

    def _load_iraf_toml_config(self) -> tuple[dict, dict]:
        toml_path = Path("parameters.toml")
        if not toml_path.exists():
            return {}, {}
        try:
            data = toml_path.read_text(encoding="utf-8")
            cfg = tomllib.loads(data)
        except Exception:
            return {}, {}

        tools = cfg.get("tools", {}) if isinstance(cfg, dict) else {}
        iraf_cfg = tools.get("iraf", {}) if isinstance(tools, dict) else {}
        if not iraf_cfg:
            iraf_cfg = cfg.get("iraf", {})
        raw_filters = iraf_cfg.get("filters", {})
        raw_aliases = iraf_cfg.get("filter_aliases", {})

        key_map = {
            "seeing_arcsec": "seeing",
            "seeing": "seeing",
            "sigma": "sigma",
            "threshold": "threshold",
            "sharplo": "sharplo",
            "datamin": "datamin",
            "nsigma": "nsigma",
            "ratio": "ratio",
            "theta": "theta",
            "sharphi": "sharphi",
            "roundlo": "roundlo",
            "roundhi": "roundhi",
            "calgorithm": "calgorithm",
            "cbox_mult": "cbox_mult",
            "cthreshold": "cthreshold",
            "minsnratio": "minsnratio",
            "cmaxiter": "cmaxiter",
            "maxshift": "maxshift",
            "clean": "clean",
            "rclean": "rclean",
            "rclip": "rclip",
            "kclean": "kclean",
            "salgorithm": "salgorithm",
            "annulus_mult": "annulus_mult",
            "dannulus_mult": "dannulus_mult",
            "skyvalue": "skyvalue",
            "smaxiter": "smaxiter",
            "sloclip": "sloclip",
            "shiclip": "shiclip",
            "snreject": "snreject",
            "sloreject": "sloreject",
            "shireject": "shireject",
            "khist": "khist",
            "binsize": "binsize",
            "smooth": "smooth",
            "rgrow": "rgrow",
            "aperture_mult": "aperture_mult",
            "zmag": "zmag",
            "mkapert": "mkapert",
            "scale": "scale",
            "emission": "emission",
            "datamax": "datamax",
            "noise": "noise",
            "readnoise": "readnoise",
            "epadu": "epadu",
            "exposure": "exposure",
            "itime": "itime",
            "sigma_ref": "sigma_ref",
        }

        filter_params = {}
        if isinstance(raw_filters, dict):
            for fkey, fvals in raw_filters.items():
                if not isinstance(fvals, dict):
                    continue
                key = normalize_filter_key(fkey)
                params = {}
                for pkey, pval in fvals.items():
                    pname = str(pkey).strip().lower().replace("-", "_").replace(" ", "_")
                    pname = key_map.get(pname, pname)
                    params[pname] = pval
                filter_params[key] = params

        filter_aliases = {}
        if isinstance(raw_aliases, dict):
            for akey, aval in raw_aliases.items():
                if aval is None:
                    continue
                filter_aliases[str(akey).strip().lower()] = normalize_filter_key(aval)

        return filter_params, filter_aliases

    def _overlay_image_dir(self) -> Path:
        return Path(self.cmp_image_edit.text())

    def _overlay_reload(self):
        img_dir = self._overlay_image_dir()
        if not img_dir.exists():
            self.overlay_file_list = []
            self.overlay_keys = []
            self.overlay_image_map = {}
            self.overlay_key_to_index = {}
            self.overlay_file_combo.clear()
            self.overlay_index_spin.setRange(0, 0)
            self.overlay_image_data = None
            self._overlay_render_empty(f"Image dir not found: {img_dir}")
            self._cmp_log(f"Overlay: image dir not found: {img_dir}")
            return

        files = sorted([p.name for p in img_dir.glob("*.fit*")])
        self.overlay_file_list = list(files)
        self.overlay_keys = []
        self.overlay_image_map = {}
        self.overlay_key_to_index = {}
        self.overlay_filter_cache = {}
        self._overlay_normalized_cache = None
        self.overlay_last_image_dir = img_dir

        for idx, fname in enumerate(self.overlay_file_list):
            key = _normalize_frame_key(Path(fname).stem)
            self.overlay_keys.append(key)
            if key not in self.overlay_image_map:
                self.overlay_image_map[key] = img_dir / fname
            if key not in self.overlay_key_to_index:
                self.overlay_key_to_index[key] = idx

        self.overlay_file_combo.blockSignals(True)
        self.overlay_file_combo.clear()
        self.overlay_file_combo.addItems(self.overlay_file_list)
        self.overlay_file_combo.blockSignals(False)

        if self.overlay_file_list:
            self.overlay_index_spin.setRange(0, max(0, len(self.overlay_file_list) - 1))
            idx = min(self.overlay_current_index, len(self.overlay_file_list) - 1)
            self.overlay_index_spin.setValue(idx)
        else:
            self.overlay_index_spin.setRange(0, 0)
            self.overlay_image_data = None
            self._overlay_render_empty("No FITS files found.")

        self._overlay_refresh_maps()
        self._cmp_log(
            f"Overlay: frames={len(self.overlay_file_list)} | image_dir={img_dir}"
        )

    def _overlay_refresh_maps(self):
        apex_dir = Path(self.cmp_apex_edit.text())
        iraf_dir = Path(self.cmp_iraf_edit.text())

        self.overlay_apex_map = {}
        self.overlay_iraf_map = {}

        if apex_dir.exists():
            for p in _iter_apex_photometry_files(apex_dir):
                key = _normalize_frame_key(p.stem)
                self.overlay_apex_map.setdefault(key, p)

        if iraf_dir.exists():
            for p in iraf_dir.rglob("*.coo"):
                key = _normalize_frame_key(p.stem)
                self.overlay_iraf_map.setdefault(key, p)
            if not self.overlay_iraf_map:
                for p in iraf_dir.rglob("*.txt"):
                    key = _normalize_frame_key(p.stem)
                    self.overlay_iraf_map.setdefault(key, p)

        self._cmp_log(
            f"Overlay: apex={len(self.overlay_apex_map)} | iraf={len(self.overlay_iraf_map)}"
        )

    def _overlay_on_file_changed(self, index):
        if index < 0 or index >= len(self.overlay_file_list):
            return
        self.overlay_current_index = index
        self.overlay_index_spin.blockSignals(True)
        self.overlay_index_spin.setValue(index)
        self.overlay_index_spin.blockSignals(False)
        self._overlay_load_current()

    def _overlay_on_index_changed(self, index):
        if index < 0 or index >= len(self.overlay_file_list):
            return
        self.overlay_current_index = index
        self.overlay_file_combo.blockSignals(True)
        self.overlay_file_combo.setCurrentIndex(index)
        self.overlay_file_combo.blockSignals(False)
        self._overlay_load_current()

    def _overlay_load_current(self):
        if not self.overlay_file_list:
            return
        img_dir = self._overlay_image_dir()
        fname = self.overlay_file_list[self.overlay_current_index]
        fpath = img_dir / fname
        self.overlay_xlim_original = None
        self.overlay_ylim_original = None
        if not fpath.exists():
            self._overlay_render_empty(f"Missing file: {fpath}")
            self._cmp_log(f"Overlay: missing file {fpath}")
            return
        try:
            with fits.open(fpath) as hdul:
                data = hdul[0].data
                header = hdul[0].header
        except Exception as e:
            self._overlay_render_empty(f"Failed to load: {fname}")
            self._cmp_log(f"Overlay: failed to load {fname}: {e}")
            return

        if data is None:
            self._overlay_render_empty(f"No image data: {fname}")
            self._cmp_log(f"Overlay: no data in {fname}")
            return

        if data.ndim > 2:
            data = data[0]

        self.overlay_image_data = data.astype(np.float32)
        self.overlay_header = header
        self._overlay_normalized_cache = None
        self._overlay_imshow_obj = None
        self._overlay_reset_stretch_plot_values()

        self._overlay_display()
        if self.overlay_viewer is not None:
            self.overlay_viewer.setFocus()

    def _overlay_display(self):
        if self.overlay_image_data is None:
            return
        if self.overlay_viewer is None:
            return

        fname = self.overlay_file_list[self.overlay_current_index]
        key = _normalize_frame_key(Path(fname).stem)

        apex_xy = self._overlay_get_apex_xy(key)
        iraf_xy = self._overlay_get_iraf_xy(key)

        apex_n = 0
        iraf_n = 0
        markers: list[OverlayMarker] = []

        if self.overlay_show_apex.isChecked() and apex_xy.size:
            apex_n = apex_xy.shape[0]
            markers.extend(
                OverlayMarker(
                    col=float(x), row=float(y), radius=5.0,
                    color=QColor(0, 200, 83, 210), line_width=1.0,
                )
                for x, y in apex_xy
            )

        if self.overlay_show_iraf.isChecked() and iraf_xy.size:
            iraf_n = iraf_xy.shape[0]
            markers.extend(
                OverlayMarker(
                    col=float(x), row=float(y), radius=4.0,
                    color=QColor(255, 82, 82, 220), rejected=True,
                    line_width=1.2,
                )
                for x, y in iraf_xy
            )

        self.overlay_viewer.set_data_auto_stf(self.overlay_image_data)
        self.overlay_viewer.set_overlay_markers(markers)
        self.overlay_viewer.fit_in_view()

        filt = self._overlay_get_filter(fname)
        self.overlay_status.setText(
            f"Frame: {fname} | Filter: {filt} | APEX: {apex_n} | IRAF: {iraf_n}"
        )

        self._cmp_log(
            f"Overlay: {fname} | filter={filt} | apex={apex_n} | iraf={iraf_n}"
        )

    def _overlay_render_empty(self, message: str):
        if self.overlay_viewer is not None:
            self.overlay_viewer.clear_overlay_markers()
            self.overlay_viewer.set_data_auto_stf(np.zeros((2, 2), dtype=np.float32))
            self.overlay_viewer.fit_in_view()
        self._overlay_imshow_obj = None
        self.overlay_status.setText(message)
        self.overlay_xlim_original = None
        self.overlay_ylim_original = None

    def _overlay_get_apex_xy(self, key: str) -> np.ndarray:
        path = self.overlay_apex_map.get(key)
        if path is None:
            return np.zeros((0, 2), dtype=float)
        try:
            df = _read_apex_tsv(path)
        except Exception:
            return np.zeros((0, 2), dtype=float)
        x_col = _pick_first(df.columns, ["xcenter", "x", "x_init"])
        y_col = _pick_first(df.columns, ["ycenter", "y", "y_init"])
        if x_col is None or y_col is None:
            return np.zeros((0, 2), dtype=float)
        xy = df[[x_col, y_col]].to_numpy(float)
        xy = xy[np.isfinite(xy).all(axis=1)]
        return xy

    def _overlay_shift_for_frame(self, key: str) -> tuple[float, float]:
        for row in self.frame_rows:
            if row.get("frame") == key:
                sx = row.get("best_shift_x", np.nan)
                sy = row.get("best_shift_y", np.nan)
                if np.isfinite(sx):
                    sx = sx - BASE_IRAF_SHIFT
                else:
                    sx = 0.0
                if np.isfinite(sy):
                    sy = sy - BASE_IRAF_SHIFT
                else:
                    sy = 0.0
                return float(sx), float(sy)
        return 0.0, 0.0

    def _filter_for_frame_key(self, frame_key: str) -> str | None:
        path = self.overlay_image_map.get(frame_key)
        if path is None:
            return None
        fname = path.name
        return self._overlay_get_filter(fname)

    def _overlay_get_iraf_xy(self, key: str) -> np.ndarray:
        path = self.overlay_iraf_map.get(key)
        if path is None:
            return np.zeros((0, 2), dtype=float)
        try:
            if path.suffix.lower() == ".coo":
                df = _read_iraf_coo(path)
            else:
                df = _read_iraf_txt(path)
        except Exception:
            return np.zeros((0, 2), dtype=float)
        if "x" not in df.columns or "y" not in df.columns:
            return np.zeros((0, 2), dtype=float)
        xy = df[["x", "y"]].to_numpy(float)
        xy = xy[np.isfinite(xy).all(axis=1)]
        shift_x, shift_y = self._overlay_shift_for_frame(key)
        if shift_x or shift_y:
            xy[:, 0] += shift_x
            xy[:, 1] += shift_y
        return xy

    def _overlay_get_filter(self, fname: str) -> str:
        if fname in self.overlay_filter_cache:
            return self.overlay_filter_cache[fname]
        img_dir = self._overlay_image_dir()
        fpath = img_dir / fname
        filt = ""
        try:
            header = fits.getheader(fpath)
            filt = normalize_filter_key(header.get("FILTER", ""))
        except Exception:
            filt = ""
        if not filt:
            filt = "unknown"
        self.overlay_filter_cache[fname] = filt
        return filt

    def _overlay_navigate(self, direction: int):
        if not self.overlay_file_list:
            return
        new_index = (self.overlay_current_index + direction) % len(self.overlay_file_list)
        self.overlay_index_spin.setValue(new_index)

    def _overlay_cycle_filter(self):
        if not self.overlay_file_list:
            return
        current = self.overlay_file_list[self.overlay_current_index]
        current_filter = self._overlay_get_filter(current)
        filters = {}
        for idx, fname in enumerate(self.overlay_file_list):
            filt = self._overlay_get_filter(fname)
            filters.setdefault(filt, []).append(idx)
        if len(filters) <= 1:
            self._cmp_log("Overlay: only one filter available.")
            return
        filter_names = sorted(filters.keys())
        if current_filter in filter_names:
            cur_idx = filter_names.index(current_filter)
        else:
            cur_idx = -1
        next_filter = filter_names[(cur_idx + 1) % len(filter_names)]
        next_idx = filters[next_filter][0]
        self.overlay_index_spin.setValue(next_idx)

    def _overlay_on_stretch_changed(self, index):
        self._overlay_normalized_cache = None
        self._overlay_reset_stretch_plot_values()
        self._overlay_display()

    def _overlay_update_stretch_label(self, value):
        if hasattr(self, "overlay_stretch_value"):
            self.overlay_stretch_value.setText(str(value))

    def _overlay_update_black_label(self, value):
        if hasattr(self, "overlay_black_value"):
            self.overlay_black_value.setText(str(value))

    def _overlay_reset_stretch(self):
        self._overlay_normalized_cache = None
        self._overlay_reset_stretch_plot_values()
        if self.overlay_viewer is not None and self.overlay_image_data is not None:
            self.overlay_viewer.set_data_auto_stf(self.overlay_image_data)

    def _overlay_reset_zoom(self):
        if self.overlay_viewer is not None:
            self.overlay_viewer.fit_in_view()

    def _overlay_redisplay(self):
        self._overlay_display()

    def _overlay_open_stretch_plot(self):
        """Open stretch plot window showing histogram with draggable min/max markers"""
        if self.overlay_image_data is None:
            QMessageBox.warning(self, "Warning", "Load an image first")
            return

        if self.overlay_stretch_plot_dialog is not None and self.overlay_stretch_plot_dialog.isVisible():
            self.overlay_stretch_plot_dialog.raise_()
            self.overlay_stretch_plot_dialog.activateWindow()
            self._overlay_update_stretch_plot()
            return

        self.overlay_stretch_plot_dialog = FittedDialog(self)
        self.overlay_stretch_plot_dialog.setWindowTitle("2D Plot - Stretch Control")
        self.overlay_stretch_plot_dialog.resize(500, 250)

        layout = QVBoxLayout(self.overlay_stretch_plot_dialog)

        self.overlay_stretch_plot_info_label = QLabel("Drag min/max markers to adjust stretch")
        self.overlay_stretch_plot_info_label.setProperty("role", "info")
        layout.addWidget(self.overlay_stretch_plot_info_label)

        self.overlay_stretch_plot_fig = Figure(figsize=(6, 2.5))
        self.overlay_stretch_plot_canvas = FigureCanvas(self.overlay_stretch_plot_fig)
        self.overlay_stretch_plot_ax = self.overlay_stretch_plot_fig.add_subplot(111)
        self.overlay_stretch_plot_fig.subplots_adjust(left=0.1, right=0.95, bottom=0.15, top=0.9)

        self.overlay_stretch_plot_canvas.mpl_connect('button_press_event', self._overlay_on_stretch_plot_press)
        self.overlay_stretch_plot_canvas.mpl_connect('motion_notify_event', self._overlay_on_stretch_plot_motion)
        self.overlay_stretch_plot_canvas.mpl_connect('button_release_event', self._overlay_on_stretch_plot_release)

        layout.addWidget(tame_canvas(self.overlay_stretch_plot_canvas, min_h=140), 1)

        hint_label = QLabel("Click and drag < > markers to adjust min/max | Changes apply in real-time")
        hint_label.setProperty("role", "caption")
        layout.addWidget(hint_label)

        self.overlay_stretch_plot_dialog.show()
        self._overlay_update_stretch_plot()

    def _overlay_update_stretch_plot(self):
        """Update the stretch plot histogram and markers"""
        if self.overlay_stretch_plot_ax is None or self.overlay_image_data is None:
            return

        ax = self.overlay_stretch_plot_ax
        ax.clear()

        data = self.overlay_image_data.copy()
        finite_mask = np.isfinite(data)
        if not finite_mask.any():
            return

        flat = data[finite_mask].flatten()
        p_low, p_high = np.percentile(flat, [1, 99])
        display_data = flat[(flat >= p_low) & (flat <= p_high)]
        if len(display_data) == 0:
            display_data = flat

        self._overlay_stretch_data_range = (float(p_low), float(p_high))

        if self._overlay_stretch_vmin is None or self._overlay_stretch_vmax is None:
            _, median_val, std_val = sigma_clipped_stats(flat, sigma=3.0, maxiters=5)
            vmin = max(np.min(flat), median_val - 2.8 * std_val)
            vmax = min(np.max(flat), np.percentile(flat, 99.9))

            if vmax <= vmin:
                vmin = np.min(flat)
                vmax = np.max(flat)

            self._overlay_stretch_vmin = float(vmin)
            self._overlay_stretch_vmax = float(vmax)

        ax.hist(display_data, bins=128, color='#3a6ea5', edgecolor='none', alpha=0.7)
        ax.set_xlim(p_low, p_high)

        vmin = self._overlay_stretch_vmin
        vmax = self._overlay_stretch_vmax

        vmin_display = max(p_low, min(p_high, vmin))
        vmax_display = max(p_low, min(p_high, vmax))

        self._overlay_stretch_marker_min_line = ax.axvline(
            vmin_display, color='#FF5722', linewidth=2, linestyle='-', label=f"Min: {vmin:.1f}"
        )
        self._overlay_stretch_marker_max_line = ax.axvline(
            vmax_display, color='#4CAF50', linewidth=2, linestyle='-', label=f"Max: {vmax:.1f}"
        )

        y_max = ax.get_ylim()[1]
        ax.text(vmin_display, y_max * 0.95, '<', color='#FF5722', fontsize=14,
                ha='center', va='top', fontweight='bold')
        ax.text(vmax_display, y_max * 0.95, '>', color='#4CAF50', fontsize=14,
                ha='center', va='top', fontweight='bold')

        ax.set_xlabel('Pixel Value')
        ax.set_ylabel('Count')
        ax.set_title('Image Histogram')
        ax.legend(loc='upper right', fontsize=8)

        if self.overlay_stretch_plot_info_label:
            self.overlay_stretch_plot_info_label.setText(
                f"Manual range | Min: {vmin:.2f} | Max: {vmax:.2f}"
            )

        self.overlay_stretch_plot_canvas.draw_idle()

    def _overlay_on_stretch_plot_press(self, event):
        """Handle mouse press on stretch plot"""
        if event.inaxes != self.overlay_stretch_plot_ax or event.xdata is None:
            return
        if self._overlay_stretch_vmin is None or self._overlay_stretch_vmax is None:
            return

        x = event.xdata
        dist_to_min = abs(x - self._overlay_stretch_vmin)
        dist_to_max = abs(x - self._overlay_stretch_vmax)
        self._overlay_stretch_drag_target = "min" if dist_to_min < dist_to_max else "max"
        self._overlay_stretch_dragging = True

    def _overlay_on_stretch_plot_motion(self, event):
        """Handle mouse motion on stretch plot (dragging)"""
        if not self._overlay_stretch_dragging or event.xdata is None:
            return

        x = event.xdata
        if self._overlay_stretch_drag_target == "min":
            new_val = min(x, self._overlay_stretch_vmax - 1)
            self._overlay_stretch_vmin = new_val
        else:
            new_val = max(x, self._overlay_stretch_vmin + 1)
            self._overlay_stretch_vmax = new_val

        self._overlay_update_stretch_plot()
        self._overlay_apply_custom_stretch()

    def _overlay_on_stretch_plot_release(self, event):
        """Handle mouse release on stretch plot"""
        self._overlay_stretch_dragging = False
        self._overlay_stretch_drag_target = None

    def _overlay_apply_custom_stretch(self):
        """Apply custom vmin/vmax stretch to the overlay image"""
        if self.overlay_image_data is None:
            return
        if self._overlay_stretch_vmin is None or self._overlay_stretch_vmax is None:
            return

        vmin = self._overlay_stretch_vmin
        vmax = self._overlay_stretch_vmax
        if vmax <= vmin:
            vmax = vmin + 1

        if self.overlay_viewer is not None:
            self.overlay_viewer.set_stretch_mode("linear")
            self.overlay_viewer.set_linear_range(vmin, vmax)

    def _overlay_reset_stretch_plot_values(self):
        """Reset stretch plot values when changing image or stretch mode"""
        self._overlay_stretch_vmin = None
        self._overlay_stretch_vmax = None
        if self.overlay_stretch_plot_dialog and self.overlay_stretch_plot_dialog.isVisible():
            self._overlay_update_stretch_plot()

    def _overlay_normalize_image(self):
        if self.overlay_image_data is None:
            return None

        cache_key = id(self.overlay_image_data)
        if self._overlay_normalized_cache is not None:
            if self._overlay_normalized_cache[0] == cache_key:
                return self._overlay_normalized_cache[1].copy()

        data = self.overlay_image_data
        finite = np.isfinite(data)
        if not finite.any():
            return np.zeros_like(data)

        _, median_val, std_val = sigma_clipped_stats(data[finite], sigma=3.0, maxiters=5)
        vmin = max(np.min(data[finite]), median_val - 2.8 * std_val)
        vmax = min(np.max(data[finite]), np.percentile(data[finite], 99.9))

        if vmax <= vmin:
            vmin = float(np.min(data[finite]))
            vmax = float(np.max(data[finite]))

        normalized = (data - vmin) / (vmax - vmin + 1e-10)
        normalized = np.clip(normalized, 0, 1)

        self._overlay_normalized_cache = (cache_key, normalized)
        return normalized.copy()

    def _overlay_calculate_zscale(self):
        finite = np.isfinite(self.overlay_image_data)
        if not finite.any():
            return 0.0, 1.0
        data = self.overlay_image_data[finite]
        _, median_val, std_val = sigma_clipped_stats(data, sigma=3.0, maxiters=5)
        vmin = float(median_val - 2.8 * std_val)
        vmax_percentile = np.percentile(data, 99.5)
        vmax_sigma = median_val + 6.0 * std_val
        vmax = float(min(vmax_percentile, vmax_sigma))
        if vmax <= vmin:
            vmin = float(np.min(data))
            vmax = float(np.max(data))
        return vmin, vmax

    def _overlay_apply_stretch(self, data):
        return data

    def _overlay_stretch_auto_siril(self, data, intensity):
        finite = data[np.isfinite(data)]
        if finite.size == 0:
            return data
        median_val = np.median(finite)
        mad = np.median(np.abs(finite - median_val))
        sigma = mad * MAD_TO_SIGMA
        shadows = max(0.0, median_val - 2.8 * sigma)
        highlights = 1.0
        stretched = (data - shadows) / (highlights - shadows + 1e-10)
        stretched = np.clip(stretched, 0, 1)
        midtone = 0.15 + (1.0 - intensity) * 0.35
        return self._overlay_mtf_function(stretched, midtone)

    def _overlay_stretch_asinh(self, data, intensity):
        beta = 1.0 + intensity * 15.0
        stretched = np.arcsinh(data * beta) / np.arcsinh(beta)
        return np.clip(stretched, 0, 1)

    def _overlay_stretch_mtf(self, data, intensity):
        midtone = 0.05 + (1.0 - intensity) * 0.45
        return self._overlay_mtf_function(data, midtone)

    def _overlay_mtf_function(self, data, midtone):
        m = np.clip(midtone, 0.001, 0.999)
        result = np.zeros_like(data)
        mask = data > 0
        result[mask] = (m - 1) * data[mask] / ((2 * m - 1) * data[mask] - m)
        result[data == 0] = 0
        result[data == 1] = 1
        return np.clip(result, 0, 1)

    def _overlay_stretch_histogram_eq(self, data):
        finite = data[np.isfinite(data)]
        if finite.size == 0:
            return data
        hist, bin_edges = np.histogram(finite.flatten(), bins=65536, range=(0, 1))
        cdf = hist.cumsum()
        cdf = cdf / cdf[-1]
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        stretched = np.interp(data, bin_centers, cdf)
        return np.clip(stretched, 0, 1)

    def _overlay_stretch_log(self, data, intensity):
        a = 100 + intensity * 900
        stretched = np.log(1 + a * data) / np.log(1 + a)
        return np.clip(stretched, 0, 1)

    def _overlay_stretch_sqrt(self, data, intensity):
        power = 0.2 + (1.0 - intensity) * 0.8
        stretched = np.power(data, power)
        return np.clip(stretched, 0, 1)

    def _overlay_on_scroll(self, event):
        return

    def _overlay_on_button_press(self, event):
        return

    def _overlay_on_button_release(self, event):
        return

    def _overlay_on_motion(self, event):
        return

    # ========================================================================
    # Run Photometry
    # ========================================================================
    def run_photometry(self):
        self._apply_params()

        data_dir = Path(self.data_edit.text())
        output_dir = Path(self.out_edit.text())

        if not data_dir.exists():
            self._set_environment_status("data", "fail", "Directory missing")
            self._set_environment_status("images", "fail", "Not checked")
            QMessageBox.warning(self, "Error", f"Data directory not found:\n{data_dir}")
            return
        if not data_dir.is_dir():
            self._set_environment_status("data", "fail", "Not a directory")
            self._set_environment_status("images", "fail", "Not checked")
            QMessageBox.warning(self, "Error", f"Data path is not a directory:\n{data_dir}")
            return

        pattern = self.pattern_edit.text().strip() or "*.fit*"
        try:
            image_files = sorted(data_dir.glob(pattern))
            image_count = len(image_files)
        except Exception as exc:
            self._set_environment_status("images", "fail", f"Pattern error: {exc}")
            QMessageBox.warning(self, "Error", f"Invalid file pattern:\n{pattern}\n\n{exc}")
            return
        self._set_environment_status("data", "ok", "Directory exists")
        if image_count == 0:
            self._set_environment_status("images", "fail", f"No match: {pattern}")
            QMessageBox.warning(
                self,
                "No Images",
                f"No images found in:\n{data_dir}\n\nPattern: {pattern}",
            )
            return
        self._set_environment_status("images", "ok", f"{image_count} file(s)")

        try:
            output_dir.mkdir(parents=True, exist_ok=True)
        except Exception as exc:
            self._set_environment_status("output", "fail", "Not writable")
            QMessageBox.warning(self, "Error", f"Output directory is not writable:\n{output_dir}\n\n{exc}")
            return
        self._set_environment_status("output", "ok", "Writable")

        use_wsl_scratch = (
            sys.platform == "win32"
            and hasattr(self, "use_wsl_scratch_check")
            and self.use_wsl_scratch_check.isChecked()
        )
        skip_existing = self.skip_existing_check.isChecked()
        pending_files = image_files
        if skip_existing:
            pending_files = [
                path for path in image_files
                if not (output_dir / f"{path.stem}.txt").exists()
            ]
        pending_count = len(pending_files)
        skipped_count = image_count - pending_count
        input_bytes = self._sum_file_sizes(pending_files)

        if pending_count == 0:
            QMessageBox.information(
                self,
                "Nothing To Run",
                f"All {image_count} matching image(s) already have IRAF .txt outputs.\n\n"
                "Disable 'Skip already processed files' to rerun.",
            )
            return

        if use_wsl_scratch and not self._confirm_wsl_scratch_run(pending_count, input_bytes):
            return

        if hasattr(self, "run_bar"):
            self.run_bar.set_running(True)
        else:
            self.run_btn.setEnabled(False)
            self.stop_btn.setEnabled(True)
        self.progress_bar.setVisible(True)
        self.progress_bar.setValue(0)
        self._progress_started_at = time.monotonic()
        self.progress_label.setText(f"0/{pending_count}: Starting | elapsed 00:00 | ETA estimating")
        self.log_text.clear()
        if self.run_log_text is not None:
            self.run_log_text.clear()
        if use_wsl_scratch:
            check = getattr(self, "_scratch_last_check", {})
            free = check.get("free_bytes")
            required = check.get("required_bytes", self._scratch_required_bytes(input_bytes))
            self._on_log(
                "[SCRATCH] Enabled: "
                f"copy {pending_count} pending FITS ({self._format_bytes(input_bytes)}) "
                "to WSL /tmp/apex_iraf per frame before IRAF; "
                f"required~{self._format_bytes(required)}, "
                f"free={self._format_bytes(free) if free is not None else 'unknown'}; "
                f"final output={output_dir}; skipped existing={skipped_count}"
            )

        self._auto_save_params()
        filter_params, filter_aliases = self._load_iraf_toml_config()
        param_defaults = self._build_iraf_param_defaults()
        use_apex_frames = bool(self.use_apex_step4_check.isChecked())
        apex_frame_params = self._load_apex_frame_params() if use_apex_frames else {}
        if apex_frame_params:
            self._on_log(f"[APEX] Using Step 4 frame overrides for {len(apex_frame_params)} frame keys")

        self.worker = IRAFPhotometryWorker(
            data_dir=data_dir,
            output_dir=output_dir,
            file_pattern=pattern,
            params=self.iraf_params,
            auto_sigma=self.auto_sigma_check.isChecked(),
            skip_existing=self.skip_existing_check.isChecked(),
            filter_params=filter_params,
            filter_aliases=filter_aliases,
            param_defaults=param_defaults,
            apex_frame_params=apex_frame_params,
            use_apex_frame_params=use_apex_frames,
            use_wsl_scratch=use_wsl_scratch,
            scratch_root="/tmp/apex_iraf",
        )

        self.worker.progress.connect(self._on_progress)
        self.worker.log.connect(self._on_log)
        self.worker.finished.connect(self._on_finished)
        self.worker.error.connect(self._on_error)

        self.worker.start()

    def stop_photometry(self):
        if self.worker:
            if hasattr(self, "run_bar"):
                self.run_bar.set_stopping()
            self.worker.stop()

    @staticmethod
    def _format_progress_duration(seconds: float) -> str:
        seconds = max(0, int(seconds))
        hours, remainder = divmod(seconds, 3600)
        minutes, secs = divmod(remainder, 60)
        if hours:
            return f"{hours:d}:{minutes:02d}:{secs:02d}"
        return f"{minutes:02d}:{secs:02d}"

    @staticmethod
    def _timestamp_log_message(msg: str) -> str:
        stamp = time.strftime("%Y-%m-%d %H:%M:%S")
        lines = str(msg).splitlines() or [""]
        return "\n".join(f"[{stamp}] {line}" if line else "" for line in lines)

    def _on_progress(self, current, total, message):
        total = max(0, int(total or 0))
        current = max(0, int(current or 0))
        if total:
            current = min(current, total)

        self.progress_bar.setMaximum(total)
        self.progress_bar.setValue(current)

        started_at = getattr(self, "_progress_started_at", None)
        elapsed = max(0.0, time.monotonic() - started_at) if started_at else 0.0
        elapsed_text = self._format_progress_duration(elapsed)
        eta_text = "estimating"
        if started_at and current > 0 and total > 0:
            seconds_per_frame = elapsed / current
            eta_text = self._format_progress_duration(seconds_per_frame * (total - current))

        self.progress_label.setText(
            f"{current}/{total}: {message} | elapsed {elapsed_text} | ETA {eta_text}"
        )

    def show_run_log_window(self):
        if self.run_log_window is None:
            self.run_log_window = QWidget(None, Qt.Window)
            self.run_log_window.setWindowTitle("IRAF Photometry Log")
            self.run_log_window.resize(900, 500)
            layout = QVBoxLayout(self.run_log_window)
            self.run_log_text = QTextEdit(self.run_log_window)
            self.run_log_text.setReadOnly(True)
            self.run_log_text.setObjectName("Log")
            layout.addWidget(self.run_log_text)
        if self.run_log_text is not None and hasattr(self, "log_text"):
            self.run_log_text.setPlainText(self.log_text.toPlainText())
            scrollbar = self.run_log_text.verticalScrollBar()
            scrollbar.setValue(scrollbar.maximum())
        self.run_log_window.show()
        self.run_log_window.raise_()
        self.run_log_window.activateWindow()

    def _on_log(self, msg):
        text = self._timestamp_log_message(msg)
        self.log_text.append(text)
        scrollbar = self.log_text.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())
        if self.run_log_text is not None:
            self.run_log_text.append(text)
            log_scrollbar = self.run_log_text.verticalScrollBar()
            log_scrollbar.setValue(log_scrollbar.maximum())

    def show_comparison_log(self):
        if hasattr(self, "tabs"):
            self.tabs.setCurrentIndex(2)
        if hasattr(self, "cmp_tabs"):
            self.cmp_tabs.setCurrentIndex(1)

    def _sync_comparison_after_iraf_run(self, output_dir: str):
        if not output_dir:
            return
        out_path = Path(output_dir)
        if hasattr(self, "cmp_iraf_edit"):
            self.cmp_iraf_edit.setText(str(out_path))
        if hasattr(self, "cmp_image_edit"):
            self.cmp_image_edit.setText(self.data_edit.text())
        self._refresh_path_summaries()
        txt_count = 0
        coo_count = 0
        try:
            txt_count = sum(1 for _ in out_path.glob("*.txt"))
            coo_count = sum(1 for _ in out_path.glob("*.coo"))
        except Exception:
            pass
        self._on_log(
            f"[COMPARISON] IRAF result path synced: {out_path} "
            f"({txt_count} txt, {coo_count} coo). Scratch files are not used for comparison."
        )
        if hasattr(self, "overlay_file_combo"):
            self._overlay_refresh_maps()

    def _on_finished(self, result):
        if hasattr(self, "run_bar"):
            self.run_bar.set_running(False)
        else:
            self.run_btn.setEnabled(True)
            self.stop_btn.setEnabled(False)
        self.progress_bar.setVisible(False)
        n = len(result.get("results", []))
        started_at = getattr(self, "_progress_started_at", None)
        elapsed = time.monotonic() - started_at if started_at else 0.0
        self._sync_comparison_after_iraf_run(result.get("output_dir", ""))
        if result.get("stopped"):
            self.progress_label.setText(
                f"Stopped: {n} images completed | elapsed {self._format_progress_duration(elapsed)}"
            )
            self._on_log(
                f"[STOP] User stopped IRAF photometry. Completed frames: {n}. "
                f"Output kept in: {result.get('output_dir', '')}"
            )
            return
        self.progress_label.setText(
            f"Completed: {n} images processed | elapsed {self._format_progress_duration(elapsed)}"
        )
        QMessageBox.information(self, "Complete",
            f"IRAF photometry completed.\n{n} images processed.\nOutput: {result.get('output_dir', '')}")

    def _on_error(self, msg):
        if hasattr(self, "run_bar"):
            self.run_bar.set_running(False)
        else:
            self.run_btn.setEnabled(True)
            self.stop_btn.setEnabled(False)
        self.progress_bar.setVisible(False)
        self._on_log(f"[ERROR] {msg}")
        QMessageBox.critical(self, "Error", msg)

    def _cmp_log(self, msg: str):
        if not hasattr(self, "cmp_log_text") or self.cmp_log_text is None:
            return
        self.cmp_log_text.append(self._timestamp_log_message(msg))
        scrollbar = self.cmp_log_text.verticalScrollBar()
        scrollbar.setValue(scrollbar.maximum())

    # ========================================================================
    # Comparison
    # ========================================================================
    def run_comparison(self):
        apex_dir = Path(self.cmp_apex_edit.text())
        iraf_dir = Path(self.cmp_iraf_edit.text())
        tol = self.cmp_tol.value()
        cmp_zmag = float(self.pp_zmag.value()) if hasattr(self, "pp_zmag") else float(self.iraf_params.zmag)

        if hasattr(self, "cmp_log_text") and self.cmp_log_text is not None:
            self.cmp_log_text.clear()
        self._cmp_log("Starting comparison")
        self._cmp_log(f"APEX dir: {apex_dir}")
        self._cmp_log(f"IRAF dir: {iraf_dir}")
        self._cmp_log(f"Tolerance: {tol:.2f} px")
        self._cmp_log(
            f"IRAF zmag: {cmp_zmag:.2f}; APEX mag_inst is shifted by +zmag only when raw instrumental scale is detected."
        )

        if not apex_dir.exists():
            self._cmp_log(f"ERROR: APEX dir not found: {apex_dir}")
            QMessageBox.warning(self, "Error", f"APEX dir not found: {apex_dir}")
            return
        if not iraf_dir.exists():
            self._cmp_log(f"ERROR: IRAF dir not found: {iraf_dir}")
            QMessageBox.warning(self, "Error", f"IRAF dir not found: {iraf_dir}")
            return

        self._overlay_refresh_maps()

        # Collect files
        apex_map = {}
        for p in _iter_apex_photometry_files(apex_dir):
            key = _normalize_frame_key(p.stem)
            apex_map.setdefault(key, p)

        iraf_map = {}
        for p in iraf_dir.rglob("*.txt"):
            key = _normalize_frame_key(p.stem)
            iraf_map.setdefault(key, p)

        self._cmp_log(f"APEX files: {len(apex_map)}")
        self._cmp_log(f"IRAF files: {len(iraf_map)}")

        frames = sorted(set(apex_map) & set(iraf_map))
        if not frames:
            self._cmp_log("ERROR: No matching frames found.")
            QMessageBox.warning(self, "Error", "No matching frames found.")
            return
        self._cmp_log(f"Matched frames: {len(frames)}")

        self.frame_rows = []
        self.frame_matches = {}
        all_matches = []

        def _match_with_iraf(iraf_df, apex_df, x_col, y_col, mag_col):
            axy = apex_df[[x_col, y_col]].to_numpy(float)
            ixy = iraf_df[["x", "y"]].to_numpy(float)

            if axy.size == 0 or ixy.size == 0:
                return pd.DataFrame()

            tree = cKDTree(axy)
            dist, idx = tree.query(ixy, distance_upper_bound=tol)
            mask = np.isfinite(dist) & (dist <= tol)

            if not np.any(mask):
                return pd.DataFrame()

            iraf_mag = iraf_df.loc[mask, "mag"].to_numpy(float)
            apex_mag_raw = apex_df.loc[idx[mask], mag_col].to_numpy(float)
            raw_dmag = apex_mag_raw - iraf_mag
            apex_zmag_offset = 0.0
            if mag_col == "mag_inst" and np.isfinite(cmp_zmag) and cmp_zmag != 0.0:
                finite_raw = raw_dmag[np.isfinite(raw_dmag)]
                if finite_raw.size:
                    raw_med = float(np.nanmedian(finite_raw))
                    shifted_med = raw_med + cmp_zmag
                    if np.isfinite(shifted_med) and abs(shifted_med) < abs(raw_med):
                        apex_zmag_offset = cmp_zmag
            apex_mag = apex_mag_raw + apex_zmag_offset

            match = pd.DataFrame({
                "iraf_x": iraf_df.loc[mask, "x"].to_numpy(),
                "iraf_y": iraf_df.loc[mask, "y"].to_numpy(),
                "iraf_mag": iraf_mag,
                "apex_x": apex_df.loc[idx[mask], x_col].to_numpy(),
                "apex_y": apex_df.loc[idx[mask], y_col].to_numpy(),
                "apex_mag": apex_mag,
                "apex_mag_raw": apex_mag_raw,
                "apex_zmag_offset": np.full(len(apex_mag), apex_zmag_offset, dtype=float),
                "dist_px": dist[mask],
            })
            match["dx"] = match["apex_x"] - match["iraf_x"]
            match["dy"] = match["apex_y"] - match["iraf_y"]
            match["raw_dmag"] = raw_dmag
            match["dmag"] = match["apex_mag"] - match["iraf_mag"]
            return match

        for frame in frames:
            apex_df = _read_apex_tsv(apex_map[frame])
            iraf_df = _read_iraf_txt(iraf_map[frame])

            mag_col = _pick_first(apex_df.columns, ["mag_inst", "mag", "mag_raw"])
            x_col = _pick_first(apex_df.columns, ["xcenter", "x", "x_init"])
            y_col = _pick_first(apex_df.columns, ["ycenter", "y", "y_init"])

            n_apex_total = len(apex_df)
            n_iraf_total = len(iraf_df)

            if not all([mag_col, x_col, y_col]):
                self._cmp_log(
                    f"{frame}: missing columns (mag={mag_col}, x={x_col}, y={y_col})"
                )
                continue

            match = _match_with_iraf(iraf_df, apex_df, x_col, y_col, mag_col)
            if not match.empty:
                dx_med = float(np.nanmedian(match["dx"]))
                dy_med = float(np.nanmedian(match["dy"]))
            else:
                dx_med = np.nan
                dy_med = np.nan

            shift_x = _auto_axis_shift(dx_med)
            shift_y = _auto_axis_shift(dy_med)
            if (shift_x != 0.0) or (shift_y != 0.0):
                iraf_adj = iraf_df.copy()
                iraf_adj["x"] = iraf_adj["x"] + shift_x
                iraf_adj["y"] = iraf_adj["y"] + shift_y
                match = _match_with_iraf(iraf_adj, apex_df, x_col, y_col, mag_col)

            if match.empty:
                self._cmp_log(f"{frame}: matched 0 (no pairs within {tol:.2f}px)")
                self.frame_rows.append({
                    "frame": frame,
                    "n": 0,
                    "dmag_med": np.nan,
                    "dmag_std": np.nan,
                    "resid_std": np.nan,
                    "dx_med": np.nan,
                    "dy_med": np.nan,
                    "dist_med": np.nan,
                    "dist_p95": np.nan,
                    "frac_within_tol": 0.0,
                    "best_shift_x": np.nan,
                    "best_shift_y": np.nan,
                    "n_iraf_total": n_iraf_total,
                    "n_apex_total": n_apex_total,
                    "apex_zmag_offset": np.nan,
                })
                self.frame_matches[frame] = pd.DataFrame()
                continue

            dmag_med = float(np.nanmedian(match["dmag"]))
            dmag_std = float(np.nanstd(match["dmag"]))
            dx_med = float(np.nanmedian(match["dx"]))
            dy_med = float(np.nanmedian(match["dy"]))
            apex_zmag_offset = float(np.nanmedian(match["apex_zmag_offset"]))

            dist_vals = match["dist_px"].to_numpy(float)
            dist_med = float(np.nanmedian(dist_vals)) if dist_vals.size else np.nan
            dist_p95 = float(np.nanpercentile(dist_vals, 95)) if dist_vals.size else np.nan
            frac_within = float(np.mean(dist_vals <= tol)) if dist_vals.size else 0.0

            best_shift_x = BASE_IRAF_SHIFT + shift_x
            best_shift_y = BASE_IRAF_SHIFT + shift_y

            self._cmp_log(
                f"{frame}: matched {len(match)} | dmag={dmag_med:.4f}±{dmag_std:.4f} "
                f"dx={dx_med:.3f} dy={dy_med:.3f} "
                f"dist_med={dist_med:.3f} p95={dist_p95:.3f} "
                f"shift=({best_shift_x:.1f},{best_shift_y:.1f}) "
                f"n_iraf={n_iraf_total} n_apex={n_apex_total}"
            )

            self.frame_rows.append({
                "frame": frame,
                "n": len(match),
                "dmag_med": dmag_med,
                "dmag_std": dmag_std,
                "resid_std": dmag_std,
                "dx_med": dx_med,
                "dy_med": dy_med,
                "dist_med": dist_med,
                "dist_p95": dist_p95,
                "frac_within_tol": frac_within,
                "best_shift_x": best_shift_x,
                "best_shift_y": best_shift_y,
                "n_iraf_total": n_iraf_total,
                "n_apex_total": n_apex_total,
                "apex_zmag_offset": apex_zmag_offset,
            })
            self.frame_matches[frame] = match
            all_matches.append(match.assign(frame=frame))

        self.matched_all = pd.concat(all_matches, ignore_index=True) if all_matches else pd.DataFrame()

        # Update table
        self.cmp_table.setRowCount(len(self.frame_rows))
        for i, row in enumerate(self.frame_rows):
            items = [
                row["frame"],
                str(row["n"]),
                f"{row['dmag_med']:.4f}" if np.isfinite(row["dmag_med"]) else "nan",
                f"{row.get('resid_std', row['dmag_std']):.4f}" if np.isfinite(row.get("resid_std", row["dmag_std"])) else "nan",
                f"{row['dx_med']:.3f}" if np.isfinite(row["dx_med"]) else "nan",
                f"{row['dy_med']:.3f}" if np.isfinite(row["dy_med"]) else "nan",
                f"{row['dist_med']:.3f}" if np.isfinite(row["dist_med"]) else "nan",
                f"{row['dist_p95']:.3f}" if np.isfinite(row["dist_p95"]) else "nan",
                f"{row['frac_within_tol']:.3f}",
                f"{row['best_shift_x']:.1f}" if np.isfinite(row["best_shift_x"]) else "nan",
                f"{row['best_shift_y']:.1f}" if np.isfinite(row["best_shift_y"]) else "nan",
                str(int(row["n_iraf_total"])),
                str(int(row["n_apex_total"])),
            ]
            for col, text in enumerate(items):
                self.cmp_table.setItem(i, col, QTableWidgetItem(text))
        self.cmp_table.resizeColumnsToContents()

        # Update summary
        if not self.matched_all.empty:
            total = len(self.matched_all)
            med = np.nanmedian(self.matched_all["dmag"])
            std = np.nanstd(self.matched_all["dmag"])
            if "apex_zmag_offset" in self.matched_all:
                zp_offsets = self.matched_all["apex_zmag_offset"].to_numpy(float)
                if zp_offsets.size and np.any(np.isfinite(zp_offsets)):
                    zp_offset = float(np.nanmedian(zp_offsets))
                    if abs(zp_offset) > 1e-9:
                        self._cmp_log(
                            f"Applied APEX mag_inst +{zp_offset:.2f} before comparison to match IRAF zmag scale."
                        )
            self.cmp_summary.setText(
                f"Total: {total} matched stars | dmag median: {med:.4f} | residual std: {std:.4f}"
            )
            self._cmp_log(
                f"Total matched: {total} | dmag median: {med:.4f} | residual std: {std:.4f}"
            )

            frame_filters = {}
            for row in self.frame_rows:
                frame_key = row.get("frame")
                if frame_key is None:
                    continue
                filt = self._filter_for_frame_key(frame_key)
                if filt:
                    frame_filters[frame_key] = filt

            if frame_filters:
                df = self.matched_all.copy()
                df["filter"] = df["frame"].map(frame_filters)
                self._cmp_log("Filter summary:")
                for filt, grp in df.groupby("filter"):
                    dmag_med = float(np.nanmedian(grp["dmag"]))
                    dmag_std = float(np.nanstd(grp["dmag"]))
                    dist_med = float(np.nanmedian(grp["dist_px"])) if "dist_px" in grp else np.nan
                    dist_p95 = float(np.nanpercentile(grp["dist_px"], 95)) if "dist_px" in grp else np.nan
                    frac_within = float(np.mean(grp["dist_px"] <= tol)) if "dist_px" in grp else np.nan
                    self._cmp_log(
                        f"  {filt}: n={len(grp)} | dmag={dmag_med:.4f}±{dmag_std:.4f} "
                        f"dist_med={dist_med:.3f} p95={dist_p95:.3f} frac<=tol={frac_within:.3f}"
                    )
        else:
            self.cmp_summary.setText("No matches found.")
            self._cmp_log("No matches found.")

        if self.cmp_table.rowCount() > 0:
            self.cmp_table.selectRow(0)

    def _plot_comparison(self):
        items = self.cmp_table.selectedItems()
        if not items:
            return

        frame = self.cmp_table.item(items[0].row(), 0).text()
        match = self.frame_matches.get(frame, pd.DataFrame())

        self._overlay_set_frame_key(frame)

        self.cmp_fig.clear()

        if match.empty:
            ax = self.cmp_fig.add_subplot(111)
            ax.text(0.5, 0.5, f"{frame}: No matches", ha="center", va="center")
            self.cmp_canvas.draw_idle()
            return

        ax1 = self.cmp_fig.add_subplot(221)
        ax2 = self.cmp_fig.add_subplot(222)
        ax3 = self.cmp_fig.add_subplot(223)
        ax4 = self.cmp_fig.add_subplot(224)

        dmag = match["dmag"].astype(float)
        med = float(np.nanmedian(dmag))
        resid = dmag - med
        scatter = float(np.nanstd(resid))
        zp_offset = 0.0
        if "apex_zmag_offset" in match:
            zp_vals = match["apex_zmag_offset"].to_numpy(float)
            if zp_vals.size and np.any(np.isfinite(zp_vals)):
                zp_offset = float(np.nanmedian(zp_vals))
        apex_label = "APEX mag_inst + zmag" if abs(zp_offset) > 1e-9 else "APEX mag"

        # Median-removed dmag vs mag
        ax1.scatter(match["apex_mag"], resid, s=10, alpha=0.6)
        ax1.axhline(0, color="red", ls="--")
        ax1.axhline(scatter, color="gray", ls=":", lw=0.8)
        ax1.axhline(-scatter, color="gray", ls=":", lw=0.8)
        ax1.set_xlabel(apex_label)
        ax1.set_ylabel("dmag residual")
        ax1.set_title(f"Residual vs mag (med={med:.4f})")

        # Residual dmag hist
        ax2.hist(resid.dropna(), bins=30, alpha=0.7, edgecolor="black")
        ax2.axvline(0, color="red", ls="--")
        ax2.axvline(scatter, color="gray", ls=":", lw=0.8)
        ax2.axvline(-scatter, color="gray", ls=":", lw=0.8)
        ax2.set_xlabel("dmag residual")
        ax2.set_title(f"resid std={scatter:.4f}")

        # dx vs dy
        ax3.scatter(match["dx"], match["dy"], s=10, alpha=0.6)
        ax3.axhline(0, color="gray", lw=0.5)
        ax3.axvline(0, color="gray", lw=0.5)
        ax3.set_xlabel("dx [px]")
        ax3.set_ylabel("dy [px]")
        ax3.set_title("Position offset")
        ax3.set_aspect("equal")

        # 1:1 after removing the median magnitude offset.
        apex_aligned = match["apex_mag"] - med
        ax4.scatter(match["iraf_mag"], apex_aligned, s=10, alpha=0.6)
        lims = [min(match["iraf_mag"].min(), apex_aligned.min()),
                max(match["iraf_mag"].max(), apex_aligned.max())]
        ax4.plot(lims, lims, "r--")
        ax4.set_xlabel("IRAF mag")
        ax4.set_ylabel(f"{apex_label} - median dmag")
        ax4.set_title("1:1 (median offset removed)")

        self.cmp_fig.suptitle(f"{frame} (N={len(match)})")
        self.cmp_fig.tight_layout()
        self.cmp_canvas.draw_idle()

    def _overlay_set_frame_key(self, frame_key: str):
        if not hasattr(self, "overlay_index_spin") or not frame_key:
            return
        img_dir = self._overlay_image_dir()
        if self.overlay_last_image_dir is None or self.overlay_last_image_dir != img_dir:
            self._overlay_reload()
        if not self.overlay_key_to_index:
            return
        idx = self.overlay_key_to_index.get(frame_key)
        if idx is None:
            self._overlay_render_empty(f"No matching image for frame: {frame_key}")
            self._cmp_log(f"Overlay: no image for frame {frame_key}")
            return
        if idx != self.overlay_current_index:
            self.overlay_index_spin.setValue(idx)
        else:
            self._overlay_load_current()

    def export_comparison(self):
        if self.matched_all is None or self.matched_all.empty:
            self._cmp_log("Export skipped: no comparison data available.")
            QMessageBox.information(self, "Export", "No data to export.")
            return

        out_dir = Path(self.cmp_apex_edit.text()) / "iraf_comparison"
        out_dir.mkdir(parents=True, exist_ok=True)

        self.matched_all.to_csv(out_dir / "iraf_compare_all.csv", index=False)
        pd.DataFrame(self.frame_rows).to_csv(out_dir / "iraf_compare_summary.csv", index=False)

        self._cmp_log(f"Exported CSVs to: {out_dir}")
        QMessageBox.information(self, "Export", f"Saved to:\n{out_dir}")

    # ========================================================================
    # Utilities
    # ========================================================================
    def _browse_dir(self, line_edit):
        path = QFileDialog.getExistingDirectory(self, "Select Directory", line_edit.text())
        if path:
            line_edit.setText(path)
            if line_edit in (self.data_edit, self.out_edit):
                self._on_photometry_path_settings_changed()
            elif line_edit in (self.cmp_apex_edit, self.cmp_iraf_edit):
                self._on_comparison_path_settings_changed()
            elif line_edit is self.cmp_image_edit:
                self._on_comparison_image_path_changed()
