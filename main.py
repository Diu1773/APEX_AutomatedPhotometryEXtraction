#!/usr/bin/env python3
"""
APEX — Automated Photometry EXtraction
Launcher: choose between Cluster CMD photometry or Light Curve analysis.
"""

import sys
import os
import warnings
from pathlib import Path

from apex import __version__
from apex.utils.ssl_certificates import configure_ssl_certificates

warnings.filterwarnings("ignore", message=".*datfix.*MJD-OBS.*", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning, module="astropy")
try:
    from astropy.wcs import FITSFixedWarning
    warnings.filterwarnings("ignore", category=FITSFixedWarning)
except ImportError:
    pass
warnings.filterwarnings("ignore", message=".*tight_layout.*", category=UserWarning)

from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QFrame, QMessageBox
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont, QFontDatabase, QIcon, QPixmap, QPainter
from apex.gui.tools.registry import iter_tool_modules

_FROZEN = getattr(sys, "frozen", False)
_SSL_CERT_OK, _SSL_CERT_DETAIL = configure_ssl_certificates()

if _FROZEN:
    _APP_DIR = Path(sys.executable).parent
    _RESOURCES = Path(sys._MEIPASS) / "apex" / "resources"
else:
    _APP_DIR = Path(__file__).parent
    _RESOURCES = _APP_DIR / "apex" / "resources"


_SMOKE_IMPORTS = (
    "apex.gui.workflow.ui_helpers",
    "apex.gui.main_window",
    "apex.gui.workflow.target_resolver",
    "apex.gui.workflow.step6_ref_build",
    "apex.gui.workflow.step7_forced_aperture_phot",
    "apex.gui.workflow.cmd.step8_psf_photometry",
    "astroquery.gaia",
    "astroquery.simbad",
    "astroquery.utils.tap.core",
    "certifi",
) + tuple(iter_tool_modules())


def _svg_to_pixmap(svg_path: Path, size: int) -> QPixmap:
    # Rasterize at the display's device-pixel ratio so the vector logo stays
    # crisp under HiDPI scaling instead of being upscaled from 1×.
    app = QApplication.instance()
    screen = app.primaryScreen() if app is not None else None
    dpr = float(screen.devicePixelRatio()) if screen is not None else 1.0
    side = max(1, int(round(size * dpr)))
    px = QPixmap(side, side)
    px.fill(Qt.transparent)
    try:
        from PyQt5.QtSvg import QSvgRenderer
        r = QSvgRenderer(str(svg_path))
        p = QPainter(px)
        r.render(p)
        p.end()
    except Exception:
        pass
    px.setDevicePixelRatio(dpr)
    return px


def _make_app_icon() -> QIcon:
    svg = _RESOURCES / "logo_base.svg"
    if not svg.exists():
        return QIcon()
    px = _svg_to_pixmap(svg, 256)
    icon = QIcon()
    for s in (16, 24, 32, 48, 64, 128, 256):
        icon.addPixmap(px.scaled(s, s, Qt.KeepAspectRatio, Qt.SmoothTransformation))
    return icon


def _ensure_parameters() -> str | None:
    """Create the workspace JSON from the bundled example on first run."""
    params = _APP_DIR / "apex_config.json"
    if params.exists():
        return None
    if _FROZEN:
        example = Path(sys._MEIPASS) / "parameters.example.toml"
    else:
        example = _APP_DIR / "parameters.example.toml"
    if not example.exists():
        return f"Default parameter file is missing:\n{example}"
    if example.exists():
        try:
            from apex.config.config_io import load_config_data, save_config_data
            data, _ = load_config_data(example)
            if not save_config_data(params, data):
                raise OSError("write failed")
        except Exception as exc:
            return (
                "Could not create runtime parameters.toml.\n\n"
                f"Target: {params}\n"
                f"Source: {example}\n\n"
                f"{exc}"
            )
    return None


def _run_smoke() -> int:
    """Import the build-critical modules and exit without opening the GUI."""
    import importlib
    import traceback

    failures = []
    for module_name in _SMOKE_IMPORTS:
        try:
            importlib.import_module(module_name)
        except Exception:
            failures.append((module_name, traceback.format_exc()))
    if failures:
        for module_name, tb in failures:
            print(f"[ERROR] Smoke import failed: {module_name}")
            print(tb)
        return 1
    print("[OK] APEX smoke imports passed.")
    return 0


def _configure_app_fonts(app: QApplication) -> None:
    candidate_font_files = [
        Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts" / "malgun.ttf",
        Path("/mnt/c/Windows/Fonts/malgun.ttf"),
    ]
    candidate_families = [
        "Malgun Gothic", "맑은 고딕", "NanumGothic", "Nanum Gothic",
        "Noto Sans CJK KR", "Noto Sans KR", "AppleGothic", "Arial Unicode MS",
    ]
    resolved_family = None
    try:
        import matplotlib
        from matplotlib import font_manager
        for font_path in candidate_font_files:
            if not font_path.exists():
                continue
            try:
                font_manager.fontManager.addfont(str(font_path))
                resolved_family = font_manager.FontProperties(fname=str(font_path)).get_name()
                break
            except Exception:
                continue
        if resolved_family is None:
            available = {f.name.lower(): f.name for f in font_manager.fontManager.ttflist}
            for family in candidate_families:
                if family.lower() in available:
                    resolved_family = available[family.lower()]
                    break
        if resolved_family:
            matplotlib.rcParams["font.family"] = [resolved_family, "DejaVu Sans"]
            matplotlib.rcParams["axes.unicode_minus"] = False
            matplotlib.rcParams["mathtext.fontset"] = "dejavusans"
            matplotlib.rcParams["mathtext.default"] = "regular"
            matplotlib.rcParams["axes.formatter.use_mathtext"] = False
    except Exception:
        pass
    try:
        db = QFontDatabase()
        qt_families = {family.lower(): family for family in db.families()}
        qt_family = None
        if resolved_family and resolved_family.lower() in qt_families:
            qt_family = qt_families[resolved_family.lower()]
        else:
            for family in candidate_families:
                if family.lower() in qt_families:
                    qt_family = qt_families[family.lower()]
                    break
        if qt_family:
            app.setFont(QFont(qt_family, 9))
    except Exception:
        pass


def _launch(launcher: "LauncherWindow", mode: str) -> None:
    launcher.hide()
    try:
        from apex.gui.main_window import MainWindowWorkflow
        window = MainWindowWorkflow(mode=mode)
        window.setAttribute(Qt.WA_DeleteOnClose)
        window.destroyed.connect(launcher.show)
        launcher._mode_window = window
        window.show()
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        QMessageBox.critical(
            launcher,
            f"APEX {mode.upper()} — 실행 오류",
            f"모드 실행 중 오류가 발생했습니다:\n{e}\n\n{tb}",
        )
        launcher.show()


class ModeButton(QPushButton):
    def __init__(self, title: str, subtitle: str, color: str, hover: str, pressed: str, parent=None):
        super().__init__(parent)
        self.setText(f"{title}\n{subtitle}")
        self.setMinimumHeight(110)
        self.setFont(QFont("", 14, QFont.Bold))
        self.setStyleSheet(f"""
            QPushButton {{
                background-color: {color};
                color: white;
                border-radius: 10px;
                padding: 16px;
                text-align: center;
                line-height: 1.4;
            }}
            QPushButton:hover {{ background-color: {hover}; }}
            QPushButton:pressed {{ background-color: {pressed}; }}
        """)


class LauncherWindow(QWidget):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("APEX — Automated Photometry EXtraction")
        self.setMinimumSize(560, 380)
        self.setMaximumSize(700, 460)
        self._mode_window = None
        self._setup_ui()

    def _setup_ui(self):
        root = QVBoxLayout(self)
        root.setSpacing(14)
        root.setContentsMargins(48, 28, 48, 28)

        svg = _RESOURCES / "logo_base.svg"
        if svg.exists():
            logo_px = _svg_to_pixmap(svg, 160)
            logo_lbl = QLabel()
            logo_lbl.setPixmap(logo_px)
            logo_lbl.setAlignment(Qt.AlignCenter)
            logo_lbl.setStyleSheet("background: transparent;")
            root.addWidget(logo_lbl)
        else:
            title = QLabel("APEX")
            title.setAlignment(Qt.AlignCenter)
            title.setFont(QFont("", 38, QFont.Bold))
            root.addWidget(title)

        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setStyleSheet("color: #cccccc;")
        root.addSpacing(4)
        root.addWidget(sep)
        root.addSpacing(4)

        btn_cmd = ModeButton(
            "성단측광", "Cluster CMD Photometry",
            "#1565C0", "#1976D2", "#0D47A1",
        )
        btn_cmd.clicked.connect(lambda: _launch(self, "cmd"))

        btn_lc = ModeButton(
            "시계열분석", "Light Curve Analysis",
            "#2E7D32", "#388E3C", "#1B5E20",
        )
        btn_lc.clicked.connect(lambda: _launch(self, "lc"))

        btn_row = QHBoxLayout()
        btn_row.setSpacing(20)
        btn_row.addWidget(btn_cmd)
        btn_row.addWidget(btn_lc)
        root.addLayout(btn_row)

        footer = QLabel(f"APEX v{__version__}")
        footer.setAlignment(Qt.AlignCenter)
        footer.setStyleSheet("color: #aaaaaa; font-size: 9px;")
        root.addSpacing(4)
        root.addWidget(footer)

    def closeEvent(self, event):
        app = QApplication.instance()
        if app is not None:
            app.quit()
        event.accept()


def main() -> int:
    os.chdir(_APP_DIR)
    if "--smoke" in sys.argv:
        return _run_smoke()

    # Render natively at the display's scale factor (Windows 125%/150% etc.).
    # Without this Qt draws at 1× and Windows bitmap-upscales the whole window,
    # which is what makes text fuzzy and the SVG logo look pixelated.
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
    app = QApplication(sys.argv)
    app.setQuitOnLastWindowClosed(False)
    app.setApplicationName("APEX")
    app.setOrganizationName("APEX Project")
    param_error = _ensure_parameters()
    if param_error:
        QMessageBox.critical(None, "Startup Error", param_error)
        return 1
    _configure_app_fonts(app)
    from apex.gui.theme import apply_theme
    apply_theme(app)
    app.setWindowIcon(_make_app_icon())

    try:
        window = LauncherWindow()
        window.show()
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setWindowTitle("Startup Error")
        msg.setText(str(e))
        msg.setDetailedText(tb)
        msg.exec_()
        return 1

    return app.exec_()


if __name__ == "__main__":
    sys.exit(main())
