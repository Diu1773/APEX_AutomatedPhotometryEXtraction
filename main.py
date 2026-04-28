#!/usr/bin/env python3
"""
APEX — Automated Photometry EXtraction
Launcher: choose between Cluster CMD photometry or Light Curve analysis.
"""

import sys
import os
import subprocess
import warnings
from pathlib import Path

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
from PyQt5.QtGui import QFont, QFontDatabase

_HERE = Path(__file__).parent

# APEX mode entry points
_CMD_MAIN = _HERE / "apex" / "cmd" / "main.py"
_LC_MAIN  = _HERE / "apex" / "lightcurve" / "main.py"

# Legacy fallbacks (used when APEX mode entry point is missing)
_LEGACY_CMD = (
    _HERE.parent
    / "Aperture_Photometry_KNUEMAO"
    / "AAPKC_GUI_project"
    / "main.py"
)
_LEGACY_LC = (
    _HERE.parent
    / "Aperture_Photometry_KNUEMAO"
    / "AAPKL_GUI_project"
    / "main.py"
)

_SUBPROC_KWARGS: dict = {}
if sys.platform == "win32":
    _SUBPROC_KWARGS["creationflags"] = subprocess.CREATE_NO_WINDOW
    _SUBPROC_KWARGS["stdin"] = subprocess.DEVNULL


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


def _launch(apex_script: Path, legacy_script: Path) -> None:
    script = apex_script if apex_script.exists() else legacy_script
    if not script.exists():
        QMessageBox.critical(
            None, "실행 오류",
            f"프로젝트를 찾을 수 없습니다:\n{script}\n\n"
            "경로가 변경된 경우 main.py의 경로 설정을 확인하세요."
        )
        return
    subprocess.Popen(
        [sys.executable, str(script)],
        cwd=str(script.parent),
        **_SUBPROC_KWARGS,
    )


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
        self._setup_ui()

    def _setup_ui(self):
        root = QVBoxLayout(self)
        root.setSpacing(18)
        root.setContentsMargins(48, 36, 48, 32)

        # ── Header ──────────────────────────────────────────────────────────
        title = QLabel("APEX")
        title.setAlignment(Qt.AlignCenter)
        title.setFont(QFont("", 38, QFont.Bold))

        subtitle = QLabel("Automated Photometry EXtraction")
        subtitle.setAlignment(Qt.AlignCenter)
        f = QFont("", 11)
        f.setItalic(True)
        subtitle.setFont(f)
        subtitle.setStyleSheet("color: #555555;")

        root.addWidget(title)
        root.addWidget(subtitle)

        sep = QFrame()
        sep.setFrameShape(QFrame.HLine)
        sep.setStyleSheet("color: #cccccc;")
        root.addSpacing(4)
        root.addWidget(sep)
        root.addSpacing(4)

        # ── Mode buttons ────────────────────────────────────────────────────
        btn_cmd = ModeButton(
            "성단측광", "Cluster CMD Photometry",
            "#1565C0", "#1976D2", "#0D47A1",
        )
        btn_cmd.clicked.connect(lambda: _launch(_CMD_MAIN, _LEGACY_CMD))

        btn_lc = ModeButton(
            "시계열분석", "Light Curve Analysis",
            "#2E7D32", "#388E3C", "#1B5E20",
        )
        btn_lc.clicked.connect(lambda: _launch(_LC_MAIN, _LEGACY_LC))

        btn_row = QHBoxLayout()
        btn_row.setSpacing(20)
        btn_row.addWidget(btn_cmd)
        btn_row.addWidget(btn_lc)
        root.addLayout(btn_row)

        # ── Footer ──────────────────────────────────────────────────────────
        footer = QLabel("APEX v0.1.0")
        footer.setAlignment(Qt.AlignCenter)
        footer.setStyleSheet("color: #aaaaaa; font-size: 9px;")
        root.addSpacing(4)
        root.addWidget(footer)


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("APEX")
    app.setOrganizationName("APEX Project")
    _configure_app_fonts(app)
    os.chdir(_HERE)

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
