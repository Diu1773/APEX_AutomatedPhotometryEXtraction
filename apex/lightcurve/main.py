#!/usr/bin/env python3
"""APEX LC — Light Curve Analysis entry point."""

import sys
import os
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

# Ensure project root is on sys.path so `import apex` works
_HERE = Path(__file__).resolve().parent          # apex/lightcurve/
_ROOT = _HERE.parent.parent                       # Automated_Photometry_EXtraction/
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from PyQt5.QtWidgets import QApplication, QMessageBox
from PyQt5.QtGui import QFont, QFontDatabase


def _configure_fonts(app: QApplication) -> None:
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
        for fp in candidate_font_files:
            if not fp.exists():
                continue
            try:
                font_manager.fontManager.addfont(str(fp))
                resolved_family = font_manager.FontProperties(fname=str(fp)).get_name()
                break
            except Exception:
                continue
        if resolved_family is None:
            available = {f.name.lower(): f.name for f in font_manager.fontManager.ttflist}
            for fam in candidate_families:
                if fam.lower() in available:
                    resolved_family = available[fam.lower()]
                    break
        if resolved_family:
            matplotlib.rcParams.update({
                "font.family": [resolved_family, "DejaVu Sans"],
                "axes.unicode_minus": False,
            })
    except Exception:
        pass
    try:
        db = QFontDatabase()
        qt_fams = {f.lower(): f for f in db.families()}
        qt_fam = resolved_family and qt_fams.get(resolved_family.lower())
        if qt_fam is None:
            for fam in candidate_families:
                qt_fam = qt_fams.get(fam.lower())
                if qt_fam:
                    break
        if qt_fam:
            app.setFont(QFont(qt_fam, 9))
    except Exception:
        pass


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("APEX LC")
    app.setOrganizationName("APEX Project")
    _configure_fonts(app)
    os.chdir(_HERE)

    try:
        from apex.lightcurve.gui.main_window_workflow import MainWindowWorkflow
        window = MainWindowWorkflow()
        window.show()
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"Failed to start APEX LC: {e}")
        print(tb)
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setWindowTitle("APEX LC — Startup Error")
        msg.setText(str(e))
        msg.setDetailedText(tb)
        msg.exec_()
        return 1

    return app.exec_()


if __name__ == "__main__":
    sys.exit(main())
