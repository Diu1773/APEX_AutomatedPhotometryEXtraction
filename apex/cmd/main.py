#!/usr/bin/env python3
"""APEX CMD entry point."""
import sys, os, warnings
from pathlib import Path

warnings.filterwarnings("ignore", message=".*datfix.*MJD-OBS.*", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning, module="astropy")
try:
    from astropy.wcs import FITSFixedWarning
    warnings.filterwarnings("ignore", category=FITSFixedWarning)
except ImportError:
    pass
warnings.filterwarnings("ignore", message=".*tight_layout.*", category=UserWarning)

_HERE = Path(__file__).resolve().parent   # apex/cmd/
_ROOT = _HERE.parent.parent               # Automated_Photometry_EXtraction/
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

os.chdir(_ROOT)  # parameters.toml lives at project root

from PyQt5.QtWidgets import QApplication, QMessageBox
from PyQt5.QtCore import Qt
from apex.utils.app_setup import configure_fonts


def main() -> int:
    # Native HiDPI rendering (see main.py): avoids fuzzy text / pixelated logo
    # under Windows display scaling.
    QApplication.setAttribute(Qt.AA_EnableHighDpiScaling, True)
    QApplication.setAttribute(Qt.AA_UseHighDpiPixmaps, True)
    app = QApplication(sys.argv)
    app.setApplicationName("APEX CMD")
    app.setOrganizationName("APEX Project")
    configure_fonts(app)
    from apex.gui.theme import apply_theme
    apply_theme(app)
    try:
        from apex.gui.main_window import _load_icon
        app.setWindowIcon(_load_icon("cmd"))
    except Exception:
        pass
    try:
        from apex.gui.main_window import MainWindowWorkflow
        window = MainWindowWorkflow(mode="cmd")
        window.show()
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"Failed to start APEX CMD: {e}\n{tb}")
        msg = QMessageBox()
        msg.setIcon(QMessageBox.Critical)
        msg.setWindowTitle("APEX CMD — Startup Error")
        msg.setText(str(e))
        msg.setDetailedText(tb)
        msg.exec_()
        return 1
    return app.exec_()

if __name__ == "__main__":
    sys.exit(main())
