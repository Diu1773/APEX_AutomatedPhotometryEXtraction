#!/usr/bin/env python3
"""APEX CMD — Cluster Photometry entry point."""
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

_HERE = Path(__file__).resolve().parent
_ROOT = _HERE.parent.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from PyQt5.QtWidgets import QApplication, QMessageBox
from apex.utils.app_setup import configure_fonts


def main() -> int:
    app = QApplication(sys.argv)
    app.setApplicationName("APEX CMD")
    app.setOrganizationName("APEX Project")
    configure_fonts(app)
    os.chdir(_HERE)
    try:
        from apex.gui.main_window import MainWindowWorkflow
        window = MainWindowWorkflow(mode="cmd")
        window.show()
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        print(f"Failed to start APEX CMD: {{e}}")
        print(tb)
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
