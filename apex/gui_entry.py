"""GUI launcher entry point usable as a console script (``apex-gui``).

The full launcher bootstrap (font configuration, application icon, first-run
parameter creation) currently lives in the repository-root ``main.py``. This
module loads that file and delegates to it so there is exactly one launcher
implementation. When the launcher is eventually relocated into the package
(e.g. ``apex/app.py``), this shim should import it directly instead.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from typing import Optional


def _load_root_main():
    """Import the repository-root main.py as a module, regardless of cwd."""
    if getattr(sys, "frozen", False):
        root = Path(sys.executable).parent
    else:
        root = Path(__file__).resolve().parent.parent
    main_py = root / "main.py"
    if not main_py.exists():
        raise FileNotFoundError(f"Launcher not found: {main_py}")
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    spec = importlib.util.spec_from_file_location("apex_root_main", main_py)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load launcher spec from {main_py}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main(mode: Optional[str] = None) -> int:
    """Launch the GUI.

    Args:
        mode: ``"cmd"`` or ``"lc"`` to open a mode directly; ``None`` shows the
              launcher window so the user can choose.
    """
    root_main = _load_root_main()

    if mode is None:
        return int(root_main.main())

    # Direct-mode launch reuses the launcher's bootstrap helpers.
    from PyQt5.QtWidgets import QApplication
    from PyQt5.QtCore import Qt

    app = QApplication(sys.argv)
    app.setQuitOnLastWindowClosed(True)
    app.setApplicationName("APEX")
    app.setOrganizationName("APEX Project")

    param_error = root_main._ensure_parameters()
    if param_error:
        from PyQt5.QtWidgets import QMessageBox

        QMessageBox.critical(None, "Startup Error", param_error)
        return 1
    root_main._configure_app_fonts(app)
    app.setWindowIcon(root_main._make_app_icon())

    from apex.gui.main_window import MainWindowWorkflow

    window = MainWindowWorkflow(mode=mode)
    window.setAttribute(Qt.WA_DeleteOnClose)
    window.show()
    return int(app.exec_())


if __name__ == "__main__":
    raise SystemExit(main())
