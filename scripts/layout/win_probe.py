"""Report one APEX window's minimum size. Run one per subprocess so a window
that pops a modal or blocks can't take the whole scan down.

    python win_probe.py <mode> <module_path> <ClassName>
"""
import sys
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))

MODE, MODULE, CLASS = sys.argv[1], sys.argv[2], sys.argv[3]

from PyQt5.QtWidgets import QApplication, QMessageBox, QFileDialog

BASE = (Path(r"E:\APEX_validation\reprocess\NGC6811") if MODE == "cmd"
        else Path(r"E:\APEX_validation\reprocess\YZBoo_2n"))
SCI = BASE / "sci"
_FRAMES = sorted(p.name for p in SCI.glob("*.fit*"))


class _FMProxy:
    filenames = _FRAMES
    path_map = {n: SCI / n for n in _FRAMES}
    night_assignments: dict = {}
    excluded_nights: set = set()
    excluded_files: set = set()
    def get_file_path(self, f): return SCI / f
    def get_file_list(self): return list(_FRAMES)


class _MainStub:
    step_names = ["step"] * 12
    mode = MODE
    def on_step_completed(self, *a, **k): pass
    def open_step(self, *a, **k): pass
    def update_step_buttons(self, *a, **k): pass


app = QApplication.instance() or QApplication([])
from apex.utils.app_setup import configure_fonts
from apex.gui.theme import apply_theme
configure_fonts(app)
apply_theme(app, "apex-light")

# Nothing may block: a modal or a file dialog on a half-built window would hang.
for _n in ("critical", "warning", "information", "question"):
    setattr(QMessageBox, _n, staticmethod(lambda *a, **k: None))
for _n in ("getOpenFileName", "getSaveFileName", "getExistingDirectory"):
    setattr(QFileDialog, _n, staticmethod(lambda *a, **k: ("", "")))

if MODE == "cmd":
    from apex.config.parameters_cmd import read_params
else:
    from apex.config.parameters_lc import read_params

import importlib
mod = importlib.import_module(MODULE)
cls = getattr(mod, CLASS)

params = read_params(str(BASE / "parameters.toml"))
params.P.result_dir = BASE / "result"

from apex.core.project_state import ProjectState
win = cls(params, _FMProxy(), ProjectState(BASE / "result"), _MainStub())
win.show()
app.processEvents()
m = win.minimumSizeHint()
print(f"RESULT\t{CLASS}\t{m.width()}\t{m.height()}")
win.close()
