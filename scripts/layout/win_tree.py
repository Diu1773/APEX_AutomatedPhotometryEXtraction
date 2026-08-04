"""Show which children drive one window's minimum size.

    python win_tree.py <mode> <module_path> <ClassName> [min_h]
"""
import sys
from pathlib import Path

REPO = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction")
sys.path.insert(0, str(REPO))

MODE, MODULE, CLASS = sys.argv[1], sys.argv[2], sys.argv[3]
MIN_H = int(sys.argv[4]) if len(sys.argv) > 4 else 300

from PyQt5.QtWidgets import (QApplication, QMessageBox, QFileDialog, QWidget,
                             QTabWidget, QScrollArea)

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
for _n in ("critical", "warning", "information", "question"):
    setattr(QMessageBox, _n, staticmethod(lambda *a, **k: None))
for _n in ("getOpenFileName", "getSaveFileName", "getExistingDirectory"):
    setattr(QFileDialog, _n, staticmethod(lambda *a, **k: ("", "")))

if MODE == "cmd":
    from apex.config.parameters_cmd import read_params
else:
    from apex.config.parameters_lc import read_params

import importlib
cls = getattr(importlib.import_module(MODULE), CLASS)
params = read_params(str(BASE / "parameters.toml"))
params.P.result_dir = BASE / "result"

from apex.core.project_state import ProjectState
win = cls(params, _FMProxy(), ProjectState(BASE / "result"), _MainStub())
win.show(); app.processEvents()
m = win.minimumSizeHint()
print(f"{CLASS}: minimumSizeHint = {m.width()} x {m.height()}\n")


def label(w):
    for attr in ("title", "text", "windowTitle"):
        if hasattr(w, attr):
            try:
                t = str(getattr(w, attr)())[:40].replace("\n", " ")
            except Exception:
                t = ""
            if t:
                return t
    return w.objectName() or ""


def tree(w, depth=0, limit=5):
    if depth > limit:
        return
    for c in w.children():
        if not isinstance(c, QWidget):
            continue
        h = c.minimumSizeHint().height()
        wd = c.minimumSizeHint().width()
        if h < MIN_H and wd < 700:
            continue
        mark = ""
        if isinstance(c, QScrollArea):
            mark = "  [scrolled]"
        print(f"{'  ' * depth}{c.__class__.__name__:16} {wd:4d}x{h:<4d} {label(c)}{mark}")
        if isinstance(c, QTabWidget):
            for i in range(c.count()):
                pg = c.widget(i)
                s = "  [scrolled]" if isinstance(pg, QScrollArea) else ""
                print(f"{'  ' * (depth+1)}tab{i} {c.tabText(i)[:24]:24} "
                      f"{pg.minimumSizeHint().width():4d}x{pg.minimumSizeHint().height():<4d}{s}")
        tree(c, depth + 1, limit)


tree(win)
win.close()
