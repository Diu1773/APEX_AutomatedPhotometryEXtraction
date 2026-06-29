"""Capture Tools-menu window screenshots fully offscreen (no window on screen).

This renders each standalone analysis tool window to a PNG without ever showing
anything on the user's display, by running Qt with the *offscreen* platform
plugin. Two environment settings are essential:

    QT_QPA_PLATFORM=offscreen     # render to memory, not the screen
    QT_QPA_FONTDIR=C:/Windows/Fonts   # offscreen has an EMPTY font DB otherwise,
                                       # so every Qt label/button would be blank

The tool windows are populated from whatever workspace ``parameters.toml`` points
at (``[io] result_dir``). Plot-only content (matplotlib) renders even without
this font fix; native Qt widget text does not — hence QT_QPA_FONTDIR.

Usage (PowerShell):

    $env:QT_QPA_PLATFORM="offscreen"; $env:QT_QPA_FONTDIR="C:/Windows/Fonts"
    python scripts/capture_tool_screenshots.py

Output: ui_screenshots/tool_*.png
"""
from __future__ import annotations

import os
import sys
import time
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_QPA_FONTDIR", "C:/Windows/Fonts")
os.environ.setdefault("QT_OPENGL", "software")

ROOT = Path(__file__).resolve().parent.parent
os.chdir(ROOT)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from PyQt5.QtWidgets import QApplication, QMessageBox, QPushButton  # noqa: E402
from PyQt5.QtGui import QFont  # noqa: E402

# Never let a modal dialog block this headless run.
QMessageBox.exec_ = lambda *a, **k: None  # type: ignore[assignment]
for _m in ("information", "warning", "critical", "question", "about"):
    setattr(QMessageBox, _m, staticmethod(lambda *a, **k: None))

OUT = ROOT / "ui_screenshots"


def main() -> int:
    app = QApplication(sys.argv)
    app.setFont(QFont("Segoe UI", 9))
    from apex.gui.theme import apply_theme
    apply_theme(app)
    from apex.gui.main_window import MainWindowWorkflow
    win = MainWindowWorkflow(mode="cmd")
    print("result_dir =", win.params.P.result_dir)

    def pump(n: int = 30) -> None:
        for _ in range(n):
            app.processEvents()

    def wait(sec: float) -> None:
        t0 = time.time()
        while time.time() - t0 < sec:
            app.processEvents()
            time.sleep(0.05)

    def click(wobj, text: str) -> bool:
        for b in wobj.findChildren(QPushButton):
            if b.text().strip() == text:
                b.click()
                return True
        return False

    def cap(attr: str, fname: str, w: int, h: int) -> None:
        wobj = getattr(win, attr, None)
        if wobj is None:
            print("NO_WINDOW", attr)
            return
        wobj.setMinimumSize(0, 0)
        wobj.resize(w, h)
        pump()
        wobj.resize(w, h)
        pump()
        ok = wobj.grab().save(str(OUT / fname))
        print("SAVED" if ok else "FAIL", fname, wobj.width(), "x", wobj.height())

    # (launcher method, window attr, filename, width, height, optional driver)
    win.open_qa_report(); pump(30)
    if click(win.qa_window, "Generate QA Report"):
        wait(30)
    cap("qa_window", "tool_qa_report.png", 1280, 880)

    win.open_gaia_3d_viewer(); pump(40); cap("gaia_3d_window", "tool_gaia_3d.png", 1200, 860)
    win.open_extinction_tool(); pump(40); cap("extinction_window", "tool_extinction.png", 1320, 900)
    win.open_iraf_tool(); pump(40); cap("iraf_window", "tool_iraf.png", 1200, 860)
    win.open_cluster_structure_tool(); pump(40); cap("cluster_structure_window", "tool_cluster_structure.png", 1200, 860)
    win.open_airmass_debug_tool(); pump(40); cap("airmass_debug_window", "tool_airmass_debug.png", 1100, 820)
    win.open_variable_star_tool(); pump(40); cap("varstar_window", "tool_variable_star.png", 1280, 860)
    win.open_transit_tool(); pump(40); cap("transit_window", "tool_transit.png", 1240, 1000)
    win.open_eb_tool(); pump(40); cap("eb_window", "tool_eb.png", 1240, 980)
    win.open_multi_night_merger(); pump(40); cap("merger_window", "tool_merger.png", 1200, 860)

    print("DONE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
