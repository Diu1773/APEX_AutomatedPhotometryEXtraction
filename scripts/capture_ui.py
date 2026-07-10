"""Capture PNG screenshots of every APEX step window for UI/UX review.

Outputs to ui_screenshots/ at the repo root. Skips steps that fail to
construct (logs the exception) so a single broken window doesn't abort
the run.
"""
from __future__ import annotations

import os
import re
import sys
import traceback
from pathlib import Path

# Offscreen render (no window flash on the user's display) — same recipe as
# capture_tool_screenshots.py. FONTDIR is required or all Qt text is blank.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
os.environ.setdefault("QT_QPA_FONTDIR", "C:/Windows/Fonts")
os.environ.setdefault("QT_OPENGL", "software")

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
os.chdir(REPO)

from PyQt5.QtWidgets import QApplication, QMessageBox
from PyQt5.QtCore import Qt

# Never let a modal dialog block (or crash) this headless run — same recipe
# as capture_tool_screenshots.py.
QMessageBox.exec_ = lambda *a, **k: None  # type: ignore[assignment]
for _m in ("information", "warning", "critical", "question", "about"):
    setattr(QMessageBox, _m, staticmethod(lambda *a, **k: None))

OUT = REPO / "ui_screenshots"
OUT.mkdir(exist_ok=True)


def _slug(text: str) -> str:
    s = text.lower().replace("&", "and").replace("/", "_")
    s = re.sub(r"[^a-z0-9]+", "_", s).strip("_")
    return s or "window"


def _pump(n: int = 6) -> None:
    app = QApplication.instance()
    for _ in range(n):
        app.processEvents()


def _grab(widget, path: Path) -> None:
    widget.show()
    widget.raise_()
    _pump(8)
    pix = widget.grab()
    pix.save(str(path), "PNG")
    print(f"[ok]   {path.name}  ({pix.width()}x{pix.height()})")


def _capture_mode(mode: str, param_file: str | None = None) -> None:
    print(f"\n=== mode: {mode}  (params: {param_file or 'parameters.toml'}) ===")
    from apex.gui.main_window import MainWindowWorkflow

    mw = MainWindowWorkflow(mode=mode, param_file=param_file)
    mw.resize(1100, 900)
    mw.show()
    _pump(10)
    _grab(mw, OUT / f"00_main_{mode}.png")

    # Optional filter: APEX_CAPTURE_STEPS="8,9,11" captures only those step
    # numbers (1-based). A native crash in one step window (e.g. an image
    # viewer under the offscreen platform) kills the whole process, so the
    # caller can loop steps in separate processes.
    only_steps = None
    steps_env = os.environ.get("APEX_CAPTURE_STEPS", "").strip()
    if steps_env:
        only_steps = {int(s) for s in steps_env.split(",") if s.strip()}

    for idx, name in enumerate(mw.step_names):
        if only_steps is not None and (idx + 1) not in only_steps:
            continue
        slug = _slug(name)
        out = OUT / f"{mode}_step{idx + 1:02d}_{slug}.png"
        try:
            win = mw._open_step_window(idx)
            if win is None:
                print(f"[skip] {mode} step {idx + 1} ({name}) — no window")
                continue
            try:
                win.resize(1400, 900)
            except Exception:
                pass
            mw.current_step_window = win
            _grab(win, out)
        except Exception as e:
            print(f"[FAIL] {mode} step{idx + 1} ({name}): {type(e).__name__}: {e}")
            traceback.print_exc(limit=3)
        finally:
            try:
                if mw.current_step_window is not None:
                    mw.current_step_window.close()
                    mw.current_step_window = None
            except Exception:
                pass
            _pump(4)

    mw.close()
    _pump(4)


def main() -> int:
    app = QApplication.instance() or QApplication(sys.argv)
    app.setAttribute(Qt.AA_DontUseNativeMenuBar, True)
    # Match the real entry points: without the global QSS the capture shows a
    # different (unthemed) UI than users actually see.
    from apex.gui.theme import apply_theme
    apply_theme(app)
    mode_params = {
        "cmd": "parameters_capture_cmd.toml" if (REPO / "parameters_capture_cmd.toml").exists() else None,
        "lc": None,
    }
    only = os.environ.get("APEX_CAPTURE_MODE")
    modes = (only,) if only in ("cmd", "lc") else ("cmd", "lc")
    for mode in modes:
        try:
            _capture_mode(mode, mode_params.get(mode))
        except Exception as e:
            print(f"[MODE FAIL] {mode}: {type(e).__name__}: {e}")
            traceback.print_exc()
    print(f"\n→ saved to {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
