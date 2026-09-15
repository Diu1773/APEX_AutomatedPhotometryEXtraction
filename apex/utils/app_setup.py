"""Shared PyQt5 / matplotlib font and application setup."""
from __future__ import annotations

import os
from pathlib import Path

from PyQt5.QtWidgets import QApplication
from PyQt5.QtGui import QFont, QFontDatabase

_CANDIDATE_FONT_FILES = [
    Path(os.environ.get("WINDIR", "C:/Windows")) / "Fonts" / "malgun.ttf",
    Path("/mnt/c/Windows/Fonts/malgun.ttf"),
]
_CANDIDATE_FAMILIES = [
    "Malgun Gothic", "맑은 고딕", "NanumGothic", "Nanum Gothic",
    "Noto Sans CJK KR", "Noto Sans KR", "AppleGothic", "Arial Unicode MS",
]


def configure_fonts(app: QApplication) -> None:
    """Configure Korean fonts for matplotlib and Qt."""
    resolved_family = None
    try:
        import matplotlib
        from matplotlib import font_manager
        for fp in _CANDIDATE_FONT_FILES:
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
            for fam in _CANDIDATE_FAMILIES:
                if fam.lower() in available:
                    resolved_family = available[fam.lower()]
                    break
        if resolved_family:
            matplotlib.rcParams.update({
                "font.family": [resolved_family, "DejaVu Sans"],
                "axes.unicode_minus": False,
                "mathtext.fontset": "dejavusans",
                "mathtext.default": "regular",
                "axes.formatter.use_mathtext": False,
            })
    except Exception:
        pass
    try:
        db = QFontDatabase()
        qt_fams = {f.lower(): f for f in db.families()}
        qt_fam = resolved_family and qt_fams.get(resolved_family.lower())
        if qt_fam is None:
            for fam in _CANDIDATE_FAMILIES:
                qt_fam = qt_fams.get(fam.lower())
                if qt_fam:
                    break
        if qt_fam:
            app.setFont(QFont(qt_fam, 9))
    except Exception:
        pass


def install_slot_exception_guard() -> None:
    """Keep an unhandled exception in a Qt slot from killing the app silently.

    PyQt5 calls ``qFatal()`` on an unhandled Python exception raised inside a
    slot, so the process just vanishes — no traceback, no log line, and
    ``faulthandler`` sees nothing because it is not a signal. Two real bugs hid
    behind that in 2026-09-15 (a missing import and two leftover lines, both
    plain ``NameError``), and Windows reported them as 0xC0000409.

    Setting ``sys.excepthook`` is enough: PyQt calls it and skips the abort.
    """
    import logging
    import sys
    import traceback

    def _hook(kind, value, tb):
        # Record only. Putting a QMessageBox here kills the process just as the
        # abort did — the hook runs while the slot is still on the stack, and a
        # modal from there does not survive (measured 2026-09-15).
        text = "".join(traceback.format_exception(kind, value, tb))
        logging.getLogger("apex").error("unhandled exception\n%s", text)
        sys.stderr.write(text)

    sys.excepthook = _hook
