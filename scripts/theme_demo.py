"""Render the SAME tool panel in two styles and save a before/after PNG.

  LEFT  = current APEX look (inline setStyleSheet, material-2014 palette,
          nested group boxes, emoji status, courier log)
  RIGHT = apex/gui/theme.py global QSS (tokens, flat cards, one accent,
          typographic hierarchy, quiet log) — note: almost no per-widget
          styling; the look comes from the global stylesheet.

Output: ui_screenshots/theme_before_after.png
Run:    python scripts/theme_demo.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
os.chdir(REPO)

from PyQt5.QtWidgets import (
    QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton,
    QGroupBox, QFormLayout, QLineEdit, QComboBox, QTextEdit, QFrame,
)
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont

from apex.gui.theme import apply_theme, refresh

OUT = REPO / "ui_screenshots"
OUT.mkdir(exist_ok=True)


# ── BEFORE: the current style, faithfully reproduced ────────────────────────
def build_before() -> QWidget:
    w = QWidget()
    w.setStyleSheet("background: #ECECEC;")  # default Qt window grey
    root = QVBoxLayout(w)
    root.setContentsMargins(10, 10, 10, 10)
    root.setSpacing(8)

    title = QLabel("Extinction (Airmass Fit)")
    title.setFont(QFont("Arial", 16, QFont.Bold))
    title.setAlignment(Qt.AlignCenter)
    root.addWidget(title)

    info = QLabel("Run per-star Bouguer extinction fitting on the cached Step 7 table.")
    info.setStyleSheet("QLabel { background-color: #E3F2FD; padding: 8px; border-radius: 5px; }")
    root.addWidget(info)

    settings = QGroupBox("Tool")
    form = QFormLayout(settings)
    form.addRow("Filter:", QComboBox())
    form.addRow("SNR min:", QLineEdit("20"))
    form.addRow("Workspace:", QLineEdit("data/example/result"))
    root.addWidget(settings)

    status = QGroupBox("Status")
    sl = QVBoxLayout(status)
    ok = QLabel("✓ Step 7 loaded")
    ok.setStyleSheet("QLabel { color: #4CAF50; font-weight: bold; }")
    warn = QLabel("⚠ 3 frames missing airmass")
    warn.setStyleSheet("QLabel { color: #FF9800; font-weight: bold; }")
    sl.addWidget(ok)
    sl.addWidget(warn)
    root.addWidget(status)

    btn_row = QHBoxLayout()
    run = QPushButton("Run Fit")
    run.setStyleSheet("""
        QPushButton { background-color: #4CAF50; color: white; font-weight: bold;
                      border: 2px solid #45a049; border-radius: 5px; padding: 8px 16px; }
        QPushButton:hover { background-color: #45a049; }
    """)
    stop = QPushButton("Stop")
    stop.setStyleSheet("""
        QPushButton { background-color: #f44336; color: white; font-weight: bold;
                      border: 2px solid #d32f2f; border-radius: 5px; padding: 8px 16px; }
    """)
    btn_row.addWidget(run)
    btn_row.addWidget(stop)
    btn_row.addStretch()
    root.addLayout(btn_row)

    log_group = QGroupBox("Log")
    lg = QVBoxLayout(log_group)
    log = QTextEdit()
    log.setReadOnly(True)
    log.setStyleSheet("QTextEdit { font-family: monospace; font-size: 9pt; }")
    log.setPlainText("[FIT] loading photometry_index.csv\n[FIT] 142 stars, 38 frames\n[FIT] k1 = 0.184 +/- 0.012")
    log.setFixedHeight(90)
    lg.addWidget(log)
    root.addWidget(log_group)
    root.addStretch()
    return w


# ── AFTER: same content, styled almost entirely by the global QSS ───────────
def _card(title_text: str) -> tuple[QFrame, QVBoxLayout]:
    card = QFrame()
    card.setObjectName("Card")
    lay = QVBoxLayout(card)
    lay.setContentsMargins(16, 14, 16, 14)
    lay.setSpacing(10)
    if title_text:
        h = QLabel(title_text)
        h.setProperty("role", "subtitle")
        lay.addWidget(h)
    return card, lay


def build_after() -> QWidget:
    w = QWidget()
    root = QVBoxLayout(w)
    root.setContentsMargins(24, 24, 24, 24)
    root.setSpacing(16)

    title = QLabel("Extinction · Airmass Fit")
    title.setProperty("role", "title")
    root.addWidget(title)
    sub = QLabel("Per-star Bouguer fit over the cached Step 7 photometry table.")
    sub.setProperty("role", "subtitle")
    root.addWidget(sub)

    settings, sl = _card("Parameters")
    form = QFormLayout()
    form.setSpacing(10)
    form.addRow("Filter", QComboBox())
    snr = QLineEdit("20")
    form.addRow("SNR min", snr)
    form.addRow("Workspace", QLineEdit("data/example/result"))
    sl.addLayout(form)
    root.addWidget(settings)

    status, stl = _card("Status")
    ok = QLabel("Step 7 loaded — 142 stars, 38 frames")
    ok.setProperty("status", "ok")
    warn = QLabel("3 frames missing airmass")
    warn.setProperty("status", "warn")
    stl.addWidget(ok)
    stl.addWidget(warn)
    root.addWidget(status)

    btn_row = QHBoxLayout()
    btn_row.setSpacing(8)
    run = QPushButton("Run fit")
    run.setProperty("variant", "primary")
    stop = QPushButton("Stop")
    stop.setProperty("variant", "ghost")
    btn_row.addWidget(run)
    btn_row.addWidget(stop)
    btn_row.addStretch()
    root.addLayout(btn_row)

    log_card, lcl = _card("Log")
    log = QTextEdit()
    log.setObjectName("Log")
    log.setReadOnly(True)
    log.setPlainText("loading photometry_index.csv\n142 stars, 38 frames\nk1 = 0.184 ± 0.012")
    log.setFixedHeight(90)
    lcl.addWidget(log)
    root.addWidget(log_card)
    root.addStretch()

    for lbl in (title, sub, ok, warn):
        refresh(lbl)
    return w


def main() -> int:
    app = QApplication.instance() or QApplication(sys.argv)

    before = build_before()
    before.resize(420, 540)
    before.show()
    for _ in range(8):
        app.processEvents()
    before.grab().save(str(OUT / "theme_before.png"), "PNG")

    # apply the global theme only for the "after" capture
    apply_theme(app)
    after = build_after()
    after.resize(420, 540)
    after.show()
    for _ in range(8):
        app.processEvents()
    after.grab().save(str(OUT / "theme_after.png"), "PNG")

    print(f"[ok] saved theme_before.png / theme_after.png -> {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
