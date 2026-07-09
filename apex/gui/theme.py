"""APEX design tokens + global Qt stylesheet.

The point of this module: stop hard-coding hex colours and borders on every
widget. Define the palette/spacing/typography once here, ship a single global
QSS, and let windows opt into roles via objectName / dynamic properties.

Usage:
    from apex.gui.theme import apply_theme
    apply_theme(QApplication.instance())

Role hooks used by the QSS (set with setProperty or objectName):
    QLabel  property "role" in {"title", "subtitle", "caption", "info"}
    QPushButton property "variant" in {"primary", "danger", "ghost"}
    QFrame  objectName "Card"           → flat surface, no heavy border
    QLabel  property "status" in {"ok", "warn", "error", "idle"}
"""

from __future__ import annotations


# ── Design tokens ──────────────────────────────────────────────────────────
class Tokens:
    # Surfaces — the canvas (BG) is deliberately a few shades below white so
    # white cards/panels/inputs read as raised "cards on a canvas" and their
    # edges stay visible even with thin borders.
    BG          = "#E9ECF1"   # window background (canvas)
    SURFACE     = "#FFFFFF"   # cards / panels
    SURFACE_ALT = "#F1F3F5"   # subtle fills (inputs, log)

    # Hairlines — strong enough to actually read as boundaries on the canvas.
    BORDER      = "#CDD4DE"
    BORDER_STRONG = "#AEB8C6"

    # Text
    TEXT        = "#1F2933"   # primary
    TEXT_SUB    = "#5B6573"   # secondary
    TEXT_MUTED  = "#98A2B3"   # captions / disabled

    # Single calm accent (replaces the rainbow of material primaries)
    ACCENT      = "#3A66DB"
    ACCENT_HOVER = "#2F56C0"
    ACCENT_PRESS = "#274AA6"
    ACCENT_SOFT = "#EAF0FE"   # accent-tinted background for info chips

    # Semantic, desaturated vs. material 2014
    OK          = "#2E9E5B"
    WARN        = "#C77A12"
    ERROR       = "#D24343"

    # Muted variant fills (for disabled primary/danger so they don't look active)
    ACCENT_MUTED = "#C2CCEA"
    ERROR_MUTED  = "#E6C5C5"

    # Geometry
    RADIUS      = 8
    RADIUS_SM   = 6

    # Spacing scale
    S1, S2, S3, S4, S5 = 4, 8, 12, 16, 24

    # Button height scale — one source of truth for sizing (the other half of
    # "hierarchy": the main action row is tall, utility buttons are compact).
    H_ACTION  = 38   # bottom action row: Run / Stop / Previous / Next
    H_BUTTON  = 32   # standard buttons (Browse, Export, Fit, …)
    H_COMPACT = 28   # header cluster: Parameters / Log / 가이드

    # Layout rhythm: everything (margins, spacing) snaps to this 8px grid so
    # buttons and rows line up as one set — gaps of 8/12/16, never 5/6/10.
    GRID        = 8
    GAP         = 8    # default spacing between sibling controls
    MARGIN      = 16   # default window/panel content margin

    # Typography
    FONT_STACK  = '"Segoe UI", "Malgun Gothic", "Noto Sans KR", Arial, sans-serif'
    FS_TITLE    = 20
    FS_SUBTITLE = 13
    FS_BASE     = 13
    FS_CAPTION  = 11


# ── Icon glyphs ──────────────────────────────────────────────────────────────
# APEX labels its action buttons with glyphs, not a bundled vector icon set.
# Bare emoji (⚙ 📜 🔒 …) render as multicolor OS emoji on Windows and break the
# flat monochrome look. Appending U+FE0E (VARIATION SELECTOR-15 / "text
# presentation") asks the font for the flat symbol glyph instead — and it is a
# zero-width, no-op when unsupported, so it can never make things worse.
#
# This is the single source for button glyphs: reference ICON["params"] etc.
# instead of pasting an emoji literal, so the icon set stays one consistent
# family (the keyline-grid idea applied to a glyph-based UI).
_FE0E = chr(0xFE0E)  # VARIATION SELECTOR-15: forces flat/text glyph presentation


def _mono(glyph: str) -> str:
    """Force a flat (monochrome) presentation of an emoji-capable *glyph*."""
    return glyph + _FE0E


ICON = {
    "params": _mono("⚙"),
    "log":    _mono("📜"),
    "guide":  "ⓘ",          # circled-i is monochrome by default
    "input":  _mono("📂"),
    "output": _mono("💾"),
    "locked": _mono("🔒"),
    "done":   "✓",
    "todo":   "○",
    "exit":   "✕",
    "prev":   "←",
    "next":   "→",
}


def global_qss(t: type[Tokens] = Tokens) -> str:
    """Return the application-wide stylesheet built from the tokens."""
    return f"""
    /* NOTE: no font-family / base font-size here on purpose. The app font is
       owned by configure_fonts()/_configure_app_fonts() (resolves Malgun
       Gothic etc. per system); setting it in QSS would override that and force
       a different family/size. The theme only governs color + shape. */
    QWidget {{
        background: {t.BG};
        color: {t.TEXT};
    }}

    /* Labels/checks are transparent so they don't paint a BG-coloured
       rectangle over the white card surface they sit on. */
    QLabel, QCheckBox, QRadioButton {{ background: transparent; }}

    /* Cards replace boxy nested QGroupBoxes: a flat surface + hairline + radius */
    QFrame#Card {{
        background: {t.SURFACE};
        border: 1px solid {t.BORDER};
        border-radius: {t.RADIUS}px;
    }}

    /* Tame the group box: no double frame, label sits as a quiet section header */
    QGroupBox {{
        background: {t.SURFACE};
        border: 1px solid {t.BORDER};
        border-radius: {t.RADIUS}px;
        margin-top: 14px;
        padding: {t.S4}px;
        font-weight: 600;
    }}
    QGroupBox::title {{
        subcontrol-origin: margin;
        left: {t.S3}px;
        padding: 0 4px;
        color: {t.TEXT_SUB};
    }}

    /* Typography roles */
    QLabel[role="title"]    {{ font-size: {t.FS_TITLE}px; font-weight: 700; color: {t.TEXT}; }}
    QLabel[role="subtitle"] {{ font-size: {t.FS_SUBTITLE}px; color: {t.TEXT_SUB}; }}
    QLabel[role="caption"]  {{ font-size: {t.FS_CAPTION}px; color: {t.TEXT_MUTED}; }}
    QLabel[role="info"] {{
        background: {t.ACCENT_SOFT};
        color: {t.ACCENT_PRESS};
        border-radius: {t.RADIUS_SM}px;
        padding: {t.S2}px {t.S3}px;
    }}

    /* Status pills */
    QLabel[status="ok"]    {{ color: {t.OK};    font-weight: 600; }}
    QLabel[status="warn"]  {{ color: {t.WARN};  font-weight: 600; }}
    QLabel[status="error"] {{ color: {t.ERROR}; font-weight: 600; }}
    QLabel[status="idle"]  {{ color: {t.TEXT_MUTED}; }}

    /* Buttons — ONE hierarchy, applied everywhere:
         (no variant)      → secondary / neutral, the quiet default
         variant=primary   → the single main action (Run, Next, Save)
         variant=danger    → Stop / destructive (Reset, Delete)
         variant=ghost     → tertiary inline (Log, 가이드, links)
       Set the role with style_button(btn, "primary") or
       btn.setProperty("variant", "primary"); never hand-paint hex on a button. */
    QPushButton {{
        background: {t.SURFACE};
        color: {t.TEXT};
        border: 1px solid {t.BORDER_STRONG};
        border-radius: {t.RADIUS_SM}px;
        padding: 6px 14px;
        font-weight: 500;
    }}
    QPushButton:hover  {{ background: {t.SURFACE_ALT}; border-color: {t.TEXT_MUTED}; }}
    QPushButton:pressed {{ background: {t.BORDER}; }}
    QPushButton:disabled {{ background: {t.SURFACE}; color: {t.TEXT_MUTED}; border-color: {t.BORDER}; }}

    QPushButton[variant="primary"] {{
        background: {t.ACCENT}; color: #FFFFFF; border: none; font-weight: 600;
    }}
    QPushButton[variant="primary"]:hover   {{ background: {t.ACCENT_HOVER}; }}
    QPushButton[variant="primary"]:pressed {{ background: {t.ACCENT_PRESS}; }}
    QPushButton[variant="primary"]:disabled {{ background: {t.ACCENT_MUTED}; color: #FFFFFF; }}

    QPushButton[variant="danger"] {{
        background: {t.ERROR}; color: #FFFFFF; border: none; font-weight: 600;
    }}
    QPushButton[variant="danger"]:hover   {{ background: #BE3A3A; }}
    QPushButton[variant="danger"]:pressed {{ background: #A83333; }}
    QPushButton[variant="danger"]:disabled {{ background: {t.ERROR_MUTED}; color: #FFFFFF; }}

    QPushButton[variant="success"] {{
        background: {t.OK}; color: #FFFFFF; border: none; font-weight: 600;
    }}
    QPushButton[variant="success"]:hover {{ background: #2A8E51; }}
    QPushButton[variant="success"]:pressed {{ background: #247D47; }}
    QPushButton[variant="success"]:disabled {{ background: #BFE0CB; color: #FFFFFF; }}

    QPushButton[variant="ghost"] {{
        background: transparent; border: none; color: {t.ACCENT};
        font-weight: 600; padding: 6px 10px;
    }}
    QPushButton[variant="ghost"]:hover   {{ color: {t.ACCENT_HOVER}; background: {t.ACCENT_SOFT}; }}
    QPushButton[variant="ghost"]:disabled {{ color: {t.TEXT_MUTED}; background: transparent; }}

    /* Inputs */
    QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {{
        background: {t.SURFACE};
        border: 1px solid {t.BORDER_STRONG};
        border-radius: {t.RADIUS_SM}px;
        padding: 6px 8px;
        selection-background-color: {t.ACCENT_SOFT};
    }}
    QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {{
        border: 1px solid {t.ACCENT};
    }}

    /* Tabs — flat underline, no chunky 3D pane */
    QTabWidget::pane {{ border: none; border-top: 1px solid {t.BORDER}; }}
    QTabBar::tab {{
        background: transparent; color: {t.TEXT_SUB};
        padding: {t.S2}px {t.S4}px; border: none;
        border-bottom: 2px solid transparent;
    }}
    QTabBar::tab:selected {{ color: {t.ACCENT}; border-bottom: 2px solid {t.ACCENT}; }}
    QTabBar::tab:hover    {{ color: {t.TEXT}; }}

    /* Tables — a white surface with readable gridlines + a divided header,
       so dense data reads as a structured grid instead of a flat wash. */
    QTableWidget, QTableView, QTreeWidget, QTreeView, QListWidget, QListView {{
        background: {t.SURFACE};
        alternate-background-color: {t.SURFACE_ALT};
        border: 1px solid {t.BORDER_STRONG};
        border-radius: {t.RADIUS_SM}px;
        gridline-color: {t.BORDER};
        selection-background-color: {t.ACCENT_SOFT};
        selection-color: {t.TEXT};
    }}
    QHeaderView::section {{
        background: {t.SURFACE_ALT};
        color: {t.TEXT_SUB};
        padding: 4px 8px;
        border: none;
        border-right: 1px solid {t.BORDER};
        border-bottom: 1px solid {t.BORDER_STRONG};
        font-weight: 600;
    }}
    QTableCornerButton::section {{ background: {t.SURFACE_ALT}; border: none; }}

    /* Standalone separators read as a real divider, not a ghost line. */
    QFrame[frameShape="4"], QFrame[frameShape="5"] {{ color: {t.BORDER_STRONG}; }}

    /* Log / console: a quiet surface, not a green terminal */
    QTextEdit#Log, QPlainTextEdit#Log {{
        background: {t.SURFACE_ALT};
        border: 1px solid {t.BORDER};
        border-radius: {t.RADIUS_SM}px;
        color: {t.TEXT_SUB};
        font-family: "Cascadia Mono", "Consolas", monospace;
        font-size: {t.FS_CAPTION}px;
    }}
    """


def apply_theme(app) -> None:
    """Install the global stylesheet on a QApplication."""
    if app is not None:
        app.setStyleSheet(global_qss())


def refresh(widget) -> None:
    """Re-polish a widget after changing a dynamic property at runtime."""
    style = widget.style()
    style.unpolish(widget)
    style.polish(widget)
    widget.update()


def style_button(btn, variant: str | None = None, *, height: int | None = None):
    """Apply the standard button role + height to *btn*.

    This is the single way APEX assigns a button its place in the hierarchy:

        style_button(btn, "primary", height=Tokens.H_ACTION)   # main action
        style_button(btn, "danger",  height=Tokens.H_ACTION)   # Stop / destructive
        style_button(btn, "ghost",   height=Tokens.H_COMPACT)  # Log / 가이드
        style_button(btn,            height=Tokens.H_BUTTON)   # neutral default

    Pass ``variant=None`` to leave the button as the neutral default. Never
    hand-paint hex via setStyleSheet — that's exactly the inconsistency this
    replaces. Returns the button for chaining.
    """
    if variant:
        btn.setProperty("variant", variant)
    else:
        # Clear any stale role so re-styling a button is predictable.
        btn.setProperty("variant", None)
    # Drop any inline stylesheet so the themed QSS actually wins.
    if btn.styleSheet():
        btn.setStyleSheet("")
    if height is not None:
        btn.setMinimumHeight(int(height))
    refresh(btn)
    return btn
