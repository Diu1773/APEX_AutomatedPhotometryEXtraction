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
    # Surfaces
    BG          = "#F6F7F9"   # window background
    SURFACE     = "#FFFFFF"   # cards / panels
    SURFACE_ALT = "#F1F3F5"   # subtle fills (inputs, log)

    # Hairlines (used sparingly — prefer spacing over borders)
    BORDER      = "#E4E7EC"
    BORDER_STRONG = "#D0D5DD"

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

    # Geometry
    RADIUS      = 8
    RADIUS_SM   = 6

    # Spacing scale
    S1, S2, S3, S4, S5 = 4, 8, 12, 16, 24

    # Typography
    FONT_STACK  = '"Segoe UI", "Malgun Gothic", "Noto Sans KR", Arial, sans-serif'
    FS_TITLE    = 20
    FS_SUBTITLE = 13
    FS_BASE     = 13
    FS_CAPTION  = 11


def global_qss(t: type[Tokens] = Tokens) -> str:
    """Return the application-wide stylesheet built from the tokens."""
    return f"""
    QWidget {{
        background: {t.BG};
        color: {t.TEXT};
        font-family: {t.FONT_STACK};
        font-size: {t.FS_BASE}px;
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

    /* Buttons — one accent, thin/no border, soft radius, real hover states */
    QPushButton {{
        background: {t.SURFACE};
        color: {t.TEXT};
        border: 1px solid {t.BORDER_STRONG};
        border-radius: {t.RADIUS_SM}px;
        padding: 7px 14px;
    }}
    QPushButton:hover  {{ background: {t.SURFACE_ALT}; }}
    QPushButton:pressed {{ background: {t.BORDER}; }}
    QPushButton:disabled {{ color: {t.TEXT_MUTED}; border-color: {t.BORDER}; }}

    QPushButton[variant="primary"] {{
        background: {t.ACCENT}; color: #FFFFFF; border: none; font-weight: 600;
    }}
    QPushButton[variant="primary"]:hover   {{ background: {t.ACCENT_HOVER}; }}
    QPushButton[variant="primary"]:pressed {{ background: {t.ACCENT_PRESS}; }}

    QPushButton[variant="danger"] {{
        background: {t.SURFACE}; color: {t.ERROR}; border: 1px solid {t.ERROR};
    }}
    QPushButton[variant="danger"]:hover {{ background: #FBECEC; }}

    QPushButton[variant="ghost"] {{
        background: transparent; border: none; color: {t.ACCENT}; font-weight: 600;
    }}
    QPushButton[variant="ghost"]:hover {{ color: {t.ACCENT_HOVER}; }}

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
