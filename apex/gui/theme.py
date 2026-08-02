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

from pathlib import Path


# ── Design tokens ──────────────────────────────────────────────────────────
class Tokens:
    # Surfaces — the canvas (BG) is deliberately a few shades off the surface
    # so cards/panels/inputs read as raised "cards on a canvas" and their
    # edges stay visible even with thin borders.
    #
    # NOTE: the colour attributes below are MUTABLE — apply_theme() overwrites
    # them from a named preset in PALETTES (AstralImage-style theme system).
    # The values written here are the "apex-light" preset (the historical
    # default look).
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

    # Single calm accent (replaces the rainbow of material primaries).
    # ACCENT is the *fill* colour (primary buttons, selected underline bars);
    # ACCENT_TEXT is the accent used *as text/outline on surfaces* — identical
    # on light themes, brighter on dark themes where a fill-strength blue is
    # too dim to read as text.
    ACCENT      = "#3A66DB"
    ACCENT_HOVER = "#2F56C0"
    ACCENT_PRESS = "#274AA6"
    ACCENT_SOFT = "#EAF0FE"   # accent-tinted background for info chips
    ACCENT_TEXT = "#3A66DB"

    # Semantic, desaturated vs. material 2014
    OK          = "#2E9E5B"
    WARN        = "#C77A12"
    ERROR       = "#D24343"
    # Soft tinted backgrounds for status banners (status="warn"/"error" cards)
    OK_SOFT     = "#E5F4EB"
    WARN_SOFT   = "#FBF1E2"
    ERROR_SOFT  = "#FBE9E9"

    # Muted variant fills (for disabled primary/danger so they don't look active)
    ACCENT_MUTED = "#C2CCEA"
    ERROR_MUTED  = "#E6C5C5"

    # Matplotlib plot surface (applied to rcParams by apply_theme, so figures
    # created after the theme is set match it — dark themes get grey plots)
    PLOT_BG      = "#FFFFFF"
    PLOT_AXES_BG = "#FFFFFF"
    PLOT_FG      = "#1F2933"
    PLOT_GRID    = "#D0D5DC"

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


# ── Theme presets (AstralImage-style) ──────────────────────────────────────
# Each preset overrides the colour tokens above; geometry/typography never
# change per theme. "apex-light" is the historical default; the dark presets
# are adapted from AstralImage's owner-approved palettes (charcoal / aurora /
# midnight) so the two apps share a family look.
#
# Keys listed in every preset — add a colour token above? add it here too,
# or dark themes will keep the stale light value after a switch.
_LIGHT = {
    "BG": "#E9ECF1", "SURFACE": "#FFFFFF", "SURFACE_ALT": "#F1F3F5",
    "BORDER": "#CDD4DE", "BORDER_STRONG": "#AEB8C6",
    "TEXT": "#1F2933", "TEXT_SUB": "#5B6573", "TEXT_MUTED": "#98A2B3",
    "ACCENT": "#3A66DB", "ACCENT_HOVER": "#2F56C0", "ACCENT_PRESS": "#274AA6",
    "ACCENT_SOFT": "#EAF0FE", "ACCENT_TEXT": "#3A66DB",
    "OK": "#2E9E5B", "WARN": "#C77A12", "ERROR": "#D24343",
    "OK_SOFT": "#E5F4EB", "WARN_SOFT": "#FBF1E2", "ERROR_SOFT": "#FBE9E9",
    "ACCENT_MUTED": "#C2CCEA", "ERROR_MUTED": "#E6C5C5",
    # Geometry participates in presets. Squared corners are the app standard
    # (instrument look, owner-approved 2026-07-11); apex-light keeps the
    # original soft radii as the legacy option.
    "RADIUS": 8, "RADIUS_SM": 6,
    # Matplotlib plot surface — light themes keep publication-white plots;
    # dark themes restyle new figures to a matching grey (never white cards).
    "PLOT_BG": "#FFFFFF", "PLOT_AXES_BG": "#FFFFFF",
    "PLOT_FG": "#1F2933", "PLOT_GRID": "#D0D5DC",
}

PALETTES: dict[str, dict[str, str]] = {
    "apex-light": dict(_LIGHT),
    # Lab Gray — instrument/legacy-Qt feel: blue-cast-free neutral greys,
    # hard borders, squared corners, one restrained slate accent. For eyes
    # that read the classic utilitarian look as "scientific".
    "lab": {
        "BG": "#DEDEDE", "SURFACE": "#F2F2F2", "SURFACE_ALT": "#E6E6E6",
        "BORDER": "#A6A6A6", "BORDER_STRONG": "#7F7F7F",
        "TEXT": "#141414", "TEXT_SUB": "#454545", "TEXT_MUTED": "#7E7E7E",
        "ACCENT": "#2F6FB0", "ACCENT_HOVER": "#3A80C6", "ACCENT_PRESS": "#275E96",
        "ACCENT_SOFT": "#D6DEE6", "ACCENT_TEXT": "#1F5187",
        "OK": "#1F7442", "WARN": "#9A6300", "ERROR": "#A83A31",
        "OK_SOFT": "#DCE8E0", "WARN_SOFT": "#EDE4D1", "ERROR_SOFT": "#EAD9D7",
        "ACCENT_MUTED": "#AFC2D5", "ERROR_MUTED": "#D5BAB6",
        "RADIUS": 2, "RADIUS_SM": 1,
        "PLOT_BG": "#FFFFFF", "PLOT_AXES_BG": "#FFFFFF",
        "PLOT_FG": "#141414", "PLOT_GRID": "#C8C8C8",
    },
    # Charcoal — neutral grey dark (AstralImage default; PixInsight-style).
    "charcoal": {
        "BG": "#212121", "SURFACE": "#2B2B2B", "SURFACE_ALT": "#333333",
        "BORDER": "#3D3D3D", "BORDER_STRONG": "#4E4E4E",
        "TEXT": "#EDEDED", "TEXT_SUB": "#B4B4B4", "TEXT_MUTED": "#8A8A8A",
        "ACCENT": "#2F78AD", "ACCENT_HOVER": "#3A8AC4", "ACCENT_PRESS": "#27638F",
        "ACCENT_SOFT": "#243A4A", "ACCENT_TEXT": "#67B7FF",
        "OK": "#2E9E5B", "WARN": "#C08A28", "ERROR": "#C75048",
        "OK_SOFT": "#22352A", "WARN_SOFT": "#38301D", "ERROR_SOFT": "#3A2626",
        "ACCENT_MUTED": "#2B4356", "ERROR_MUTED": "#5A3535",
        "RADIUS": 2, "RADIUS_SM": 1,
        "PLOT_BG": "#2B2B2B", "PLOT_AXES_BG": "#333333",
        "PLOT_FG": "#D6D6D6", "PLOT_GRID": "#4A4A4A",
    },
    # Aurora — deep cool slate + teal-cyan accent (AstralImage signature).
    "aurora": {
        "BG": "#151A21", "SURFACE": "#1D242D", "SURFACE_ALT": "#232B35",
        "BORDER": "#2D3742", "BORDER_STRONG": "#3B4753",
        "TEXT": "#E6EBF1", "TEXT_SUB": "#AEB9C6", "TEXT_MUTED": "#77828F",
        "ACCENT": "#2A7480", "ACCENT_HOVER": "#338B99", "ACCENT_PRESS": "#22606A",
        "ACCENT_SOFT": "#1C3238", "ACCENT_TEXT": "#4FC8CF",
        "OK": "#3FA56C", "WARN": "#C69B3D", "ERROR": "#CC5F55",
        "OK_SOFT": "#1C3226", "WARN_SOFT": "#332B18", "ERROR_SOFT": "#38221F",
        "ACCENT_MUTED": "#24444B", "ERROR_MUTED": "#54322E",
        "RADIUS": 2, "RADIUS_SM": 1,
        "PLOT_BG": "#1D242D", "PLOT_AXES_BG": "#232B35",
        "PLOT_FG": "#CCD4DE", "PLOT_GRID": "#3B4753",
    },
    # Midnight Navy — deep rich blue.
    "midnight": {
        "BG": "#0C1730", "SURFACE": "#12233F", "SURFACE_ALT": "#182C4D",
        "BORDER": "#1F3355", "BORDER_STRONG": "#2C4470",
        "TEXT": "#CDD8EE", "TEXT_SUB": "#9DABC9", "TEXT_MUTED": "#66748F",
        "ACCENT": "#23558F", "ACCENT_HOVER": "#2C67AC", "ACCENT_PRESS": "#1C4574",
        "ACCENT_SOFT": "#16294D", "ACCENT_TEXT": "#5A93E0",
        "OK": "#4CB87C", "WARN": "#CCAE52", "ERROR": "#CC5F55",
        "OK_SOFT": "#143832", "WARN_SOFT": "#31301F", "ERROR_SOFT": "#3A2430",
        "ACCENT_MUTED": "#1B3455", "ERROR_MUTED": "#4E3048",
        "RADIUS": 2, "RADIUS_SM": 1,
        "PLOT_BG": "#12233F", "PLOT_AXES_BG": "#182C4D",
        "PLOT_FG": "#CDD8EE", "PLOT_GRID": "#2C4470",
    },
}

# (key, menu label) — drives the theme menu; add a preset to PALETTES and
# list it here to expose it.
THEME_PRESETS: tuple[tuple[str, str], ...] = (
    ("apex-light", "APEX Light"),
    ("lab",        "Lab Gray"),
    ("charcoal",   "Charcoal"),
    ("aurora",     "Aurora"),
    ("midnight",   "Midnight Navy"),
)

DEFAULT_THEME = "lab"  # squared instrument look (owner pick, 2026-07-11)

_THEME_FILE = Path.home() / ".apex" / "theme.txt"


def current_theme() -> str:
    """Name of the persisted (or default) theme preset."""
    try:
        name = _THEME_FILE.read_text(encoding="utf-8").strip()
    except OSError:
        return DEFAULT_THEME
    return name if name in PALETTES else DEFAULT_THEME


def _persist_theme(name: str) -> None:
    try:
        _THEME_FILE.parent.mkdir(parents=True, exist_ok=True)
        _THEME_FILE.write_text(name, encoding="utf-8")
    except OSError:
        pass  # persistence is best-effort; the session still gets the theme


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
        color: {t.ACCENT_TEXT};
        border-radius: {t.RADIUS_SM}px;
        padding: {t.S2}px {t.S3}px;
    }}

    /* Status pills */
    QLabel[status="ok"]    {{ color: {t.OK};    font-weight: 600; }}
    QLabel[status="warn"]  {{ color: {t.WARN};  font-weight: 600; }}
    QLabel[status="error"] {{ color: {t.ERROR}; font-weight: 600; }}
    QLabel[status="idle"]  {{ color: {t.TEXT_MUTED}; }}

    /* Status banners — tinted card versions of the pills, for the step-top
       "do X first" notices that used to hand-paint pink/blue rectangles. */
    QLabel[banner="ok"], QLabel[banner="warn"], QLabel[banner="error"] {{
        border-radius: {t.RADIUS_SM}px;
        padding: {t.S2}px {t.S3}px;
        font-weight: 600;
    }}
    QLabel[banner="ok"]    {{ background: {t.OK_SOFT};    color: {t.OK}; }}
    QLabel[banner="warn"]  {{ background: {t.WARN_SOFT};  color: {t.WARN}; }}
    QLabel[banner="error"] {{ background: {t.ERROR_SOFT}; color: {t.ERROR}; }}

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
    /* Checkable neutral buttons (Flip X/Y, toggles) need a visible on-state */
    QPushButton:checked {{
        background: {t.ACCENT_SOFT}; border-color: {t.ACCENT}; color: {t.ACCENT_TEXT};
    }}
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
        background: transparent; border: none; color: {t.ACCENT_TEXT};
        font-weight: 600; padding: 6px 10px;
    }}
    QPushButton[variant="ghost"]:hover   {{ color: {t.ACCENT_TEXT}; background: {t.ACCENT_SOFT}; }}
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
        border: 1px solid {t.ACCENT_TEXT};
    }}

    /* Menus — follow the theme instead of the OS default (matters on dark) */
    QMenuBar {{ background: {t.BG}; }}
    QMenuBar::item {{ background: transparent; padding: 4px 10px; }}
    QMenuBar::item:selected {{ background: {t.ACCENT_SOFT}; border-radius: 4px; }}
    QMenu {{ background: {t.SURFACE}; border: 1px solid {t.BORDER_STRONG}; }}
    QMenu::item {{ padding: 5px 24px 5px 12px; }}
    QMenu::item:selected {{ background: {t.ACCENT_SOFT}; }}
    QMenu::separator {{ height: 1px; background: {t.BORDER}; margin: 4px 8px; }}

    /* Tabs — flat underline, no chunky 3D pane */
    QTabWidget::pane {{ border: none; border-top: 1px solid {t.BORDER}; }}
    QTabBar::tab {{
        background: transparent; color: {t.TEXT_SUB};
        padding: {t.S2}px {t.S4}px; border: none;
        border-bottom: 2px solid transparent;
    }}
    QTabBar::tab:selected {{ color: {t.ACCENT_TEXT}; border-bottom: 2px solid {t.ACCENT_TEXT}; }}
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


def _apply_mpl_theme() -> None:
    """Point matplotlib rcParams at the active plot tokens.

    Dark presets get grey plot surfaces instead of white cards. Only figures
    created *after* this call pick it up — step/tool windows build their
    canvases on open, so a theme switch applies to the next window opened.
    Saved PNGs match the screen (savefig.* follows the same tokens).
    """
    try:
        import matplotlib as mpl
    except ImportError:  # headless installs without plotting
        return
    t = Tokens
    mpl.rcParams.update({
        # Hangul in plot text (axis labels, empty-state notices) rendered as
        # tofu boxes: the default DejaVu Sans has no Hangul glyphs. Matplotlib
        # >= 3.6 falls back per glyph through this list, so DejaVu stays first
        # — Latin text and every existing paper figure keep their exact face —
        # and Korean-capable fonts only fill the glyphs DejaVu lacks.
        "font.sans-serif": ["DejaVu Sans", "Malgun Gothic", "Noto Sans KR",
                            "Segoe UI", "Arial"],
        "figure.facecolor": t.PLOT_BG,
        "figure.edgecolor": t.PLOT_BG,
        "savefig.facecolor": t.PLOT_BG,
        "savefig.edgecolor": t.PLOT_BG,
        "axes.facecolor": t.PLOT_AXES_BG,
        "axes.edgecolor": t.PLOT_FG,
        "axes.labelcolor": t.PLOT_FG,
        "axes.titlecolor": t.PLOT_FG,
        "xtick.color": t.PLOT_FG,
        "ytick.color": t.PLOT_FG,
        "text.color": t.PLOT_FG,
        "grid.color": t.PLOT_GRID,
        "legend.facecolor": t.PLOT_AXES_BG,
        "legend.edgecolor": t.PLOT_GRID,
    })


def apply_theme(app, name: str | None = None) -> str:
    """Load a theme preset into Tokens and (re)install the global stylesheet.

    ``name=None`` uses the persisted preset (``~/.apex/theme.txt``) so all
    three entry points pick up the user's choice. Returns the applied name.
    """
    theme = name if name in PALETTES else current_theme()
    for key, value in PALETTES[theme].items():
        setattr(Tokens, key, value)
    _apply_mpl_theme()
    if app is not None:
        app.setStyleSheet(global_qss())
    return theme


def set_theme(app, name: str) -> str:
    """Apply *and persist* a theme preset (the theme-menu entry point).

    Widgets styled purely by the global QSS restyle instantly; code that
    hand-painted colours must re-run its styling (e.g. the main window's
    step buttons) — callers refresh those explicitly.
    """
    theme = apply_theme(app, name)
    _persist_theme(theme)
    return theme


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
