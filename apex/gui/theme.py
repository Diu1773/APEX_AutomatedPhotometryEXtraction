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
    TEXT_MUTED  = "#8692A6"   # captions / disabled

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
    OK          = "#247A46"
    WARN        = "#985D0E"
    ERROR       = "#C73030"
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
    "TEXT": "#1F2933", "TEXT_SUB": "#5B6573", "TEXT_MUTED": "#8692A6",
    "ACCENT": "#3A66DB", "ACCENT_HOVER": "#2F56C0", "ACCENT_PRESS": "#274AA6",
    "ACCENT_SOFT": "#EAF0FE", "ACCENT_TEXT": "#3A66DB",
    "OK": "#247A46", "WARN": "#985D0E", "ERROR": "#C73030",
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
        "OK": "#1F7442", "WARN": "#8B5900", "ERROR": "#A83A31",
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
        "OK": "#33B267", "WARN": "#C8902A", "ERROR": "#D57D77",
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
        "OK": "#40A96E", "WARN": "#C69B3D", "ERROR": "#D3756D",
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
        "OK": "#4CB87C", "WARN": "#CCAE52", "ERROR": "#D67C75",
        "OK_SOFT": "#143832", "WARN_SOFT": "#31301F", "ERROR_SOFT": "#3A2430",
        "ACCENT_MUTED": "#1B3455", "ERROR_MUTED": "#4E3048",
        "RADIUS": 2, "RADIUS_SM": 1,
        "PLOT_BG": "#12233F", "PLOT_AXES_BG": "#182C4D",
        "PLOT_FG": "#CDD8EE", "PLOT_GRID": "#2C4470",
    },

    # ---- Industry-standard dark presets -------------------------------------
    # Surfaces, text and accents are the published values of each design; only
    # the tokens APEX needs but the original doesn't define (SOFT/MUTED fills)
    # are derived, by blending the accent or status colour into the surface.
    #
    # Three ERROR values are lifted a few percent in lightness off their
    # published hex (VS Code #F14C4C, One Dark #E06C75, Gruvbox #FB4934).
    # Those are editor *syntax* colours, tuned against each editor's own
    # background; on APEX's panel surface they land at 3.8-4.4 contrast, and
    # APEX also uses ERROR as label text, not just as a token tint. Hue and
    # saturation are untouched, so the designs still read as themselves.

    # Chrome dark mode (Google). Neutral greys, no blue cast in the surfaces —
    # the blue lives only in the accent. RADIUS 4 matches Chrome's rounding.
    "chrome-dark": {
        "BG": "#202124", "SURFACE": "#292A2D", "SURFACE_ALT": "#35363A",
        "BORDER": "#3C4043", "BORDER_STRONG": "#5F6368",
        "TEXT": "#E8EAED", "TEXT_SUB": "#BDC1C6", "TEXT_MUTED": "#9AA0A6",
        "ACCENT": "#3B6FCC", "ACCENT_HOVER": "#4C82E0", "ACCENT_PRESS": "#2F5CAB",
        "ACCENT_SOFT": "#26313F", "ACCENT_TEXT": "#8AB4F8",
        "OK": "#81C995", "WARN": "#FDD663", "ERROR": "#F28B82",
        "OK_SOFT": "#25342A", "WARN_SOFT": "#3A3222", "ERROR_SOFT": "#3E2A28",
        "ACCENT_MUTED": "#2E3B4D", "ERROR_MUTED": "#4A3330",
        "RADIUS": 4, "RADIUS_SM": 2,
        "PLOT_BG": "#292A2D", "PLOT_AXES_BG": "#35363A",
        "PLOT_FG": "#E8EAED", "PLOT_GRID": "#4A4D51",
    },

    # VS Code Dark+ (Microsoft). The most familiar developer dark: near-black
    # neutral surfaces with the signature #007ACC chrome blue.
    "vscode-dark": {
        "BG": "#1E1E1E", "SURFACE": "#252526", "SURFACE_ALT": "#2D2D30",
        "BORDER": "#3E3E42", "BORDER_STRONG": "#565659",
        "TEXT": "#D4D4D4", "TEXT_SUB": "#BBBBBB", "TEXT_MUTED": "#858585",
        "ACCENT": "#0E639C", "ACCENT_HOVER": "#1177BB", "ACCENT_PRESS": "#0A4C78",
        "ACCENT_SOFT": "#1B2C38", "ACCENT_TEXT": "#3794FF",
        "OK": "#4EC9B0", "WARN": "#CCA700", "ERROR": "#F25A5A",
        "OK_SOFT": "#1D3330", "WARN_SOFT": "#332E1A", "ERROR_SOFT": "#2C1A1A",
        "ACCENT_MUTED": "#28394A", "ERROR_MUTED": "#4A2E2E",
        "RADIUS": 2, "RADIUS_SM": 1,
        "PLOT_BG": "#252526", "PLOT_AXES_BG": "#2D2D30",
        "PLOT_FG": "#D4D4D4", "PLOT_GRID": "#4A4A4E",
    },

    # One Dark (Atom). Slightly cool blue-grey surfaces, high-legibility text.
    "one-dark": {
        "BG": "#21252B", "SURFACE": "#282C34", "SURFACE_ALT": "#2F343D",
        "BORDER": "#3E4451", "BORDER_STRONG": "#4F5666",
        "TEXT": "#ABB2BF", "TEXT_SUB": "#9199A6", "TEXT_MUTED": "#727986",
        "ACCENT": "#3D6FB5", "ACCENT_HOVER": "#4B84D4", "ACCENT_PRESS": "#325C96",
        "ACCENT_SOFT": "#25303F", "ACCENT_TEXT": "#61AFEF",
        "OK": "#98C379", "WARN": "#E5C07B", "ERROR": "#E2747D",
        "OK_SOFT": "#28331F", "WARN_SOFT": "#332E1F", "ERROR_SOFT": "#33242A",
        "ACCENT_MUTED": "#2C3A4B", "ERROR_MUTED": "#432E33",
        "RADIUS": 3, "RADIUS_SM": 2,
        "PLOT_BG": "#282C34", "PLOT_AXES_BG": "#2F343D",
        "PLOT_FG": "#ABB2BF", "PLOT_GRID": "#3E4451",
    },

    # Dracula (official spec). Violet-leaning surfaces; the accent is purple,
    # so filled buttons take dark ink automatically via ink_on().
    "dracula": {
        "BG": "#21222C", "SURFACE": "#282A36", "SURFACE_ALT": "#343746",
        "BORDER": "#44475A", "BORDER_STRONG": "#5A5D75",
        "TEXT": "#F8F8F2", "TEXT_SUB": "#D6D6D0", "TEXT_MUTED": "#6272A4",
        "ACCENT": "#6D50B8", "ACCENT_HOVER": "#8062D0", "ACCENT_PRESS": "#5A4199",
        "ACCENT_SOFT": "#2E2A44", "ACCENT_TEXT": "#BD93F9",
        "OK": "#50FA7B", "WARN": "#F1FA8C", "ERROR": "#FF5555",
        "OK_SOFT": "#23372B", "WARN_SOFT": "#36371F", "ERROR_SOFT": "#37232E",
        "ACCENT_MUTED": "#37324F", "ERROR_MUTED": "#4A2F3A",
        "RADIUS": 4, "RADIUS_SM": 2,
        "PLOT_BG": "#282A36", "PLOT_AXES_BG": "#343746",
        "PLOT_FG": "#F8F8F2", "PLOT_GRID": "#44475A",
    },

    # Gruvbox Dark (retro groove). Warm brown-grey surfaces and an amber
    # accent — the one preset here with no blue anywhere.
    "gruvbox-dark": {
        "BG": "#282828", "SURFACE": "#32302F", "SURFACE_ALT": "#3C3836",
        "BORDER": "#504945", "BORDER_STRONG": "#665C54",
        "TEXT": "#EBDBB2", "TEXT_SUB": "#D5C4A1", "TEXT_MUTED": "#A89984",
        "ACCENT": "#AF7A18", "ACCENT_HOVER": "#C88F22", "ACCENT_PRESS": "#8F6412",
        "ACCENT_SOFT": "#3B3223", "ACCENT_TEXT": "#FABD2F",
        "OK": "#B8BB26", "WARN": "#FE8019", "ERROR": "#FC6D5C",
        "OK_SOFT": "#333520", "WARN_SOFT": "#3B3123", "ERROR_SOFT": "#2D1E1B",
        "ACCENT_MUTED": "#463B29", "ERROR_MUTED": "#4C332D",
        "RADIUS": 2, "RADIUS_SM": 1,
        "PLOT_BG": "#32302F", "PLOT_AXES_BG": "#3C3836",
        "PLOT_FG": "#EBDBB2", "PLOT_GRID": "#504945",
    },
}

# (key, menu label) — drives the theme menu; add a preset to PALETTES and
# list it here to expose it.
THEME_PRESETS: tuple[tuple[str, str], ...] = (
    ("apex-light",   "APEX Light"),
    ("lab",          "Lab Gray"),
    ("charcoal",     "Charcoal"),
    ("aurora",       "Aurora"),
    ("midnight",     "Midnight Navy"),
    ("chrome-dark",  "Chrome Dark"),
    ("vscode-dark",  "VS Code Dark+"),
    ("one-dark",     "One Dark"),
    ("dracula",      "Dracula"),
    ("gruvbox-dark", "Gruvbox Dark"),
)

# Presets that reproduce a published third-party design rather than an APEX
# one. The theme menu draws a separator before the first of these so the two
# groups don't read as one long undifferentiated list.
STANDARD_THEMES: frozenset[str] = frozenset({
    "chrome-dark", "vscode-dark", "one-dark", "dracula", "gruvbox-dark",
})

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

    /* Tick boxes, radio dots and slider handles used to be left to the native
       Windows style, so they painted the *OS* accent (blue) in every preset —
       fine while every APEX dark theme was navy, wrong the moment a palette
       isn't blue. Drive them from the palette instead. */
    QCheckBox::indicator, QRadioButton::indicator {{
        width: 14px; height: 14px;
        border: 1px solid {t.BORDER_STRONG};
        background: {t.SURFACE};
    }}
    QCheckBox::indicator {{ border-radius: {t.RADIUS_SM}px; }}
    QRadioButton::indicator {{ border-radius: 8px; }}
    QCheckBox::indicator:hover, QRadioButton::indicator:hover {{
        border-color: {t.ACCENT};
    }}
    /* A filled box is the checked state — APEX ships no tick glyph resource,
       and pointing at a missing one would blank the indicator entirely. */
    QCheckBox::indicator:checked {{
        background: {t.ACCENT}; border-color: {t.ACCENT};
    }}
    QRadioButton::indicator:checked {{
        background: {t.ACCENT}; border-color: {t.ACCENT};
    }}
    QCheckBox::indicator:disabled, QRadioButton::indicator:disabled {{
        background: {t.SURFACE_ALT}; border-color: {t.BORDER};
    }}
    QCheckBox::indicator:checked:disabled, QRadioButton::indicator:checked:disabled {{
        background: {t.ACCENT_MUTED}; border-color: {t.ACCENT_MUTED};
    }}

    QSlider::groove:horizontal {{
        height: 4px; background: {t.BORDER}; border-radius: 2px;
    }}
    QSlider::sub-page:horizontal {{ background: {t.ACCENT}; border-radius: 2px; }}
    QSlider::handle:horizontal {{
        background: {t.ACCENT}; border: none;
        width: 12px; margin: -5px 0; border-radius: {t.RADIUS_SM + 1}px;
    }}
    QSlider::handle:horizontal:hover {{ background: {t.ACCENT_HOVER}; }}
    QSlider::handle:horizontal:disabled {{ background: {t.BORDER_STRONG}; }}

    /* Same story as the indicators: unthemed, QProgressBar drew the native
       near-white trough, so every dark preset had a glaring white bar across
       its run row. Individual windows had started pasting their own fix. */
    QProgressBar {{
        background: {t.SURFACE_ALT};
        border: 1px solid {t.BORDER};
        border-radius: {t.RADIUS_SM}px;
        text-align: center;
        color: {t.TEXT};
    }}
    QProgressBar::chunk {{
        background: {t.ACCENT};
        border-radius: {t.RADIUS_SM}px;
    }}

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
        background: {t.ACCENT}; color: {ink_on(t.ACCENT)};
        border: none; font-weight: 600;
    }}
    QPushButton[variant="primary"]:hover   {{ background: {t.ACCENT_HOVER}; }}
    QPushButton[variant="primary"]:pressed {{ background: {t.ACCENT_PRESS}; }}
    QPushButton[variant="primary"]:disabled {{
        background: {t.ACCENT_MUTED}; color: {ink_on(t.ACCENT_MUTED)};
    }}

    QPushButton[variant="danger"] {{
        background: {t.ERROR}; color: {ink_on(t.ERROR)};
        border: none; font-weight: 600;
    }}
    QPushButton[variant="danger"]:hover   {{ background: {shade(t.ERROR, 1.12)}; }}
    QPushButton[variant="danger"]:pressed {{ background: {shade(t.ERROR, 0.86)}; }}
    QPushButton[variant="danger"]:disabled {{
        background: {t.ERROR_MUTED}; color: {ink_on(t.ERROR_MUTED)};
    }}

    QPushButton[variant="success"] {{
        background: {t.OK}; color: {ink_on(t.OK)};
        border: none; font-weight: 600;
    }}
    QPushButton[variant="success"]:hover {{ background: {shade(t.OK, 1.12)}; }}
    QPushButton[variant="success"]:pressed {{ background: {shade(t.OK, 0.86)}; }}
    QPushButton[variant="success"]:disabled {{
        background: {mix(t.OK, t.SURFACE, 0.62)}; color: {ink_on(mix(t.OK, t.SURFACE, 0.62))};
    }}

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


def _srgb_luminance(rgb) -> float:
    """WCAG relative luminance of an (r, g, b) 0-255 tuple."""
    def _lin(c: float) -> float:
        c /= 255.0
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (_lin(v) for v in rgb)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _contrast(a, b) -> float:
    la, lb = _srgb_luminance(a), _srgb_luminance(b)
    hi, lo = max(la, lb), min(la, lb)
    return (hi + 0.05) / (lo + 0.05)


def _parse_hex(value: str):
    """``#rgb`` / ``#rrggbb`` -> (r, g, b), or None if it isn't a colour."""
    value = str(value or "").strip().lstrip("#")
    if len(value) == 3:
        value = "".join(c * 2 for c in value)
    if len(value) != 6:
        return None
    try:
        return tuple(int(value[i:i + 2], 16) for i in (0, 2, 4))
    except ValueError:
        return None


def _to_hex(rgb) -> str:
    return "#{:02X}{:02X}{:02X}".format(
        *(max(0, min(255, int(round(v)))) for v in rgb))


def shade(color: str, factor: float) -> str:
    """Lighten (``factor`` > 1) or darken (< 1) *color* for interaction states.

    A design system derives a control's hover/pressed colours from its base
    rather than pinning a second hex, so every palette gets a coherent set for
    free — the old QSS hard-coded ``#BE3A3A``/``#2A8E51`` for those states,
    which stayed put no matter which preset was loaded. Lightening blends
    toward white so a near-black base still moves.
    """
    rgb = _parse_hex(color)
    if rgb is None:
        return color
    r, g, b = rgb
    if factor >= 1.0:
        t = min(1.0, factor - 1.0)
        return _to_hex((r + (255 - r) * t, g + (255 - g) * t, b + (255 - b) * t))
    return _to_hex((r * factor, g * factor, b * factor))


def mix(color: str, other: str, t: float) -> str:
    """Blend *color* toward *other* by fraction *t* (0..1)."""
    a, b = _parse_hex(color), _parse_hex(other)
    if a is None or b is None:
        return color
    return _to_hex(tuple(a[i] + (b[i] - a[i]) * t for i in range(3)))


def ink_on(background: str) -> str:
    """Readable text colour for a *filled* control on ``background``.

    Chrome's own filled buttons put dark text on their light blue; pinning
    white — as the QSS used to — fails contrast on any light accent, so a
    palette with a pale accent would ship unreadable button labels.
    """
    rgb = _parse_hex(background)
    if rgb is None:
        return "#FFFFFF"
    return "#FFFFFF" if _contrast((255, 255, 255), rgb) >= 4.5 else "#111418"


def readable_on(color: str, background: str | None = None,
                min_ratio: float = 4.5) -> str:
    """Nudge ``color`` until it is legible on ``background``, keeping its hue.

    Data-bound text (a legend tinted with its own marker colour) has to stay
    recognisably that colour AND stay readable — pure gold reads fine on a dark
    canvas and vanishes on a white one, which is why such labels used to carry
    a second, hand-darkened hex that then drifted from the marker. This blends
    the colour toward black or white (whichever raises contrast) in small steps
    until it clears ``min_ratio``, so one source of truth serves every theme.
    """
    rgb = _parse_hex(color)
    bg = _parse_hex(background if background is not None else Tokens.BG)
    if rgb is None or bg is None:
        return color
    # Blend toward whichever end is further from the background.
    target = (0, 0, 0) if _srgb_luminance(bg) > 0.5 else (255, 255, 255)
    best = rgb
    for step in range(21):                      # up to 100% in 5% steps
        if _contrast(best, bg) >= min_ratio:
            break
        t = step / 20.0
        best = tuple(round(c + (target[i] - c) * t) for i, c in enumerate(rgb))
    return "#{:02X}{:02X}{:02X}".format(*best)


def mono_note_style() -> str:
    """Stylesheet for a monospace note/summary QLabel (paths, fit results).

    Built from the live Tokens so every preset keeps its own surface and
    border; no QSS property exists for a mono QLabel, hence the one sanctioned
    f-string. Windows re-apply it on reopen after a theme switch.
    """
    t = Tokens
    return (
        f"QLabel {{ background: {t.SURFACE_ALT}; color: {t.TEXT}; "
        f"padding: {t.GAP}px; border: 1px solid {t.BORDER}; "
        f"border-radius: {t.RADIUS_SM}px; "
        f"font-family: 'Cascadia Mono', 'Consolas', monospace; "
        f"font-size: {t.FS_CAPTION}px; }}"
    )


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
