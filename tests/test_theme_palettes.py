"""Guard: every theme preset must be complete and actually readable.

A palette is 27 hand-written hex values. Two failure modes have real cost:

* a **missing key** — ``apply_theme`` only overwrites the tokens a preset
  lists, so a forgotten key silently keeps the *previous* theme's value. Switch
  from light to a dark preset that forgot ``PLOT_BG`` and you get white plots
  on a dark window.
* an **unreadable pair** — status text sitting on its own tinted background, or
  a filled button's label, can fall below the contrast where it stops being
  legible. That is exactly what happened before the migration: hand-painted
  cards kept light backgrounds on the dark themes.

These tests pin both so a new preset can't ship broken.
"""
from __future__ import annotations

import pytest

pytest.importorskip("PyQt5")

from apex.gui.theme import (  # noqa: E402
    PALETTES, THEME_PRESETS, _contrast, _parse_hex, ink_on, mix, shade,
)

# apex-light is the reference: it defines every token a preset may override.
REFERENCE_KEYS = set(PALETTES["apex-light"])

# Every preset is held to the same bar (owner decision, 2026-08-05). The five
# original palettes were corrected to clear it by hue-preserving lightness
# moves, so a new preset can't be added below the standard either.
STRICT = tuple(sorted(PALETTES))


def _cr(a: str, b: str) -> float:
    return _contrast(_parse_hex(a), _parse_hex(b))


@pytest.mark.parametrize("name", sorted(PALETTES))
def test_palette_defines_every_token(name):
    keys = set(PALETTES[name])
    assert keys == REFERENCE_KEYS, (
        f"{name} differs from the reference token set; missing "
        f"{sorted(REFERENCE_KEYS - keys)}, unexpected {sorted(keys - REFERENCE_KEYS)}. "
        "A missing key keeps the previously applied theme's value."
    )


@pytest.mark.parametrize("name", sorted(PALETTES))
def test_palette_colours_parse(name):
    for key, value in PALETTES[name].items():
        if key.startswith("RADIUS"):
            assert isinstance(value, int), f"{name}.{key} must be an int"
            continue
        assert _parse_hex(value) is not None, f"{name}.{key} = {value!r} is not a colour"


@pytest.mark.parametrize("name", sorted(PALETTES))
def test_body_text_is_readable_on_its_surfaces(name):
    p = PALETTES[name]
    for surface in ("BG", "SURFACE", "SURFACE_ALT"):
        ratio = _cr(p["TEXT"], p[surface])
        assert ratio >= 4.5, f"{name}: TEXT on {surface} is {ratio:.2f}, below 4.5"


@pytest.mark.parametrize("name", sorted(PALETTES))
def test_every_preset_is_reachable_from_the_menu(name):
    assert name in {key for key, _ in THEME_PRESETS}, (
        f"{name} exists in PALETTES but no menu entry exposes it"
    )


@pytest.mark.parametrize("name", STRICT)
def test_status_text_is_readable_on_its_own_tint(name):
    """Banners pair a status colour with its *_SOFT background."""
    p = PALETTES[name]
    for fg, bg in (("OK", "OK_SOFT"), ("WARN", "WARN_SOFT"), ("ERROR", "ERROR_SOFT")):
        ratio = _cr(p[fg], p[bg])
        assert ratio >= 4.5, f"{name}: {fg} on {bg} is {ratio:.2f}, below 4.5"


@pytest.mark.parametrize("name", STRICT)
def test_filled_button_labels_are_readable(name):
    """ink_on picks the label colour, so this also pins ink_on's threshold."""
    p = PALETTES[name]
    for fill in ("ACCENT", "OK", "ERROR", "ACCENT_MUTED", "ERROR_MUTED"):
        ratio = _cr(ink_on(p[fill]), p[fill])
        assert ratio >= 4.5, f"{name}: label on {fill} is {ratio:.2f}, below 4.5"


@pytest.mark.parametrize("name", STRICT)
def test_accent_text_reads_on_the_surface(name):
    p = PALETTES[name]
    ratio = _cr(p["ACCENT_TEXT"], p["SURFACE"])
    assert ratio >= 4.5, f"{name}: ACCENT_TEXT on SURFACE is {ratio:.2f}, below 4.5"


@pytest.mark.parametrize("name", STRICT)
def test_status_colours_read_as_labels_on_a_panel(name):
    """``QLabel[status="ok"|"warn"|"error"]`` puts the status colour straight
    on a panel, with no tint behind it — a separate pairing from the banner."""
    p = PALETTES[name]
    for key in ("OK", "WARN", "ERROR"):
        ratio = _cr(p[key], p["SURFACE"])
        assert ratio >= 4.5, f"{name}: {key} on SURFACE is {ratio:.2f}, below 4.5"


@pytest.mark.parametrize("name", STRICT)
def test_captions_clear_the_large_text_bar(name):
    p = PALETTES[name]
    ratio = _cr(p["TEXT_MUTED"], p["SURFACE"])
    assert ratio >= 3.0, f"{name}: TEXT_MUTED on SURFACE is {ratio:.2f}, below 3.0"


@pytest.mark.parametrize("name", STRICT)
def test_banner_tints_stay_subtle(name):
    """A status banner is a wash behind its text, not a solid block: if the
    tint drifts too far from the panel it reads as a filled bar instead. This
    caught a solver that 'fixed' the light theme's pink card by blacking it."""
    p = PALETTES[name]
    for key in ("OK_SOFT", "WARN_SOFT", "ERROR_SOFT"):
        ratio = _cr(p[key], p["SURFACE"])
        assert 1.0 <= ratio <= 2.0, (
            f"{name}: {key} sits at {ratio:.2f} against SURFACE; a banner tint "
            "should stay within 2.0 of the panel it covers"
        )


def test_shade_lightens_and_darkens():
    assert shade("#808080", 1.5) != "#808080"
    assert _parse_hex(shade("#808080", 1.5)) > _parse_hex("#808080")
    assert _parse_hex(shade("#808080", 0.5)) < _parse_hex("#808080")


def test_shade_moves_a_near_black_base():
    """Multiplying a near-black by >1 would barely move it; lightening has to
    blend toward white so hover states stay visible on dark accents."""
    assert _cr(shade("#0A0A0A", 1.2), "#0A0A0A") > 1.2


def test_shade_stays_in_range():
    for factor in (0.0, 0.5, 1.0, 1.5, 4.0):
        for base in ("#000000", "#FFFFFF", "#3B6FCC"):
            rgb = _parse_hex(shade(base, factor))
            assert rgb is not None and all(0 <= c <= 255 for c in rgb)


def test_mix_interpolates_between_endpoints():
    assert mix("#000000", "#FFFFFF", 0.0) == "#000000"
    assert mix("#000000", "#FFFFFF", 1.0) == "#FFFFFF"
    assert mix("#000000", "#FFFFFF", 0.5) == "#808080"


def test_ink_on_flips_to_dark_for_light_fills():
    assert ink_on("#1A237E") == "#FFFFFF"      # deep blue keeps white
    assert ink_on("#FABD2F") != "#FFFFFF"      # gruvbox amber needs dark ink
    assert ink_on("#8AB4F8") != "#FFFFFF"      # chrome's light blue too


def test_bad_input_is_passed_through_not_crashed():
    assert shade("not-a-colour", 1.2) == "not-a-colour"
    assert mix("nope", "#FFFFFF", 0.5) == "nope"
    assert ink_on("nope") == "#FFFFFF"
