"""Contrast guard for data-bound label colours (apex.gui.theme.readable_on).

A legend label tinted with its own marker colour has to stay recognisably that
colour AND stay legible. Before this existed, such labels carried a second,
hand-darkened hex (step 8's check label was #8A6A00 against a #FFD700 marker)
which then drifted the moment a user recoloured the overlay.
"""

from __future__ import annotations

import pytest

from apex.gui.theme import _contrast, readable_on


def _rgb(value: str):
    value = value.lstrip("#")
    return tuple(int(value[i:i + 2], 16) for i in (0, 2, 4))


LIGHT_BG = "#E9ECF1"     # apex-light canvas
DARK_BG = "#0C1730"      # midnight canvas


def test_gold_is_darkened_on_a_light_canvas():
    """Pure gold is unreadable on white — the historical reason for the
    hand-picked second colour."""
    assert _contrast(_rgb("#FFD700"), _rgb(LIGHT_BG)) < 4.5
    out = readable_on("#FFD700", LIGHT_BG)
    assert out != "#FFD700"
    assert _contrast(_rgb(out), _rgb(LIGHT_BG)) >= 4.5


def test_gold_is_left_alone_on_a_dark_canvas():
    assert readable_on("#FFD700", DARK_BG) == "#FFD700"


def test_dark_red_is_lightened_on_a_dark_canvas():
    assert _contrast(_rgb("#C62828"), _rgb(DARK_BG)) < 4.5
    out = readable_on("#C62828", DARK_BG)
    assert out != "#C62828"
    assert _contrast(_rgb(out), _rgb(DARK_BG)) >= 4.5


def test_already_legible_colour_is_untouched():
    """No gratuitous shifting — a colour that passes keeps its exact hex."""
    assert readable_on("#C62828", LIGHT_BG) == "#C62828"


@pytest.mark.parametrize("color", ["#FFD700", "#C62828", "#D32F2F", "#4CAF50",
                                   "#00BCD4", "#FFD54F", "#FF9800"])
@pytest.mark.parametrize("bg", [LIGHT_BG, DARK_BG, "#FFFFFF", "#000000"])
def test_every_overlay_colour_becomes_legible(color, bg):
    """The full step-8 overlay palette on the extremes of both theme families."""
    assert _contrast(_rgb(readable_on(color, bg)), _rgb(bg)) >= 4.5


def test_bad_input_is_returned_unchanged():
    assert readable_on("", LIGHT_BG) == ""
    assert readable_on("not-a-colour", LIGHT_BG) == "not-a-colour"
    assert readable_on("#FFD700", "nonsense") == "#FFD700"


def test_short_hex_is_accepted():
    out = readable_on("#FD0", LIGHT_BG)          # == #FFDD00
    assert out.startswith("#") and len(out) == 7
    assert _contrast(_rgb(out), _rgb(LIGHT_BG)) >= 4.5
