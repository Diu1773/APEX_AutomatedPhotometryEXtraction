"""Every key the mode presets are read for has to exist in every preset.

The Step 8 parameter dialog offers normal / crowded / faint. Picking one calls
`_apply_mode_to_widgets`, which reads the preset with `p["key"]` — a plain
subscript, so a key the preset does not define raises KeyError and the dialog
dies on a click.

That is exactly what happened: `c0ddf95` added grouper-budget widgets and the
lines that load them from a preset, but never added the keys to the three
preset dicts. Nothing failed at import, no test covered the click, and the
dialog was left unable to switch modes. Adding a widget is a two-place edit and
this pins the second place.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from apex.gui.workflow.cmd.step8_psf_photometry import _PSF_MODE_PRESETS

SOURCE = Path(__file__).absolute().parents[1] / "apex" / "gui" / "workflow" / "cmd" / "step8_psf_photometry.py"


def _keys_the_dialog_reads() -> set[str]:
    text = SOURCE.read_text(encoding="utf-8")
    start = text.index("def _apply_mode_to_widgets(mode_key):")
    end = text.index("_epsf_only_widgets = [", start)
    return set(re.findall(r'p\["(\w+)"\]', text[start:end]))


@pytest.mark.parametrize("mode", sorted(_PSF_MODE_PRESETS))
def test_a_preset_answers_every_key_the_dialog_asks_for(mode):
    missing = sorted(_keys_the_dialog_reads() - set(_PSF_MODE_PRESETS[mode]))
    assert not missing, f"{mode} 프리셋에 없는 키: {missing}"


def test_the_presets_all_carry_the_same_keys():
    """A key in one mode only is a setting that changes by accident."""
    shapes = {mode: frozenset(values) for mode, values in _PSF_MODE_PRESETS.items()}
    assert len(set(shapes.values())) == 1, {
        mode: sorted(keys ^ next(iter(shapes.values())))
        for mode, keys in shapes.items()
    }


def test_no_mode_turns_the_profile_error_term_on():
    """It is instrument-dependent, not a crowding dial — see the audit README."""
    for mode, values in _PSF_MODE_PRESETS.items():
        assert values["psf_profile_error_frac"] == 0.0, mode
