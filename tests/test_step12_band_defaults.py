"""Which magnitude a CMD opens against, when the user has not said.

Step 12 populates its band combos from whatever calibrated columns the Step 10
table has, then fell back to the first entry. On Johnson data that is B, so a
fresh workspace opened plotting "B vs B−V" — readable, but not what a CMD is:
the magnitude axis conventionally carries the redder half of the colour.

A saved choice still wins. Changing target used to be a restart, so a default
that quietly overrode a stored selection would be the more annoying bug of the
two — these pin both halves.
"""

from __future__ import annotations

import pandas as pd
import pytest

pytest.importorskip("PyQt5")

from PyQt5.QtWidgets import QApplication, QComboBox     # noqa: E402

from apex.gui.workflow.cmd.step12_isochrone_model import (  # noqa: E402
    IsochroneModelWindow,
)


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


class _Bare:
    """Just the two combos and the state the method touches."""

    def __init__(self, color_items, mag_items, pending_color=None, pending_mag=None):
        self.color_combo = QComboBox()
        self.color_combo.addItems(color_items)
        self.mag_combo = QComboBox()
        self.mag_combo.addItems(mag_items)
        self._pending_band_color_text = pending_color
        self._pending_band_mag_text = pending_mag


def _repopulate(host, df):
    IsochroneModelWindow._repopulate_band_combos(host, df)
    return host.color_combo.currentText(), host.mag_combo.currentText()


def _johnson_table():
    return pd.DataFrame({
        "mag_std_B": [15.0, 16.0], "mag_std_err_B": [0.01, 0.01],
        "mag_std_V": [14.5, 15.4], "mag_std_err_V": [0.01, 0.01],
        "mag_std_R": [14.2, 15.0], "mag_std_err_R": [0.01, 0.01],
    })


def _sdss_table():
    return pd.DataFrame({
        "mag_std_g": [15.0, 16.0], "mag_std_err_g": [0.01, 0.01],
        "mag_std_r": [14.5, 15.4], "mag_std_err_r": [0.01, 0.01],
        "mag_std_i": [14.2, 15.0], "mag_std_err_i": [0.01, 0.01],
    })


def test_johnson_opens_on_the_redder_half(app):
    """B-V should plot against V, not B."""
    host = _Bare(["g-r"], ["g", "r", "i"])       # stale defaults, nothing saved
    color, mag = _repopulate(host, _johnson_table())
    assert color.endswith("-" + mag), f"{color} 인데 밝기축이 {mag}"
    assert mag == "V"


def test_sdss_opens_on_the_redder_half(app):
    host = _Bare(["B-V"], ["B", "V", "R"])
    color, mag = _repopulate(host, _sdss_table())
    assert color.endswith("-" + mag), f"{color} 인데 밝기축이 {mag}"


def test_a_saved_choice_is_honoured(app):
    """The workspace remembers; a convention must not overrule it."""
    host = _Bare(["g-r"], ["g", "r", "i"],
                 pending_color="B-V", pending_mag="B")
    color, mag = _repopulate(host, _johnson_table())
    assert (color, mag) == ("B-V", "B")


def test_a_saved_band_that_no_longer_exists_falls_back(app):
    """A gri workspace opened with a stored Johnson band must not go blank."""
    host = _Bare(["g-r"], ["g", "r", "i"], pending_mag="V")
    color, mag = _repopulate(host, _sdss_table())
    assert mag in {"g", "r", "i"}
    assert color.endswith("-" + mag)
