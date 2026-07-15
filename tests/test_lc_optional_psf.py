from __future__ import annotations

import pytest


pytest.importorskip("PyQt5")

from apex.gui.main_window import _migrate_lc_optional_psf_state


class _ProjectState:
    def __init__(self):
        self.state = {
            "completed_steps": list(range(9)),
            "current_step": 8,
            "step_data": {},
        }
        self.save_count = 0

    def save(self):
        self.save_count += 1


def test_lc_optional_psf_migration_shifts_legacy_progress_once():
    state = _ProjectState()

    assert _migrate_lc_optional_psf_state(state) is True
    assert state.state["completed_steps"] == list(range(10))
    assert state.state["current_step"] == 9
    assert state.state["step_data"]["psf_photometry"]["skip_psf"] is True
    assert state.save_count == 1

    assert _migrate_lc_optional_psf_state(state) is False
    assert state.state["completed_steps"] == list(range(10))
    assert state.save_count == 1
