import json

import pytest

from apex.core.project_state import ProjectState, STATE_SCHEMA_VERSION


def test_load_without_steps_preserves_completed_steps(tmp_path):
    state_file = tmp_path / "project_state.json"
    state_file.write_text(
        json.dumps({
            "project_name": "state-load-test",
            "created": "2026-05-02T00:00:00",
            "last_modified": "2026-05-02T00:00:00",
            "current_step": 3,
            "completed_steps": [0, 1, 2, 3],
            "step_data": {"source_detection": {"detection_complete": True}},
        }),
        encoding="utf-8",
    )

    state = ProjectState(tmp_path)

    assert state.state["completed_steps"] == [0, 1, 2, 3]
    assert state.state["current_step"] == 3


def test_assign_steps_clamps_loaded_state_after_preserving_it(tmp_path):
    state_file = tmp_path / "project_state.json"
    state_file.write_text(
        json.dumps({
            "project_name": "state-clamp-test",
            "created": "2026-05-02T00:00:00",
            "last_modified": "2026-05-02T00:00:00",
            "current_step": 8,
            "completed_steps": [0, 1, 2, 8],
            "step_data": {},
        }),
        encoding="utf-8",
    )

    state = ProjectState(tmp_path)
    state.assign_steps(["File Selection", "Image Crop", "Sky Preview & QC"])

    assert state.state["completed_steps"] == [0, 1, 2]
    assert state.state["current_step"] == 2


# ── Detector-calibration off-chain state (schema v2 migration) ────────────────

def test_legacy_state_migrates_to_calibration_skipped(tmp_path):
    # A pre-Step-0 project (no version, no calibration key) must migrate without
    # shifting any progress: completed_steps/current_step stay identical, and
    # calibration is marked "skipped" (a resolved state, no nagging).
    state_file = tmp_path / "project_state.json"
    state_file.write_text(
        json.dumps({
            "project_name": "legacy",
            "created": "2026-05-02T00:00:00",
            "last_modified": "2026-05-02T00:00:00",
            "current_step": 3,
            "completed_steps": [0, 1, 2, 3],
            "step_data": {},
        }),
        encoding="utf-8",
    )
    state = ProjectState(tmp_path)
    assert state.calibration_status() == "skipped"
    assert state.state["state_schema_version"] == STATE_SCHEMA_VERSION
    # progress untouched (Scheme B: no index shift)
    assert state.state["completed_steps"] == [0, 1, 2, 3]
    assert state.state["current_step"] == 3


def test_fresh_state_calibration_not_run(tmp_path):
    state = ProjectState(tmp_path)
    assert state.calibration_status() == "not_run"
    assert not state.is_calibration_done()


def test_mark_calibration_roundtrip_and_persist(tmp_path):
    state = ProjectState(tmp_path)
    state.mark_calibration("done")
    assert state.is_calibration_done()
    # persisted + reloaded
    reloaded = ProjectState(tmp_path)
    assert reloaded.calibration_status() == "done"
    with pytest.raises(ValueError):
        state.mark_calibration("bogus")


def test_migration_idempotent_on_v2(tmp_path):
    state = ProjectState(tmp_path)
    state.mark_calibration("skipped")
    # reload twice; status and version stay put
    again = ProjectState(tmp_path)
    assert again.calibration_status() == "skipped"
    assert again.state["state_schema_version"] == STATE_SCHEMA_VERSION
