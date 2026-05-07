import json

from apex.core.project_state import ProjectState


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
