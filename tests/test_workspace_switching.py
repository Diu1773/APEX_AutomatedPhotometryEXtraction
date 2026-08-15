"""Switching workspace from the menu, without carrying the old one along.

Changing target used to mean editing a config by hand or relaunching with
``--params``, so several things that were shared across the whole application
were never noticed to be shared. A File > Open Workspace makes them visible at
once: progress marks, open step/tool windows, and the window title all belonged
to the process rather than to the target.

These pin the three that would silently mislead — an untouched cluster showing
the previous one's completed steps, a step window still writing into the
previous target's directories, and a title bar that cannot tell two workspaces
apart.
"""

from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

from apex.gui.main_window import MainWindowWorkflow


def _workspace(tmp_path, name: str, target: str, completed=None) -> Path:
    root = tmp_path / name
    (root / "sci").mkdir(parents=True)
    (root / "result").mkdir(parents=True)
    config = root / "apex_config.json"
    config.write_text(json.dumps({
        "io": {"data_dir": str(root / "sci"),
               "result_dir": str(root / "result"), "cache_dir": "cache"},
        "target": {"name": target},
    }), encoding="utf-8")
    if completed is not None:
        (root / "result" / "project_state.json").write_text(
            json.dumps({"completed_steps": completed,
                        "current_step": (max(completed) + 1) if completed else 0}),
            encoding="utf-8")
    return config


def _fake(config: Path, mode: str = "cmd") -> SimpleNamespace:
    data = json.loads(config.read_text(encoding="utf-8"))
    return SimpleNamespace(
        mode=mode,
        params=SimpleNamespace(
            param_file=config,
            P=SimpleNamespace(result_dir=data["io"]["result_dir"]),
        ),
    )


def test_progress_belongs_to_the_workspace(tmp_path):
    """The bug this closes: an untouched cluster reporting 9/12 done."""
    busy = _workspace(tmp_path, "M13_ws", "M13", completed=[0, 1, 2, 3, 4])
    fresh = _workspace(tmp_path, "M67_ws", "M67")

    busy_state = MainWindowWorkflow._make_project_state(_fake(busy))
    fresh_state = MainWindowWorkflow._make_project_state(_fake(fresh))

    assert sorted(busy_state.state.get("completed_steps", [])) == [0, 1, 2, 3, 4]
    assert fresh_state.state.get("completed_steps", []) == []


def test_state_lives_next_to_the_results_it_describes(tmp_path):
    config = _workspace(tmp_path, "M13_ws", "M13")
    state = MainWindowWorkflow._make_project_state(_fake(config))
    assert state.state_file.parent == config.parent / "result"


def test_a_workspace_without_a_copy_inherits_the_legacy_shared_file(tmp_path, monkeypatch):
    """Existing users keep their progress the first time they open a workspace."""
    legacy_root = tmp_path / "pkg"
    legacy = legacy_root / ".state" / "cmd"
    legacy.mkdir(parents=True)
    (legacy / "project_state.json").write_text(
        json.dumps({"completed_steps": [0, 1], "current_step": 2}), encoding="utf-8")
    monkeypatch.setattr("apex.gui.main_window.__file__",
                        str(legacy_root / "gui" / "main_window.py"))

    config = _workspace(tmp_path, "M13_ws", "M13")
    state = MainWindowWorkflow._make_project_state(_fake(config), migrate_legacy=True)
    assert sorted(state.state.get("completed_steps", [])) == [0, 1]


def test_migration_never_overwrites_a_workspace_that_has_its_own(tmp_path, monkeypatch):
    legacy_root = tmp_path / "pkg"
    legacy = legacy_root / ".state" / "cmd"
    legacy.mkdir(parents=True)
    (legacy / "project_state.json").write_text(
        json.dumps({"completed_steps": [0, 1, 2, 3, 4, 5, 6]}), encoding="utf-8")
    monkeypatch.setattr("apex.gui.main_window.__file__",
                        str(legacy_root / "gui" / "main_window.py"))

    config = _workspace(tmp_path, "M13_ws", "M13", completed=[0])
    state = MainWindowWorkflow._make_project_state(_fake(config), migrate_legacy=True)
    assert sorted(state.state.get("completed_steps", [])) == [0]


def test_every_child_window_is_closed_and_dropped(tmp_path):
    """A window kept open would go on using the previous target's paths."""
    from PyQt5.QtWidgets import QApplication, QWidget

    app = QApplication.instance() or QApplication([])
    assert app is not None

    class Fake:
        pass

    holder = Fake()
    closed = []

    class Tracked(QWidget):
        def close(self):
            closed.append(self)
            return super().close()

    holder.current_step_window = Tracked()
    holder.qa_window = Tracked()
    holder.a_new_tool_window = Tracked()      # a tool nobody updated a list for
    holder.not_a_window = object()

    MainWindowWorkflow._close_child_windows(holder)

    assert len(closed) == 3
    assert holder.current_step_window is None
    assert holder.qa_window is None
    assert holder.a_new_tool_window is None
    assert holder.not_a_window is not None    # untouched


def test_the_title_names_the_target_and_the_folder(tmp_path):
    config = _workspace(tmp_path, "NGC6811_ws", "NGC 6811")
    title = MainWindowWorkflow._workspace_title(_fake(config), config)
    assert "NGC 6811" in title and "NGC6811_ws" in title


def test_the_title_survives_a_config_it_cannot_read(tmp_path):
    broken = tmp_path / "broken_ws" / "apex_config.json"
    broken.parent.mkdir()
    broken.write_text("{ not json", encoding="utf-8")
    title = MainWindowWorkflow._workspace_title(_fake_broken(), broken)
    assert "broken_ws" in title


def _fake_broken() -> SimpleNamespace:
    return SimpleNamespace(mode="cmd", params=SimpleNamespace(param_file=""))
