from __future__ import annotations

from types import SimpleNamespace

from apex.gui.main_window import MainWindowWorkflow


def test_failed_step_load_keeps_previous_window_and_state(monkeypatch):
    class PreviousWindow:
        closed = False

        def close(self):
            self.closed = True

    class Status:
        def showMessage(self, _message):
            pass

        def clearMessage(self):
            pass

    current_steps = []
    project_state = SimpleNamespace(
        is_step_accessible=lambda _index: True,
        set_current_step=current_steps.append,
    )
    previous = PreviousWindow()
    logs = []
    window = SimpleNamespace(
        project_state=project_state,
        step_names=["one", "two"],
        current_step_window=previous,
        append_log=logs.append,
        statusBar=lambda: Status(),
    )

    def fail_open(_index):
        raise RuntimeError("broken constructor")

    window._open_step_window = fail_open
    monkeypatch.setattr(
        "apex.gui.main_window.QApplication.setOverrideCursor", lambda *_args: None
    )
    monkeypatch.setattr(
        "apex.gui.main_window.QApplication.restoreOverrideCursor", lambda: None
    )
    monkeypatch.setattr(
        "apex.gui.main_window.QApplication.processEvents", lambda *_args: None
    )
    monkeypatch.setattr(
        "apex.gui.main_window.QMessageBox.critical", lambda *_args: None
    )

    MainWindowWorkflow.open_step(window, 1)

    assert not previous.closed
    assert current_steps == []
    assert any("broken constructor" in message for message in logs)
