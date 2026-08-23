"""The window has to put what it did on the directory's record.

Most of the result directories behind the paper were made in windows, not by
`apex run`. A ledger over sixteen of them on 2026-08-23 found a parameter record
in three — all three left by headless runs, two of them that same morning by
accident. So the window wrote products for months and never wrote down what made
them.

`StepWindowBase._finalize_valid_step` is the one place a step is persisted and
marked complete, so it is where the record is written. Two things make that
awkward, and both are what these tests hold:

* finalising is not running. Closing a window on an already-valid step finalises
  it again, so a naive hook would bury the real runs under identical lines.
* the GUI numbers its windows from zero and every other record numbers steps
  from one. If the conversion is dropped the two records disagree about which
  step made a file, and a journal that lies about that is worse than none.

The method is called directly on a stand-in rather than through a constructed
window: it reads four attributes and no Qt state, so building a QApplication
would test the harness instead of the behaviour.
"""

from __future__ import annotations

import os
from types import SimpleNamespace

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from apex.gui.workflow.step_window_base import StepWindowBase  # noqa: E402
from apex.utils import run_journal as journal  # noqa: E402


def _window(result_dir, step_index=3, **settings):
    """Everything `_record_in_journal` touches, and nothing else."""
    P = SimpleNamespace(result_dir=str(result_dir), **settings)
    return SimpleNamespace(params=SimpleNamespace(P=P, param_file=None),
                           step_index=step_index, step_name="Detection")


def _record(stub):
    StepWindowBase._record_in_journal(stub)


def test_the_window_puts_the_step_on_the_record(tmp_path):
    _record(_window(tmp_path, detect_thresh=5.0))

    steps = journal.latest_steps(tmp_path)
    assert list(steps) == [4], "a window at index 3 is step 4 everywhere else"
    line = steps[4]
    assert line["source"] == "gui" and line["title"] == "Detection"
    assert line["settings"]["detect_thresh"] == 5.0


def test_the_values_are_marked_as_the_workspace_not_as_what_the_step_read():
    """A weaker claim than the headless one, and it has to be legible as such.

    Headless stands a recorder in front of the namespace and knows which
    settings the step asked for. The window does not, so it records the
    workspace as it stood — and a reader must be able to tell which of the two
    kinds of statement they are looking at.
    """
    import inspect

    source = inspect.getsource(StepWindowBase._record_in_journal)
    assert 'settings_scope="workspace"' in source


def test_closing_a_window_again_does_not_repeat_the_line(tmp_path):
    """Finalising is not running; visits must not look like runs."""
    stub = _window(tmp_path, detect_thresh=5.0)
    for _ in range(5):
        _record(stub)

    lines = [l for l in journal.read(tmp_path) if l.get("event") == "step"]
    assert len(lines) == 1


def test_changing_a_setting_does_add_a_line(tmp_path):
    """The point of dropping repeats is to make the changes visible."""
    stub = _window(tmp_path, detect_thresh=5.0)
    _record(stub)
    stub.params.P.detect_thresh = 7.5
    _record(stub)

    lines = [l for l in journal.read(tmp_path) if l.get("event") == "step"]
    assert [l["settings"]["detect_thresh"] for l in lines] == [5.0, 7.5]


def test_a_headless_line_is_not_treated_as_a_repeat(tmp_path):
    """The two sources make different claims, so one cannot silence the other.

    A headless line records what the step read; a window line records the
    workspace. Identical-looking values do not mean the same thing, and dropping
    the window's line because a headless run happened to match would delete the
    record that the window ran at all.
    """
    journal.record_step(tmp_path, "r", index=4, key="detect", status="ok",
                        settings={"detect_thresh": 5.0}, source="headless")
    _record(_window(tmp_path, detect_thresh=5.0))

    lines = [l for l in journal.read(tmp_path) if l.get("event") == "step"]
    assert [l["source"] for l in lines] == ["headless", "gui"]


def test_a_workspace_with_no_result_dir_writes_nothing_and_does_not_raise(tmp_path):
    stub = SimpleNamespace(params=SimpleNamespace(P=SimpleNamespace()),
                           step_index=0, step_name="Scan")
    _record(stub)                       # must not raise
    assert not journal.journal_path(tmp_path).exists()


def test_a_failing_record_never_reaches_the_window(tmp_path):
    """The window's job is the step. Bookkeeping does not get to interrupt it."""

    class Hostile:
        result_dir = str(tmp_path)

        def __getattr__(self, name):
            raise RuntimeError("no")

    stub = SimpleNamespace(params=SimpleNamespace(P=Hostile()),
                           step_index=0, step_name="Scan")
    _record(stub)                       # must not raise


def test_finalising_a_step_is_what_writes_the_line():
    """Pins the hook to the single place a step is persisted and completed.

    There are three callers — Next, the manual complete action, and closing the
    window — and they must not each grow their own copy.
    """
    import inspect

    source = inspect.getsource(StepWindowBase._finalize_valid_step)
    assert "self._record_in_journal()" in source
    assert source.index("self.persist_params()") < source.index("self._record_in_journal()")
