"""A result directory has to be able to say what was done to it.

The failure these tests hold shut is specific. `pipeline_run.json` is written
with `write_text`, so on 2026-08-23 a `--steps 10` run over M13 replaced its
seven-step record with a one-step one: the products of steps 1-7 stayed on disk
and the statement that they had run did not. Ten of the twelve directories
behind the paper's numbers held no parameter record at all.

So the properties worth testing are not "a file is written" but:

* a later partial run does not erase an earlier full one,
* a run killed halfway leaves what it finished,
* the record holds parameter *values*, not just which names were touched,
* recording never breaks the run it is recording.
"""

from __future__ import annotations

import json
from types import SimpleNamespace

import pytest

from apex.utils import run_journal as journal


class _P:
    """Stands in for a parameter namespace: attribute access, nothing else."""

    def __init__(self, **values):
        self.__dict__.update(values)


# ── the file ───────────────────────────────────────────────────────────────

def test_lines_accumulate_and_never_replace_each_other(tmp_path):
    journal.append(tmp_path, "step", {"run": "a", "index": 1})
    journal.append(tmp_path, "step", {"run": "a", "index": 2})
    journal.append(tmp_path, "step", {"run": "b", "index": 1})

    lines = list(journal.read(tmp_path))
    assert [(l["run"], l["index"]) for l in lines] == [("a", 1), ("a", 2), ("b", 1)]


def test_a_later_partial_run_does_not_erase_the_earlier_full_one(tmp_path):
    """The M13 case, in miniature.

    A run of steps 1-3, then a run of step 3 alone. Afterwards the directory
    must still be able to say that steps 1 and 2 happened, and to say that
    step 3's product came from the second run.
    """
    for index in (1, 2, 3):
        journal.record_step(tmp_path, "first", index=index, key=f"s{index}",
                            status="ok", settings={"thr": 5.0})
    journal.record_step(tmp_path, "second", index=3, key="s3",
                        status="ok", settings={"thr": 9.0})

    newest = journal.latest_steps(tmp_path)
    assert sorted(newest) == [1, 2, 3]
    assert newest[1]["run"] == "first"
    assert newest[3]["run"] == "second"
    assert newest[3]["settings"]["thr"] == 9.0


def test_a_killed_run_leaves_the_steps_it_finished(tmp_path):
    """No `run_end` line, and the finished steps are still readable.

    This is the case a manifest assembled at the end loses entirely, and it is
    the case that matters: a twelve-hour reprocess that dies at hour eleven.
    """
    journal.record_run_start(tmp_path, "r", mode="cmd", plan=[1, 2, 3])
    journal.record_step(tmp_path, "r", index=1, key="s1", status="ok")
    journal.record_step(tmp_path, "r", index=2, key="s2", status="ok")

    run = journal.history(tmp_path)[0]
    assert [s["index"] for s in run["steps"]] == [1, 2]
    assert run["ended"] is None and run["success"] is None


def test_a_half_written_last_line_is_reported_not_dropped(tmp_path):
    """A truncated tail is the normal shape of a kill, so it must be visible.

    A reader that silently skips unparseable lines under-reports exactly the
    runs worth looking at.
    """
    journal.record_step(tmp_path, "r", index=1, key="s1", status="ok")
    with journal.journal_path(tmp_path).open("a", encoding="utf-8") as handle:
        handle.write('{"event": "step", "run": "r", "ind')

    lines = list(journal.read(tmp_path))
    assert len(lines) == 2
    assert lines[1]["event"] == "unreadable" and lines[1]["line"] == 2


# ── what a line holds ──────────────────────────────────────────────────────

def test_the_record_holds_values_not_only_names():
    """The gap that made the old record unusable.

    `parameters_used.json` lists which step read which setting, and the settings
    as they stood when the run *ended*. For anything produced before the config
    was next edited those are two different claims.
    """
    params = SimpleNamespace(P=_P(aperture_radius=0.8, detect_thresh=5.0, unused=1))
    got = journal.settings_snapshot(params, {"aperture_radius", "detect_thresh"})
    assert got == {"aperture_radius": 0.8, "detect_thresh": 5.0}


def test_an_unreadable_setting_is_recorded_as_such_not_omitted():
    """A step that reached through `_raw` leaves a hole; the hole is the point."""

    class Hostile:
        def __getattr__(self, name):
            raise RuntimeError("gone")

    got = journal.settings_snapshot(SimpleNamespace(P=Hostile()), {"thr"})
    assert got == {"thr": "<unreadable>"}


def test_values_json_cannot_carry_are_spelled_out(tmp_path):
    """NaN would round-trip as `null`, which reads as "never set"."""
    params = SimpleNamespace(P=_P(bad=float("nan"), fine=2.5))
    snap = journal.settings_snapshot(params, ["bad", "fine"])
    journal.append(tmp_path, "step", {"settings": snap})

    written = json.loads(journal.journal_path(tmp_path).read_text(encoding="utf-8"))
    assert written["settings"] == {"bad": "nan", "fine": 2.5}


def test_the_config_is_pinned_by_digest_not_by_path(tmp_path):
    """The path is not evidence — this session edited every config three times."""
    cfg = tmp_path / "apex_config.json"
    cfg.write_text('{"a": 1}', encoding="utf-8")
    first = journal.config_fingerprint(SimpleNamespace(param_file=str(cfg)))
    cfg.write_text('{"a": 2}', encoding="utf-8")
    second = journal.config_fingerprint(SimpleNamespace(param_file=str(cfg)))

    assert first["path"] == second["path"]
    assert first["sha256"] and first["sha256"] != second["sha256"]


def test_a_missing_config_is_named_with_a_null_digest():
    got = journal.config_fingerprint(SimpleNamespace(param_file="nowhere.json"))
    assert got["path"].endswith("nowhere.json") and got["sha256"] is None


def test_a_note_records_what_the_pipeline_cannot_know(tmp_path):
    """`result_nocr` means something; no file has ever said what."""
    journal.record_note(tmp_path, "cosmic-ray rejection off, for the CR ablation")
    note = list(journal.read(tmp_path))[0]
    assert note["event"] == "note" and "ablation" in note["text"]


# ── it must not break a run ────────────────────────────────────────────────

def test_a_journal_that_cannot_be_written_does_not_raise(tmp_path):
    """A twelve-hour reprocess must not die reporting that it can't be described."""
    blocked = tmp_path / "wall"
    blocked.write_text("not a directory", encoding="utf-8")
    assert journal.append(blocked / "inside", "step", {"index": 1}) is False


def test_an_unserialisable_value_does_not_raise(tmp_path):
    class Odd:
        def __repr__(self):
            return "<odd>"

    assert journal.append(tmp_path, "step", {"settings": {"x": Odd()}}) is True
    assert json.loads(journal.journal_path(tmp_path).read_text(encoding="utf-8")
                      )["settings"]["x"] == "<odd>"


def test_reading_a_directory_with_no_journal_is_empty_not_an_error(tmp_path):
    assert list(journal.read(tmp_path)) == []
    assert journal.history(tmp_path) == []
    assert journal.latest_steps(tmp_path) == {}


# ── folded views ───────────────────────────────────────────────────────────

def test_history_folds_a_run_into_one_entry(tmp_path):
    journal.record_run_start(tmp_path, "r", mode="cmd", plan=[1, 2], force=True)
    journal.record_step(tmp_path, "r", index=1, key="s1", status="ok",
                        duration_s=1.234, outputs=[tmp_path / "out.tsv"])
    journal.record_step(tmp_path, "r", index=2, key="s2", status="skipped")
    journal.record_run_end(tmp_path, "r", success=True)

    runs = journal.history(tmp_path)
    assert len(runs) == 1
    run = runs[0]
    assert run["mode"] == "cmd" and run["force"] is True
    assert run["success"] is True and run["ended"]
    assert [s["status"] for s in run["steps"]] == ["ok", "skipped"]
    assert run["steps"][0]["duration_s"] == pytest.approx(1.234)
    assert run["steps"][0]["outputs"] == [str(tmp_path / "out.tsv")]


def test_skipped_and_blocked_steps_are_on_the_record_too(tmp_path):
    """"Step 4 was skipped because it was already complete" explains a directory.

    A record holding only the steps that executed would leave it out, and a
    reader would conclude step 4 never ran at all.
    """
    journal.record_step(tmp_path, "r", index=4, key="s4", status="skipped",
                        message="already complete")
    journal.record_step(tmp_path, "r", index=5, key="s5", status="blocked",
                        message="missing inputs: cat.tsv")

    statuses = {l["index"]: l["status"] for l in journal.read(tmp_path)}
    assert statuses == {4: "skipped", 5: "blocked"}


# ── how a reader sees it ───────────────────────────────────────────────────

def test_one_desktop_sitting_is_one_entry(tmp_path):
    """Not one entry per step.

    The window mints no plan, so an id per finalised step made every step look
    like a separate run that had been interrupted — twelve "no end" lines that
    told a reader nothing. One id for the process groups the sitting, which is
    what actually happened.
    """
    first, second = journal.session_id(), journal.session_id()
    assert first == second

    journal.record_step(tmp_path, first, index=1, key="", status="ok", source="gui")
    journal.record_step(tmp_path, first, index=2, key="", status="ok", source="gui")

    runs = journal.history(tmp_path)
    assert len(runs) == 1
    assert [s["index"] for s in runs[0]["steps"]] == [1, 2]


def test_only_an_announced_run_can_be_called_interrupted(tmp_path):
    """`announced` is what separates "a plan died" from "someone used the app".

    A headless run declares its plan up front, so a missing end means it was cut
    short. A desktop sitting never declares one, and reporting it as interrupted
    states the opposite of what happened.
    """
    journal.record_run_start(tmp_path, "plan", mode="cmd", plan=[1])
    journal.record_step(tmp_path, "plan", index=1, key="scan", status="ok")
    journal.record_step(tmp_path, "sitting", index=2, key="", status="ok", source="gui")

    by_id = {r["run"]: r for r in journal.history(tmp_path)}
    assert by_id["plan"]["announced"] is True
    assert by_id["sitting"]["announced"] is False


def test_the_last_line_of_an_entry_is_tracked(tmp_path):
    """A sitting has no end line, so its span is first line to last."""
    journal.record_step(tmp_path, "s", index=1, key="", status="ok", source="gui")
    journal.record_step(tmp_path, "s", index=2, key="", status="ok", source="gui")

    run = journal.history(tmp_path)[0]
    assert run["last"] >= run["started"] and run["ended"] is None
