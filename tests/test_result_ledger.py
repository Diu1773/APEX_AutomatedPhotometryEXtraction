"""What `scripts/result_ledger.py` is willing to call provenance.

The ledger's whole job is to say which result directories can explain
themselves. It said so from three signals — a journal, a parameter record, a
run manifest — and two of them could be satisfied by something that is not a
record of a run at all:

  * `apex journal <dir> --note "…"` writes a journal line, and the ledger
    counted journal *entries*. The one command meant for labelling the folders
    nobody can explain was therefore also the way to make them look explained,
    and labelling `result_nocr` and `result_pre20260807` was the next planned
    use of it.
  * a `parameters_used.json` that would not parse was recorded as `-1`, which
    is true, so an unreadable file counted as a parameter record.

These pin both, and the reader that used to stop the whole sweep on one
hand-edited CSV.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO))

from apex.utils import run_journal as journal  # noqa: E402


def _load():
    spec = importlib.util.spec_from_file_location(
        "result_ledger", REPO / "scripts" / "result_ledger.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules["result_ledger"] = module
    spec.loader.exec_module(module)
    return module


ledger = _load()


def _bare_result_dir(tmp_path: Path) -> Path:
    """A directory with products and nothing that says what made them."""
    d = tmp_path / "M67" / "result_nocr"
    (d / "step1_file_selection").mkdir(parents=True)
    return d


# ── a note is not a run ────────────────────────────────────────────────────

def test_a_note_does_not_make_a_directory_explain_itself(tmp_path):
    d = _bare_result_dir(tmp_path)
    assert ledger.provenance_verdict(ledger.describe(d)) == "미상"

    journal.record_note(d, "우주선 제거를 끄고 돌린 판본")

    entry = ledger.describe(d)
    assert entry["journal_runs"] == 0
    assert entry["journal_entries"] == 1
    assert entry["notes"] == ["우주선 제거를 끄고 돌린 판본"]
    assert ledger.provenance_verdict(entry) == "메모만"


def test_a_run_does(tmp_path):
    d = _bare_result_dir(tmp_path)
    journal.record_step(d, "run-A", index=7, key="forcedphot", status="ok",
                        settings={"apcorr_small_scale": 0.8})

    entry = ledger.describe(d)
    assert entry["journal_runs"] == 1
    assert ledger.provenance_verdict(entry) == "저널"


def test_a_note_beside_a_run_does_not_hide_the_run(tmp_path):
    d = _bare_result_dir(tmp_path)
    journal.record_step(d, "run-A", index=7, key="forcedphot", status="ok")
    journal.record_note(d, "논문 Fig 6 이 인용하는 실행")

    entry = ledger.describe(d)
    assert entry["journal_runs"] == 1
    assert entry["journal_entries"] == 2
    assert ledger.provenance_verdict(entry) == "저널"


# ── an unreadable record is not a record ───────────────────────────────────

def test_a_parameter_record_that_will_not_parse_is_not_evidence(tmp_path):
    d = _bare_result_dir(tmp_path)
    (d / "parameters_used.json").write_text("{ truncated", encoding="utf-8")

    entry = ledger.describe(d)
    assert entry["parameters"] is None
    assert entry["parameters_unreadable"] is True
    assert ledger.provenance_verdict(entry) == "미상"


def test_a_parameter_record_that_parses_is(tmp_path):
    d = _bare_result_dir(tmp_path)
    (d / "parameters_used.json").write_text(
        json.dumps({"settings": {"apcorr_small_scale": 0.8}}), encoding="utf-8")

    entry = ledger.describe(d)
    assert entry["parameters"] == 1
    assert ledger.provenance_verdict(entry) == "파라미터만"


def test_a_manifest_that_will_not_parse_is_not_a_step_list(tmp_path):
    """`["<읽을 수 없음>"]` has length 1, which used to read as one step."""
    d = _bare_result_dir(tmp_path)
    (d / "pipeline_run.json").write_text("{ truncated", encoding="utf-8")

    entry = ledger.describe(d)
    assert entry["manifest_readable"] is False
    assert ledger.provenance_verdict(entry) == "미상"


# ── one bad file must not stop the sweep ───────────────────────────────────

def test_a_headers_csv_in_another_encoding_does_not_stop_the_ledger(tmp_path):
    """Only OSError was caught here, while every JSON reader beside it caught
    its parse errors too. A cp949 file raised straight out of describe()."""
    d = _bare_result_dir(tmp_path)
    (d / "step1_file_selection" / "headers.csv").write_bytes(
        "file,filter,대상\nx.fit,V,황소자리\n".encode("cp949"))

    seen, filters = ledger.read_headers(d)
    assert seen is None and filters == []

    entry = ledger.describe(d)          # must not raise
    assert entry["frames_seen"] is None


def test_a_readable_headers_csv_still_reports(tmp_path):
    d = _bare_result_dir(tmp_path)
    (d / "step1_file_selection" / "headers.csv").write_text(
        "file,filter\na.fit,V\nb.fit,B\n", encoding="utf-8")

    seen, filters = ledger.read_headers(d)
    assert seen == 2
    assert filters == ["B", "V"]


# ── the summary the paper reads ────────────────────────────────────────────

def test_the_summary_counts_a_noted_directory_as_still_needing_a_rerun(tmp_path, capsys):
    noted = _bare_result_dir(tmp_path)
    journal.record_note(noted, "판본 설명")
    ran = tmp_path / "M67" / "result"
    (ran / "step1_file_selection").mkdir(parents=True)
    journal.record_step(ran, "run-A", index=7, key="forcedphot", status="ok")

    entries = [ledger.describe(p) for p in ledger.find_result_dirs(tmp_path)]
    assert len(entries) == 2

    verdicts = sorted(ledger.provenance_verdict(e) for e in entries)
    assert verdicts == ["메모만", "저널"]
