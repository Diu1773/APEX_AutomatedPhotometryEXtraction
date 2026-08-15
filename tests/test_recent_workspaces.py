"""Recent-workspace list: it has to survive the restart it exists to save.

The menu it feeds replaces "edit the config by hand to change cluster", so the
failure modes that matter are the quiet ones — an entry that points at a file
that is gone, the same workspace listed twice under two spellings, or a list
that resets because the home directory was not writable.
"""

from __future__ import annotations

import json

from apex.config.recent_workspaces import (
    MAX_ENTRIES, display_label, forget, load_recent, remember,
)


def _config(tmp_path, name: str):
    workspace = tmp_path / name
    workspace.mkdir()
    config = workspace / "apex_config.json"
    config.write_text("{}", encoding="utf-8")
    return config


def test_most_recent_comes_first(tmp_path):
    store = tmp_path / "recent.json"
    first = _config(tmp_path, "M13")
    second = _config(tmp_path, "NGC6811")
    remember(first, store)
    remember(second, store)
    assert load_recent(store) == [second, first]


def test_reopening_moves_it_up_instead_of_duplicating(tmp_path):
    store = tmp_path / "recent.json"
    first = _config(tmp_path, "M13")
    second = _config(tmp_path, "NGC6811")
    remember(first, store)
    remember(second, store)
    remember(first, store)
    assert load_recent(store) == [first, second]


def test_a_workspace_that_vanished_is_not_offered(tmp_path):
    store = tmp_path / "recent.json"
    kept = _config(tmp_path, "M13")
    gone = _config(tmp_path, "M67")
    remember(kept, store)
    remember(gone, store)
    gone.unlink()
    assert load_recent(store) == [kept]


def test_the_list_is_bounded(tmp_path):
    store = tmp_path / "recent.json"
    for index in range(MAX_ENTRIES + 5):
        remember(_config(tmp_path, f"target{index}"), store)
    assert len(load_recent(store)) == MAX_ENTRIES


def test_forget_removes_one_entry(tmp_path):
    store = tmp_path / "recent.json"
    first = _config(tmp_path, "M13")
    second = _config(tmp_path, "NGC6811")
    remember(first, store)
    remember(second, store)
    assert forget(second, store) == [first]
    assert load_recent(store) == [first]


def test_a_corrupt_store_reads_as_empty_rather_than_raising(tmp_path):
    store = tmp_path / "recent.json"
    store.write_text("{ this is not json", encoding="utf-8")
    assert load_recent(store) == []


def test_an_unwritable_store_still_returns_the_new_order(tmp_path):
    """Losing the history must not stop the workspace from opening."""
    store = tmp_path / "nodir" / "sub" / "recent.json"
    store.parent.mkdir(parents=True)
    store.parent.chmod(0o555)
    try:
        config = _config(tmp_path, "M13")
        assert remember(config, store) == [config]
    finally:
        store.parent.chmod(0o755)


def test_the_label_names_the_target_not_the_file(tmp_path):
    """Every workspace's config has the same file name; the folder differs."""
    config = _config(tmp_path, "M13")
    label = display_label(config)
    assert "M13" in label and "apex_config.json" in label


def test_the_stored_form_is_a_plain_list_of_paths(tmp_path):
    store = tmp_path / "recent.json"
    config = _config(tmp_path, "M13")
    remember(config, store)
    data = json.loads(store.read_text(encoding="utf-8"))
    assert data["workspaces"] == [str(config.absolute())]
