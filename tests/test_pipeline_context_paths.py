"""`--result-dir` must take the run's caches with it.

`read_params` resolves `cache_dir` against the config's `result_dir`, so an
override applied afterwards used to leave the detection/WCS caches pointing at
the config's tree.  Two runs with different `--result-dir` then shared one
cache and could reuse each other's per-file results — the "stale cache mixed
into new output" failure the performance plan forbids outright.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from apex.pipeline.context import RunContext


@pytest.fixture
def workspace(tmp_path):
    """A minimal CMD config whose io paths point inside tmp_path."""
    configured = tmp_path / "configured"
    (configured / "sci").mkdir(parents=True)
    config = tmp_path / "apex_config.json"
    config.write_text(json.dumps({
        "io": {
            "data_dir": str(configured / "sci"),
            "result_dir": str(configured / "result"),
            "cache_dir": "cache",
        },
        "target": {"name": "TEST"},
    }, ensure_ascii=False), encoding="utf-8")
    return config, configured


def test_result_dir_override_moves_the_cache(workspace, tmp_path):
    config, configured = workspace
    elsewhere = tmp_path / "elsewhere"

    ctx = RunContext.build("cmd", config, result_dir=elsewhere)

    assert ctx.result_dir == elsewhere
    cache = Path(ctx.params.P.cache_dir)
    assert cache == elsewhere / "cache"
    # …and specifically NOT left behind in the configured tree.
    assert not cache.is_relative_to(configured)


def test_without_an_override_the_configured_cache_is_kept(workspace, tmp_path):
    config, configured = workspace

    ctx = RunContext.build("cmd", config)

    assert Path(ctx.params.P.cache_dir) == configured / "result" / "cache"


def test_two_overridden_runs_do_not_share_a_cache(workspace, tmp_path):
    config, _ = workspace

    a = RunContext.build("cmd", config, result_dir=tmp_path / "run_a")
    b = RunContext.build("cmd", config, result_dir=tmp_path / "run_b")

    assert Path(a.params.P.cache_dir) != Path(b.params.P.cache_dir)


def test_absolute_cache_outside_the_result_tree_is_left_alone(tmp_path):
    """An absolute cache_dir is a deliberate choice, not an accident."""
    shared = tmp_path / "shared_cache"
    config = tmp_path / "apex_config.json"
    (tmp_path / "sci").mkdir()
    config.write_text(json.dumps({
        "io": {
            "data_dir": str(tmp_path / "sci"),
            "result_dir": str(tmp_path / "result"),
            "cache_dir": str(shared),
        },
        "target": {"name": "TEST"},
    }, ensure_ascii=False), encoding="utf-8")

    ctx = RunContext.build("cmd", config, result_dir=tmp_path / "elsewhere")

    assert Path(ctx.params.P.cache_dir) == shared
