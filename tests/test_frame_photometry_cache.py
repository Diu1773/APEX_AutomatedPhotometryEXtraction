"""The light-curve frame cache must survive a change of result directory (O3).

The builder previously kept one directory at a time and cleared everything
when it changed, so a multi-night workspace — which interleaves frames from
several result dirs — re-read the lot on every switch. Measured: three stars
over 124 frames cost 372 frame reads instead of 124, amplification exactly
linear in stars (benchmark/perf/20260809/b3_lc_load.json).

Bounding by bytes rather than entries matters for the same reason the read
amplification did: a frame table scales with the master catalog, so a count
that is safe at 1,357 stars is not safe at 5,000.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from apex.utils.photometry_loader import FramePhotometryCache


def _table(rows: int = 100, cols: int = 8) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    return pd.DataFrame({f"c{i}": rng.normal(size=rows) for i in range(cols)})


@pytest.fixture
def cache():
    return FramePhotometryCache(budget_mb=64)


def test_a_stored_table_comes_back(cache):
    df = _table()
    cache.put("/night1", "f001.fits", df)
    assert cache.get("/night1", "f001.fits") is df
    assert cache.hits == 1


def test_a_miss_is_none_not_an_error(cache):
    assert cache.get("/night1", "missing.fits") is None
    assert cache.misses == 1


def test_two_result_dirs_coexist(cache):
    """The regression: alternating directories used to clear the cache."""
    a, b = _table(), _table()
    cache.put("/night1", "f001.fits", a)
    cache.put("/night2", "f001.fits", b)

    # Same frame name, different night — both must still be there.
    assert cache.get("/night1", "f001.fits") is a
    assert cache.get("/night2", "f001.fits") is b
    assert len(cache) == 2


def test_alternating_directories_never_re_read(cache):
    """Interleave two nights the way a merged workspace does."""
    for night in ("/night1", "/night2"):
        for i in range(20):
            cache.put(night, f"f{i:03d}.fits", _table(rows=10))

    misses_before = cache.misses
    for i in range(20):                       # alternate, twice over
        assert cache.get("/night1", f"f{i:03d}.fits") is not None
        assert cache.get("/night2", f"f{i:03d}.fits") is not None
    assert cache.misses == misses_before


def test_paths_are_normalised(cache):
    df = _table()
    cache.put("/night1/", "f001.fits", df)
    assert cache.get("/night1", "f001.fits") is df


def test_contains_takes_a_pair(cache):
    cache.put("/night1", "f001.fits", _table())
    assert ("/night1", "f001.fits") in cache
    assert ("/night2", "f001.fits") not in cache


# ── the bound ──────────────────────────────────────────────────────────────

def test_the_budget_is_bytes_not_entries():
    """Same entry count, different row count — only the big one evicts."""
    small = FramePhotometryCache(budget_mb=1)
    for i in range(10):
        small.put("/night1", f"f{i:03d}.fits", _table(rows=100))
    assert small.evictions == 0

    big = FramePhotometryCache(budget_mb=1)
    for i in range(10):
        big.put("/night1", f"f{i:03d}.fits", _table(rows=20_000))
    assert big.evictions > 0
    # The budget is respected unless a single table exceeds it on its own —
    # the cache always keeps one entry so a caller is never handed nothing.
    assert big.stats()["bytes"] <= big.budget_bytes or len(big) == 1


def test_eviction_is_least_recently_used():
    # Budget sized from a real table so exactly three fit and a fourth evicts.
    one = _table(rows=5_000)
    entry_mb = one.memory_usage(deep=True).sum() / 1e6
    cache = FramePhotometryCache(budget_mb=entry_mb * 3.5)
    for i in range(3):
        cache.put("/n", f"f{i}.fits", _table(rows=5_000))
    assert cache.evictions == 0, "three entries should fit"

    cache.get("/n", "f0.fits")          # f0 becomes most recent
    cache.put("/n", "new.fits", _table(rows=5_000))

    assert cache.get("/n", "f0.fits") is not None, "MRU entry was evicted"
    assert cache.get("/n", "f1.fits") is None, "LRU entry should have gone"


def test_the_cache_never_empties_itself_completely():
    """One table larger than the whole budget must still be usable."""
    cache = FramePhotometryCache(budget_mb=1)
    huge = _table(rows=200_000)
    cache.put("/n", "huge.fits", huge)
    assert cache.get("/n", "huge.fits") is huge


def test_replacing_a_key_does_not_double_count():
    cache = FramePhotometryCache(budget_mb=64)
    cache.put("/n", "f.fits", _table(rows=1_000))
    first = cache.stats()["bytes"]
    cache.put("/n", "f.fits", _table(rows=1_000))
    assert cache.stats()["bytes"] == first
    assert len(cache) == 1


def test_clear_resets_the_contents_but_keeps_the_counters(cache):
    cache.put("/n", "f.fits", _table())
    cache.get("/n", "f.fits")
    cache.clear()
    assert len(cache) == 0 and cache.stats()["bytes"] == 0
    assert cache.hits == 1          # history is diagnostics, not state


def test_the_default_budget_is_derived_and_bounded():
    cache = FramePhotometryCache()
    budget_mb = cache.budget_bytes / 1e6
    assert FramePhotometryCache._BUDGET_MIN_MB <= budget_mb <= FramePhotometryCache._BUDGET_MAX_MB
