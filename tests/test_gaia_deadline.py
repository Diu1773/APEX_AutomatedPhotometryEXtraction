"""A timed-out Gaia query must not keep the process alive.

`tap_query_with_deadline` frees the caller after `deadline_s`, but freeing the
caller is not enough: the abandoned query has to be on a thread the
interpreter is willing to leave behind. An earlier version used a
`ThreadPoolExecutor` and dropped it with `shutdown(wait=False)` — futures'
workers are non-daemon and `concurrent.futures` installs an `atexit` hook that
joins them, so a benchmark run whose pipeline finished in 226 s took 2,424 s
to exit, waiting on a query nobody was reading.

The last test here is the one that matters: it starts a real interpreter,
abandons a slow query, and requires the process to be gone long before the
query would have finished.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
import threading
import time

import pytest

from apex.utils.gaia_catalog_service import tap_query_with_deadline


class _Service:
    """Stands in for astroquery's Gaia/TAP service."""

    def __init__(self, delay: float, result="rows", error=None):
        self.delay, self.result, self.error = delay, result, error

    def launch_job_async(self, adql, **kwargs):
        time.sleep(self.delay)
        if self.error is not None:
            raise self.error
        service = self

        class _Job:
            def get_results(self):
                return service.result

        return _Job()


def test_a_fast_query_returns_its_rows():
    assert tap_query_with_deadline(_Service(0.0), "SELECT 1",
                                   deadline_s=5) == "rows"


def test_a_slow_query_raises_timeout_at_the_deadline():
    started = time.perf_counter()
    with pytest.raises(TimeoutError, match="timed out"):
        tap_query_with_deadline(_Service(10.0), "SELECT 1", deadline_s=0.3)
    # The caller is freed at the deadline, not when the server answers.
    assert time.perf_counter() - started < 3.0


def test_the_error_message_maps_to_a_retryable_failure():
    """`gaia_failure_reason` keys off "timed out" so callers try another mirror."""
    from apex.utils.gaia_catalog_service import gaia_failure_reason

    with pytest.raises(TimeoutError) as caught:
        tap_query_with_deadline(_Service(10.0), "SELECT 1", deadline_s=0.2)
    assert gaia_failure_reason(caught.value) == "timeout"


def test_a_query_error_reaches_the_caller():
    boom = ValueError("503 server_down")
    with pytest.raises(ValueError, match="server_down"):
        tap_query_with_deadline(_Service(0.0, error=boom), "SELECT 1",
                                deadline_s=5)


def test_the_worker_thread_is_a_daemon():
    before = {t.ident for t in threading.enumerate()}
    with pytest.raises(TimeoutError):
        tap_query_with_deadline(_Service(5.0), "SELECT 1", deadline_s=0.2)
    leaked = [t for t in threading.enumerate()
              if t.ident not in before and t.name == "apex-tap-query"]
    assert leaked, "the abandoned query thread should still be running"
    assert all(t.daemon for t in leaked), \
        "a non-daemon thread blocks interpreter exit (atexit joins it)"


def test_the_process_exits_without_waiting_for_the_query():
    """The regression itself: abandon a 60 s query, exit in well under that."""
    program = textwrap.dedent("""
        import sys, time
        sys.path.insert(0, %r)
        from apex.utils.gaia_catalog_service import tap_query_with_deadline

        class S:
            def launch_job_async(self, adql, **kw):
                time.sleep(60)
                raise AssertionError("should never be reached")

        try:
            tap_query_with_deadline(S(), "SELECT 1", deadline_s=0.5)
        except TimeoutError:
            pass
        print("done", flush=True)
    """) % str(__import__("pathlib").Path(__file__).resolve().parents[1])

    started = time.perf_counter()
    proc = subprocess.run([sys.executable, "-c", program],
                          capture_output=True, text=True, timeout=45)
    elapsed = time.perf_counter() - started

    assert proc.returncode == 0, proc.stderr
    assert "done" in proc.stdout
    assert elapsed < 20, (
        f"process took {elapsed:.1f}s to exit; it waited for the abandoned "
        f"query instead of leaving it behind")
