"""Wall-time and peak-RSS instrumentation for benchmarks.

Every benchmark in the performance plan (docs/audit/APEX_PERF_DEV_PLAN.md)
records the same machine-readable envelope: wall time, peak RSS including
child processes, and an environment snapshot (git commit, package versions,
worker settings).  Without the child-process part a Step 0 run driven through
``subprocess`` would report the parent's few tens of MB and miss the actual
working set — which is the number the streaming optimisation is judged on.

``psutil`` is an optional dependency (``pip install apex-photometry[bench]``).
Without it, wall time and the environment snapshot still work; RSS fields are
``None`` rather than wrong.

Qt-free.
"""

from __future__ import annotations

import json
import os
import platform
import subprocess
import sys
import threading
import time
from contextlib import contextmanager
from pathlib import Path

try:
    import psutil
except Exception:  # pragma: no cover - optional dependency
    psutil = None

HAS_PSUTIL = psutil is not None

_POLL_S = 0.1
_PACKAGES = ("numpy", "scipy", "astropy", "photutils", "sep", "pandas",
             "astroscrappy", "bottleneck", "psutil")


def environment_snapshot() -> dict:
    """Git commit, interpreter, CPU and package versions — once per record."""
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10,
            cwd=Path(__file__).resolve().parents[2],
        ).stdout.strip() or None
    except Exception:  # pragma: no cover - git absent
        commit = None

    versions = {}
    for name in _PACKAGES:
        try:
            versions[name] = __import__(name).__version__
        except Exception:
            versions[name] = None

    return {
        "git_commit": commit,
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "cpu_logical": os.cpu_count(),
        "apex_max_workers_env": os.environ.get("APEX_MAX_WORKERS"),
        "packages": versions,
    }


class _RssPoller(threading.Thread):
    """Samples RSS of a process tree until stopped; keeps the maximum."""

    def __init__(self, pid: int, poll_s: float = _POLL_S):
        super().__init__(daemon=True)
        self._proc = psutil.Process(pid)
        self._poll_s = poll_s
        self._halt = threading.Event()
        self.peak_self = 0
        self.peak_total = 0

    def run(self) -> None:  # pragma: no cover - timing-dependent
        while not self._halt.is_set():
            try:
                rss = self._proc.memory_info().rss
                total = rss
                for child in self._proc.children(recursive=True):
                    try:
                        total += child.memory_info().rss
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        continue
                self.peak_self = max(self.peak_self, rss)
                self.peak_total = max(self.peak_total, total)
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                break
            self._halt.wait(self._poll_s)

    def stop(self) -> None:
        self._halt.set()
        self.join(timeout=2.0)


@contextmanager
def measure(label: str, extra: dict | None = None, poll_s: float = _POLL_S):
    """Measure the current process (and any children it spawns) around a block.

    Yields the metrics dict so the block can add its own fields; ``wall_s`` and
    the RSS peaks are filled in on exit.
    """
    metrics: dict = {
        "label": label,
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "env": environment_snapshot(),
    }
    if extra:
        metrics.update(extra)

    poller = None
    if HAS_PSUTIL:
        poller = _RssPoller(os.getpid(), poll_s)
        poller.start()

    t0 = time.perf_counter()
    try:
        yield metrics
    finally:
        metrics["wall_s"] = round(time.perf_counter() - t0, 3)
        if poller is not None:
            poller.stop()
            metrics["peak_rss_mb"] = round(poller.peak_self / 1e6, 1)
            metrics["peak_rss_total_mb"] = round(poller.peak_total / 1e6, 1)
        else:  # pragma: no cover - psutil absent
            metrics["peak_rss_mb"] = None
            metrics["peak_rss_total_mb"] = None


def measure_command(cmd: list, *, label: str, env: dict | None = None,
                    cwd: str | None = None, poll_s: float = _POLL_S) -> tuple[int, dict]:
    """Run a command and measure its whole process tree.

    Returns ``(returncode, metrics)``.  Output is inherited so long pipeline
    runs stay observable in the launcher's log.
    """
    metrics: dict = {
        "label": label,
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "cmd": [str(c) for c in cmd],
        "env": environment_snapshot(),
    }
    t0 = time.perf_counter()
    proc = subprocess.Popen(cmd, env=env, cwd=cwd)

    poller = None
    if HAS_PSUTIL:
        try:
            poller = _RssPoller(proc.pid, poll_s)
            poller.start()
        except Exception:  # pragma: no cover - race with a fast child
            poller = None

    returncode = proc.wait()
    metrics["wall_s"] = round(time.perf_counter() - t0, 3)
    metrics["returncode"] = returncode
    if poller is not None:
        poller.stop()
        # For a subprocess the child tree IS the workload; self==the child.
        metrics["peak_rss_mb"] = round(poller.peak_self / 1e6, 1)
        metrics["peak_rss_total_mb"] = round(poller.peak_total / 1e6, 1)
    else:
        metrics["peak_rss_mb"] = None
        metrics["peak_rss_total_mb"] = None
    return returncode, metrics


def save_metrics(metrics: dict, path: Path) -> Path:
    """Write one metrics record as pretty JSON (UTF-8), creating parents."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metrics, ensure_ascii=False, indent=1),
                    encoding="utf-8")
    return path
