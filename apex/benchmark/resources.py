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


def machine_state() -> dict:
    """How busy the machine was — the half of reproducibility that git can't fix.

    A step-7 run measured 1.74x slower than the same code a day earlier, with
    identical inputs and an identical step-4 time.  The cause was a game plus
    its anti-cheat holding ~1.5 GB and burning ~4.7 CPU-hours in the
    background: free RAM had fallen from ~12 GB to 3.5 GB, so the frames that
    step 4 leaves in the OS page cache were evicted before step 7 re-read them.
    Nothing in the recorded envelope could have shown that after the fact.

    Deliberately aggregate — free RAM and load are what invalidate a
    measurement; *what* the user was running is not the benchmark's business.
    """
    if not HAS_PSUTIL:  # pragma: no cover - optional extra
        return {}
    try:
        vm = psutil.virtual_memory()
        # A non-blocking first call returns 0.0; a short interval is the price
        # of a number that means anything.
        cpu_pct = psutil.cpu_percent(interval=0.3)
        heavy = sum(1 for p in psutil.process_iter(["memory_info"])
                    if (p.info["memory_info"] or None)
                    and p.info["memory_info"].rss > 200e6)
        return {
            "ram_total_mb": round(vm.total / 1e6),
            "ram_available_mb": round(vm.available / 1e6),
            "ram_available_pct": round(100.0 * vm.available / vm.total, 1),
            "cpu_percent": cpu_pct,
            "processes_over_200mb": heavy,
        }
    except Exception:  # pragma: no cover - psutil quirks must not break a run
        return {}


def environment_snapshot() -> dict:
    """Git commit, interpreter, CPU, package versions and machine load."""
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
        "machine": machine_state(),
    }


# Below this, a measurement is competing for RAM with something else and its
# absolute numbers should not be quoted (see machine_state).
IDLE_RAM_PCT = 40.0


def warn_if_busy(printer=print) -> dict:
    """Print a warning when the machine is too loaded to trust wall times."""
    state = machine_state()
    pct = state.get("ram_available_pct")
    if pct is not None and pct < IDLE_RAM_PCT:
        printer(f"WARNING: only {pct:.0f}% RAM free "
                f"({state['ram_available_mb']:,} MB of {state['ram_total_mb']:,} MB), "
                f"{state['processes_over_200mb']} processes over 200 MB. "
                f"Absolute wall times from this run are not comparable across "
                f"days — re-measure on an idle machine before quoting them.")
    return state


class _RssPoller(threading.Thread):
    """Samples RSS of a process tree until stopped; keeps the maximum."""

    def __init__(self, pid: int, poll_s: float = _POLL_S):
        super().__init__(daemon=True)
        self._proc = psutil.Process(pid)
        self._poll_s = poll_s
        self._halt = threading.Event()
        self.peak_self = 0
        self.peak_total = 0
        self.peak_uss = 0
        self.peak_uss_total = 0

    def run(self) -> None:  # pragma: no cover - timing-dependent
        while not self._halt.is_set():
            try:
                rss = self._proc.memory_info().rss
                total = rss
                try:
                    children = self._proc.children(recursive=True)
                except psutil.Error:
                    # A transient failure enumerating children must not end the
                    # watch: an earlier version broke out of the loop here and
                    # reported a process's start-up size (5 MB) as its peak.
                    children = []
                # RSS counts resident pages of memory-mapped files, which the OS
                # can evict on demand — a streaming reader that mmaps its inputs
                # looks like it uses more memory than one that allocates real
                # arrays. USS excludes those shared/file-backed pages, so it is
                # the number that reflects actual OOM pressure.
                try:
                    uss = self._proc.memory_full_info().uss
                except psutil.Error:
                    uss = 0
                uss_total = uss
                for child in children:
                    try:
                        total += child.memory_info().rss
                        uss_total += child.memory_full_info().uss
                    except psutil.Error:
                        continue
                self.peak_self = max(self.peak_self, rss)
                self.peak_total = max(self.peak_total, total)
                self.peak_uss = max(self.peak_uss, uss)
                self.peak_uss_total = max(self.peak_uss_total, uss_total)
            except psutil.NoSuchProcess:
                break                     # the process really is gone
            except psutil.Error:
                pass                      # transient; keep watching
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
            metrics["peak_uss_mb"] = round(poller.peak_uss / 1e6, 1) or None
        else:  # pragma: no cover - psutil absent
            metrics["peak_rss_mb"] = None
            metrics["peak_rss_total_mb"] = None
            metrics["peak_uss_mb"] = None


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
        # Report the whole tree: a venv's python.exe on Windows is a ~5 MB
        # launcher stub that execs the real interpreter as a child, so the
        # parent's own numbers describe the stub, not the workload.
        metrics["peak_rss_mb"] = round(poller.peak_total / 1e6, 1)
        metrics["peak_rss_total_mb"] = round(poller.peak_total / 1e6, 1)
        metrics["peak_rss_launcher_mb"] = round(poller.peak_self / 1e6, 1)
        metrics["peak_uss_mb"] = round(poller.peak_uss_total / 1e6, 1) or None
    else:
        metrics["peak_rss_mb"] = None
        metrics["peak_rss_total_mb"] = None
        metrics["peak_uss_mb"] = None
    return returncode, metrics


def save_metrics(metrics: dict, path: Path) -> Path:
    """Write one metrics record as pretty JSON (UTF-8), creating parents."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metrics, ensure_ascii=False, indent=1),
                    encoding="utf-8")
    return path
