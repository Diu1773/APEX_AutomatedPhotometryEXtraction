"""The measurements that need an idle machine, run back to back.

Everything here was already tried on a loaded machine and came back
unresolvable: detect's parallel gain sat inside its own noise, and the
three-policy comparison would have landed the same way. Both are cheap to
repeat and expensive to misread, so they wait for a quiet machine rather than
being quoted from a busy one.

Each stage checkpoints per run, so an interrupted night keeps what it finished.

Usage:  python -X utf8 benchmark/overnight.py
"""

from __future__ import annotations

import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from apex.benchmark import resources  # noqa: E402

STAGES = [
    ("three worker policies, steps 1-7, 2 repeats",
     ["benchmark/o2_arms.py", "arms", "--repeats", "2"]),
    ("detect/wcs sweep, steps 1-5, 3 repeats",
     ["benchmark/o2_arms.py", "sweep", "--workers", "1,2,4,8", "--repeats", "3"]),
]


def main() -> int:
    print(f"=== overnight batch started {time.strftime('%Y-%m-%d %H:%M:%S')} ===",
          flush=True)
    state = resources.machine_state()
    print(f"machine: {state.get('ram_available_mb', 0):,} MB free "
          f"({state.get('ram_available_pct', 0):.0f}%), "
          f"{state.get('processes_over_200mb', 0)} processes over 200 MB",
          flush=True)
    resources.warn_if_busy()

    failed = []
    for label, argv in STAGES:
        print(f"\n=== {label} — {time.strftime('%H:%M:%S')} ===", flush=True)
        started = time.perf_counter()
        rc = subprocess.call([sys.executable, "-X", "utf8", "-u", *argv],
                             cwd=str(REPO))
        mins = (time.perf_counter() - started) / 60.0
        print(f"=== {label}: rc={rc} in {mins:.0f} min ===", flush=True)
        if rc != 0:
            failed.append(label)

    print(f"\n=== overnight batch done {time.strftime('%H:%M:%S')} ===", flush=True)
    if failed:
        print("failed stages: " + "; ".join(failed), flush=True)
    final = resources.machine_state()
    print(f"machine at end: {final.get('ram_available_mb', 0):,} MB free",
          flush=True)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
