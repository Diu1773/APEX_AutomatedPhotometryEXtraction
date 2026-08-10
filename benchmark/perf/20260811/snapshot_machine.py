"""Record machine state around the Phase 3 clean run.

The clean run is the single source of the paper's performance numbers, so the
conditions it ran under have to be recorded with it — a wall time measured
while twelve other processes hold 13 GB is not the same number as one measured
idle, and the audit trail must say which it was.

Usage: snapshot_machine.py <label>
"""

import json
import shutil
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))

from apex.benchmark.resources import machine_state  # noqa: E402

label = sys.argv[1] if len(sys.argv) > 1 else "snapshot"
state = machine_state()
state["label"] = label
state["time"] = time.strftime("%Y-%m-%d %H:%M:%S")
state["e_free_gb"] = round(shutil.disk_usage("E:\\").free / 1e9, 1)
state["c_free_gb"] = round(shutil.disk_usage("C:\\").free / 1e9, 1)

out = Path(__file__).parent / f"machine_{label}.json"
out.write_text(json.dumps(state, indent=1), encoding="utf-8")
print(json.dumps(state, indent=1))
