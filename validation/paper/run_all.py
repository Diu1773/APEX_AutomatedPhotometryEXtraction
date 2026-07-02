"""Run the full APEX paper-figure suite in order.

Usage (always with the deploy venv interpreter):
    .venv-deploy/Scripts/python.exe validation/paper/run_all.py
    .venv-deploy/Scripts/python.exe validation/paper/run_all.py --only 1 5
    .venv-deploy/Scripts/python.exe validation/paper/run_all.py --fast   # skip fig3 sweeps
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

PAPER = Path(__file__).resolve().parent

STEPS: list[tuple[str, str]] = [
    ("data", "_make_canonical_data.py"),   # canonical injection dataset (~1 min)
    ("1", "fig1_completeness.py"),
    ("2", "fig2_error_model.py"),
    ("3", "fig3_parameter_sweep.py"),      # slow: ~19 benchmark runs (~10 min)
    ("4", "fig4_crosscheck_sep.py"),
    ("5", "fig5_crosscheck_iraf.py"),
    ("6", "fig6_qc_validation.py"),
    ("assemble", "assemble_figures.py"),
]

SLOW = {"3"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="*", default=None,
                        help="run only these step keys (e.g. --only 2 5)")
    parser.add_argument("--fast", action="store_true", help="skip slow steps (fig3)")
    args = parser.parse_args()

    selected = [
        (key, script) for key, script in STEPS
        if (args.only is None or key in set(args.only))
        and not (args.fast and key in SLOW)
    ]
    # canonical data is required by fig1; keep it unless explicitly excluded
    if args.only is None and not (PAPER / "data" / "canonical_summary.json").exists():
        pass  # "data" step is already first in STEPS

    failures: list[str] = []
    for key, script in selected:
        path = PAPER / script
        if not path.exists():
            print(f"[skip] {script} (not present)")
            continue
        print(f"[run ] {script} ...", flush=True)
        t0 = time.perf_counter()
        proc = subprocess.run([sys.executable, str(path)], cwd=str(PAPER.parents[1]))
        dt = time.perf_counter() - t0
        status = "ok" if proc.returncode == 0 else f"FAIL rc={proc.returncode}"
        print(f"[{status:>4}] {script} ({dt:.1f}s)")
        if proc.returncode != 0:
            failures.append(script)

    if failures:
        print("FAILED:", ", ".join(failures))
        return 1
    print("all selected steps completed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
