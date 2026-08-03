"""Run the final 17-figure APEX paper suite in manuscript order.

Usage (always with the deploy venv interpreter):
    .venv-deploy/Scripts/python.exe validation/paper/run_all.py
    .venv-deploy/Scripts/python.exe validation/paper/run_all.py --only 1 5 17
    .venv-deploy/Scripts/python.exe validation/paper/run_all.py --fast   # skip parameter sweep
"""

from __future__ import annotations

import argparse
import subprocess
import sys
import time
from pathlib import Path

PAPER = Path(__file__).absolute().parent

STEPS: list[tuple[str, str]] = [
    ("data", "_make_canonical_data.py"),   # canonical injection dataset (~1 min)
    ("1", "fig_architecture.py"),
    ("2", "fig_calibration_step0.py"),
    ("3", "fig11_detector.py"),
    ("4", "fig12_preproc_crosscheck.py"),
    ("5", "fig13_cross_instrument.py"),
    ("6", "fig6_qc_validation.py"),
    ("7", "fig_completeness_realvssynth.py"),
    ("8", "fig_detection_threshold.py"),
    ("9", "fig_wcs_engines.py"),
    ("10", "fig2_error_model.py"),
    ("11", "fig3_parameter_sweep.py"),     # slow: benchmark sweeps (~10 min)
    ("12", "fig_photometry_crosschecks.py"),
    ("13", "fig_psf_validation.py"),
    ("14", "fig9_crowded_field.py"),
    ("15", "fig_external_validation.py"),
    ("16", "fig_timeseries_validation.py"),
    ("17", "fig_lc_yzboo.py"),
    ("render", "render_preview.py"),
]

SLOW = {"11"}
# Steps that require retained validation products or the external data volume
# and therefore cannot run in a source-only checkout.
NEEDS_DATA_VOLUME = {"2", "3", "4", "5", "7", "8", "9", "12", "13", "14", "15", "17"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="*", default=None,
                        help="run only these step keys (e.g. --only 2 5)")
    parser.add_argument("--fast", action="store_true", help="skip slow steps (Figure 11 parameter sweep)")
    args = parser.parse_args()

    selected = [
        (key, script) for key, script in STEPS
        if (args.only is None or key in set(args.only))
        and not (args.fast and key in SLOW)
    ]
    # Canonical data is required by several synthetic tests; keep it on full runs.
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
        proc = subprocess.run(
            [sys.executable, "-X", "utf8", str(path)],
            cwd=str(PAPER.parents[1]),
        )
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
