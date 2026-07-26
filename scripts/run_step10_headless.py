"""Run Step-10 (CMD zeropoint calibration) headless, without the GUI.

Reads parameters.toml (or --params), backs up the key cmd_zeropoint products
into cmd_zeropoint/_pre_run_backup/, then executes the production
ZeropointCalibrationWorker synchronously and prints its log lines.

    .venv-deploy/Scripts/python.exe scripts/run_step10_headless.py
    .venv-deploy/Scripts/python.exe scripts/run_step10_headless.py --params other.toml
"""

from __future__ import annotations

import argparse
import shutil
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

BACKUP_NAMES = (
    "median_by_ID_filter_wide_cmd.csv",
    "median_by_ID_filter_wide.csv",
    "median_by_ID_filter_wide_raw.csv",
    "zp_fit_coefficients.csv",
    "gaia_sdss_calibrator_by_ID.csv",
    "gaia_cmd_comparison_summary.csv",
    "gaia_cmd_snr_sweep.csv",
    "gaia_cmd_drift_by_mag.csv",
    "zp_qc_summary.csv",
)


def _console_text(value) -> str:
    encoding = getattr(sys.stdout, "encoding", None) or "utf-8"
    return str(value).encode(encoding, errors="backslashreplace").decode(encoding)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--params", default=str(REPO / "parameters.toml"))
    parser.add_argument("--no-backup", action="store_true")
    args = parser.parse_args()

    from PyQt5.QtCore import QCoreApplication

    app = QCoreApplication.instance() or QCoreApplication(sys.argv)

    from apex.config.parameters_cmd import read_params
    from apex.gui.workflow.cmd.step10_zeropoint_calibration import (
        ZeropointCalibrationWorker,
    )

    params = read_params(args.params)
    P = params.P
    zp_dir = Path(P.result_dir) / "cmd_zeropoint"

    if not args.no_backup and zp_dir.exists():
        bak = zp_dir / "_pre_run_backup"
        bak.mkdir(exist_ok=True)
        for name in BACKUP_NAMES:
            src = zp_dir / name
            if src.exists():
                shutil.copyfile(src, bak / name)
        print(f"[backup] {bak}")

    worker = ZeropointCalibrationWorker(params, P.data_dir, P.result_dir, P.cache_dir)
    worker._log = lambda message: print(
        f"[LOG] {_console_text(message)}",
        flush=True,
    )
    t0 = time.perf_counter()
    worker.run()  # synchronous
    ok = dict(worker.last_summary)
    if worker.last_error:
        print(f"[ERROR] {_console_text(worker.last_error)}", flush=True)
    success = bool(ok.get("ok", False))
    print(f"[done] elapsed {time.perf_counter() - t0:.1f}s  ok={success}")
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
