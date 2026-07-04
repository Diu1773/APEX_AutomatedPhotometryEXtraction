"""Run Step-7 (forced aperture photometry) headless, without the GUI.

Reads parameters.toml (or --params), re-processes every frame already indexed
under <result_dir>/step7_forced_phot/ (or, if that directory is empty, every
frame under the configured data_dir), and executes the production
ForcedPhotWorker synchronously.

    .venv-deploy/Scripts/python.exe scripts/run_step7_headless.py
    .venv-deploy/Scripts/python.exe scripts/run_step7_headless.py --params parameters_M5.toml
"""

from __future__ import annotations

import argparse
import glob
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--params", default=str(REPO / "parameters.toml"))
    args = parser.parse_args()

    from PyQt5.QtCore import QCoreApplication

    QCoreApplication(sys.argv)

    from apex.config.parameters_cmd import read_params
    from apex.gui.workflow.step7_forced_aperture_phot import ForcedPhotWorker

    params = read_params(args.params)
    P = params.P
    result_dir = Path(P.result_dir)

    tsvs = sorted(glob.glob(str(result_dir / "step7_forced_phot" / "photometry_*.tsv")))
    if tsvs:
        file_list = [Path(t).name[len("photometry_"):-len(".tsv")] for t in tsvs]
    else:
        file_list = sorted(p.name for p in Path(P.data_dir).glob("*.fit*"))
    print(f"frames: {len(file_list)}"
          f"  (first: {file_list[0] if file_list else '-'})")
    print(f"annulus_scale={P.fitsky_annulus_scale}  dannulus_scale={P.fitsky_dannulus_scale}")

    worker = ForcedPhotWorker(
        params=params,
        data_dir=Path(P.data_dir),
        result_dir=result_dir,
        cache_dir=Path(P.cache_dir) if Path(P.cache_dir).is_absolute() else result_dir / P.cache_dir,
        file_list=file_list,
    )
    worker.log.connect(lambda m: print(f"[LOG] {m}") if any(
        k in str(m).lower() for k in ("error", "fail", "complete", "saved", "apcorr")
    ) else None)
    worker.error.connect(lambda m: print(f"[ERROR] {m}"))
    done: dict = {}
    worker.finished.connect(lambda d: done.update(d if isinstance(d, dict) else {"_": d}))

    t0 = time.perf_counter()
    worker.run()  # synchronous
    print(f"[done] elapsed {time.perf_counter() - t0:.1f}s  keys={sorted(done)[:6]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
