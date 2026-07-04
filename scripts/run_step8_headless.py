"""Run Step-8 (PSF photometry) headless, without the GUI.

Reads parameters.toml (or --params) and re-processes every frame already
indexed under <result_dir>/step7_forced_phot/, executing the production
Step6PSFWorker (CMD Step-8 PSF photometry) synchronously.

    .venv-deploy/Scripts/python.exe scripts/run_step8_headless.py
    .venv-deploy/Scripts/python.exe scripts/run_step8_headless.py --params parameters_M5.toml
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
    from apex.gui.workflow.cmd.step8_psf_photometry import Step6PSFWorker

    params = read_params(args.params)
    P = params.P
    result_dir = Path(P.result_dir)

    tsvs = sorted(glob.glob(str(result_dir / "step7_forced_phot" / "photometry_*.tsv")))
    file_list = [Path(t).name[len("photometry_"):-len(".tsv")] for t in tsvs]
    print(f"frames: {len(file_list)}")

    worker = Step6PSFWorker(
        file_list=file_list,
        params=params,
        data_dir=Path(P.data_dir),
        result_dir=result_dir,
        cache_dir=Path(P.cache_dir) if Path(P.cache_dir).is_absolute() else result_dir / P.cache_dir,
        use_cropped=False,
    )
    worker.log.connect(lambda m: print(f"[LOG] {m}", flush=True) if any(
        k in str(m).lower() for k in ("error", "fail", "complete", "saved", "[")
    ) else None)
    worker.error.connect(lambda m: print(f"[ERROR] {m}", flush=True))
    done: dict = {}
    worker.finished.connect(lambda d: done.update(d if isinstance(d, dict) else {"_": d}))

    t0 = time.perf_counter()
    worker.run()
    print(f"[done] elapsed {time.perf_counter() - t0:.1f}s  keys={sorted(done)[:8]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
