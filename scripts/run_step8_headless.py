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


def _console_text(value) -> str:
    encoding = getattr(sys.stdout, "encoding", None) or "utf-8"
    return str(value).encode(encoding, errors="backslashreplace").decode(encoding)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--params", default=str(REPO / "parameters.toml"))
    args = parser.parse_args()

    from apex.config.parameters_cmd import read_params
    from apex.analysis.cmd.psf_photometry_runner import (
        PsfPhotometryRunner,
        build_psf_output_signature,
        export_psf_qc_products,
        write_psf_output_signature,
    )

    params = read_params(args.params)
    P = params.P
    result_dir = Path(P.result_dir)

    tsvs = sorted(glob.glob(str(result_dir / "step7_forced_phot" / "photometry_*.tsv")))
    file_list = [Path(t).name[len("photometry_"):-len(".tsv")] for t in tsvs]
    print(f"frames: {len(file_list)}")

    worker = PsfPhotometryRunner(
        file_list=file_list,
        params=params,
        data_dir=Path(P.data_dir),
        result_dir=result_dir,
        cache_dir=Path(P.cache_dir) if Path(P.cache_dir).is_absolute() else result_dir / P.cache_dir,
        use_cropped=False,
    )
    def _print_interesting(message):
        if any(k in str(message).lower()
               for k in ("error", "fail", "complete", "saved", "[")):
            print(f"[LOG] {_console_text(message)}", flush=True)

    worker.on_log.subscribe(_print_interesting)
    worker.on_error.subscribe(
        lambda stage, m: print(f"[ERROR] {stage}: {_console_text(m)}", flush=True))
    done: dict = {}
    worker.on_finished.subscribe(lambda d: done.update(d if isinstance(d, dict) else {"_": d}))

    t0 = time.perf_counter()
    worker.run()
    elapsed = time.perf_counter() - t0
    complete = (
        len(file_list) > 0
        and int(done.get("processed", 0)) == len(file_list)
        and int(done.get("stopped", 0)) == 0
    )
    if complete:
        signature = build_psf_output_signature(
            params,
            file_list,
            use_cropped=False,
            cache_dir=(
                Path(P.cache_dir)
                if Path(P.cache_dir).is_absolute()
                else result_dir / P.cache_dir
            ),
        )
        signature_path = write_psf_output_signature(result_dir, signature)
        print(f"[signature] {signature_path}")
    else:
        signature_path = result_dir / "cmd_psf" / "psf_output_signature.json"
        signature_path.unlink(missing_ok=True)
        print(
            "[signature] not written: "
            f"processed={done.get('processed', 0)}/{len(file_list)} "
            f"stopped={done.get('stopped', 0)}"
        )
    # QC 산출물 — 창이 만드는 것과 같은 함수를 쓴다. Step 10 은 워커 안에서
    # 내보내므로 헤드리스도 그림이 나오는데 Step 8 은 창에만 있어서 헤드리스
    # 검증에 PSF 리포트가 하나도 안 남았다.
    try:
        qc_paths = export_psf_qc_products(
            result_dir / "cmd_psf", params=params, result_dir=result_dir
        )
        for p in qc_paths:
            print(f"[qc] {p.name}")
        if not qc_paths:
            print("[qc] 산출물 없음 (입력 표가 비었다)")
    except Exception as exc:
        print(f"[qc] 실패: {_console_text(exc)}")

    print(f"[done] elapsed {elapsed:.1f}s  keys={sorted(done)[:8]}")
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
