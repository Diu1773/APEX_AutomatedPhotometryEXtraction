"""Re-run Step 8 on an already-injected trial with one setting changed.

Asking whether a fit setting helps is only answerable if nothing else moves.
Re-running the whole artificial-star pipeline would re-detect and re-inject,
so a difference in the answer could come from a different set of stars rather
than from the setting. This copies the trial's finished detection and forced
catalogue, changes the one key under test, and runs Step 8 alone. Everything
upstream is byte-identical by construction.

    run_apex_variant.py --run-dir <run> --trial 1 --tag prof \
        --parameter-file <base apex_config.json> \
        --set psf.profile_error_frac=0.05

The output lands in ``<run>/apex<tag>_trial<N>/`` — the layout
``compare_three_engines.py`` already reads.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).absolute().parents[2]
PYTHON = str(REPO / ".venv-deploy" / "Scripts" / "python.exe")

# Copied from the trial so Step 8 starts from the same place it did there.
# ``cache`` is deliberately not copied: its entries are keyed by result_dir.
REUSED_SUBDIRS = ("step4_detection", "step7_forced_phot")


def _parse_override(text: str) -> tuple[list[str], object]:
    """``psf.profile_error_frac=0.05`` -> (['psf', 'profile_error_frac'], 0.05)."""
    if "=" not in text:
        raise SystemExit(f"--set 은 key=value 형식이어야 한다: {text!r}")
    key, _, raw = text.partition("=")
    path = [part for part in key.split(".") if part]
    if not path:
        raise SystemExit(f"--set 의 키가 비었다: {text!r}")
    try:
        value = json.loads(raw)
    except json.JSONDecodeError:
        value = raw
    return path, value


def _apply(config: dict, path: list[str], value: object) -> None:
    node = config
    for part in path[:-1]:
        child = node.get(part)
        if not isinstance(child, dict):
            child = {}
            node[part] = child
        node = child
    node[path[-1]] = value


def _prepare_workspace(run_dir: Path, trial: int, tag: str,
                       parameter_file: Path, overrides: list[str],
                       force: bool) -> tuple[Path, Path, str]:
    source = run_dir / f"trial_{trial:04d}"
    if not source.is_dir():
        raise SystemExit(f"주입 회차가 없다: {source}")
    frames = sorted((source / "data").glob("*.fit*"))
    if len(frames) != 1:
        raise SystemExit(f"data 에 프레임이 정확히 하나여야 한다: {frames}")

    workspace = run_dir / f"apex{tag}_trial{trial}"
    if workspace.exists():
        if not force:
            raise SystemExit(f"이미 있다 (덮어쓰려면 --force): {workspace}")
        shutil.rmtree(workspace)
    (workspace / "data").mkdir(parents=True)
    shutil.copy2(frames[0], workspace / "data" / frames[0].name)

    result = workspace / "result"
    result.mkdir()
    for name in REUSED_SUBDIRS:
        origin = source / "result" / name
        if not origin.is_dir():
            raise SystemExit(f"회차에 {name} 이 없다: {origin}")
        shutil.copytree(origin, result / name)

    config = json.loads(parameter_file.read_text(encoding="utf-8"))
    for text in overrides:
        _apply(config, *_parse_override(text))
    _apply(config, ["io", "data_dir"], str(workspace / "data"))
    _apply(config, ["io", "result_dir"], str(result))
    target = workspace / "apex_config.json"
    target.write_text(json.dumps(config, indent=2, ensure_ascii=False),
                      encoding="utf-8")

    # The workspace-identity accident of 2026-08-06 in miniature: a config that
    # still points at the source workspace would write Step 8 output back over
    # the run being compared against. Read it back rather than trust the write.
    written = json.loads(target.read_text(encoding="utf-8"))
    for key in ("data_dir", "result_dir"):
        if Path(written["io"][key]).absolute() != (
                workspace / ("data" if key == "data_dir" else "result")).absolute():
            raise SystemExit(f"io.{key} 가 변형 작업공간을 가리키지 않는다")
    return workspace, target, frames[0].name


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, type=Path,
                    help="인공별 주입 실행 폴더 (summary.json 이 있는 곳)")
    ap.add_argument("--trial", type=int, default=1)
    ap.add_argument("--tag", required=True,
                    help="변형 이름. 결과는 apex<tag>_trial<N>/ 에 쌓인다")
    ap.add_argument("--parameter-file", required=True, type=Path,
                    help="기준 apex_config.json. 여기에 --set 을 얹는다")
    ap.add_argument("--set", dest="overrides", action="append", default=[],
                    metavar="KEY=VALUE",
                    help="점 표기 설정 덮어쓰기. 여러 번 쓸 수 있다")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    run_dir = args.run_dir.absolute()
    meta = json.loads((run_dir / "summary.json").read_text(encoding="utf-8"))["metadata"]
    workspace, config, frame = _prepare_workspace(
        run_dir, args.trial, args.tag, args.parameter_file.absolute(),
        args.overrides, args.force)

    command = [
        PYTHON, "-X", "utf8", str(REPO / "validation" / "run_real_gui_psf.py"),
        "--data-dir", str(workspace / "data"),
        "--result-dir", str(workspace / "result"),
        "--parameter-file", str(config),
        "--fit-engine", "apex_iterative",
        "--fitter-max-iter", "8",
        "--fit-shape-fwhm-mult", str(meta["fit_shape_fwhm_mult"]),
        "--fit-window-mode", meta["fit_window_mode"],
        "--fit-encircled-energy", str(meta["fit_encircled_energy"]),
        "--postfit-qfit-noise-max", str(meta["postfit_qfit_noise_max"]),
        "--residual-passes", "2",
        "--epsf-max-stars", "0",
        "--epsf-contamination", "on",
        "--flux-scale", meta["flux_scale"],
        "--use-grouper", meta["use_grouper"],
        "--grouper-max-size", str(meta["grouper_max_size"]),
        "--grouper-radius-fwhm", str(meta["grouper_radius_fwhm"]),
        "--forced-match-radius-fwhm", str(meta["forced_match_radius_fwhm"]),
        "--core-cut", "off",
        "--core-center-mode", "auto",
        "--skip-step4",                      # detection copied from the trial
        frame,
    ]
    print(f"작업공간 {workspace}\n  변경 {args.overrides or ['(없음)']}", flush=True)
    started = time.perf_counter()
    log = workspace / "run.log"
    with log.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(command, cwd=str(REPO), text=True,
                                   stdout=handle, stderr=subprocess.STDOUT)
    print(f"step8 {'완료' if completed.returncode == 0 else '실패'} · "
          f"{time.perf_counter() - started:.0f}s · 로그 {log}", flush=True)
    return completed.returncode


if __name__ == "__main__":
    raise SystemExit(main())
