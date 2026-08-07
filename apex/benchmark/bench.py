"""Performance benchmarks over fixed fixtures (``apex bench …``).

Implements the measurement half of docs/audit/APEX_PERF_DEV_PLAN.md:

* ``step0``  — B1: Step 0 calibration of an N-frame raw subset, in-process,
  so peak RSS is the calibration engine's own working set.
* ``sweep``  — B2: shared chain (Steps 1–7) at several worker counts, each in
  a fresh result dir, with a bit-identity digest of the Step 7 tables.
* ``repro``  — B4: the same run repeated at a fixed worker count; any drift
  between repeats is the pipeline's reproducibility noise floor, which parity
  tolerances must sit above.
* ``parity`` — delegates to :mod:`apex.benchmark.parity` (T0.5).

Fixtures are the plan's fixed set; scratch output lands on the validation
drive, metrics JSON lands in ``benchmark/perf/<date>/`` for committing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import sys
import time
from pathlib import Path

from apex.benchmark import resources

REPO = Path(__file__).resolve().parents[2]
VENV_PY = sys.executable

# ── fixed fixtures (docs/audit/APEX_PERF_DEV_PLAN.md T0.3) ──────────────────
RAW_M67 = Path(r"E:\observe_raw_Analysis\M67_20260208")
BIAS_POOL = Path(r"E:\bias")
DARK_POOL = Path(r"E:\darks")
SCI_NGC6811 = Path(r"E:\APEX_validation\reprocess\NGC6811\sci")
# The clean-run config is JSON (the config authority since 2026-08-06); the
# dated suffix keeps it from becoming the directory's implicit authority while
# still being explicitly passable to --config. The first overnight batch failed
# here by pointing at a TOML that a silenced `cp 2>/dev/null` never created.
CFG_NGC6811 = Path(r"E:\APEX_validation\reprocess\NGC6811\apex_config_20260807.json")
SCRATCH_ROOT = Path(r"E:\APEX_validation\bench")
TZ_OFFSET_HOURS = 9.0


def _out_root(args) -> Path:
    # getattr, not args.out_root: the option is declared with SUPPRESS so it
    # can appear on either side of the subcommand (see _global_options).
    out_root = getattr(args, "out_root", None)
    if out_root:
        return Path(out_root)
    # benchmark/runs is gitignored (bulky science-run artifacts); the perf
    # envelopes are small machine-readable records the plan requires to be
    # committed, so they live in a tracked sibling instead.
    return REPO / "benchmark" / "perf" / time.strftime("%Y%m%d")


def _fresh_dir(path: Path) -> Path:
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)
    return path


def _digest_step7(result_dir: Path) -> dict[str, str]:
    """SHA-256 per Step 7 photometry table — the bit-identity fingerprint."""
    digests: dict[str, str] = {}
    for p in sorted((result_dir / "step7_forced_phot").glob("photometry_*.tsv")):
        digests[p.name] = hashlib.sha256(p.read_bytes()).hexdigest()
    return digests


def _sum_detections(result_dir: Path) -> int | None:
    total, seen = 0, False
    for p in (result_dir / "step4_detection").glob("detect_*.json"):
        try:
            n = json.loads(p.read_text(encoding="utf-8")).get("n_sources")
        except Exception:
            continue
        if isinstance(n, (int, float)):
            total += int(n)
            seen = True
    return total if seen else None


# ── B1: Step 0 RSS vs frame count ───────────────────────────────────────────

def cmd_step0(args) -> int:
    from collections import Counter

    from apex.analysis import calibration_scan as scan
    from apex.analysis.calibration import CalibrationOptions
    from apex.analysis.calibration_run import ALL_NIGHTS, run_calibration

    n = int(args.frames)
    out = _fresh_dir(SCRATCH_ROOT / f"step0_M67_n{n:02d}")

    frames = []
    for root in (RAW_M67, BIAS_POOL, DARK_POOL):
        frames.extend(scan.scan_folder(str(root), tz_offset_hours=TZ_OFFSET_HOURS,
                                       warn=lambda m: None))
    lights = sorted((f for f in frames if f.ftype == "light"),
                    key=lambda f: str(f.path))
    if len(lights) < n:
        print(f"[bench] only {len(lights)} lights available (< {n})")
        return 2
    keep = {str(f.path) for f in lights[:n]}
    subset = [f for f in frames if f.ftype != "light" or str(f.path) in keep]
    print(f"[bench step0] n={n}  scanned={dict(Counter(f.ftype for f in subset))}")

    opts = CalibrationOptions(combine_method="median", pedestal_mode="none")
    with resources.measure(f"step0_M67_n{n:02d}",
                           extra={"n_lights": n, "scratch": str(out)}) as m:
        summary = run_calibration(subset, ALL_NIGHTS, out, opts, log=lambda s: None)
    m["n_calibrated"] = summary.get("n_calibrated")

    path = resources.save_metrics(m, _out_root(args) / f"step0_M67_n{n:02d}.json")
    print(f"[bench step0] n={n}  wall={m['wall_s']:.1f}s  "
          f"peak_rss={m['peak_rss_mb']} MB  -> {path}")
    return 0


# ── shared runner for sweep / repro ─────────────────────────────────────────

def _run_chain(result_dir: Path, workers: int, steps: str, label: str) -> dict:
    _fresh_dir(result_dir)
    cmd = [VENV_PY, "-X", "utf8", "-m", "apex.cli", "run", "--mode", "cmd",
           "--config", str(CFG_NGC6811), "--steps", steps,
           "--data-dir", str(SCI_NGC6811), "--result-dir", str(result_dir),
           "--force"]
    env = {**os.environ, "APEX_MAX_WORKERS": str(workers)}
    rc, metrics = resources.measure_command(cmd, label=label, env=env, cwd=str(REPO))
    metrics["workers"] = workers
    metrics["steps"] = steps
    if rc == 0:
        metrics["step7_digests"] = _digest_step7(result_dir)
        metrics["n_detect_sum"] = _sum_detections(result_dir)
    return metrics


def _identity_line(reference: dict, current: dict) -> str:
    if not reference or not current:
        return "digest=—"
    if reference == current:
        return "digest=IDENTICAL"
    differing = sum(1 for k in reference if current.get(k) != reference[k])
    missing = len(set(reference) ^ set(current))
    return f"digest=DIFFERS ({differing} files differ, {missing} missing)"


def cmd_sweep(args) -> int:
    workers = [int(w) for w in str(args.workers).split(",")]
    rows, reference = [], None
    for w in workers:
        metrics = _run_chain(SCRATCH_ROOT / f"sweep_w{w:02d}", w, args.steps,
                             f"sweep_w{w:02d}")
        if metrics.get("returncode") != 0:
            print(f"[bench sweep] w={w} FAILED rc={metrics.get('returncode')}")
            rows.append(metrics)
            continue
        if reference is None:
            reference = metrics.get("step7_digests")
        line = _identity_line(reference, metrics.get("step7_digests"))
        print(f"[bench sweep] w={w:>2}  wall={metrics['wall_s']:7.1f}s  "
              f"peak_rss_total={metrics['peak_rss_total_mb']} MB  "
              f"n_detect={metrics.get('n_detect_sum')}  {line}")
        rows.append(metrics)

    out = resources.save_metrics(
        {"label": "sweep", "workers": workers, "runs": rows},
        _out_root(args) / "sweep_NGC6811.json")
    print(f"[bench sweep] saved -> {out}")
    return 0


def cmd_repro(args) -> int:
    repeat, w = int(args.repeat), int(args.workers)
    rows, reference = [], None
    for i in range(repeat):
        metrics = _run_chain(SCRATCH_ROOT / f"repro_r{i}", w, args.steps,
                             f"repro_r{i}_w{w:02d}")
        if metrics.get("returncode") != 0:
            print(f"[bench repro] run {i} FAILED")
            rows.append(metrics)
            continue
        if reference is None:
            reference = metrics.get("step7_digests")
        line = _identity_line(reference, metrics.get("step7_digests"))
        print(f"[bench repro] run {i}  wall={metrics['wall_s']:7.1f}s  "
              f"n_detect={metrics.get('n_detect_sum')}  {line}")
        rows.append(metrics)

    out = resources.save_metrics(
        {"label": "repro", "workers": w, "repeat": repeat, "runs": rows},
        _out_root(args) / "repro_NGC6811.json")
    print(f"[bench repro] saved -> {out}")
    return 0


def cmd_parity(args) -> int:
    from apex.benchmark import parity

    return parity.main(args.parity_args)


# ── CLI ─────────────────────────────────────────────────────────────────────

SUBCOMMANDS = ("step0", "sweep", "repro", "parity")


def _global_options() -> argparse.ArgumentParser:
    """Options usable on either side of the subcommand.

    ``default=SUPPRESS`` is what makes that work: the same option lives on the
    top-level parser and on every subparser, and without SUPPRESS the
    subparser's ``None`` default would overwrite a value the top-level parser
    had already stored.  Read the result with ``getattr(args, …, None)``.
    """
    shared = argparse.ArgumentParser(add_help=False)
    shared.add_argument("--out-root", default=argparse.SUPPRESS,
                        help="metrics output dir (default benchmark/perf/<date>)")
    return shared


def build_parser() -> argparse.ArgumentParser:
    shared = _global_options()
    parser = argparse.ArgumentParser(
        prog="apex bench", parents=[shared],
        description="Performance benchmarks over fixed fixtures.")
    sub = parser.add_subparsers(dest="bench_cmd", required=True)

    p0 = sub.add_parser("step0", parents=[shared],
                        help="B1: Step 0 peak RSS on an N-frame M67 subset.")
    p0.add_argument("--frames", required=True, type=int)
    p0.set_defaults(func=cmd_step0)

    ps = sub.add_parser("sweep", parents=[shared],
                        help="B2: worker sweep over Steps 1-7 (NGC 6811).")
    ps.add_argument("--workers", default="1,2,4,8,12,16")
    ps.add_argument("--steps", default="1-7")
    ps.set_defaults(func=cmd_sweep)

    pr = sub.add_parser("repro", parents=[shared],
                        help="B4: repeat the same run; drift = noise floor.")
    pr.add_argument("--repeat", default=3, type=int)
    pr.add_argument("--workers", default=12, type=int)
    pr.add_argument("--steps", default="1-7")
    pr.set_defaults(func=cmd_repro)

    pp = sub.add_parser("parity", parents=[shared],
                        help="T0.5: parity gate (see apex.benchmark.parity).")
    pp.add_argument("parity_args", nargs=argparse.REMAINDER)
    pp.set_defaults(func=cmd_parity)
    return parser


def split_parity(argv: list) -> tuple:
    """Split ``argv`` at a ``parity`` subcommand: ``(head, tail)``.

    ``parity`` forwards a free-form argument list on to
    :mod:`apex.benchmark.parity`, and argparse's REMAINDER cannot absorb a
    leading option token (bpo-17050) — so the split has to happen before the
    parser sees it.  Anything that is a value of a global option (``--out-root
    parity``) must not be mistaken for the subcommand, hence the skip.

    Returns ``(None, None)`` when this is not a parity invocation.
    """
    takes_value = {"--out-root", "-o"}
    i = 0
    while i < len(argv):
        token = argv[i]
        if token in takes_value:
            i += 2
            continue
        if token.startswith("-"):
            i += 1
            continue
        return (argv[:i], argv[i + 1:]) if token == "parity" else (None, None)
    return (None, None)


def main(argv=None) -> int:
    argv = list(argv) if argv is not None else sys.argv[1:]
    head, tail = split_parity(argv)
    if tail is not None:
        from apex.benchmark import parity

        # Still parse the head so a misspelled global flag is an error rather
        # than being silently handed to the parity parser.
        build_parser().parse_args([*head, "parity"])
        return parity.main(tail)
    resources.warn_if_busy()
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
