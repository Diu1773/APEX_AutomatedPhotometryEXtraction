"""O2 A/B: one global worker count vs the measured per-stage ceilings.

The single-run comparison this replaces was confounded twice over.  Both
policies wrote into the *same* result directory, and — until the fix in
`apex/pipeline/context.py` — `--result-dir` did not move `cache_dir`, so the
detection and WCS caches were shared by every run regardless.  Whichever
policy ran second inherited the other's per-file cache and looked faster for
free.

So: a fresh result directory per run (which now carries its own cache), one
discarded warm-up run so nobody pays the cold read of the input frames, and
the two policies alternated across repeats so a slow disk or a background
process cannot land on one policy only.

Usage:  python -X utf8 benchmark/o2_ab.py [--repeats 2]
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from apex.benchmark import resources  # noqa: E402

CFG = Path(r"E:\APEX_validation\reprocess\NGC6811\apex_config_20260807.json")
SCI = Path(r"E:\APEX_validation\reprocess\NGC6811\sci")
SCRATCH = Path(r"E:\APEX_validation\bench")
STEPS = "1-7"

# "old" is the production default before O2: one count for the whole pipeline.
# "new" is get_parallel_workers(stage=...) with the measured ceilings, which is
# what an unset APEX_MAX_WORKERS now selects.
POLICIES = (("old_global12", "12"), ("new_stagecaps", None))


def run(label: str, workers: str | None, out: Path) -> dict:
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)

    env = {k: v for k, v in os.environ.items() if k != "APEX_MAX_WORKERS"}
    if workers:
        env["APEX_MAX_WORKERS"] = workers

    cmd = [sys.executable, "-X", "utf8", "-m", "apex.cli", "run",
           "--mode", "cmd", "--config", str(CFG), "--steps", STEPS,
           "--data-dir", str(SCI), "--result-dir", str(out), "--force"]
    rc, metrics = resources.measure_command(cmd, label=label, env=env, cwd=str(REPO))
    metrics["policy"] = label
    metrics["result_dir"] = str(out)

    manifest = out / "pipeline_run.json"
    if manifest.exists():
        steps = json.loads(manifest.read_text(encoding="utf-8"))["steps"]
        metrics["per_step"] = {s["key"]: s["duration_s"] for s in steps}
        metrics["all_ok"] = all(s["status"] == "ok" for s in steps)
    else:
        metrics["per_step"], metrics["all_ok"] = {}, False

    # The cache must have landed inside this run's own directory.
    metrics["own_cache"] = (out / "cache").exists()
    return metrics


def _line(m: dict, rep: int) -> str:
    ps = m["per_step"]
    return (f"{m['policy']:14s} #{rep}  wall={m['wall_s']:7.1f}s  "
            f"RSS={m['peak_rss_mb']:>6,.0f} USS={m['peak_uss_mb'] or 0:>6,.0f} MB  "
            f"detect={ps.get('detect', 0):5.0f} wcs={ps.get('wcs', 0):5.0f} "
            f"phot={ps.get('forcedphot', 0):6.0f}  "
            f"{'ok' if m['all_ok'] else 'FAILED'} cache={'own' if m['own_cache'] else 'SHARED!'}")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repeats", type=int, default=2)
    ap.add_argument("--out", default=str(REPO / "benchmark" / "perf" / "20260808"
                                         / "o2_ab_NGC6811.json"))
    args = ap.parse_args()

    print("warm-up (discarded) …", flush=True)
    run("warmup", None, SCRATCH / "ab_warm")

    runs = []
    for rep in range(args.repeats):
        order = POLICIES if rep % 2 == 0 else POLICIES[::-1]
        for label, workers in order:
            m = run(label, workers, SCRATCH / f"ab_{label}_{rep}")
            m["repeat"] = rep
            runs.append(m)
            print(_line(m, rep), flush=True)

    print()
    summary = {}
    for label, _ in POLICIES:
        walls = [m["wall_s"] for m in runs if m["policy"] == label]
        summary[label] = {
            "n": len(walls),
            "wall_s_median": round(statistics.median(walls), 1),
            "wall_s_min": round(min(walls), 1),
            "wall_s_max": round(max(walls), 1),
            "peak_uss_mb_max": max((m["peak_uss_mb"] or 0)
                                   for m in runs if m["policy"] == label),
        }
        s = summary[label]
        print(f"{label:14s} median {s['wall_s_median']:7.1f}s  "
              f"[{s['wall_s_min']:.1f}, {s['wall_s_max']:.1f}]  "
              f"USS max {s['peak_uss_mb_max']:,.0f} MB")

    old, new = (summary[p]["wall_s_median"] for p, _ in POLICIES)
    print(f"\nspeed-up (median): {old / new:.2f}x")

    resources.save_metrics(
        {"label": "o2_ab_NGC6811", "steps": STEPS, "config": str(CFG),
         "summary": summary, "runs": runs}, Path(args.out))
    print(f"saved -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
