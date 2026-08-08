"""O2, measured properly: three worker policies, plus a clean detect/wcs sweep.

The first A/B could not see what it claimed to. The NGC 6811 workspace config
carries ``[parallel] max_workers = 1``, and a stage ceiling is a *maximum* —
``min(1, 4) = 1`` — so the "per-stage caps" arm ran everything serially. What
it actually compared was global-1 against global-12.

Three arms, then:

* ``global12``  — ``APEX_MAX_WORKERS=12``; the env knob bypasses the ceilings
  by design, so this is the pre-O2 production behaviour.
* ``stagecaps`` — a config copy with ``max_workers = 12`` and no env override,
  which is the only way the ceilings actually bind: detect 4, wcs 4,
  forcedphot 1.
* ``serial``    — the workspace config as-is (``max_workers = 1``).

``sweep`` re-measures detect and wcs alone (Steps 1–5) across worker counts.
B2's wcs row was invalidated by the shared-cache bug (see
benchmark/perf/20260807/RESULTS.md); running only the cheap half of the
pipeline buys that curve back for ~25 min instead of ~2 h.

Usage:
    python -X utf8 benchmark/o2_arms.py arms  [--repeats 2]
    python -X utf8 benchmark/o2_arms.py sweep [--workers 1,2,4,8,12,16]
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
OUT_DIR = REPO / "benchmark" / "perf" / "20260808"


def config_with_max_workers(value: int) -> Path:
    """A copy of the workspace config whose only change is max_workers."""
    data = json.loads(CFG.read_text(encoding="utf-8"))
    data.setdefault("parallel", {})["max_workers"] = int(value)
    out = SCRATCH / f"cfg_maxworkers{value}.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(data, ensure_ascii=False, indent=1), encoding="utf-8")
    return out


def run(label: str, *, config: Path, env_workers: str | None, steps: str,
        out: Path) -> dict:
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)

    env = {k: v for k, v in os.environ.items() if k != "APEX_MAX_WORKERS"}
    if env_workers:
        env["APEX_MAX_WORKERS"] = env_workers

    cmd = [sys.executable, "-X", "utf8", "-m", "apex.cli", "run",
           "--mode", "cmd", "--config", str(config), "--steps", steps,
           "--data-dir", str(SCI), "--result-dir", str(out), "--force"]
    rc, m = resources.measure_command(cmd, label=label, env=env, cwd=str(REPO))
    m["policy"] = label

    manifest = out / "pipeline_run.json"
    if manifest.exists():
        steps_done = json.loads(manifest.read_text(encoding="utf-8"))["steps"]
        m["per_step"] = {s["key"]: s["duration_s"] for s in steps_done}
        m["all_ok"] = all(s["status"] == "ok" for s in steps_done)
        m["messages"] = {s["key"]: s.get("message", "") for s in steps_done}
    else:
        m["per_step"], m["all_ok"], m["messages"] = {}, False, {}
    m["own_cache"] = (out / "cache").exists()
    return m


def checkpoint(path: Path, payload: dict) -> None:
    """Write the batch after every run, not only at the end.

    A scheduled shutdown killed a sweep at 05:20 that had already finished
    three of six worker counts; the JSON was written last, so those results
    only survived because the console log happened to be on disk.
    """
    resources.save_metrics(payload, path)


def _row(m: dict) -> str:
    p = m["per_step"]
    return (f"{m['policy']:12s} wall={m['wall_s']:7.1f}s  "
            f"detect={p.get('detect', 0):6.1f} wcs={p.get('wcs', 0):6.1f} "
            f"phot={p.get('forcedphot', 0):7.1f}  "
            f"USS={m['peak_uss_mb'] or 0:>6,.0f} MB  "
            f"{'ok' if m['all_ok'] else 'FAILED'}")


def cmd_arms(args) -> int:
    arms = (
        ("global12", CFG, "12"),
        ("stagecaps", config_with_max_workers(12), None),
        ("serial", CFG, None),          # config already pins max_workers = 1
    )
    print("warm-up (discarded) …", flush=True)
    run("warmup", config=CFG, env_workers=None, steps="1-7",
        out=SCRATCH / "arm_warm")

    runs = []
    for rep in range(args.repeats):
        order = arms if rep % 2 == 0 else arms[::-1]
        for label, cfg, env_w in order:
            m = run(label, config=cfg, env_workers=env_w, steps="1-7",
                    out=SCRATCH / f"arm_{label}_{rep}")
            m["repeat"] = rep
            runs.append(m)
            print(f"#{rep} " + _row(m), flush=True)
            checkpoint(OUT_DIR / "o2_arms_NGC6811.json",
                       {"label": "o2_arms_NGC6811", "complete": False,
                        "runs": runs})

    print()
    summary = {}
    for label, _, _ in arms:
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
        print(f"{label:12s} median {s['wall_s_median']:7.1f}s "
              f"[{s['wall_s_min']:.1f}, {s['wall_s_max']:.1f}]  "
              f"USS max {s['peak_uss_mb_max']:,.0f} MB")

    resources.save_metrics({"label": "o2_arms_NGC6811", "complete": True,
                            "summary": summary, "runs": runs},
                           OUT_DIR / "o2_arms_NGC6811.json")
    print(f"saved -> {OUT_DIR / 'o2_arms_NGC6811.json'}")
    return 0


def cmd_sweep(args) -> int:
    """B2′: detect and wcs across worker counts, each with its own cache.

    Repeats matter here more than extra worker counts: a single pass put wcs
    between 99.6 s and 172.5 s with no trend, which cannot distinguish "no
    parallel gain" from "noisy".  The order is reversed on odd repeats so a
    drifting machine cannot favour one end of the range.
    """
    counts = [int(w) for w in args.workers.split(",")]
    out_path = OUT_DIR / "b2prime_detect_wcs.json"
    print("warm-up (discarded) …", flush=True)
    run("warmup", config=CFG, env_workers="4", steps="1-5",
        out=SCRATCH / "sweep2_warm")

    runs = []
    for rep in range(args.repeats):
        for w in (counts if rep % 2 == 0 else counts[::-1]):
            m = run(f"w{w:02d}", config=CFG, env_workers=str(w), steps="1-5",
                    out=SCRATCH / f"sweep2_w{w:02d}")
            m["workers"], m["repeat"] = w, rep
            runs.append(m)
            p = m["per_step"]
            print(f"#{rep} w={w:<3d} detect={p.get('detect', 0):6.1f}s  "
                  f"wcs={p.get('wcs', 0):6.1f}s  wall={m['wall_s']:7.1f}s  "
                  f"USS={m['peak_uss_mb'] or 0:>6,.0f} MB  "
                  f"{'ok' if m['all_ok'] else 'FAILED'} "
                  f"cache={'own' if m['own_cache'] else 'SHARED!'}", flush=True)
            checkpoint(out_path, {"label": "b2prime_detect_wcs_NGC6811",
                                  "steps": "1-5", "complete": False,
                                  "runs": runs})

    print()
    print(f"{'w':>3} {'detect med':>11} {'range':>15} {'wcs med':>9} {'range':>15}")
    summary = {}
    for w in counts:
        det = sorted(m["per_step"].get("detect", 0)
                     for m in runs if m["workers"] == w)
        wcs = sorted(m["per_step"].get("wcs", 0)
                     for m in runs if m["workers"] == w)
        summary[w] = {"n": len(det),
                      "detect_median": round(statistics.median(det), 1),
                      "detect_min": round(det[0], 1), "detect_max": round(det[-1], 1),
                      "wcs_median": round(statistics.median(wcs), 1),
                      "wcs_min": round(wcs[0], 1), "wcs_max": round(wcs[-1], 1)}
        s = summary[w]
        det_range = "[{:.1f}, {:.1f}]".format(s["detect_min"], s["detect_max"])
        wcs_range = "[{:.1f}, {:.1f}]".format(s["wcs_min"], s["wcs_max"])
        print(f"{w:>3} {s['detect_median']:>11.1f} {det_range:>15} "
              f"{s['wcs_median']:>9.1f} {wcs_range:>15}")

    resources.save_metrics({"label": "b2prime_detect_wcs_NGC6811",
                            "steps": "1-5", "complete": True,
                            "summary": {str(k): v for k, v in summary.items()},
                            "runs": runs}, out_path)
    print(f"saved -> {out_path}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("arms")
    a.add_argument("--repeats", type=int, default=2)
    a.set_defaults(func=cmd_arms)
    s = sub.add_parser("sweep")
    s.add_argument("--workers", default="1,2,4,8,12,16")
    s.add_argument("--repeats", type=int, default=1)
    s.set_defaults(func=cmd_sweep)
    resources.warn_if_busy()
    args = ap.parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
