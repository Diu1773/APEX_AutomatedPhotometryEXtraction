"""Does the qfit cut help the thing users actually look at?

The artificial-star benchmark says the cut removes measurements that are 2-5x
more wrong than the ones it keeps, on three instruments. That is a good reason
to believe the statistic carries information. It is *not* a reason to switch it
on by default, because the benchmark's metric is scatter against injected truth
and what a user reads is a zero point, a colour term and a CMD.

The threshold sweep made the question sharper rather than easier: scatter
improves smoothly all the way from qfit ≤ 10 down to qfit ≤ 1, with no knee.
There is no correct value waiting to be discovered — only a trade to be priced.
So price it on the product: run Step 10 on the same night with the cut at
several thresholds and read off what changes in the zero point.

Method: Step 10 reads Step 8's per-frame tables and keeps `flags_psf == 0`.
This script writes filtered copies of those tables into a scratch workspace,
one set per threshold, runs the headless Step 10 on each, and collects the
per-frame zero-point scatter and the number of calibration stars. Nothing in
the real workspace is modified.
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).absolute().parent
REPO = HERE.parents[1]
PYTHON = str(REPO / ".venv-deploy" / "Scripts" / "python.exe")


def filtered_workspace(source: Path, dest: Path, threshold: float) -> int:
    """Copy a workspace, dropping Step 8 rows above the qfit threshold."""
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True)
    for item in source.iterdir():
        if item.name == "result":
            continue
        (shutil.copytree if item.is_dir() else shutil.copy2)(item, dest / item.name)

    result_src, result_dst = source / "result", dest / "result"
    result_dst.mkdir()
    for item in result_src.iterdir():
        if item.name.startswith("cmd_zeropoint"):
            continue          # Step 10 rebuilds these
        (shutil.copytree if item.is_dir() else shutil.copy2)(item, result_dst / item.name)

    dropped = 0
    for tsv in (result_dst / "cmd_psf").glob("photometry_*.tsv"):
        table = pd.read_csv(tsv, sep="\t")
        if "qfit" not in table.columns:
            continue
        qfit = pd.to_numeric(table["qfit"], errors="coerce")
        keep = ~np.isfinite(qfit) | (qfit <= threshold)
        dropped += int((~keep).sum())
        table[keep].to_csv(tsv, sep="\t", index=False)

    # Repoint the copied workspace at itself. This is the whole safety of the
    # experiment: Step 10 takes its output directory from `io.result_dir` in the
    # config, so a copy that still names the source writes its results **into
    # the source**. That happened on the first attempt (2026-08-15) — the run
    # succeeded, reported "완료", and silently replaced the real M13 zero point.
    # It was recoverable only because Step 10 is deterministic and the Step 8
    # inputs were untouched. Hence: rewrite, then verify, then refuse.
    config = dest / "apex_config.json"
    if not config.exists():
        raise SystemExit(f"{dest}: apex_config.json 이 없다 — 경로 재지정 불가")
    data = json.loads(config.read_text(encoding="utf-8"))
    io = data.setdefault("io", {})
    for key in ("data_dir", "result_dir"):
        value = io.get(key)
        if isinstance(value, str) and value:
            io[key] = str(dest / Path(value).name)
    config.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")

    check = json.loads(config.read_text(encoding="utf-8")).get("io", {})
    for key in ("data_dir", "result_dir"):
        target = str(check.get(key, ""))
        if str(source) in target or str(dest) not in target:
            raise SystemExit(
                f"{dest}: io.{key} 가 여전히 원본을 가리킨다 ({target}) — 중단")
    return dropped


def read_zeropoint(work: Path) -> dict:
    """Per-band zero point health, from Step 10's own per-frame table.

    `zp_scatter` is the scatter of the calibration stars *within* one frame, so
    it reads the effect of the cut directly rather than through frame-to-frame
    variation. `n_ref` is what the cut costs in calibration stars.
    """
    path = work / "result" / "cmd_zeropoint" / "frame_zeropoint.csv"
    if not path.exists():
        return {}
    table = pd.read_csv(path)
    out: dict = {}
    for band, group in table.groupby("filter"):
        scatter = pd.to_numeric(group["zp_scatter"], errors="coerce").dropna()
        n_ref = pd.to_numeric(group["n_ref"], errors="coerce").dropna()
        outlier = pd.to_numeric(group.get("outlier_fraction"), errors="coerce").dropna()
        if len(scatter):
            out[f"zp_scatter_{band}_mmag"] = float(np.median(scatter)) * 1000
        if len(n_ref):
            out[f"n_ref_{band}"] = float(np.median(n_ref))
        if len(outlier):
            out[f"outlier_frac_{band}"] = float(np.median(outlier))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", default=r"E:\APEX_validation\phase3\M13")
    ap.add_argument("--scratch", default=r"E:\APEX_validation\psf_engines\qfit_ab")
    ap.add_argument("--thresholds", type=float, nargs="+",
                    default=[1e9, 6.0, 3.0, 2.0, 1.0])
    args = ap.parse_args()

    source = Path(args.source)
    scratch = Path(args.scratch)
    scratch.mkdir(parents=True, exist_ok=True)
    rows = []
    for threshold in args.thresholds:
        tag = "off" if threshold > 1e8 else f"{threshold:g}"
        work = scratch / f"qfit_{tag}"
        print(f"\n=== qfit <= {tag} ===", flush=True)
        dropped = filtered_workspace(source, work, threshold)
        print(f"  step8 행 {dropped}개 제거", flush=True)
        params = work / "apex_config.json"
        started = time.perf_counter()
        log = work / "step10.log"
        with log.open("w", encoding="utf-8", errors="replace") as handle:
            result = subprocess.run(
                [PYTHON, "-X", "utf8", str(REPO / "scripts" / "run_step10_headless.py"),
                 "--params", str(params), "--no-backup"],
                stdout=handle, stderr=subprocess.STDOUT, text=True)
        elapsed = time.perf_counter() - started
        ok = result.returncode == 0
        print(f"  step10 {'완료' if ok else '실패'} ({elapsed / 60:.1f}분)", flush=True)
        if not ok:
            print(log.read_text(encoding="utf-8", errors="replace")[-1200:], flush=True)
        row = {"threshold": tag, "rows_dropped": dropped, "ok": ok}
        row.update(read_zeropoint(work))
        rows.append(row)
        for key, value in sorted(row.items()):
            if key.startswith(("zp_scatter", "n_ref", "outlier_frac")):
                print(f"    {key} = {value:.1f}", flush=True)

    table = pd.DataFrame(rows)
    out = scratch / "qfit_gate_product_ab.csv"
    table.to_csv(out, index=False)
    print(f"\n{table.to_string(index=False)}")
    print(f"\nsaved -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
