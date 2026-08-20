"""Does the candidate table shrink when there are more nights?

D-014 closed on "do not pick a period; lay the candidates out". That only helps
if the table reflects the situation — a run with more baseline should offer
fewer candidates, because more baseline is exactly what kills the aliases. If it
always returns eight rows regardless, the table is decoration.

Controlled: one object, one reduction, one ensemble. The only thing that changes
is how much of the curve the step is allowed to see. YZ Boo was observed on two
nights, so this runs night 1 alone, night 2 alone, and both.
"""
import json
import logging
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).absolute().parent.parent))

import numpy as np
import pandas as pd

from apex.analysis.light_curve.target_config import LcTarget, write_selection
from apex.pipeline.base import StepStatus
from apex.pipeline.context import RunContext
from apex.pipeline.steps.lc_period import LcPeriodStep
from apex.utils.step_paths_lc import step9_lc_dir, step11_period_dir

SOURCE = Path(r"E:\APEX_validation\reprocess\YZBoo_2n\result")
TARGET_ID = 153
LIT = 0.104092          # YZ Boo, literature


def _params(result_dir):
    from types import SimpleNamespace

    return SimpleNamespace(P=SimpleNamespace(
        result_dir=str(result_dir), data_dir=str(result_dir),
        lc_target_id=TARGET_ID, lc_target_name="YZ Boo", lc_comparison_ids="",
        lc_comparison_mode="manual", lc_comparison_count=6, lc_filter="",
        lc_period_methods="ls,pdm",
        period_min_days=0.01, period_max_days=10.0,
    ))


def run_subset(curve: pd.DataFrame, label: str) -> dict:
    """Run Step 11 on this subset in a throwaway workspace."""
    workspace = Path(tempfile.mkdtemp(prefix="apex_nights_"))
    try:
        out = Path(step9_lc_dir(workspace))
        out.mkdir(parents=True, exist_ok=True)
        curve.to_csv(out / f"lightcurve_ID{TARGET_ID}_raw.csv", index=False)
        write_selection(workspace, LcTarget(target_id=TARGET_ID), [1, 2, 3])

        ctx = RunContext(mode="lc", params=_params(workspace),
                         result_dir=workspace, data_dir=workspace,
                         logger=logging.getLogger("nights"))
        result = LcPeriodStep().run(ctx)
        if result.status is not StepStatus.OK:
            return {"label": label, "error": result.message}

        table = step11_period_dir(workspace) / f"period_candidates_all_ID{TARGET_ID}.csv"
        analysis = step11_period_dir(workspace) / f"period_analysis_all_ID{TARGET_ID}.json"
        if not analysis.exists():                       # single band -> named file
            found = sorted(step11_period_dir(workspace).glob("period_analysis_*.json"))
            analysis = found[0] if found else None
        body = json.loads(analysis.read_text(encoding="utf-8")) if analysis else {}
        alias = body.get("alias_analysis") or {}
        rows = pd.read_csv(table) if table.exists() else pd.DataFrame()
        ranked = rows[rows["source"] == "alias candidate"] if not rows.empty else rows

        best, gap = None, None
        if not ranked.empty:
            ordered = ranked.sort_values("rank")
            best = float(ordered.iloc[0]["period_days"])
            # Distance to the runner-up, in delta BIC. This — not the number of
            # rows — is what says how sure the run is: more baseline resolves
            # MORE aliases, so the table gets longer as the answer gets better.
            if len(ordered) > 1:
                gap = float(ordered.iloc[1]["strength"]) - float(ordered.iloc[0]["strength"])
        return {
            "label": label,
            "n_points": int(len(curve)),
            "n_nights": int(alias.get("n_nights", 0)),
            "baseline_d": float(alias.get("baseline_days", float("nan"))),
            "status": str(alias.get("status", "")),
            "n_candidates": int(len(ranked)),
            "rank1_period": best,
            "rank1_err_pct": (best - LIT) / LIT * 100 if best else None,
            "bic_gap": gap,
        }
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def main() -> None:
    curve = pd.read_csv(SOURCE / f"lc_lightcurve/lightcurve_ID{TARGET_ID}_raw.csv")
    jd = curve["BJD_TDB"].to_numpy(float)
    # Split on the real gap, not on a stored night_id — the window's saved run
    # labelled both nights 0, so the column cannot be trusted as the boundary.
    order = np.argsort(jd)
    gaps = np.diff(jd[order])
    cut = jd[order][int(np.argmax(gaps))] + float(gaps.max()) / 2.0
    night1, night2 = curve[curve["BJD_TDB"] <= cut], curve[curve["BJD_TDB"] > cut]

    print(f"  YZ Boo · 총 {len(curve)} 점 · 실제 간극 {gaps.max():.3f} 일에서 분리")
    print(f"  밤 1: {len(night1)} 점 · 밤 2: {len(night2)} 점\n")
    print(f"  {'경우':<14}{'점':>5}{'밤':>4}{'기준선(일)':>11}{'판정':>13}"
          f"{'후보 수':>8}{'2위와 ΔBIC':>12}{'rank1 (일)':>12}{'문헌대비':>10}")
    print("  " + "-" * 92)
    for label, subset in (("밤 1 만", night1), ("밤 2 만", night2),
                          ("두 밤 전부", curve)):
        r = run_subset(subset, label)
        if "error" in r:
            print(f"  {label:<14} 실패: {r['error'][:60]}")
            continue
        err = f"{r['rank1_err_pct']:+.2f}%" if r["rank1_err_pct"] is not None else "—"
        per = f"{r['rank1_period']:.6f}" if r["rank1_period"] else "—"
        gap = f"{r['bic_gap']:.1f}" if r["bic_gap"] is not None else "—"
        print(f"  {label:<14}{r['n_points']:>5}{r['n_nights']:>4}"
              f"{r['baseline_d']:>11.3f}{r['status']:>13}"
              f"{r['n_candidates']:>8}{gap:>12}{per:>12}{err:>10}")


if __name__ == "__main__":
    main()
