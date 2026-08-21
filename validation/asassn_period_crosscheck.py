"""Does APEX's period analysis recover published periods on archive data?

The photometry is already checked — against IRAF, Gaia, PS1 and AIPPI. What is
not is the period stage, and the reason is data: YZ Boo has two nights, the
longest covering 1.23 cycles, so its aliases genuinely cannot be resolved. That
is a property of the observing run, not of the code, and no amount of local
reprocessing changes it.

Archive tables skip the parts already verified — detection, WCS, ID matching,
forced photometry — and hand the period stage a long, densely sampled series
with a published period beside it. ASAS-SN's variable star database is the
one used here: it serves per-star V-band tables and carries `period` and
`variable_type` for each entry, so every run has something to be wrong against.

Detached eclipsing binaries (`EA`) are the targets. They have a deep primary and
a shallower secondary minimum, so folding at half the true period puts both on
top of each other and fits nearly as well. That is where a period search is
most likely to be confidently wrong, and where APEX's candidate table has to
show both readings rather than pick one. Verified so far only on a synthetic
signal (`tests/test_lc_period_step.py`).

Cross-check only — nothing here feeds a figure or the manuscript. ASAS-SN's own
`period` is the comparison, and it is a catalogue value, not ground truth: it
carries its own alias risk, which is part of what the comparison shows.
"""
from __future__ import annotations

import csv
import io
import json
import logging
import re
import shutil
import sys
import tempfile
import urllib.request
from pathlib import Path

sys.path.insert(0, str(Path(__file__).absolute().parent.parent))

import numpy as np
import pandas as pd

from apex.analysis.light_curve.target_config import LcTarget, write_selection
from apex.pipeline.base import StepStatus
from apex.pipeline.context import RunContext
from apex.pipeline.steps.lc_period import LcPeriodStep
from apex.utils.step_paths_lc import step9_lc_dir, step11_period_dir

BASE = "https://asas-sn.osu.edu"
CATALOG = f"{BASE}/variables.csv?action=index&controller=variables"
UA = {"User-Agent": "Mozilla/5.0 (APEX validation; period cross-check)"}
CACHE = Path(__file__).absolute().parent / "asassn_cache"
TARGET_ID = 1


def _get(url: str, timeout: int = 90) -> bytes:
    return urllib.request.urlopen(
        urllib.request.Request(url, headers=UA), timeout=timeout).read()


def find_targets(want_type: str = "EA", n: int = 5, max_rows: int = 60000,
                 min_prob: float = 0.98, max_vmag: float = 13.0,
                 period_range: tuple[float, float] = (0.3, 5.0),
                 min_amplitude: float = 0.3) -> list[dict]:
    """Catalogue entries worth testing, asked for by type.

    The endpoint takes `variable_type` and honours it — a thousand rows of
    exactly that class, verified for RRAB, DSCT, HADS and EW. The first version
    streamed the unfiltered table and stopped after a fixed number of rows, so
    the rarer classes simply never appeared: RRAB and HADS came back empty from
    a scan that had 54 and 8 of them further down. A sample of zero is not a
    result about the rule, it is a result about the scan.
    """
    found, scanned = [], 0
    req = urllib.request.Request(f"{CATALOG}&variable_type={want_type}",
                                 headers=UA)
    with urllib.request.urlopen(req, timeout=120) as resp:
        reader = csv.DictReader(io.TextIOWrapper(resp, encoding="utf-8",
                                                 errors="replace"))
        for row in reader:
            scanned += 1
            if scanned > max_rows or len(found) >= n:
                break
            # Trust but check: a filter that silently stopped filtering would
            # quietly fill a pulsator sample with eclipsing binaries.
            if (row.get("variable_type") or "").strip() != want_type:
                continue
            try:
                prob = float(row.get("class_probability") or 0)
                vmag = float(row.get("mean_vmag") or 99)
                period = float(row.get("period") or 0)
                amp = float(row.get("amplitude") or 0)
            except ValueError:
                continue
            if (prob >= min_prob and vmag <= max_vmag and amp >= min_amplitude
                    and period_range[0] <= period <= period_range[1]):
                found.append({
                    "source_id": str(row["source_id"]),
                    "name": str(row["asassn_name"]),
                    "catalog_period": period,
                    "vmag": vmag,
                    "amplitude": amp,
                    "variable_type": want_type,
                })
    print(f"  카탈로그 {scanned:,} 행에서 {want_type} {len(found)} 개")
    return found


def fetch_light_curve(source_id: str) -> pd.DataFrame:
    """The star's V-band table, cached so a rerun does not re-fetch.

    The obvious-looking `/variables/<uuid>.csv` returns a header and no rows for
    every star tried; the page's own "Export V-band Data" link is the one that
    carries photometry.
    """
    CACHE.mkdir(parents=True, exist_ok=True)
    cached = CACHE / f"asassn_{source_id}.csv"
    if not cached.exists():
        page = _get(f"{BASE}/variables/{source_id}", timeout=60).decode("utf-8", "replace")
        m = re.search(r'href="(/variables/[0-9a-f-]{36})/star_data/export', page)
        if not m:
            return pd.DataFrame()
        cached.write_bytes(_get(f"{BASE}{m.group(1)}/star_data/export?type=csv"))
    frame = pd.read_csv(cached)
    frame.columns = [c.strip().lower().replace(" ", "_") for c in frame.columns]
    return frame


def to_apex_curve(frame: pd.DataFrame) -> pd.DataFrame:
    """The archive table in the shape Step 9 would have written.

    `diff_mag_raw` is the star's own magnitude with its median removed: ASAS-SN
    publishes calibrated magnitudes, so there is no ensemble to subtract and
    pretending otherwise would invent a comparison that was never measured.
    """
    out = pd.DataFrame({
        "file": [f"asassn_{i:05d}" for i in range(len(frame))],
        "filter": "V",
        "BJD_TDB": pd.to_numeric(frame["hjd"], errors="coerce"),
        "mag": pd.to_numeric(frame["mag"], errors="coerce"),
        "mag_err": pd.to_numeric(frame.get("mag_err"), errors="coerce"),
    })
    out = out[np.isfinite(out["BJD_TDB"]) & np.isfinite(out["mag"])].copy()
    out["mag_err"] = out["mag_err"].fillna(0.02)
    out["diff_mag_raw"] = out["mag"] - out["mag"].median()
    out["diff_err"] = out["mag_err"]
    return out.sort_values("BJD_TDB").reset_index(drop=True)


def _params(result_dir, period_max: float):
    from types import SimpleNamespace

    return SimpleNamespace(P=SimpleNamespace(
        result_dir=str(result_dir), data_dir=str(result_dir),
        lc_target_id=TARGET_ID, lc_target_name="", lc_comparison_ids="",
        lc_comparison_mode="manual", lc_comparison_count=3, lc_filter="",
        lc_period_methods="ls,pdm",
        period_min_days=0.05, period_max_days=float(period_max),
    ))


def run_period(curve: pd.DataFrame, period_max: float) -> dict:
    """Step 11 on this curve, in a throwaway workspace."""
    workspace = Path(tempfile.mkdtemp(prefix="apex_asassn_"))
    try:
        out = Path(step9_lc_dir(workspace))
        out.mkdir(parents=True, exist_ok=True)
        curve.to_csv(out / f"lightcurve_ID{TARGET_ID}_raw.csv", index=False)
        write_selection(workspace, LcTarget(target_id=TARGET_ID), [2, 3, 4])

        ctx = RunContext(mode="lc", params=_params(workspace, period_max),
                         result_dir=workspace, data_dir=workspace,
                         logger=logging.getLogger("asassn"))
        result = LcPeriodStep().run(ctx)
        if result.status is not StepStatus.OK:
            return {"error": result.message}

        folder = Path(step11_period_dir(workspace))
        analysis = next(iter(sorted(folder.glob("period_analysis_*.json"))), None)
        body = json.loads(analysis.read_text(encoding="utf-8")) if analysis else {}
        alias = body.get("alias_analysis") or {}
        table = next(iter(sorted(folder.glob("period_candidates_*.csv"))), None)
        rows = pd.read_csv(table) if table is not None else pd.DataFrame()
        return {
            "results": {k: (v or {}).get("best_period")
                        for k, v in (body.get("results") or {}).items()},
            "status": str(alias.get("status", "")),
            "adopted": alias.get("adopted_period"),
            "n_nights": int(alias.get("n_nights", 0)),
            "baseline_days": float(alias.get("baseline_days", float("nan"))),
            "candidates": rows,
        }
    finally:
        shutil.rmtree(workspace, ignore_errors=True)


def _closest(candidates: pd.DataFrame, truth: float) -> tuple[float, int]:
    """The candidate nearest the published period, and where it ranks."""
    if candidates.empty or "period_days" not in candidates.columns:
        return float("nan"), -1
    ranked = candidates[candidates["source"] == "alias candidate"].copy()
    if ranked.empty:
        ranked = candidates.copy()
    idx = (ranked["period_days"] - truth).abs().idxmin()
    row = ranked.loc[idx]
    rank = int(row["rank"]) if str(row.get("rank", "")).strip() not in ("", "nan") else -1
    return float(row["period_days"]), rank


def main() -> None:
    targets = find_targets()
    if not targets:
        print("  후보 없음")
        return
    print(f"\n  {'이름':<32}{'점':>5}{'밤':>4}{'기준선':>8}{'카탈로그 P':>11}"
          f"{'APEX 채택':>11}{'오차%':>8}{'판정':>13}{'카탈로그P 표에':>12}")
    print("  " + "-" * 102)

    for t in targets:
        frame = fetch_light_curve(t["source_id"])
        if frame.empty or "hjd" not in frame.columns:
            print(f"  {t['name']:<32} 광곡선 없음")
            continue
        curve = to_apex_curve(frame)
        truth = t["catalog_period"]
        out = run_period(curve, period_max=min(10.0, truth * 4))
        if "error" in out:
            print(f"  {t['name']:<32} 실패: {out['error'][:44]}")
            continue

        adopted = out["adopted"]
        near, rank = _closest(out["candidates"], truth)
        # The question is whether the *catalogue* period is on the table at all.
        # Asking whether P/2 is there answers nothing: when the run adopts P/2
        # it is trivially there, which is what the first version of this column
        # measured and why every row said "yes".
        periods = (out["candidates"]["period_days"].to_numpy(float)
                   if not out["candidates"].empty else np.array([]))
        has_half = bool(periods.size and
                        (np.abs(periods - truth) / truth < 0.03).any())
        err = (adopted - truth) / truth * 100 if adopted else float("nan")
        print(f"  {t['name']:<32}{len(curve):>5}{out['n_nights']:>4}"
              f"{out['baseline_days']:>8.0f}{truth:>11.5f}"
              f"{(adopted or float('nan')):>11.5f}{err:>8.2f}"
              f"{out['status']:>13}{('있음' if has_half else '없음'):>12}")
        if rank > 0:
            print(f"      └ 카탈로그 주기에 가장 가까운 후보: {near:.5f} d "
                  f"({(near - truth) / truth * 100:+.2f} %), 표에서 {rank} 위")

    print(f"\n  광곡선 캐시: {CACHE}")
    print("  ASAS-SN 카탈로그 주기는 참값이 아니라 대조값이다 — 그쪽도 별칭을 탈 수 있다.")


if __name__ == "__main__":
    main()


# ── the broad run: several kinds, several of each ──────────────────────────
#
# One kind cannot answer the question. Eclipsing binaries *should* be doubled
# relative to the periodogram peak; pulsators should not, and a rule that
# doubles them is worse than no rule. Both directions are checked here.

KINDS: tuple[tuple[str, str, bool], ...] = (
    ("EA", "분리형 식쌍성", True),
    ("EB", "베타 라이라형 식쌍성", True),
    ("EW", "접촉형 식쌍성", True),
    ("RRAB", "RR Lyrae 기본진동", False),
    ("RRC", "RR Lyrae 1차 배진동", False),
    ("HADS", "고진폭 방패자리 델타", False),
)
# DSCT (ordinary delta Scuti) is absent on purpose: at V<=13.5 and amplitude
# >=0.15 the catalogue has none in a thousand rows — they are faint and shallow.
# HADS is the same pulsation with a large enough amplitude to measure, so the
# pulsator side is represented without loosening the cuts for one class only.
"""Type, a name to print, and whether its period should be twice the peak.

Eclipsing types dip twice per orbit so their periodogram peak sits at half the
catalogue period. Pulsators brighten and fade once per cycle, so the peak *is*
the period and doubling it would be an error.
"""


def survey(per_kind: int = 8, max_rows: int = 60000) -> "pd.DataFrame":
    """Run every kind through the period step and score both directions."""
    from apex.analysis.light_curve.period_harmonic_service import resolve_harmonic
    from astropy.timeseries import LombScargle

    rows = []
    for kind, label, should_double in KINDS:
        found = find_targets(want_type=kind, n=per_kind, max_rows=max_rows,
                             min_prob=0.95, max_vmag=13.5,
                             period_range=(0.15, 8.0), min_amplitude=0.15)
        for target in found:
            frame = fetch_light_curve(target["source_id"])
            if frame.empty or "hjd" not in frame.columns:
                continue
            curve = to_apex_curve(frame)
            if len(curve) < 60:
                continue
            t = curve["BJD_TDB"].to_numpy(float)
            y = curve["diff_mag_raw"].to_numpy(float)
            e = curve["diff_err"].to_numpy(float)
            truth = float(target["catalog_period"])
            pmax = min(12.0, truth * 4)
            freq, power = LombScargle(t, y, e).autopower(
                minimum_frequency=1.0 / pmax, maximum_frequency=1.0 / 0.05,
                samples_per_peak=12)
            peak = float(1.0 / freq[int(np.argmax(power))])
            verdict = resolve_harmonic(t, y, peak, mag_err=e,
                                       min_period=0.05, max_period=pmax)
            rows.append({
                "kind": kind, "label": label, "should_double": should_double,
                "name": target["name"], "n": len(curve),
                "catalog_period": truth, "peak_period": peak,
                "adopted": verdict.adopted_period, "factor": verdict.factor,
                "peak_err_pct": (peak - truth) / truth * 100,
                "adopted_err_pct": (verdict.adopted_period - truth) / truth * 100,
                "coherence": next((r.get("subcycle_coherence") for r in verdict.candidates
                                   if abs(r["factor"] - 2) < 0.01), float("nan")),
                "subcycle_sigma": next((r.get("subcycle_sigma") for r in verdict.candidates
                                        if abs(r["factor"] - 2) < 0.01), float("nan")),
            })
    return pd.DataFrame(rows)


def report_survey(table: "pd.DataFrame") -> None:
    if table.empty:
        print("  표본 없음")
        return
    near = lambda col: table[col].abs() < 1.0
    print(f"\n  {'부류':<22}{'표본':>5}{'주기도표 정답':>13}{'배음 후 정답':>13}"
          f"{'배 늘림':>9}{'맥동성 오탐':>12}{'아직 짧음':>10}")
    print("  " + "-" * 78)
    for kind, label, should_double in KINDS:
        part = table[table["kind"] == kind]
        if part.empty:
            continue
        doubled = (part["factor"] > 1.5).sum()
        # Two different errors, and lumping them hid one. A pulsator that got
        # doubled is the rule misfiring; an eclipsing binary that is still short
        # after doubling was simply not doubled enough (its periodogram peak sat
        # at a quarter of the period).
        if should_double:
            wrong = int(((part["factor"] > 1.5)
                         & (part["adopted_err_pct"].abs() >= 1.0)
                         & (part["adopted_err_pct"] < -1.0)).sum())
        else:
            wrong = int((part["factor"] > 1.5).sum())
        short = int(((part["factor"] > 1.5) & (part["adopted_err_pct"] < -1.0)).sum())
        print(f"  {kind + ' ' + label:<22}{len(part):>5}"
              f"{int((part['peak_err_pct'].abs() < 1.0).sum()):>13}"
              f"{int((part['adopted_err_pct'].abs() < 1.0).sum()):>13}"
              f"{int(doubled):>9}{(wrong if not should_double else 0):>12}{short:>10}")
    ecl = table[table["should_double"]]
    pul = table[~table["should_double"]]
    print("  " + "-" * 78)
    print(f"  식쌍성  {len(ecl):>3} 개 중 카탈로그 1% 이내: "
          f"주기도표 {int((ecl['peak_err_pct'].abs() < 1.0).sum())} → "
          f"배음 후 {int((ecl['adopted_err_pct'].abs() < 1.0).sum())}")
    print(f"  맥동성  {len(pul):>3} 개 중 카탈로그 1% 이내: "
          f"주기도표 {int((pul['peak_err_pct'].abs() < 1.0).sum())} → "
          f"배음 후 {int((pul['adopted_err_pct'].abs() < 1.0).sum())}"
          f"   (잘못 늘린 것 {int((pul['factor'] > 1.5).sum())} 개)")
