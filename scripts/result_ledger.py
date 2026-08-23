"""What is known about each result directory, and what is not.

Twelve directories hold the numbers behind the paper. Their names — `result`,
`result_psf`, `result_nocr`, `result_pre20260807` — are the only statement
anywhere of how they differ, and the names do not agree between clusters. On
2026-08-23 a check found that ten of the twelve carried no record of the
parameters that made them, and that the two which did had been written that
morning by accident.

This walks a root and says, per directory, what can be established from files
that are actually there:

  products      which steps left output, and whether PSF photometry ran
  frames        how many the run saw, and how many are still on disk
  history       runs in `apex_journal.jsonl` (appended, complete)
  manifest      steps in `pipeline_run.json` (rewritten, last run only)
  parameters    whether `parameters_used.json` exists, and whether the config
                beside it still hashes to what was recorded
  stray         records no code maintains any more

Nothing is inferred. A directory with no journal and no parameter record is
printed as `미상` — unknown — because that is the true state, and filling it in
from the config file sitting there now would assert something no evidence
supports. `apex_config.json` is edited; the products are not re-made.

Read-only. Run it as:

    .venv-deploy/Scripts/python.exe -X utf8 scripts/result_ledger.py
    .venv-deploy/Scripts/python.exe -X utf8 scripts/result_ledger.py --detail M13
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from apex.utils import run_journal as journal  # noqa: E402

DEFAULT_ROOT = Path("E:/APEX_validation/reprocess")

# A directory is a result directory if it holds any of these. `step1_*` alone is
# enough: a run that got no further still produced something worth explaining.
MARKERS = ("step1_file_selection", "pipeline_run.json", journal.JOURNAL_NAME)

# Product directories, in pipeline order, with the short label used in the table.
PRODUCTS = [
    ("step1_file_selection", "scan"),
    ("step2_crop", "crop"),
    ("step4_detection", "detect"),
    ("step5_wcs", "wcs"),
    ("step6_refbuild", "refbuild"),
    ("step7_forced_phot", "aperture"),
    ("cmd_psf", "PSF"),
    ("cmd_zeropoint", "zeropoint"),
    ("cmd_isochrone", "isochrone"),
    ("cmd_plot", "cmd-plot"),
    ("lc_selection", "lc-select"),
    ("lc_lightcurve", "lc-curve"),
    ("lc_period", "lc-period"),
]

# Written by one-off scripts that no longer exist in the tree. Their contents
# can contradict the directory they sit in — M67/result_psf's names `result` as
# its result_dir and a config format removed in August — so they are listed as
# strays rather than read as provenance.
STRAY_RECORDS = ("real_gui_run.json",)


def sha256_of(path: Path) -> str | None:
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return None


def count_frames(directory: Path | None) -> int | None:
    if directory is None or not directory.is_dir():
        return None
    return sum(1 for p in directory.iterdir()
               if p.suffix.lower() in (".fit", ".fits", ".fts"))


def read_headers(result_dir: Path) -> tuple[int | None, list[str]]:
    """How many frames the run saw, and in which filters.

    From the run's own `headers.csv` rather than from the data directory, so it
    still answers after the science frames have been cleaned up.
    """
    path = result_dir / "step1_file_selection" / "headers.csv"
    if not path.exists():
        return None, []
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
    except OSError:
        return None, []
    field = next((f for f in ("filter", "FILTER", "filter_key") if rows and f in rows[0]), None)
    filters = sorted({str(r[field]).strip() for r in rows if r.get(field)}) if field else []
    return len(rows), filters


def describe(result_dir: Path) -> dict:
    """Everything establishable about one directory."""
    out: dict = {"path": result_dir, "name": result_dir.name,
                 "workspace": result_dir.parent.name}

    out["products"] = [label for sub, label in PRODUCTS if (result_dir / sub).is_dir()]
    out["frames_seen"], out["filters"] = read_headers(result_dir)

    runs = journal.history(result_dir)
    out["runs"] = runs
    out["journal_runs"] = len(runs)
    out["journal_steps"] = sorted(journal.latest_steps(result_dir))
    out["notes"] = [n["text"] for r in runs for n in r["notes"]]

    manifest = result_dir / "pipeline_run.json"
    out["manifest_steps"] = []
    out["manifest_started"] = None
    out["manifest_packages"] = 0
    if manifest.exists():
        try:
            data = json.loads(manifest.read_text(encoding="utf-8"))
            out["manifest_steps"] = [s.get("index") for s in data.get("steps") or []]
            out["manifest_started"] = data.get("started")
            out["manifest_packages"] = len((data.get("environment") or {}).get("packages") or {})
        except (OSError, json.JSONDecodeError):
            out["manifest_steps"] = ["<읽을 수 없음>"]

    record = result_dir / "parameters_used.json"
    out["parameters"] = None
    out["config_match"] = None
    if record.exists():
        try:
            data = json.loads(record.read_text(encoding="utf-8"))
            out["parameters"] = len(data.get("settings") or {})
            cfg = (data.get("config") or {})
            recorded, path = cfg.get("sha256"), cfg.get("path")
            if recorded and path and Path(path).exists():
                out["config_match"] = (sha256_of(Path(path)) == recorded)
        except (OSError, json.JSONDecodeError):
            out["parameters"] = -1

    # Inputs. The config beside the workspace is the only pointer to them, and
    # it may have been edited since — so this reports what is there now and
    # never claims it is what the run used.
    out["data_dir"] = None
    out["frames_now"] = None
    config = result_dir.parent / "apex_config.json"
    if config.exists():
        try:
            io = json.loads(config.read_text(encoding="utf-8")).get("io") or {}
            if io.get("data_dir"):
                out["data_dir"] = Path(io["data_dir"])
                out["frames_now"] = count_frames(out["data_dir"])
        except (OSError, json.JSONDecodeError):
            pass

    out["strays"] = [n for n in STRAY_RECORDS if (result_dir / n).exists()]
    return out


def provenance_verdict(entry: dict) -> str:
    """One phrase for "can this directory explain itself"."""
    if entry["journal_runs"]:
        return "저널"
    if entry["parameters"]:
        return "파라미터만"
    if entry["manifest_steps"]:
        return "스텝목록만"
    return "미상"


def find_result_dirs(root: Path) -> list[Path]:
    found = []
    for workspace in sorted(p for p in root.iterdir() if p.is_dir()):
        for child in sorted(p for p in workspace.iterdir() if p.is_dir()):
            if any((child / m).exists() for m in MARKERS):
                found.append(child)
    return found


def print_table(entries: list[dict]) -> None:
    head = (f"  {'워크스페이스':<10} {'폴더':<20} {'측광':<9} {'필터':<9} "
            f"{'프레임':>10} {'출처':<10} {'설정대조':<9} 산출물")
    print(head)
    print("  " + "─" * (len(head) + 18))
    for e in entries:
        phot = "PSF" if "PSF" in e["products"] else ("구경" if "aperture" in e["products"] else "—")
        filters = ",".join(e["filters"]) or "—"
        seen, now = e["frames_seen"], e["frames_now"]
        frames = f"{seen if seen is not None else '?'}"
        if now is not None and seen is not None and now != seen:
            frames += f"→{now}"
        match = {True: "일치", False: "달라짐", None: "—"}[e["config_match"]]
        extra = [p for p in e["products"]
                 if p not in ("scan", "crop", "detect", "wcs", "refbuild", "aperture", "PSF")]
        print(f"  {e['workspace']:<10} {e['name']:<20} {phot:<9} {filters:<9} "
              f"{frames:>10} {provenance_verdict(e):<10} {match:<9} {','.join(extra) or '—'}")


def print_detail(entry: dict) -> None:
    print(f"\n  ── {entry['workspace']} / {entry['name']} " + "─" * 40)
    print(f"     경로       {entry['path']}")
    print(f"     산출물     {', '.join(entry['products']) or '없음'}")
    seen = entry["frames_seen"]
    print(f"     프레임     실행이 본 것 {seen if seen is not None else '기록 없음'}"
          f" · 지금 data_dir 에 {entry['frames_now'] if entry['frames_now'] is not None else '?'}")
    if entry["data_dir"]:
        print(f"                {entry['data_dir']}")

    if entry["runs"]:
        print(f"     저널       실행 {len(entry['runs'])}회")
        for run in entry["runs"]:
            steps = ",".join(str(s.get("index")) for s in run["steps"])
            print(f"       {run['started']}  {run['mode'] or '?':<4} "
                  f"스텝 [{steps}] 성공={run['success']}")
            cfg = run["config"] or {}
            if cfg.get("path"):
                print(f"         설정 {cfg['path']}  {str(cfg.get('sha256'))[:12]}")
    else:
        print("     저널       없음 — 이 폴더는 자기 역사를 못 말한다")

    if entry["manifest_steps"]:
        print(f"     매니페스트 스텝 {entry['manifest_steps']} "
              f"({entry['manifest_started']}) · 패키지 {entry['manifest_packages']}개")
        if entry["journal_steps"] and len(entry["manifest_steps"]) < len(entry["journal_steps"]):
            print("                ↑ 저널보다 적다 — 나중 부분 실행이 덮어썼다")

    if entry["parameters"]:
        match = {True: "일치", False: "달라짐 — 기록된 뒤 설정을 고쳤다",
                 None: "대조 불가 — 그 설정 파일이 지금 없다"}[entry["config_match"]]
        print(f"     파라미터   {entry['parameters']}개 기록 · 설정 대조 {match}")
    else:
        print("     파라미터   없음")

    for note in entry["notes"]:
        print(f"     메모       {note}")
    for stray in entry["strays"]:
        print(f"     잔재       {stray} — 이 파일을 쓰는 코드가 없다")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--root", type=Path, default=DEFAULT_ROOT,
                    help=f"워크스페이스 루트 (기본 {DEFAULT_ROOT})")
    ap.add_argument("--detail", nargs="*", default=None,
                    help="이름이 맞는 워크스페이스·폴더를 자세히 (인자 없으면 전부)")
    args = ap.parse_args()

    if not args.root.is_dir():
        print(f"  루트가 없다: {args.root}")
        return 1

    entries = [describe(p) for p in find_result_dirs(args.root)]
    if not entries:
        print(f"  결과 디렉터리를 못 찾았다: {args.root}")
        return 1

    print(f"\n  {args.root}  —  결과 디렉터리 {len(entries)}개\n")
    print_table(entries)

    unknown = [e for e in entries if provenance_verdict(e) == "미상"]
    partial = [e for e in entries if provenance_verdict(e) in ("파라미터만", "스텝목록만")]
    print(f"\n  출처: 저널 {len(entries) - len(unknown) - len(partial)}개 · "
          f"부분 {len(partial)}개 · 미상 {len(unknown)}개")
    if unknown:
        print("  미상 — 어떤 설정으로 만들었는지 남은 게 없다. 소급 복원은 불가하며,")
        print("        근거로 쓰려면 다시 돌려야 한다:")
        for e in unknown:
            print(f"    {e['workspace']}/{e['name']}")

    if args.detail is not None:
        wanted = [w.lower() for w in args.detail]
        for entry in entries:
            if not wanted or any(w in f"{entry['workspace']}/{entry['name']}".lower()
                                 for w in wanted):
                print_detail(entry)
    print()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
