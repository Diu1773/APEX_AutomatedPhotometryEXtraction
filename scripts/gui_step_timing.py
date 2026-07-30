"""실제 GUI 스텝 창을 하나씩 열어 로드 시간과 상태를 잰다.

앱을 사람이 클릭하며 도는 대신, 같은 창 클래스를 같은 순서로 생성해서
「열리는가 · 얼마나 걸리는가 · 표가 채워지는가 · 화면에 들어가는가」를 기록한다.
사용자가 실제로 보는 것과 같은 코드경로다(오프스크린 렌더만 다름).

    .venv-deploy/Scripts/python.exe -X utf8 scripts/gui_step_timing.py \
        --mode cmd --params <parameters.toml> --out validation/lc_gui_run/cmd_steps.json

주의: 창 생성 전에 apex/.state/<mode>/project_state.json 의 file_selection 을
이 파라미터 파일의 경로로 맞춘다. 안 그러면 이전 세션의 프로젝트가
params.P.result_dir 을 덮어써서 엉뚱한 자료를 잰다(main_window.py 의
_bootstrap_file_selection_state).
"""
from __future__ import annotations

import argparse
import json
import os
import time
import traceback
from pathlib import Path

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import sys

REPO = Path(__file__).absolute().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))


def _sync_state(mode: str, data_dir: str, result_dir: str) -> None:
    """실사용 재현 — GUI 는 Step 1 에서 고른 폴더를 상태에 저장해 두고 복원한다."""
    state_path = REPO / "apex" / ".state" / mode / "project_state.json"
    state_path.parent.mkdir(parents=True, exist_ok=True)
    doc = {}
    if state_path.exists():
        try:
            doc = json.loads(state_path.read_text(encoding="utf-8"))
        except Exception:
            doc = {}
    step_data = doc.setdefault("step_data", {})
    fs = step_data.setdefault("file_selection", {})
    fs["data_dir"] = data_dir
    fs["result_dir"] = result_dir
    state_path.write_text(json.dumps(doc, indent=2), encoding="utf-8")


def _table_shape(win) -> str:
    """창 안의 첫 QTableWidget 모양 — 표가 비었는지 한눈에 본다."""
    from PyQt5.QtWidgets import QTableWidget

    tables = win.findChildren(QTableWidget)
    if not tables:
        return "-"
    parts = []
    for t in tables[:2]:
        filled = sum(
            1
            for c in range(t.columnCount())
            if t.item(0, c) is not None and t.item(0, c).text().strip()
        ) if t.rowCount() else 0
        parts.append(f"{t.rowCount()}x{t.columnCount()}(1행채움 {filled})")
    return " | ".join(parts)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=("cmd", "lc"), required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", default="")
    ap.add_argument("--steps", default="", help="예: 0-7 (기본: 전부)")
    args = ap.parse_args()

    from PyQt5.QtWidgets import QApplication

    app = QApplication.instance() or QApplication(sys.argv)

    from apex.config.parameters_cmd import Parameters as PCmd
    from apex.config.parameters_lc import Parameters as PLc

    params = (PCmd if args.mode == "cmd" else PLc)(args.params)
    _sync_state(args.mode, str(params.P.data_dir), str(params.P.result_dir))

    from apex.gui.main_window import MainWindowWorkflow

    t0 = time.perf_counter()
    mw = MainWindowWorkflow(mode=args.mode, param_file=args.params)
    boot = time.perf_counter() - t0
    print(f"[boot] 메인창 {boot:.2f}s | result_dir={mw.params.P.result_dir}", flush=True)

    n_steps = len(mw.step_names)
    if args.steps:
        lo, _, hi = args.steps.partition("-")
        idx = range(int(lo), int(hi or lo) + 1)
    else:
        idx = range(n_steps)

    rows = []
    for i in idx:
        if i >= n_steps:
            break
        name = mw.step_names[i]
        rec = {"index": i, "name": name}
        t1 = time.perf_counter()
        try:
            win = mw._open_step_window(i)
            app.processEvents()
            rec["open_s"] = round(time.perf_counter() - t1, 2)
            if win is None:
                rec["status"] = "창 없음(None 반환)"
            else:
                rec["class"] = type(win).__name__
                # 디스크 산출물 복원 — 실제 창이 열릴 때 하는 일
                t2 = time.perf_counter()
                for meth in ("_load_from_disk", "load_existing", "_restore_outputs"):
                    fn = getattr(win, meth, None)
                    if callable(fn):
                        try:
                            fn()
                        except Exception as exc:
                            rec.setdefault("load_warn", str(exc)[:120])
                        break
                for meth in ("update_frame_table", "_refresh_table", "refresh_table"):
                    fn = getattr(win, meth, None)
                    if callable(fn):
                        try:
                            fn()
                        except Exception as exc:
                            rec.setdefault("table_warn", str(exc)[:120])
                        break
                app.processEvents()
                rec["load_s"] = round(time.perf_counter() - t2, 2)
                rec["table"] = _table_shape(win)
                sz = win.size()
                rec["size"] = f"{sz.width()}x{sz.height()}"
                tabs = getattr(win, "tabs", None)
                if tabs is not None and hasattr(tabs, "count"):
                    rec["tabs"] = [tabs.tabText(k) for k in range(tabs.count())]
                rec["status"] = "OK"
                try:
                    win.close()
                except Exception:
                    pass
        except Exception as exc:
            rec["open_s"] = round(time.perf_counter() - t1, 2)
            rec["status"] = f"예외: {type(exc).__name__}: {exc}"[:200]
            rec["trace"] = traceback.format_exc()[-600:]
        rows.append(rec)
        print(
            f"  step{i:>2} {name:<24} {rec.get('status','')[:60]:<60} "
            f"open {rec.get('open_s','-')}s load {rec.get('load_s','-')}s "
            f"| {rec.get('table','-')} | {rec.get('size','-')}",
            flush=True,
        )

    out = {
        "mode": args.mode,
        "params": args.params,
        "result_dir": str(mw.params.P.result_dir),
        "boot_s": round(boot, 2),
        "steps": rows,
    }
    if args.out:
        p = Path(args.out)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(out, indent=2, ensure_ascii=False), encoding="utf-8")
        print(f"[saved] {p}")
    ok = sum(1 for r in rows if r.get("status") == "OK")
    print(f"[done] {ok}/{len(rows)} 스텝 창이 열렸다")
    return 0 if ok == len(rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
