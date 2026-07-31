"""도구·스텝 창을 실제로 열어 실행하고 화면을 PNG 로 남긴다.

사람이 메뉴를 눌러 여는 것과 같은 진입점(main_window 의 launcher 메서드)을 쓴다.
오프스크린이라 화면은 안 뜨지만 위젯 렌더링은 그대로이므로 `grab()` 이 실제
화면과 같은 그림을 준다.

    python run_tools_capture.py --mode cmd --params <toml> --out <dir>

각 창마다 기록하는 것: 열림/예외 · 소요 · 창 크기 · 표 행수 · 떠오른 다이얼로그
· 캔버스 유무 · PNG 경로. 문제는 삼키지 않고 전부 JSON 에 남긴다.
"""
from __future__ import annotations

import argparse
import json
import os
import time
import traceback
from pathlib import Path

import sys

# 기본은 실제 윈도우 플랫폼이다. offscreen 으로 잡으면 시스템 폰트를 못 읽어
# 캡처에 **글자가 하나도 안 나온다**(색 블록만 찍힌다 — 실측). 창은 show() 하지
# 않고 grab() 만 하므로 화면을 오래 가리지 않는다.
if "--offscreen" in sys.argv:
    os.environ["QT_QPA_PLATFORM"] = "offscreen"

REPO = Path(__file__).absolute().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

DIALOGS: list[str] = []


def _patch_dialogs():
    """모달을 띄우지 않고 내용만 모은다.

    오프스크린에서 모달을 exec 하면 중첩 이벤트루프를 돌다 죽는다(실측: segfault).
    사용자에게 뜨는 메시지 자체가 기록 대상이므로 버리지 않고 모은다.
    """
    from PyQt5.QtWidgets import QMessageBox, QFileDialog, QInputDialog

    def _mk(kind):
        def _f(*a, **k):
            title = str(a[1]) if len(a) > 1 else ""
            body = str(a[2])[:300] if len(a) > 2 else ""
            DIALOGS.append(f"{kind}[{title}] {body}")
            return QMessageBox.Ok if kind != "question" else QMessageBox.No
        return staticmethod(_f)

    for name in ("critical", "warning", "information", "about"):
        setattr(QMessageBox, name, _mk(name))
    setattr(QMessageBox, "question", _mk("question"))
    # 파일 대화상자는 사용자 입력을 기다리므로 무조건 취소로 막는다
    for name in ("getOpenFileName", "getSaveFileName"):
        setattr(QFileDialog, name, staticmethod(lambda *a, **k: ("", "")))
    setattr(QFileDialog, "getOpenFileNames", staticmethod(lambda *a, **k: ([], "")))
    setattr(QFileDialog, "getExistingDirectory", staticmethod(lambda *a, **k: ""))
    setattr(QInputDialog, "getText", staticmethod(lambda *a, **k: ("", False)))
    setattr(QInputDialog, "getItem", staticmethod(lambda *a, **k: ("", False)))


def _pump(app, seconds: float):
    from PyQt5.QtCore import QEventLoop, QTimer

    loop = QEventLoop()
    QTimer.singleShot(int(seconds * 1000), loop.quit)
    loop.exec_()
    app.processEvents()


# 도구별 「주 실행 버튼」 — 사람이 창을 열고 실제로 누르는 그 버튼이다.
# 여기 없는 도구는 열기만 해도 내용이 차거나(뷰어), 입력이 더 필요한 것이다.
RUN_BUTTONS = {
    "step12_isochrone": ("Run MCMC Auto-Fit",),
    "qa_report": ("Generate QA Report",),
    "iraf_photometry": ("Run Comparison",),
    "extinction_fit": ("Run Fit",),
    "gaia_3d_viewer": ("Plot",),
    "variable_star": ("Compute Periodogram", "Analyze", "Run"),
    "transit": ("Fit Transit", "Run Fit", "Run"),
    "eclipsing_binary": ("Run Fit", "Analyze", "Run"),
    "multi_night_merger": ("Merge", "Run Merge"),
    "airmass_debug": ("Scan", "Reload", "Run"),
}


def _click_run(win, tid: str, app, timeout_s: float) -> dict:
    """주 실행 버튼을 눌러 결과가 나올 때까지 기다린다."""
    from PyQt5.QtWidgets import QPushButton

    wanted = RUN_BUTTONS.get(tid, ())
    if not wanted:
        return {"run": "실행 버튼 미지정"}
    target = None
    for name in wanted:
        for b in win.findChildren(QPushButton):
            if b.text().replace("\n", " ").strip() == name and b.isEnabled():
                target = b
                break
        if target is not None:
            break
    if target is None:
        return {"run": f"버튼 없음/비활성: {wanted[0]}"}

    t0 = time.perf_counter()
    try:
        target.click()
    except Exception as exc:
        return {"run": f"클릭 실패: {type(exc).__name__}: {exc}"[:150]}
    app.processEvents()
    # 워커가 도는 동안 기다린다. 버튼이 다시 활성화되면 끝난 것으로 본다.
    deadline = time.perf_counter() + timeout_s
    while time.perf_counter() < deadline:
        _pump(app, 2.0)
        try:
            if target.isEnabled() and time.perf_counter() - t0 > 4.0:
                break
        except RuntimeError:      # 위젯이 파괴된 경우
            break
    return {"run": f"{wanted[0]} 클릭", "run_s": round(time.perf_counter() - t0, 1)}


def _inspect(win) -> dict:
    from PyQt5.QtWidgets import QTableWidget, QPushButton
    try:
        from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FC
    except Exception:
        FC = None

    info: dict = {}
    sz = win.size()
    info["size"] = f"{sz.width()}x{sz.height()}"
    tables = win.findChildren(QTableWidget)
    if tables:
        info["tables"] = [f"{t.rowCount()}x{t.columnCount()}" for t in tables[:3]]
    if FC is not None:
        info["canvases"] = len(win.findChildren(FC))
    btns = [b.text().replace("\n", " ") for b in win.findChildren(QPushButton) if b.text().strip()]
    info["buttons"] = btns[:12]
    return info


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mode", choices=("cmd", "lc"), required=True)
    ap.add_argument("--params", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--settle", type=float, default=6.0, help="창마다 기다리는 초")
    ap.add_argument("--run", action="store_true",
                    help="각 도구의 주 실행 버튼을 눌러 결과까지 만든다")
    ap.add_argument("--run-timeout", type=float, default=180.0,
                    help="실행 버튼을 누른 뒤 기다리는 최대 초")
    ap.add_argument("--offscreen", action="store_true", help="(폰트가 빠진다 — 진단용)")
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    from PyQt5.QtWidgets import QApplication

    app = QApplication.instance() or QApplication(sys.argv)
    # 실제 앱은 세 진입점 모두 이걸 부른다 — 안 부르면 캡처가 실사용과 달라진다
    from apex.gui.theme import apply_theme
    apply_theme(app)
    _patch_dialogs()

    from apex.config.parameters_cmd import Parameters as PCmd
    from apex.config.parameters_lc import Parameters as PLc

    params = (PCmd if args.mode == "cmd" else PLc)(args.params)

    # 실사용 재현 — GUI 는 Step 1 에서 고른 폴더를 상태에 저장해 두고 복원한다
    state_path = REPO / "apex" / ".state" / args.mode / "project_state.json"
    state_path.parent.mkdir(parents=True, exist_ok=True)
    doc = {}
    if state_path.exists():
        try:
            doc = json.loads(state_path.read_text(encoding="utf-8"))
        except Exception:
            doc = {}
    fs = doc.setdefault("step_data", {}).setdefault("file_selection", {})
    fs["data_dir"] = str(params.P.data_dir)
    fs["result_dir"] = str(params.P.result_dir)
    state_path.write_text(json.dumps(doc, indent=2), encoding="utf-8")

    from apex.gui.main_window import MainWindowWorkflow
    from apex.gui.tools.registry import iter_tools_for_mode

    mw = MainWindowWorkflow(mode=args.mode, param_file=args.params)
    print(f"[boot] {args.mode} | result_dir={mw.params.P.result_dir}", flush=True)

    targets: list[tuple[str, str, callable]] = []
    if args.mode == "cmd":
        targets.append(("step12_isochrone", "Step 12 Isochrone Model",
                        lambda: mw._open_step_window(11)))
    for spec in iter_tools_for_mode(args.mode):
        fn = getattr(mw, spec.launcher, None)
        if fn is None:
            targets.append((spec.id, spec.label, None))
            continue
        targets.append((spec.id, spec.label, fn))

    records = []
    for tid, label, fn in targets:
        DIALOGS.clear()
        rec = {"id": tid, "label": label, "mode": args.mode}
        if fn is None:
            rec["status"] = "launcher 없음"
            records.append(rec)
            print(f"  {tid:<22} launcher 없음", flush=True)
            continue
        t0 = time.perf_counter()
        win = None
        # launcher 는 창을 반환하지 않고 self.qa_window / self.iraf_window 처럼
        # 도구마다 다른 속성에 넣는다. 호출 전후로 **main_window 의 새 속성**을
        # 먼저 보고, 없을 때만 최상위 창으로 넘어간다. 최상위 창만 보면 숨은
        # 헬퍼 위젯을 잡아 100x30 짜리 껍데기를 캡처하게 된다(실측).
        from PyQt5.QtWidgets import QWidget

        attrs_before = dict(vars(mw))
        before = set(id(w) for w in app.topLevelWidgets())
        try:
            win = fn()
            app.processEvents()
            _pump(app, args.settle)
            if win is None:
                cands = [v for k, v in vars(mw).items()
                         if k not in attrs_before and isinstance(v, QWidget)]
                cands += [v for k, v in vars(mw).items()
                          if isinstance(v, QWidget) and v is not attrs_before.get(k)
                          and v is not mw and k.endswith(("window", "tool", "viewer", "dialog"))]
                # 실제 내용이 있는 창을 고른다 — 껍데기보다 큰 것
                cands = [c for c in cands if c is not mw]
                if cands:
                    win = max(cands, key=lambda w: w.size().width() * w.size().height())
            if win is None:
                fresh = [w for w in app.topLevelWidgets()
                         if id(w) not in before and w is not mw]
                if fresh:
                    win = max(fresh, key=lambda w: w.size().width() * w.size().height())
            rec["open_s"] = round(time.perf_counter() - t0, 2)
            if win is None:
                rec["status"] = "창 참조를 못 얻음"
            else:
                rec["class"] = type(win).__name__
                if args.run:
                    DIALOGS.clear()
                    rec.update(_click_run(win, tid, app, args.run_timeout))
                    if DIALOGS:
                        rec["run_dialogs"] = list(DIALOGS)
                rec.update(_inspect(win))
                png = out_dir / f"{args.mode}_{tid}.png"
                try:
                    win.grab().save(str(png))
                    rec["png"] = png.name
                except Exception as exc:
                    rec["png_error"] = str(exc)[:150]
                rec["status"] = "OK"
        except Exception as exc:
            rec["open_s"] = round(time.perf_counter() - t0, 2)
            rec["status"] = f"예외: {type(exc).__name__}: {exc}"[:200]
            rec["trace"] = traceback.format_exc()[-800:]
        if DIALOGS:
            rec["dialogs"] = list(DIALOGS)
        # 도구가 params 를 건드리면 뒤따라 여는 창이 엉뚱한 곳을 본다.
        # 실제로 variable_star 가 result_dir 대신 sci/result 를 보고 있었다.
        now_result = str(mw.params.P.result_dir)
        if now_result != str(params.P.result_dir):
            rec["result_dir_changed"] = now_result
            print(f"      * result_dir 이 바뀌었다 -> {now_result}", flush=True)
        records.append(rec)
        print(f"  {tid:<22} {rec.get('status','')[:60]:<60} {rec.get('open_s','-')}s "
              f"{rec.get('size','-')} tables={rec.get('tables','-')} "
              f"canvas={rec.get('canvases','-')}", flush=True)
        if DIALOGS:
            for d in DIALOGS[:3]:
                print(f"      ! {d[:160]}", flush=True)
        try:
            if win is not None:
                win.close()
                app.processEvents()
        except Exception:
            pass

    out_json = out_dir / f"{args.mode}_tools.json"
    out_json.write_text(json.dumps(records, indent=2, ensure_ascii=False), encoding="utf-8")
    ok = sum(1 for r in records if r.get("status") == "OK")
    print(f"[done] {ok}/{len(records)} 창이 열렸다 -> {out_json}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
