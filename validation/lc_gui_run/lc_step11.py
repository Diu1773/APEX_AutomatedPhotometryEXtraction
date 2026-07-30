"""LC Step 11 — 주기 분석 (실제 창의 실행 경로).

    python lc_step11.py <parameters.toml>
"""
from __future__ import annotations
import os, sys, time, warnings
warnings.filterwarnings("ignore")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO = r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction"
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from PyQt5.QtCore import QEventLoop, QTimer
from PyQt5.QtWidgets import QApplication


def _pump(app, seconds: float):
    loop = QEventLoop()
    QTimer.singleShot(int(seconds * 1000), loop.quit)
    loop.exec_()
    app.processEvents()


def main() -> int:
    param_file = sys.argv[1]
    app = QApplication.instance() or QApplication(sys.argv)
    from apex.gui.main_window import MainWindowWorkflow

    t0 = time.perf_counter()
    mw = MainWindowWorkflow(mode="lc", param_file=param_file)
    w = mw._open_step_window(11)
    app.processEvents()
    print(f"[step11] 창 {time.perf_counter()-t0:.1f}s | {type(w).__name__}", flush=True)

    for name in ("_auto_load_target_id", "_load_from_disk", "load_state"):
        fn = getattr(w, name, None)
        if callable(fn):
            try:
                fn(); app.processEvents()
                print(f"[step11] {name} OK", flush=True)
            except Exception as exc:
                print(f"[step11] {name}: {type(exc).__name__}: {str(exc)[:80]}", flush=True)

    print(f"[step11] lc_path={getattr(w, '_auto_lc_path', None)}", flush=True)

    # 완료 신호를 잡기 위해 결과 훅을 감싼다
    state = {"done": False, "payload": None}
    loop = QEventLoop()

    def _finish(payload=None):
        if payload is not None:
            state["payload"] = payload
        if not state["done"]:
            state["done"] = True
            loop.quit()

    for hook in ("_on_analysis_finished", "_on_period_finished", "_on_worker_finished"):
        orig = getattr(w, hook, None)
        if callable(orig):
            def make(o):
                def wrapped(payload=None, *a, **k):
                    try:
                        return o(payload, *a, **k) if payload is not None else o(*a, **k)
                    finally:
                        _finish(payload)
                return wrapped
            setattr(w, hook, make(orig))
            print(f"[step11] 훅 연결: {hook}", flush=True)

    t1 = time.perf_counter()
    try:
        w._run_analysis()
    except Exception as exc:
        import traceback
        print(f"[step11] _run_analysis 실패: {type(exc).__name__}: {str(exc)[:150]}", flush=True)
        print(traceback.format_exc()[-600:], flush=True)
        return 3

    QTimer.singleShot(600_000, _finish)
    if not state["done"]:
        loop.exec_()
    app.processEvents()
    _pump(app, 3)
    print(f"[step11] 분석 {time.perf_counter()-t1:.1f}s", flush=True)

    p = state["payload"]
    if isinstance(p, dict):
        for k in sorted(p):
            v = p[k]
            if not isinstance(v, (list, tuple, dict)) or len(str(v)) < 90:
                print(f"    {k}: {str(v)[:110]}", flush=True)

    for attr in ("best_period", "_best_period", "result", "_last_result"):
        v = getattr(w, attr, None)
        if v is not None and not callable(v):
            print(f"[step11] {attr}: {str(v)[:160]}", flush=True)

    try:
        w.save_state(); app.processEvents(); _pump(app, 3)
        print("[step11] 저장", flush=True)
    except Exception as exc:
        print(f"[step11] 저장 실패: {str(exc)[:90]}", flush=True)
    print(f"[step11] 총 {time.perf_counter()-t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
