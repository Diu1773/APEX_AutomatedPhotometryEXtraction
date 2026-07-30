"""LC Step 9 — 비교성 QC + 광곡선 생성 (실제 창의 실행 경로).

    python lc_step9.py <parameters.toml>
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
    w = mw._open_step_window(9)
    app.processEvents()
    print(f"[step9] 창 {time.perf_counter()-t0:.1f}s | {type(w).__name__}", flush=True)

    for name in ("_load_from_disk", "load_selections", "_auto_load_ids", "load_state"):
        fn = getattr(w, name, None)
        if callable(fn):
            try:
                fn(); app.processEvents()
                print(f"[step9] {name} OK", flush=True)
            except Exception as exc:
                print(f"[step9] {name} 실패: {type(exc).__name__}: {str(exc)[:90]}", flush=True)

    tgt = getattr(w, "target_id", None)
    comps = getattr(w, "comparison_ids", None) or getattr(w, "active_comp_ids", None)
    print(f"[step9] 타깃={tgt} 비교성={comps}", flush=True)

    # 비교성 QC
    t1 = time.perf_counter()
    try:
        w.run_comp_qc(); app.processEvents(); _pump(app, 3)
        print(f"[step9] run_comp_qc {time.perf_counter()-t1:.1f}s", flush=True)
    except Exception as exc:
        print(f"[step9] run_comp_qc 실패: {type(exc).__name__}: {str(exc)[:120]}", flush=True)

    # 광곡선 생성 — 입력칸에 ID 를 직접 넣고 runtime_mode 로 동기 실행한다.
    target = sys.argv[2] if len(sys.argv) > 2 else "115"
    comps_arg = sys.argv[3] if len(sys.argv) > 3 else "121,36,194"
    w.target_edit.setText(target)
    w.comp_edit.setText(comps_arg)
    app.processEvents()
    print(f"[step9] 입력: 타깃 {target} / 비교성 {comps_arg} | datasets={len(getattr(w,'datasets',[]) or [])}", flush=True)

    t2 = time.perf_counter()
    try:
        w.runtime_mode = True
        res = w.build_light_curve()
        app.processEvents()
        _pump(app, 5)
        if isinstance(res, dict):
            keys = sorted(res)[:10]
            print(f"[step9] build 결과 keys={keys}", flush=True)
            for k in ("n_points", "rms", "path", "output", "csv"):
                if k in res:
                    print(f"    {k}: {str(res[k])[:120]}", flush=True)
        print(f"[step9] build_light_curve {time.perf_counter()-t2:.1f}s", flush=True)
    except Exception as exc:
        import traceback
        print(f"[step9] build_light_curve 실패: {type(exc).__name__}: {str(exc)[:160]}", flush=True)
        print(traceback.format_exc()[-700:], flush=True)

    for name in ("save_state", "save_frame_excludes"):
        fn = getattr(w, name, None)
        if callable(fn):
            try:
                fn(); app.processEvents()
            except Exception as exc:
                print(f"[step9] {name} 실패: {str(exc)[:80]}", flush=True)
    _pump(app, 4)
    print(f"[step9] 총 {time.perf_counter()-t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
