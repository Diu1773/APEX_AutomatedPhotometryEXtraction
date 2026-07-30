"""LC Step 8 — 타깃 + 비교성을 명시 지정 (자동선택이 안 붙어 수동 경로 사용).

    python lc_step8_manual.py <parameters.toml> <target_id> <comp_id,comp_id,...>
"""
from __future__ import annotations
import os, sys, time, warnings
warnings.filterwarnings("ignore")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

REPO = r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction"
if REPO not in sys.path:
    sys.path.insert(0, REPO)

from PyQt5.QtWidgets import QApplication


def _row_of(tbl, wanted: str):
    for i in range(tbl.rowCount()):
        it = tbl.item(i, 0)
        if it and it.text().strip() == wanted:
            return i
    return None


def main() -> int:
    param_file = sys.argv[1]
    target_id = sys.argv[2]
    comps = [c.strip() for c in sys.argv[3].split(",") if c.strip()]

    app = QApplication.instance() or QApplication(sys.argv)
    from apex.gui.main_window import MainWindowWorkflow

    t0 = time.perf_counter()
    mw = MainWindowWorkflow(mode="lc", param_file=param_file)
    w = mw._open_step_window(8)
    app.processEvents()
    w.load_master_catalogs()
    app.processEvents()
    tbl = w.master_table
    print(f"[step8] 창+카탈로그 {time.perf_counter()-t0:.1f}s | 별 {tbl.rowCount()}개", flush=True)

    if _row_of(tbl, target_id) is None:
        print(f"[step8] 타깃 ID={target_id} 없음")
        return 2

    # 선택 저장은 지연 실행이고 **마지막 상태**를 쓴다. 그래서 순서가 중요하다:
    # 비교성을 먼저 다 찍고 타깃을 마지막에 잡은 뒤 루프를 돌려야 둘 다 남는다.
    ok = []
    for cid in comps:
        r = _row_of(tbl, cid)
        if r is None:
            print(f"[step8]   비교성 {cid} 없음"); continue
        tbl.selectRow(r); tbl.setCurrentCell(r, 0); app.processEvents()
        try:
            w.toggle_comparison_selected(); app.processEvents()
            ok.append(cid)
        except Exception as exc:
            print(f"[step8]   비교성 {cid} 실패: {type(exc).__name__}: {exc}")
    print(f"[step8] 비교성 {len(ok)}개 지정: {ok}", flush=True)

    r = _row_of(tbl, target_id)
    tbl.selectRow(r); tbl.setCurrentCell(r, 0); app.processEvents()
    w.set_target_selected(); app.processEvents()
    print(f"[step8] 타깃 = {target_id}", flush=True)

    roles = {}
    for i in range(tbl.rowCount()):
        it = tbl.item(i, 7)
        k = (it.text().strip() if it else "") or "-"
        roles[k] = roles.get(k, 0) + 1
    print(f"[step8] 역할 분포: {roles}", flush=True)

    try:
        w.save_master_catalog(log_action="headless step8 manual")
    except Exception as exc:
        print(f"[step8] 마스터 저장 실패: {type(exc).__name__}: {exc}", flush=True)

    # 선택 저장은 _queue_selection_save() 로 지연 실행된다 — 이벤트 루프를
    # 돌려 타이머가 발화하게 두지 않으면 comparison_ids 가 빈 채로 남는다.
    for fn_name in ("_flush_selection_save", "_save_selections", "save_selections"):
        fn = getattr(w, fn_name, None)
        if callable(fn):
            try:
                fn()
                print(f"[step8] 선택 즉시 저장({fn_name})", flush=True)
                break
            except Exception:
                pass
    else:
        from PyQt5.QtCore import QEventLoop, QTimer
        loop = QEventLoop()
        QTimer.singleShot(6000, loop.quit)
        loop.exec_()
        app.processEvents()
        print("[step8] 지연 저장 대기 6s", flush=True)
    print(f"[step8] 총 {time.perf_counter()-t0:.1f}s", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
