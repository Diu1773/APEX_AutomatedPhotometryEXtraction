"""파라미터 대화상자가 가로로 넘치는 자리를 **재서** 찾는다.

Step 10 의 파라미터 창에 가로 스크롤바가 생긴다. 가로 스크롤바는 읽는 사람이
값을 보려고 좌우로 밀어야 한다는 뜻이라 그 자체로 결함이다.

**어느 위젯이 넘치는지 눈으로 짐작하지 않는다.** 대화상자를 실제로 띄우고
(모달이라 `exec_` 을 가로챈다) 스크롤 안쪽 내용의 최소 너비와 보이는 너비를 재고,
한 줄씩 내려가며 **가장 넓은 위젯을 이름과 함께** 뽑는다. 그러면 고칠 자리가
하나로 좁혀진다.

    python -X utf8 validation/gui_tools/measure_param_dialog.py --step 10

`--params` 를 주면 그 워크스페이스로 연다. 안 주면 저장소 기본 설정을 쓴다.
`--shot` 을 주면 대화상자를 그림으로도 남긴다.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

REPO = Path(__file__).absolute().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

#: 창 번호 -> main_window 의 스텝 인덱스 (창 번호보다 하나 작다)
CAPTURED: list = []


def _patch_exec():
    """모달을 띄우지 않고 붙잡는다.

    `exec_` 은 중첩 이벤트루프를 돌며 사용자를 기다리므로 그대로 두면 스크립트가
    멈춘다. 대화상자는 이미 다 만들어진 뒤에 `exec_` 이 불리므로, 여기서
    가로채면 **사람이 보는 것과 같은 위젯 구성**을 그대로 잴 수 있다.
    """
    from PyQt5.QtWidgets import QDialog

    def _fake(self, *a, **k):
        CAPTURED.append(self)
        return QDialog.Rejected

    QDialog.exec_ = _fake
    QDialog.exec = _fake


def _pump(app, seconds: float) -> None:
    """이벤트 루프를 실제로 돌린다 — `processEvents` 한 번으로는 안 앉는다."""
    from PyQt5.QtCore import QEventLoop, QTimer

    loop = QEventLoop()
    QTimer.singleShot(int(seconds * 1000), loop.quit)
    loop.exec_()
    app.processEvents()


def _widget_label(w) -> str:
    """위젯을 사람이 알아볼 수 있게 한 줄로."""
    from PyQt5.QtWidgets import (QCheckBox, QComboBox, QLabel, QLineEdit,
                                 QPushButton)

    kind = type(w).__name__
    text = ""
    if isinstance(w, (QLabel, QPushButton, QCheckBox)):
        text = w.text().replace("\n", " ")[:70]
    elif isinstance(w, QLineEdit):
        text = (w.text() or w.placeholderText()).replace("\n", " ")[:70]
    elif isinstance(w, QComboBox):
        text = " | ".join(w.itemText(i) for i in range(min(w.count(), 3)))[:70]
    return f"{kind}: {text}" if text else kind


def measure(dialog) -> dict:
    """스크롤 안쪽이 얼마나 넓은지, 무엇이 넓게 만드는지."""
    from PyQt5.QtWidgets import QScrollArea, QWidget

    out: dict = {"title": dialog.windowTitle(),
                 "dialog_w": dialog.width(), "dialog_h": dialog.height()}
    areas = dialog.findChildren(QScrollArea)
    if not areas:
        out["note"] = "스크롤 영역이 없다"
        return out
    area = areas[0]
    content = area.widget()
    vp = area.viewport()
    out["viewport_w"] = vp.width()
    out["content_min_w"] = content.minimumSizeHint().width()
    out["content_hint_w"] = content.sizeHint().width()
    out["overflow"] = out["content_min_w"] - out["viewport_w"]
    out["h_scrollbar_visible"] = bool(area.horizontalScrollBar().isVisible())

    # 가장 넓은 것부터. 자식을 품은 컨테이너는 그 자식이 원인이므로, 잎에
    # 가까운 것만 남기려고 **자식이 없는 위젯**을 따로 모은다.
    rows = []
    for w in content.findChildren(QWidget):
        try:
            mw = w.minimumSizeHint().width()
            hw = w.sizeHint().width()
        except Exception:  # noqa: BLE001
            continue
        if max(mw, hw) < 200:
            continue
        rows.append({
            "min_w": int(mw), "hint_w": int(hw),
            "leaf": not w.findChildren(QWidget),
            "what": _widget_label(w),
        })
    rows.sort(key=lambda r: -max(r["min_w"], r["hint_w"]))
    out["widest"] = rows[:14]
    out["widest_leaves"] = [r for r in rows if r["leaf"]][:10]

    # **줄 단위로 갈라야 고칠 자리가 나온다.** 폼 레이아웃의 최소 너비는
    # 「이름 칸의 최소 + 값 칸의 최소 + 사이 간격」이므로, 어느 줄이 그 합을
    # 키우는지 보면 무엇을 줄여야 하는지 바로 나온다.
    from PyQt5.QtWidgets import QFormLayout

    def _section_of(w) -> str:
        """그 위젯을 품은 접이 섹션의 제목."""
        from apex.gui.workflow.ui_helpers import CollapsibleSection
        node = w
        while node is not None:
            if isinstance(node, CollapsibleSection):
                return node.title() or node.toggle_button.text()
            node = node.parentWidget()
        return "(섹션 밖)"

    all_rows = []
    for form in content.findChildren(QFormLayout):
        rws = []
        for i in range(form.rowCount()):
            lab = form.itemAt(i, QFormLayout.LabelRole)
            fld = form.itemAt(i, QFormLayout.FieldRole)
            span = form.itemAt(i, QFormLayout.SpanningRole)
            lw = lab.widget().minimumSizeHint().width() if lab and lab.widget() else 0
            fw = fld.widget().minimumSizeHint().width() if fld and fld.widget() else 0
            sw = span.widget().minimumSizeHint().width() if span and span.widget() else 0
            total = max(lw + fw + form.horizontalSpacing(), sw)
            name = ""
            if lab and lab.widget():
                name = _widget_label(lab.widget())
            elif span and span.widget():
                name = "(가로지름) " + _widget_label(span.widget())
            what = _widget_label(fld.widget()) if fld and fld.widget() else ""
            holder = ((lab.widget() if lab else None)
                      or (fld.widget() if fld else None)
                      or (span.widget() if span else None))
            rws.append({"total": int(total), "label_w": int(lw),
                        "field_w": int(fw), "span_w": int(sw),
                        "label": name, "field": what,
                        "section": _section_of(holder) if holder else "?"})
        all_rows += rws
    # 폼의 최소 너비는 그 안에서 가장 넓은 줄이 정하므로 전부 모아 정렬한다.
    all_rows.sort(key=lambda r: -r["total"])
    out["widest_rows"] = all_rows[:10]
    return out


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--step", type=int, default=10, help="창 번호 (기본 10)")
    ap.add_argument("--mode", default="cmd", choices=["cmd", "lc"])
    ap.add_argument("--params", default="", help="apex_config.json 경로")
    ap.add_argument("--shot", default="", help="그림으로 남길 경로")
    ap.add_argument("--section", default="", help="이 이름의 섹션만 따로 그림으로")
    ap.add_argument("--collapsed", dest="expand", action="store_false",
                    help="접힌 채로 잰다 (기본은 전부 펼쳐서 잰다)")
    ap.add_argument("--no-autofit", dest="autofit", action="store_false",
                    help="폭 자동 맞추기를 끄고 잰다 — 고침 전후를 같은 자리에서 견줄 때")
    a = ap.parse_args(argv)

    from PyQt5.QtWidgets import QApplication, QPushButton

    from apex.gui.theme import apply_theme
    from apex.utils.app_setup import configure_fonts

    app = QApplication.instance() or QApplication(sys.argv[:1])
    # **폰트를 안 잡으면 굴림으로 떨어져 너비가 실제와 달라진다.** 가로 넘침을
    # 재는 일이라 이 한 줄이 결과를 바꾼다.
    configure_fonts(app)
    apply_theme(app)
    _patch_exec()
    if not a.autofit:
        # **대조군이다.** 고침 전후를 다른 날 다른 상태에서 견주면 무엇이 바뀐
        # 것인지 못 가르므로, 같은 실행 조건에서 끄고 켜 본다.
        from apex.gui.workflow import ui_helpers
        ui_helpers.fit_parameter_dialog_width = lambda _dlg: 0
        print("(폭 자동 맞추기를 껐다 — 대조군)")

    from apex.gui.main_window import MainWindowWorkflow

    mw = MainWindowWorkflow(mode=a.mode, param_file=(a.params or None))
    win = mw._open_step_window(a.step - 1)
    if win is None:
        print(f"{a.step} 번 창을 못 열었다")
        return 1
    app.processEvents()

    # 사람이 누르는 그 버튼을 누른다 — 헤더의 Parameters.
    target = None
    for b in win.findChildren(QPushButton):
        if "param" in b.text().lower() or "파라" in b.text():
            target = b
            break
    if target is None:
        print("Parameters 버튼을 못 찾았다")
        return 1
    target.click()
    app.processEvents()

    if not CAPTURED:
        print("대화상자가 안 잡혔다")
        return 1
    dlg = CAPTURED[-1]
    # **접힌 채로 재면 안 넘친다.** 사용자가 실제로 부딪히는 것은 섹션을 펼친
    # 상태이므로, 기본은 전부 펼쳐서 잰다.
    if a.expand:
        from apex.gui.workflow.ui_helpers import CollapsibleSection
        n = 0
        for sec in dlg.findChildren(CollapsibleSection):
            sec.set_expanded(True)
            n += 1
        print(f"(섹션 {n} 개를 펼쳐서 잰다)")
    # **레이아웃이 앉기 전에 재면 숫자가 앞뒤가 안 맞는다** — 창 폭보다 보이는
    # 너비가 넓게 나온다. 실제로 띄우고 이벤트 루프를 한두 번 돌린 뒤에 잰다
    # (`memory/project_gui_render_harness.md`).
    dlg.show()
    _pump(app, 1.5)
    info = measure(dlg)
    info["expanded"] = bool(a.expand)

    print(f"창: {info.get('title')}  {info.get('dialog_w')}x{info.get('dialog_h')}")
    print(f"보이는 너비 {info.get('viewport_w')} · 내용 최소 너비 "
          f"{info.get('content_min_w')} · 넘침 {info.get('overflow')} px")
    print(f"가로 스크롤바: {info.get('h_scrollbar_visible')}")
    print()
    print("가장 넓은 것 (자식이 없는 위젯만):")
    for r in info.get("widest_leaves", []):
        print(f"  {r['min_w']:>5} / {r['hint_w']:>5}  {r['what']}")
    print()
    print("컨테이너까지 포함:")
    for r in info.get("widest", []):
        mark = " " if r["leaf"] else "▣"
        print(f"  {mark} {r['min_w']:>5} / {r['hint_w']:>5}  {r['what']}")

    print()
    print("가장 넓은 폼 줄 (이름 칸 + 값 칸):")
    for r in info.get("widest_rows", []):
        print(f"  {r['total']:>5} = {r['label_w']:>4} + {r['field_w']:>4}"
              f"{' (가로지름 ' + str(r['span_w']) + ')' if r['span_w'] else ''}"
              f"   [{r.get('section', '?')[:26]}] {r['label'][:40]}"
              f"  |  {r['field'][:24]}")

    if a.shot:
        p = Path(a.shot)
        p.parent.mkdir(parents=True, exist_ok=True)
        dlg.grab().save(str(p))
        print(f"\n그림: {p}")
        # 한 섹션만 따로도 남긴다 — 스크롤 아래에 있는 섹션은 창 그림에 안 잡힌다.
        if a.section:
            from apex.gui.workflow.ui_helpers import CollapsibleSection
            for sec in dlg.findChildren(CollapsibleSection):
                if a.section.lower() in (sec.title() or "").lower():
                    q = p.with_name(p.stem + "_section" + p.suffix)
                    sec.grab().save(str(q))
                    print(f"섹션 그림: {q}  ({sec.width()}x{sec.height()})")
                    break

    out = REPO / "validation/gui_tools/param_dialog_measure.json"
    out.write_text(json.dumps(info, ensure_ascii=False, indent=2),
                   encoding="utf-8")
    print(f"썼다: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
