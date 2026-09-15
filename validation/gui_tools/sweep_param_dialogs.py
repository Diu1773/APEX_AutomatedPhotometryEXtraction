"""파라미터 창 전부를 재서 가로로 넘치는 곳을 찾는다 (2026-09-15).

Step 10 의 파라미터 창이 12 px 넘쳤고 원인은 폼 레이아웃의 성질이었다
(`Main/FAILURES.md` F-325) — **이름 칸과 값 칸의 너비는 줄마다 공유되므로, 가장 긴
이름과 가장 넓은 값이 서로 다른 줄에 있어도 둘이 더해진다.** 그 모양이 다른 창에도
있는지는 안 봤다.

**눈으로 짐작하지 않고 창마다 실제로 띄워서 잰다.** 재는 일 자체는
`measure_param_dialog.py` 가 이미 하므로, 여기서는 그것을 창마다 **따로 띄워서**
돌리고 결과를 한 표로 모은다.

창마다 프로세스를 새로 띄우는 이유는 둘이다. `QDialog.exec_` 을 가로채는 덧칠이
한 프로세스 안에 남고, Qt 위젯도 창 사이에 상태를 물려받아서 **같은 프로세스에서
이어 재면 앞 창이 뒤 창의 숫자를 바꾼다.**

## 반드시 실제 화면에서 잰다

**넘침은 화면 크기에 딸린 값이다.** 창은 `clamp_to_screen` 으로 화면에 맞춰 줄어들고,
줄어든 폭보다 내용이 넓으면 넘친 것으로 잡힌다. 그래서 화면이 좁으면 멀쩡한 창도
넘친 것으로 나온다.

    offscreen 플랫폼의 가상 화면   800 × 600
    이 기계의 실제 화면            2560 × 1440

2026-09-15 에 `QT_QPA_PLATFORM=offscreen` 으로 한 번 돌렸다가 창 셋이 249~1169 px
넘치는 것으로 나왔는데, 전부 창이 634 px 로 잘려서 생긴 허위였다. 그래서 **기본은
실제 화면**이고, 잰 화면 크기를 결과에 같이 적는다(`--offscreen` 은 일부러 좁은
화면을 시험할 때만 쓴다).

## 덮는 범위 (2026-09-15, CMD 기준)

열두 창 가운데 **다섯만 실제로 잰다.** 나머지를 못 재는 이유가 서로 다르므로 갈라
둔다 — 「넘치는 창 0 개」는 **잰 다섯에 대한 말**이지 열둘 전부가 아니다.

    잰 창          4 · 6 · 8 · 9 · 10
    버튼이 없다    1 · 2 · 7 · 11 · 12   — 파라미터 창 자체가 없어 잴 것이 없다
    앱이 죽는다    3 · 5                 — **도구가 아니라 앱 쪽 결함이다**

**앞서 이 자리에 「앱이 아니라 도구의 한계」라고 적었던 것은 틀렸다.** 처음에는
`QDialog.exec_` 을 가로친 탓으로 보고 재는 시점을 옮기고, 다시 진짜 모달로 띄우는
길까지 만들어 봤지만 세 방법 모두 같은 자리에서 죽었다. 그래서 **아무것도 덧칠하지
않고 Parameters 만 눌러 보니**, 3 번은 「Photometry Parameters」가 뜬 뒤 닫는 순간에,
5 번은 대화상자가 뜨기도 전에 프로세스가 `0xC0000409`(스택 버퍼 오버런)로 죽었다.
같은 방식으로 4 번은 끝까지 멀쩡하다. 재현은
`Main/FAILURES.md` 의 해당 항목에 적어 두었고, **이 둘은 못 본 채로 남아 있다.**

실행:
    python -X utf8 validation/gui_tools/sweep_param_dialogs.py
    python -X utf8 validation/gui_tools/sweep_param_dialogs.py --mode lc
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).absolute().parents[2]
TOOL = REPO / "validation/gui_tools/measure_param_dialog.py"
OUT = REPO / "validation/gui_tools/sweep_param_dialogs.json"

#: 재는 창 번호. CMD·LC 둘 다 1~12 이고, 파라미터 버튼이 없는 창은 스스로 걸러진다.
STEPS = range(1, 13)
PAT_SIZE = re.compile(r"보이는 너비 (\d+) · 내용 최소 너비 (\d+) · 넘침 (-?\d+) px")
PAT_TITLE = re.compile(r"^창: (.*?)\s+(\d+)x(\d+)", re.M)
PAT_BAR = re.compile(r"가로 스크롤바: (\w+)")


def _screen_size(offscreen: bool) -> tuple[int, int]:
    """이 측정이 어떤 화면에 견준 것인지 — 넘침 숫자의 뜻이 여기에 달려 있다."""
    env = dict(os.environ)
    if offscreen:
        env["QT_QPA_PLATFORM"] = "offscreen"
    else:
        env.pop("QT_QPA_PLATFORM", None)
    code = ("import sys;from PyQt5.QtWidgets import QApplication;"
            "a=QApplication(sys.argv[:1]);g=a.primaryScreen().availableGeometry();"
            "print(g.width(),g.height())")
    p = subprocess.run([sys.executable, "-c", code], capture_output=True,
                       text=True, env=env, timeout=120)
    try:
        w, h = (p.stdout or "").strip().split()[-2:]
        return int(w), int(h)
    except (ValueError, IndexError):
        return 0, 0


def measure_one(step: int, mode: str, params: str,
                offscreen: bool = False) -> dict:
    """창 하나를 새 프로세스에서 재고 숫자만 뽑아 온다."""
    env = dict(os.environ)
    if offscreen:
        env["QT_QPA_PLATFORM"] = "offscreen"
    else:
        # **실제 화면에서 재야 한다.** offscreen 은 800×600 을 보고하므로 창이
        # 그만큼 잘리고, 멀쩡한 창도 넘친 것으로 나온다(위 설명).
        env.pop("QT_QPA_PLATFORM", None)
    cmd = [sys.executable, "-X", "utf8", str(TOOL), "--step", str(step),
           "--mode", mode]
    if params:
        cmd += ["--params", params]
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, env=env,
                           timeout=420, encoding="utf-8", errors="replace")
    except subprocess.TimeoutExpired:
        return dict(step=step, mode=mode, status="시간초과")
    out = (p.stdout or "") + (p.stderr or "")
    m = PAT_SIZE.search(out)
    if not m:
        # 파라미터 버튼이 없거나 창을 못 여는 것은 결함이 아니다 — 다만 **무엇이었는지
        # 남겨야** 재지 못한 창을 「넘치지 않는 창」과 헷갈리지 않는다.
        #
        # **도구가 내는 말을 이름으로 찾는다.** 여기서는 stdout 뒤에 stderr 를 이어
        # 붙이므로 「마지막 줄」은 언제나 경고 쪽이고, 실제로 그렇게 집었다가 사유
        # 자리에 설정 파일 경로가 찍혔다.
        known = ("Parameters 버튼을 못 찾았다", "번 창을 못 열었다",
                 "대화상자가 안 잡혔다")
        note = next((ln.strip() for ln in out.splitlines()
                     if any(k in ln for k in known)), "")
        if not note:
            tail = [ln.strip() for ln in out.splitlines() if ln.strip()]
            note = tail[-1][:90] if tail else "출력 없음"
        return dict(step=step, mode=mode, status="못 쟀다", note=note[:90])
    title = PAT_TITLE.search(out)
    bar = PAT_BAR.search(out)
    widest = [ln.strip() for ln in out.splitlines() if ln.strip().startswith("-")][:3]
    return dict(step=step, mode=mode, status="쟀다",
                title=(title.group(1) if title else ""),
                dialog_w=int(title.group(2)) if title else None,
                viewport_w=int(m.group(1)), content_min_w=int(m.group(2)),
                overflow=int(m.group(3)),
                h_scrollbar=(bar.group(1) if bar else "?"),
                widest=widest)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="파라미터 창을 전부 재서 넘침을 찾는다")
    ap.add_argument("--mode", default="cmd", choices=["cmd", "lc", "both"])
    ap.add_argument("--params", default="", help="apex_config.json 경로")
    ap.add_argument("--offscreen", action="store_true",
                    help="일부러 좁은 가상 화면(800×600)에서 잰다 — 평소엔 쓰지 않는다")
    a = ap.parse_args(argv)
    modes = ["cmd", "lc"] if a.mode == "both" else [a.mode]

    # **잰 화면을 먼저 적는다.** 넘침이 화면 크기에 딸린 값이라, 화면을 안 적으면
    # 나중에 이 숫자가 무엇에 견준 것인지 알 수 없다.
    screen = _screen_size(a.offscreen)
    print(f"잰 화면: {screen[0]} × {screen[1]}"
          + ("  ← 일부러 좁게 잡았다" if a.offscreen else ""))
    if screen[0] < 1280:
        print("  **주의: 화면이 좁아 멀쩡한 창도 넘친 것으로 나온다.**")
    print()

    rows = []
    for mode in modes:
        print(f"=== {mode.upper()} 모드 — 창 {STEPS.start}~{STEPS.stop - 1}")
        print(f"{'창':>4}{'상태':>8}{'창 폭':>8}{'보이는':>8}{'내용 최소':>10}"
              f"{'넘침':>8}  가로 스크롤바")
        print("-" * 64)
        for step in STEPS:
            r = measure_one(step, mode, a.params, a.offscreen)
            rows.append(r)
            if r["status"] != "쟀다":
                print(f"{step:>4}{r['status']:>8}   {r.get('note', '')}")
                continue
            mark = "  ← 넘친다" if r["overflow"] > 0 else ""
            print(f"{step:>4}{'ok':>8}{r['dialog_w'] or 0:>8}{r['viewport_w']:>8}"
                  f"{r['content_min_w']:>10}{r['overflow']:>8}  "
                  f"{r['h_scrollbar']}{mark}")
        print()

    bad = [r for r in rows if r["status"] == "쟀다" and r["overflow"] > 0]
    print(f"넘치는 창 {len(bad)} 개")
    for r in bad:
        print(f"  {r['mode']} {r['step']} 번 ({r.get('title', '')}) — "
              f"{r['overflow']} px")
        for w in r.get("widest", []):
            print(f"      {w}")
    OUT.write_text(json.dumps(dict(screen=list(screen),
                                   offscreen=bool(a.offscreen), rows=rows),
                              ensure_ascii=False, indent=2, default=str),
                   encoding="utf-8")
    print(f"\n썼다: {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
