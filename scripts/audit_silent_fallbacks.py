"""값을 만들면서 아무 기록도 안 남기는 갈래를 센다.

사장님 지시(2026-09-14, `Main/OPERATOR.md` C-237): *"fallback들 있으면 항상 무슨
로직으로 들어가는지 디버깅관리도 해야돼"*.

**조용한 fallback 은 나중에 아무도 못 찾는다.** 계산이 원래 길로 못 가서 다른 값을
쓰고 이어졌다면, 결과를 보는 사람은 **자기가 무엇을 보고 있는지 모른다.** 색 열이
없어 색항을 0 으로 두고 보정한 등급은 색보정을 한 등급과 숫자만 보고는 가를 수
없다.

## 무엇을 세나

`except` 블록 가운데 **값을 만들어 내는 것**만 센다 — 대입이 있거나 값을 돌려주는
것. 계산이 그 값으로 이어지기 때문이다. 창을 닫다 난 예외처럼 값을 안 만드는 것은
세지 않는다. 그리고 그 안에서 **기록하거나 다시 던지면** 조용하지 않으므로 뺀다.

    세는 것     except 안에 대입/return 이 있고, 로그·예외가 없다
    안 세는 것  값을 안 만드는 정리 코드 · 로그를 남기는 것 · 다시 던지는 것

**이 수를 0 으로 만드는 것이 목표가 아니다.** 어떤 자리는 정말로 조용해도 되고
(기록이 실행을 깨뜨리면 안 되는 자리), 어떤 자리는 화면 꾸미기다. 목표는 **수가
슬그머니 늘지 않게 지켜보는 것**이고, 계산 층부터 줄여 가는 것이다.

실행:
    python -X utf8 scripts/audit_silent_fallbacks.py
    python -X utf8 scripts/audit_silent_fallbacks.py --layer analysis --show
"""
from __future__ import annotations

import argparse
import ast
import json
from collections import Counter
from pathlib import Path

REPO = Path(__file__).absolute().parents[1]

#: 이름에 이것이 들어간 호출이 있으면 「적었다」로 본다.
LOGGY = ("log", "warn", "error", "print", "info", "exception", "debug",
         "critical", "note_fallback")

#: 화면(gui)과 계산을 가른다 — 줄여야 할 곳은 계산 쪽이 먼저다.
CALC_PREFIXES = ("apex/analysis", "apex/core", "apex/utils", "apex/pipeline",
                 "apex/config", "apex/cmd", "apex/lightcurve")


def _records_something(handler: ast.ExceptHandler) -> bool:
    """그 갈래가 기록을 남기거나 예외를 다시 던지나."""
    for node in ast.walk(handler):
        if isinstance(node, ast.Raise):
            return True
        if isinstance(node, ast.Call):
            func = node.func
            name = (func.attr if isinstance(func, ast.Attribute)
                    else func.id if isinstance(func, ast.Name) else "")
            if any(t in name.lower() for t in LOGGY):
                return True
    return False


def _produces_a_value(handler: ast.ExceptHandler) -> bool:
    """그 갈래가 값을 만들어 계산을 잇나."""
    for node in ast.walk(handler):
        if isinstance(node, (ast.Assign, ast.AnnAssign, ast.AugAssign)):
            return True
        if isinstance(node, ast.Return) and node.value is not None:
            return True
    return False


def scan(root: Path = REPO / "apex") -> list[dict]:
    rows: list[dict] = []
    for path in sorted(root.rglob("*.py")):
        rel = path.relative_to(REPO).as_posix()
        try:
            tree = ast.parse(path.read_text(encoding="utf-8"))
        except (OSError, SyntaxError):
            continue
        for node in ast.walk(tree):
            if not isinstance(node, ast.ExceptHandler):
                continue
            if not _produces_a_value(node) or _records_something(node):
                continue
            rows.append({
                "file": rel, "line": node.lineno,
                "layer": "calc" if rel.startswith(CALC_PREFIXES) else "gui",
            })
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="조용한 fallback 을 센다")
    ap.add_argument("--layer", choices=["calc", "gui", "all"], default="all")
    ap.add_argument("--show", action="store_true", help="자리를 하나씩 찍는다")
    ap.add_argument("--json", default="", help="이 경로에 목록을 쓴다")
    a = ap.parse_args(argv)

    rows = scan()
    if a.layer != "all":
        rows = [r for r in rows if r["layer"] == a.layer]

    by_layer = Counter(r["layer"] for r in rows)
    by_file = Counter(r["file"] for r in rows)
    print(f"값을 만들면서 조용한 갈래 {len(rows)} 곳 · 파일 {len(by_file)} 개")
    print(f"  계산 {by_layer.get('calc', 0)} · 화면 {by_layer.get('gui', 0)}")
    print()
    print("많은 파일부터:")
    for f, n in by_file.most_common(15):
        print(f"  {n:>4}  {f}")

    if a.show:
        print()
        for r in rows:
            print(f"  {r['file']}:{r['line']}")

    if a.json:
        Path(a.json).write_text(
            json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
        print(f"\n썼다: {a.json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
