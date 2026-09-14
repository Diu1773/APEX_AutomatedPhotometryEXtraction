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

## 둘로 가른다

    없음을 표시      `None`·`NaN`·빈 표·사유를 담은 값. 아래로 흘러가도 **없음이
                    그대로 드러나므로** 조용해도 해롭지 않다.
    그럴듯한 값      `0`·`False`·기본값. **진짜 결과와 겉보기가 같다** — 색항 0 은
                    진짜 색항 0 과 숫자만 보고는 가를 수 없다. 이쪽이 고칠 것이다.

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


#: 이 이름의 자리에 값을 넣는 것도 「적었다」로 본다. 로그로 흘려보내지 않고
#: **산출물에 사유를 담는** 방식이고, 오히려 더 오래 남는다.
NOTE_TARGETS = ("note", "reason", "status", "error", "warning", "message")


def _assigns_a_reason(handler: ast.ExceptHandler) -> bool:
    """사유를 담을 자리에 값을 넣나 — `out["qc_note"] = ...` 같은 것."""
    for node in ast.walk(handler):
        targets: list = []
        if isinstance(node, ast.Assign):
            targets = list(node.targets)
        elif isinstance(node, (ast.AnnAssign, ast.AugAssign)):
            targets = [node.target]
        for t in targets:
            name = ast.unparse(t).lower()
            if any(k in name for k in NOTE_TARGETS):
                return True
    return False


def _records_something(handler: ast.ExceptHandler) -> bool:
    """그 갈래가 기록을 남기거나 예외를 다시 던지나."""
    if _assigns_a_reason(handler):
        return True
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


def _substituted_values(handler: ast.ExceptHandler) -> list[ast.expr]:
    """그 갈래가 만들어 내보내는 값들. 비었으면 계산을 안 잇는다."""
    out: list[ast.expr] = []
    for node in ast.walk(handler):
        if isinstance(node, ast.Assign):
            out.append(node.value)
        elif isinstance(node, (ast.AnnAssign, ast.AugAssign)) and node.value is not None:
            out.append(node.value)
        elif isinstance(node, ast.Return) and node.value is not None:
            out.append(node.value)
    return out


#: 「없다」를 뜻하는 값. 이런 값은 아래로 흘러가도 없음이 그대로 드러나므로
#: 조용해도 해롭지 않다. 해로운 것은 **그럴듯한 값**이다 — 색항 0 은 진짜
#: 색항 0 과 겉보기가 같아서 숫자만 보고는 가를 수 없다.
_ABSENCE = {"None", "{}", "[]", "()", "''", '""', "pd.DataFrame()", "set()",
            "dict()", "list()", "tuple()"}


def _marks_absence(value: ast.expr) -> bool:
    src = ast.unparse(value)
    if src in _ABSENCE or "nan" in src.lower():
        return True
    # 사유를 담은 값도 조용하지 않다 — 읽는 쪽이 무슨 일이 있었는지 알 수 있다.
    lowered = src.lower()
    return any(t in lowered for t in ("'error'", '"error"', "str(exc)", "str(e)",
                                      "reason", "_error:", "failed"))


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
            if _records_something(node):
                continue
            values = _substituted_values(node)
            if not values:
                continue
            rows.append({
                "file": rel, "line": node.lineno,
                "layer": "calc" if rel.startswith(CALC_PREFIXES) else "gui",
                "kind": ("absence" if all(_marks_absence(v) for v in values)
                         else "plausible"),
            })
    return rows


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="조용한 fallback 을 센다")
    ap.add_argument("--layer", choices=["calc", "gui", "all"], default="all")
    ap.add_argument("--kind", choices=["absence", "plausible", "all"],
                    default="all", help="plausible 이 고쳐야 할 쪽이다")
    ap.add_argument("--show", action="store_true", help="자리를 하나씩 찍는다")
    ap.add_argument("--json", default="", help="이 경로에 목록을 쓴다")
    a = ap.parse_args(argv)

    rows = scan()
    everything = list(rows)
    if a.layer != "all":
        rows = [r for r in rows if r["layer"] == a.layer]
    if a.kind != "all":
        rows = [r for r in rows if r["kind"] == a.kind]

    split = Counter((r["layer"], r["kind"]) for r in everything)
    by_file = Counter(r["file"] for r in rows)
    print(f"값을 만들면서 조용한 갈래 {len(everything)} 곳")
    print("  없음을 표시 (해롭지 않다) · 그럴듯한 값으로 갈아치움 (고쳐야 한다)")
    for layer in ("calc", "gui"):
        print(f"    {layer:<5} {split.get((layer, 'absence'), 0):>4}"
              f" · {split.get((layer, 'plausible'), 0):>4}")
    print()
    print(f"골라 본 것 {len(rows)} 곳 · 파일 {len(by_file)} 개")
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
