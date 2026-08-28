"""APEX 산출물이 디스크를 어디에 쓰고 있는지 센다. **아무것도 지우지 않는다.**

    python -X utf8 scripts/disk_survey.py                 # 기본 뿌리 전부
    python -X utf8 scripts/disk_survey.py E:/APEX_validation/reprocess
    python -X utf8 scripts/disk_survey.py --depth 2       # 더 잘게

왜 지우지 않는가. 2026-08-20 에 중복 삭제로 30 GB 를 회수했는데 **안전 검사가 두 번
모자랐고 둘 다 되돌려서 구했다**(F-061·F-062).

- `validation/` 은 `E:\APEX_validation_output` 을 가리키는 링크다. 레포 안처럼
  보이는 경로가 실제로는 외장 드라이브이고, **git 이 추적하는 파일이 「삭제 가능」
  구역 안에 있었다.**
- 도달성 검사가 `validation/` 만 훑어서 **`tests/` 가 이름으로 적어 둔 보정 프레임을
  지웠다.** 그 테스트는 실패가 아니라 **스킵**했고, 총계만 보면 통과처럼 보였다.

그래서 이 도구는 판단을 사람에게 넘긴다. 「지워도 되는 것」이 아니라
**「이만큼이 재생성 가능해 보이고, 이것들은 무슨 일이 있어도 건드리면 안 된다」**를
낸다.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

DEFAULT_ROOTS = [
    Path("E:/APEX_validation/reprocess"),
    Path("E:/APEX_validation_output"),
]

# 이름만으로 재생성 가능하다고 말할 수 있는 것. 나머지는 전부 「판단 필요」로 둔다.
REGENERABLE_DIRS = {"cache", "__pycache__", ".pytest_cache"}

# 산출물의 출처를 말해 주는 파일. 작고, 지우면 폴더가 자기를 설명하지 못한다.
RECORD_NAMES = {"apex_journal.jsonl", "pipeline_run.json", "parameters_used.json"}

PRODUCT_SUFFIXES = {".csv", ".tsv", ".ecsv", ".png", ".pdf", ".svg", ".md"}
BULK_SUFFIXES = {".fits", ".fit", ".fts", ".npy", ".npz"}


def human(n: int) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.1f} {unit}" if unit != "B" else f"{n} B"
        n /= 1024.0
    return f"{n:.1f} TB"


def tracked_files() -> set[Path]:
    """git 이 추적하는 파일의 실제 경로. 링크를 따라 해석한다.

    `validation/` 이 외장 드라이브를 가리키는 링크라, 레포 상대경로만 모으면
    E 드라이브 쪽 실제 경로와 대조가 안 된다 — F-061 이 정확히 그 구멍이었다.
    """
    try:
        out = subprocess.run(
            ["git", "ls-files", "-z"], cwd=REPO, capture_output=True, check=True
        ).stdout.decode("utf-8", "replace")
    except Exception as exc:  # git 이 없거나 저장소가 아니면 보호를 포기하지 않는다
        print(f"!! git ls-files 실패({exc}) — 보호 목록을 만들 수 없다. 중단한다.")
        raise SystemExit(2)
    paths = set()
    for rel in out.split("\0"):
        if not rel:
            continue
        try:
            paths.add(Path(os.path.realpath(REPO / rel)))
        except OSError:
            pass
    return paths


def referenced_names() -> set[str]:
    """레포 어딘가가 **이름으로** 적어 둔 파일들.

    `tests/` 만이 아니라 레포 전체를 훑는다. F-062 는 `validation/` 만 훑어서
    `tests/` 가 적어 둔 보정 프레임을 놓친 사고였다.
    """
    names: set[str] = set()
    exts = ("*.py", "*.md", "*.json", "*.toml", "*.cfg", "*.txt")
    skip = {".git", ".venv", ".venv-deploy", "node_modules", "__pycache__"}
    for pattern in exts:
        for path in REPO.rglob(pattern):
            if any(part in skip for part in path.parts):
                continue
            try:
                text = path.read_text(encoding="utf-8", errors="ignore")
            except OSError:
                continue
            for token in text.replace("\\", "/").split("/"):
                token = token.strip("\"'`(), \n\t")
                if "." in token and 4 < len(token) < 120:
                    names.add(token)
    return names


def classify(path: Path, parts: tuple[str, ...]) -> str:
    if any(p in REGENERABLE_DIRS for p in parts):
        return "재생성 가능(캐시)"
    if path.name in RECORD_NAMES:
        return "출처 기록"
    suffix = path.suffix.lower()
    if suffix in PRODUCT_SUFFIXES:
        return "표·그림"
    if suffix in BULK_SUFFIXES:
        return "영상·배열"
    return "그 밖"


def survey(root: Path, depth: int, tracked: set[Path], named: set[str]) -> None:
    if not root.exists():
        print(f"\n=== {root} — 없음 ===")
        return

    print(f"\n=== {root} ===")
    by_group: dict[str, int] = defaultdict(int)
    by_kind: dict[str, int] = defaultdict(int)
    kind_files: dict[str, int] = defaultdict(int)
    protected_tracked = 0
    protected_named: set[str] = set()
    total = 0
    n_files = 0

    root_parts = len(root.parts)
    for dirpath, dirnames, filenames in os.walk(root, followlinks=False):
        here = Path(dirpath)
        rel_parts = here.parts[root_parts:]
        group = "/".join(rel_parts[:depth]) or "."
        for name in filenames:
            f = here / name
            try:
                size = f.stat().st_size
            except OSError:
                continue
            total += size
            n_files += 1
            by_group[group] += size
            kind = classify(f, rel_parts)
            by_kind[kind] += size
            kind_files[kind] += 1
            try:
                if Path(os.path.realpath(f)) in tracked:
                    protected_tracked += 1
            except OSError:
                pass
            if name in named:
                protected_named.add(name)

    print(f"  파일 {n_files:,} 개 · {human(total)}")

    print(f"\n  {'무엇':<20} {'크기':>12} {'파일 수':>10}")
    for kind, size in sorted(by_kind.items(), key=lambda kv: -kv[1]):
        print(f"  {kind:<20} {human(size):>12} {kind_files[kind]:>10,}")

    print(f"\n  큰 폴더 (깊이 {depth})")
    for group, size in sorted(by_group.items(), key=lambda kv: -kv[1])[:15]:
        print(f"  {human(size):>12}  {group}")

    print("\n  건드리면 안 되는 것")
    print(f"  {'git 이 추적하는 파일':<28} {protected_tracked:>6,} 개")
    print(f"  {'레포가 이름으로 적어 둔 파일':<28} {len(protected_named):>6,} 개")
    if protected_tracked:
        print("  ! 이 뿌리 안에 git 추적 파일이 있다 — 링크로 이어진 레포 구역이다."
              " 통째 삭제 금지 (F-061)")

    cache = by_kind.get("재생성 가능(캐시)", 0)
    if cache:
        print(f"\n  이름만으로 재생성 가능하다고 말할 수 있는 것: {human(cache)}"
              f" ({cache / total * 100:.1f} %)" if total else "")
        print("  나머지는 이 도구가 판정하지 않는다. 재실행 비용과 논문 인용 여부를"
              " 아는 사람이 정한다.")


def main() -> None:
    ap = argparse.ArgumentParser(description="APEX 산출물 디스크 조사 (읽기 전용)")
    ap.add_argument("roots", nargs="*", type=Path, default=None)
    ap.add_argument("--depth", type=int, default=1, help="폴더를 몇 단계까지 묶을지")
    args = ap.parse_args()

    roots = args.roots or DEFAULT_ROOTS
    print("이 도구는 아무것도 지우지 않는다. 세기만 한다.")
    tracked = tracked_files()
    print(f"git 추적 파일 {len(tracked):,} 개를 보호 목록에 올렸다.")
    named = referenced_names()
    print(f"레포가 이름으로 적어 둔 파일 이름 {len(named):,} 개를 모았다.")

    for root in roots:
        survey(Path(root), args.depth, tracked, named)


if __name__ == "__main__":
    main()
