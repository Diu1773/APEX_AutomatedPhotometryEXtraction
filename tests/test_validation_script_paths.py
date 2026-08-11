"""`validation/` is a junction, so scripts under it must not call `resolve()`.

On this machine `validation/` is a directory junction onto the E: drive. A
script that computes its repository root with `Path(__file__).resolve()`
follows the link and lands outside the repository — `E:\\` instead of the
checkout — so the `sys.path` entry it then inserts does not contain the `apex`
package and every import of it fails.

The failure is quiet in the worst way: scripts that only import their siblings
keep working, because Python puts the script's own directory on `sys.path`
regardless. Only the ones that import `apex` break, and they break at run time,
not at collection. On 2026-08-12 twelve scripts were in this state, including
the three that regenerate paper figures 14, 15 and 16, and nobody had noticed.

`Path(__file__).absolute()` does not follow links and is the fix. This test
keeps it fixed.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO = Path(__file__).absolute().parents[1]
VALIDATION = REPO / "validation"

# Matches a root-directory assignment built from __file__ via resolve(), which
# is the exact construct that follows the junction.
ROOT_VIA_RESOLVE = re.compile(
    r"Path\(__file__\)\.resolve\(\)\s*\.parents\[\d+\]"
)


# Throwaway probes live here and are not tracked, so a fresh one must not turn
# the suite red for everybody.
SCRATCH = {"_scratch"}


def _scripts() -> list[Path]:
    if not VALIDATION.exists():
        return []
    return sorted(
        p for p in VALIDATION.rglob("*.py")
        if p.is_file() and not SCRATCH & set(p.relative_to(VALIDATION).parts)
    )


@pytest.mark.skipif(not VALIDATION.exists(), reason="validation/ 없음")
def test_no_validation_script_derives_its_root_through_resolve():
    offenders = []
    for script in _scripts():
        try:
            text = script.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        if ROOT_VIA_RESOLVE.search(text):
            offenders.append(script.relative_to(REPO).as_posix())

    assert not offenders, (
        "validation/ 은 정션이므로 resolve() 는 레포 밖(E:\\)으로 나간다. "
        "absolute() 로 바꿀 것:\n  " + "\n  ".join(offenders)
    )


@pytest.mark.skipif(not VALIDATION.exists(), reason="validation/ 없음")
def test_absolute_reaches_the_repo_but_resolve_does_not():
    """Pin the behaviour the rule rests on, so the rule is not cargo-cult.

    If `validation/` ever stops being a junction the two calls agree and this
    test says so by skipping rather than by silently passing.
    """
    probe = VALIDATION / "_junction_probe.py"
    if probe.exists():  # pragma: no cover - defensive
        pytest.skip("probe 이름이 이미 쓰이고 있다")

    absolute_root = (VALIDATION / "x.py").absolute().parents[1]
    resolved_root = (VALIDATION / "x.py").resolve().parents[1]
    if absolute_root == resolved_root:
        pytest.skip("validation/ 이 이 체크아웃에서는 정션이 아니다")

    assert absolute_root == REPO
    assert (absolute_root / "apex").is_dir()
    assert not (resolved_root / "apex").is_dir(), (
        "resolve() 가 레포 밖으로 나가는데도 그곳에 apex 가 있다 — "
        "이 테스트의 전제를 다시 확인할 것"
    )
