"""Pytest configuration.

Lets the test suite run in a HEADLESS environment (no PyQt5 installed), which
is how the cross-platform `test` CI job validates the headless core install.
When PyQt5 is unavailable, the handful of tests that import real GUI worker /
window classes at module load are skipped from collection; everything else
(the Qt-free analysis, pipeline, utils, and config tests) still runs.

With PyQt5 present (local dev, the `test-gui` CI job), nothing is ignored and
the full suite runs.
"""

from __future__ import annotations

import importlib.util

_HAS_PYQT5 = importlib.util.find_spec("PyQt5") is not None

# Test modules whose MODULE-LEVEL imports pull in real GUI classes (which import
# PyQt5). Verified individually: these fail to import without PyQt5; all other
# test modules import cleanly headless.
#
# Two left this list on 2026-08-16 when Steps 8 and 10 moved their calculation
# to apex.analysis: those tests exercise the photometry itself, which is the
# part that should never have needed a widget toolkit.
_GUI_DEPENDENT_TESTS = [
    "test_iso_cache_worker.py",
    "test_isochrone_fitter_v2.py",
    "test_lc_night_classification.py",
    "test_step6_union_master.py",
    "test_variable_star_phase_plot.py",
]

collect_ignore: list = [] if _HAS_PYQT5 else list(_GUI_DEPENDENT_TESTS)


# ---------------------------------------------------------------------------
# 저장소 밖(외장 디스크)을 가져다 쓰는 시험
# ---------------------------------------------------------------------------
#
# `validation/` 은 외장 디스크 `E:` 를 가리키는 디렉터리 접합이다. **디스크가
# 빠지면 그 폴더에서 가져다 쓰는 시험 모듈이 수집 단계에서 `ModuleNotFoundError`
# 로 죽고, 그 바람에 나머지 시험까지 통째로 멈춘다.**
#
# 루트 `conftest.py` 가 폴더 자체를 빼지만 **시험 파일은 거기서 못 뺀다** —
# pytest 는 어떤 경로를 건너뛸지 정할 때 그 경로 바로 위 폴더의 conftest 만 보기
# 때문이다. 그래서 이 짝이 필요하다.
#
# 조용히 빼지 않는다. 무엇을 왜 뺐는지 보고서 머리에 적는다.

from pathlib import Path as _Path  # noqa: E402

_REPO = _Path(__file__).absolute().parents[1]
_TESTS = _Path(__file__).absolute().parent

#: 열리지 않는 저장소 루트 폴더 이름. 정상이면 빈 목록이다.
_MISSING_TREES = [
    name for name in ("validation",)
    if not _Path(_REPO / name).exists()
]


def _modules_importing(name: str) -> list[str]:
    """그 폴더에서 가져다 쓰는 시험 파일 이름.

    가려내는 꼴은 둘 — `from <name> import ...` 과 `sys.path` 에 그 폴더를 끼워
    넣는 줄. 본문에 이름이 스치기만 한 것(주석·경로 문자열)은 세지 않는다.
    이 저장소에서 이름이 스치는 시험은 스물이 넘고 정말로 가져다 쓰는 것은 셋이다.
    """
    hits: list[str] = []
    for path in sorted(_TESTS.glob("test_*.py")):
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        imports = f"from {name} import" in text or f"from {name}." in text
        on_path = any("sys.path" in line and f'"{name}"' in line
                      for line in text.splitlines())
        if imports or on_path:
            hits.append(path.name)
    return hits


_ORPHANED: list[str] = []
for _name in _MISSING_TREES:
    _ORPHANED += _modules_importing(_name)
collect_ignore += _ORPHANED


def pytest_report_header(config) -> list[str]:
    """무엇을 왜 뺐는지 적는다 — 조용히 빼면 안 돌았다는 것도 모른다."""
    if not _ORPHANED:
        return []
    # 머리말은 영문으로 — 콘솔 코드페이지에 따라 한글이 깨진다(루트 conftest 의
    # 같은 훅에 이유를 적어 두었다).
    return [
        "[tests/conftest] skipped tests that import from a missing path: "
        + ", ".join(sorted(_ORPHANED)),
        "[tests/conftest] missing: " + ", ".join(_MISSING_TREES)
        + " - plug the external drive in and they come back",
    ]
