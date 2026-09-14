"""저장소 루트의 pytest 설정 — 끊긴 접합(junction) 때문에 수집이 죽는 것을 막는다.

`validation/` 은 외장 디스크 `E:` 를 가리키는 **디렉터리 접합**이다(20 GB 를 E 로
옮긴 결정, `Main/DECISIONS.md`). **E 가 빠져 있으면 그 항목은 남아 있는데 열리지
않는다.** 윈도에서 pytest 는 루트의 자식마다 `os.path.samefile` 을 부르므로
(`_pytest/main.py`, 짧은 경로 대응 #11895) 거기서 이렇게 죽는다.

    FileNotFoundError: [WinError 3] 지정된 경로를 찾을 수 없습니다:
        '...\\Automated_Photometry_EXtraction\\validation'
    Interrupted: 1 error during collection

**시험을 한 개도 못 돌린다. 그리고 이 메시지는 이유를 말해 주지 않는다** — 읽는
사람은 저장소가 깨진 줄 안다. 실제로는 외장 디스크를 꽂으면 그만이다.

여기서 하는 일은 둘이다.

    끊긴 접합을 수집에서 뺀다        그래야 나머지 시험이 돈다
    **무엇을 왜 뺐는지 적는다**      그래야 읽는 사람이 디스크를 꽂는다

**정상일 때는 아무것도 안 한다.** 접합이 살아 있으면 목록이 비므로 동작이 그대로다.
"""

from __future__ import annotations

import os
from pathlib import Path

_ROOT = Path(__file__).absolute().parent

#: 끊긴 채 발견된 항목. 아래 헤더 훅이 이것을 보고 한 줄 적는다.
_DANGLING: list[str] = []


def _dangling_entries() -> list[str]:
    """이름은 있는데 열리지 않는 루트 항목.

    끊긴 접합·심링크가 이렇다 — `os.listdir` 에는 나오지만 `stat` 이 실패한다.
    숨김 폴더와 `.git` 은 애초에 수집 대상이 아니므로 건너뛴다.
    """
    found: list[str] = []
    try:
        names = os.listdir(_ROOT)
    except OSError:                      # 루트를 못 읽으면 여기서 할 일이 없다
        return found
    for name in names:
        if name.startswith("."):
            continue
        try:
            (_ROOT / name).stat()
        except OSError:
            found.append(name)
    return found


_DANGLING = _dangling_entries()

#: pytest 가 수집에서 제외할 경로. 정상일 때는 빈 목록이다.
collect_ignore = list(_DANGLING)

# **시험 파일을 빼는 일은 여기서 못 한다.** pytest 는 어떤 경로를 건너뛸지 정할 때
# **그 경로의 바로 위 폴더에 있는 conftest** 의 `collect_ignore` 만 본다
# (`_pytest/main.py` 의 `_getconftest_pathlist(..., path=collection_path.parent)`).
# 그래서 `tests/` 안의 모듈은 `tests/conftest.py` 가 뺀다 — 거기에 짝이 있다.


def pytest_report_header(config) -> list[str]:
    """수집에서 뺀 것이 있으면 보고서 머리에 이유와 함께 적는다.

    조용히 빼면 「왜 저 폴더의 시험이 안 돌았지」를 물을 방법이 없다.
    """
    if not _DANGLING:
        return []
    # **머리말은 영문으로 적는다.** 한글은 콘솔 코드페이지에 따라 깨지거나
    # `외...` 로 escape 되어 나온다(2026-09-14 실측). 이 줄은 사람이 상황을
    # 알아채라고 있는 것이라 안 깨지는 쪽이 먼저다. 저장소의 다른 로그 줄도
    # 영문이다 (`[ZP QC] bright->faint drift: ...`).
    joined = ", ".join(sorted(_DANGLING))
    return [
        f"[conftest] skipped unreadable path(s): {joined}",
        "[conftest] this is a directory junction to an external drive that is "
        "not connected - plug it in and the path comes back",
    ]
