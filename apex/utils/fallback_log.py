"""갈래를 하나 골랐으면 골랐다고 적는다 — 한 가지 꼴로.

**조용한 fallback 은 나중에 아무도 못 찾는다.** 계산이 원래 길로 못 가서 다른 길로
갔다면, 그 결과를 보는 사람은 **자기가 무엇을 보고 있는지 모른다.** 색 열이 없어
색항을 0 으로 두고 보정한 등급은 색보정을 한 등급과 겉보기가 같다. 숫자만 보고는
가를 방법이 없다.

이 모듈이 하는 일은 **꼴을 하나로 맞추는 것**뿐이다. 줄이 전부 같은 머리말로
시작하면 어떤 로그에서든 한 줄로 뽑을 수 있다.

    [fallback] zeropoint.color_term: 0 (color column 'color_g_r' is missing)

    grep "\\[fallback\\]" run.log

## 어디까지 이 모듈로 적나

**계산 결과를 바꾸는 갈래만.** 창을 닫다 난 예외나 그림 꾸미기가 실패한 것은
여기 해당하지 않는다. 판별은 이 물음 하나다 — *이 갈래로 갔다는 것을 모르면
결과를 잘못 읽게 되는가.* 그렇다면 적는다.

## 로그만으로는 모자란 자리가 있다

**로그는 그 실행을 지켜본 사람에게만 존재한다.** 나중에 산출물 파일만 열어 보는
사람에게는 안 남는다. 그래서 **산출물에 빈칸이나 뜻이 달라진 값이 생기는 갈래는
표에도 사유를 적는다** (예: `zp_qc_summary.csv` 의 `drift_note`).
이 모듈은 그 표를 대신하지 않는다.
"""

from __future__ import annotations

from typing import Callable, Optional

#: 로그에서 한 줄로 뽑을 수 있게 하는 머리말.
PREFIX = "[fallback]"


def format_fallback(where: str, chose: str, why: str) -> str:
    """한 줄로 만든다. 로그를 안 쓰는 곳에서도 이 꼴을 쓴다.

    Parameters
    ----------
    where : 어느 계산인지. ``모듈.항목`` 꼴이 읽기 좋다 (``zeropoint.color_term``).
    chose : **무엇을 골랐는지.** 「실패했다」가 아니라 「무엇으로 갔다」를 적는다 —
        읽는 사람이 알아야 하는 것은 지금 손에 든 값이 무엇이냐이기 때문이다.
    why : 왜 원래 길로 못 갔는지. 없는 열 이름·모자란 개수처럼 **확인할 수 있는
        사실**을 적는다.
    """
    return f"{PREFIX} {where}: {chose} ({why})"


def note_fallback(log: Optional[Callable[[str], None]],
                  where: str, chose: str, why: str) -> str:
    """갈래를 적고 그 줄을 돌려준다.

    `log` 가 ``None`` 이면 적지 않고 줄만 돌려준다 — 부르는 쪽이 표에 넣거나
    모아 두었다가 한꺼번에 낼 수 있게.
    """
    line = format_fallback(where, chose, why)
    if log is not None:
        try:
            log(line)
        except Exception:  # noqa: BLE001
            # **적다가 터져서 계산이 멈추면 안 된다.** 이 모듈의 존재 이유가
            # 「기록이 실행을 깨뜨리면 안 된다」이므로 여기만은 조용히 넘어간다.
            pass
    return line
