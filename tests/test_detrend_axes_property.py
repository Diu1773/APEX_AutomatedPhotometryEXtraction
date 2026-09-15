"""읽기 전용 속성 위에 창이 대입하려다 LC 11 번이 통째로 안 열렸다.

`DetrendNightMergeWindow` 는 `DetrendRunner` 를 상속하는데, 그 실행기가
`ax_raw`·`ax_corr`·`ax_diag` 를 **읽기 전용 property** 로 두었다. 창은 자기 캔버스에
세 축을 만들어 `self.ax_raw = ...` 로 대입하므로 여기서 바로 터진다.

    AttributeError: property 'ax_raw' of 'DetrendNightMergeWindow' object has no setter
    step10_detrend_merge.py:893  in setup_step_ui

**창 하나가 못 뜨는 정도가 아니라 그 단계 전체가 막힌다** — `__init__` 이
`setup_step_ui()` 를 부르므로 창 생성 자체가 실패한다(2026-09-15 실측).

고침은 세 속성에 setter 를 붙여 **대입한 축이 이긴다**로 만드는 것이다. 읽기만 하는
헤드리스 경로는 예전대로 필요할 때 축을 만들어 쓰고, 창이 대입하면 그 축을 쓴다.

여기서 지키는 것은 둘이다 — **대입이 되는가**, 그리고 **대입하지 않았을 때 헤드리스가
여전히 제 축을 만드는가**.
"""
from __future__ import annotations

import pytest


class _Axes:
    """축 대신 쓰는 표식. 무엇이 어디에 들어갔는지만 보면 된다."""

    def __init__(self, name: str):
        self.name = name


@pytest.fixture()
def runner():
    """실행기만 세운다 — 전체 초기화 없이 축 속성만 시험한다."""
    from apex.analysis.light_curve.detrend_runner import DetrendRunner

    return DetrendRunner.__new__(DetrendRunner)


def test_assigning_one_axis_sticks(runner):
    runner.ax_raw = _Axes("raw")
    assert isinstance(runner.ax_raw, _Axes)
    assert runner.ax_raw.name == "raw"


def test_assigning_all_three_keeps_them_apart(runner):
    """창은 셋을 잇달아 대입한다 — 뒤의 대입이 앞의 것을 지우면 안 된다."""
    runner.ax_raw = _Axes("raw")
    runner.ax_corr = _Axes("corr")
    runner.ax_diag = _Axes("diag")
    assert [runner.ax_raw.name, runner.ax_corr.name, runner.ax_diag.name] == [
        "raw", "corr", "diag"]


def test_reading_without_assigning_still_builds_axes(runner):
    """헤드리스 경로는 대입 없이 읽기만 한다 — 그 길이 막히면 배치 실행이 죽는다."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    runner._plot_figure = lambda: Figure()
    assert runner.ax_raw is not None
    assert runner.ax_corr is not None
    assert runner.ax_diag is not None
    assert runner.ax_raw is not runner.ax_corr


def test_one_assignment_does_not_erase_the_other_two(runner):
    """하나만 대입하고 나머지를 읽으면, 대입한 것은 그대로고 둘은 새로 만들어진다."""
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    from matplotlib.figure import Figure

    runner._plot_figure = lambda: Figure()
    mine = _Axes("mine")
    runner.ax_raw = mine
    assert runner.ax_corr is not None
    assert runner.ax_raw is mine, "다른 축을 읽었더니 대입한 축이 날아갔다"
