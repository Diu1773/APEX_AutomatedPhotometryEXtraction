"""갈래를 골랐으면 골랐다고 적는다 — 꼴이 한 가지여야 찾을 수 있다.

사장님 교정(2026-09-14, C-237): *"fallback들 있으면 항상 무슨로직으로 들어가는지
디버깅관리도 해야돼; 알지 ?"*

**조용한 fallback 은 나중에 아무도 못 찾는다.** 색 열이 없어 색항을 0 으로 두고
보정한 등급은 색보정을 한 등급과 **숫자만 보고는 가를 수 없다.** 그래서 두 군데에
남긴다 — 로그에 한 줄, 그리고 **산출물 표에 한 칸.** 로그는 그 실행을 지켜본
사람에게만 존재하지만 표는 나중에 파일만 여는 사람에게도 남기 때문이다.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO))

from apex.utils.fallback_log import (  # noqa: E402
    PREFIX, format_fallback, note_fallback,
)


def test_the_line_can_be_grepped_out_of_any_log():
    """머리말이 한 가지라 어떤 로그에서든 한 줄로 뽑힌다."""
    line = format_fallback("zeropoint.color_term[g]", "0", "column missing")
    assert line.startswith(PREFIX)
    assert "zeropoint.color_term[g]" in line
    assert "column missing" in line


def test_it_says_what_was_chosen_not_just_that_something_failed():
    """「실패했다」가 아니라 「무엇으로 갔다」를 적는다.

    읽는 사람이 알아야 하는 것은 **지금 손에 든 값이 무엇이냐**다.
    """
    line = format_fallback("wcs", "header WCS", "internal solver found 0/30")
    assert "header WCS" in line


def test_it_writes_the_line_to_the_given_log():
    lines: list[str] = []
    note_fallback(lines.append, "a.b", "default", "why")
    assert len(lines) == 1 and lines[0].startswith(PREFIX)


def test_without_a_log_it_only_returns_the_line():
    """부르는 쪽이 표에 넣거나 모았다가 낼 수 있게 줄만 돌려준다."""
    assert note_fallback(None, "a.b", "default", "why").startswith(PREFIX)


def test_a_broken_log_does_not_stop_the_calculation():
    """**적다가 터져서 계산이 멈추면 안 된다.**

    이 모듈의 존재 이유가 「기록이 실행을 깨뜨리면 안 된다」인데, 기록하는 줄이
    스스로 예외를 내면 그 원칙을 정확히 뒤집는다 (F-323 과 같은 모양).
    """
    def _explode(_line):
        raise RuntimeError("로그가 터졌다")

    assert note_fallback(_explode, "a.b", "default", "why").startswith(PREFIX)


# ---------------------------------------------------------------------------
# 실제로 붙인 자리 — 색항을 못 붙이고 0 으로 간 경우
# ---------------------------------------------------------------------------


def _coeff() -> pd.DataFrame:
    return pd.DataFrame([{"filter": "g", "zp": -5.0, "ct": 0.3, "ct2": 0.0,
                          "color_col": "g_r", "N": 200, "scatter_rms": 0.02},
                         {"filter": "r", "zp": -5.1, "ct": 0.2, "ct2": 0.0,
                          "color_col": "g_r", "N": 200, "scatter_rms": 0.02}])


def test_the_summary_says_which_filters_lost_their_colour_term():
    """색항을 못 붙인 필터가 표에서 구별된다."""
    from apex.analysis.cmd.zeropoint_runner import build_zp_qc_summary

    summary = build_zp_qc_summary(_coeff(), None, None, None, None, 20.0, {"g"})
    got = dict(zip(summary["filter"], summary["color_term_applied"]))
    assert got["g"] is False or got["g"] == False  # noqa: E712
    assert got["r"] is True or got["r"] == True    # noqa: E712


def test_by_default_every_filter_counts_as_having_its_colour_term():
    """안 넘기면 예전처럼 전부 붙인 것으로 본다 — 없던 경고를 만들지 않는다."""
    from apex.analysis.cmd.zeropoint_runner import build_zp_qc_summary

    summary = build_zp_qc_summary(_coeff(), None, None, None)
    assert bool(summary["color_term_applied"].all())


@pytest.mark.parametrize("skipped", [set(), {"g"}, {"g", "r"}])
def test_the_column_is_always_present(skipped):
    """열이 상황에 따라 생겼다 없어졌다 하면 읽는 쪽이 갈라져야 한다."""
    from apex.analysis.cmd.zeropoint_runner import build_zp_qc_summary

    summary = build_zp_qc_summary(_coeff(), None, None, None, None, 20.0, skipped)
    assert "color_term_applied" in summary.columns
    assert int((~summary["color_term_applied"].astype(bool)).sum()) == len(skipped)


def test_a_filter_with_no_colour_axis_is_not_flagged():
    """색축이 "none" 인 필터는 원래 안 붙이는 것이라 결함이 아니다."""
    from apex.analysis.cmd.zeropoint_runner import build_zp_qc_summary

    co = _coeff()
    co.loc[0, "color_col"] = "none"
    summary = build_zp_qc_summary(co, None, None, None, None, 20.0, set())
    assert bool(summary["color_term_applied"].all())
    assert not np.any(summary["filter"].isna())
