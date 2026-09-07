"""Step 0 창에서 필터 이름을 그 자리에서 고칠 수 있는가.

**왜 창이어야 하나.** 관측소마다 FILTER 헤더를 자기 방식으로 적는다. 이름을
코드에 박으면 새 관측소마다 우리 코드를 고쳐야 하고, 설정 파일로 빼기만 하면
자료를 앞에 둔 사람에게 JSON 을 열게 하는 셈이다. 스캔이 헤더에 무엇이 있는지
이미 알고 있으니 그것을 보여 주고 그 자리에서 고치게 한다
(사용자 지적, 2026-09-08 · `Main/OPERATOR.md` C-201).

여기서 보는 것 넷.

1. 스캔한 헤더 값이 표에 그대로 뜨는가 (많은 것부터).
2. 표에서 고치면 판정이 바로 바뀌는가.
3. 그것이 설정 파일에 남아 다시 열어도 살아 있는가.
4. **빈 칸으로 지우면 설정에서도 없어지는가** — 새 표만 쓰면 파일에 옛 줄이
   남아 다음에 열 때 되살아난다.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pytest

pytest.importorskip("PyQt5")  # GUI 없는 CI 에서는 건너뛴다

from apex.analysis.calibration_scan import FrameInfo  # noqa: E402
from apex.utils.astro_utils import (  # noqa: E402
    normalize_filter_name,
    register_filter_aliases,
)

REPO = Path(__file__).absolute().parents[1]


@pytest.fixture(autouse=True)
def _clean_aliases():
    register_filter_aliases({})
    yield
    register_filter_aliases({})


@pytest.fixture
def window(tmp_path, qapp):
    from apex.config.parameters_cmd import read_params
    from apex.core.project_state import ProjectState
    from apex.gui.workflow.step0_detector_calibration import (
        DetectorCalibrationWindow,
    )

    shutil.copy(REPO / "parameters.example.json", tmp_path / "apex_config.json")
    params = read_params(str(tmp_path / "apex_config.json"))
    win = DetectorCalibrationWindow(params, ProjectState(tmp_path))
    win._frames = (
        [FrameInfo(path=f"/x/g{i}.fits.fz", ftype="light", exp=10.0,
                   filt="gp", night="20210317", name=f"g{i}.fits.fz")
         for i in range(60)]
        + [FrameInfo(path=f"/x/v{i}.fits.fz", ftype="light", exp=30.0,
                     filt="Bessell-V", night="20210317", name=f"v{i}.fits.fz")
           for i in range(8)]
    )
    win._rebuild_filter_table()
    yield win, tmp_path
    win.close()


@pytest.fixture(scope="session")
def qapp():
    from PyQt5.QtWidgets import QApplication
    import os
    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    app = QApplication.instance() or QApplication([])
    return app


def _row_for(table, raw):
    for r in range(table.rowCount()):
        if table.item(r, 0).text() == raw:
            return r
    raise AssertionError(f"{raw!r} not in the table")


def test_the_header_values_are_shown_most_common_first(window):
    win, _ = window
    t = win.filter_table
    assert t.rowCount() == 2
    assert t.item(0, 0).text() == "gp"          # 60 장
    assert t.item(0, 1).text() == "60"
    assert t.item(1, 0).text() == "Bessell-V"   # 8 장
    assert win.filter_box.isVisible() or t.rowCount() > 0


def test_a_name_apex_knows_is_prefilled(window):
    win, _ = window
    assert win.filter_table.item(_row_for(win.filter_table, "gp"), 2).text() == "g"


def test_a_name_apex_does_not_know_passes_through(window):
    win, _ = window
    row = _row_for(win.filter_table, "Bessell-V")
    assert win.filter_table.item(row, 2).text() == "Bessell-V"


def test_editing_the_key_changes_the_verdict_at_once(window):
    win, _ = window
    assert normalize_filter_name("Bessell-V") == "Bessell-V"
    win.filter_table.item(_row_for(win.filter_table, "Bessell-V"), 2).setText("V")
    assert normalize_filter_name("Bessell-V") == "V"


def test_the_edit_reaches_the_config_file(window):
    win, work = window
    win.filter_table.item(_row_for(win.filter_table, "Bessell-V"), 2).setText("V")
    saved = json.loads((work / "apex_config.json").read_text(encoding="utf-8"))
    assert saved["filters"]["aliases"]["Bessell-V"] == "V"


def test_it_survives_reopening_the_workspace(window):
    from apex.config.parameters_cmd import read_params

    win, work = window
    win.filter_table.item(_row_for(win.filter_table, "Bessell-V"), 2).setText("V")
    register_filter_aliases({})                      # 세션을 새로 여는 셈
    assert normalize_filter_name("Bessell-V") == "Bessell-V"
    read_params(str(work / "apex_config.json"))      # 설정을 다시 읽으면
    assert normalize_filter_name("Bessell-V") == "V"


def test_clearing_the_cell_removes_it_from_the_file_too(window):
    """새 표만 쓰면 파일에 옛 줄이 남아 다음에 열 때 되살아난다."""
    win, work = window
    row = _row_for(win.filter_table, "Bessell-V")
    win.filter_table.item(row, 2).setText("V")
    assert json.loads((work / "apex_config.json").read_text(
        encoding="utf-8"))["filters"]["aliases"] == {"Bessell-V": "V"}

    win.filter_table.item(row, 2).setText("")
    assert normalize_filter_name("Bessell-V") == "Bessell-V"
    assert json.loads((work / "apex_config.json").read_text(
        encoding="utf-8"))["filters"]["aliases"] == {}


def test_a_users_key_beats_the_built_in_table(window):
    """관측소가 우리 표와 다른 뜻으로 쓰면 그쪽이 옳다."""
    win, _ = window
    assert normalize_filter_name("gp") == "g"
    win.filter_table.item(_row_for(win.filter_table, "gp"), 2).setText("V")
    assert normalize_filter_name("gp") == "V"


def test_the_table_grows_with_its_rows_but_stops(window):
    """필터가 많은 자료에서 이 표가 창을 다 먹으면 안 된다."""
    win, _ = window
    two_rows = win.filter_table.height()

    win._frames = [
        FrameInfo(path=f"/x/{k}{i}.fits.fz", ftype="light", exp=10.0,
                  filt=f"F{k:02d}", night="20210317", name=f"{k}{i}.fits.fz")
        for k in range(12) for i in range(3)
    ]
    win._rebuild_filter_table()
    assert win.filter_table.rowCount() == 12
    many_rows = win.filter_table.height()
    assert many_rows > two_rows
    row_h = win.filter_table.verticalHeader().defaultSectionSize()
    assert many_rows <= two_rows + 5 * row_h + 4, "여섯 줄에서 멈춰야 한다"
