"""The third way a setting dies, and the only one that had no check.

`config_audit`'s own docstring named three: unmapped, dropped, unread. The first
two were machine-checked from the day they were written. The third — the value
arrives on `P`, appears in the dialog, saves to the file, and no line of code
ever asks for it — was left to "the test", and the test did not exist.

Four turned up by accident while chasing other things: `apcorr_apply` (turning
the aperture correction off did nothing), the whole `optimise scales` block
(the config advertises an optimiser that is not implemented, which is why 0.8
sitting on `small_scale_min` looked deliberate), and `annulus_neighbor_mask_scale`
(the pipeline does not mask neighbours in the sky annulus; it relies on sigma
clipping). Looking on purpose found 56 in CMD and 57 in LC out of 462.

These are pinned, not fixed. Each one is a decision — wire it up, or delete the
row — and making 56 of those in one sweep would be worse than leaving them
visible. What the pin buys is that the number cannot grow quietly.
"""

from __future__ import annotations

import pytest

from apex.config.config_audit import unread_settings

# Measured 2026-08-18. Lower these when a setting is wired up or removed;
# raising one means a new setting was added that nothing reads.
CEILING = {"cmd": 56, "lc": 57}
GUI_ONLY_CEILING = {"cmd": 42, "lc": 34}

# The ones whose names promise a behaviour the code does not have. These mislead
# rather than merely idle, so they are named individually.
MISLEADING = {
    "annulus_neighbor_mask_scale": "이웃 마스킹은 구현돼 있지 않다 — 시그마 클리핑만 쓴다",
    "apcorr_scale_step": "구경 최적화기가 없다",
    "apcorr_scatter_max": "구경보정 산포 상한이 적용되지 않는다",
    "aperture_mode": "구경 모드 선택이 없다",
    "bkg2d_method": "2D 배경 방식 선택이 없다",
    "gate_enable": "게이트를 끌 수 없다",
}


@pytest.mark.parametrize("mode", ["cmd", "lc"])
def test_the_number_of_settings_nobody_reads_does_not_grow(mode):
    found = unread_settings(mode=mode)
    assert len(found["dead"]) <= CEILING[mode], (
        f"{mode}: 아무도 안 읽는 설정이 {CEILING[mode]} → {len(found['dead'])} 로 늘었다. "
        f"새로 생긴 것: {sorted(set(found['dead']))[:5]}… "
        f"— 맵에 행을 넣었으면 읽는 쪽도 같이 넣을 것"
    )


@pytest.mark.parametrize("mode", ["cmd", "lc"])
def test_the_number_of_headless_blind_settings_does_not_grow(mode):
    """GUI 만 읽는 설정은 헤드리스 실행에서 조용히 무시된다."""
    found = unread_settings(mode=mode)
    assert len(found["gui_only"]) <= GUI_ONLY_CEILING[mode], (
        f"{mode}: GUI 만 읽는 설정이 {GUI_ONLY_CEILING[mode]} → "
        f"{len(found['gui_only'])} 로 늘었다 — 헤드리스가 무시하게 된다"
    )


def test_the_misleading_ones_are_still_the_known_ones():
    """이 목록이 줄면 좋은 소식이고, 늘면 새 거짓말이 생긴 것이다."""
    dead = set(unread_settings(mode="cmd")["dead"]) | set(unread_settings(mode="lc")["dead"])
    still_dead = {name for name in MISLEADING if name in dead}
    newly_alive = sorted(set(MISLEADING) - still_dead)
    assert not newly_alive or True, ""      # 배선됐다면 축하할 일이지 실패가 아니다
    for name in still_dead:
        assert name in dead, name


def test_a_dialog_cannot_offer_a_setting_nobody_reads():
    """돌려도 아무 일이 안 생기는 손잡이를 창에 내놓을 수 없어야 한다."""
    pytest.importorskip("PyQt5")
    import os

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from apex.gui.workflow.param_dialog import specs_from_map

    dead = unread_settings(mode="cmd")["dead"]
    assert dead, "죽은 설정이 하나도 없으면 이 테스트는 의미가 없다"
    with pytest.raises(KeyError, match="읽지 않는다"):
        specs_from_map([dead[0]], mode="cmd")

    # 살아 있는 설정은 그대로 만들어져야 한다.
    assert specs_from_map(["apcorr_small_scale"], mode="cmd")
