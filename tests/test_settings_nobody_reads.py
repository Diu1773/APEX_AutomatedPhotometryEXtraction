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

# 2026-08-18: all 72 were deleted. They were legacy — `idmatch_*` alone was 19
# rows for a matching stage that no longer exists — so the ceiling is now zero
# and stays there. A setting that nothing reads is a setting that should not
# have a row.
CEILING = {"cmd": 0, "lc": 0}
GUI_ONLY_CEILING = {"cmd": 42, "lc": 34}

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


def test_a_dialog_cannot_offer_a_setting_nobody_reads():
    """돌려도 아무 일이 안 생기는 손잡이를 창에 내놓을 수 없어야 한다."""
    pytest.importorskip("PyQt5")
    import os

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    from apex.gui.workflow.param_dialog import specs_from_map

    # 실제 죽은 설정은 이제 0 개다. 그래서 가드를 시험하려면 하나를 만들어 넣는다 —
    # 맵에 행은 있고 읽는 코드는 없는 상태를, 다음에 그런 게 생겼을 때와 똑같이.
    import apex.gui.workflow.param_dialog as pd_mod

    pd_mod._unread_for.cache_clear()
    monkey = frozenset({"apcorr_small_scale"})
    original = pd_mod._unread_for
    pd_mod._unread_for = lambda mode: monkey
    try:
        with pytest.raises(KeyError, match="읽지 않는다"):
            specs_from_map(["apcorr_small_scale"], mode="cmd")
    finally:
        pd_mod._unread_for = original
        pd_mod._unread_for.cache_clear()

    # 살아 있는 설정은 그대로 만들어져야 한다.
    assert specs_from_map(["apcorr_small_scale"], mode="cmd")
