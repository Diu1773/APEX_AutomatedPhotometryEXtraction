"""한 장도 못 풀었는데 초록으로 끝나면 안 된다.

**Step 5 는 몇 장을 풀었든 `OK` 를 돌려주고 있었다.** kb26 의 BANZAI 프레임
서른 장을 돌렸더니 `0/30 frames solved` 라고 적으면서 상태는 `ok` 였고, 파이프라인
요약도 초록이었다. 사람이 로그를 읽으면 알 수 있지만 **자동 실행은 못 알아챈다.**

그런데 그 자리에서 실제로 멈췄어야 했느냐 하면 그것도 아니다. 이 프레임들은
관측소가 붙여 놓은 측성 해를 헤더에 갖고 있고, 6 단계는 그것을 읽어 쓴다
(`refbuild.py` 의 `w.has_celestial` 갈래). 실제로 마스터 목록 2,500 별이
제대로 나왔다. **즉 「한 장도 못 풀었다」와 「이어서 못 간다」는 다른 말이다.**

그래서 갈라야 할 것은 이것이다.

    한 장도 못 풀었고 헤더에도 측성 해가 없다   → 이어서 갈 수 없다. FAILED
    한 장도 못 풀었지만 헤더에 측성 해가 있다   → 갈 수 있다. OK 로 두되
                                                 무엇을 쓰는지 분명히 적는다
    몇 장이라도 풀었다                          → 지금까지대로 OK

`is_complete` 가 요구하는 `wcs_solve_summary.csv` 는 한 장도 못 풀면 아예 안
쓰인다. **스스로 「완료 아님」이라고 판정할 산출물을 안 내고서 OK 를 돌려주던
셈이다.**
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO))

from apex.pipeline.base import StepStatus  # noqa: E402
from apex.pipeline.steps.wcs import WcsStep  # noqa: E402
from apex.utils.step_paths import step1_dir, step4_dir, step5_wcs_dir  # noqa: E402


class _Log:
    def info(self, *a, **k):
        pass

    def warning(self, *a, **k):
        pass

    def exception(self, *a, **k):
        pass


class _Ctx:
    """러너가 넘기는 것 가운데 이 단계가 실제로 쓰는 것만 담은 껍데기."""

    def __init__(self, data_dir: Path, result_dir: Path, cache_dir: Path):
        self.data_dir = data_dir
        self.result_dir = result_dir
        self.logger = _Log()

        class _P:
            pass

        p = _P()
        p.cache_dir = str(cache_dir)

        class _Params:
            pass

        self.params = _Params()
        self.params.P = p


#: 하늘 좌표가 제대로 붙은 헤더 (M67 언저리, 0.58 초각/화소)
_WCS_KEYS = dict(CTYPE1="RA---TAN", CTYPE2="DEC--TAN",
                 CRVAL1=132.8242, CRVAL2=11.7997,
                 CRPIX1=32.0, CRPIX2=32.0,
                 CD1_1=-0.000161, CD1_2=0.0, CD2_1=0.0, CD2_2=0.000161)


def _frame(path: Path, with_wcs: bool) -> None:
    hdu = fits.PrimaryHDU(np.zeros((64, 64), dtype=np.float32))
    hdu.header["EXPTIME"] = 40.0
    hdu.header["FILTER"] = "rp"
    if with_wcs:
        for k, v in _WCS_KEYS.items():
            hdu.header[k] = v
    hdu.writeto(path)


def _workspace(tmp_path: Path, with_wcs: bool) -> _Ctx:
    data = tmp_path / "data"
    result = tmp_path / "result"
    cache = tmp_path / "cache"
    for d in (data, result, cache):
        d.mkdir(parents=True, exist_ok=True)
    names = ["f1.fits", "f2.fits"]
    for n in names:
        _frame(data / n, with_wcs)
    step1_dir(result).mkdir(parents=True, exist_ok=True)
    (step1_dir(result) / "selection.json").write_text(
        json.dumps({"filenames": names}), encoding="utf-8")
    step4_dir(result).mkdir(parents=True, exist_ok=True)
    return _Ctx(data, result, cache)


def _run_with_no_solutions(ctx: _Ctx, monkeypatch) -> object:
    """풀이 엔진이 한 장도 못 푼 상황을 만든다."""
    import apex.analysis.wcs_solve as ws

    monkeypatch.setattr(ws, "run_wcs_solve",
                        lambda *a, **k: {"total": 2, "ok": 0, "wcs_qc_pass": 0})
    monkeypatch.setattr(ws, "resolve_wcs_engine", lambda *_: "internal")
    step5_wcs_dir(ctx.result_dir).mkdir(parents=True, exist_ok=True)
    return WcsStep().run(ctx)


def test_nothing_solved_and_no_header_wcs_is_a_failure(tmp_path, monkeypatch):
    """이어서 갈 수 없으면 실패로 알린다 — 지금은 초록으로 끝난다."""
    ctx = _workspace(tmp_path, with_wcs=False)
    res = _run_with_no_solutions(ctx, monkeypatch)
    assert res.status == StepStatus.FAILED, (
        f"한 장도 못 풀고 헤더에도 좌표가 없는데 상태가 {res.status!r} 였다")


def test_nothing_solved_but_header_wcs_present_still_proceeds(tmp_path, monkeypatch):
    """헤더에 측성 해가 있으면 이어서 간다 — 관측소가 이미 푼 프레임의 길이다."""
    ctx = _workspace(tmp_path, with_wcs=True)
    res = _run_with_no_solutions(ctx, monkeypatch)
    assert res.status == StepStatus.OK
    assert "헤더" in res.message or "header" in res.message.lower(), (
        f"무엇을 쓰는지 안 적었다: {res.message!r}")


def test_the_message_still_says_how_many_were_solved(tmp_path, monkeypatch):
    """몇 장을 풀었는지는 예전처럼 그대로 적는다."""
    ctx = _workspace(tmp_path, with_wcs=True)
    res = _run_with_no_solutions(ctx, monkeypatch)
    assert "0/2" in res.message


def test_some_solved_is_unchanged(tmp_path, monkeypatch):
    """한 장이라도 풀면 지금까지와 같다."""
    import apex.analysis.wcs_solve as ws

    ctx = _workspace(tmp_path, with_wcs=False)
    monkeypatch.setattr(ws, "run_wcs_solve",
                        lambda *a, **k: {"total": 2, "ok": 1, "wcs_qc_pass": 1})
    monkeypatch.setattr(ws, "resolve_wcs_engine", lambda *_: "internal")
    step5_wcs_dir(ctx.result_dir).mkdir(parents=True, exist_ok=True)
    res = WcsStep().run(ctx)
    assert res.status == StepStatus.OK
    assert "1/2" in res.message


@pytest.mark.parametrize("with_wcs", [True, False])
def test_missing_selection_is_still_blocked(tmp_path, with_wcs):
    """앞 단계 산출물이 없으면 예전처럼 BLOCKED 다."""
    ctx = _workspace(tmp_path, with_wcs=with_wcs)
    (step1_dir(ctx.result_dir) / "selection.json").unlink()
    assert WcsStep().run(ctx).status == StepStatus.BLOCKED
