"""헤더에서 만든 설정이 맞는가 — 특히 시야가 돌아 있을 때.

**`CD1_1` 하나로 화소 크기를 재면 안 된다.** CD 행렬은 화소 크기와 회전을 함께
담으므로 `CD1_1 = 화소크기 · cos(각)` 이다. 시야가 90 도 돌아 있으면 `CD1_1` 이
거의 0 이 되고, 거기에 3600 을 곱한 값은 화소 크기가 아니다.

실제 자료에서 그렇게 났다. LCO kb26 의 BANZAI 프레임은 CD 행렬이

    CD1_1 = 2.11e-06    CD1_2 =  1.59e-04
    CD2_1 = -1.59e-04   CD2_2 =  2.11e-06

이라 `|CD1_1|·3600 = 0.0076` 초각이 나오는데 실제 화소는 0.573 초각이다.
**75 배 틀린다.** 제대로 재려면 한 열의 길이를 쓴다.

    화소 크기 = hypot(CD1_1, CD2_1) · 3600

이것이 두 군데를 망친다. **하나**, 헤더에 `PIXSCALE` 이 있는 기기에서는 「둘이
2 % 넘게 다르다」는 거짓 경고가 뜬다. **둘**, `PIXSCALE` 이 없는 기기에서는 그
0.0076 이 설정에 그대로 들어가 측성 해(사진에 하늘 좌표를 붙이는 것)를 깨뜨린다.
둘째가 진짜 고장이다.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import make_lco_config  # noqa: E402

#: 0.573 초각/화소로 90 도 돌아 있는 시야. LCO kb26 의 BANZAI 프레임 값이다.
ROTATED = dict(CD1_1=2.11439029347e-06, CD1_2=0.000159051924498,
               CD2_1=-0.000159051924498, CD2_2=2.11439029347e-06)

#: 같은 화소 크기인데 안 돌아 있는 시야.
UNROTATED = dict(CD1_1=-0.000159065, CD1_2=0.0, CD2_1=0.0, CD2_2=0.000159065)


def _job(tmp_path: Path, extra: dict) -> Path:
    """과학 프레임 한 장을 담은 최소 작업 폴더."""
    d = tmp_path / "inputs" / "science" / "rp"
    d.mkdir(parents=True)
    hdu = fits.PrimaryHDU(np.zeros((8, 8), dtype=np.float32))
    h = hdu.header
    h["OBJECT"] = "M67"
    h["CRVAL1"] = 132.8242
    h["CRVAL2"] = 11.7997
    h["GAIN"] = 1.0
    h["RDNOISE"] = 14.5
    h["SATURATE"] = 102400.0
    h["CCDSUM"] = "1 1"
    for k, v in extra.items():
        h[k] = v
    hdu.writeto(d / "frame-e91.fits")
    return tmp_path


def _scale(job: Path) -> float:
    return float(make_lco_config.build(job)["match"]["pixel_scale_arcsec"])


def test_rotated_field_without_pixscale_still_gets_the_right_scale(tmp_path):
    """**이것이 진짜 고장이다.** PIXSCALE 이 없으면 CD 행렬이 유일한 근거다."""
    got = _scale(_job(tmp_path, ROTATED))
    assert got == pytest.approx(0.5726, abs=0.002), (
        f"돌아 있는 시야에서 화소 크기를 {got} 로 읽었다 — 0.573 이어야 한다")


def test_unrotated_field_gets_the_same_scale(tmp_path):
    """안 돌아 있어도 값이 같아야 한다 — 회전은 화소 크기를 안 바꾼다."""
    assert _scale(_job(tmp_path, UNROTATED)) == pytest.approx(0.5726, abs=0.002)


def test_pixscale_in_the_header_wins_and_raises_no_false_warning(tmp_path, capsys):
    """헤더의 PIXSCALE 이 있으면 그것을 쓰고, 돌아 있다고 경고하지 않는다."""
    job = _job(tmp_path, dict(ROTATED, PIXSCALE=0.58))
    assert _scale(job) == pytest.approx(0.58)
    assert "2 % 넘게 다르다" not in capsys.readouterr().out


def test_a_genuine_mismatch_still_warns(tmp_path, capsys):
    """진짜로 어긋나면 경고는 그대로 떠야 한다 — 경고를 없애는 것이 아니다."""
    job = _job(tmp_path, dict(ROTATED, PIXSCALE=0.40))
    make_lco_config.build(job)
    assert "2 % 넘게 다르다" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# 이미 보정된 프레임
# ---------------------------------------------------------------------------
#
# **아카이브는 원본과 보정본을 같은 이름꼴로 내놓는다.** LCO 는 `e00` 이 원본,
# `e91` 이 관측소 파이프라인(BANZAI)이 보정한 것이다. 보정본에는 bias·dark·flat
# 을 다시 빼면 안 되고, 오버스캔은 이미 잘려 나갔는데 `BIASSEC` 키워드는 헤더에
# 그대로 남아 있다. 그대로 두면 설정이 「오버스캔 36 열을 오른쪽에서 빼고 자른다」
# 로 만들어져 **없는 영역을 자른다.**
#
# 알아보는 근거는 `RLEVEL` 이다. 0 이 원본, 91 이 보정본이다.


def test_a_reduced_frame_turns_calibration_and_overscan_off(tmp_path):
    """RLEVEL 91 이면 다시 보정하지 않고 오버스캔도 건드리지 않는다."""
    job = _job(tmp_path, dict(ROTATED, PIXSCALE=0.58, RLEVEL=91,
                              BIASSEC="[3100:3135,1:2048]"))
    data = make_lco_config.build(job)
    assert data["calibration"]["enabled"] is False
    assert data["calibration"]["overscan"]["enable"] is False


def test_a_raw_frame_still_calibrates_and_trims(tmp_path):
    """원본은 그대로다 — 알아보는 것이지 끄는 것이 아니다."""
    job = _job(tmp_path, dict(ROTATED, PIXSCALE=0.58, RLEVEL=0,
                              BIASSEC="[3100:3135,1:2048]"))
    data = make_lco_config.build(job)
    assert data["calibration"]["enabled"] is True
    assert data["calibration"]["overscan"]["enable"] is True
    assert data["calibration"]["overscan"]["width"] == 36


def test_no_rlevel_keyword_is_treated_as_raw(tmp_path):
    """RLEVEL 이 없는 기기는 예전처럼 원본으로 본다."""
    job = _job(tmp_path, dict(ROTATED, PIXSCALE=0.58,
                              BIASSEC="[3100:3135,1:2048]"))
    assert make_lco_config.build(job)["calibration"]["enabled"] is True
