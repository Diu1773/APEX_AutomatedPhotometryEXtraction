"""받기가 한 장에서 걸려도 나머지를 마저 받는가.

**수습하는 코드가 터지면 받기 전체가 죽는다.** 한 장이 실패하면 `download` 은
남은 `.part` 를 지우고 다음 장으로 넘어가게 돼 있는데, **그 지우는 줄이 스스로
예외를 낼 수 있다.** 윈도에서 다른 프로세스가 그 파일을 아직 붙들고 있으면
`unlink` 이 `PermissionError` 를 던지고, 그것은 이미 `except` 안이라 아무도 안
받는다. 그러면 남은 스물네 장을 못 받은 채 받기가 끝난다.

2026-09-14 에 실제로 그렇게 났다. BANZAI 프레임 서른 장 가운데 여섯 장을 받은
자리에서 죽었고, 로그 마지막 줄이 받기 실패가 아니라 `tmp.unlink` 의
`PermissionError` 였다.

**이 모듈의 원칙은 「한 장이 실패해도 나머지는 계속 받는다」이고**(같은 함수의
주석) 그 원칙을 깨뜨린 것이 원칙을 지키려던 줄이었다.
"""
from __future__ import annotations

import sys
import urllib.error
from pathlib import Path

import pytest

REPO = Path(__file__).absolute().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import fetch_lco  # noqa: E402


class _Resp:
    """읽으면 내용을 한 번 주고 끝나는 가짜 응답."""

    def __init__(self, body: bytes):
        self._body = body

    def read(self, _n):
        body, self._body = self._body, b""
        return body

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def _rows(n: int) -> list[dict]:
    return [{"filename": f"frame{i}.fits.fz", "url": f"http://example/{i}"}
            for i in range(n)]


def test_one_failed_frame_does_not_stop_the_rest(tmp_path, monkeypatch):
    """가운데 한 장이 실패해도 뒤의 장들이 받아진다."""
    monkeypatch.setattr(fetch_lco, "RETRIES", 1)
    monkeypatch.setattr(fetch_lco, "_sleep", lambda *_: None)

    def fake_open(req, timeout=None):
        if req.full_url.endswith("/1"):
            raise urllib.error.URLError("끊김")
        return _Resp(b"data")

    monkeypatch.setattr(fetch_lco.urllib.request, "urlopen", fake_open)
    fetch_lco.download(_rows(3), tmp_path)

    got = sorted(p.name for p in tmp_path.glob("*.fits.fz"))
    assert got == ["frame0.fits.fz", "frame2.fits.fz"]


def test_a_locked_part_file_does_not_stop_the_rest(tmp_path, monkeypatch):
    """**수습하다 터져도** 뒤의 장들이 받아진다.

    실패한 장의 `.part` 를 지우려는데 그 파일이 다른 프로세스에 잡혀 있는
    상황이다. 윈도에서 실제로 나는 일이고, 고치기 전에는 여기서 받기가 통째로
    끝났다.
    """
    monkeypatch.setattr(fetch_lco, "RETRIES", 1)
    monkeypatch.setattr(fetch_lco, "_sleep", lambda *_: None)

    def fake_open(req, timeout=None):
        if req.full_url.endswith("/1"):
            raise urllib.error.URLError("끊김")
        return _Resp(b"data")

    monkeypatch.setattr(fetch_lco.urllib.request, "urlopen", fake_open)

    real_unlink = Path.unlink

    def locked_unlink(self, missing_ok=False):
        if self.name == "frame1.fits.fz.part":
            raise PermissionError(
                32, "다른 프로세스가 파일을 사용 중이기 때문에 접근할 수 없습니다")
        return real_unlink(self, missing_ok=missing_ok)

    monkeypatch.setattr(Path, "unlink", locked_unlink)
    fetch_lco.download(_rows(3), tmp_path)

    got = sorted(p.name for p in tmp_path.glob("*.fits.fz"))
    assert got == ["frame0.fits.fz", "frame2.fits.fz"], (
        "잠긴 .part 를 지우려다 받기 전체가 멈췄다")


def test_already_downloaded_frames_are_skipped(tmp_path, monkeypatch):
    """이미 받은 장은 건너뛴다 — 끊긴 뒤 다시 돌릴 수 있어야 한다."""
    (tmp_path / "frame0.fits.fz").write_bytes(b"old")
    calls: list[str] = []

    def fake_open(req, timeout=None):
        calls.append(req.full_url)
        return _Resp(b"data")

    monkeypatch.setattr(fetch_lco.urllib.request, "urlopen", fake_open)
    fetch_lco.download(_rows(2), tmp_path)

    assert calls == ["http://example/1"]
    assert (tmp_path / "frame0.fits.fz").read_bytes() == b"old"


def test_download_returns_the_byte_count(tmp_path, monkeypatch):
    """받은 양을 돌려준다 — 부르는 쪽이 속도를 적는다."""
    monkeypatch.setattr(fetch_lco.urllib.request, "urlopen",
                        lambda req, timeout=None: _Resp(b"12345"))
    total = fetch_lco.download(_rows(2), tmp_path)
    assert total == 10


@pytest.mark.parametrize("name", ["frame0.fits.fz"])
def test_a_stale_part_from_a_crashed_run_is_replaced(tmp_path, monkeypatch, name):
    """앞서 죽은 실행이 남긴 `.part` 가 있어도 새로 받아 덮는다."""
    (tmp_path / f"{name}.part").write_bytes(b"leftover")
    monkeypatch.setattr(fetch_lco.urllib.request, "urlopen",
                        lambda req, timeout=None: _Resp(b"data"))
    fetch_lco.download(_rows(1), tmp_path)

    assert (tmp_path / name).read_bytes() == b"data"
    assert not (tmp_path / f"{name}.part").exists()
