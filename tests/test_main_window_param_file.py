"""파라미터 파일을 명시해서 연 창이 이전 프로젝트 상태에 밀리지 않는지.

`apex/.state/<mode>/project_state.json` 은 모드당 하나뿐이라 마지막에 GUI 로 연
프로젝트가 남는다. 그 상태를 그대로 복원하면 `--config` 로 지정한 워크스페이스
대신 엉뚱한 자료를 읽는다(실제로 검증 스크립트가 364장 대신 29장을 쟀다).

평소 GUI 실행(파라미터 파일을 안 넘기는 경우)에는 마지막에 고른 폴더를 되살리는
편의가 그대로여야 하므로, 두 방향을 모두 고정한다.
"""
from __future__ import annotations

from pathlib import Path
from types import SimpleNamespace

from apex.gui.main_window import MainWindowWorkflow


def _make_window(tmp_path: Path, *, explicit: bool, saved: dict) -> SimpleNamespace:
    config_data = tmp_path / "config_data"
    config_result = tmp_path / "config_result"
    for p in (config_data, config_result):
        p.mkdir(parents=True, exist_ok=True)

    file_manager = SimpleNamespace(
        ref_filename=None,
        set_multi_night_dirs=lambda *_a, **_k: None,
        clear_multi_night_dirs=lambda: None,
    )
    return SimpleNamespace(
        mode="lc",
        _explicit_param_file=explicit,
        params=SimpleNamespace(P=SimpleNamespace(
            data_dir=config_data,
            result_dir=config_result,
            cache_dir=config_result / "cache",
            filename_prefix="",
        )),
        project_state=SimpleNamespace(get_step_data=lambda _key: saved),
        file_manager=file_manager,
        _register_state_mirror=lambda: None,
    )


def test_explicit_param_file_beats_saved_project_state(tmp_path):
    """--config 로 연 경우 저장된 다른 프로젝트가 경로를 덮어쓰면 안 된다."""
    other = tmp_path / "other_project"
    (other / "data").mkdir(parents=True, exist_ok=True)
    window = _make_window(
        tmp_path,
        explicit=True,
        saved={
            "data_dir": str(other / "data"),
            "result_dir": str(other / "result"),
            "filename_prefix": "other_",
        },
    )

    MainWindowWorkflow._bootstrap_file_selection_state(
        window, respect_explicit_param=True
    )

    assert window.params.P.data_dir == tmp_path / "config_data"
    assert window.params.P.result_dir == tmp_path / "config_result"


def test_default_launch_still_restores_last_folder(tmp_path):
    """평소 실행(파라미터 파일 미지정)에서는 마지막에 고른 폴더를 되살린다."""
    saved_data = tmp_path / "saved" / "data"
    saved_result = tmp_path / "saved" / "result"
    saved_data.mkdir(parents=True, exist_ok=True)
    window = _make_window(
        tmp_path,
        explicit=False,
        saved={
            "data_dir": str(saved_data),
            "result_dir": str(saved_result),
        },
    )

    MainWindowWorkflow._bootstrap_file_selection_state(
        window, respect_explicit_param=True
    )

    assert window.params.P.data_dir == saved_data
    assert window.params.P.result_dir == saved_result


def test_loading_a_project_moves_params_even_when_explicit(tmp_path):
    """프로젝트 불러오기는 일부러 params 를 그 프로젝트로 옮기는 경로다.

    생성자에서만 respect_explicit_param 을 켜므로, 여기서는 명시 여부와 무관하게
    저장된 경로가 그대로 적용돼야 한다.
    """
    saved_data = tmp_path / "loaded" / "data"
    saved_result = tmp_path / "loaded" / "result"
    saved_data.mkdir(parents=True, exist_ok=True)
    window = _make_window(
        tmp_path,
        explicit=True,
        saved={
            "data_dir": str(saved_data),
            "result_dir": str(saved_result),
        },
    )

    MainWindowWorkflow._bootstrap_file_selection_state(window)   # 기본값 = False

    assert window.params.P.data_dir == saved_data
    assert window.params.P.result_dir == saved_result
