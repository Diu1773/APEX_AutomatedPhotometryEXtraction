# Repository Guidelines

## Project Structure & Module Organization

APEX is a Python 3.10+ PyQt5 desktop application for astronomical photometry. The root `main.py` launches either mode. Mode entry points live in `apex/cmd/main.py` for cluster CMD work and `apex/lightcurve/main.py` for light-curve analysis. Shared application code is under `apex/`: `core/` contains state and file-management primitives, `config/` holds TOML-backed parameter models, `utils/` contains path, cache, I/O, logging, and astronomy helpers, `analysis/` contains science services, and `gui/` contains windows, workflow steps, tools, and widgets. SVG assets are in `apex/resources/`. Runtime settings are read from `parameters.toml`.

## Build, Test, and Development Commands

- `python3 main.py`: start the launcher and choose CMD or LC mode.
- `python3 apex/cmd/main.py`: run CMD mode directly.
- `python3 apex/lightcurve/main.py`: run LC mode directly.
- `python3 -m compileall apex main.py`: syntax-check the package after edits.
- `python3 -m pytest tests`: run tests once a `tests/` suite is added.

The README references `requirements.txt`, but this checkout does not currently include it. Keep dependency changes documented in `README.md` and include PyQt5, astropy, photutils, numpy, pandas, scipy, matplotlib, astroquery, and local astrometry.net requirements.

## Coding Style & Naming Conventions

Use 4-space indentation and standard Python naming: `snake_case` for functions, methods, and modules; `PascalCase` for classes; uppercase for constants such as `_CMD_MAIN`. Prefer `pathlib.Path` for filesystem work and keep workflow step files named by step and purpose, for example `step7_forced_aperture_phot.py`. Keep GUI changes consistent with the existing PyQt5 patterns in `apex/gui/main_window.py` and `apex/gui/workflow/step_window_base.py`.

## Testing Guidelines

There is no committed test suite yet. For new logic, add focused pytest tests under `tests/` using names like `test_step_paths_lc.py` or `test_period_analysis_service.py`. Prefer unit tests for `apex/utils/`, `apex/core/`, and `apex/analysis/`; GUI-heavy changes should at least pass `compileall` and include a short manual validation note.

## Commit & Pull Request Guidelines

Recent history uses concise lowercase prefixes such as `feat:`, `fix:`, and `remove:`. Keep commits scoped and imperative, for example `fix: preserve lc parameter cache path`. Pull requests should summarize the mode affected, list validation commands, note dependency or `parameters.toml` changes, and include screenshots for visible GUI updates.

---

# 하네스 프로토콜 (research-os)

> 제어층: `C:\Users\bmffr\Desktop\Main` — `NOW.md`, `PORTFOLIO.yaml`, `dashboard.html`

## 세션을 시작할 때

**TRACK.md (개발) · TRACK_PAPER.md (논문) 를 먼저 읽는다.** 이 트랙의 상태에 대한 유일한 진실이다.
사용자에게 "전에 뭘 했었죠"를 묻지 않는다. 파일에 다 있다.

**`## 사용자 의견` 절을 반드시 읽는다.** 사용자가 대시보드에서 남긴 방향 지시가
거기 최신순으로 쌓인다. 그 방향과 다르게 가려면 먼저 물어본다.


## 세션을 끝내기 전에 (필수)

1. **오라클 실행**
   ```bash
   .venv-deploy/Scripts/python.exe -m pytest tests/ -q
   ```
   통과 기준: 614 passed, 0 failed (약 7분). 실패하면 완료가 아니다 — 롤백하거나 원인을 TRACK.md에 남긴다.
2. **커밋 + 푸쉬.** 미푸쉬 커밋을 남기고 끝내지 않는다.
3. **TRACK.md 갱신** — `## 지금` / `## 다음 3개` / `## 함정` 세 절을 다시 쓴다.
   판단이 필요한 건 `## 사용자 판단 필요`에. 채팅에만 쓰면 세션과 함께 사라진다.

Stop hook (`Main/scripts/hook_track_freshness.py`) 이 이걸 검사한다. 차단하지는 않지만 경고한다.

## Codex 역할 — 검토자

이 하네스에서 Codex는 **구현자가 아니라 검토자**다. Claude가 만든 변경을 독립적으로 본다.
작성자와 검토자가 다른 공급자면 같은 오류를 공유할 확률이 줄고, Claude 사용 한도도 아낀다.

검토할 때 보는 것:

- 오라클이 실제로 돌았는가, 통과 기준을 만족했는가
- 회귀 — 기존 동작이 조용히 바뀌지 않았는가
- 수치 · 단위 · 시간대 · 좌표계 일관성
- 시계열 데이터라면 **미래 정보 누수**
- 하드코딩된 가변 값 (기기명 · 필터 · 경로 · 스텝값)

검토 결과는 PR 코멘트 또는 TRACK.md의 `## 함정`에 남긴다. **구현하지 말고 문제만 보고한다.**

## 이 레포의 함정

- `validation/` 은 20 GB다. 절대 커밋하지 않는다.
- 실행은 `run.bat` 또는 `.venv-deploy\Scripts\python`. 시스템 python 아님.
- `step10_zeropoint_calibration.py` 의 pandas `DataFrameGroupBy.apply` 는 폐기 예정이다.
  pandas를 올리면 영점보정 결과가 조용히 달라질 수 있다.
