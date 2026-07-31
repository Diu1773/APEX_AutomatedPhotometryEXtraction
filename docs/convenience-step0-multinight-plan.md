# 편의성 개선 계획 — Step 0 자동 분류·매칭 + 멀티나잇 전반

> 2026-07-31 구조 검토 세션에서 도출. 절차: **이 계획 → 2차 검토(사용자) → Opus 구현**.
> 구현 세션은 이 문서만 읽고 착수할 수 있어야 한다. 검토 대상 코드:
> `apex/analysis/calibration_scan.py`, `apex/gui/workflow/step0_detector_calibration.py`,
> `apex/gui/workflow/lc/step1_night_setup.py`, `apex/core/file_manager.py`,
> `apex/analysis/merge/{workspace_scan,id_match,workspace_build}.py`,
> `apex/gui/tools/multi_night_merger.py`.

## 사용자 지시 (2026-07-31, 확정 방향)

1. **flat 매칭의 "하루"는 자정이 아니라 해가 떠 있는 시간대를 경계로 잘라야 한다**
   — 저녁 flat과 새벽(자정 이후) flat·light 가 같은 관측밤으로 묶여야 함.
2. **dark 온도 매칭은 0.1 °C 차이도 크게 보는 사용자가 있다** — 허용오차를 고민할 것.
3. **디스크 절약은 좋으나, FITS 같은 무거운 파일 제어가 우선이다.**
   `storage_mode` 필드는 현재 `"full"` 만 기록되는 자리표시자(다른 모드 미구현)임을 확인함.
   사실관계: 머저는 **FITS 를 복사하지 않는다**(`file_path_map.json` 경로 참조만).
   복사되는 것은 프레임당 photometry TSV 전량이다.

---

## P1. 관측밤(night) 판정 통일 — 일출·일몰 기준 [사용자 지시 1]

**문제**
- Step 0: `night_from_path(경로 8자리 날짜) or night_from_dateobs(정오분할)` 로
  **경로 날짜가 우선** (`calibration_scan.py:181`). 자정 이후 촬영분 파일명에
  다음날 타임스탬프가 박히는 캡처 SW(설정에 따라 MaxIm 포함)에서는 한 밤이 두
  night 로 찢어지고 dark/flat 매칭이 갈라진다.
- LC Step 1 은 별도 정의(JD 간격 8h, `step1_night_setup.py:_classify_nights_by_jd_gap`)
  — 같은 앱에 night 정의가 두 벌.

**수정안**
- `apex/utils/night_utils.py` 신설: `observing_night(date_obs_utc, lon_east_deg) -> "YYYYMMDD"`.
  DATE-OBS(UTC)를 경도 기반 지역 태양시(UTC + lon/15h)로 옮긴 뒤 **지역 정오에서 분할**.
  정오는 (극지 제외) 항상 일출~일몰 사이이므로 "낮에 자른다"는 요구를
  위도 정보·추가 의존성 없이 만족한다. 실제 태양고도 계산은 불필요 —
  낮 시간대 어느 시각에 잘라도 밤 묶임 결과가 동일하다. (2차 검토 질문 Q1)
- Step 0: 우선순위 반전 — **DATE-OBS 정오분할 우선**, 경로 날짜는 DATE-OBS 부재 시
  폴백. 둘 다 있고 값이 다르면 스캔 로그에 경고 1줄.
- LC Step 1: 공용 유틸 기반 정오분할을 기본 분류로, JD-gap 은 DATE-OBS/JD 없는
  프레임 폴백으로 강등 (2차 검토 질문 Q2). night_id 는 지금처럼 1-based 연속
  정수로 재매핑 → **산출물 스키마 불변**. `night_gap_hours` UI/파라미터는 유지.
- JD 도 DATE-OBS 도 없는 프레임: 현재 조용히 Night 1 (`id_map.get(i, 1)`) →
  **night 0 "불명"** 으로 분리하고 night 표에 경고 행 표시. 기본 제외 여부는 Q3.
- `night_assignments.json` 에 night_id → 관측일(YYYYMMDD) 매핑을 함께 저장
  (Step 10·머저 로그에서 "Night 3 = 2026-05-14" 표기 가능).

**검증** — 기존 검증 데이터(M13·NGC6811·YZ Boo, `E:\observed_Analysis`)에서
before/after night 배정 동일성 확인(이 폴더 구조에선 결과가 안 바뀌어야 정상).
자정 넘김 파일명 시나리오는 합성 헤더가 아닌 실제 헤더 사본으로 단위테스트.

## P2. Dark 매칭 허용오차 명시화 [사용자 지시 2]

**문제** — 현재 int 반올림 버킷(1 °C 폭) + nearest 매칭. (a) 버킷 경계에서
−10.4 와 −10.6 이 다른 그룹, (b) ΔT 가 아무리 커도 조용히 매칭, (c) 노출도
nearest 폴백인데 `dark_scale=False` 면 10s light 에 300s dark 가 무경고 진입.

**수정안**
- 그룹핑(트리 표시)은 버킷 유지, **매칭은 실제 온도 nearest + `temp_match_tol_c`
  게이트**로 변경. 기본 1.0 °C, 범위 0.1~10.0 (0.1 °C 민감 사용자는 내려서 쓰게).
- tol 초과 시: 기본은 경고 로그 + 트리 매칭 미리보기에 `ΔT=2.3°C` 명시.
  `strict_temp=True`(opt-in)면 매칭 거부 → 해당 light 는 dark 스킵 경고.
- `calibration.json` 의 프레임 레코드에 매칭 ΔT·Δexp 를 기록(provenance).
- 노출 nearest 폴백 + `dark_scale=False` 조합일 때 경고 1줄 추가.

## P3. `[calibration]` TOML 연결 + 설정 영속화

**문제** — `parameters.example.toml` 의 `[calibration]` 섹션을 **파싱하는 코드가
없다**(죽은 섹션). GUI 설정은 창 인스턴스에만 살고 닫으면 초기화(전역 철칙 11 위반).
헤드리스 `pipeline/steps/calibration.py` 는 `not_implemented` 스텁.

**수정안**
- `apex/config/schema.py` 에 CalibrationConfig 추가, `CalibrationOptions` 필드와
  1:1 (+ P2 의 `temp_match_tol_c`, `strict_temp`). parameters_cmd/lc 양쪽 로드.
- Step 0 GUI: 초기값을 TOML 에서 읽고, Parameters 저장 시 param 파일에 되쓰기
  (Step 1 의 `_persist_param_file` 패턴 재사용).
- 헤드리스 스텝 구현: scan → match → run 을 GUI `_CalibrationWorker._run` 과
  **같은 코어 함수**로 (GUI-코어 동일경로 원칙). 배포 로드맵(Phase 1 이식)과 합류.

## P4. 미분류 프레임 가시화 + 수동 재분류

**문제** — `read_frame_info` 가 분류 실패 시 None 반환 → 파일이 **소리 없이
스캔 결과에서 빠진다**. 오분류·미분류를 GUI 에서 고칠 방법이 전혀 없어서 헤더가
비표준인 데이터를 만나면 Step 0 자체를 못 쓴다.

**수정안**
- 분류 실패 프레임을 `ftype="unknown"` 으로 살려 반환, 트리에 "Unclassified (N)"
  노드로 표시(Run 대상에서는 제외).
- 트리 우클릭 → 타입(bias/dark/flat/light) 재지정, 필터·night 수정.
- 오버라이드는 `step0_calibration/classification_overrides.json` 에 영속
  (같은 폴더 재스캔 시 자동 적용 — 철칙 11).

## P5. 머저 위치 매칭 greedy → 전역 최근접 배정

**문제** — `id_match.py` 의 positional 매칭이 **행 순서대로** 최근접 canonical 을
선점. 크라우딩(M5/M13급)에서 행 A 가 1.8″ 별을 먼저 가져가면 0.3″ 진짜 주인
행 B 는 신규 별로 **중복 생성**되고, 선점당한 행은 2등 후보를 시도하지 않는다.

**수정안** — (row, canonical) 후보쌍을 separation 오름차순 정렬 후 일대일 배정
(전역 greedy = 사실상 mutual nearest). 수정 전후를 `merge_id_map.csv` 의
method/sep 분포와 신규 생성 수로 회귀 비교.

## P6. 머저 입력 게이트 완화

**문제** — `workspace_scan.py:merge_ready` 가 폴더마다 Step 9 `lightcurve_*.csv`
를 요구하지만, materialize 는 Step 7 TSV + Step 8 카탈로그만 쓰고 LC 는 merged
workspace 에서 다시 만든다. **밤마다 Step 9 완주를 강제하는 헛수고.**

**수정안** — 필수 게이트는 Step 7 + Step 8 로 낮추고, Step 9 부재는 폴더 스캔
표에 경고 표시로 강등.

## P7. Merged workspace `storage_mode="reference"` [사용자 지시 3]

**현황** — FITS 는 이미 참조만. 복사되는 건 프레임당 photometry TSV 전량
(3000프레임 멀티나잇이면 원본 총합만큼 중복). `storage_mode` 는 `"full"` 고정
자리표시자. full 모드의 존재 이유는 **이식성**(merged 폴더 단독으로 완결,
원본 폴더를 옮기거나 지워도 동작)으로 추정된다.

**수정안**
- `reference` 모드 추가(opt-in, 기본은 full 유지 — Q4): per-frame TSV 복사 생략.
  매니페스트에 폴더별 원본 경로 + ID 리맵 테이블(`merge_id_map.csv` 활용)을 남기고,
  `photometry_loader.load_frame_photometry` 가 merged workspace 에서는 원본 폴더
  TSV 를 읽어 ID 를 on-the-fly 리맵.
- **주의(구현 전 필수)**: 로더를 우회해 `step7_forced_phot/` 를 직접 glob 하는
  코드가 있으면 깨진다. Step 9/10/11 + `_preload_photometry_cache` 경로 전수 확인.
- 트레이드오프 명시: reference 모드는 원본 폴더 이동·삭제 시 깨짐 — UI 에 고지.
- 규모 中. P1~P6 과 독립적이라 별도 커밋 단위로.

## P8. 소규모 수정 묶음 (반나절)

| # | 내용 | 위치 |
|---|------|------|
| a | flat 매칭을 `normalize_filter_key()` 경유로 (지금은 소문자 비교만 — "Ha"/"HA" 오탐) | `calibration_scan.py:match_flat` |
| b | `FrameInfo` 에 비닝/shape 추가, dark·flat 매칭 키에 포함, 불일치 트리 경고 (2×2 dark ↔ 1×1 light 방지) | `calibration_scan.py` |
| c | Step 0 완료 후 `data_dir` 만 바꾸고 `filename_prefix` 방치 → `pp_` 산출물과 prefix 불일치로 Step 1 에서 0개 뜨는 문제. prefix 도 함께 갱신/초기화 | `step0_detector_calibration.py:_set_data_dir` |
| d | Step 0 재실행 시 출력 폴더의 기존 `calibration.json` 감지 → 이어하기/재계산 선택 | `step0_detector_calibration.py:_on_run` |
| e | 겹치는 관측일을 가진 workspace 끼리 merge 시 경고 (현재는 같은 밤이 merged night 2개로 갈라짐) | `workspace_build.py` |
| f | selection 기본값을 base 폴더에서만 가져옴 → base→나머지 우선순위 union | `multi_night_merger.py:compute_selection_defaults` |
| g | 머저 창 수제 스타일(`#1565C0` 버튼, `#E3F2FD` 박스, Arial) → `style_button`/`Tokens` 로 교체 (다크 테마 대응, UI 백로그와 합류) | `multi_night_merger.py` |

---

## 구현 순서 권고

**P1 → P2 → P8(a–d) → P3 → P4 → P5 → P6 → P8(e–g) → P7**

P1·P2 가 사용자 지시 직결 + 과학적 정합성. P7 은 독립 트랙이라 마지막.
각 P 는 커밋 단위로 끊고, P1·P5 는 실데이터 before/after 비교를 커밋 메시지에 남긴다.

## 검증 (전 단계 공통)

```bash
.venv-deploy/Scripts/python.exe -m pytest tests/ -q   # 오라클: 625 passed 기준
```

- P1: M13·NGC6811·YZ Boo 실폴더 night 배정 before/after 동일성 (달라지면 원인 규명 후 진행)
- P2: ΔT·Δexp 가 `calibration.json` 에 남는지 + strict 모드 거부 동작 단위테스트
- P5: 실측 M5/M13 카탈로그로 merge 재실행 → 신규 생성 수 감소·sep 분포 개선 확인
- P7: full vs reference 로 만든 merged workspace 의 Step 9 산출물 비트 동일성

## 2차 검토 질문 (사용자 판단 필요)

- **Q1** 밤 경계: 지역 정오분할로 충분한가, 실제 일출·일몰(태양고도) 계산까지 갈
  것인가. **권고: 정오분할** — 밤 묶임 결과가 동일하고 의존성·위도 입력이 불필요.
- **Q2** LC Step 1 에서 JD-gap 을 폴백으로 강등할지, 정오분할과 병행 표시할지.
  **권고: 폴백 강등** (정의 단일화).
- **Q3** JD 불명 프레임(night 0)의 기본 처리: 제외 vs 포함+경고. **권고: 제외 기본.**
- **Q4** reference 모드 기본값: **권고 full 유지, reference 는 opt-in** (이식성 보존).
- **Q5** `temp_match_tol_c` 기본 1.0 °C + strict off 가 적절한가.
