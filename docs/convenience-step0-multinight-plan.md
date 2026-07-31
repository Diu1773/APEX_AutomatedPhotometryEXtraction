# 편의성 개선 계획 — Step 0 자동 분류·매칭 + 멀티나잇 전반

> 2026-07-31 구조 검토 세션에서 도출. 절차: 계획 → 2차 검토 → Opus 구현.
> **상태: P1~P8 전부 구현·검증 완료 (2026-07-31). 오라클 689 passed, 0 failed.**
> 구현 결과와 실데이터에서 드러난 사실은 아래 「구현 결과」 절에 있다.
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
  낮 시간대 어느 시각에 잘라도 밤 묶임 결과가 동일하다. (Q1 — 2차 검토로 확정)
- **[2차 검토 수정] 경도/시간대 정보가 없으면 DATE-OBS 우선이 오히려 위험하다.**
  `night_from_dateobs` 의 경도-부재 폴백은 UTC −12h 규칙(그리니치 기준)이라
  동아시아에서 한 밤을 찢는다(현재 코드 주석에도 명시). 확정 우선순위:
  1. DATE-OBS + 경도(헤더 SITELONG 등 `_LON_KEYS`) → 태양시 정오분할
  2. DATE-OBS + `params.P.site_tz_offset_hours`(TOML `[site] tz_offset_hours`,
     GUI 설정에 이미 존재) → 시민시 정오분할. **Step 0 창도 params 접근이
     가능하므로 이 폴백을 받는다** (현재 Step 0 는 헤더만 봄 — 갭).
  3. 경로 8자리 날짜 (경도·tz 둘 다 없을 때만)
  4. DATE-OBS −12h (최후 폴백)
  1·2 가 가능한데 경로 날짜와 다르면 스캔 로그에 경고 1줄.
- LC Step 1: 공용 유틸 기반 정오분할을 기본 분류로, JD-gap 은 DATE-OBS/JD 없는
  프레임 폴백으로 강등 (Q2). 경도는 헤더 테이블에 안 남으므로(2차 검토 확인)
  LC 쪽은 `site_tz_offset_hours` 를 쓴다. night_id 는 지금처럼 1-based 연속
  정수로 재매핑 → **산출물 스키마 불변**. `night_gap_hours` UI/파라미터는 유지.
  ⚠️ `site_tz_offset_hours` 기본 0.0 인 채 정오분할하면 폴백 4 와 같아지므로,
  tz 미설정 + 경도 부재면 **JD-gap 을 그대로 쓴다** (조용한 오분류 방지).
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
  1:1 (+ P2 의 `temp_match_tol_c`, `strict_temp`). parameters_cmd/lc 양쪽 로드
  (`parameter_map.py` 의 (("calibration", …), "…") 매핑 패턴 따름).
- Step 0 GUI: 초기값을 TOML 에서 읽고, Parameters 저장 시 param 파일에 되쓰기.
  **[2차 검토 확인]** `_persist_param_file` 은 Step 1 창
  (`step1_file_selection_common.py:554`)에만 있고 Step 0 의 `ToolWindowBase`
  에는 없다 — 유틸로 추출해 공용화하거나 Step 0 에 동등 메서드를 둔다.
- **[2차 검토 발견] 헤드리스 Step 0 는 이미 실전 가동 중이다** —
  `scripts/_reprocess_step0.py` 가 offscreen QApplication 으로 GUI 워커
  (`_ScanWorker`/`_CalibrationWorker`)를 직접 구동한다(패리티 목적).
  올바른 이식: `_CalibrationWorker._run` 본체를 Qt-free
  `apex/analysis/calibration_run.py` (progress 콜백 인자)로 추출하고
  **GUI 워커·pipeline 스텝·재처리 스크립트 3자가 그 한 함수를 호출** —
  GUI-코어 동일경로 원칙을 지키면서 QApplication 의존을 없앤다.

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

**수정안 — [2차 검토로 설계 확정]**
- `reference` 모드 추가(opt-in, 기본은 full 유지 — Q4): per-frame TSV 복사 생략.
- **로더 중앙화가 성립한다** (호출부 전수조사 완료): TSV 경로 결정은
  `photometry_loader._resolve_photometry_path` 단일 지점이고, LC 소비자 전원
  (step8 selection·step9 builder·step10 detrend 의 `photometry_source_service`·
  `photometric_qc`·extinction_fit·qa_report)이 `load_frame_photometry` 를 경유한다.
  **step9 builder 는 `photometry_index.csv` 의 `path` 컬럼을 쓰지 않음**을 확인.
- 구현: merged workspace 에서 파일이 없으면 `merge_manifest.json` 의
  folder_tag→원본 result_dir 매핑으로 `F01_xxx__원본명` 을 파싱해 원본 TSV 를
  읽고, `merge_id_map.csv` 의 (folder_tag, filter, local source_id →
  merged source_id) 리맵을 적용(모듈 캐시). 이후의 sid_map(merged 카탈로그
  기반 ID 부여)은 기존 코드 그대로 동작한다.
- 알려진 열화: `qa_report.py` 는 `photometry_*.tsv` 를 직접 glob 하는 경로가
  일부 있어 reference workspace 에서 빈 결과가 된다 — 해당 도구에
  "reference merged workspace 미지원" 안내 1줄. CMD 전용 glob 은 무관.
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

## 2차 검토 결과 (2026-07-31) — Q1~Q5 확정

2차 검토(코드 대조 전수 확인)를 거쳐 권고안대로 확정한다. 구현 세션은 아래를
기본값으로 진행하고, 사용자가 반대 의견을 주면 그때 수정한다.

- **Q1 확정: 지역 정오분할.** 태양고도 계산 불필요 — 낮 어느 시각에 잘라도
  결과 동일. 단, 경도·tz 둘 다 없으면 정오분할 자체를 포기하는 폴백 체인이
  필수 (P1 의 [2차 검토 수정] 항목 — UTC −12h 는 동아시아에서 밤을 찢는다).
- **Q2 확정: JD-gap 폴백 강등.** 단 tz 미설정 + 경도 부재 환경에서는 JD-gap 유지.
- **Q3 확정: night 0(시각 불명) 기본 제외 + 표 경고행.**
- **Q4 확정: full 기본 유지, reference 는 opt-in.** full 의 존재 이유(이식성)
  확인. 로더 중앙화 성립도 확인(P7) — 구현 리스크 낮음.
- **Q5 확정: `temp_match_tol_c` 기본 1.0 °C, `strict_temp` 기본 off.**
  0.1 °C 민감 사용자는 tol 을 내리거나 strict 를 켠다 (GUI decimals=1, 하한 0.1).

### 2차 검토에서 나온 계획 수정·보강 (구현 시 필독)

1. **P1 우선순위 체인 수정** — "무조건 DATE-OBS 우선" 은 위험. 경도/tz 가용성에
   따른 4단 폴백으로 확정 (P1 본문). Step 0 에 `site_tz_offset_hours` 폴백 추가.
2. **P3 경로 구체화** — 헤드리스 Step 0 는 `scripts/_reprocess_step0.py` 로 이미
   가동 중(GUI 워커를 offscreen 구동). `_CalibrationWorker._run` 을 Qt-free 로
   추출해 GUI·pipeline·스크립트 3자 단일 경로화. `_persist_param_file` 은
   Step 1 전용이라 공용 유틸로 추출 필요.
3. **P7 설계 확정** — 로더 단일 지점(`_resolve_photometry_path`) 확인,
   step9 는 index 의 `path` 컬럼 미사용, qa_report 직접 glob 만 알려진 열화.
4. **P8c 참고** — prefix 초기화 패턴은 `step1_file_selection.py:485`
   (`_persist_param_file(io_updates={"filename_prefix": ""})`) 재사용.
5. **LC 헤더 테이블에는 SITELONG 이 안 남는다** — LC 정오분할의 경도 출처는
   params 뿐. 헤더 스캔에 경도 키를 추가로 보존하는 확장은 선택 사항(비채택 —
   tz 파라미터로 충분).

---

# 구현 결과 (2026-07-31)

전 항목 구현·검증 완료. 오라클 **689 passed, 0 failed** (8분 7초, 기존 625 → +64).

| # | 커밋 | 내용 | 실데이터 영향 |
|---|------|------|---------------|
| P1 | `451e659` | 관측밤 판정 통일 (`apex/utils/night_utils.py`) | YZ Boo 128/879장 재배정 — 전부 설명됨, 보정 결과 불변 |
| P2 | `ebad05b` | dark 온도 허용오차 설정화 | **YZ Boo 283장이 ΔT 5.0 °C — 숨겨져 있던 결함 발견** |
| P8a–d | `c56da6b` | 필터 정규화·비닝 검사·prefix·재실행 확인 | 매칭 손실 0건 |
| P3 | `bffaac4` | `[calibration]` TOML 연결 + Qt-free 코어 | M13 15/15 프레임 비트 동일 |
| P4 | `3d7d94a` | 미분류 프레임 가시화·수동 재분류 | GUI 비포/애프터 확인 |
| P5 | `37f3617` | 머저 위치매칭 전역 최근접 배정 | **영향 없음** (실측 최근접 이웃 3.97" > 허용 2") |
| P6 | `6e6f29e` | Step 9 완주 강제 제거 | 밤마다 Step 9 불필요 |
| P8e–g | `dcc8ef3` | 밤 중복 경고·selection union·테마 | 다크 테마 캡처 확인 |
| P7 | `b69a974` | merged workspace 원본참조 모드 | full/reference 로드 결과 동일 |

## 실데이터에서 새로 드러난 것 (후속 판단 필요)

### YZ Boo 의 dark 가 관측 조건과 안 맞는다 — P2 가 드러냄

`E:\observe_raw_Analysis\YZbootis` 전체에 dark 가 **30 s / −10 °C 한 세트
(2025-04-28)** 뿐이다. 기존 UI 는 어느 프레임에나 `dark 30s/-10C` 로만
표시해서 이 사실이 보이지 않았다.

- **night 20250430 의 283장은 −5 °C 에서 촬영** → ΔT **5.0 °C** (최대 5.29).
  암전류는 약 6 °C 마다 2배이므로 빼는 dark 가 1.8배가량 과소하다.
  이 283장에는 **TRACK 에 기록된 LC 완주용 g 필터 97장이 포함**된다.
- 2026 nights 는 40/50/60/120 s 노출인데 dark 는 30 s → Δexp 최대 90 s.
  `dark_scale` 기본값이 켜져 있어 비율 보정되지만 30→120 s 는 4배 외삽이다.
- 2026 년 light 에 2025-04-28 dark 를 쓰고 있다 (약 11개월 차이).

차등측광이 상수 오프셋을 상당 부분 상쇄하므로 **변광 검출 자체가 뒤집힐
사안은 아니지만**, 논문에 YZ Boo 를 쓸 경우 이 조건은 명시해야 한다.
지금은 실행할 때마다 로그·`calibration.json`·트리에 ΔT 가 남는다.

**원인을 추적하니 자료가 아니라 코드였다 (2026-07-31, `8cd2f7f`).**
`E:\darks\dark-30s-5` 에 30 s/−5 °C dark 10장(2024-10-04)이 **있었다**.
`group_for_night` 이 같은 밤 dark 가 하나라도 있으면 전역 풀을 통째로
배제해서, YZ Boo 폴더의 30 s/−10 °C dark 가 라이브러리를 막고 있었다.
수정 후 같은 조건에서 ΔT 5.01 → **0.16 °C**. 20260325 의 120 s light 31장도
40 s dark(Δexp 80 s) 대신 120 s dark(Δexp 0)를 쓴다.

**실측한 비용** — 벌크 암전류는 무시 가능하다(30 s 에 median 0.1 ADU
= 0.07 e-, read noise 2.1 e-). 문제는 **hot pixel 마스크**다. APEX 는 이걸
마스터 다크에서 만드는데, −10 °C 다크로 만든 마스크가 −5 °C 에서 실제로
hot 인 픽셀의 **29.4%(34,022 px)를 놓친다**(median 11.3 · p90 16.5 ·
max 899 ADU). 온도가 노후화보다 중요하다: 같은 온도 6개월 차 마스터끼리
차이 std 7.29 ADU, 5 °C 차 1개월 차는 11.03 ADU — **온도 맞은 오래된 다크가
온도 틀린 최신 다크보다 낫다.**

**온도 보간은 권하지 않는다.** (a) 보간할 대상인 벌크 암전류가 이미
무시 가능하고, (b) 실제로 중요한 hot pixel 은 개별 결함(활성화 에너지가
제각각, RTS 픽셀은 무작위로 튄다)이라 매끄러운 보간 법칙이 가장 못 맞히는
부분이다. E: 라이브러리가 −15/−10/−5/0/+5 °C 를 덮으므로 최근접 실측
다크로도 최대 2.5 °C 안에 든다.

**사용자 판단 필요**: YZ Boo 를 새 코드로 재처리할지 (TRACK D-006).

### 온도 허용오차 기본값의 근거

정상 자료(M13·NGC6811)의 ΔT 는 중앙값 0.05~0.11 °C, 최대 0.23 °C 다.
0.1 °C 허용오차는 **센서 온도조절 잡음 수준**이라 789장 중 398장이 플래그된다.
기본 1.0 °C 가 타당하고, 0.1 °C 까지 내리는 것은 사용자 선택으로 남겼다.

### P5 는 현재 자료를 개선하지 않는다

M13(1,347개)·M5(1,406개) 실측 마스터 카탈로그의 최근접 이웃 거리는
**최소 3.97" · 중앙값 10.4"**, 3" 안쪽 쌍이 0% 다. 기본 허용오차 2" 의
두 배 이상이라 경합이 일어나지 않는다. UI 최대치 5" + 2" 지터로 밀어붙여도
매칭·신규 건수가 같고 평균 separation 만 0.014~0.020" 줄었다.
**도달 가능한 결함을 막는 방어적 수정**이지 개선이 아니다.

## 검증 방법 (재현용)

- night 배정 before/after · dark ΔT 분포 · geometry 필터 손실 · P5 old/new 비교:
  세션 중 스크래치 스크립트로 수행 (실데이터 `E:\observe_raw_Analysis`,
  `E:\observed_Analysis`). 결론은 위 표와 커밋 메시지에 기록됨
- P3 패리티: M13 15장을 새 코어로 재보정해
  `E:\APEX_validation\reprocess\M13\sci` (옛 GUI 워커 산출물)와 픽셀 비교
- P4·P8g: 실제 창을 `apex-light`/`midnight` 테마로 띄워 캡처
