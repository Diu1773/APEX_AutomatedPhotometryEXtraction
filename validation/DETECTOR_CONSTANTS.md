# 검출기 상수 — 측정 근거와 재현 절차

**카메라** Moravian C3-61000 (Sony IMX455), 2×2 average binning, 4788×3194
**측정일** 2026-08-06 · **적용 커밋** `4edf2e0`

이 문서는 논문·설정·그림이 쓰는 gain 과 읽기잡음이 **어떤 코드로, 어떤 자료에서,
어떤 조건으로** 나왔는지를 남긴다. 이 값이 다시 필요해지면 아래 명령 하나로
재현되며, 재측정할 이유가 없다.

## 확정값

| 상수 | 값 | 비고 |
|---|---|---|
| gain (저장 화소) | **0.680 ± 0.008 e⁻/ADU** | 통계 ±0.0023, 필터·연도 계통 포함 ±0.008 |
| 읽기잡음 (저장 화소) | **2.35 e⁻** (3.4512 ADU) | bias 24 쌍, 산포 0.0064 ADU (0.2 %) |
| gain (광자 화소) | 0.170 e⁻/ADU | 저장 화소 값 ÷ 4 (2×2 average) |
| 읽기잡음 (광자 화소) | 1.17 e⁻ | 저장 화소 값 ÷ 2 |
| 암전류 | 0.0068 e⁻/s | 10–480 s 사다리, R² 0.997 |
| 헤더 `EGAIN` | 0.04952 | **측정값의 1/13.7 — 쓰면 안 된다** |

설정 필드: `[instrument] gain_e_per_adu = 0.68`, `rdnoise_e = 2.35`,
`noise_use_fits_header = false`.

## 어떤 코드로

`apex/analysis/detector_ptc.py` — Janesick 광자전달법. 2026-08-06 이전에는 이
계산이 외부 저장소(`AstralImage/core/camera_calib.py`)에 있어 논문 검증을
재현할 수 없었다. 자체 구현으로 옮기면서 **기존 구현과 7개 산출값이 비트 동일**
함을 확인했다(gain·읽기잡음·화소환산 2종·기울기·절편·쌍수, 최대 상대차 0).

- 읽기잡음: `RN = sqrt(var(bias_a − bias_b) / 2)`, 여러 쌍의 중앙값
- gain: `var((flat_a − flat_b)/2) = S/g + RN²` 의 기울기. 쌍 차분이 PRNU 와
  비네팅을 지우므로 flat 의 고정 패턴이 분산에 들어가지 않는다
- 모든 통계는 4σ 클리핑 — 우주선·핫픽셀·먼지가 분산에 들어가지 않는다

## 어떤 자료로

| | |
|---|---|
| bias | `E:\bias` 의 `bias-*.fit` 102 장 → 앞 24 쌍 |
| flat | 동일 기기 3,454 장에서 신호 20,000 ADU 이상을 5,000 ADU 구간별로 고르게 추림 → 181 장 → **46 쌍 채택** |
| flat 출처 | `E:\observe_raw_Analysis`, `E:\observe_DSY`, `E:\observed_Analysis` 의 flat 폴더 전부 |
| 기기 동일성 | `INSTRUME`·`XBINNING`·`NAXIS1/2` 로 확인. 1×1(20 장)·4×4(10 장)는 제외 |
| 신호 구간 | **20,334 – 56,192 ADU** (지렛대 64 %) |
| 필터 | B·V·R·g·r·i·z·Hα·OIII·SII 등 5 종이 채택 쌍에 등장 |
| 관측 기간 | 2024 – 2026 |
| 통계 영역 | 각 프레임 중앙 600×600 화소 |

## 어떤 조건으로

**절편을 읽기잡음² 로 고정한다.** 측광에 쓸 만큼 밝은 flat 은 신호가 원점에서
멀어 자유 적합의 절편(= RN²)을 구속하지 못한다. 실제로 20 k 이상만 쓰면 자유
절편이 **+1,178 ADU²** 로 발산한다(물리값 11.9). 읽기잡음은 bias 쌍에서 독립적으로
재므로 그 값을 절편에 넣고 기울기만 적합한다. 이렇게 하면 gain 오차가 4 배
좋아진다(±0.0141 → ±0.0035).

**포화 징후는 없다.** 상한 컷을 30 k/35 k/40 k/45 k/50 k/60 k 로 바꿔가며 적합하면
gain 이 0.670–0.679 사이에서 흔들릴 뿐 단조 증가하지 않는다. 포화가 시작되면
분산이 꺾여 gain 이 계속 커져야 하므로, 56 k 까지 선형 영역이라는 뜻이다.
광자 화소 gain 0.170 에서 16-bit ADC 가 11.1 ke⁻ 에 포화하는데 이는 제조사가
명시한 전 우물 51.4 ke⁻ 의 21.6 % 이므로, **물리적 우물이 차기 전에 ADC 가 먼저
잘려** 응답이 끝까지 선형이다.

**계통 산포가 통계오차보다 크다.**

| 갈래 | gain |
|---|---|
| 필터 g / i / r / Hα | 0.683 / 0.668 / 0.675 / 0.678 |
| 연도 2024 / 2025 / 2026 | 0.675 / 0.685 / 0.668 |

그래서 논문에 쓸 오차는 ±0.0023(통계)이 아니라 **±0.008**이다.

## 외부 정합

| | gain | 읽기잡음 |
|---|---|---|
| 이 측정 (광자 화소) | 0.170 e⁻/ADU | 1.17 e⁻ |
| Alarcón et al. 2023 (동일 IMX455, 고이득 모드) | — | 1.42 e⁻ |
| Moravian 공식 사양 | **공표하지 않음** | 최대 3.5 e⁻ |

제조사가 gain 을 공표하지 않는 이유는 읽기 모드마다 값이 달라지기 때문이다.
사양표에는 읽기잡음 상한과 전 우물만 있다. 우리 읽기잡음이 공식 상한 아래이고
같은 센서의 독립 측정과 같은 수준이며, ADC 포화가 전 우물의 21.6 % 라는 것까지
셋이 **고이득 읽기 모드** 하나로 일관되게 설명된다.

## 재현

```bash
.venv-deploy/Scripts/python.exe -X utf8 -m apex.benchmark.detector_characterize
```

기본 flat 폴더 두 개만 쓰므로 신호 구간이 20,238–23,790 ADU 로 좁고 gain 이
0.705 로 나온다 — **지렛대 부족 경고가 함께 뜬다.** 확정값을 재현하려면 flat
폴더를 넓게 준다.

```bash
.venv-deploy/Scripts/python.exe -X utf8 -m apex.benchmark.detector_characterize \
    --flats "E:\observe_raw_Analysis" "E:\observe_DSY" --signal-floor 20000
```

GUI 에서는 **Tools → Detector Characterisation (Gain / Read Noise)**. bias 폴더와
flat 폴더를 지정하면 같은 계산부를 돌리고, 헤더 `EGAIN` 과 10 % 넘게 어긋나면
경고하며, 측정값을 작업공간 설정에 바로 쓴다.

## 미해결

**채택 쌍 70 개 중 24 개가 버려진다** (분산이 0 이거나 비유한). 버려진 쌍의 신호는
27,415 – 64,342 ADU 에 퍼져 있고 중앙값이 51,011 이라 포화만으로는 설명되지 않는다.
같은 파일명이 여러 경로에 있는 경우(40 건)를 의심해 6 건을 열어 봤으나 내용이 모두
달라 중복도 아니었다. **원인 미규명.**

결과에는 영향이 없다고 판단한다 — 상한 컷을 30 k 까지 낮추면 탈락이 거의 없는데도
gain 이 0.679 로 같기 때문이다. 고신호 탈락이 편향을 만들었다면 이 값이 달라져야 한다.

## 이 값을 쓰는 곳

값을 바꿀 때 같이 고쳐야 하는 자리다.

| 자리 | 필드 |
|---|---|
| `parameters.toml` (재처리 씨앗 — `scripts/reprocess_batch.py` 가 읽는다) | `gain_e_per_adu`·`rdnoise_e`·`epadu`·`readnoise` |
| `apex_config.json` (레포 작업공간, 미추적) | 같음 |
| `parameters.example.toml` / `.json` (새 작업공간 템플릿) | 같음 |
| `parameters_M13.toml`·`parameters_M5.toml`·`parameters_capture_cmd.toml` | 같음 |
| `E:\APEX_validation\reprocess\<대상>\parameters.toml` 5 개 | 같음 |
| `validation/paper/fig1_completeness_snr.py`·`fig_completeness_realvssynth.py`·`make_injection_cutouts.py` | `GAIN` |
| `validation/analyze_m13_psf_asymmetry.py` | 나눗셈 상수 |
| `docs/manual/01-getting-started.md`·`06-parameters-reference.md` | 표의 기본값 |
| 그림 캡션 `fig_completeness_realvssynth.md`·`fig_qc_depth.md`·`fig11_detector.md` | 본문 수치 |

**GUI 의 IRAF 도구는 더 이상 하드코딩하지 않는다** — `IRAFParameters(instrument)`
가 작업공간의 `[instrument]` 에서 읽는다. IRAF 대조는 양쪽이 같은 상수를 받아야
의미가 있기 때문이다.
