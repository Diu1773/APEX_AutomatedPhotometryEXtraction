# PSF 엔진 대조 — IRAF DAOPHOT ALLSTAR (2026-08-12)

`APEX_ENGINE_SCORECARD.md` 의 A 축("기존 엔진과 실측 대조했는가")에서 step 8 은
✕ 다. 이유는 좁다 — APEX 에는 이미 2,591 줄짜리 IRAF 대조 하네스가 있지만
(`apex/benchmark/iraf_crosscheck*.py`) **그 하네스는 `phot`, 즉 IRAF 의 *조리개*
태스크만 부른다.** daophot 패키지를 올려놓고 `psf`·`allstar` 는 한 번도 쓰지
않는다. **Step 7 은 IRAF 와 맞대 봤고 Step 8 은 한 번도 없었다.**

이 디렉터리가 그 구멍을 메운다.

- `daophot_allstar.py` — WSL PyRAF 로 `phot` → `pstselect` → `psf` → `allstar`
  를 돌리고 평평한 표를 낸다
- `compare_recovery.py` — 같은 인공별 진리값에 두 엔진을 채점한다

## 왜 인공별인가

두 엔진 중 어느 쪽도 진리가 아니다. 서로 비교하면 **다르다**까지만 말할 수 있고
**누가 맞다**는 말할 수 없다. 알려진 밝기의 별을 실제 프레임에 심으면 각 엔진이
낸 오차가 의견이 아니라 사실이 된다. 주입은 APEX 의 기존 자산
(`apex/benchmark/psf_artificial_stars.py`)을 쓴다.

## 공정성 — 두 번 틀렸고 두 번 고쳤다

**첫째, `phot` 의 재중심을 켜서 DAOPHOT 을 불리하게 만들었다.** 좌표를 주는
강제측광인데 `calgorithm="centroid"` 를 걸었더니, 0.75 FWHM 떨어진 이웃이 있는
희미한 별의 중심이 이웃으로 끌려갔다. **주입별 25 개 중 7 개만 살아남는 동안
실제 별은 1,599 개 중 1,575 개가 통과**해서 발각됐다. `calgorithm="none"` 으로
고쳤다(ALLSTAR 는 적합 중에 여전히 위치를 미세조정한다 — APEX 가 그러듯이).
회수 8 → 11 개.

**둘째, 서로 다른 등급 체계를 빼고 있었다.** ALLSTAR 등급은 `phot` 이 쓴 구경의
영점을 물려받고 이 사슬은 구경보정을 걸지 않는다. APEX 의 ePSF 유량은 전체
유량이다. 밝고 한산한 별에서 잰 상수 오프셋이 **+0.296 mag** 으로 나왔는데
이는 1·FWHM 구경보정의 크기와 일치한다 — 즉 그 차이는 엔진이 아니라 구경이었다.
이제 엔진마다 상수를 제거한 뒤에만 bias 를 말한다(원값도 함께 기록).

**셋째, "회수"의 뜻이 서로 달랐다.** APEX 벤치마크는 주입 좌표를 step 7 강제
카탈로그에 넣으므로 step 8 은 모든 위치에서 값을 낸다 — `psf_recovered` 는
25/25 로 **정의상 1.00** 이다. ALLSTAR 는 못 믿을 별을 버린다. 이걸 완전도
비교로 적으면 거짓이 된다. `--apex-mode gated` 가 APEX 자신의 post-fit 정책
(`flags_psf==0`, `snr_psf≥3`, `qfit≤3`, `reduced_chi2≤25`)을 적용해 **버리는
엔진과 같은 질문**으로 맞춘다.

세 질문의 답이 다 다르다:

| 질문 | APEX |
|---|---|
| step 8 이 값을 냈나 | 25/25 |
| **자체 품질 임계를 통과했나** (ALLSTAR 대응) | **22/25** |
| 위치를 안 알려줘도 찾았나 | 7/25 |

세 번째는 **검출** 질문이고 ALLSTAR 에게는 묻지 않았다(좌표를 줬다).

## 방법 파라미터를 맞추지 않는 이유

검출기 상수(gain 0.68 e⁻/ADU · 읽기잡음 2.35 e⁻ · 유효구간 0.1–55,000 ADU ·
노출)는 **자료의 성질**이라 APEX 설정에서 그대로 가져와 두 엔진에 같이 준다.
반면 PSF 반경·적합 반경·해석함수는 **각 엔진 저자의 선택**이다. 이걸 APEX 에
맞추면 DAOPHOT 을 APEX 의 답 쪽으로 조율하는 것이 된다. 그래서 DAOPHOT 표준
지침(Massey & Davis 1992)을 따른다 — `psfrad = 4·FWHM+1`, `fitrad = FWHM`,
`function = auto`, `varorder = 0`.

## 명시할 교란요인

- **적합창이 다르다.** APEX 는 포위에너지 90 % 자동창(이 프레임에서 21 px,
  약 1.7·FWHM)이고 DAOPHOT 표준은 1·FWHM 이다. `--fitrad-fwhm 1.7` 로 같은
  사슬을 다시 돌려 **차이가 엔진 때문인지 창 때문인지 가려야 한다. 아직 안 했다.**
- **varorder = 0 으로 고정했다.** APEX 가 프레임당 ePSF 하나이므로 대응이지만,
  DAOPHOT 은 varorder = 1·2 로 시야 변화를 줄 수 있고 그 경우 더 나을 수 있다.
- **속도는 같은 조건이 아니다.** APEX 는 Windows 네이티브, DAOPHOT 은 WSL 경유다.
- **`psf` 태스크가 시간의 대부분을 쓴다**(535 s 중 449 s). 이 값은 기준성 60 개
  기준이며 `maxnpsf` 에 민감하다.

## 현재 상태 — 하네스는 돌고, 결론은 아직 없다

2026-08-12 기준 실행: M13 `pp_messier13-0005-B.fit` (Moravian C3-61000 ·
2026-05-15 · 60 s · FWHM 7.05 px), 좌표 1,599 개, 주입별 **25 개**.

**25 개로는 구간별 통계를 못 낸다.** 여러 칸의 산포가 별 2~3 개에서 나온다.
아래는 하네스가 끝까지 돈다는 증거이지 결과가 아니다.

| 엔진 | 유지 | bias(영점 제거) | 산포 |
|---|---:|---:|---:|
| APEX step8 (gated) | 22/25 (0.88) | +0.075 | 0.149 |
| IRAF ALLSTAR | 11/25 (0.44) | +0.007 | 0.047 |

읽히는 방향은 **"APEX 는 더 많이 답하고, ALLSTAR 는 답한 것에 대해 더
정확하다"** 이다. 이건 APEX 에 유리한 그림이 아니며, 그래서 더더욱 표본을
늘려서 확인해야 한다.

**다음**: 주입을 수백 개로 늘려(`--injections`) 밝기×혼잡 칸마다 통계가 서게
한 뒤 다시 채점한다. 그리고 `--fitrad-fwhm 1.7` 감도 확인.

## 재현

```bash
python validation/psf_engines/daophot_allstar.py \
  --frame <injected.fit> --positions <step7_photometry_*.tsv> \
  --output <daophot.csv> --fwhm <px>

python validation/psf_engines/compare_recovery.py \
  --truth <truth.csv> --apex-recovery <recovery.csv> --daophot <daophot.csv>
```
