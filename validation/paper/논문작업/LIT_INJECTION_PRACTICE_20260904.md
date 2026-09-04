# 인공별 주입을 남들은 무엇에 쓰는가 — 진행 중

**묻는 것**: 측광 소프트웨어 논문들은 「값이 옳다」를 무엇으로 보이는가.
특히 인공별 주입(artificial star injection)을 **「찾았나」(완전도·한계등급)**에만
쓰는가, **「제대로 쟀나」(등급 정확도)**에도 쓰는가.

이 답에 따라 APEX 가 뽑은 등급 회수 편차 표(11~18 등급에서 ±11 mmag)의 자리가
정해진다 — 표준을 따르는 것인지, 표준을 넘는 것인지, 표준이 요구하는 것을
빠뜨린 것인지.

**진행 방식**: 한 편씩 원문을 열어 확인하고 그때그때 여기 적는다.
(2026-09-04 에 22 개 병렬로 돌리다 전량 손실 — `Main/FAILURES.md` F-251)

---

## 1. AutoPhOT — Brennan & Fraser (2022), A&A 667, A62

**확인 경로**: arXiv:2201.02635 전문(ar5iv). A&A 사이트는 403 으로 막힘.
**APEX 와의 관계**: 가장 가까운 유례. 자동 측광 파이프라인, 최근, 소프트웨어 논문.

검증 절은 §7 Testing and validation (§7.1 Testing of photometry packages,
§7.2 Performance).

**인공별 주입 → 한계등급에만 쓴다.**

> *"Finally, the most rigorous limiting magnitude is determined though injecting
> and recovering artificial sources."*

한계등급의 정의도 회수 여부로만 되어 있다 — 개별 측정이 β < 0.75 면 잃은 것으로
보고, **80 % 가 잃어지는 등급**을 한계로 잡는다. 즉 「찾았나 못 찾았나」다.

**등급 정확도 → DAOPHOT 대조로 보인다.**

> *"AutoPhOT's ability to recover source fluxes is consistent with commonly used
> software e.g. DAOPHOT, using both aperture and PSF photometry."*

> *"The PSF fitting package from AutoPhOT can match the recovered instrumental
> magnitude from DAOPHOT, even at faint magnitudes where the flux from the source
> becomes comparable to the sky background."*

Figure 15 가 그 비교인데, **본문에 mmag 수치가 명시되어 있지 않다.**

**외부 카탈로그 대조는 안 한다.** 대신 각 이미지 안의 sequence star 로 영점을
맞춘다.

### 이 한 편이 말하는 것

주입한 별의 **등급을 되찾는 정확도는 보고하지 않는다.** 주입은 깊이를 재는
도구이고, 정확도는 다른 구현과의 일치로 답한다. **APEX 가 뽑은 표는 AutoPhOT 이
안 한 것이다.**

다만 한 편이므로 이것만으로 표준을 말할 수 없다. 계속 확인한다.

---

## 2. LSST 측광 알고리즘 비교 — Becker et al. (2007), PASP 119, 1462

**확인 경로**: IOPscience 원문. doi:10.1086/524710.
저자: Becker, Silvestri, Owen, Ivezić, Lupton.

LSST 요구사항을 놓고 측광 알고리즘들(Photo · DAOPHOT · ALLFRAME)을 비교한 논문이다.
**주입을 안 쓰고 실자료 대조만 한다.** 그리고 그 한계를 저자들이 직접 적어 놨다.

> *"Because the absolute 'truth' is not known here, these comparisons are by
> necessity relative."*

참값을 모르니 비교가 상대적일 수밖에 없다는 것이다. 그래서 Photo 의 결과를
참으로 **가정**하고 시작한다.

> *"start with the assumption that Photo's star-galaxy classification is 'truth.'"*

**이 문장이 사용자가 지적한 것과 같은 말이다** — 「다른기기에서도 서겠지, 근데
그거랑 값이 옳은가는 다르지」. 대조는 일치를 재지 정확성을 재지 못한다는 것을
PASP 논문이 자기 방법의 한계로 명시했다.

방법은 「여러 밤에 걸쳐 같은 천체를 재고 알고리즘 사이의 계통차를 본다」이다.

> *"Our analyses are designed to ascertain the level of systematics inherent to
> each photometry algorithm by comparing the measured properties of objects on
> multiple nights."*

**수치**: 밝은 쪽에서 구경 측광의 계통 바닥이 **0.007 등급**이고, LSST 요구가
**0.005 등급**이다. PSF 측광은 2~3 배 모자란다. 중심 위치는 PSF FWHM 의 1/100.

### 이 한 편이 말하는 것

**허용 한계의 문헌 근거가 하나 생겼다** — LSST 가 요구하는 측광 계통 오차가
0.005 등급이고, 2007 년에 실제로 달성된 것이 구경 측광 0.007 등급이다.
APEX 가 잰 ±11 mmag 는 이 두 값과 같은 자릿수다.

그리고 **대조만으로는 안 된다는 근거**가 이 논문 자체에 있다.

---

## 3. 붐비는 시야의 표준 관행 — 인공별로 계통 오차까지 잰다

여러 논문이 같은 말을 한다(출처 확정 진행 중).

> *"Artificial star tests have been the preferred method to assess the robustness
> of crowded-field photometry (Stetson & Harris 1988)."*

> *"The difference in input and recovered flux for ASTs provides a more realistic
> accounting of photometric uncertainties than the Poisson noise reported by the
> crowded-field photometric process alone."*

> *"This process delivers robust estimates of photometric performance including
> completeness, random errors, and systematic errors."*

**여기서 갈린다.** AutoPhOT 은 주입을 한계등급에만 썼는데, 붐비는 시야
(구상성단·근접은하 별 분해) 쪽에서는 **넣은 등급과 나온 등급의 차**를 재서
완전도·우연오차·**계통오차** 셋을 다 뽑는 것이 표준이라고 적혀 있다.
「gold standard」라는 표현까지 쓴다.

**APEX 가 뽑은 표가 정확히 이것이다.** 다만 위 인용의 출처를 아직 확정하지
않았으므로 다음 항목에서 원문을 잡는다.

### 출처 확정 — Jang (2023), MNRAS 521, 1532–1546

**확인 경로**: Oxford Academic 원문. doi:10.1093/mnras/stad619.
제목이 그대로 *"Tests of photometry: the case of the NGC 3370 ACS field"* 다.

**이 논문이 APEX 가 한 것과 똑같은 것을 한다.**

> *"Systematic errors (ΔF814W = F814W_input − F814W_output) are derived as a
> function of input F814W magnitudes."*

넣은 등급의 함수로 계통 오차를 낸다 — APEX 의 표가 정확히 이 모양이다.
흩어짐도 같은 방식으로 낸다.

> *"the standard deviation of the magnitude difference between the injected and
> recovered stars ... as a proxy for the statistical errors"*

그리고 Stetson & Harris (1988) 을 이 관행의 시작으로 인용한다.

> *"Artificial star tests have been the preferred method to assess robustness of
> crowded-field photometry"*

**보고한 수치** (DOLPHOT 과 DAOPHOT 을 같은 자료로 비교)

| 양 | 값 |
|---|---|
| 회수율 50 % 지점의 계통 오차 | DOLPHOT 약 0.05 등급 · DAOPHOT 약 0.03 등급(안쪽 영역) |
| 우연 오차 | DAOPHOT 이 DOLPHOT 보다 20~30 % 작다 |
| 어두운 끝의 등급 의존 편차 | 최대 0.15 등급 |
| **소프트웨어가 스스로 보고한 오차가 실제 오차를 과소평가하는 배수** | **1.25 ~ 2.68** |

**마지막 줄이 세 번째 확인이다.** Merline & Howell (1995) 이 CCD equation 이
S/N 을 과대평가한다고 예측했고, PhoPS 가 밝은 별에서 정규화 잔차 폭 2.069 를
쟀고, 이제 Jang 이 HST 자료에서 1.25~2.68 을 쟀다. **서로 다른 세 경로가 같은
방향을 가리킨다.**

### 관행의 시작 — Stetson & Harris (1988), AJ 96, 909

*"CCD Photometry of the Globular Cluster M92"*, bibcode `1988AJ.....96..909S`.

M92 의 깊은 프레임에 250 개짜리 인공별 묶음 여섯 벌을 원본 CCD 이미지 일곱 장에
더해 합성 프레임 42 장(각 3,000 개 이상)을 만들었다. 이렇게 얻은 오차 추정은
광자 통계와 별 혼잡을 **둘 다** 반영한다.

---

# 답 — APEX 의 방법은 표준이다. 다만 어느 쪽 표준인지가 갈린다.

문헌이 두 갈래로 나뉜다.

**첫째, 파이프라인 소프트웨어 논문.** AutoPhOT 이 대표다. 인공별 주입은
**한계등급을 정하는 데만** 쓰고, 「값이 옳은가」는 DAOPHOT 같은 기존 구현과
맞춰서 답한다.

**둘째, 분해된 별 측광과 붐비는 시야.** Stetson & Harris (1988) 이 시작했고
Jang (2023) 이 지금도 그대로 한다. 인공별을 넣고 **넣은 등급과 나온 등급의 차를
등급의 함수로** 낸다. 여기서 나오는 것이 완전도 하나가 아니라 완전도·우연
오차·계통 오차 셋이다. 이 갈래에서는 이것을 gold standard 라고 부른다.

**APEX 가 한 것은 둘째 갈래의 방법이고, 쓰고 있는 것은 첫째 갈래의 논문이다.**

그래서 답은 이렇다. **표준을 빠뜨린 것이 아니다.** 오히려 소프트웨어 논문
갈래에서 보면 남들이 안 하는 것을 했고, 성단 측광 갈래에서 보면 그 분야가
40 년 가까이 표준으로 삼아 온 것을 했다. 성단 CMD 를 내는 도구이므로 둘째
갈래의 기준을 따르는 것이 오히려 맞다.

## 대조만으로는 안 되는 이유가 문헌에 있다

Becker et al. (2007) 이 LSST 를 위해 측광 알고리즘들을 비교하면서 자기 방법의
한계를 이렇게 적었다.

> *"Because the absolute 'truth' is not known here, these comparisons are by
> necessity relative."*

**사용자가 지적한 것과 같은 말이다** — 「다른기기에서도 서겠지, 근데 그거랑
값이 옳은가는 다르지」. 인용할 수 있는 근거가 생겼다.

## 허용 한계의 문헌 근거

| 출처 | 양 | 값 |
|---|---|---|
| Becker et al. (2007) | LSST 측광 계통 오차 요구 | 0.005 등급 |
| Becker et al. (2007) | 2007 년 실제 달성(구경 측광, 밝은 끝) | 0.007 등급 |
| Jang (2023) | 회수율 50 % 에서 계통 오차 (DAOPHOT, 안쪽) | 0.03 등급 |
| Jang (2023) | 회수율 50 % 에서 계통 오차 (DOLPHOT) | 0.05 등급 |
| Jang (2023) | 어두운 끝의 등급 의존 편차 | 최대 0.15 등급 |
| **APEX (2026-09-03 측정)** | **11~18 등급의 계통 오차** | **0.011 등급 이내** |

**직접 비교는 조심해야 한다.** Jang 은 HST ACS 로 먼 은하의 별을 분해한 것이고
APEX 는 지상 망원경으로 성단을 본 것이라 혼잡도와 S/N 이 다르다. 다만 같은 양을
같은 방법으로 쟀으므로 자릿수 비교는 성립한다.

## 아직 안 본 것

SExtractor · photutils · DAOPHOT 원논문(1987) · STDWeb · 표준성야 관행 ·
저널 요구사항 · 인접 분야. 22 개 병렬 조사가 전량 실패해서 위 다섯 편만
직접 확인했다.

---

# 정정 — 허용 한계 표가 서로 다른 것들을 한 칸에 넣고 있었다

사용자 지적: *"허용한계들이 왤케 제각각임"*

맞다. 위의 허용 한계 표는 **세 종류를 한 열에 세웠다.**

**첫째, 곡선 위의 다른 지점이다.** Becker et al. 의 0.007 은 **밝은 끝의 계통
바닥**이고, Jang 의 0.03~0.05 는 **회수율 50 % 지점**, 즉 가장 나쁜 끝이다.
APEX 자신의 표에서도 편차가 11 등급 +0.0007 에서 18 등급 −0.0110 으로 열다섯 배
변한다. 문헌 값들이 갈리는 이유가 같다 — **어디서 쟀느냐가 다르다.**

**둘째, 양의 종류가 다르다.**

| 값 | 무엇인가 | 참값이 있나 |
|---|---|---|
| LSST 0.005 | **요구 사양**. 달성한 값이 아니라 달성해야 할 값 | — |
| Becker 0.007 | 알고리즘끼리의 **상대 계통차** | 없음 |
| Jang 0.03 · 0.05 | 주입한 참값 대비 **편차** | 있음 |
| APEX | 주입한 참값 대비 **편차** | 있음 |

**Jang 과 APEX 만 같은 종류다.** 앞의 둘은 나란히 놓을 수 없는데 놓았다.

**셋째, APEX 값이 뭉개진 값이었다.** 등급을 다 합쳐 냈더니 어두운 쪽 편차가
희석됐다.

## 같은 지점에서 다시 쟀다

프레임마다 자기 m50(회수율 50 % 등급)이 기록돼 있으므로, 그 언저리 ±0.25 등급의
주입별만 골라 다시 냈다. **Jang 이 잰 것과 같은 지점이다.**

| 프레임 | m50 | 50 % 지점의 계통 편차 | n |
|---|---|---|---|
| M13R_sharp | 15.78 | −0.0062 | 180 |
| M13V | 14.84 | −0.0006 | 137 |
| M67g_broad | 14.87 | +0.0041 | 163 |
| M67r_mid | 15.79 | −0.0029 | 171 |
| NGC6811R | 15.61 | −0.0053 | 145 |
| **NGC6811R_broad** | 14.21 | **+0.0334** | 156 |

프레임 중앙값 −0.0017, **절대값 최대 0.0334.**

**앞에서 「APEX 0.011 이 Jang 의 0.03~0.05 안쪽」이라고 적은 것은 틀렸다.**
같은 지점에서 재면 최악 프레임이 0.033 으로 **Jang 의 DAOPHOT 값과 같은
자리**다. 열 배 좋은 것이 아니었다.

`NGC6811R_broad` 한 장이 다른 다섯 장보다 여섯 배 나쁘다. 앞서 등급 뭉친
표에서도 이 프레임의 편차가 가장 컸다(+0.0092). **한 프레임의 문제이므로 원인을
따로 봐야 한다.**

밝은 끝은 비교가 안 된다. 주입 등급 범위가 대부분 프레임에서 m50 −3 등급까지
안 내려가서, 그 구간에 별이 있는 것이 M13V 한 장뿐이다(거기서는 편차 0.0000).
**Becker 의 밝은 끝 0.007 과 견주려면 더 밝은 별을 주입해야 한다.**

M67i 는 `completeness_fit.json` 이 없어 이 표에서 빠졌다(일곱 장 중 여섯 장).

---

# 문헌만 — 각 논문이 무엇으로 검증했나 (APEX 이야기 없이)

| 논문 | 검증 방법 |
|---|---|
| **Stetson & Harris (1988)**, AJ 96, 909 | M92 프레임 7 장에 250 개짜리 인공별 6 벌을 더해 합성 프레임 42 장을 만들고 되찾음 |
| **Stetson (1987)**, PASP 99, 191 — DAOPHOT | 주입 도구 `ADDSTAR` 를 소프트웨어에 넣어 배포. 표준 사용법은 「넣은 등급 대 되찾은 등급 + 잃은 것」 |
| **Bertin & Arnouts (1996)**, A&AS 117, 393 — SExtractor | **완전 합성 이미지**를 만들어 시험. 하늘 밝기 + Poisson 잡음 + Moffat 별 + 은하. 19 등급 별에 동반성을 거리별로 붙여 등급 오차를 잼 |
| **Dolphin (2000)**, PASP 112, 1383 — HSTphot | 같은 시야 반복 관측(IC 1613 의 F555W 2400 s 8 장 합성)과 DoPHOT 대조. **인공별 시험은 안 한다** |
| **Becker et al. (2007)**, PASP 119, 1462 | **주입 없음.** 여러 밤의 실자료로 알고리즘끼리 비교 |
| **Hu et al. (2010)**, PASP, doi:10.1086/658162 | 인공별을 이미지가 아니라 카탈로그에 더하는 계산 절약법 |
| **Brennan & Fraser (2022)**, A&A 667, A62 — AutoPhOT | 주입은 **한계등급 계산에만**. 정확도는 DAOPHOT 대조 |
| **Jang (2023)**, MNRAS 521, 1532 | 주입한 등급의 함수로 계통 오차를 냄. 같은 자료로 DOLPHOT 과 DAOPHOT 을 비교 |

## 방법은 넷이고 논문마다 조합이 다르다

1. **실제 이미지에 인공별 주입** — Stetson & Harris (1988) 이래. DAOPHOT 이
   `ADDSTAR` 로 도구를 제공. Dolphin, Jang 이 씀.
2. **완전 합성 이미지** — SExtractor 가 씀. 하늘부터 은하까지 다 만듦.
3. **다른 구현과 대조** — Becker(이것만), AutoPhOT(정확도용), Dolphin(외부),
   Jang(DOLPHOT vs DAOPHOT).
4. **같은 시야 반복 관측** — Dolphin(내부 일관성), Becker(밤 사이 계통차).

## 주입을 «무엇에» 쓰는지가 갈린다

- **한계등급만** — AutoPhOT
- **완전도 + 계통 편차 + 오차 크기** — Dolphin (2000), Jang (2023),
  그리고 Stetson & Harris 이래의 붐비는 시야 관행

## 두 논문이 같은 것을 발견했다

코드가 스스로 보고하는 오차가 실제보다 작다.

- **(정정) Dolphin (2000) 은 이 발견의 출처가 아니다.** 검색 조각을 원문 확인 없이 옮겼던 것이고, 원문은 인공별 시험 자체를 안 한다. 「25~50 %」는 DOLPHOT 을 쓴 후대 논문의 값으로 보이며 출처 미확정이다
- **Jang (2023)**: 코드 오차가 실제를 **1.25~2.68 배 과소평가**한다

---

# 왜 그 방법을 골랐나 — 논문이 밝힌 이유

앞의 표는 **무엇을 했나**만 적었다. 사용자 지적대로 논문들은 **왜 그것을
골랐는지**도 적어 놓는데, 그 이유가 방법 자체보다 중요하다.

## 인공별을 넣는 이유 — 코드가 매기는 오차는 혼잡을 못 담는다

Jang (2023) 이 명시적이다.

> *"artificial star experiments are critical to properly estimate photometric
> errors; using the errors computed directly from photometry codes underestimates
> the true errors, especially in crowded fields"*

> *"the internally computed errors are smaller than the errors determined from
> the artificial star experiments in all cases"*

즉 **주입은 「정확도를 자랑하려고」 하는 것이 아니라 「코드가 보고하는 오차를
믿을 수 없어서」 한다.** 소프트웨어가 계산하는 오차는 Poisson 통계에서 나오는데,
별이 겹쳐서 생기는 오차는 거기 안 들어간다. 넣어 보는 것 말고는 잴 방법이 없다.

Stetson & Harris (1988) 이 이 관행을 시작한 이유도 같다 — 광자 통계와 별 혼잡을
**둘 다** 반영하는 오차 추정을 얻으려고.

## 대조를 함께 하는 이유 — 주입은 모델 PSF 를 쓰기 때문에 못 보는 것이 있다

**이것이 이번 조사에서 가장 중요한 문장이다.**

> *"One example is the error associated with the use of incomplete PSF models...
> It is not possible to measure this bias from artificial star tests alone, since
> the artificial stars are injected using the model PSF, which differs from the
> real sources"*

인공별은 **코드가 가진 모델 PSF 로 그려서 넣는다.** 그러니 그 모델이 실제 별과
다를 때 생기는 편차는 주입으로 절대 안 보인다. 넣을 때 쓴 모양으로 다시 재니
잘 맞을 수밖에 없다. **주입에도 자기만의 순환이 있다.**

그 편차를 드러내는 것은 **독립적으로 처리한 다른 구현**뿐이다.

### 이것이 내가 앞서 적은 것을 뒤집는다

이 문서 앞부분과 `MAG_ACCURACY_20260903.md` 에서 나는 이렇게 적었다 —
「대조 상대가 없다는 것이 이 값의 힘이다. 참값을 우리가 넣었으므로 둘이 같은
방향으로 틀려서 일치하는 문제가 성립하지 않는다.」

**절반만 맞다.** 배경·잡음·혼잡에 대해서는 맞지만, **PSF 모델에 대해서는 주입이
바로 그 순환에 빠진다.** 그러니 주입이 대조보다 강한 것이 아니라 **둘이 서로의
사각지대를 덮는다.** Jang 이 둘 다 한 이유가 이것이다.

## 참값이 없을 때만 대조하는 경우

Becker et al. (2007) 은 주입을 아예 안 하는데, 그 이유를 이렇게 적었다.

> *"Because the absolute 'truth' is not known here, these comparisons are by
> necessity relative."*

이들이 다룬 것은 실제 하늘의 천체들이라 참값이 없었다. **참값을 만들 수 있으면
주입하고, 없으면 대조한다**는 선택 기준이 여기서 나온다.

## 연구 자체의 동기

Jang (2023) 이 이 시험을 한 이유는 소프트웨어 평가가 아니라 **거리 측정**이다.

> *"the mean photometric error for a single TRGB star at ∼20 Mpc is
> σ_F814W ∼ 0.15 mag in HST imaging...These individual errors are much larger
> than the final error of the Hubble constant"*

개별 별의 오차가 최종 결과의 오차보다 훨씬 크니, 그 오차가 어떻게 줄어드는지를
정확히 알아야 한다는 것이다.

## 아직 이유를 못 찾은 것

- **SExtractor** 가 완전 합성 이미지를 고른 이유 (A&AS PDF 가 403 으로 막힘)
- **Dolphin (2000)** 은 반복 관측과 DoPHOT 대조를 쓰면서 **이유를 명시하지 않는다**
- **AutoPhOT** 이 주입을 한계등급에만 쓴 이유
