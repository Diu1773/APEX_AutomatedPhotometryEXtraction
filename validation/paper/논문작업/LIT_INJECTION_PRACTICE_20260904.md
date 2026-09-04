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
