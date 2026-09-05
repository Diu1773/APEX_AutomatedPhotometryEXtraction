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

**확인 등급.** [미확인] 이 붙은 줄은 원문을 열지 않고 검색 요약으로 적은 것이다.
2026-09-04 에 Dolphin (2000) 에서 검색 요약이 실제 논문과 다른 것을 확인했으므로
(Main/FAILURES.md F-253), [미확인] 줄은 확인 전까지 근거로 쓰지 않는다.

검색이 틀리는 방식이 특정된다. 소프트웨어 이름과 원논문 저자·연도가 다를 때
(DOLPHOT 대 Dolphin 2000, DAOPHOT 대 Stetson 1987) 검색은 **그 소프트웨어를 쓴
후대 논문들의 내용을 원논문의 것처럼 합쳐서** 돌려준다.

| 논문 | 검증 방법 |
|---|---|
| **Stetson & Harris (1988)**, AJ 96, 909 | M92 깊은 이미지 7 장에 인공별 1,500 개를 여섯 벌로 나눠 넣어 합성 프레임 42 장을 만들고, 되찾은 952 개의 (관측 − 입력) 차로 **계통 오차를 보정하고 우연 오차를 추정** — **[원문 확인]** |
| **Stetson**, *User’s Manual for DAOPHOT II* (2006 Apr 21판) | ADDSTAR 로 합성별을 넣고 되찾아 **star-finding efficiency 와 photometric accuracy 를 둘 다** 추정 — **[원문 확인]** |
| **Bertin & Arnouts (1996)**, A&AS 117, 393–404 — SExtractor | B 등급 10~27 의 별과 은하를 넣은 **완전 합성 이미지** 600 장(512×512)으로 **deblending · photometry · star/galaxy separation 셋을 시험** — **[원문 확인]** |
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


---

## ADDSTAR 의 목적 — Stetson 본인의 매뉴얼에서 (2026-09-04 원문 확인)

앞서 이 줄을 검색 요약으로 적고 [미확인] 을 달아 뒀다. PDF 를 직접 열어
확인했다.

**출처**: *User’s Manual for DAOPHOT II: The Next Generation*,
Peter B. Stetson, Dominion Astrophysical Observatory / Herzberg Institute of
Astrophysics. 이 판은 **2006 April 21**. 총 75 쪽, ADDSTAR 는 §XX (52 쪽).

전문 그대로다.

> *"This routine is used to add synthetic stars, either placed at random by the
> computer, or in accordance with positions and magnitudes specified by you, to
> your picture. They can then be found by FIND, reduced by PHOTOMETRY and the
> rest, and the star-finding efficiency and the photometric accuracy can be
> estimated by comparing the output data for these stars to what was put in."*

**두 가지를 다 하라고 적혀 있다** — star-finding efficiency(찾는 효율, 즉
완전도)와 **photometric accuracy(측광 정확도)**. 방법은 「넣은 것과 나온 것을
견주어서」다.

즉 **측광 소프트웨어를 만든 사람이 자기 도구의 설명서에, 주입한 별의 등급을
되찾는 정확도를 재라고 써 놓았다.** 이 관행의 출처를 더 거슬러 올라갈 곳이
없다.

곁들여, 같은 절이 재현성도 다룬다. 난수 씨앗을 사용자가 지정하게 해서
**어느 컴퓨터에서든 똑같은 인공 이미지가 만들어지도록** 했다.

### 다만 1987 년 논문은 아직 확인 못 했다

이것은 **DAOPHOT II 매뉴얼**이고, Stetson (1987) PASP 99, 191 은 DAOPHOT
Classic 논문이다. IOPscience 는 그 논문의 초록만 공개한다. 초록에는 FIND ·
PHOT · GROUP · NSTAR 만 나오고 **ADDSTAR 도 검증 절도 언급되지 않는다.**
그러므로 「1987 년 논문이 ADDSTAR 를 기술한다」고는 아직 쓸 수 없다.


---

## 관행의 원본 — Stetson & Harris (1988) §b) Artificial-Star Tests (2026-09-04 원문 확인)

ADS 스캔 PDF(79 쪽, 텍스트층 있음)를 직접 열었다. 이 절이 인공별 시험이 무엇을
위한 것인지, 그리고 **어떻게 해야 유효한지**를 다 적어 놓았다.

### 무엇을 위해 하는가

> *"the differences between the “observed” magnitudes and colors derived for the
> artificial stars and their known input values will be used to calibrate the
> systematic errors and to estimate the random errors that invariably accompany
> attempts to perform deep photometry in such a crowded field."*

**계통 오차를 보정하고 우연 오차를 추정하는 것**이 목적이라고 명시한다.
완전도만이 아니다.

### 왜 필요한가 — 코드의 내부 오차가 너무 작아서

> *"it will be shown from the artificial-star tests ... that these internal error
> estimates are still systematically too small, and accordingly we multiply this
> figure by the factor of 1.19 obtained below."*

**1988 년에 이미 같은 결론이 나와 있었다.** 소프트웨어가 계산한 오차가 계통적으로
작아서 **1.19 배**를 곱해 쓴다. Jang (2023) 의 1.25~2.9 배와 같은 이야기이고,
35 년 사이에 값만 달라졌다.

### 유효한 시험의 요건 셋 — 이게 이 논문의 진짜 값어치다

**첫째, 원본과 합성 프레임을 똑같이 처리해야 한다.**

> *"The original and artificial frames must be reduced identically for the
> comparisons to be valid"*

**둘째, 처리하는 동안 어느 것이 넣은 별인지 몰라야 한다.**

> *"A second requirement for a valid artificial-star test is that the reductions
> be performed without knowing which stars in the synthetic frames are added and
> which are real."*

그래서 이들은 profile-fitting 을 다 끝내고 표준 BV 계로 변환까지 마친 **뒤에야**
어느 검출이 넣은 별인지 대조했다.

**셋째, 넣는 양이 혼잡을 바꾸면 안 된다.**

> *"The principal difference is that the artificial frames contain of order 10%
> more stars than the original ones."*

원본보다 10 % 정도만 늘렸다. 많이 넣으면 재려는 그 혼잡 자체가 달라진다.

### 어떻게 묶어서 보고했나

되찾은 952 개를 **관측 등급으로**(입력 등급이 아니라) 정렬해 79 개씩 열한 묶음과
83 개 한 묶음으로 나눴다. 그 결과가 Table VII 「Artificial-star (observed − input)
comparisons」다.

> *"The actual data for the 952 recovered artificial stars were sorted by observed
> (not input) magnitude and divided into 11 groups of 79 and one group of 83."*

**입력이 아니라 관측 등급으로 묶는다는 점이 중요하다.** 어느 쪽으로 묶느냐에 따라
답이 달라진다.

### 결과의 방향

> *"stars at all observed magnitude levels tend, on average, to have been measured
> too bright—although for the brighter artificial stars the effect is only on the
> order of a few millimagnitudes, at fainter levels it becomes quite important."*

모든 등급에서 평균적으로 **실제보다 밝게** 측정되고, 밝은 쪽에서는 몇 mmag 이지만
어두워질수록 커진다.

### 저자들이 스스로 밝힌 함정

넣은 인공별의 광도함수가 F_input = 24.6 에서 잘려 있어서, 가장 어두운 묶음의
편차가 **과소평가되었을 수 있다**고 적었다. 더 어두운 별이 그 묶음으로 흘러들어올
자리가 없기 때문이다.

> *"it is possible that the effect of fainter artificial stars being scattered into
> the last bin is somewhat underestimated, due to the unphysical truncation of
> their luminosity function."*


---

## SExtractor 가 합성 이미지를 고른 이유 (2026-09-04 원문 확인)

A&AS PDF 는 403 이고 ADS 스캔은 JBIG2 순수 이미지라 텍스트층이 없었다.
PyMuPDF 로 12 쪽을 이미지로 렌더해 눈으로 읽었다.

**서지**: Bertin, E. & Arnouts, S., *SExtractor: Software for source extraction*,
Astron. Astrophys. Suppl. Ser. **117**, 393–404 (1996). 접수 1995-07-18,
게재확정 1995-08-17. Institut d’Astrophysique de Paris · ESO.

### 무엇을 시험했나 — 셋이다

부록 A 첫 문단이 명시한다.

> *"The simulated images we have used to test deblending, photometry and
> star/galaxy separation contain galaxies and stars with B magnitudes ranging
> from 10 to 27."*

**분리(deblending) · 측광 · 별은하 판별 셋을 다 합성 이미지로 시험했다.**
신경망 훈련용만이 아니다.

### 왜 합성인가 — 개수가 필요했고, 속도와 현실성을 맞바꿨다

> *"In order to simulate the large number of images needed for the neural network
> training, we have tried to find a compromise between realism and speed.
> Our concern was not to build a cosmological tool, but simply a fast code
> capable of producing convincing sky images."*

**우주론 도구를 만들려는 게 아니라 그럴듯한 하늘 이미지를 빨리 찍어내는 코드가
필요했다**고 스스로 적었다. 512×512 짜리 **600 장**을 만들어 훈련했고,
각 이미지를 8 가지 검출 문턱으로 돌려 **약 100 만 개** 항목의 목록을 얻었다.

### 합성이 실제보다 «더 어렵게» 되어 있다

> *"the crowding in the simulated images is higher than what one sees on real
> images of the field, allowing for the presence of many “difficult cases”
> (close double stars, truncated profiles, etc...) that the neural network
> classifier will have to deal with."*

혼잡도를 실제보다 **일부러 높였다.** 가까운 이중성이나 잘린 프로파일 같은 어려운
경우를 훈련에 넣기 위해서다.

### 무엇을 어떻게 만들었나

은하는 Schechter (1976) 광도함수에서 절대등급을 뽑고, 구형 성분은 de Vaucouleurs
법칙, 원반은 지수 프로파일을 쓴다. Hubble type 을 −5 에서 +10 사이에서 무작위로
고르고, 적색이동에 놓고 거리로 어둡게 한 뒤, **Moffat (1969) 함수로 PSF 를
씌운다.** 화소마다 3×3 으로 과표본화해 만들고 마지막 화소 크기로만 합성곱한다.

별은 **은하와 같은 등급-개수 분포**를 주어서, 훈련에 들어오는 어떤 패턴이든
별일 확률과 은하일 확률이 각각 50 % 가 되게 했다.

PSF 파라미터는 무작위다 — seeing FWHM 0.025~5.5 초각, Moffat β 2~4.
화소 크기는 항상 FWHM 의 0.7 배보다 작게 두어 표본화가 깨지지 않게 했다.

### 저자들이 밝힌 한계와 전망

> *"real data unavoidably differ a bit from simulated ones"*

합성과 실제가 다를 수밖에 없다는 것을 인정한다. 그러면서 이 접근을 확장할 수
있다고 적었는데, 확장 대상으로 **측정 과정 자체**를 꼽았다.

> *"could be advantageously extended to ... even to the measurement process itself
> (optimal determination of positions, magnitudes, etc.)"*

### 확인하지 못한 것

앞서 검색 요약으로 「하늘 밝기와 Poisson 잡음을 Gaussian 으로 더했다」고 적었는데,
**부록 A(402~403 쪽)에서 그 서술을 찾지 못했다.** 부록은 은하 모형과 PSF 합성곱
까지만 다룬다. 잡음 모형이 본문 다른 절에 있는지 아직 확인 못 했으므로
이 항목은 쓰지 않는다.


---

## 네 번째 방법 — 같은 별을 여러 밤 다시 재기 (2026-09-04 원문 확인)

**Lechapelain (2026)**, *Multi-night photometric repeatability of STDWeb: an
empirical error budget from a 14-night campaign on a single field*,
arXiv:2608.10017 (2026-08-08 투고).

앞서 「Dolphin (2000) 은 반복 관측을 쓰면서 이유를 명시하지 않는다」고 적었다.
이 논문이 그 이유를 정면으로 적어 놓았다.

> *"The formal magnitude uncertainty they report describes a single reduction of
> a single frame; campaign science instead requires the scatter of the same star,
> same field, night after night."*

**파이프라인이 내놓는 오차는 한 프레임을 한 번 처리한 것에 대한 값인데, 실제
연구는 같은 별을 밤마다 다시 잰 흩어짐을 필요로 한다.** 그래서 반복 관측으로
잰다.

### 방법과 결과

Einstein Probe 트리거 시야를 2026 년 7 월 **14 밤**에 걸쳐 관측한 **157 장**을
STDWeb 으로 똑같이 처리하고, 변하지 않는 비교별 **101 개**(10.5 < G < 14.5)를
분석했다.

| 양 | 값 |
|---|---|
| 파이프라인이 보고한 오차의 중앙값 | 8.1 mmag (밝은 절반은 6.7) |
| 하룻밤 안의 실제 흩어짐 | 보고값과 일치 (χ_within ≈ 1) |
| 밤마다 생기는 별별 영점 어긋남 | 7~9 mmag |
| 전 캠페인 흩어짐 | 11~12 mmag |
| **보고 오차 대 실제 오차** | **χ_camp = 1.5~1.7 (보고값이 실제의 약 절반)** |

**하룻밤 안에서는 맞는데 밤을 건너면 안 맞는다**는 것이 이 논문의 핵심이다.

### 원인과 처방

원인 전체가 **시각에 따라 달라지는 색항(epoch-dependent colour term)** 이었다.
카탈로그 계 등급 `m_sys = mag_calib + 색항 × (BP−RP)` 로 바꾸면 밤 항이
3.4 mmag 로 줄고, χ_camp 가 1.0~1.1 로 돌아오고, Gaia 기준 별별 재현성이
**43 → 9.5 mmag**, 밤별 영점의 최대-최소가 **130 → 11 mmag** 로 준다.

그리고 색을 광도곡선 자체에서 적합해도 Gaia BP−RP 를 0.027 등급(r = 0.98)으로
재현하므로, **외부 색 카탈로그 없이도 이 보정이 된다**고 적었다.

### 「보고 오차가 실제보다 작다」의 네 번째 확인

| 출처 | 값 |
|---|---|
| Merline & Howell (1995) | CCD equation 이 S/N 을 과대평가한다 (예측) |
| Stetson & Harris (1988) | 내부 오차가 작아 **1.19 배**를 곱해 씀 |
| Jang (2023) | 내부 오차가 실제를 **1.25~2.9 배** 과소평가 |
| Lechapelain (2026) | 캠페인 규모에서 **χ = 1.5~1.7** |

방법이 서로 다르다 — 하나는 이론, 둘은 인공별 주입, 하나는 반복 관측이다.
**서로 다른 네 경로가 같은 방향을 가리킨다.**

---

## 저널이 요구하는 것 — AAS 정책 원문 (2026-09-04 확인)

**출처**: *Policy Statement on Software*, AAS Journals, 2024 년 2 월 갱신.

> *"Such articles should contain a description of the software, its novel
> features and its intended use."*

> *"Such articles need not include research results produced using the software,
> although including examples of applications can be helpful."*

**검증이나 정확성 증명을 요구하는 조항이 없다.** 요구하는 것은 소프트웨어의
설명, 새로운 기능, 의도된 용도 셋이다. 연구 결과는 **넣지 않아도 된다**고
명시한다.

공개 방식은 권고이지 의무가 아니다 — 오픈소스 라이선스와 Zenodo/FigShare 아카이빙을
권하지만, *"any articles which provide a clear statement on how to access the code
– for example, by contacting the author – are acceptable"* 이라고 적혀 있다.

### 이것이 내가 앞서 한 말과 어긋난다

이 세션에서 나는 「PASP 소프트웨어 논문은 astrophysical use 와 실제 연구 결과의
예를 요구한다」고 적었고, 그것을 근거로 「과정 경험 중점만으로는 PASP 를 통과하지
못한다」고 말했다. **AAS 정책 원문은 정반대를 말한다.**

PASP 자체 안내문에 별도 조항이 있을 가능성은 남아 있으나,
`journals.aas.org/pasp-author-instructions/` 가 404 라 아직 확인하지 못했다.
**확인 전까지 위의 「PASP 가 연구 결과를 요구한다」는 주장은 근거 없음으로 둔다.**


---

## 정정 — 방금 인용한 AAS 정책은 PASP 에 적용되지 않는다 (2026-09-04)

바로 위에서 AAS 소프트웨어 정책을 인용하면서 「저널이 요구하는 것」이라고 적었다.
**그 정책이 PASP 를 규율하지 않는다.**

AAS 의 Scope Statements 페이지가 다루는 저널은 여섯이다 — **AJ · ApJ · ApJL ·
ApJS · PSJ · RNAAS.** **PASP 는 그 목록에 없다.**

IOPscience 의 PASP 저널 페이지가 관계를 명시한다.

> *"the technical journal of the Astronomical Society of the Pacific (ASP)"*

**PASP 는 IOP 가 태평양천문학회(ASP)를 위해 내는 저널이고, AAS(미국천문학회)의
저널이 아니다.** 두 학회가 다르다. 그러므로 AAS 의 소프트웨어 정책 원문
(「검증 요구 없음」·「연구 결과 넣지 않아도 됨」)은 **PASP 투고의 근거가 될 수
없다.** 위 절은 AAS 저널에 투고할 때의 기준으로만 읽어야 한다.

### PASP 자체가 밝힌 것 — 여기까지는 확인했다

IOPscience 의 PASP 저널 소개에서 논문 범주를 확인했다.

> *"Astronomical Software, Data Analysis, and Techniques"* — *"Original research
> that describes the software, data analysis, and research techniques used in all
> astrophysical contexts."*

**범주가 존재하고 그 정의가 「소프트웨어·자료분석·연구기법을 기술하는 원저
연구」다.** 다만 이 페이지는 검증이나 응용 사례를 요구하는지 말하지 않는다.

### 확인하지 못한 것 — 그리고 이 세션에서 내가 한 주장

PASP 의 실제 투고 안내문(`iopscience.iop.org/1538-3873/page/instructions_for_authors`)
은 **오늘 Radware 봇 검사(CAPTCHA)에 막혀 열지 못했다.** CAPTCHA 는 풀지 않는다.

이 세션에서 나는 「PASP 소프트웨어 논문은 astrophysical use 와 실제 연구 결과의
예를 요구한다」고 말했고, 그 근거로 2026-08-27 에 그 안내문을 읽은 것을 들었다.
그때 같은 문서에서 「익명 심사는 선택」을 확인했고 그 확인은 지금도 유효하다.
**그러나 「연구 결과의 예를 요구한다」는 부분은 오늘 다시 확인하지 못했다.**

**그러므로 이 주장은 「2026-08-27 에 읽었으나 재확인 불가」 상태로 둔다.**
논문 전략의 근거로 쓰기 전에 안내문을 다시 열어야 한다.


---

## photutils 는 심사받은 논문이 없다 (2026-09-04 원문 확인)

공식 인용 안내(`photutils.readthedocs.io/en/latest/getting_started/citation.html`)
가 이렇게 적는다.

> *"This research made use of Photutils, an Astropy package for detection and
> photometry of astronomical sources (Bradley et al. <YEAR>)."*
> *"...where (Bradley et al. <YEAR>) is a citation to the Zenodo record of the
> Photutils version that was used."*

**심사받은 저널 논문이 없고 Zenodo 기록만 있다.** 개념 DOI 는
`10.5281/zenodo.596036` 이고, 판마다 별도 DOI 가 붙는다(3.0.0 은
`10.5281/zenodo.19636730`, 2026-04-17 공개).

그러므로 **파이썬에서 가장 널리 쓰이는 측광 패키지에 「검증 절」이라는 것이
아예 없다.** 대조할 상대로 삼을 발표된 정확도 수치가 없다는 뜻이기도 하다.

---

## 인접 분야 — 의료영상은 같은 구조를 규격으로 써 놓았다 (2026-09-04 원문 확인)

**Choi, Kim, Ko, Cho, Jang, Ahn, Kim & Kim (2024)**, *Whole Process of
Standardization of Diffusion-Weighted Imaging: Phantom Validation and Clinical
Application According to the QIBA Profile*, **Diagnostics** 14(6), 583,
doi:10.3390/diagnostics14060583.

### 두 축을 요구한다

> *"imaging biomarkers for treatment response assessment should be validated for
> both accuracy (i.e., how close the ADCs are to the true values) and precision
> (i.e., how close the ADCs are between repeatable measurements)."*

**정확도(참값에 얼마나 가까운가)와 정밀도(다시 재면 얼마나 같은가) 둘 다**를
검증해야 한다고 적는다. 앞에서 문헌이 넷으로 나뉜다고 정리했는데, 그 넷 중
주입 계열이 정확도이고 반복 관측 계열이 정밀도다. **의료영상은 이 둘을 나란히
요구 조항으로 써 놓았다.**

### 참값은 팬텀에서 온다

실제 조직의 참 ADC 값은 알 수 없으므로, **알려진 값을 담은 물리 팬텀**을 찍어서
잰다.

> *"The QIBA developed the QIBA diffusion phantom to validate the accuracy and
> repeatability of DWI acquisition and ADC measurement."*

**인공별 주입과 같은 자리다** — 참값을 아는 대상을 일부러 만들어 넣고 되찾는다.

### 허용 한계가 «측정 전에» 문서로 정해져 있다

이것이 천문 측광에 없는 것이다. QIBA Profile 이 정한 값들이다.

| 양 | 허용 한계 |
|---|---|
| ADC 편차(bias) | ≤ 3.6 % |
| 단기 재현성 wCV | ≤ 0.5 % |
| 장기 재현성 wCV | < 2.2 % |
| 선형성 | R² > 0.90, 기울기 0.95~1.05 |
| b 값 의존성 | ≤ 2 % |
| 무작위 측정 오차 | ≤ 2 % |
| 신호대잡음비 | ≥ 45 |

**profile 이라는 문서가 「무엇을 재고 어느 선을 넘지 말아야 하는가」를 미리
정해 놓는다.** 측정한 뒤에 그 값을 보고 기준을 정하는 것이 아니다.

### 확인하지 못한 것

검색 요약에 QIBA 가 물리 팬텀 외에 **digital reference object(DRO)** 도 쓴다는
서술이 있었는데, **이 논문은 DRO 를 다루지 않는다.** DRO 는 QIBA Profile 원문
(Radiology, doi:10.1148/radiol.233055)에 있을 것으로 보이나 그 쪽은 403 이라
확인하지 못했다. 그러므로 DRO 는 쓰지 않는다.
