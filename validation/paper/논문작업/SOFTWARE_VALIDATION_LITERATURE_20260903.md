# 과학 소프트웨어가 옳다는 것을 어떻게 검증하는가 — 문헌 조사

이 조사의 대상은 성단 변수의 허용 한계가 아니라 **검증 방법론 자체**다. 소프트웨어
공학 쪽에 확립된 분야가 있고, 용어와 인용 근거가 이미 있다. 지금까지 APEX 가 한
검증에는 그 이름이 붙어 있지 않았다.

Research OS 서고에는 이 주제가 없다(172 건 훑어 match 0, SCAN 경로가 외부 검색을
요구함). 아래는 전부 외부에서 찾은 것이다.

---

## 우리가 겪는 어려움에는 이름이 있다 — the oracle problem

**test oracle problem** 은 「선택한 시험 입력에 대해 기대 출력이 무엇인지 정할 수
없거나, 실제 출력이 기대 출력과 맞는지 판정할 수 없는 어려움」을 가리킨다.

APEX 가 정확히 이 상황이다. 어떤 별의 참 등급을 모르므로 출력이 맞는지 직접 확인할
수 없다. 문헌값과 비교하는 것도 그 문헌값 자체가 다른 측정의 결과이므로 참값이
아니다.

**이 이름을 쓰면 3 절이 「우리가 이런저런 비교를 했다」에서 「oracle 이 없는 상황을
이렇게 다뤘다」로 바뀐다.** 후자가 심사에서 훨씬 강하다.

관련 조사: Kanewala & Bieman, *"Testing scientific software: A systematic literature
review"*, Information and Software Technology 56(10), 2014. **서지 확인 필요.**
그 조사는 과학 소프트웨어 시험의 어려움을 둘로 나눈다 — oracle problem 같은
소프트웨어 자체의 성질에서 오는 것과, 과학자와 소프트웨어 공학 공동체의 문화
차이에서 오는 것.

---

## oracle 이 없을 때의 표준 답 — metamorphic testing

출력이 맞는지 직접 못 보는 대신, **여러 번 실행한 결과 사이에 반드시 성립해야 하는
관계**를 검사한다. 그 관계를 metamorphic relation(MR)이라 부른다.

> *"Metamorphic relations (MRs) are necessary properties of the intended
> functionality of the software, and must involve multiple executions of the
> software. Any inconsistency (after taking rounding errors into consideration)
> indicates a failure of the program."*

1998~2015 년 조사에서 metamorphic testing 적용 사례의 **4 분의 1 이상이 과학
소프트웨어**였다.

### APEX 는 이미 하나를 하고 있다

**GUI 로 돌린 결과와 command line 으로 돌린 결과가 같아야 한다**는 것이 바로
metamorphic relation 이다. 같은 입력을 다른 실행 경로로 통과시켰을 때 출력이 같아야
한다는 성질이기 때문이다. 지금은 이걸 「재현성」이라 부르고 있는데, **정확한 이름은
metamorphic testing 이고 그렇게 부르면 인용 근거가 생긴다.**

### 싸게 추가할 수 있는 관계들

측광 파이프라인에는 반드시 성립해야 하는 관계가 여럿 있고, 대부분 참값을 몰라도
검사할 수 있다.

- **화소값을 전부 k 배** 하면 instrumental magnitude 가 −2.5 log₁₀(k) 만큼 이동하고
  **색지수는 안 바뀐다**
- **영상을 회전하거나 뒤집으면** 등급은 그대로이고 WCS 만 따라 바뀐다
- **배경에 상수를 더하면** sky subtraction 뒤 flux 가 안 바뀐다
- **입력 프레임의 순서를 섞으면** master source list 와 등급이 같아야 한다
- **한 프레임을 복제해 넣으면** 가중치가 예측대로 바뀐다

**셋째 것이 특히 값어치 있다.** 앞서 「bias·dark 의 몇 DN 은 sky subtraction 이
지운다」고 계산으로 주장했는데, 이 관계를 검사하면 그 주장이 실측이 된다.

**넷째 것은 이미 필요하다고 지적받았다.** ARS 심사에서 「입력 순서를 섞어 다시
돌리는 것을 3.12 에 포함할 것」이 나왔다. 그게 metamorphic relation 이다.

---

## 검증을 세 갈래로 나누는 표준 — Oberkampf & Roy

Oberkampf & Roy, *Verification and Validation in Scientific Computing*,
Cambridge University Press, 2010. 이 분야의 표준 교과서다.

세 갈래로 나눈다.

| 갈래 | 묻는 것 |
|---|---|
| **code verification** | 코드가 의도한 수식을 올바로 구현했는가 |
| **solution verification** | 이번 실행의 수치가 목적에 비추어 충분히 정확한가 |
| **validation** | 모형이 실제와 맞는가 |

**APEX 의 3 절이 지금 이 셋을 섞고 있다.** 나누면 이렇게 된다.

| APEX 절 | 갈래 |
|---|---|
| 3.3 라이브러리 호출 검증, ccdproc 비트동일 | code verification |
| 3.2 오차 예산, 허용 한계 판정 | solution verification |
| 4 절 문헌값 재현 | validation |

**이 구분이 「비트동일의 증거력이 0」이라는 심사 지적도 정리해 준다.** 비트동일은
code verification 이고, 그것이 solution verification 이나 validation 을 대신하지
않는다. 심사자가 지적한 것이 정확히 그 혼동이었다.

### method of manufactured solutions

정답을 아는 문제를 일부러 만들어 넣고 코드가 그것을 되찾는지 보는 방법이다.
미분방정식 코드의 code verification 표준 기법이다.

**APEX 의 대응물이 artificial star injection 이다.** 이미 완전도 측정에 쓰고 있고,
같은 것을 「정답을 아는 입력」으로 더 넓게 쓸 수 있다.

---

## 이 논문을 쓰는 이유가 되는 문헌 — Hatton

Hatton 등이 지구과학용 소프트웨어 여러 개를 조사해서, **각각이 그럴듯하지만 서로
본질적으로 다른 결과를 냈다**는 것을 보였다.

> *"Hatton et al. found that several software systems written for geoscientists
> produced reasonable yet essentially different results."*

**같은 계산을 하는 독립 구현들이 서로 다른 답을 낸다는 것이 실증된 것이다.**
그리고 그 조사가 이 분야의 검증 논의를 촉발했다.

서론이 「새로 만든 도구에는 빌려 올 신뢰가 없다」로 시작하는데, 이 문헌은 그보다
센 말을 할 수 있게 해 준다 — **오래 쓰인 도구도 서로 다른 답을 낸다.**

같은 조사에서 나온 것들도 적어 둔다. 체계적인 시험이 없으면 프로그램을 죽이지
않으면서 출력만 바꾸는 결함이 남고, 실제로 지진 자료 처리에서 정밀도가 손실되고
**소프트웨어 결함 때문에 발표한 연구를 철회한 사례**가 있다.

**서지 확인 필요** — Hatton 의 원 논문(*"The T-experiments"* 계열)과 Hatton & Roberts
의 지구과학 소프트웨어 조사.

---

## 3.1 을 어떻게 다시 쓸 수 있나

지금 3.1 은 「무엇을 얼마나 재는지 정하는 기준」이고 근거가 우리 판단뿐이다.
위 문헌을 쓰면 이렇게 된다.

1. **이 소프트웨어는 oracle problem 아래 있다.** 별의 참 등급을 모른다.
2. **그래서 검증을 세 갈래로 나눈다** — code verification, solution verification,
   validation (Oberkampf & Roy 2010).
3. **code verification 은 metamorphic relation 과 독립 구현 대조로 한다.**
   비트동일은 여기 속하고, 여기까지만 보인다.
4. **solution verification 은 오차 예산으로 한다.** 허용 한계는 같은 방법을 쓰는
   선행 연구에서 가져온다(Oliveira+2013).
5. **validation 은 4 절에서 문헌값과 대조해 한다.**

**이렇게 쓰면 3 절의 구조가 우리 취향이 아니라 확립된 틀이 된다.**

---

## 아직 안 한 것

**서지 확인 셋** — Kanewala & Bieman 2014, Hatton 의 원 논문, 그리고 metamorphic
testing 의 원 출처(Chen 등). 지금은 검색 요약만 봤고 **DOI 조작 금지 규약에 따라
확인 전에는 원고에 안 넣는다.**

**천문 쪽에 이 틀을 쓴 사례가 있는지** 아직 안 찾았다. 있으면 그 논문이 우리보다
앞선 것이고, 없으면 「소프트웨어 공학의 확립된 틀을 이 분야에 처음 가져왔다」가
된다. **어느 쪽인지가 이 논문의 신규성 주장에 직접 영향을 준다.**

**「한 프로그램에 다 넣으면 새로 생기는 위험」도 이 문헌군에서 찾아야 한다.**
소프트웨어 공학 쪽에 monolithic 대 modular 의 결함률 비교가 있을 수 있다.
