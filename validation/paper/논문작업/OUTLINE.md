# APEX 논문 목차 (v2 구조, 2026-07-11)

기존 초안이 조각조각 쌓여 "왜 갑자기?"가 많아 → **구조부터 재확정.**
핵심 = 관통 논리선(spine) + 각 절이 그 선의 한 마디가 되게 여는 것.

## 관통 논리선 (spine)
> capable한 사람이 도구 벽(IRAF)에 막혀 진짜 과학을 못 한다 → APEX가 그 벽을 없앤다
> (가이드 raw→science) → 쉽고 AI로 만들었으니 "믿을 수 있냐"가 남는다 → 검증이 답이다
> (각 컴포넌트를 APEX가 안 만든 것에 대조) → "쉬운데 신뢰가능한" 과학이 된다.

모든 절이 이 선에 연결. §1끝(AI→검증 필요) ↔ §2끝(자작은 소수) ↔ §3시작(그래서 자작+산출물만 대조) 세 지점이 손잡음.

---

## §1 Introduction (~900w) — 문제와 기여
- 1.1 훅: "착착 돌아가는 SW"(HOPS) ↔ IRAF 벽 (진입장벽 서사, 저자경험 1-2문장 절제)
- 1.2 도구 지형·2축 갭 (표: IRAF/photutils/AutoPhOT/AstroImageJ/HOPS vs APEX)
- 1.3 APEX 한 문단 + 기여 4개 + "쉽고 AI개발 → 검증이 신뢰의 전제"(§3 예고)
- 1.4 로드맵

## §2 Design & implementation (~1400w) — APEX가 뭔지 (검증 해석 전제)
- 2.1 아키텍처 — GUI/Qt-free 코어 분리 = "검증한 게 곧 실행분"
- 2.2 **설계철학: 신뢰 라이브러리 래퍼 + 자작 코어 (Table 1)**  ← 구조의 핵심, §3 예비
    - 래퍼: 검출(photutils/sep)·구경/PSF측광(photutils)·WCS변환(astropy)·주기 LS/BLS(astropy)·BJD(astropy)
    - 자작: apcorr(photutils COG primitive 위 워크플로)·PDM(astropy 부재)·SYSREM(라이브러리 부재)·(옵션)quad솔버(의존성-free)
- 2.3 공유 측광체인 Step 0–7 (보정식·SEP검출·WCS·마스터·구경측광+apcorr·mag_inst·ZP)
- 2.4 모드분기 (CMD: PSF→ZP→CMD→이소크론 / LC: LC빌드→SYSREM→주기)
- 2.5 재현성 (TOML·JSON상태·캐시무효화·헤드리스)

## §3 Validation (~2800w) — 작동하는가
[여는 문단: 접근성·AI개발 → 신뢰 증명 필요. novelty=통합 → 산출물 정확성 + 자작코어만,
 각각 APEX가 안 만든 것(독립엔진·독립카탈로그·주입진실)에 대조. §2와 같은 파이프라인 순서.]
- 3.1 접근·데이터 (정직성 원칙 · 래퍼=인용+end-to-end / 자작=직접 · 단일기기 명시)
- 3.2 검출기 보정 (Step0) — 합성 inject-recover·ccdproc 비트일치·PTC·LCO 교차기기
- 3.3 검출·완전도 (Step4, photutils/sep) — erf 검출모델(Masci/Kashyap/AutoPhOT)·SEP 3.2σ
- 3.4 천체측정 (Step5) — Gaia 잔차 + 내장 quad솔버 vs ASTAP/astnet 일치  [신규]
- 3.5 구경측광·오차모델 (Step7, photutils+자작 apcorr) — pull·sep·IRAF(T2 파라미터표)·apcorr복원·민감도
- 3.6 PSF·crowded (Step8, photutils) — 구경↔PSF 일치 (M5·M13)
- 3.7 표준화·CMD (Step10-11) — Gaia·PS1 독립대조 + 정직성: BP faint bias(반증 사례)
- 3.8 시계열 코어 (LC, 자작 PDM·SYSREM) — 주입신호 복원·PDM↔LS·BJD sanity  [신규]
- 3.9 프레임 QC — 주입결함→판정 + 투명도 blind spot 정직보고

## §4 End-to-end science demonstrations (~900w) — 진짜 과학 산출
- 4.1 성단 CMD 갤러리 — all-APEX 5성단 (단일카메라 명시)
- 4.2 변광성 주기 — YZ Boo(문헌일치) / AE UMa(future work)
- 4.3 워크드 예제 — raw→CMD A4 1p, raw→주기 A4 1p (부록)

## §5 Discussion & limitations (~1000w)
- 5.1 확립된 것 · 5.2 scope-out(이소크론·LC) · 5.3 일반성 한계(단일기기·멀티앰프·구조배경·sub-res·apcorr 전역스칼라) · 5.4 실용권고 + 도구비교표(T5)

## §6 Conclusion (~350w) — 접근성 thesis 회귀 + 정직한 out-of-scope

Back matter: Data&Code · AI disclosure · CRediT · 부록 A4예제

---

## 표
- T1 스텝 × 라이브러리/자작 × 검증방식 (§2.2) ← 구조 핵심
- T2 IRAF/DAOPHOT 파라미터 1:1 매칭 (§3.5)
- T3 데이터 출처 (사이트+기기+기간)
- T4 검증결과 요약 (컴포넌트당 지표 한 줄)
- T5 도구 비교 매트릭스 (§5.4)

## 자작 4코어 정당화 (§2.2 각 한 줄, 담백하게 — 팩트체크 완료)
- **apcorr**: photutils `CurveOfGrowth` primitive는 있으나 프레임별 자동 apcorr 워크플로(참조별 선택→robust값→전소스 적용)는 없어 조립
- **PDM**: astropy.timeseries는 LombScargle·BoxLeastSquares만 제공, PDM 부재 → 구현(비정현파 변광용, 벡터화)
- **SYSREM**: 표준 라이브러리 부재 → 구현 (차등측광 계통제거 표준, Tamuz+2005)
- **quad 솔버**: 의존성-free 옵션. 기본은 ASTAP/astnet(외부); 내장솔버는 그 외부해를 레퍼런스로 검증
- (보정 산술도 자작이나 단순 사칙연산 + ccdproc로 직접검증됨 — 알고리즘 코어로 미분류)
