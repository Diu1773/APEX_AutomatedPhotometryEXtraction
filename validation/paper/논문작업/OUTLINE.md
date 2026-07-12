# APEX 논문 — 전체 목차 + 세부절 설계 (v4, 2026-07-11)

레퍼런스(AutoPhOT / AstroImageJ) 구조에 근거해 확정. 목차 레벨에서 각 절의
**역할 · 담을 내용 · 그림/표 · (검증절은) 검증대상·직접성**까지 설계. 프로즈는 이 설계 확정 후.

### 레퍼런스 근거로 확정한 구조 결정 (2026-07-11)
- **§3 = 9 상위절(주제그룹) + 하위절** — AutoPhOT·AIJ 둘 다 주제별로 묶음, 16개 flat 나열 안 함.
- **QC·완전도는 서사 배치 허용** — AutoPhOT가 완전도(limiting mag)를 검출 소속인데도 측광 뒤 독립절로 늦게 둠. → 프레임 QC 끝에 "정직 사례로 마무리" 유지.
- **도구 비교표 = 큰 매트릭스 빼고 산문 + 소규모 대조** — AIJ無·AutoPhOT 산문. §1은 산문 포지셔닝.
- **§2 설계 / §3 검증 분리 유지** — 우리 thesis(AI로 쉽게 만든 걸 검증으로 신뢰)는 검증 집중이 유리. AutoPhOT는 thesis가 파이프라인이라 섞음.
- 출처: [AutoPhOT 2201.02635](https://arxiv.org/abs/2201.02635), [AstroImageJ 1701.04817](https://arxiv.org/abs/1701.04817)

---

## 관통 논리선 (spine)
> capable한 사람이 도구 벽(IRAF)에 막혀 진짜 과학을 못 한다 → APEX가 그 벽을 없앤다(가이드 raw→science)
> → 쉽고 AI로 만들었으니 "믿을 수 있냐"가 남는다 → **검증이 답이다**(각 컴포넌트를 APEX가 안 만든 것에 대조)
> → "쉬운데 신뢰가능한" 과학이 된다.

세 접합점: §1끝(AI→검증 필요) ↔ §2끝(자작은 소수) ↔ §3시작(그래서 자작+산출물만 독립대조).

---

## §1 Introduction (~900w)

| 세부절 | 역할 | 담을 내용 | 그림/표 |
|---|---|---|---|
| 1.1 훅 | 기대-현실 간극 | HOPS로 transit 광곡선 한나절 → "측광=버튼하나" 기대 / 소형망원경 두 작업(CMD·변광성)은 아직 아님 | — |
| 1.2 도구 지형·갭 | "왜 기존 도구로 안 되나" | ①명령행(DAOPHOT/SExtractor) 진입문턱 ②파이썬(photutils) 스크립팅 문턱 → 유능한 사람도 막힘 ③그래픽 도구(AIJ·HOPS·MuniWin·VStar) 범위 좁음 + 공용 측광체인 부재 | **산문 포지셔닝** (표 아님) |
| 1.3 APEX+기여+pivot | 해법 + 신뢰 축 | APEX가 벽 없앰(가이드 raw→science, Qt-free 코어=검증분=실행분) / **pivot: 쉬움이 오류를 가림 + AI개발 → 검증이 신뢰의 전제** / 기여 (i)통합 종단 파이프라인 (ii)다층 재현가능 검증 | — |
| 1.4 로드맵 | 길잡이 | §2설계(+래퍼/자작 구분) → §3 파이프라인순 검증 → §4 과학 → §5 한계 → §6 결론 | — |

[결정·확정] 저자 경험: 1인칭 없이 본문 일반화 + AI disclosure에서 저자 개발 명시.

## §2 Design & implementation (~1400w)

| 세부절 | 역할 | 담을 내용 | 그림/표 |
|---|---|---|---|
| 2.1 아키텍처 | 검증해석 전제 | GUI(PyQt5)/Qt-free 코어(`apex.analysis`+`apex.pipeline`) 분리 = "검증한 게 곧 실행분" · 헤드리스 러너(멱등·JSON manifest) | — |
| 2.2 **설계철학: 래퍼 + 자작 4코어** | 구조의 핵 / §3 예비 | 라이브러리로 덮는 단계(검출·구경/PSF측광·WCS변환·주기LS/BLS·BJD=photutils/astropy)와 **자작 4코어**(apcorr·PDM·SYSREM·quad솔버). 이 구분이 §3 검증전략을 조직 | **T_split** (스텝×라이브러리/자작) |
| 2.3 공용 측광체인 Step 0–7 | APEX가 뭘 하나 | Step0 보정식 · **SEP 검출**(seg 아님) · WCS(ASTAP→astnet, Gaia게이트) · 마스터목록 · 구경측광+annulus · **apcorr(growth curve)** · mag_inst 계수율 · ZP(Gaia앵커·품질컷) | — |
| 2.4 모드 분기 | 갈래 | CMD: PSF(8)→ZP(10)→CMD(11)→이소크론(12) / LC: 선택(8)→광곡선(9)→SYSREM(10)→주기(11) | — |
| 2.5 재현성·설정 | 재현 근거 | TOML 단일설정 · JSON 상태영속 · 경로헬퍼 · 캐시 시그니처·무효화 · 헤드리스 | — |

[결정·확정] 사실 정정: 구경 기본값 0.8→**실제 1.0×FWHM** 정합 · 워커 "75%"→현 캡 반영 · Step0 raw/precalibrated 스킵가드 문구.
[결정·확정] T_split(§2.2 설계관점) 와 T1(§3.1 검증커버리지) 둘 다 유지 — 역할 다름.

## §3 Validation (~3000w) — 9 상위절 + 하위절, 파이프라인 순서

**여는 문단(3.1)**: 접근성·AI개발→신뢰 증명 필요. novelty=통합 → 산출물 정확성+자작코어만, 각각 APEX가 안 만든 것(독립엔진·독립카탈로그·주입진실)에 대조. §2와 같은 파이프라인 순서.

| 상위절 | 하위절 | 검증대상 | 근거·방법 | 그림/표 | 직접성 |
|---|---|---|---|---|---|
| **3.1** 접근·데이터 | — | — | 정직성 원칙·래퍼=인용+end-to-end/자작=직접·단일기기 명시 | **T1**·**T3** | — |
| **3.2** 검출기 보정 (Step 0) | 3.2.1 합성 inject-recover | Step0 | offset−0.016DN/MAD3.18 + CR/hotpx | Fig_cal | 직접 |
| | 3.2.2 검출기 특성화 | Step0 | PTC gain0.681·RN2.35·dark0.0077 (헤더EGAIN 16×틀림) | Fig_ptc | 직접 |
| | 3.2.3 ccdproc 교차 | Step0 | 픽셀 비트동일 Δ=0 | Fig_ccdproc | 직접 |
| | 3.2.4 교차기기 | Step0 | LCO QHY600·Sinistro vs BANZAI | Fig_lco | 직접 |
| **3.3** 검출 완전도 (Step 4) | — | Step4 | erf 검출모델(Masci/Kashyap/AutoPhOT)·SEP 3.2σ·m50=17.56 | Fig_compl | 직접 |
| **3.4** 천체측정 해 (Step 5) ★STUB | — | Step5 | Gaia 잔차 게이트 + 내장 quad솔버 vs ASTAP/astnet | **Fig_wcs ★신규** | 직접 |
| **3.5** 구경 측광 정확도 (Step 7) | 3.5.1 오차 모델 | Step7 | pull 단위정규(std1.014,N3404)·σ=1.0857/SNR | Fig_pull | 직접 |
| | 3.5.2 파라미터·관측조건 민감도 | Step7 | 구경/하늘/시상 스윕·문턱robust | Fig_sens | 직접 |
| | 3.5.3 SEP 합성 교차 | Step7 | 독립엔진 SEP vs APEX, MAD0.006 | Fig_sep | 직접 |
| | 3.5.4 IRAF 실데이터 교차 | Step7 | IRAF/DAOPHOT vs APEX, MAD0.009 | Fig_iraf·**T2** | 직접(최강) |
| **3.6** 밀집장 PSF (Step 8) | — | Step8 | 구경↔PSF 내부일치, M5·M13 코어, null | Fig_crowd | 간접(내부) |
| **3.7** 표준화와 CMD (Step 10-11) | 3.7.1 참조카탈로그(정직사례) | Step10 | PS1 독립대조 → B밴드 편차=Gaia BP 결함 | Fig_ps1 | 직접 |
| | 3.7.2 CMD 재현 | Step10-11 | Gaia·PS1 독립계 능선 19mmag | Fig_cmd | 직접(산출물) |
| **3.8** 시계열 코어 (LC) ★STUB | — | LC 10-11 | 주입신호 복원·PDM↔LS·BJD sanity | **Fig_ts ★신규** | 직접 |
| **3.9** 프레임 QC | — | Step3 | 주입결함→판정 + 투명도 blind spot 정직보고 (서사상 끝에 배치) | Fig_qc | 직접 |

## §4 End-to-end science (~900w)

| 세부절 | 역할 | 담을 내용 | 그림 |
|---|---|---|---|
| 4.1 성단 CMD 갤러리 | 통합 산출 실증 | all-APEX 5성단 raw→CMD (단일카메라 명시) | Fig_gallery ★재처리 완주 필요 |
| 4.2 변광성 주기 | LC 종단 실증 | YZ Boo(문헌 P=0.104092d 일치) / AE UMa=future work | Fig_yzboo ★LC재실행 필요 |
| 4.3 워크드 예제 | 재현 안내 | raw→CMD 1p, raw→주기 1p (부록) | 부록 |

## §5 Discussion & limitations (~1000w)

| 세부절 | 담을 내용 |
|---|---|
| 5.1 확립된 것 | 검출~측광~QC가 잘보정된 기기처럼 작동 (§3 종합) |
| 5.2 scope-out | 이소크론 물리량 복원(축퇴) · LC는 §4까지만 · 소광/매칭 부분검증 |
| 5.3 일반성 한계 | **단일기기(측광)** · 구조배경 · sub-res 겹침 · apcorr 전역스칼라 · 멀티앰프 |
| 5.4 실용권고 + 위치짓기 | 구경1.0-1.2FWHM · QC 2단계 · faint검증은 G/PS1 · 밀집장 내부진단 / **AutoPhOT식 소규모 대조**(산문 위주, 큰 매트릭스 아님) |

## §6 Conclusion (~350w)
접근성 thesis 회귀 + 정직한 out-of-scope. 측광 주장, 이소크론/구조배경/sub-res/미검증기기 제외 명시.

Back matter: Data&Code · Software · AI disclosure · CRediT · 부록 A4예제

---

## 표 목록 (레퍼런스 근거로 축소)
- **T1** 검증 커버리지 (컴포넌트×checked-by×직접성) — §3.1 · *존재*
- **T2** IRAF/DAOPHOT 파라미터 1:1 매칭 — §3.5.4 · *존재*
- **T3** 데이터 출처 (사이트+기기+기간) — §3.1(or 5.3)
- **T4** 검증결과 요약 (컴포넌트당 지표 한 줄) — §3말 or §5.1
- **T_split** 스텝×라이브러리/자작 — §2.2
- ~~T_tools / T5 대형 기능 매트릭스~~ → **폐기.** §1 산문 포지셔닝 + §5.4 소규모 대조(AutoPhOT식)

## 자작 4코어 정당화 (§2.2, 담백하게 — 팩트체크 완료)
- **apcorr**: photutils `CurveOfGrowth` primitive는 있으나 프레임별 자동 apcorr 워크플로(참조별 선택→robust값→전소스 적용) 부재 → 조립
- **PDM**: `astropy.timeseries`는 LS·BLS만, PDM 부재 → 구현(비정현파 변광·벡터화)
- **SYSREM**: 표준 라이브러리 부재 → 구현 (Tamuz+2005)
- **quad 솔버**: 의존성-free 옵션. 기본은 ASTAP/astnet(외부); 내장솔버는 그 외부해를 레퍼런스로 검증
- (보정 산술도 자작이나 단순 사칙연산 + ccdproc 직접검증 → 알고리즘 코어로 미분류)

## ★ 미완성 자산 (프로즈 전 필요)
- **Fig_wcs**(§3.4 천체측정) — 데이터 있음, 지금 생성 가능
- **Fig_ts**(§3.8 시계열) — LC 재실행 필요(YZ Boo)
- **Fig_gallery**(§4.1) — 5성단 재처리 완주 필요 (M3 stuck)
- **Fig_yzboo**(§4.2) — LC 재실행 필요

## 그림 번호
현재 MANUSCRIPT 그림번호(Fig1=완전도…Fig13=교차기기)는 **옛 순서**. 프로즈 확정 후 파이프라인 순서로 일괄 재번호.

## MANUSCRIPT 반영 대기 (프로즈 단계)
- 현재 MANUSCRIPT §3는 flat 16세부절 — 이 v4의 **9상위절+하위절 계층**으로 regroup 필요(3.2.x·3.5.x·3.7.x 하위 헤더 도입, 절번호 재조정).
- 국문본 §2 이하 조리개→구경, 검출 SEP화, 절번호 갱신 잔존.
