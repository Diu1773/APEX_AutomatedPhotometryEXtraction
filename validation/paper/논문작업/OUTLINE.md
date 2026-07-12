# APEX 논문 — 전체 목차 + 세부절 설계 (v3, 2026-07-11)

목차 레벨에서 각 세부절의 **역할 · 담을 내용 · 그림/표 · (검증절은) 검증대상·직접성**까지 설계.
프로즈는 이 설계가 확정된 뒤. 현재 MANUSCRIPT 실체(§3=16세부절, 파이프라인 순서, stub 2개) 반영.

기호: ★STUB=아직 그림/데이터 없음(자리표시) · [ref 필요]=외부 도구능력 확인 후 확정 · [결정]=설계 판단 대기.

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
| 1.2 도구 지형·갭 | "왜 기존 도구로 안 되나" | ①명령행(DAOPHOT/SExtractor) 진입문턱 ②파이썬(photutils) 스크립팅 문턱 → 유능한 사람도 막힘 ③그래픽 도구(AIJ·HOPS·MuniWin·VStar) 범위 좁음 + 공용 측광체인 부재 | **T_tools** [ref 필요] |
| 1.3 APEX+기여+pivot | 해법 + 신뢰 축 | APEX가 벽 없앰(가이드 raw→science, Qt-free 코어=검증분=실행분) / **pivot: 쉬움이 오류를 가림 + AI개발 → 검증이 신뢰의 전제** / 기여 (i)통합 종단 파이프라인 (ii)다층 재현가능 검증 | — |
| 1.4 로드맵 | 길잡이 | §2설계(+래퍼/자작 구분) → §3 파이프라인순 검증 → §4 과학 → §5 한계 → §6 결론 | — |

[결정] 저자 경험 노출: 1인칭 없이 본문 일반화 + AI disclosure에서 저자 개발 명시 (추천).

## §2 Design & implementation (~1400w)

| 세부절 | 역할 | 담을 내용 | 그림/표 |
|---|---|---|---|
| 2.1 아키텍처 | 검증해석 전제 | GUI(PyQt5)/Qt-free 코어(`apex.analysis`+`apex.pipeline`) 분리 = "검증한 게 곧 실행분" · 헤드리스 러너(멱등·JSON manifest) | — |
| 2.2 **설계철학: 래퍼 + 자작 4코어** | 구조의 핵 / §3 예비 | 라이브러리로 덮는 단계(검출·구경/PSF측광·WCS변환·주기LS/BLS·BJD=photutils/astropy)와 **자작 4코어**(apcorr·PDM·SYSREM·quad솔버). 이 구분이 §3 검증전략을 조직 | **T_split** (스텝×라이브러리/자작) |
| 2.3 공용 측광체인 Step 0–7 | APEX가 뭘 하나 | Step0 보정식 · **SEP 검출**(seg 아님) · WCS(ASTAP→astnet, Gaia게이트) · 마스터목록 · 구경측광+annulus · **apcorr(growth curve)** · mag_inst 계수율 · ZP(Gaia앵커·품질컷) | — |
| 2.4 모드 분기 | 갈래 | CMD: PSF(8)→ZP(10)→CMD(11)→이소크론(12) / LC: 선택(8)→광곡선(9)→SYSREM(10)→주기(11) | — |
| 2.5 재현성·설정 | 재현 근거 | TOML 단일설정 · JSON 상태영속 · 경로헬퍼 · 캐시 시그니처·무효화 · 헤드리스 | — |

[결정] 사실 정정: 구경 기본값 0.8→**실제 1.0×FWHM** 정합 · 워커 "75%"→현 캡 반영 · Step0 raw/precalibrated 스킵가드 문구.
[결정] **T_split(§2.2) vs Table1(§3.1) 관계** — 겹침. T_split=설계관점(라이브러리/자작), T1=검증커버리지(직접/간접). 둘 유지 or 통합? (추천: 유지, 역할 다름)

## §3 Validation (~3000w, 16세부절, 파이프라인 순서)

**여는 문단(3.1)**: 접근성·AI개발→신뢰 증명 필요. novelty=통합 → 산출물 정확성+자작코어만, 각각 APEX가 안 만든 것(독립엔진·독립카탈로그·주입진실)에 대조. §2와 같은 파이프라인 순서.

| 세부절 | 검증대상(스텝) | 근거·방법 | 그림/표 | 직접성 |
|---|---|---|---|---|
| **3.1** 접근·데이터 | — | 정직성 원칙 · 래퍼=인용+end-to-end / 자작=직접 · 단일기기 명시 | **T1**(커버리지)·**T3**(데이터출처) | — |
| **3.2** 검출기보정 | Step0 | 합성 inject-recover(offset−0.016DN/MAD3.18) + CR/hotpx 제거 | Fig_cal | 직접 |
| **3.3** 검출기특성화 | Step0 | PTC gain 0.681e⁻/ADU·RN2.35·dark0.0077 (헤더EGAIN 16×틀림) | Fig_ptc | 직접 |
| **3.4** ccdproc 교차 | Step0 | 픽셀 비트동일(Δ=0), 전파이프 5e-4DN | Fig_ccdproc | 직접 |
| **3.5** 교차기기 | Step0 | LCO QHY600(+0.06e⁻)·Sinistro 4앰프 vs BANZAI | Fig_lco | 직접 |
| **3.6** 검출 완전도 | Step4 | erf 검출모델(Masci/Kashyap/AutoPhOT)·SEP 3.2σ·m50=17.56 | Fig_compl | 직접 |
| **3.7** 천체측정 해 ★STUB | Step5 | Gaia 잔차 게이트 + 내장 quad솔버 vs ASTAP/astnet 일치 | **Fig_wcs ★신규** | 직접 |
| **3.8** 오차모델 | Step7 | pull 단위정규(std1.014,N3404)·σ=1.0857/SNR | Fig_pull | 직접 |
| **3.9** 파라미터·관측조건 민감도 | Step7 | 구경/하늘/시상 스윕·문턱robust | Fig_sens | 직접 |
| **3.10** SEP 합성교차 | Step7 | 독립엔진 SEP vs APEX, MAD0.006 | Fig_sep | 직접 |
| **3.11** IRAF 실데이터교차 | Step7 | IRAF/DAOPHOT vs APEX, MAD0.009(NGC6811 V) | Fig_iraf·**T2**(파라미터매칭) | 직접(최강) |
| **3.12** 밀집장 GC | Step8 | 구경↔PSF 내부일치, M5·M13 코어, null | Fig_crowd | 간접(내부) |
| **3.13** 참조카탈로그 | Step10 | PS1 독립대조 → B밴드 편차=Gaia BP 결함(정직사례) | Fig_ps1 | 직접 |
| **3.14** CMD 재현 | Step10-11 | Gaia·PS1 독립계 능선 19mmag 일치 | Fig_cmd | 직접(산출물) |
| **3.15** 시계열 코어 ★STUB | LC 10-11 | 주입신호 복원·PDM↔LS·BJD sanity | **Fig_ts ★신규** | 직접 |
| **3.16** 프레임 QC | Step3 | 주입결함→판정 + 투명도 blind spot 정직보고 | Fig_qc | 직접 |

## §4 End-to-end science (~900w)

| 세부절 | 역할 | 담을 내용 | 그림 |
|---|---|---|---|
| 4.1 성단 CMD 갤러리 | 통합 산출 실증 | all-APEX 5성단 raw→CMD (단일카메라 명시) | Fig_gallery ★재처리 완주 필요 |
| 4.2 변광성 주기 | LC 종단 실증 | YZ Boo(문헌 P=0.104092d 일치) / AE UMa=future work | Fig_yzboo ★LC재실행 필요 |
| 4.3 워크드 예제 | 재현 안내 | raw→CMD 1p, raw→주기 1p (부록) | 부록 |

## §5 Discussion & limitations (~1000w)

| 세부절 | 담을 내용 |
|---|---|
| 5.1 확립된 것 | 검출~측광~QC가 잘보정된 기기처럼 작동 (§3.6–3.16 종합) |
| 5.2 scope-out | 이소크론 물리량 복원(축퇴) · LC는 §4까지만 · 소광/매칭 부분검증 |
| 5.3 일반성 한계 | **단일기기(측광)** · 구조배경 · sub-res 겹침 · apcorr 전역스칼라 · 멀티앰프 |
| 5.4 실용권고 | 구경1.0-1.2FWHM · QC 2단계 · faint검증은 G/PS1 · 밀집장 내부진단 + **T5 도구비교** [ref 필요] |

## §6 Conclusion (~350w)
접근성 thesis 회귀 + 정직한 out-of-scope. 측광 주장, 이소크론/구조배경/sub-res/미검증기기 제외 명시.

Back matter: Data&Code · Software · AI disclosure · CRediT · 부록 A4예제

---

## 표 목록
- **T1** 검증 커버리지 (컴포넌트×checked-by×직접성) — §3.1 · *존재*
- **T2** IRAF/DAOPHOT 파라미터 1:1 매칭 — §3.11 · *존재*
- **T3** 데이터 출처 (사이트+기기+기간) — §3.1(or 5.3)
- **T4** 검증결과 요약 (컴포넌트당 지표 한 줄) — §3말 or §5.1
- **T5** 도구 비교 매트릭스 — §5.4 · **[ref 필요]**
- **T_tools** §1 포지셔닝 표 — **[ref 필요]** (T5와 통합 검토)
- **T_split** 스텝×라이브러리/자작 — §2.2

## 자작 4코어 정당화 (§2.2, 담백하게 — 팩트체크 완료)
- **apcorr**: photutils `CurveOfGrowth` primitive는 있으나 프레임별 자동 apcorr 워크플로(참조별 선택→robust값→전소스 적용) 부재 → 조립
- **PDM**: `astropy.timeseries`는 LS·BLS만, PDM 부재 → 구현(비정현파 변광·벡터화)
- **SYSREM**: 표준 라이브러리 부재 → 구현 (Tamuz+2005)
- **quad 솔버**: 의존성-free 옵션. 기본은 ASTAP/astnet(외부); 내장솔버는 그 외부해를 레퍼런스로 검증
- (보정 산술도 자작이나 단순 사칙연산 + ccdproc 직접검증 → 알고리즘 코어로 미분류)

## ★ 미완성 자산 (프로즈 전 필요)
- **Fig_wcs**(§3.7 천체측정) — 데이터 있음, 지금 생성 가능
- **Fig_ts**(§3.15 시계열) — LC 재실행 필요(YZ Boo)
- **Fig_gallery**(§4.1) — 5성단 재처리 완주 필요 (M3 stuck)
- **Fig_yzboo**(§4.2) — LC 재실행 필요
- **T5/T_tools** — 각 도구 능력 레퍼런스 확인 필요

## 그림 번호
현재 MANUSCRIPT 그림번호(Fig1=완전도…Fig13=교차기기)는 **옛 순서** — §3 재정렬로 어긋남.
프로즈 확정 후 파이프라인 순서로 Fig 일괄 재번호 필요.
