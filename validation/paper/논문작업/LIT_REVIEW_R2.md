# APEX 논문 — 레퍼런스 리서치 R2 (병렬 4축, 2026-07-18)

축1(도구지형)·축3(컴포넌트 실측검증관례)·축4(저널관례) 완료. **축2(합성vs실측 검증방법론 — 최우선)는 세션 한도로 실패, 재실행 대기.**
전 항목 원문/DOI 실존 확인. 확인 실패는 정직 표기.

---

## ★A. 지금 원고대로 쓰면 **오류**가 되는 것 (최우선 수정)

| # | 현재 상태 | 사실 | 조치 |
|---|---|---|---|
| A1 | SEP를 `bertin1996`으로 인용(§3.10·Software) | SEP 원전은 Barbary 2016 JOSS 1(6),58 doi:10.21105/joss.00058 | `barbary2016`으로 교체 (R1에서 이미 지적) |
| A2 | AIJ를 "단일시야·주기분석 없음"류로 서술 | **AIJ 6.0.0.00(2025-10-18)에 BLS/TLS/Lomb-Scargle 페리오도그램 추가**(Greg Srdoc) | "주기분석 부재" 문구 금지. → "테이블 기반 페리오도그램은 추가됐으나 다중밤 병합·앙상블 디트렌딩 가이드는 없음" |
| A3 | HOPS를 저널 논문처럼 인용(`tsiaras2019`) | **refereed 논문 부재.** EPSC-DPS 2019 초록(2019EPSC...13.1594T) + GitHub만 | 학회초록+소프트웨어로 표기 정정 |
| A4 | "Munipack(Hroch)" 단일 표기 | **별개 두 패키지**: Munipack(Hroch, CLI, ascl:1402.006) ↔ C-Munipack/Muniwin(Motl, GUI, SourceForge) | 분리 인용. GUI 논거=Motl, CLI CMD=Hroch |
| A5 | ASTROPOP을 엔진 교차검증 사례로 볼 여지 | **ASTROPOP의 IRAF(PCCDPACK) 비교는 편광 전용.** 측광은 SkyMaker 모의+표준성야 카탈로그 대조 | 엔진교차 사례로 인용 금지. 영점·표준화 사례로만 |
| A6 | PASP를 AAS 계열로 상정 | **PASP는 ASP 소유·IOP 발행.** AI 정책은 AAS가 아니라 **IOP 정책** 적용 | 백업 타겟 서술 정정 |

## ★B. 갭 논증 — 성립하나 **정밀화 필수** (축1 적대 점검 결과)

**결론: "raw→성단 CMD / raw→다중밤 LC를 가이드 GUI로" 갭은 반례 없이 성립.** 단 아래 3개를 선제 인용해야 리뷰어 반격을 막음:

1. **VaST** (Sokolovsky & Lebedev 2018, A&C 22,28; arXiv:1702.07715) — **우리 LC 모드와 기능적으로 가장 가까운 기존 도구.** 다중밤 + 프레임간 등급보정 + APASS/UCAC5 절대보정 자동화. 단 **터미널 실행 + PGPLOT 키보드 조작**, Windows는 WSL 필요, CMD 없음. → 갭을 "다중밤 자동화는 존재하나 가이드형 GUI가 아니다"로 정밀화.
2. **STDWeb** (Karpov 2025, Acta Polytechnica 65(1),50 doi:10.14311/AP.2025.65.0050) — STDPipe 위의 **웹 GUI**. "파이프라인형=CLI 전제"의 부분 반례. → "GUI화 흐름이 시작됐으나 transient 영역 한정"으로 프레임하면 오히려 APEX 동기 강화.
3. **Munipack(Hroch) M67 CMD 튜토리얼** (munipack.physics.muni.cz/cmd.html) — raw→CMD의 가장 근접 기존 경로. 단 **CLI 명령열 + TOPCAT 2단 조합**. → "가능하지만 가이드 GUI 아님"의 직접 증거.

**추가 발굴(채택 권고)**: EXOTIC(Zellem 2020 PASP 132,054401 — 시민과학, 단일 transit), Siril(Richard 2024 JOSS doi:10.21105/joss.07242 — GUI, 단일세션 LC까지), stellarphot(Zenodo 10.5281/zenodo.10679636 — reduced 입력, 노트북 위젯), Peranso(Paunzen & Vanmunster 2016 AN 337,239 — refereed, VStar류), APT(Laher 2012 PASP 124,737 — GUI 단일이미지 탐색).
**교육 축(§1 훅 보강)**: astrosource(Fitzgerald 2021 JOSS 6,2641 — 측광 *이후*만), OSS 파이프라인(Fitzgerald 2018 RTSRE — 서버측 완전자동=블랙박스), Makali'i(2006ASPC..351..544H), SalsaJ(2011IAUS..260..715F). → **교육 진영의 해법은 "완전 자동화(블랙박스)"이지 "가이드형 학습 경로"가 아님** = APEX 차별점과 정확히 대비.

**반례 정면 탐색 결과**: "GUI로 raw→성단 CMD" 현대 도구 **발견 실패 = 갭 성립**.

## ★C. 컴포넌트별 실측 검증 관례 (축3) — 우리 그림 사양의 근거

| 컴포넌트 | 채택할 실측 검증 관례 | 근거 레퍼런스 (실존 확인) |
|---|---|---|
| 천체측정(§3.7 stub) | Gaia 교차 ΔRA/ΔDec 산포+히스토그램, RMS(mas/″) vs 등급, **solve 실패율** | Ofek 2019 PASP 131,054504 (>5×10⁴ 실이미지, ~14 mas, 실패율 ≲2×10⁻⁵) · Masci 2019 ZTF (45–85 mas) · Lang 2010(성공률 관례) |
| 독립엔진 교차(§3.11) | 동일 실프레임 Δmag vs mag + σ대역, **계통차 발견 시 원인 귀속까지** | **Schechter 1993 PASP 105,1342 (DoPHOT vs DAOPHOT — 이 장르의 원조, faint/밀집 계통차를 sky fitting에 귀속)** · prose 2022(26 TESS 실관측 vs AIJ) · WFC3 ISR 2017-10(IRAF 다각형근사 ~0.1%) |
| 시계열 SYSREM(§3.15 stub) | **rms(unbinned/binned) vs 등급 전/후 + red noise 수치 + basis의 물리 원인 귀속** | **Collier Cameron 2006 MNRAS 373,799 Fig2/3 (red noise 2.5→1.5 mmag)** · NGTS 2018 Fig14(fractional RMS vs 등급 + 노이즈모델) |
| 주기(PDM, §3.15/§4) | **문헌 주기 회복률** vs 등급/관측수/클래스 | **Graham 2013 MNRAS 434,3423 (11 알고리즘, 실광곡선 ~67k)** → 우리 YZ Boo·AE UMa 문헌주기 재현 계획이 정확히 이 관례 |
| 구경보정(§2.3·§3) | 성장곡선(Δmag vs 반경) + apcorr 적용 후 **잔차 산포 vs 등급/위치** | Stetson 1990 · **Sirianni 2005 PASP 117,1049(ACS encircled energy)** · Bosch 2018 HSC(별 기반+공간변화 모델) |
| 프레임 QC(§3.16) | **인간 라벨 대비 일치율** (유일하게 확립된 틀) | Teimoorinia 2020 AJ 159,170 (CFHT 60k 라벨, 97% 일치) · **소규모 파이프라인 논문엔 QC 검증 관례 부재 — 본문에 명시해도 무방** |

**중요**: photutils의 refereed 자체검증 논문 **부재 확인** → WFC3 ISR 2017-10로 대체 인용.

## ★D. 저널 관례·제출 요건 (축4)

**RASTI 실게재 도구논문 실측 규범** (PyTICS 2025 RASTI 4,rzaf021 / DRUID 2025 RASTI 4,rzaf006 / GOTO 2026 RASTI 5,rzag042):
- 분량 12–18쪽, **그림 14–16개, 표 1–2** → **우리 13개는 규범 정중앙, 공식 상한 없음** (2–3개 여유)
- **검증 절이 본문의 40–50%**, 패턴 = 시뮬 1 + 실데이터 ≥2 도메인 + **경쟁 도구 교차비교** → 우리 구도(IRAF·Gaia·PS1·BANZAI)가 정확히 부합
- 비교 도구를 Methods 안 소절로 소개하는 관례(DRUID §3.2–3.3)
- 전수 결과·부속 도구는 **부록**으로

**제출 체크리스트(RASTI)**: `rasti` 클래스 2단 PDF ≤10MB · 초록 **≤250단어** · running head ≤45자 · 절 구조 고정(첫=Introduction, 끝=Conclusions) · endmatter 순서 **Ack→Data Availability→CoI→References→부록** · **그림마다 alt text**(신규 접근성 요건) · ≥400dpi · **소프트웨어 공개 + 버그리포트 절차 명시**(GitHub Issues로 충족) · AI 공개는 **Methods 또는 Acknowledgements + 커버레터**(커버레터 누락 주의)
- **Zenodo DOI는 RASTI 비강제**(PyTICS·GOTO 모두 GitHub만으로 게재). 단 무료이므로 릴리스 태그+Zenodo 발급 권장. PASP는 명시 권고.
- 비용: RASTI APC £1,339(재확인 요) + **waiver 제도 존재** / PASP $129/조판쪽 or OA $3,490 → **RASTI 1지망 재확인**

**AI 공개 문구 보강 필수**: IOP(PASP) 정책이 "생성형 AI는 원자료·결과·플롯을 생성/조작할 수 없음"을 명시 → 우리 공개문에 **"AI는 코드 작성을 보조했고, 모든 수치·그림은 그 코드가 실데이터에서 산출했으며 저자가 검증했다"**는 인과 구분 문장 추가.

## ★E. 미해결 — 축2 재실행 대기

**사용자 최우선 질문**: "합성 프레임에 인공별을 주입해 검증 = 자기참조 순환 아닌가. 표준은 **실제 관측 프레임에 주입**하는 것 아닌가?"
→ DAOPHOT ADDSTAR 계보·HST·DES가 실프레임 주입인지, 합성 검증의 역할 한정(수치 구현 정확성 vs 최종 정확도) 명문화 사례, referee의 simulation-only 비판 사례를 확인해야 §3.6·§3.10의 운명(실프레임 주입 전환 / 위상 격하 / 실측 앵커 추가)이 결정됨. **재실행 필요.**

## 확인 실패 (정직 표기)
HOPS refereed·ASCL 부재 / PP·prose의 raw 보정 포함 여부 본문 미확인 / MetroPSF·stellarphot·Afterglow·ASTAP refereed 서지 미확보 / RASTI 심사 익명모델·APC 정확금액 / PASP 이중익명 의무 여부 / Sirianni 2005 성단명·Schechter 1993 그림번호
