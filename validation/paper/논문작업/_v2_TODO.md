# 논문작업 — 초안 지도 & v2 반영 TODO

이 폴더 = **APEX 논문의 활성 초안 홈.** 상위 설계는 `../PAPER_PLAN.md v2`.

## 파일
- `MANUSCRIPT.md` — **영문 정본 초안** (§1–§6 + 백매터 완성, 2026-07-05 집필 + 07-09 Step0 반영).
  거의 완성 상태 — v2 델타만 반영하면 됨.
- `MANUSCRIPT_ko.md` — 국문 병렬본 (§1–§2 완료, §3–§6 미러 미완).
- `references.bib` (23건) · `REFERENCES.md` — 인용. (원본 복사본; `.tex`는 상위에 유지)

> **판단**: 기존 초안은 스텁이 아니라 near-submission 정본이라 삭제하지 않고 **활성 초안으로 승격**.
> 예전 것을 "참고만"이 아니라 이걸 직접 다듬는 게 맞음. PAPER.md/PAPER_PLAN.md/FIGURES.md는 상위에 원천으로 유지.

---

## ✅ 이미 v2-정렬된 것 (기존 초안이 잘 되어 있음)
- 포지셔닝: Abstract·§1이 이미 "접근성(알고리즘 아님)" + HOPS/IRAF/photutils 지형 = v2와 일치
- Step0 raw→science: §2.2·§3.11–3.14·Abstract에 이미 반영 ("raw 미지원" 옛 caveat 없음)
- 정직성 스파인: §3.1·§3.8·Table 1(direct/indirect/not-yet)·§5.2 다 있음
- 측광 방법: §2.2에 구경+apcorr+PSF, apcorr 성장곡선 서술 있음

## 🔧 반영 완료 (2026-07-10)
- **AIPPI 검증 제거** — Table 1 / §3.11 / §Software / §2.2 4곳. AIPPI는 provenance로만 강등
- **§3.11 실수치**: 합성 inject–recover(calibrated offset −0.016 DN·scatter 3.18 DN≈잡음바닥·
  bias/dark ~1 DN RMS·flat vignette 21→2 DN) + ccdproc(§3.13)로 AIPPI 대체. (커밋 e0f538e)
- **§3.2 완전도**: **erf 검출모델**로 확정 — AutoPhOT(brennan2022) 선례 확인 결과 걔네가 로지스틱
  아니라 error function β=½[1−erf(z/√2)] 사용, Masci 2011+Kashyap 2010 인용. 내 SNR-임계(7.4)
  물리와 정확히 일치 → 로지스틱 폐기, **Stetson1987(주입법)+Masci2011+Kashyap2010(erf)+brennan2022
  (선례)** 인용. Fig1 양패널 erf 재적합. bib에 kashyap2010·masci2011 추가. (커밋 53e5b4a)
- **§3.6 IRAF 매칭 + T2 표**: all-APEX NGC6811 V, 구경 1.0×FWHM·annulus 6-9 매칭 재실행
  (N=498, MAD 0.0092, r=0.99984). NGC457(g,278)과 짝지어 "단일프레임" 한계 해소. (커밋 17620df)

## 🔲 남은 v2 델타 (우선순위)

### A. 그림 재생성 연동 (본문 수치·문장 대기)
1. **§3.11 Fig10 패널**: 옛 (c)=AIPPI 등가성 → **합성 inject–recover 또는 ccdproc 패널로 교체**.
   본문은 이미 교체했으나 **합성 복원 정확 수치**(잔차 median 등)를 실측으로 채워야 함
   (calibration_validate.py 산출값 확인 → 삽입).
2. **§3.6 IRAF (Fig5→F6)**: NGC457 1프레임 → **NGC6811 all-APEX + 파라미터 매칭 재실행**으로 교체/보강.
   MAD 0.0097(현재 구경 0.8≠1.0 불일치) → 매칭 후 재측정. **파라미터 매칭 표(T2) 신설.**
3. **§3.9 CMD (Fig8→F15)**: all-APEX NGC6811로 재생성 (현 ridge 16mmag). PS1 위치 재매칭.
4. **§3.2 완전도 (Fig1→F2)**: 로지스틱 = "경험적 요약 + 인용"(Fleming 1995 대체) 라벨.
   **SNR-임계 물리 뒷받침 한 문장** 추가(전이 SNR≈7.4; 로지스틱=잡음-임계 누적가우시안 근사).
   밝은쪽 결손 ~4% 정직 명시(이미 문장 있음, 유지).

### B. 신규 검증 절 (자작 4코어 커버리지 완성)
5. **§3.3 천체측정 신설**: Step5 자작 quad 솔버 — **Gaia 잔차 분포 + 내부↔astrometry.net 일치**(F7).
   현재 Table 1에서 WCS="self-reported"인데, 이걸 direct로 승격. (지표 이미 산출됨)
6. **§3.7-ish 시계열 코어 신설**: SYSREM 주입복원 + PDM↔LS 일치(F17). LC 재처리 후.

### C. 실증·데이터
7. **§3.1 데이터 문장**: "previously reduced APEX result trees" → **"raw에서 APEX Step 0부터 전 과정 재환원"**
   으로 갱신. 대상 목록에 **M67·M3 추가**(5성단), 단일카메라 명시 유지.
8. **§4 (LC 실증)**: placeholder → **YZ Boo P=0.104092d 주기 복원**(F18). LC 재처리(--lc) 필요.
9. **§5.1 CMD 갤러리**: 5성단 all-APEX CMD(F16) 언급 추가.

### D. 마감
10. **국문 §3–§6 미러** (MANUSCRIPT_ko.md) — v2 반영 후 일괄.
11. **인용 추가**: Fleming1995(완전도)·Stellingwerf1978·VanderPlas2018·Kovács2002·Eastman2010 등
    lit-review에서 실검증(bibcode). references.bib 증분.
12. **Table 1 갱신**: WCS direct 승격, 시계열 코어 행 추가.
13. **F1 아키텍처 다이어그램 + GUI 스크린샷** (§2).
14. **백매터**: CRediT(교수 공저 확정 후) · Funding · Data availability에 reprocess 경로.

---

## 현재 데이터 (2026-07-10)
- all-APEX CMD: NGC6811·M67·M13 완주 / M3·M5 배치 진행 중
- LC(AE UMa·YZ Boo): Step0 대기 (E: 236GB 확보, `--lc`)
- 검증자산: Fig1–13 + fig_apex_cmd(16mmag) + fig_apex_iraf(0.0097, 매칭 재실행 예정)
