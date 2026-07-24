# -*- coding: utf-8 -*-
"""APEX 논문 지식그래프 — 정본 데이터(knowledge_graph.json) + HTML 렌더 생성.

- knowledge_graph.json  : 노드/엣지/카드 정본. Claude가 Grep으로 조회하는 RAG-인덱스.
  (노드 1줄 = json 1줄 — grep 친화. 카드에 경로·수치·근거 포함.)
- 지식그래프.html        : KMTNet 기상수치모델의 캔버스 엔진을 재사용해 렌더한 뷰어.

갱신 절차: 이 파일의 NODES/EDGES 수정 → 실행 → json+html 재생성 → 커밋.
    .venv-deploy\\Scripts\\python.exe "validation/paper/논문작업/build_knowledge_graph.py"
"""
from __future__ import annotations
import json, re
from pathlib import Path

HERE = Path(__file__).parent
ENGINE_SRC = Path(r"C:\Users\bmffr\Desktop\Me\2026-1-천문연_kmtnet\기상수치모델\references\지식그래프.html")
OUT_JSON = HERE / "knowledge_graph.json"
OUT_HTML = HERE / "지식그래프.html"

CATS = {
    1: "재현 그림",
    2: "레퍼런스 논문",
    3: "데이터 자산",
    4: "확정 결정·원칙",
    5: "조사·진단",
    6: "보류·미결",
}

def N(id, label, cat, card, title=None, type="n"):
    return {"id": id, "label": label, "title": title or label,
            "type": type, "cat": cat, "card": card}

NODES = [
    # ── 허브 ──
    N("hub_paper", "RASTI/PASP 논문", 0, """# APEX 도구 논문
RASTI/PASP급. 원고 `validation/paper/MANUSCRIPT_ko.md`(+영문). 구조 v4 = `OUTLINE.md`.
- 그림 헌장: [[dec_figure_spec|FIGURE_SPEC]] · 재현 현황: `REPRODUCTION_STATUS.md`
- 리뷰 R1: Major Revision(`REVIEW_ROUND1.md`) · 서지 40소스: `LIT_REVIEW.md`+`LIT_REVIEW_R2.md`
- 원칙: [[dec_repro_def|재현=방법/법칙]] · [[dec_vv_ladder|verification→validation 사다리]]""", type="hub"),
    N("hub_repro", "레퍼런스 재현 (/goal)", 0, """# 레퍼런스 검증-프레임 재현
사용자 /goal: "각 주요 논문 검증프레임을 할 수 있는 만큼 재현". 정본 = `REPRODUCTION_STATUS.md`.
- 완료 4: [[fig_completeness|완전도]] [[fig_astrometry|측성]] [[fig_iraf|IRAF교차]] [[fig_precision|정밀도floor]]
- 보류: [[pend_lc_period|LC주기]] [[pend_zp_resid|광도보정잔차]] [[pend_cog|CoG]] [[pend_bouguer|소광]]""", type="hub"),

    # ── 재현 그림 (cat1) ──
    N("fig_completeness", "F5 완전도(실측주입)", 1, """# F5 — 검출 완전도, 실측-프레임 주입
`fig_completeness_realvssynth.py` → `figures/fig_completeness_realvssynth.png`
- 재현: [[ref_haynes|Haynes 2002]] · [[ref_balrog|DES Balrog]] · [[ref_daophot|DAOPHOT ADDSTAR]] · [[ref_autophot|AutoPhOT App.D]]
- 실측 3+4프레임: M67 i **m50 17.65** / NGC6811 R **15.64** / M13 V **14.90** (+법칙검증용 4프레임 추가)
- (b) S/N 붕괴 패널: 7프레임 전부 peak-S/N축에서 단일 erf로 붕괴, **S/N50=4.05±0.18** (3.2σ 매치드필터와 일치) — [[dec_m13_ok|M13 이상치 아님]]
- 하단: 실주입별 컷아웃 사다리 m13→16 (`make_injection_cutouts.npz`, inject_flux_catalog 정규경로)
- 합성은 verification rung: S/N 기준 약간 보수적([[dec_synth_ok|2차 정정]]) · 조사: [[inv_completeness|완전도 조사노트]]"""),
    N("fig_astrometry", "F6 측성(Gaia잔차)", 1, """# F6 — 측성 검증
`fig_astrometry.py` · 재현: [[ref_ofek_masci|Ofek 2019 / Masci 2019(ZTF)]]
- 실측 66프레임(M13 15·M67 30·NGC6811 21) `step5_wcs/frame_wcs_qc.csv`
- Gaia DR3 잔차 RMS **중앙값 0.258″**, 66/66 solved, 피팅 없음(파이프라인이 실제 배포한 WCS)
- (b)패널은 별수-안정성으로 대체됨 — 원문은 노출/시간축. 개선 여지(사용자 지적)"""),
    N("fig_iraf", "F10 독립엔진(IRAF)", 1, """# F10 — 독립엔진 교차검증
`fig_iraf_crosscheck.py` · 재현: [[ref_schechter|Schechter 1993]] · AutoPhOT Fig14
- NGC6811 V 전-APEX 리덕션, 고정좌표 499별: APEX vs IRAF phot(DAOPHOT)
- **MAD 9.7 mmag · RMS 18.7 · r=0.99989**, binned median faint까지 평평(계통 없음)
- 데이터: `benchmark/runs/ngc6811_iraf_allapex_v1/phot_fixed_coords/fixed_comparison.csv`
- 기기의존도 최저(같은 프레임, 두 코드) — 재현 4종 중 가장 깨끗"""),
    N("fig_precision", "정밀도 floor(예비)", 1, """# 정밀도 floor — RMS vs mag
`fig_precision_floor.py` · 귀속: [[ref_honeycutt|Honeycutt 1992]] (앙상블 방법검증)
- 실측 M67 r 10프레임 1073별, Honeycutt 앙상블 ZP: **floor g/r/i = 5.4/5.2/6.7 mmag**
- 경험RMS = 보고 mag_err (photon영역 일치 = 에러모델 정직) + bright end 계통 floor 노출
- ⚠️ **10에폭 = under-powered**: 단일밤·binning 불가. floor값 자체는 에폭수 무감(5ep 5.18±0.33)
- ⚠️ [[pend_precision_fix|CC/Kovács 귀속 정정 필요]] — 원문 의도=detrending 전후 비교(SYSREM/TFA), 미재현"""),

    # ── 레퍼런스 (cat2) ──
    N("ref_haynes", "Haynes 2002", 2, """# Haynes 2002 (MNRAS 334,262)
합성 vs 실측 프레임 주입 **0.4–0.5 mag 계통 차** + 순환성 경고("merely confirm the starting hypothesis").
→ §3.6을 실측-프레임 주입으로 전환시킨 근거. APEX 합성은 S/N 기준 오히려 약간 보수적
([[dec_synth_ok|2차 정정]] — m50 우연 일치 철회, FWHM 트렌드로 설명)."""),
    N("ref_balrog", "DES Balrog", 2, """# Suchyta+ (DES Balrog)
"real pixel-level CCD images"에 소스 주입 — 실측 주입의 서베이 표준. 컷아웃 예시 그림 관행의 출처 중 하나.
HSC SynPipe도 동일 계열."""),
    N("ref_daophot", "Stetson (DAOPHOT)", 2, """# Stetson 1987 — DAOPHOT ADDSTAR
인공별 실측 주입의 원조 루틴. crowded-field 측광 표준."""),
    N("ref_autophot", "AutoPhOT (A&A 667 A62)", 2, """# Brennan & Fraser 2022 — AutoPhOT
QC·벤치마크 준거 논문. App.D=주입-회수(S/N erf), Fig12=주입 컷아웃, Fig14=IRAF 대조.
[[fig_completeness|F5]]·[[fig_iraf|F10]] 둘 다 이 논문 형식 재현."""),
    N("ref_ofek_masci", "Ofek 2019 · Masci 2019", 2, """# Ofek 2019 (PASP 131) · Masci 2019 (ZTF)
서베이 측성 검증 보고 형식: 프레임별 Gaia 잔차 RMS + solve 신뢰도. [[fig_astrometry|F6]]이 재현.
PP(Mommert 2017) 0.3″을 참조선으로 병기."""),
    N("ref_schechter", "Schechter 1993", 2, """# Schechter+ 1993 (PASP 105,1342) — DoPHOT
신규코드 vs 기성코드(DAOPHOT), 같은 별·고정좌표, Δmag-vs-mag 평평성으로 계통 귀속.
[[fig_iraf|F10]]이 재현."""),
    N("ref_honeycutt", "Honeycutt 1992", 2, """# Honeycutt 1992 (PASP 104,435)
비교성 지정 없는 **앙상블 차등측광**(모든 별+프레임ZP 동시해). RMS-vs-mag은 앙상블이
광자한계에 도달함을 보이는 방법검증. [[fig_precision|floor 그림]]의 올바른 단독 귀속."""),
    N("ref_cc_kovacs", "CC 2006 · Kovács 2005", 2, """# Collier Cameron 2006 (WASP) · Kovács 2005 (TFA)
**detrending 논문** — 의도는 SYSREM/TFA 전후 floor 하강 비교(수천 에폭, 트랜짓 감도).
10에폭 단일밤으론 재현 불가 → [[pend_lc_period|LC 재현]]에서 SYSREM 전후로 제대로.
현 floor 그림에서 귀속 제거 예정([[pend_precision_fix|정정 미결]])."""),
    N("ref_period", "Stellingwerf · Graham", 2, """# Stellingwerf 1978 (PDM) · Graham 2013 (주기복원 비교)
주기 복원 검증 준거. AE UMa(P=0.086017d)·YZ Boo(P=0.104092d) 문헌 확정값 보유.
LC 전처리 deferred라 보류: [[pend_lc_period|LC 주기복원]]."""),

    # ── 데이터 자산 (cat3) ──
    N("data_reprocess", "reprocess 3성단", 3, """# E:\\APEX_validation\\reprocess\\
전-APEX 헤드리스 재처리(raw→science). step7 forced phot: M13 15fr(BVR) · M67 30fr(gri) ·
NGC6811 21fr(BVR). M5는 step7 없음. WCS QC 66프레임 전부 solved.
과학프레임 위치: `{target}/sci/pp_*.fit` (M13만 calibrated/ 사용례 있음)."""),
    N("data_injections", "실측주입 7런", 3, """# validation/paper/data_realframe_*
run_artificial_star_suite(reference_frame=실측) 60trial×50별. stars.csv는 git 추적.
- 대표3: M67i(m50 17.65) · NGC6811R(15.64) · M13V(14.90)
- 법칙용4(2026-07-24): M67r_mid · M67g_broad · NGC6811R_broad · M13R_sharp — [[done_lawtest|법칙 검증 완료]]
- 컷아웃: `data_realframe_M13V/injection_cutouts.npz`"""),
    N("data_iraf", "IRAF 교차 런", 3, """# benchmark/runs/
`ngc6811_iraf_allapex_v1`(499별, 본검증) · `ngc457_iraf_crosscheck_g0016_v1`(278별, 구버전).
PyRAF는 WSL. fixed_comparison.csv 소비."""),
    N("data_synth", "합성 21k 주입", 3, """# validation/paper/data/artificial_star/
합성 참조프레임(FWHM 3.4px·sky 150) 21,000주입. m50=17.59. verification rung 전용
([[dec_vv_ladder|사다리 프레이밍]]). 실측과 혼용 금지."""),
    N("data_lc", "LC 데이터(잠김)", 3, """# E:\\observed_Analysis (LC 전처리 대기)
AE UMa·YZ Boo 등 변광성 시계열. **사용자 수동 삭제 대기로 100GB 전처리 deferred** —
[[pend_lc_period|LC 주기복원]]과 CC/Kovács 제대로 된 재현의 선결조건."""),

    # ── 결정·원칙 (cat4) ──
    N("dec_repro_def", "재현=방법/법칙", 4, """# "재현"의 정의 (2026-07-23 확정)
레퍼런스들은 전부 다른 기기(CCD 서베이 등), APEX=IMX455 CMOS 소형망원경.
- **비교 불가**: floor 절대값·한계등급·광자곡선 위치 → "WASP만큼 정밀"류 문구 금지
- **재현 대상**: 검증 프로토콜(주입·앙상블·독립엔진) + 보편 법칙(RMS모양, m50∝σ·FWHM²)
- 수치는 항상 "APEX(C3-61000 CMOS)에서" 로 보고. `REPRODUCTION_STATUS.md`에 명문화."""),
    N("dec_synth_ok", "합성=보수적(정정)", 4, """# 합성 프레임 판정 — 2차 정정 (2026-07-24)
~~"M67 실측 m50=합성이라 낙관 아님"~~ ← **철회**: m50 일치는 노이즈↑·PSF샤프↑ 상쇄의
파라미터 우연(사용자 의심 적중). 직접 계산하면 합성 S/N50=**5.19**로 실측 법칙(4.05±0.18)
위에 없음. **새 체계**: S/N50이 커널 FWHM에 단조 의존(ρ=-0.982; 샤프할수록 높음, 기전=
최소 픽셀면적 요구). 합성(3.4px)도 이 트렌드 위 → **S/N 기준 합성은 약간 보수적**.
상세: COMPLETENESS_REALFRAME_INVESTIGATION.md 후속정정 절."""),
    N("dec_m13_ok", "M13 이상치 아님", 4, """# M13 V 프레임은 그냥 얕은 프레임
m50 격차 완전 분해: 하늘밝기 1.92 + seeing 0.82 = 2.75 mag (관측 2.75, Δ0.00).
법칙 m50 = C − 2.5logσ − 5logFWHM 이 3프레임 Δ≤0.04로 성립(지수 고정).
peak-detection 절벽 메커니즘: mag15 별 peak 84e < 3.2σ 임계 96e."""),
    N("dec_vv_ladder", "V&V 사다리", 4, """# verification → validation 프레이밍
Oberkampf&Trucano 2002 · Portillo 2020 3-rung. 합성=기계 수치정확성(verification),
실측=실성능(validation). §3 전반과 [[fig_completeness|F5]] 범례에 적용."""),
    N("dec_figure_spec", "FIGURE_SPEC 헌장", 4, """# 그림 전면 재작성 헌장 (`FIGURE_SPEC.md`)
10개 공통병 진단(패널제목=변수명·내부ID 노출·suptitle 박힘·AIPPI 패널 등).
규칙: message-first 제목 · 내부 데이터셋ID 금지 · **AIPPI는 검증 그림에서 제외**(자체툴=순환) ·
F1-17 파이프라인 순서 번호."""),
    N("dec_citation_fixes", "인용 오류 6건", 4, """# 확정 인용 정정 목록 (LIT_REVIEW_R2 ★A)
1) SEP→**barbary2016** (bertin1996 아님) 2) AIJ 6.0은 주기도 **있음** 3) HOPS=EPSC 학회초록
4) Munipack(Hroch,CLI)≠Muniwin(Motl,GUI) 5) ASTROPOP IRAF비교=편광측정 한정 6) PASP=IOP.
원고 반영 미완 → 다음 원고 패스에서."""),

    # ── 조사·진단 (cat5) ──
    N("inv_completeness", "완전도 조사노트", 5, """# COMPLETENESS_REALFRAME_INVESTIGATION.md
2.7mag 격차 전체 추적 기록: flux스케일 무죄 → PSF폭=실측seeing 확인 → peak-detection 절벽 →
파이프라인 자체 검출롤오프와 독립 일치 → 다중프레임 확장(법칙). [[dec_m13_ok|결론]]·[[dec_synth_ok|반증]]의 근거문서."""),
    N("inv_detector", "검출기 실측", 5, """# C3-61000/IMX455 특성 실측
gain **0.689 e/ADU**(PTC) — 헤더 EGAIN 0.0495는 14× 오류(펌웨어 미구현). RN 2.1e·dark 0.0077e/s.
주입 flux 전자환산의 기반. 메모리 project_detector_characterization."""),
    N("inv_bfilter", "B필터 판결(종결)", 5, """# B필터 faint 편차 — PS1 대조 종결
Gaia BP faint 감광 유죄(ref +0.022), APEX 무죄(vs PS1 +0.008 평평).
교훈: faint 검증에 BP기반 참조 금지 — [[pend_zp_resid|광도보정 잔차]] 설계 시 필수 반영."""),

    # ── 보류·미결 (cat6) ──
    N("done_lawtest", "법칙 검증(완료)", 5, """# 7프레임 검출법칙 검증 — 완료 (2026-07-24)
m50 = C − 2.5logσ − 5logFWHM, **지수 이론고정·C만 fit: 잔차 RMS 103 mmag** (동적범위 3.4mag, N=7).
- 추가 4런: M67r_mid·M67g_broad·NGC6811R_broad·M13R_sharp (off-diagonal로 σ·FWHM 항 분리)
- **교훈**: FWHM은 frame_stats(밝은별)가 아니라 **각 런의 주입커널**(empirical_psf.fits 반높이면적)로
  — soft 프레임에서 9.30 vs 7.23px 괴리가 0.4mag 이탈 일으켰음(추적해 해소)
- 이 법칙이 [[dec_repro_def|기기무관 재현]]의 실체 — 절대 깊이가 아니라 법칙이 검증 대상"""),
    N("pend_precision_fix", "floor 귀속 정정", 6, """# 정밀도 floor 그림 정정 (미승인 제안)
CC/Kovács 귀속 제거 → Honeycutt 단독 + "single-night N=10 예비" 라벨.
확정판은 [[pend_lc_period|LC]]에서 SYSREM 전후 비교로. **사용자 승인 대기.**"""),
    N("pend_lc_period", "LC 주기복원", 6, """# LC 재현 (Stellingwerf/Graham + CC/Kovács 진짜 재현)
선결: 사용자 observed_Analysis 수동삭제 → LC 전처리(100GB) → step1-11.
목표: AE UMa·YZ Boo 주기 문헌 대조 + SYSREM 전후 RMS-vs-mag(수백 에폭)."""),
    N("pend_zp_resid", "광도보정 잔차", 6, """# ZP·색항·잔차 그림 (Padmanabhan/Schlafly식)
데이터는 cmd_zeropoint/에 있으나 **Gaia-CMD 보정계=논란영역**(이소크론·B필터 이력).
어느 레퍼런스로 보정할지 사용자 결정 필요. [[inv_bfilter|BP 금지 교훈]] 적용."""),
    N("pend_cog", "Curve-of-growth", 6, """# 구경보정 CoG (Stetson/Howell)
apcorr_summary는 프레임당 단일값뿐 — 다중구경 재실행 필요."""),
    N("pend_bouguer", "소광 Bouguer", 6, """# zp vs airmass (Hardie 1962)
frame_zeropoint.csv에 airmass 있음. 단 ZP 레퍼런스 얽힘 + airmass span 미확인."""),
    N("pend_fig_a_layout", "fig(a) 배치 미결", 6, """# fig (a) 등급공간 완전도의 논문 배치 — 사용자 판단 대기
관행 4계열 조사 완료 → `FIG_A_PRESENTATION_SURVEY.md`. 옵션: A=현상유지(도구논문식 예시),
B=66프레임 예측 m50 분포 추가(Balrog식), **C★=predicted-m50 QC 게이트 구현 후 '예측 vs 실현
깊이' 운영검증 그림(Kessler/ATLAS식, 백로그 task_5bb4af3e 선행)**. 실별 6% 예측실증도 C에 포함.
사용자: "좀 더 생각해봐야" — A/B/C 미결."""),
    N("pend_manuscript", "원고 §3.6 갱신", 6, """# 원고 반영 대기
§3.6을 다중프레임 실측주입+법칙으로 갱신 · §3.7(측성)·§3.15(시계열) 스텁 채우기 ·
[[dec_citation_fixes|인용 6건]] 반영 · Data Availability 그림번호 깨짐 수정."""),
]

EDGES = [
    # hubs
    ("hub_paper", "hub_repro"),
    ("hub_paper", "dec_figure_spec"), ("hub_paper", "dec_citation_fixes"),
    ("hub_paper", "dec_vv_ladder"), ("hub_paper", "pend_manuscript"),
    ("hub_repro", "fig_completeness"), ("hub_repro", "fig_astrometry"),
    ("hub_repro", "fig_iraf"), ("hub_repro", "fig_precision"),
    ("hub_repro", "dec_repro_def"),
    # F5
    ("fig_completeness", "ref_haynes"), ("fig_completeness", "ref_balrog"),
    ("fig_completeness", "ref_daophot"), ("fig_completeness", "ref_autophot"),
    ("fig_completeness", "data_injections"), ("fig_completeness", "data_synth"),
    ("fig_completeness", "inv_completeness"), ("fig_completeness", "dec_synth_ok"),
    ("fig_completeness", "dec_m13_ok"), ("fig_completeness", "done_lawtest"),
    # F6
    ("fig_astrometry", "ref_ofek_masci"), ("fig_astrometry", "data_reprocess"),
    # F10
    ("fig_iraf", "ref_schechter"), ("fig_iraf", "ref_autophot"), ("fig_iraf", "data_iraf"),
    # precision
    ("fig_precision", "ref_honeycutt"), ("fig_precision", "ref_cc_kovacs"),
    ("fig_precision", "data_reprocess"), ("fig_precision", "pend_precision_fix"),
    # investigations / decisions
    ("inv_completeness", "dec_m13_ok"), ("inv_completeness", "dec_synth_ok"),
    ("dec_vv_ladder", "data_synth"), ("inv_detector", "data_injections"),
    ("dec_repro_def", "fig_precision"), ("dec_repro_def", "fig_completeness"),
    # pending chain
    ("pend_lc_period", "ref_cc_kovacs"), ("pend_lc_period", "ref_period"),
    ("pend_lc_period", "data_lc"), ("pend_precision_fix", "ref_cc_kovacs"),
    ("pend_zp_resid", "inv_bfilter"),
    ("done_lawtest", "data_injections"), ("pend_manuscript", "dec_citation_fixes"),
    ("hub_repro", "pend_zp_resid"), ("hub_repro", "pend_cog"), ("hub_repro", "pend_bouguer"),
    ("hub_repro", "pend_lc_period"), ("fig_completeness", "pend_fig_a_layout"),
]


def main():
    ids = {n["id"] for n in NODES}
    bad = [(a, b) for a, b in EDGES if a not in ids or b not in ids]
    assert not bad, f"unknown edge ids: {bad}"
    # card 내부링크 무결성
    for n in NODES:
        for m in re.finditer(r"\[\[([a-z0-9_]+)\|", n["card"]):
            assert m.group(1) in ids, f"{n['id']} links unknown [[{m.group(1)}]]"

    G = {"cats": {str(k): v for k, v in CATS.items()},
         "nodes": NODES, "edges": [list(e) for e in EDGES]}

    # json 정본 — 노드 1줄씩(grep 친화)
    lines = ['{"cats": ' + json.dumps(G["cats"], ensure_ascii=False) + ',', ' "nodes": [']
    lines += [" " + json.dumps(n, ensure_ascii=False) + ("," if i < len(NODES) - 1 else "")
              for i, n in enumerate(NODES)]
    lines += [" ],", ' "edges": ' + json.dumps(G["edges"]) + "}"]
    OUT_JSON.write_text("\n".join(lines), encoding="utf-8")

    # html — KMTNet 엔진 재사용, G 블롭·제목만 교체
    src = ENGINE_SRC.read_text(encoding="utf-8")
    blob = json.dumps(G, ensure_ascii=False)
    out = re.sub(r"const G = \{.*?\};\n", lambda m: "const G = " + blob + ";\n",
                 src, count=1, flags=re.S)
    out = out.replace("<title>지식그래프 — KMTNet 시잉·기상 예보 레퍼런스</title>",
                      "<title>지식그래프 — APEX 논문·재현·검증</title>")
    out = out.replace("레퍼런스 지식그래프 · 노드", "APEX 논문 지식그래프 · 노드")
    OUT_HTML.write_text(out, encoding="utf-8")
    print(f"nodes={len(NODES)} edges={len(EDGES)}")
    print("wrote", OUT_JSON.name, OUT_JSON.stat().st_size, "B /",
          OUT_HTML.name, OUT_HTML.stat().st_size, "B")


if __name__ == "__main__":
    main()
