#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Render MANUSCRIPT markdown -> academic-paper HTML, self-contained.
Math via sub/sup/unicode, citations from references.bib, figures embedded as base64."""
import re, html, base64
from pathlib import Path

ROOT = Path(__file__).resolve().parent   # 워크트리에서도 자기 트리를 읽는다
SRC = ROOT / "논문작업" / "MANUSCRIPT_ko.md"
BIB = ROOT / "references.bib"
FIGDIR = ROOT / "figures"
OUT = ROOT / "MANUSCRIPT_ko_preview.html"   # 정본 자리(2026-07-09부터 이 경로)
OUT_ARTIFACT = ROOT / "MANUSCRIPT_ko_artifact.html"   # 아티팩트 게시용(doctype 없음)

# ---------- figures ----------
# 본문 「그림 N」 -> 실제 파일. 번호는 파이프라인 순서다.
# 파일명이 fig<숫자>_ 패턴이 아닌 새 그림들(실측 주입 완전도 등)도 여기서 직접 잡는다 —
# 예전 자동 매칭은 그 패턴만 찾아서 새 그림을 통째로 놓치고 옛 판을 붙이고 있었다.
FIGFILES = {
    2:  "fig_calibration_step0.png",           # 3.2 검출기 보정(0단계)
    1:  "fig_architecture.png",              # 2.1 작업 흐름·계층
    3:  "fig11_detector.png",                # 3.2 검출기 특성화(PTC)
    4:  "fig12_preproc_crosscheck.png",      # 3.2 ccdproc 교차
    5:  "fig13_cross_instrument.png",        # 3.2 교차기기 보정
    6:  "fig6_qc_validation.png",            # 3.3 프레임 QC
    7:  "fig_detection_threshold.png",       # 3.4 검출 문턱·헛검출 오염
    8:  "fig_completeness_realvssynth.png",  # 3.5 완전도(실측 주입)
    9:  "fig_wcs_engines.png",               # 3.6 astrometric solution
    10:  "fig2_error_model.png",              # 3.7 오차 모형
    11: "fig3_parameter_sweep.png",          # 3.8 민감도
    12: "fig_photometry_crosschecks.png",    # 3.9 SEP + IRAF
    13: "fig_psf_validation.png",            # 3.10 PSF 측광
    14: "fig9_crowded_field.png",            # 3.11 밀집장
    15: "fig_external_validation.png",       # 3.12 PS1 잔차 + CMD
    16: "fig_timeseries_validation.png",     # 3.13 시계열 (SYSREM·PDM)
    17: "fig_lc_yzboo.png",                  # 4    과학 적용(YZ Boo)
    # 옛 그림 2(fig10_calibration.png, 레거시 4패널)는 2026-08-03 사용자 지시로
    # 뺐다. 3.2 의 보정 수치는 본문이 담고, 대체 그림 제작은 fig 세션 몫이다.
}
FIGMAP = {k: FIGDIR / v for k, v in FIGFILES.items() if (FIGDIR / v).exists()}
_missing = [f"{k}:{v}" for k, v in FIGFILES.items() if not (FIGDIR / v).exists()]
if _missing:
    print("[warn] 그림 파일 없음:", ", ".join(_missing))
def fig_uri(p):
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()

# descriptive captions, keyed by pipeline figure number
CAPTIONS = {
 2: r"검출기 보정(0단계): 적용하는 프레임과 그 효과. \textbf{(a–c)} 그 밤의 보정 노출로 APEX가 만든 마스터 bias(중앙값 512 DN, 8장)·dark(60초에 1.0 DN, 8장)·flat(중앙값 1.000, 5장). \textbf{(d)} flat의 가로 프로파일 — 시야를 가로지르는 12\% 감도 기울기가 나눠 내는 대상이다. \textbf{(e–g)} 같은 NGC 6811 $B$ 60초 프레임을 연산마다 본 것: raw, bias·dark를 뺀 뒤, flat으로 나눈 뒤. (f)와 (g)는 같은 회색조라 차이가 flat 보정 몫이다. \textbf{(h)} 각자의 중앙값으로 정규화한 하늘 프로파일 — 좌우 진폭이 12.9\%에서 1.6\%로 줄어든다. 전부 실측 프레임이다(Moravian C3-61000, $2\times2$, 2026-06-11 밤). 이 연산들의 수치 검증은 따로다: 참값 복원(본문), 검출기 상수(그림 3), 세 경로의 절대 산출값 표(그림 4), 카메라 두 대 추가 재현(그림 5).",
 1: r"APEX의 작업 흐름과 계층. 0–7단계는 두 모드가 공유하며, 측광이 끝난 뒤에야 CMD 모드와 LC 모드로 갈라진다. 각 단계는 관측자에게 **결정 하나**를 요구하고(어떤 보정 프레임을 쓸지, 이 프레임을 받아들일지, 검출 문턱을 어디에 둘지, 몇 장에서 보여야 별로 인정할지, 구경을 얼마로 할지) 정해진 경로에 검사 가능한 산출물을 남긴다. 계산은 Qt를 부르지 않는 핵심부에 있고 그래픽 계층은 그것을 부르기만 하므로, 3절의 화면 없는 검증이 곧 화면이 돌리는 코드를 시험한다. 각 단계 아래의 붉은 번호는 그 단계를 검증하는 절이다.",
 3: r"같은 2026-06-11 보정 자료를 APEX, Python `ccdproc`, IRAF `ccdproc`으로 각각 줄여 얻은 검출기 상수. **(a)** 표는 gain $0.681\pm0.014$ e$^-$/ADU, 읽기잡음 $2.35$ e$^-$, 암전류 $0.0077$ e$^-$/s를 도구별로 나란히 보인다. **(b)** 물리량별 막대는 서로 다른 단위를 한 축에 섞지 않고 각각의 실제 값으로 비교한다. 세 결과가 표시 정밀도에서 겹치는 것은 상수 추정의 일치 결과이며, 실제 보정 영상의 절대 산출값은 그림 4의 표에 따로 제시한다. 헤더·제조사·실험실 gain은 이 비교에 넣지 않았다.",
 4: r"같은 NGC 6811 $B$ 60초 자료를 세 경로로 줄인 **절대 산출값**. 각 칸은 전체 프레임의 중앙값 $\pm$ robust $\sigma$($1.4826\times$MAD)이며, APEX에서 뺀 차이값을 그리지 않았다. 행은 master bias, 60초 master dark, 중앙값 1로 정규화한 master flat, bias·dark·flat을 모두 적용한 과학 영상이다. 열은 APEX, 독립적인 Python `ccdproc`, PyRAF로 실행한 IRAF `ccdproc`이다. Python 값은 표시 정밀도에서 APEX와 일치하며, 별도의 화소 잔차 감사 자료는 그림의 수치가 아니다. 2026-06-11 밤, Moravian C3-61000($2\times2$), bias 8·dark 8·flat 5; 우주선·핫픽셀 수리는 산술 비교에서 제외했다.",
 5: r"두 LCO 카메라의 raw를 APEX로 보정해 독립 파이프라인 BANZAI 산출물과 비교. QHY600 CMOS는 전체가 균일한 $+0.06$ e$^-$, 4-앰프 Sinistro CCD는 $\approx0.3\%$ 일치(사분면 패턴은 앰프별 조립의 차이). 보정 산술이 검출기를 넘어 일반화됨을 보인다.",
 6: r"주입 결함으로 만든 44-프레임 밤에서 자동 프레임 QC. 정상 24장을 오탐 없이 통과, 나쁜 시상·밝은 하늘·거짓 헤더 프레임을 모두 검출. 회색 투명도 손실만 놓치는데(균일 변화라 영상 통계에 안 잡힘), 이것이 2단계 측광-QC의 근거다.",
 7: r"검출 문턱값에 따른 헛검출 오염을, 외부 목록 없이 프레임 자신에게서 잰 것. 배경을 뺀 영상의 부호를 뒤집어 같은 검출기를 다시 돌리면 잡음에서 기원한 헛검출이 그대로 세어진다 \citep{serra2012, molino2014}. **(a)** 구상성단 둘·산개성단 둘의 실측 단일 노출 다섯 장($B$·$R$·$g'$). 빈 역삼각형은 상한이다 — 부호를 뒤집은 영상에서 검출이 하나도 안 나온 경우로, $1/N_+$에 찍고 잇는 선에서는 뺐다(이으면 오염이 아니라 $N_+$가 만든 기울기가 그려진다). 큰 빈 기호는 각 프레임 자신의 하한, 곧 오염이 5\% 아래로 유지되는 가장 낮은 문턱값이며 **그 값이 $1.5,\,1.5,\,1.8,\,2.0,\,2.2$로 프레임마다 다르다.** M13의 두 곡선은 **같은 성단을 같은 밤에 찍은 두 필터**라 이 갈림을 대상 탓으로 돌릴 수 없음을 보인다. 기본값 $3.2\sigma$(세로선)에서는 다섯 장 모두 오염 2\% 이하다. 음영은 다섯이 함께 무너지는 $1.5$–$1.2\sigma$ 구간으로, 잡음이 최소 연결 면적에 걸쳐 문턱을 넘을 확률은 문턱값만의 함수이므로 프레임의 별 개수와 무관하다. **(b)** Gaia DR3 대조 검증. 가로축은 Gaia와 짝지어지지 않은 실측 헛검출 수, 세로축은 목록 없이 낸 추정값이며, 44점이 네 자릿수에 걸쳐 2배(음영) 안에 든다. 무너지는 구간에서 2.4\%로 맞고, 그 위에서는 추정이 실측보다 크다 — 게이트로서는 안전한 방향이고, Gaia 한계보다 어두운 진짜 별이 가로축에서 헛검출로 세어지는 탓도 있다.",
 8: r"실측 프레임에 인공별을 주입해 잰 검출 완전도. **(a)** 관측 조건을 대표하는 실측 단일 노출 세 장 — 어두운 하늘·양호한 시상(M67 $i$), 밝은 하늘·양호한 시상(NGC 6811 $R$), 밝은 하늘·불량 시상(M13 $V$) — 의 주입 등급 대비 회수율. 점은 구간별 회수율(Wilson 95\% 이항 구간)이고 **등급 공간에는 어떤 함수도 맞추지 않는다**. 50\% 깊이 $m_{50}=17.65,\ 15.64,\ 14.90$ 은 곡선이 0.5를 지나는 지점을 읽은 값이며, 깊이는 방법이 아니라 그 프레임의 하늘 밝기와 시상이 정한다. 회색 파선은 합성 verification 프레임. 아래 스트립은 가장 얕은 프레임에 실제로 주입된 별들의 컷아웃이다. **(b)** 각 별을 기대 peak-화소 S/N으로 재표현하면 깊이가 3.4등급에 걸쳐 벌어졌던 **일곱 프레임이 단일 곡선으로 붕괴**한다. 합동 표본의 오차함수 피팅은 $\mathrm{S/N}_{50}=4.0$, 프레임별 독립 판독은 $4.05\pm0.18$ 로 일치한다.",
 9: r"같은 M13 단일 노출 여덟 장(각 60초, Moravian C3-61000, 2026-05-15; 가로축은 노출번호+필터)을 엔진별 사본 트리에서 세 solve 엔진으로 푼 비교. 4단계 검출 목록은 세 트리가 공유하고, 해는 각 사본의 FITS 헤더에 기록된다. \textbf{(a)} 프레임마다 5×5 화소 격자를 세 엔진의 WCS로 하늘 좌표에 투영해 잰 엔진 쌍 사이 각거리(격자 중앙값). 여덟 장 가운데 여섯(R·V 다섯 장과 B 한 장)에서 세 쌍 모두 0.2–0.4″에 들고, 나머지 B 두 장에서는 내장 해가 외부 두 해로부터 0.8–1.0″ 벌어지되 외부 두 해끼리는 0.06–0.58″로 가깝다(점선은 1화소 = 0.395″). \textbf{(b)} 파이프라인 QC 코드와 독립적으로, 각 엔진의 헤더 WCS로 Gaia DR3 별을 화소에 투영해 검출과 2″ 반경 최근접 매칭한 잔차 RMS. 프레임 중앙값은 내장 0.77 px, astrometry.net 0.93 px, ASTAP 1.10 px다. 내장 솔버는 채점에 쓰인 것과 같은 검출·Gaia 쌍에 해를 맞추므로 이 지표가 유리한 것은 부분적으로 구성 탓이다. 수용 게이트는 내장·ASTAP 여덟 장 전부와 astrometry.net 일곱 장을 통과시켰고(한 장은 99퍼센타일 잔차로 검토 표시), 검출 목록이 없던 예비 실행은 여덟 장 전부 기각하였다.",
 10: r"합성 몬테카를로로 확인한 측광 오차 모형. 잡음 모형은 합성 프레임 생성기와 같고, 측정은 APEX의 측광 루틴 \texttt{phot\_vectorized}가 한다(격자 주입 12등급 단계 $\times$ 36별 $\times$ 12실현, 시드 20260702, 고정 구경 $1.0\times$FWHM). \textbf{(a)} 평균 편향은 $m=17.0$까지 5 mmag 안(실측 2 mmag 이하)이고, 더 어두운 쪽은 잡음 바닥 근처의 비선형 등급 변환이 키운다($m=20.0$ 러닝 중앙값 $-227$ mmag). \textbf{(b)} 구간별 RMS는 CCD 관계식 $\sigma_m=1.0857/\mathrm{SNR}$을 1.9 dex에 걸쳐 따른다(11.6 mmag $\to$ 0.97 등급). \textbf{(c)} SNR $>5$의 pull은 단위 가우스다(평균 $-0.148$, 표준편차 $1.014$, $N=3404$). \textbf{(d)} 보고 $\sigma$ 대 실측 RMS는 $y=x$를 따르고 중앙값 비 $0.98$. 부족은 $m=19.0$에서 17%로 가장 크며, 가장 어두운 $m=20.0$ 구간은 양수 플럭스 선택과 얽혀 3%로 되돌아온다.",
 11: r"파이프라인 파라미터와 관측 조건 훑기. 합성 프레임(APEX 생성기, 시드 20260702, $1024^2$ px, 별 435개, gain 1.5 e$^-$/ADU, 읽기잡음 5 e$^-$) 위에서 점마다 30회 시행 $\times$ 40개 인공별을 생산 코드 경로(경험 PSF 주입, 4단계 검출, 강제 구경 측광)로 재고, 오차막대는 시행 클러스터 부트스트랩의 95% 구간이다. \textbf{(a)} 구경 반지름(FWHM 6.0 px 고정 프레임): 산포는 $0.80\times$FWHM에서 최소 0.051등급이고 $2.5\times$FWHM에서 0.181등급까지 오른다. 회색 띠는 생산 하한 4 px 아래라 도달할 수 없는 반지름이고, 빈 기호는 $0.5\times$FWHM 요청이 하한으로 잘린 점이다. 기본 배율 0.8(파선)이 실측 최소점에 놓인다. \textbf{(b)} 위치 RMSE는 전 구경에서 0.255 px로 동일하다. 중심은 검출에서 오고 측광 구경을 쓰지 않는다. \textbf{(c,d)} 하늘 50$\to$1000 ADU에서 깊이 $m_{50}$ 18.38$\to$17.06등급, 산포 36$\to$88 mmag. \textbf{(e,f)} 시상 3.0$\to$6.0 px에서 깊이 18.07$\to$16.96등급, 산포 42$\to$60 mmag. 그림에 없는 동반 훑기에서 검출 문턱 2.0$\to$6.0$\sigma$는 깊이를 18.36$\to$17.32등급으로 1.04등급 바꾸고, 헛검출은 어느 문턱에서도 시행 30회 합계 1~2건이다.",
 12: r"합성 참값에서 APEX와 독립 구현 엔진 SEP \citep{barbary2016}의 일치. 합성 프레임(생성기와 같은 잡음 모형, 시드 20260702, $800^2$ px, 고립별 150개 중 SNR$>10$ 생존 95개, FWHM 3.5 px, 하늘 150 ADU)에서 두 엔진이 같은 중심에 같은 구경($1.0\times$FWHM $=3.5$ px)을 적분한다. 영점을 한 번 정렬한 뒤 MAD $0.006$등급, RMS $0.010$등급, Pearson $r=0.99995$, 별의 95%가 0.02등급 안에서 맞는다. 같은 화소에 두 독립 엔진이 같은 플럭스를 낸다.",
 13: r"APEX가 raw에서부터 전부 줄인 NGC 6811 $V$ 단일 노출(Moravian C3-61000, 2026-06-11, 30 s)의 실제 별 499개를, 표 3에 맞춘 파라미터로 APEX 강제 측광과 IRAF `phot`(DAOPHOT)이 **같은 고정 좌표에서** 각각 잰 결과: MAD $9.7$ mmag, $r=0.99989$, 구간 중앙값 잔차가 어두운 쪽까지 평평하다 \citep{schechter1993}. 양쪽 모두 재중심을 껐으므로 잔차가 중심 잡기의 차이를 흡수할 수 없다. 이 불일치는 두 코드가 스스로 보고한 형식 오차의 제곱합근($27.8$ mmag)의 3분의 1 남짓이다.",
 14: r"PSF 측광(8단계)과 강제 구경 측광(7단계)의 별 단위 일치. 두 표를 검출 식별자로 짝짓고, 양쪽 모두 신호대잡음비 20을 넘고 혼잡 신뢰불가·포화 플래그가 없는 별에서 $\Delta = m_{\rm PSF}-m_{\rm ap}$ 를 계산해 프레임별 중앙값(EPSF 정규화 오프셋, 절대 보정의 프레임 영점이 흡수)을 뺐다. **(a)** 카메라 세 대(Moravian C3-61000 CMOS, LCO QHY600 CMOS, LCO Sinistro CCD 4-앰프)·자료 여섯 벌·단일 노출 67장·별 90,201개의 분포와 MAD(7–40 mmag). 표시 범위 ±0.2 등급 밖의 별은 성긴 장에서 0.8–2.3%, M13에서 8.9%, 은하가 시야에 있는 NGC 5985 장에서 12.9%로, 구경 쪽이 이웃 별빛·배경 구조를 담는 방향이며 3.11절에서 다룬다. **(b)** QHY600 풀프레임 M45 $2.0^\circ\times1.3^\circ$ 19장의 반경 의존성. 반경별 중앙값은 모서리 구간을 제외하면 14 mmag 폭 안에 있고, 안쪽(정규화 반경 0.3 미만)과 바깥(0.7 초과)의 차이는 4 mmag, 마지막 모서리 구간만 $+28$ mmag(별 단위 MAD와 같은 크기)다. 회색 점은 무작위 12,000개 표본.",
 15: r"두 구상성단(M5·M13) 코어에서 APEX의 두 측광법 — 강제 구경 대 PSF — 의 내부 일치. 최근접 이웃 거리에 대한 중앙 차이가 평평($\pm0.02$–$0.04$등급, 분해 한계 $\sim10$ px까지). Gaia에 의존하지 않는 내부 일관성 시험.",
 16: r"NGC 6811을 Gaia와 독립인 Pan-STARRS 1 \citep{chambers2016}에 교차대조. $B$에서 Gaia 변환 참조가 어두운 쪽으로 $+0.022$등급 흐르지만 APEX 자체는 PS1 대비 평평($+0.010$) — 어두운 쪽 편차는 APEX가 아니라 Gaia BP의 알려진 결함 \citep{riello2021}.",
 17: r"NGC 6811의 Johnson 색-등급도($V$ 대 $B-V$, 별 1921개). APEX 지상 측광이 Gaia 변환 우주 기반 참조와 주계열 능선 $19$ mmag로 일치하며, 독립 PS1 계에서도 같은 형태. CMD 산출물 자체를 검증(이소크론 맞추기는 별개).",
 18: r"자체 구현한 시계열 모듈 둘을 주입한 참값으로 확인한 결과. 검증은 LC 모드의 단계 창이 호출하는 것과 같은 진입점을 쓰고 고정 시드에서 만든 합성 계열만 쓴다. **(a)** 위상 분산 최소화(PDM)의 주기 복원. 4절과 같은 관측 배치(하룻밤 5.2시간·80점) 위에 기본 진동과 배진동을 더한 비대칭 신호를 얹고, 주기 0.06–0.20일과 잡음 5–40 mmag의 각 칸마다 잡음 실현 12개를 돌려 $|P_{\rm rec}-P_{\rm true}|/P_{\rm true}$의 중앙값을 표시하였다. 잡음 20 mmag 이하에서 중앙값 1.4\%이고, 오차를 키우는 것은 잡음보다 주기다. 기선 0.217일이 주기 0.20일을 한 번밖에 담지 못하기 때문이며, 알고리즘이 아니라 기선의 한계다. **(b)** 같은 계열에 라이브러리의 Lomb–Scargle \citep{vanderplas2018} 을 함께 돌린 결과. 288개 계열이 일대일 선을 따르므로 자체 구현한 PDM이 다른 문제를 풀고 있지 않음을 확인한다. 어긋나는 점은 기선이 한 주기를 겨우 담는 0.16–0.22일 구간에 몰려 있다. **(c)** SYSREM \citep{tamuz2005}. 별 60개·프레임 120장에 프레임별 투과율과 별별 감도의 곱으로 계통 성분을 주입하고 한 별에만 주기 0.104092일·진폭 0.21등급의 변광을 넣었다. APEX는 성분을 **비교성만으로** 풀고 대상별에 적용한다(파란 점, 진폭 109\% 유지). 대상별을 성분 추출에 넣으면 비교성 산포는 똑같이 줄지만 변광 진폭은 0.4\%만 남는다(보라 x). 추세 제거를 산포 감소만으로 평가하면 안 되는 이유다.",
 19: r"LC 모드의 종단(end-to-end) 과학 산출물: 고진폭 $\delta$ Scuti 별 YZ Boötis를 raw 프레임에서 APEX만으로 줄였다. **(a)** 하룻밤(2026-03-28, $r$ 밴드, 5.2시간에 걸친 80점)을 문헌 주기로 접으면 출판된 톱니 곡선이 재현되고 진폭은 pk-pk $0.39$ 등급이다. **(b)** Lomb–Scargle 주기도. 단일밤 최고 봉우리는 $0.1046$ 일로 문헌값 $0.10409$ 일(파선)과 0.5\% 차이지만, 하루 간격의 두 밤을 병합하면 최고 봉우리가 $+1$ 주기/일 alias인 $0.0946$ 일로 옮겨가고 참 주기는 부봉우리로만 남는다. 이는 관측 창(window)의 성질이지 파이프라인의 것이 아니며, 숨기지 않고 보인다."
}

# Final 17-figure sequence.  Keep these captions factual and let the prose carry
# interpretation; each image itself contains its data provenance.
CAPTIONS.update({
  1: r"APEX의 작업 흐름과 소프트웨어 계층. 0–7단계는 CMD와 LC 모드가 공유하며 측광 뒤에 분기한다. 그래픽 계층은 Qt와 분리된 계산 핵심부를 호출하므로, 3절의 화면 없는 시험은 GUI가 사용하는 것과 같은 계산 경로를 거친다. 상자 아래 번호는 해당 단계의 검증 또는 적용 절이다.",
  2: r"0단계 보정 자료와 보정 효과. **(a)** master bias 영상과 **(b)** 그 화소값 분포·중앙값·산포, **(c)** master dark 영상과 **(d)** 그 분포·상위 백분위, **(e)** master flat 영상과 **(f)** flat 가로 프로파일을 짝지어 보인다. **(g)** raw NGC 6811 $B$ 과학 영상, **(h)** bias·dark 제거와 flat 나눗셈을 끝낸 영상이다. (h)의 작은 삽입도에서 보정 전후 하늘 프로파일을 각 중앙값으로 정규화해 비교한다. master는 bias 8장, 60초 dark 8장, flat 5장으로 만들었고, 전부 Moravian C3-61000의 2026-06-11 실측 프레임이다. 하늘 배경 좌우 진폭은 12.9\%에서 1.6\%로 줄었다.",
  3: r"광자전달곡선으로 측정한 검출기 상수. **(a)** flat 쌍의 신호-분산 관계에서 gain $0.681\pm0.014$ e$^-$/ADU를 얻었고, **(b)** dark 노출 사다리의 기울기에서 암전류 $0.0077$ e$^-$/s를 얻었다. (b) 아래 잔차는 선형 적합의 $R^2=0.998$을 보인다. 읽기잡음은 $2.35$ e$^-$이며 FITS 헤더 EGAIN은 측정값보다 약 16배 작다.",
  4: r"같은 NGC 6811 $B$ 60초 자료를 세 경로로 줄인 **절대 산출값**. 각 칸은 전체 프레임의 중앙값 $\pm$ robust $\sigma$($1.4826\times$MAD)이며, APEX에서 뺀 차이값을 그리지 않았다. 행은 master bias, 60초 master dark, 중앙값 1로 정규화한 master flat, bias·dark·flat을 모두 적용한 과학 영상이다. 열은 APEX, 독립적인 Python `ccdproc`, PyRAF로 실행한 IRAF `ccdproc`이다. Python 값은 bias $512\pm1.483$, dark $1\pm2.224$, flat $1\pm0.04428$, 완전 보정 $633.5\pm31.39$ DN으로 APEX와 같은 수치가 표시되며, IRAF 값은 $512\pm1.483$, $1\pm2.224$, $1.001\pm0.04431$, $633\pm31.37$ DN이다. 별도의 화소 잔차 감사 자료는 그림의 수치가 아니다. 2026-06-11 밤, Moravian C3-61000($2\times2$), bias 8·dark 8·flat 5; 우주선·핫픽셀 수리는 산술 비교에서 제외했다.",
  5: r"서로 다른 두 LCO 카메라에서 APEX 보정 영상과 BANZAI 산출물을 비교한 결과. **위 행**은 0.4 m의 QHY600 단일 증폭기 CMOS와 Proxima Cen 장이며, **아래 행**은 1 m의 Sinistro 네 증폭기 CCD와 NGC 5985 장이다. QHY600 차이는 거의 균일한 $+0.06$ e$^-$이고, Sinistro 화소값은 약 0.3\% 범위에서 일치하지만 증폭기 사분면 구조가 남는다.",
  6: r"주입 결함으로 구성한 44장 자료에서 자동 프레임 QC의 판정. **(a)** 주입한 정상·불량 종류와 QC 판정의 혼동행렬: 정상 24장은 모두 PASS, 나쁜 시상은 FAIL, 밝은 하늘과 헤더와 실제 읽기잡음이 다른 영상은 REVIEW로 분류됐다. **(b)** QC가 사용하는 진단 평면에서 균일한 0.7등급 투과율 손실은 영상 통계만으로 분리되지 않는다.",
 7: r"검출 문턱값에 따른 헛검출 오염. **(a)** 실제 영상의 부호를 뒤집은 검출로 얻은 프레임별 안전 하한과 기본값 $3.2\sigma$의 오염률(다섯 영상 모두 2\% 이하). **(b)** Gaia DR3와의 독립 대조에서 목록 없이 추정한 헛검출 수가 44개 자료점에서 실측값의 2배 이내에 드는지 보인다. 부호반전은 대칭 잡음을 세는 상한이며 우주선·핫픽셀 같은 양의 결함은 별도 단계의 대상이다.",
 8: r"실측 영상 일곱 장에 경험적 PSF로 인공별을 주입한 검출 완전도. **(a)** 세 대표 단일 노출의 등급 공간 회수율과 $m_{50}$ (14.90–17.65등급); 점은 Wilson 95\% 구간이며 회색 파선은 합성 verification이다. 아래 스트립은 가장 얕은 M13 $V$ 프레임의 주입 컷아웃이다. **(b)** 이 그림의 주 결과로, 모든 주입별을 기대 peak-pixel S/N으로 바꾼 일곱 프레임이 한 곡선으로 모이며 합동 오차함수 적합은 $\mathrm{S/N}_{50}=4.0$, 프레임별 판독은 $4.05\pm0.18$이다.",
  9: r"M13 60초 영상 8장에서 built-in, ASTAP, astrometry.net 해를 비교했다. **(a)** 동일한 5×5 화소 격자를 각 WCS로 투영했을 때 세 해 사이의 각거리. **(b)** 같은 Gaia DR3 대조 절차로 얻은 검출 위치 잔차. 내장 quad 솔버가 기본 기준이고 ASTAP와 astrometry.net은 선택한 외부 엔진의 교차검사다.",
 10: r"합성 몬테카를로 자료로 점검한 강제 구경측광 오차 모형. **(a)** 평균 편향은 $m\le17$에서 5 mmag 이하이고 어두운 쪽에서 등급 변환 비선형성이 커진다. **(b)** 구간별 RMS는 CCD 관계식 $1.0857/\mathrm{SNR}$를 따른다. **(c)** S/N $>5$ 표본의 pull 평균과 표준편차는 각각 $-0.148$, 1.014다. **(d)** 보고 오차 대 경험 RMS의 중앙비는 0.98이며 $m=19$에서 부족이 17\%로 가장 크다. 이 그림이 측광 오차 모형의 주 검증이고 그림 8-b는 검출 완전도의 주 검증이다.",
 11: r"파라미터와 관측 조건에 대한 민감도. **(a)** FWHM 6 px에서 구경 반지름의 오차 MAD와 위치 RMSE: 최소는 $0.80\times$FWHM, 위치 RMSE 중앙값은 0.255 px이다. **(b)** 하늘 배경 50→1000 ADU에서 $m_{50}$은 18.38→17.06등급, MAD는 36→88 mmag이다. **(c)** FWHM 3→6 px에서 $m_{50}$은 18.07→16.96등급, MAD는 42→60 mmag이다.",
 12: r"독립 측광 엔진과의 대조. **(a)** 합성 영상에서 같은 중심과 구경을 쓴 APEX–SEP 영점 정렬 잔차(MAD 6.0 mmag, RMS 9.7 mmag). **(b)** APEX로 보정한 NGC 6811 $V$ 영상의 고정 중심 499개에서 APEX–IRAF/DAOPHOT 잔차(MAD 9.7 mmag). 각 패널의 실선은 보정 함수가 아닌 구간별 중앙값이다.",
 13: r"세 카메라, 여섯 자료군, 67장, 90,201개 별에서 PSF 측광과 강제 구경측광의 차이. **(a)** 세 카메라와 여섯 자료군의 별 단위 잔차 분포(MAD 7–40 mmag). **(b)** QHY600 M45 광시야에서 정규화 반경에 따른 잔차; 마지막 모서리 구간만 중앙부보다 $+28$ mmag이다.",
 14: r"구상성단 M5와 M13 중심부에서 PSF 측광과 강제 구경측광의 차이를 최근접별 거리로 나타냈다. **(a)** 두 성단의 코어 밀도 대비 배경 밀도. **(b)** 거리 구간별 측광 차이; 약 10 px의 분해 하한부터 수백 px까지 중앙값은 0–0.04등급 이내다. 외부 절대 정확도나 10 px 아래 미분해 겹침은 검증하지 않는다.",
 15: r"NGC 6811의 외부 목록 및 CMD 대조. **(a)** $B$에서 PS1과 Gaia 변환 참조에 대한 밝기별 잔차. **(b)** $V$에서 같은 대조. **(c)** APEX와 Gaia 변환값으로 만든 Johnson CMD 1,917개의 주계열 ridgeline RMS 19 mmag. 이 그림은 CMD 산출물을 검증하며 등시선 모수의 유일성을 검증하지 않는다.",
 16: r"합성 시계열로 점검한 PDM과 SYSREM. **(a)** PDM의 주기 복원(0.06–0.20일, 잡음 5–40 mmag). **(b)** 같은 288개 계열에서 자체 PDM과 Lomb–Scargle 결과의 일치. **(c)** SYSREM: 비교성만으로 공통 성분을 구하면 주입 변광 진폭의 109\%를 유지하지만 변광성을 성분 추출에 넣으면 0.4\%만 남는다.",
 17: r"YZ Boo의 LC 모드 종단 산출물. **(a)** 2026-03-28 $r$ 대역 80장, 5.2시간 자료를 문헌 주기 0.10409일로 위상 접은 광도곡선. **(b)** Lomb–Scargle 주기도: 단일 밤 최고 봉우리는 0.1046일로 문헌값과 0.5\% 차이지만 여러 밤을 합치면 관측창의 $+1$ cycle/day alias가 가장 높아진다."
})

# Figure 3 is a cross-tool detector-constant comparison.  Keep this override
# after the legacy caption table above so the preview cannot silently show the
# former single-pipeline PTC description.
CAPTIONS[3] = r"같은 2026-06-11 보정 자료를 APEX, Python ccdproc, IRAF ccdproc으로 각각 줄여 얻은 검출기 상수. **(a)** 표는 gain $0.681\pm0.014$ e$^-$/ADU, 읽기잡음 $2.35$ e$^-$, 암전류 $0.0077$ e$^-$/s를 도구별로 나란히 보인다. **(b)** 물리량별 막대는 서로 다른 단위를 한 축에 섞지 않고 각각의 실제 값으로 비교한다. 세 결과가 표시 정밀도에서 겹치는 것은 상수 추정의 일치 결과이며, 실제 보정 영상의 절대 산출값은 그림 4의 표에 따로 제시한다. 헤더·제조사·실험실 gain은 이 비교에 넣지 않았다."

# ---------- bibliography ----------
def parse_bib(text):
    labels = {}
    for m in re.finditer(r"@\w+\s*\{\s*([^,]+),(.*?)\n\}", text, re.S):
        key = m.group(1).strip(); body = m.group(2)
        au = re.search(r"author\s*=\s*[{\"](.+?)[}\"]\s*,?\s*\n", body, re.S)
        yr = re.search(r"year\s*=\s*[{\"]?\s*(\d{4})", body)
        year = yr.group(1) if yr else ""
        if au:
            authors = re.split(r"\s+and\s+", au.group(1).strip())
            def last(a):
                a = a.strip().strip("{}")
                if "," in a: return a.split(",")[0].strip().strip("{}")
                if a.endswith("Collaboration"): return a
                parts = a.split(); return parts[-1] if parts else a
            l0 = last(authors[0])
            if len(authors) == 1: lab = f"{l0} {year}"
            elif len(authors) == 2: lab = f"{l0} & {last(authors[1])} {year}"
            else: lab = f"{l0} et al. {year}"
        else:
            lab = f"{key} {year}".strip()
        labels[key] = lab.strip()
    return labels
LAB = parse_bib(BIB.read_text(encoding="utf-8")) if BIB.exists() else {}

# bib 의 LaTeX 악센트를 실제 글자로. 안 풀면 "Alarc{\'o}n et al. 2023" 이 그대로 찍힌다.
_ACC = {"'": "\u0301", "`": "\u0300", '"': "\u0308", "^": "\u0302", "~": "\u0303",
        "c": "\u0327", "v": "\u030C", "u": "\u0306", "=": "\u0304", ".": "\u0307"}
def de_tex(s):
    import unicodedata
    def one(m):
        acc, ch = m.group(1), m.group(2)
        return unicodedata.normalize("NFC", ch + _ACC[acc]) if acc in _ACC else ch
    s = re.sub(r"\{?\\([\"'`^~=.]|[cvu])\s*\{(\w)\}\}?", one, s)   # {\'o} · \'{o}
    s = re.sub(r"\\([\"'`^~=.])\s*(\w)", one, s)                   # \'o
    s = s.replace(r"\ss", "ß").replace(r"\o", "ø").replace(r"\aa", "å")
    return s.replace("{", "").replace("}", "")

LAB = {k: de_tex(v) for k, v in LAB.items()}
def cite_label(k): return LAB.get(k.strip(), k.strip())

def repl_citations(s):
    def paren(m):
        pre, post, keys = m.group(1), m.group(2), m.group(3)
        inner = "; ".join(cite_label(k) for k in keys.split(","))
        if pre: inner = f"{pre} {inner}"
        if post: inner = f"{inner}, {post}"
        return f"({inner})"
    def plain(m): return "; ".join(cite_label(k) for k in m.group(1).split(","))
    s = re.sub(r"\\citep(?:\[([^\]]*)\])?(?:\[([^\]]*)\])?\{([^}]+)\}", paren, s)
    s = re.sub(r"\\cite(?:t|alt)\{([^}]+)\}", plain, s)
    return s

# ---------- math ----------
GREEK = {r"\\simeq":"≃", r"\\sigma":"σ", r"\\tau":"τ", r"\\theta":"θ", r"\\mu":"μ", r"\\delta":"δ",
         r"\\Delta":"Δ", r"\\alpha":"α", r"\\beta":"β", r"\\chi":"χ", r"\\pi":"π",
         r"\\approx":"≈", r"\\times":"×", r"\\sim":"∼", r"\\pm":"±", r"\\geq":"≥",
         r"\\leq":"≤", r"\\ge":"≥", r"\\le":"≤", r"\\cdot":"·", r"\\to":"→",
         r"\\infty":"∞", r"\\propto":"∝", r"\\equiv":"≡"}
def render_math(x):
    for pat, rep in GREEK.items(): x = re.sub(pat, rep, x)
    x = re.sub(r"\\(?:mathrm|rm|text|mathbf|bf|it)\s*\{([^{}]*)\}", r"\1", x)
    x = re.sub(r"\\(?:mathrm|rm|text|mathbf|bf|it)\b\s*", "", x)
    x = re.sub(r"\\(log|max|min|sec|exp|sin|cos|ln|arcsin)", r"\1", x)
    x = x.replace("\\,", " ").replace("\\;", " ").replace("\\!", "").replace("\\ast", "*")
    x = re.sub(r"\\frac\s*\{([^{}]*)\}\s*\{([^{}]*)\}", r"(\1)/(\2)", x)
    x = re.sub(r"\\sqrt\s*\{([^{}]*)\}", r"√(\1)", x)
    x = re.sub(r"\^\{([^{}]*)\}", r"<sup>\1</sup>", x)
    x = re.sub(r"_\{([^{}]*)\}", r"<sub>\1</sub>", x)
    x = re.sub(r"\^(\*|[-+]|\w)", lambda m: f"<sup>{'−' if m.group(1)=='-' else m.group(1)}</sup>", x)
    x = re.sub(r"_(\w)", r"<sub>\1</sub>", x)
    x = x.replace("{", "").replace("}", "").replace("$", "")
    x = re.sub(r"\\[a-zA-Z]+", "", x)
    x = re.sub(r"  +", " ", x).strip()
    return f'<span class="math">{x}</span>'

def inline(s):
    holds = []
    def hold(t): holds.append(t); return f"\0{len(holds)-1}\0"
    s = re.sub(r"\\texttt\{([^}]*)\}", r"`\1`", s)
    s = re.sub(r"`([^`]+)`", lambda m: hold(f'<code>{html.escape(m.group(1))}</code>'), s)
    s = re.sub(r"\$([^$]+)\$", lambda m: hold(render_math(m.group(1))), s)
    s = repl_citations(s)
    s = s.replace("\\%", "%").replace("\\&", "&").replace("\\_", "_").replace("\\#", "#").replace("\\$", "$")
    s = html.escape(s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"(?<!\*)\*(?!\*)([^*]+)\*(?!\*)", r"<em>\1</em>", s)
    # [텍스트](대상) — 파일 링크는 걸 곳이 없으니 텍스트만 남긴다
    s = re.sub(r"\[([^\]\n]+)\]\(([^)\s]+)\)",
               lambda m: m.group(1) if not m.group(2).startswith(("http://", "https://"))
               else f'<a href="{m.group(2)}">{m.group(1)}</a>', s)
    s = re.sub(r"\0(\d+)\0", lambda m: holds[int(m.group(1))], s)
    return s

# ---------- parse ----------
lines = [ln for ln in SRC.read_text(encoding="utf-8").splitlines()
         if not re.match(r"\s*<!--.*-->\s*$", ln)]
out = []
para = []
sub_title = ""
sub_figs = []   # figure numbers referenced in current subsection (in order)
emitted = set()

def emit_fig_numbers(numbers):
    """Emit figures at the paragraph that introduces them.

    The previous renderer accumulated all references until the subsection
    ended.  That made a numerical sentence in §3.2 and its evidence figure
    drift apart by several pages.  Keep the old subsection fallback for
    unusual references, but place ordinary paragraph references immediately.
    """
    if not re.match(r"^\d+(\.\d+)?[\.\s]", sub_title):
        return
    title = re.sub(r"^\d+(\.\d+)?\.?\s*", "", sub_title)
    for n in numbers:
        if n in FIGMAP and n not in emitted:
            emitted.add(n)
            body = inline(CAPTIONS[n]) if n in CAPTIONS else html.escape(title)
            cap = f"<b>그림 {n}.</b> {body}"
            out.append(f'<figure><img alt="그림 {n}" src="{fig_uri(FIGMAP[n])}">'
                       f'<figcaption>{cap}</figcaption></figure>')

def flush_para(buf):
    global sub_figs
    if buf:
        txt = " ".join(buf).strip()
        if txt:
            cls = ""
            if txt.startswith("*[") or txt.startswith("["): cls = ' class="pending"'
            elif re.match(r"\*?영문 제출본", txt): cls = ' class="docnote"'
            out.append(f"<p{cls}>{inline(txt)}</p>")
            refs = []
            for mm in re.finditer(r"그림\s*(\d+)", txt):
                n = int(mm.group(1))
                if n not in refs:
                    refs.append(n)
            if refs:
                emit_fig_numbers(refs)
                sub_figs = [n for n in sub_figs if n not in refs]
    return []

def emit_figs():
    global sub_figs
    # 번호가 붙은 절이면 그림을 받는다. 하위절(3.6)뿐 아니라 최상위 절(4. 과학 적용)도 —
    # 예전 조건은 "N.M" 만 허용해서 §4 의 그림이 통째로 빠졌다.
    if re.match(r"^\d+(\.\d+)?[\.\s]", sub_title):
        title = re.sub(r"^\d+(\.\d+)?\.?\s*", "", sub_title)
        emit_fig_numbers(sub_figs)
    sub_figs = []

i = 0
while i < len(lines):
    ln = lines[i]
    if re.match(r"^#\s+", ln):
        para = flush_para(para); emit_figs()
        out.append(f'<h1 class="title">{inline(ln[2:].strip())}</h1>')
    elif re.match(r"^##\s+", ln):
        para = flush_para(para); emit_figs()
        t = ln[3:].strip()
        out.append(f'<h2{" id=abstract" if t.startswith("초록") else ""}>{inline(t)}</h2>')
        sub_title = t
    elif re.match(r"^###\s+", ln):
        para = flush_para(para); emit_figs()
        t = ln[4:].strip()
        out.append(f"<h3>{inline(t)}</h3>")
        sub_title = t
    elif re.match(r"^####\s+", ln):
        para = flush_para(para)
        out.append(f"<h4>{inline(ln[5:].strip())}</h4>")
    elif ln.strip() == "---":
        para = flush_para(para)
    elif ln.strip().startswith("|"):
        para = flush_para(para)
        tbl = []
        while i < len(lines) and lines[i].strip().startswith("|"):
            tbl.append(lines[i]); i += 1
        i -= 1
        rows = [[c.strip() for c in r.strip().strip("|").split("|")] for r in tbl]
        body_rows = [r for r in rows if not all(re.match(r"^:?-{2,}:?$", c) for c in r)]
        head, data = body_rows[0], body_rows[1:]
        h = "".join(f"<th>{inline(c)}</th>" for c in head)
        b = "".join("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>" for r in data)
        out.append(f'<div class="tw"><table><thead><tr>{h}</tr></thead><tbody>{b}</tbody></table></div>')
    elif ln.strip().startswith(">"):
        para = flush_para(para)
        q = []
        while i < len(lines) and lines[i].strip().startswith(">"):
            q.append(lines[i].strip()[1:].strip()); i += 1
        i -= 1
        out.append(f"<blockquote>{inline(' '.join(q))}</blockquote>")
    elif ln.strip() == "":
        para = flush_para(para)
    else:
        para.append(ln.strip())
        for mm in re.finditer(r"그림 (\d+)", ln):
            n = int(mm.group(1))
            if n not in sub_figs: sub_figs.append(n)
    i += 1
flush_para(para); emit_figs()
BODY = "\n".join(out)

# ---------- 논문 체재로 조립: 표지 · 목차 · 절 단위 페이지 ----------
_m = re.search(r'<h1 class="title">(.*?)</h1>', BODY, re.S)
DOC_TITLE = re.sub(r"<[^>]+>", "", _m.group(1)).strip() if _m else "APEX"
BODY = re.sub(r'<h1 class="title">.*?</h1>', "", BODY, count=1, flags=re.S)
_note = re.search(r'<p class="docnote">(.*?)</p>', BODY, re.S)
DOC_NOTE = _note.group(1).strip() if _note else ""
BODY = re.sub(r'<p class="docnote">.*?</p>', "", BODY, count=1, flags=re.S)

# 목차 — 본문의 절 제목과 그림 목록에서 만든다
_toc = []
for _h in re.finditer(r'<(h2|h3)[^>]*>(.*?)</\1>', BODY, re.S):
    _t = re.sub(r"<[^>]+>", "", _h.group(2)).strip()
    if _t and not _t.startswith("그림"):
        _toc.append((_h.group(1), _t))
def _plain(s):
    """태그를 벗기고 엔티티를 한 번 되돌린다 — inline() 이 이미 이스케이프했으므로
    여기서 다시 escape 하면 '&amp;' 가 화면에 그대로 보인다."""
    return html.unescape(re.sub(r"<[^>]+>", "", s)).strip()

# 쪽번호 칸(.pn)은 조판 전에 미리 넣어 둔다. 조판이 끝난 뒤 칸을 추가하면
# 줄바꿈이 달라져 이미 배치한 지면이 넘칠 수 있다.
_toc_html = "".join(
    f'<li class="{lvl} go"><span class="tx">{html.escape(_plain(t))}</span>'
    f'<span class="pn"></span></li>' for lvl, t in _toc)
_fig_html = "".join(
    f'<li class="go"><span class="tx">'
    f'{html.escape("그림 " + str(n) + ". " + _plain(inline(CAPTIONS[n]))[:110])}'
    f'</span><span class="pn"></span></li>'
    for n in sorted(FIGMAP) if n in CAPTIONS)

# ── 표제면: A&A 667, A62 (AutoPhOT 논문) 판면을 참고하되, 투고 후보는 A&C로 표시 ──
#    저널 머리 → 제목 → 저자 → 소속 → 접수일 → 초록(전폭) → Key words → 2단 본문
_abs = re.search(r'<h2[^>]*\bid=["\']?abstract["\']?[^>]*>.*?</h2>(.*?)(?=<h2)', BODY, re.S)
ABS_HTML, KW_HTML = "", ""
if _abs:
    ABS_HTML = _abs.group(1).strip()
    BODY = BODY.replace(_abs.group(0), "", 1)
    # 원고가 이미 가진 핵심어 줄을 그대로 쓴다 (따로 지어내지 않는다)
    _kw = re.search(r'<p>\s*<strong>핵심어[^<]*</strong>(.*?)</p>', ABS_HTML, re.S)
    if _kw:
        KW_HTML = f'<p class="kw"><b>핵심어.</b>{_kw.group(1)}</p>'
        ABS_HTML = ABS_HTML.replace(_kw.group(0), "", 1).strip()

TITLEBLOCK = f'''<div class="titleblock">
  <div class="jrnl">
    <div class="jl">Astronomy and Computing<br><span>투고 준비 원고 · 국문 검토용</span></div>
    <div class="jr">A&amp;C</div>
  </div>
  <h1 class="ptitle">{html.escape(DOC_TITLE)}</h1>
  <p class="pauth">저자 미기재 (투고 전 확정)</p>
  <p class="paff">한국천문연구원 · 소형망원경 측광 파이프라인<br>
     <span class="pmail">2026erpcosmos@gmail.com</span></p>
  <p class="pdate">{DOC_NOTE}</p>
</div>
<div class="absblock">
  <div class="abshead">초록</div>
  {ABS_HTML}
  {KW_HTML}
</div>'''

# 목차는 A&A 판면에 없다. 검토 편의를 위해 기본은 켜 두되 한 줄로 끌 수 있게 한다.
WANT_TOC = True
TOCBLOCK = (f'''<div class="tocblock">
  <div class="toc-h">차례</div>
  <ol class="toc-list">{_toc_html}</ol>
  <div class="toc-h">그림</div>
  <ol class="toc-figs">{_fig_html}</ol>
</div>''' if WANT_TOC else "")

BODY = TITLEBLOCK + TOCBLOCK + f'<div class="flow">{BODY}</div>'

# ============================================================================
#  판면 — A&A 667, A62 (2022) 를 자로 삼는다.
#  A4 210×297mm = 96dpi 에서 794×1123px. 좌우 여백 17mm, 위 20mm, 아래 16mm,
#  단 간격 6mm, 단 너비 ≈ 85mm. 본문 9pt(12px) Times, 제목은 산세리프.
# ============================================================================
CSS = r"""
:root{
  color-scheme:light;
  --pw:794px; --ph:1123px;          /* A4 @96dpi */
  --mx:64px; --mt:76px; --mb:60px;  /* 판면 여백 */
  --gut:24px;                        /* 단 사이 */
  --ink:#111; --muted:#4b5057; --link:#1a3d8f;
  --serif:"Times New Roman",Times,"Noto Serif KR","Batang",serif;
  --sans:"Helvetica Neue",Arial,"Malgun Gothic","Noto Sans KR",sans-serif;
}
*{box-sizing:border-box}
html,body{margin:0;background:#8f9298;}
body{font-family:var(--serif);font-size:11.2px;line-height:1.32;color:var(--ink);
  font-variant-numeric:lining-nums;-webkit-font-smoothing:antialiased;}

/* ── 지면 ── */
#src{display:none;}
#stage{overflow:hidden;}
/* 판면 모드: PDF 뷰어처럼 A4 쪽을 세로로 이어 붙인다. 휠은 자연스럽게
   현재 쪽의 끝에서 다음 쪽으로 이어지고, ‹›·PageUp/Down은 해당 쪽으로 이동한다. */
body.paged{overflow:hidden;}
body.paged #stage{display:block;overflow:auto;min-height:100vh;}
body.paged #sizer{display:block;}
body.paged #book{display:block;}
body.paged .page{display:block;}
body.paged .page.active{display:block;}
#sizer{position:relative;}
#book{transform-origin:top left;width:var(--pw);}
.page{position:relative;width:var(--pw);height:var(--ph);background:#fff;
  margin:0 0 14px;box-shadow:0 2px 10px rgba(0,0,0,.35);overflow:hidden;
  break-after:page;page-break-after:always;break-inside:avoid;page-break-inside:avoid;}
.page:last-child{break-after:auto;page-break-after:auto;}
.page.singlecol .cols{display:block!important;flex:none!important;}
.page.singlecol .col{display:block!important;width:100%!important;flex:none!important;overflow:visible;}
.page.singlecol .col + .col{display:none!important;}
.pinner{position:absolute;left:var(--mx);right:var(--mx);top:var(--mt);bottom:var(--mb);
  display:flex;flex-direction:column;}
.span:not(:empty){margin-bottom:10px;}
.span.bot:not(:empty){margin:10px 0 0;}
.cols{flex:1 1 auto;min-height:0;display:flex;gap:var(--gut);}
.col{flex:1 1 0;min-width:0;overflow:hidden;}
.run{position:absolute;left:var(--mx);right:var(--mx);top:34px;font-size:10px;
  color:#111;text-align:center;}
.folio{position:absolute;left:var(--mx);right:var(--mx);bottom:32px;font-size:10px;color:#111;}
.folio.r{text-align:right;} .folio.l{text-align:left;}

/* ── 본문 조판 ── */
/* 좁은 단에서 한글 양끝맞춤은 어절 사이가 벌어지므로, 한글은 어디서나 줄바꿈을
   허용(word-break:normal)해 균일하게 채운다. 국문 학술지 조판 관례다. */
p{margin:0;text-align:justify;word-break:normal;overflow-wrap:break-word;text-indent:1.1em;}
p.pcont{text-indent:0;}                 /* 단 경계에서 쪼개진 문단의 뒤토막 */
p.psplit{text-align-last:justify;}      /* 앞토막은 마지막 줄도 양끝맞춤 */
h2 + p, h3 + p, h4 + p, figure + p, .tw + p, blockquote + p{text-indent:0;}
/* 제목·캡션은 어절 중간에서 끊기면 안 된다. 본문 단만 음절 단위 줄바꿈을 허용한다. */
h2,h3,h4,.ptitle,.abshead,.toc-h,figcaption,thead th{word-break:keep-all;}
h2{font-family:var(--sans);font-size:12.4px;font-weight:700;line-height:1.25;
  margin:11px 0 4px;text-indent:0;}
h3{font-family:var(--sans);font-size:11.2px;font-weight:700;line-height:1.25;
  margin:9px 0 3px;text-indent:0;}
h4{font-family:var(--sans);font-size:11.2px;font-weight:700;font-style:italic;
  margin:7px 0 2px;text-indent:0;}
.col > *:first-child{margin-top:0;}
strong{font-weight:700;}
em{font-style:italic;}
code{font-family:"Courier New",monospace;font-size:.92em;}
.math{font-family:"Cambria Math","Times New Roman",serif;white-space:nowrap;}
.math sub,.math sup{font-size:.72em;}
blockquote{margin:5px 0 5px 10px;padding-left:8px;border-left:1.5px solid #999;
  color:var(--muted);font-size:11.4px;text-align:justify;}
p.pending{color:var(--muted);font-style:italic;text-indent:0;}
a{color:var(--link);text-decoration:none;}

/* ── 표제면 ── */
.jrnl{display:flex;justify-content:space-between;align-items:flex-start;
  font-size:10.5px;line-height:1.35;margin:0 0 26px;}
.jl span{color:var(--link);}
.jr{font-family:var(--sans);font-weight:700;font-size:15px;letter-spacing:.02em;
  border-bottom:2px solid var(--ink);padding-bottom:1px;}
.ptitle{font-family:var(--sans);font-size:22px;line-height:1.28;font-weight:700;
  text-align:center;margin:0 0 16px;text-wrap:balance;}
.pauth{font-size:13px;text-align:center;text-indent:0;margin:0 0 14px;}
.paff{font-size:10.5px;text-align:center;text-indent:0;line-height:1.45;margin:0 0 12px;}
.pmail{font-family:"Courier New",monospace;font-size:9.6px;}
.pdate{font-size:10px;text-align:center;text-indent:0;color:var(--muted);
  font-style:italic;margin:0 0 18px;}
.absblock{margin:0 44px 6px;}
.abshead{font-family:var(--sans);font-size:10.5px;font-weight:700;text-align:center;
  letter-spacing:.14em;margin:0 0 5px;}
.absblock p{font-size:11px;line-height:1.4;text-align:justify;text-indent:0;}
.absblock .kw{margin-top:7px;font-size:10.6px;}

/* ── 차례 ── */
.tocblock{font-size:11px;}
.toc-h{font-family:var(--sans);font-size:12.4px;font-weight:700;margin:0 0 5px;
  padding-bottom:2px;border-bottom:1px solid #111;}
.toc-list,.toc-figs{list-style:none;padding:0;margin:0 0 12px;}
.toc-list li{display:flex;gap:.5em;align-items:baseline;margin:1px 0;}
.toc-list li.h2{font-weight:700;margin-top:5px;}
.toc-list li.h3{padding-left:1.1em;}
.toc-figs li{display:flex;gap:.5em;align-items:baseline;margin:1px 0;font-size:10.4px;}
.tx{flex:1 1 auto;min-width:0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;}
.pn{flex:0 0 auto;font-variant-numeric:tabular-nums;color:#333;}
.go{cursor:pointer;}
.go:hover .tx{color:var(--link);}

/* ── 그림·표 (전폭) ── */
figure{margin:0 0 9px;text-align:center;}
figure img{max-width:100%;height:auto;}
figcaption{font-size:10.6px;line-height:1.34;text-align:justify;margin:4px 0 0;text-indent:0;}
.tw{margin:0 0 9px;}
table{border-collapse:collapse;width:100%;font-size:9.8px;line-height:1.28;
  font-variant-numeric:tabular-nums;}
thead th{border-top:1.1px solid #111;border-bottom:.7px solid #111;text-align:left;
  padding:2.6px 4px;font-weight:700;vertical-align:bottom;}
tbody td{padding:2.2px 4px;vertical-align:top;}
tbody tr:last-child td{border-bottom:1.1px solid #111;}

/* ── 읽기 모드 — A4 2단은 좁은 화면에서 못 읽는다(375px 폰에서 실효 5.5px).
   판면을 접고 한 단으로 흘려 읽는다. PDF 뷰어의 reflow 에 해당한다. ── */
body.reflow #stage{display:none;}
body.reflow #src{display:block;max-width:40rem;margin:0 auto;background:#fff;
  padding:1.4rem 1.15rem 5rem;font-size:16px;line-height:1.78;}
body.reflow #src .tocblock{display:none;}
body.reflow #src p{text-align:left;text-indent:0;margin:0 0 .95rem;word-break:keep-all;}
body.reflow #src h2{font-size:1.5rem;line-height:1.3;margin:1.8rem 0 .6rem;}
body.reflow #src h3{font-size:1.16rem;margin:1.4rem 0 .4rem;}
body.reflow #src h4{font-size:1.02rem;margin:1.1rem 0 .3rem;}
body.reflow #src figure{margin:1.4rem 0;}
body.reflow #src figure img{max-height:none;width:100%;}
body.reflow #src figcaption{font-size:.86rem;line-height:1.5;}
body.reflow #src table{font-size:.84rem;}
body.reflow #src .tw{overflow-x:auto;}
body.reflow #src .titleblock .ptitle{font-size:1.6rem;}
body.reflow #src .absblock{margin:0 0 1.2rem;}
body.reflow #src .jrnl{display:block;margin-bottom:1.4rem;}
body.reflow #src .jr{display:inline-block;margin-top:.45rem;font-size:13px;}
body.reflow #src .absblock p{font-size:.98rem;}
body.reflow #zout,body.reflow #zin,body.reflow #zfit,
body.reflow #zoomind,body.reflow #pgind,body.reflow #toctop,
body.reflow #prev,body.reflow #next{display:none;}

/* ── 메모 ── */
mark.apexnote{background:#fff2a8;color:inherit;padding:0 .05em;border-bottom:1.5px solid #d9b100;}
mark.apexnote.flash{animation:apexflash 1.2s ease-out;}
@keyframes apexflash{0%{background:#ffd23f}100%{background:#fff2a8}}
#notebtn{position:absolute;z-index:40;display:none;font-family:var(--sans);font-size:12px;
  background:#1f2329;color:#fff;border:0;border-radius:4px;padding:6px 10px;cursor:pointer;
  box-shadow:0 2px 8px rgba(0,0,0,.35);}
#notepanel{position:fixed;right:0;top:0;bottom:0;width:min(23rem,88vw);z-index:41;
  background:#fbfbfc;border-left:1px solid #cfd3da;display:none;flex-direction:column;
  font-family:var(--sans);box-shadow:-2px 0 14px rgba(0,0,0,.18);}
#notepanel.on{display:flex;}
#notepanel header{display:flex;align-items:center;gap:.5rem;padding:.7rem .8rem;
  border-bottom:1px solid #e2e5ea;font-size:13px;font-weight:700;}
#notepanel header .sp{flex:1;}
#notepanel header button{font:inherit;font-size:11.5px;font-weight:400;padding:.3rem .55rem;
  border:1px solid #c4c8ce;background:#fff;border-radius:3px;cursor:pointer;}
#notelist{flex:1;overflow-y:auto;padding:.5rem .6rem 2rem;}
.noteitem{border:1px solid #e2e5ea;border-radius:4px;padding:.55rem .6rem;margin:0 0 .5rem;
  background:#fff;font-size:12.5px;line-height:1.5;}
.noteitem .when{font-size:10.5px;color:#7b828c;}
.noteitem .quote{margin:.3rem 0;padding-left:.5rem;border-left:3px solid #ffd23f;color:#333;
  font-size:11.5px;max-height:4.4em;overflow:hidden;}
.noteitem .body{white-space:pre-wrap;color:#111;}
.noteitem .row{display:flex;gap:.35rem;margin-top:.45rem;}
.noteitem .row button{font:inherit;font-size:11px;padding:.22rem .5rem;border:1px solid #c4c8ce;
  background:#fff;border-radius:3px;cursor:pointer;}
.noteitem .row .del{color:#7a1010;}
#noteempty{color:#7b828c;font-size:12px;padding:1rem .3rem;line-height:1.6;}
#noteexport{position:fixed;inset:0;z-index:42;background:rgba(20,22,26,.55);display:none;
  align-items:center;justify-content:center;padding:1rem;}
#noteexport.on{display:flex;}
#noteexport .box{background:#fff;border-radius:6px;width:min(46rem,96vw);max-height:86vh;
  display:flex;flex-direction:column;padding:.9rem;}
#noteexport textarea{flex:1;min-height:44vh;font:12px/1.5 "Courier New",monospace;
  border:1px solid #cfd3da;border-radius:4px;padding:.6rem;resize:none;}
#noteexport .row{display:flex;gap:.4rem;justify-content:flex-end;margin-top:.6rem;}
#noteexport button{font:inherit;font-size:12px;padding:.4rem .8rem;border:1px solid #c4c8ce;
  background:#fff;border-radius:3px;cursor:pointer;}

#notecompose,#notepop{position:absolute;z-index:44;background:#fff;border:1px solid #c8ccd3;
  border-radius:6px;box-shadow:0 6px 22px rgba(0,0,0,.22);padding:.6rem;display:none;
  width:min(20rem,86vw);font-family:var(--sans);}
#notecompose.on,#notepop.on{display:block;}
#notecompose .q,#notepop .q{font-size:11px;color:#6b7078;border-left:3px solid #ffd23f;
  padding-left:.45rem;margin:0 0 .45rem;max-height:3.2em;overflow:hidden;line-height:1.45;}
#notecompose textarea{width:100%;box-sizing:border-box;min-height:5.2em;font:13px/1.55 inherit;
  border:1px solid #cfd3da;border-radius:4px;padding:.45rem;resize:vertical;}
#notecompose .row,#notepop .row{display:flex;gap:.35rem;justify-content:flex-end;margin-top:.45rem;}
#notecompose button,#notepop button{font:inherit;font-size:12px;padding:.32rem .7rem;
  border:1px solid #c4c8ce;background:#fff;border-radius:3px;cursor:pointer;}
#notecompose .save{background:#1f2329;color:#fff;border-color:#1f2329;}
#notepop .body{font-size:12.5px;line-height:1.55;white-space:pre-wrap;margin:0 0 .2rem;}
#notepop .when{font-size:10.5px;color:#7b828c;margin-bottom:.35rem;}
#notepop .del{color:#7a1010;}
mark.apexnote{cursor:pointer;}
#shortcuts{position:fixed;inset:0;z-index:45;display:none;align-items:center;justify-content:center;
  background:rgba(20,22,26,.55);padding:1rem;font-family:var(--sans);}
#shortcuts.on{display:flex;}
#shortcuts .box{width:min(32rem,94vw);background:#fff;border-radius:6px;box-shadow:0 8px 28px rgba(0,0,0,.28);
  padding:1rem 1.1rem;}
#shortcuts header{display:flex;align-items:center;gap:.6rem;margin-bottom:.65rem;}
#shortcuts h2{font-family:var(--sans);font-size:15px;margin:0;}
#shortcuts header .sp{flex:1;}
#shortcuts header button{font:inherit;font-size:11.5px;padding:.3rem .55rem;border:1px solid #c4c8ce;
  background:#fff;border-radius:3px;cursor:pointer;}
#shortcuts dl{display:grid;grid-template-columns:10rem 1fr;gap:.35rem .8rem;margin:0;font-size:12px;line-height:1.45;}
#shortcuts dt{font-weight:700;}
#shortcuts kbd{display:inline-block;min-width:1.4em;padding:.08rem .3rem;border:1px solid #c4c8ce;
  border-bottom-width:2px;border-radius:3px;background:#f6f7f8;font:11px/1.25 var(--sans);text-align:center;}
@media print{ #notebtn,#notepanel,#noteexport,#notecompose,#notepop{display:none!important;} mark.apexnote{background:none;} }

/* ── 뷰어 크롬 ── */
#hud{position:fixed;right:12px;bottom:12px;z-index:20;display:flex;gap:6px;
  align-items:center;background:rgba(28,30,34,.9);color:#eceef1;border-radius:4px;
  padding:6px 9px;font-family:var(--sans);font-size:12px;}
#hud button{font:inherit;background:#3a3e45;color:#eceef1;border:0;border-radius:3px;
  padding:4px 8px;cursor:pointer;}
#hud button:hover{background:#4b505a;}
#hud .nav{font-size:15px;line-height:1;padding:3px 8px;}
#pgind{font-variant-numeric:tabular-nums;min-width:4.6em;text-align:center;}
#zoomind{font-variant-numeric:tabular-nums;min-width:3.2em;text-align:center;font-size:11px;color:#c9ced6;}
#loading{position:fixed;inset:0;display:flex;align-items:center;justify-content:center;
  background:#8f9298;color:#f0f1f3;font-family:var(--sans);font-size:13px;z-index:30;}

@media print{
  @page{size:A4;margin:0;}
  html,body{background:#fff;}
  #hud,#loading,#shortcuts{display:none!important;}
  #stage{overflow:visible;height:auto!important;}
  #sizer{width:auto!important;height:auto!important;margin:0!important;}
  #book{transform:none!important;margin:0!important;}
  .page{box-shadow:none;margin:0;break-after:page;page-break-after:always;
    break-inside:avoid;page-break-inside:avoid;}
  .page:last-child{break-after:auto;page-break-after:auto;}
}
"""

JS = r"""
(function(){
  var PW=794;
  var src=document.getElementById('src'), book=document.getElementById('book'),
      stage=document.getElementById('stage'), load=document.getElementById('loading'),
      ind=document.getElementById('pgind'), sizer=document.getElementById('sizer'),
      zind=document.getElementById('zoomind');
  var pages=[], pending=[], scale=1;

  function el(t,c){ var e=document.createElement(t); if(c) e.className=c; return e; }
  function isHead(n){ return n && /^H[234]$/.test(n.tagName); }
  // Tables remain float-like because a wide table can consume the whole page.
  // Figures are page-level floats.  Keeping a figure inside a text column makes
  // the two columns unequal and leaves a large hole beside the image; the
  // source paragraph and figure may therefore land on adjacent pages.
  var INLINE_FIGURES={};
  function isInlineFigure(n){
    if (!n || n.tagName!=='FIGURE') return false;
    var im=n.querySelector('img'), m=im && /그림\s*(\d+)/.exec(im.alt||'');
    return !!(m && INLINE_FIGURES[parseInt(m[1],10)]);
  }
  function isSpanBlock(n){
    return n && (n.classList.contains('tw') ||
                 (n.tagName==='FIGURE' && !isInlineFigure(n)));
  }

  /* 판면 치수는 CSS 변수가 아니라 여기서 정하고 인라인으로 박는다.
     아티팩트처럼 호스트가 자기 스타일시트를 얹는 환경에서 :root 변수나 클래스
     규칙이 밀리면 지면 높이가 안 먹고 단이 무한정 늘어난다(2026-07-28 실제로 발생:
     로컬 23쪽 / 아티팩트 4쪽, 오른쪽 단이 빈 채 왼쪽 단만 넘침). 인라인 선언은
     호스트 스타일시트보다 우선하므로 조판이 환경에 상관없이 같아진다. */
  var GEO={PW:794, PH:1123, MX:64, MT:76, MB:60, GUT:24};
  function newPage(){
    var pg=el('div','page');
    pg.innerHTML='<div class="run"></div><div class="pinner">'+
      '<div class="span"></div><div class="cols"><div class="col"></div>'+
      '<div class="col"></div></div><div class="span bot"></div></div>'+
      '<div class="folio"></div>';
    book.appendChild(pg);
    var inner=pg.querySelector('.pinner'), cols=pg.querySelector('.cols'),
        cs=pg.querySelectorAll('.col'),
        colW=(GEO.PW-2*GEO.MX-GEO.GUT)/2, innerH=GEO.PH-GEO.MT-GEO.MB;
    pg.style.cssText='position:relative;width:'+GEO.PW+'px;height:'+GEO.PH+'px;'+
      'box-sizing:border-box;overflow:hidden;background:#fff;margin:0 0 14px;'+
      'box-shadow:0 2px 10px rgba(0,0,0,.35);flex:none;';
    inner.style.cssText='position:absolute;left:'+GEO.MX+'px;top:'+GEO.MT+'px;'+
      'width:'+(GEO.PW-2*GEO.MX)+'px;height:'+innerH+'px;box-sizing:border-box;'+
      'display:flex;flex-direction:column;';
    cols.style.cssText='flex:1 1 auto;min-height:0;display:flex;gap:'+GEO.GUT+'px;';
    for (var i=0;i<cs.length;i++){
      cs[i].style.cssText='width:'+colW+'px;flex:0 0 '+colW+'px;min-width:0;'+
        'overflow:hidden;box-sizing:border-box;';
    }
    var P={el:pg, span:pg.querySelector('.span'), spanB:pg.querySelector('.span.bot'),
           cols:[cs[0],cs[1]], ci:0, inner:inner};
    pages.push(P); return P;
  }
  function fits(c){ return c.scrollHeight<=c.clientHeight+1; }
  /* 단 안에서 내용이 실제로 어디까지 찼는지. overflow:hidden 이라 scrollHeight 로는
     못 잰다(내용이 짧아도 clientHeight 밑으로 안 내려간다). */
  function contentBottom(c){
    var k=c.lastElementChild;
    if (!k) return 0;
    return k.getBoundingClientRect().bottom - c.getBoundingClientRect().top;
  }

  /* 전폭 요소(그림·표)는 지면 머리와 발치의 두 자리에 앉힌다 (A&A 식 상·하 float).
     자리가 지면당 하나뿐이면 3절처럼 본문 9쪽이 그림 17장을 참조하는 구간에서
     줄이 5-7쪽씩 밀린다(2026-08-03 실측 — 사용자: "3절 fig들은 왜 없고").
     두 자리 합쳐 판면의 72% 까지 쓰고, 남는 가운데 띠에 본문이 흐른다. */
  var draining=false;
  function spanCap(P, sp){
    // Readability first: the larger caption and figure budgets keep panel
    // labels legible. This intentionally trades a few extra pages for size.
    if (draining) return P.inner.clientHeight*0.66;
    var other=(sp===P.span ? P.spanB : P.span).offsetHeight;
    var candidate=sp.firstElementChild || pending[0];
    var hasCols=P.cols.some(function(c){ return c.children.length>0; });
    /* A table alone may use most of a sheet; when text is also present, keep
       at least 36% of the page for the two columns so they do not collapse. */
    var tableOnly=candidate && candidate.classList.contains('tw') && !hasCols && !other;
    var total=P.inner.clientHeight*(tableOnly ? 0.82 : 0.64);
    var room=total-other;
    var base=P.inner.clientHeight*(tableOnly ? 0.82 : 0.52);
    return Math.min(base, room);
  }
  /* 반복 0.9 곱 축소는 쓰지 않는다. giveBack 으로 큐에 돌아갔다가 다시 앉을 때
     이미 줄어든 크기에서 또 줄어 그림이 우표만 해진다(2026-08-03: 그림 3 이 89px).
     캡션을 뺀 목표 높이로 한 번에 맞추고, 너무 작아질 자리면 아예 안 앉힌다. */
  function sizeTo(node, sp, cap){
    var im=node.querySelector('img'); if(!im) return true;
    im.style.maxHeight=''; im.style.width='';        // 자연 크기에서 다시 계산
    var imH=im.getBoundingClientRect().height||400;
    var chrome=sp.offsetHeight-imH;                  // 캡션·여백
    var target=cap-chrome;
    if (imH<=target) return true;
    if ((target<140 || target<imH*0.72) && !draining) return false;
    im.style.maxHeight=Math.max(90, Math.round(target))+'px';
    im.style.width='auto';
    return true;
  }
  function flushOne(P, sp){
    while (pending.length){
      if (P._floatKind) return;
      var cap=spanCap(P, sp);      // 들어올 때의 밀린 정도로 한 번만 정한다
      var f=pending[0], wasEmpty=!sp.firstChild;
      sp.appendChild(f);
      if (sp.offsetHeight>cap){
        if (wasEmpty && sizeTo(f,sp,cap)){
          pending.shift(); P._floatKind=f.classList.contains('tw')?'table':'figure';
        }
        else { sp.removeChild(f); break; }
      } else {
        pending.shift(); P._floatKind=f.classList.contains('tw')?'table':'figure';
      }
    }
  }
  function flushSpan(P){
    flushOne(P, P.span);
    if (pending.length) flushOne(P, P.spanB);
  }
  function fresh(){ var P=newPage(); flushSpan(P); return P; }
  function advance(P){ if (P.ci===0){ P.ci=1; return P; } return fresh(); }

  function splitList(list, P){
    var tag=list.tagName.toLowerCase(), box=el(tag,list.className);
    P.cols[P.ci].appendChild(box);
    var kids=[].slice.call(list.children);
    for (var i=0;i<kids.length;i++){
      box.appendChild(kids[i]);
      if (!fits(P.cols[P.ci])){
        box.removeChild(kids[i]);
        if (!box.children.length){ box.appendChild(kids[i]); continue; }
        P=advance(P); box=el(tag,list.className);
        P.cols[P.ci].appendChild(box); box.appendChild(kids[i]);
      }
    }
    return P;
  }

  /* ── 문단 쪼개기 ─────────────────────────────────────────────
     문단을 원자 블록으로 다루면 단 끝 자투리가 늘 버려지고, 큰 문단이 빈 단에도
     안 들어가면 giveBackSpan 이 발치 그림을 도로 뺏는다(2026-08-03: 8쪽 발치가
     그래서 비었고 그림 줄이 지면당 한 장으로 후퇴). 진짜 조판처럼 문단을 단
     경계에서 자른다. 인라인 마크업은 Range.extractContents 가 보존한다. */
  function truncCloneHeight(orig, c, k){
    var cl=orig.cloneNode(true);
    c.appendChild(cl);
    var cnt=0, w=document.createTreeWalker(cl, NodeFilter.SHOW_TEXT, null), t, cut=null;
    outer: while((t=w.nextNode())){
      var re=/\S+\s*/g, m;
      while((m=re.exec(t.nodeValue))){
        cnt++;
        if (cnt===k){ cut={nd:t, off:m.index+m[0].length}; break outer; }
      }
    }
    if (cut){
      var r=document.createRange();
      r.setStart(cut.nd, cut.off); r.setEnd(cl, cl.childNodes.length);
      r.deleteContents();
    }
    var h=cl.getBoundingClientRect().height;
    c.removeChild(cl);
    return h;
  }
  function splitPara(node, c, avail){
    if (avail < 40) return null;                         // 두 줄도 안 되는 자투리
    var nW=(node.textContent.match(/\S+/g)||[]).length;
    if (nW < 14) return null;
    var lo=4, hi=nW-4, best=0;
    while (lo<=hi){
      var mid=(lo+hi)>>1;
      if (truncCloneHeight(node, c, mid) <= avail){ best=mid; lo=mid+1; }
      else hi=mid-1;
    }
    if (best < 4) return null;
    var cnt=0, w=document.createTreeWalker(node, NodeFilter.SHOW_TEXT, null), t, cut=null;
    outer: while((t=w.nextNode())){
      var re=/\S+\s*/g, m;
      while((m=re.exec(t.nodeValue))){
        cnt++;
        if (cnt===best){ cut={nd:t, off:m.index+m[0].length}; break outer; }
      }
    }
    if (!cut) return null;
    var r=document.createRange();
    r.setStart(cut.nd, cut.off); r.setEnd(node, node.childNodes.length);
    var frag=r.extractContents();
    var b=document.createElement('p');
    b.className=((node.className||'')+' pcont').trim();
    b.appendChild(frag);
    node.classList.add('psplit');                        // 앞토막: 마지막 줄도 양끝맞춤
    return b;
  }

  /* 빈 단인데도 블록이 안 들어가면, 지면 머리(그림·표)가 자리를 너무 먹은 것이다.
     마지막 것부터 pending 앞으로 되돌려 단을 넓혀 본다. 앞에 넣으므로 차례는 그대로다.
     되돌릴 것이 없으면 false — 그때는 어쩔 수 없이 넘치게 둔다. */
  function giveBackSpan(node, P, c){
    // 발치부터 되돌린다. 발치의 그림 번호가 머리보다 크므로, 발치를 먼저
    // unshift 하고 머리를 그 앞에 넣어야 pending 이 오름차순으로 남는다.
    var sps=[P.spanB, P.span];
    for (var si=0; si<sps.length; si++){
      var sp=sps[si];
      while (sp.children.length){
        var back=sp.lastElementChild;
        sp.removeChild(back);
        var bim=back.querySelector('img');           // 다음 자리에서 새로 계산하도록
        if (bim){ bim.style.maxHeight=''; bim.style.width=''; }
        pending.unshift(back);
        c.appendChild(node);
        if (fits(c)) return true;
        c.removeChild(node);
      }
    }
    return false;
  }
  /* 블록을 지금 단에 앉힌다. 안 들어가면 지면 머리를 되돌려 단을 넓혀 본다.
     되돌릴 것도 없으면 넘치게 두는 수밖에 없다 — 이 함수를 거치지 않고 그냥
     appendChild 하면 검사 없이 넘친다(2026-07-30: 10쪽 오른쪽 단 +141px). */
  function accept(node, P){
    var c=P.cols[P.ci];
    c.appendChild(node);
    if (fits(c)) return P;
    c.removeChild(node);
    if (giveBackSpan(node,P,c)) return P;
    c.appendChild(node);
    return P;
  }
  function onlyHeads(c){
    for (var i=0;i<c.children.length;i++) if (!isHead(c.children[i])) return false;
    return true;
  }

  function place(node, P){
    // A wide table belongs after the preceding text, not after the next
    // subsection.  Put it at the foot of the current sheet when there is room;
    // otherwise flush it on a fresh sheet before continuing with the heading.
    if (isSpanBlock(node)){
      if (P._floatKind){ pending.push(node); return P; }
      var bothColsUsed=P.ci===1 || P.cols[1].children.length>0;
      // 짧은 절 뒤에 빈 공간이 남아도 두 단이 모두 찰 때까지 기다리면,
      // 다음 그림·표가 새 쪽으로 밀려 지면 하단이 크게 비게 된다. 제목만
      // 있는 단은 제외하고, 실제 본문이 한 단이라도 있으면 하단 float를
      // 시도한다. 두 단이 찬 경우의 기존 동작도 그대로 유지한다.
      var hasBody=P.cols.some(function(c){
        return [].some.call(c.children,function(ch){ return !isHead(ch); });
      });
      var canBottom=bothColsUsed || hasBody;
      if (node.classList.contains('tw') && !pending.length && !P._bottomUsed && canBottom){
        var tsp=P.spanB;
        tsp.appendChild(node);
        if (tsp.offsetHeight<=spanCap(P,tsp) && P.cols.every(fits)){
          P._bottomUsed=true; P._floatKind='table';
          return P;
        }
        tsp.removeChild(node);
      }
      if (node.tagName==='FIGURE' && !pending.length && !P._bottomUsed && canBottom){
        var fsp=P.spanB;
        fsp.appendChild(node);
        var fcap=spanCap(P,fsp);
        if ((fsp.offsetHeight<=fcap || sizeTo(node,fsp,fcap) && fsp.offsetHeight<=fcap) &&
            P.cols.every(fits)){
          P._bottomUsed=true; P._floatKind='figure';
          return P;
        }
        fsp.removeChild(node);
      }
      pending.push(node); return P;
    }
    // A bottom float occupies the visual end of this A4 sheet.  The next
    // source block must therefore begin on the following sheet, otherwise it
    // would be painted above the float by the two-column flex layout.
    if (P._bottomUsed) P=fresh();
    // Do not defer a wide table across a section heading: that reverses the
    // table and heading in the rendered reading order.
    if (pending.length && isHead(node)) P=fresh();
    if (isInlineFigure(node)){
      // Finish any wide table that is waiting.  Prefer the bottom of the
      // current sheet so a numeric paragraph and its evidence are adjacent;
      // if it does not fit, start a fresh sheet and place the figure at the top.
      if (pending.length) P=fresh();
      var c0=P.cols[P.ci];
      if (c0.children.length && !P.spanB.firstChild){
        var bsp=P.spanB;
        bsp.appendChild(node);
        var bcap=spanCap(P,bsp);
        if ((bsp.offsetHeight<=bcap || sizeTo(node,bsp,bcap) && bsp.offsetHeight<=bcap) &&
            P.cols.every(fits)){
          P._bottomUsed=true;
          return P;
        }
        bsp.removeChild(node);
      }
      if (c0.children.length || P.span.firstChild || P.spanB.firstChild) P=fresh();
      P.span.appendChild(node);
      sizeTo(node, P.span, spanCap(P, P.span));
      return P;
    }
    var c=P.cols[P.ci];
    c.appendChild(node);
    if (fits(c)){
      // 제목을 넣고 남는 자리가 두 줄도 안 되면 본문이 못 따라오므로 제목만 남는다.
      // 넘치기를 기다리지 말고 그 자리에서 다음 단으로 넘긴다.
      //
      // 남은 높이는 scrollHeight 로 재면 안 된다. overflow:hidden 인 요소의
      // scrollHeight 는 내용이 짧아도 clientHeight 밑으로 안 내려가서 차이가 늘 0 이
      // 되고, 그러면 모든 제목이 밀려난다(2026-07-29: 24쪽 -> 27쪽, 잘린 제목 1 -> 5).
      if (isHead(node) && c.children.length>1 && c.clientHeight-contentBottom(c) < 48){
        c.removeChild(node);
        // 바로 앞이 또 제목이면 같이 끌고 간다. 안 그러면 h2 뒤의 h3 만 넘어가며
        // h2 가 단 끝에 홀로 남는다(2026-08-03: 17쪽 "5. 논의").
        var chain=[];
        while (c.children.length && isHead(c.lastElementChild)){
          chain.unshift(c.lastElementChild); c.removeChild(c.lastElementChild);
        }
        var Pn=advance(P);
        chain.forEach(function(h){ Pn.cols[Pn.ci].appendChild(h); });
        Pn.cols[Pn.ci].appendChild(node); return Pn;
      }
      return P;
    }
    c.removeChild(node);
    // 문단이면 단 끝 자투리를 채우도록 쪼갠다. 나머지는 다음 자리로 흘린다.
    if (node.tagName==='P'){
      var availH=c.clientHeight - contentBottom(c) - 3;
      var rest=splitPara(node, c, availH);
      if (rest){
        c.appendChild(node);
        return place(rest, P);
      }
    }
    var listy=(node.tagName==='OL'||node.tagName==='UL')&&node.children.length>1;
    // 이 단에 있는 것이 제목뿐(또는 빈 단)이면 옮길 데가 없다. 지면 머리를 되돌려
    // 제목과 본문을 한자리에 앉힌다. 안 그러면 제목만 남고 본문이 다음 단으로 간다.
    if (onlyHeads(c) && giveBackSpan(node,P,c)) return P;
    if (!c.firstChild){                       // 빈 단인데도 안 들어감
      if (listy) return splitList(node,P);
      c.appendChild(node);
      return advance(P);
    }
    // 단 끝에 남은 제목은 **연속된 것을 전부** 다음 단으로 옮긴다. 하나만 옮기면
    // h2 바로 뒤에 h3 가 붙은 자리에서 h2 가 홀로 남는다
    // (2026-07-29 실제 발생: 4쪽 오른쪽 단 끝의 "2. 설계와 구현").
    var orphans=[];
    while (c.children.length>1 && isHead(c.lastElementChild)){
      orphans.unshift(c.lastElementChild);
      c.removeChild(c.lastElementChild);
    }
    var P2=advance(P);
    orphans.forEach(function(h){ P2.cols[P2.ci].appendChild(h); });
    if (listy) return splitList(node,P2);
    var c2=P2.cols[P2.ci];
    c2.appendChild(node);
    if (!fits(c2)){
      c2.removeChild(node);
      if (onlyHeads(c2) && giveBackSpan(node,P2,c2)) return P2;
      if (!c2.firstChild){ c2.appendChild(node); return advance(P2); }
      return accept(node, advance(P2));
    }
    return P2;
  }

  /* 차례는 한 지면에 본문과 같은 2단으로 짠다(2026-08-03 사용자 지시: 목차 1페이지).
     넘치면 오른 단으로, 두 단이 다 차야만 새 지면을 연다. */
  function newTocPage(){
    var P=newPage();
    P.el.classList.add('p-toc');
    return P;
  }
  function layoutToc(toc){
    var P=newTocPage(), ci=0, col=P.cols[0];
    function ok(){ return col.scrollHeight<=col.clientHeight+1; }
    function nextCol(){
      if (ci===0){ ci=1; col=P.cols[1]; }
      else { P=newTocPage(); ci=0; col=P.cols[0]; }
    }
    [].slice.call(toc.children).forEach(function(node){
      // 그림 목록은 제목과 항목이 서로 다른 단으로 갈라지면 목차의
      // 시각적 계층이 사라진다. 새 단(또는 새 목차 페이지)에서 함께
      // 시작하게 해 오른쪽 위에 마지막 항목만 고립되는 것을 막는다.
      var isFigureHeading=node.tagName==='DIV' &&
        node.classList.contains('toc-h') && node.nextElementSibling &&
        node.nextElementSibling.classList.contains('toc-figs');
      if (isFigureHeading && col.children.length) nextCol();
      var c=node.cloneNode(true);
      col.appendChild(c);
      if (ok()) return;
      col.removeChild(c);
      if (c.tagName==='OL'){
        var box=document.createElement('ol'); box.className=c.className;
        col.appendChild(box);
        [].slice.call(c.children).forEach(function(li){
          box.appendChild(li);
          if (!ok()){
            box.removeChild(li);
            if (!box.children.length) col.removeChild(box);
            nextCol();
            box=document.createElement('ol'); box.className=c.className;
            col.appendChild(box); box.appendChild(li);
          }
        });
        if (!box.children.length) col.removeChild(box);
      } else {
        nextCol(); col.appendChild(c);
      }
    });
  }

  var reflowInit=false;
  function build(){
    // 읽기 모드가 켜져 있으면 #stage 가 display:none 이라 단 높이가 0 으로 잡힌다.
    // 조판하는 동안만 판면을 보이게 한다.
    // (2026-07-30: 그림이 늦게 실려 재조판될 때 이 상태로 돌아 24쪽이 4쪽이 됐다)
    var wasReflow=document.body.classList.contains('reflow');
    if (wasReflow) document.body.classList.remove('reflow');
    // 진단 배너는 #book 밖에 붙으므로 여기서 직접 지운다. 안 지우면 그림이 실리기
    // 전 첫 조판에서 뜬 배너가 정상 재조판 뒤에도 남는다.
    var oldbn=document.getElementById('typoerr'); if (oldbn) oldbn.remove();
    book.innerHTML=''; pages=[]; pending=[]; draining=false;
    src.style.display='none';
    stage.style.cssText='overflow:hidden;';
    book.style.cssText='transform-origin:top left;width:'+GEO.PW+'px;'+
      'display:block;margin:0;padding:0;';
    // 체재(2026-08-03 사용자 지시): 1쪽 제목+안내 · 2쪽 차례 한 지면 ·
    // 3쪽부터 초록을 지면 머리에 얹고 바로 2단 본문 (A&A/AutoPhOT 식).
    // 1) 표제면 — 제목과 안내만.
    var T0=newPage(); T0.el.classList.add('p-title');
    var tb=src.querySelector('.titleblock');
    if (tb) T0.span.appendChild(tb.cloneNode(true));
    // 2) 차례 — 한 지면에 2단으로
    var toc=src.querySelector('.tocblock');
    if (toc) layoutToc(toc);
    // 3) 본문 — 초록을 첫 지면 머리(전폭)에 얹고 그 아래에서 2단이 시작된다
    var P=newPage();
    var abs=src.querySelector('.absblock');
    if (abs) P.span.appendChild(abs.cloneNode(true));
    var flow=src.querySelector('.flow');
    [].slice.call(flow.children).forEach(function(k){ P=place(k.cloneNode(true), P); });
    draining=true;                                              // 남은 그림 마무리
    while (pending.length){
      var n0=pending.length, Q=fresh();
      if (pending.length===n0){
        var f=pending.shift(); Q.span.appendChild(f); sizeTo(f,Q.span,spanCap(Q,Q.span));
      }
    }
    var N=pages.length;
    pages.forEach(function(p,i){
      if (i>0) p.el.querySelector('.run').textContent=
        'APEX — GUI 기반 구경·PSF 측광 파이프라인 (국문 원고)';
      var f=p.el.querySelector('.folio');
      f.className='folio '+((i+1)%2 ? 'r' : 'l');
      f.textContent='APEX, '+(i+1)+' / '+N+' 쪽';
    });
    /* 그림·표만 놓인 전용 쪽은 전폭으로 보여 주되, 본문이 한 단이라도
       있으면 3쪽 이후의 2단 구조를 유지한다. 그림이 있다는 이유만으로
       둘째 단을 숨기면 짧은 절이 전폭 한 단으로 늘어나 지면 흐름이 깨진다. */
    pages.forEach(function(p){
      if (p.el.querySelector('figure,.tw') && p.cols[0] && !p.cols[0].children.length &&
          p.cols[1] && !p.cols[1].children.length){
        p.el.classList.add('singlecol');
      }
    });
    linkToc();
    load.style.display='none';
    if (!reflowInit){ reflowInit=true; applyInitialReflow(); }  // 조판 뒤에 적용
    else if (wasReflow) setReflow(true);                        // 재조판이면 원래대로
    renderNotes(); markAll();
    // 조판 자체 점검. 넘치는 단이 있으면 지면 높이가 안 먹은 것이다.
    var over=0;
    pages.forEach(function(p){
      if(!p.cols) return;
      p.cols.forEach(function(c){
        if (c.offsetParent!==null && c.scrollHeight>c.clientHeight+2) over++;
      });
    });
    // 흘려보낸 블록이 다 앉았는지도 센다. 조판 함수가 다른 함수에 가려지면 쪽수만으로는
    // 안 잡힌다(2026-07-30: 메모 말풍선의 place() 가 조판 place() 를 덮어 본문 0블록).
    var want=flow.children.length, got=0;
    pages.forEach(function(p){
      if (p.cols) p.cols.forEach(function(c){ got+=c.children.length; });
      got+=p.span.children.length;
      if (p.spanB) got+=p.spanB.children.length;
    });
    if (over || pages.length<6 || got<want){
      var m=document.createElement('div'); m.id='typoerr';
      m.style.cssText='position:fixed;left:8px;bottom:8px;z-index:30;background:#7a1010;'+
        'color:#fff;font:12px/1.4 sans-serif;padding:6px 9px;border-radius:3px;';
      m.textContent='조판 오류: '+pages.length+'쪽 · 넘치는 단 '+over+'개 · 블록 '+got+'/'+want;
      document.body.appendChild(m);
      if (window.console) console.warn('[APEX] 조판 오류 '+m.textContent);
    }
    setPaged(true);
    var hp=/[#&?]p=(\d+)/.exec(location.hash||'');     // #p=7 로 그 쪽을 바로 연다
    if (hp){ var pi=parseInt(hp[1],10)-1;
      if (pages[pi]) goto(pi); }
  }

  function linkToc(){
    var where={};
    pages.forEach(function(p,i){
      [].slice.call(p.el.querySelectorAll('h2,h3')).forEach(function(h){
        var t=h.textContent.trim(); if(!(t in where)) where[t]=i; });
      [].slice.call(p.el.querySelectorAll('figcaption')).forEach(function(c){
        var m=/그림\s*(\d+)/.exec(c.textContent);
        if(m && !('f'+m[1] in where)) where['f'+m[1]]=i; });
    });
    function mark(li,i){
      if (i===undefined) return;
      li.querySelector('.pn').textContent=i+1;
      li.classList.add('go');
      li.addEventListener('click',function(){ goto(i); });
    }
    [].slice.call(book.querySelectorAll('.toc-list li')).forEach(function(li){
      mark(li, where[li.querySelector('.tx').textContent.trim()]); });
    [].slice.call(book.querySelectorAll('.toc-figs li')).forEach(function(li){
      var m=/그림\s*(\d+)/.exec(li.textContent); mark(li, m?where['f'+m[1]]:undefined); });
  }

  /* 확대율. 1 = 창 높이와 폭에 맞는 판면. 사용자가 확대·축소할 수 있다.
     transform 은 레이아웃에 영향을 주지 않으므로, 실제 크기를 갖는 #sizer 로
     스크롤 영역을 만들고 그 안에서 #book 을 시각적으로만 확대한다. */
  var zoom=1, currentPage=0;
  function applyPageVisibility(){
    if (!pages.length) return;
    pages.forEach(function(p,i){
      p.el.classList.toggle('active', i===currentPage);
    });
    if (ind) ind.textContent=(currentPage+1)+' / '+pages.length;
  }
  function setPaged(on){
    document.body.classList.toggle('paged', !!on);
    if (on){
      applyPageVisibility();
      if (stage) stage.scrollTo(0,0);
    } else {
      pages.forEach(function(p){ p.el.classList.remove('active'); });
    }
    fit();
  }
  function setZoom(z){
    /* 고배율에서도 판면 좌우를 충분히 탐색할 수 있게 상한을 600%로 둔다. */
    zoom=Math.max(0.5, Math.min(6, z));
    fit();
    try { localStorage.setItem('apexPaperZoom', String(zoom)); } catch(e){}
  }
  function fit(){
    if (document.body.classList.contains('reflow')) return;
    var w=stage.clientWidth||document.documentElement.clientWidth||PW;
    if (document.body.classList.contains('paged')){
      var h=window.innerHeight||GEO.PH;
      var pageScale=Math.min(1,(w-40)/PW,(h-38)/GEO.PH);
      scale=Math.max(0.5,pageScale*zoom);
      /* 판면 모드도 모든 쪽을 같은 스크롤 영역에 둔다. 페이지를 숨겼다
         다시 보이는 방식은 확대 페이지 끝에서 다음 쪽으로 이어지지 않는다. */
      book.style.height='auto';
      var bookH=Math.max(GEO.PH, book.scrollHeight);
      var pcw=Math.ceil(PW*scale), pch=Math.ceil(bookH*scale);
      sizer.style.width=pcw+'px';
      sizer.style.height=pch+'px';
      sizer.style.margin='18px auto 28px';
      book.style.transform='scale('+scale+')';
      book.style.width=PW+'px';
      book.style.height='auto';
      stage.style.height=Math.max(420,h)+'px';
      stage.style.overflowX=(pcw>w ? 'auto' : 'hidden');
      stage.style.overflowY='auto';
      if (zind) zind.textContent=Math.round(scale*100)+'%';
      applyPageVisibility();
      onScroll();
      return;
    }
    var base=Math.min(1,(w-20)/PW);
    scale=base*zoom;
    var cw=Math.ceil(PW*scale), ch=Math.ceil(book.scrollHeight*scale);
    sizer.style.width=cw+'px';
    sizer.style.height=ch+'px';
    sizer.style.margin='0 auto';
    book.style.transform='scale('+scale+')';
    stage.style.overflowX = cw>w ? 'auto' : 'hidden';
    stage.style.height=(ch+2)+'px';
    if (zind) zind.textContent=Math.round(scale*100)+'%';
    onScroll();
  }
  function goto(i){
    var p=pages[i]; if(!p) return;
    if (document.body.classList.contains('paged')){
      currentPage=i;
      fit();
      var top=(sizer.offsetTop||18)+(p.el.offsetTop*scale)-8;
      if (stage) stage.scrollTo({top:Math.max(0,top), behavior:'smooth'});
      if (ind) ind.textContent=(currentPage+1)+' / '+pages.length;
      return;
    }
    window.scrollTo({top:Math.max(0,p.el.offsetTop*scale-6), behavior:'smooth'});
  }
  function onScroll(){
    if (!pages.length) return;
    if (document.body.classList.contains('paged')){
      var base=sizer.offsetTop||0;
      var y=((stage.scrollTop||0)-base)/scale + 40, cur=0;
      for (var j=0;j<pages.length;j++){
        if (pages[j].el.offsetTop<=y) cur=j; else break;
      }
      currentPage=cur;
      if (ind) ind.textContent=(currentPage+1)+' / '+pages.length;
      return;
    }
    var y=(window.scrollY||window.pageYOffset||0)/scale + 40, cur=1;
    for (var i=0;i<pages.length;i++){ if (pages[i].el.offsetTop<=y) cur=i+1; else break; }
    ind.textContent=cur+' / '+pages.length;
  }
  window.addEventListener('scroll', onScroll, {passive:true});
  stage.addEventListener('scroll', onScroll, {passive:true});
  var t=null;
  window.addEventListener('resize', function(){ clearTimeout(t); t=setTimeout(fit,150); });
  function on(id, fn){ var e=document.getElementById(id); if (e) e.addEventListener('click', fn); }
  on('prev', function(){ goto(currentPage-1); });
  on('next', function(){ goto(currentPage+1); });
  var zoomWheelLock=false;
  document.addEventListener('wheel', function(ev){
    if (!document.body.classList.contains('paged') ||
        document.body.classList.contains('reflow')) return;
    /* PDF 뷰어 관례: Ctrl/Cmd+휠은 현재 쪽의 확대율을 바꾼다. 브라우저의
       기본 페이지 확대가 실행되지 않도록 먼저 이벤트를 소비한다. */
    if (ev.ctrlKey || ev.metaKey){
      ev.preventDefault();
      if (zoomWheelLock || Math.abs(ev.deltaY)<1) return;
      zoomWheelLock=true;
      setZoom(zoom*(ev.deltaY<0 ? 1.1 : 1/1.1));
      window.setTimeout(function(){ zoomWheelLock=false; },70);
      return;
    }
    /* Shift+휠은 세로 페이지 이동 대신 가로 패닝으로 쓴다. deltaX가
       없는 마우스 휠도 같은 동작을 하도록 deltaY를 가로값으로 대체한다. */
    if (ev.shiftKey){
      ev.preventDefault();
      var dx=ev.deltaX || ev.deltaY;
      if (Math.abs(dx)>0) stage.scrollLeft+=Math.max(-160,Math.min(160,dx*.65));
      return;
    }
    /* 일반 휠은 stage의 자연스러운 연속 스크롤에 맡긴다. */
  }, {passive:false});
  document.addEventListener('keydown', function(ev){
    var tag=(ev.target && ev.target.tagName)||'';
    if (tag==='INPUT' || tag==='TEXTAREA' || tag==='SELECT' || tag==='BUTTON' ||
        (ev.target && ev.target.isContentEditable)) return;
    var k=ev.key, mod=ev.ctrlKey||ev.metaKey;
    if (k==='?' || (mod && k==='/')) { ev.preventDefault(); toggleShortcuts(); return; }
    if (k==='Escape'){
      toggleShortcuts(false);
      if (typeof closeCompose==='function') closeCompose();
      if (typeof npop!=='undefined' && npop) npop.classList.remove('on');
      if (typeof togglePanel==='function') togglePanel(false);
      return;
    }
    if (mod && (k==='+' || k==='=' || k==='Add')) { ev.preventDefault(); setZoom(zoom*1.25); return; }
    if (mod && (k==='-' || k==='_' || k==='Subtract')) { ev.preventDefault(); setZoom(zoom/1.25); return; }
    if (mod && k==='0') { ev.preventDefault(); setZoom(1); return; }
    if (k==='+' || k==='=' || k==='Add') { ev.preventDefault(); setZoom(zoom*1.25); return; }
    if (k==='-' || k==='_' || k==='Subtract') { ev.preventDefault(); setZoom(zoom/1.25); return; }
    if (k==='0') { ev.preventDefault(); setZoom(1); return; }
    if (k==='m' || k==='M') { ev.preventDefault(); togglePanel(); return; }
    if (k==='r' || k==='R'){
      ev.preventDefault();
      setReflow(!document.body.classList.contains('reflow'));
      return;
    }
    if (document.body.classList.contains('reflow')) return;
    if (ev.shiftKey && (k==='ArrowLeft' || k==='ArrowRight')){
      ev.preventDefault();
      stage.scrollLeft=Math.max(0, Math.min(stage.scrollWidth-stage.clientWidth,
        stage.scrollLeft+(k==='ArrowRight' ? 80 : -80)));
      return;
    }
    if (k==='a' || k==='A' || k==='d' || k==='D'){
      ev.preventDefault();
      stage.scrollLeft=Math.max(0, Math.min(stage.scrollWidth-stage.clientWidth,
        stage.scrollLeft+((k==='d'||k==='D') ? 80 : -80)));
      return;
    }
    if (k==='Home') { ev.preventDefault(); goto(0); return; }
    if (k==='End') { ev.preventDefault(); goto(pages.length-1); return; }
    if (k==='ArrowLeft' || k==='PageUp' || k==='[') { ev.preventDefault(); goto(currentPage-1); return; }
    if (k==='ArrowRight' || k==='PageDown' || k===' ' || k===']') { ev.preventDefault(); goto(currentPage+1); return; }
  });
  on('zin',  function(){ setZoom(zoom*1.25); });
  on('zout', function(){ setZoom(zoom/1.25); });
  on('zfit', function(){ setZoom(1); });
  var shortcutPanel=document.getElementById('shortcuts');
  function toggleShortcuts(on){
    if (!shortcutPanel) return;
    var show=on===undefined ? !shortcutPanel.classList.contains('on') : !!on;
    shortcutPanel.classList.toggle('on', show);
    shortcutPanel.setAttribute('aria-hidden', show ? 'false' : 'true');
  }
  on('helpbtn', function(){ toggleShortcuts(); });
  on('shortcut-close', function(){ toggleShortcuts(false); });
  if (shortcutPanel) shortcutPanel.addEventListener('click', function(e){
    if (e.target===shortcutPanel) toggleShortcuts(false);
  });
  var rb=document.getElementById('rflow');
  function setReflow(on){
    document.body.classList.toggle('reflow', on);
    src.style.display = on ? 'block' : 'none';   // build() 가 건 인라인을 덮는다
    if (rb) rb.textContent = on ? '판면' : '읽기';
    try { localStorage.setItem('apexPaperReflow', on ? '1' : '0'); } catch(e){}
    if (!on) fit();
    markAll();
  }
  if (rb) rb.addEventListener('click', function(){
    setReflow(!document.body.classList.contains('reflow'));
  });
  // 좁은 화면은 읽기 모드로 시작한다. 사용자가 고른 값이 있으면 그것을 따른다.
  // 적용은 조판이 끝난 뒤에 한다(applyInitialReflow).
  function applyInitialReflow(){
    var saved = null;
    try { saved = localStorage.getItem('apexPaperReflow'); } catch(e){}
    setReflow(saved === null ? (window.innerWidth < 640) : saved === '1');
  }
  try {
    var z=parseFloat(localStorage.getItem('apexPaperZoom'));
    if (z && z>0.4 && z<6.2) zoom=z;
  } catch(e){}
  document.getElementById('top').addEventListener('click', function(){ goto(0); });
  document.getElementById('toctop').addEventListener('click', function(){ goto(1); });


  /* ══ 메모 ══════════════════════════════════════════════════════════
     인용문 문자열로 위치를 잡는다. DOM 경로로 잡으면 조판이 다시 될 때마다
     끊어진다. 원고가 바뀌어 인용문을 못 찾으면 표시만 못 하고 메모는 남는다. */
  var NKEY='apexPaperNotes';
  function nLoad(){ try{ return JSON.parse(localStorage.getItem(NKEY)||'[]'); }catch(e){ return []; } }
  function nSave(v){ try{ localStorage.setItem(NKEY, JSON.stringify(v)); }catch(e){} }
  var notes=nLoad();

  var nbtn=document.getElementById('notebtn'),
      npanel=document.getElementById('notepanel'),
      nlist=document.getElementById('notelist'),
      nempty=document.getElementById('noteempty'),
      nexp=document.getElementById('noteexport'),
      nta=document.getElementById('noteta'),
      ncount=document.getElementById('notecount');
  var pendingQuote='';

  function fmtTime(ts){
    var d=new Date(ts);
    function p(n){ return (n<10?'0':'')+n; }
    return d.getFullYear()+'-'+p(d.getMonth()+1)+'-'+p(d.getDate())+' '+p(d.getHours())+':'+p(d.getMinutes());
  }
  function activeRoot(){
    return document.body.classList.contains('reflow') ? src : book;
  }

  /* 선택이 생기면 그 자리에 단추를 띄운다 */
  document.addEventListener('mouseup', function(){
    setTimeout(function(){
      var sel=window.getSelection();
      if (!sel || sel.isCollapsed){ nbtn.style.display='none'; return; }
      var t=sel.toString().trim();
      if (t.length<2){ nbtn.style.display='none'; return; }
      var root=activeRoot();
      if (!root.contains(sel.anchorNode)){ nbtn.style.display='none'; return; }
      var r=sel.getRangeAt(0).getBoundingClientRect();
      pendingQuote=t.replace(/\s+/g,' ').slice(0,300);
      nbtn.style.display='block';
      nbtn.style.left=Math.max(8, r.left+window.scrollX)+'px';
      nbtn.style.top=(r.bottom+window.scrollY+6)+'px';
    }, 10);
  });
  document.addEventListener('scroll', function(){ nbtn.style.display='none'; }, {passive:true});

  /* 샌드박스 iframe 에서는 prompt/alert/confirm 이 막힌다. 전부 인라인 UI 로 한다. */
  var ncomp=document.getElementById('notecompose'),
      ncta=ncomp.querySelector('textarea'),
      ncq=ncomp.querySelector('.q'),
      npop=document.getElementById('notepop');
  var editingId=null;

  /* 조판기의 place(node,P) 와 이름이 겹치면 안 된다. 함수 선언은 나중 것이 이겨서
     본문 조판이 통째로 이 함수로 넘어간다(2026-07-30: 24쪽 -> 4쪽, 본문 0블록). */
  function placePop(el, x, y){
    el.classList.add('on');
    var w=el.offsetWidth, vw=document.documentElement.clientWidth;
    el.style.left=Math.max(8, Math.min(x, vw-w-8))+'px';
    el.style.top=(y+8)+'px';
  }
  function openCompose(quote, x, y, existing){
    npop.classList.remove('on');
    editingId = existing ? existing.id : null;
    ncq.textContent = quote;
    ncta.value = existing ? (existing.body||'') : '';
    placePop(ncomp, x, y);
    ncta.focus();
  }
  function closeCompose(){ ncomp.classList.remove('on'); editingId=null; }

  nbtn.addEventListener('click', function(){
    var r=nbtn.getBoundingClientRect();
    nbtn.style.display='none';
    openCompose(pendingQuote, r.left+window.scrollX, r.bottom+window.scrollY);
  });
  ncomp.querySelector('.cancel').addEventListener('click', closeCompose);
  ncomp.querySelector('.save').addEventListener('click', saveCompose);
  ncta.addEventListener('keydown', function(e){
    if (e.key==='Enter' && !e.shiftKey){ e.preventDefault(); saveCompose(); }
    if (e.key==='Escape'){ closeCompose(); }
  });
  function saveCompose(){
    var body=ncta.value.trim();
    if (editingId){
      notes.forEach(function(n){ if (n.id===editingId) n.body=body; });
    } else {
      if (!body && !ncq.textContent) { closeCompose(); return; }
      notes.push({ id:Date.now()+'-'+Math.random().toString(36).slice(2,7),
                   quote:ncq.textContent, body:body, at:Date.now() });
    }
    /* 저장 직후에는 본문과 작성창만 갱신한다. 메모 목록은 사용자가
       메모 버튼이나 M 단축키로 열었을 때만 표시한다. */
    nSave(notes); closeCompose(); renderNotes(); markAll();
  }

  /* 형광펜을 누르면 말풍선 */
  document.addEventListener('click', function(e){
    var m=e.target.closest && e.target.closest('mark.apexnote');
    if (!m){
      if (!e.target.closest || !e.target.closest('#notepop')) npop.classList.remove('on');
      if (e.target.closest && !e.target.closest('#notecompose') && e.target!==nbtn) {
        /* 작성 중이면 유지 */
      }
      return;
    }
    var n=notes.filter(function(x){ return x.id===m.dataset.nid; })[0];
    if (!n) return;
    npop.querySelector('.when').textContent=fmtTime(n.at);
    npop.querySelector('.q').textContent=n.quote;
    npop.querySelector('.body').textContent=n.body||'(메모 없음)';
    npop.dataset.nid=n.id;
    var r=m.getBoundingClientRect();
    placePop(npop, r.left+window.scrollX, r.bottom+window.scrollY);
  });
  npop.querySelector('.close').addEventListener('click', function(){ npop.classList.remove('on'); });
  npop.querySelector('.edit').addEventListener('click', function(){
    var n=notes.filter(function(x){ return x.id===npop.dataset.nid; })[0];
    if (!n) return;
    var r=npop.getBoundingClientRect();
    openCompose(n.quote, r.left+window.scrollX, r.top+window.scrollY, n);
  });
  npop.querySelector('.del').addEventListener('click', function(){
    notes=notes.filter(function(x){ return x.id!==npop.dataset.nid; });
    nSave(notes); npop.classList.remove('on'); renderNotes(); markAll();
  });

  /* 인용문을 본문에서 찾아 표시한다 */
  function clearMarks(){
    [].slice.call(document.querySelectorAll('mark.apexnote')).forEach(function(m){
      var p=m.parentNode; while(m.firstChild) p.insertBefore(m.firstChild,m);
      p.removeChild(m); p.normalize();
    });
  }
  /* 인용문을 블록(문단·제목·항목·캡션) 단위로 정규화해 찾는다.
     - 예전에는 한 텍스트 노드 안에서만 찾아서, 기울임·코드 조각을 낀 선택과
       제목 선택이 표시되지 않았다. 제목은 같은 문자열이 차례에 먼저 나와
       차례 쪽에 표시가 붙기도 했다(2026-08-03) — 차례·러닝헤드는 건너뛴다.
     - 블록 안에서 텍스트 노드들을 이어붙여 공백 정규화 문자열과
       (노드, 오프셋) 매핑을 만들고, 걸치는 노드마다 mark 를 하나씩 감싼다. */
  function mkMark(n){
    var m=document.createElement('mark');
    m.className='apexnote'; m.dataset.nid=n.id; m.title=n.body||'(메모 없음)';
    return m;
  }
  function markAll(){
    clearMarks();
    var root=activeRoot();
    var blocks=[].slice.call(root.querySelectorAll('p,h2,h3,h4,li,figcaption,td,caption'));
    notes.forEach(function(n){
      var q=(n.quote||'').replace(/\s+/g,' ').trim();
      if (q.length<2) return;
      for (var bi=0;bi<blocks.length;bi++){
        var b=blocks[bi];
        if (b.closest('.p-toc') || b.closest('.tocblock') ||
            b.closest('.run') || b.closest('.folio')) continue;
        var w=document.createTreeWalker(b, NodeFilter.SHOW_TEXT, null), t, nodes=[];
        while ((t=w.nextNode())) nodes.push(t);
        if (!nodes.length) continue;
        var norm='', map=[], lastSp=true;
        for (var k=0;k<nodes.length;k++){
          var sv=nodes[k].nodeValue;
          for (var i=0;i<sv.length;i++){
            if (/\s/.test(sv[i])){
              if (lastSp) continue;
              norm+=' '; map.push({nd:nodes[k], off:i}); lastSp=true;
            } else { norm+=sv[i]; map.push({nd:nodes[k], off:i}); lastSp=false; }
          }
        }
        var j=norm.indexOf(q);
        if (j<0) continue;
        var s0=map[j], e0=map[j+q.length-1];
        try {
          var rg=document.createRange();
          rg.setStart(s0.nd, s0.off); rg.setEnd(e0.nd, e0.off+1);
          if (rg.startContainer===rg.endContainer){
            rg.surroundContents(mkMark(n));
          } else {
            var w2=document.createTreeWalker(b, NodeFilter.SHOW_TEXT, null), nd, seg=[];
            while ((nd=w2.nextNode())) if (rg.intersectsNode(nd)) seg.push(nd);
            seg.forEach(function(tn){
              var r2=document.createRange();
              r2.setStart(tn, tn===rg.startContainer ? rg.startOffset : 0);
              r2.setEnd(tn, tn===rg.endContainer ? rg.endOffset : tn.nodeValue.length);
              if (r2.toString().trim().length){
                try { r2.surroundContents(mkMark(n)); } catch(e){}
              }
            });
          }
        } catch(e){}
        break;
      }
    });
  }

  function renderNotes(){
    ncount.textContent=notes.length ? '메모 '+notes.length : '메모';
    nlist.innerHTML='';
    nempty.style.display = notes.length ? 'none' : 'block';
    notes.slice().sort(function(a,b){ return b.at-a.at; }).forEach(function(n){
      var d=document.createElement('div'); d.className='noteitem';
      var w=document.createElement('div'); w.className='when'; w.textContent=fmtTime(n.at);
      var q=document.createElement('div'); q.className='quote'; q.textContent=n.quote;
      var b=document.createElement('div'); b.className='body'; b.textContent=n.body||'(메모 없음)';
      var row=document.createElement('div'); row.className='row';
      var go=document.createElement('button'); go.textContent='이동';
      go.addEventListener('click', function(){ jumpTo(n); });
      var ed=document.createElement('button'); ed.textContent='수정';
      ed.addEventListener('click', function(){
        var r=ed.getBoundingClientRect();
        openCompose(n.quote, Math.max(8, r.left+window.scrollX-260), r.bottom+window.scrollY, n);
      });
      var del=document.createElement('button'); del.className='del'; del.textContent='삭제';
      del.addEventListener('click', function(){
        notes=notes.filter(function(x){ return x.id!==n.id; });
        nSave(notes); renderNotes(); markAll();
      });
      row.appendChild(go); row.appendChild(ed); row.appendChild(del);
      d.appendChild(w); d.appendChild(q); d.appendChild(b); d.appendChild(row);
      nlist.appendChild(d);
    });
  }
  function jumpTo(n){
    var m=document.querySelector('mark.apexnote[data-nid="'+n.id+'"]');
    if (!m){ nempty.style.display='block';
      nempty.textContent='본문에서 이 부분을 찾지 못했습니다. 원고가 수정되었을 수 있습니다.';
      return; }
    m.scrollIntoView({block:'center', behavior:'smooth'});
    m.classList.remove('flash'); void m.offsetWidth; m.classList.add('flash');
  }
  function togglePanel(on){
    npanel.classList.toggle('on', on===undefined ? !npanel.classList.contains('on') : on);
  }

  function exportText(){
    if (!notes.length) return '메모가 없습니다.';
    var out=['# 원고 메모 ('+notes.length+'건, '+fmtTime(Date.now())+' 내보냄)',''];
    notes.slice().sort(function(a,b){ return a.at-b.at; }).forEach(function(n,i){
      out.push('## '+(i+1)+'. '+fmtTime(n.at));
      out.push('> '+n.quote);
      out.push('');
      out.push(n.body||'(메모 없음)');
      out.push('');
    });
    return out.join('\n');
  }
  document.getElementById('nexport').addEventListener('click', function(){
    nta.value=exportText(); nexp.classList.add('on'); nta.focus(); nta.select();
  });
  document.getElementById('ncopy').addEventListener('click', function(){
    nta.select();
    try { document.execCommand('copy'); } catch(e){}
    if (navigator.clipboard) { navigator.clipboard.writeText(nta.value).catch(function(){}); }
    document.getElementById('ncopy').textContent='복사됨';
    setTimeout(function(){ document.getElementById('ncopy').textContent='복사'; }, 1400);
  });
  document.getElementById('nclose').addEventListener('click', function(){ nexp.classList.remove('on'); });
  var clearArm=false, clearBtn=document.getElementById('nclear');
  clearBtn.addEventListener('click', function(){
    if (!notes.length) return;
    if (!clearArm){                      // confirm 이 막히므로 두 번 눌러 확인한다
      clearArm=true; clearBtn.textContent='한 번 더';
      setTimeout(function(){ clearArm=false; clearBtn.textContent='전체 삭제'; }, 3000);
      return;
    }
    clearArm=false; clearBtn.textContent='전체 삭제';
    notes=[]; nSave(notes); renderNotes(); markAll();
  });
  document.getElementById('npanelclose').addEventListener('click', function(){ togglePanel(false); });
  document.getElementById('nopen').addEventListener('click', function(){ togglePanel(); });

  /* 그림이 아직 안 실린 상태에서 조판하면 높이를 0으로 재서 쪽수가 틀어진다.
     한 번 짜고, 남은 그림이 다 실리면 한 번만 다시 짠다. */
  function boot(){
    build();
    var late=[].slice.call(src.querySelectorAll('img')).filter(function(i){ return !i.complete; });
    if (!late.length) return;
    var left=late.length, done=function(){ if (--left<=0) build(); };
    late.forEach(function(i){
      i.addEventListener('load', done, {once:true});
      i.addEventListener('error', done, {once:true});
    });
  }
  if (document.readyState==='complete') boot();
  else window.addEventListener('load', boot);
})();
"""

HTML = f"""<meta charset="utf-8">
<title>APEX — GUI 기반 구경·PSF 측광 파이프라인 (국문 원고)</title>
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>{CSS}</style>
<div id="src">
{BODY}
</div>
<div id="stage"><div id="sizer"><div id="book"></div></div></div>
<button id="notebtn" type="button">메모</button>
<div id="notecompose">
  <p class="q"></p>
  <textarea placeholder="메모를 적으세요 (Enter 저장 · Shift+Enter 줄바꿈)"></textarea>
  <div class="row"><button type="button" class="cancel">취소</button>
    <button type="button" class="save">저장</button></div>
</div>
<div id="notepop">
  <div class="when"></div>
  <p class="q"></p>
  <p class="body"></p>
  <div class="row"><button type="button" class="edit">수정</button>
    <button type="button" class="del">삭제</button>
    <button type="button" class="close">닫기</button></div>
</div>

<aside id="notepanel">
  <header><span id="notecount">메모</span><span class="sp"></span>
    <button id="nexport" type="button">내보내기</button>
    <button id="nclear" type="button">전체 삭제</button>
    <button id="npanelclose" type="button">닫기</button>
  </header>
  <div id="notelist"></div>
  <div id="noteempty">본문에서 문장을 드래그하면 「메모」 단추가 뜹니다.<br><br>
    메모는 이 브라우저에 저장되며, 원고를 다시 올려도 남습니다.<br>
    다 적은 뒤 <b>내보내기</b>로 복사해 대화에 붙여 주세요.</div>
</aside>
<div id="noteexport"><div class="box">
  <textarea id="noteta" readonly></textarea>
  <div class="row">
    <button id="ncopy" type="button">복사</button>
    <button id="nclose" type="button">닫기</button>
  </div>
</div></div>
<div id="shortcuts" aria-hidden="true">
  <div class="box" role="dialog" aria-modal="true" aria-labelledby="shortcut-title">
    <header><h2 id="shortcut-title">단축키</h2><span class="sp"></span>
      <button id="shortcut-close" type="button">닫기</button></header>
    <dl>
      <dt><kbd>+</kbd> / <kbd>=</kbd></dt><dd>확대</dd>
      <dt><kbd>-</kbd></dt><dd>축소</dd>
      <dt><kbd>Ctrl</kbd> + 휠</dt><dd>휠 위/아래로 확대·축소</dd>
      <dt><kbd>Shift</kbd> + 휠</dt><dd>좌우 패닝</dd>
      <dt><kbd>Shift</kbd>+<kbd>←</kbd>/<kbd>→</kbd> 또는 <kbd>A</kbd>/<kbd>D</kbd></dt><dd>좌우 80px 이동</dd>
      <dt><kbd>0</kbd> 또는 <kbd>Ctrl</kbd>+<kbd>0</kbd></dt><dd>맞춤 배율</dd>
      <dt><kbd>←</kbd> <kbd>→</kbd></dt><dd>이전/다음 페이지</dd>
      <dt><kbd>Home</kbd> / <kbd>End</kbd></dt><dd>첫 페이지/마지막 페이지</dd>
      <dt><kbd>M</kbd></dt><dd>메모 목록 열기/닫기</dd>
      <dt><kbd>R</kbd></dt><dd>읽기 모드 전환</dd>
      <dt><kbd>Esc</kbd></dt><dd>현재 팝업 닫기</dd>
      <dt><kbd>?</kbd></dt><dd>이 단축키 안내</dd>
    </dl>
  </div>
</div>
<div id="loading">판면을 조판하는 중…</div>
<div id="hud">
  <button id="top" type="button">처음</button>
  <button id="toctop" type="button">차례</button>
  <button id="prev" class="nav" type="button" title="이전 페이지">‹</button>
  <span id="pgind">– / –</span>
  <button id="next" class="nav" type="button" title="다음 페이지">›</button>
  <button id="zout" type="button" title="축소">&minus;</button>
  <span id="zoomind">100%</span>
  <button id="zin" type="button" title="확대">+</button>
  <button id="zfit" type="button" title="폭 맞춤">맞춤</button>
  <button id="rflow" type="button" title="좁은 화면에서 읽기">읽기</button>
  <button id="nopen" type="button" title="메모">메모</button>
  <button id="helpbtn" type="button" title="단축키 (?)">?</button>
</div>
<script>{JS}</script>"""

# 단독 파일은 doctype 이 있어야 한다. 없으면 quirks 모드로 떨어져 판면 계산이 틀린다.
# 아티팩트용 사본은 게시 때 skeleton 이 씌워지므로 뺀다.
OUT.write_text("<!doctype html>\n" + HTML, encoding="utf-8")
OUT_ARTIFACT.write_text(HTML, encoding="utf-8")
print("wrote", OUT, round(len(HTML) / 1e6, 2), "MB | citations", len(LAB), "| figs", len(FIGMAP))
print("wrote", OUT_ARTIFACT.name, "(아티팩트용 — doctype 없음)")
