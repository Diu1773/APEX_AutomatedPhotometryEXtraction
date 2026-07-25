# Figure — Operational depth QC: every frame's depth monitored, and predicted before detection runs

**Figure X.** Per-frame detection depth as an operational quality-control quantity,
for all 66 single exposures of the three clusters reduced raw→science by APEX
(cf. the survey convention of per-image depth monitoring; Kessler et al. 2015 Fig. 7).
**(a)** Distribution of the *realized* 50% detection depth — **each entry is one
single exposure**, not a stack or a per-target average — coloured by target
(M67, NGC 6811, M13). Realized depth is measured from real master-catalogue stars,
not from injections: each star's "true" magnitude is the median of its count-rate
instrumental magnitudes over all frames of the same filter, which is exposure-time
invariant and removes the single-frame Eddington bias that would otherwise inflate
the depth of the frame being tested; the frame's `detected_flag` roll-off through
50% is then read off directly. The spread — 14.3 to 17.6 mag, 3.3 mag — is the
observing conditions, not the method: it tracks sky brightness and seeing across
nights and filters, and is exactly the quantity a survey needs logged per frame.
**(b)** The same frames with the *predicted* depth on the abscissa. The prediction
uses **only the frame's background noise and PSF peak fraction** — no injections, no
knowledge of what was actually detected — through the injection-calibrated peak-S/N
law of panel (b) of the completeness figure,
$m_{50} = \mathrm{ZP} - 2.5\log_{10}\!\left(\mathrm{S/N_{50}}\,\sigma_e / p_{\rm peak}\right)$
with $\mathrm{S/N_{50}} = 4.05$, evaluated by the same function the pipeline's Step-7
QC gate calls (`apex.analysis.detection_limit.predict_frame_m50`), so the figure and
the shipped gate cannot drift apart. Over the 57 frames that pass the circularity
guard the agreement is **RMS 78 mmag** about the 1:1 line, worst frame 206 mmag, with
a mean offset of **+38 mmag** (realized fainter than predicted). That offset is a
known, stated bias and not noise: in crowded fields a real catalogue star can be
matched to a neighbouring detection, so the realized roll-off sits slightly faint.
The grey band is the $\pm0.5$ mag tolerance of the shipped QC gate
(`depth_qc_tolerance_mag`); every valid frame lies well inside it, which is what
makes the gate usable — a frame that *does* leave the band is anomalous for a reason
the noise and PSF cannot explain (focus drift, cloud, tracking, a calibration
defect), not merely a poor-seeing frame.

**Circularity guard (open symbols).** The realized depth of the *deepest* frame of a
filter is not measurable this way: the faint tail of that filter's master catalogue
exists only because those frames detected it, so the frame is scored against stars it
selected itself and its depth inflates. Frames whose depth is within 0.7 mag of the
filter master's own 90th-percentile limit are therefore flagged `depth_valid = False`
and excluded from the RMS (9 of 66: all five M13 $R$, three M13 $V$, one NGC 6811 $V$
— M13 is the most vulnerable because it was observed only at 60 s, so its master is
shallow). The guard is not cosmetic: the excluded M13 $R$ frame at realized 17.26 vs
predicted 15.59 sits visibly off the 1:1 line in panel (b), a $+1.45$ mag circular
inflation that would otherwise have been read as a prediction failure. Independently,
on the six injection-calibrated frames that pass the guard (heavy black outlines),
realized $-$ injected $= +0.11$ mag mean, 0.15 mag RMS — consistent with the same
blend-matching bias.

**Data.** 66 single exposures (30–480 s) of M13, M67 and NGC 6811 in six filters
($B$, $V$, $R$, $g$, $r$, $i$), one camera (Moravian C3-61000, gain 0.689 $e^-$/ADU
from photon-transfer measurement), reduced raw→science entirely by APEX; 57 pass the
circularity guard (M67 30, NGC 6811 20, M13 7), FWHM 5.2–9.3 px, 1139–2001 catalogue
stars per frame. Prediction reads no pixel data — background noise comes from the
Step-7 photometry table and the PSF peak fraction from the frame's own bright clean
stars. Per-frame values: `validation/paper/data_qc_depth/qc_depth_summary.csv`.

---

## 국문 (MANUSCRIPT_ko용)

**그림 X.** 운영 품질관리(QC) 양으로서의 프레임별 검출 깊이. APEX가 raw→science
전 과정을 처리한 세 성단의 단일 노출 66장 전체를 대상으로 한다(프레임별 깊이 상시
기록은 서베이의 표준 관행 — Kessler et al. 2015 Fig. 7 참조).
**(a)** *실현(realized)* 50% 검출 깊이의 분포로, **각 항목은 스택이나 대상별 평균이
아닌 단일 노출 1장**이며 대상별로 색을 구분했다(M67, NGC 6811, M13). 실현 깊이는
주입이 아니라 **실제 마스터 카탈로그 별**에서 측정한다: 각 별의 "참" 등급은 동일 필터
전 프레임에서의 count-rate 기기등급 중앙값으로, 노출시간에 불변이며 검사 대상 프레임의
깊이를 부풀릴 단일 프레임 Eddington 편향을 제거한다. 이어서 그 프레임의
`detected_flag`가 50%를 지나는 지점을 직접 읽는다. 14.3–17.6등급(3.3등급)의 산포는
방법이 아니라 관측 조건이다 — 밤과 필터에 따른 하늘 밝기·시상을 그대로 따라가며,
서베이가 프레임마다 기록해야 하는 바로 그 양이다.
**(b)** 같은 프레임들을 가로축 *예측* 깊이에 대해 나타낸 것. 예측은 **그 프레임의 배경
잡음과 PSF peak fraction만** 사용하며, 주입도 실제 검출 정보도 쓰지 않는다. 완전도
그림 (b)의 주입 보정된 peak-S/N 법칙
$m_{50} = \mathrm{ZP} - 2.5\log_{10}\!\left(\mathrm{S/N_{50}}\,\sigma_e / p_{\rm peak}\right)$,
$\mathrm{S/N_{50}} = 4.05$ 를 쓰되, 파이프라인 Step 7의 QC 게이트가 호출하는 것과
**동일한 함수**(`apex.analysis.detection_limit.predict_frame_m50`)로 계산해 그림과
배포 게이트가 어긋날 수 없게 했다. 순환 guard를 통과한 57프레임에서 1:1선 대비
**RMS 78 mmag**, 최악 프레임 206 mmag, 평균 편차 **+38 mmag**(실현이 예측보다 어두움)
이다. 이 편차는 잡음이 아니라 **알려진 편향이며 명시한다**: 혼잡 영역에서 실제 카탈로그
별이 이웃 검출에 매칭될 수 있어 실현 롤오프가 약간 어두운 쪽에 놓인다. 회색 띠는 배포된
QC 게이트의 허용 범위 $\pm0.5$등급(`depth_qc_tolerance_mag`)으로, 유효 프레임이 모두
그 안에 넉넉히 들어온다는 점이 게이트를 쓸 만하게 만든다 — 이 띠를 벗어나는 프레임은
단지 시상이 나쁜 프레임이 아니라, 잡음과 PSF로 설명되지 않는 이유(초점 이탈·구름·추적
불량·보정 결함)가 있는 프레임이다.

**순환 guard (빈 기호).** 한 필터에서 *가장 깊은* 프레임의 실현 깊이는 이 방식으로
측정할 수 없다. 그 필터 마스터 카탈로그의 어두운 꼬리는 바로 그 프레임들이 검출했기
때문에 존재하므로, 프레임이 **자기가 선택한 별로 채점**되어 깊이가 부풀려진다. 따라서
필터 마스터 자체의 90분위 한계 대비 여유가 0.7등급 이내인 프레임은
`depth_valid = False`로 표시해 RMS에서 제외했다(66중 9: M13 $R$ 전 5장, M13 $V$ 3장,
NGC 6811 $V$ 1장 — M13은 60초 노출로만 관측돼 마스터가 얕아 특히 취약하다). 이 guard는
형식적인 장치가 아니다: 제외된 M13 $R$ 프레임은 실현 17.26 vs 예측 15.59로 패널 (b)의
1:1선에서 눈에 띄게 벗어나며, guard가 없었다면 $+1.45$등급의 순환 부풀림이 예측 실패로
오독됐을 것이다. 독립적으로, guard를 통과한 주입 보정 프레임 6장(굵은 검은 테두리)에서
실현 $-$ 주입 $= +0.11$등급(평균), RMS 0.15등급으로 동일한 블렌드 매칭 편향과 일관된다.

**데이터.** M13·M67·NGC 6811의 단일 노출 66장(30–480초), 6개 필터($B$, $V$, $R$,
$g$, $r$, $i$), 단일 카메라(Moravian C3-61000, gain 0.689 $e^-$/ADU, 광자전달
실측), APEX가 raw→science 전 과정 처리. 순환 guard 통과 57장(M67 30, NGC 6811 20,
M13 7), FWHM 5.2–9.3 px, 프레임당 카탈로그 별 1139–2001개. 예측은 픽셀 데이터를 전혀
읽지 않는다 — 배경 잡음은 Step 7 측광 테이블에서, PSF peak fraction은 그 프레임의 밝고
깨끗한 별에서 가져온다. 프레임별 값:
`validation/paper/data_qc_depth/qc_depth_summary.csv`.
