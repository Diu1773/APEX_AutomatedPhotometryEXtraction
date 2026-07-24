# Figure — Detection completeness by real-frame injection, and its collapse onto a single detection law

*(supersedes `fig1_completeness.md`, which described the earlier synthetic-only version)*

**Figure X.** Detection completeness of the APEX pipeline measured by artificial-star
injection into real science frames (DAOPHOT `ADDSTAR` heritage; cf. DES *Balrog*,
AutoPhOT App. D).
**(a)** Recovery fraction versus injected magnitude (count-rate instrumental system,
ZP = 25; 0.25 mag bins, ≥15 stars per bin) for three representative fully-APEX-reduced
frames spanning the observed range of conditions: a dark-sky, sharp frame (M67 $i$,
background noise $\sigma = 5\,e^-$ px$^{-1}$, FWHM 5.2 px), a bright-sky, sharp frame
(NGC 6811 $R$, $\sigma = 30\,e^-$, 5.3 px), and a bright-sky, soft frame (M13 $V$,
$\sigma = 30\,e^-$, 7.4 px). Stars are injected with the empirical PSF measured from
each frame, with source Poisson noise, through the identical detection path as the
science reduction; injections confounded by a pre-existing source are excluded.
Vertical dotted lines mark the 50% completeness depths, read off where each binned
curve crosses 0.5 (no functional fit in magnitude space): $m_{50} = 17.65$, $15.64$,
and $14.90$, respectively. The 50% depth is a property of the frame — set by sky
brightness and seeing — not of the method; the grey dashed curve shows the synthetic
verification frame (known truth; $\sigma \simeq 8\,e^-$, FWHM 3.4 px, $m_{50}=17.59$),
which the best real frame matches. The bottom strip shows six of the actual injected
stars in the shallowest frame (M13 $V$) on a shared greyscale stretch, straddling the
transition: stars brighter than $m_{50}$ are recovered (blue borders), fainter ones
vanish into the sky noise (orange).
**(b)** The same injections for all seven real frames used in this test — spanning a
factor 11 in background noise ($\sigma = 5$–$58\,e^-$) and 1.7 in seeing (FWHM
5.2–9.0 px), i.e. 3.4 mag in depth — re-expressed as a function of each star's
expected peak-pixel signal-to-noise ratio, $\mathrm{S/N_{peak}} = F\,p_{\rm peak} /
\sigma$, where $F$ is the injected flux in electrons, $p_{\rm peak}$ the peak fraction
of that frame's injection PSF, and $\sigma$ that frame's background noise. All seven
curves collapse onto a single completeness law. The black dashed curve is a
three-parameter error-function fit to the pooled sample,
$C = \tfrac{A}{2}\left[1 + \mathrm{erf}\!\left((\log \mathrm{S/N} - \log
\mathrm{S/N_{50}})/(\sqrt{2}\,w)\right)\right]$, giving $\mathrm{S/N_{50}} = 4.0$
($A = 0.99$, $w = 0.07$ dex); the per-frame read-off values agree to
$\mathrm{S/N_{50}} = 4.05 \pm 0.18$ (frame-to-frame scatter, 4%), consistent with the
pipeline's 3.2$\sigma$ matched-filter detection threshold. Absolute depth is therefore
a frame property, while the quantity the pipeline controls — the S/N at which sources
are recovered — is invariant across frame conditions; magnitude-space completeness
curves are this single law translated by each frame's $\sigma \cdot \mathrm{FWHM}^2$.

**Data.** Seven frames of three clusters (M13, M67, NGC 6811), single camera
(Moravian C3-61000, gain 0.689 $e^-$/ADU from photon-transfer measurement), reduced
raw→science entirely by APEX; 60 trials × 50 stars per frame (21,000 total across the
family, 3,000 per frame), seeds fixed for reproducibility. Per-frame numbers:

| frame | $\sigma$ ($e^-$/px) | FWHM (px) | $m_{50}$ | $\mathrm{S/N_{50}}$ |
|---|---|---|---|---|
| M67 $i$ | 5.1 | 5.2 | 17.65 | 4.29 |
| M67 $r$ | 17.0 | 6.9 | 15.79 | 3.97 |
| M67 $g$ | 24.8 | 9.0 | 14.91 | 3.82 |
| M13 $R$ | 22.4 | 5.3 | 15.81 | 4.15 |
| M13 $V$ | 30.1 | 7.4 | 14.90 | 3.92 |
| NGC 6811 $R$ | 30.3 | 5.3 | 15.64 | 4.29 |
| NGC 6811 $R$ (soft) | 58.3 | 7.2 | 14.24 | 3.94 |

---

## 국문 (MANUSCRIPT_ko용)

**그림 X.** 실측 과학 프레임에 인공별을 주입해 측정한 APEX 파이프라인의 검출
완전도(detection completeness) (DAOPHOT `ADDSTAR` 계보; DES *Balrog*, AutoPhOT App. D 참조).
**(a)** 관측 조건 범위를 대표하는 세 실측 프레임 — 어두운 하늘·양호한 시상(M67 $i$),
밝은 하늘·양호한 시상(NGC 6811 $R$), 밝은 하늘·불량 시상(M13 $V$) — 의 주입 등급 대비
회수율(count-rate 기기등급, ZP = 25; 0.25등급 bin). 각 프레임에서 측정한 경험적
PSF로 광자 잡음을 포함해 주입하며, 과학 리덕션과 동일한 검출 경로를 통과한다. 점선은
50% 완전도 깊이로, 등급 공간의 함수 피팅 없이 곡선이 0.5를 지나는 지점을 읽은 값이다:
각각 $m_{50}=17.65,\ 15.64,\ 14.90$. **깊이는 하늘 밝기와 시상이 정하는 프레임 속성**
이며, 가장 좋은 실측 프레임은 합성 verification 프레임(회색 파선, $m_{50}=17.59$)과
같은 깊이에 도달한다. 하단 스트립: 가장 얕은 프레임(M13 $V$)에 실제로 주입된 별 6개
(공유 스트레치) — $m_{50}$보다 밝으면 회수(파란 테두리), 어두우면 하늘 잡음에 소실(주황).
**(b)** 이 검증의 핵심: 일곱 실측 프레임 전체(배경 잡음 11배, 시상 1.7배, 깊이 3.4등급
스팬)의 주입별을 기대 peak-화소 S/N으로 재표현하면 **모든 곡선이 단일 완전도 법칙으로
붕괴**한다. 합동 표본의 오차함수(erf) 피팅은 $\mathrm{S/N_{50}}=4.0$, 프레임별 판독값은
$4.05\pm0.18$(프레임 간 4%)로, 파이프라인의 3.2$\sigma$ 매치드필터 검출 임계와
일치한다. 즉 절대 깊이는 프레임 속성이고, 파이프라인이 통제하는 양 — 소스가 회수되는
S/N — 은 프레임 조건과 무관하게 불변이다. 등급 공간의 완전도 곡선들은 이 단일 법칙이
각 프레임의 $\sigma \cdot \mathrm{FWHM}^2$만큼 평행이동된 것이다.
