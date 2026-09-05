# APEX 지금 상태 — 여기서 시작한다

**갱신 2026-09-05.** 오래 떠나 있다가 돌아왔으면 이 문서 하나만 읽으면 된다.
왜 만드는지는 [RESEARCH_FRAME.md](RESEARCH_FRAME.md), 단계 번호와 출력 경로는
[PLAN.md](PLAN.md), 쓰는 법은 [README.md](README.md) 에 있다.

## 지금 어디까지 왔나

**소프트웨어는 끝났다.** 원본 사진에서 성단 색-등급도와 여러 밤 광도곡선까지 한
프로그램으로 간다. 성단 다섯 개와 변광성 두 개에서 처음부터 끝까지 통과했고,
시험 **1,426 개**가 세 운영체제와 파이썬 세 판에서 돈다(2026-08-27 기준).

**남은 것은 논문이다.** 투고처는 PASP 로 정해져 있다(2026-08-24).

---

## 목표 셋

[RESEARCH_FRAME.md](RESEARCH_FRAME.md) 의 「아직 답하지 못한 것」이 그대로 목표다.
여기에는 **무엇을 하면 끝난 것으로 치는지**만 덧붙인다.

### 하나 — 단계의 오차가 최종 결과에 얼마나 남는가

**지금**: 넷에 참값 시험이 있다. 검출기 보정(합성 raw 로 참 프레임 복원, 계통
1 DN 미만) · 강제측광(주입 10,993 개, 50 % 회수 지점 편차 중앙 −0.002, 최대
0.033 등급) · 구경보정(주입 18,000 개, 45 프레임) · 검출 완전도(7 프레임,
S/N_50 = 4.05 ± 0.18).

**정정 (2026-09-06)**: 앞서 이 자리에 「주기 탐색 · SYSREM · 등시선에 참값 시험이
없다」고 적었는데 **틀렸다. 셋 다 이미 있다.**

| 무엇 | 어디 | 무엇을 판정하나 |
|---|---|---|
| 주기 탐색 | `tests/test_lc_period_step.py` 의 `test_it_recovers_a_known_period` | 0.1 일을 넣고 0.002 일 이내로 되찾는가 |
| SYSREM | `tests/test_sysrem.py` 의 `test_sysrem_recovers_single_systematic` | 주입한 계통 성분과 상관 0.95 를 넘는가 |
| 등시선 | `tests/test_isochrone_mcmc.py` 의 `test_synthetic_recovery_truth_in_credible_interval` | 참 나이가 신용구간 안이고 로그나이 편차가 0.15 미만인가 |

이것들이 단위 시험이라 `docs/audit/APEX_ENGINE_SCORECARD.md` 에는 안 잡힌다.

**그러면 진짜 빈 곳은 둘이다.**

**하나, 저 시험들은 통과와 실패만 낸다.** 「0.002 일 이내」는 작동을 증명하지만
「주기 정확도가 얼마다」를 말하지 않는다. 논문에 쓰려면 조건을 쓸어 가며 잰 곡선이
있어야 한다.

**둘, 단계의 오차가 최종에 얼마나 남는지를 잰 것이 없다.** 이것이 RESEARCH_FRAME 이
실제로 묻는 것이다.

**끝난 것으로 치는 조건**: 등급 오차를 정해진 크기로 얹었을 때 최종 성단 파라미터와
주기가 얼마나 흔들리는지를 표로 낸다.

**착수함 (2026-09-06)**: CMD 갈래를 먼저 잰다.
`validation/error_propagation/run_isochrone_propagation.py` 가 참 파라미터로 합성
성단을 만들고 측광 오차 0.005~0.080 등급을 여섯 단계로 얹어 열 번씩 다시 적합한다.
오차값은 APEX 가 인공별로 실제로 잰 흩어짐 범위에 맞췄다. 파이프라인 전체가 아니라
**이소크론 적합 층 하나만** 돈다.

### 둘 — 다른 카메라에서도 되는가

**지금**: 검출기 보정만 세 기기다(Moravian C3-61000 · QHY600 · Sinistro).
**측광부터 끝까지는 Moravian 한 대뿐이다.** `docs/audit/APEX_ENGINE_SCORECARD.md`
의 범용성 열이 Step 0 을 빼면 전부 비어 있다.

**끝난 것으로 치는 조건**: 다른 기기 자료로 Step 1–7 과 CMD 또는 LC 한 갈래를
끝까지 돌리고, 등급 회수 편차를 위와 같은 방식으로 낸다.

**길이 열렸다 (2026-09-05).** 갤럭시북이 놀고 있어서 내려받기와 다른 기기 대조를
그쪽에서 돌릴 수 있다. 그리고 **LCO 공개 아카이브에 M67 raw 프레임이 있다.**

APEX 는 이미 M67 을 Moravian C3-61000 으로 처리해 두었다
(`validation/paper/data_realframe_M67g_broad` · `M67r_mid` · `M67i`).
LCO 의 `sq30`(0.4 m, QHY600 계열)이 같은 M67 을 **V 와 rp** 로 찍어 두었고,
`rp` 가 APEX 의 `r` 과 겹친다. **같은 성단 · 같은 필터 · 다른 기기**이므로
표준성야보다 직접적인 대조다.

**받는 도구를 만들어 두었다** — `scripts/fetch_lco.py`. 인증이 필요 없고 표준
라이브러리만 쓰므로 갤럭시북에 레포만 있으면 그대로 돈다.

**가장 좋은 후보는 MuSCAT3 다.** 2 m 망원경(ogg 2m0a)에 붙은 `ep02`~`ep05` 가
**gp · rp · ip · zs 네 밴드를 동시에** 찍는다. APEX 의 Moravian M67 이 g · r · i
이므로 **세 밴드가 겹친다.** 2021-03 에 M67 을 찍은 raw 가 밴드마다 131 장 있다
(proposal `MuSCAT Commissioning`).

보정 프레임도 다 있다 — `ep02` 기준 BIAS 31,360 · DARK 9,694 · SKYFLAT 17,580 이
공개되어 있고, 그 밤의 master bias 도 이름이 붙어 있다
(`ogg2m001-ep02-20210322-bias-MUSCAT_SLOW-bin1x1.fits.fz`).

무엇이 있는지 먼저 본다.

    python -X utf8 scripts/fetch_lco.py --target M67 --level 0 --limit 300 --list

과학 프레임을 받는다(밴드마다 한 번씩).

    python -X utf8 scripts/fetch_lco.py --target M67 --instrument ep02         --level 0 --start 2021-03-01 --end 2021-04-01 --limit 40         --out E:/APEX_validation/external/LCO_M67_muscat3/rp

보정 프레임을 받는다.

    python -X utf8 scripts/fetch_lco.py --instrument ep02 --config-type BIAS         --level 91 --start 2021-03-20 --end 2021-03-23 --limit 5         --out E:/APEX_validation/external/LCO_M67_muscat3/cal

`--level 0` 이 raw 이고 `91` 이 BANZAI 가 처리한 것이다. **둘 다 받아 두면** APEX 의
Step 0 을 BANZAI 와 화소 단위로 견줄 수 있다(Fig 13 이 그 방식이고 QHY600 에서
0.06 e- 로 맞았다).

**크기 어림**: MuSCAT3 는 2048×2048 이라 raw 한 장이 8 MB 안팎이다. 세 밴드 40 장씩
이면 약 1 GB 이고 보정까지 더해도 1.5 GB 안쪽이다.

**하는 순서**: raw 와 보정을 받는다 → Step 0 으로 보정한다 → Step 1–7 을 돌린다 →
CMD 10 까지 간다 → 등급 회수 편차를 Moravian M67 과 같은 방식으로 낸다.

**미리 볼 것 하나**: 2 m 망원경의 10 초 노출은 M67 의 밝은 별을 포화시킬 수 있다.
60 초짜리도 밴드마다 10 장 있으니 둘 다 받아 보고 포화 상태를 먼저 확인한다.

**남는 자리**: `E:\APEX_validation\external\stetson\` 이 비어 있다. Stetson
표준성야는 여러 망원경의 발표 등급이 있어 절대 영점까지 물을 수 있으므로, LCO
대조가 끝난 뒤의 다음 단계로 둔다.

### 셋 — 남이 받아서 쓰는가

**지금**: 답 없음. 이것이 궁극 질문의 본체다.

**끝난 것으로 치는 조건**: 아직 정해지지 않았다. 사용자가 정해야 한다.

---

## 저절로 도는 것

**시험** — GitHub Actions 가 `main` 에 밀 때마다 돌고, **매주 월요일에도 돈다**
(`.github/workflows/tests.yml`). 세 운영체제 × 파이썬 3.10/3.11/3.12 조합이다.
실패하면 GitHub 이 알린다.

그게 전부다. 아래 둘은 저절로 못 돈다.

## 사람이나 세션이 있어야 하는 것

**벤치마크** — `benchmark/overnight.py` 가 단계마다 저장하며 밤새 돌도록 되어
있지만, **E 드라이브의 실제 관측 자료와 놀고 있는 기계**가 필요하다. CI 에서는
못 돈다. 기계가 빌 때 이 한 줄을 돌린다.

    .venv-deploy/Scripts/python.exe -X utf8 benchmark/overnight.py

**문헌 조사** — **자동화하면 안 된다.** 2026-09-04~05 에 검색 요약을 믿었다가
세 번 틀렸다(`Main/FAILURES.md` F-253, `Main/OPERATOR.md` C-168). 검색은 소프트웨어
이름과 원논문을 섞어서, 후대 논문의 내용을 원논문의 것처럼 돌려준다. **논문 X 가
Y 를 했다고 적으려면 X 를 직접 열어야 한다.** 검색은 논문을 찾는 데만 쓴다.

자동화할 수 있는 것은 판단이 아니라 **수집**이다 — arXiv 감시 목록을 만들어
후보를 쌓아 두고, 읽는 것은 사람이나 세션이 한다. 아직 안 만들었다.

---

## 오늘(2026-09-05) 문헌조사가 답한 것

`validation/paper/논문작업/LIT_INJECTION_PRACTICE_20260904.md` 에 원문 인용과 함께
다 있다. 요지만 옮기면 이렇다.

**목적이 APEX 와 같은 도구들은 검증을 거의 안 한다.** AstroImageJ (2017) 는
IRAF 와 맞췄다고 한 문장 적고 수치를 안 냈고, AutoPhOT (2022) 는 인공별을
한계등급 계산에만 썼고, STDWeb (2024) 은 검증 절이 없고, PhoPS (2026) 는 Gaia 로
측성만 잰다.

**반면 성단 측광 쪽은 1988 년부터 인공별로 계통 오차를 잰다.** Stetson & Harris
(1988) 이 시작했고 Jang (2023) 이 지금도 같은 것을 한다. Stetson 은 DAOPHOT 에
`ADDSTAR` 를 넣고 매뉴얼에 「넣은 것과 나온 것을 견주어 photometric accuracy 를
추정하라」고 적었다.

**APEX 는 그 둘 사이에 있다** — 목적은 첫째 갈래이고 산출물은 둘째 갈래다.

**IRAF 를 기준으로 삼아도 된다.** 다만 그것이 주는 것은 등급의 정당성이고,
**오차막대는 안 따라온다.** DAOPHOT 의 보고 오차를 두고 네 논문이 갈리며 하나는
방향까지 반대다.

---

## 문서 지도 (전체 156 개)

| 어디 | 몇 개 | 무엇 |
|---|---|---|
| 루트 | 10 | 진입·계약·상태. **이 문서 · README · PLAN · RESEARCH_FRAME · ARCHITECTURE** |
| `docs/audit/` | 22 | 엔진 성적표 · 선행연구 · 감사 기록 |
| `docs/` | 18 | 설계 문서 |
| `docs/manual/` | 8 | 한국어 실사용 매뉴얼(스크린샷 포함) |
| `validation/paper/논문작업/` | 39 | **논문 작업 문서.** 목차 · 원고 · 심사 · 문헌 |
| `validation/paper/captions/` | 20 | 그림 캡션 |
| `validation/*/` | 30 | 검증 실행 기록(엔진 대조 · 인공별 · 기기 교차) |

**세션 상태**는 `TRACK.md`(소프트웨어)와 `TRACK_PAPER.md`(논문)에 있다. 둘 다
길어서 관련된 절만 읽으면 된다.

---

## 돌아왔을 때 첫 명령

    git log --oneline -15
    .venv-deploy/Scripts/python.exe -m pytest tests -q

시험이 1,426 개 통과하면 소프트웨어는 그대로다. 그 다음은 위 「목표 셋」에서
고른다.

## 아직 안 고친 문서

`ARCHITECTURE.md`(2026-06-11) · `GOOGLE_SITE_CONTENT.md`(2026-06-30) ·
`CHANGELOG.md`(2026-07-24) 가 낡았다. `CHANGELOG.md` 의 `[Unreleased]` 에는
2026-07-24 이후 여섯 주치가 안 들어가 있다.
