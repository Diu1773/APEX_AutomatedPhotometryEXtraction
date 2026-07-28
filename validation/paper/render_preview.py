#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Render MANUSCRIPT markdown -> academic-paper HTML, self-contained.
Math via sub/sup/unicode, citations from references.bib, figures embedded as base64."""
import re, html, base64
from pathlib import Path

ROOT = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction\validation\paper")
SRC = ROOT / "논문작업" / "MANUSCRIPT_ko.md"
BIB = ROOT / "references.bib"
FIGDIR = ROOT / "figures"
OUT = ROOT / "MANUSCRIPT_ko_preview.html"   # 정본 자리(2026-07-09부터 이 경로)

# ---------- figures ----------
# 본문 「그림 N」 -> 실제 파일. 번호는 파이프라인 순서다.
# 파일명이 fig<숫자>_ 패턴이 아닌 새 그림들(실측 주입 완전도 등)도 여기서 직접 잡는다 —
# 예전 자동 매칭은 그 패턴만 찾아서 새 그림을 통째로 놓치고 옛 판을 붙이고 있었다.
FIGFILES = {
    1:  "fig_architecture.png",              # 2.1 설계 의도 — 작업 흐름·계층
    2:  "fig10_calibration.png",             # 3.2 검출기 보정
    3:  "fig11_detector.png",                # 3.3 검출기 특성화
    4:  "fig12_preproc_crosscheck.png",      # 3.4 ccdproc 교차
    5:  "fig13_cross_instrument.png",        # 3.5 교차기기 보정
    6:  "fig_completeness_realvssynth.png",  # 3.6 완전도(실측 주입)
    7:  "fig2_error_model.png",              # 3.8 오차 모형
    8:  "fig3_parameter_sweep.png",          # 3.9 민감도
    9:  "fig4_crosscheck_sep.png",           # 3.10 합성 독립 엔진
    10: "fig_iraf_crosscheck.png",           # 3.11 실측 독립 엔진(NGC6811 499별)
    11: "fig9_crowded_field.png",            # 3.12 밀집장
    12: "fig7_reference_crosscheck.png",     # 3.13 참조 목록
    13: "fig8_cmd_reproduction.png",         # 3.14 CMD 재현
    14: "fig6_qc_validation.png",            # 3.16 프레임 QC
    15: "fig_lc_yzboo.png",                  # 4   과학 적용(YZ Boo)
}
FIGMAP = {k: FIGDIR / v for k, v in FIGFILES.items() if (FIGDIR / v).exists()}
_missing = [f"{k}:{v}" for k, v in FIGFILES.items() if not (FIGDIR / v).exists()]
if _missing:
    print("[warn] 그림 파일 없음:", ", ".join(_missing))
def fig_uri(p):
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()

# descriptive captions, keyed by pipeline figure number
CAPTIONS = {
 1: r"APEX의 작업 흐름과 계층. 0–7단계는 두 모드가 공유하며, 측광이 끝난 뒤에야 CMD 모드와 LC 모드로 갈라진다. 각 단계는 관측자에게 **결정 하나**를 요구하고(어떤 보정 프레임을 쓸지, 이 프레임을 받아들일지, 검출 문턱을 어디에 둘지, 몇 장에서 보여야 별로 인정할지, 구경을 얼마로 할지) 정해진 경로에 검사 가능한 산출물을 남긴다. 계산은 Qt를 부르지 않는 핵심부에 있고 그래픽 계층은 그것을 부르기만 하므로, 3절의 화면 없는 검증이 곧 화면이 돌리는 코드를 시험한다. 각 단계 아래의 붉은 번호는 그 단계를 검증하는 절이다.",
 2: r"합성 bias·dark·flat을 주입해 되돌린 검출기 보정 검증. 실제 M13 프레임의 보정 전후로 비네팅과 고정패턴이 제거되고, L.A.Cosmic \citep{vandokkum2001}으로 우주선을 지운다. 보정 잔차는 offset $-0.016$ DN, MAD $3.18$ DN으로 주입 잡음 바닥과 일치.",
 3: r"광자전달곡선(photon-transfer curve)으로 자료에서 직접 잰 검출기 특성: gain $0.681$ e$^-$/ADU, 읽기잡음 $2.35$ e$^-$, 암전류 $0.0077$ e$^-$/s($R^2=0.998$). 헤더의 EGAIN 값은 $\approx16\times$ 틀리므로 쓰지 않는다.",
 4: r"각 보정 단계를 표준 파이썬 패키지 astropy \texttt{ccdproc}과 픽셀 대 픽셀로 비교. 마스터 bias·dark와 적용 세 단계가 비트동일($\Delta=0$), 전체 파이프라인도 $5\times10^{-4}$ DN 이내 — 읽기잡음보다 네 자릿수 이상 아래.",
 5: r"두 LCO 카메라의 raw를 APEX로 보정해 독립 파이프라인 BANZAI 산출물과 비교. QHY600 CMOS는 전체가 균일한 $+0.06$ e$^-$, 4-앰프 Sinistro CCD는 $\approx0.3\%$ 일치(사분면 패턴은 앰프별 조립의 차이). 보정 산술이 검출기를 넘어 일반화됨을 보인다.",
 6: r"실측 프레임에 인공별을 주입해 잰 검출 완전도. **(a)** 관측 조건을 대표하는 실측 단일 노출 세 장 — 어두운 하늘·양호한 시상(M67 $i$), 밝은 하늘·양호한 시상(NGC 6811 $R$), 밝은 하늘·불량 시상(M13 $V$) — 의 주입 등급 대비 회수율. 점은 구간별 회수율(Wilson 95\% 이항 구간)이고 **등급 공간에는 어떤 함수도 맞추지 않는다**. 50\% 깊이 $m_{50}=17.65,\ 15.64,\ 14.90$ 은 곡선이 0.5를 지나는 지점을 읽은 값이며, 깊이는 방법이 아니라 그 프레임의 하늘 밝기와 시상이 정한다. 회색 파선은 합성 verification 프레임. 아래 스트립은 가장 얕은 프레임에 실제로 주입된 별들의 컷아웃이다. **(b)** 각 별을 기대 peak-화소 S/N으로 재표현하면 깊이가 3.4등급에 걸쳐 벌어졌던 **일곱 프레임이 단일 곡선으로 붕괴**한다. 합동 표본의 오차함수 피팅은 $\mathrm{S/N}_{50}=4.0$, 프레임별 독립 판독은 $4.05\pm0.18$ 로 일치한다.",
 7: r"주입–되찾기로 검증한 측광 오차 모형. 보고한 오차가 실측 산포와 일치하며(pull 표준편차 $1.014$, $N=3404$), 구간별 RMS가 $\sigma_m=1.0857/\mathrm{SNR}$을 두 자릿수에 걸쳐 따른다.",
 8: r"파이프라인·관측 조건 스윕. 측광 산포는 구경 $1.2\times$FWHM에서 최소($0.058$등급)이고, 하늘밝기와 시상이 나빠지면 깊이가 얕아진다. 검출은 $2$–$6\sigma$ 문턱에 무관($<0.01$등급).",
 9: r"합성 참값에서 APEX와 독립 구현 엔진 Barbary SEP \citep{bertin1996}의 일치(별 95개): MAD $0.006$등급, Pearson $r=0.99995$. 같은 화소에 두 독립 엔진이 같은 플럭스를 낸다.",
 10: r"APEX가 raw에서부터 전부 줄인 NGC 6811 $V$ 프레임의 실제 별 499개를, 표 2에 맞춘 파라미터로 APEX 강제 측광과 IRAF `phot`(DAOPHOT)이 **같은 고정 좌표에서** 각각 잰 결과: MAD $9.7$ mmag, $r=0.99989$, 구간 중앙값 잔차가 어두운 쪽까지 평평하다 \citep{schechter1993}. 양쪽 모두 재중심을 껐으므로 잔차가 중심 잡기의 차이를 흡수할 수 없다. 이 불일치는 두 코드가 스스로 보고한 형식 오차의 제곱합근($27.8$ mmag)의 3분의 1 남짓이다.",
 11: r"두 구상성단(M5·M13) 코어에서 APEX의 두 측광법 — 강제 구경 대 PSF — 의 내부 일치. 최근접 이웃 거리에 대한 중앙 차이가 평평($\pm0.02$–$0.04$등급, 분해 한계 $\sim10$ px까지). Gaia에 의존하지 않는 내부 일관성 시험.",
 12: r"NGC 6811을 Gaia와 독립인 Pan-STARRS 1 \citep{chambers2016}에 교차대조. $B$에서 Gaia 변환 참조가 어두운 쪽으로 $+0.022$등급 흐르지만 APEX 자체는 PS1 대비 평평($+0.010$) — 어두운 쪽 편차는 APEX가 아니라 Gaia BP의 알려진 결함 \citep{riello2021}.",
 13: r"NGC 6811의 Johnson 색-등급도($V$ 대 $B-V$, 별 1921개). APEX 지상 측광이 Gaia 변환 우주 기반 참조와 주계열 능선 $19$ mmag로 일치하며, 독립 PS1 계에서도 같은 형태. CMD 산출물 자체를 검증(이소크론 맞추기는 별개).",
 14: r"주입 결함으로 만든 44-프레임 밤에서 자동 프레임 QC. 정상 24장을 오탐 없이 통과, 나쁜 시상·밝은 하늘·거짓 헤더 프레임을 모두 검출. 회색 투명도 손실만 놓치는데(균일 변화라 영상 통계에 안 잡힘), 이것이 2단계 측광-QC의 근거다.",
 15: r"LC 모드의 종단(end-to-end) 과학 산출물: 고진폭 $\delta$ Scuti 별 YZ Boötis를 raw 프레임에서 APEX만으로 줄였다. **(a)** 하룻밤(2026-03-28, $r$ 밴드, 5.2시간에 걸친 80점)을 문헌 주기로 접으면 출판된 톱니 곡선이 재현되고 진폭은 pk-pk $0.39$ 등급이다. **(b)** Lomb–Scargle 주기도. 단일밤 최고 봉우리는 $0.1046$ 일로 문헌값 $0.10409$ 일(파선)과 0.5\% 차이지만, 하루 간격의 두 밤을 병합하면 최고 봉우리가 $+1$ 주기/일 alias인 $0.0946$ 일로 옮겨가고 참 주기는 부봉우리로만 남는다. 이는 관측 창(window)의 성질이지 파이프라인의 것이 아니며, 숨기지 않고 보인다."
}

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
    x = re.sub(r"\^(\*|\w)", r"<sup>\1</sup>", x)
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
    s = re.sub(r"\0(\d+)\0", lambda m: holds[int(m.group(1))], s)
    return s

# ---------- parse ----------
lines = [ln for ln in SRC.read_text(encoding="utf-8").splitlines()
         if not re.match(r"\s*<!--.*-->\s*$", ln)]
out = []
para = []
sub_title = ""
sub_figs = []   # figure numbers referenced in current subsection (in order)

def flush_para(buf):
    if buf:
        txt = " ".join(buf).strip()
        if txt:
            cls = ""
            if txt.startswith("*[") or txt.startswith("["): cls = ' class="pending"'
            elif re.match(r"\*?영문 제출본", txt): cls = ' class="docnote"'
            out.append(f"<p{cls}>{inline(txt)}</p>")
    return []

emitted = set()
def emit_figs():
    global sub_figs
    # 번호가 붙은 절이면 그림을 받는다. 하위절(3.6)뿐 아니라 최상위 절(4. 과학 적용)도 —
    # 예전 조건은 "N.M" 만 허용해서 §4 의 그림이 통째로 빠졌다.
    if re.match(r"^\d+(\.\d+)?[\.\s]", sub_title):
        title = re.sub(r"^\d+(\.\d+)?\.?\s*", "", sub_title)
        for n in sub_figs:
            if n in FIGMAP and n not in emitted:
                emitted.add(n)
                body = inline(CAPTIONS[n]) if n in CAPTIONS else html.escape(title)
                cap = f"<b>그림 {n}.</b> {body}"
                out.append(f'<figure><img alt="그림 {n}" src="{fig_uri(FIGMAP[n])}">'
                           f'<figcaption>{cap}</figcaption></figure>')
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

CSS = r"""
:root{ color-scheme:light; --ink:#000; --muted:#2c2c2c; --rule:#000; --code:#eef0f2; }
*{box-sizing:border-box}
html,body{background:#fff;}
/* 한글 본문 폰트: Noto Serif KR(본명조 계열)을 1순위로 둔다. 예전 스택은 미설치
   Nanum Myeongjo 를 먼저 찾다가 Batang 으로 떨어져 화면에서 구식으로 보였다. */
body{margin:0;color:var(--ink);
  font-family:"Georgia","Times New Roman",Times,"Noto Serif KR",
              "Apple SD Gothic Neo","Malgun Gothic",serif;
  line-height:1.78;font-size:17px;text-rendering:optimizeLegibility;
  word-break:keep-all;overflow-wrap:break-word;
  -webkit-font-smoothing:antialiased;}
.wrap{max-width:44rem;margin:0 auto;padding:3.2rem 1.5rem 6rem;}

h1.title{font-size:1.62rem;line-height:1.3;text-align:center;text-wrap:balance;
  font-weight:700;margin:0 0 1rem;}
.docnote{font-size:.82rem;color:var(--muted);text-align:center;font-style:italic;
  margin:0 auto 2.2rem;max-width:33rem;line-height:1.5;}

h2{font-size:1.15rem;font-weight:700;margin:1.9rem 0 .55rem;text-wrap:balance;}
#abstract{text-align:center;font-size:1.05rem;margin:1.6rem 0 .5rem;}
h3{font-size:1.02rem;font-weight:700;margin:1.4rem 0 .4rem;font-style:italic;text-wrap:balance;}
h4{font-size:.97rem;font-weight:700;font-style:italic;margin:1.1rem 0 .3rem;}

/* 한글은 양끝맞춤하면 어절 간격이 벌어져 읽기 나빠진다. 왼쪽 정렬 + 문단 간격. */
p{margin:0 0 .95rem;text-align:left;}
h1+p,h2+p,h3+p,h4+p,blockquote+p,figure+p,.tw+p,p.docnote,p.pending{text-indent:0;}
p.pending{color:var(--muted);font-style:italic;}
strong{font-weight:700;}
code{font-family:"Courier New",ui-monospace,monospace;background:var(--code);
  padding:.02em .28em;font-size:.86em;}
.math{font-family:"Cambria Math","Times New Roman",serif;white-space:nowrap;}
.math sub,.math sup{font-size:.74em;}

/* abstract: centered heading, narrowed justified block (LaTeX style) */
h2#abstract + p{margin:.2rem 1.6rem 0;font-size:.92rem;line-height:1.5;text-indent:0;}

blockquote{margin:.9rem 1.6rem;padding-left:.9rem;border-left:2px solid #999;
  color:var(--muted);font-size:.95rem;}

.tw{overflow-x:auto;margin:1rem 0;}
table{border-collapse:collapse;width:100%;font-size:.8rem;
  font-variant-numeric:tabular-nums;line-height:1.4;margin:0 auto;}
thead th{border-top:1.3px solid #000;border-bottom:.8px solid #000;
  text-align:left;padding:.4rem .55rem;font-weight:700;vertical-align:bottom;}
tbody td{padding:.34rem .55rem;vertical-align:top;}
tbody tr:last-child td{border-bottom:1.3px solid #000;}

figure{margin:1.6rem 0;text-align:center;}
figure img{max-width:100%;height:auto;background:#fff;}
figcaption{font-size:.82rem;color:var(--ink);margin:.5rem auto 0;max-width:38rem;
  text-align:left;line-height:1.45;}
figcaption b{font-weight:700;}

::selection{background:#ccd8ee;}
a{color:#7a1010;text-decoration:none;}
"""
HTML = f"""<title>APEX — 검증된 그래픽 측광 파이프라인 (국문 프리뷰)</title>
<style>{CSS}</style>
<div class="wrap">
{BODY}
</div>"""
OUT.write_text(HTML, encoding="utf-8")
print("wrote", OUT, round(len(HTML)/1e6,2), "MB | citations", len(LAB), "| figs", len(FIGMAP))
