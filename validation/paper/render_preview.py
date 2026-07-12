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
OUT = Path(r"C:\Users\bmffr\AppData\Local\Temp\claude\C--Users-bmffr-Desktop-Result-Automated-Photometry-EXtraction\671e5e9e-3ddb-4c65-9711-7a874bcba688\scratchpad\paper_ko.html")

# ---------- figures ----------
# manuscript figure number (pipeline order) -> source fig file number (old order)
NEWFIG = {1:10, 2:11, 3:12, 4:13, 5:1, 6:2, 7:3, 8:4, 9:5, 10:9, 11:7, 12:8, 13:6}
_filemap = {}
for p in sorted(FIGDIR.glob("fig*.png")):
    m = re.match(r"fig(\d+)_", p.name)
    if m:
        _filemap[int(m.group(1))] = p
FIGMAP = {k: _filemap[v] for k, v in NEWFIG.items() if v in _filemap}
def fig_uri(p):
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()

# descriptive captions, keyed by pipeline figure number
CAPTIONS = {
 1: r"합성 bias·dark·flat을 주입해 되돌린 검출기 보정 검증. 실제 M13 프레임의 보정 전후로 비네팅과 고정패턴이 제거되고, L.A.Cosmic \citep{vandokkum2001}으로 우주선을 지운다. 보정 잔차는 offset $-0.016$ DN, MAD $3.18$ DN으로 주입 잡음 바닥과 일치.",
 2: r"광자전달곡선(photon-transfer curve)으로 자료에서 직접 잰 검출기 특성: gain $0.681$ e$^-$/ADU, 읽기잡음 $2.35$ e$^-$, 암전류 $0.0077$ e$^-$/s($R^2=0.998$). 헤더의 EGAIN 값은 $\approx16\times$ 틀리므로 쓰지 않는다.",
 3: r"각 보정 단계를 표준 파이썬 패키지 astropy \texttt{ccdproc}과 픽셀 대 픽셀로 비교. 마스터 bias·dark와 적용 세 단계가 비트동일($\Delta=0$), 전체 파이프라인도 $5\times10^{-4}$ DN 이내 — 읽기잡음보다 네 자릿수 이상 아래.",
 4: r"두 LCO 카메라의 raw를 APEX로 보정해 독립 파이프라인 BANZAI 산출물과 비교. QHY600 CMOS는 전체가 균일한 $+0.06$ e$^-$, 4-앰프 Sinistro CCD는 $\approx0.3\%$ 일치(사분면 패턴은 앰프별 조립의 차이). 보정 산술이 검출기를 넘어 일반화됨을 보인다.",
 5: r"(a) 등급에 대한 완전도 — 인공별을 주입해 실제 검출기로 되찾은 비율(점, Wilson 95\%). 레퍼런스(AutoPhOT·Masci)를 따라 등급 공간엔 곡선을 맞추지 않고 되찾음이 각 비율을 지나는 등급을 읽어 낸다: $m_{50}=17.59$, $m_{90}=17.39$, $m_{10}=17.91$등급. (b) 같은 되찾음을 S/N으로 그리면 오차함수(erf) 검출확률 모형에 모이며 50\%가 S/N $\approx7.4$ — 검출을 지배하는 것은 등급이 아니라 S/N이다 \citep{masci2011, kashyap2010, brennan2022}.",
 6: r"주입–되찾기로 검증한 측광 오차 모형. 보고한 오차가 실측 산포와 일치하며(pull 표준편차 $1.014$, $N=3404$), 구간별 RMS가 $\sigma_m=1.0857/\mathrm{SNR}$을 두 자릿수에 걸쳐 따른다.",
 7: r"파이프라인·관측 조건 스윕. 측광 산포는 구경 $1.2\times$FWHM에서 최소($0.058$등급)이고, 하늘밝기와 시상이 나빠지면 깊이가 얕아진다. 검출은 $2$–$6\sigma$ 문턱에 무관($<0.01$등급).",
 8: r"합성 참값에서 APEX와 독립 구현 엔진 Barbary SEP \citep{bertin1996}의 일치(별 95개): MAD $0.006$등급, Pearson $r=0.99995$. 같은 화소에 두 독립 엔진이 같은 플럭스를 낸다.",
 9: r"완전-APEX로 처리한 NGC 6811 $V$ 프레임에서 APEX 강제 측광과 IRAF/DAOPHOT \citep{stetson1987}의 일치(별 498개): MAD $0.009$등급, $r=0.99984$. 같은 좌표·같은 파라미터라 잔차는 두 플럭스 적분기의 차이만 격리한다.",
 10: r"두 구상성단(M5·M13) 코어에서 APEX의 두 측광법 — 강제 구경 대 PSF — 의 내부 일치. 최근접 이웃 거리에 대한 중앙 차이가 평평($\pm0.02$–$0.04$등급, 분해 한계 $\sim10$ px까지). Gaia에 의존하지 않는 내부 일관성 시험.",
 11: r"NGC 6811을 Gaia와 독립인 Pan-STARRS 1 \citep{chambers2016}에 교차대조. $B$에서 Gaia 변환 참조가 어두운 쪽으로 $+0.022$등급 흐르지만 APEX 자체는 PS1 대비 평평($+0.010$) — 어두운 쪽 편차는 APEX가 아니라 Gaia BP의 알려진 결함 \citep{riello2021}.",
 12: r"NGC 6811의 Johnson 색-등급도($V$ 대 $B-V$, 별 1921개). APEX 지상 측광이 Gaia 변환 우주 기반 참조와 주계열 능선 $19$ mmag로 일치하며, 독립 PS1 계에서도 같은 형태. CMD 산출물 자체를 검증(이소크론 맞추기는 별개).",
 13: r"주입 결함으로 만든 44-프레임 밤에서 자동 프레임 QC. 정상 24장을 오탐 없이 통과, 나쁜 시상·밝은 하늘·거짓 헤더 프레임을 모두 검출. 회색 투명도 손실만 놓치는데(균일 변화라 영상 통계에 안 잡힘), 이것이 2단계 측광-QC의 근거다.",
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
GREEK = {r"\\sigma":"σ", r"\\tau":"τ", r"\\theta":"θ", r"\\mu":"μ", r"\\delta":"δ",
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
    if re.match(r"^\d+\.\d", sub_title):   # only numbered §3 subsections carry figures
        title = re.sub(r"^\d+\.\d+\s*", "", sub_title)
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
body{margin:0;color:var(--ink);
  font-family:"Georgia","Times New Roman",Times,"Nanum Myeongjo","Batang","Apple SD Gothic Neo",serif;
  line-height:1.66;font-size:17.5px;text-rendering:optimizeLegibility;
  -webkit-font-smoothing:antialiased;}
.wrap{max-width:45rem;margin:0 auto;padding:3.4rem 1.7rem 6rem;}

h1.title{font-size:1.62rem;line-height:1.3;text-align:center;text-wrap:balance;
  font-weight:700;margin:0 0 1rem;}
.docnote{font-size:.82rem;color:var(--muted);text-align:center;font-style:italic;
  margin:0 auto 2.2rem;max-width:33rem;line-height:1.5;}

h2{font-size:1.15rem;font-weight:700;margin:1.9rem 0 .55rem;text-wrap:balance;}
#abstract{text-align:center;font-size:1.05rem;margin:1.6rem 0 .5rem;}
h3{font-size:1.02rem;font-weight:700;margin:1.4rem 0 .4rem;font-style:italic;text-wrap:balance;}
h4{font-size:.97rem;font-weight:700;font-style:italic;margin:1.1rem 0 .3rem;}

p{margin:0;text-align:justify;text-justify:inter-word;text-indent:1.6em;}
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
