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
    s = re.sub(r"`([^`]+)`", lambda m: hold(f'<code>{html.escape(m.group(1))}</code>'), s)
    s = re.sub(r"\$([^$]+)\$", lambda m: hold(render_math(m.group(1))), s)
    s = repl_citations(s)
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
                cap = f"<b>그림 {n}.</b> {html.escape(title)}"
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
