# -*- coding: utf-8 -*-
"""Render MANUSCRIPT.tex -> self-contained journal-style HTML (arXiv-white).
Bounded LaTeX subset: sections, figures, \\citep/\\ref, inline math, emphasis."""
import re, html, base64
from pathlib import Path

ROOT = Path(r"C:\Users\bmffr\Desktop\Result\Automated_Photometry_EXtraction\validation\paper")
SRC = ROOT / "MANUSCRIPT.tex"
BIB = ROOT / "references.bib"
FIGDIR = ROOT / "figures"
OUT = Path(__file__).parent / "paper_draft.html"

raw = SRC.read_text(encoding="utf-8")

# ---------- strip comments (%, but not \%) and preamble ----------
lines = []
for ln in raw.split("\n"):
    ln = re.sub(r"(?<!\\)%.*$", "", ln)
    lines.append(ln)
tex = "\n".join(lines)

title = re.search(r"\\title(?:\[[^\]]*\])?\{(.+?)\}", tex, re.S).group(1).strip()
abstract = re.search(r"\\begin\{abstract\}(.*?)\\end\{abstract\}", tex, re.S).group(1).strip()
body = tex[re.search(r"\\end\{abstract\}", tex).end():]
body = body[:re.search(r"\\bibliography\{", body).start()] if re.search(r"\\bibliography\{", body) else body

# ---------- bib ----------
def parse_bib(text):
    d = {}
    for m in re.finditer(r"@\w+\s*\{\s*([^,]+),(.*?)\n\}", text, re.S):
        key = m.group(1).strip(); bd = m.group(2)
        au = re.search(r"author\s*=\s*[{\"](.+?)[}\"]\s*,?\s*\n", bd, re.S)
        yr = re.search(r"year\s*=\s*[{\"]?\s*(\d{4})", bd)
        ti = re.search(r"title\s*=\s*[{\"](.+?)[}\"]\s*,?\s*\n", bd, re.S)
        jr = re.search(r"journal\s*=\s*[{\"](.+?)[}\"]", bd, re.S)
        vol = re.search(r"volume\s*=\s*[{\"]?([^,}\"]+)", bd)
        year = yr.group(1) if yr else ""
        def last(a):
            a = a.strip().strip("{}")
            if "," in a: return a.split(",")[0].strip().strip("{}")
            if a.endswith("Collaboration"): return a
            p = a.split(); return p[-1] if p else a
        authors = re.split(r"\s+and\s+", au.group(1).strip()) if au else ["?"]
        l0 = last(authors[0])
        if len(authors) == 1: short = f"{l0} {year}"
        elif len(authors) == 2: short = f"{l0} & {last(authors[1])} {year}"
        else: short = f"{l0} et al. {year}"
        full_au = ", ".join(last(a) for a in authors[:3]) + (" et al." if len(authors) > 3 else "")
        ti_t = re.sub(r"[{}]", "", ti.group(1)).strip() if ti else ""
        jr_t = re.sub(r"[{}]", "", jr.group(1)).strip() if jr else ""
        d[key] = dict(short=short, sort=(l0.lower(), year), full=full_au, year=year,
                      title=ti_t, journal=jr_t, vol=vol.group(1) if vol else "")
    return d
bib = parse_bib(BIB.read_text(encoding="utf-8"))
cited = set()

# ---------- pass 1: number sections & figures, build label map ----------
labelnum = {}
sec = 0; sub = 0; fignum = 0
tokens = re.finditer(
    r"\\section\{(?P<s>[^}]*)\}(?:\s*\\label\{(?P<sl>[^}]*)\})?"
    r"|\\subsection\{(?P<ss>[^}]*)\}(?:\s*\\label\{(?P<ssl>[^}]*)\})?"
    r"|\\begin\{figure\*?\}(?P<fig>.*?)\\end\{figure\*?\}", body, re.S)
for m in tokens:
    if m.group("s") is not None:
        sec += 1; sub = 0
        if m.group("sl"): labelnum[m.group("sl")] = str(sec)
    elif m.group("ss") is not None:
        sub += 1
        if m.group("ssl"): labelnum[m.group("ssl")] = f"{sec}.{sub}"
    elif m.group("fig") is not None:
        fignum += 1
        fl = re.search(r"\\label\{([^}]*)\}", m.group("fig"))
        if fl: labelnum[fl.group(1)] = str(fignum)

# ---------- inline transforms ----------
def math(s):
    s = s.replace(r"\times", "×").replace(r"\pm", "±").replace(r"\approx", "≈")
    s = s.replace(r"\sim", "~").replace(r"\sigma", "σ").replace(r"\Delta", "Δ")
    s = s.replace(r"\delta", "δ").replace(r"\mu", "µ").replace(r"\lesssim", "≲")
    s = s.replace(r"\tau", "τ").replace(r"\odot", "☉").replace(r"\propto", "∝")
    s = s.replace(r"\mathrm", "").replace(r"\rm", "").replace(r"\,", "\u2009")
    s = re.sub(r"\\log_\{?10\}?", "log₁₀", s); s = s.replace(r"\log", "log")
    s = re.sub(r"\^\{([^}]*)\}", lambda m: "<sup>" + m.group(1) + "</sup>", s)
    s = re.sub(r"\^(\w)", lambda m: "<sup>" + m.group(1) + "</sup>", s)
    s = re.sub(r"_\{([^}]*)\}", lambda m: "<sub>" + m.group(1) + "</sub>", s)
    s = re.sub(r"_(\w)", lambda m: "<sub>" + m.group(1) + "</sub>", s)
    s = s.replace("{", "").replace("}", "")
    return "<span class='m'>" + s + "</span>"

def cite(keys, style):
    out = []
    for k in [x.strip() for x in keys.split(",")]:
        cited.add(k)
        b = bib.get(k)
        s = b["short"] if b else k
        if style == "t":  # \citet -> Author (Year)
            s = re.sub(r"\s(\d{4})$", r" (\1)", s)
        out.append(s)
    joined = "; ".join(out)
    return "(" + joined + ")" if style == "p" else joined

def inline(s):
    s = re.sub(r"\$(.+?)\$", lambda m: "\x00" + math(m.group(1)) + "\x01", s)
    s = html.escape(s)
    s = s.replace("\x00", "").replace("\x01", "")
    # restore math spans (escape put &lt; on them)
    s = s.replace("&lt;span class='m'&gt;", "<span class='m'>").replace("&lt;/span&gt;", "</span>")
    s = s.replace("&lt;sup&gt;", "<sup>").replace("&lt;/sup&gt;", "</sup>")
    s = s.replace("&lt;sub&gt;", "<sub>").replace("&lt;/sub&gt;", "</sub>")
    s = re.sub(r"\\citep\{([^}]*)\}", lambda m: cite(m.group(1), "p"), s)
    s = re.sub(r"\\citealt\{([^}]*)\}", lambda m: cite(m.group(1), "a"), s)
    s = re.sub(r"\\citealp\{([^}]*)\}", lambda m: cite(m.group(1), "a"), s)
    s = re.sub(r"\\citet\{([^}]*)\}", lambda m: cite(m.group(1), "t"), s)
    s = re.sub(r"\\ref\{([^}]*)\}", lambda m: labelnum.get(m.group(1), "?"), s)
    s = re.sub(r"\\emph\{([^}]*)\}", r"<em>\1</em>", s)
    s = re.sub(r"\\textit\{([^}]*)\}", r"<em>\1</em>", s)
    s = re.sub(r"\\texttt\{([^}]*)\}", r"<code>\1</code>", s)
    s = re.sub(r"\\textsc\{([^}]*)\}", r"<span class='sc'>\1</span>", s)
    s = re.sub(r"\\thanks\{[^}]*\}", "", s)
    s = s.replace(r"\%", "%").replace(r"\&", "&amp;").replace(r"\_", "_").replace(r"\#", "#")
    s = s.replace("---", "—").replace("--", "–").replace(r"\ ", " ").replace("~", "\u00a0")
    s = re.sub(r"``([^']*)''", r"“\1”", s)
    s = re.sub(r"\\cf\.?|\\,", " ", s)
    s = re.sub(r"\\[a-zA-Z]+\*?", "", s)  # drop any stray commands
    return s

# ---------- pass 2: render body ----------
def fig_uri(name):
    p = FIGDIR / name
    if not p.exists(): return None
    return "data:image/png;base64," + base64.b64encode(p.read_bytes()).decode()

parts = []
fignum = 0
pos = 0
struct = re.compile(
    r"\\section\{(?P<s>[^}]*)\}(?:\s*\\label\{[^}]*\})?"
    r"|\\subsection\{(?P<ss>[^}]*)\}(?:\s*\\label\{[^}]*\})?"
    r"|\\paragraph\*?\{(?P<pg>[^}]*)\}"
    r"|\\begin\{figure\*?\}(?P<fig>.*?)\\end\{figure\*?\}", re.S)
sec = 0; sub = 0
def flush(text):
    text = text.strip()
    if not text: return
    for para in re.split(r"\n\s*\n", text):
        para = para.strip()
        if para:
            parts.append("<p>" + inline(para) + "</p>")
for m in struct.finditer(body):
    flush(body[pos:m.start()]); pos = m.end()
    if m.group("s") is not None:
        sec += 1; sub = 0
        parts.append(f"<h2><span class='num'>{sec}</span>{inline(m.group('s'))}</h2>")
    elif m.group("ss") is not None:
        sub += 1
        parts.append(f"<h3><span class='num'>{sec}.{sub}</span>{inline(m.group('ss'))}</h3>")
    elif m.group("pg") is not None:
        parts.append(f"<p class='runin'><strong>{inline(m.group('pg'))}</strong> ")
    elif m.group("fig") is not None:
        fignum += 1
        blk = m.group("fig")
        inc = re.search(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]*)\}", blk)
        cap = re.search(r"\\caption\{(.*?)\}\s*(?:\\label|\\end|$)", blk, re.S)
        name = inc.group(1) if inc else ""
        if not name.endswith(".png"): name += ".png"
        uri = fig_uri(name)
        img = f"<img src='{uri}' alt='Figure {fignum}'>" if uri else f"<div class='missing'>[{name} not found]</div>"
        capt = inline(cap.group(1)) if cap else ""
        parts.append(f"<figure>{img}<figcaption><b>Figure {fignum}.</b> {capt}</figcaption></figure>")
flush(body[pos:])

# ---------- references ----------
def delatex(s):
    s = s.replace(r"\&", "&").replace(r"\delta", "δ")
    s = re.sub(r"\\[Hvcuo'`\"^~=.]\{?([A-Za-z])\}?", r"\1", s)  # accents -> base letter
    s = re.sub(r"\\[a-zA-Z]+", "", s)
    return s.replace("{", "").replace("}", "").strip()
refs = sorted(cited, key=lambda k: bib.get(k, {}).get("sort", (k, "")))
reflis = []
for k in refs:
    b = bib.get(k)
    if not b: reflis.append(f"<li>{k}</li>"); continue
    j = f" <i>{html.escape(delatex(b['journal']))}</i>" + (f" {html.escape(b['vol'])}" if b['vol'] else "") if b['journal'] else ""
    reflis.append(f"<li>{html.escape(delatex(b['full']))} ({b['year']}). {html.escape(delatex(b['title']))}.{j}</li>")

body_html = "\n".join(parts)
abstract_html = inline(re.sub(r"\s+", " ", abstract))
title_html = inline(title)

CSS = """
<style>
:root{--bg:#ffffff;--ink:#1a1a1a;--mut:#5a5f66;--line:#e2e4e8;--accent:#7a2020;
 --serif:"Iowan Old Style",Georgia,"Times New Roman",serif;
 --sans:system-ui,-apple-system,"Segoe UI",sans-serif;--mono:ui-monospace,"Cascadia Code",Consolas,monospace;}
@media (prefers-color-scheme:dark){:root{--bg:#12141a;--ink:#e6e8ec;--mut:#9aa1ab;--line:#2a2e37;--accent:#e0928c;}}
:root[data-theme="light"]{--bg:#ffffff;--ink:#1a1a1a;--mut:#5a5f66;--line:#e2e4e8;--accent:#7a2020;}
:root[data-theme="dark"]{--bg:#12141a;--ink:#e6e8ec;--mut:#9aa1ab;--line:#2a2e37;--accent:#e0928c;}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);font-family:var(--serif);line-height:1.62;}
article{max-width:760px;margin:0 auto;padding:48px 26px 96px;}
h1{font-size:1.72rem;line-height:1.22;font-weight:600;margin:0 0 .5rem;text-wrap:balance;letter-spacing:-.01em;}
.byline{font-family:var(--sans);font-size:.82rem;color:var(--mut);margin:0 0 1.8rem;padding-bottom:1.2rem;border-bottom:1px solid var(--line);}
.abstract{background:color-mix(in srgb,var(--ink) 3%,transparent);border:1px solid var(--line);border-radius:8px;padding:16px 20px;margin:0 0 2.2rem;}
.abstract h4{font-family:var(--sans);text-transform:uppercase;letter-spacing:.08em;font-size:.72rem;color:var(--accent);margin:0 0 .5rem;}
.abstract p{margin:0;font-size:.94rem;}
h2{font-size:1.24rem;font-weight:600;margin:2.4rem 0 .7rem;padding-top:.4rem;text-wrap:balance;}
h3{font-size:1.04rem;font-weight:600;margin:1.7rem 0 .5rem;color:var(--ink);text-wrap:balance;}
.num{font-family:var(--sans);color:var(--accent);font-weight:600;margin-right:.6em;font-size:.9em;}
h3 .num{color:var(--mut);}
p{margin:0 0 .9rem;}
p.runin{margin-top:1.1rem;}
code{font-family:var(--mono);font-size:.86em;background:color-mix(in srgb,var(--ink) 6%,transparent);padding:.05em .35em;border-radius:4px;}
.sc{font-variant:small-caps;letter-spacing:.03em;}
.m{font-family:var(--serif);font-style:italic;white-space:nowrap;}
.m sub,.m sup{font-style:normal;}
figure{margin:1.8rem 0;text-align:center;}
figure img{max-width:100%;height:auto;border:1px solid var(--line);border-radius:6px;background:#fff;}
figcaption{font-family:var(--sans);font-size:.8rem;line-height:1.5;color:var(--mut);text-align:left;margin-top:.6rem;}
figcaption b{color:var(--ink);}
.refh{font-size:1.1rem;margin-top:3rem;border-top:1px solid var(--line);padding-top:1.4rem;}
ol.refs{font-family:var(--sans);font-size:.82rem;line-height:1.5;color:var(--mut);padding-left:1.6rem;}
ol.refs li{margin:.35rem 0;}
ol.refs i{color:var(--ink);}
sup,sub{font-size:.72em;line-height:0;}
</style>
"""
OUT.write_text(CSS + f"""
<article>
<h1>{title_html}</h1>
<p class="byline">APEX validation manuscript · working draft (pre-submission, un-typeset) · figures embedded</p>
<div class="abstract"><h4>Abstract</h4><p>{abstract_html}</p></div>
{body_html}
<h2 class="refh">References</h2><ol class="refs">{''.join(reflis)}</ol>
</article>""", encoding="utf-8")
print("figures:", fignum, "| cited:", len(cited), "| sections:", sec, "| out:", OUT)
