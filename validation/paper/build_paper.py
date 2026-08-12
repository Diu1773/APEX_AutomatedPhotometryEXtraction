#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Build the paper review artifact from the canonical Korean manuscript.

The canonical source remains ``논문작업/MANUSCRIPT_ko.md``.  This command creates
three reproducible outputs:

* ``latex/review.tex`` and ``latex/generated/*.tex`` — the real TeX source used
  when XeLaTeX/LuaLaTeX (or another configured engine) is available;
* ``build/review/APEX_review.pdf`` — the page-faithful review PDF; and
* ``build/review/build_status.json`` — the engine, source hash and page count.

The bundled environment does not necessarily contain a TeX distribution.  In
that case PyMuPDF's Story engine is used as a deterministic review fallback;
the command never labels that fallback as an A&A submission PDF.  Installing
TeX is deliberately left to the user because it is a machine-level dependency.
"""

from __future__ import annotations

import argparse
import base64
import hashlib
import html
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

try:
    import fitz  # PyMuPDF, bundled in .venv-deploy
except ImportError:  # pragma: no cover - the CLI reports a useful error below
    fitz = None


# ``validation/`` is a junction in this checkout.  Keep generated paths on the
# user-facing worktree instead of silently switching to the junction target.
ROOT = Path(__file__).absolute().parent
SOURCE = ROOT / "논문작업" / "MANUSCRIPT_ko.md"
PREVIEW = ROOT / "MANUSCRIPT_ko_preview.html"
LATEX = ROOT / "latex"
GENERATED = LATEX / "generated"
DEFAULT_OUT = ROOT / "build" / "review"

# Keep this mapping in sync with render_preview.py.  The preview source already
# contains the base64 image, so the PDF fallback does not need to read data*/.
FIG_FILES = {
    1: "fig_architecture.png",
    2: "fig_calibration_step0.png",
    3: "fig11_detector.png",
    4: "fig12_preproc_crosscheck.png",
    5: "fig13_cross_instrument.png",
    6: "fig6_qc_validation.png",
    7: "fig_detection_threshold.png",
    8: "fig_completeness_realvssynth.png",
    9: "fig_wcs_engines.png",
    10: "fig2_error_model.png",
    11: "fig3_parameter_sweep.png",
    12: "fig_photometry_crosschecks.png",
    13: "fig_psf_validation.png",
    14: "fig9_crowded_field.png",
    15: "fig_external_validation.png",
    16: "fig_timeseries_validation.png",
    17: "fig_lc_yzboo.png",
}


def _strip_tags(value: str) -> str:
    value = re.sub(r"<br\s*/?>", " ", value, flags=re.I)
    value = re.sub(r"<[^>]+>", "", value)
    return re.sub(r"\s+", " ", html.unescape(value)).strip()


def _inner(html_doc: str, element_id: str) -> str:
    match = re.search(
        rf'<div[^>]*(?:\bid=["\']{re.escape(element_id)}["\']|\bclass=["\'][^"\']*\b{re.escape(element_id)}\b[^"\']*["\'])[^>]*>',
        html_doc,
        re.I,
    )
    if not match:
        raise RuntimeError(f"preview fragment not found: {element_id}")
    start = match.end()
    depth = 1
    tag_re = re.compile(r"</?div\b[^>]*>", re.I)
    for tag in tag_re.finditer(html_doc, start):
        raw = tag.group(0)
        if raw.startswith("</"):
            depth -= 1
            if depth == 0:
                return html_doc[start:tag.start()]
        elif not raw.rstrip().endswith("/>"):
            depth += 1
    raise RuntimeError(f"unterminated preview fragment: {element_id}")


def _source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def refresh_preview() -> None:
    """Regenerate the canonical HTML source before extracting its fragments."""

    subprocess.run(
        [sys.executable, "-X", "utf8", str(ROOT / "render_preview.py")],
        cwd=ROOT,
        check=True,
    )


def read_fragments(preview_path: Path) -> dict[str, str]:
    doc = preview_path.read_text(encoding="utf-8")
    src = _inner(doc, "src")
    return {
        "src": src,
        "title": _inner(src, "titleblock"),
        "abstract": _inner(src, "absblock"),
        "toc": _inner(src, "tocblock"),
        "flow": _inner(src, "flow"),
    }


def _protect(value: str, pattern: str, transform, store: list[str]) -> str:
    def repl(match: re.Match[str]) -> str:
        store.append(transform(match.group(0)))
        return f"\x00{len(store) - 1}\x00"

    return re.sub(pattern, repl, value, flags=re.S)


def tex_inline(value: str) -> str:
    """Convert the Markdown inline subset used by the manuscript to TeX."""

    tokens: list[str] = []
    value = value.replace("\n", " ")
    value = _protect(value, r"`[^`]*`", lambda m: r"\texttt{\detokenize{" + m[1:-1] + "}}", tokens)
    value = _protect(value, r"\$[^$]+\$", lambda m: m, tokens)
    value = _protect(
        value,
        r"\\(?:citep|citet|citealt|citealp|ref)\{[^{}]*\}",
        lambda m: m,
        tokens,
    )
    value = _protect(
        value,
        r"\\(?:textbf|textit|emph|texttt|mathrm|mathbf|mathit)\{[^{}]*\}",
        lambda m: m,
        tokens,
    )
    value = re.sub(r"\*\*([^*]+)\*\*", r"\\textbf{\1}", value)
    value = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"\\emph{\1}", value)
    value = value.replace("&", r"\&").replace("%", r"\%").replace("#", r"\#")
    value = value.replace("_", r"\_")
    value = value.replace("~", r"\textasciitilde{}")
    value = re.sub(r"\x00(\d+)\x00", lambda m: tokens[int(m.group(1))], value)
    return value


def markdown_blocks(source: str) -> tuple[str, str, list[tuple[str, str]]]:
    """Return title, abstract and body blocks from canonical Markdown.

    A block is ``(kind, payload)`` where kind is heading, paragraph, table,
    figure-ref, or raw.  Tables and figures are kept as separate blocks so the
    fallback PDF can give them a full-width review sheet and never split a
    caption across columns.
    """

    lines = source.splitlines()
    title = ""
    abstract: list[str] = []
    body: list[str] = []
    mode = "before"
    for line in lines:
        if line.startswith("# ") and not title:
            title = line[2:].strip()
            continue
        if line.startswith("<!--"):
            continue
        if line.startswith("*영문 제출본"):
            continue
        if re.match(r"^##\s+초록\s*$", line):
            mode = "abstract"
            continue
        if mode == "abstract" and re.match(r"^##\s+", line):
            mode = "body"
        if mode == "abstract":
            abstract.append(line)
        elif mode == "body":
            body.append(line)
    blocks: list[tuple[str, str]] = []
    buf: list[str] = []

    def flush() -> None:
        if buf and " ".join(x.strip() for x in buf).strip():
            blocks.append(("paragraph", " ".join(x.strip() for x in buf).strip()))
        buf.clear()

    i = 0
    while i < len(body):
        line = body[i]
        if re.match(r"^#{2,4}\s+", line):
            flush()
            level = len(line) - len(line.lstrip("#"))
            blocks.append((f"h{level - 1}", line[level:].strip()))
            i += 1
            continue
        if line.strip().startswith("|") and i + 1 < len(body) and re.match(r"^\s*\|?\s*:?-{2,}", body[i + 1]):
            flush()
            rows = [line]
            i += 1
            while i < len(body) and body[i].strip().startswith("|"):
                rows.append(body[i])
                i += 1
            blocks.append(("table", "\n".join(rows)))
            continue
        if not line.strip():
            flush()
        else:
            buf.append(line)
        i += 1
    flush()
    abstract_text = "\n".join(x for x in abstract if x.strip() and not x.startswith("---"))
    return title, abstract_text, blocks


def _table_tex(raw: str) -> str:
    rows = []
    for line in raw.splitlines():
        cells = [x.strip() for x in line.strip().strip("|").split("|")]
        if all(re.fullmatch(r":?-{2,}:?", x.replace(" ", "")) for x in cells):
            continue
        rows.append(cells)
    if not rows:
        return ""
    n = max(len(x) for x in rows)
    cols = "X" * n
    out = [r"\begin{table*}[!t]", r"\centering\small", rf"\begin{{tabularx}}{{\textwidth}}{{{cols}}}", r"\toprule"]
    for ri, row in enumerate(rows):
        row = row + [""] * (n - len(row))
        out.append(" & ".join(tex_inline(x) for x in row) + r" \\")
        if ri == 0:
            out.append(r"\midrule")
    out += [r"\bottomrule", r"\end{tabularx}", r"\end{table*}"]
    return "\n".join(out)


def _figure_tex(number: int, caption: str) -> str:
    name = FIG_FILES.get(number)
    if not name:
        return ""
    cap = tex_inline(caption)
    return "\n".join(
        [
            r"\begin{figure*}[!t]",
            r"\centering",
            rf"\includegraphics[width=\textwidth]{{../figures/{name}}}",
            rf"\caption[{tex_inline(f'그림 {number}')}]{{{cap}}}",
            rf"\label{{fig:{number}}}",
            r"\end{figure*}",
        ]
    )


def generate_latex(source: str, fragments: dict[str, str], captions: dict[int, str]) -> Path:
    title, abstract, blocks = markdown_blocks(source)
    GENERATED.mkdir(parents=True, exist_ok=True)
    body_lines: list[str] = []
    seen_figs: set[int] = set()
    for kind, payload in blocks:
        if kind.startswith("h"):
            level = int(kind[1:])
            command = {1: "section", 2: "subsection", 3: "subsubsection"}.get(level, "paragraph")
            body_lines.append(rf"\{command}{{{tex_inline(payload)}}}")
            continue
        if kind == "table":
            body_lines.append(_table_tex(payload))
            continue
        body_lines.append(tex_inline(payload))
        body_lines.append("")
        for match in re.finditer(r"그림\s*(\d+)", payload):
            number = int(match.group(1))
            if number not in seen_figs and number in FIG_FILES:
                body_lines.append(_figure_tex(number, captions.get(number, f"그림 {number}.")))
                seen_figs.add(number)
    (GENERATED / "abstract.tex").write_text(
        "\n".join(rf"\paragraph*{{}}{tex_inline(x.strip())}" for x in abstract.splitlines() if x.strip() and not x.startswith("**핵심어"))
        + "\n\\medskip\n"
        + tex_inline(next((x.strip() for x in abstract.splitlines() if x.strip().startswith("**핵심어")), "")),
        encoding="utf-8",
    )
    (GENERATED / "body.tex").write_text("\n\n".join(body_lines), encoding="utf-8")
    tex = rf"""% Generated by build_paper.py. Do not edit generated/*.tex by hand.
\documentclass[10pt,a4paper]{{article}}
\usepackage{{iftex}}
\ifXeTeX
  \usepackage{{fontspec}}
  \usepackage{{xeCJK}}
  \setmainfont{{TeX Gyre Termes}}
  \setCJKmainfont{{Malgun Gothic}}
\else\ifLuaTeX
  \usepackage{{fontspec}}
  \usepackage{{luatexja-fontspec}}
  \setmainfont{{TeX Gyre Termes}}
  \setmainjfont{{Malgun Gothic}}
\else
  \PackageError{{APEX}}{{Korean review requires XeLaTeX or LuaLaTeX}}{{Use the PyMuPDF fallback or install a Unicode TeX engine.}}
\fi\fi
\usepackage[inner=18mm,outer=18mm,top=17mm,bottom=18mm,headheight=14pt]{{geometry}}
\usepackage{{graphicx,booktabs,tabularx,array,amsmath,natbib,multicol,fancyhdr,lastpage,hyperref,placeins}}
\graphicspath{{{{../figures/}}}}
\hypersetup{{colorlinks=true,linkcolor=blue,citecolor=blue,urlcolor=blue}}
\setlength{{\columnsep}}{{6mm}}
\setlength{{\parindent}}{{1em}}
\setlength{{\parskip}}{{0pt}}
\pagestyle{{fancy}}
\fancyhf{{}}
\fancyhead[C]{{APEX — raw 영상부터 CMD·광도 곡선까지의 통합 측광 파이프라인}}
\fancyfoot[C]{{APEX · \thepage\ / \pageref{{LastPage}}}}
\renewcommand{{\headrulewidth}}{{0.3pt}}
\renewcommand{{\footrulewidth}}{{0pt}}
\title{{{tex_inline(title)}}}
\author{{저자 미기재 (투고 전 확정)}}
\date{{국문 검토본 · A\&A Section 15 목표}}
\begin{{document}}
\pagenumbering{{gobble}}
\maketitle
\thispagestyle{{empty}}
\clearpage
\renewcommand{{\contentsname}}{{차례}}
\renewcommand{{\listfigurename}}{{그림 목록}}
\tableofcontents
\listoffigures
\listoftables
\clearpage
\pagenumbering{{arabic}}
\setcounter{{page}}{{3}}
\begin{{abstract}}
\input{{generated/abstract.tex}}
\end{{abstract}}
\begin{{multicols}}{{2}}
\input{{generated/body.tex}}
\end{{multicols}}
\bibliographystyle{{plainnat}}
\bibliography{{../논문작업/references}}
\end{{document}}
"""
    tex_path = LATEX / "review.tex"
    tex_path.write_text(tex, encoding="utf-8")
    return tex_path


PDF_CSS = r"""
@page { size: A4; }
* { box-sizing: border-box; }
body { margin: 0; color: #111; font-family: "Malgun Gothic", "Noto Serif CJK KR", serif;
       font-size: 9.2pt; line-height: 1.27; }
h1 { font-size: 22pt; line-height: 1.22; text-align: center; margin: 0 0 14pt; }
h2 { font-size: 13pt; line-height: 1.2; margin: 9pt 0 4pt; break-after: avoid; }
h3 { font-size: 10.6pt; line-height: 1.2; margin: 7pt 0 3pt; break-after: avoid; }
h4 { font-size: 10pt; margin: 5pt 0 2pt; }
p { margin: 0 0 5pt; text-align: justify; }
code { font-family: "Courier New", monospace; font-size: 0.88em; }
.math { font-family: "Times New Roman", serif; font-style: italic; }
.titleblock { text-align: center; margin: 25pt 0 10pt; }
.jrnl { font-size: 9pt; text-align: left; margin-bottom: 28pt; }
.jr { float: right; font-weight: bold; font-size: 14pt; }
.pauth, .paff, .pdate { text-align: center; margin: 4pt 0; }
.pdate { font-size: 8pt; color: #555; }
.absblock { margin: 0 18pt 8pt; }
.abshead { text-align: center; font-weight: bold; letter-spacing: 1pt; margin-bottom: 4pt; }
.absblock p { text-align: justify; }
.kw { font-size: 8.5pt; }
.flow { margin-top: 3pt; }
figure { margin: 4pt 0 8pt; text-align: center; }
figure img { max-width: 100%; height: auto; }
figcaption { font-size: 8pt; line-height: 1.18; text-align: justify; margin-top: 3pt; }
.tw { margin: 4pt 0 8pt; }
table { border-collapse: collapse; width: 100%; font-size: 8pt; line-height: 1.17; }
th { border-top: 0.7pt solid #111; border-bottom: 0.5pt solid #111; text-align: left; padding: 2pt 3pt; }
td { padding: 1.5pt 3pt; vertical-align: top; }
tr:last-child td { border-bottom: 0.7pt solid #111; }
.toc-title { font-size: 14pt; font-weight: bold; margin-bottom: 6pt; border-bottom: .7pt solid #111; }
.toc-table { font-size: 7.2pt; line-height: 1.15; }
.toc-table td { padding: 1.1pt 2pt; border: 0; }
.toc-table td:last-child { width: 28pt; text-align: right; }
.cover-note { font-size: 8pt; color: #555; text-align: center; margin-top: 20pt; }
"""


def _story(html_fragment: str, extra_css: str = ""):
    if fitz is None:
        raise RuntimeError("PyMuPDF is required for the review fallback")
    return fitz.Story(html=html_fragment, user_css=PDF_CSS + extra_css)


MEDIA = fitz.Rect(0, 0, 595.276, 841.89) if fitz else None
MARGIN = 42.0


def render_full_page(html_fragment: str, path: Path) -> None:
    writer = fitz.DocumentWriter(str(path))
    story = _story(html_fragment)
    def rectfn(_page: int, _filled):
        return MEDIA, fitz.Rect(MARGIN, MARGIN, MEDIA.width - MARGIN, MEDIA.height - MARGIN), fitz.Identity
    story.write(writer, rectfn)
    writer.close()


def split_flow(flow: str) -> list[tuple[str, str]]:
    pattern = re.compile(r"(<h[2-4][^>]*>.*?</h[2-4]>|<p[^>]*>.*?</p>|<figure[^>]*>.*?</figure>|<div class=\"tw\">.*?</div>)", re.S | re.I)
    blocks = [("text", x) for x in pattern.findall(flow)]
    # The preview flow is block-only; retain any unparsed text as prose.
    if not blocks:
        return [("text", flow)]
    out: list[tuple[str, str]] = []
    pos = 0
    for match in pattern.finditer(flow):
        if flow[pos:match.start()].strip():
            out.append(("text", flow[pos:match.start()]))
        raw = match.group(1)
        normalized = raw.lstrip().lower()
        kind = "figure" if normalized.startswith("<figure") else "table" if "class=\"tw\"" in normalized else "text"
        out.append((kind, raw))
        pos = match.end()
    if flow[pos:].strip():
        out.append(("text", flow[pos:]))
    return out


def render_body_pdf(abstract_html: str, flow_html: str, path: Path) -> None:
    """Render body with an abstract + two columns and conservative floats.

    A float is first offered the unused lower part of the current page when the
    preceding prose is short enough.  If it cannot fit, it is promoted to its
    own page.  This keeps the fallback readable without creating the large
    blank pages that a strict "one float = one page" policy produces.
    """

    writer = fitz.DocumentWriter(str(path))
    chunks = split_flow(flow_html)
    text_chunks: list[str] = []
    # We put prose into a Story until a figure/table boundary, then give the
    # float a full-width page. This is intentionally conservative for review:
    # no caption can be clipped or force an empty column.
    for kind, raw in chunks:
        if kind == "text":
            text_chunks.append(raw)
            continue
        if text_chunks:
            text_raw = "".join(text_chunks)
            first_abstract = abstract_html if path.stat().st_size == 0 else ""
            if not first_abstract and len(text_raw) <= 14000:
                _draw_text_with_float(writer, text_raw, raw)
            else:
                _draw_text_story(writer, text_raw, include_abstract=first_abstract)
                if raw:
                    _draw_full_story(writer, raw)
            abstract_html = ""
            text_chunks.clear()
        else:
            _draw_full_story(writer, raw)
    if text_chunks:
        _draw_text_story(writer, "".join(text_chunks), include_abstract=abstract_html)
    writer.close()


def _draw_full_story(writer, raw: str) -> None:
    story = _story(raw)
    more = True
    while more:
        dev = writer.begin_page(MEDIA)
        more, _filled = story.place(fitz.Rect(MARGIN, MARGIN, MEDIA.width - MARGIN, MEDIA.height - MARGIN))
        story.draw(dev)
        writer.end_page()


def _draw_text_with_float(writer, text_raw: str, float_raw: str) -> None:
    """Put a short prose run in two columns and try its float below it."""

    story = _story(text_raw)
    float_story = _story(float_raw)
    more_text = True
    first = True
    more_float = True
    while more_text or (first and more_float):
        dev = writer.begin_page(MEDIA)
        filled = []
        columns = [
            fitz.Rect(MARGIN, MARGIN, MEDIA.width / 2 - 8, MEDIA.height - MARGIN),
            fitz.Rect(MEDIA.width / 2 + 8, MARGIN, MEDIA.width - MARGIN, MEDIA.height - MARGIN),
        ]
        for rect in columns:
            if not more_text:
                break
            more_text, used = story.place(rect)
            story.draw(dev)
            filled.append(fitz.Rect(used))
        if first and not more_text and filled:
            top = max(rect.y1 for rect in filled) + 10
            bottom = MEDIA.height - MARGIN
            if top + 40 < bottom:
                more_float, _ = float_story.place(fitz.Rect(MARGIN, top, MEDIA.width - MARGIN, bottom))
                float_story.draw(dev)
        writer.end_page()
        first = False
        if not more_text:
            break
    if more_float:
        _draw_full_story(writer, float_raw)


def _draw_text_story(writer, raw: str, include_abstract: str = "") -> None:
    story = _story(raw)
    abs_story = _story(include_abstract) if include_abstract else None
    more_text = True
    more_abs = bool(abs_story)
    first = True
    while more_text or more_abs:
        dev = writer.begin_page(MEDIA)
        drew_abstract = False
        if more_abs and abs_story:
            # Reserve the abstract rectangle first, but draw it after the body
            # so that Story's device clip cannot hide the two-column text.
            more_abs, _ = abs_story.place(fitz.Rect(MARGIN, MARGIN, MEDIA.width - MARGIN, 430))
            drew_abstract = True
        if not more_abs:
            top = 450 if first and include_abstract else MARGIN
            columns = [
                fitz.Rect(MARGIN, top, MEDIA.width / 2 - 8, MEDIA.height - MARGIN),
                fitz.Rect(MEDIA.width / 2 + 8, top, MEDIA.width - MARGIN, MEDIA.height - MARGIN),
            ]
            for rect in columns:
                if not more_text:
                    break
                more_text, _ = story.place(rect)
                story.draw(dev)
        if drew_abstract and abs_story:
            # Draw body first.  Story.draw leaves the device clipped to its
            # placement rectangle; drawing the full-width abstract first would
            # clip the following two columns in some PyMuPDF versions.
            abs_story.draw(dev)
        writer.end_page()
        first = False
        if not more_abs and not more_text:
            break


def _body_page_map(path: Path, headings: list[str], figures: list[str]) -> dict[str, int]:
    doc = fitz.open(path)
    result: dict[str, int] = {}
    for key in headings + figures:
        needle = re.sub(r"\s+", " ", key).strip()
        for i, page in enumerate(doc):
            text = re.sub(r"\s+", " ", page.get_text("text")).strip()
            if needle and (needle in text or needle[:24] in text):
                result[key] = i + 1
                break
    doc.close()
    return result


def _heading_texts(flow: str) -> list[str]:
    return [_strip_tags(x) for x in re.findall(r"<h[23][^>]*>.*?</h[23]>", flow, re.S | re.I)]


def _figure_labels(flow: str) -> list[str]:
    labels = []
    for caption in re.findall(r"<figcaption[^>]*>(.*?)</figcaption>", flow, re.S | re.I):
        text = _strip_tags(caption)
        match = re.search(r"그림\s*\d+", text)
        if match:
            labels.append(match.group(0))
    return labels


def make_toc(title: str, headings: list[str], figures: list[str], page_map: dict[str, int], body_offset: int) -> str:
    rows = []
    for item in headings:
        page = page_map.get(item)
        rows.append(f"<tr><td>{html.escape(item)}</td><td>{'' if page is None else page + body_offset}</td></tr>")
    figure_rows = []
    for item in figures:
        page = page_map.get(item)
        figure_rows.append(f"<tr><td>{html.escape(item)}</td><td>{'' if page is None else page + body_offset}</td></tr>")
    return (
        f"<h1 class='toc-title'>차례 · {html.escape(title)}</h1>"
        "<table class='toc-table'>" + "".join(rows) + "</table>"
        "<h2 class='toc-title'>그림</h2><table class='toc-table'>" + "".join(figure_rows) + "</table>"
    )


def add_footers(path: Path) -> int:
    doc = fitz.open(path)
    total = doc.page_count
    for i, page in enumerate(doc):
        page.insert_textbox(
            fitz.Rect(36, MEDIA.height - 30, MEDIA.width - 36, MEDIA.height - 14),
            f"APEX · {i + 1} / {total}",
            fontname="helv",
            fontsize=8,
            color=(0.18, 0.18, 0.18),
            align=1,
        )
    tmp = path.with_suffix(".numbered.pdf")
    doc.save(tmp, garbage=4, deflate=True)
    doc.close()
    tmp.replace(path)
    return total


def merge_pdfs(output: Path, *parts: Path) -> None:
    out = fitz.open()
    for part in parts:
        with fitz.open(part) as src:
            out.insert_pdf(src)
    out.save(output, garbage=4, deflate=True)
    out.close()


def build_review_pdf(fragments: dict[str, str], title: str, output: Path) -> int:
    output.parent.mkdir(parents=True, exist_ok=True)
    tmp = output.parent / "_parts"
    tmp.mkdir(exist_ok=True)
    body = tmp / "body.pdf"
    cover = tmp / "cover.pdf"
    toc = tmp / "toc.pdf"
    cover_html = fragments["title"] + "<p class='cover-note'>국문 검토본 · A&A Section 15 목표<br>정본: 논문작업/MANUSCRIPT_ko.md<br>이 PDF는 build_paper.py가 생성한 페이지 검토본이다.</p>"
    render_full_page(cover_html, cover)
    render_body_pdf(fragments["abstract"], fragments["flow"], body)
    headings = _heading_texts(fragments["flow"])
    figures = _figure_labels(fragments["flow"])
    page_map = _body_page_map(body, headings, figures)
    toc_html = make_toc(title, headings, figures, page_map, body_offset=2)
    render_full_page(toc_html, toc)
    merge_pdfs(output, cover, toc, body)
    pages = add_footers(output)
    # The three intermediate PDFs can be hundreds of megabytes because the
    # preview embeds raster figures.  Keep only the final review artifact.
    shutil.rmtree(tmp, ignore_errors=True)
    return pages


def try_tex(tex_path: Path, output: Path, engine: str | None) -> tuple[str, str | None]:
    candidates = [engine] if engine else ["xelatex", "lualatex", "tectonic", "pdflatex"]
    selected = next((x for x in candidates if x and shutil.which(x)), None)
    if not selected:
        return "pymupdf-fallback", "No TeX engine found on PATH"
    outdir = output.parent / "tex"
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [selected]
    if selected != "tectonic":
        cmd += ["-interaction=nonstopmode", "-halt-on-error", f"-output-directory={outdir}"]
    else:
        cmd += ["--outdir", str(outdir)]
    cmd += [str(tex_path)]
    try:
        for _ in range(2):
            subprocess.run(cmd, cwd=LATEX, check=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        produced = outdir / (tex_path.stem + ".pdf")
        if not produced.exists():
            return "pymupdf-fallback", f"{selected} did not create {produced.name}"
        shutil.copy2(produced, output)
        return selected, None
    except (OSError, subprocess.CalledProcessError) as exc:
        return "pymupdf-fallback", f"{selected} failed: {exc}"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--preview", type=Path, default=PREVIEW)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--engine", help="xelatex, lualatex, tectonic or pdflatex")
    parser.add_argument("--no-refresh-preview", action="store_true")
    args = parser.parse_args(argv)
    source = args.source if args.source.is_absolute() else ROOT / args.source
    preview = args.preview if args.preview.is_absolute() else ROOT / args.preview
    if not args.no_refresh_preview:
        refresh_preview()
    fragments = read_fragments(preview)
    title, abstract, _blocks = markdown_blocks(source.read_text(encoding="utf-8"))
    captions = {
        i: _strip_tags(c)
        for i, c in enumerate(re.findall(r"<figcaption[^>]*>(.*?)</figcaption>", fragments["flow"], re.S | re.I), 1)
    }
    tex_path = generate_latex(source.read_text(encoding="utf-8"), fragments, captions)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = args.out_dir / "APEX_review.pdf"
    engine, reason = try_tex(tex_path, pdf_path, args.engine)
    if engine == "pymupdf-fallback":
        pages = build_review_pdf(fragments, title, pdf_path)
    else:
        pages = fitz.open(pdf_path).page_count if fitz else None
    status = {
        "source": str(source),
        "source_sha256": _source_hash(source),
        "engine": engine,
        "reason": reason,
        "pdf": str(pdf_path),
        "pages": pages,
        "latex": str(tex_path),
        "layout": {"paper": "A4", "cover": 1, "toc": 2, "abstract_body": 3, "body_columns": 2},
    }
    (args.out_dir / "build_status.json").write_text(json.dumps(status, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(status, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
