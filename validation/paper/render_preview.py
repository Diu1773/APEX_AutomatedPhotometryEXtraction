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
OUT_ARTIFACT = ROOT / "MANUSCRIPT_ko_artifact.html"   # 아티팩트 게시용(doctype 없음)

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

# 절(<h2>) 마다 새 페이지
_parts = [p for p in re.split(r'(?=<h2)', BODY) if p.strip()]
_sheets = "\n".join(f'<section class="sheet">{p}</section>' for p in _parts)

COVER = f'''<section class="sheet cover">
  <div class="cover-inner">
    <p class="cover-kicker">RAS Techniques and Instruments · 투고 준비 원고</p>
    <h1 class="cover-title">{html.escape(DOC_TITLE)}</h1>
    <p class="cover-authors">저자 미기재 (투고 전 확정)</p>
    <p class="cover-affil">한국천문연구원 · 소형망원경 측광 파이프라인</p>
    <div class="cover-meta">
      <span>국문 원고</span><span>그림 {len(FIGMAP)}점</span><span>인용 {len(LAB)}건</span>
    </div>
    <p class="cover-note">{DOC_NOTE}</p>
  </div>
</section>'''

TOC = f'''<section class="sheet toc">
  <h2 class="toc-h">목차</h2>
  <ol class="toc-list">{_toc_html}</ol>
  <h2 class="toc-h">그림 목록</h2>
  <ol class="toc-figs">{_fig_html}</ol>
</section>'''

BODY = COVER + TOC + _sheets


CSS = r"""
:root{ color-scheme:light; --ink:#111; --muted:#4a4f56; --code:#eef0f2;
       --page-w:46rem; --page-h:1040px; }
*{box-sizing:border-box}
html,body{margin:0;background:#e6e7ea;}
/* 한글 본문 폰트: Noto Serif KR(본명조 계열)을 1순위로 둔다. 예전 스택은 미설치
   Nanum Myeongjo 를 먼저 찾다가 Batang 으로 떨어져 화면에서 구식으로 보였다. */
/* 숫자가 많은 문서라 lining figure 가 필수다. Georgia 는 old-style figure 가 기본이라
   "Step 0" 이 "Step o" 로 보인다 — Cambria/Times 를 앞에 둔다. */
body{color:var(--ink);
  font-family:"Cambria","Times New Roman",Times,"Noto Serif KR",
              "Apple SD Gothic Neo","Malgun Gothic",serif;
  font-size:15px;line-height:1.70;text-rendering:optimizeLegibility;
  font-variant-numeric:lining-nums;font-feature-settings:"lnum" 1;
  word-break:keep-all;overflow-wrap:break-word;
  -webkit-font-smoothing:antialiased;}

/* ── 위계: 본문(15px) < h4 < h3(17.4px) < h2(24px) < 표지 제목 ── */
h2{font-size:1.60rem;line-height:1.28;font-weight:700;margin:0 0 1rem;
   letter-spacing:-.012em;text-wrap:balance;}
h3{font-size:1.16rem;line-height:1.42;font-weight:700;margin:1.55rem 0 .42rem;
   text-wrap:balance;}
h4{font-size:1.00rem;font-weight:700;font-style:italic;margin:1.1rem 0 .28rem;}

/* 한글은 양끝맞춤하면 어절 간격이 벌어져 읽기 나빠진다. 왼쪽 정렬 + 문단 간격. */
p{margin:0 0 .82rem;text-align:left;}
p.pending{color:var(--muted);font-style:italic;}
strong{font-weight:700;}
code{font-family:"Courier New",ui-monospace,monospace;background:var(--code);
  padding:.02em .28em;font-size:.88em;}
.math{font-family:"Cambria Math","Times New Roman",serif;white-space:nowrap;}
.math sub,.math sup{font-size:.74em;}

/* 초록: 제목은 작고 자간 넓게, 본문은 괘선 사이 좁은 단 (논문 관례) */
h2#abstract{font-size:.94rem;letter-spacing:.22em;text-align:center;font-weight:700;
  margin:0 0 1.1rem;}
h2#abstract + p{font-size:.94rem;line-height:1.66;margin:0 .6rem .9rem;
  padding:1.1rem .2rem;border-top:1.4px solid #111;border-bottom:1.4px solid #111;}

blockquote{margin:.9rem 1.4rem;padding-left:.9rem;border-left:2px solid #999;
  color:var(--muted);font-size:.95rem;}

.tw{overflow-x:auto;margin:1rem 0;}
table{border-collapse:collapse;width:100%;font-size:.79rem;
  font-variant-numeric:tabular-nums;line-height:1.38;margin:0 auto;}
thead th{border-top:1.3px solid #000;border-bottom:.8px solid #000;
  text-align:left;padding:.4rem .55rem;font-weight:700;vertical-align:bottom;}
tbody td{padding:.32rem .55rem;vertical-align:top;}
tbody tr:last-child td{border-bottom:1.3px solid #000;}

figure{margin:1.3rem 0;text-align:center;}
figure img{max-width:100%;height:auto;background:#fff;}
figcaption{font-size:.79rem;color:var(--ink);margin:.5rem auto 0;max-width:38rem;
  text-align:left;line-height:1.45;}
figcaption b{font-weight:700;}
::selection{background:#ccd8ee;}
a{color:#7a1010;text-decoration:none;}

/* ── 지면(page) — JS 조판이 여기에 내용을 흘려 넣는다 ── */
#src{display:none;}
#book{padding:.75rem 0 .75rem;}
/* 조판 중에는 모든 지면을 펼쳐 둔다 — display:none 이면 높이가 0 이라 측정이 안 된다 */
#book.measuring .page{display:block;visibility:hidden;position:absolute;top:0;left:-99999px;}
.page{position:relative;background:#fff;width:min(var(--page-w),94vw);
  height:var(--page-h);margin:0 auto;padding:3.1rem 3rem 2.7rem;
  box-shadow:0 1px 5px rgba(0,0,0,.17);border-radius:2px;display:none;}
.page.on{display:block;}
.page-body{height:100%;overflow:hidden;}
.folio{position:absolute;left:0;right:0;bottom:1.05rem;text-align:center;
  font-size:.72rem;color:#9299a1;letter-spacing:.1em;font-variant-numeric:tabular-nums;}
.p-cover .page-body{display:flex;align-items:center;}
@media (max-width:680px){ .page{padding:1.7rem 1.3rem 2.3rem;} }

/* ── 표지 ── */
/* p{text-align:left} 가 상속을 이기므로 표지 문단은 따로 가운데로 되돌린다 */
.cover-inner{width:100%;text-align:center;}
.cover-inner p{text-align:center;}
.cover-kicker{font-size:.72rem;letter-spacing:.16em;color:#5a5f66;margin:0 0 2.2rem;
  text-transform:uppercase;}
.cover-title{font-size:1.92rem;line-height:1.36;font-weight:700;margin:0 0 2rem;
  letter-spacing:-.015em;text-wrap:balance;}
.cover-authors{font-size:1rem;margin:0 0 .3rem;}
.cover-affil{font-size:.86rem;color:#5a5f66;margin:0 0 2.4rem;}
.cover-meta{display:flex;gap:1.3rem;justify-content:center;font-size:.76rem;
  color:#5a5f66;border-top:1px solid #dfe1e6;border-bottom:1px solid #dfe1e6;
  padding:.66rem 0;margin:0 auto 1.6rem;max-width:26rem;}
.cover-note{font-size:.75rem;color:#7a7f86;font-style:italic;margin:0;}

/* ── 목차 ── */
.toc-h{font-size:1.16rem;margin:0 0 .7rem;padding-bottom:.35rem;
  border-bottom:1px solid #dfe1e6;font-weight:700;}
.toc-list,.toc-figs{list-style:none;padding:0;margin:0 0 1.6rem;font-size:.88rem;}
.toc-list li.h2{font-weight:700;margin:.62rem 0 .16rem;}
.toc-list li.h3{padding-left:1.3rem;color:#3a3f46;margin:.1rem 0;}
.toc-figs li{margin:.16rem 0;color:#3a3f46;font-size:.81rem;}
/* 목차 항목: float 로 쪽번호를 붙이면 두 줄짜리 항목에서 아래 항목 번호와 겹친다.
   flex 로 본문/쪽번호를 나눈다. */
.go{cursor:pointer;display:flex;align-items:baseline;gap:.6rem;}
.go .tx{flex:1 1 auto;min-width:0;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;}
.go .pn{flex:0 0 auto;color:#9299a1;font-size:.82em;font-variant-numeric:tabular-nums;}
.go:hover .tx{color:#7a1010;text-decoration:underline;}

/* ── 넘김 막대 ── */
#nav{position:fixed;left:0;right:0;bottom:0;z-index:9;display:flex;height:48px;
  align-items:center;justify-content:center;gap:.5rem;padding:0 .5rem;
  background:rgba(243,243,246,.97);border-top:1px solid #d4d7dd;}
#nav button{font:inherit;font-size:.84rem;line-height:1;padding:.5rem .85rem;
  border:1px solid #c2c6cd;background:#fff;color:#20242a;border-radius:3px;
  cursor:pointer;-webkit-tap-highlight-color:transparent;}
#nav button:disabled{opacity:.34;cursor:default;}
#pgind{font-size:.82rem;color:#3a3f46;min-width:6.2rem;text-align:center;
  font-variant-numeric:tabular-nums;}
#loading{max-width:24rem;margin:22vh auto;text-align:center;color:#6b7078;
  font-size:.9rem;}

@media print{
  html,body{background:#fff;}
  #nav,#loading{display:none!important;}
  #book{padding:0;}
  .page{display:block!important;width:auto;max-width:none;height:auto;
    box-shadow:none;margin:0;padding:0 0 1.2rem;page-break-after:always;}
  .page-body{height:auto;overflow:visible;}
  .folio{position:static;margin-top:.6rem;}
}
"""

JS = r"""
(function(){
  var src=document.getElementById('src'), book=document.getElementById('book'),
      nav=document.getElementById('nav'), ind=document.getElementById('pgind'),
      bPrev=document.getElementById('prev'), bNext=document.getElementById('next'),
      bToc=document.getElementById('toc'), load=document.getElementById('loading');
  var pages=[], cur=0, tocPage=0;

  function metrics(){
    var pr=document.createElement('div'); pr.className='page';
    pr.style.cssText='visibility:hidden;display:block;position:absolute;';
    book.appendChild(pr);
    var cs=getComputedStyle(pr),
        padV=parseFloat(cs.paddingTop)+parseFloat(cs.paddingBottom),
        w=pr.clientWidth;
    book.removeChild(pr);
    // 지면은 항상 화면 안에 다 들어와야 한다(넘김이지 스크롤이 아니다).
    // 화면이 충분히 높으면 A4 비율까지만 키운다.
    // 넘김 막대 높이는 CSS 에서 48px 로 고정돼 있다. offsetHeight 를 쓰면 폰트 로딩 전
    // 값이 잡혀 지면이 화면보다 작게 조판된다. 뷰포트 높이는 innerHeight 를 쓴다 —
    // quirks 모드에서 documentElement.clientHeight 는 뷰포트가 아니라 내용 높이를 준다.
    var avail=(window.innerHeight||document.documentElement.clientHeight) - 48 - 22;
    var h=Math.max(380, Math.min(Math.round(w*1.414), avail));
    return {h:h, body:Math.max(240, h-padV)};
  }
  function newPage(cls){
    var pg=document.createElement('div'); pg.className='page'+(cls?' '+cls:'');
    var pb=document.createElement('div'); pb.className='page-body';
    pg.appendChild(pb); book.appendChild(pg);
    var o={el:pg, body:pb}; pages.push(o); return o;
  }
  function isHead(n){ return n && /^H[234]$/.test(n.tagName); }

  function splitList(list, page, M){
    var tag=list.tagName.toLowerCase(),
        box=document.createElement(tag); box.className=list.className;
    page.body.appendChild(box);
    var kids=[].slice.call(list.children);
    for (var i=0;i<kids.length;i++){
      box.appendChild(kids[i]);
      if (page.body.scrollHeight>M.body){
        box.removeChild(kids[i]);
        if (!box.children.length){ box.appendChild(kids[i]); continue; }
        page=newPage();
        box=document.createElement(tag); box.className=list.className;
        page.body.appendChild(box); box.appendChild(kids[i]);
      }
    }
    return page;
  }
  // 그림이 남은 자리에 안 들어가면 지면 절반이 비어 버린다. LaTeX 의 float 처럼
  // 뒤로 미뤄 두고 본문을 마저 흘린 뒤, 다음 지면 머리에 앉힌다.
  var pending=[];
  // 캡션까지 합쳐 지면보다 큰 그림은 어느 지면에도 안 들어가 큐를 영영 막는다.
  // 빈 지면에서는 그림을 줄여서라도 반드시 앉힌다.
  function shrinkToFit(fig, page, M){
    var im=fig.querySelector('img'); if(!im) return;
    for (var k=0; k<14 && page.body.scrollHeight>M.body; k++){
      var h=parseFloat(im.style.maxHeight)||(M.body*0.66);
      im.style.maxHeight=Math.round(h*0.88)+'px';
    }
  }
  function flushPending(page, M){
    while (pending.length){
      var empty=!page.body.firstChild, fig=pending[0];
      page.body.appendChild(fig);
      if (page.body.scrollHeight<=M.body){ pending.shift(); }
      else if (empty){ shrinkToFit(fig, page, M); pending.shift(); }
      else { page.body.removeChild(fig); break; }
    }
    return page;
  }
  function place(node, page, M){
    page.body.appendChild(node);
    if (page.body.scrollHeight<=M.body) return page;
    page.body.removeChild(node);
    var listy = (node.tagName==='OL'||node.tagName==='UL') && node.children.length>1;
    if (!page.body.firstChild){                 // 빈 지면인데도 안 들어감
      if (listy) return splitList(node, page, M);
      page.body.appendChild(node);
      if (node.tagName==='FIGURE') shrinkToFit(node, page, M);
      return flushPending(newPage(), M);
    }
    if (node.tagName==='FIGURE'){ pending.push(node); return page; }
    var last=page.body.lastElementChild, np=newPage();
    if (isHead(last) && page.body.children.length>1) np.body.appendChild(last); // 제목 고아 방지
    flushPending(np, M);
    if (listy) return splitList(node, np, M);
    np.body.appendChild(node);
    if (np.body.scrollHeight>M.body && np.body.children.length>1){
      np.body.removeChild(node);
      var np2=flushPending(newPage(), M); np2.body.appendChild(node); return np2;
    }
    return np;
  }

  function build(){
    book.innerHTML=''; pages=[]; pending=[];
    book.classList.add('measuring');
    var M=metrics();
    document.documentElement.style.setProperty('--page-h', M.h+'px');
    var secs=[].slice.call(src.querySelectorAll('.sheet')), page=null;
    secs.forEach(function(sec){
      var cover=sec.classList.contains('cover');
      if (!page || page.body.firstChild){ page=newPage(cover?'p-cover':''); if(!cover) flushPending(page,M); }
      else if (cover) page.el.classList.add('p-cover');
      if (sec.classList.contains('toc')) tocPage=pages.length-1;
      [].slice.call(sec.children).forEach(function(k){
        var c=k.cloneNode(true);
        c.querySelectorAll && [].slice.call(c.querySelectorAll('img')).forEach(function(im){
          im.style.maxHeight=Math.round(M.body*0.66)+'px'; im.style.width='auto';
        });
        page=place(c, page, M);
      });
    });
    while (pending.length){                       // 미뤄 둔 그림 마무리
      var n0=pending.length, pg=flushPending(newPage(), M);
      if (pending.length===n0) pg.body.appendChild(pending.shift());
    }
    pages.forEach(function(p,i){
      var f=document.createElement('div'); f.className='folio';
      f.textContent = i===0 ? '' : (i+1)+' / '+pages.length;
      p.el.appendChild(f);
    });
    book.classList.remove('measuring');
    linkToc();
    var hp=/[#&?]p=(\d+)/.exec(location.hash||'');
    if (hp) cur=parseInt(hp[1],10)-1;
    if (cur>=pages.length) cur=pages.length-1;
    show(cur);
    load.style.display='none';
  }

  function linkToc(){
    var where={};
    pages.forEach(function(p,i){
      [].slice.call(p.body.querySelectorAll('h2,h3')).forEach(function(h){
        var t=h.textContent.trim(); if(!(t in where)) where[t]=i;
      });
      [].slice.call(p.body.querySelectorAll('figcaption')).forEach(function(c){
        var m=/그림\s*(\d+)/.exec(c.textContent); if(m && !('fig'+m[1] in where)) where['fig'+m[1]]=i;
      });
    });
    pages.forEach(function(p){
      [].slice.call(p.body.querySelectorAll('.toc-list li.go')).forEach(function(li){
        mark(li, where[li.querySelector('.tx').textContent.trim()]);
      });
      [].slice.call(p.body.querySelectorAll('.toc-figs li.go')).forEach(function(li){
        var m=/그림\s*(\d+)/.exec(li.textContent);
        mark(li, m?where['fig'+m[1]]:undefined);
      });
    });
    function mark(li,i){
      if (i===undefined){ li.classList.remove('go'); return; }
      li.querySelector('.pn').textContent=i+1;
      li.addEventListener('click', function(){ show(i); });
    }
  }

  function show(i){
    if (!pages.length) return;
    cur=Math.max(0,Math.min(pages.length-1,i));
    pages.forEach(function(p,k){ p.el.classList.toggle('on', k===cur); });
    ind.textContent=(cur+1)+' / '+pages.length;
    bPrev.disabled=(cur===0); bNext.disabled=(cur===pages.length-1);
    window.scrollTo(0,0);
  }
  bPrev.addEventListener('click',function(){show(cur-1);});
  bNext.addEventListener('click',function(){show(cur+1);});
  bToc.addEventListener('click',function(){show(tocPage);});
  document.addEventListener('keydown',function(e){
    if (e.key==='ArrowRight'||e.key==='PageDown'||e.key===' ') { show(cur+1); e.preventDefault(); }
    else if (e.key==='ArrowLeft'||e.key==='PageUp') { show(cur-1); e.preventDefault(); }
    else if (e.key==='Home') show(0);
    else if (e.key==='End') show(pages.length-1);
  });
  var x0=null,y0=null;
  book.addEventListener('touchstart',function(e){
    x0=e.changedTouches[0].clientX; y0=e.changedTouches[0].clientY;},{passive:true});
  book.addEventListener('touchend',function(e){
    if(x0===null) return;
    var dx=e.changedTouches[0].clientX-x0, dy=e.changedTouches[0].clientY-y0;
    if (Math.abs(dx)>45 && Math.abs(dx)>Math.abs(dy)*1.4) show(cur + (dx<0?1:-1));
    x0=null;},{passive:true});

  var t=null;
  window.addEventListener('resize',function(){
    clearTimeout(t); t=setTimeout(function(){ load.style.display=''; build(); },220);});
  window.addEventListener('beforeprint',function(){
    pages.forEach(function(p){p.el.classList.add('on');});});
  window.addEventListener('afterprint',function(){ show(cur); });

  // 첫 조판이 끝난 뒤 가용 높이가 달라졌으면(초기 레이아웃이 덜 자리잡은 경우) 한 번만 다시 짠다
  var settled=false;
  function buildOnce(){
    build();
    if (settled) return;
    requestAnimationFrame(function(){
      settled=true;
      var want=metrics().h,
          got=parseFloat(getComputedStyle(document.documentElement).getPropertyValue('--page-h'));
      if (Math.abs(want-got)>6) build();
    });
  }
  if (document.readyState==='complete') buildOnce();
  else window.addEventListener('load', buildOnce);
})();
"""

HTML = f"""<meta charset="utf-8">
<title>APEX — 검증된 그래픽 측광 파이프라인 (국문 프리뷰)</title>
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>{CSS}</style>
<div id="src">
{BODY}
</div>
<div id="loading">지면을 조판하는 중…</div>
<div id="book"></div>
<nav id="nav">
  <button id="prev" type="button">‹ 이전</button>
  <span id="pgind">– / –</span>
  <button id="next" type="button">다음 ›</button>
  <button id="toc" type="button">목차</button>
</nav>
<script>{JS}</script>"""
# 단독 파일은 doctype 이 있어야 한다. 없으면 quirks 모드로 떨어져 뷰포트 높이 계산이
# 틀리고 지면이 화면을 못 채운다. 아티팩트용 사본은 게시 때 skeleton 이 씌워지므로 뺀다.
OUT.write_text("<!doctype html>\n" + HTML, encoding="utf-8")
OUT_ARTIFACT.write_text(HTML, encoding="utf-8")
print("wrote", OUT, round(len(HTML)/1e6, 2), "MB | citations", len(LAB), "| figs", len(FIGMAP))
print("wrote", OUT_ARTIFACT.name, "(아티팩트용 — doctype 없음)")
