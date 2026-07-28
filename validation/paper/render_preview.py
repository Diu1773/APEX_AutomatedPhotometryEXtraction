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

# ── 표제면: A&A 667, A62 (AutoPhOT 논문) 판면을 그대로 따른다 ──
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
    <div class="jl">RAS Techniques and Instruments<br><span>투고 준비 원고 · 국문 검토용</span></div>
    <div class="jr">RASTI</div>
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
body{font-family:var(--serif);font-size:12.2px;line-height:1.36;color:var(--ink);
  font-variant-numeric:lining-nums;-webkit-font-smoothing:antialiased;}

/* ── 지면 ── */
#src{display:none;}
#stage{overflow:hidden;}
#book{transform-origin:top left;width:var(--pw);}
.page{position:relative;width:var(--pw);height:var(--ph);background:#fff;
  margin:0 0 14px;box-shadow:0 2px 10px rgba(0,0,0,.35);overflow:hidden;}
.pinner{position:absolute;left:var(--mx);right:var(--mx);top:var(--mt);bottom:var(--mb);
  display:flex;flex-direction:column;}
.span:not(:empty){margin-bottom:10px;}
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
h2 + p, h3 + p, h4 + p, figure + p, .tw + p, blockquote + p{text-indent:0;}
/* 제목·캡션은 어절 중간에서 끊기면 안 된다. 본문 단만 음절 단위 줄바꿈을 허용한다. */
h2,h3,h4,.ptitle,.abshead,.toc-h,figcaption,thead th{word-break:keep-all;}
h2{font-family:var(--sans);font-size:12.8px;font-weight:700;line-height:1.25;
  margin:11px 0 4px;text-indent:0;}
h3{font-family:var(--sans);font-size:11.6px;font-weight:700;line-height:1.25;
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
figcaption{font-size:10.2px;line-height:1.36;text-align:justify;margin:4px 0 0;text-indent:0;}
.tw{margin:0 0 9px;}
table{border-collapse:collapse;width:100%;font-size:10px;line-height:1.3;
  font-variant-numeric:tabular-nums;}
thead th{border-top:1.1px solid #111;border-bottom:.7px solid #111;text-align:left;
  padding:2.6px 4px;font-weight:700;vertical-align:bottom;}
tbody td{padding:2.2px 4px;vertical-align:top;}
tbody tr:last-child td{border-bottom:1.1px solid #111;}

/* ── 뷰어 크롬 ── */
#hud{position:fixed;right:12px;bottom:12px;z-index:20;display:flex;gap:6px;
  align-items:center;background:rgba(28,30,34,.9);color:#eceef1;border-radius:4px;
  padding:6px 9px;font-family:var(--sans);font-size:12px;}
#hud button{font:inherit;background:#3a3e45;color:#eceef1;border:0;border-radius:3px;
  padding:4px 8px;cursor:pointer;}
#hud button:hover{background:#4b505a;}
#pgind{font-variant-numeric:tabular-nums;min-width:4.6em;text-align:center;}
#loading{position:fixed;inset:0;display:flex;align-items:center;justify-content:center;
  background:#8f9298;color:#f0f1f3;font-family:var(--sans);font-size:13px;z-index:30;}

@media print{
  @page{size:A4;margin:0;}
  html,body{background:#fff;}
  #hud,#loading{display:none!important;}
  #stage{overflow:visible;height:auto!important;}
  #book{transform:none!important;margin:0!important;}
  .page{box-shadow:none;margin:0;page-break-after:always;}
}
"""

JS = r"""
(function(){
  var PW=794;
  var src=document.getElementById('src'), book=document.getElementById('book'),
      stage=document.getElementById('stage'), load=document.getElementById('loading'),
      ind=document.getElementById('pgind');
  var pages=[], pending=[], scale=1;

  function el(t,c){ var e=document.createElement(t); if(c) e.className=c; return e; }
  function isHead(n){ return n && /^H[234]$/.test(n.tagName); }
  function isSpanBlock(n){ return n && (n.tagName==='FIGURE' || n.classList.contains('tw')); }

  function newPage(){
    var pg=el('div','page');
    pg.innerHTML='<div class="run"></div><div class="pinner">'+
      '<div class="span"></div><div class="cols"><div class="col"></div>'+
      '<div class="col"></div></div></div><div class="folio"></div>';
    book.appendChild(pg);
    var cs=pg.querySelectorAll('.col');
    var P={el:pg, span:pg.querySelector('.span'), cols:[cs[0],cs[1]], ci:0,
           inner:pg.querySelector('.pinner')};
    pages.push(P); return P;
  }
  function fits(c){ return c.scrollHeight<=c.clientHeight+1; }

  /* 전폭 요소(그림·표)는 지면 머리에 앉힌다. 판면의 58% 를 넘지 않게 줄인다. */
  function spanCap(P){ return P.inner.clientHeight*0.58; }
  function shrink(node, P){
    var im=node.querySelector('img'); if(!im) return;
    for (var k=0;k<16 && P.span.offsetHeight>spanCap(P);k++){
      var h=im.getBoundingClientRect().height||400;
      im.style.maxHeight=Math.round(h*0.9)+'px'; im.style.width='auto';
    }
  }
  function flushSpan(P){
    while (pending.length){
      var f=pending[0], wasEmpty=!P.span.firstChild;
      P.span.appendChild(f);
      if (P.span.offsetHeight>spanCap(P)){
        if (wasEmpty){ shrink(f,P); pending.shift(); }
        else { P.span.removeChild(f); break; }
      } else { pending.shift(); }
    }
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

  function place(node, P){
    if (isSpanBlock(node)){ pending.push(node); return P; }
    var c=P.cols[P.ci];
    c.appendChild(node);
    if (fits(c)) return P;
    c.removeChild(node);
    var listy=(node.tagName==='OL'||node.tagName==='UL')&&node.children.length>1;
    if (!c.firstChild){                       // 빈 단인데도 안 들어감
      if (listy) return splitList(node,P);
      c.appendChild(node);
      return advance(P);
    }
    var last=c.lastElementChild, orphan=null;
    if (isHead(last) && c.children.length>1){ orphan=last; c.removeChild(last); }
    var P2=advance(P);
    if (orphan) P2.cols[P2.ci].appendChild(orphan);
    if (listy) return splitList(node,P2);
    var c2=P2.cols[P2.ci];
    c2.appendChild(node);
    if (!fits(c2) && c2.children.length>1){
      c2.removeChild(node);
      var P3=advance(P2); P3.cols[P3.ci].appendChild(node); return P3;
    }
    return P2;
  }

  function build(){
    book.innerHTML=''; pages=[]; pending=[];
    var P=newPage();
    ['.titleblock','.absblock'].forEach(function(sel){          // 표제면은 전폭
      var n=src.querySelector(sel); if(n) P.span.appendChild(n.cloneNode(true));
    });
    var flow=src.querySelector('.flow');
    [].slice.call(flow.children).forEach(function(k){ P=place(k.cloneNode(true), P); });
    while (pending.length){                                     // 남은 그림 마무리
      var n0=pending.length, Q=fresh();
      if (pending.length===n0){
        var f=pending.shift(); Q.span.appendChild(f); shrink(f,Q);
      }
    }
    var toc=src.querySelector('.tocblock');                     // 차례는 2쪽에 1단 전폭
    if (toc){
      var T=el('div','page');
      T.innerHTML='<div class="run"></div><div class="pinner"><div class="span"></div>'+
        '<div class="cols"><div class="col"></div></div></div><div class="folio"></div>';
      T.querySelector('.col').appendChild(toc.cloneNode(true));
      book.insertBefore(T, pages[0].el.nextSibling);
      pages.splice(1,0,{el:T});
    }
    var N=pages.length;
    pages.forEach(function(p,i){
      if (i>0) p.el.querySelector('.run').textContent=
        'APEX — GUI 기반 구경·PSF 측광 파이프라인 (국문 원고)';
      var f=p.el.querySelector('.folio');
      f.className='folio '+((i+1)%2 ? 'r' : 'l');
      f.textContent='APEX, '+(i+1)+' / '+N+' 쪽';
    });
    linkToc();
    fit();
    load.style.display='none';
    var hp=/[#&?]p=(\d+)/.exec(location.hash||'');     // #p=7 로 그 쪽을 바로 연다
    if (hp){ var pg=pages[parseInt(hp[1],10)-1];
      if (pg) window.scrollTo(0, Math.max(0, pg.el.offsetTop*scale-6)); }
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

  function fit(){
    var w=stage.clientWidth;
    scale=Math.min(1,(w-20)/PW);
    book.style.transform='scale('+scale+')';
    book.style.marginLeft=Math.max(0,(w-PW*scale)/2)+'px';
    stage.style.height=Math.ceil(book.scrollHeight*scale)+'px';
    onScroll();
  }
  function goto(i){
    var p=pages[i]; if(!p) return;
    window.scrollTo({top:Math.max(0,p.el.offsetTop*scale-6), behavior:'smooth'});
  }
  function onScroll(){
    if (!pages.length) return;
    var y=(window.scrollY||window.pageYOffset||0)/scale + 40, cur=1;
    for (var i=0;i<pages.length;i++){ if (pages[i].el.offsetTop<=y) cur=i+1; else break; }
    ind.textContent=cur+' / '+pages.length;
  }
  window.addEventListener('scroll', onScroll, {passive:true});
  var t=null;
  window.addEventListener('resize', function(){ clearTimeout(t); t=setTimeout(fit,150); });
  document.getElementById('top').addEventListener('click', function(){ goto(0); });
  document.getElementById('toctop').addEventListener('click', function(){ goto(1); });

  if (document.readyState==='complete') build();
  else window.addEventListener('load', build);
})();
"""

HTML = f"""<meta charset="utf-8">
<title>APEX — GUI 기반 구경·PSF 측광 파이프라인 (국문 원고)</title>
<meta name="viewport" content="width=device-width,initial-scale=1">
<style>{CSS}</style>
<div id="src">
{BODY}
</div>
<div id="stage"><div id="book"></div></div>
<div id="loading">판면을 조판하는 중…</div>
<div id="hud">
  <button id="top" type="button">처음</button>
  <button id="toctop" type="button">차례</button>
  <span id="pgind">– / –</span>
</div>
<script>{JS}</script>"""

# 단독 파일은 doctype 이 있어야 한다. 없으면 quirks 모드로 떨어져 판면 계산이 틀린다.
# 아티팩트용 사본은 게시 때 skeleton 이 씌워지므로 뺀다.
OUT.write_text("<!doctype html>\n" + HTML, encoding="utf-8")
OUT_ARTIFACT.write_text(HTML, encoding="utf-8")
print("wrote", OUT, round(len(HTML) / 1e6, 2), "MB | citations", len(LAB), "| figs", len(FIGMAP))
print("wrote", OUT_ARTIFACT.name, "(아티팩트용 — doctype 없음)")
