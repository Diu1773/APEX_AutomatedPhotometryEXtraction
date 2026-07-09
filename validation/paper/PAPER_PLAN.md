# APEX Paper — Chapter Plan (ars-plan output)

*Working spec for turning `validation/paper/PAPER.md` into a submittable RASTI paper.
Produced 2026-07-05 via Socratic plan mode. This drives ars-lit-review → ars-outline → drafting.*

---

## Locked decisions

| Item | Value | Rationale |
|---|---|---|
| **Venue** | **RASTI** (RAS Techniques & Instruments) | Purpose-built for software/technique/instrument validation; open access; no algorithmic-novelty requirement; AUTOPHOT-adjacent. PASP = stretch alternative. |
| **Contribution frame** | **Integrated tool + rigorous validation** | AUTOPHOT model. Not a new algorithm — an accessible end-to-end GUI pipeline, validated to research grade with a *reproducible* validation suite. |
| **Structure** | Techniques-paper 6 sections (not strict IMRaD) | Intro → Design → Validation → Science app → Discussion → Conclusion. |
| **LC (science app)** | Reserve §4, aim to fill with **YZ Boo** (clean, no prewhitening) | AE UMa = explicit future work (needs prewhitening dev). LC mode not yet validated → section stays provisional. |
| **AI-assisted development** | **Balanced frame** | One honest sentence in Intro motivation (AI lowered the barrier to *building* tools, not just using them) + full AI-disclosure statement. Paper's claim rests on tool+validation, not on how it was built. |

## Thesis (one sentence)

> A GUI photometry pipeline built from standard algorithms — one that takes a non-expert
> from raw frames to a science product (cluster CMD, variable-star light curve) without
> scripting — can pass a multi-layer, independent, *reproducible* validation and reach
> research-grade accuracy.

## Honest scope (state once, explicitly, up front and in Discussion)

- **Validated & publishable now:** shared Steps 1–7 measurement chain (detection → forced aperture/PSF photometry), automatic frame QC, the CMD *product* diagram (Fig 8), crowded-field behaviour (Fig 9). All *measurement-level* claims.
- **Forthcoming (reserved §4):** LC-mode science application (YZ Boo period recovery). LC mode still being fixed (period/prewhitening); section written only when the result is solid.
- **Explicitly out of scope:** automated isochrone-parameter recovery (degenerate inverse problem — §Discussion), sub-resolution blending finer than this dataset resolves, cross-instrument generality (single camera: Moravian C3-61000).

---

## Chapter Plan

### §1 Introduction  (~900 words · NEW writing)

- **Purpose:** motivate APEX from the accessibility gap (author's lived experience, depersonalized) and declare the tool+validation contribution.
- **Core argument:** existing photometry tooling forces a choice between legacy power with a brutal onboarding cliff (IRAF/DAOPHOT) and GUI convenience that is either scripting-bound (photutils) or scoped to transits/single-field differential work (HOPS, AstroImageJ); nothing guides a non-expert from raw frames to a *cluster CMD* or *multi-night light curve*. APEX fills that gap and earns trust through validation.
- **Evidence / beats:** Hook (the "software that just works" experience, generalized from HOPS transit analysis) → tooling landscape and its two-axis gap → gap statement (1 para) → APEX in one paragraph + contribution statement + one sentence on AI-assisted development lowering the build barrier → paper roadmap.
- **Tools cited (landscape, not all personally used):** IRAF/DAOPHOT (Stetson 1987), SExtractor (Bertin & Arnouts 1996), photutils, HOPS, AstroImageJ (Collins et al. 2017), MuniWin/Munipack, VStar (AAVSO). *Author used HOPS + IRAF/DS9 through to results; attempted VStar/MuniWin-class tools but found the UI prohibitive. photutils/AIJ cited as landscape, not claimed as personal pain points.*
- **Risks to defend:** (1) "why another tool vs AIJ/photutils?" → niche must be crisp (scope + integration + reproducible headless core); (2) skepticism of AI-assisted scientific code → the entire §3 validation is the answer.

### §2 APEX: Design & Implementation  (~1200 words · NEW writing — the current gap)

- **Purpose:** tell the reader what APEX *is* before validating it.
- **Depth policy:** standard algorithms cited + summarized (DAOPHOT/SExtractor/astropy lineage); detail only APEX-specific choices — empirical PSF extracted from the frame, forced aperture photometry at fixed coordinates, sky-annulus handling, and the Qt-free / headless reproducible core.
- **Content:** architecture overview (GUI + Qt-free core split); the shared measurement chain Steps 1–7 (file selection → crop → sky preview → detection → WCS plate solve → master catalog → forced aperture photometry) **in detail** — this is the validated spine; then CMD and LC branches **each in one short, provisional paragraph** (LC not yet validated, CMD interpretation still settling → describe, don't overclaim).
- **Figures:** (F-arch) workflow/step diagram; (F-gui) one representative GUI screenshot (from `ui_screenshots/`).
- **Risk:** scope discipline — do not describe mode features whose validation isn't in §3.

### §3 Validation  (~2500 words · REORGANIZE existing PAPER.md §2–3)

- **Purpose:** the evidence core. Port Figs 1–9 almost verbatim from current `PAPER.md`, re-sectioned under this heading.
- **Structure:** Data & Methods (synthetic self-contained + real-data cross-checks + scope-discipline note) → results in the existing order: completeness (F1) · error model (F2) · parameter/observing sweeps (F3) · sep cross-check (F4) · IRAF/DAOPHOT cross-check (F5) · frame QC + the honest transparency blind spot (F6) · Gaia-BP-vs-PS1 reference artifact (F7) · CMD reproduction across systems (F8) · crowded-field M5+M13 replication (F9).
- **Design principle (keep, it's the paper's spine):** every check is against something APEX did not produce (independent engine / independent catalog / injected truth); no free parameter is fit to the data being judged; the one anomaly investigated (B-band drift) resolved to a *reference* artifact, not an APEX error — reported precisely because a validation that only confirms is not a validation.
- **Status:** material complete, scripts committed, reproducible via `run_all.py`.

### §4 Science application: LC mode  (~800 words · RESERVED, provisional)

- **Purpose:** demonstrate the pipeline end-to-end on a real science result (the AUTOPHOT-style "science application" section).
- **Plan:** YZ Boötis (HADS δ Sct, P = 0.104092 d; Yang et al. 2018) — clean, fundamental-mode + harmonics, recoverable by LS/PDM without prewhitening. Re-run LC Steps 1–11 with current code; compare recovered period/amplitude to literature.
- **Explicit future work in text:** AE UMa (dual-mode, needs prewhitening — feature not yet implemented). The prior-code 8% period discrepancy (0.079150 vs 0.086017 d) is the honest motivating anecdote for why prewhitening matters — mention as future work, do not present as a result.
- **Gate:** section is written only once the YZ Boo result is solid. Until then it is a placeholder in the outline.

### §5 Discussion & Limitations  (~1000 words · REORGANIZE existing §4)

- What is established (measurement chain behaves as a well-calibrated instrument).
- Out of scope: isochrone-parameter recovery (degeneracy); LC mode until §4 lands.
- Limits on generality: **single instrument (verified from FITS headers, not assumed)** — the genuine independence axes are the cross-check *methods*, not the choice of cluster; sub-resolution crowding untested; M13 broken-Gaia-match cautionary sub-case.
- Practical recommendations (aperture ~1.0–1.2 FWHM, two-stage QC, faint-end reference choice, internal aperture-vs-PSF crowding diagnostic).

### §6 Conclusion  (~350 words · REORGANIZE existing §5)

- Answer the thesis; restate the measurement-level claim, deliberately scoped away from isochrone fitting, LC (until §4), sub-resolution blending, and un-re-verified instruments/clusters.

### Back matter

- Data & Code Availability (repo, MIT license, `run_all.py`, external data volume note).
- AI-usage disclosure (balanced frame — generate via ars-disclosure for RASTI).
- Author Contributions (CRediT — solo author), Conflict of Interest, Funding.
- Appendix A — Reproducibility table (exists in current PAPER.md).

---

## INSIGHT Collection

- `[INSIGHT: thesis]` — standard-algorithm accessible GUI pipeline + multi-layer reproducible validation = research-grade, trustworthy. Contribution is tool+validation, not algorithm.
- `[INSIGHT: motivation]` — author lived the three barriers: (1) the "software that just works" ideal (HOPS transit analysis), (2) legacy onboarding cliff + doc rot + parameter opacity (IRAF, personally, through to results), (3) alternative GUI tools too hard (VStar/MuniWin). Gap = guided GUI from raw frames to cluster CMD / multi-night LC. AI-assisted dev made a non-expert single-author build feasible.
- `[INSIGHT: niche]` — vs IRAF (onboarding), vs photutils (scripting), vs HOPS/AIJ (transit/single-field scope), vs MuniWin/VStar (narrow/hard). APEX = integrated detection→WCS→master→forced-phot→CMD/LC in one GUI + reproducible headless core.
- `[INSIGHT: honesty_spine]` — the validation only earns trust because it checks against what APEX didn't produce and reports the disconfirming case (B-band → Gaia BP artifact). This is also the reviewer-facing answer to "can AI-built science code be trusted?"
- `[INSIGHT: scope_discipline]` — validated (measurement chain, QC, CMD product, crowding) vs forthcoming (LC/YZ Boo) vs out-of-scope (isochrone fitting, sub-resolution blending, cross-instrument). Stated up front and in Discussion.

## Next steps (skill sequence)

1. **ars-lit-review** — pin exact citations for the tooling landscape + validation references (Stetson 1987, Bertin & Arnouts 1996, Collins et al. 2017, Riello et al. 2021, Tamuz et al. 2005, Stellingwerf 1978, Yang et al. 2018, Brennan & Fraser 2022, etc.). Verify each (no fabricated refs).
2. **ars-outline** — expand this Chapter Plan into a paragraph-level outline + evidence map.
3. **Drafting** (Opus) — §1 and §2 are net-new; §3/§5/§6 are re-editing existing PAPER.md prose into RASTI voice; §4 waits on YZ Boo.
4. ars-citation-check → ars-disclosure (RASTI) → ars-format-convert (LaTeX + .bib) → ars-reviewer (mock peer review).
