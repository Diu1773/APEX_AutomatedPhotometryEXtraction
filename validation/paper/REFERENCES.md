# APEX paper — Annotated bibliography (verified)

*ars-lit-review output, 2026-07-05. All 21 references confirmed against NASA ADS / arXiv / CrossRef.
BibTeX in [`references.bib`](references.bib). Zero fabricated DOIs. Corrections vs prior memory noted inline.*

## Verification status: 23/23 confirmed (21 initial + `lang2010`, `astap` added for §2)

> §2 additions: **`lang2010`** — Astrometry.net (Lang et al. 2010, AJ 139, 1782; DOI 10.1088/0004-6256/139/5/1782) VERIFIED, the fallback plate solver. **`astap`** — ASTAP (Kleijn), the primary plate solver: software only, no refereed paper or ASCL record found; cite the URL/GitHub and record the version used.


| # | Key | Reference | Status |
|---|---|---|---|
| 1 | `stetson1987` | DAOPHOT — PASP 99, 191 | ✅ VERIFIED |
| 2 | `bertin1996` | SExtractor — A&AS 117, 393 | ✅ VERIFIED |
| 3 | `collins2017` | AstroImageJ — AJ 153, 77 | ✅ VERIFIED |
| 4 | `tsiaras2019` | HOPS — EPSC-DPS2019-1594 | ⚠️ conference abstract only (no paper) |
| 5 | `hroch2014` | Munipack/MuniWin — ASCL:1402.006 | ⚠️ software record (no paper) |
| 6 | `benn2012` | VStar — JAAVSO 40, 852 | ✅ VERIFIED |
| 7 | `photutils` | photutils — Zenodo | ⚠️ update to installed-version DOI |
| 8 | `brennan2022` | AutoPhOT — A&A 667, A62 | ✅ VERIFIED |
| 9 | `riello2021` | Gaia EDR3 photometric validation — A&A 649, **A3** | ✅ VERIFIED (not A1!) |
| 10 | `pancino2022` | Gaia standard stars — A&A 664, A109 | ✅ VERIFIED |
| 11 | `tamuz2005` | SYSREM — MNRAS 356, 1466 | ✅ VERIFIED |
| 12 | `stellingwerf1978` | PDM — ApJ 224, 953 | ✅ VERIFIED (no DOI, pre-DOI era) |
| 13 | `naylor2006` | τ² CMD fitting — MNRAS 373, 1251 | ✅ VERIFIED (**not** completeness-likelihood) |
| 14 | `yang2018` | YZ Boo — RAA 18, 2 | ✅ VERIFIED (P=0.104091579 d) |
| 15 | `zhou2001` | AE UMa — A&A 374, 235 | ✅ VERIFIED (single author) |
| 16 | `paparo2001` | AE UMa — A&A 368, 880 | ⚠️ verify full author list |
| 17 | `gaia_dr3_2023` | Gaia DR3 summary — A&A 674, A1 | ✅ VERIFIED |
| 18 | `chambers2016` | PS1 Surveys — arXiv:1612.05560 | ⚠️ arXiv-only preprint |
| 19 | `astropy2013` | Astropy — A&A 558, A33 | ✅ VERIFIED |
| 20 | `astropy2018` | Astropy v2.0 — AJ 156, 123 | ✅ VERIFIED |
| 21 | `astropy2022` | Astropy v5.0 — ApJ 935, 167 | ✅ VERIFIED |

## §3.12 additions (2026-08-01) — detection-threshold contamination

Added for the sign-flipped contamination measurement (Section~\ref{sec:threshold}).
Every DOI below was resolved and the landing page checked against the title;
one initial guess (`10.1086/339315` for `hopkins2002`) turned out to be a
different paper and was discarded rather than shipped.

| # | Key | Reference | Status |
|---|---|---|---|
| 22 | `serra2012` | Using Negative Detections to Estimate Source-Finder Reliability — PASA 29, 296 | ✅ VERIFIED (10.1071/AS11065) |
| 23 | `serra2015` | SoFiA: a flexible source finder for 3D spectral line data — MNRAS 448, 1922 | ✅ VERIFIED (10.1093/mnras/stv079) |
| 24 | `molino2014` | The ALHAMBRA Survey: Bayesian photometric redshifts — MNRAS 441, 2891 | ✅ VERIFIED (10.1093/mnras/stu387) |
| 25 | `hopkins2002` | A New Source Detection Algorithm Using the False-Discovery Rate — AJ 123, 1086 | ✅ VERIFIED (10.1086/338316) |

- **`serra2012`** — the method's origin: negative detections give each positive
  detection a probability of being real, on the assumption that the noise is
  symmetric and real sources have positive flux. APEX uses the same ratio
  $(N_+-N_-)/N_+$ but as one number per frame, not per detection.
- **`serra2015` (SoFiA)** — the implementation that made it standard (WALLABY).
  Cited to show the method is established practice, and to mark what APEX does
  *not* borrow: SoFiA's kernel-density estimate over a 3-D parameter space needs
  hundreds of negative detections, which a single 2-D optical frame does not
  supply (NGC 6811 yields none at $3.2\sigma$).
- **`molino2014` (ALHAMBRA)** — the optical precedent, and the closest prior
  work: SExtractor run on the inverted image to pick a threshold at 3 per cent
  contamination. The difference APEX reports is that the safe floor moves frame
  to frame, so one threshold chosen once is not sufficient.
- **`hopkins2002` (FDR)** — the alternative way to set a threshold statistically
  (Benjamini--Hochberg multiple-testing control, implemented in SFIND 2.0).
  Cited as the road not taken: it operates per pixel and does not compose with
  a minimum-connected-area criterion.

## Annotations & section mapping

### Tooling landscape → §1 Introduction (contrast set that defines the niche)
- **`stetson1987` / `bertin1996`** — the standard detection/photometry algorithms APEX builds on (DAOPHOT crowded-field, SExtractor extraction). Cited to establish "standard algorithms, not novel"; also the lineage APEX's forced-aperture + empirical-PSF choices descend from.
- **`collins2017` (AstroImageJ)** — the primary GUI comparator: powerful for single-field differential/transit photometry, but not aimed at cluster-CMD chains or multi-night LC building. Defines the scope gap.
- **`tsiaras2019` (HOPS)** — "HOlomon Photometric Software", a Python transit-photometry GUI used by ~70% of ExoClock citizen-science participants. **This is the author's origin-experience tool** (the "software that just works" that motivated APEX). Transit-focused → motivates extending GUI convenience to cluster/variable-star work. *Cite abstract + GitHub; no refereed paper exists.*
- **`hroch2014` (Munipack/MuniWin)`, `benn2012` (VStar)** — variable-star tools the author attempted but found hard to use; represent the "alternative GUI still has a barrier" point.
- **`photutils`** — the Python-library alternative (flexible but scripting-bound); also an actual APEX dependency → belongs in Data/Code + Introduction.

### Validation methodology → §3 Validation
- **`brennan2022` (AutoPhOT)** — the model paper: an automated photometry pipeline validated and published in A&A. APEX's paper follows this "tool + rigorous validation" template. Cite as precedent for the paper's own form.
- **`riello2021`** — documents Gaia BP faint-source flux systematics; the *reference-artifact* explanation for §3.7's B-band drift. **Must cite A3, not the A1 summary paper.**
- **`pancino2022`** — Gaia→Johnson-Kron-Cousins standard-star transformations; supports the Gaia-reference construction used in §3.7–3.8.
- **`naylor2006`** — τ² max-likelihood CMD fitting. Relevant to §4.2/Discussion (why isochrone-parameter recovery is a separate, degenerate problem). **Correction: this is the τ² fitting statistic paper, not a completeness-likelihood paper as an earlier note assumed.**

### LC methods → §2 Design (LC branch) + §4 Science application
- **`tamuz2005` (SYSREM)`, `stellingwerf1978` (PDM)** — the detrending and period-search algorithms in APEX's LC mode. Cited when describing/using LC mode (YZ Boo).

### Science application → §4
- **`yang2018`** — literature YZ Boo period (0.104091579 d) that APEX's recovery is compared against.
- **`zhou2001` / `paparo2001`** — AE UMa dual-mode literature; cited in the *future-work* framing (AE UMa needs prewhitening, deferred).

### Infrastructure → §2 Design + Data/Code Availability
- **`gaia_dr3_2023`** — Gaia DR3 reference catalog (calibration, cross-match).
- **`chambers2016`** — Pan-STARRS1, the independent reference in §3.7–3.8. arXiv-only.
- **`astropy2013/2018/2022`** — core dependency (WCS, FITS, coordinates); all three cited per Astropy policy.

## ⚠️ Fix before submission (5 items)
1. **`riello2021`** — ensure the manuscript cites **A3** (photometric validation), never A1 (summary).
2. **`photutils`** — swap the concept DOI for the version-specific Zenodo DOI of the installed release.
3. **`paparo2001`** — pull the complete/ordered author list directly from ADS (`2001A%26A...368..880P`).
4. **`tsiaras2019` (HOPS) + `hroch2014` (Munipack)** — cite as software/abstract with GitHub/ASCL; do not invent a journal ref or DOI.
5. **`chambers2016`** — cite as arXiv e-print (no journal volume/page).

## Coverage / gaps
- **Well covered:** photometry algorithms, GUI-tool comparators, Gaia/PS1 references, LC period methods, δ Sct science targets, the AutoPhOT precedent.
- **Added for §2:** Astrometry.net (`lang2010`) and ASTAP (`astap`) — plate solvers, now cited in §2.3.
- **Used in §3 (verify/optionally add):** SEP (`sep`) — §3.5 currently cites `bertin1996` (SExtractor, the shared algorithm); optionally add the sep software ref (Barbary 2016, JOSS 1, 58 — **verify before adding**). IRAF/DAOPHOT — §3.6 cites `stetson1987` (DAOPHOT); optionally add an IRAF software ref (Tody 1986/1993 — verify) if §3.6 leans on IRAF itself.
- **Possible additions during drafting:** an aperture-photometry-theory / CCD-equation reference (e.g. Howell) for §3.2's error model; a δ Scuti / HADS review for §4 context. Add only if the drafted text actually leans on them.
