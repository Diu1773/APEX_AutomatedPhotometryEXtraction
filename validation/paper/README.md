# APEX paper-grade validation figures

Publication-quality validation figure set for APEX. Every figure is generated
by its own `figN_*.py` script from **self-contained experiments** (no external
data, fixed seeds) or from committed benchmark outputs, using the production
APEX code paths (real Step-4 detector, real `phot_vectorized`).

| Fig | Script | What it proves | Data |
| --- | --- | --- | --- |
| 1 | `fig1_completeness.py` | Completeness function, m50 depth + bootstrap CI | canonical injection run (`data/`) |
| 2 | `fig2_error_model.py` | Reported σ is honest: pull ~ N(0,1), CCD-equation floor | self-generated Monte-Carlo |
| 3 | `fig3_parameter_sweep.py` | Aperture optimum; depth vs sky / seeing trends | self-generated sweeps |
| 4 | `fig4_crosscheck_sep.py` | APEX ≡ independent `sep` C engine (synthetic truth) | self-generated frame |
| 5 | `fig5_crosscheck_iraf.py` | APEX ≡ IRAF/DAOPHOT on real NGC 457 g-band data | `benchmark/runs/ngc457_iraf_crosscheck_g0016_v1/` (committed) |
| 6 | `fig6_qc_validation.py` | Frame-QC decisions vs injected frame defects | self-generated synthetic night |
| 7 | `fig7_reference_crosscheck.py` | Faint drift is a Gaia BP artifact, not an APEX error | NGC 6811 reduction + PS1 (VizieR, cached) |
| 8 | `fig8_cmd_reproduction.py` | APEX/Gaia/PS1 CMDs agree (19 mmag ridgeline RMS) | NGC 6811 reduction + PS1 cache |
| 9 | `fig9_crowded_field.py` | No crowding-dependent bias in a real globular core (M5) | M5 reduction (re-run Steps 7/8/10) + NGC 6811 |

Shared infrastructure:

- `apex_paper_style.py` — one publication style for all figures (Okabe-Ito
  colorblind-safe palette, serif/CM math, 300 dpi, vector PDF + PNG on every save).
- `_make_canonical_data.py` — regenerates the canonical artificial-star dataset
  (840 injections / 12 trials, ~1 min) under `data/`.
- `assemble_figures.py` — writes `FIGURES.md` (figure + caption index) and
  `figures/overview_contact_sheet.png`.
- `run_all.py` — runs everything above in the right order.

## How to run

Always use the deploy venv interpreter (all deps installed):

```bash
.venv-deploy/Scripts/python.exe validation/paper/run_all.py            # everything
.venv-deploy/Scripts/python.exe validation/paper/run_all.py --only 2 5 # subset
.venv-deploy/Scripts/python.exe validation/paper/run_all.py --fast     # skip slow fig3 sweeps
```

Approximate runtimes: fig1/4/5 seconds; fig2 ~1 min; fig6 ~3-5 min;
fig3 ~10 min (19 full benchmark runs); fig7/8 ~1 min (PS1 query cached after
first run); fig9 seconds (consumes an already-reduced M5 tree); canonical
data ~1 min.

## Rules / gotchas

- **Windows MAX_PATH (260 chars):** benchmark runs write deeply nested cache
  paths. Figure outputs stay under the repo; heavy *run* directories go to a
  SHORT root (`C:\Users\<user>\AppData\Local\Temp\apx_*`). Never point a run
  at a deep scratch path.
- `data/` and `_sweep_ref.fits` are **regenerable and gitignored** — only
  scripts, captions, and final figures are committed. Rerun
  `_make_canonical_data.py` after cloning to rebuild `data/`.
- Fig 5 consumes the committed IRAF cross-check
  (`benchmark/runs/ngc457_iraf_crosscheck_g0016_v1/phot_fixed_coords/fixed_comparison.csv`);
  regenerating that requires PyRAF in WSL (`apex/benchmark/iraf_crosscheck.py`).
- Figs 5, 7, 8, 9 need the external data volume (`E:\observed_Analysis`),
  not just this repo — they are not reproducible from a bare checkout.
- Fig 9 (M5) requires Steps 7/8/10 to have been re-run against the current
  code with `parameters_M5.toml` (gitignored, uncommitted — see
  `scripts/run_step7_headless.py` / `run_step8_headless.py` /
  `run_step10_headless.py --params parameters_M5.toml`); it does not
  regenerate that reduction itself, only reads it.
- Captions live in `captions/figN_*.md` — paper-style "Figure N." text plus
  the exact numbers from the run that produced the committed figure.
