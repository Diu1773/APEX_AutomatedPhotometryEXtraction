# Figure — per-step preprocessing cross-check vs ccdproc

**(a)–(g)** |APEX − ccdproc| difference map for every calibration stage
(8×8 max-pooled over the full 3194×4788 frame):
master bias/dark/flat construction, bias/dark/flat application, and the full
end-to-end pipeline. Six stages are **bit-identical** (Δ = 0 at every pixel);
only the full chain shows float32 rounding at max|Δ| =
8.6e-04 DN (robust σ =
2.6e-05 DN). **(h)** That worst
disagreement against the detector read noise (3.45 DN) and sky shot
noise (41 DN) — more than three orders of magnitude below any
real noise term. Inputs: 8 bias, 8 darks (60 s),
5 flats, one 60 s NGC 6811 $B$ light (Moravian C3-61000, 2×2,
night 2026-06-11); reference astropy ccdproc 2.5.1. The cosmetic
(L.A.Cosmic + hot-pixel) stage is disabled here — it repairs ~1% of pixels by
design and is validated separately by injection; this figure isolates the
bias/dark/flat arithmetic. Generator: `calib_crosscheck_ngc6811.py` (input
file names recorded in the JSON).
