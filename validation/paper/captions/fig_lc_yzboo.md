# Figure (LC) — YZ Boötis: period reproduction and multi-night aliasing

**Figure.** APEX LC-mode reproduction of the high-amplitude δ Scuti (HADS)
star **YZ Boötis**, whose literature period is
$P = 0.104091579$ d (Yang et al. 2018, RAA 18, 2 = arXiv:1709.08798).
**(a)** A single good night (2026-03-28, $r$ band, $N = 80$ over 5.2 h)
folded on the *literature* period reproduces the published high-amplitude
sawtooth — fast rise, slow decline — with a peak-to-peak amplitude of
$\approx 0.39$ mag ($r$; consistent with the literature $V$ peak-to-peak
$\sim 0.42$ mag, since $r$-band amplitudes of δ Sct stars are slightly
lower than $V$). An unconstrained Lomb-Scargle search of that single night
recovers $P = 0.1046$ d, within **0.5%** of the literature value — the
residual is the expected period resolution of a 5.2 h ($\approx 2$-cycle)
baseline. **(b)** Lomb-Scargle periodograms: the single night (blue) has one
clean broad peak at the true period, whereas a 2-night merge with a $\sim$1-day
gap (orange) produces a dense alias comb whose *tallest* peak jumps to
$0.0946$ d — the $+1$ cycle-day$^{-1}$ alias of the true frequency — while
the true period (dashed line) survives only as a secondary peak.

**Honest scope.** This reproduces the known period and pulse shape of a real
variable star from APEX photometry, validating the LC-mode measurement +
period-analysis chain on a single night. It also exposes a real limitation,
reported rather than hidden: the current period analysis takes the single
tallest periodogram peak, so on multi-night ground-based data with strong
1-day sampling aliases it can select an alias instead of the true period —
exactly the "$\sim$0.95/0.095 d" behaviour seen in practice. This is the
standard spectral-window aliasing of single-site time series, not a
photometric error (the true period is present among the top peaks, and the
single-night solution disambiguates it); robustly resolving it on merged
multi-night data needs alias-aware peak selection or pre-whitening, which is
tracked as a separate development item (see also AE UMa, a double-mode HADS
whose beat requires pre-whitening not yet implemented).
