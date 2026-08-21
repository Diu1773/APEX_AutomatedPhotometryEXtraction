"""Which harmonic of a periodogram peak is the star's actual period.

A Lomb-Scargle periodogram fits one sinusoid, so it reports the frequency with
the most power — which is not always the star's period. An eclipsing binary
completes one orbit per period P but dips twice, at primary and secondary
eclipse, so the strongest single sinusoid sits at 2/P and the periodogram says
P/2. Fold at P/2 and the two eclipses land on top of each other; the curve looks
tidy and the answer is half of the truth.

Measured on ASAS-SN's variable star database, five detached binaries (`EA`),
900-1900 day baselines, 143-467 points each (2026-08-21):

    catalogue P   1.92079  1.30085  1.38502  3.33455  3.16715
    APEX adopted  0.96041  0.65042  0.69252  3.33302  1.58362
                   −50 %    −50 %    −50 %   −0.05 %   −50 %

Four of five, exactly half, and two of those were reported RESOLVED. Worse: the
candidate table did not contain the catalogue period in three cases, because
`build_alias_candidates` builds its list from the periodogram peak and the
sampling-window offsets — never from the harmonics. A table that is meant to
carry the answer for a person to pick did not have it.

The remedy has two parts, and both are needed.

**Put the harmonics on the candidate table.** P/3, P/2, P, 2P, 3P. Without this
the true period is not merely mis-ranked, it is absent — which it was for three
of the five stars above.

**Decide by folding, and only when three separate things agree.** Each closes a
different way of being wrong, and each was added because the one before it let
something through:

1. *The sub-cycles differ.* Fold at 2P, cut the cycle in half, compare the two
   halves bin by bin. A secondary eclipse shows up here; a curve that really
   repeats at P does not. Measured: 10-69 sigma on real binaries.

2. *The difference repeats across the run.* A second, unrelated pulsation
   frequency also makes sub-cycles differ — on a synthetic two-mode signal it
   reached 249 sigma, higher than any real binary here — but it drifts, so the
   difference profile from the first half of the observations does not match the
   second half. An eclipse falls at the same phase forever.

3. *The curve is flat between dips.* This is what keeps a pulsator's correct
   period from being doubled. Coherence alone looked like enough on six stars
   and was not on twenty-three: an RR Lyrae passed at 0.52, just over the line.

    flat fraction    EA 0.65-0.90   EB 0.42-0.69   EW 0.34-0.50   RRc 0.30-0.42

The last row is the honest cost, and it is not a threshold that wants tuning.
Contact binaries are tidally distorted and vary continuously, so folded at half
their orbital period they are the same shape as a first-overtone RR Lyrae folded
at its own period — not similar, the same. None of the twelve measured here is
doubled automatically.

Two ways out were measured and neither works (2026-08-21):

* *Colour.* Gaia BP-RP over ~780 stars per class: EW spans 0.56-1.44, the
  pulsators 0.40-1.37. Almost complete overlap.
* *Fourier shape.* The amplitude ratio R21 separates the classes beautifully
  **at the catalogue period** — EW 8.5-137 against pulsators 0.05-0.50 — and
  that is circular. At the periodogram peak, which is the only period available
  when the question is asked, an EW folded at P/2 and an RRc folded at P give
  the same R21, because they are the same curve.

What does separate them is absolute magnitude: RR Lyrae are standard candles and
contact binaries are not. That needs a parallax and a classifier, which is a
different piece of software from a period search. Until then their harmonics stay
on the candidate table for a person to choose from, which is what the table is
for.

Across 70 ASAS-SN stars of six kinds (2026-08-21):

    eclipsing   EA 11  EB 12  EW 12    catalogue period matched  0/35 -> 10/35
    pulsating   RRab 12  RRc 12  HADS 11   33/35 -> 33/35, doubled 0/35

Every adoption carried a sub-cycle difference of 8.8-1150 sigma and a coherence
of 0.75-0.99, well clear of the thresholds; nothing was adopted marginally. The
pulsator sample matters as much as the binary one — the whole risk of a rule
like this is that it "fixes" periods that were already right, and thirty-five
untouched pulsators is the evidence that it does not.

Three things were tried and removed rather than kept at a lower confidence: a x3
factor (invented, and it tripled a star whose doubling test had correctly
failed), a concentration measure meant to tell an eclipse from a second
frequency (0.78 against 0.78-0.98 — no separation), and iterating the doubling
to reach periodogram peaks that sit at P/4 (it overshot to x8). What this cannot
justify from the data, it does not do.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

DEFAULT_MULTIPLES: tuple[int, ...] = (2,)
"""Harmonics this may *adopt*, as opposed to merely list.

Doubling is the case with a mechanism behind it: two eclipses per orbit. Three
was included at first for "rarer triple-dip geometries", which was a guess and
behaved like one — on ASAS-SN J225021, whose two eclipses are indistinguishable
and whose doubling test correctly failed, the x3 fold passed both tests and
tripled the period. A factor this rule cannot justify from the data is a factor
it should not apply.

`build_alias_candidates` still puts x2 and x3 on the candidate table. Listing a
possibility and adopting it are different acts, and the table is for a person to
decide from.
"""

DEFAULT_BINS = 20
"""Fallback when the bin count is not chosen from the data."""

BIN_LADDER: tuple[int, ...] = (10, 20, 40, 80, 160)
"""Bin counts to choose from, coarsest first.

The count decides the answer, not just its precision. A detached binary's
eclipse can be a few per cent of the cycle wide, so at twenty bins the eclipse
and its flat surroundings share a bin and folding at half the period barely
changes the dispersion. Measured on ASAS-SN J234910 (449 points):

    20 bins   theta(P/2)=0.1601  theta(P)=0.2015   → half looks better
    50 bins   theta(P/2)=0.1355  theta(P)=0.0684   → the truth looks better

Same star, same data, opposite conclusion. Twenty bins was hiding a 17.7-sigma
difference between the two eclipse depths.
"""

MIN_POINTS_PER_BIN = 3
"""A bin with one or two points has no meaningful internal scatter."""

TARGET_POINTS_PER_BIN = 5
"""Occupancy the chosen bin count aims for.

Finer bins resolve a narrow eclipse but thin the counts, and a bin with two
points contributes almost no degrees of freedom — which flatters whichever
period happens to scatter its points into the emptiest bins.
"""


def choose_bins(n_points: int, ladder: tuple[int, ...] = BIN_LADDER) -> int:
    """The finest binning this many points can fill, coarsest as a floor."""
    usable = [n for n in ladder if n * TARGET_POINTS_PER_BIN <= int(n_points)]
    return max(usable) if usable else min(ladder)


@dataclass
class HarmonicVerdict:
    """What the fold test concluded, and the evidence for it."""

    adopted_period: float
    base_period: float
    factor: float
    reason: str
    candidates: list[dict] = field(default_factory=list)

    @property
    def changed(self) -> bool:
        return bool(np.isfinite(self.factor) and abs(self.factor - 1.0) > 1e-9)

    def as_dict(self) -> dict:
        return {
            "adopted_period": self.adopted_period,
            "base_period": self.base_period,
            "factor": self.factor,
            "reason": self.reason,
            "candidates": list(self.candidates),
        }


def _clean(time, mag, mag_err=None):
    t = np.asarray(time, dtype=float)
    y = np.asarray(mag, dtype=float)
    ok = np.isfinite(t) & np.isfinite(y)
    e = None
    if mag_err is not None:
        e = np.asarray(mag_err, dtype=float)
        ok &= np.isfinite(e) & (e > 0)
    return (t[ok], y[ok], (e[ok] if e is not None else None))


def fold_profile(time, mag, period: float, n_bins: Optional[int] = None,
                 mag_err=None) -> dict:
    """Binned mean and scatter of the light curve folded at `period`."""
    t, y, e = _clean(time, mag, mag_err)
    n_bins = int(n_bins) if n_bins else choose_bins(t.size)
    if t.size < n_bins * MIN_POINTS_PER_BIN or not np.isfinite(period) or period <= 0:
        return {"bins": np.array([]), "mean": np.array([]),
                "count": np.array([]), "sem": np.array([])}
    phase = np.mod(t / float(period), 1.0)
    index = np.minimum((phase * n_bins).astype(int), n_bins - 1)
    count = np.bincount(index, minlength=n_bins).astype(float)
    total = np.bincount(index, weights=y, minlength=n_bins)
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = np.where(count > 0, total / np.maximum(count, 1), np.nan)
        sq = np.bincount(index, weights=y * y, minlength=n_bins)
        var = np.where(count > 1, sq / np.maximum(count, 1) - mean ** 2, np.nan)
        var = np.clip(var, 0.0, None)
        sem = np.where(count > 1, np.sqrt(var / np.maximum(count, 1)), np.nan)
    return {"bins": (np.arange(n_bins) + 0.5) / n_bins,
            "mean": mean, "count": count, "sem": sem}


def phase_dispersion(time, mag, period: float,
                     n_bins: Optional[int] = None) -> float:
    """Stellingwerf (1978) theta: within-bin variance over total variance.

    Below 1 the fold explains structure; at or above 1 it explains nothing.
    Bins with fewer than two points carry no degrees of freedom and are skipped,
    which is what keeps a fold with mostly-empty bins from scoring well.
    """
    t, y, _ = _clean(time, mag)
    n_bins = int(n_bins) if n_bins else choose_bins(t.size)
    if t.size < n_bins * MIN_POINTS_PER_BIN or not np.isfinite(period) or period <= 0:
        return float("nan")
    total_var = float(np.var(y, ddof=1))
    if not np.isfinite(total_var) or total_var <= 0:
        return float("nan")
    phase = np.mod(t / float(period), 1.0)
    index = np.minimum((phase * n_bins).astype(int), n_bins - 1)
    count = np.bincount(index, minlength=n_bins).astype(float)
    total = np.bincount(index, weights=y, minlength=n_bins)
    sq = np.bincount(index, weights=y * y, minlength=n_bins)
    usable = count >= 2
    if not usable.any():
        return float("nan")
    ss = sq[usable] - total[usable] ** 2 / count[usable]
    dof = count[usable] - 1.0
    if dof.sum() <= 0:
        return float("nan")
    return float((ss.sum() / dof.sum()) / total_var)


def _removed_half_cycle_contrast() -> None:
    """Gone. It promised "how unequal the two eclipses are" and never delivered.

    First it compared the two deepest *bins*, which on a single-dip fold is the
    eclipse against the flat part — 21 sigma for a star whose eclipses agree to
    1.2. Rewritten to find local minima, it picked noise bumps. Rewritten again
    to roll-and-halve, it came out *smaller* at half the period than at the
    period, which is the opposite of its whole premise.

    Nothing in the decision used it; it was a number in the output that a reader
    would have trusted. Three attempts is enough to say the construction was
    wrong rather than the tuning.
    """




def harmonic_periods(base_period: float,
                     multiples: Sequence[int] = DEFAULT_MULTIPLES,
                     min_period: float = 0.0,
                     max_period: float = float("inf")) -> list[float]:
    """`base_period` scaled by each harmonic factor, inside the search window."""
    if not np.isfinite(base_period) or base_period <= 0:
        return []
    factors = [1.0]
    for m in multiples:
        m = int(m)
        if m > 1:
            factors.extend((float(m), 1.0 / m))
    out = []
    for factor in factors:
        period = float(base_period) * factor
        if min_period <= period <= max_period and period > 0:
            if not any(abs(period - kept) <= 1e-9 for kept in out):
                out.append(period)
    return out


def resolve_harmonic(
    time,
    mag,
    base_period: float,
    *,
    mag_err=None,
    multiples: Sequence[int] = DEFAULT_MULTIPLES,
    min_period: float = 0.0,
    max_period: float = float("inf"),
    n_bins: Optional[int] = None,
    sigma_threshold: float = 3.0,
    coherence_threshold: float = 0.5,
    flat_threshold: float = 0.60,
) -> HarmonicVerdict:
    """Adopt `m x base_period` when its sub-cycles demonstrably differ.

    The shorter period is the null hypothesis and it keeps the answer unless
    the longer fold shows structure the shorter one cannot hold. `m` is tried
    smallest first, so a star whose halves differ is doubled and not tripled.

    Doubling is applied at most once. A periodogram peak can sit at a quarter
    of the period — ASAS-SN J021350, catalogue 2.21434 d, peak 0.55358 d — and
    one doubling leaves that answer still halved. Iterating the test to reach it
    was tried and overshot to eight times the peak on the same star, so it is
    not here: an extension that cannot be shown to stop in the right place is
    worse than the shortfall it was meant to fix. Such a star keeps its
    unresolved harmonic on the candidate table.

    `sigma_threshold` is on `(chi2 - dof) / sqrt(2 dof)` comparing the
    sub-cycles bin by bin. Three sigma is where the ASAS-SN detached binaries
    measured here separate cleanly from the synthetic multi-mode signals that
    must not be doubled.
    """
    t, y, e = _clean(time, mag, mag_err)
    bins = int(n_bins) if n_bins else choose_bins(t.size)
    periods = harmonic_periods(base_period, multiples, min_period, max_period)
    if t.size == 0 or not periods:
        return HarmonicVerdict(float(base_period), float(base_period), 1.0,
                               "no usable points or no harmonic in range")

    rows = []
    for period in periods:
        factor = period / float(base_period)
        row = {
            "period": period,
            "factor": factor,
            "n_bins": bins,
            "theta": phase_dispersion(t, y, period, bins),
        }
        if factor > 1.0 + 1e-9:
            m = int(round(factor))
            row.update({f"subcycle_{k}": v for k, v in
                        subcycle_difference(t, y, period, m, bins, e).items()})
            row["subcycle_coherence"] = subcycle_coherence(t, y, period, m, e)
        row["flat_fraction"] = flat_fraction(t, y, period, bins, e)
        rows.append(row)

    for row in sorted((r for r in rows if r["factor"] > 1.0 + 1e-9),
                      key=lambda r: r["factor"]):
        sigma = float(row.get("subcycle_sigma", float("nan")))
        coherence = float(row.get("subcycle_coherence", float("nan")))
        flat = float(row.get("flat_fraction", float("nan")))
        # Three conditions, each closing a way of being wrong. Sigma: the
        # sub-cycles differ. Coherence: the difference repeats across the run,
        # so it is an eclipse and not a second unrelated frequency. Flatness:
        # the curve has flat stretches between dips, so it is an eclipsing
        # star at all.
        #
        # The last is what keeps a pulsator's period from being doubled. On 23
        # ASAS-SN stars, coherence alone let an RRc through at 0.52 — the gap
        # that looked clean on six stars was not there on twenty-three. Flatness
        # separates detached binaries (0.65-0.90) from RRc (0.30-0.42) with room
        # to spare. The cost is contact binaries: they vary continuously, sit at
        # 0.34-0.50, and are not auto-doubled. That is a real limit, not a
        # tuning failure — an EW light curve and an RRc light curve are the same
        # shape, and telling them apart needs colour or Fourier decomposition.
        if (np.isfinite(sigma) and sigma >= sigma_threshold
                and np.isfinite(coherence) and coherence >= coherence_threshold
                and np.isfinite(flat) and flat >= flat_threshold):
            m = int(round(row["factor"]))
            return HarmonicVerdict(
                adopted_period=row["period"],
                base_period=float(base_period),
                factor=row["factor"],
                reason=(f"the curve is flat between dips ({flat:.2f}) and the "
                        f"{m} sub-cycles of the x{m} fold differ at "
                        f"{sigma:.1f} sigma, repeating across "
                        f"the run (coherence {coherence:.2f}, largest gap "
                        f"{row.get('subcycle_max_gap', float('nan')):.3f} mag), so the "
                        f"periodogram peak is a harmonic"),
                candidates=rows,
            )

    base_row = next((r for r in rows if abs(r["factor"] - 1.0) <= 1e-9), rows[0])
    longer = [r for r in rows if r["factor"] > 1.0 + 1e-9]
    best = max((float(r.get("subcycle_sigma", float("nan"))) for r in longer),
               default=float("nan"))
    coh = max((float(r.get("subcycle_coherence", float("nan"))) for r in longer),
              default=float("nan"))
    if not np.isfinite(best):
        reason = "no harmonic fold splits into differing sub-cycles"
    elif best < sigma_threshold:
        reason = (f"the harmonic folds repeat themselves (best {best:.1f} sigma, "
                  f"below {sigma_threshold:g}), so the periodogram period stands")
    else:
        reason = (f"a harmonic fold does split ({best:.1f} sigma) but the split "
                  f"does not repeat across the run (coherence {coh:.2f}, below "
                  f"{coherence_threshold:g}) — that is what a second, unrelated "
                  f"frequency looks like, so the periodogram period stands and "
                  f"the harmonic remains a candidate")
    return HarmonicVerdict(
        adopted_period=float(base_row["period"]),
        base_period=float(base_period),
        factor=1.0,
        reason=reason,
        candidates=rows,
    )


def subcycle_difference(time, mag, period: float, m: int,
                        n_bins: Optional[int] = None, mag_err=None) -> dict:
    """Do the `m` sub-cycles of a fold at `period` differ from each other?

    This is the question the harmonic choice actually turns on, and it compares
    a fold with *itself* rather than with another fold — so it is not biased by
    how many bins a longer period gets, and it does not thin the counts the way
    matching phase width does.

    Fold at `period`, cut the cycle into `m` equal sub-cycles, and compare their
    binned profiles. If the star really repeats at `period / m`, the sub-cycles
    are copies and the difference is noise. If it repeats at `period` — an
    eclipsing binary, whose primary and secondary differ — the sub-cycles do
    not match, and that mismatch is the evidence for the longer period.

    Returns chi-square against the bins' own standard errors, its degrees of
    freedom, and the equivalent sigma. Sigma uses the normal approximation
    `(chi2 - dof) / sqrt(2 dof)`, which is what makes "3 sigma" mean the usual
    thing here.
    """
    t, y, e = _clean(time, mag, mag_err)
    n_bins = int(n_bins) if n_bins else choose_bins(t.size)
    m = int(m)
    blank = {"chi2": float("nan"), "dof": 0, "sigma": float("nan"),
             "max_gap": float("nan"), "concentration": float("nan")}
    if m < 2 or not np.isfinite(period) or period <= 0 or t.size == 0:
        return blank
    per_cycle = n_bins // m
    if per_cycle < 3:
        return blank

    profile = fold_profile(t, y, period, per_cycle * m, e)
    mean, sem, count = profile["mean"], profile["sem"], profile["count"]
    if mean.size != per_cycle * m:
        return blank

    blocks = mean.reshape(m, per_cycle)
    errs = sem.reshape(m, per_cycle)
    counts = count.reshape(m, per_cycle)
    usable = np.isfinite(blocks) & np.isfinite(errs) & (errs > 0) & (counts >= 2)
    # A phase bin is only comparable where every sub-cycle has points in it.
    columns = np.flatnonzero(usable.all(axis=0))
    if columns.size < 3:
        return blank

    chi2, dof, gaps = 0.0, 0, []
    per_bin = np.zeros(columns.size, dtype=float)
    for a in range(m - 1):
        for b in range(a + 1, m):
            diff = blocks[a, columns] - blocks[b, columns]
            err = np.hypot(errs[a, columns], errs[b, columns])
            contribution = (diff / err) ** 2
            per_bin += contribution
            chi2 += float(np.sum(contribution))
            dof += int(columns.size)
            gaps.append(float(np.max(np.abs(diff))))
    if dof <= 0:
        return blank
    sigma = (chi2 - dof) / np.sqrt(2.0 * dof)

    # Where the difference sits, not just how big it is. A secondary eclipse is
    # a few phase bins deep and flat elsewhere, so most of the chi-square comes
    # from a handful of bins. A second, unrelated pulsation frequency also makes
    # the sub-cycles differ — it has a different phase in each one — but that
    # difference is spread over the whole cycle. Without this, a two-mode signal
    # doubles its period at 249 sigma.
    top = max(1, int(np.ceil(0.2 * per_bin.size)))
    order = np.argsort(per_bin)[::-1]
    concentration = (float(per_bin[order[:top]].sum()) / float(per_bin.sum())
                     if per_bin.sum() > 0 else float("nan"))
    return {"chi2": float(chi2), "dof": int(dof), "sigma": float(sigma),
            "max_gap": float(max(gaps)) if gaps else float("nan"),
            "concentration": concentration}


def subcycle_coherence(time, mag, period: float, m: int = 2,
                       mag_err=None) -> float:
    """Is the sub-cycle difference the same in the first half of the data as the second?

    Sub-cycles differ for two quite different reasons and one fold cannot tell
    them apart: a secondary eclipse (the period really is `period`) or a second,
    unrelated pulsation frequency (it is not). Both raise the chi-square — on a
    synthetic two-mode signal the difference reached 249 sigma, higher than any
    of the real binaries measured here.

    Time tells them apart. An eclipse falls at the same phase every cycle, so
    the difference profile measured on the first half of the observations
    matches the one from the second half. An unrelated frequency drifts out of
    phase, so the two profiles are uncorrelated. Measured (2026-08-21):

        synthetic two-mode          -0.23
        ASAS-SN J225021 (no real secondary)   0.01
        ASAS-SN J230059              0.29
        ASAS-SN J234910              0.97
        ASAS-SN J182310              1.00
        ASAS-SN J032449              0.98

    J230059 is a real doubling that this misses — 191 points split in half
    leaves too few per bin for its shallow secondary to survive. It stays on the
    candidate table; it is simply not adopted automatically.
    """
    t, y, e = _clean(time, mag, mag_err)
    if t.size < 40 or not np.isfinite(period) or period <= 0 or int(m) < 2:
        return float("nan")
    m = int(m)
    order = np.argsort(t)
    t, y = t[order], y[order]
    e = e[order] if e is not None else None
    cut = t.size // 2

    profiles = []
    for part in (slice(0, cut), slice(cut, None)):
        bins = (choose_bins(t[part].size) // m) * m
        if bins < 2 * m:
            return float("nan")
        folded = fold_profile(t[part], y[part], period, bins,
                              (e[part] if e is not None else None))
        mean = folded["mean"]
        if mean.size != bins:
            return float("nan")
        blocks = mean.reshape(m, -1)
        profiles.append(blocks[0] - blocks[1] if m == 2
                        else blocks[0] - blocks[1:].mean(axis=0))
    first, second = profiles
    if first.size != second.size:
        return float("nan")
    ok = np.isfinite(first) & np.isfinite(second)
    if ok.sum() < 4:
        return float("nan")
    if np.std(first[ok]) == 0 or np.std(second[ok]) == 0:
        return float("nan")
    return float(np.corrcoef(first[ok], second[ok])[0, 1])


def flat_fraction(time, mag, period: float, n_bins: Optional[int] = None,
                  mag_err=None, level: float = 0.25) -> float:
    """Fraction of the folded cycle spent near maximum brightness.

    An eclipsing binary sits at constant brightness between eclipses, so most
    of its folded cycle is flat and bright. A pulsator is always on its way up
    or down. Measured on ASAS-SN (2026-08-21):

        EA  detached   0.65 - 0.90
        EB  beta Lyr   0.42 - 0.69
        EW  contact    0.34 - 0.50
        RRc pulsator   0.30 - 0.42

    Detached binaries separate from pulsators with room to spare. Contact
    binaries do not, and that is physics rather than a measurement problem:
    tidally distorted stars vary continuously, so an EW light curve and an RRc
    light curve are the same shape. Telling those two apart needs colour or
    Fourier decomposition, not one fold.
    """
    t, y, e = _clean(time, mag, mag_err)
    n_bins = int(n_bins) if n_bins else choose_bins(t.size)
    profile = fold_profile(t, y, period, n_bins, e)
    mean = profile["mean"]
    mean = mean[np.isfinite(mean)]
    if mean.size < 6:
        return float("nan")
    bright, faint = np.percentile(mean, 5), np.percentile(mean, 95)
    if faint - bright <= 0:
        return float("nan")
    return float(np.mean((mean - bright) / (faint - bright) < float(level)))
