"""Paper-quality figures for CMD isochrone fitting (Qt-free, reusable).

Two publication-style figures are produced from an MCMC isochrone fit:

* :func:`cmd_figure` — the colour-magnitude diagram with the best-fit isochrone
  overlaid and a parameter annotation box.
* :func:`corner_figure` — the posterior corner plot (uses the ``corner`` package
  when available, with a self-contained fallback grid otherwise).

Both return a Matplotlib ``Figure`` so a GUI can embed it via
``FigureCanvasQTAgg`` *and* can be written straight to disk with
:func:`save_figure` (or the ``save_*`` convenience wrappers). This module imports
Matplotlib lazily and never imports Qt, so it is safe to use from the headless
pipeline, the CLI, the validation runner, and the GUI alike.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Sequence

import numpy as np

__all__ = [
    "cmd_figure",
    "segmented_isochrone_line",
    "corner_figure",
    "save_figure",
    "save_cmd_plot",
    "save_corner_plot",
]


# Per-parameter display helpers: label + value formatting for the annotation box.
_PARAM_LABELS = {
    "log_age": "log(age)",
    "age_gyr": "age",
    "mh": "[M/H]",
    "dm": "(m−M)₀",
    "e_color": "E(color)",
    "e_bv": "E(B−V)",
}


def _lazy_plt():
    import matplotlib

    matplotlib.use("Agg")  # safe default; a Qt canvas overrides this when embedding
    import matplotlib.pyplot as plt

    return plt


def segmented_isochrone_line(
    color: np.ndarray,
    magnitude: np.ndarray,
    *,
    min_jump: float = 0.5,
    jump_factor: float = 8.0,
    local_window: int = 8,
):
    """Insert NaN breaks where adjacent isochrone rows are discontinuous.

    Isochrone tables can append remnant / separately-computed evolutionary blocks
    after the visible stellar sequence; connecting those rows with a normal line
    draws long diagonals across the CMD. Splitting the *display* line at large
    jumps (relative to the local median step) keeps the curve interpolated within
    each evolutionary segment without joining the endpoints across gaps. Mirrors
    the CMD-viewer's behaviour. Pure numpy (Qt-free).
    """
    x = np.asarray(color, dtype=float)
    y = np.asarray(magnitude, dtype=float)
    if x.shape != y.shape or x.ndim != 1 or len(x) < 2:
        return x.copy(), y.copy()
    steps = np.hypot(np.diff(x), np.diff(y))
    break_after = ~(np.isfinite(x[:-1]) & np.isfinite(y[:-1])
                    & np.isfinite(x[1:]) & np.isfinite(y[1:]))
    for index, step in enumerate(steps):
        if break_after[index] or not np.isfinite(step) or step <= min_jump:
            continue
        start = max(0, index - local_window)
        stop = min(len(steps), index + local_window + 1)
        nearby = np.concatenate((steps[start:index], steps[index + 1:stop]))
        nearby = nearby[np.isfinite(nearby) & (nearby > 0)]
        local_step = float(np.median(nearby)) if len(nearby) else 0.0
        if local_step == 0.0 or step > jump_factor * local_step:
            break_after[index] = True
    if not break_after.any():
        return x.copy(), y.copy()
    extra = int(break_after.sum())
    line_x = np.empty(len(x) + extra, dtype=float)
    line_y = np.empty(len(y) + extra, dtype=float)
    j = 0
    for i in range(len(x)):
        line_x[j] = x[i]; line_y[j] = y[i]; j += 1
        if i < len(break_after) and break_after[i]:
            line_x[j] = np.nan; line_y[j] = np.nan; j += 1
    return line_x, line_y


def cmd_figure(
    obs_color: np.ndarray,
    obs_mag: np.ndarray,
    iso_color: np.ndarray,
    iso_mag: np.ndarray,
    *,
    color_label: str = "g−r",
    mag_label: str = "g",
    title: str = "CMD with best-fit isochrone",
    annotations: Optional[dict] = None,
    member_mask: Optional[np.ndarray] = None,
):
    """Return a Matplotlib ``Figure`` of the CMD with the isochrone overlaid.

    Parameters
    ----------
    obs_color, obs_mag :
        Observed colour and magnitude of every plotted star.
    iso_color, iso_mag :
        The best-fit isochrone track (evolutionary order), drawn as a line.
    color_label, mag_label :
        Axis labels (e.g. ``"g-r"`` / ``"g"``).
    title :
        Figure title.
    annotations :
        Optional ``{key: (median, lo, hi)}`` or ``{key: value}`` of fitted
        parameters (keys from :data:`_PARAM_LABELS`) rendered in a corner box.
    member_mask :
        Optional boolean mask; ``True`` stars are drawn as members (darker), the
        rest as faint field stars, so membership selection is visible.
    """
    plt = _lazy_plt()
    obs_color = np.asarray(obs_color, float)
    obs_mag = np.asarray(obs_mag, float)

    # Publication style: monochrome, minimal chrome.
    fig, ax = plt.subplots(figsize=(5.6, 6.6))

    if member_mask is not None and np.any(~np.asarray(member_mask, bool)):
        m = np.asarray(member_mask, bool)
        ax.scatter(obs_color[~m], obs_mag[~m], s=5, alpha=0.20, c="0.72",
                   linewidths=0, label="Field", zorder=1, rasterized=True)
        ax.scatter(obs_color[m], obs_mag[m], s=7, alpha=0.65, c="0.30",
                   linewidths=0, label="Members", zorder=2, rasterized=True)
    else:
        ax.scatter(obs_color, obs_mag, s=6, alpha=0.55, c="0.45",
                   linewidths=0, label="Observed", zorder=2, rasterized=True)

    iso_color = np.asarray(iso_color, float)
    iso_mag = np.asarray(iso_mag, float)
    # Interpolated line, segmented at large jumps so the giant tip / remnant
    # blocks don't connect across the CMD (matches the CMD-viewer behaviour).
    line_x, line_y = segmented_isochrone_line(iso_color, iso_mag)
    ax.plot(line_x, line_y, "-", lw=1.6, c="black", solid_capstyle="round",
            label="Best-fit isochrone", zorder=3)

    xlo, xhi = np.nanpercentile(obs_color, [1, 99])
    ylo, yhi = np.nanpercentile(obs_mag, [1, 99])
    xpad = max(0.05, 0.08 * (xhi - xlo))
    ypad = max(0.20, 0.08 * (yhi - ylo))
    ax.set_xlim(xlo - xpad, xhi + xpad)
    ax.set_ylim(yhi + ypad, ylo - ypad)  # magnitudes increase downward
    ax.set_xlabel(color_label, fontsize=12)
    ax.set_ylabel(mag_label, fontsize=12)
    ax.tick_params(direction="in", top=True, right=True)
    ax.legend(loc="upper right", frameon=False, fontsize=9, handletextpad=0.4)
    if title:
        ax.set_title(title, fontsize=11)

    if annotations:
        lines = []
        for key, val in annotations.items():
            lab = _PARAM_LABELS.get(key, key)
            if isinstance(val, (tuple, list)) and len(val) == 3:
                med, lo, hi = val
                lines.append(f"{lab} = {med:.3f}$^{{+{hi - med:.3f}}}_{{-{med - lo:.3f}}}$")
            elif val is not None and np.isfinite(val):
                lines.append(f"{lab} = {val:.3f}")
        if lines:
            # lower-left (faint/blue corner is usually empty); plain box, no colour.
            ax.text(0.04, 0.04, "\n".join(lines), transform=ax.transAxes,
                    fontsize=9, va="bottom", ha="left",
                    bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.6", lw=0.8))
    fig.tight_layout()
    return fig


def corner_figure(
    flat_chain: np.ndarray,
    labels: Sequence[str],
    truths: Optional[Sequence[float]] = None,
):
    """Return a posterior corner ``Figure`` (uses ``corner`` if installed)."""
    plt = _lazy_plt()
    flat_chain = np.asarray(flat_chain, float)
    labels = list(labels)
    try:
        import corner  # type: ignore

        fig = corner.corner(
            flat_chain,
            labels=labels,
            truths=list(truths) if truths is not None else None,
            quantiles=[0.16, 0.5, 0.84],
            show_titles=True,
            title_fmt=".3f",
        )
        return fig
    except Exception:
        pass

    # Self-contained fallback: pairwise scatter + marginal histograms.
    n = flat_chain.shape[1]
    fig, axes = plt.subplots(n, n, figsize=(2.4 * n, 2.4 * n))
    axes = np.atleast_2d(axes)
    for i in range(n):
        for j in range(n):
            ax = axes[i, j]
            if i == j:
                ax.hist(flat_chain[:, i], bins=40, color="#33bb66", alpha=0.85)
                if truths is not None:
                    ax.axvline(truths[i], color="#d62728", lw=1)
            elif j < i:
                ax.scatter(flat_chain[:, j], flat_chain[:, i], s=2, alpha=0.15,
                           c="#333333", linewidths=0)
            else:
                ax.axis("off")
            if i == n - 1 and ax.axison:
                ax.set_xlabel(labels[j])
            if j == 0 and ax.axison:
                ax.set_ylabel(labels[i])
    fig.tight_layout()
    return fig


def save_figure(fig, output, *, dpi: int = 160) -> Path:
    """Write *fig* to *output* (creating parent dirs) and close it. Returns path."""
    plt = _lazy_plt()
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    return output


def save_cmd_plot(output, obs_color, obs_mag, iso_color, iso_mag, *,
                  color_label="g−r", mag_label="g",
                  title="CMD with best-fit isochrone",
                  annotations=None, member_mask=None, dpi: int = 160) -> Path:
    """Convenience: build the CMD figure and save it in one call."""
    fig = cmd_figure(obs_color, obs_mag, iso_color, iso_mag,
                     color_label=color_label, mag_label=mag_label, title=title,
                     annotations=annotations, member_mask=member_mask)
    return save_figure(fig, output, dpi=dpi)


def save_corner_plot(output, flat_chain, labels, truths=None, *, dpi: int = 130) -> Path:
    """Convenience: build the corner figure and save it in one call."""
    fig = corner_figure(flat_chain, labels, truths=truths)
    return save_figure(fig, output, dpi=dpi)
