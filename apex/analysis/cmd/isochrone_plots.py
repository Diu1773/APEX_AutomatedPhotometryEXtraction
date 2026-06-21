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

    fig, ax = plt.subplots(figsize=(6.6, 7.8))

    if member_mask is not None and np.any(~np.asarray(member_mask, bool)):
        m = np.asarray(member_mask, bool)
        ax.scatter(obs_color[~m], obs_mag[~m], s=4, alpha=0.18, c="#b8b8b8",
                   linewidths=0, label="Field / non-member", zorder=1)
        ax.scatter(obs_color[m], obs_mag[m], s=6, alpha=0.55, c="#444444",
                   linewidths=0, label="Members", zorder=2)
    else:
        ax.scatter(obs_color, obs_mag, s=5, alpha=0.40, c="#555555",
                   linewidths=0, label="Observed", zorder=2)

    iso_color = np.asarray(iso_color, float)
    iso_mag = np.asarray(iso_mag, float)
    fin = np.isfinite(iso_color) & np.isfinite(iso_mag)
    ax.plot(iso_color[fin], iso_mag[fin], "-", lw=1.8, c="#d62728",
            alpha=0.95, label="Best-fit isochrone", zorder=3)

    xlo, xhi = np.nanpercentile(obs_color, [1, 99])
    ylo, yhi = np.nanpercentile(obs_mag, [1, 99])
    xpad = max(0.05, 0.08 * (xhi - xlo))
    ypad = max(0.20, 0.08 * (yhi - ylo))
    ax.set_xlim(xlo - xpad, xhi + xpad)
    ax.set_ylim(yhi + ypad, ylo - ypad)  # magnitudes increase downward
    ax.set_xlabel(color_label, fontsize=12)
    ax.set_ylabel(mag_label, fontsize=12)
    ax.grid(True, ls=":", alpha=0.30)
    ax.legend(loc="upper right", framealpha=0.9, fontsize=9)
    ax.set_title(title, fontsize=12)

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
            ax.text(0.03, 0.03, "\n".join(lines), transform=ax.transAxes,
                    fontsize=9.5, va="bottom", ha="left",
                    bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#999999",
                              alpha=0.9))
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
