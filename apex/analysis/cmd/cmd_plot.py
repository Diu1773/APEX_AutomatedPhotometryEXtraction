"""Step 11 CMD figure for a batch run.

The desktop CMD viewer is an interactive instrument — region-of-interest
selection, parallax and proper-motion sliders, quality masks, view switching,
a Teff colour bar. None of that survives a translation to a batch run, and
pretending otherwise would produce a worse viewer rather than a figure.

What a batch run needs is the picture: the calibrated colour-magnitude diagram
of everything that passed, from the same table the viewer loads
(`median_by_ID_filter_wide_cmd.csv`), on the same axes, with the same
magnitude convention. That is what belongs in a paper and what a run should
leave behind. The viewer stays where it is, for looking.

Which colour and which magnitude are not guessed. They come from the config
(`isochrone.colors`, `isochrone.mag_band`) when it says, because that is the
pair the user already declared for the fit and the two figures should agree.
Failing that, the widest colour baseline available in the table is used and
named in the caption, so nobody has to reverse-engineer the axes.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from apex.utils.step_paths_cmd import step10_zp_dir, step11_cmd_dir

CMD_TABLE = "median_by_ID_filter_wide_cmd.csv"


def _bands(df: pd.DataFrame, prefix: str = "mag_std_") -> list[str]:
    """Photometric bands present as calibrated magnitudes, error columns aside."""
    out = []
    for col in df.columns:
        if not col.startswith(prefix) or col.startswith(prefix + "err"):
            continue
        band = col[len(prefix):]
        if band and pd.to_numeric(df[col], errors="coerce").notna().any():
            out.append(band)
    return out


def choose_axes(df: pd.DataFrame, params=None) -> tuple[str, str, str, str]:
    """(blue band, red band, magnitude band, why) for the colour-magnitude axes.

    The config wins when it names a pair, so this figure and the isochrone fit
    are drawn on the same axes. Otherwise pick the widest baseline available and
    say so — an unexplained axis choice is how two figures of one dataset end up
    looking like two datasets.
    """
    prefix = "mag_std_"
    available = _bands(df, prefix)
    if not available:
        prefix = "mag_inst_"
        available = _bands(df, prefix)
    if len(available) < 2:
        return "", "", "", "not enough bands"

    P = getattr(params, "P", None)
    if P is not None:
        raw = str(getattr(P, "iso_colors", "") or "").strip()
        first = raw.replace(";", ",").split(",")[0].strip()
        if "-" in first:
            blue, _, red = first.partition("-")
            blue, red = blue.strip(), red.strip()
            if blue in available and red in available:
                mag = str(getattr(P, "iso_mag_band", "") or "").strip()
                if mag not in available:
                    mag = red
                return blue, red, mag, "from isochrone.colors"

    # Widest baseline: bluest and reddest of what is here, in the usual order.
    order = ["u", "U", "B", "g", "V", "b", "v", "r", "R", "i", "I", "z", "Y"]
    ranked = sorted(available, key=lambda b: order.index(b) if b in order else 99)
    blue, red = ranked[0], ranked[-1]
    return blue, red, red, "widest available baseline"


def draw_cmd(fig, df: pd.DataFrame, params=None) -> bool:
    """Colour-magnitude diagram of the calibrated table."""
    fig.clear()
    if df is None or df.empty:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No CMD table", transform=ax.transAxes,
                ha="center", va="center", color="gray")
        ax.set_xticks([]); ax.set_yticks([])
        return False

    blue, red, mag_band, why = choose_axes(df, params)
    if not blue:
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "Fewer than two calibrated bands",
                transform=ax.transAxes, ha="center", va="center", color="gray")
        ax.set_xticks([]); ax.set_yticks([])
        return False

    prefix = "mag_std_" if f"mag_std_{blue}" in df.columns else "mag_inst_"
    b = pd.to_numeric(df.get(f"{prefix}{blue}"), errors="coerce")
    r = pd.to_numeric(df.get(f"{prefix}{red}"), errors="coerce")
    m = pd.to_numeric(df.get(f"{prefix}{mag_band}"), errors="coerce")
    colour = b - r

    ok = np.isfinite(colour) & np.isfinite(m)
    if not ok.any():
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, "No finite colour/magnitude pairs",
                transform=ax.transAxes, ha="center", va="center", color="gray")
        return False

    err_col = f"{prefix}err_{mag_band}"
    err = pd.to_numeric(df.get(err_col), errors="coerce") if err_col in df.columns else None

    ax = fig.add_subplot(111)
    if err is not None and np.isfinite(err[ok]).any():
        # Colour by measurement error rather than a flat scatter: the shape of a
        # CMD is set by the faint end, and the faint end is where the error is.
        scatter = ax.scatter(colour[ok], m[ok], c=err[ok], s=7, alpha=0.75,
                             cmap="viridis", linewidths=0, rasterized=True)
        bar = fig.colorbar(scatter, ax=ax, pad=0.02)
        bar.set_label(f"{prefix}err_{mag_band} (mag)", fontsize=9)
    else:
        ax.scatter(colour[ok], m[ok], s=7, alpha=0.75, color="#212121",
                   linewidths=0, rasterized=True)

    ax.invert_yaxis()                       # brighter is up, as a CMD is read
    ax.set_xlabel(f"{blue} - {red}", fontsize=10)
    ax.set_ylabel(mag_band, fontsize=10)
    ax.grid(True, alpha=0.25)
    system = "standard" if prefix == "mag_std_" else "instrumental"
    ax.set_title(
        f"CMD | {int(ok.sum())} stars | {system} | axes {why}", fontsize=10)
    fig.tight_layout()
    return True


def export_cmd_plot(result_dir, params=None) -> list[Path]:
    """Write `step11_cmd.png` from the Step 10 table. Returns what was written."""
    from matplotlib.figure import Figure

    table = Path(step10_zp_dir(result_dir)) / CMD_TABLE
    if not table.exists():
        return []
    try:
        df = pd.read_csv(table)
    except Exception:                               # noqa: BLE001
        return []

    out_dir = Path(step11_cmd_dir(result_dir))
    out_dir.mkdir(parents=True, exist_ok=True)
    fig = Figure(figsize=(7.0, 7.6), dpi=120)
    if not draw_cmd(fig, df, params):
        return []
    path = out_dir / "step11_cmd.png"
    fig.savefig(path, dpi=160, bbox_inches="tight")
    return [path]
