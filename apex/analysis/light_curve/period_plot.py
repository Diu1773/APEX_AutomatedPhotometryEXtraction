"""The period summary figure — three panels the window drew and nobody else could.

Light curve, periodogram, phase fold. It is the picture that says whether a
period is real: a clean fold at the peak period, or a smear that says the peak
was an alias of the sampling window.

It lived on the window because that is where the canvas was, and the only thing
it truly needed the canvas for was its width — the layout collapses to a
stacked column when the pane is narrow. A saved figure has no pane, so the
batch path fixes a wide canvas and gets the full two-by-two grid, while the
window keeps its responsive behaviour by overriding three hooks.

The rest is arithmetic on arrays the analysis already produced, which is why it
moves cleanly: `RawLightCurveBuilder` had the same shape.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


from apex.analysis.light_curve.period_analysis_service import compute_ls
from apex.analysis.light_curve.period_io_service import ALL_FILTER_KEY
from apex.utils.common_helpers import normalize_filter_key as _normalize_filter_key
from pathlib import Path


_FILTER_PLOT_COLORS = {
    "u": "#7E57C2",
    "b": "#1E88E5",
    "v": "#43A047",
    "g": "#2A9D68",
    "r": "#D64B45",
    "i": "#7259C7",
    "z": "#8D6E63",
}


def _load_check_star_for_plot(result_dir: Path, filt: str | None = None):
    """Load check star CSV from lc_lightcurve/ for plotting. Returns (check_id, df_or_None)."""
    try:
        from apex.analysis.light_curve.check_star_io import load_check_star_csv
        check_id, df = load_check_star_csv(result_dir, filt=filt)
        return check_id, (df if not df.empty else None)
    except Exception:
        return None, None

class _DefaultPlotColors:
    """The palette a batch run draws with.

    Same values as the light preset in `apex.gui.theme`. They are duplicated
    rather than imported at module scope because `analysis` must not depend on
    `gui` — `_colors()` below picks up the live, possibly re-themed tokens
    whenever the GUI is installed, so a window's figure still follows the
    user's chosen theme.
    """

    PLOT_AXES_BG = "#FFFFFF"
    PLOT_FG = "#1F2933"
    PLOT_GRID = "#D0D5DC"
    TEXT_MUTED = "#8692A6"
    WARN = "#B26A00"


def _colors():
    """The window's live palette when there is one, the defaults otherwise."""
    try:
        from apex.gui.theme import Tokens as _live      # noqa: PLC0415
        return _live
    except Exception:                                   # noqa: BLE001
        return _DefaultPlotColors

class PeriodSummaryPlotter:
    """Draws the three-panel period summary.

    A user supplies `results` (from `run_period_analysis`), `lc_data` (from
    `load_period_lightcurve_csv`), and optionally `alias_analysis`,
    `check_star_data`, `params`. The window supplies them from its own state
    and overrides the three canvas hooks; a batch run uses the defaults below.
    """

    results: dict | None = None
    lc_data: dict | None = None
    alias_analysis: dict | None = None
    multimode_diagnostic: dict | None = None
    _summary_layout_compact = False
    _summary_layout_stacked = False

    # -- what the window reads from widgets ---------------------------------
    #
    # The figure needs four things a person would have typed: the search
    # window, how finely it was sampled, and whether to mark aliases. The
    # window reads them off its spin boxes; a batch run is handed them.

    search_min_period: float = 0.01
    search_max_period: float = 10.0
    search_samples_per_peak: int = 10
    show_alias_marks: bool = True
    params = None
    _check_star_override = None

    # Written on first use, read on the line before — the window set both in
    # its constructor, so a plotter built any other way has to declare them.
    # Immutable defaults only: a class-level dict would be shared by every
    # instance and hand back another run's check star.
    _check_star_plot_cache_key = None
    _check_star_plot_cache = None

    def _search_window(self) -> tuple[float, float, int]:
        return (float(self.search_min_period), float(self.search_max_period),
                int(self.search_samples_per_peak))

    def _show_alias_marks(self) -> bool:
        return bool(self.show_alias_marks)

    def _load_check_star(self, requested_filter):
        """The check star to overplot, if this run has one."""
        if self._check_star_override is not None:
            return self._check_star_override
        result_dir = getattr(getattr(self.params, "P", None), "result_dir", None)
        if result_dir is None:
            return (None, None)
        return _load_check_star_for_plot(Path(result_dir), requested_filter)

    @property
    def _check_star_ls_cache(self) -> dict:
        # Per instance, not per class — a class-level dict would be shared by
        # every plotter ever made and hand back another run's periodogram.
        if "_ls_cache" not in self.__dict__:
            self.__dict__["_ls_cache"] = {}
        return self.__dict__["_ls_cache"]

    @_check_star_ls_cache.setter
    def _check_star_ls_cache(self, value: dict) -> None:
        # The window assigns this in its constructor and a test replaces it
        # outright. A getter-only property turned both into AttributeError.
        self.__dict__["_ls_cache"] = value

    # -- the canvas seam -----------------------------------------------------

    def _summary_figure(self):
        """The figure to draw on. The window returns its live canvas figure."""
        if getattr(self, "_headless_figure", None) is None:
            from matplotlib.figure import Figure
            self._headless_figure = Figure(figsize=(11.0, 7.0), dpi=160)
        return self._headless_figure

    def _summary_canvas_width(self) -> int:
        """Pixel width the layout branches on. Wide enough for the full grid."""
        return 1760

    def _summary_redraw(self) -> None:
        """The window schedules a repaint; a saved figure needs nothing."""


    @staticmethod
    def _summary_uses_compact_layout(canvas_width: int) -> bool:
        """Use shorter labels and fixed spacing on constrained widths."""
        return int(canvas_width) < 900

    @staticmethod
    def _summary_uses_stacked_layout(canvas_width: int) -> bool:
        """Stack plots only while the controls are explicitly open."""
        return int(canvas_width) < 620

    def _summary_data_type(self) -> str:
        requested = str((self.alias_analysis or {}).get("input_series", "")).lower()
        if requested == "corrected":
            requested = "corr"
        if requested in {"raw", "corr"}:
            return requested
        has_corr = self.lc_data is not None and self.lc_data.get("mag_corr") is not None
        preserves_baseline = bool(
            self.lc_data.get("correction_preserves_nightly_baseline", True)
            if self.lc_data
            else True
        )
        return "corr" if has_corr and preserves_baseline else "raw"

    def _summary_period(self) -> float:
        adopted = float((self.alias_analysis or {}).get("adopted_period", np.nan))
        if np.isfinite(adopted) and adopted > 0:
            return adopted
        for method in ("ls", "pdm", "bls"):
            result = self._summary_method_result(method)
            period = float(result.get("best_period", np.nan)) if result else np.nan
            if np.isfinite(period) and period > 0:
                return period
        return np.nan

    def _summary_method_result(self, method: str) -> dict | None:
        dtype = self._summary_data_type()
        result = self.results.get(f"{dtype}_{method}")
        if isinstance(result, dict) and "error" not in result:
            return result
        fallback = self.results.get(f"raw_{method}")
        if isinstance(fallback, dict) and "error" not in fallback:
            return fallback
        return None

    @staticmethod
    def _median_center_for_plot(values: np.ndarray, filters: np.ndarray) -> np.ndarray:
        centered = np.asarray(values, dtype=float).copy()
        labels = np.asarray(filters, dtype=str)
        for label in np.unique(labels):
            mask = labels == label
            finite = mask & np.isfinite(centered)
            if np.any(finite):
                centered[mask] -= float(np.nanmedian(centered[finite]))
        return centered

    def _phase_plot_check_time_column(self, check_df: pd.DataFrame) -> str | None:
        generic_cols = ["BJD_TDB", "BJD", "bjd", "JD", "jd", "HJD", "hjd", "time"]
        base_time_col = str(self.lc_data.get("col_time", "") or "").strip() if self.lc_data else ""
        if base_time_col:
            if base_time_col in check_df.columns:
                return base_time_col
            for col in check_df.columns:
                if str(col).strip().lower() == base_time_col.lower():
                    return str(col)
        return next((c for c in generic_cols if c in check_df.columns), None)

    def _check_star_ls_result(self, check: dict) -> dict:
        key = (
            self._check_star_plot_cache_key,
            *self._search_window(),
        )
        cached = self._check_star_ls_cache.get(key)
        if cached is not None:
            return cached
        result = compute_ls(
            check["time"],
            check["mag"],
            check["mag_err"],
            "check",
            key[1],
            key[2],
            key[3],
        )
        self._check_star_ls_cache[key] = result
        return result

    def _check_star_plot_data(self) -> dict | None:
        if self.lc_data is None:
            return None
        selected_filter = str(self.lc_data.get("filter", ""))
        requested_filter = None if selected_filter == ALL_FILTER_KEY else selected_filter
        cache_key = (
            str(self.lc_data.get("source_file", "")),
            int(self.lc_data.get("target_id", 0)),
            selected_filter,
            self._summary_data_type(),
        )
        if cache_key == self._check_star_plot_cache_key:
            return self._check_star_plot_cache
        check_id, frame = self._load_check_star(requested_filter)
        if frame is None or frame.empty:
            self._check_star_plot_cache_key = cache_key
            self._check_star_plot_cache = None
            return None

        time_col = self._phase_plot_check_time_column(frame)
        dtype = self._summary_data_type()
        if dtype == "corr":
            mag_candidates = ["diff_mag", "diff_mag_corr", "diff_mag_raw", "mag"]
        else:
            mag_candidates = ["diff_mag_raw", "diff_mag", "mag"]
        mag_col = next((col for col in mag_candidates if col in frame.columns), None)
        if time_col is None or mag_col is None:
            self._check_star_plot_cache_key = cache_key
            self._check_star_plot_cache = None
            return None

        time = pd.to_numeric(frame[time_col], errors="coerce").to_numpy(float)
        if time_col == "rel_time_hr":
            time = time / 24.0
        mag = pd.to_numeric(frame[mag_col], errors="coerce").to_numpy(float)
        err_col = next(
            (col for col in ("diff_err", "diff_err_corr", "mag_err", "err") if col in frame.columns),
            None,
        )
        mag_err = (
            pd.to_numeric(frame[err_col], errors="coerce").to_numpy(float)
            if err_col
            else None
        )
        if "filter" in frame.columns:
            filters = (
                frame["filter"]
                .astype(str)
                .map(lambda value: _normalize_filter_key(value) or str(value).strip())
                .to_numpy(dtype=str)
            )
        else:
            filters = np.full(len(frame), requested_filter or "check", dtype=str)
        mag = self._median_center_for_plot(mag, filters)
        valid = np.isfinite(time) & np.isfinite(mag)
        if mag_err is not None:
            valid &= np.isfinite(mag_err) & (mag_err > 0)
        if not np.any(valid):
            self._check_star_plot_cache_key = cache_key
            self._check_star_plot_cache = None
            return None
        payload = {
            "check_id": check_id,
            "time": time[valid],
            "mag": mag[valid],
            "mag_err": mag_err[valid] if mag_err is not None else None,
            "filters": filters[valid],
        }
        self._check_star_plot_cache_key = cache_key
        self._check_star_plot_cache = payload
        return payload

    @staticmethod
    def _two_harmonic_phase_model(
        phase: np.ndarray,
        mag: np.ndarray,
        mag_err: np.ndarray | None,
        phase_grid: np.ndarray,
    ) -> np.ndarray | None:
        phase = np.asarray(phase, dtype=float)
        mag = np.asarray(mag, dtype=float)
        valid = np.isfinite(phase) & np.isfinite(mag)
        err = None if mag_err is None else np.asarray(mag_err, dtype=float)
        if err is not None:
            valid &= np.isfinite(err) & (err > 0)
        if np.count_nonzero(valid) < 8:
            return None

        q = phase[valid]
        y = mag[valid]
        design = np.column_stack(
            [
                np.ones(len(q)),
                np.sin(2.0 * np.pi * q),
                np.cos(2.0 * np.pi * q),
                np.sin(4.0 * np.pi * q),
                np.cos(4.0 * np.pi * q),
            ]
        )
        if err is not None:
            root_weight = 1.0 / err[valid]
            design_fit = design * root_weight[:, None]
            y_fit = y * root_weight
        else:
            design_fit = design
            y_fit = y
        try:
            coefficients, *_ = np.linalg.lstsq(design_fit, y_fit, rcond=None)
        except np.linalg.LinAlgError:
            return None

        grid = np.asarray(phase_grid, dtype=float) % 1.0
        grid_design = np.column_stack(
            [
                np.ones(len(grid)),
                np.sin(2.0 * np.pi * grid),
                np.cos(2.0 * np.pi * grid),
                np.sin(4.0 * np.pi * grid),
                np.cos(4.0 * np.pi * grid),
            ]
        )
        return grid_design @ coefficients

    def _update_summary_plot(self):
        fig = self._summary_figure()
        fig.clear()
        if not self.results or self.lc_data is None:
            self._summary_redraw()
            return

        canvas_width = self._summary_canvas_width()
        compact = self._summary_uses_compact_layout(canvas_width)
        stacked = self._summary_uses_stacked_layout(canvas_width)
        self._summary_layout_compact = compact
        self._summary_layout_stacked = stacked
        fig.set_layout_engine(None if compact else "constrained")
        if stacked:
            grid = fig.add_gridspec(3, 1, height_ratios=(1.05, 1.0, 1.0))
            ax_lc = fig.add_subplot(grid[0, 0])
            ax_period = fig.add_subplot(grid[1, 0])
            ax_phase = fig.add_subplot(grid[2, 0])
        else:
            grid = fig.add_gridspec(
                2, 2, height_ratios=(1.02, 1.0), width_ratios=(1.08, 0.92)
            )
            ax_lc = fig.add_subplot(grid[0, :])
            ax_period = fig.add_subplot(grid[1, 0])
            ax_phase = fig.add_subplot(grid[1, 1])

        dtype = self._summary_data_type()
        time = np.asarray(self.lc_data.get("time", []), dtype=float)
        raw_mag = self.lc_data.get("mag_corr") if dtype == "corr" else self.lc_data.get("mag_raw")
        if raw_mag is None:
            raw_mag = self.lc_data.get("mag_raw")
            dtype = "raw"
        mag = np.asarray(raw_mag, dtype=float)
        mag_err_raw = self.lc_data.get("mag_err")
        mag_err = None if mag_err_raw is None else np.asarray(mag_err_raw, dtype=float)
        filters = np.asarray(
            self.lc_data.get(
                "filter_values",
                np.full(len(time), str(self.lc_data.get("filter", "data"))),
            ),
            dtype=str,
        )
        mag = self._median_center_for_plot(mag, filters)
        valid = np.isfinite(time) & np.isfinite(mag)
        if mag_err is not None:
            valid &= np.isfinite(mag_err) & (mag_err > 0)
        time = time[valid]
        mag = mag[valid]
        filters = filters[valid]
        mag_err = mag_err[valid] if mag_err is not None else None
        if len(time) == 0:
            self._summary_redraw()
            return

        t0 = float(np.nanmin(time))
        baseline = float(np.nanmax(time) - t0)
        fallback_colors = ["#2A9D68", "#D64B45", "#7259C7", "#1E88E5", "#C77A12"]
        band_order = {band: index for index, band in enumerate(("u", "b", "v", "g", "r", "i", "z"))}
        filter_order = sorted(
            set(filters.tolist()),
            key=lambda label: (band_order.get(str(label).lower(), len(band_order)), str(label)),
        )
        filter_colors = {
            label: _FILTER_PLOT_COLORS.get(
                str(label).lower(), fallback_colors[idx % len(fallback_colors)]
            )
            for idx, label in enumerate(filter_order)
        }

        # Scatter only: connecting sequential multi-band points suggests an
        # interpolation that was never observed and exaggerates frame glitches.
        for label in filter_order:
            mask = filters == label
            ax_lc.scatter(
                (time[mask] - t0) * 24.0,
                mag[mask],
                s=18,
                color=filter_colors[label],
                alpha=0.82,
                edgecolors=_colors().PLOT_AXES_BG,
                linewidths=0.25,
                label=(
                    f"{label} (n={np.count_nonzero(mask)})"
                    if compact
                    else f"{label} target (n={np.count_nonzero(mask)})"
                ),
                zorder=3,
            )

        check = self._check_star_plot_data()
        if check is not None:
            ax_lc.scatter(
                (check["time"] - t0) * 24.0,
                check["mag"],
                s=10,
                marker="x",
                color=_colors().TEXT_MUTED,
                alpha=0.4,
                linewidths=0.7,
                label=(
                    f"check (n={len(check['time'])})"
                    if compact
                    else f"check star (n={len(check['time'])})"
                ),
                zorder=2,
            )
        ax_lc.axhline(0.0, color=_colors().PLOT_GRID, lw=0.8)
        ax_lc.set_title(
            "Observed differential light curve"
            if compact
            else "Observed differential light curve - per-filter median removed",
            loc="left",
        )
        ax_lc.set_xlabel("Hours from first exposure")
        ax_lc.set_ylabel(
            "Diff. mag"
            if compact
            else "Median-centered differential magnitude [mag]"
        )
        ax_lc.invert_yaxis()
        ax_lc.grid(True, alpha=0.25)
        ax_lc.legend(
            loc="upper right",
            ncol=(
                min(2, max(1, len(filter_order) + 1))
                if compact
                else min(4, max(1, len(filter_order) + 1))
            ),
            fontsize=7 if compact else 8,
        )

        period_colors = {"ls": "#1E88E5", "pdm": "#C77A12", "bls": "#7259C7"}
        period_labels = {"ls": "Lomb-Scargle", "pdm": "PDM (1-theta)", "bls": "BLS (normalized)"}
        period_lines = 0
        period_notes: list[str] = []
        for method in ("ls", "pdm", "bls"):
            result = self._summary_method_result(method)
            if result is None:
                continue
            if "frequency" in result:
                periods = 1.0 / np.asarray(result["frequency"], dtype=float)
            else:
                periods = np.asarray(result.get("trial_periods", []), dtype=float)
            power = np.asarray(result.get("power", []), dtype=float)
            finite = np.isfinite(periods) & np.isfinite(power) & (periods > 0)
            if not np.any(finite):
                continue
            periods = periods[finite]
            power = power[finite]
            order = np.argsort(periods)
            if method == "bls" and np.nanmax(np.abs(power)) > 0:
                power = power / np.nanmax(np.abs(power))
            ax_period.plot(
                periods[order], power[order], color=period_colors[method], lw=1.2,
                label=period_labels[method],
            )
            period_notes.append(f"{method.upper()} {float(result['best_period']):.6f} d")
            period_lines += 1

        if check is not None and self._summary_method_result("ls") is not None:
            check_ls = self._check_star_ls_result(check)
            if "error" not in check_ls:
                periods = 1.0 / np.asarray(check_ls["frequency"], dtype=float)
                order = np.argsort(periods)
                ax_period.plot(
                    periods[order], np.asarray(check_ls["power"])[order],
                    color=_colors().TEXT_MUTED, lw=0.9, ls=":", label="check-star LS",
                )

        adopted = self._summary_period()
        if np.isfinite(adopted):
            ax_period.axvline(
                adopted, color=_colors().PLOT_FG, lw=1.2,
                label=f"adopted {adopted:.6f} d",
            )
        if self._show_alias_marks():
            for candidate in (self.alias_analysis or {}).get("candidates", [])[1:6]:
                candidate_period = float(candidate.get("period", np.nan))
                if np.isfinite(candidate_period):
                    ax_period.axvline(
                        candidate_period, color=_colors().WARN, lw=0.8, ls="-.", alpha=0.75
                    )
        min_period, max_period, _ = self._search_window()
        ax_period.set_xlim(min_period, max_period)
        if max_period / max(min_period, 1e-12) >= 30.0:
            ax_period.set_xscale("log")
        ax_period.set_title(f"Period search - {dtype}", loc="left")
        ax_period.set_xlabel("Trial period [days]")
        ax_period.set_ylabel("Statistic" if compact else "Periodogram statistic")
        ax_period.grid(True, alpha=0.25)
        if period_lines:
            ax_period.legend(
                loc="upper right", ncol=2 if compact else 1,
                fontsize=7 if compact else 8,
            )
        if period_notes and not compact:
            ax_period.text(
                0.02, 0.04,
                "\n".join(period_notes),
                transform=ax_period.transAxes, va="bottom",
                fontsize=8,
            )

        if np.isfinite(adopted) and adopted > 0:
            phase = ((time - t0) / adopted) % 1.0
            phase_grid = np.linspace(0.0, 2.0, 400)
            for label in filter_order:
                mask = filters == label
                for shift in (0.0, 1.0):
                    ax_phase.scatter(
                        phase[mask] + shift, mag[mask], s=16,
                        color=filter_colors[label], alpha=0.78,
                        edgecolors=_colors().PLOT_AXES_BG, linewidths=0.25,
                        label=(
                            label if compact and shift == 0
                            else f"{label} data" if shift == 0
                            else None
                        ),
                    )
                model = self._two_harmonic_phase_model(
                    phase[mask], mag[mask], mag_err[mask] if mag_err is not None else None,
                    phase_grid,
                )
                if model is not None:
                    ax_phase.plot(
                        phase_grid, model, color=filter_colors[label], lw=1.4,
                        label=None if compact else f"{label} fit",
                    )
            ax_phase.set_xlim(0.0, 2.0)
            ax_phase.invert_yaxis()
            ax_phase.set_title(
                f"Phase | P={adopted:.6f} d"
                if compact
                else f"Phase-folded at {adopted:.9f} d",
                loc="left",
            )
            ax_phase.set_xlabel("Phase")
            ax_phase.set_ylabel(
                "Diff. mag"
                if compact
                else "Median-centered differential magnitude [mag]"
            )
            ax_phase.grid(True, alpha=0.25)
            ax_phase.legend(
                loc="upper right", ncol=3 if compact else 2,
                fontsize=7 if compact else 8,
            )
            status = str((self.alias_analysis or {}).get("status", "")).upper()
            mode_status = str((self.multimode_diagnostic or {}).get("status", "")).upper()
            annotation = "\n".join(
                line for line in (f"alias: {status}" if status else "", f"mode: {mode_status}" if mode_status else "")
                if line
            )
            if annotation and not compact:
                ax_phase.text(0.02, 0.04, annotation, transform=ax_phase.transAxes, va="bottom", fontsize=8)
        else:
            ax_phase.text(0.5, 0.5, "No valid period", ha="center", va="center", transform=ax_phase.transAxes)

        target_id = int(self.lc_data.get("target_id", 0))
        cycles = baseline / adopted if np.isfinite(adopted) and adopted > 0 else np.nan
        cycle_text = f" | {cycles:.2f} cycles" if np.isfinite(cycles) else ""
        if compact:
            fig.suptitle(
                f"ID {target_id} | {baseline * 24.0:.2f} h{cycle_text}", fontsize=10
            )
        else:
            fig.suptitle(
                f"Target ID {target_id} - period analysis | "
                f"{baseline * 24.0:.2f} h{cycle_text}"
            )
        if compact:
            if stacked:
                fig.subplots_adjust(
                    left=0.16, right=0.98, bottom=0.09, top=0.90, hspace=0.80
                )
            else:
                fig.subplots_adjust(
                    left=0.10, right=0.98, bottom=0.14, top=0.87,
                    wspace=0.35, hspace=0.85,
                )
        self._summary_redraw()



def save_period_summary_figure(out_path, *, results, lc_data,
                               alias_analysis=None, check_star_data=None,
                               params=None, width_px: int = 1760,
                               search_window=None):
    """Draw the summary for one filter and write it to `out_path`.

    This is what the pipeline step calls. It is the same drawing the window
    performs — `PeriodSummaryPlotter._update_summary_plot` — with a fixed-size
    figure standing in for the canvas.
    """
    plotter = PeriodSummaryPlotter()
    plotter.results = results
    plotter.lc_data = lc_data
    plotter.alias_analysis = alias_analysis
    plotter._check_star_override = check_star_data
    plotter.params = params
    plotter._summary_canvas_width = lambda: int(width_px)
    if search_window is not None:
        plotter.search_min_period, plotter.search_max_period = search_window[:2]
        if len(search_window) > 2:
            plotter.search_samples_per_peak = int(search_window[2])

    plotter._update_summary_plot()
    fig = plotter._summary_figure()
    if not fig.axes:
        return None
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=160, bbox_inches="tight")
    return out_path
