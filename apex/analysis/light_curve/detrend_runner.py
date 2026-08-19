"""Detrending — taking out what the sky did so what the star did is left.

Three corrections live here: a per-night zero-point offset, a colour-dependent
term, and a global ensemble, with SYSREM as a fourth. Each takes the raw
differential light curve and returns a corrected one plus the fit that produced
it, and none of them needs a window.

They lived in one anyway, and this looked like the hard port because the
calculation appeared to read its inputs from widgets and write its results to
them. Half of that was true. The reads were already funnelled through a single
method — `_sync_state_from_controls` copied every spin box into a plain
attribute — so the calculation had been reading attributes all along, and only
four widget reads bypassed it. The writes are all presentation.

What a user supplies: the settings below (whose defaults are the window's own
starting values), `raw_df`, and the paths. What comes back is what the window
gets — `fit_and_apply` fills `corr_df`, `fit_params` and `summary_text`.
"""

from __future__ import annotations


from apex.analysis.light_curve.period_plot import _colors, _load_check_star_for_plot
from apex.analysis.light_curve.detrend_output_service import annotate_step10_output, build_detrend_summary_report_text, build_global_summary_text, write_step10_current_meta
from apex.analysis.light_curve.global_ensemble import solve_global_ensemble
from apex.analysis.light_curve.photometry_source_service import load_lightcurve_frame_photometry, resolve_lightcurve_photometry_source
from apex.utils.astro_utils import compute_airmass_from_jd_array, compute_airmass_from_header, compute_bjd_tdb_array, is_reasonable_airmass
from apex.utils.common_helpers import safe_float as _safe_float, normalize_filter_key as _normalize_filter_key, parse_jd as _parse_jd
from apex.utils.io_utils import (
    coerce_int64_source_id,
    read_csv_int64_source_id,
    load_headers_table as _load_headers_table_util,
    load_night_assignments as _load_night_assignments_util,
)
from apex.utils.photometry_provenance import build_photometry_provenance, format_photometry_provenance, summarize_photometry_table
from apex.utils.qc_utils import filter_frame_df_by_qc, load_frame_excludes, should_use_frame_quality_qc
from apex.utils.step_paths import forced_phot_input_dir
from apex.utils.step_paths_lc import step1_dir, step8_selection_dir, step9_lc_dir, step10_detrend_dir, step10_current_lc_path, step10_current_params_path, step10_current_summary_path, step10_current_plot_path, step10_current_meta_path, step10_current_global_zp_path, step10_current_global_mean_path, step10_current_global_diag_path, step10_history_dir, load_detrend_preference
from astropy.io import fits
from astropy.time import Time
from datetime import datetime
from pathlib import Path
from types import MappingProxyType
import json
import numpy as np
import pandas as pd
import re



def _load_headers_table(result_dir: Path) -> pd.DataFrame:
    return _load_headers_table_util(result_dir)


def _load_night_assignments_from_disk(result_dir: Path) -> dict[str, int]:
    """Load filename -> night_id from step1/night_assignments.json."""
    return _load_night_assignments_util(result_dir)


def _load_step9_comparison_ids_by_filter(result_dir: Path) -> dict[str, set[int]]:
    """Load authoritative comparison IDs for each Step 9 filter selection."""
    selection_dir = step8_selection_dir(result_dir)
    if not selection_dir.exists():
        return {}

    output: dict[str, set[int]] = {}
    for selection_path in sorted(selection_dir.glob("selection_*.json")):
        try:
            data = json.loads(selection_path.read_text(encoding="utf-8"))
        except Exception:
            continue
        filter_key = _normalize_filter_key(
            data.get("filter") or selection_path.stem.replace("selection_", "")
        )
        if not filter_key:
            continue
        comp_ids = {
            int(value)
            for value in data.get("comparison_ids", [])
            if value is not None
        }
        comp_source_ids = [
            int(value)
            for value in data.get("comparison_source_ids", [])
            if value is not None
        ]
        if comp_source_ids:
            source_to_id = _load_step8_source_to_id_map(result_dir, filter_key)
            resolved_ids = {
                int(source_to_id[source_id])
                for source_id in comp_source_ids
                if source_id in source_to_id
            }
            if resolved_ids:
                comp_ids = resolved_ids
        output[filter_key] = comp_ids
    return output


def _load_step9_selection_ids(result_dir: Path) -> tuple[int | None, list[int]]:
    """Load target/comp IDs from Step 9 selections, merging per-filter comp sets."""
    s9 = step8_selection_dir(result_dir)
    if not s9.exists():
        return None, []

    target_ids: set[int] = set()
    comp_ids: set[int] = set()

    for sp in sorted(s9.glob("selection_*.json")):
        try:
            data = json.loads(sp.read_text(encoding="utf-8"))
        except Exception:
            continue

        flt = sp.stem.replace("selection_", "")
        sid_map = _load_step8_source_to_id_map(result_dir, flt)

        tid = data.get("target_id")
        target_sid = data.get("target_source_id")
        if tid is None and target_sid is not None and int(target_sid) in sid_map:
            tid = int(sid_map[int(target_sid)])
        if tid is not None:
            target_ids.add(int(tid))

        cids = [int(x) for x in data.get("comparison_ids", []) if x is not None]
        comp_sids = [int(x) for x in data.get("comparison_source_ids", []) if x is not None]
        if not cids and comp_sids:
            cids = sorted({
                int(sid_map[int(sid)])
                for sid in comp_sids
                if int(sid) in sid_map
            })
        comp_ids.update(int(cid) for cid in cids)

    if len(target_ids) == 1:
        return next(iter(target_ids)), sorted(comp_ids)
    return None, sorted(comp_ids)


def _parse_color_expr(expr: str | None) -> tuple[str, str] | None:
    """Parse ``"B-V"``/``"g_r"`` into canonical bands (case preserved via normalize)."""
    if not expr:
        return None
    s = str(expr).strip().replace(" ", "")
    if "-" in s:
        parts = [p for p in s.split("-") if p]
    elif "_" in s:
        parts = [p for p in s.split("_") if p]
    else:
        return None
    if len(parts) != 2:
        return None
    a = _normalize_filter_key(parts[0])
    b = _normalize_filter_key(parts[1])
    if not a or not b:
        return None
    return a, b

def _load_step8_source_to_id_map(result_dir: Path, flt: str | None = None) -> dict[int, int]:
    """Load source_id -> final ID map from Step 9 outputs."""
    step9_out = step8_selection_dir(result_dir)
    if not step9_out.exists():
        return {}
    key = _normalize_filter_key(flt or "")
    candidates: list[tuple[Path, str]] = []
    if key:
        candidates.extend(
            [
                (step9_out / f"master_catalog_{key}.tsv", "\t"),
                (step9_out / f"id_mapping_{key}.csv", ","),
            ]
        )
    candidates.extend(
        [(p, "\t") for p in sorted(step9_out.glob("master_catalog_*.tsv"))]
    )
    candidates.extend(
        [(p, ",") for p in sorted(step9_out.glob("id_mapping_*.csv"))]
    )

    mapping: dict[int, int] = {}
    for path, sep in candidates:
        if not path.exists():
            continue
        try:
            df = read_csv_int64_source_id(path, sep=sep)
        except Exception:
            continue
        if "det_uid" in df.columns and "ID" not in df.columns:
            df = df.rename(columns={"det_uid": "ID"})
        if not {"source_id", "ID"} <= set(df.columns):
            continue
        sid_vals = coerce_int64_source_id(df["source_id"])
        id_vals = pd.to_numeric(df["ID"], errors="coerce").astype("Int64")
        for sid_val, id_val in zip(sid_vals, id_vals):
            if pd.isna(sid_val) or pd.isna(id_val):
                continue
            sid_int = int(sid_val)
            if sid_int not in mapping:
                mapping[sid_int] = int(id_val)
    return mapping

# Shared read-only stand-ins for state the window builds in its constructor.
# Read-only on purpose: an in-place write would otherwise reach every instance,
# and the mapping proxy turns that mistake into an error instead of a haunting.
_EMPTY_FRAME = pd.DataFrame()
_EMPTY_MAP = MappingProxyType({})

class _Silent:
    """Stands in for a widget the batch path has no equivalent of.

    Every write the calculation makes to a label, a table or a progress bar is
    presentation. Rather than edit fifty call sites — and risk changing the
    calculation while moving it, which is how `_build_ensemble_series` got
    renamed out from under its caller — the runner hands those names something
    that accepts any call and does nothing. The window's real widgets are
    instance attributes and shadow these.

    Reads are deliberately not covered here: the four values the calculation
    needs from the user come through the accessors below, so a silent stand-in
    can never feed the science a wrong number.
    """

    def __getattr__(self, _name):
        return lambda *args, **kwargs: None

    def __bool__(self):
        return False


class DetrendRunner:
    """The detrend calculation, minus the window.

    Settings are class attributes so a batch runner sets the handful it cares
    about and inherits the rest. The defaults are the window's own starting
    values — what its spin boxes read before anyone touches them.
    """

    # -- which correction, and the clipping every one of them shares ---------
    mode: str = "offset"                  # offset | color | global | sysrem
    sigma_clip: bool = True
    clip_sigma: float = 3.0
    clip_iters: int = 3

    # -- phase folding: the display axis, so it reaches the saved plot only --
    phase_period: float = 0.0
    phase_t0: float = 0.0
    phase_cycles: float = 2.0
    x_axis_mode: str = "time"

    # -- colour mode ---------------------------------------------------------
    color_by: str = ""

    # -- global ensemble -----------------------------------------------------
    global_min_comps: int = 3
    global_sigma: float = 3.0
    global_iters: int = 3
    global_rms_pct: float = 80.0
    global_rms_threshold: float = 0.05
    global_frame_sigma: float = 5.0
    global_gauge: str = ""
    global_robust: bool = True
    global_interp_missing: bool = False
    global_normalize: bool = False
    global_rescale_errors: bool = False
    use_global_k2: bool = False

    # -- SYSREM --------------------------------------------------------------
    sysrem_iter: int = 3
    sysrem_apply: int = 1

    # -- the run itself ------------------------------------------------------
    target_id: int = 0
    filter_selection: str = "All"
    selected_dates = None                 # None = every night present
    runtime_mode: bool = False

    # -- presentation the batch path has no equivalent of --------------------
    log_text = _Silent()
    recommendation_label = _Silent()
    photometry_source_label = _Silent()
    busy_status_label = _Silent()
    busy_progress_bar = _Silent()
    analysis_text = _Silent()
    color_status_label = _Silent()
    result_table = _Silent()
    id_info_label = _Silent()

    def _update_results_table(self, *args, **kwargs) -> None:
        """The window fills a table; a batch run has the DataFrame."""

    def _populate_date_list(self, *args, **kwargs) -> None:
        """The window fills its night list; a batch run was told the nights."""

    def _refresh_filter_combo(self, *args, **kwargs) -> None:
        """The window fills its filter box; a batch run was told the filter."""

    def _update_comp_label(self, *args, **kwargs) -> None:
        """Comparison-star caption."""

    def _tell_user(self, level: str, title: str, message: str) -> None:
        """Say something the window would have put in a dialog.

        Every call site already sits behind `if not silent:` or `if update_ui:`,
        so a batch run reaches none of them today — but "today" is the kind of
        thing that changes, and a calculation that reaches for `QMessageBox`
        cannot run where Qt is not installed. So the runner owns the channel and
        the window subscribes, the same trade `PsfPhotometryRunner` made when it
        stopped calling `.emit()`.

        The default writes to the log rather than raising, because the caller
        already decided what happens next — these sites return False or raise on
        their own branch, and turning a notice into an exception here would
        change control flow while moving code.
        """
        self.log("[%s] %s: %s" % (level.upper(), title, message))

    def log(self, message: str) -> None:
        print(message)

    # -- what one method fills and another reads ------------------------------
    #
    # The window builds these in its constructor, so the calculation could
    # assume they exist. A runner built any other way cannot, and the failure
    # is an AttributeError partway through a load. `None` rather than an empty
    # frame or list: a class-level mutable would be shared by every instance,
    # and the code already treats absent as absent.

    # Each default matches how the attribute is consumed, which is why they are
    # not all `None`: `comp_*` are iterated, `delta_c_map` is `.items()`-ed,
    # `global_mean_df` is asked for `.empty`. Every one of them is replaced
    # wholesale rather than mutated in place — checked — so sharing an empty
    # value across instances cannot leak between runs.
    raw_df = _EMPTY_FRAME
    corrected_df = _EMPTY_FRAME
    params_df = _EMPTY_FRAME
    global_mean_df = _EMPTY_FRAME
    summary_text = ""
    comp_active_ids = ()
    comp_candidate_ids = ()
    delta_c_map = _EMPTY_MAP
    global_diagnostics = None

    # -- which panels the saved figure shows ---------------------------------
    #
    # The window defaults to `corr`, one panel, because its canvas shares a
    # screen with the controls. A file has no such constraint, and the two
    # panels the single view hides are the ones that make the figure worth
    # keeping: the raw curve it was corrected from, and Delta-mag against
    # airmass, which is where a bad night shows itself. So a batch run defaults
    # to `all` and the window keeps its own.
    _plot_view_mode: str = "all"

    def _apply_plot_view(self, redraw: bool = True) -> None:
        """Re-lay the raw/corr/diag axes for the selected view.

        The three axes objects stay alive (50+ call sites draw into them);
        only their subplotspec/visibility changes. In single view the hidden
        axes share the visible one's spec so tight_layout math stays sane.
        """
        from matplotlib.gridspec import GridSpec
        fig = self._plot_figure()
        order = ("raw", "corr", "diag")
        axmap = {"raw": self.ax_raw, "corr": self.ax_corr, "diag": self.ax_diag}
        mode = getattr(self, "_plot_view_mode", "corr")
        if mode == "all":
            gs = GridSpec(3, 1, figure=fig)
            for i, key in enumerate(order):
                axmap[key].set_subplotspec(gs[i])
                axmap[key].set_visible(True)
        else:
            gs = GridSpec(1, 1, figure=fig)
            for key in order:
                axmap[key].set_subplotspec(gs[0])
                axmap[key].set_visible(key == mode)
        if redraw:
            try:
                fig.tight_layout()
            except Exception:
                pass
            self._plot_redraw()

    # -- three axes the window builds on its canvas --------------------------

    @property
    def ax_raw(self):
        return self._axes()[0]

    @property
    def ax_corr(self):
        return self._axes()[1]

    @property
    def ax_diag(self):
        return self._axes()[2]

    def _axes(self):
        """Raw / corrected / diagnostic, stacked — the window's 311-312-313."""
        if getattr(self, "_headless_axes", None) is None:
            fig = self._plot_figure()
            self._headless_axes = (fig.add_subplot(311), fig.add_subplot(312),
                                   fig.add_subplot(313))
        return self._headless_axes

    # -- presentation the batch path has no equivalent of --------------------

    def _refresh_style(self, _widget=None) -> None:
        """Re-apply the Qt stylesheet after a property change. Nothing to do here."""

    def _queue_background_fit(self, *args, **kwargs) -> None:
        """The window runs the fit on a worker thread; a batch run is the thread."""

    def _rebuild_color_map_controls(self, *args, **kwargs) -> None:
        """Colour-mapping combo boxes."""

    def _sync_mode_controls_from_state(self, *args, **kwargs) -> None:
        """Push `self.mode` back onto the radio buttons."""

    def _set_busy_state(self, *args, **kwargs) -> None:
        """Busy cursor and progress bar."""

    def _set_busy_message(self, *args, **kwargs) -> None:
        """Busy caption."""

    # -- remaining widgets, write-only from here ------------------------------
    mode_offset = _Silent()
    mode_color = _Silent()
    mode_global = _Silent()
    mode_sysrem = _Silent()
    target_edit = _Silent()
    plot_canvas = _Silent()

    # -- state a window would have been constructed with ----------------------
    params = None
    project_state = None
    file_manager = None
    main_window = None
    datasets = ()
    color_map_by_filter = None

    # -- the canvas: the window has one, a batch run makes one ---------------
    #
    # This is where the four detrend figures come from. They were never drawn
    # headless because the only figure in reach belonged to a `FigureCanvas`,
    # and a batch run has no canvas — not because the drawing needed Qt. It
    # does not: `_update_plots` and `_plot_global_diagnostics` are matplotlib
    # and nothing else.

    plot_width_in: float = 11.0
    plot_height_in: float = 7.0

    def _plot_figure(self):
        if getattr(self, "_headless_fig", None) is None:
            from matplotlib.figure import Figure
            self._headless_fig = Figure(
                figsize=(self.plot_width_in, self.plot_height_in), dpi=150)
        return self._headless_fig

    def _ensure_plot_drawn(self) -> None:
        """Draw the figure about to be saved.

        A batch run has had no UI pass to draw it. Overridden to nothing in the
        window, whose canvas is already current.
        """
        self._apply_plot_view(redraw=False)
        self._update_plots()

    def _plot_redraw(self) -> None:
        """The window schedules a repaint; a figure about to be saved needs none."""

    # -- the four inputs the window reads off widgets ------------------------

    def _target_id_text(self) -> str:
        return str(self.target_id or "").strip()

    def _filter_selection(self) -> str:
        return str(self.filter_selection or "All")

    def _use_global_k2(self) -> bool:
        return bool(self.use_global_k2)

    def _sync_state_from_controls(self) -> None:
        """The window copies its spin boxes into the attributes above.

        A batch run set them directly, so there is nothing to copy. This is the
        seam that already existed: the calculation reads attributes, not
        widgets, and always did.
        """

    def _append_background_log(self, msg: str) -> None:
        try:
            self.log_text.append(str(msg))
        except RuntimeError:
            pass

    def _flush_busy_log_buffer(self) -> None:
        if not self._busy_log_buffer:
            return
        self.log_text.append("\n".join(self._busy_log_buffer))
        self._busy_log_buffer = []

    def _update_analysis_panel(self):
        """Update analysis panel with data statistics and recommendations."""
        if self.raw_df.empty:
            self.analysis_text.setText("Step 10 결과가 있으면 자동 로드되어 분석 결과가 표시됩니다.")
            self.recommendation_label.setText("")
            self.recommendation_label.setVisible(False)
            return

        if self.mode == "global":
            n_points = len(self.raw_df)
            n_dates = self.raw_df["date"].nunique() if "date" in self.raw_df.columns else 0
            filters = sorted({_normalize_filter_key(f) for f in self.raw_df.get("filter", []) if str(f).strip()})
            self.analysis_text.setText(
                f"<b>Global Ensemble:</b> {n_points}점, {n_dates}일, 필터: {', '.join(filters) or 'N/A'}"
            )
            self.recommendation_label.setText("")
            self.recommendation_label.setVisible(False)
            return
        if self.mode == "sysrem":
            n_points = len(self.raw_df)
            n_dates = self.raw_df["date"].nunique() if "date" in self.raw_df.columns else 0
            filters = sorted({_normalize_filter_key(f) for f in self.raw_df.get("filter", []) if str(f).strip()})
            self.analysis_text.setText(
                f"<b>SYSREM:</b> {n_points}점, {n_dates}일, 필터: {', '.join(filters) or 'N/A'}"
            )
            self.recommendation_label.setText("")
            self.recommendation_label.setVisible(False)
            return

        lines = []

        # 1. Data summary
        n_points = len(self.raw_df)
        n_dates = self.raw_df["date"].nunique() if "date" in self.raw_df.columns else 0
        filters = sorted({_normalize_filter_key(f) for f in self.raw_df.get("filter", []) if str(f).strip()})
        lines.append(f"<b>데이터:</b> {n_points}점, {n_dates}일, 필터: {', '.join(filters) or 'N/A'}")

        # 2. Airmass range
        airmass = self.raw_df["airmass"].to_numpy(float)
        airmass = airmass[np.isfinite(airmass)]
        if airmass.size > 0:
            x_min, x_max = np.min(airmass), np.max(airmass)
            x_range = x_max - x_min
            lines.append(f"<b>Airmass:</b> {x_min:.3f} ~ {x_max:.3f} (ΔX = {x_range:.3f})")
        else:
            x_range = 0.0
            lines.append("<b>Airmass:</b> 데이터 없음")

        # 3. Color index difference (ΔC)
        delta_c_values = []
        delta_c_info = []
        for fkey, dc in self.delta_c_map.items():
            if np.isfinite(dc):
                delta_c_values.append(abs(dc))
                delta_c_info.append(f"{fkey or 'all'}: {dc:+.3f}")

        if delta_c_values:
            max_dc = max(delta_c_values)
            lines.append(f"<b>|ΔC| (Target-Comp):</b> {', '.join(delta_c_info)}")
        else:
            max_dc = 0.0
            lines.append("<b>|ΔC|:</b> 계산 불가 (color index 설정 필요)")

        # 4. Expected color term effect
        if delta_c_values and airmass.size > 0:
            # Typical k'' ~ 0.02-0.05 mag/mag for broad-band filters
            k2_typical = 0.03
            effect_max = k2_typical * max_dc * x_range
            lines.append(f"<b>예상 색항 효과:</b> ~{effect_max:.4f} mag (k''≈0.03 가정)")

        # 5. Raw scatter
        raw_mag = self.raw_df["diff_mag_raw"].to_numpy(float)
        raw_mag = raw_mag[np.isfinite(raw_mag)]
        if raw_mag.size > 0:
            raw_rms = np.std(raw_mag)
            lines.append(f"<b>Raw RMS:</b> {raw_rms:.4f} mag")

        self.analysis_text.setText("<br>".join(lines))

        # Generate recommendation
        recommendation = self._generate_recommendation(max_dc, x_range, n_dates)
        if recommendation:
            self.recommendation_label.setText(recommendation)
            self.recommendation_label.setVisible(True)
        else:
            self.recommendation_label.setVisible(False)

    def _generate_recommendation(self, max_dc: float, x_range: float, n_dates: int) -> str:
        """Generate mode recommendation based on data characteristics."""
        reasons = []

        # Decision logic based on astronomical data science principles
        use_color_mode = False
        color_mode_warning = False

        # 1. Check airmass range FIRST - critical for k'' fitting stability
        airmass_sufficient = x_range >= 0.3

        # 2. Check color index difference
        if max_dc >= 0.5:
            reasons.append(f"|ΔC| = {max_dc:.2f} ≥ 0.5: 색차가 매우 큼")
            if airmass_sufficient:
                use_color_mode = True
            else:
                color_mode_warning = True
                reasons.append(f"  → 단, ΔX = {x_range:.2f} < 0.3: k'' 피팅 불안정 우려")
        elif max_dc >= 0.3:
            reasons.append(f"|ΔC| = {max_dc:.2f} ≥ 0.3: 색차가 상당함")
            if airmass_sufficient:
                use_color_mode = True
            else:
                color_mode_warning = True
        elif max_dc > 0:
            reasons.append(f"|ΔC| = {max_dc:.2f} < 0.3: 색차가 작음")

        # 3. Check airmass range details
        if x_range >= 0.5:
            reasons.append(f"ΔX = {x_range:.2f} ≥ 0.5: airmass 범위 충분")
        elif x_range >= 0.3:
            reasons.append(f"ΔX = {x_range:.2f}: airmass 범위 적절")
        else:
            reasons.append(f"ΔX = {x_range:.2f} < 0.3: airmass 범위 좁음 (k'' 피팅 불안정)")
            use_color_mode = False  # Override - not enough airmass range

        # 4. Multi-night consideration
        if n_dates > 1:
            reasons.append(f"{n_dates}일 데이터: 밤별 ZP₀ 보정 필수")

        # Build recommendation message
        # The two "권장" cases used to differ only by orange vs amber, which the
        # wording already carries; both map to the theme's warn banner.
        if use_color_mode:
            mode_text = "⚠️ <b>Color 모드 권장</b>"
            banner = "warn"
        elif color_mode_warning:
            mode_text = "⚠️ <b>Offset 모드 권장</b> (색차 있으나 ΔX 부족)"
            banner = "warn"
        else:
            mode_text = "✓ <b>Offset 모드 적합</b>"
            banner = "ok"

        self.recommendation_label.setProperty("banner", banner)
        self._refresh_style(self.recommendation_label)

        return f"{mode_text}<br>{'<br>'.join('• ' + r for r in reasons)}"

    def _parse_target_id_from_name(self, name: str) -> int | None:
        m = re.search(r"lightcurve_ID(\d+)_raw\.csv", name)
        if not m:
            return None
        try:
            return int(m.group(1))
        except Exception:
            return None

    def _apply_frame_excludes(self, df: pd.DataFrame, result_dir: Path, label: str) -> pd.DataFrame:
        if df.empty or "file" not in df.columns:
            return df
        exclude_map = load_frame_excludes(result_dir, exclude_dir=step9_lc_dir(result_dir))
        if not exclude_map:
            return df
        before = len(df)
        df = df[~df["file"].astype(str).isin(exclude_map.keys())]
        removed = before - len(df)
        if removed > 0:
            self.log(f"[Frame QC] {label}: {before} → {len(df)} (excluded {removed})")
        return df

    def _raw_lightcurve_candidates(self, result_dir: Path, target_id: int) -> list[Path]:
        step10_path = step9_lc_dir(result_dir)
        return [
            step10_path / f"lightcurve_ID{target_id}_raw.csv",
            step10_path / f"lightcurve_combined_ID{target_id}_raw.csv",
        ]

    def _rebuild_step10_raw_outputs(self, target_id: int) -> bool:
        if not self.datasets:
            return False

        self._load_comp_selection()
        comp_ids = self.comp_active_ids or self.comp_candidate_ids
        comp_ids = [int(x) for x in comp_ids if str(x).strip()]
        if not comp_ids:
            self.log("[LOAD] Step 10 raw rebuild skipped: comparison IDs not found")
            return False

        self.log(f"[LOAD] Step 10 raw CSV missing. Rebuilding for target ID {target_id}...")
        try:
            from .step9_lightcurve_builder import LightCurveBuilderWindow

            builder = LightCurveBuilderWindow(
                self.params,
                self.file_manager,
                self.project_state,
                self.main_window,
            )
            builder.hide()
            builder.show_log_window = lambda: None
            builder.datasets = [(label, Path(path)) for label, path in self.datasets]
            builder.target_edit.setText(str(target_id))
            builder.comp_edit.setText(",".join(str(i) for i in comp_ids))
            builder.comp_candidate_ids = list(comp_ids)
            builder._update_comp_ids_from_input()
            builder.build_light_curve()
            builder.close()
            builder.deleteLater()
        except Exception as e:
            self.log(f"[LOAD] Step 10 raw rebuild failed: {e}")
            return False

        for _label, result_dir in self.datasets:
            result_dir = Path(result_dir)
            if any(path.exists() for path in self._raw_lightcurve_candidates(result_dir, target_id)):
                return True
        return False

    def load_raw_data(
        self,
        silent: bool = False,
        *,
        rebuild_missing: bool = True,
        allow_fits_airmass: bool = True,
    ) -> bool:
        if not self.datasets:
            if not silent:
                self._tell_user("info", "Detrend", "데이터셋이 없습니다.")
            return False
        target_text = self._target_id_text()
        if not target_text:
            base_dir = step9_lc_dir(self.params.P.result_dir)
            raw_paths = list(base_dir.glob("lightcurve_ID*_raw.csv"))
            if not raw_paths and self.datasets:
                raw_paths = list(step9_lc_dir(self.datasets[0][1]).glob("lightcurve_ID*_raw.csv"))
            if raw_paths:
                target_id = self._parse_target_id_from_name(raw_paths[0].name)
                if target_id is not None:
                    self.target_edit.setText(str(target_id))
                    target_text = str(target_id)
                    self.log(f"[LOAD] Target ID from raw filename: {target_id}")
        if not target_text:
            if not silent:
                self._tell_user("info", "Detrend", "대상 ID가 필요합니다.")
            return False
        target_id = int(target_text)

        raw_frames = []
        def _collect_raw_frames() -> list[pd.DataFrame]:
            frames: list[pd.DataFrame] = []
            for label, result_dir in self.datasets:
                result_dir = Path(result_dir)
                step10_path = step9_lc_dir(result_dir)
                raw_path = next((cand for cand in self._raw_lightcurve_candidates(result_dir, target_id) if cand.exists()), None)

                if raw_path is None:
                    self.log(f"[LOAD] Missing raw in {step10_path}")
                    self.log(f"  Tried: lightcurve_ID{target_id}_raw.csv, lightcurve_combined_ID{target_id}_raw.csv")
                    if step10_path.exists():
                        available = list(step10_path.glob("lightcurve_*.csv"))
                        if available:
                            self.log(f"  Available: {[f.name for f in available[:5]]}")
                    continue

                try:
                    df = pd.read_csv(raw_path)
                    self.log(f"[LOAD] Loaded: {raw_path.name} ({len(df)} rows)")
                except Exception as e:
                    self.log(f"[LOAD] Failed to read {raw_path.name}: {e}")
                    continue
                df = df.copy()
                df["dataset"] = label
                df = self._apply_frame_excludes(df, result_dir, str(label))
                frames.append(df)
            return frames

        raw_frames = _collect_raw_frames()
        if (
            not raw_frames
            and rebuild_missing
            and self._rebuild_step10_raw_outputs(target_id)
        ):
            raw_frames = _collect_raw_frames()

        if not raw_frames:
            step10_path = step9_lc_dir(self.params.P.result_dir)
            msg = f"Raw 데이터를 찾지 못했습니다.\n\n경로: {step10_path}\n\nStep 10에서 먼저 'Build Light Curve'를 실행하세요."
            self.log(f"[LOAD] {msg.replace(chr(10), ' ')}")
            if not silent:
                self._tell_user("info", "Detrend", msg)
            return False

        self.raw_df = pd.concat(raw_frames, ignore_index=True)
        self.photometry_source_label.setText(
            format_photometry_provenance(
                summarize_photometry_table(self.raw_df)
            )
        )
        if "diff_mag_raw" not in self.raw_df.columns and "diff_mag" in self.raw_df.columns:
            self.raw_df["diff_mag_raw"] = self.raw_df["diff_mag"]
        if "diff_mag_raw" not in self.raw_df.columns:
            self.raw_df["diff_mag_raw"] = np.nan
        if "diff_err" not in self.raw_df.columns:
            self.raw_df["diff_err"] = np.nan
        self.raw_df["diff_mag_raw"] = pd.to_numeric(self.raw_df["diff_mag_raw"], errors="coerce")
        self.raw_df["airmass"] = pd.to_numeric(self.raw_df.get("airmass", np.nan), errors="coerce")
        self.raw_df["JD"] = pd.to_numeric(self.raw_df.get("JD", np.nan), errors="coerce")
        # BJD_TDB가 있으면 JD 컬럼을 덮어쓰기 (위상/시간 분석 모두 BJD_TDB 기준)
        if "BJD_TDB" in self.raw_df.columns:
            bjd = pd.to_numeric(self.raw_df["BJD_TDB"], errors="coerce")
            valid_bjd = bjd.notna() & np.isfinite(bjd)
            if valid_bjd.any():
                self.raw_df.loc[valid_bjd, "JD"] = bjd[valid_bjd]
        if "date" not in self.raw_df.columns:
            self.raw_df["date"] = "unknown"
        self._fill_date_from_jd()
        self._fill_night_id()
        self._fill_airmass_from_headers(
            allow_fits_fallback=allow_fits_airmass
        )

        self._load_comp_selection()
        self._populate_date_list()
        self._refresh_filter_combo(self.raw_df.get("filter", pd.Series([], dtype=str)).astype(str).tolist())
        self._rebuild_color_map_controls(self.raw_df.get("filter", pd.Series([], dtype=str)).astype(str).tolist())
        self.log(f"[LOAD] Raw points: {len(self.raw_df)}")

        self._refresh_delta_c_map()
        self._log_color_index_info()
        self._update_color_mode_enabled()
        self.corrected_df = pd.DataFrame()
        self.params_df = pd.DataFrame()
        self.global_mean_df = pd.DataFrame()
        self.global_diagnostics = {}
        self._load_saved_current_result(silent=True)
        if not self.runtime_mode:
            self._update_results_table()
            self._update_plots()
            self._update_analysis_panel()
        return True

    def _load_saved_current_result(self, silent: bool = True) -> bool:
        target_text = self._target_id_text()
        if not target_text:
            return False
        try:
            target_id = int(target_text)
        except Exception:
            return False

        result_root = self._find_saved_current_result_root(target_id)
        if result_root is None:
            return False
        lc_path = step10_current_lc_path(result_root, target_id)

        try:
            corrected_df = pd.read_csv(lc_path)
        except Exception as e:
            self.log(f"[LOAD] Failed to read current Step11 result: {e}")
            if not silent:
                self._tell_user("warn", "Detrend", f"Current 결과 로드 실패:\n{e}")
            return False

        if corrected_df.empty:
            self.log(f"[LOAD] Current Step11 result is empty: {lc_path.name}")
            return False

        mode = ""
        if "correction_mode" in corrected_df.columns:
            vals = corrected_df["correction_mode"].dropna().astype(str).str.strip().str.lower()
            if not vals.empty:
                mode = vals.iloc[0]
        if mode not in ("global", "color", "offset", "sysrem"):
            try:
                mode = str(load_detrend_preference(result_root, target_id=target_id) or "").strip().lower()
            except Exception:
                mode = ""
        if mode not in ("global", "color", "offset", "sysrem"):
            mode = self.mode if self.mode in ("global", "color", "offset", "sysrem") else "offset"

        for col in ("JD", "jd", "BJD_TDB", "airmass", "diff_mag_raw", "diff_mag_corr", "diff_err", "diff_err_corr", "residual", "fit_value"):
            if col in corrected_df.columns:
                corrected_df[col] = pd.to_numeric(corrected_df[col], errors="coerce")
        if "JD" not in corrected_df.columns and "jd" in corrected_df.columns:
            corrected_df["JD"] = corrected_df["jd"]
        if "date" not in corrected_df.columns and not self.raw_df.empty and "file" in corrected_df.columns and "file" in self.raw_df.columns:
            date_map = self.raw_df.groupby(self.raw_df["file"].astype(str))["date"].first()
            corrected_df["date"] = corrected_df["file"].astype(str).map(date_map).fillna("unknown")
        if "night_id" not in corrected_df.columns and not self.raw_df.empty and "file" in corrected_df.columns and "file" in self.raw_df.columns:
            night_map = self.raw_df.groupby(self.raw_df["file"].astype(str))["night_id"].first() if "night_id" in self.raw_df.columns else None
            if night_map is not None:
                corrected_df["night_id"] = corrected_df["file"].astype(str).map(night_map)

        params_df = pd.DataFrame()
        global_mean_df = pd.DataFrame()
        global_diagnostics: dict = {}
        params_path = step10_current_params_path(result_root, target_id)

        if mode == "global":
            zp_path = step10_current_global_zp_path(result_root, target_id)
            params_candidate = zp_path if zp_path.exists() else params_path
            if params_candidate.exists():
                try:
                    params_df = pd.read_csv(params_candidate)
                except Exception as e:
                    self.log(f"[LOAD] Failed to read global ZP/current params: {e}")
            mean_path = step10_current_global_mean_path(result_root, target_id)
            if mean_path.exists():
                try:
                    global_mean_df = pd.read_csv(mean_path)
                except Exception as e:
                    self.log(f"[LOAD] Failed to read global mean: {e}")
            diag_path = step10_current_global_diag_path(result_root, target_id)
            if diag_path.exists():
                try:
                    global_diagnostics = json.loads(diag_path.read_text(encoding="utf-8"))
                except Exception as e:
                    self.log(f"[LOAD] Failed to read global diagnostics: {e}")
        elif params_path.exists():
            try:
                params_df = pd.read_csv(params_path)
            except Exception as e:
                self.log(f"[LOAD] Failed to read current params: {e}")

        self.corrected_df = corrected_df
        self.params_df = params_df
        self.global_mean_df = global_mean_df
        self.global_diagnostics = global_diagnostics
        self.mode = mode
        self._sync_mode_controls_from_state()
        self.log(
            f"[LOAD] Restored Step11 current result: {lc_path.name} "
            f"(mode={mode}, params={len(self.params_df)})"
        )
        return True

    def _saved_current_result_roots(self) -> list[Path]:
        roots: list[Path] = []
        seen: set[str] = set()
        for root in [Path(self.params.P.result_dir), *[Path(path) for _, path in self.datasets]]:
            try:
                key = str(root.resolve())
            except Exception:
                key = str(root)
            if key in seen:
                continue
            seen.add(key)
            roots.append(root)
        return roots

    def _find_saved_current_result_root(self, target_id: int) -> Path | None:
        for root in self._saved_current_result_roots():
            if step10_current_lc_path(root, target_id).exists():
                return root
        return None

    def _load_global_ensemble_df(
        self,
        target_id_override: int | None = None,
        comp_ids_override: list[int] | None = None,
    ) -> pd.DataFrame:
        if not self.datasets:
            raise RuntimeError("No datasets available")

        if target_id_override is not None:
            target_id = int(target_id_override)
        else:
            target_text = self._target_id_text()
            if not target_text:
                raise RuntimeError("Target ID is required")
            target_id = int(target_text)

        if comp_ids_override is not None:
            comp_ids = [int(c) for c in comp_ids_override if str(c).strip() and int(c) != target_id]
        else:
            if not self.comp_active_ids and not self.comp_candidate_ids:
                self._load_comp_selection()
            comp_ids = self.comp_active_ids or self.comp_candidate_ids
            comp_ids = [int(c) for c in comp_ids if str(c).strip() and int(c) != target_id]
        if not comp_ids:
            raise RuntimeError("Comparison IDs not found")

        current_root = Path(self.params.P.result_dir).resolve()
        source_by_root: dict[str, dict] = {}
        for _, dataset_root in self.datasets:
            resolved = Path(dataset_root).resolve()
            state = self.project_state if resolved == current_root else None
            source_by_root[str(resolved)] = resolve_lightcurve_photometry_source(
                resolved, state
            )
        use_psf = bool(source_by_root) and all(
            info["source"] == "psf" for info in source_by_root.values()
        )
        if not use_psf:
            for root_key in source_by_root:
                root = Path(root_key)
                aperture_dir = forced_phot_input_dir(root)
                source_by_root[root_key] = {
                    **build_photometry_provenance("aperture", "mag", "mag_err"),
                    "directory": aperture_dir,
                    "index_path": aperture_dir / "photometry_index.csv",
                    "reason": "Using one aperture source across all datasets",
                }
        self.active_photometry_source = "psf" if use_psf else "aperture"
        if source_by_root:
            first_source = next(iter(source_by_root.values()))
            self.photometry_source_label.setText(
                format_photometry_provenance(first_source)
            )
            self.photometry_source_label.setToolTip(
                str(first_source.get("reason", ""))
            )
        self.log(
            f"[GLOBAL] Photometry source: {self.active_photometry_source.upper()}"
        )

        rows = []
        for label, result_dir in self.datasets:
            result_dir = Path(result_dir)
            source_info = source_by_root[str(result_dir.resolve())]
            filter_comp_ids = _load_step9_comparison_ids_by_filter(result_dir)
            idx_path = Path(source_info["index_path"])
            if not idx_path.exists():
                self.log(f"[GLOBAL] photometry_index.csv missing: {result_dir}")
                continue

            try:
                idx = pd.read_csv(idx_path)
            except Exception as e:
                self.log(f"[GLOBAL] Failed to read index: {e}")
                continue
            if "file" not in idx.columns:
                self.log("[GLOBAL] photometry_index.csv missing 'file'")
                continue
            if source_info["source"] == "psf":
                allowed = set(source_info.get("frames", []))
                idx = idx[
                    idx["file"].astype(str).map(lambda value: Path(value).name).isin(allowed)
                ].reset_index(drop=True)

            # Apply step4 QC pass filter (frame_quality.csv) if the parameter is enabled
            use_qc = should_use_frame_quality_qc(
                result_dir,
                self.params.P,
                "phot_use_qc_pass_only",
                default=False,
            )
            if use_qc:
                idx, qc_info = filter_frame_df_by_qc(result_dir, idx, file_col="file", require_qc=True)
                if qc_info.get("applied"):
                    self.log(f"[GLOBAL] Frame QC: {qc_info['kept']}/{qc_info['total']} frames kept")

            # Apply manual frame exclusions (D-key in step9)
            # LC excludes live in lc_lightcurve/; falls back to result_dir/ for compat.
            exclude_map = load_frame_excludes(result_dir, exclude_dir=step9_lc_dir(result_dir))
            if exclude_map:
                before = len(idx)
                idx = idx[~idx["file"].astype(str).isin(exclude_map.keys())].reset_index(drop=True)
                if before != len(idx):
                    self.log(f"[GLOBAL] Frame excludes: removed {before - len(idx)} frame(s)")

            headers_df = _load_headers_table(result_dir)
            jd_map = {}
            filt_map = {}
            if not headers_df.empty and "Filename" in headers_df.columns:
                if "JD" in headers_df.columns:
                    jd_map = dict(
                        zip(
                            headers_df["Filename"].astype(str),
                            pd.to_numeric(headers_df["JD"], errors="coerce"),
                        )
                    )
                elif "DATE-OBS" in headers_df.columns:
                    jd_map = dict(
                        zip(
                            headers_df["Filename"].astype(str),
                            headers_df["DATE-OBS"].astype(str),
                        )
                    )
                for col in ("FILTER", "filter"):
                    if col in headers_df.columns:
                        filt_map = dict(zip(headers_df["Filename"].astype(str), headers_df[col].astype(str)))
                        break

            for _, row in idx.iterrows():
                fname = str(row.get("file", "")).strip()
                if not fname:
                    continue
                phot = load_lightcurve_frame_photometry(
                    result_dir,
                    fname,
                    source_info,
                    str(row.get("filter", "") or ""),
                )
                if phot is None or phot.empty:
                    continue

                id_col = None
                for cand in ("ID", "id", "det_uid"):
                    if cand in phot.columns:
                        id_col = cand
                        break
                if not id_col:
                    continue
                if id_col == "det_uid":
                    phot = phot.rename(columns={"det_uid": "ID"})
                    id_col = "ID"

                filt_val = ""
                if "FILTER" in phot.columns:
                    filt_val = str(phot["FILTER"].iloc[0]).strip()
                if not filt_val:
                    filt_val = str(row.get("filter", "") or "")
                if not filt_val and fname in filt_map:
                    filt_val = str(filt_map.get(fname, "") or "")
                filt_key = _normalize_filter_key(filt_val)

                phot[id_col] = pd.to_numeric(phot[id_col], errors="coerce").astype("Int64")
                frame_comp_ids = comp_ids
                if filt_key in filter_comp_ids:
                    selected_ids = filter_comp_ids[filt_key]
                    frame_comp_ids = [
                        comp_id for comp_id in comp_ids if comp_id in selected_ids
                    ]
                wanted = [target_id] + frame_comp_ids
                phot = phot[phot[id_col].isin(wanted)].copy()
                if phot.empty:
                    continue

                date_obs = jd_map.get(fname)
                jd_val = _safe_float(date_obs)
                if not np.isfinite(jd_val):
                    jd_val = _parse_jd(date_obs) if date_obs else np.nan

                for _, r in phot.iterrows():
                    sid = int(r[id_col]) if pd.notna(r[id_col]) else None
                    if sid is None:
                        continue
                    mag = _safe_float(r.get("mag"))
                    err = _safe_float(r.get("mag_err"))
                    if not np.isfinite(mag):
                        continue
                    time_id = f"{label}:{fname}" if label else fname
                    rows.append(
                        dict(
                            time_id=time_id,
                            jd=jd_val,
                            filter=filt_key,
                            star_id=sid,
                            mag_inst=mag,
                            err=err,
                            file=fname,
                            dataset=label,
                            photometry_source=source_info["source"],
                            mag_input_column=source_info["mag_column"],
                            mag_error_input_column=source_info["mag_error_column"],
                        )
                    )

        if not rows:
            raise RuntimeError("No photometry rows found for global ensemble")
        return pd.DataFrame(rows)

    def _fill_date_from_jd(self) -> None:
        if self.raw_df.empty or "JD" not in self.raw_df.columns:
            return
        date_col = self.raw_df.get("date", pd.Series(["unknown"] * len(self.raw_df))).astype(str)
        jd = self.raw_df["JD"].to_numpy(float)
        fill_mask = np.isfinite(jd) & date_col.astype(str).str.strip().str.lower().isin(["", "nan", "none", "unknown"])
        if not np.any(fill_mask):
            self._fill_date_from_filename()
            return
        try:
            times = Time(jd[fill_mask], format="jd").to_datetime()
            date_vals = [t.strftime("%Y-%m-%d") for t in times]
            date_col = date_col.to_numpy(object)
            date_col[fill_mask] = date_vals
            self.raw_df["date"] = date_col
            self._fill_date_from_filename()
        except Exception:
            return

    def _fill_date_from_filename(self) -> None:
        if self.raw_df.empty or "file" not in self.raw_df.columns:
            return
        date_col = self.raw_df.get("date", pd.Series(["unknown"] * len(self.raw_df))).astype(str).to_numpy(object)
        files = self.raw_df["file"].astype(str).tolist()
        for i, fname in enumerate(files):
            if str(date_col[i]).strip().lower() not in ["", "nan", "none", "unknown"]:
                continue
            m = re.search(r"(\d{4}-\d{2}-\d{2})", fname)
            if m:
                date_col[i] = m.group(1)
                continue
            m = re.search(r"(\d{8})", fname)
            if m:
                d = m.group(1)
                date_col[i] = f"{d[:4]}-{d[4:6]}-{d[6:8]}"
        self.raw_df["date"] = date_col

    def _fill_night_id(self) -> None:
        """Fill night_id column from file_manager or disk, then overwrite date with 'Night N' labels."""
        if self.raw_df.empty:
            return

        def _apply_night_labels(nids: list) -> bool:
            if not any(n > 0 for n in nids):
                return False
            date_col = self.raw_df["date"].astype(str).to_numpy(object)
            for i, nid in enumerate(nids):
                if nid > 0:
                    date_col[i] = f"Night {nid}"
            self.raw_df["date"] = date_col
            self.raw_df["night_id"] = nids
            n_labeled = sum(1 for n in nids if n > 0)
            self.log(f"[NIGHT] {n_labeled}/{len(nids)} rows → 'Night N' 레이블 적용")
            return True

        # Case 1: night_id column already present with valid values (from step10 CSV)
        if "night_id" in self.raw_df.columns:
            nids = pd.to_numeric(self.raw_df["night_id"], errors="coerce").fillna(0).astype(int).tolist()
            if _apply_night_labels(nids):
                return
            self.log("[NIGHT] night_id 컬럼 있으나 모두 0 → JSON fallback 시도")

        # Case 2: load night_assignments from file_manager or disk JSON
        if "file" not in self.raw_df.columns:
            self.log("[NIGHT] 'file' 컬럼 없음 → Night 레이블 불가")
            return

        night_id_map: dict[str, int] = {}
        fm = getattr(self, "file_manager", None)
        if fm is not None:
            night_id_map = {Path(k).name: v for k, v in getattr(fm, "night_assignments", {}).items()}
            if night_id_map:
                self.log(f"[NIGHT] file_manager에서 {len(night_id_map)}개 로드")

        if not night_id_map:
            for _label, rdir in self.datasets:
                m = _load_night_assignments_from_disk(Path(rdir))
                if m:
                    night_id_map.update({Path(k).name: v for k, v in m.items()})
                    self.log(f"[NIGHT] disk JSON에서 {len(m)}개 로드: {rdir}")

        if not night_id_map:
            m = _load_night_assignments_from_disk(Path(self.params.P.result_dir))
            if m:
                night_id_map.update({Path(k).name: v for k, v in m.items()})
                self.log(f"[NIGHT] result_dir JSON에서 {len(m)}개 로드")

        if not night_id_map:
            self.log("[NIGHT] night_assignments.json 없음 → 원본 날짜 유지")
            return

        files = self.raw_df["file"].astype(str).tolist()
        nids = [night_id_map.get(Path(fn).name, 0) for fn in files]
        n_matched = sum(1 for n in nids if n > 0)
        self.log(f"[NIGHT] 파일명 매칭: {n_matched}/{len(nids)} (샘플 file: {Path(files[0]).name if files else '-'})")
        _apply_night_labels(nids)

    def _fill_airmass_from_headers(
        self, *, allow_fits_fallback: bool = True
    ) -> None:
        if self.raw_df.empty or "file" not in self.raw_df.columns:
            return
        airmass_series = pd.to_numeric(self.raw_df.get("airmass", np.nan), errors="coerce")
        refill_mask = airmass_series.isna() | ~airmass_series.map(is_reasonable_airmass)
        if int(refill_mask.sum()) == 0:
            return
        dataset_map = {label: Path(path) for label, path in self.datasets} if self.datasets else {}
        filled_from_headers = 0
        filled_from_header_times = 0
        filled_from_fits = 0
        total_missing = int(refill_mask.sum())
        site_lat = _safe_float(getattr(self.params.P, "site_lat_deg", np.nan), np.nan)
        site_lon = _safe_float(getattr(self.params.P, "site_lon_deg", np.nan), np.nan)
        site_alt = _safe_float(getattr(self.params.P, "site_alt_m", 0.0), 0.0)
        site_tz = _safe_float(getattr(self.params.P, "site_tz_offset_hours", 0.0), 0.0)
        target_ra = _safe_float(getattr(self.params.P, "target_ra_deg", np.nan), np.nan)
        target_dec = _safe_float(getattr(self.params.P, "target_dec_deg", np.nan), np.nan)
        formula = getattr(self.params.P, "airmass_formula", None)

        def _filename_map(names: pd.Series, values: pd.Series) -> dict[str, float]:
            mapped: dict[str, float] = {}
            for fname, value in zip(names.astype(str), values):
                if not np.isfinite(value):
                    continue
                mapped[str(fname)] = float(value)
                mapped[Path(str(fname)).name] = float(value)
            return mapped

        def _map_rows(files: pd.Series, mapping: dict[str, float]) -> pd.Series:
            return files.astype(str).map(
                lambda value: mapping.get(value, mapping.get(Path(value).name, np.nan))
            )

        if not dataset_map:
            dataset_map = {"": Path(self.params.P.result_dir)}
            self.raw_df["dataset"] = self.raw_df.get("dataset", "")
        for label, result_dir in dataset_map.items():
            if "dataset" in self.raw_df.columns:
                sel = (self.raw_df["dataset"].astype(str) == str(label)) & refill_mask
            else:
                sel = refill_mask
            if not np.any(sel):
                continue
            headers_path = step1_dir(result_dir) / "headers.csv"
            if not headers_path.exists():
                headers_path = Path(result_dir) / "headers.csv"
            if not headers_path.exists():
                continue
            try:
                hdf = pd.read_csv(headers_path)
            except Exception:
                continue
            if "Filename" not in hdf.columns:
                continue
            col_air = None
            for cand in ("AIRMASS", "airmass"):
                if cand in hdf.columns:
                    col_air = cand
                    break
            if col_air is not None:
                header_airmass = pd.to_numeric(hdf[col_air], errors="coerce")
                valid_header = header_airmass.map(is_reasonable_airmass)
                amap = _filename_map(
                    hdf.loc[valid_header, "Filename"],
                    header_airmass.loc[valid_header],
                )
                vals = _map_rows(self.raw_df.loc[sel, "file"], amap)
                n_fill = int(vals.notna().sum())
                if n_fill > 0:
                    self.raw_df.loc[sel, "airmass"] = vals
                    filled_from_headers += n_fill

            current_airmass = pd.to_numeric(self.raw_df.get("airmass", np.nan), errors="coerce")
            remaining = sel & (
                current_airmass.isna()
                | ~current_airmass.map(is_reasonable_airmass)
            )
            can_compute = all(
                np.isfinite(value)
                for value in (site_lat, site_lon, target_ra, target_dec)
            )
            if not np.any(remaining) or not can_compute:
                continue

            jd_values = pd.Series(np.nan, index=hdf.index, dtype=float)
            if "JD" in hdf.columns:
                jd_values = pd.to_numeric(hdf["JD"], errors="coerce")
            if jd_values.notna().sum() == 0 and "DATE-OBS" in hdf.columns:
                date_mask = hdf["DATE-OBS"].notna()
                try:
                    jd_values.loc[date_mask] = Time(
                        hdf.loc[date_mask, "DATE-OBS"].astype(str).tolist(),
                        scale="utc",
                    ).jd
                except Exception:
                    pass
            computed = compute_airmass_from_jd_array(
                jd_values.to_numpy(dtype=float),
                target_ra,
                target_dec,
                site_lat,
                site_lon,
                site_alt,
                formula=formula,
            )
            valid_computed = np.array(
                [is_reasonable_airmass(value, max_airmass=12.0) for value in computed],
                dtype=bool,
            )
            cmap = _filename_map(
                hdf.loc[valid_computed, "Filename"],
                pd.Series(computed[valid_computed], index=hdf.index[valid_computed]),
            )
            vals = _map_rows(self.raw_df.loc[remaining, "file"], cmap)
            n_fill = int(vals.notna().sum())
            if n_fill > 0:
                self.raw_df.loc[remaining, "airmass"] = vals
                filled_from_header_times += n_fill

        airmass_series = pd.to_numeric(self.raw_df.get("airmass", np.nan), errors="coerce")
        refill_mask = airmass_series.isna() | ~airmass_series.map(is_reasonable_airmass)
        if allow_fits_fallback and np.any(refill_mask):
            if np.isfinite(site_lat) and np.isfinite(site_lon):
                fits_cache: dict[str, float | None] = {}
                for row_idx in self.raw_df.index[refill_mask]:
                    fname = str(self.raw_df.at[row_idx, "file"])
                    dataset_key = (
                        str(self.raw_df.at[row_idx, "dataset"])
                        if "dataset" in self.raw_df.columns
                        else ""
                    )
                    cache_key = f"{dataset_key}:{Path(fname).name}"
                    if cache_key in fits_cache:
                        cached_airmass = fits_cache[cache_key]
                        if cached_airmass is not None:
                            self.raw_df.at[row_idx, "airmass"] = cached_airmass
                            filled_from_fits += 1
                        continue
                    try:
                        src_path = Path(self.params.get_file_path(fname))
                    except Exception:
                        src_path = None
                    if src_path is None or not src_path.exists():
                        fits_cache[cache_key] = None
                        continue
                    try:
                        with fits.open(src_path) as hdul:
                            hdr = hdul[0].header
                        info = compute_airmass_from_header(
                            hdr,
                            site_lat,
                            site_lon,
                            site_alt,
                            site_tz,
                            formula=formula,
                        )
                    except Exception:
                        fits_cache[cache_key] = None
                        continue
                    airmass_val = float(info.get("airmass", np.nan))
                    airmass_source = str(info.get("airmass_source", "") or "")
                    if not is_reasonable_airmass(airmass_val, max_airmass=12.0):
                        fits_cache[cache_key] = None
                        continue
                    if airmass_source == "header_suspicious":
                        fits_cache[cache_key] = None
                        continue
                    fits_cache[cache_key] = airmass_val
                    self.raw_df.at[row_idx, "airmass"] = airmass_val
                    filled_from_fits += 1

        filled = filled_from_headers + filled_from_header_times + filled_from_fits
        if filled > 0:
            msg = f"[LOAD] Filled airmass: {filled}/{total_missing}"
            details = []
            if filled_from_headers > 0:
                details.append(f"headers={filled_from_headers}")
            if filled_from_header_times > 0:
                details.append(f"header_times={filled_from_header_times}")
            if filled_from_fits > 0:
                details.append(f"computed={filled_from_fits}")
            if details:
                msg += f" ({', '.join(details)})"
            self.log(msg)

    def _load_comp_selection(self) -> None:
        def _load_from_step9() -> bool:
            rd = Path(self.params.P.result_dir)
            tid, cids = _load_step9_selection_ids(rd)
            if tid is not None and not self._target_id_text():
                self.target_edit.setText(str(int(tid)))
            if not cids:
                return False
            self.comp_active_ids = list(cids)
            self.comp_candidate_ids = list(cids)
            return True

        if not self.datasets:
            return
        base_dir = step9_lc_dir(self.params.P.result_dir)
        sel_path = base_dir / "comp_selection.json"
        if not sel_path.exists() and self.datasets:
            sel_path = step9_lc_dir(self.datasets[0][1]) / "comp_selection.json"
        if not sel_path.exists():
            if _load_from_step9():
                self._update_comp_label()
            return
        try:
            data = json.loads(sel_path.read_text(encoding="utf-8"))
            self.comp_active_ids = [int(x) for x in data.get("comp_active_ids", []) if str(x).strip()]
            self.comp_candidate_ids = [int(x) for x in data.get("comp_candidate_ids", []) if str(x).strip()]
            if not self.comp_active_ids and self.comp_candidate_ids:
                self.comp_active_ids = list(self.comp_candidate_ids)
            if not self.comp_active_ids and not self.comp_candidate_ids:
                _load_from_step9()
            self._update_comp_label()
        except Exception:
            return

    def _compute_delta_c_map(self, df: pd.DataFrame) -> dict[str, float]:
        if df.empty:
            return {}
        if "color_index" not in df.columns or "color_index_ref" not in df.columns:
            return {}
        delta = pd.to_numeric(df["color_index"], errors="coerce") - pd.to_numeric(
            df["color_index_ref"], errors="coerce"
        )
        if "filter" in df.columns:
            filters = df["filter"].astype(str).map(_normalize_filter_key)
        else:
            filters = pd.Series([""] * len(df))
        out: dict[str, float] = {}
        temp = pd.DataFrame({"filter": filters, "delta_c": delta})
        for fkey, sub in temp.groupby("filter"):
            vals = pd.to_numeric(sub["delta_c"], errors="coerce").to_numpy(float)
            vals = vals[np.isfinite(vals)]
            if vals.size == 0:
                continue
            out[str(fkey)] = float(np.nanmedian(vals))
        if "" not in out and "filter" not in df.columns:
            vals = pd.to_numeric(delta, errors="coerce").to_numpy(float)
            vals = vals[np.isfinite(vals)]
            if vals.size:
                out[""] = float(np.nanmedian(vals))
        return out

    def _load_color_median_table(self, result_dir: Path) -> pd.DataFrame:
        candidates = [
            result_dir / "median_by_ID_filter_wide_cmd.csv",
            result_dir / "median_by_ID_filter_wide.csv",
        ]
        path = next((p for p in candidates if p.exists()), None)
        if path is None:
            return pd.DataFrame()
        try:
            return pd.read_csv(path)
        except Exception:
            return pd.DataFrame()

    def _compute_delta_c_map_from_median(
        self,
        result_dir: Path,
        target_id: int,
        comp_ids: list[int],
        mapping: dict[str, str],
    ) -> dict[str, float]:
        if not mapping:
            return {}
        df = self._load_color_median_table(result_dir)
        if df.empty or "ID" not in df.columns:
            return {}
        ids = pd.to_numeric(df["ID"], errors="coerce").astype("Int64")
        out: dict[str, float] = {}
        for fkey, expr in mapping.items():
            bands = _parse_color_expr(expr)
            if not bands:
                continue
            col_a = None
            col_b = None
            for prefix in ("mag_cal_", "mag_std_", "mag_inst_"):
                cand_a = f"{prefix}{bands[0]}"
                cand_b = f"{prefix}{bands[1]}"
                if cand_a in df.columns and cand_b in df.columns:
                    col_a = cand_a
                    col_b = cand_b
                    break
            if col_a is None or col_b is None:
                continue
            color = pd.to_numeric(df[col_a], errors="coerce") - pd.to_numeric(df[col_b], errors="coerce")
            cmap: dict[int, float] = {}
            for id_val, c in zip(ids.to_numpy(), color.to_numpy(float)):
                if pd.isna(id_val):
                    continue
                cmap[int(id_val)] = float(c) if np.isfinite(c) else np.nan
            target_color = cmap.get(target_id)
            if target_color is None or not np.isfinite(target_color):
                continue
            comp_colors = [cmap.get(int(cid), np.nan) for cid in comp_ids]
            comp_mean = float(np.nanmean(comp_colors)) if comp_colors else np.nan
            if not np.isfinite(comp_mean):
                continue
            out[_normalize_filter_key(fkey)] = float(target_color - comp_mean)
        return out

    def _compute_delta_c_map_from_raw(self, df: pd.DataFrame, mapping: dict[str, str]) -> dict[str, float]:
        if df.empty or not mapping:
            return {}
        if "filter" not in df.columns or "mag" not in df.columns or "comp_avg" not in df.columns:
            return {}
        data = df.copy()
        data["filter_key"] = data["filter"].astype(str).map(_normalize_filter_key)
        data["mag"] = pd.to_numeric(data["mag"], errors="coerce")
        data["comp_avg"] = pd.to_numeric(data["comp_avg"], errors="coerce")
        med_target = data.groupby("filter_key")["mag"].median()
        med_comp = data.groupby("filter_key")["comp_avg"].median()
        out: dict[str, float] = {}
        for fkey, expr in mapping.items():
            bands = _parse_color_expr(expr)
            if not bands:
                continue
            b1, b2 = bands
            if b1 not in med_target.index or b2 not in med_target.index:
                continue
            t1 = _safe_float(med_target.get(b1))
            t2 = _safe_float(med_target.get(b2))
            c1 = _safe_float(med_comp.get(b1))
            c2 = _safe_float(med_comp.get(b2))
            if not all(np.isfinite(v) for v in (t1, t2, c1, c2)):
                continue
            out[_normalize_filter_key(fkey)] = float((t1 - t2) - (c1 - c2))
        return out

    def _refresh_delta_c_map(self) -> None:
        result_dir = self.datasets[0][1] if self.datasets else self.params.P.result_dir
        target_text = self._target_id_text()
        mapping = self.color_map_by_filter
        if mapping:
            delta_map = self._compute_delta_c_map_from_raw(self.raw_df, mapping)
        else:
            delta_map = {}
        if not delta_map:
            if target_text:
                target_id = int(target_text)
                comp_ids = self.comp_active_ids if self.comp_active_ids else self.comp_candidate_ids
                if comp_ids and mapping:
                    delta_map = self._compute_delta_c_map_from_median(
                        Path(result_dir), target_id, comp_ids, mapping
                    )
            if not delta_map:
                delta_map = self._compute_delta_c_map(self.raw_df)
        self.delta_c_map = delta_map

    def _update_color_mode_enabled(self) -> bool:
        has_color = False
        if self.delta_c_map:
            has_color = any(np.isfinite(v) for v in self.delta_c_map.values())
        self.mode_color.setEnabled(has_color)
        if not has_color:
            self.color_status_label.setText("ⓘ Color mode 사용 불가: ΔC 데이터 없음")
        else:
            self.color_status_label.setText("")
        return has_color

    def _delta_c_for_filter(self, fkey: str) -> float:
        if not self.delta_c_map:
            return np.nan
        key = _normalize_filter_key(fkey)
        if key in self.delta_c_map:
            return self.delta_c_map[key]
        if "" in self.delta_c_map:
            return self.delta_c_map[""]
        return np.nan

    def _log_color_index_info(self) -> None:
        mapping = self.color_map_by_filter or getattr(self.params.P, "lightcurve_color_index_by_filter", {}) or {}
        if mapping:
            pairs = ", ".join(f"{k}:{v}" for k, v in mapping.items())
            self.log(f"[COLOR] color_index_by_filter = {pairs}")
        if self.delta_c_map:
            for fkey in sorted(self.delta_c_map):
                self.log(f"[COLOR] ΔC median {fkey or 'all'} = {self.delta_c_map[fkey]:.5f}")

    def _selected_dates(self) -> set:
        """Which nights this fit covers.

        The window reads its date list widget; a batch run is handed the
        set, or None meaning every night present in the data.
        """
        if self.selected_dates is not None:
            return set(self.selected_dates)
        raw = getattr(self, 'raw_df', None)
        if raw is None or getattr(raw, 'empty', True):
            return set()
        if 'date' not in raw.columns:
            return set()
        return {str(v) for v in raw['date'].dropna().unique()}

    def _phase_reference_t0(self, *extra_frames: pd.DataFrame) -> float:
        if self.phase_t0 > 0:
            return float(self.phase_t0)
        frames: list[pd.DataFrame] = []
        for frame in (self.raw_df, self.corrected_df, *extra_frames):
            if isinstance(frame, pd.DataFrame) and not frame.empty:
                frames.append(frame)
        if not frames:
            return np.nan
        finite_chunks: list[np.ndarray] = []
        for frame in frames:
            if "JD" not in frame.columns:
                continue
            jd = pd.to_numeric(frame["JD"], errors="coerce").to_numpy(float)
            finite = jd[np.isfinite(jd)]
            if finite.size:
                finite_chunks.append(finite)
        if not finite_chunks:
            return np.nan
        return float(np.nanmin(np.concatenate(finite_chunks)))

    def _mask_by_ranges(self, df: pd.DataFrame) -> np.ndarray:
        return np.ones(len(df), dtype=bool)

    def fit_and_apply(
        self,
        update_ui: bool = True,
        save_outputs: bool = True,
        selected_dates: set[str] | None = None,
        use_global_k2: bool | None = None,
        target_id_override: int | None = None,
        comp_ids_override: list[int] | None = None,
        sync_controls: bool = True,
    ):
        if update_ui and not self.runtime_mode:
            self._queue_background_fit(
                save_outputs=save_outputs,
                selected_dates=selected_dates,
                use_global_k2=use_global_k2,
                target_id_override=target_id_override,
                comp_ids_override=comp_ids_override,
                sync_controls=sync_controls,
            )
            return
        if sync_controls:
            self._sync_state_from_controls()
        if self.mode == "global":
            self._run_global_ensemble(
                update_ui=update_ui,
                save_outputs=save_outputs,
                target_id_override=target_id_override,
                comp_ids_override=comp_ids_override,
                sync_controls=False,
            )
            return
        if self.mode == "sysrem":
            self._run_sysrem(
                update_ui=update_ui,
                save_outputs=save_outputs,
                target_id_override=target_id_override,
                comp_ids_override=comp_ids_override,
            )
            return
        if self.raw_df.empty:
            if update_ui:
                self._tell_user("info", "Detrend", "Raw 데이터가 없습니다.")
            else:
                raise RuntimeError("Raw 데이터가 없습니다.")
            return

        dates = set(selected_dates or self._selected_dates())
        if not dates:
            if update_ui:
                self._tell_user("info", "Detrend", "날짜를 하나 이상 선택하세요.")
            else:
                raise RuntimeError("날짜를 하나 이상 선택하세요.")
            return

        # Color mode validation
        if self.mode == "color":
            has_delta_c = any(np.isfinite(self._delta_c_for_filter(fkey)) for fkey in self.delta_c_map)
            if not has_delta_c:
                msg = (
                    "Color mode를 사용하려면 ΔC (색지수 차이)가 필요합니다.\n\n"
                    "해결방법:\n"
                    "• Color Index 설정에서 필터별 색지수 조합 선택\n"
                    "• Step 10에서 color_index 컬럼이 있는 데이터 생성\n\n"
                    "현재는 Offset 모드로 진행합니다."
                )
                if update_ui:
                    self._tell_user("warn", "Color Mode", msg)
                    self.color_status_label.setText("⚠ ΔC 없음 - Offset 모드 사용")
                self.mode = "offset"
                if update_ui:
                    self.mode_offset.setChecked(True)
            else:
                if update_ui:
                    self.color_status_label.setText("")
        if use_global_k2 is None:
            use_global_k2 = bool(self._use_global_k2())
        if update_ui:
            self._set_busy_state(True, "Preparing fit...")
        try:
            fit_df = self.raw_df

            # Color mode with global k'' fitting
            global_k2_by_filter: dict[str, tuple[float, float]] = {}

            if self.mode == "color" and use_global_k2:
                if update_ui:
                    self._set_busy_message("Fitting global k''...")
                self.log("[FIT] Global k'' fitting mode enabled")
                all_filters = [""]
                if "filter" in fit_df.columns:
                    all_filters = sorted({_normalize_filter_key(f) for f in fit_df["filter"].astype(str)})

                for fkey in all_filters:
                    sub = fit_df
                    if fkey and "filter" in sub.columns:
                        sub = fit_df[fit_df["filter"].astype(str).map(_normalize_filter_key) == fkey]
                    if sub.empty:
                        continue

                    delta_c_const = self._delta_c_for_filter(fkey)
                    if not np.isfinite(delta_c_const):
                        self.log(f"[FIT] Global k'' for {fkey or 'all'}: ΔC missing, skipped")
                        continue

                    date_mask = sub["date"].astype(str).isin([str(d) for d in dates])
                    sub_all = sub[date_mask]
                    if sub_all.empty:
                        continue

                    mask = self._mask_by_ranges(sub_all)
                    y = sub_all["diff_mag_raw"].to_numpy(float)
                    x_air = sub_all["airmass"].to_numpy(float)
                    x = x_air * float(delta_c_const)
                    err = sub_all.get("diff_err", pd.Series([np.nan] * len(sub_all))).to_numpy(float)
                    w = np.where(np.isfinite(err) & (err > 0), 1.0 / (err * err), 1.0)

                    base_mask = mask & np.isfinite(y) & np.isfinite(x)
                    if np.sum(base_mask) < 5:
                        continue

                    zp, k2, zp_err, k2_err, _cov, _ = self._fit_linear(y, x, w, base_mask, fit_slope=True)
                    global_k2_by_filter[fkey] = (k2, k2_err)
                    self.log(f"[FIT] Global k'' ({fkey or 'all'}): {k2:.5f} ± {k2_err:.5f}")

                    if abs(k2) > 0.15:
                        x_range = np.ptp(x_air[base_mask])
                        self.log(f"[WARNING] |k''| = {abs(k2):.3f} >> 0.05 (비정상)")
                        self.log(f"  → Airmass 범위 ΔX = {x_range:.3f} (좁으면 피팅 불안정)")
                        self.log(f"  → Offset 모드 권장 (k'' 피팅 불가)")

            if update_ui:
                self._set_busy_message("Fitting nightly groups...")
            params_rows = []
            for date_val in sorted(dates):
                sub_date = fit_df[fit_df["date"].astype(str) == str(date_val)]
                if sub_date.empty:
                    continue
                filters = [""]
                if "filter" in sub_date.columns:
                    filters = sorted({_normalize_filter_key(f) for f in sub_date["filter"].astype(str)})

                for fkey in filters:
                    sub = sub_date
                    if fkey:
                        sub = sub_date[sub_date["filter"].astype(str).map(_normalize_filter_key) == fkey]
                    if sub.empty:
                        continue

                    mask = self._mask_by_ranges(sub)
                    y = sub["diff_mag_raw"].to_numpy(float)
                    x_air = sub["airmass"].to_numpy(float)
                    err = sub.get("diff_err", pd.Series([np.nan] * len(sub))).to_numpy(float)
                    w = np.where(np.isfinite(err) & (err > 0), 1.0 / (err * err), 1.0)

                    x = x_air
                    if self.mode == "color":
                        delta_c_const = self._delta_c_for_filter(fkey)
                        if not np.isfinite(delta_c_const):
                            self.log(f"[FIT] {date_val}/{fkey or 'all'}: ΔC missing, skipped")
                            continue
                        x = x_air * float(delta_c_const)

                    base_mask = mask & np.isfinite(y)
                    if self.mode != "offset":
                        base_mask = base_mask & np.isfinite(x)
                    if not np.any(base_mask):
                        n_total = len(sub)
                        n_y = int(np.sum(np.isfinite(y)))
                        n_x = int(np.sum(np.isfinite(x)))
                        n_mask = int(np.sum(mask))
                        self.log(
                            f"[FIT] {date_val}/{fkey or 'all'}: no valid points "
                            f"(N={n_total}, y={n_y}, x={n_x}, mask={n_mask})"
                        )
                        continue

                    cov_zp_slope = 0.0
                    if self.mode == "offset":
                        zp, slope, zp_err, slope_err, _cov, used_mask = self._fit_linear(
                            y, x, w, base_mask, fit_slope=False
                        )
                    elif self.mode == "color" and use_global_k2 and fkey in global_k2_by_filter:
                        k2_global, k2_global_err = global_k2_by_filter[fkey]
                        y_adjusted = y - k2_global * x
                        zp, _, zp_err, _, _cov, used_mask = self._fit_linear(
                            y_adjusted, x, w, base_mask, fit_slope=False
                        )
                        slope = k2_global
                        slope_err = k2_global_err
                        # Sequential fit: ZP and global k'' are from separate regressions,
                        # so the ZP-k'' covariance for this nightly group is zero.
                        cov_zp_slope = 0.0
                    elif self.mode == "color":
                        zp, slope, zp_err, slope_err, cov_zp_slope, used_mask = self._fit_linear(
                            y, x, w, base_mask, fit_slope=True
                        )
                    else:
                        zp, slope, zp_err, slope_err, _cov, used_mask = self._fit_linear(
                            y, x, w, base_mask, fit_slope=False
                        )

                    y_fit = zp + slope * x
                    rms_before = np.nanstd(y[base_mask])
                    rms_after = np.nanstd((y - y_fit)[used_mask]) if np.any(used_mask) else np.nan

                    var_excess = np.nan
                    excess_ratio = np.nan
                    if np.any(used_mask) and np.isfinite(rms_after):
                        err_used = err[used_mask]
                        good_err = np.isfinite(err_used) & (err_used > 0)
                        if np.any(good_err):
                            median_var = float(np.nanmedian(err_used[good_err] ** 2))
                            var_excess = float(max(0.0, rms_after ** 2 - median_var))
                            if median_var > 0:
                                excess_ratio = float(var_excess / median_var)
                            if excess_ratio > 2.0:
                                self.log(
                                    f"[WARN] {date_val}/{fkey or 'all'}: "
                                    f"excess variance {excess_ratio:.1f}x photometric noise"
                                )

                    params_rows.append({
                        "date": date_val,
                        "filter": fkey,
                        "zp_offset": zp,
                        "zp_offset_err": zp_err,
                        "ext_slope": slope,
                        "ext_slope_err": slope_err,
                        "cov_zp_slope": cov_zp_slope,
                        "fit_slope": self.mode == "color",  # True = slope was fitted
                        "n_used": int(np.sum(used_mask)),
                        "rms_before": rms_before,
                        "rms_after": rms_after,
                        "var_excess": var_excess,
                        "excess_ratio": excess_ratio,
                        "global_k2": use_global_k2 and self.mode == "color",
                    })

            if not params_rows:
                self.log("[FIT] No fit groups. Check airmass/ΔC/Date selection.")
                if update_ui:
                    self._tell_user("info", "Detrend", "피팅할 데이터가 없습니다.")
                else:
                    raise RuntimeError("피팅할 데이터가 없습니다.")
                return

            if update_ui:
                self._set_busy_message("Applying corrections...")
            self.params_df = pd.DataFrame(params_rows)
            self.corrected_df = self._apply_params(self.raw_df, self.params_df)

            if update_ui:
                self._set_busy_message("Refreshing plots...")
                self._update_results_table()
                self._update_plots()
                self._update_analysis_panel()
            self.log(f"[FIT] Applied corrections for {len(self.params_df)} groups")
            if update_ui:
                self._log_fit_summary()

            if save_outputs:
                if update_ui:
                    self._set_busy_message("Saving outputs...")
                self._save_comprehensive_results()
        finally:
            if update_ui:
                self._set_busy_state(False)

    def _apply_bjd_to_raw_df(self, target_id: int | None = None) -> None:
        """Convert raw_df["JD"] (plain JD_UTC) to BJD_TDB and store in BJD_TDB column.

        Idempotent: if BJD_TDB is already populated (non-NaN), skips to avoid
        double-applying the light-travel correction.
        JD column is intentionally left as original JD_UTC so downstream
        consumers that expect plain JD are not surprised.
        """
        if self.raw_df.empty or "JD" not in self.raw_df.columns:
            return
        # Idempotent guard: skip if BJD_TDB already filled
        if "BJD_TDB" in self.raw_df.columns and self.raw_df["BJD_TDB"].notna().any():
            return
        site_lat = float(getattr(self.params.P, "site_lat_deg", np.nan))
        site_lon = float(getattr(self.params.P, "site_lon_deg", np.nan))
        site_alt = float(getattr(self.params.P, "site_alt_m", 0.0))
        if not (np.isfinite(site_lat) and np.isfinite(site_lon)):
            return
        # Target RA/Dec: try master_catalog first, then params
        tgt_ra, tgt_dec = np.nan, np.nan
        if target_id is not None and self.datasets:
            result_dir = Path(self.datasets[0][1])
            step9_out = step8_selection_dir(result_dir)
            for path in list(step9_out.glob("master_catalog_*.tsv")) + [step9_out / "master_catalog.tsv"]:
                if not path.exists():
                    continue
                try:
                    df_cat = read_csv_int64_source_id(path, sep="\t")
                    row = df_cat[pd.to_numeric(df_cat.get("ID", pd.Series([])), errors="coerce") == int(target_id)]
                    if not row.empty and "ra_deg" in df_cat.columns:
                        tgt_ra = float(pd.to_numeric(row["ra_deg"].values[0], errors="coerce"))
                        tgt_dec = float(pd.to_numeric(row["dec_deg"].values[0], errors="coerce"))
                        if np.isfinite(tgt_ra) and np.isfinite(tgt_dec):
                            break
                except Exception:
                    continue
        if not np.isfinite(tgt_ra):
            # `P.target.ra_deg` was the old spelling; the parameter model builds
            # no nested `target` object, so this fallback returned None every
            # time and BJD_TDB came back all-NaN on any workspace whose Step 8
            # had not run. Step 9 had the same dead lookup and was fixed on
            # 2026-08-18; this copy survived in the window until the move.
            def _cfg(*names):
                for name in names:
                    try:
                        value = float(getattr(self.params.P, name, np.nan))
                    except (TypeError, ValueError):
                        continue
                    if np.isfinite(value):
                        return value
                return np.nan

            tgt_ra = _cfg("target_ra_deg", "ra_deg")
            tgt_dec = _cfg("target_dec_deg", "dec_deg")
        if not (np.isfinite(tgt_ra) and np.isfinite(tgt_dec)):
            return
        jd_arr = self.raw_df["JD"].to_numpy(float)
        bjd_arr = compute_bjd_tdb_array(jd_arr, tgt_ra, tgt_dec, site_lat, site_lon, site_alt)
        valid = np.isfinite(bjd_arr)
        if valid.any():
            if "BJD_TDB" not in self.raw_df.columns:
                self.raw_df["BJD_TDB"] = np.nan
            self.raw_df.loc[valid, "BJD_TDB"] = bjd_arr[valid]
            # JD column is kept as original JD_UTC — do not overwrite
            delta = np.nanmedian(bjd_arr[valid] - jd_arr[valid]) * 86400
            self.log(f"[BJD] BJD_TDB computed ({valid.sum()} pts, median correction {delta:+.1f}s)")

    def _run_sysrem(
        self,
        update_ui: bool = True,
        save_outputs: bool = True,
        target_id_override: int | None = None,
        comp_ids_override: list[int] | None = None,
    ) -> None:
        """SYSREM systematic-noise removal (Tamuz, Mazeh & Zucker 2005)."""
        from apex.analysis.light_curve.sysrem import sysrem, apply_sysrem

        if update_ui:
            self._set_busy_state(True, "Loading Step 7 data for SYSREM...")
        try:
            df_all = self._load_global_ensemble_df(
                target_id_override=target_id_override,
                comp_ids_override=comp_ids_override,
            )
        except Exception as exc:
            if update_ui:
                self._tell_user("warn", "SYSREM", str(exc))
            else:
                raise
            return
        finally:
            if update_ui:
                self._set_busy_state(False)

        target_id = int(target_id_override or self._target_id_text())
        if comp_ids_override is not None:
            comp_ids = [int(c) for c in comp_ids_override if int(c) != target_id]
        else:
            comp_ids = self.comp_active_ids or self.comp_candidate_ids
            comp_ids = [int(c) for c in comp_ids if int(c) != target_id]

        if not comp_ids:
            if update_ui:
                self._tell_user("warn", "SYSREM", "비교성 ID가 없습니다.")
            return

        n_iter = int(self.sysrem_iter)
        n_apply = min(int(self.sysrem_apply), n_iter)

        filters = sorted(df_all["filter"].astype(str).map(_normalize_filter_key).unique())
        if not filters:
            filters = [""]

        corrected_rows: list[pd.DataFrame] = []
        params_rows: list[dict] = []

        for fkey in filters:
            if update_ui:
                self._set_busy_message(f"SYSREM filter={fkey or 'all'}...")

            sub = df_all
            if fkey:
                sub = df_all[df_all["filter"].astype(str).map(_normalize_filter_key) == fkey]
            if sub.empty:
                continue

            frames = sorted(sub["time_id"].astype(str).unique())
            if len(frames) < 5:
                self.log(f"[SYSREM] {fkey or 'all'}: too few frames ({len(frames)}), skipping")
                continue

            n_f = len(frames)

            # Build comp star matrix via pivot_table — avoids O(comps × rows) iterrows
            comp_data = sub[sub["star_id"].isin(comp_ids)].copy()
            comp_data["time_id"] = comp_data["time_id"].astype(str)
            comp_data["star_id"] = comp_data["star_id"].astype(int)

            mag_wide = (
                comp_data.pivot_table(
                    index="star_id", columns="time_id",
                    values="mag_inst", aggfunc="first"
                )
                .reindex(index=comp_ids, columns=frames)
            )
            mag_mat = mag_wide.to_numpy(float)   # (n_comps, n_frames)

            if "err" in comp_data.columns:
                err_wide = (
                    comp_data.pivot_table(
                        index="star_id", columns="time_id",
                        values="err", aggfunc="first"
                    )
                    .reindex(index=comp_ids, columns=frames)
                )
                err_raw = err_wide.to_numpy(float)
                # zero / negative errors → NaN so sysrem uses fallback weight
                err_mat = np.where(err_raw > 0, err_raw, np.nan)
            else:
                err_mat = np.full_like(mag_mat, np.nan)

            # Run SYSREM on comp stars
            result = sysrem(mag_mat, err_mat, n_iter=n_iter)

            rms_info = " → ".join(
                f"{c.rms_before:.4f}→{c.rms_after:.4f}" for c in result.components
            )
            self.log(f"[SYSREM] {fkey or 'all'}: RMS per iter: {rms_info}")

            for comp in result.components:
                params_rows.append({
                    "filter": fkey,
                    "iteration": comp.iteration,
                    "rms_before": comp.rms_before,
                    "rms_after": comp.rms_after,
                })

            # Apply to target — build target vectors via pivot (same pattern)
            tgt_data = sub[sub["star_id"] == target_id].copy()
            tgt_data["time_id"] = tgt_data["time_id"].astype(str)

            frame_idx = {f: i for i, f in enumerate(frames)}

            tgt_inst = np.full(n_f, np.nan)
            tgt_err  = np.full(n_f, np.nan)
            tgt_jd   = np.full(n_f, np.nan)
            tgt_tid  = np.array([""] * n_f, dtype=object)

            if not tgt_data.empty:
                tgt_pivot_mag = (
                    tgt_data.pivot_table(index="time_id", values="mag_inst", aggfunc="first")
                    .reindex(frames)
                )
                tgt_pivot_err = (
                    tgt_data.pivot_table(index="time_id", values="err", aggfunc="first")
                    .reindex(frames)
                ) if "err" in tgt_data.columns else None
                tgt_pivot_jd = (
                    tgt_data.pivot_table(index="time_id", values="jd", aggfunc="first")
                    .reindex(frames)
                ) if "jd" in tgt_data.columns else None

                tgt_inst = tgt_pivot_mag["mag_inst"].to_numpy(float)
                if tgt_pivot_err is not None:
                    tgt_err = tgt_pivot_err["err"].to_numpy(float)
                if tgt_pivot_jd is not None:
                    tgt_jd = tgt_pivot_jd["jd"].to_numpy(float)
                for i, f in enumerate(frames):
                    tgt_tid[i] = f

            tgt_corr = apply_sysrem(tgt_inst, tgt_err, result, n_components=n_apply)

            # Weighted mean of comp stars per frame (consistent with _build_ensemble_series).
            # err_mat rows where err≤0 are NaN; those positions get weight=0.
            w_mat = np.where(np.isfinite(err_mat) & (err_mat > 0), 1.0 / err_mat**2, 0.0)
            w_sum = w_mat.sum(axis=0)                              # (n_frames,)
            has_weight = w_sum > 0
            comp_mean_inst = np.full(n_f, np.nan)
            if has_weight.any():
                comp_mean_inst[has_weight] = (
                    np.nansum(w_mat[:, has_weight] * mag_mat[:, has_weight], axis=0)
                    / w_sum[has_weight]
                )
            # Fall back to unweighted nanmean for frames with no valid errors
            no_weight = ~has_weight
            if no_weight.any():
                comp_mean_inst[no_weight] = np.nanmean(mag_mat[:, no_weight], axis=0)

            diff_raw  = tgt_inst - comp_mean_inst
            diff_corr = tgt_corr - comp_mean_inst

            valid = np.isfinite(tgt_inst) & np.isfinite(comp_mean_inst)
            if not valid.any():
                self.log(f"[SYSREM] {fkey or 'all'}: no valid target frames")
                continue

            fdf = pd.DataFrame({
                "JD":            tgt_jd[valid],
                "time_id":       tgt_tid[valid],
                "filter":        fkey,
                "ID":            target_id,
                "diff_mag_raw":  diff_raw[valid],
                "diff_mag_corr": diff_corr[valid],
                "diff_err":      tgt_err[valid],
                "diff_err_corr": tgt_err[valid],   # SYSREM doesn't propagate formal errors
            })

            # Merge date/airmass from raw_df if available
            if not self.raw_df.empty and "JD" in self.raw_df.columns:
                merge_cols = [c for c in ("date", "airmass") if c in self.raw_df.columns]
                if merge_cols:
                    extra = self.raw_df[["JD"] + merge_cols].drop_duplicates("JD")
                    fdf = fdf.merge(extra, on="JD", how="left")

            corrected_rows.append(fdf)

        if not corrected_rows:
            if update_ui:
                self._tell_user("warn", "SYSREM", "유효한 결과가 없습니다.")
            return

        corrected = pd.concat(corrected_rows, ignore_index=True)
        corrected["correction_mode"] = "sysrem"
        corrected["correction_formula"] = f"SYSREM n_iter={n_iter} n_apply={n_apply}"

        self.corrected_df = corrected
        self.params_df    = pd.DataFrame(params_rows)

        if update_ui:
            self._update_results_table()
            self._update_plots()
            n_pts = len(corrected)
            self.log(f"[SYSREM] Done: {n_pts} corrected points, {len(result.components)} components")

        if save_outputs:
            self._save_comprehensive_results()

    def _run_global_ensemble(
        self,
        update_ui: bool = True,
        save_outputs: bool = True,
        target_id_override: int | None = None,
        comp_ids_override: list[int] | None = None,
        sync_controls: bool = True,
    ) -> None:
        if sync_controls:
            self._sync_state_from_controls()
        if update_ui:
            self._set_busy_state(True, "Loading Step 7 forced photometry...")
        try:
            try:
                df_global = self._load_global_ensemble_df(
                    target_id_override=target_id_override,
                    comp_ids_override=comp_ids_override,
                )
            except Exception as e:
                if update_ui:
                    self._tell_user("warn", "Global Ensemble", str(e))
                else:
                    raise
                return

            target_id = int(target_id_override) if target_id_override is not None else int(self._target_id_text())
            if comp_ids_override is not None:
                comp_ids = [int(c) for c in comp_ids_override if str(c).strip() and int(c) != target_id]
            else:
                comp_ids = self.comp_active_ids or self.comp_candidate_ids
                comp_ids = [int(c) for c in comp_ids if str(c).strip() and int(c) != target_id]
            if not comp_ids:
                if update_ui:
                    self._tell_user("warn", "Global Ensemble", "비교성 ID가 필요합니다.")
                else:
                    raise RuntimeError("비교성 ID가 필요합니다.")
                return

            self.log(
                "[GLOBAL] min_comps={mc} sigma={sg} iters={it} rms_pct={rp} rms_thr={rt} frame_sigma={fs} gauge={g}".format(
                    mc=self.global_min_comps,
                    sg=self.global_sigma,
                    it=self.global_iters,
                    rp=self.global_rms_pct,
                    rt=self.global_rms_threshold,
                    fs=self.global_frame_sigma,
                    g=self.global_gauge,
                )
            )

            if update_ui:
                self._set_busy_message("Solving global ensemble...")
            try:
                result = solve_global_ensemble(
                    df_global,
                    target_id=target_id,
                    comp_ids=comp_ids,
                    min_comps=self.global_min_comps,
                    sigma=self.global_sigma,
                    n_iter=self.global_iters,
                    gauge=self.global_gauge,
                    per_filter=True,
                    robust=self.global_robust,
                    rms_clip_pct=self.global_rms_pct,
                    rms_clip_threshold=self.global_rms_threshold if self.global_rms_threshold > 0 else None,
                    frame_sigma=self.global_frame_sigma,
                    interp_missing=self.global_interp_missing,
                    normalize_target=self.global_normalize,
                    rescale_errors=self.global_rescale_errors,
                    log=self.log,
                )
            except Exception as e:
                if update_ui:
                    self._tell_user("warn", "Global Ensemble", f"Fit failed: {e}")
                else:
                    raise
                return

            zp_df = result.get("zp_df", pd.DataFrame()).copy()
            mean_df = result.get("mean_df", pd.DataFrame()).copy()
            lc_df = result.get("lc_df", pd.DataFrame()).copy()
            diagnostics = result.get("diagnostics", {}) or {}
            self.log(
                "[GLOBAL] Solver output: input_rows={rows} frames={frames} zp_rows={zp} mean_rows={mean} lc_rows={lc}".format(
                    rows=len(df_global),
                    frames=df_global["time_id"].astype(str).nunique() if "time_id" in df_global.columns else 0,
                    zp=len(zp_df),
                    mean=len(mean_df),
                    lc=len(lc_df),
                )
            )
            for fkey, diag in diagnostics.get("filters", {}).items():
                network = diag.get("network", {})
                self.log(
                    f"[GLOBAL] filter={fkey}: connected={network.get('connected', False)}, "
                    f"components={network.get('n_components', 0)}, "
                    f"active_comps={diag.get('n_comp_final', 0)}, "
                    f"removed_comps={len(diag.get('removed_comps', []))}, "
                    f"imputed_errors={diag.get('n_errors_imputed', 0)}"
                )
            if zp_df.empty:
                msg = (
                    "Global ensemble produced no zeropoint rows.\n"
                    "Check Step 7 forced photometry / comparison-star availability / filter mapping."
                )
                if update_ui:
                    self._tell_user("warn", "Global Ensemble", msg)
                else:
                    raise RuntimeError(msg)
                return

            self.global_input_df = df_global
            self.params_df = zp_df
            self.global_mean_df = mean_df
            self.corrected_df = lc_df
            self.raw_df = self.corrected_df.copy()
            self.global_diagnostics = diagnostics

            if "JD" not in self.corrected_df.columns and "jd" in self.corrected_df.columns:
                self.corrected_df["JD"] = self.corrected_df["jd"]
            if "JD" not in self.raw_df.columns and "jd" in self.raw_df.columns:
                self.raw_df["JD"] = self.raw_df["jd"]
            if "JD" not in self.raw_df.columns:
                self.raw_df["JD"] = np.nan
            self._apply_bjd_to_raw_df(target_id)
            if "date" not in self.raw_df.columns:
                self.raw_df["date"] = "unknown"
            self._fill_date_from_jd()
            self._fill_night_id()
            self.corrected_df = self.raw_df.copy()

            if update_ui:
                self._set_busy_message("Refreshing plots...")
                self._populate_date_list()
                self._refresh_filter_combo(self.raw_df.get("filter", pd.Series([], dtype=str)).astype(str).tolist())
                self._update_results_table()
                self._update_plots()
                self._update_analysis_panel()
            self.log("[GLOBAL] Fit complete")

            if save_outputs:
                if update_ui:
                    self._set_busy_message("Saving outputs...")
                self._save_comprehensive_results()
        finally:
            if update_ui:
                self._set_busy_state(False)

    def _log_fit_summary(self):
        """Log fit summary with astronomical interpretation."""
        if self.params_df.empty:
            return

        self.log("\n" + "=" * 50)
        self.log("[SUMMARY] 피팅 결과 분석")
        self.log("=" * 50)

        # Global RMS improvement (all data combined)
        if not self.corrected_df.empty:
            raw_all = pd.to_numeric(self.corrected_df["diff_mag_raw"], errors="coerce").to_numpy(float)
            corr_all = pd.to_numeric(self.corrected_df["diff_mag_corr"], errors="coerce").to_numpy(float)
            raw_all = raw_all[np.isfinite(raw_all)]
            corr_all = corr_all[np.isfinite(corr_all)]

            if raw_all.size > 0 and corr_all.size > 0:
                global_rms_before = np.std(raw_all)
                global_rms_after = np.std(corr_all)
                global_improve = (1 - global_rms_after / global_rms_before) * 100 if global_rms_before > 0 else 0
                self.log(f"  전체 RMS: {global_rms_before:.4f} → {global_rms_after:.4f} mag ({global_improve:+.1f}%)")

                if self.mode == "offset" and abs(global_improve) < 1:
                    self.log("  → Offset 모드: 밤별 RMS는 동일, 밤간 정렬로 개선")

        # Per-night average RMS (for color mode comparison)
        rms_before = self.params_df["rms_before"].mean()
        rms_after = self.params_df["rms_after"].mean()
        if np.isfinite(rms_before) and np.isfinite(rms_after) and rms_before > 0:
            improvement = (1 - rms_after / rms_before) * 100
            self.log(f"  밤별 평균 RMS: {rms_before:.4f} → {rms_after:.4f} mag ({improvement:+.1f}%)")

            if self.mode == "color" and improvement < -5:
                self.log("  → [경고] 밤별 RMS 증가: 과적합 가능성 (데이터 부족)")

        # Mode-specific analysis
        if self.mode == "color":
            slopes = self.params_df["ext_slope"].to_numpy(float)
            slopes = slopes[np.isfinite(slopes)]
            if slopes.size > 0:
                k2_mean = np.mean(slopes)
                k2_std = np.std(slopes) if slopes.size > 1 else 0.0
                self.log(f"  k'' (2차 소광계수): {k2_mean:.4f} ± {k2_std:.4f} mag/mag")

                # Typical k'' values for reference (should be 0.02-0.05)
                if abs(k2_mean) < 0.01:
                    self.log("  → k'' ≈ 0: 색항 효과 미미 (Offset 모드로 충분)")
                elif abs(k2_mean) <= 0.05:
                    self.log("  → k'' 정상 범위 (0.02-0.05 typical for broad-band)")
                elif abs(k2_mean) <= 0.15:
                    self.log("  → k'' 다소 큼: 데이터 품질 확인 권장")
                else:
                    self.log(f"  → [경고] |k''| = {abs(k2_mean):.2f} >> 0.05 (비정상)")
                    self.log("    원인: Airmass 범위(ΔX)가 좁아 k'' 피팅 불안정")
                    self.log("    해결: Offset 모드 사용 권장")

        # ZP scatter analysis
        zps = self.params_df["zp_offset"].to_numpy(float)
        zps = zps[np.isfinite(zps)]
        if zps.size > 1:
            zp_scatter = np.std(zps)
            self.log(f"  밤별 ZP 산포: {zp_scatter:.4f} mag")
            if zp_scatter > 0.3:
                self.log("  → [경고] ZP 산포가 매우 큼")
                if self.mode == "color":
                    self.log("    원인: k'' 과적합으로 ZP₀가 보상하는 중일 수 있음")
                    self.log("    해결: Offset 모드로 다시 시도 권장")
            elif zp_scatter > 0.1:
                self.log("  → 밤별 조건 변화가 큼 (정상 범위)")

        # Night-by-night alignment verification
        self._log_alignment_verification()

        self.log("=" * 50 + "\n")

    def _log_alignment_verification(self):
        """Log per-night raw mean, ZP, and corrected mean to verify alignment."""
        if self.corrected_df.empty or self.params_df.empty:
            return

        self.log("\n  [밤간 정렬 검증]")
        self.log("  " + "-" * 46)
        self.log(f"  {'Date':<12} {'Raw평균':>10} {'ZP₀':>10} {'Corr평균':>10} {'검증':>6}")
        self.log("  " + "-" * 46)

        df = self.corrected_df
        dates = sorted(df["date"].astype(str).unique())

        alignment_ok = True
        corr_means = []

        for date_val in dates:
            date_mask = df["date"].astype(str) == date_val
            sub = df[date_mask]

            raw_vals = pd.to_numeric(sub["diff_mag_raw"], errors="coerce").to_numpy(float)
            corr_vals = pd.to_numeric(sub["diff_mag_corr"], errors="coerce").to_numpy(float)

            raw_mean = np.nanmean(raw_vals) if np.any(np.isfinite(raw_vals)) else np.nan
            corr_mean = np.nanmean(corr_vals) if np.any(np.isfinite(corr_vals)) else np.nan

            # Get ZP for this date (may have multiple filters, take mean)
            zp_rows = self.params_df[self.params_df["date"].astype(str) == date_val]
            zp_mean = zp_rows["zp_offset"].mean() if not zp_rows.empty else np.nan

            # Verify: raw_mean ≈ ZP (for offset mode, ZP = weighted mean of raw)
            if np.isfinite(raw_mean) and np.isfinite(zp_mean):
                diff = abs(raw_mean - zp_mean)
                check = "✓" if diff < 0.05 else "~"
            else:
                check = "-"

            if np.isfinite(corr_mean):
                corr_means.append(corr_mean)

            self.log(f"  {date_val:<12} {raw_mean:>10.4f} {zp_mean:>10.4f} {corr_mean:>10.4f} {check:>6}")

        self.log("  " + "-" * 46)

        # Check if corrected means are aligned
        if len(corr_means) > 1:
            corr_scatter = np.std(corr_means)
            self.log(f"  보정 후 밤간 평균 산포: {corr_scatter:.4f} mag")
            if corr_scatter < 0.02:
                self.log("  → ✓ 밤간 정렬 양호")
            elif corr_scatter < 0.05:
                self.log("  → ~ 밤간 정렬 적절")
            else:
                self.log("  → ⚠ 밤간 정렬 불완전 (필터별 차이 또는 피팅 문제)")
                alignment_ok = False

    def _fit_linear(self, y, x, w, base_mask, fit_slope: bool = True):
        """Weighted linear fit  y = zp + slope * x  with iterative sigma clipping.

        Returns
        -------
        zp, slope, zp_err, slope_err, cov_zp_slope, used_mask
          cov_zp_slope : off-diagonal covariance Cov(zp, k'') from the WLS
            solution — needed for correct error propagation in color mode.
            NaN when unavailable (offset mode, or n ≤ 2).
        """
        mask = base_mask.copy()
        zp = 0.0
        slope = 0.0
        zp_err = np.nan
        slope_err = np.nan
        cov_zp_slope = 0.0   # offset mode: ZP and slope are independent by design

        for _ in range(self.clip_iters if self.sigma_clip else 1):
            if np.sum(mask) < 2:
                break
            yy = y[mask]
            xx = x[mask]
            ww = w[mask]

            if not fit_slope:
                zp = float(np.average(yy, weights=ww))
                slope = 0.0
                if np.sum(ww) > 0:
                    zp_err = float(np.sqrt(1.0 / np.sum(ww)))
                slope_err = 0.0
                cov_zp_slope = 0.0
            else:
                A = np.vstack([np.ones_like(xx), xx]).T
                Aw = A * np.sqrt(ww[:, None])
                yw = yy * np.sqrt(ww)
                try:
                    coeff, residuals, rank, s = np.linalg.lstsq(Aw, yw, rcond=1e-10)
                    zp = float(coeff[0])
                    slope = float(coeff[1])

                    n = len(yy)
                    if n > 2:
                        resid_fit = yy - (zp + slope * xx)
                        mse = np.sum(ww * resid_fit**2) / (n - 2)
                        try:
                            # (A^T W A) via broadcasting — avoids O(N²) np.diag(ww)
                            AtWA = (A * ww[:, None]).T @ A
                            cov = mse * np.linalg.inv(AtWA)
                            zp_err = float(np.sqrt(cov[0, 0]))
                            slope_err = float(np.sqrt(cov[1, 1]))
                            cov_zp_slope = float(cov[0, 1])
                        except Exception:
                            zp_err = np.nan
                            slope_err = np.nan
                            cov_zp_slope = np.nan
                except Exception:
                    zp, slope = 0.0, 0.0
                    zp_err, slope_err = np.nan, np.nan
                    cov_zp_slope = np.nan

            resid = y - (zp + slope * x)
            sigma = np.nanstd(resid[mask])
            if not self.sigma_clip or not np.isfinite(sigma) or sigma == 0:
                break
            mask = mask & (np.abs(resid) <= self.clip_sigma * sigma)

        return zp, slope, zp_err, slope_err, cov_zp_slope, mask

    def _apply_params(self, df: pd.DataFrame, params_df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        out["diff_mag_corr"] = out["diff_mag_raw"]

        if "diff_err" not in out.columns:
            out["diff_err"] = np.nan
        out["diff_err_corr"] = out["diff_err"].copy()

        for _, row in params_df.iterrows():
            date_val = str(row["date"])
            zp = _safe_float(row.get("zp_offset"))
            zp_err = _safe_float(row.get("zp_offset_err", 0.0))
            slope = _safe_float(row.get("ext_slope"))
            slope_err = _safe_float(row.get("ext_slope_err", 0.0))

            idx = out["date"].astype(str) == date_val
            fval = _normalize_filter_key(row.get("filter", ""))
            if fval and "filter" in out.columns:
                idx &= out["filter"].astype(str).map(_normalize_filter_key) == fval

            x = out.loc[idx, "airmass"].to_numpy(float)
            if self.mode == "color":
                delta_c_const = self._delta_c_for_filter(fval)
                if not np.isfinite(delta_c_const):
                    continue
                x = x * float(delta_c_const)

            # fit_slope flag (stored in params_df) is the authoritative indicator of
            # whether a slope was actually fitted.  Avoids float equality pitfalls
            # where a legitimate near-zero fitted slope would be treated as "no slope".
            fit_slope_flag = bool(row.get("fit_slope", False))
            slope_term = slope * x if fit_slope_flag else np.zeros_like(x)
            out.loc[idx, "diff_mag_corr"] = out.loc[idx, "diff_mag_raw"] - zp - slope_term

            raw_err = out.loc[idx, "diff_err"].to_numpy(float)
            raw_var = np.where(np.isfinite(raw_err), raw_err**2, 0.0)
            zp_var    = zp_err**2    if np.isfinite(zp_err)    else 0.0
            slope_var = slope_err**2 if (fit_slope_flag and np.isfinite(slope_err)) else 0.0
            # Load ZP-k'' covariance (stored by color-mode fits; zero for offset mode)
            cov_zp_slope = _safe_float(row.get("cov_zp_slope", 0.0))
            if not np.isfinite(cov_zp_slope):
                cov_zp_slope = 0.0
            x_eff = x if fit_slope_flag else np.zeros_like(x)  # zero axis when slope not fitted
            # Full variance: σ²_raw + σ²_ZP + x²σ²_k'' + 2x·Cov(ZP,k'')
            corr_var = (
                raw_var
                + zp_var
                + x_eff**2 * slope_var
                + 2.0 * x_eff * cov_zp_slope
            )
            # Clamp negative (can occur if cov term dominates near zero airmass range)
            corr_var = np.where(corr_var >= 0, corr_var, raw_var)
            out.loc[idx, "diff_err_corr"] = np.sqrt(corr_var)

        return out

    @staticmethod
    def _stage_step10_path(path: Path) -> Path:
        return path.with_name(f"{path.stem}.__tmp__{path.suffix}")

    def _cleanup_staged_step10_outputs(self, staged_outputs: list[tuple[Path, Path]]) -> None:
        for staged_path, _ in staged_outputs:
            try:
                staged_path.unlink()
            except FileNotFoundError:
                pass
            except Exception:
                pass

    def _promote_staged_step10_outputs(self, target_id: int, staged_outputs: list[tuple[Path, Path]]) -> None:
        if not staged_outputs:
            return
        keep_names = {staged_path.name for staged_path, _ in staged_outputs}
        self._archive_step10_outputs(target_id, keep_names=keep_names)
        try:
            for staged_path, final_path in staged_outputs:
                staged_path.replace(final_path)
        except Exception:
            self._cleanup_staged_step10_outputs(staged_outputs)
            raise

    def _save_comprehensive_results(self) -> None:
        """Save comprehensive result files including formula, corrections, residuals, and summary."""
        if self.corrected_df.empty:
            return
        if self.mode == "global":
            self._save_global_results()
            return
        target_text = self._target_id_text()
        if not target_text:
            return

        target_id = int(target_text)
        out_dir = step10_detrend_dir(self.params.P.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        mode_tag = self.mode if self.mode in ("offset", "color", "sysrem") else "mode"
        if self.mode == "offset":
            formula = "Δm_corr = Δm_raw - ZP₀"
        elif self.mode == "sysrem":
            formula = str(
                self.corrected_df.get(
                    "correction_formula",
                    pd.Series(["SYSREM systematic-component correction"]),
                ).dropna().iloc[0]
                if "correction_formula" in self.corrected_df.columns
                and not self.corrected_df["correction_formula"].dropna().empty
                else "SYSREM systematic-component correction"
            )
        else:
            formula = "Δm_corr = Δm_raw - ZP₀ - k''·ΔC·X"
        lc_path = step10_current_lc_path(self.params.P.result_dir, target_id)
        params_path = step10_current_params_path(self.params.P.result_dir, target_id)
        summary_path = step10_current_summary_path(self.params.P.result_dir, target_id)
        plot_path = step10_current_plot_path(self.params.P.result_dir, target_id)
        meta_path = step10_current_meta_path(self.params.P.result_dir, target_id)
        mode_lc_path = out_dir / f"lightcurve_ID{target_id}_{mode_tag}.csv"
        mode_params_path = out_dir / f"fit_params_ID{target_id}_{mode_tag}.csv"
        mode_summary_path = out_dir / f"summary_ID{target_id}_{mode_tag}.txt"
        mode_plot_path = out_dir / f"plot_ID{target_id}_{mode_tag}.png"
        staged_outputs: list[tuple[Path, Path]] = []
        try:
            # ===== 1. Main corrected light curve with residuals =====
            out_df = self.corrected_df.copy()

            # Add residuals (diff_mag_corr - median)
            corr_vals = pd.to_numeric(out_df["diff_mag_corr"], errors="coerce")
            median_corr = np.nanmedian(corr_vals)
            out_df["residual"] = corr_vals - median_corr

            # Add phase if period is set
            if self.phase_period > 0 and "JD" in out_df.columns:
                jd = pd.to_numeric(out_df["JD"], errors="coerce").to_numpy(float)
                t0 = self._phase_reference_t0(out_df)
                if np.isfinite(t0):
                    out_df["phase"] = ((jd - t0) / self.phase_period) % 1.0

            # Add fit value (what was subtracted)
            out_df["fit_value"] = out_df["diff_mag_raw"] - out_df["diff_mag_corr"]

            # Reorder columns for clarity
            priority_cols = [
                "file", "JD", "date", "filter", "airmass",
                "diff_mag_raw", "diff_mag_corr", "fit_value", "residual",
                "diff_err", "diff_err_corr"
            ]
            if "phase" in out_df.columns:
                priority_cols.insert(3, "phase")
            other_cols = [c for c in out_df.columns if c not in priority_cols]
            out_df = out_df[[c for c in priority_cols if c in out_df.columns] + other_cols]
            out_df = annotate_step10_output(out_df, mode_tag, formula)

            staged_lc_path = self._stage_step10_path(lc_path)
            out_df.to_csv(staged_lc_path, index=False)
            staged_outputs.append((staged_lc_path, lc_path))
            staged_mode_lc_path = self._stage_step10_path(mode_lc_path)
            out_df.to_csv(staged_mode_lc_path, index=False)
            staged_outputs.append((staged_mode_lc_path, mode_lc_path))

            # ===== 2. Fit parameters by date/filter =====
            params_saved = None
            if not self.params_df.empty:
                params_out = self.params_df.copy()
                params_out["correction_mode"] = mode_tag
                params_out["formula"] = formula
                staged_params_path = self._stage_step10_path(params_path)
                params_out.to_csv(staged_params_path, index=False)
                staged_outputs.append((staged_params_path, params_path))
                params_saved = params_path
                staged_mode_params_path = self._stage_step10_path(mode_params_path)
                params_out.to_csv(staged_mode_params_path, index=False)
                staged_outputs.append((staged_mode_params_path, mode_params_path))

            # ===== 3. Summary report (text file) =====
            summary_text = build_detrend_summary_report_text(
                mode=self.mode,
                use_global_k2=self._use_global_k2(),
                delta_c_map=self.delta_c_map,
                params_df=self.params_df,
                df=out_df,
                sigma_clip=self.sigma_clip,
                clip_sigma=self.clip_sigma,
                clip_iters=self.clip_iters,
                phase_period=self.phase_period,
                phase_t0=self.phase_t0,
            )
            staged_summary_path = self._stage_step10_path(summary_path)
            staged_summary_path.write_text(summary_text, encoding="utf-8")
            staged_outputs.append((staged_summary_path, summary_path))
            staged_mode_summary_path = self._stage_step10_path(mode_summary_path)
            staged_mode_summary_path.write_text(summary_text, encoding="utf-8")
            staged_outputs.append((staged_mode_summary_path, mode_summary_path))

            # ===== 4. Save plot as PNG =====
            # Every `_update_plots()` call sits behind `if update_ui:`, so a
            # batch run reaches this with a figure nothing has drawn on — and
            # `savefig` writes a blank page without complaining. The window
            # overrides this to nothing: its canvas already holds the view the
            # user is looking at, and redrawing here would change what it saves.
            self._ensure_plot_drawn()
            staged_plot_path = self._stage_step10_path(plot_path)
            self._plot_figure().savefig(staged_plot_path, dpi=150, bbox_inches="tight")
            staged_outputs.append((staged_plot_path, plot_path))
            staged_mode_plot_path = self._stage_step10_path(mode_plot_path)
            self._plot_figure().savefig(staged_mode_plot_path, dpi=150, bbox_inches="tight")
            staged_outputs.append((staged_mode_plot_path, mode_plot_path))
            staged_meta_path = self._stage_step10_path(meta_path)
            write_step10_current_meta(
                self.params.P.result_dir,
                target_id,
                mode_tag,
                formula,
                lc_path,
                params_path=params_saved,
                summary_path=summary_path,
                plot_path=plot_path,
                meta_path=staged_meta_path,
            )
            staged_outputs.append((staged_meta_path, meta_path))
            self._promote_staged_step10_outputs(target_id, staged_outputs)

            self.log(f"[SAVE] Current light curve: {lc_path.name}")
            self.log(f"[SAVE] Mode light curve: {mode_lc_path.name}")
            if params_saved is not None:
                self.log(f"[SAVE] Current params: {params_path.name}")
                self.log(f"[SAVE] Mode params: {mode_params_path.name}")
            self.log(f"[SAVE] Current summary: {summary_path.name}")
            self.log(f"[SAVE] Mode summary: {mode_summary_path.name}")
            self.log(f"[SAVE] Current plot: {plot_path.name}")
            self.log(f"[SAVE] Mode plot: {mode_plot_path.name}")
            self.log(f"[SAVE] Current meta: {meta_path.name}")

            self.log(f"[SAVE] 모든 결과 저장 완료: {out_dir}")

        except Exception as e:
            self._cleanup_staged_step10_outputs(staged_outputs)
            self.log(f"[SAVE] Failed: {e}")
            self._tell_user("warn", "Save Error", f"저장 실패: {e}")

    def _archive_step10_outputs(self, target_id: int, keep_names: set[str] | None = None) -> None:
        keep_names = keep_names or set()
        out_dir = step10_detrend_dir(self.params.P.result_dir)
        if not out_dir.exists():
            return

        patterns = [
            f"lightcurve_ID{target_id}_*.csv",
            f"fit_params_ID{target_id}_*.csv",
            f"params_ID{target_id}_*.csv",
            f"summary_ID{target_id}_*.txt",
            f"plot_ID{target_id}_*.png",
            f"result_ID{target_id}_*.json",
            f"global_zp_ID{target_id}_*.csv",
            f"global_mean_ID{target_id}_*.csv",
            f"global_diagnostics_ID{target_id}_*.json",
            f"global_zp_ID{target_id}.csv",
            f"global_mean_ID{target_id}.csv",
            f"global_diagnostics_ID{target_id}.json",
        ]
        existing: list[Path] = []
        seen: set[Path] = set()
        for pattern in patterns:
            for path in out_dir.glob(pattern):
                if path.name in keep_names or path.parent != out_dir or not path.is_file():
                    continue
                if path not in seen:
                    seen.add(path)
                    existing.append(path)

        if not existing:
            return

        stamp = datetime.utcnow().strftime("%Y%m%d_%H%M%SZ")
        archive_dir = step10_history_dir(self.params.P.result_dir) / f"ID{target_id}" / stamp
        archive_dir.mkdir(parents=True, exist_ok=True)
        for path in existing:
            try:
                path.rename(archive_dir / path.name)
            except Exception as e:
                self.log(f"[SAVE] archive skipped for {path.name}: {e}")

    def _save_global_results(self) -> None:
        if self.corrected_df.empty:
            return
        target_text = self._target_id_text()
        if not target_text:
            return
        target_id = int(target_text)
        out_dir = step10_detrend_dir(self.params.P.result_dir)
        out_dir.mkdir(parents=True, exist_ok=True)
        formula = "Δm_corr = (mag_inst_target - Z_t) - <M_comp>"
        lc_path = step10_current_lc_path(self.params.P.result_dir, target_id)
        params_path = step10_current_params_path(self.params.P.result_dir, target_id)
        summary_path = step10_current_summary_path(self.params.P.result_dir, target_id)
        plot_path = step10_current_plot_path(self.params.P.result_dir, target_id)
        meta_path = step10_current_meta_path(self.params.P.result_dir, target_id)
        zp_path = step10_current_global_zp_path(self.params.P.result_dir, target_id)
        mean_path = step10_current_global_mean_path(self.params.P.result_dir, target_id)
        diag_path = step10_current_global_diag_path(self.params.P.result_dir, target_id)
        mode_lc_path = out_dir / f"lightcurve_ID{target_id}_global.csv"
        mode_params_path = out_dir / f"fit_params_ID{target_id}_global.csv"
        mode_summary_path = out_dir / f"summary_ID{target_id}_global.txt"
        mode_plot_path = out_dir / f"plot_ID{target_id}_global.png"
        staged_outputs: list[tuple[Path, Path]] = []
        try:
            # Rebuild the figure from the current in-memory state so the saved PNG
            # cannot lag behind the freshly solved global Z_t table.
            self._update_results_table()
            self._update_plots()
            self._plot_redraw()

            out_df = annotate_step10_output(self.corrected_df, "global", formula)
            staged_lc_path = self._stage_step10_path(lc_path)
            out_df.to_csv(staged_lc_path, index=False)
            staged_outputs.append((staged_lc_path, lc_path))
            staged_mode_lc_path = self._stage_step10_path(mode_lc_path)
            out_df.to_csv(staged_mode_lc_path, index=False)
            staged_outputs.append((staged_mode_lc_path, mode_lc_path))

            params_saved = None
            extra_files: list[Path] = []
            if not self.params_df.empty:
                params_out = self.params_df.copy()
                params_out["correction_mode"] = "global"
                params_out["formula"] = formula
                staged_params_path = self._stage_step10_path(params_path)
                params_out.to_csv(staged_params_path, index=False)
                staged_outputs.append((staged_params_path, params_path))
                params_saved = params_path
                staged_mode_params_path = self._stage_step10_path(mode_params_path)
                params_out.to_csv(staged_mode_params_path, index=False)
                staged_outputs.append((staged_mode_params_path, mode_params_path))

                staged_zp_path = self._stage_step10_path(zp_path)
                self.params_df.to_csv(staged_zp_path, index=False)
                staged_outputs.append((staged_zp_path, zp_path))
                extra_files.append(zp_path)

            if not self.global_mean_df.empty:
                staged_mean_path = self._stage_step10_path(mean_path)
                self.global_mean_df.to_csv(staged_mean_path, index=False)
                staged_outputs.append((staged_mean_path, mean_path))
                extra_files.append(mean_path)

            if self.global_diagnostics:
                staged_diag_path = self._stage_step10_path(diag_path)
                with open(staged_diag_path, "w", encoding="utf-8") as f:
                    json.dump(self.global_diagnostics, f, indent=2)
                staged_outputs.append((staged_diag_path, diag_path))
                extra_files.append(diag_path)

            summary_text = build_global_summary_text(
                target_id=target_id,
                params_df=self.params_df,
                corrected_df=self.corrected_df,
                formula=formula,
            )
            staged_summary_path = self._stage_step10_path(summary_path)
            staged_summary_path.write_text(summary_text, encoding="utf-8")
            staged_outputs.append((staged_summary_path, summary_path))
            staged_mode_summary_path = self._stage_step10_path(mode_summary_path)
            staged_mode_summary_path.write_text(summary_text, encoding="utf-8")
            staged_outputs.append((staged_mode_summary_path, mode_summary_path))
            # Every `_update_plots()` call sits behind `if update_ui:`, so a
            # batch run reaches this with a figure nothing has drawn on — and
            # `savefig` writes a blank page without complaining. The window
            # overrides this to nothing: its canvas already holds the view the
            # user is looking at, and redrawing here would change what it saves.
            self._ensure_plot_drawn()
            staged_plot_path = self._stage_step10_path(plot_path)
            self._plot_figure().savefig(staged_plot_path, dpi=150, bbox_inches="tight")
            staged_outputs.append((staged_plot_path, plot_path))
            staged_mode_plot_path = self._stage_step10_path(mode_plot_path)
            self._plot_figure().savefig(staged_mode_plot_path, dpi=150, bbox_inches="tight")
            staged_outputs.append((staged_mode_plot_path, mode_plot_path))
            staged_meta_path = self._stage_step10_path(meta_path)
            write_step10_current_meta(
                self.params.P.result_dir,
                target_id,
                "global",
                formula,
                lc_path,
                params_path=params_saved,
                summary_path=summary_path,
                plot_path=plot_path,
                extra_files=extra_files,
                meta_path=staged_meta_path,
            )
            staged_outputs.append((staged_meta_path, meta_path))
            self._promote_staged_step10_outputs(target_id, staged_outputs)

            self.log(f"[SAVE] Current light curve: {lc_path.name}")
            self.log(f"[SAVE] Mode light curve: {mode_lc_path.name}")
            if params_saved is not None:
                self.log(f"[SAVE] Current params: {params_path.name}")
                self.log(f"[SAVE] Mode params: {mode_params_path.name}")
                self.log(f"[SAVE] Global ZP: {zp_path.name}")
            if not self.global_mean_df.empty:
                self.log(f"[SAVE] Global mean: {mean_path.name}")
            if self.global_diagnostics:
                self.log(f"[SAVE] Diagnostics: {diag_path.name}")
            self.log(f"[SAVE] Current summary: {summary_path.name}")
            self.log(f"[SAVE] Mode summary: {mode_summary_path.name}")
            self.log(f"[SAVE] Current plot: {plot_path.name}")
            self.log(f"[SAVE] Mode plot: {mode_plot_path.name}")
            self.log(f"[SAVE] Current meta: {meta_path.name}")
        except Exception as e:
            self._cleanup_staged_step10_outputs(staged_outputs)
            self.log(f"[SAVE] Global failed: {e}")
            self._tell_user("warn", "Save Error", f"저장 실패: {e}")

    def _global_zp_is_loaded(self) -> bool:
        if self.params_df.empty:
            return False
        cols = {str(c) for c in self.params_df.columns}
        return "time_id" in cols and "Z" in cols

    def _get_global_zp_df(self, silent: bool = True) -> pd.DataFrame:
        """Return a plottable global-ZP table from memory or saved current outputs."""
        if self._global_zp_is_loaded():
            return self.params_df.copy()

        if self.mode == "global":
            self._restore_saved_global_payload(silent=silent)
            if self._global_zp_is_loaded():
                return self.params_df.copy()

        if self.params_df.empty:
            return pd.DataFrame()

        cols = {str(c) for c in self.params_df.columns}
        if not {"time_id", "Z"}.issubset(cols):
            return pd.DataFrame()
        return self.params_df.copy()

    def _restore_saved_global_payload(self, silent: bool = True) -> bool:
        if self.mode != "global":
            return False
        if self._global_zp_is_loaded():
            return True

        target_text = self.target_edit.text().strip()
        if not target_text:
            return False
        try:
            target_id = int(target_text)
        except Exception:
            return False

        root = self._find_saved_current_result_root(target_id)
        if root is None:
            root = Path(self.params.P.result_dir)

        zp_path = step10_current_global_zp_path(root, target_id)
        if not zp_path.exists():
            return False

        try:
            params_df = pd.read_csv(zp_path)
        except Exception as e:
            self.log(f"[LOAD] Failed to recover global ZP: {e}")
            if not silent:
                self._tell_user("warn", "Detrend", f"Global ZP 로드 실패:\n{e}")
            return False

        if params_df.empty:
            self.log(f"[LOAD] Saved global ZP is empty: {zp_path.name}")
            return False

        self.params_df = params_df

        mean_path = step10_current_global_mean_path(root, target_id)
        if mean_path.exists():
            try:
                self.global_mean_df = pd.read_csv(mean_path)
            except Exception as e:
                self.log(f"[LOAD] Failed to recover global mean: {e}")

        diag_path = step10_current_global_diag_path(root, target_id)
        if diag_path.exists():
            try:
                self.global_diagnostics = json.loads(diag_path.read_text(encoding="utf-8"))
            except Exception as e:
                self.log(f"[LOAD] Failed to recover global diagnostics: {e}")

        self.log(f"[LOAD] Recovered global ZP from saved current result: {zp_path.name} ({len(self.params_df)} rows)")
        return True

    def _filter_linestyle(self, fkey: str) -> str:
        mapping = {
            "g": "-",
            "r": "--",
            "i": ":",
            "b": "-.",
            "v": (0, (3, 1, 1, 1)),
        }
        key = _normalize_filter_key(fkey)
        return mapping.get(key, "-")

    def _filter_color(self, fkey: str) -> str:
        mapping = {
            "g": "#2ca02c",
            "r": "#d62728",
            "i": "#9467bd",
            "b": "#1f77b4",
            "v": "#2ca02c",
            "z": "#8c564b",
        }
        key = _normalize_filter_key(fkey)
        return mapping.get(key, "#1f77b4")

    def _phase_xy(self, df: pd.DataFrame, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        if self.x_axis_mode != "phase" or self.phase_period <= 0 or "JD" not in df.columns:
            x = df["JD"].to_numpy(float) if "JD" in df.columns else np.arange(len(df), dtype=float)
            return x, y
        jd = pd.to_numeric(df["JD"], errors="coerce").to_numpy(float)
        if not np.any(np.isfinite(jd)):
            return jd, y
        t0 = self._phase_reference_t0()
        if not np.isfinite(t0):
            return jd, y
        phase = ((jd - t0) / self.phase_period) % 1.0
        cycles = float(self.phase_cycles) if self.phase_cycles > 1 else 1.0
        full = int(cycles)
        frac = cycles - full
        xs: list[np.ndarray] = []
        ys: list[np.ndarray] = []
        if full < 1:
            full = 1
        for k in range(full):
            xs.append(phase + k)
            ys.append(y)
        if frac > 1e-6:
            mask = phase <= (frac + 1e-9)
            xs.append(phase[mask] + full)
            ys.append(y[mask])
        return (np.concatenate(xs), np.concatenate(ys)) if xs else (phase, y)

    def _update_plots(self):
        self.ax_raw.clear()
        self.ax_corr.clear()
        self.ax_diag.clear()

        if self.raw_df.empty:
            self._plot_redraw()
            return

        selected_dates = self._selected_dates()
        if selected_dates:
            dates = sorted(selected_dates)
        else:
            dates = sorted({str(d) for d in self.raw_df.get("date", []) if str(d)})
        date_colors = {}
        palette = [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
            "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
        ]
        for i, d in enumerate(dates):
            date_colors[d] = palette[i % len(palette)]

        filter_sel = "All"
        if hasattr(self, "filter_combo"):
            filter_sel = self._filter_selection() or "All"
        filter_key = "" if filter_sel == "All" else _normalize_filter_key(filter_sel)

        raw = self.raw_df.copy()
        if filter_key and "filter" in raw.columns:
            raw = raw[raw["filter"].astype(str).map(_normalize_filter_key) == filter_key]
        if dates:
            raw = raw[raw["date"].astype(str).isin(dates)]

        x_label = "JD"
        if self.x_axis_mode == "phase" and self.phase_period > 0 and "JD" in raw.columns:
            x_label = "Phase"
            if self.phase_cycles > 1:
                x_label = f"Phase (0-{self.phase_cycles:g})"

        if self.color_by == "Filter" and "filter" in raw.columns:
            for fkey, sub in raw.groupby(raw["filter"].astype(str).map(_normalize_filter_key)):
                y = sub["diff_mag_raw"].to_numpy(float)
                x, y = self._phase_xy(sub, y)
                m = np.isfinite(x) & np.isfinite(y)
                if np.any(m):
                    self.ax_raw.plot(
                        x[m], y[m], marker="o", linestyle="None",
                        color=self._filter_color(fkey), markersize=3, alpha=0.7, label=fkey or "all"
                    )
        else:
            for d in dates:
                sub = raw[raw["date"].astype(str) == str(d)]
                y = sub["diff_mag_raw"].to_numpy(float)
                x, y = self._phase_xy(sub, y)
                m = np.isfinite(x) & np.isfinite(y)
                if np.any(m):
                    self.ax_raw.plot(x[m], y[m], marker="o", linestyle="None", color=date_colors[d], markersize=3, alpha=0.7, label=d)

        self.ax_raw.set_title("Raw", fontsize=10)
        self.ax_raw.set_xlabel(x_label, fontsize=9)
        self.ax_raw.set_ylabel("Δmag", fontsize=9)
        self.ax_raw.invert_yaxis()
        self.ax_raw.grid(True, alpha=0.3)
        self.ax_raw.tick_params(labelsize=8)

        # Check star overlay on raw plot
        try:
            _rd = self.datasets[0][1] if self.datasets else Path(self.params.P.result_dir)
            _ck_id, _ck_df = _load_check_star_for_plot(Path(_rd), filter_key)
            if _ck_df is not None and not _ck_df.empty:
                _y_col = next((c for c in ["diff_mag_raw", "diff_mag", "mag"] if c in _ck_df.columns), None)
                if _y_col and "JD" in _ck_df.columns:
                    if filter_key and "filter" in _ck_df.columns:
                        _ck_df = _ck_df[_ck_df["filter"].astype(str).map(_normalize_filter_key) == filter_key].copy()
                    _cy = pd.to_numeric(_ck_df[_y_col], errors="coerce").to_numpy(float)
                    _cx, _cy = self._phase_xy(_ck_df, _cy)
                    _m = np.isfinite(_cx) & np.isfinite(_cy)
                    if _m.any():
                        _ck_label = f"Check ID {_ck_id}" if _ck_id is not None else "Check"
                        self.ax_raw.scatter(
                            _cx[_m], _cy[_m], s=8, color="#FFD700", alpha=0.5, zorder=2,
                            label=_ck_label, marker="^"
                        )
        except Exception:
            pass

        corr = self.corrected_df if not self.corrected_df.empty else raw
        if filter_key and "filter" in corr.columns:
            corr = corr[corr["filter"].astype(str).map(_normalize_filter_key) == filter_key]
        if dates and "date" in corr.columns:
            corr = corr[corr["date"].astype(str).isin(dates)]

        if self.color_by == "Filter" and "filter" in corr.columns:
            for fkey, sub in corr.groupby(corr["filter"].astype(str).map(_normalize_filter_key)):
                y = sub["diff_mag_corr"].to_numpy(float) if "diff_mag_corr" in sub.columns else sub["diff_mag_raw"].to_numpy(float)
                x, y = self._phase_xy(sub, y)
                m = np.isfinite(x) & np.isfinite(y)
                if np.any(m):
                    self.ax_corr.plot(
                        x[m], y[m], marker="o", linestyle="None",
                        color=self._filter_color(fkey), markersize=3, alpha=0.7, label=fkey or "all"
                    )
        else:
            if "date" in corr.columns:
                for d in dates:
                    sub = corr[corr["date"].astype(str) == str(d)]
                    y = sub["diff_mag_corr"].to_numpy(float) if "diff_mag_corr" in sub.columns else sub["diff_mag_raw"].to_numpy(float)
                    x, y = self._phase_xy(sub, y)
                    m = np.isfinite(x) & np.isfinite(y)
                    if np.any(m):
                        self.ax_corr.plot(x[m], y[m], marker="o", linestyle="None", color=date_colors[d], markersize=3, alpha=0.7, label=d)
            else:
                # No date column, plot all data
                y = corr["diff_mag_corr"].to_numpy(float) if "diff_mag_corr" in corr.columns else corr["diff_mag_raw"].to_numpy(float)
                x, y = self._phase_xy(corr, y)
                m = np.isfinite(x) & np.isfinite(y)
                if np.any(m):
                    self.ax_corr.plot(x[m], y[m], marker="o", linestyle="None", color="#1f77b4", markersize=3, alpha=0.7)

        self.ax_corr.set_title("Corrected", fontsize=10)
        self.ax_corr.set_xlabel(x_label, fontsize=9)
        self.ax_corr.set_ylabel("Δmag", fontsize=9)
        self.ax_corr.invert_yaxis()
        self.ax_corr.grid(True, alpha=0.3)
        self.ax_corr.tick_params(labelsize=8)

        if self.mode == "global":
            self._plot_global_diagnostics(dates, filter_key, date_colors)
            for ax in (self.ax_raw, self.ax_corr, self.ax_diag):
                handles, labels = ax.get_legend_handles_labels()
                if handles and len(handles) <= 8:
                    ax.legend(loc="best", fontsize=7, framealpha=0.8)
            self._plot_figure().tight_layout()
            self._plot_redraw()
            return

        # Diagnostic plot
        diag = raw
        has_airmass = "airmass" in diag.columns
        diag_x_col = "airmass" if has_airmass else ("JD" if "JD" in diag.columns else None)
        if self.color_by == "Filter" and "filter" in diag.columns:
            for fkey, sub in diag.groupby(diag["filter"].astype(str).map(_normalize_filter_key)):
                x = pd.to_numeric(sub[diag_x_col], errors="coerce").to_numpy(float) if diag_x_col else np.full(len(sub), np.nan)
                y = pd.to_numeric(sub["diff_mag_raw"], errors="coerce").to_numpy(float)
                m = np.isfinite(x) & np.isfinite(y)
                if np.any(m):
                    self.ax_diag.plot(
                        x[m], y[m], marker="o", linestyle="None",
                        color=self._filter_color(fkey), markersize=3, alpha=0.7, label=fkey or "all"
                    )
        else:
            for d in dates:
                sub = diag[diag["date"].astype(str) == str(d)]
                x = pd.to_numeric(sub[diag_x_col], errors="coerce").to_numpy(float) if diag_x_col else np.full(len(sub), np.nan)
                y = pd.to_numeric(sub["diff_mag_raw"], errors="coerce").to_numpy(float)
                m = np.isfinite(x) & np.isfinite(y)
                if np.any(m):
                    self.ax_diag.plot(x[m], y[m], marker="o", linestyle="None", color=date_colors[d], markersize=3, alpha=0.7, label=d)

        # Fit lines
        if not self.params_df.empty:
            for _, row in self.params_df.iterrows():
                date_val = str(row.get("date", ""))
                if dates and date_val not in dates:
                    continue
                fkey = _normalize_filter_key(row.get("filter", ""))
                if filter_key and fkey != filter_key:
                    continue
                sub = diag[diag["date"].astype(str) == date_val]
                if fkey and "filter" in sub.columns:
                    sub = sub[sub["filter"].astype(str).map(_normalize_filter_key) == fkey]
                xvals = pd.to_numeric(sub[diag_x_col], errors="coerce").to_numpy(float) if diag_x_col else np.array([])
                xvals = xvals[np.isfinite(xvals)]
                if xvals.size == 0:
                    continue
                xmin = float(np.nanmin(xvals))
                xmax = float(np.nanmax(xvals))
                if not np.isfinite(xmin) or not np.isfinite(xmax):
                    continue
                if xmin == xmax:
                    xmin -= 0.05
                    xmax += 0.05
                xline = np.linspace(xmin, xmax, 50)
                xfit = xline
                if self.mode == "color":
                    delta_c_const = self._delta_c_for_filter(fkey)
                    if not np.isfinite(delta_c_const):
                        continue
                    xfit = xline * float(delta_c_const)
                yline = float(row.get("zp_offset", 0.0)) + float(row.get("ext_slope", 0.0)) * xfit
                linestyle = "-" if filter_key else self._filter_linestyle(fkey)
                line_color = date_colors.get(date_val, _colors().TEXT) if self.color_by == "Date" else self._filter_color(fkey)
                self.ax_diag.plot(xline, yline, color=line_color, linestyle=linestyle, linewidth=1.5, alpha=0.9)

        diag_xlabel = diag_x_col if diag_x_col else "Index"
        self.ax_diag.set_title(f"Δmag vs {diag_xlabel} (diagnostics)", fontsize=10)
        self.ax_diag.set_xlabel(diag_xlabel, fontsize=9)
        self.ax_diag.set_ylabel("Δmag", fontsize=9)
        self.ax_diag.invert_yaxis()
        self.ax_diag.grid(True, alpha=0.3)
        self.ax_diag.tick_params(labelsize=8)

        # Legends (compact)
        for ax in (self.ax_raw, self.ax_corr, self.ax_diag):
            handles, labels = ax.get_legend_handles_labels()
            if handles and len(handles) <= 8:
                ax.legend(loc="best", fontsize=7, framealpha=0.8)

        self._plot_figure().tight_layout()
        self._plot_redraw()

    def _plot_global_diagnostics(self, dates, filter_key, date_colors) -> None:
        self.ax_diag.clear()
        diag_df = self._get_global_zp_df()
        if diag_df.empty:
            self.ax_diag.text(0.5, 0.5, "No global ZP available", ha="center", va="center")
            return

        ref_df = self.corrected_df if not self.corrected_df.empty else self.raw_df
        if not ref_df.empty and "time_id" in ref_df.columns:
            if "JD" in ref_df.columns:
                jd_map = ref_df.groupby("time_id")["JD"].median()
                diag_df["JD"] = diag_df["time_id"].map(jd_map)
            if "date" in ref_df.columns:
                date_map = ref_df.groupby("time_id")["date"].first()
                diag_df["date"] = diag_df["time_id"].map(date_map)

        if filter_key and "filter" in diag_df.columns:
            diag_df = diag_df[diag_df["filter"].astype(str).map(_normalize_filter_key) == filter_key]
        if dates and "date" in diag_df.columns:
            diag_df = diag_df[diag_df["date"].astype(str).isin(dates)]

        def _col(df, col):
            """Safely extract a numeric column as float array; NaN-filled if missing."""
            if col in df.columns:
                return pd.to_numeric(df[col], errors="coerce").to_numpy(float)
            return np.full(len(df), np.nan)

        # Prefer JD when it has usable values. Otherwise fall back to frame order,
        # because time_id is typically a filename-like string and cannot be plotted numerically.
        jd_vals = _col(diag_df, "JD")
        has_jd = np.isfinite(jd_vals).any()
        x_col = "JD"
        x_label = "JD"
        title = "Global Z_t vs JD"
        if not has_jd:
            order_map = None
            if not ref_df.empty and "time_id" in ref_df.columns:
                ordered_time_ids = pd.Index(ref_df["time_id"].astype(str)).drop_duplicates()
                order_map = {tid: i + 1 for i, tid in enumerate(ordered_time_ids)}
            if "time_id" in diag_df.columns:
                plot_x = None
                if order_map:
                    plot_x = pd.to_numeric(
                        diag_df["time_id"].astype(str).map(order_map),
                        errors="coerce",
                    )
                if plot_x is None or not np.isfinite(plot_x.to_numpy(float)).any():
                    plot_x = pd.Series(
                        np.arange(1, len(diag_df) + 1, dtype=float),
                        index=diag_df.index,
                    )
                diag_df["_plot_x"] = plot_x.to_numpy(float)
            else:
                diag_df["_plot_x"] = np.arange(1, len(diag_df) + 1, dtype=float)
            x_col = "_plot_x"
            x_label = "Frame Index"
            title = "Global Z_t vs Frame"

        if self.color_by == "Filter" and "filter" in diag_df.columns:
            for fkey, sub in diag_df.groupby(diag_df["filter"].astype(str).map(_normalize_filter_key)):
                xs = _col(sub, x_col)
                ys = _col(sub, "Z")
                m = np.isfinite(xs) & np.isfinite(ys)
                if np.any(m):
                    self.ax_diag.plot(
                        xs[m], ys[m], marker="o", linestyle="None",
                        color=self._filter_color(fkey), markersize=3, alpha=0.7, label=fkey or "all"
                    )
        elif "date" in diag_df.columns:
            for d in sorted(diag_df["date"].astype(str).unique().tolist()):
                sub = diag_df[diag_df["date"].astype(str) == str(d)]
                xs = _col(sub, x_col)
                ys = _col(sub, "Z")
                m = np.isfinite(xs) & np.isfinite(ys)
                if np.any(m):
                    self.ax_diag.plot(xs[m], ys[m], marker="o", linestyle="None", color=date_colors.get(d, _colors().TEXT), markersize=3, alpha=0.7, label=d)
        else:
            xs = _col(diag_df, x_col)
            ys = _col(diag_df, "Z")
            m = np.isfinite(xs) & np.isfinite(ys)
            if np.any(m):
                self.ax_diag.plot(xs[m], ys[m], marker="o", linestyle="None", color=_colors().TEXT, markersize=3, alpha=0.7)

        if not self.ax_diag.lines:
            self.ax_diag.text(0.5, 0.5, "No plottable global ZP rows", ha="center", va="center")

        self.ax_diag.set_title(title, fontsize=10)
        self.ax_diag.set_xlabel(x_label, fontsize=9)
        self.ax_diag.set_ylabel("Z_t", fontsize=9)
        self.ax_diag.grid(True, alpha=0.3)
        self.ax_diag.tick_params(labelsize=8)


class HeadlessDetrendRunner(DetrendRunner):
    """A batch runner for the same detrend the window performs.

    The window gathers its state across a constructor and a dozen widget
    callbacks; this gathers the equivalent in one place. There is no second
    implementation of the correction here — `load_raw_data` and
    `fit_and_apply` are inherited, which is the point of having moved them.
    """

    def __init__(self, params, result_dirs, *, logger=None, project_state=None,
                 file_manager=None, target_id: int = 0, **settings):
        self.params = params
        self.project_state = project_state
        self.file_manager = file_manager
        self.target_id = int(target_id or 0)
        self.runtime_mode = True

        dirs = [Path(d) for d in result_dirs]
        self.datasets = [(d.name, d) for d in dirs]

        # The names the calculation actually fills. Inventing `corr_df` and
        # `fit_params` here made the step report "corrected no points" on a run
        # that had corrected all 364 of them — the same mistake as asking step 6
        # for `master_sources.csv`.
        self.raw_df = pd.DataFrame()
        self.corrected_df = pd.DataFrame()
        self.params_df = pd.DataFrame()
        self.summary_text = ""
        self._logger = logger

        # Only settings the class already declares; a typo must not become a
        # silently ignored option.
        for key, value in settings.items():
            if not hasattr(DetrendRunner, key):
                raise TypeError(f"unknown detrend setting: {key}")
            setattr(self, key, value)

    def log(self, message: str) -> None:
        if self._logger is not None:
            self._logger(message)

    def save_plot(self, path) -> "Path | None":
        """Write the figure `_update_plots` drew. The window saves its canvas."""
        fig = getattr(self, "_headless_fig", None)
        if fig is None or not fig.axes:
            return None
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=150, bbox_inches="tight")
        return path

