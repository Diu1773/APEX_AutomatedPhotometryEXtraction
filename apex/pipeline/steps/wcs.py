"""Step 5 (headless): WCS plate solving.

Delegates to the Qt-free :func:`apex.analysis.wcs_solve.run_wcs_solve`, the same
orchestration the GUI workers now use. Reads the file list from Step 1's
``selection.json`` and the detections from Step 4, runs the configured solver
(default: the internal Python engine), and writes ``step5_wcs/`` (per-frame FITS
WCS headers, sidecar JSONs, and ``wcs_solve_summary.csv`` / ``frame_wcs_qc.csv``).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import List

from apex.pipeline.base import PipelineStep, StepResult, StepStatus
from apex.pipeline.context import RunContext
from apex.utils.step_paths import (
    step1_dir,
    step2_cropped_dir,
    step4_dir,
    step5_wcs_dir,
    crop_is_active,
)


def _selection_path(result_dir: Path) -> Path:
    return step1_dir(result_dir) / "selection.json"


def _resolve_use_cropped(result_dir: Path) -> bool:
    """Mirror the GUI rule: crop active AND cropped dir holds FITS files."""
    cropped_dir = step2_cropped_dir(result_dir)
    if not crop_is_active(result_dir):
        return False
    if not cropped_dir.exists():
        return False
    return bool(list(cropped_dir.glob("*.fit*")))


def _frames_carrying_a_header_wcs(
    file_list: List[str], data_dir: Path, result_dir: Path, use_cropped: bool
) -> int:
    """헤더에 하늘 좌표가 붙어 있는 프레임의 수.

    **「한 장도 못 풀었다」와 「이어서 못 간다」는 다른 말이다.** 6 단계는 5 단계의
    산출물이 없으면 프레임 헤더의 측성 해를 읽어 쓴다(`refbuild.py` 의
    `w.has_celestial` 갈래). 관측소가 이미 푼 프레임 — LCO 의 BANZAI 처리본 같은 것
    — 은 그 길로 정상적으로 끝까지 간다.

    그러니 풀이 엔진이 한 장도 못 풀었을 때 실패로 볼지 말지는 **헤더에 쓸 수 있는
    좌표가 남아 있느냐**로 갈린다. 헤더만 읽으므로 화소는 건드리지 않는다.
    """
    from astropy.wcs import WCS

    from apex.utils.io_utils import read_fits_header

    roots = [step2_cropped_dir(result_dir), data_dir] if use_cropped else [data_dir]
    n = 0
    for name in file_list:
        for root in roots:
            path = Path(root) / str(name)
            if not path.exists():
                continue
            try:
                if WCS(read_fits_header(path), relax=True).has_celestial:
                    n += 1
            except Exception:  # noqa: BLE001 - 못 읽는 헤더는 없는 것으로 친다
                pass
            break
    return n


class WcsStep(PipelineStep):
    index = 5
    key = "wcs"
    title = "WCS plate solving"
    interactive = False

    def required_inputs(self, ctx: RunContext) -> List[Path]:
        return [step4_dir(ctx.result_dir)]

    def outputs(self, ctx: RunContext) -> List[Path]:
        return [step5_wcs_dir(ctx.result_dir)]

    def is_complete(self, ctx: RunContext) -> bool:
        out = step5_wcs_dir(ctx.result_dir)
        if not out.exists():
            return False
        summary = out / "wcs_solve_summary.csv"
        return summary.exists() and summary.stat().st_size > 0

    def run(self, ctx: RunContext) -> StepResult:
        sel_path = _selection_path(ctx.result_dir)
        if not sel_path.exists():
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message=f"missing required input: {sel_path}",
            )

        selection = json.loads(sel_path.read_text(encoding="utf-8"))
        file_list = list(selection.get("filenames", []))
        if not file_list:
            return StepResult(
                index=self.index, key=self.key, status=StepStatus.BLOCKED,
                message="selection.json contains no filenames",
            )

        use_cropped = _resolve_use_cropped(ctx.result_dir)

        from apex.analysis.wcs_solve import run_wcs_solve, resolve_wcs_engine

        engine = resolve_wcs_engine(ctx.params)

        summary = run_wcs_solve(
            file_list,
            ctx.params,
            ctx.data_dir,
            ctx.result_dir,
            ctx.params.P.cache_dir,
            engine=engine,
            use_cropped=use_cropped,
            logger=ctx.logger,
        )

        out_dir = step5_wcs_dir(ctx.result_dir)
        summary = summary if isinstance(summary, dict) else {}
        n_total = int(summary.get("total", 0) or 0)
        n_ok = int(summary.get("ok", 0) or 0)
        n_qc = int(summary.get("wcs_qc_pass", 0) or 0)
        msg = (
            f"{n_ok}/{len(file_list)} frames solved ({engine}); "
            f"WCS-QC pass={n_qc}"
        )
        # Step 5 had no figure at all — not even in the window. The numbers
        # were in frame_wcs_qc.csv and nothing looked at them together.
        try:
            from apex.analysis.wcs_qc_plots import export_wcs_qc
            qc_figures = export_wcs_qc(ctx.result_dir, ctx.params)
            if qc_figures:
                msg += f'; {len(qc_figures)} QC figures'
        except Exception:  # noqa: BLE001 - QC must not fail a finished solve
            if ctx.logger is not None:
                ctx.logger.exception('Could not write Step 5 QC figure')

        # **한 장도 못 풀었는데 초록으로 끝나면 안 된다.** 여태 이 단계는 몇 장을
        # 풀었든 OK 를 돌려주었다. `0/30 frames solved` 라고 적기는 하지만 상태가
        # 초록이라 자동 실행은 못 알아챈다. 게다가 한 장도 못 풀면
        # `wcs_solve_summary.csv` 가 아예 안 쓰이므로, **스스로 `is_complete` 에서
        # 「완료 아님」이라고 판정할 산출물을 안 내고서 OK 를 돌려주던 셈이다.**
        #
        # 그렇다고 무조건 실패로 볼 수도 없다. 6 단계는 프레임 헤더의 측성 해를
        # 읽어 쓸 수 있고, 관측소가 이미 푼 프레임(LCO 의 BANZAI 처리본 등)은 그
        # 길로 끝까지 간다. 그래서 **이어서 갈 수 있느냐**로 가른다.
        if n_ok == 0 and file_list:
            n_hdr = _frames_carrying_a_header_wcs(
                file_list, Path(ctx.data_dir), ctx.result_dir, use_cropped)
            if n_hdr:
                msg += (f"; 한 장도 못 풀었다 — 헤더에 측성 해가 있는 "
                        f"{n_hdr}/{len(file_list)} 장으로 이어서 간다")
            else:
                return StepResult(
                    index=self.index, key=self.key, status=StepStatus.FAILED,
                    message=(msg + "; 헤더에도 측성 해가 없어 이어서 갈 수 없다"),
                    outputs=[str(out_dir)],
                )

        return StepResult(
            index=self.index, key=self.key, status=StepStatus.OK,
            message=msg, outputs=[str(out_dir)],
        )
