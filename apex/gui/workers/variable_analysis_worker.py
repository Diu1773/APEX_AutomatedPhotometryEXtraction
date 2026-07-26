"""Cancelable Qt worker for the headless variable-analysis service."""

from __future__ import annotations

from PyQt5.QtCore import QThread, pyqtSignal

from apex.analysis.light_curve.variable_analysis_cache import (
    load_cached_result,
    store_cached_result,
)
from apex.analysis.light_curve.variable_analysis_contract import VariableAnalysisRequest
from apex.analysis.light_curve.variable_analysis_service import run_variable_analysis


class VariableAnalysisWorker(QThread):
    progress = pyqtSignal(str, int)
    finished = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, request: VariableAnalysisRequest, parent=None):
        super().__init__(parent)
        self.request = request
        self._stop_requested = False

    def stop(self) -> None:
        self._stop_requested = True

    def _cancelled(self) -> bool:
        return bool(self._stop_requested)

    def _progress(self, stage: str, fraction: float) -> None:
        self.progress.emit(str(stage), int(round(float(fraction) * 100.0)))

    def run(self) -> None:
        try:
            self._progress("Checking analysis cache", 0.01)
            cached = load_cached_result(self.request)
            if cached is not None and not self._cancelled():
                self._progress("Loaded cached advanced analysis", 1.0)
                self.finished.emit(cached)
                return
            result = run_variable_analysis(
                self.request,
                progress_callback=self._progress,
                cancel_requested=self._cancelled,
            )
            if result.status == "COMPLETE":
                try:
                    store_cached_result(self.request, result)
                except OSError:
                    pass
            self.finished.emit(result)
        except Exception as exc:
            self.failed.emit(str(exc))
