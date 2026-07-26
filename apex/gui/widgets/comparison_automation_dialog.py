"""Progress window for automated comparison and check-star selection."""

from __future__ import annotations

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
)

from apex.gui.theme import Tokens, style_button


class ComparisonAutomationDialog(QDialog):
    cancel_requested = pyqtSignal()

    def __init__(self, filters: list[str], total_steps: int, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Comparison Automation")
        self.setWindowModality(Qt.WindowModal)
        self.setMinimumSize(620, 400)
        self.resize(700, 460)
        self._running = True

        layout = QVBoxLayout(self)
        layout.setContentsMargins(Tokens.MARGIN, Tokens.MARGIN, Tokens.MARGIN, Tokens.MARGIN)
        layout.setSpacing(Tokens.S3)

        self.title_label = QLabel("Automated comparison and check-star selection")
        self.title_label.setProperty("role", "title")
        layout.addWidget(self.title_label)

        self.scope_label = QLabel("Filters: " + ", ".join(str(value) for value in filters))
        self.scope_label.setProperty("role", "caption")
        layout.addWidget(self.scope_label)

        self.stage_label = QLabel("Preparing inputs")
        self.stage_label.setProperty("role", "subtitle")
        layout.addWidget(self.stage_label)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, max(1, int(total_steps)))
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(True)
        layout.addWidget(self.progress_bar)

        self.detail_label = QLabel("")
        self.detail_label.setProperty("role", "caption")
        self.detail_label.setWordWrap(True)
        layout.addWidget(self.detail_label)

        self.log_text = QTextEdit()
        self.log_text.setObjectName("Log")
        self.log_text.setReadOnly(True)
        layout.addWidget(self.log_text, 1)

        self.summary_label = QLabel("")
        self.summary_label.setWordWrap(True)
        layout.addWidget(self.summary_label)

        actions = QHBoxLayout()
        actions.addStretch()
        self.action_button = QPushButton("Cancel")
        style_button(self.action_button, height=Tokens.H_BUTTON)
        self.action_button.clicked.connect(self._handle_action)
        actions.addWidget(self.action_button)
        layout.addLayout(actions)

    def _handle_action(self) -> None:
        if self._running:
            self.action_button.setEnabled(False)
            self.stage_label.setText("Canceling")
            self.cancel_requested.emit()
        else:
            self.accept()

    def set_progress(
        self,
        value: int,
        total: int,
        stage: str,
        detail: str = "",
    ) -> None:
        self.progress_bar.setRange(0, max(1, int(total)))
        self.progress_bar.setValue(max(0, min(int(value), int(total))))
        self.stage_label.setText(str(stage))
        self.detail_label.setText(str(detail))
        if detail:
            self.log_text.append(f"{stage}: {detail}")

    def append_log(self, message: str) -> None:
        if message:
            self.log_text.append(str(message))

    def finish(self, summary: str, *, success: bool) -> None:
        self._running = False
        self.progress_bar.setValue(self.progress_bar.maximum())
        self.stage_label.setText("Complete" if success else "Completed with issues")
        self.summary_label.setText(str(summary))
        self.action_button.setText("Close")
        self.action_button.setEnabled(True)
        style_button(
            self.action_button,
            "primary" if success else None,
            height=Tokens.H_BUTTON,
        )

    def reject(self) -> None:
        if self._running:
            self._handle_action()
            return
        super().reject()
