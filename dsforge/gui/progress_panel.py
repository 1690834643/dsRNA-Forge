"""
进度监控面板
- 进度条
- 实时日志
- 取消按钮
"""

from PyQt6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QProgressBar,
    QTextEdit,
    QPushButton,
    QGroupBox,
)
from PyQt6.QtCore import pyqtSignal


class ProgressPanel(QGroupBox):
    """进度监控面板"""

    cancel_requested = pyqtSignal()

    def __init__(self):
        super().__init__("Progress")
        self._setup_ui()
        self.reset()

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(8)

        # 进度条
        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 100)
        layout.addWidget(self.progress_bar)

        # 状态标签
        self.status_label = QLabel("Idle")
        layout.addWidget(self.status_label)

        # 日志输出
        self.log_edit = QTextEdit()
        self.log_edit.setMaximumHeight(120)
        self.log_edit.setReadOnly(True)
        layout.addWidget(self.log_edit)

        # 统计
        self.stats_label = QLabel("")
        layout.addWidget(self.stats_label)

    def reset(self):
        """重置状态"""
        self.progress_bar.setValue(0)
        self.status_label.setText("Idle")
        self.log_edit.clear()
        self.stats_label.setText("")

    def set_progress(self, percent: float):
        """设置进度百分比"""
        self.progress_bar.setValue(int(percent))

    def set_status(self, message: str):
        """设置状态文本"""
        self.status_label.setText(message)

    def log(self, message: str):
        """添加日志"""
        self.log_edit.append(message)

    def set_stats(self, stats: str):
        """设置统计信息"""
        self.stats_label.setText(stats)
