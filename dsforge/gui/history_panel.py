"""
历史任务面板
- 显示过去的设计任务列表
- 点击重载结果到表格
"""

from PyQt6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QPushButton,
    QHeaderView,
    QAbstractItemView,
    QLabel,
    QGroupBox,
)
from PyQt6.QtCore import pyqtSignal


class HistoryPanel(QGroupBox):
    """历史任务面板"""

    task_selected = pyqtSignal(int)  # 发送 task_id

    def __init__(self):
        super().__init__("History")
        self._setup_ui()

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(8)

        self.table = QTableWidget()
        self.table.setColumnCount(5)
        self.table.setHorizontalHeaderLabels(
            ["ID", "Mode", "Target", "Status", "Created"]
        )
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.itemSelectionChanged.connect(self._on_selection_changed)
        layout.addWidget(self.table)

        btn_layout = QHBoxLayout()
        self.refresh_btn = QPushButton("Refresh")
        self.refresh_btn.clicked.connect(self.load_history)
        self.load_btn = QPushButton("Load Selected")
        self.load_btn.setEnabled(False)
        self.load_btn.clicked.connect(self._on_load)

        btn_layout.addWidget(self.refresh_btn)
        btn_layout.addWidget(self.load_btn)
        layout.addLayout(btn_layout)

        self.status_label = QLabel("No history loaded")
        layout.addWidget(self.status_label)

    def load_history(self, tasks: list = None):
        """加载历史任务"""
        if tasks is None:
            self.table.setRowCount(0)
            self.status_label.setText("Click Refresh to load from database")
            return

        self.table.setRowCount(0)
        for task in tasks:
            row = self.table.rowCount()
            self.table.insertRow(row)
            self.table.setItem(row, 0, QTableWidgetItem(str(task.get("id", ""))))
            self.table.setItem(row, 1, QTableWidgetItem(task.get("mode", "")))
            self.table.setItem(row, 2, QTableWidgetItem(task.get("target_seq_id", "") or "N/A"))
            self.table.setItem(row, 3, QTableWidgetItem(task.get("status", "")))
            self.table.setItem(row, 4, QTableWidgetItem(str(task.get("created_at", ""))))

        self.status_label.setText(f"{len(tasks)} task(s) in history")

    def _on_selection_changed(self):
        """选择变化"""
        self.load_btn.setEnabled(len(self.table.selectedItems()) > 0)

    def _on_load(self):
        """加载选中任务"""
        selected = self.table.selectedItems()
        if selected:
            task_id = int(self.table.item(selected[0].row(), 0).text())
            self.task_selected.emit(task_id)
