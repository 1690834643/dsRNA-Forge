"""
Saved transcriptome cache management dialog.
"""

from PyQt6.QtWidgets import (
    QDialog,
    QVBoxLayout,
    QHBoxLayout,
    QTableWidget,
    QTableWidgetItem,
    QPushButton,
    QInputDialog,
    QMessageBox,
    QHeaderView,
)

from dsforge.core.sequence import TranscriptomeIndex


class CacheManagerDialog(QDialog):
    """Simple cache manager for saved transcriptomes."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Manage Transcriptome Cache")
        self.resize(760, 420)
        self.entries = []
        self._setup_ui()
        self.refresh()

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        self.table = QTableWidget()
        self.table.setColumnCount(5)
        self.table.setHorizontalHeaderLabels(["Name", "Sequences", "Total nt", "Source", "Key"])
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table)

        buttons = QHBoxLayout()
        self.rename_btn = QPushButton("Rename")
        self.delete_btn = QPushButton("Delete Selected")
        self.clear_btn = QPushButton("Clear All")
        self.refresh_btn = QPushButton("Refresh")
        self.close_btn = QPushButton("Close")
        for btn in [self.rename_btn, self.delete_btn, self.clear_btn, self.refresh_btn, self.close_btn]:
            buttons.addWidget(btn)
        layout.addLayout(buttons)

        self.rename_btn.clicked.connect(self.rename_selected)
        self.delete_btn.clicked.connect(self.delete_selected)
        self.clear_btn.clicked.connect(self.clear_all)
        self.refresh_btn.clicked.connect(self.refresh)
        self.close_btn.clicked.connect(self.accept)

    def refresh(self):
        self.entries = TranscriptomeIndex.list_saved()
        self.table.setRowCount(0)
        for entry in self.entries:
            row = self.table.rowCount()
            self.table.insertRow(row)
            stats = entry.get("stats") or {}
            values = [
                entry.get("name", ""),
                str(stats.get("num_sequences", "")),
                str(stats.get("total_nt", "")),
                entry.get("source_file", ""),
                entry.get("key", ""),
            ]
            for col, value in enumerate(values):
                self.table.setItem(row, col, QTableWidgetItem(value))

    def _selected_entry(self):
        selected = self.table.selectionModel().selectedRows()
        if not selected:
            return None
        row = selected[0].row()
        if row < 0 or row >= len(self.entries):
            return None
        return self.entries[row]

    def rename_selected(self):
        entry = self._selected_entry()
        if entry is None:
            return
        new_name, ok = QInputDialog.getText(self, "Rename Transcriptome", "New name:", text=entry.get("name", ""))
        if ok and new_name.strip():
            TranscriptomeIndex.rename_saved(entry["key"], new_name.strip())
            self.refresh()

    def delete_selected(self):
        entry = self._selected_entry()
        if entry is None:
            return
        confirm = QMessageBox.question(
            self,
            "Delete Saved Transcriptome",
            f"Delete '{entry.get('name', entry.get('key'))}' and its cached indexes?",
        )
        if confirm == QMessageBox.StandardButton.Yes:
            TranscriptomeIndex.delete_saved(entry["key"], delete_cache=True)
            self.refresh()

    def clear_all(self):
        confirm = QMessageBox.question(
            self,
            "Clear Cache",
            "Delete all saved transcriptomes and local indexes?",
        )
        if confirm == QMessageBox.StandardButton.Yes:
            TranscriptomeIndex.clear_cache()
            self.refresh()
