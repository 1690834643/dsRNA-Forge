"""
转录组管理面板
- FASTA 上传
- 索引状态显示
- 转录组统计
- 目标序列选择
"""

from pathlib import Path, PureWindowsPath

from PyQt6.QtWidgets import (
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QFileDialog,
    QLineEdit,
    QTextEdit,
    QGroupBox,
    QComboBox,
)
from PyQt6.QtCore import pyqtSignal

from dsforge.core.sequence import TranscriptomeIndex
from dsforge.core.sequence import load_first_fasta_record, make_safe_sequence_id, normalize_sequence


def _path_stem_for_label(path: str) -> str:
    if "\\" in path:
        return PureWindowsPath(path).stem
    return Path(path).stem


class TranscriptPanel(QGroupBox):
    """转录组管理面板"""

    transcriptome_loaded = pyqtSignal(object, dict)  # TranscriptomeIndex, stats
    transcriptome_cleared = pyqtSignal()
    manage_cache_requested = pyqtSignal()

    def __init__(self):
        super().__init__("1. Load transcriptome and choose target")
        self.transcriptome = None
        self._all_target_ids = []
        self._uploaded_target = None
        self._saved_entries = []
        self._backgrounds = []
        self._setup_ui()

    @staticmethod
    def background_label_from_path(path: str) -> str:
        return make_safe_sequence_id(_path_stem_for_label(path), "background")

    def _setup_ui(self):
        layout = QVBoxLayout(self)
        layout.setSpacing(8)

        # 已入库转录组
        saved_layout = QHBoxLayout()
        saved_layout.addWidget(QLabel("Saved transcriptomes:"))
        self.saved_combo = QComboBox()
        self.load_saved_btn = QPushButton("Load Saved")
        self.manage_cache_btn = QPushButton("Manage Cache")
        self.load_saved_btn.clicked.connect(self._on_load_saved)
        self.manage_cache_btn.clicked.connect(self.manage_cache_requested.emit)
        saved_layout.addWidget(self.saved_combo)
        saved_layout.addWidget(self.load_saved_btn)
        saved_layout.addWidget(self.manage_cache_btn)
        layout.addLayout(saved_layout)

        # 文件选择行
        file_layout = QHBoxLayout()
        self.path_edit = QLineEdit()
        self.path_edit.setPlaceholderText("Choose transcriptome FASTA...")
        self.path_edit.setReadOnly(True)

        self.browse_btn = QPushButton("Browse...")
        self.browse_btn.clicked.connect(self._on_browse)

        file_layout.addWidget(self.path_edit)
        file_layout.addWidget(self.browse_btn)
        layout.addLayout(file_layout)

        # 目标序列来源
        source_layout = QHBoxLayout()
        source_layout.addWidget(QLabel("Target source:"))
        self.target_source_combo = QComboBox()
        self.target_source_combo.addItems([
            "Search in transcriptome",
            "Paste target sequence",
            "Upload target FASTA",
        ])
        self.target_source_combo.currentIndexChanged.connect(self._on_target_source_changed)
        source_layout.addWidget(self.target_source_combo)
        layout.addLayout(source_layout)

        # 目标序列选择
        self.search_edit = QLineEdit()
        self.search_edit.setPlaceholderText("Type a gene/transcript ID to filter the target list...")
        self.search_edit.textChanged.connect(self._filter_targets)
        layout.addWidget(self.search_edit)

        target_layout = QHBoxLayout()
        target_layout.addWidget(QLabel("Target:"))
        self.target_combo = QComboBox()
        self.target_combo.setEnabled(False)
        target_layout.addWidget(self.target_combo)
        layout.addLayout(target_layout)

        # 自定义目标序列
        self.custom_seq_edit = QTextEdit()
        self.custom_seq_edit.setMaximumHeight(90)
        self.custom_seq_edit.setPlaceholderText("Paste a target sequence here. DNA/RNA and FASTA text are both OK.")
        layout.addWidget(self.custom_seq_edit)

        target_file_layout = QHBoxLayout()
        self.target_path_edit = QLineEdit()
        self.target_path_edit.setReadOnly(True)
        self.target_path_edit.setPlaceholderText("Choose a single-target FASTA...")
        self.target_browse_btn = QPushButton("Upload Target FASTA...")
        self.target_browse_btn.clicked.connect(self._on_target_browse)
        target_file_layout.addWidget(self.target_path_edit)
        target_file_layout.addWidget(self.target_browse_btn)
        layout.addLayout(target_file_layout)

        background_layout = QHBoxLayout()
        self.background_label = QLabel("Extra backgrounds: 0")
        self.add_background_btn = QPushButton("Add Background FASTA...")
        self.clear_background_btn = QPushButton("Clear Backgrounds")
        self.add_background_btn.clicked.connect(self._on_add_background)
        self.clear_background_btn.clicked.connect(self._on_clear_backgrounds)
        background_layout.addWidget(self.background_label)
        background_layout.addWidget(self.add_background_btn)
        background_layout.addWidget(self.clear_background_btn)
        layout.addLayout(background_layout)

        # 统计信息
        self.stats_label = QLabel("No transcriptome loaded")
        layout.addWidget(self.stats_label)

        # 日志/状态显示
        self.log_edit = QTextEdit()
        self.log_edit.setMaximumHeight(80)
        self.log_edit.setReadOnly(True)
        self.log_edit.setPlaceholderText("Loading messages will appear here.")
        layout.addWidget(self.log_edit)

        # 操作按钮
        btn_layout = QHBoxLayout()
        self.load_btn = QPushButton("Load & Index")
        self.load_btn.setEnabled(False)
        self.load_btn.clicked.connect(self._on_load)

        self.clear_btn = QPushButton("Clear")
        self.clear_btn.clicked.connect(self._on_clear)

        btn_layout.addWidget(self.load_btn)
        btn_layout.addWidget(self.clear_btn)
        layout.addLayout(btn_layout)
        self._refresh_saved_transcriptomes()
        self._on_target_source_changed(0)

    def _on_browse(self):
        """浏览文件"""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Transcriptome FASTA",
            "",
            "FASTA Files (*.fa *.fasta *.fna);;All Files (*)",
        )
        if path:
            self.path_edit.setText(path)
            self.load_btn.setEnabled(True)
            self._log(f"Selected: {path}")

    def _on_load(self):
        """加载并索引"""
        path = self.path_edit.text()
        if not path:
            return

        self._log(f"Loading transcriptome from {path}...")
        try:
            self.transcriptome = TranscriptomeIndex()
            self.transcriptome.load_fasta(path)
            stats = self.transcriptome.get_stats()
            self.set_transcriptome(self.transcriptome)

            self._log(f"Loaded {stats['num_sequences']} sequences, {stats['total_nt']} nt")
            self._refresh_saved_transcriptomes()
            self.transcriptome_loaded.emit(self.transcriptome, stats)

        except Exception as e:
            self._log(f"ERROR: {e}")
            self.transcriptome = None
            self.transcriptome_cleared.emit()

    def _on_clear(self):
        """清除"""
        self.transcriptome = None
        self._all_target_ids = []
        self._uploaded_target = None
        self._backgrounds = []
        self.path_edit.clear()
        self.load_btn.setEnabled(False)
        self.target_combo.clear()
        self.target_combo.setEnabled(False)
        self.search_edit.clear()
        self.custom_seq_edit.clear()
        self.target_path_edit.clear()
        self.stats_label.setText("No transcriptome loaded")
        self.log_edit.clear()
        self._update_background_label()
        self._refresh_saved_transcriptomes()
        self.transcriptome_cleared.emit()

    def get_target_id(self):
        """获取选中的目标序列 ID"""
        selection = self.get_target_selection()
        if selection is None:
            return None
        return selection["id"]

    def get_target_selection(self):
        """Return selected target metadata for MainWindow."""
        source = self.target_source_combo.currentIndex()
        if source == 0:
            if self.target_combo.count() == 0:
                return None
            return {
                "source": "transcriptome",
                "id": self.target_combo.currentText(),
                "sequence": None,
            }

        if source == 1:
            text = self.custom_seq_edit.toPlainText().strip()
            if not text:
                return None
            if text.startswith(">"):
                import tempfile
                from pathlib import Path

                with tempfile.NamedTemporaryFile("w", delete=False, encoding="utf-8", suffix=".fa") as tmp:
                    tmp.write(text)
                    tmp_path = tmp.name
                try:
                    record_id, sequence = load_first_fasta_record(tmp_path)
                finally:
                    Path(tmp_path).unlink(missing_ok=True)
            else:
                record_id = "custom_target"
                sequence = normalize_sequence(text)
            return {
                "source": "custom",
                "id": make_safe_sequence_id(record_id, "custom_target"),
                "sequence": sequence,
            }

        if self._uploaded_target is None:
            return None
        record_id, sequence = self._uploaded_target
        return {
            "source": "custom",
            "id": make_safe_sequence_id(record_id, "custom_target"),
            "sequence": sequence,
        }

    def add_background_index(self, label: str, transcriptome: TranscriptomeIndex):
        """Add an already loaded extra off-target background."""
        self._backgrounds.append((make_safe_sequence_id(label, "background"), transcriptome))
        self._update_background_label()

    def get_background_indexes(self):
        """Return loaded extra off-target backgrounds."""
        return list(self._backgrounds)

    def _update_background_label(self):
        self.background_label.setText(f"Extra backgrounds: {len(self._backgrounds)}")

    def set_transcriptome(self, transcriptome: TranscriptomeIndex):
        """Populate panel state from an already-loaded transcriptome."""
        self.transcriptome = transcriptome
        stats = transcriptome.get_stats()
        self._all_target_ids = transcriptome.list_ids()
        self._filter_targets(self.search_edit.text())
        self.target_combo.setEnabled(bool(self._all_target_ids))
        self.path_edit.setText(transcriptome.source_file or "")
        self.load_btn.setEnabled(False)
        self.stats_label.setText(
            f"Sequences: {stats['num_sequences']} | "
            f"Total: {stats['total_nt']} nt | "
            f"Avg: {stats['avg_length']:.0f} nt | "
            f"GC: {stats['gc_content']:.1f}%"
        )

    def apply_target_selection(self, target: dict):
        """Restore target controls from a saved project target payload."""
        if not isinstance(target, dict):
            return
        if target.get("source") == "transcriptome":
            self.target_source_combo.setCurrentIndex(0)
            target_id = target.get("id", "")
            if self.search_edit.text():
                self.search_edit.clear()
            idx = self.target_combo.findText(target_id)
            if idx >= 0:
                self.target_combo.setCurrentIndex(idx)
            return
        if target.get("source") == "custom":
            self.target_source_combo.setCurrentIndex(1)
            self.custom_seq_edit.setPlainText(target.get("sequence") or "")

    def _on_add_background(self):
        """Load an extra transcriptome/background for off-target scanning."""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Extra Off-target Background FASTA",
            "",
            "FASTA Files (*.fa *.fasta *.fna);;All Files (*)",
        )
        if not path:
            return
        try:
            background = TranscriptomeIndex()
            background.load_fasta(path)
            label = self.background_label_from_path(path)
            self.add_background_index(label, background)
            stats = background.get_stats()
            self._log(f"Added background {label}: {stats.get('num_sequences', 0)} sequences")
        except Exception as e:
            self._log(f"ERROR loading background: {e}")

    def _on_clear_backgrounds(self):
        self._backgrounds = []
        self._update_background_label()
        self._log("Cleared extra off-target backgrounds")

    def _filter_targets(self, query: str):
        """Filter target IDs by user search text."""
        query = query.strip().lower()
        current = self.target_combo.currentText()
        self.target_combo.clear()
        matches = [seq_id for seq_id in self._all_target_ids if query in seq_id.lower()]
        self.target_combo.addItems(matches)
        if current in matches:
            self.target_combo.setCurrentText(current)
        self.target_combo.setEnabled(bool(matches))

    def _on_target_browse(self):
        """Load one target sequence from a FASTA file."""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select Target FASTA",
            "",
            "FASTA Files (*.fa *.fasta *.fna);;All Files (*)",
        )
        if not path:
            return
        try:
            self._uploaded_target = load_first_fasta_record(path)
            self.target_path_edit.setText(path)
            self._log(f"Loaded target FASTA: {self._uploaded_target[0]}")
        except Exception as e:
            self._uploaded_target = None
            self._log(f"ERROR loading target FASTA: {e}")

    def _refresh_saved_transcriptomes(self):
        """Refresh saved transcriptome selector from the local manifest."""
        self._saved_entries = TranscriptomeIndex.list_saved()
        current_key = self.saved_combo.currentData() if hasattr(self, "saved_combo") else None
        self.saved_combo.clear()
        for entry in self._saved_entries:
            stats = entry.get("stats") or {}
            label = (
                f"{entry.get('name', entry.get('key', 'transcriptome'))} "
                f"({stats.get('num_sequences', '?')} seqs)"
            )
            self.saved_combo.addItem(label, entry.get("key"))
        if current_key:
            idx = self.saved_combo.findData(current_key)
            if idx >= 0:
                self.saved_combo.setCurrentIndex(idx)
        has_saved = self.saved_combo.count() > 0
        self.saved_combo.setEnabled(has_saved)
        self.load_saved_btn.setEnabled(has_saved)

    def _on_load_saved(self):
        """Load a transcriptome that was previously indexed and saved."""
        key = self.saved_combo.currentData()
        if not key:
            return
        try:
            self._log(f"Loading saved transcriptome: {self.saved_combo.currentText()}...")
            self.transcriptome = TranscriptomeIndex.load_saved(key)
            stats = self.transcriptome.get_stats()
            self.set_transcriptome(self.transcriptome)
            self._log(f"Loaded saved transcriptome with {stats['num_sequences']} sequences")
            self.transcriptome_loaded.emit(self.transcriptome, stats)
        except Exception as e:
            self._log(f"ERROR loading saved transcriptome: {e}")

    def _on_target_source_changed(self, index: int):
        """Show only controls relevant to the selected target source."""
        from_transcriptome = index == 0
        paste = index == 1
        upload = index == 2
        self.search_edit.setVisible(from_transcriptome)
        self.target_combo.setVisible(from_transcriptome)
        self.custom_seq_edit.setVisible(paste)
        self.target_path_edit.setVisible(upload)
        self.target_browse_btn.setVisible(upload)

    def _log(self, message: str):
        """添加日志"""
        self.log_edit.append(message)
