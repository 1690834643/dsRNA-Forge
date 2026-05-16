#!/usr/bin/env python3
"""
无头 GUI 组件测试
验证 GUI 组件能正确初始化（无需显示服务器）
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# 使用 offscreen 平台
os.environ["QT_QPA_PLATFORM"] = "offscreen"

from PyQt6.QtWidgets import QApplication, QFileDialog, QMessageBox
from PyQt6.QtCore import Qt

app = QApplication(sys.argv)

print("=" * 60)
print("GUI Component Test (Headless)")
print("=" * 60)

# Test 1: MainWindow
print("\nTest 1: MainWindow...")
from dsforge.gui.main_window import MainWindow
from dsforge.core.sequence import TranscriptomeIndex
window = MainWindow()
assert window is not None
loaded_for_clear = TranscriptomeIndex()
loaded_for_clear.sequences = {"gene_clear": "AUGCAUGCAUGCAUGCAUGCAUGC"}
loaded_for_clear._compute_stats()
window._on_transcriptome_loaded(loaded_for_clear, loaded_for_clear.get_stats())
assert window.transcriptome is loaded_for_clear
assert window.config_panel.start_btn.isEnabled()
window.transcript_panel._on_clear()
assert window.transcriptome is None
assert not window.config_panel.start_btn.isEnabled()
print("  ✓ MainWindow created")

# Test 2: TranscriptPanel
print("\nTest 2: TranscriptPanel...")
from dsforge.gui.transcript_panel import TranscriptPanel
panel = TranscriptPanel()
assert panel is not None
assert panel.saved_combo is not None
assert panel.load_saved_btn is not None
assert panel.manage_cache_btn is not None
background = TranscriptomeIndex()
background.sequences = {"bg_gene": "UUUUUUUUUUUUUUUUUUUU"}
background._compute_stats()
panel.add_background_index("host", background)
assert panel.get_background_indexes()[0][0] == "host"
assert TranscriptPanel.background_label_from_path(r"C:\Users\me\host transcriptome.fa") == "host_transcriptome"
loaded = TranscriptomeIndex()
loaded.sequences = {"gene_alpha": "AUGCAUGCAUGC", "gene_beta": "GGGGCCCCAAAA"}
loaded.source_file = "saved.fa"
loaded.source_hash = "abc"
loaded._compute_stats()
panel.set_transcriptome(loaded)
assert panel.target_combo.count() == 2
panel.search_edit.setText("alpha")
panel.apply_target_selection({"source": "transcriptome", "id": "gene_beta"})
assert panel.target_combo.currentText() == "gene_beta"
panel._all_target_ids = ["gene_alpha", "gene_beta", "other"]
panel._filter_targets("")
panel.target_combo.setEnabled(True)
panel.search_edit.setText("beta")
assert panel.target_combo.count() == 1
assert panel.target_combo.currentText() == "gene_beta"
panel.target_source_combo.setCurrentIndex(1)
panel.custom_seq_edit.setPlainText("atgc atgc")
selection = panel.get_target_selection()
assert selection["id"] == "custom_target"
assert selection["sequence"] == "AUGCAUGC"
print("  ✓ TranscriptPanel created")

# Test 3: ConfigPanel
print("\nTest 3: ConfigPanel...")
from dsforge.gui.config_panel import ConfigPanel
panel = ConfigPanel()
config = panel.get_config()
assert config["mode"] == "siRNA"  # 默认第一个选项
assert config["length"]["min"] == 21
assert config["length"]["max"] == 21
assert "sgRNA" in [panel.mode_combo.itemText(i) for i in range(panel.mode_combo.count())][-1]
panel.mode_combo.setCurrentIndex(3)
assert panel.get_config()["mode"] == "sgRNA"
panel.apply_config({
    "mode": "DsiRNA",
    "preset": "strict",
    "length_min": 27,
    "length_max": 27,
    "gc_min": 36,
    "gc_max": 48,
    "enabled_rules": ["consensus", "jagla"],
    "n_cores": 1,
})
restored = panel.get_config()
assert restored["mode"] == "DsiRNA"
assert restored["preset"] == "strict"
assert restored["gc"]["min"] == 36
assert restored["rules"]["jagla"] is True
assert restored["rules"]["reynolds"] is False
assert config["preset"] == "balanced"
panel.preset_combo.setCurrentText("Relaxed - rescue difficult targets")
assert panel.get_config()["preset"] == "relaxed"
assert panel.advanced_group.isHidden()
panel.advanced_toggle.setChecked(True)
assert not panel.advanced_group.isHidden()
print(f"  ✓ ConfigPanel created, default mode: {config['mode']}")

# Test 4: ProgressPanel
print("\nTest 4: ProgressPanel...")
from dsforge.gui.progress_panel import ProgressPanel
panel = ProgressPanel()
panel.set_progress(50.0)
panel.set_status("Testing...")
assert panel.progress_bar.value() == 50
print("  ✓ ProgressPanel works")

# Test 5: ResultsPanel
print("\nTest 5: ResultsPanel...")
from dsforge.gui.results_panel import ResultsPanel
panel = ResultsPanel()
panel.add_result(
    1,
    "AUGC",
    "0-4",
    85.5,
    True,
    cluster_size=3,
    risk_level="medium",
    risk_score=30,
    top_targets="gene_x",
    validation_direction="优先验证 gene_x",
)
assert panel.table.rowCount() == 1
assert panel.table.columnCount() >= 10
assert panel.export_report_btn is not None
assert panel.export_primers_btn is not None
panel.load_results([{
    "rank": 1,
    "sequence": "AUGC",
    "position": "0-4",
    "consensus_score": 85.5,
    "passed": True,
    "explanation": {"summary": "推荐理由：测试", "risk_notes": ["低风险"]},
    "validation_hits": [{"target_id": "gene_x", "validation_action": "qPCR"}],
    "region_map": "candidate [####]",
    "primers": {"t7_forward_primer": {"sequence": "TAATACGACTCACTATAGGGATGC"}},
}, {
    "rank": 2,
    "sequence": "GGGGGGGGGGGGGGGGGGGG",
    "position": "25-45",
    "consensus_score": 75.0,
    "passed": True,
    "sgrna": {
        "pam": "AGG",
        "genomic_pam": "AGG",
        "strand": "+",
        "cut_site": 42,
    },
    "explanation": {"summary": "推荐理由：sgRNA second"},
    "validation_hits": [],
    "primers": {
        "sgrna_cloning_oligos": {
            "px330_forward_oligo": {"sequence": "CACCGGGGGGGGGGGGGGGGGGGG"},
            "px330_reverse_oligo": {"sequence": "AAACCCCCCCCCCCCCCCCCCCCCC"},
        }
    },
}])
panel.table.sortItems(0, Qt.SortOrder.DescendingOrder)
panel.table.selectRow(0)
panel._show_selected_detail()
detail_text = panel.detail_edit.toPlainText()
assert "sgRNA second" in detail_text
assert "PAM AGG" in detail_text
assert "cut site 42" in detail_text
assert "CACCG" in detail_text
print("  ✓ ResultsPanel works")

# Test 5b: export errors are shown instead of escaping the slot
print("\nTest 5b: Export error handling...")
from dsforge.gui.main_window import MainWindow
from dsforge.controller import exporter as exporter_module
window = MainWindow()
window.current_results = [{"rank": 1, "sequence": "AUGC"}]
messages = []
original_get_save = QFileDialog.getSaveFileName
original_export_csv = exporter_module.ResultExporter.export_csv
original_critical = QMessageBox.critical

def fake_get_save(*args, **kwargs):
    return ("blocked.csv", "CSV Files (*.csv)")

def fake_export_csv(self, results, path):
    raise PermissionError("blocked path")

def fake_critical(parent, title, message):
    messages.append((title, message))

QFileDialog.getSaveFileName = fake_get_save
exporter_module.ResultExporter.export_csv = fake_export_csv
QMessageBox.critical = fake_critical
try:
    window._on_export_csv()
finally:
    QFileDialog.getSaveFileName = original_get_save
    exporter_module.ResultExporter.export_csv = original_export_csv
    QMessageBox.critical = original_critical

assert messages and "Export Failed" in messages[0][0]
assert "blocked path" in messages[0][1]
print("  ✓ Export error handling works")

# Test 6: HistoryPanel
print("\nTest 6: HistoryPanel...")
from dsforge.gui.history_panel import HistoryPanel
panel = HistoryPanel()
panel.load_history([
    {"id": 1, "mode": "siRNA", "target_seq_id": "gene1", "status": "completed", "created_at": "2026-05-15"},
])
assert panel.table.rowCount() == 1
print("  ✓ HistoryPanel works")

# Test 7: Workers
print("\nTest 7: DesignTaskWorker...")
from dsforge.gui.workers import DesignTaskWorker, WorkerSignals
from dsforge.core.sequence import TranscriptomeIndex
from dsforge.controller.design_task import DesignConfig

transcriptome = TranscriptomeIndex()
transcriptome.sequences = {"test": "AUGCAUGCAUGC"}
transcriptome._compute_stats()

config = DesignConfig(mode="siRNA")
worker = DesignTaskWorker(transcriptome, "test", config)
assert worker.signals is not None
print("  ✓ DesignTaskWorker created")

print("\n" + "=" * 60)
print("ALL GUI TESTS PASSED ✓")
print("=" * 60)
