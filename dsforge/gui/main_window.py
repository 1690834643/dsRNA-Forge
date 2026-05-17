"""
dsRNA-Forge 主窗口
PyQt6 GUI 入口
"""

from PyQt6.QtWidgets import (
    QMainWindow,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QTabWidget,
    QLabel,
    QSplitter,
    QMessageBox,
    QFileDialog,
)
from PyQt6.QtCore import Qt, QThreadPool
from dataclasses import fields

from dsforge.gui.transcript_panel import TranscriptPanel
from dsforge.gui.config_panel import ConfigPanel
from dsforge.gui.progress_panel import ProgressPanel
from dsforge.gui.results_panel import ResultsPanel
from dsforge.gui.history_panel import HistoryPanel
from dsforge.gui.workers import DesignTaskWorker
from dsforge.gui.cache_dialog import CacheManagerDialog
from dsforge.controller.design_task import DesignConfig
from dsforge.database.manager import DatabaseManager
from dsforge.core.diagnostics import diagnose_design_outcome
from dsforge.core.presets import apply_preset_to_config
from dsforge.core.project import load_project_file, save_project_file
from dsforge.core.sequence import (
    TranscriptomeIndex,
    clone_with_custom_sequence,
    make_safe_sequence_id,
    merge_background_transcriptomes,
)


def _design_config_from_payload(payload):
    if not isinstance(payload, dict):
        return None
    allowed = {field.name for field in fields(DesignConfig)}
    values = {key: value for key, value in payload.items() if key in allowed}
    try:
        return DesignConfig(**values)
    except Exception:
        return None


class MainWindow(QMainWindow):
    """主窗口"""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("dsRNA-Forge v0.1.2")
        self.setMinimumSize(1400, 900)

        self.thread_pool = QThreadPool()
        self.db = DatabaseManager()
        self.current_worker = None
        self.transcriptome = None
        self.current_results = []
        self.current_design_config = None
        self.current_target_selection = None

        self._setup_ui()
        self._connect_signals()
        self._load_history()

    def _setup_ui(self):
        """设置界面布局"""
        central = QWidget()
        self.setCentralWidget(central)

        main_layout = QHBoxLayout(central)
        main_layout.setContentsMargins(10, 10, 10, 10)
        main_layout.setSpacing(10)

        # 左侧：配置 + 转录组（上下分割）
        left_splitter = QSplitter(Qt.Orientation.Vertical)

        self.transcript_panel = TranscriptPanel()
        left_splitter.addWidget(self.transcript_panel)

        self.config_panel = ConfigPanel()
        left_splitter.addWidget(self.config_panel)

        left_splitter.setSizes([350, 450])

        # 右侧：结果 + 进度 + 历史（上下分割）
        right_splitter = QSplitter(Qt.Orientation.Vertical)

        # 结果和历史用标签页
        self.tabs = QTabWidget()
        self.results_panel = ResultsPanel()
        self.history_panel = HistoryPanel()
        self.tabs.addTab(self.results_panel, "Results")
        self.tabs.addTab(self.history_panel, "History")
        right_splitter.addWidget(self.tabs)

        self.progress_panel = ProgressPanel()
        right_splitter.addWidget(self.progress_panel)

        right_splitter.setSizes([600, 200])

        # 主分割器：左右
        main_splitter = QSplitter(Qt.Orientation.Horizontal)
        main_splitter.addWidget(left_splitter)
        main_splitter.addWidget(right_splitter)
        main_splitter.setSizes([450, 950])

        main_layout.addWidget(main_splitter)

        self.statusBar().showMessage("Ready — Load a transcriptome to begin")

    def _connect_signals(self):
        """连接信号"""
        # 转录组加载完成
        self.transcript_panel.transcriptome_loaded.connect(self._on_transcriptome_loaded)
        self.transcript_panel.transcriptome_cleared.connect(self._on_transcriptome_cleared)
        self.transcript_panel.manage_cache_requested.connect(self._on_manage_cache)

        # 开始设计
        self.config_panel.start_btn.clicked.connect(self._on_start_design)
        self.config_panel.cancel_btn.clicked.connect(self._on_cancel_design)

        # 导出按钮
        self.results_panel.export_csv_btn.clicked.connect(self._on_export_csv)
        self.results_panel.export_excel_btn.clicked.connect(self._on_export_excel)
        self.results_panel.export_fasta_btn.clicked.connect(self._on_export_fasta)
        self.results_panel.export_report_btn.clicked.connect(self._on_export_report)
        self.results_panel.export_primers_btn.clicked.connect(self._on_export_primers)
        self.results_panel.save_project_btn.clicked.connect(self._on_save_project)
        self.results_panel.open_project_btn.clicked.connect(self._on_open_project)

        # 历史任务加载
        self.history_panel.task_selected.connect(self._on_load_history_task)
        self.history_panel.refresh_btn.clicked.connect(self._load_history)

    def _load_history(self):
        """加载历史任务"""
        tasks = self.db.list_tasks(limit=50)
        self.history_panel.load_history(tasks)

    def _on_transcriptome_loaded(self, transcriptome, stats):
        """转录组加载完成"""
        self.transcriptome = transcriptome
        self.statusBar().showMessage(
            f"Transcriptome loaded: {stats['num_sequences']} sequences, "
            f"{stats['total_nt']} nt, GC={stats['gc_content']:.1f}%"
        )
        self.config_panel.start_btn.setEnabled(True)

    def _on_transcriptome_cleared(self):
        """Keep main-window state in sync when the transcriptome panel is cleared."""
        self.transcriptome = None
        self.current_target_selection = None
        self.config_panel.start_btn.setEnabled(False)
        self.statusBar().showMessage("Ready — Load a transcriptome to begin")

    def _on_start_design(self):
        """开始设计"""
        if self.transcriptome is None:
            QMessageBox.warning(self, "No Transcriptome", "Please load a transcriptome first.")
            return

        # 获取配置
        gui_config = self.config_panel.get_config()
        config = DesignConfig(
            mode=gui_config["mode"],
            length_min=gui_config["length"]["min"],
            length_max=gui_config["length"]["max"],
            gc_min=gui_config["gc"]["min"],
            gc_max=gui_config["gc"]["max"],
            enabled_rules=[name for name, enabled in gui_config["rules"].items() if enabled],
            n_cores=gui_config["cores"],
            preset=gui_config.get("preset", "balanced"),
        )
        config = apply_preset_to_config(config, gui_config.get("preset", "balanced"))

        # 获取目标序列
        target_selection = self.transcript_panel.get_target_selection()
        if target_selection is None:
            QMessageBox.warning(self, "No Target", "Please select, paste, or upload a target sequence.")
            return

        design_transcriptome = self.transcriptome
        target_id = target_selection["id"]
        if target_selection["source"] == "custom":
            target_seq = target_selection.get("sequence", "")
            if not target_seq:
                QMessageBox.warning(self, "No Target Sequence", "Please paste or upload a valid target sequence.")
                return
            target_id = make_safe_sequence_id(target_id, "custom_target")
            if target_id in self.transcriptome.sequences:
                target_id = make_safe_sequence_id(f"{target_id}_custom", "custom_target")
            design_transcriptome = clone_with_custom_sequence(
                self.transcriptome,
                target_id,
                target_seq,
            )

        backgrounds = self.transcript_panel.get_background_indexes()
        if backgrounds:
            design_transcriptome = merge_background_transcriptomes(design_transcriptome, backgrounds)

        self.current_design_config = config
        self.current_target_selection = target_selection

        # 重置进度和结果
        self.progress_panel.reset()
        self.results_panel.clear_results()
        self.config_panel.start_btn.setEnabled(False)
        self.config_panel.cancel_btn.setEnabled(True)

        # 创建工作线程
        self.current_worker = DesignTaskWorker(
            transcriptome=design_transcriptome,
            target_seq_id=target_id,
            config=config,
        )
        self.current_worker.signals.progress.connect(self._on_progress)
        self.current_worker.signals.result.connect(self._on_result)
        self.current_worker.signals.error.connect(self._on_error)
        self.current_worker.signals.finished.connect(self._on_finished)

        self.thread_pool.start(self.current_worker)
        bg_text = f" with {len(backgrounds)} extra background(s)" if backgrounds else ""
        self.statusBar().showMessage(f"Running {config.mode} design on '{target_id}'{bg_text}...")

    def _on_cancel_design(self):
        """取消设计"""
        if self.current_worker:
            self.current_worker.cancel()
            self.progress_panel.set_status("Cancelling...")

    def _on_progress(self, step: str, percent: float):
        """接收进度更新"""
        self.progress_panel.set_status(step)
        self.progress_panel.set_progress(percent)
        self.progress_panel.log(step)

    def _on_result(self, result: dict):
        """接收结果"""
        self.current_results = result.get("results", [])
        self.progress_panel.set_status("Loading results into table...")
        self.progress_panel.set_progress(99)
        self.results_panel.load_results(self.current_results)
        self.progress_panel.set_progress(100)
        self.progress_panel.set_status("Design complete")

        summary = result.get("summary", {})
        if "raw_candidates" in summary:
            stats_text = (
                f"Raw: {summary.get('raw_candidates', 0)}, "
                f"Recommended: {summary.get('nonredundant_candidates', summary.get('total_candidates', 0))}, "
                f"Passed: {summary.get('passed_candidates', 0)}"
            )
        else:
            stats_text = (
                f"Total: {summary.get('total_candidates', 0)}, "
                f"Passed: {summary.get('passed_candidates', 0)}"
            )
        self.progress_panel.set_stats(
            stats_text
        )

        self.statusBar().showMessage(
            f"Design complete: {summary.get('nonredundant_candidates', summary.get('total_candidates', 0))} recommendations, "
            f"{summary.get('passed_candidates', 0)} passed"
        )

        if summary.get("total_candidates", 0) == 0 or summary.get("passed_candidates", 0) == 0:
            target_seq = result.get("target_seq", "")
            if not target_seq and hasattr(self, "current_worker") and self.current_worker:
                target_seq = self.current_worker.transcriptome.get_sequence(self.current_worker.target_seq_id) or ""
            diagnosis = diagnose_design_outcome(
                target_seq=target_seq,
                config=result.get("config"),
                total_candidates=summary.get("total_candidates", 0),
                passed_candidates=summary.get("passed_candidates", 0),
            )
            message = "可能原因：\n- " + "\n- ".join(diagnosis["reasons"])
            message += "\n\n建议：\n- " + "\n- ".join(diagnosis["suggestions"])
            if getattr(result.get("config"), "preset", "balanced") != "relaxed":
                message += "\n- 可切换到 Relaxed - rescue difficult targets 后重跑，以放宽 GC/脱靶阈值。"
            self.progress_panel.log(message)
            QMessageBox.information(self, "No passing design found", message)

        # 刷新历史
        self._load_history()

    def _on_error(self, error_msg: str):
        """接收错误"""
        self.progress_panel.log(f"ERROR: {error_msg}")
        QMessageBox.critical(self, "Design Error", error_msg[:500])

    def _on_finished(self):
        """任务完成"""
        self.config_panel.start_btn.setEnabled(True)
        self.config_panel.cancel_btn.setEnabled(False)
        self.current_worker = None

    def _on_load_history_task(self, task_id: int):
        """从历史加载任务结果"""
        results = self.db.get_results(task_id)
        if results:
            self.current_results = results
            self.results_panel.load_results(results)
            self.statusBar().showMessage(f"Loaded {len(results)} results from task {task_id}")
        else:
            QMessageBox.information(self, "No Results", f"Task {task_id} has no results.")

    def _on_export_csv(self):
        """导出 CSV"""
        if not hasattr(self, 'current_results') or not self.current_results:
            QMessageBox.information(self, "No Results", "No results to export.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export CSV", "results.csv", "CSV Files (*.csv)")
        if path:
            from dsforge.controller.exporter import ResultExporter
            self._run_export(lambda: ResultExporter().export_csv(self.current_results, path), f"Exported to {path}")

    def _run_export(self, export_func, success_message: str):
        try:
            export_func()
            self.statusBar().showMessage(success_message)
        except ImportError as e:
            QMessageBox.warning(self, "Missing Dependency", str(e))
        except Exception as e:
            QMessageBox.critical(self, "Export Failed", str(e))

    def _on_export_report(self):
        """导出实验验证报告"""
        if not self.current_results:
            QMessageBox.information(self, "No Results", "No results to export.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export Validation Report", "validation_report.xlsx", "Excel Files (*.xlsx)")
        if path:
            from dsforge.controller.exporter import ResultExporter
            self._run_export(
                lambda: ResultExporter().export_validation_report(self.current_results, path),
                f"Exported validation report to {path}",
            )

    def _on_export_primers(self):
        """导出 T7 引物订购表"""
        if not self.current_results:
            QMessageBox.information(self, "No Results", "No results to export.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export Primers", "t7_primers.csv", "CSV Files (*.csv)")
        if path:
            from dsforge.controller.exporter import ResultExporter
            self._run_export(
                lambda: ResultExporter().export_primer_order_csv(self.current_results, path),
                f"Exported primers to {path}",
            )

    def _project_payload(self):
        transcriptome_meta = {}
        if self.transcriptome is not None:
            transcriptome_meta = {
                "source_file": self.transcriptome.source_file,
                "cache_key": self.transcriptome.source_hash,
                "cache_path": self.transcriptome.cache_path,
                "stats": self.transcriptome.get_stats(),
            }
        return {
            "transcriptome": transcriptome_meta,
            "target": self.current_target_selection or {},
            "config": self.current_design_config or {},
            "results": self.current_results,
        }

    def _on_save_project(self):
        """保存 dsRNA-Forge 项目文件"""
        try:
            path, _ = QFileDialog.getSaveFileName(
                self,
                "Save dsRNA-Forge Project",
                "experiment.dsforge_project",
                "dsRNA-Forge Project (*.dsforge_project);;JSON Files (*.json)",
            )
            if path:
                save_project_file(path, self._project_payload())
                self.statusBar().showMessage(f"Saved project to {path}")
        except Exception as e:
            QMessageBox.critical(self, "Save Project Failed", str(e))

    def _on_open_project(self):
        """打开 dsRNA-Forge 项目文件"""
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Open dsRNA-Forge Project",
            "",
            "dsRNA-Forge Project (*.dsforge_project);;JSON Files (*.json)",
        )
        if not path:
            return
        try:
            project = load_project_file(path)
            cache_key = (project.get("transcriptome") or {}).get("cache_key")
            if cache_key:
                try:
                    self.transcriptome = TranscriptomeIndex.load_saved(cache_key)
                    self.transcript_panel.set_transcriptome(self.transcriptome)
                    self.config_panel.start_btn.setEnabled(True)
                except Exception:
                    pass
            config_payload = project.get("config") or {}
            self.config_panel.apply_config(config_payload)
            self.current_design_config = _design_config_from_payload(config_payload)
            self.current_results = project.get("results") or []
            self.current_target_selection = project.get("target") or {}
            self.transcript_panel.apply_target_selection(self.current_target_selection)
            self.results_panel.load_results(self.current_results)
            self.tabs.setCurrentWidget(self.results_panel)
            self.statusBar().showMessage(f"Opened project {path}")
        except Exception as e:
            QMessageBox.critical(self, "Open Project Failed", str(e))

    def _on_manage_cache(self):
        """打开缓存管理窗口"""
        dialog = CacheManagerDialog(self)
        dialog.exec()
        self.transcript_panel._refresh_saved_transcriptomes()

    def _on_export_excel(self):
        """导出 Excel"""
        if not hasattr(self, 'current_results') or not self.current_results:
            QMessageBox.information(self, "No Results", "No results to export.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export Excel", "results.xlsx", "Excel Files (*.xlsx)")
        if path:
            from dsforge.controller.exporter import ResultExporter
            self._run_export(lambda: ResultExporter().export_excel(self.current_results, path), f"Exported to {path}")

    def _on_export_fasta(self):
        """导出 FASTA"""
        if not hasattr(self, 'current_results') or not self.current_results:
            QMessageBox.information(self, "No Results", "No results to export.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export FASTA", "results.fa", "FASTA Files (*.fa *.fasta)")
        if path:
            from dsforge.controller.exporter import ResultExporter
            self._run_export(lambda: ResultExporter().export_fasta(self.current_results, path), f"Exported to {path}")
