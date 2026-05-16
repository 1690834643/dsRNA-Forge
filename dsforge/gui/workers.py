"""
GUI 工作线程
QThreadPool + QRunnable 封装设计任务
"""

from PyQt6.QtCore import QRunnable, pyqtSignal, QObject

from dsforge.controller.design_task import DesignCancelled, DesignTask, DesignConfig
from dsforge.controller.design_task_parallel import ParallelDesignTask
from dsforge.core.sequence import TranscriptomeIndex


def create_design_task(db, config: DesignConfig):
    """Create the execution backend that matches the GUI core-count setting."""
    if config.n_cores and config.n_cores > 1:
        return ParallelDesignTask(
            db_manager=db,
            n_cores=config.n_cores,
            batch_size=config.batch_size,
        )
    return DesignTask(db_manager=db)


class WorkerSignals(QObject):
    """工作线程信号"""
    progress = pyqtSignal(str, float)  # step, percent
    result = pyqtSignal(dict)          # 最终结果
    error = pyqtSignal(str)            # 错误信息
    finished = pyqtSignal()            # 完成信号


class DesignTaskWorker(QRunnable):
    """设计任务工作线程"""

    def __init__(
        self,
        transcriptome: TranscriptomeIndex,
        target_seq_id: str,
        config: DesignConfig,
        db_path: str = None,
    ):
        super().__init__()
        self.transcriptome = transcriptome
        self.target_seq_id = target_seq_id
        self.config = config
        self.db_path = db_path
        self.signals = WorkerSignals()
        self._cancelled = False

    def run(self):
        """执行任务"""
        try:
            from dsforge.database.manager import DatabaseManager

            db = DatabaseManager(self.db_path)
            task = create_design_task(db, self.config)

            def progress_cb(step: str, percent: float):
                if self._cancelled:
                    raise DesignCancelled("Design task was cancelled by the user")
                self.signals.progress.emit(step, percent)

            if isinstance(task, ParallelDesignTask):
                result = task.run_parallel(
                    transcriptome=self.transcriptome,
                    target_seq_id=self.target_seq_id,
                    config=self.config,
                    progress_callback=progress_cb,
                )
            else:
                result = task.run(
                    transcriptome=self.transcriptome,
                    target_seq_id=self.target_seq_id,
                    config=self.config,
                    progress_callback=progress_cb,
                )

            if not self._cancelled:
                self.signals.result.emit(result)

        except DesignCancelled:
            pass
        except Exception as e:
            import traceback
            error_msg = f"{e}\n{traceback.format_exc()}"
            self.signals.error.emit(error_msg)
        finally:
            self.signals.finished.emit()

    def cancel(self):
        """Request cancellation; the running design stops at the next progress checkpoint."""
        self._cancelled = True
