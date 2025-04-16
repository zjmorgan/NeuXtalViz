import sys
import traceback

from qtpy.QtCore import QRunnable, QThreadPool, Signal, QObject, Slot


class WorkerSignals(QObject):
    finished = Signal()
    error = Signal(tuple)
    progress = Signal(str, int)
    result = Signal(object)


class Worker(QRunnable):
    def __init__(self, task, *args, **kwargs):
        super().__init__()

        self.signals = WorkerSignals()
        self.task = task
        self.args = args
        self.kwargs = kwargs

        self.kwargs["progress"] = self.emit_progress

    @Slot()
    def run(self):
        try:
            result = self.task(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit(result)
        finally:
            self.signals.finished.emit()

    def emit_progress(self, message, progress):
        self.signals.progress.emit(message, progress)

    def connect_result(self, process):
        self.signals.result.connect(process)

    def connect_finished(self, process):
        self.signals.finished.connect(process)

    def connect_progress(self, process):
        self.signals.progress.connect(process)


class ThreadPool(QThreadPool):
    def __init__(self):
        super().__init__()

    def start_worker_pool(self, worker):
        self.start(worker)
