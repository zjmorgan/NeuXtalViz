import sys
import traceback

from PyQt5.QtCore import QRunnable, QThreadPool, pyqtSignal, QObject, pyqtSlot

class WorkerSignals(QObject):

    finished = pyqtSignal()
    error = pyqtSignal(tuple)
    progress = pyqtSignal(str, int)
    result = pyqtSignal(object)

class Worker(QRunnable):
    def __init__(self, task, *args, **kwargs):
        super().__init__()
        self.signals = WorkerSignals()
        self.task = task
        self.args = args
        self.kwargs = kwargs
        
        self.kwargs['progress'] = self.signals.progress

    @pyqtSlot()
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

class ThreadPool(QThreadPool):

    def __init__(self):

        super().__init__()
