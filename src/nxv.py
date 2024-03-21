import sys
from qtpy.QtWidgets import QApplication, QMainWindow

from mantid.kernel import Logger
from mantidqt.gui_helper import set_matplotlib_backend

set_matplotlib_backend()

from NeuXtalViz._version import __version__  
from NeuXtalViz.views.main_window import MainWindow

logger = Logger('NeuXtalViz')

class NeuXtalViz(QMainWindow):

    __instance = None

    def __new__(cls):
        if NeuXtalViz.__instance is None:
            NeuXtalViz.__instance = QMainWindow.__new__(cls)  
        return NeuXtalViz.__instance

    def __init__(self, parent=None):
        super().__init__(parent)
        logger.information(f'NeuXtalViz {__version__}')

        self.setWindowTitle(f'NeuXtalViz {__version__}')
        self.main_window = MainWindow(self)
        self.setCentralWidget(self.main_window)

def gui():
    app = QApplication(sys.argv)
    window = NeuXtalViz()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    gui()