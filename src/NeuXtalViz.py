import os
import sys

os.environ['QT_API'] = 'pyqt5'

from qtpy.QtWidgets import (QApplication,
                            QMainWindow,
                            QWidget, 
                            QAction,
                            QStackedWidget,
                            QVBoxLayout)

from mantid.kernel import Logger
from mantidqt.gui_helper import set_matplotlib_backend

set_matplotlib_backend()

from NeuXtalViz._version import __version__  

logger = Logger('NeuXtalViz')

import pyvista
pyvista.set_plot_theme('document')

# from NeuXtalViz.views.reciprocal_space_viewer import ReciprocalSpaceViewerView
# from NeuXtalViz.models.reciprocal_space_viewer import ReciprocalSpaceViewerModel
# from NeuXtalViz.presenters.reciprocal_space_viewer import ReciprocalSpaceViewer

from NeuXtalViz.views.crystal_structure_tools import CrystalStructureView
from NeuXtalViz.models.crystal_structure_tools import CrystalStructureModel
from NeuXtalViz.presenters.crystal_structure_tools import CrystalStructure

class NeuXtalViz(QMainWindow):

    __instance = None

    def __new__(cls):
        if NeuXtalViz.__instance is None:
            NeuXtalViz.__instance = QMainWindow.__new__(cls)  
        return NeuXtalViz.__instance

    def __init__(self, parent=None):
        super().__init__(parent)
        logger.information('NeuXtalViz {}'.format(__version__))

        self.setWindowTitle('NeuXtalViz {}'.format(__version__))
        self.resize(1200, 900)
        app_menu = self.menuBar().addMenu('Applications')

        main_window = QWidget(self)
        self.setCentralWidget(main_window)

        layout = QVBoxLayout(main_window)

        self.stack = QStackedWidget()

        cs_action = QAction('Crystal Structure', self)
        cs_action.triggered.connect(lambda: self.stack.setCurrentIndex(0))

        cs_view = CrystalStructureView(self)
        cs_model = CrystalStructureModel()
        self.cs = CrystalStructure(cs_view, cs_model)
        self.stack.addWidget(cs_view)

        layout.addWidget(self.stack)

        app_menu.addAction(cs_action)

def gui():
    app = QApplication(sys.argv)
    window = NeuXtalViz()
    window.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    gui()