import os
os.environ['QT_API'] = 'pyqt5'

from qtpy.QtWidgets import (QWidget, QAction, QStackedWidget, QVBoxLayout)

import pyvista
pyvista.set_plot_theme('document')

# from NeuXtalViz.views.reciprocal_space_viewer import ReciprocalSpaceViewerView
# from NeuXtalViz.models.reciprocal_space_viewer import ReciprocalSpaceViewerModel
# from NeuXtalViz.presenters.reciprocal_space_viewer import ReciprocalSpaceViewer

from NeuXtalViz.views.crystal_structure_tools import CrystalStructureView
from NeuXtalViz.models.crystal_structure_tools import CrystalStructureModel
from NeuXtalViz.presenters.crystal_structure_tools import CrystalStructure

class MainWindow(QWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        # app_menu = self.menuBar().addMenu('Applications')

        # self.stack = QStackedWidget()

        # layout = QVBoxLayout()

        # # rsv_view = ReciprocalSpaceViewerView(self)
        # # rsv_model = ReciprocalSpaceViewerModel()
        # # self.rsv = ReciprocalSpaceViewer(rsv_view, rsv_model)
        # # layout.addWidget(rsv_view)

        # # spi_view = SatellitePeakIndexerView(self)
        # # spi_model = SatellitePeakIndexerModel()
        # # self.spi = SatellitePeakIndexer(spi_view, spi_model)
        # # self.tabs.addTab(spi_view, 'SatellitePeakIndexer')

        # # rss_view = ReciprocalSpaceSlicerView(self)
        # # rss_model = ReciprocalSpaceSlicerModel()
        # # self.rss = ReciprocalSpaceSlicer(rss_view, rss_model)
        # # self.tabs.addTab(rss_view, 'ReciprocalSpaceSlicer')

        # cs_action = QAction('Crystal Structure', self)
        # cs_action.triggered.connect(lambda: self.stack.setCurrentIndex(0))

        # cs_view = CrystalStructureView(self)
        # cs_model = CrystalStructureModel()
        # self.cs = CrystalStructure(cs_view, cs_model)
        # self.stack.addWidget.addWidget(cs_view)

        # app_menu.addAction(cs_action)

        # # s_view = SampleView(self)
        # # s_model = SampleModel()
        # # self.s = Sample(s_view, s_model)

        # self.addWidget(self.stack)