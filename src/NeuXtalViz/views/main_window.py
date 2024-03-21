import os
os.environ['QT_API'] = 'pyqt5'

from qtpy.QtWidgets import (QHBoxLayout,
                            QVBoxLayout,
                            QWidget,
                            QTabWidget,
                            QPushButton)

import pyvista
pyvista.set_plot_theme('document')

class MainWindow(QWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        # self.tabs = QTabWidget()

        # rsv_view = ReciprocalSpaceViewerView(self)
        # rsv_model = ReciprocalSpaceViewerModel()
        # self.rsv = ReciprocalSpaceViewer(rsv_view, rsv_model)
        # self.tabs.addTab(rsv_view, 'ReciprocalSpaceViewer')

        # spi_view = SatellitePeakIndexerView(self)
        # spi_model = SatellitePeakIndexerModel()
        # self.spi = SatellitePeakIndexer(spi_view, spi_model)
        # self.tabs.addTab(spi_view, 'SatellitePeakIndexer')

        # co_view = CoverageOptimizerView(self)
        # co_model = CoverageOptimizerModel()
        # self.co = CoverageOptimizer(co_view, co_model)
        # self.tabs.addTab(co_view, 'CoverageOptimizer')

        # rss_view = ReciprocalSpaceSlicerView(self)
        # rss_model = ReciprocalSpaceSlicerModel()
        # self.rss = ReciprocalSpaceSlicer(rss_view, rss_model)
        # self.tabs.addTab(rss_view, 'ReciprocalSpaceSlicer')

        # cs_view = CrystalStructureView(self)
        # cs_model = CrystalStructureModel()
        # self.cs = CrystalStructure(cs_view, cs_model)
        # self.tabs.addTab(cs_view, 'CrystalStructure')

        # s_view = SampleView(self)
        # s_model = SampleModel()
        # self.s = Sample(s_view, s_model)

        layout = QVBoxLayout()
        #layout.addWidget(self.s)

        self.setLayout(layout)