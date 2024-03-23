from qtpy.QtWidgets import (QWidget,
                            QFrame,
                            QGridLayout,
                            QHBoxLayout,
                            QVBoxLayout,
                            QPushButton,
                            QCheckBox,
                            QComboBox,
                            QLineEdit,
                            QLabel)

from qtpy.QtGui import QDoubleValidator
from PyQt5.QtCore import Qt

import numpy as np
import pyvista as pv

from NeuXtalViz.views.base_view import NeuXtalVizWidget

class ReciprocalSpaceViewerView(NeuXtalVizWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        layout = QVBoxLayout()
        widget = NeuXtalVizWidget()
        layout.addWidget(widget)

        self.setLayout(layout)

    def add_peaks(self, peak_dict):

        self.plotter.clear_actors()

        transforms = peak_dict['transforms']
        intensities = peak_dict['intensities']
        numbers = peak_dict['numbers']

        sphere = pv.Icosphere(radius=1, nsub=1)

        geoms, self.indexing = [], {}
        for i, (T, I, no) in enumerate(zip(transforms, intensities, numbers)):
            ellipsoid = sphere.copy().transform(T)
            ellipsoid['scalars'] = np.full(sphere.n_cells, I)
            geoms.append(ellipsoid)
            self.indexing[i] = no

        multiblock = pv.MultiBlock(geoms)

        _, mapper = self.plotter.add_composite(multiblock,
                                               scalars='scalars',
                                               log_scale=True,
                                               smooth_shading=True,
                                               scalar_bar_args={'title':
                                                                'Intensity'})

        self.mapper = mapper

        self.plotter.enable_block_picking(callback=self.highlight,
                                          side='left')
        self.plotter.enable_block_picking(callback=self.highlight,
                                          side='right')

        self.plotter.add_camera_orientation_widget()
        self.plotter.enable_depth_peeling()

        self.change_proj()

    def highlight(self, index, dataset):

        color = self.mapper.block_attr[index].color

        if color == 'pink':
            color, select = None, False
        else:
            color, select = 'pink', True

        self.mapper.block_attr[index].color = color

        print('peak_no = {}'.format(self.indexing[index]))

        return self.indexing[index], select