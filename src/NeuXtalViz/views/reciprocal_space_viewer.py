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
from pyvistaqt import QtInteractor

class ReciprocalSpaceViewerView(QWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.proj_box = QCheckBox('Parallel Projection', self)
        self.proj_box.clicked.connect(self.change_proj)

        self.reset_button = QPushButton('Reset View', self)
        self.reset_button.clicked.connect(self.reset_view)

        self.view_combo = QComboBox(self)
        self.view_combo.addItem('[uvw]')
        self.view_combo.addItem('[hkl]')

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(-100, 100, 5, notation=notation)

        self.axis1_line = QLineEdit()
        self.axis2_line = QLineEdit()
        self.axis3_line = QLineEdit()

        self.axis1_line.setValidator(validator)
        self.axis2_line.setValidator(validator)
        self.axis3_line.setValidator(validator)

        self.manual_button = QPushButton('View Plane', self)

        self.px_button = QPushButton('+Qx', self)
        self.py_button = QPushButton('+Qy', self)
        self.pz_button = QPushButton('+Qz', self)

        self.mx_button = QPushButton('-Qx', self)
        self.my_button = QPushButton('-Qy', self)
        self.mz_button = QPushButton('-Qz', self)

        self.a_star_button = QPushButton('a*', self)
        self.b_star_button = QPushButton('b*', self)
        self.c_star_button = QPushButton('c*', self)

        self.a_button = QPushButton('a', self)
        self.b_button = QPushButton('b', self)
        self.c_button = QPushButton('c', self)

        self.frame = QFrame()

        self.plotter = QtInteractor(self.frame)

        layout = QVBoxLayout()
        camera_layout = QGridLayout()
        plot_layout = QHBoxLayout()

        self.axis1_label = QLabel('h', self)
        self.axis2_label = QLabel('k', self)
        self.axis3_label = QLabel('l', self)

        camera_layout.addWidget(self.proj_box, 0, 0)
        camera_layout.addWidget(self.reset_button, 1, 0)

        camera_layout.addWidget(self.px_button, 0, 1)
        camera_layout.addWidget(self.py_button, 0, 2)
        camera_layout.addWidget(self.pz_button, 0, 3)
        camera_layout.addWidget(self.a_star_button, 0, 4)
        camera_layout.addWidget(self.b_star_button, 0, 5)
        camera_layout.addWidget(self.c_star_button, 0, 6)

        camera_layout.addWidget(self.mx_button, 1, 1)
        camera_layout.addWidget(self.my_button, 1, 2)
        camera_layout.addWidget(self.mz_button, 1, 3)
        camera_layout.addWidget(self.a_button, 1, 4)
        camera_layout.addWidget(self.b_button, 1, 5)
        camera_layout.addWidget(self.c_button, 1, 6)

        camera_layout.addWidget(self.axis1_label, 0, 7, Qt.AlignCenter)
        camera_layout.addWidget(self.axis2_label, 0, 8, Qt.AlignCenter)
        camera_layout.addWidget(self.axis3_label, 0, 9, Qt.AlignCenter)

        camera_layout.addWidget(self.axis1_line, 1, 7)
        camera_layout.addWidget(self.axis2_line, 1, 8)
        camera_layout.addWidget(self.axis3_line, 1, 9)

        camera_layout.addWidget(self.view_combo, 0, 10)
        camera_layout.addWidget(self.manual_button, 1, 10)

        plot_layout.addWidget(self.plotter.interactor)

        layout.addLayout(camera_layout)
        layout.addLayout(plot_layout)

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

    def change_proj(self):

        if self.proj_box.isChecked():
            self.plotter.enable_parallel_projection()
        else:
            self.plotter.disable_parallel_projection()

    def reset_view(self):

        self.plotter.reset_camera()
        self.plotter.view_isometric()

    def set_transform(self, T):

        if T is not None:

            b = pv._vtk.vtkMatrix4x4()
            for i in range(3):
                for j in range(3):
                    b.SetElement(i,j,T[i,j])

            actor = self.plotter.add_axes(xlabel='a*',
                                          ylabel='b*',
                                          zlabel='c*')
            actor.SetUserMatrix(b)

    def view_xy(self):

        self.plotter.view_xy()

    def view_yz(self):

        self.plotter.view_yz()

    def view_zx(self):

        self.plotter.view_zx()

    def view_yx(self):

        self.plotter.view_yx()

    def view_zy(self):

        self.plotter.view_zy()

    def view_xz(self):

        self.plotter.view_xz()

    def view_vector(self, vecs):

        if len(vecs) == 2:
            self.plotter.view_vector(*vecs)
        else:
            self.plotter.view_vector(vecs)

    def get_manual_indices(self):

        axes_type = self.view_combo.currentText()

        axes = [self.axis1_line, self.axis2_line, self.axis3_line]
        valid_axes = all([axis.hasAcceptableInput() for axis in axes])

        if valid_axes:

            axis1 = float(self.axis1_line.text())
            axis2 = float(self.axis2_line.text())
            axis3 = float(self.axis3_line.text())

            ind = np.array([axis1,axis2,axis3])

            return axes_type, ind