from qtpy.QtWidgets import (QWidget,
                            QFrame,
                            QGridLayout,
                            QHBoxLayout,
                            QVBoxLayout,
                            QPushButton,
                            QLabel,
                            QCheckBox,
                            QComboBox,
                            QLineEdit)

from qtpy.QtGui import QDoubleValidator
from PyQt5.QtCore import Qt

import numpy as np
import pyvista as pv

from pyvistaqt import QtInteractor

class ReciprocalSpaceSlicerView(QWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.proj_box = QCheckBox('Parallel Projection', self)
        self.proj_box.clicked.connect(self.change_proj)

        self.reset_button = QPushButton('Reset View', self)
        self.reset_button.clicked.connect(self.reset_view)

        self.view_combo = QComboBox(self)
        self.view_combo.addItem('[hkl]')
        self.view_combo.addItem('[uvw]')

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(-100, 100, 5, notation=notation)

        self.axis1_line = QLineEdit()
        self.axis2_line = QLineEdit()
        self.axis3_line = QLineEdit()

        self.axis1_line.setValidator(validator)
        self.axis2_line.setValidator(validator)
        self.axis3_line.setValidator(validator)

        self.axis1_label = QLabel('h', self)
        self.axis2_label = QLabel('k', self)
        self.axis3_label = QLabel('l', self)

        self.manual_button = QPushButton('View Axis', self)

        self.a_star_button = QPushButton('a*', self)
        self.b_star_button = QPushButton('b*', self)
        self.c_star_button = QPushButton('c*', self)

        self.a_button = QPushButton('a', self)
        self.b_button = QPushButton('b', self)
        self.c_button = QPushButton('c', self)

        self.frame = QFrame()

        self.plotter = QtInteractor(self.frame)

        self.clim_combo = QComboBox(self)
        self.clim_combo.addItem('3-sig')
        self.clim_combo.addItem('1.5-IQR')

        self.replot_button = QPushButton('Redraw', self)

        self.slice_button = QPushButton('Slice Plane', self)

        self.slice_combo = QComboBox(self)
        self.slice_combo.addItem('[hkl]')
        self.slice_combo.addItem('[uvw]')

        layout = QVBoxLayout()
        camera_layout = QGridLayout()
        plot_layout = QHBoxLayout()
        slice_layout = QGridLayout()

        camera_layout.addWidget(self.proj_box, 0, 0)
        camera_layout.addWidget(self.reset_button, 1, 0)
        camera_layout.addWidget(self.a_star_button, 0, 1)
        camera_layout.addWidget(self.b_star_button, 0, 2)
        camera_layout.addWidget(self.c_star_button, 0, 3)
        camera_layout.addWidget(self.a_button, 1, 1)
        camera_layout.addWidget(self.b_button, 1, 2)
        camera_layout.addWidget(self.c_button, 1, 3)
        camera_layout.addWidget(self.axis1_label, 0, 4, Qt.AlignCenter)
        camera_layout.addWidget(self.axis2_label, 0, 5, Qt.AlignCenter)
        camera_layout.addWidget(self.axis3_label, 0, 6, Qt.AlignCenter)
        camera_layout.addWidget(self.axis1_line, 1, 4)
        camera_layout.addWidget(self.axis2_line, 1, 5)
        camera_layout.addWidget(self.axis3_line, 1, 6)
        camera_layout.addWidget(self.view_combo, 0, 7)
        camera_layout.addWidget(self.manual_button, 1, 7)

        plot_layout.addWidget(self.plotter.interactor)

        self.normal1_line = QLineEdit()
        self.normal2_line = QLineEdit()
        self.normal3_line = QLineEdit()
        self.scalar_line = QLineEdit()

        self.normal1_line.setValidator(validator)
        self.normal2_line.setValidator(validator)
        self.normal3_line.setValidator(validator)
        self.scalar_line.setValidator(validator)

        self.normal1_label = QLabel('h', self)
        self.normal2_label = QLabel('k', self)
        self.normal3_label = QLabel('l', self)
        scalar_label = QLabel('value', self)

        slice_layout.addWidget(self.clim_combo, 0, 0)
        slice_layout.addWidget(self.normal1_label, 0, 1, Qt.AlignCenter)
        slice_layout.addWidget(self.normal2_label, 0, 2, Qt.AlignCenter)
        slice_layout.addWidget(self.normal3_label, 0, 3, Qt.AlignCenter)
        slice_layout.addWidget(scalar_label, 0, 4, Qt.AlignCenter)
        slice_layout.addWidget(self.slice_combo, 0, 5)
        slice_layout.addWidget(self.replot_button, 1, 0)
        slice_layout.addWidget(self.normal1_line, 1, 1)
        slice_layout.addWidget(self.normal2_line, 1, 2)
        slice_layout.addWidget(self.normal3_line, 1, 3)
        slice_layout.addWidget(self.scalar_line, 1, 4)
        slice_layout.addWidget(self.slice_button, 1, 5)

        layout.addLayout(camera_layout)
        layout.addLayout(plot_layout)
        layout.addLayout(slice_layout)

        self.setLayout(layout)

    def add_histo(self, histo_dict, clim, normal, origin):

        self.plotter.clear()        
        self.plotter.clear_actors()

        signal = histo_dict['signal']
        labels = histo_dict['labels']

        min_lim = histo_dict['min_lim']
        max_lim = histo_dict['max_lim']
        spacing = histo_dict['spacing']

        P = histo_dict['projection']
        T = histo_dict['transform']
        S = histo_dict['scales']

        grid = pv.ImageData()

        grid.dimensions = np.array(signal.shape)+1

        grid.origin = min_lim
        grid.spacing = spacing
        
        min_bnd = min_lim*S
        max_bnd = max_lim*S

        bounds = np.array([[min_bnd[i], max_bnd[i]] for i in [0,1,2]])
        limits = np.array([[min_lim[i], max_lim[i]] for i in [0,1,2]])

        a = pv._vtk.vtkMatrix3x3()
        b = pv._vtk.vtkMatrix4x4()
        for i in range(3):
            for j in range(3):
                a.SetElement(i,j,T[i,j])
                b.SetElement(i,j,P[i,j])

        trans = np.log10(signal)
        if clim is not None:
            trans[trans < clim[0]] = clim[0]
            trans[trans > clim[1]] = clim[1]

        grid.cell_data['scalars'] = trans.flatten(order='F')

        normal = np.array(normal).astype(float)
        #normal /= np.linalg.norm(normal)

        origin = -normal*origin

        self.clip = self.plotter.add_volume_clip_plane(grid,
                                                       opacity='linear',
                                                       log_scale=False,
                                                       clim=clim,
                                                       normal=normal,
                                                       origin=origin,
                                                       origin_translation=True,
                                                       show_scalar_bar=False,
                                                       normal_rotation=False,
                                                       user_matrix=b)

        prop = self.clip.GetOutlineProperty()
        prop.SetOpacity(0)

        actor = self.plotter.show_grid(xtitle=labels[0],
                                       ytitle=labels[1],
                                       ztitle=labels[2],
                                       font_size=8,
                                       minor_ticks=True)

        actor.SetAxisBaseForX(*T[:,0])
        actor.SetAxisBaseForY(*T[:,1])
        actor.SetAxisBaseForZ(*T[:,2])

        actor.bounds = bounds.ravel()
        actor.SetXAxisRange(limits[0])
        actor.SetYAxisRange(limits[1])
        actor.SetZAxisRange(limits[2])

        axis0_args = *limits[0], actor.n_xlabels, actor.x_label_format
        axis1_args = *limits[1], actor.n_ylabels, actor.y_label_format
        axis2_args = *limits[2], actor.n_zlabels, actor.z_label_format

        axis0_label = pv.plotting.cube_axes_actor.make_axis_labels(*axis0_args)
        axis1_label = pv.plotting.cube_axes_actor.make_axis_labels(*axis1_args)
        axis2_label = pv.plotting.cube_axes_actor.make_axis_labels(*axis2_args)

        actor.SetAxisLabels(0, axis0_label)
        actor.SetAxisLabels(1, axis1_label)
        actor.SetAxisLabels(2, axis2_label)

        self.change_proj()

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

    def view_vector(self, vecs):

        if len(vecs) == 2:
            vec = np.cross(vecs[0],vecs[1])
            self.plotter.view_vector(vecs[0],vec)
        else:
            self.plotter.view_vector(vecs)

    def update_axis_labels(self):

        axes_type = self.view_combo.currentText()

        if axes_type == '[hkl]':
            self.axis1_label.setText('h')
            self.axis2_label.setText('k')
            self.axis3_label.setText('l')
        else:
            self.axis1_label.setText('u')
            self.axis2_label.setText('v')
            self.axis3_label.setText('w')

    def update_normal_labels(self):

        axes_type = self.slice_combo.currentText()

        if axes_type == '[hkl]':
            self.normal1_label.setText('h')
            self.normal2_label.setText('k')
            self.normal3_label.setText('l')
        else:
            self.normal1_label.setText('u')
            self.normal2_label.setText('v')
            self.normal3_label.setText('w')

    def get_manual_axis_indices(self):

        axes_type = self.view_combo.currentText()

        axes = [self.axis1_line, self.axis2_line, self.axis3_line]
        valid_axes = all([axis.hasAcceptableInput() for axis in axes])

        if valid_axes:

            axis1 = float(self.axis1_line.text())
            axis2 = float(self.axis2_line.text())
            axis3 = float(self.axis3_line.text())

            ind = np.array([axis1,axis2,axis3])

            return axes_type, ind

    def get_manual_normal_indices(self):

        axes_type = self.slice_combo.currentText()

        axes = [self.normal1_line, self.normal2_line, self.normal3_line]
        valid_axes = all([axis.hasAcceptableInput() for axis in axes])

        if valid_axes:

            norm1 = float(self.normal1_line.text())
            norm2 = float(self.normal2_line.text())
            norm3 = float(self.normal3_line.text())

            ind = np.array([norm1,norm2,norm3])

            return axes_type, ind

    def get_manual_normal_value(self):

        if self.scalar_line.hasAcceptableInput():

            return float(self.scalar_line.text())

    def get_clim_clip_type(self):
        
        return self.clim_combo.currentText()