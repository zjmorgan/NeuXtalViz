from qtpy.QtWidgets import (QWidget,
                            QFrame,
                            QGridLayout,
                            QHBoxLayout,
                            QVBoxLayout,
                            QPushButton,
                            QLabel,
                            QTabWidget,
                            QCheckBox,
                            QComboBox,
                            QLineEdit,
                            QFileDialog)

from qtpy.QtGui import QDoubleValidator
from PyQt5.QtCore import Qt

import numpy as np
import pyvista as pv

from NeuXtalViz.views.base_view import NeuXtalVizWidget

class VolumeSlicerView(NeuXtalVizWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.tab_widget = QTabWidget(self)

        self.slicer_tab()

        self.layout().addWidget(self.tab_widget)

    def slicer_tab(self):

        slice_tab = QWidget()
        self.tab_widget.addTab(slice_tab, 'Slicer')

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(-100, 100, 5, notation=notation)

        slice_layout = QVBoxLayout()
        plane_layout = QHBoxLayout()
        draw_layout = QHBoxLayout()

        self.clim_combo = QComboBox(self)
        self.clim_combo.addItem('Min/Max')
        self.clim_combo.addItem('μ±3×σ')
        self.clim_combo.addItem('Q₃/Q₁±1.5×IQR')

        self.load_NXS_button = QPushButton('Load NXS', self)
        self.slice_button = QPushButton('Slice Plane', self)

        draw_layout.addWidget(self.slice_button)
        draw_layout.addWidget(self.load_NXS_button)

        self.slice_combo = QComboBox(self)
        self.slice_combo.addItem('Axis 1/2')
        self.slice_combo.addItem('Axis 1/3')
        self.slice_combo.addItem('Axis 2/3')
        self.slice_combo.setCurrentIndex(0)

        scalar_label = QLabel('Slice:', self)

        self.scalar_line = QLineEdit('0.0')

        self.scalar_line.setValidator(validator)

        plane_layout.addWidget(self.clim_combo)
        plane_layout.addWidget(self.slice_combo)
        plane_layout.addWidget(scalar_label)
        plane_layout.addWidget(self.scalar_line)

        slice_layout.addLayout(draw_layout)
        slice_layout.addLayout(plane_layout)
        slice_layout.addStretch(1)

        slice_tab.setLayout(slice_layout)

    def connect_load_NXS(self, load_NXS):

        self.load_NXS_button.clicked.connect(load_NXS)

    def connect_slice(self, slice_plane):

        self.slice_button.clicked.connect(slice_plane)

    def load_NXS_file_dialog(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        filename, _ = QFileDialog.getOpenFileName(self,
                                                  'Load NXS file',
                                                  '',
                                                  'NXS files (*.nxs)',
                                                  options=options)

        return filename

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

        normal /= np.linalg.norm(normal)

        origin = -normal*origin

        self.clip = self.plotter.add_volume_clip_plane(grid,
                                                       opacity='linear',
                                                       log_scale=False,
                                                       clim=clim,
                                                       normal=normal,
                                                       origin=origin,
                                                       origin_translation=True,
                                                       show_scalar_bar=False,
                                                       normal_rotation=True,
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

        self.reset_view()

    def get_normal_value(self):

        if self.scalar_line.hasAcceptableInput():

            return float(self.scalar_line.text())

    def get_clim_clip_type(self):

        return self.clim_combo.currentText()

    def get_slice(self):

        return self.slice_combo.currentText()