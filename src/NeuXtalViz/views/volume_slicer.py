from qtpy.QtWidgets import (QWidget,
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

import numpy as np
import pyvista as pv

from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.transforms import Affine2D

from NeuXtalViz.views.base_view import NeuXtalVizWidget

class VolumeSlicerView(NeuXtalVizWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.tab_widget = QTabWidget(self)

        self.slicer_tab()

        self.layout().addWidget(self.tab_widget, stretch=1)

    def slicer_tab(self):

        slice_tab = QWidget()
        self.tab_widget.addTab(slice_tab, 'Slicer')

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(-100, 100, 5, notation=notation)

        plots_layout = QVBoxLayout()
        slice_params_layout = QHBoxLayout()
        cut_params_layout = QHBoxLayout()
        draw_layout = QHBoxLayout()

        self.clim_combo = QComboBox(self)
        self.clim_combo.addItem('Min/Max')
        self.clim_combo.addItem('μ±3×σ')
        self.clim_combo.addItem('Q₃/Q₁±1.5×IQR')

        self.load_NXS_button = QPushButton('Load NXS', self)
        self.redraw_button = QPushButton('Redraw Volume', self)
        self.slice_button = QPushButton('Slice Plane', self)
        self.cut_button = QPushButton('Cut Line', self)

        draw_layout.addWidget(self.redraw_button)
        draw_layout.addWidget(self.clim_combo)
        draw_layout.addWidget(self.load_NXS_button)

        self.slice_combo = QComboBox(self)
        self.slice_combo.addItem('Axis 1/2')
        self.slice_combo.addItem('Axis 1/3')
        self.slice_combo.addItem('Axis 2/3')
        self.slice_combo.setCurrentIndex(0)

        self.cut_combo = QComboBox(self)
        self.cut_combo.addItem('Axis 1')
        self.cut_combo.addItem('Axis 2')
        self.cut_combo.setCurrentIndex(0)

        slice_label = QLabel('Slice:', self)
        cut_label = QLabel('Cut:', self)

        self.slice_line = QLineEdit('0.0')
        self.slice_line.setValidator(validator)

        self.cut_line = QLineEdit('0.0')
        self.cut_line.setValidator(validator)

        validator = QDoubleValidator(0.0001, 100, 5, notation=notation)

        slice_thickness_label = QLabel('Thickness:', self)
        cut_thickness_label = QLabel('Thickness:', self)

        self.slice_thickness_line = QLineEdit('0.01')
        self.cut_thickness_line = QLineEdit('0.01')

        self.slice_thickness_line.setValidator(validator)
        self.cut_thickness_line.setValidator(validator)

        self.slice_check = QCheckBox('Log', self)
        self.cut_check = QCheckBox('Log', self)

        self.slice_check.setChecked(True)
        self.cut_check.setChecked(True)

        slice_params_layout.addWidget(self.slice_button)
        slice_params_layout.addWidget(self.slice_combo)
        slice_params_layout.addWidget(slice_label)
        slice_params_layout.addWidget(self.slice_line)
        slice_params_layout.addWidget(slice_thickness_label,)
        slice_params_layout.addWidget(self.slice_thickness_line)
        slice_params_layout.addWidget(self.slice_check)

        cut_params_layout.addWidget(self.cut_button)
        cut_params_layout.addWidget(self.cut_combo)
        cut_params_layout.addWidget(cut_label)
        cut_params_layout.addWidget(self.cut_line)
        cut_params_layout.addWidget(cut_thickness_label)
        cut_params_layout.addWidget(self.cut_thickness_line)
        cut_params_layout.addWidget(self.cut_check)

        plots_layout.addLayout(draw_layout)

        self.canvas_slice = FigureCanvas(Figure(constrained_layout=True))
        self.canvas_cut = FigureCanvas(Figure(constrained_layout=True))

        plots_layout.addWidget(NavigationToolbar2QT(self.canvas_slice, self))
        plots_layout.addWidget(self.canvas_slice)
        plots_layout.addLayout(slice_params_layout)
        plots_layout.addWidget(NavigationToolbar2QT(self.canvas_cut, self))
        plots_layout.addWidget(self.canvas_cut)
        plots_layout.addLayout(cut_params_layout)

        self.fig_slice = self.canvas_slice.figure
        self.fig_cut = self.canvas_cut.figure

        self.ax_slice = self.fig_slice.subplots(1, 1)
        self.ax_cut = self.fig_cut.subplots(1, 1)

        self.cb = None

        slice_tab.setLayout(plots_layout)

    def connect_load_NXS(self, load_NXS):

        self.load_NXS_button.clicked.connect(load_NXS)

    def connect_redraw(self, redraw_data):

        self.redraw_button.clicked.connect(redraw_data)

    def connect_slice(self, slice_plane):

        self.slice_button.clicked.connect(slice_plane)

    def connect_cut(self, cut_line):

        self.cut_button.clicked.connect(cut_line)

    def load_NXS_file_dialog(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        filename, _ = QFileDialog.getOpenFileName(self,
                                                  'Load NXS file',
                                                  '',
                                                  'NXS files (*.nxs)',
                                                  options=options)

        return filename

    def add_histo(self, histo_dict, normal, origin):

        self.plotter.clear_plane_widgets()
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

        grid.cell_data['scalars'] = signal.flatten(order='F')

        normal /= np.linalg.norm(normal)

        origin = -normal*origin

        clim = [np.nanmin(signal), np.nanmax(signal)]

        self.clip = self.plotter.add_volume_clip_plane(grid,
                                                       opacity='linear',
                                                       log_scale=False,
                                                       clim=clim,
                                                       normal=normal,
                                                       origin=origin,
                                                       origin_translation=True,
                                                       show_scalar_bar=False,
                                                       normal_rotation=True,
                                                       mapper='smart',
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

    def add_slice(self, slice_dict):

        x = slice_dict['x']
        y = slice_dict['y']

        labels = slice_dict['labels']
        title = slice_dict['title']
        signal = slice_dict['signal']

        scale = 'log' if self.slice_check.isChecked() else 'linear'

        T = slice_dict['transform']
        aspect = slice_dict['aspect']

        transform = Affine2D(T)+self.ax_slice.transData
        self.transform = transform

        self.xlim = np.array([x.min(), x.max()])
        self.ylim = np.array([y.min(), y.max()])

        if self.cb is not None:
            self.cb.remove()

        self.ax_slice.clear()

        im = self.ax_slice.pcolormesh(x,
                                      y,
                                      signal,
                                      norm=scale,
                                      shading='flat',
                                      transform=transform)

        self.ax_slice.set_aspect(aspect)
        self.ax_slice.set_xlabel(labels[0])
        self.ax_slice.set_ylabel(labels[1])
        self.ax_slice.set_title(title)
        self.ax_slice.minorticks_on()

        self.ax_slice.xaxis.get_major_locator().set_params(integer=True)
        self.ax_slice.yaxis.get_major_locator().set_params(integer=True)

        self.cb = self.fig_slice.colorbar(im, ax=self.ax_slice)
        self.cb.minorticks_on()

        self.canvas_slice.draw_idle()
        self.canvas_slice.flush_events()

    def add_cut(self, cut_dict):

        x = cut_dict['x']
        y = cut_dict['y']
        e = cut_dict['e']

        val = cut_dict['value']

        label = cut_dict['label']
        title = cut_dict['title']

        scale = 'log' if self.cut_check.isChecked() else 'linear'

        line_cut = self.get_cut()

        lines = self.ax_slice.get_lines()
        for line in lines:
            line.remove()

        xlim = self.xlim
        ylim = self.ylim

        if line_cut == 'Axis 2':
            line = [val,val], ylim
        else:
            line = xlim, [val,val]

        self.ax_slice.plot(*line, 'w--', linewidth=1, transform=self.transform)

        self.ax_cut.clear()

        self.ax_cut.errorbar(x, y, e)
        self.ax_cut.set_xlabel(label)
        self.ax_cut.set_yscale(scale)
        self.ax_cut.set_title(title)
        self.ax_cut.minorticks_on()

        self.ax_cut.xaxis.get_major_locator().set_params(integer=True)

        self.canvas_cut.draw_idle()
        self.canvas_cut.flush_events()

        self.canvas_slice.draw_idle()
        self.canvas_slice.flush_events()

    def get_slice_value(self):

        if self.slice_line.hasAcceptableInput():

            return float(self.slice_line.text())

    def get_cut_value(self):

        if self.cut_line.hasAcceptableInput():

            return float(self.cut_line.text())

    def get_slice_thickness(self):

        if self.slice_thickness_line.hasAcceptableInput():

            return float(self.slice_thickness_line.text())

    def get_cut_thickness(self):

        if self.cut_thickness_line.hasAcceptableInput():

            return float(self.cut_thickness_line.text())

    def get_clim_clip_type(self):

        return self.clim_combo.currentText()

    def get_slice(self):

        return self.slice_combo.currentText()

    def get_cut(self):

        return self.cut_combo.currentText()