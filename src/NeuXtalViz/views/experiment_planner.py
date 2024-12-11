import sys

from qtpy.QtWidgets import (QWidget,
                            QTableWidget,
                            QTableWidgetItem,
                            QHeaderView,
                            QFrame,
                            QGridLayout,
                            QHBoxLayout,
                            QVBoxLayout,
                            QFormLayout,
                            QPushButton,
                            QCheckBox,
                            QComboBox,
                            QLabel,
                            QLineEdit,
                            QTabWidget,
                            QFileDialog)

from qtpy.QtGui import QDoubleValidator, QIntValidator
from PyQt5.QtCore import Qt

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure

from NeuXtalViz.views.base_view import NeuXtalVizWidget

class ExperimentView(NeuXtalVizWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.tab_widget = QTabWidget(self)

        self.coverage_tab()

        self.layout().addWidget(self.tab_widget, stretch=1)

    def coverage_tab(self):

        cov_tab = QWidget()
        self.tab_widget.addTab(cov_tab, 'Coverage')

        coverage_layout = QVBoxLayout()

        self.instrument_combo = QComboBox(self)
        self.instrument_combo.addItem('TOPAZ')
        self.instrument_combo.addItem('MANDI')
        self.instrument_combo.addItem('CORELLI')
        self.instrument_combo.addItem('SNAP')
        self.instrument_combo.addItem('WAND²')
        self.instrument_combo.addItem('DEMAND')

        self.mode_combo = QComboBox(self)

        self.crystal_combo = QComboBox(self)
        self.point_group_combo = QComboBox(self)
        self.lattice_centering_combo = QComboBox(self)

        self.crystal_combo.addItem('Triclinic')
        self.crystal_combo.addItem('Monoclinic')
        self.crystal_combo.addItem('Orthorhombic')
        self.crystal_combo.addItem('Tetragonal')
        self.crystal_combo.addItem('Trigonal/Rhombohedral')
        self.crystal_combo.addItem('Trigonal/Hexagonal')
        self.crystal_combo.addItem('Hexagonal')
        self.crystal_combo.addItem('Cubic')

        self.load_UB_button = QPushButton('Load UB', self)
        self.optimize_button = QPushButton('Optimize Coverage', self)

        self.wl_min_line = QLineEdit('0.4')
        self.wl_max_line = QLineEdit('3.5')

        self.d_min_line = QLineEdit('0.7')

        wl_label = QLabel('λ:')
        d_min_label = QLabel('d(min):')

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0.2, 10, 5, notation=notation)

        self.wl_min_line.setValidator(validator)
        self.wl_max_line.setValidator(validator)

        validator = QDoubleValidator(0.4, 10, 5, notation=notation)

        self.d_min_line.setValidator(validator)

        angstrom_label = QLabel('Å')

        settings_label = QLabel('Settings')
        self.settings_line = QLineEdit('20')

        validator = QIntValidator(1, 1000)

        self.settings_line.setValidator(validator)

        resize = QHeaderView.Stretch

        self.goniometer_table = QTableWidget()

        self.goniometer_table.setRowCount(0)
        self.goniometer_table.setColumnCount(3)

        labels = ['Motor', 'Min', 'Max']

        self.goniometer_table.horizontalHeader().setStretchLastSection(True)
        self.goniometer_table.horizontalHeader().setSectionResizeMode(resize)
        self.goniometer_table.setHorizontalHeaderLabels(labels)

        self.motor_table = QTableWidget()
        self.motor_table = QTableWidget()

        self.motor_table.setRowCount(0)
        self.motor_table.setColumnCount(2)

        labels = ['Motor', 'Value']

        self.motor_table.horizontalHeader().setStretchLastSection(True)
        self.motor_table.horizontalHeader().setSectionResizeMode(resize)
        self.motor_table.setHorizontalHeaderLabels(labels)

        optimize_layout = QHBoxLayout()

        optimize_layout.addWidget(self.load_UB_button)
        optimize_layout.addWidget(self.instrument_combo)
        optimize_layout.addWidget(self.mode_combo)
        optimize_layout.addWidget(self.optimize_button)

        settings_layout = QHBoxLayout()

        settings_layout.addWidget(settings_label)
        settings_layout.addWidget(self.settings_line)
        settings_layout.addWidget(self.crystal_combo)
        settings_layout.addWidget(self.point_group_combo)
        settings_layout.addWidget(self.lattice_centering_combo)

        params_layout = QHBoxLayout()

        params_layout.addWidget(wl_label)
        params_layout.addWidget(self.wl_min_line)
        params_layout.addWidget(self.wl_max_line)
        params_layout.addWidget(d_min_label)
        params_layout.addWidget(self.d_min_line)
        params_layout.addWidget(angstrom_label)

        result_layout = QVBoxLayout()

        values_tab = QTabWidget()

        goniometer_tab = QWidget()
        motor_tab = QWidget()

        goniometer_layout = QVBoxLayout()
        motor_layout = QVBoxLayout()

        goniometer_layout.addWidget(self.goniometer_table)
        motor_layout.addWidget(self.motor_table)

        goniometer_tab.setLayout(goniometer_layout)
        motor_tab.setLayout(motor_layout)

        values_tab.addTab(goniometer_tab, 'Goniometers')
        values_tab.addTab(motor_tab, 'Motors')

        result_layout.addWidget(values_tab)

        self.canvas = FigureCanvas(Figure(tight_layout=True))

        result_layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        result_layout.addWidget(self.canvas)

        fig = self.canvas.figure

        self.ax = fig.subplots(1, 1)
        self.ax.set_xlim(0,1)
        self.ax.set_ylim(0,100)
        self.ax.minorticks_on()
        self.ax.set_xlabel('Iteration [#]')
        self.ax.set_ylabel('Coverage [%]')

        coverage_layout.addLayout(optimize_layout)
        coverage_layout.addLayout(settings_layout)
        coverage_layout.addLayout(params_layout)
        coverage_layout.addLayout(result_layout)

        cov_tab.setLayout(coverage_layout)

    def connect_switch_crystal(self, switch_crystal):

        self.crystal_combo.activated.connect(switch_crystal)

    def connect_switch_point_group(self, switch_group):

        self.point_group_combo.activated.connect(switch_group)

    def connect_switch_instrument(self, switch_instrument):

        self.instrument_combo.activated.connect(switch_instrument)

    def connect_update_goniometer(self, update_goniometer):

        self.mode_combo.activated.connect(update_goniometer)

    def connect_optimize(self, optimize):

        self.optimize_button.clicked.connect(optimize)

    def connect_load_UB(self, load_UB):

        self.load_UB_button.clicked.connect(load_UB)

    def connect_wavelength(self, update_wavelength):

        self.wl_min_line.editingFinished.connect(update_wavelength)

    def load_UB_file_dialog(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filename, _ = file_dialog.getOpenFileName(self,
                                                  'Load UB file',
                                                  '',
                                                  'UB files (*.mat)',
                                                  options=options)

        return filename

    def get_crystal_system(self):

        return self.crystal_combo.currentText()

    def set_point_groups(self, groups):

        self.point_group_combo.clear()
        for group in groups:
            self.point_group_combo.addItem(group)

    def get_point_group(self):

        return self.point_group_combo.currentText()

    def set_lattice_centerings(self, centerings):

        self.lattice_centering_combo.clear()
        for centering in centerings:
            self.lattice_centering_combo.addItem(centering)

    def get_lattice_centering(self):

        return self.lattice_centering_combo.currentText()

    def get_mode(self):

        return self.mode_combo.currentText()

    def set_modes(self, modes):

        self.mode_combo.clear()
        for mode in modes:
            self.mode_combo.addItem(mode)

    def get_mode(self):

        return self.mode_combo.currentText()

    def set_wavelength(self, wavelength):

        if type(wavelength) is list:
            self.wl_min_line.setText(str(wavelength[0]))
            self.wl_max_line.setText(str(wavelength[1]))
            self.wl_max_line.setEnabled(True)
        else:
            self.wl_min_line.setText(str(wavelength))
            self.wl_max_line.setText(str(wavelength))
            self.wl_max_line.setEnabled(False)

    def get_wavelength(self):

        params = self.wl_min_line, self.wl_max_line

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:

            return [float(param.text()) for param in params]

    def update_wavelength(self, lamda_min):

        if not self.wl_max_line.isEnabled():
            self.wl_max_line.setText(str(lamda_min))

    def update_tables(self, goniometers, motors):

        self.goniometer_table.setRowCount(0)
        self.goniometer_table.setRowCount(len(goniometers))

        for row, gon in enumerate(goniometers):
            angle, amin, amax = gon
            amin, amax = str(amin), str(amax)
            self.goniometer_table.setItem(row, 0, QTableWidgetItem(angle))
            self.goniometer_table.setItem(row, 1, QTableWidgetItem(amin))
            self.goniometer_table.setItem(row, 2, QTableWidgetItem(amax))

        self.motor_table.setRowCount(0)
        self.motor_table.setRowCount(len(motors))

        for row, mot in enumerate(motors):
            setting, val = mot
            val = str(val)
            self.motor_table.setItem(row, 0, QTableWidgetItem(setting))
            self.motor_table.setItem(row, 1, QTableWidgetItem(val))

        for row in range(self.goniometer_table.rowCount()):
            item = self.goniometer_table.item(row, 0)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)

        for row in range(self.motor_table.rowCount()):
            item = self.motor_table.item(row, 0)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)

    def get_instrument(self):

        return self.instrument_combo.currentText()

    def get_laue_symmetry(self):

        return self.laue_combo.currentText()

    def get_motors(self):

        logs = {}
        for row in range(self.motor_table.rowCount()):
            setting = self.motor_table.item(row, 0).text()
            logs[setting] = float(self.motor_table.item(row, 1).text())

        return logs

    # def get_cluster_parameters(self):

    #     params = [self.param_eps_line, self.param_min_line]
    #     valid_params = all([param.hasAcceptableInput() for param in params])

    #     if valid_params:

    #         return float(self.param_eps_line.text()), \
    #                 int(self.param_min_line.text())

    def add_coverage(self, coverage_dict):

        self.plotter.clear_actors()

        signal = coverage_dict['signal']
        spacing = coverage_dict['spacing']
        scalars = coverage_dict['scalars']
        limits = coverage_dict['limits']

        grid = pv.ImageData(dimensions=signal.shape)
        grid['scalars'] = signal.flatten(order='F')

        grid.spacing = spacing
        grid.origin = limits[0::2]

        actor = self.plotter.add_volume(grid,
                                        scalars=scalars,
                                        mapper='smart')
        labels = ['Qx', 'Qy', 'Qz']

        actor = self.plotter.show_grid(xtitle=labels[0],
                                       ytitle=labels[1],
                                       ztitle=labels[2],
                                       font_size=8,
                                       minor_ticks=True)

        actor.SetXAxisRange(limits[0:2])
        actor.SetYAxisRange(limits[2:4])
        actor.SetZAxisRange(limits[4:6])

        axis0_args = *limits[0:2], actor.n_xlabels, actor.x_label_format
        axis1_args = *limits[2:4], actor.n_ylabels, actor.y_label_format
        axis2_args = *limits[4:6], actor.n_zlabels, actor.z_label_format

        axis0_label = pv.plotting.cube_axes_actor.make_axis_labels(*axis0_args)
        axis1_label = pv.plotting.cube_axes_actor.make_axis_labels(*axis1_args)
        axis2_label = pv.plotting.cube_axes_actor.make_axis_labels(*axis2_args)

        actor.SetAxisLabels(0, axis0_label)
        actor.SetAxisLabels(1, axis1_label)
        actor.SetAxisLabels(2, axis2_label)

        # points = coverage_dict['points']
        # labels = coverage_dict['labels']

        # for point in points:

        #     mesh = pv.Line(pointa=(0,0,0),
        #                    pointb=point,
        #                    resolution=1)

        #     self.plotter.add_mesh(mesh,
        #                           color='k',
        #                           style='wireframe',
        #                           render_lines_as_tubes=True)

        # self.plotter.add_point_labels(points,
        #                               labels,
        #                               always_visible=True,
        #                               point_size=10,
        #                               render_points_as_spheres=True)

        self.reset_view()
