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

from PyQt5.QtWidgets import QApplication, QMainWindow

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

        self.layout().addWidget(self.tab_widget)

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
        self.mode_combo.addItem('Ambient')
        self.mode_combo.addItem('Cryogenic')

        self.laue_combo = QComboBox(self)
        self.laue_combo.addItem('None')
        self.laue_combo.addItem('-1')
        self.laue_combo.addItem('2/m')
        self.laue_combo.addItem('mmm')
        self.laue_combo.addItem('4/m')
        self.laue_combo.addItem('4/mmm')
        self.laue_combo.addItem('-3')
        self.laue_combo.addItem('-3m')
        self.laue_combo.addItem('6/m')
        self.laue_combo.addItem('6/mmm')
        self.laue_combo.addItem('m-3')
        self.laue_combo.addItem('m-3m')

        self.load_UB_button = QPushButton('Load UB', self)
        self.optimize_button = QPushButton('Optimize Coverage', self)

        self.wl_min_line = QLineEdit('0.3')
        self.wl_max_line = QLineEdit('3.5')

        wl_label = QLabel('λ:')

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0.2, 10, 5, notation=notation)

        self.wl_min_line.setValidator(validator)
        self.wl_max_line.setValidator(validator)

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
        optimize_layout.addWidget(self.optimize_button)

        settings_layout = QHBoxLayout()

        settings_layout.addWidget(settings_label)
        settings_layout.addWidget(self.settings_line)
        settings_layout.addWidget(self.laue_combo)

        params_layout = QHBoxLayout()

        params_layout.addWidget(wl_label)
        params_layout.addWidget(self.wl_min_line)
        params_layout.addWidget(self.wl_max_line)
        params_layout.addWidget(self.mode_combo)

        result_layout = QVBoxLayout()

        result_layout.addWidget(self.goniometer_table)
        result_layout.addWidget(self.motor_table)

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

    def connect_switch_instrument(self, switch_instrument):

        self.instrument_combo.activated.connect(switch_instrument)

    def connect_optimize(self, optimize):

        self.optimize_button.clicked.connect(optimize)

    def connect_load_UB(self, load_UB):

        self.load_UB_button.clicked.connect(load_UB)

    def connect_wavelength(self, update_wavelength):

        self.wl_min_line.editingFinished.connect(update_wavelength)

    def load_UB_file_dialog(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        filename, _ = QFileDialog.getOpenFileName(self,
                                                  'Load UB file',
                                                  '',
                                                  'UB files (*.mat)',
                                                  options=options)

        return filename

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

        points = coverage_dict['points']
        labels = coverage_dict['labels']

        for point in points:

            mesh = pv.Line(pointa=(0,0,0),
                           pointb=point,
                           resolution=1)

            self.plotter.add_mesh(mesh,
                                  color='k',
                                  style='wireframe',
                                  render_lines_as_tubes=True)

        self.plotter.add_point_labels(points,
                                      labels,
                                      always_visible=True,
                                      point_size=10,
                                      render_points_as_spheres=True)

        self.reset_view()

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.setWindowTitle('Testing Look')

        widget = ExperimentView()
        self.setCentralWidget(widget)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())