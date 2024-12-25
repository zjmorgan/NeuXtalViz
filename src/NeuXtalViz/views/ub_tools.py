import re

import numpy as np
import pyvista as pv

from qtpy.QtWidgets import (
    QWidget,
    QTableWidget,
    QTableWidgetItem,
    QHeaderView,
    QLineEdit,
    QLabel,
    QComboBox,
    QPushButton,
    QCheckBox,
    QHBoxLayout,
    QVBoxLayout,
    QGridLayout,
    QTabWidget,
    QFileDialog,
)

from qtpy.QtGui import QDoubleValidator, QIntValidator
from qtpy.QtCore import Qt, Signal

from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.transforms import Affine2D
from matplotlib.ticker import FormatStrFormatter

# from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.axisartist import Axes, GridHelperCurveLinear
from mpl_toolkits.axisartist.grid_finder import ExtremeFinderSimple, MaxNLocator

from NeuXtalViz.views.base_view import NeuXtalVizWidget

cmaps = {
    "Sequential": "viridis",
    "Binary": "binary",
    "Diverging": "bwr",
    "Rainbow": "turbo",
}


class UBView(NeuXtalVizWidget):
    roi_ready = Signal()
    index_ready = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)

        self.tab_widget = QTabWidget(self)

        self.parameters_tab()
        self.table_tab()
        self.verify_tab()

        self.layout().addWidget(self.tab_widget, stretch=1)

        self.last_highlight = None

        self.x_min, self.x_max = None, None
        self.y_min, self.y_max = None, None

    def parameters_tab(self):
        ub_peaks_tab = QWidget()
        self.tab_widget.addTab(ub_peaks_tab, "Parameters")

        ub_layout = QVBoxLayout()

        self.save_q_button = QPushButton("Save Q", self)
        self.load_q_button = QPushButton("Load Q", self)

        self.save_peaks_button = QPushButton("Save Peaks", self)
        self.load_peaks_button = QPushButton("Load Peaks", self)

        self.save_ub_button = QPushButton("Save UB", self)
        self.load_ub_button = QPushButton("Load UB", self)

        convert_io_layout = QHBoxLayout()

        convert_io_layout.addStretch(1)
        convert_io_layout.addWidget(self.save_q_button)
        convert_io_layout.addWidget(self.load_q_button)

        peaks_io_layout = QHBoxLayout()

        peaks_io_layout.addStretch(1)
        peaks_io_layout.addWidget(self.save_peaks_button)
        peaks_io_layout.addWidget(self.load_peaks_button)

        ub_io_layout = QHBoxLayout()

        ub_io_layout.addStretch(1)
        ub_io_layout.addWidget(self.save_ub_button)
        ub_io_layout.addWidget(self.load_ub_button)

        convert_tab = self.__init_convert_tab()
        peaks_tab = self.__init_peaks_tab()
        ub_tab = self.__init_ub_tab()
        values_tab = self.__init_values_tab()

        ub_layout.addWidget(convert_tab)
        ub_layout.addLayout(convert_io_layout)

        ub_layout.addWidget(peaks_tab)
        ub_layout.addLayout(peaks_io_layout)

        ub_layout.addWidget(ub_tab)
        ub_layout.addLayout(ub_io_layout)

        ub_layout.addWidget(values_tab)

        ub_peaks_tab.setLayout(ub_layout)

    def __init_values_tab(self):
        values_tab = QTabWidget()

        parameters_tab = QWidget()
        orientation_tab = QWidget()
        satellite_tab = QWidget()

        self.a_line = QLineEdit()
        self.b_line = QLineEdit()
        self.c_line = QLineEdit()

        self.alpha_line = QLineEdit()
        self.beta_line = QLineEdit()
        self.gamma_line = QLineEdit()

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0.1, 1000, 4, notation=notation)

        self.a_line.setValidator(validator)
        self.b_line.setValidator(validator)
        self.c_line.setValidator(validator)

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(10, 170, 4, notation=notation)

        self.alpha_line.setValidator(validator)
        self.beta_line.setValidator(validator)
        self.gamma_line.setValidator(validator)

        a_label = QLabel("a:")
        b_label = QLabel("b:")
        c_label = QLabel("c:")

        alpha_label = QLabel("α:")
        beta_label = QLabel("β:")
        gamma_label = QLabel("γ:")

        angstrom_label = QLabel("Å")
        degree_label = QLabel("°")

        parameters_layout = QGridLayout()

        parameters_layout.addWidget(a_label, 0, 0)
        parameters_layout.addWidget(self.a_line, 0, 1)
        parameters_layout.addWidget(b_label, 0, 2)
        parameters_layout.addWidget(self.b_line, 0, 3)
        parameters_layout.addWidget(c_label, 0, 4)
        parameters_layout.addWidget(self.c_line, 0, 5)
        parameters_layout.addWidget(angstrom_label, 0, 6)
        parameters_layout.addWidget(alpha_label, 1, 0)
        parameters_layout.addWidget(self.alpha_line, 1, 1)
        parameters_layout.addWidget(beta_label, 1, 2)
        parameters_layout.addWidget(self.beta_line, 1, 3)
        parameters_layout.addWidget(gamma_label, 1, 4)
        parameters_layout.addWidget(self.gamma_line, 1, 5)
        parameters_layout.addWidget(degree_label, 1, 6)

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(-5, 5, 4, notation=notation)

        self.dh1_line = QLineEdit("0.0")
        self.dk1_line = QLineEdit("0.0")
        self.dl1_line = QLineEdit("0.0")

        self.dh2_line = QLineEdit("0.0")
        self.dk2_line = QLineEdit("0.0")
        self.dl2_line = QLineEdit("0.0")

        self.dh3_line = QLineEdit("0.0")
        self.dk3_line = QLineEdit("0.0")
        self.dl3_line = QLineEdit("0.0")

        self.dh1_line.setValidator(validator)
        self.dk1_line.setValidator(validator)
        self.dl1_line.setValidator(validator)

        self.dh2_line.setValidator(validator)
        self.dk2_line.setValidator(validator)
        self.dl2_line.setValidator(validator)

        self.dh3_line.setValidator(validator)
        self.dk3_line.setValidator(validator)
        self.dl3_line.setValidator(validator)

        self.max_order_line = QLineEdit("0")

        mod_vec1_label = QLabel("1:")
        mod_vec2_label = QLabel("2:")
        mod_vec3_label = QLabel("3:")

        dh_label = QLabel("Δh")
        dk_label = QLabel("Δk")
        dl_label = QLabel("Δl")

        max_order_label = QLabel("Max Order")

        self.cross_box = QCheckBox("Cross Terms", self)
        self.cross_box.setChecked(False)

        satellite_layout = QGridLayout()

        satellite_layout.addWidget(dh_label, 0, 1, Qt.AlignCenter)
        satellite_layout.addWidget(dk_label, 0, 2, Qt.AlignCenter)
        satellite_layout.addWidget(dl_label, 0, 3, Qt.AlignCenter)
        satellite_layout.addWidget(max_order_label, 0, 4, Qt.AlignCenter)
        satellite_layout.addWidget(mod_vec1_label, 1, 0)
        satellite_layout.addWidget(self.dh1_line, 1, 1)
        satellite_layout.addWidget(self.dk1_line, 1, 2)
        satellite_layout.addWidget(self.dl1_line, 1, 3)
        satellite_layout.addWidget(self.max_order_line, 1, 4)
        satellite_layout.addWidget(mod_vec2_label, 2, 0)
        satellite_layout.addWidget(self.dh2_line, 2, 1)
        satellite_layout.addWidget(self.dk2_line, 2, 2)
        satellite_layout.addWidget(self.dl2_line, 2, 3)
        satellite_layout.addWidget(self.cross_box, 2, 4)
        satellite_layout.addWidget(mod_vec3_label, 3, 0)
        satellite_layout.addWidget(self.dh3_line, 3, 1)
        satellite_layout.addWidget(self.dk3_line, 3, 2)
        satellite_layout.addWidget(self.dl3_line, 3, 3)

        x_label = QLabel("x:")
        y_label = QLabel("y:")
        z_label = QLabel("z:")

        a_star_label = QLabel("a*")
        b_star_label = QLabel("b*")
        c_star_label = QLabel("c*")

        self.uh_line = QLineEdit()
        self.uk_line = QLineEdit()
        self.ul_line = QLineEdit()

        self.vh_line = QLineEdit()
        self.vk_line = QLineEdit()
        self.vl_line = QLineEdit()

        self.wh_line = QLineEdit()
        self.wk_line = QLineEdit()
        self.wl_line = QLineEdit()

        self.uh_line.setReadOnly(False)
        self.uk_line.setReadOnly(False)
        self.ul_line.setReadOnly(False)

        self.vh_line.setReadOnly(False)
        self.vk_line.setReadOnly(False)
        self.vl_line.setReadOnly(False)

        self.wh_line.setReadOnly(False)
        self.wk_line.setReadOnly(False)
        self.wl_line.setReadOnly(False)

        orientation_layout = QGridLayout()

        orientation_layout.addWidget(a_star_label, 0, 1, Qt.AlignCenter)
        orientation_layout.addWidget(b_star_label, 0, 2, Qt.AlignCenter)
        orientation_layout.addWidget(c_star_label, 0, 3, Qt.AlignCenter)
        orientation_layout.addWidget(y_label, 1, 0)
        orientation_layout.addWidget(self.wh_line, 1, 1)
        orientation_layout.addWidget(self.wk_line, 1, 2)
        orientation_layout.addWidget(self.wl_line, 1, 3)
        orientation_layout.addWidget(z_label, 2, 0)
        orientation_layout.addWidget(self.uh_line, 2, 1)
        orientation_layout.addWidget(self.uk_line, 2, 2)
        orientation_layout.addWidget(self.ul_line, 2, 3)
        orientation_layout.addWidget(x_label, 3, 0)
        orientation_layout.addWidget(self.vh_line, 3, 1)
        orientation_layout.addWidget(self.vk_line, 3, 2)
        orientation_layout.addWidget(self.vl_line, 3, 3)

        lattice_layout = QVBoxLayout()

        lattice_layout.addLayout(parameters_layout)
        lattice_layout.addStretch(1)

        parameters_tab.setLayout(lattice_layout)
        orientation_tab.setLayout(orientation_layout)
        satellite_tab.setLayout(satellite_layout)

        values_tab.addTab(parameters_tab, "Lattice Parameters")
        values_tab.addTab(orientation_tab, "Sample Orientation")
        values_tab.addTab(satellite_tab, "Modulation Parameters")

        return values_tab

    def __init_convert_tab(self):
        convert_tab = QTabWidget()

        convert_to_q_tab = QWidget()
        convert_to_q_tab_layout = QVBoxLayout()

        experiment_params_layout = QHBoxLayout()
        run_params_layout = QHBoxLayout()
        wavelength_params_layout = QHBoxLayout()
        instrument_params_layout = QGridLayout()

        self.instrument_combo = QComboBox(self)
        self.instrument_combo.addItem("TOPAZ")
        self.instrument_combo.addItem("MANDI")
        self.instrument_combo.addItem("CORELLI")
        self.instrument_combo.addItem("SNAP")
        self.instrument_combo.addItem("WAND²")
        self.instrument_combo.addItem("DEMAND")

        ipts_label = QLabel("IPTS:")
        exp_label = QLabel("Experiment:")
        run_label = QLabel("Runs:")
        filter_time_label = QLabel("Time Stop [s]:")
        angstrom_label = QLabel("Å")

        validator = QIntValidator(1, 1000000000, self)

        self.runs_line = QLineEdit("")

        self.ipts_line = QLineEdit("")
        self.ipts_line.setValidator(validator)

        self.exp_line = QLineEdit("")
        self.exp_line.setValidator(validator)

        self.cal_line = QLineEdit("")
        self.tube_line = QLineEdit("")

        self.wl_min_line = QLineEdit("0.3")
        self.wl_max_line = QLineEdit("3.5")

        wl_label = QLabel("λ:")

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0.2, 10, 5, notation=notation)

        self.wl_min_line.setValidator(validator)
        self.wl_max_line.setValidator(validator)

        validator = QIntValidator(1, 1000, self)

        self.filter_time_line = QLineEdit("300")
        self.filter_time_line.setValidator(validator)

        self.cal_browse_button = QPushButton("Detector", self)
        self.tube_browse_button = QPushButton("Tube", self)

        experiment_params_layout.addWidget(self.instrument_combo)
        experiment_params_layout.addWidget(ipts_label)
        experiment_params_layout.addWidget(self.ipts_line)
        experiment_params_layout.addWidget(exp_label)
        experiment_params_layout.addWidget(self.exp_line)

        run_params_layout.addWidget(run_label)
        run_params_layout.addWidget(self.runs_line)

        wavelength_params_layout.addWidget(wl_label)
        wavelength_params_layout.addWidget(self.wl_min_line)
        wavelength_params_layout.addWidget(self.wl_max_line)
        wavelength_params_layout.addWidget(angstrom_label)

        instrument_params_layout.addWidget(self.cal_line, 1, 0)
        instrument_params_layout.addWidget(self.cal_browse_button, 1, 1)
        instrument_params_layout.addWidget(self.tube_line, 2, 0)
        instrument_params_layout.addWidget(self.tube_browse_button, 2, 1)

        self.convert_to_q_button = QPushButton("Convert", self)

        self.lorentz_box = QCheckBox("Lorentz Correction", self)
        self.lorentz_box.setChecked(True)

        convert_to_q_action_layout = QHBoxLayout()
        convert_to_q_action_layout.addWidget(self.convert_to_q_button)
        convert_to_q_action_layout.addWidget(self.lorentz_box)
        convert_to_q_action_layout.addWidget(filter_time_label)
        convert_to_q_action_layout.addWidget(self.filter_time_line)
        convert_to_q_action_layout.addStretch(1)
        convert_to_q_action_layout.addLayout(wavelength_params_layout)

        convert_to_q_tab_layout.addLayout(experiment_params_layout)
        convert_to_q_tab_layout.addLayout(run_params_layout)
        convert_to_q_tab_layout.addLayout(wavelength_params_layout)
        convert_to_q_tab_layout.addLayout(instrument_params_layout)
        convert_to_q_tab_layout.addStretch(1)
        convert_to_q_tab_layout.addLayout(convert_to_q_action_layout)

        convert_to_q_tab.setLayout(convert_to_q_tab_layout)

        convert_tab.addTab(convert_to_q_tab, "Convert To Q")

        return convert_tab

    def __init_peaks_tab(self):
        peaks_tab = QTabWidget()

        find_tab = QWidget()
        find_tab_layout = QVBoxLayout()

        max_peaks_label = QLabel("Max Peaks:")
        min_distance_label = QLabel("Min Distance:")
        density_threshold_label = QLabel("Min Density:")
        find_edge_label = QLabel("Edge Pixels:")
        distance_unit_label = QLabel("Å⁻¹")

        validator = QIntValidator(10, 1000, self)

        self.max_peaks_line = QLineEdit("100")
        self.max_peaks_line.setValidator(validator)

        validator = QIntValidator(1, 100000, self)

        self.density_threshold_line = QLineEdit("100")
        self.density_threshold_line.setValidator(validator)

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0.01, 10, 4, notation=notation)

        self.min_distance_line = QLineEdit("0.2")
        self.min_distance_line.setValidator(validator)

        validator = QIntValidator(0, 64, self)

        self.find_edge_line = QLineEdit("0")
        self.find_edge_line.setValidator(validator)

        find_params_layout = QGridLayout()

        find_params_layout.addWidget(max_peaks_label, 0, 0)
        find_params_layout.addWidget(self.max_peaks_line, 0, 1)
        find_params_layout.addWidget(density_threshold_label, 0, 2)
        find_params_layout.addWidget(self.density_threshold_line, 0, 3)
        find_params_layout.addWidget(min_distance_label, 1, 0)
        find_params_layout.addWidget(self.min_distance_line, 1, 1)
        find_params_layout.addWidget(distance_unit_label, 1, 2)
        find_params_layout.addWidget(find_edge_label, 2, 0)
        find_params_layout.addWidget(self.find_edge_line, 2, 1)

        self.find_button = QPushButton("Find", self)

        find_action_layout = QHBoxLayout()
        find_action_layout.addWidget(self.find_button)
        find_action_layout.addStretch(1)

        find_tab_layout.addLayout(find_params_layout)
        find_tab_layout.addStretch(1)
        find_tab_layout.addLayout(find_action_layout)

        find_tab.setLayout(find_tab_layout)

        index_tab = QWidget()
        index_tab_layout = QVBoxLayout()

        index_tolerance_label = QLabel("Tolerance:")

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0.01, 1, 5, notation=notation)

        self.index_sat_box = QCheckBox("Satellite", self)
        self.index_sat_box.setChecked(False)

        self.index_tolerance_line = QLineEdit("0.1")
        self.index_tolerance_line.setValidator(validator)

        self.index_sat_tolerance_line = QLineEdit("0.1")
        self.index_sat_tolerance_line.setValidator(validator)

        index_params_layout = QGridLayout()

        index_params_layout.addWidget(index_tolerance_label, 0, 0)
        index_params_layout.addWidget(self.index_tolerance_line, 0, 1)
        index_params_layout.addWidget(self.index_sat_tolerance_line, 0, 2)
        index_params_layout.addWidget(self.index_sat_box, 1, 2)

        self.round_box = QCheckBox("Round hkl", self)
        self.round_box.setChecked(True)

        self.index_button = QPushButton("Index", self)

        index_action_layout = QHBoxLayout()
        index_action_layout.addWidget(self.index_button)
        index_action_layout.addWidget(self.round_box)
        index_action_layout.addStretch(1)

        index_tab_layout.addLayout(index_params_layout)
        index_tab_layout.addStretch(1)
        index_tab_layout.addLayout(index_action_layout)

        index_tab.setLayout(index_tab_layout)

        centering_label = QLabel("Centering:")

        self.centering_combo = QComboBox(self)
        self.centering_combo.addItem("P")
        self.centering_combo.addItem("I")
        self.centering_combo.addItem("F")
        self.centering_combo.addItem("Robv")
        self.centering_combo.addItem("Rrev")
        self.centering_combo.addItem("A")
        self.centering_combo.addItem("B")
        self.centering_combo.addItem("C")
        self.centering_combo.addItem("H")

        min_d_unit_label = QLabel("Å")

        min_d_label = QLabel("Min d-spacing:")
        predict_edge_label = QLabel("Edge Pixels:")

        self.predict_sat_box = QCheckBox("Satellite", self)
        self.predict_sat_box.setChecked(False)

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0.4, 100, 3, notation=notation)

        self.min_d_line = QLineEdit("0.7")
        self.min_d_line.setValidator(validator)

        self.min_sat_d_line = QLineEdit("1.0")
        self.min_sat_d_line.setValidator(validator)

        validator = QIntValidator(0, 64, self)

        self.predict_edge_line = QLineEdit("0")
        self.predict_edge_line.setValidator(validator)

        predict_tab = QWidget()
        predict_tab_layout = QVBoxLayout()

        predict_params_layout = QGridLayout()

        predict_params_layout.addWidget(centering_label, 0, 0)
        predict_params_layout.addWidget(self.centering_combo, 0, 1)
        predict_params_layout.addWidget(min_d_label, 1, 0)
        predict_params_layout.addWidget(self.min_d_line, 1, 1)
        predict_params_layout.addWidget(self.min_sat_d_line, 1, 2)
        predict_params_layout.addWidget(min_d_unit_label, 1, 3)
        predict_params_layout.addWidget(predict_edge_label, 2, 0)
        predict_params_layout.addWidget(self.predict_edge_line, 2, 1)
        predict_params_layout.addWidget(self.predict_sat_box, 2, 2)

        self.predict_button = QPushButton("Predict", self)

        predict_action_layout = QHBoxLayout()
        predict_action_layout.addWidget(self.predict_button)
        predict_action_layout.addStretch(1)

        predict_tab_layout.addLayout(predict_params_layout)
        predict_tab_layout.addStretch(1)
        predict_tab_layout.addLayout(predict_action_layout)

        predict_tab.setLayout(predict_tab_layout)

        self.centroid_box = QCheckBox("Centroid", self)
        self.centroid_box.setChecked(True)

        self.adaptive_box = QCheckBox("Adaptive Envelope", self)
        self.adaptive_box.setChecked(True)

        radius_label = QLabel("Radius:")
        inner_label = QLabel("Inner Factor:")
        outer_label = QLabel("Outer Factor:")
        radius_unit_label = QLabel("Å⁻¹")

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0, 1, 3, notation=notation)

        self.radius_line = QLineEdit("0.25")
        self.radius_line.setValidator(validator)

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(1, 3, 3, notation=notation)

        self.inner_line = QLineEdit("1.5")
        self.inner_line.setValidator(validator)

        self.outer_line = QLineEdit("2")
        self.outer_line.setValidator(validator)

        integrate_tab = QWidget()
        integrate_tab_layout = QVBoxLayout()

        integrate_params_layout = QGridLayout()

        integrate_params_layout.addWidget(radius_label, 0, 0)
        integrate_params_layout.addWidget(self.radius_line, 0, 1)
        integrate_params_layout.addWidget(radius_unit_label, 0, 2)
        integrate_params_layout.addWidget(inner_label, 2, 0)
        integrate_params_layout.addWidget(self.inner_line, 2, 1)
        integrate_params_layout.addWidget(outer_label, 2, 2)
        integrate_params_layout.addWidget(self.outer_line, 2, 3)

        self.integrate_button = QPushButton("Integrate", self)

        integrate_action_layout = QHBoxLayout()
        integrate_action_layout.addWidget(self.integrate_button)
        integrate_action_layout.addWidget(self.centroid_box)
        integrate_action_layout.addWidget(self.adaptive_box)
        integrate_action_layout.addStretch(1)

        integrate_tab_layout.addLayout(integrate_params_layout)
        integrate_tab_layout.addStretch(1)
        integrate_tab_layout.addLayout(integrate_action_layout)

        integrate_tab.setLayout(integrate_tab_layout)

        self.filter_combo = QComboBox(self)
        self.filter_combo.addItem("I/σ")
        self.filter_combo.addItem("d")
        self.filter_combo.addItem("λ")
        self.filter_combo.addItem("Q")
        self.filter_combo.addItem("h^2+k^2+l^2")
        self.filter_combo.addItem("m^2+n^2+p^2")
        self.filter_combo.addItem("Run #")

        self.comparison_combo = QComboBox(self)
        self.comparison_combo.addItem(">")
        self.comparison_combo.addItem("<")
        self.comparison_combo.addItem(">=")
        self.comparison_combo.addItem("<=")
        self.comparison_combo.addItem("=")
        self.comparison_combo.addItem("!=")

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(-1e6, 1e6, 3, notation=notation)

        self.filter_line = QLineEdit("0")
        self.filter_line.setValidator(validator)

        filter_tab = QWidget()
        filter_tab_layout = QVBoxLayout()

        filter_params_layout = QHBoxLayout()

        filter_params_layout.addWidget(self.filter_combo)
        filter_params_layout.addWidget(self.comparison_combo)
        filter_params_layout.addWidget(self.filter_line)

        self.filter_button = QPushButton("Filter", self)

        filter_action_layout = QHBoxLayout()
        filter_action_layout.addWidget(self.filter_button)
        filter_action_layout.addStretch(1)

        filter_tab_layout.addLayout(filter_params_layout)
        filter_tab_layout.addStretch(1)
        filter_tab_layout.addLayout(filter_action_layout)

        filter_tab.setLayout(filter_tab_layout)

        peaks_tab.addTab(find_tab, "Find Peaks")
        peaks_tab.addTab(index_tab, "Index Peaks")
        peaks_tab.addTab(predict_tab, "Predict Peaks")
        peaks_tab.addTab(integrate_tab, "Integrate Peaks")
        peaks_tab.addTab(filter_tab, "Filter Peaks")

        return peaks_tab

    def __init_ub_tab(self):
        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0.01, 1, 5, notation=notation)

        ub_tab = QTabWidget()

        calculate_tolerance_label = QLabel("Tolerance:")

        self.calculate_tolerance_line = QLineEdit("0.1")
        self.calculate_tolerance_line.setValidator(validator)

        max_scalar_error_label = QLabel("Max Scalar Error:")

        self.max_scalar_error_line = QLineEdit("0.2")
        self.max_scalar_error_line.setValidator(validator)

        calculate_tab = QWidget()
        calculate_tab_layout = QVBoxLayout()

        calculate_params_layout = QGridLayout()

        calculate_params_layout.addWidget(calculate_tolerance_label, 0, 0)
        calculate_params_layout.addWidget(self.calculate_tolerance_line, 0, 1)
        calculate_params_layout.addWidget(max_scalar_error_label, 0, 2)
        calculate_params_layout.addWidget(self.max_scalar_error_line, 0, 3)

        self.conventional_button = QPushButton("Conventional", self)
        self.niggli_button = QPushButton("Niggli", self)
        self.select_button = QPushButton("Select", self)

        self.form_line = QLineEdit("")
        self.form_line.setReadOnly(True)

        form_label = QLabel("Form:")

        min_const_label = QLabel("Min(a,b,c) [Å]:")
        max_const_label = QLabel("Max(a,b,c) [Å]:")

        self.min_const_line = QLineEdit("5")
        self.max_const_line = QLineEdit("15")

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0.1, 1000, 4, notation=notation)

        self.min_const_line.setValidator(validator)
        self.max_const_line.setValidator(validator)

        const_layout = QHBoxLayout()
        const_layout.addWidget(min_const_label)
        const_layout.addWidget(self.min_const_line)
        const_layout.addWidget(max_const_label)
        const_layout.addWidget(self.max_const_line)

        calculate_action_layout = QHBoxLayout()
        calculate_action_layout.addWidget(self.conventional_button)
        calculate_action_layout.addStretch(1)
        calculate_action_layout.addLayout(const_layout)
        calculate_action_layout.addWidget(self.niggli_button)
        calculate_action_layout.addWidget(form_label)
        calculate_action_layout.addWidget(self.form_line)
        calculate_action_layout.addWidget(self.select_button)

        stretch = QHeaderView.Stretch

        self.cell_table = QTableWidget()
        self.cell_table.setRowCount(0)
        self.cell_table.setColumnCount(9)

        header = ["Error", "Bravais", "a", "b", "c", "α", "β", "γ", "V"]

        self.cell_table.horizontalHeader().setSectionResizeMode(stretch)
        self.cell_table.setHorizontalHeaderLabels(header)
        self.cell_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.cell_table.setSelectionBehavior(QTableWidget.SelectRows)
        # self.cell_table.setSortingEnabled(True)

        calculate_tab_layout.addLayout(calculate_params_layout)
        calculate_tab_layout.addWidget(self.cell_table)
        calculate_tab_layout.addLayout(calculate_action_layout)

        calculate_tab.setLayout(calculate_tab_layout)

        transform_tolerance_label = QLabel("Tolerance:")

        self.transform_tolerance_line = QLineEdit("0.1")
        self.transform_tolerance_line.setValidator(validator)

        self.lattice_combo = QComboBox(self)
        self.lattice_combo.addItem("Triclinic")
        self.lattice_combo.addItem("Monoclinic")
        self.lattice_combo.addItem("Orthorhombic")
        self.lattice_combo.addItem("Tetragonal")
        self.lattice_combo.addItem("Rhombohedral")
        self.lattice_combo.addItem("Hexagonal")
        self.lattice_combo.addItem("Cubic")

        self.symmetry_combo = QComboBox(self)
        self.symmetry_combo.addItem("x,y,z")
        self.symmetry_combo.addItem("-x,-y,-z")

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(-10, 10, 5, notation=notation)

        self.T11_line = QLineEdit("1")
        self.T12_line = QLineEdit("0")
        self.T13_line = QLineEdit("0")

        self.T21_line = QLineEdit("0")
        self.T22_line = QLineEdit("1")
        self.T23_line = QLineEdit("0")

        self.T31_line = QLineEdit("0")
        self.T32_line = QLineEdit("0")
        self.T33_line = QLineEdit("1")

        self.T11_line.setValidator(validator)
        self.T12_line.setValidator(validator)
        self.T13_line.setValidator(validator)

        self.T21_line.setValidator(validator)
        self.T22_line.setValidator(validator)
        self.T23_line.setValidator(validator)

        self.T31_line.setValidator(validator)
        self.T32_line.setValidator(validator)
        self.T33_line.setValidator(validator)

        hp_label = QLabel("h′:")
        kp_label = QLabel("k′:")
        lp_label = QLabel("l′:")

        h_label = QLabel("h")
        k_label = QLabel("k")
        l_label = QLabel("l")

        transform_tab = QWidget()
        transform_tab_layout = QVBoxLayout()

        transform_params_layout = QHBoxLayout()

        transform_params_layout.addWidget(transform_tolerance_label)
        transform_params_layout.addWidget(self.transform_tolerance_line)
        transform_params_layout.addWidget(self.lattice_combo)
        transform_params_layout.addWidget(self.symmetry_combo)

        transform_matrix_layout = QGridLayout()

        transform_matrix_layout.addWidget(h_label, 0, 1, Qt.AlignCenter)
        transform_matrix_layout.addWidget(k_label, 0, 2, Qt.AlignCenter)
        transform_matrix_layout.addWidget(l_label, 0, 3, Qt.AlignCenter)
        transform_matrix_layout.addWidget(hp_label, 1, 0)
        transform_matrix_layout.addWidget(self.T11_line, 1, 1)
        transform_matrix_layout.addWidget(self.T12_line, 1, 2)
        transform_matrix_layout.addWidget(self.T13_line, 1, 3)
        transform_matrix_layout.addWidget(kp_label, 2, 0)
        transform_matrix_layout.addWidget(self.T21_line, 2, 1)
        transform_matrix_layout.addWidget(self.T22_line, 2, 2)
        transform_matrix_layout.addWidget(self.T23_line, 2, 3)
        transform_matrix_layout.addWidget(lp_label, 3, 0)
        transform_matrix_layout.addWidget(self.T31_line, 3, 1)
        transform_matrix_layout.addWidget(self.T32_line, 3, 2)
        transform_matrix_layout.addWidget(self.T33_line, 3, 3)

        self.transform_button = QPushButton("Transform", self)

        transform_action_layout = QHBoxLayout()
        transform_action_layout.addWidget(self.transform_button)
        transform_action_layout.addStretch(1)
        transform_action_layout.addLayout(transform_params_layout)

        transform_tab_layout.addLayout(transform_matrix_layout)
        transform_tab_layout.addStretch(1)
        transform_tab_layout.addLayout(transform_action_layout)

        transform_tab.setLayout(transform_tab_layout)

        refine_tolerance_label = QLabel("Tolerance:")

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0.01, 1, 5, notation=notation)

        self.refine_tolerance_line = QLineEdit("0.1")
        self.refine_tolerance_line.setValidator(validator)

        self.optimize_combo = QComboBox(self)
        self.optimize_combo.addItem("Unconstrained")
        self.optimize_combo.addItem("Constrained")
        self.optimize_combo.addItem("Triclinic")
        self.optimize_combo.addItem("Monoclinic")
        self.optimize_combo.addItem("Orthorhombic")
        self.optimize_combo.addItem("Tetragonal")
        self.optimize_combo.addItem("Rhombohedral")
        self.optimize_combo.addItem("Hexagonal")
        self.optimize_combo.addItem("Cubic")

        refine_tab = QWidget()
        refine_tab_layout = QVBoxLayout()

        refine_params_layout = QHBoxLayout()
        refine_params_layout.addWidget(refine_tolerance_label)
        refine_params_layout.addWidget(self.refine_tolerance_line)
        refine_params_layout.addWidget(self.optimize_combo)

        self.refine_button = QPushButton("Refine", self)

        refine_action_layout = QHBoxLayout()
        refine_action_layout.addWidget(self.refine_button)
        refine_action_layout.addStretch(1)

        refine_tab_layout.addLayout(refine_params_layout)
        refine_tab_layout.addStretch(1)
        refine_tab_layout.addLayout(refine_action_layout)

        refine_tab.setLayout(refine_tab_layout)

        ub_tab.addTab(calculate_tab, "Calculate UB")
        ub_tab.addTab(transform_tab, "Transform UB")
        ub_tab.addTab(refine_tab, "Refine UB")

        return ub_tab

    def table_tab(self):
        peaks_table_tab = QWidget()
        self.tab_widget.addTab(peaks_table_tab, "Peaks")

        peaks_layout = QVBoxLayout()

        calculator_layout = QGridLayout()

        h_label = QLabel("h", self)
        k_label = QLabel("k", self)
        l_label = QLabel("l", self)

        peak_1_label = QLabel("1:", self)
        peak_2_label = QLabel("2:", self)

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(-100, 100, 5, notation=notation)

        self.h1_line = QLineEdit()
        self.k1_line = QLineEdit()
        self.l1_line = QLineEdit()

        self.h2_line = QLineEdit()
        self.k2_line = QLineEdit()
        self.l2_line = QLineEdit()

        self.h1_line.setValidator(validator)
        self.k1_line.setValidator(validator)
        self.l1_line.setValidator(validator)

        self.h2_line.setValidator(validator)
        self.k2_line.setValidator(validator)
        self.l2_line.setValidator(validator)

        d_label = QLabel("d [Å]", self)

        phi_label = QLabel("φ [°]", self)

        self.d1_line = QLineEdit()
        self.d2_line = QLineEdit()
        self.phi_line = QLineEdit()

        self.d1_line.setEnabled(False)
        self.d2_line.setEnabled(False)
        self.phi_line.setEnabled(False)

        self.calculate = QPushButton("Calculate", self)

        calculator_layout.addWidget(h_label, 0, 1, Qt.AlignCenter)
        calculator_layout.addWidget(k_label, 0, 2, Qt.AlignCenter)
        calculator_layout.addWidget(l_label, 0, 3, Qt.AlignCenter)
        calculator_layout.addWidget(d_label, 0, 4, Qt.AlignCenter)
        calculator_layout.addWidget(phi_label, 0, 5, Qt.AlignCenter)

        calculator_layout.addWidget(peak_1_label, 1, 0)
        calculator_layout.addWidget(self.h1_line, 1, 1)
        calculator_layout.addWidget(self.k1_line, 1, 2)
        calculator_layout.addWidget(self.l1_line, 1, 3)
        calculator_layout.addWidget(self.d1_line, 1, 4)
        calculator_layout.addWidget(self.phi_line, 1, 5)

        calculator_layout.addWidget(peak_2_label, 2, 0)
        calculator_layout.addWidget(self.h2_line, 2, 1)
        calculator_layout.addWidget(self.k2_line, 2, 2)
        calculator_layout.addWidget(self.l2_line, 2, 3)
        calculator_layout.addWidget(self.d2_line, 2, 4)
        calculator_layout.addWidget(self.calculate, 2, 5)

        stretch = QHeaderView.Stretch

        self.peaks_table = QTableWidget()
        self.peaks_table.setRowCount(0)
        self.peaks_table.setColumnCount(7)

        header = ["h", "k", "l", "d", "λ", "I", "I/σ"]

        self.peaks_table.horizontalHeader().setSectionResizeMode(stretch)
        self.peaks_table.setHorizontalHeaderLabels(header)
        self.peaks_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.peaks_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.peaks_table.setSortingEnabled(True)

        extended_info = QGridLayout()

        d_label = QLabel("d [Å]:", self)
        lambda_label = QLabel("λ [Å]:", self)

        run_label = QLabel("Run #", self)
        bank_label = QLabel("Bank #", self)
        row_label = QLabel("Row #", self)
        col_label = QLabel("Col #", self)

        self.d_line = QLineEdit()
        self.lambda_line = QLineEdit()
        self.run_line = QLineEdit()
        self.bank_line = QLineEdit()
        self.row_line = QLineEdit()
        self.col_line = QLineEdit()

        self.d_line.setEnabled(False)
        self.lambda_line.setEnabled(False)
        self.run_line.setEnabled(False)
        self.bank_line.setEnabled(False)
        self.row_line.setEnabled(False)
        self.col_line.setEnabled(False)

        self.intensity_line = QLineEdit()
        self.sigma_line = QLineEdit()

        self.intensity_line.setEnabled(False)
        self.sigma_line.setEnabled(False)

        intensity_label = QLabel("I: ", self)
        sigma_label = QLabel("± σ:", self)

        extended_info.addWidget(intensity_label, 0, 0)
        extended_info.addWidget(self.intensity_line, 0, 1)
        extended_info.addWidget(sigma_label, 0, 2)
        extended_info.addWidget(self.sigma_line, 0, 3)

        extended_info.addWidget(d_label, 0, 4)
        extended_info.addWidget(self.d_line, 0, 5)
        extended_info.addWidget(lambda_label, 0, 6)
        extended_info.addWidget(self.lambda_line, 0, 7)

        extended_info.addWidget(run_label, 1, 0)
        extended_info.addWidget(self.run_line, 1, 1)
        extended_info.addWidget(bank_label, 1, 2)
        extended_info.addWidget(self.bank_line, 1, 3)
        extended_info.addWidget(row_label, 1, 4)
        extended_info.addWidget(self.row_line, 1, 5)
        extended_info.addWidget(col_label, 1, 6)
        extended_info.addWidget(self.col_line, 1, 7)

        hkl_info = QHBoxLayout()
        peak_info = QGridLayout()

        left_label = QLabel("(", self)
        left_comma_label = QLabel(",", self)
        right_comma_label = QLabel(",", self)
        right_label = QLabel(")", self)

        index_label = QLabel("Indexed:", self)
        total_label = QLabel("Total:", self)

        self.index_line = QLineEdit("0")
        self.total_line = QLineEdit("0")

        self.index_line.setEnabled(False)
        self.total_line.setEnabled(False)

        int_h_label = QLabel("h", self)
        int_k_label = QLabel("k", self)
        int_l_label = QLabel("l", self)

        int_m_label = QLabel("m", self)
        int_n_label = QLabel("n", self)
        int_p_label = QLabel("p", self)

        self.h_line = QLineEdit()
        self.k_line = QLineEdit()
        self.l_line = QLineEdit()

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(-100, 100, 5, notation=notation)

        self.h_line.setValidator(validator)
        self.k_line.setValidator(validator)
        self.l_line.setValidator(validator)

        validator = QIntValidator(-1000000000, 1000000000, self)

        self.int_h_line = QLineEdit()
        self.int_k_line = QLineEdit()
        self.int_l_line = QLineEdit()

        self.int_h_line.setValidator(validator)
        self.int_k_line.setValidator(validator)
        self.int_l_line.setValidator(validator)

        self.int_m_line = QLineEdit()
        self.int_n_line = QLineEdit()
        self.int_p_line = QLineEdit()

        self.int_m_line.setValidator(validator)
        self.int_n_line.setValidator(validator)
        self.int_p_line.setValidator(validator)

        hkl_info.addWidget(left_label)
        hkl_info.addWidget(self.h_line)
        hkl_info.addWidget(left_comma_label)
        hkl_info.addWidget(self.k_line)
        hkl_info.addWidget(right_comma_label)
        hkl_info.addWidget(self.l_line)
        hkl_info.addWidget(right_label)
        hkl_info.addStretch(1)
        hkl_info.addWidget(index_label)
        hkl_info.addWidget(self.index_line)
        hkl_info.addWidget(total_label)
        hkl_info.addWidget(self.total_line)

        peak_info.addWidget(int_h_label, 0, 0, Qt.AlignCenter)
        peak_info.addWidget(int_k_label, 0, 1, Qt.AlignCenter)
        peak_info.addWidget(int_l_label, 0, 2, Qt.AlignCenter)
        peak_info.addWidget(int_m_label, 0, 3, Qt.AlignCenter)
        peak_info.addWidget(int_n_label, 0, 4, Qt.AlignCenter)
        peak_info.addWidget(int_p_label, 0, 5, Qt.AlignCenter)

        peak_info.addWidget(self.int_h_line, 1, 0)
        peak_info.addWidget(self.int_k_line, 1, 1)
        peak_info.addWidget(self.int_l_line, 1, 2)
        peak_info.addWidget(self.int_m_line, 1, 3)
        peak_info.addWidget(self.int_n_line, 1, 4)
        peak_info.addWidget(self.int_p_line, 1, 5)

        peaks_layout.addLayout(calculator_layout)
        peaks_layout.addWidget(self.peaks_table)
        peaks_layout.addLayout(hkl_info)
        peaks_layout.addLayout(peak_info)
        peaks_layout.addLayout(extended_info)

        peaks_table_tab.setLayout(peaks_layout)

    def verify_tab(self):
        inspect_verify_tab = QTabWidget()
        self.tab_widget.addTab(inspect_verify_tab, "Views")

        inspect_tab = self.__init_inspect_tab()
        verify_tab = self.__init_verify_tab()

        inspect_verify_tab.addTab(inspect_tab, "Slice View")
        inspect_verify_tab.addTab(verify_tab, "Detector View")

    def __init_inspect_tab(self):
        convert_to_hkl_tab = QWidget()
        convert_to_hkl_tab_layout = QVBoxLayout()

        convert_to_hkl_params_layout = QGridLayout()

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(-10, 10, 5, notation=notation)

        self.U1_line = QLineEdit("1")
        self.U2_line = QLineEdit("0")
        self.U3_line = QLineEdit("0")

        self.V1_line = QLineEdit("0")
        self.V2_line = QLineEdit("1")
        self.V3_line = QLineEdit("0")

        self.W1_line = QLineEdit("0")
        self.W2_line = QLineEdit("0")
        self.W3_line = QLineEdit("1")

        self.U1_line.setValidator(validator)
        self.U2_line.setValidator(validator)
        self.U3_line.setValidator(validator)

        self.V1_line.setValidator(validator)
        self.V2_line.setValidator(validator)
        self.V3_line.setValidator(validator)

        self.W1_line.setValidator(validator)
        self.W2_line.setValidator(validator)
        self.W3_line.setValidator(validator)

        ax1_label = QLabel("1:")
        ax2_label = QLabel("2:")
        ax3_label = QLabel("3:")

        h_label = QLabel("h")
        k_label = QLabel("k")
        l_label = QLabel("l")

        convert_to_hkl_params_layout.addWidget(h_label, 0, 1, Qt.AlignCenter)
        convert_to_hkl_params_layout.addWidget(k_label, 0, 2, Qt.AlignCenter)
        convert_to_hkl_params_layout.addWidget(l_label, 0, 3, Qt.AlignCenter)
        convert_to_hkl_params_layout.addWidget(ax1_label, 1, 0, Qt.AlignCenter)
        convert_to_hkl_params_layout.addWidget(ax2_label, 2, 0, Qt.AlignCenter)
        convert_to_hkl_params_layout.addWidget(ax3_label, 3, 0, Qt.AlignCenter)

        convert_to_hkl_params_layout.addWidget(self.U1_line, 1, 1)
        convert_to_hkl_params_layout.addWidget(self.V1_line, 2, 1)
        convert_to_hkl_params_layout.addWidget(self.W1_line, 3, 1)

        convert_to_hkl_params_layout.addWidget(self.U2_line, 1, 2)
        convert_to_hkl_params_layout.addWidget(self.V2_line, 2, 2)
        convert_to_hkl_params_layout.addWidget(self.W2_line, 3, 2)

        convert_to_hkl_params_layout.addWidget(self.U3_line, 1, 3)
        convert_to_hkl_params_layout.addWidget(self.V3_line, 2, 3)
        convert_to_hkl_params_layout.addWidget(self.W3_line, 3, 3)

        self.convert_to_hkl_button = QPushButton("Convert", self)

        self.clim_combo = QComboBox(self)
        self.clim_combo.addItem("Min/Max")
        self.clim_combo.addItem("μ±3×σ")
        self.clim_combo.addItem("Q₃/Q₁±1.5×IQR")
        self.clim_combo.setCurrentIndex(1)

        self.cbar_combo = QComboBox(self)
        self.cbar_combo.addItem("Sequential")
        self.cbar_combo.addItem("Rainbow")
        self.cbar_combo.addItem("Binary")
        self.cbar_combo.addItem("Diverging")
        self.cbar_combo.setCurrentIndex(2)

        self.slice_combo = QComboBox(self)
        self.slice_combo.addItem("Axis 1/2")
        self.slice_combo.addItem("Axis 1/3")
        self.slice_combo.addItem("Axis 2/3")
        self.slice_combo.setCurrentIndex(0)

        slice_label = QLabel("Slice:", self)

        self.slice_line = QLineEdit("0.0")
        self.slice_line.setValidator(validator)

        validator = QDoubleValidator(0.0001, 100, 5, notation=notation)

        slice_thickness_label = QLabel("Thickness:", self)

        self.slice_thickness_line = QLineEdit("0.1")
        self.slice_thickness_line.setValidator(validator)

        validator = QDoubleValidator(0.01, 0.5, 5, notation=notation)

        slice_width_label = QLabel("Width:", self)

        self.slice_width_line = QLineEdit("0.05")
        self.slice_width_line.setValidator(validator)

        self.slice_scale_combo = QComboBox(self)
        self.slice_scale_combo.addItem("Linear")
        self.slice_scale_combo.addItem("Log")

        convert_to_hkl_action_layout = QHBoxLayout()
        convert_to_hkl_action_layout.addWidget(self.convert_to_hkl_button)
        convert_to_hkl_action_layout.addWidget(self.slice_combo)
        convert_to_hkl_action_layout.addWidget(slice_label)
        convert_to_hkl_action_layout.addWidget(self.slice_line)
        convert_to_hkl_action_layout.addWidget(slice_thickness_label)
        convert_to_hkl_action_layout.addWidget(self.slice_thickness_line)
        convert_to_hkl_action_layout.addWidget(slice_width_label)
        convert_to_hkl_action_layout.addWidget(self.slice_width_line)

        convert_to_hkl_view_layout = QHBoxLayout()
        convert_to_hkl_view_layout.addWidget(self.cbar_combo)
        convert_to_hkl_view_layout.addWidget(self.clim_combo)
        convert_to_hkl_view_layout.addWidget(self.slice_scale_combo)

        convert_to_hkl_tab_layout.addLayout(convert_to_hkl_params_layout)
        convert_to_hkl_tab_layout.addStretch(1)
        convert_to_hkl_tab_layout.addLayout(convert_to_hkl_action_layout)

        self.canvas_slice = FigureCanvas(Figure(figsize=[12.8, 12.8]))

        self.ax_xint = None
        self.ax_yint = None
        self.cb_slice = None
        self.cb_inst = None

        self.fig_slice = self.canvas_slice.figure

        self.ax_slice = self.fig_slice.subplots(1, 1)

        slice_layout = QVBoxLayout()

        slice_layout.addWidget(NavigationToolbar2QT(self.canvas_slice, self))
        slice_layout.addWidget(self.canvas_slice)

        convert_to_hkl_tab_layout.addLayout(slice_layout)
        convert_to_hkl_tab_layout.addLayout(convert_to_hkl_view_layout)

        convert_to_hkl_tab.setLayout(convert_to_hkl_tab_layout)

        return convert_to_hkl_tab

    def __init_verify_tab(self):
        instrument_tab = QWidget()
        instrument_tab_layout = QVBoxLayout()

        notation = QDoubleValidator.StandardNotation

        self.data_combo = QComboBox(self)

        d_min_label = QLabel("d(min):", self)
        d_max_label = QLabel("d(max):", self)

        validator = QDoubleValidator(0, float("inf"), 5, notation=notation)

        self.d_min_line = QLineEdit("0")
        self.d_min_line.setValidator(validator)

        self.d_max_line = QLineEdit("inf")
        self.d_max_line.setValidator(validator)

        data_layout = QHBoxLayout()
        data_layout.addWidget(self.data_combo)
        data_layout.addStretch(1)
        data_layout.addWidget(d_min_label)
        data_layout.addWidget(self.d_min_line)
        data_layout.addWidget(d_max_label)
        data_layout.addWidget(self.d_max_line)

        vertical_label = QLabel("Vertical Angle:", self)
        horizontal_label = QLabel("Horizontal Angle:", self)

        vertical_roi_label = QLabel("ROI:", self)
        horizontal_roi_label = QLabel("ROI:", self)

        validator = QDoubleValidator(-180, 180, 5, notation=notation)

        self.vertical_line = QLineEdit("0")
        self.vertical_line.setValidator(validator)

        self.horizontal_line = QLineEdit("0")
        self.horizontal_line.setValidator(validator)

        validator = QDoubleValidator(0, 180, 5, notation=notation)

        self.vertical_roi_line = QLineEdit("5")
        self.vertical_roi_line.setValidator(validator)

        self.horizontal_roi_line = QLineEdit("5")
        self.horizontal_roi_line.setValidator(validator)

        angle_layout = QHBoxLayout()
        angle_layout.addWidget(horizontal_label)
        angle_layout.addWidget(self.horizontal_line)
        angle_layout.addWidget(horizontal_roi_label)
        angle_layout.addWidget(self.horizontal_roi_line)
        angle_layout.addWidget(vertical_label)
        angle_layout.addWidget(self.vertical_line)
        angle_layout.addWidget(vertical_roi_label)
        angle_layout.addWidget(self.vertical_roi_line)

        self.add_peak_button = QPushButton("Add Peak", self)

        self.diffraction_label = QLabel("Axis:", self)

        validator = QDoubleValidator(
            -float("inf"), float("inf"), 5, notation=notation
        )

        self.diffraction_line = QLineEdit("0")
        self.diffraction_line.setValidator(validator)

        peak_layout = QHBoxLayout()
        peak_layout.addWidget(self.diffraction_label)
        peak_layout.addWidget(self.diffraction_line)
        peak_layout.addStretch(1)
        peak_layout.addWidget(self.add_peak_button)

        self.canvas_inst = FigureCanvas(Figure(constrained_layout=True))
        self.canvas_scan = FigureCanvas(Figure(constrained_layout=True))

        self.fig_inst = self.canvas_inst.figure
        self.fig_scan = self.canvas_scan.figure

        self.ax_inst = self.fig_inst.subplots(1, 1)
        self.ax_scan = self.fig_scan.subplots(1, 1)

        view_layout = QVBoxLayout()

        view_layout.addLayout(data_layout)
        view_layout.addWidget(NavigationToolbar2QT(self.canvas_inst, self))
        view_layout.addWidget(self.canvas_inst)

        view_layout.addLayout(angle_layout)
        view_layout.addWidget(NavigationToolbar2QT(self.canvas_scan, self))
        view_layout.addWidget(self.canvas_scan)

        view_layout.addLayout(peak_layout)

        instrument_tab_layout.addLayout(view_layout)

        instrument_tab.setLayout(instrument_tab_layout)

        return instrument_tab

    def connect_h_index(self, update_index):
        self.h_line.editingFinished.connect(update_index)

    def connect_k_index(self, update_index):
        self.k_line.editingFinished.connect(update_index)

    def connect_l_index(self, update_index):
        self.l_line.editingFinished.connect(update_index)

    def connect_integer_h_index(self, update_index):
        self.int_h_line.editingFinished.connect(update_index)

    def connect_integer_k_index(self, update_index):
        self.int_k_line.editingFinished.connect(update_index)

    def connect_integer_l_index(self, update_index):
        self.int_l_line.editingFinished.connect(update_index)

    def connect_integer_m_index(self, update_index):
        self.int_m_line.editingFinished.connect(update_index)

    def connect_integer_n_index(self, update_index):
        self.int_n_line.editingFinished.connect(update_index)

    def connect_integer_p_index(self, update_index):
        self.int_p_line.editingFinished.connect(update_index)

    def connect_data_combo(self, update_inst_data):
        self.data_combo.currentIndexChanged.connect(update_inst_data)

    def connect_add_peak(self, add_peak):
        self.add_peak_button.clicked.connect(add_peak)

    def connect_diffraction(self, update_inst_data):
        self.diffraction_line.editingFinished.connect(update_inst_data)

    def connect_d_min(self, update_inst_data):
        self.d_min_line.editingFinished.connect(update_inst_data)

    def connect_d_max(self, update_inst_data):
        self.d_max_line.editingFinished.connect(update_inst_data)

    def connect_horizontal(self, update_inst_data):
        self.horizontal_line.editingFinished.connect(update_inst_data)

    def connect_vertical(self, update_inst_data):
        self.vertical_line.editingFinished.connect(update_inst_data)

    def connect_horizontal_roi(self, update_inst_data):
        self.horizontal_roi_line.editingFinished.connect(update_inst_data)

    def connect_vertical_roi(self, update_inst_data):
        self.vertical_roi_line.editingFinished.connect(update_inst_data)

    def connect_convert_to_hkl(self, convert_to_hkl):
        self.convert_to_hkl_button.clicked.connect(convert_to_hkl)

    def connect_browse_calibration(self, load_detector_cal):
        self.cal_browse_button.clicked.connect(load_detector_cal)

    def connect_browse_tube(self, load_tube_cal):
        self.tube_browse_button.clicked.connect(load_tube_cal)

    def connect_convert_Q(self, convert_Q):
        self.convert_to_q_button.clicked.connect(convert_Q)

    def connect_find_peaks(self, find_peaks):
        self.find_button.clicked.connect(find_peaks)

    def connect_find_conventional(self, find_conventional):
        self.conventional_button.clicked.connect(find_conventional)

    def connect_find_niggli(self, find_niggli):
        self.niggli_button.clicked.connect(find_niggli)

    def connect_select_form(self, select_form):
        self.select_button.clicked.connect(select_form)

    def connect_convert_HKL(self, convert_HKL):
        self.convert_to_hkl_button.clicked.connect(convert_HKL)

    def connect_switch_instrument(self, switch_instrument):
        self.instrument_combo.activated.connect(switch_instrument)

    def connect_wavelength(self, update_wavelength):
        self.wl_min_line.editingFinished.connect(update_wavelength)

    def connect_load_Q(self, load_Q):
        self.load_q_button.clicked.connect(load_Q)

    def connect_save_Q(self, save_Q):
        self.save_q_button.clicked.connect(save_Q)

    def connect_load_peaks(self, load_peaks):
        self.load_peaks_button.clicked.connect(load_peaks)

    def connect_save_peaks(self, save_peaks):
        self.save_peaks_button.clicked.connect(save_peaks)

    def connect_load_UB(self, load_UB):
        self.load_ub_button.clicked.connect(load_UB)

    def connect_save_UB(self, save_UB):
        self.save_ub_button.clicked.connect(save_UB)

    def connect_lattice_transform(self, lattice_transform):
        self.lattice_combo.currentIndexChanged.connect(lattice_transform)

    def connect_symmetry_transform(self, symmetry_transform):
        self.symmetry_combo.currentIndexChanged.connect(symmetry_transform)

    def connect_transform_UB(self, transform_UB):
        self.transform_button.clicked.connect(transform_UB)

    def connect_optimize_UB(self, optimize_UB):
        self.refine_button.clicked.connect(optimize_UB)

    def connect_index_peaks(self, index_peaks):
        self.index_button.clicked.connect(index_peaks)

    def connect_predict_peaks(self, predict_peaks):
        self.predict_button.clicked.connect(predict_peaks)

    def connect_filter_peaks(self, filter_peaks):
        self.filter_button.clicked.connect(filter_peaks)

    def connect_integrate_peaks(self, integrate_peaks):
        self.integrate_button.clicked.connect(integrate_peaks)

    def connect_calculate_peaks(self, calculate_peaks):
        self.calculate.clicked.connect(calculate_peaks)

    def connect_peak_row_highligter(self, highlight_row):
        self.peaks_table.itemSelectionChanged.connect(highlight_row)

    def connect_cell_row_highligter(self, highlight_row):
        self.cell_table.itemSelectionChanged.connect(highlight_row)

    def connect_select_cell(self, select_cell):
        self.select_button.clicked.connect(select_cell)

    def load_detector_cal_dialog(self, path=""):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        file_filters = "Calibration files (*.DetCal *.detcal *.xml)"

        filename, _ = file_dialog.getOpenFileName(
            self, "Load calibration file", path, file_filters, options=options
        )

        return filename

    def load_tube_cal_dialog(self, path=""):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        file_filters = "Tube files (*.h5 *.nxs)"

        filename, _ = file_dialog.getOpenFileName(
            self, "Load calibration file", path, file_filters, options=options
        )

        return filename

    def load_Q_file_dialog(self, path=""):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filename, _ = file_dialog.getOpenFileName(
            self, "Load Q file", path, "Q files (*.nxs)", options=options
        )

        return filename

    def save_Q_file_dialog(self, path=""):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filename, _ = file_dialog.getSaveFileName(
            self, "Save Q file", path, "Q files (*.nxs)", options=options
        )

        if filename is not None:
            if not filename.endswith(".nxs"):
                filename += ".nxs"

        return filename

    def load_peaks_file_dialog(self, path=""):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filename, _ = QFileDialog.getOpenFileName(
            self,
            "Load peaks file",
            path,
            "Peaks files (*.nxs)",
            options=options,
        )

        return filename

    def save_peaks_file_dialog(self, path=""):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filename, _ = file_dialog.getSaveFileName(
            self,
            "Save peaks file",
            path,
            "Peaks files (*.nxs)",
            options=options,
        )

        if filename is not None:
            if not filename.endswith(".nxs"):
                filename += ".nxs"

        return filename

    def load_UB_file_dialog(self, path=""):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filename, _ = file_dialog.getOpenFileName(
            self, "Load UB file", path, "UB files (*.mat)", options=options
        )

        return filename

    def save_UB_file_dialog(self, path=""):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filename, _ = file_dialog.getSaveFileName(
            self, "Save UB file", path, "UB files (*.mat)", options=options
        )

        if filename is not None:
            if not filename.endswith(".mat"):
                filename += ".mat"

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

    def get_instrument(self):
        return self.instrument_combo.currentText()

    def update_diffraction_label(self, mono):
        text = "Wavelength:" if not mono else "Angle:"

        self.diffraction_label.setText(text)

    def clear_run_info(self, filepath):
        self.exp_line.setText("")
        # self.run_line.setText('')
        # self.ipts_line.setText('')
        self.cal_line.setText("")
        self.tube_line.setText("")

        if "exp" in filepath:
            self.exp_line.setEnabled(True)
        else:
            self.exp_line.setEnabled(False)

        if "SNS" in filepath:
            self.filter_time_line.setEnabled(True)
            self.filter_time_line.setText("300")
            self.tube_line.setEnabled(False)
            self.tube_browse_button.setEnabled(False)
            if "CORELLI" in filepath:
                self.tube_line.setEnabled(True)
                self.tube_browse_button.setEnabled(True)
        else:
            self.filter_time_line.setEnabled(False)
            self.filter_time_line.setText("")
            self.cal_line.setEnabled(False)
            self.cal_browse_button.setEnabled(False)
            self.tube_line.setEnabled(False)
            self.tube_browse_button.setEnabled(False)

    def get_tube_calibration(self):
        return self.tube_line.text()

    def get_detector_calibration(self):
        return self.cal_line.text()

    def set_tube_calibration(self, filename):
        return self.tube_line.setText(filename)

    def set_detector_calibration(self, filename):
        return self.cal_line.setText(filename)

    def get_IPTS(self):
        if self.ipts_line.hasAcceptableInput():
            return self.ipts_line.text()

    def get_experiment(self):
        if self.exp_line.hasAcceptableInput():
            return self.exp_line.text()

    def runs_string_to_list(self, runs_str):
        """
        Convert runs string to list using regex validation.
        Return None for invalid formats.

        Parameters
        ----------
        runs_str : str
            Condensed notation for run numbers.

        Returns
        -------
        runs : list or None
            Integer run numbers or None if the input is invalid.

        """

        pattern = r"^(\d+(:\d+)?)(,\d+(:\d+)?)*$"

        if not re.match(pattern, runs_str):
            return None

        ranges = runs_str.split(",")
        runs = []
        for part in ranges:
            if ":" in part:
                start, end = map(int, part.split(":"))
                if start > end:
                    return None
                runs.extend(range(start, end + 1))
            else:
                runs.append(int(part))
        return runs

    def get_runs(self):
        run_str = self.runs_line.text()

        return self.runs_string_to_list(run_str)

    def get_lorentz(self):
        return self.lorentz_box.isChecked()

    def get_time_stop(self):
        if self.filter_time_line.hasAcceptableInput():
            return self.filter_time_line.text()

    def add_Q_viz(self, Q_dict):
        self.clear_scene()

        signal = Q_dict.get("signal")
        x = Q_dict.get("x")
        y = Q_dict.get("y")
        z = Q_dict.get("z")

        if all([elem is not None for elem in [signal, x, y, z]]):
            points = np.column_stack([x, y, z])

            point_cloud = pv.PolyData(points)
            point_cloud["scalars"] = signal

            self.plotter.add_mesh(
                point_cloud,
                scalars="scalars",
                cmap="binary",
                show_scalar_bar=False,
                opacity=1,
                log_scale=True,
                point_size=1,
                smooth_shading=False,
                culling=False,
                emissive=False,
                style="points_gaussian",
            )

        transforms = Q_dict.get("transforms")
        intensities = Q_dict.get("intensities")
        indexings = Q_dict.get("indexings")
        numbers = Q_dict.get("numbers")

        params = [transforms, intensities, indexings, numbers]

        integrate = np.any(intensities)

        if all([elem is not None for elem in params]) and len(numbers) > 0:
            sphere = pv.Icosphere(radius=1, nsub=1)

            geoms, self.indexing = [], {}
            for i, (T, I, ind, no) in enumerate(zip(*params)):
                ellipsoid = sphere.copy().transform(T)
                color = I if integrate else ind
                ellipsoid["scalars"] = np.full(sphere.n_cells, color)
                geoms.append(ellipsoid)
                self.indexing[i] = i

            multiblock = pv.MultiBlock(geoms)

            mu = np.nanmean(intensities)
            sigma = np.nanstd(intensities)

            cmap = "viridis" if integrate else ["lightblue", "lightgreen"]
            n_colors = 256 if integrate else 2
            clim = [mu - 3 * sigma, mu + 3 * sigma] if integrate else [0, 1]

            _, mapper = self.plotter.add_composite(
                multiblock,
                scalars="scalars",
                color=None,
                log_scale=False,
                style="surface",
                cmap=cmap,
                clim=clim,
                n_colors=n_colors,
                show_scalar_bar=False,
                smooth_shading=True,
            )

            self.mapper = mapper

            self.plotter.enable_block_picking(
                callback=self.highlight, side="left"
            )
            self.plotter.enable_block_picking(
                callback=self.highlight, side="right"
            )

            self.last_highlight = None

        self.reset_scene()

    def highlight(self, index, dataset):
        if self.last_highlight is not None:
            self.mapper.block_attr[self.last_highlight].color = None

        color = self.mapper.block_attr[index].color

        self.peaks_table.clearSelection()

        if color == "pink":
            color = None
        else:
            color = "pink"
            self.last_highlight = index

        self.mapper.block_attr[index].color = color

        ind = self.indexing[index - 1]

        if color == "pink":
            selected = self.peaks_table.selectedIndexes()
            if selected:
                selected_row = selected[0].row()
                if selected_row == ind:
                    return
            self.peaks_table.selectRow(ind)

    def highlight_peak(self, index):
        if self.last_highlight is not None:
            self.mapper.block_attr[self.last_highlight].color = None

        self.mapper.block_attr[index].color = "pink"
        self.last_highlight = index

    def set_sample_directions(self, params):
        v, w, u = params

        self.uh_line.setText("{}".format(u[0]))
        self.uk_line.setText("{}".format(u[1]))
        self.ul_line.setText("{}".format(u[2]))

        self.vh_line.setText("{}".format(v[0]))
        self.vk_line.setText("{}".format(v[1]))
        self.vl_line.setText("{}".format(v[2]))

        self.wh_line.setText("{}".format(w[0]))
        self.wk_line.setText("{}".format(w[1]))
        self.wl_line.setText("{}".format(w[2]))

    def set_lattice_constants(self, params):
        self.a_line.setText("{:.4f}".format(params[0]))
        self.b_line.setText("{:.4f}".format(params[1]))
        self.c_line.setText("{:.4f}".format(params[2]))

        self.alpha_line.setText("{:.4f}".format(params[3]))
        self.beta_line.setText("{:.4f}".format(params[4]))
        self.gamma_line.setText("{:.4f}".format(params[5]))

    def get_lattice_constants(self):
        params = (
            self.a_line,
            self.b_line,
            self.c_line,
            self.alpha_line,
            self.beta_line,
            self.gamma_line,
        )

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:
            return [float(param.text()) for param in params]

    def get_min_max_constants(self):
        params = self.min_const_line, self.max_const_line

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:
            return [float(param.text()) for param in params]

    def get_find_peaks_parameters(self):
        params = self.density_threshold_line, self.max_peaks_line

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:
            return [int(param.text()) for param in params]

    def get_find_peaks_distance(self):
        param = self.min_distance_line

        if param.hasAcceptableInput():
            return float(param.text())

    def get_find_peaks_edge(self):
        param = self.find_edge_line

        if param.hasAcceptableInput():
            return int(param.text())

    def get_calculate_UB_tol(self):
        param = self.calculate_tolerance_line

        if param.hasAcceptableInput():
            return float(param.text())

    def get_lattice_transform(self):
        return self.lattice_combo.currentText()

    def get_symmetry_symbol(self):
        return self.symmetry_combo.currentText()

    def update_symmetry_symbols(self, symbols):
        self.symmetry_combo.clear()
        for symbol in symbols:
            self.symmetry_combo.addItem(symbol)

    def get_transform_matrix(self):
        params = (
            self.T11_line,
            self.T12_line,
            self.T13_line,
            self.T21_line,
            self.T22_line,
            self.T23_line,
            self.T31_line,
            self.T32_line,
            self.T33_line,
        )
        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:
            params = [float(param.text()) for param in params]

            return params

    def get_projection_matrix(self):
        params = (
            self.U1_line,
            self.U2_line,
            self.U3_line,
            self.V1_line,
            self.V2_line,
            self.V3_line,
            self.W1_line,
            self.W2_line,
            self.W3_line,
        )
        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:
            params = [float(param.text()) for param in params]

            return params

    def set_transform_matrix(self, params):
        self.T11_line.setText("{:.0f}".format(params[0][0]))
        self.T12_line.setText("{:.0f}".format(params[0][1]))
        self.T13_line.setText("{:.0f}".format(params[0][2]))

        self.T21_line.setText("{:.0f}".format(params[1][0]))
        self.T22_line.setText("{:.0f}".format(params[1][1]))
        self.T23_line.setText("{:.0f}".format(params[1][2]))

        self.T31_line.setText("{:.0f}".format(params[2][0]))
        self.T32_line.setText("{:.0f}".format(params[2][1]))
        self.T33_line.setText("{:.0f}".format(params[2][2]))

    def get_transform_UB_tol(self):
        param = self.transform_tolerance_line

        if param.hasAcceptableInput():
            return float(param.text())

    def get_refine_UB_option(self):
        return self.optimize_combo.currentText()

    def get_refine_UB_tol(self):
        param = self.refine_tolerance_line

        if param.hasAcceptableInput():
            return float(param.text())

    def get_modulatation_offsets(self):
        params = (
            self.dh1_line,
            self.dk1_line,
            self.dl1_line,
            self.dh2_line,
            self.dk2_line,
            self.dl2_line,
            self.dh3_line,
            self.dk3_line,
            self.dl3_line,
        )

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:
            params = [float(param.text()) for param in params]

            return params

    def set_modulatation_offsets(self, params):
        self.dh1_line.setText("{:.3f}".format(params[0][0]))
        self.dk1_line.setText("{:.3f}".format(params[0][1]))
        self.dl1_line.setText("{:.3f}".format(params[0][2]))

        self.dh2_line.setText("{:.3f}".format(params[1][0]))
        self.dk2_line.setText("{:.3f}".format(params[1][1]))
        self.dl2_line.setText("{:.3f}".format(params[1][2]))

        self.dh3_line.setText("{:.3f}".format(params[2][0]))
        self.dk3_line.setText("{:.3f}".format(params[2][1]))
        self.dl3_line.setText("{:.3f}".format(params[2][2]))

    def get_max_order_cross_terms(self):
        param = self.max_order_line

        if param.hasAcceptableInput():
            return int(param.text()), self.cross_box.isChecked()

    def get_max_scalar_error(self):
        param = self.max_scalar_error_line

        if param.hasAcceptableInput():
            return float(param.text())

    def get_form_number(self):
        form = self.form_line.text()

        if form != "":
            return int(form)

    def get_index_peaks_parameters(self):
        param = self.index_tolerance_line

        sat_param = self.index_sat_tolerance_line

        if param.hasAcceptableInput():
            tol = float(param.text())

            if sat_param.hasAcceptableInput():
                sat_tol = float(sat_param.text())
            else:
                sat_tol = tol

            return tol, sat_tol

    def get_index_satellite_peaks(self):
        return self.index_sat_box.isChecked()

    def get_index_peaks_round(self):
        return self.round_box.isChecked()

    def get_predict_peaks_centering(self):
        return self.centering_combo.currentText()

    def get_predict_peaks_parameters(self):
        param = self.min_d_line

        sat_param = self.min_sat_d_line

        if param.hasAcceptableInput():
            d_min = float(param.text())

            if sat_param.hasAcceptableInput():
                sat_d_min = float(sat_param.text())
            else:
                sat_d_min = d_min

            return d_min, sat_d_min

    def get_predict_peaks_edge(self):
        param = self.predict_edge_line

        if param.hasAcceptableInput():
            return int(param.text())

    def get_predict_satellite_peaks(self):
        return self.predict_sat_box.isChecked()

    def get_integrate_peaks_parameters(self):
        params = self.radius_line, self.inner_line, self.outer_line

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:
            params = [float(param.text()) for param in params]

            return params

    def get_centroid(self):
        return self.centroid_box.isChecked()

    def get_ellipsoid(self):
        return self.adaptive_box.isChecked()

    def get_filter_variable(self):
        return self.filter_combo.currentText()

    def get_filter_comparison(self):
        return self.comparison_combo.currentText()

    def get_filter_value(self):
        param = self.filter_line

        if param.hasAcceptableInput():
            return float(param.text())

    def update_peaks_table(self, peaks):
        self.peaks_table.clearSelection()
        self.peaks_table.setRowCount(0)
        self.peaks_table.setRowCount(len(peaks))

        ind, tot = 0, 0
        for row, peak in enumerate(peaks):
            self.set_peak(row, peak)
            ind += peak[-2]
            tot += 1

        self.index_line.setText("{}".format(ind))
        self.total_line.setText("{}".format(tot))

    def set_peak(self, row, peak):
        hkl, d, lamda, intens, sig_noise, *_ = peak
        h, k, l = hkl
        h = "{:.3f}".format(h)
        k = "{:.3f}".format(k)
        l = "{:.3f}".format(l)
        d = "{:.4f}".format(d)
        lamda = "{:.4f}".format(lamda)
        intens = "{:.2e}".format(intens)
        sig_noise = "{:.2f}".format(sig_noise)
        self.peaks_table.setItem(row, 0, QTableWidgetItem(h))
        self.peaks_table.setItem(row, 1, QTableWidgetItem(k))
        self.peaks_table.setItem(row, 2, QTableWidgetItem(l))
        self.peaks_table.setItem(row, 3, QTableWidgetItem(d))
        self.peaks_table.setItem(row, 4, QTableWidgetItem(lamda))
        self.peaks_table.setItem(row, 5, QTableWidgetItem(intens))
        self.peaks_table.setItem(row, 6, QTableWidgetItem(sig_noise))

    def clear_niggli_info(self):
        self.cell_table.clearSelection()
        self.cell_table.setRowCount(0)
        self.form_line.setText("")

    def update_cell_table(self, cells):
        self.cell_table.clearSelection()
        self.cell_table.setRowCount(0)
        self.cell_table.setRowCount(len(cells))

        for row, cell in enumerate(cells):
            self.set_cell(row, cell)

    def set_cell(self, row, cell):
        form, error, bl, params = cell
        a, b, c, alpha, beta, gamma, vol = params
        error = "{:.4f}".format(error)
        bravais = " ".join(bl)
        a = "{:.2f}".format(a)
        b = "{:.2f}".format(b)
        c = "{:.2f}".format(c)
        alpha = "{:.1f}".format(alpha)
        beta = "{:.1f}".format(beta)
        gamma = "{:.1f}".format(gamma)
        vol = "{:.0f}".format(vol)
        self.cell_table.setVerticalHeaderItem(row, QTableWidgetItem(str(form)))
        self.cell_table.setItem(row, 0, QTableWidgetItem(error))
        self.cell_table.setItem(row, 1, QTableWidgetItem(bravais))
        self.cell_table.setItem(row, 2, QTableWidgetItem(a))
        self.cell_table.setItem(row, 3, QTableWidgetItem(b))
        self.cell_table.setItem(row, 4, QTableWidgetItem(c))
        self.cell_table.setItem(row, 5, QTableWidgetItem(alpha))
        self.cell_table.setItem(row, 6, QTableWidgetItem(beta))
        self.cell_table.setItem(row, 7, QTableWidgetItem(gamma))
        self.cell_table.setItem(row, 8, QTableWidgetItem(vol))

    def get_form(self):
        row = self.cell_table.currentRow()
        if row is not None:
            item = int(self.cell_table.verticalHeaderItem(row).text())
            return item

    def set_cell_form(self, form):
        self.form_line.setText(str(form))

    def get_peak(self):
        row = self.peaks_table.currentRow()
        if row is not None:
            return row

    def set_peak_info(self, peak):
        (
            hkl,
            d,
            lamda,
            intens,
            signal_noise,
            sigma,
            int_hkl,
            int_mnp,
            run,
            bank,
            row,
            col,
            ind,
            Q,
        ) = peak

        self.set_indices(hkl, int_hkl, int_mnp)

        self.intensity_line.setText("{:.2e}".format(intens))
        self.sigma_line.setText("{:.2e}".format(sigma))

        self.lambda_line.setText("{:.4f}".format(lamda))
        self.d_line.setText("{:.4f}".format(d))

        self.run_line.setText(str(run))
        self.bank_line.setText(str(bank))
        self.row_line.setText(str(row))
        self.col_line.setText(str(col))

    def update_table_index(self, row, hkl):
        h, k, l = hkl
        h = "{:.3f}".format(h)
        k = "{:.3f}".format(k)
        l = "{:.3f}".format(l)
        self.peaks_table.setItem(row, 0, QTableWidgetItem(h))
        self.peaks_table.setItem(row, 1, QTableWidgetItem(k))
        self.peaks_table.setItem(row, 2, QTableWidgetItem(l))

    def set_indices(self, hkl, int_hkl, int_mnp):
        H, K, L = hkl

        h, k, l = int_hkl
        m, n, p = int_mnp

        self.h_line.blockSignals(True)
        self.k_line.blockSignals(True)
        self.l_line.blockSignals(True)

        self.int_h_line.blockSignals(True)
        self.int_k_line.blockSignals(True)
        self.int_l_line.blockSignals(True)

        self.int_m_line.blockSignals(True)
        self.int_n_line.blockSignals(True)
        self.int_p_line.blockSignals(True)

        self.h_line.setText("{:.3f}".format(H))
        self.k_line.setText("{:.3f}".format(K))
        self.l_line.setText("{:.3f}".format(L))

        self.int_h_line.setText("{:.0f}".format(h))
        self.int_k_line.setText("{:.0f}".format(k))
        self.int_l_line.setText("{:.0f}".format(l))

        self.int_m_line.setText("{:.0f}".format(m))
        self.int_n_line.setText("{:.0f}".format(n))
        self.int_p_line.setText("{:.0f}".format(p))

        self.h_line.blockSignals(False)
        self.k_line.blockSignals(False)
        self.l_line.blockSignals(False)

        self.int_h_line.blockSignals(False)
        self.int_k_line.blockSignals(False)
        self.int_l_line.blockSignals(False)

        self.int_m_line.blockSignals(False)
        self.int_n_line.blockSignals(False)
        self.int_p_line.blockSignals(False)

    def get_indices(self):
        params_hkl = self.h_line, self.k_line, self.l_line
        params_int_hkl = self.int_h_line, self.int_k_line, self.int_l_line
        params_int_mnp = self.int_m_line, self.int_n_line, self.int_p_line

        params = params_hkl + params_int_hkl + params_int_mnp

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:
            hkl = [float(param.text()) for param in params_hkl]
            int_hkl = [int(param.text()) for param in params_int_hkl]
            int_mnp = [int(param.text()) for param in params_int_mnp]
            return hkl, int_hkl, int_mnp

    def connect_hand_index_peak(self, reindex):
        self.index_ready.connect(reindex)

    def get_input_hkls(self):
        params_1 = self.h1_line, self.k1_line, self.l1_line
        params_2 = self.h2_line, self.k2_line, self.l2_line

        valid_params = all([param.hasAcceptableInput() for param in params_1])

        if valid_params:
            params_1 = [float(param.text()) for param in params_1]
        else:
            params_1 = None

        valid_params = all([param.hasAcceptableInput() for param in params_2])

        if valid_params:
            params_2 = [float(param.text()) for param in params_2]
        else:
            params_2 = None

        return params_1, params_2

    def set_d_phi(self, d_1, d_2, phi_12):
        if d_1 is not None:
            self.d1_line.setText("{:.4f}".format(d_1))
        else:
            self.d1_line.setText("")
        if d_2 is not None:
            self.d2_line.setText("{:.4f}".format(d_2))
        else:
            self.d2_line.setText("")
        if phi_12 is not None:
            self.phi_line.setText("{:.4f}".format(phi_12))
        else:
            self.phi_line.setText("")

    def get_data_combo(self):
        return self.data_combo.currentIndex()

    def get_diffraction(self):
        if self.diffraction_line.hasAcceptableInput():
            return float(self.diffraction_line.text())

    def set_diffraction(self, val):
        self.diffraction_line.setText(str(round(val, 3)))

    def get_d_min(self):
        if self.d_min_line.hasAcceptableInput():
            return float(self.d_min_line.text())

    def get_d_max(self):
        text = self.d_max_line.text()

        if self.d_max_line.hasAcceptableInput() or text == "inf":
            return float(text)

    def get_horizontal(self):
        if self.horizontal_line.hasAcceptableInput():
            return float(self.horizontal_line.text())

    def get_vertical(self):
        if self.vertical_line.hasAcceptableInput():
            return float(self.vertical_line.text())

    def set_horizontal(self, val):
        self.horizontal_line.setText(str(round(val, 2)))

    def set_vertical(self, val):
        self.vertical_line.setText(str(round(val, 2)))

    def get_horizontal_roi(self):
        if self.horizontal_roi_line.hasAcceptableInput():
            return float(self.horizontal_roi_line.text())

    def get_vertical_roi(
        self,
    ):
        if self.vertical_roi_line.hasAcceptableInput():
            return float(self.vertical_roi_line.text())

    def update_instrument_view(self, inst_view, norm="log"):
        gamma = inst_view["gamma"]
        nu = inst_view["nu"]
        counts = inst_view["counts"]

        if self.cb_inst is not None:
            self.cb_inst.remove()
            self.cb_inst = None

        self.ax_inst.clear()

        self.im = self.ax_inst.scatter(
            gamma, nu, c=counts, marker="o", norm=norm, rasterized=True
        )

        self.ax_inst.set_aspect(1)
        self.ax_inst.minorticks_on()

        self.ax_inst.set_xlabel(r"$\gamma$")
        self.ax_inst.set_ylabel(r"$\nu$")

        fmt_str_form = FormatStrFormatter(r"$%d^\circ$")

        self.ax_inst.xaxis.set_major_formatter(fmt_str_form)
        self.ax_inst.yaxis.set_major_formatter(fmt_str_form)

        # self.cb_inst = self.fig_inst.colorbar(self.im, ax=self.ax_inst)
        # self.cb_inst.minorticks_on()

        self.canvas_inst.draw_idle()
        self.canvas_inst.flush_events()

    def update_roi_view(self, roi_view):
        horz = roi_view["horz"]
        vert = roi_view["vert"]
        horz_roi = roi_view["horz_roi"]
        vert_roi = roi_view["vert_roi"]

        for line in self.ax_inst.lines:
            line.remove()

        self.ax_inst.axvline(x=horz - horz_roi, color="k", linestyle="--")
        self.ax_inst.axvline(x=horz + horz_roi, color="k", linestyle="--")

        self.ax_inst.axhline(y=vert - vert_roi, color="k", linestyle="--")
        self.ax_inst.axhline(y=vert + vert_roi, color="k", linestyle="--")

        self.canvas_inst.draw_idle()
        self.canvas_inst.flush_events()

        self.inst_roi = {"roi": (horz_roi, vert_roi)}

        self.fig_inst.canvas.mpl_connect(
            "button_press_event", self.on_press_inst
        )

    def update_scan_view(self, roi_view):
        x = roi_view["x"]
        y = roi_view["y"]
        val = roi_view["val"]
        label = roi_view["label"]

        self.ax_scan.clear()

        self.ax_scan.errorbar(x, y, yerr=np.sqrt(y), fmt="o", color="C0")
        self.ax_scan.plot(x, y, color="C1")
        # self.ax_scan.set_yscale('log')
        self.line_scan = self.ax_scan.axvline(x=val, color="k", linestyle="--")
        self.ax_scan.minorticks_on()

        if label == "wavelength":
            xlabel = r"$\lambda$ [Å]"
        else:
            xlabel = r"$\vartheta$ [°]"

        self.ax_scan.set_xlabel(xlabel)

        self.canvas_scan.draw_idle()
        self.canvas_scan.flush_events()

        self.fig_scan.canvas.mpl_connect(
            "button_press_event", self.on_press_scan
        )

    def on_press_scan(self, event):
        if (
            event.inaxes == self.ax_scan
            and self.fig_scan.canvas.toolbar.mode == ""
        ):
            val = event.xdata

            self.diffraction_line.blockSignals(True)

            self.set_diffraction(val)

            self.diffraction_line.blockSignals(False)

            self.line_scan.set_xdata([val])

            self.canvas_scan.draw_idle()
            self.canvas_scan.flush_events()

    def on_press_inst(self, event):
        if (
            event.inaxes == self.ax_inst
            and self.fig_inst.canvas.toolbar.mode == ""
        ):
            for line in self.ax_inst.lines:
                line.remove()

            horz_roi, vert_roi = self.inst_roi["roi"]

            horz, vert = event.xdata, event.ydata

            self.horizontal_line.blockSignals(True)
            self.vertical_line.blockSignals(True)

            self.set_horizontal(horz)
            self.set_vertical(vert)

            self.horizontal_line.blockSignals(False)
            self.vertical_line.blockSignals(False)

            self.ax_inst.axvline(x=horz - horz_roi, color="k", linestyle="--")
            self.ax_inst.axvline(x=horz + horz_roi, color="k", linestyle="--")

            self.ax_inst.axhline(y=vert - vert_roi, color="k", linestyle="--")
            self.ax_inst.axhline(y=vert + vert_roi, color="k", linestyle="--")

            self.canvas_inst.draw_idle()
            self.canvas_inst.flush_events()

            self.roi_ready.emit()

    def connect_roi_ready(self, replot):
        self.roi_ready.connect(replot)

    def get_slice_value(self):
        if self.slice_line.hasAcceptableInput():
            return float(self.slice_line.text())

    def get_slice_thickness(self):
        if self.slice_thickness_line.hasAcceptableInput():
            return float(self.slice_thickness_line.text())

    def get_slice_width(self):
        if self.slice_width_line.hasAcceptableInput():
            return float(self.slice_width_line.text())

    def get_clim_clip_type(self):
        return self.clim_combo.currentText()

    def get_slice(self):
        return self.slice_combo.currentText()

    def get_slice_scale(self):
        return self.slice_scale_combo.currentText().lower()

    def get_colormap(self):
        return self.cbar_combo.currentText()

    def __format_axis_coord(self, x, y):
        x, y, _ = np.dot(self.T_inv, [x, y, 1])
        return "x={:.3f}, y={:.3f}".format(x, y)

    def update_slice(self, slice_dict):
        cmap = cmaps[self.get_colormap()]

        x = slice_dict["x"]
        y = slice_dict["y"]

        labels = slice_dict["labels"]
        title = slice_dict["title"]
        signal = slice_dict["signal"]
        clip = slice_dict["clip"]

        scale = self.get_slice_scale()

        vmin = np.nanmin(clip)
        vmax = np.nanmax(clip)

        if scale == "log" and np.isclose(vmin, 0):
            vmin = np.nanmin(signal[signal > 0])

        if np.isclose(vmax, vmin) or not np.isfinite([vmin, vmax]).all():
            vmin, vmax = (0.1, 1) if scale == "log" else (0, 1)

        T = slice_dict["transform"]
        aspect = slice_dict["aspect"]

        transform = Affine2D(T)

        self.T_inv = np.linalg.inv(T)

        self.ax_slice.format_coord = self.__format_axis_coord

        if self.cb_slice is not None:
            self.cb_slice.remove()
            self.cb_slice = None

        self.ax_slice.remove()

        # if self.ax_xint:
        #     self.ax_xint.remove()
        # if self.ax_yint:
        #     self.ax_yint.remove()

        extreme_finder = ExtremeFinderSimple(20, 20)

        grid_locator1 = MaxNLocator(nbins=10)
        grid_locator2 = MaxNLocator(nbins=10)

        grid_locator1.set_params(integer=True)
        grid_locator2.set_params(integer=True)

        grid_helper = GridHelperCurveLinear(
            transform,
            extreme_finder=extreme_finder,
            grid_locator1=grid_locator1,
            grid_locator2=grid_locator2,
        )

        self.ax_slice = self.fig_slice.add_subplot(
            1, 1, 1, axes_class=Axes, grid_helper=grid_helper
        )

        # self.ax_slice.set_xlabel(labels[0])
        # self.ax_slice.set_ylabel(labels[1])
        self.ax_slice.set_aspect(aspect)

        # divider = make_axes_locatable(self.ax_slice)

        # self.ax_yint = divider.append_axes('right',
        #                                    '10%',
        #                                    pad=0.15,
        #                                    sharey=self.ax_slice)

        # self.ax_xint = divider.append_axes('top',
        #                                    '10%',
        #                                    pad=0.15,
        #                                    sharex=self.ax_slice)

        trans = transform + self.ax_slice.transData

        im = self.ax_slice.pcolormesh(
            x,
            y,
            clip,
            norm=scale,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            shading="flat",
            transform=trans,
            rasterized=True,
        )

        self.ax_slice.set_xlabel(labels[0])
        self.ax_slice.set_ylabel(labels[1])

        # self.ax_slice.set_xticks([])
        # self.ax_slice.set_yticks([])

        # xlim = self.ax_slice.get_xlim()
        # ylim = self.ax_slice.get_ylim()

        # ascale = (ylim[1] - ylim[0]) / (xlim[1] - xlim[0]) * aspect

        # xstart = 1+0.05 if ascale > 1 else 1+0.05*ascale
        # ystart = 1+0.05 if ascale < 1 else 1+0.05*ascale

        # xwidth = 0.1 if ascale < 1 else 0.1 * ascale
        # ywidth = 0.1 if ascale > 1 else 0.1 / ascale

        # self.ax_xint = self.ax_slice.inset_axes(
        #     [0, 0 - ywidth, 1, ywidth], sharex=self.ax_slice
        # )

        # self.ax_yint = self.ax_slice.inset_axes(
        #     [0 - xwidth, 0, xwidth, 1], sharey=self.ax_slice
        # )

        # xint = signal.sum(axis=0)
        # yint = signal.sum(axis=1)
        # sigx = np.sqrt(xint)
        # sigy = np.sqrt(yint)

        # self.ax_xint.errorbar(
        #     0.5 * (x[1:] + x[:-1]),
        #     xint,
        #     yerr=sigx,
        #     fmt=".",
        #     linestyle="-",
        #     color="C0",
        # )

        # self.ax_yint.errorbar(
        #     yint,
        #     0.5 * (y[1:] + y[:-1]),
        #     xerr=sigy,
        #     fmt=".",
        #     linestyle="-",
        #     color="C1",
        # )

        # self.ax_xint.minorticks_on()
        # self.ax_yint.minorticks_on()

        # self.ax_xint.xaxis.get_major_locator().set_params(integer=True)
        # self.ax_yint.yaxis.get_major_locator().set_params(integer=True)

        # self.ax_xint.set_xlabel(labels[0])
        # self.ax_yint.set_ylabel(labels[1])

        # self.ax_xint.yaxis.tick_right()
        # self.ax_yint.xaxis.tick_top()

        # self.ax_xint.ticklabel_format(style="sci", axis="y", scilimits=(0, 0))
        # self.ax_yint.ticklabel_format(style="sci", axis="x", scilimits=(0, 0))

        # self.ax_xint.grid(True)
        # self.ax_yint.grid(True)

        # self.ax_xint.set_xticks([])
        # self.ax_yint.set_yticks([])

        self.im = im
        self.vmin, self.vmax = self.im.norm.vmin, self.im.norm.vmax

        self.ax_slice.set_title(title)
        self.ax_slice.grid(True)

        # ax = [self.ax_slice, self.ax_xint, self.ax_yint]

        # cax = self.ax_yint.inset_axes([1.1, 0, 0.25, 1])

        self.cb_slice = self.fig_slice.colorbar(self.im, ax=self.ax_slice)
        self.cb_slice.minorticks_on()

        # self.fig_slice.tight_layout()

        self.canvas_slice.draw_idle()
        self.canvas_slice.flush_events()
