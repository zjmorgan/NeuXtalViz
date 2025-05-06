import sys

from qtpy.QtWidgets import (
    QWidget,
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
    QFileDialog,
)

from qtpy.QtGui import QDoubleValidator, QIntValidator
from qtpy.QtCore import Qt, Signal

import numpy as np

# import matplotlib.pyplot as plt
import pyvista as pv

from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.ticker import FormatStrFormatter

from NeuXtalViz.views.base_view import NeuXtalVizWidget


class ExperimentView(NeuXtalVizWidget):
    roi_ready = Signal()
    viz_ready = Signal()

    def __init__(self, parent=None):
        super().__init__(parent)

        self.tab_widget = QTabWidget(self)

        self.coverage_tab()
        self.peak_tab()

        self.layout().addWidget(self.tab_widget, stretch=1)

    def coverage_tab(self):
        cov_tab = QWidget()
        self.tab_widget.addTab(cov_tab, "Coverage")

        coverage_layout = QVBoxLayout()

        self.instrument_combo = QComboBox(self)
        self.instrument_combo.addItem("TOPAZ")
        self.instrument_combo.addItem("MANDI")
        self.instrument_combo.addItem("CORELLI")
        self.instrument_combo.addItem("SNAP")
        self.instrument_combo.addItem("WAND²")
        self.instrument_combo.addItem("DEMAND")

        self.mode_combo = QComboBox(self)

        self.crystal_combo = QComboBox(self)
        self.point_group_combo = QComboBox(self)
        self.lattice_centering_combo = QComboBox(self)

        self.crystal_combo.addItem("Triclinic")
        self.crystal_combo.addItem("Monoclinic")
        self.crystal_combo.addItem("Orthorhombic")
        self.crystal_combo.addItem("Tetragonal")
        self.crystal_combo.addItem("Trigonal/Rhombohedral")
        self.crystal_combo.addItem("Trigonal/Hexagonal")
        self.crystal_combo.addItem("Hexagonal")
        self.crystal_combo.addItem("Cubic")

        self.load_UB_button = QPushButton("Load UB", self)
        self.optimize_button = QPushButton("Optimize Coverage", self)
        self.delete_button = QPushButton("Delete Highlighted", self)
        self.highlight_button = QPushButton("Highlight All", self)

        self.count_combo = QComboBox(self)
        self.update_button = QPushButton("Update Highlighted", self)
        self.count_line = QLineEdit("1.0")
        self.title_line = QLineEdit("Scan Title")

        notation = QDoubleValidator.StandardNotation
        validator = QDoubleValidator(0.001, 10000, 5, notation=notation)

        self.count_line.setValidator(validator)

        self.save_plan_button = QPushButton("Save CSV", self)
        self.save_experiment_button = QPushButton("Save Experiment", self)
        self.load_experiment_button = QPushButton("Load Experiment", self)

        self.wl_min_line = QLineEdit("0.4")
        self.wl_max_line = QLineEdit("3.5")

        self.d_min_line = QLineEdit("0.7")

        wl_label = QLabel("λ:")
        d_min_label = QLabel("d(min):")

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0.2, 10, 5, notation=notation)

        self.wl_min_line.setValidator(validator)
        self.wl_max_line.setValidator(validator)

        validator = QDoubleValidator(0.4, 10, 5, notation=notation)

        self.d_min_line.setValidator(validator)

        angstrom_label = QLabel("Å")

        settings_label = QLabel("Settings")
        self.settings_line = QLineEdit("20")

        validator = QIntValidator(1, 1000)

        self.settings_line.setValidator(validator)

        resize = QHeaderView.Stretch

        self.goniometer_table = QTableWidget()

        self.goniometer_table.setRowCount(0)
        self.goniometer_table.setColumnCount(3)

        labels = ["Motor", "Min", "Max"]

        self.goniometer_table.horizontalHeader().setStretchLastSection(True)
        self.goniometer_table.horizontalHeader().setSectionResizeMode(resize)
        self.goniometer_table.setHorizontalHeaderLabels(labels)

        self.motor_table = QTableWidget()

        self.motor_table.setRowCount(0)
        self.motor_table.setColumnCount(2)

        self.plan_table = QTableWidget()

        labels = ["Motor", "Value"]

        self.motor_table.horizontalHeader().setStretchLastSection(True)
        self.motor_table.horizontalHeader().setSectionResizeMode(resize)
        self.motor_table.setHorizontalHeaderLabels(labels)

        settings_layout = QHBoxLayout()

        settings_layout.addWidget(self.load_UB_button)
        settings_layout.addWidget(self.crystal_combo)
        settings_layout.addWidget(self.point_group_combo)
        settings_layout.addWidget(self.lattice_centering_combo)

        params_layout = QHBoxLayout()

        params_layout.addWidget(self.instrument_combo)
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
        plan_tab = QWidget()

        goniometer_layout = QVBoxLayout()
        motor_layout = QVBoxLayout()
        plan_layout = QVBoxLayout()

        mode_layout = QHBoxLayout()
        mode_layout.addWidget(self.mode_combo)
        mode_layout.addStretch(1)

        planning_layout = QHBoxLayout()
        planning_layout.addWidget(self.title_line)
        planning_layout.addWidget(self.count_combo)
        planning_layout.addWidget(self.count_line)
        planning_layout.addWidget(self.update_button)
        planning_layout.addStretch(1)
        planning_layout.addWidget(settings_label)
        planning_layout.addWidget(self.settings_line)
        planning_layout.addWidget(self.optimize_button)

        save_layout = QHBoxLayout()
        save_layout.addWidget(self.delete_button)
        save_layout.addWidget(self.highlight_button)
        save_layout.addStretch(1)
        save_layout.addWidget(self.save_plan_button)
        save_layout.addWidget(self.save_experiment_button)
        save_layout.addWidget(self.load_experiment_button)

        goniometer_layout.addLayout(mode_layout)
        goniometer_layout.addWidget(self.goniometer_table)
        motor_layout.addWidget(self.motor_table)
        plan_layout.addLayout(planning_layout)
        plan_layout.addWidget(self.plan_table)
        plan_layout.addLayout(save_layout)

        goniometer_tab.setLayout(goniometer_layout)
        motor_tab.setLayout(motor_layout)
        plan_tab.setLayout(plan_layout)

        values_tab.addTab(goniometer_tab, "Goniometers")
        values_tab.addTab(motor_tab, "Motors")
        values_tab.addTab(plan_tab, "Plan")

        result_layout.addWidget(values_tab)

        self.canvas_cov = FigureCanvas(
            Figure(constrained_layout=True, figsize=(6.4, 4.8))
        )

        result_layout.addWidget(NavigationToolbar2QT(self.canvas_cov, self))
        result_layout.addWidget(self.canvas_cov)

        fig = self.canvas_cov.figure

        self.ax_cov = fig.subplots(3, 1, sharex=True)
        self.ax_cov[2].set_xlabel("Resolution [Å]")
        self.ax_cov[0].set_ylabel("Completeness [%]")
        self.ax_cov[1].set_ylabel("Multiplicity")
        self.ax_cov[2].set_ylabel("Reflections")

        coverage_layout.addLayout(settings_layout)
        coverage_layout.addLayout(params_layout)
        coverage_layout.addLayout(result_layout)

        cov_tab.setLayout(coverage_layout)

    def peak_tab(self):
        inst_tab = QWidget()
        self.tab_widget.addTab(inst_tab, "Peak")

        peak_layout = QVBoxLayout()

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

        self.calculate_single_button = QPushButton("Individual Peak", self)
        self.calculate_double_button = QPushButton("Simultaneous Peaks", self)

        calculator_layout.addWidget(h_label, 0, 1, Qt.AlignCenter)
        calculator_layout.addWidget(k_label, 0, 2, Qt.AlignCenter)
        calculator_layout.addWidget(l_label, 0, 3, Qt.AlignCenter)

        calculator_layout.addWidget(peak_1_label, 1, 0)
        calculator_layout.addWidget(self.h1_line, 1, 1)
        calculator_layout.addWidget(self.k1_line, 1, 2)
        calculator_layout.addWidget(self.l1_line, 1, 3)
        calculator_layout.addWidget(self.calculate_single_button, 1, 4)

        calculator_layout.addWidget(peak_2_label, 2, 0)
        calculator_layout.addWidget(self.h2_line, 2, 1)
        calculator_layout.addWidget(self.k2_line, 2, 2)
        calculator_layout.addWidget(self.l2_line, 2, 3)
        calculator_layout.addWidget(self.calculate_double_button, 2, 4)

        peak_layout.addLayout(calculator_layout)

        self.canvas_inst = FigureCanvas(Figure(constrained_layout=True))

        peak_layout.addWidget(NavigationToolbar2QT(self.canvas_inst, self))
        peak_layout.addWidget(self.canvas_inst)

        self.fig_inst = self.canvas_inst.figure
        self.ax_inst = self.fig_inst.subplots(1, 1)
        self.ax_inst.clear()
        self.ax_inst.invert_xaxis()

        self.cb_inst = None
        self.cb_inst_alt = None

        orientation_layout = QHBoxLayout()

        self.add_button = QPushButton("Add Orientation", self)

        self.angles_line = QLineEdit()
        self.horizontal_line = QLineEdit()
        self.vertical_line = QLineEdit()

        self.angles_line.setEnabled(False)
        self.horizontal_line.setEnabled(False)
        self.vertical_line.setEnabled(False)

        self.angles_combo = QComboBox(self)

        orientation_layout.addWidget(self.angles_combo)
        orientation_layout.addWidget(self.angles_line)
        orientation_layout.addWidget(self.horizontal_line)
        orientation_layout.addWidget(self.vertical_line)
        orientation_layout.addWidget(self.add_button)

        peak_layout.addLayout(orientation_layout)

        stretch = QHeaderView.Stretch

        self.peaks_table = QTableWidget()
        self.peaks_table.setRowCount(0)
        self.peaks_table.setColumnCount(5)

        header = ["h", "k", "l", "d", "λ"]

        self.peaks_table.horizontalHeader().setSectionResizeMode(stretch)
        self.peaks_table.setHorizontalHeaderLabels(header)
        self.peaks_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.peaks_table.setSelectionBehavior(QTableWidget.SelectRows)
        self.peaks_table.setSortingEnabled(True)

        peak_layout.addWidget(self.peaks_table)

        inst_tab.setLayout(peak_layout)

    def connect_peak_table(self, update_table):
        self.angles_combo.activated.connect(update_table)

    def connect_add_orientation(self, add_orientation):
        self.add_button.clicked.connect(add_orientation)

    def connect_delete_angles(self, delete_angles):
        self.delete_button.clicked.connect(delete_angles)

    def connect_highlight_angles(self, highlight_angles):
        self.highlight_button.clicked.connect(highlight_angles)

    def connect_calculate_single(self, calculate_single):
        self.calculate_single_button.clicked.connect(calculate_single)

    def connect_calculate_double(self, calculate_double):
        self.calculate_double_button.clicked.connect(calculate_double)

    def connect_switch_crystal(self, switch_crystal):
        self.crystal_combo.activated.connect(switch_crystal)

    def connect_switch_point_group(self, switch_group):
        self.point_group_combo.activated.connect(switch_group)

    def connect_switch_lattice_centering(self, switch_centering):
        self.lattice_centering_combo.activated.connect(switch_centering)

    def connect_switch_instrument(self, switch_instrument):
        self.instrument_combo.activated.connect(switch_instrument)

    def connect_update_goniometer(self, update_goniometer):
        self.mode_combo.activated.connect(update_goniometer)

    def connect_optimize(self, optimize):
        self.optimize_button.clicked.connect(optimize)

    def connect_load_UB(self, load_UB):
        self.load_UB_button.clicked.connect(load_UB)

    def connect_save_CSV(self, save_CSV):
        self.save_plan_button.clicked.connect(save_CSV)

    def connect_load_experiment(self, load_experiment):
        self.load_experiment_button.clicked.connect(load_experiment)

    def connect_save_experiment(self, save_experiment):
        self.save_experiment_button.clicked.connect(save_experiment)

    def connect_wavelength(self, update_wavelength):
        self.wl_min_line.editingFinished.connect(update_wavelength)

    def connect_update(self, update):
        self.update_button.clicked.connect(update)

    def load_UB_file_dialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filename, _ = file_dialog.getOpenFileName(
            self, "Load UB file", "", "UB files (*.mat)", options=options
        )

        return filename

    def save_CSV_file_dialog(self, path=""):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filename, _ = file_dialog.getSaveFileName(
            self,
            "Save peaks file",
            path,
            "Experiment files (*.csv)",
            options=options,
        )

        if filename is not None:
            if not filename.endswith(".csv"):
                filename += ".csv"

        return filename

    def load_experiment_file_dialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filename, _ = file_dialog.getOpenFileName(
            self,
            "Load experiment file",
            "",
            "Experiment files (*.nxs)",
            options=options,
        )

        return filename

    def save_experiment_file_dialog(self, path=""):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filename, _ = file_dialog.getSaveFileName(
            self,
            "Save experiment file",
            path,
            "Experiment files (*.nxs)",
            options=options,
        )

        if filename is not None:
            if not filename.endswith(".nxs"):
                filename += ".nxs"

        return filename

    def get_d_min(self):
        if self.d_min_line.hasAcceptableInput():
            return float(self.d_min_line.text())

    def set_d_min(self, d_min):
        self.d_min_line.setText(str(d_min))

    def get_crystal_system(self):
        return self.crystal_combo.currentText()

    def set_crystal_system(self, crystal_system):
        index = self.crystal_combo.findText(crystal_system)
        if index >= 0:
            self.crystal_combo.blockSignals(True)
            self.crystal_combo.setCurrentIndex(index)
            self.crystal_combo.blockSignals(False)

    def set_point_groups(self, groups):
        self.point_group_combo.clear()
        for group in groups:
            self.point_group_combo.addItem(group)

    def get_point_group(self):
        return self.point_group_combo.currentText()

    def set_point_group(self, point_group):
        index = self.point_group_combo.findText(point_group)
        if index >= 0:
            self.point_group_combo.blockSignals(True)
            self.point_group_combo.setCurrentIndex(index)
            self.point_group_combo.blockSignals(False)

    def set_lattice_centerings(self, centerings):
        self.lattice_centering_combo.clear()
        for centering in centerings:
            self.lattice_centering_combo.addItem(centering)

    def get_lattice_centering(self):
        return self.lattice_centering_combo.currentText()

    def set_lattice_centering(self, lattice_centering):
        index = self.lattice_centering_combo.findText(lattice_centering)
        if index >= 0:
            self.lattice_centering_combo.blockSignals(True)
            self.lattice_centering_combo.setCurrentIndex(index)
            self.lattice_centering_combo.blockSignals(False)

    def get_mode(self):
        return self.mode_combo.currentText()

    def set_mode(self, mode):
        index = self.mode_combo.findText(mode)
        if index >= 0:
            self.mode_combo.blockSignals(True)
            self.mode_combo.setCurrentIndex(index)
            self.mode_combo.blockSignals(False)

    def set_modes(self, modes):
        self.mode_combo.clear()
        for mode in modes:
            self.mode_combo.addItem(mode)

    def set_counting_options(self, options):
        self.count_combo.clear()
        for option in options:
            self.count_combo.addItem(option)

    def get_counting_options(self):
        return [
            self.count_combo.itemText(i)
            for i in range(self.count_combo.count())
        ]

    def get_counting_index(self):
        return self.count_combo.currentIndex()

    def get_count_value(self):
        if self.count_line.hasAcceptableInput():
            return float(self.count_line.text())

    def set_peak_list(self, rows):
        self.angles_combo.blockSignals(True)
        self.angles_combo.clear()
        self.angles_combo.addItem("0")
        for row in range(rows):
            self.angles_combo.addItem((str(row + 1)))
        self.angles_combo.blockSignals(False)

    def get_peak_list(self):
        val = self.angles_combo.currentText()
        if val is not None:
            if val.isdigit():
                return int(val) - 1

    def set_wavelength(self, wavelength):
        self.wl_min_line.blockSignals(True)
        self.wl_max_line.blockSignals(True)
        if type(wavelength) is list:
            self.wl_min_line.setText(str(wavelength[0]))
            self.wl_max_line.setText(str(wavelength[1]))
            self.wl_max_line.setEnabled(True)
        else:
            self.wl_min_line.setText(str(wavelength))
            self.wl_max_line.setText(str(wavelength))
            self.wl_max_line.setEnabled(False)
        self.wl_min_line.blockSignals(False)
        self.wl_max_line.blockSignals(False)

    def get_wavelength(self):
        params = self.wl_min_line, self.wl_max_line

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:
            return [float(param.text()) for param in params]

    def update_wavelength(self, lamda_min):
        if not self.wl_max_line.isEnabled():
            self.wl_max_line.setText(str(lamda_min))

    def update_tables(self, title, goniometers, motors):
        self.goniometer_table.clearContents()
        self.goniometer_table.setRowCount(0)
        self.goniometer_table.setRowCount(len(goniometers))

        free = []
        for row, gon in enumerate(goniometers):
            angle, amin, amax = gon
            amin, amax = str(amin), str(amax)
            self.goniometer_table.setItem(row, 0, QTableWidgetItem(angle))
            self.goniometer_table.setItem(row, 1, QTableWidgetItem(amin))
            self.goniometer_table.setItem(row, 2, QTableWidgetItem(amax))
            item = self.goniometer_table.item(row, 0)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            if float(amin) == float(amax):
                for j in [1, 2]:
                    item = self.goniometer_table.item(row, j)
                    item.setFlags(item.flags() & ~Qt.ItemIsEditable)
            else:
                free.append(angle)

        self.motor_table.setRowCount(0)
        self.motor_table.setRowCount(len(motors))

        for row, mot in enumerate(motors):
            setting, val = mot
            val = str(val)
            self.motor_table.setItem(row, 0, QTableWidgetItem(setting))
            self.motor_table.setItem(row, 1, QTableWidgetItem(val))
            item = self.motor_table.item(row, 0)
            item.setFlags(item.flags() & ~Qt.ItemIsEditable)

        self.plan_table.blockSignals(True)
        self.plan_table.clearContents()
        self.plan_table.setRowCount(0)
        self.plan_table.setColumnCount(0)
        self.plan_table.setColumnCount(len(free) + 5)

        labels = [title] + free + ["Comment", "Wait For", "Value", "Use"]

        resize = QHeaderView.Stretch

        self.plan_table.horizontalHeader().setStretchLastSection(True)
        self.plan_table.horizontalHeader().setSectionResizeMode(resize)
        self.plan_table.setHorizontalHeaderLabels(labels)

        self.plan_table.itemChanged.connect(self.handle_item_changed)
        self.plan_table.blockSignals(False)

    def delete_angles(self):
        self.plan_table.blockSignals(True)

        rows = set(index.row() for index in self.plan_table.selectedIndexes())
        for row in sorted(rows, reverse=True):
            self.plan_table.removeRow(row)

        self.set_peak_list(self.get_number_of_orientations())
        self.plan_table.blockSignals(False)

        return rows

    def highlight_angles(self):
        self.plan_table.setSelectionBehavior(self.plan_table.SelectRows)
        self.plan_table.selectAll()

    def get_title(self):
        return self.title_line.text()

    def update_counting(self):
        self.plan_table.blockSignals(True)

        title = self.get_title()
        index = self.get_counting_index()
        value = self.get_count_value()

        col = self.plan_table.columnCount() - 3

        rows = set(index.row() for index in self.plan_table.selectedIndexes())
        for row in rows:
            if title is not None:
                item = QTableWidgetItem(title)
                self.plan_table.setItem(row, 0, item)
            if index is not None:
                widget = self.plan_table.cellWidget(row, col)
                if isinstance(widget, QComboBox):
                    widget.setCurrentIndex(index)
            if value is not None:
                item = QTableWidgetItem("{:.3f}".format(value))
                self.plan_table.setItem(row, col + 1, item)

        self.plan_table.blockSignals(False)

        return rows

    def get_all_angles(self):
        rows = self.goniometer_table.rowCount()

        angles = [
            self.goniometer_table.item(row, 0).text() for row in range(rows)
        ]

        return angles

    def get_free_angles(self):
        cols = self.plan_table.columnCount() - 5

        angles = [
            self.plan_table.horizontalHeaderItem(i + 1).text()
            for i in range(cols)
        ]

        return angles

    def get_number_of_orientations(self):
        return self.plan_table.rowCount()

    def get_orientations_to_use(self):
        col = self.plan_table.columnCount() - 1

        use = []
        for row in range(self.get_number_of_orientations()):
            item = self.plan_table.item(row, col)
            use.append(item.checkState() == Qt.Checked)

        return use

    def get_all_titles(self):
        title = []
        for row in range(self.get_number_of_orientations()):
            item = self.plan_table.item(row, 0).text()
            title.append(item)

        return title

    def get_all_values(self):
        col = self.plan_table.columnCount() - 2

        value = []
        for row in range(self.get_number_of_orientations()):
            item = self.plan_table.item(row, col).text()
            item = float(item) if item.replace(".", "").isnumeric() else 0.0
            value.append(item)

        return value

    def get_all_countings(self):
        col = self.plan_table.columnCount() - 3

        count = []
        for row in range(self.get_number_of_orientations()):
            widget = self.plan_table.cellWidget(row, col)
            if isinstance(widget, QComboBox):
                count.append(widget.currentText())
            else:
                count.append("")

        return count

    def get_all_comments(self):
        col = self.plan_table.columnCount() - 4

        comment = []
        for row in range(self.get_number_of_orientations()):
            comment.append(self.plan_table.item(row, col).text())

        return comment

    def get_all_settings(self):
        settings = []
        for row in range(self.get_number_of_orientations()):
            setting = self.get_angle_setting(row)
            settings.append(setting)

        return settings

    def get_settings(self):
        if self.settings_line.hasAcceptableInput():
            return int(self.settings_line.text())

    def get_optimized_settings(self):
        col = self.plan_table.columnCount() - 5

        opt = []
        for row in range(self.get_number_of_orientations()):
            item = self.plan_table.item(row, col + 1)
            opt.append(item.text() == "CrystalPlan")

        return opt

    def get_angle_setting(self, row):
        cols = self.plan_table.columnCount() - 5

        setting = []
        for col in range(cols):
            setting.append(float(self.plan_table.item(row, col + 1).text()))

        return setting

    def add_orientation(self, title, comment, angles):
        row = self.get_number_of_orientations()
        self.plan_table.blockSignals(True)
        self.plan_table.setSortingEnabled(False)
        self.plan_table.setRowCount(row + 1)

        col = 0

        item = QTableWidgetItem(title)
        self.plan_table.setItem(row, col, item)
        col += 1

        for angle in angles:
            item = QTableWidgetItem("{:.1f}".format(angle))
            self.plan_table.setItem(row, col, item)
            col += 1

        self.plan_table.setItem(row, col, QTableWidgetItem(comment))
        col += 1

        combobox = QComboBox()
        options = self.get_counting_options()
        for option in options:
            combobox.addItem(option)
        index = self.get_counting_index()
        if index is not None:
            combobox.setCurrentIndex(index)
        self.plan_table.setCellWidget(row, col, combobox)
        col += 1

        val = self.get_count_value()
        if val is not None:
            item = QTableWidgetItem("{:.3f}".format(val))
            self.plan_table.setItem(row, col, item)
        col += 1

        flags = Qt.ItemIsUserCheckable | Qt.ItemIsEnabled

        checkbox = QTableWidgetItem("")
        checkbox.setText("")
        checkbox.setFlags(flags)
        checkbox.setCheckState(Qt.Checked)
        self.plan_table.setItem(row, col, checkbox)

        self.set_peak_list(self.get_number_of_orientations())
        self.plan_table.blockSignals(False)
        self.plan_table.setSortingEnabled(True)

    def add_settings(self, titles, settings, comments, counts, values, use):
        self.plan_table.setUpdatesEnabled(False)
        self.plan_table.setSortingEnabled(False)
        self.plan_table.blockSignals(True)
        self.plan_table.clearContents()
        self.plan_table.setRowCount(len(use))

        for row, angles in enumerate(settings):
            col = 0

            item = QTableWidgetItem(titles[row])
            self.plan_table.setItem(row, col, item)
            col += 1

            for angle in angles:
                item = QTableWidgetItem("{:.1f}".format(angle))
                self.plan_table.setItem(row, col, item)
                col += 1

            self.plan_table.setItem(row, col, QTableWidgetItem(comments[row]))
            col += 1

            combobox = QComboBox()
            options = self.get_counting_options()
            for option in options:
                combobox.addItem(option)
            index = options.index(counts[row])
            combobox.setCurrentIndex(index)
            self.plan_table.setCellWidget(row, col, combobox)
            col += 1

            item = QTableWidgetItem("{:.3f}".format(values[row]))
            self.plan_table.setItem(row, col, item)
            col += 1

            flags = Qt.ItemIsUserCheckable | Qt.ItemIsEnabled

            checkbox = QTableWidgetItem("")
            checkbox.setText("")
            checkbox.setFlags(flags)
            checkbox.setCheckState(Qt.Checked if use[row] else Qt.Unchecked)
            self.plan_table.setItem(row, col, checkbox)

        self.plan_table.setUpdatesEnabled(True)
        self.plan_table.blockSignals(False)
        self.plan_table.setSortingEnabled(True)

        self.set_peak_list(self.get_number_of_orientations())

    def handle_item_changed(self, item):
        self.plan_table.blockSignals(True)

        col = item.column()

        if col == self.plan_table.columnCount() - 1:
            self.viz_ready.emit()

        self.plan_table.blockSignals(False)

    def connect_viz_ready(self, visualize):
        self.viz_ready.connect(visualize)

    def get_instrument(self):
        return self.instrument_combo.currentText()

    def set_instrument(self, instrument):
        index = self.instrument_combo.findText(instrument)
        if index >= 0:
            self.instrument_combo.blockSignals(True)
            self.instrument_combo.setCurrentIndex(index)
            self.instrument_combo.blockSignals(False)

    def get_motors(self):
        logs = {}
        for row in range(self.motor_table.rowCount()):
            setting = self.motor_table.item(row, 0).text()
            logs[setting] = float(self.motor_table.item(row, 1).text())

        return logs

    def set_motors(self, values):
        for row, value in enumerate(values):
            self.motor_table.setItem(row, 1, self.set_item_value(str(value)))

    def get_goniometer_limits(self):
        limits = []
        for row in range(self.goniometer_table.rowCount()):
            amin = float(self.goniometer_table.item(row, 1).text())
            amax = float(self.goniometer_table.item(row, 2).text())
            limits.append([amin, amax])

        return limits

    def set_goniometer_limits(self, limits):
        for row, limit in enumerate(limits):
            amin, amax = str(limit[0]), str(limit[1])
            self.goniometer_table.setItem(row, 1, self.set_item_value(amin))
            self.goniometer_table.setItem(row, 2, self.set_item_value(amax))

    def add_peaks(self, peak_dict):
        self.plotter.clear_actors()

        coords = np.array(peak_dict["coords"])
        colors = np.array(peak_dict["colors"])
        # sizes = np.array(peak_dict['sizes'])

        points = pv.PolyData(coords)
        points["colors"] = colors
        # points['sizes'] = 5*sizes

        self.plotter.add_mesh(
            points,
            scalars="colors",
            rgb=True,
            smooth_shading=True,
            point_size=10,
            render_points_as_spheres=True,
        )

        self.plotter.enable_depth_peeling()
        # self.plotter.add_axes_at_origin()

        coords = np.array(peak_dict["axis_coords"])
        colors = np.array(peak_dict["axis_colors"])

        for i in range(3):
            arrow = pv.Arrow([0, 0, 0], coords[i], scale="auto")
            self.plotter.add_mesh(arrow, color=colors[i], smooth_shading=True)

        radius = 0.2 * np.sqrt(np.min(np.sum(coords**2, axis=1)))
        sphere = pv.Sphere(radius=radius)

        self.plotter.add_mesh(sphere, color="w", smooth_shading=True)

        Q_max = 2 * np.pi / peak_dict["axis_limit"]

        mesh = pv.Line(
            pointa=(-Q_max, 0, 0), pointb=(Q_max, 0, 0), resolution=1
        )

        self.plotter.add_mesh(
            mesh, color="k", style="wireframe", render_lines_as_tubes=True
        )

        mesh = pv.Line(
            pointa=(0, -Q_max, 0), pointb=(0, Q_max, 0), resolution=1
        )

        self.plotter.add_mesh(
            mesh, color="k", style="wireframe", render_lines_as_tubes=True
        )

        mesh = pv.Line(
            pointa=(0, 0, -Q_max), pointb=(0, 0, Q_max), resolution=1
        )

        self.plotter.add_mesh(
            mesh, color="k", style="wireframe", render_lines_as_tubes=True
        )

        self.reset_view()

    def update_peaks_table(self, peaks):
        self.peaks_table.clearSelection()
        self.peaks_table.setSortingEnabled(False)
        self.peaks_table.setRowCount(0)
        self.peaks_table.setRowCount(len(peaks))

        for row, peak in enumerate(peaks):
            self.set_peak(row, peak)

        self.peaks_table.setSortingEnabled(True)

    def set_peak(self, row, peak):
        h, k, l, d, lamda = peak
        h = "{:.3f}".format(h)
        k = "{:.3f}".format(k)
        l = "{:.3f}".format(l)
        d = "{:.4f}".format(d)
        lamda = "{:.4f}".format(lamda)
        self.peaks_table.setItem(row, 0, self.set_item_value(h))
        self.peaks_table.setItem(row, 1, self.set_item_value(k))
        self.peaks_table.setItem(row, 2, self.set_item_value(l))
        self.peaks_table.setItem(row, 3, self.set_item_value(d))
        self.peaks_table.setItem(row, 4, self.set_item_value(lamda))

    def set_item_value(self, value):
        item = QTableWidgetItem()
        item.setData(Qt.DisplayRole, float(value))
        return item

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

    def plot_statistics(self, shel, comp, mult, refl):
        self.ax_cov[0].clear()
        self.ax_cov[1].clear()
        self.ax_cov[2].clear()

        self.ax_cov[0].bar(shel, comp, color="C0")
        self.ax_cov[1].bar(shel, mult, color="C1")
        self.ax_cov[2].bar(shel, refl, color="C2")

        self.ax_cov[0].set_ylim(0, 100)

        self.ax_cov[0].minorticks_on()
        self.ax_cov[1].minorticks_on()
        self.ax_cov[2].minorticks_on()

        self.ax_cov[2].set_xlabel("Resolution [Å]")
        self.ax_cov[0].set_ylabel("Completeness [%]")
        self.ax_cov[1].set_ylabel("Redundancy")
        self.ax_cov[2].set_ylabel("Reflections")

        self.canvas_cov.draw_idle()
        self.canvas_cov.flush_events()

    def plot_instrument(self, gamma_inst, nu_inst, gamma, nu, lamda):
        if self.cb_inst is not None:
            self.cb_inst.remove()
            self.cb_inst = None

        if self.cb_inst_alt is not None:
            self.cb_inst_alt.remove()
            self.cb_inst_alt = None

        self.ax_inst.clear()
        self.ax_inst.invert_xaxis()

        self.ax_inst.scatter(
            gamma_inst, nu_inst, color="lightgray", marker="o", rasterized=True
        )

        self.im = self.ax_inst.scatter(
            gamma, nu, c=lamda, marker="o", rasterized=True
        )

        self.ax_inst.set_aspect(1)
        self.ax_inst.minorticks_on()

        self.ax_inst.set_xlabel(r"$\gamma$")
        self.ax_inst.set_ylabel(r"$\nu$")

        fmt_str_form = FormatStrFormatter(r"$%d^\circ$")

        self.ax_inst.xaxis.set_major_formatter(fmt_str_form)
        self.ax_inst.yaxis.set_major_formatter(fmt_str_form)

        if len(lamda) > 0:
            self.cb_inst = self.fig_inst.colorbar(
                self.im, ax=self.ax_inst, orientation="horizontal"
            )
            self.cb_inst.minorticks_on()
            self.cb_inst.ax.set_xlabel(r"$\lambda$ [Å]")

        self.fig_inst.canvas.mpl_connect(
            "button_press_event", self.on_press_inst
        )

        self.canvas_inst.draw_idle()
        self.canvas_inst.flush_events()

    def plot_instrument_alternate(
        self,
        gamma_inst,
        nu_inst,
        gamma_1,
        nu_1,
        lamda_1,
        gamma_2,
        nu_2,
        lamda_2,
    ):
        if self.cb_inst is not None:
            self.cb_inst.remove()
            self.cb_inst = None

        if self.cb_inst_alt is not None:
            self.cb_inst_alt.remove()
            self.cb_inst_alt = None

        self.ax_inst.clear()
        self.ax_inst.invert_xaxis()

        self.ax_inst.scatter(
            gamma_inst, nu_inst, color="lightgray", marker="o", rasterized=True
        )

        self.im = self.ax_inst.scatter(
            gamma_1, nu_1, c=lamda_1, marker="o", cmap="GnBu", rasterized=True
        )

        self.im_alt = self.ax_inst.scatter(
            gamma_2, nu_2, c=lamda_2, marker="o", cmap="RdPu", rasterized=True
        )

        self.ax_inst.set_aspect(1)
        self.ax_inst.minorticks_on()

        self.ax_inst.set_xlabel(r"$\gamma$")
        self.ax_inst.set_ylabel(r"$\nu$")

        fmt_str_form = FormatStrFormatter(r"$%d^\circ$")

        self.ax_inst.xaxis.set_major_formatter(fmt_str_form)
        self.ax_inst.yaxis.set_major_formatter(fmt_str_form)

        if len(lamda_2) > 0:
            self.cb_inst_alt = self.fig_inst.colorbar(
                self.im_alt, ax=self.ax_inst, orientation="horizontal"
            )
            self.cb_inst_alt.minorticks_on()

        if len(lamda_1) > 0:
            self.cb_inst = self.fig_inst.colorbar(
                self.im, ax=self.ax_inst, orientation="horizontal"
            )
            self.cb_inst.minorticks_on()

        if len(lamda_2) > 0:
            self.cb_inst_alt.ax.set_xlabel(r"$\lambda$ [Å]")
        elif len(lamda_1) > 0:
            self.cb_inst.ax.set_xlabel(r"$\lambda$ [Å]")

        self.fig_inst.canvas.mpl_connect(
            "button_press_event", self.on_press_inst
        )

        self.canvas_inst.draw_idle()
        self.canvas_inst.flush_events()

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

    def set_angles(self, values):
        ang = "(" + ", ".join(np.array(values).astype(str)) + ")"

        self.angles_line.setText(ang)

    def get_angles(self):
        ang = self.angles_line.text()
        ang = ang.strip("(").strip(")").split(",")

        return [float(val) for val in ang if val != ""]

    def on_press_inst(self, event):
        if (
            event.inaxes == self.ax_inst
            and self.fig_inst.canvas.toolbar.mode == ""
        ):
            horz, vert = event.xdata, event.ydata
            self.set_horizontal(horz)
            self.set_vertical(vert)

            self.roi_ready.emit()

    def update_inst(self):
        for line in self.ax_inst.lines:
            line.remove()

        horz, vert = self.get_horizontal(), self.get_vertical()

        self.ax_inst.axvline(x=horz, color="k", linestyle="--")
        self.ax_inst.axhline(y=vert, color="k", linestyle="--")

        self.canvas_inst.draw_idle()
        self.canvas_inst.flush_events()

    def connect_roi_ready(self, lookup):
        self.roi_ready.connect(lookup)
