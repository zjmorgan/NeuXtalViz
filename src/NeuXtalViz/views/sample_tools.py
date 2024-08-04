import sys
import numpy as np

from qtpy.QtWidgets import (QWidget,
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
                            QFileDialog)

from qtpy.QtGui import QDoubleValidator, QIntValidator, QRegExpValidator
from PyQt5.QtCore import Qt, QRegExp

import pyvista as pv

from NeuXtalViz.views.base_view import NeuXtalVizWidget

class SampleView(NeuXtalVizWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.tab_widget = QTabWidget(self)

        self.sample_tab()

        self.layout().addWidget(self.tab_widget, stretch=1)

    def sample_tab(self):

        samp_tab = QWidget()
        self.tab_widget.addTab(samp_tab, 'Sample')

        self.sample_combo = QComboBox(self)
        self.sample_combo.addItem('Sphere')
        self.sample_combo.addItem('Cylinder')
        self.sample_combo.addItem('Plate')

        self.add_sample_button = QPushButton('Add Sample', self)
        self.load_UB_button = QPushButton('Load UB', self)

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0, 100, 5, notation=notation)

        param1_label = QLabel('Width', self)
        param2_label = QLabel('Height', self)
        param3_label = QLabel('Thickness', self)

        unit_label = QLabel('cm', self)

        self.param1_line = QLineEdit('0.50')
        self.param2_line = QLineEdit('0.50')
        self.param3_line = QLineEdit('0.50')

        self.param2_line.setDisabled(True)
        self.param3_line.setDisabled(True)

        self.param1_line.setValidator(validator)
        self.param2_line.setValidator(validator)
        self.param3_line.setValidator(validator)

        param_layout = QHBoxLayout()
        generate_layout = QHBoxLayout()

        generate_layout.addWidget(self.sample_combo)
        generate_layout.addWidget(self.load_UB_button)
        generate_layout.addWidget(self.add_sample_button)

        param_layout.addWidget(param1_label)
        param_layout.addWidget(self.param1_line)

        param_layout.addWidget(param2_label)
        param_layout.addWidget(self.param2_line)

        param_layout.addWidget(param3_label)
        param_layout.addWidget(self.param3_line)

        param_layout.addWidget(unit_label)

        material_layout = QHBoxLayout()

        self.chem_line = QLineEdit()
        self.Z_line = QLineEdit()
        self.V_line = QLineEdit()

        exp = '-^((\(\d+[A-Z][a-z]?\)|[DT]|[A-Z][a-z]?)(\d+(\.\d+)?)?)'+\
              '(-((\(\d+[A-Z][a-z]?\)|[DT]|[A-Z][a-z]?)(\d+(\.\d+)?)?))*$'

        regexp = QRegExp(exp)
        validator = QRegExpValidator(regexp)

        validator = QIntValidator(1, 10000, self)

        self.Z_line.setValidator(validator)

        validator = QDoubleValidator(0, 100000, 4, notation=notation)

        self.V_line.setValidator(validator)

        Z_label = QLabel('Z')
        V_label = QLabel('Ω')
        uc_vol_label = QLabel('Å^3')

        material_layout.addWidget(self.chem_line)
        material_layout.addWidget(Z_label)
        material_layout.addWidget(self.Z_line)
        material_layout.addWidget(V_label)
        material_layout.addWidget(self.V_line)
        material_layout.addWidget(uc_vol_label)

        cryst_layout = QGridLayout()

        scattering_label = QLabel('Scattering', self)
        absorption_label = QLabel('Absorption', self)

        sigma_label = QLabel('σ', self)
        mu_label = QLabel('μ', self)

        sigma_unit_label = QLabel('barn', self)
        mu_unit_label = QLabel('1/cm', self)

        self.sigma_a_line = QLineEdit()
        self.sigma_s_line = QLineEdit()

        self.mu_a_line = QLineEdit()
        self.mu_s_line = QLineEdit()

        self.sigma_a_line.setEnabled(False)
        self.sigma_s_line.setEnabled(False)

        self.mu_a_line.setEnabled(False)
        self.mu_s_line.setEnabled(False)

        cryst_layout.addWidget(scattering_label, 0, 1, Qt.AlignCenter)
        cryst_layout.addWidget(absorption_label, 0, 2, Qt.AlignCenter)

        cryst_layout.addWidget(sigma_label, 1, 0)
        cryst_layout.addWidget(self.sigma_a_line, 1, 1)
        cryst_layout.addWidget(self.sigma_s_line, 1, 2)
        cryst_layout.addWidget(sigma_unit_label, 1, 3)

        cryst_layout.addWidget(mu_label, 2, 0)
        cryst_layout.addWidget(self.mu_a_line, 2, 1)
        cryst_layout.addWidget(self.mu_s_line, 2, 2)
        cryst_layout.addWidget(mu_unit_label, 2, 3)

        N_label = QLabel('N', self)
        M_label = QLabel('M', self)
        n_label = QLabel('n', self)
        rho_label = QLabel('rho', self)
        v_label = QLabel('V', self)
        m_label = QLabel('m', self)

        N_unit_label = QLabel('atoms', self)
        M_unit_label = QLabel('g/mol', self)
        n_unit_label = QLabel('1/Å^3', self)
        rho_unit_label = QLabel('g/cm^3', self)
        v_unit_label = QLabel('cm^3', self)
        m_unit_label = QLabel('g', self)

        self.N_line = QLineEdit()
        self.M_line = QLineEdit()
        self.n_line = QLineEdit()
        self.rho_line = QLineEdit()
        self.v_line = QLineEdit()
        self.m_line = QLineEdit()

        self.N_line.setEnabled(False)
        self.M_line.setEnabled(False)
        self.n_line.setEnabled(False)
        self.rho_line.setEnabled(False)
        self.v_line.setEnabled(False)
        self.m_line.setEnabled(False)

        cryst_layout.addWidget(N_label, 3, 0)
        cryst_layout.addWidget(self.N_line, 3, 1)
        cryst_layout.addWidget(N_unit_label, 3, 2)

        cryst_layout.addWidget(M_label, 4, 0)
        cryst_layout.addWidget(self.M_line, 4, 1)
        cryst_layout.addWidget(M_unit_label, 4, 2)

        cryst_layout.addWidget(n_label, 5, 0)
        cryst_layout.addWidget(self.n_line, 5, 1)
        cryst_layout.addWidget(n_unit_label, 5, 2)

        cryst_layout.addWidget(rho_label, 6, 0)
        cryst_layout.addWidget(self.rho_line, 6, 1)
        cryst_layout.addWidget(rho_unit_label, 6, 2)

        cryst_layout.addWidget(v_label, 7, 0)
        cryst_layout.addWidget(self.v_line, 7, 1)
        cryst_layout.addWidget(v_unit_label, 7, 2)

        cryst_layout.addWidget(m_label, 8, 0)
        cryst_layout.addWidget(self.m_line, 8, 1)
        cryst_layout.addWidget(m_unit_label, 8, 2)

        sample_layout = QVBoxLayout()

        a_star_label = QLabel('a*', self)
        b_star_label = QLabel('b*', self)
        c_star_label = QLabel('c*', self)

        face_index_label = QLabel('Face Index', self)
        u_label = QLabel('Along Thickness:', self)
        v_label = QLabel('In-plane Lateral:', self)

        self.hu_line = QLineEdit('0')
        self.ku_line = QLineEdit('0')
        self.lu_line = QLineEdit('1')

        self.hv_line = QLineEdit('1')
        self.kv_line = QLineEdit('0')
        self.lv_line = QLineEdit('0')

        orient_layout = QGridLayout()

        orient_layout.addWidget(face_index_label, 0, 0, Qt.AlignCenter)
        orient_layout.addWidget(a_star_label, 0, 1, Qt.AlignCenter)
        orient_layout.addWidget(b_star_label, 0, 2, Qt.AlignCenter)
        orient_layout.addWidget(c_star_label, 0, 3, Qt.AlignCenter)

        orient_layout.addWidget(u_label, 1, 0)
        orient_layout.addWidget(self.hu_line, 1, 1)
        orient_layout.addWidget(self.ku_line, 1, 2)
        orient_layout.addWidget(self.lu_line, 1, 3)
        orient_layout.addWidget(v_label, 2, 0)
        orient_layout.addWidget(self.hv_line, 2, 1)
        orient_layout.addWidget(self.kv_line, 2, 2)
        orient_layout.addWidget(self.lv_line, 2, 3)

        stretch = QHeaderView.Stretch

        self.gon_table = QTableWidget()

        self.gon_table.setRowCount(3)
        self.gon_table.setColumnCount(6)

        labels = ['name','x','y','z','sense','angle']

        self.gon_table.horizontalHeader().setSectionResizeMode(stretch)
        self.gon_table.setHorizontalHeaderLabels(labels)
        self.gon_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.gon_table.setSelectionBehavior(QTableWidget.SelectRows)

        data = [['ω', 0, 1, 0, 1, 0.0],
                ['χ', 0, 0, 1, 1, 0.0],
                ['φ', 0, 1, 0, 1, 0.0]]

        for row, gon in enumerate(data):
            for col, val in enumerate(gon):
                item = QTableWidgetItem(str(val))
                item.setTextAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
                self.gon_table.setItem(row, col, item)

        goniometer_layout = QHBoxLayout()

        self.name_line = QLineEdit()
        self.x_line = QLineEdit()
        self.y_line = QLineEdit()
        self.z_line = QLineEdit()
        self.sense_line = QLineEdit()
        self.angle_line = QLineEdit()

        self.name_line.setReadOnly(True)

        validator = QIntValidator(-1, 1, self)

        self.x_line.setValidator(validator)
        self.y_line.setValidator(validator)
        self.z_line.setValidator(validator)

        regexp = QRegExp('^-1$|^1$')
        validator = QRegExpValidator(regexp)

        self.sense_line.setValidator(validator)

        validator = QDoubleValidator(-360, 360, 1, notation=notation)

        self.angle_line.setValidator(validator)

        goniometer_layout.addWidget(self.name_line)
        goniometer_layout.addWidget(self.x_line)
        goniometer_layout.addWidget(self.y_line)
        goniometer_layout.addWidget(self.z_line)
        goniometer_layout.addWidget(self.sense_line)
        goniometer_layout.addWidget(self.angle_line)

        sample_layout.addLayout(generate_layout)
        sample_layout.addLayout(param_layout)
        sample_layout.addWidget(self.gon_table)
        sample_layout.addLayout(goniometer_layout)
        sample_layout.addLayout(orient_layout)
        sample_layout.addLayout(material_layout)
        sample_layout.addLayout(cryst_layout)

        samp_tab.setLayout(sample_layout)

    def connect_sample_parameters(self, update_parameters):

        self.sample_combo.activated.connect(update_parameters)

        self.param1_line.editingFinished.connect(update_parameters)
        self.param2_line.editingFinished.connect(update_parameters)
        self.param3_line.editingFinished.connect(update_parameters)

    def connect_row_highligter(self, highlight_row):

        self.gon_table.itemSelectionChanged.connect(highlight_row)

    def connect_load_UB(self, load_UB):

        self.load_UB_button.clicked.connect(load_UB)

    def connect_goniometer_table(self, set_gonioneter_table):

        self.name_line.editingFinished.connect(set_gonioneter_table)
        self.x_line.editingFinished.connect(set_gonioneter_table)
        self.y_line.editingFinished.connect(set_gonioneter_table)
        self.z_line.editingFinished.connect(set_gonioneter_table)
        self.sense_line.editingFinished.connect(set_gonioneter_table)
        self.angle_line.editingFinished.connect(set_gonioneter_table)

    def connect_add_sample(self, add_sample):

        self.add_sample_button.clicked.connect(add_sample)

    def load_UB_file_dialog(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        filename, _ = QFileDialog.getOpenFileName(self,
                                                  'Load UB file',
                                                  '',
                                                  'UB files (*.mat)',
                                                  options=options)

        return filename

    def get_sample_shape(self):

        return self.sample_combo.currentText()

    def set_sample_constants(self, params):

        self.param1_line.setText('{:.2f}'.format(params[0]))
        self.param2_line.setText('{:.2f}'.format(params[1]))
        self.param3_line.setText('{:.2f}'.format(params[2]))

    def get_sample_constants(self):

        params = self.param1_line, self.param2_line, self.param3_line

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:

            return [float(param.text()) for param in params]

    def constrain_size(self, const):

        params = self.param1_line, self.param2_line, self.param3_line

        for fixed, param in zip(const, params):
            param.setDisabled(fixed)

    def set_goniometer(self, row, goniometer):

        for col, val in enumerate(goniometer):
            item = QTableWidgetItem(str(val))
            item.setTextAlignment(Qt.AlignHCenter | Qt.AlignVCenter)
            self.gon_table.setItem(row, col, item)

    def get_goniometer(self):

        row = self.gon_table.currentRow()
        if row is not None:
            return self.get_goniometer_angle(row)

    def set_angle(self, goniometer):

        self.name_line.setText(goniometer[0])
        self.x_line.setText(str(goniometer[1]))
        self.y_line.setText(str(goniometer[2]))
        self.z_line.setText(str(goniometer[3]))
        self.sense_line.setText(str(goniometer[4]))
        self.angle_line.setText(str(goniometer[5]))

    def get_goniometer_angle(self, row):

        name = self.gon_table.item(row, 0).text()
        x = self.gon_table.item(row, 1).text()
        y = self.gon_table.item(row, 2).text()
        z = self.gon_table.item(row, 3).text()
        sense = self.gon_table.item(row, 4).text()
        angle = self.gon_table.item(row, 5).text()
        axis = [int(val) for val in [x, y, z, sense]]
        goniometer = [name, *axis, float(angle)]

        return goniometer

    def get_goniometers(self):

        n = self.gon_table.rowCount()

        goniometers = []
        for row in range(n):
            goniometer = self.get_goniometer_angle(row)
            goniometers.append(goniometer)

        return goniometers

    def set_goniometer_table(self):

        row = self.gon_table.currentRow()

        params = self.name_line, self.x_line, self.y_line, self.z_line, \
                 self.sense_line, self.angle_line

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params and row:

            goniometer = [params[0].text(), \
                          *[int(param.text()) for param in params[1:-1]],
                          float(params[-1].text())]

            self.set_goniometer(row, goniometer)

    def set_unit_cell_volume(self, vol):

        self.V_line.setText('{:.4f}'.format(vol))

    def get_material_paremters(self):

        params = self.chem_line, self.Z_line, self.V_line,

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:

            vals = [params[0].text(),
                    float(params[1].text()),
                    float(params[2].text())]

            return vals

    def add_sample(self, sample_mesh):

        self.plotter.clear_actors()

        triangles = []
        for triangle in sample_mesh:
            triangles.append(pv.Triangle(triangle))

        multiblock = pv.MultiBlock(triangles)

        _, mapper = self.plotter.add_composite(multiblock,
                                               smooth_shading=True)

        self.plotter.add_legend_scale(corner_offset_factor=2,
                                      bottom_border_offset=50,
                                      top_border_offset=50,
                                      left_border_offset=100,
                                      right_border_offset=100,
                                      legend_visibility=True,
                                      xy_label_mode=False)

        self.plotter.add_axes_at_origin()

        self.reset_view()

    def set_absortion_parameters(self, abs_dict):

        self.sigma_a_line.setText('{:.4f}'.format(abs_dict['sigma_a']))
        self.sigma_s_line.setText('{:.4f}'.format(abs_dict['sigma_s']))

        self.mu_a_line.setText('{:.4f}'.format(abs_dict['mu_a']))
        self.mu_s_line.setText('{:.4f}'.format(abs_dict['mu_s']))

        self.N_line.setText('{:.4f}'.format(abs_dict['N']))
        self.M_line.setText('{:.4f}'.format(abs_dict['M']))
        self.n_line.setText('{:.4f}'.format(abs_dict['n']))
        self.rho_line.setText('{:.4f}'.format(abs_dict['rho']))
        self.v_line.setText('{:.4f}'.format(abs_dict['V']))
        self.m_line.setText('{:.4f}'.format(abs_dict['m']))

    def get_face_indexing(self):

        params = self.hu_line, self.ku_line, self.lu_line, \
                 self.hv_line, self.kv_line, self.lv_line

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:

            vals = [float(param.text()) for param in params]

            return vals[0:3], vals[3:6]