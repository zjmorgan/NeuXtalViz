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
                            QFrame,
                            QFileDialog)

from PyQt5.QtWidgets import QApplication, QMainWindow

from qtpy.QtGui import QDoubleValidator, QIntValidator
from PyQt5.QtCore import Qt

import pyvista as pv
from pyvistaqt import QtInteractor

import matplotlib.colors

from garnet.config.atoms import colors, radii

class CrystalStructureView(QWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        structure_layout = self.__init_structure()
        viewer_layout = self.__init_viewer()
        factors_layout = self.__init_factors()

        vert_sep_left = QFrame()
        vert_sep_right = QFrame()

        vert_sep_left.setFrameShape(QFrame.VLine)
        vert_sep_right.setFrameShape(QFrame.VLine)

        layout = QHBoxLayout()

        layout.addLayout(structure_layout)
        layout.addWidget(vert_sep_left)
        layout.addLayout(viewer_layout)
        layout.addWidget(vert_sep_right)
        layout.addLayout(factors_layout)

        self.setLayout(layout)

    def __init_structure(self):

        structure_layout = QVBoxLayout()

        crystal_layout = QHBoxLayout()
        parameters_layout = QGridLayout()

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

        a_label = QLabel('a')
        b_label = QLabel('b')
        c_label = QLabel('c')

        alpha_label = QLabel('α')
        beta_label = QLabel('β')
        gamma_label = QLabel('γ')

        angstrom_label = QLabel('Å')
        degree_label = QLabel('°')

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

        self.crystal_system_combo = QComboBox(self)
        self.crystal_system_combo.addItem('Triclinic')
        self.crystal_system_combo.addItem('Monoclinic')
        self.crystal_system_combo.addItem('Orthorhombic')
        self.crystal_system_combo.addItem('Tetragonal')
        self.crystal_system_combo.addItem('Trigonal')
        self.crystal_system_combo.addItem('Hexagonal')
        self.crystal_system_combo.addItem('Cubic')

        self.space_group_combo = QComboBox(self)
        self.setting_combo = QComboBox(self)

        self.crystal_system_combo.setEnabled(False)
        self.space_group_combo.setEnabled(False)
        self.setting_combo.setEnabled(False)

        self.load_CIF_button = QPushButton('Load CIF', self)

        crystal_layout.addWidget(self.crystal_system_combo)
        crystal_layout.addWidget(self.space_group_combo)
        crystal_layout.addWidget(self.setting_combo)
        crystal_layout.addWidget(self.load_CIF_button)

        structure_layout.addLayout(crystal_layout)
        structure_layout.addLayout(parameters_layout)

        stretch = QHeaderView.Stretch

        self.atm_table = QTableWidget()

        self.atm_table.setRowCount(0)
        self.atm_table.setColumnCount(6)

        self.atm_table.horizontalHeader().setSectionResizeMode(stretch)
        self.atm_table.setHorizontalHeaderLabels(['atm','x','y','z','occ','U'])
        self.atm_table.setEditTriggers(QTableWidget.NoEditTriggers)
        self.atm_table.setSelectionBehavior(QTableWidget.SelectRows)

        structure_layout.addWidget(self.atm_table)

        scatterer_layout = QHBoxLayout()

        self.atm_button = QPushButton('', self)

        self.x_line = QLineEdit()
        self.y_line = QLineEdit()
        self.z_line = QLineEdit()
        self.occ_line = QLineEdit()
        self.Uiso_line = QLineEdit()

        validator = QDoubleValidator(-1, 1, 4, notation=notation)

        self.x_line.setValidator(validator)
        self.y_line.setValidator(validator)
        self.z_line.setValidator(validator)

        validator = QDoubleValidator(0, 1, 4, notation=notation)

        self.occ_line.setValidator(validator)

        validator = QDoubleValidator(0, 100, 4, notation=notation)

        self.Uiso_line.setValidator(validator)

        scatterer_layout.addWidget(self.atm_button)
        scatterer_layout.addWidget(self.x_line)
        scatterer_layout.addWidget(self.y_line)
        scatterer_layout.addWidget(self.z_line)
        scatterer_layout.addWidget(self.occ_line)
        scatterer_layout.addWidget(self.Uiso_line)

        sample_layout = QGridLayout()

        self.chem_line = QLineEdit()
        self.Z_line = QLineEdit()
        self.V_line = QLineEdit()

        self.chem_line.setReadOnly(True)
        self.Z_line.setReadOnly(True)
        self.V_line.setReadOnly(True)

        Z_label = QLabel('Z')
        V_label = QLabel('Ω')
        uc_vol_label = QLabel('Å^3')

        sample_layout.addWidget(self.chem_line, 0, 0, 1, 5)
        sample_layout.addWidget(Z_label, 1, 0)
        sample_layout.addWidget(self.Z_line, 1, 1)
        sample_layout.addWidget(V_label, 1, 2)
        sample_layout.addWidget(self.V_line, 1, 3)
        sample_layout.addWidget(uc_vol_label, 1, 4)

        structure_layout.addLayout(scatterer_layout)
        structure_layout.addLayout(sample_layout)

        return structure_layout

    def __init_factors(self):

        factors_layout = QVBoxLayout()

        calculate_layout = QHBoxLayout()

        dmin_label = QLabel('d(min)')
        angstrom_label = QLabel('Å')

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0.1, 1000, 4, notation=notation)

        self.dmin_line = QLineEdit()
        self.dmin_line.setValidator(validator)

        self.calculate_button = QPushButton('Calculate', self)

        calculate_layout.addWidget(dmin_label)
        calculate_layout.addWidget(self.dmin_line)
        calculate_layout.addWidget(angstrom_label)
        calculate_layout.addStretch(1)
        calculate_layout.addWidget(self.calculate_button)

        stretch = QHeaderView.Stretch

        self.f2_table = QTableWidget()

        self.f2_table.setRowCount(0)
        self.f2_table.setColumnCount(5)

        self.f2_table.horizontalHeader().setSectionResizeMode(stretch)
        self.f2_table.setHorizontalHeaderLabels(['h','k','l','d','F²'])
        self.f2_table.setEditTriggers(QTableWidget.NoEditTriggers)

        indivdual_layout = QHBoxLayout()

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(-100, 100, 5, notation=notation)

        self.h_line = QLineEdit()
        self.k_line = QLineEdit()
        self.l_line = QLineEdit()

        self.h_line.setValidator(validator)
        self.k_line.setValidator(validator)
        self.l_line.setValidator(validator)

        self.individual_button = QPushButton('Calculate', self)

        hkl_label = QLabel('hkl')

        indivdual_layout.addWidget(hkl_label)
        indivdual_layout.addWidget(self.h_line)
        indivdual_layout.addWidget(self.k_line)
        indivdual_layout.addWidget(self.l_line)
        indivdual_layout.addStretch(1)
        indivdual_layout.addWidget(self.individual_button)

        factors_layout.addLayout(calculate_layout)
        factors_layout.addWidget(self.f2_table)
        factors_layout.addLayout(indivdual_layout)

        return factors_layout

    def __init_viewer(self):

        self.proj_box = QCheckBox('Parallel Projection', self)

        self.reset_button = QPushButton('Reset View', self)

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

        viewer_layout = QVBoxLayout()
        camera_layout = QGridLayout()

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

        viewer_layout.addLayout(camera_layout)
        viewer_layout.addWidget(self.plotter.interactor)

        return viewer_layout

    def set_transform(self, T):

        if T is not None:

            a = pv._vtk.vtkMatrix4x4()
            for i in range(3):
                for j in range(3):
                    a.SetElement(i,j,T[i,j])

            actor = self.plotter.add_axes(xlabel='a',
                                          ylabel='b',
                                          zlabel='c')
            actor.SetUserMatrix(a)

    def draw_cell(self, A):

        T = np.eye(4)
        T[:3,:3] = A

        mesh = pv.Box(bounds=(0,1,0,1,0,1), level=0, quads=True)
        mesh.transform(T, inplace=True)

        self.plotter.add_mesh(mesh,
                              color='k',
                              style='wireframe',
                              render_lines_as_tubes=True)

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

    def load_CIF_file_dialog(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        filename, _ = QFileDialog.getOpenFileName(self,
                                                  'Load CIF file',
                                                  '',
                                                  'CIf files (*.cif)',
                                                  options=options)

        return filename

    def get_crystal_system(self):

        return self.crystal_system_combo.currentText()

    def set_crystal_system(self, crystal_system):

        index = self.crystal_system_combo.findText(crystal_system)
        if index >= 0:
            self.crystal_system_combo.setCurrentIndex(index)

    def update_space_groups(self, nos):

        self.space_group_combo.clear()
        for no in nos:
            self.space_group_combo.addItem(no)

    def get_space_group(self):

        return self.space_group_combo.currentText()

    def set_space_group(self, space_group):

        index = self.space_group_combo.findText(space_group)
        if index >= 0:
            self.space_group_combo.setCurrentIndex(index)

    def update_settings(self, settings):

        self.setting_combo.clear()
        for setting in settings:
            self.setting_combo.addItem(setting)

    def get_setting(self):

        return self.setting_combo.currentText()

    def set_setting(self, setting):

        index = self.setting_combo.findText(setting)
        if index >= 0:
            self.setting_combo.setCurrentIndex(index)

    def set_lattice_constants(self, params):

        self.a_line.setText('{:.4f}'.format(params[0]))
        self.b_line.setText('{:.4f}'.format(params[1]))
        self.c_line.setText('{:.4f}'.format(params[2]))

        self.alpha_line.setText('{:.4f}'.format(params[3]))
        self.beta_line.setText('{:.4f}'.format(params[4]))
        self.gamma_line.setText('{:.4f}'.format(params[5]))

    def get_lattice_constants(self):

        params = self.a_line, self.b_line, self.c_line, \
                 self.alpha_line, self.beta_line, self.gamma_line

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:

            return [float(param.text()) for param in params]

    def set_unit_cell_volume(self, vol):

        self.V_line.setText('{:.4f}'.format(vol))

    def set_scatterers(self, scatterers):

        self.atm_table.setRowCount(0)
        self.atm_table.setRowCount(len(scatterers))

        for row, scatterer in enumerate(scatterers):
            self.set_scatterer(row, scatterer)

    def set_scatterer(self, row, scatterer):

        atm, *xyz, occ, Uiso = scatterer
        xyz = ['{:.4f}'.format(val) for val in xyz]
        occ = '{:.4f}'.format(occ)
        Uiso = '{:.4f}'.format(Uiso)
        self.atm_table.setItem(row, 0, QTableWidgetItem(atm))
        self.atm_table.setItem(row, 1, QTableWidgetItem(xyz[0]))
        self.atm_table.setItem(row, 2, QTableWidgetItem(xyz[1]))
        self.atm_table.setItem(row, 3, QTableWidgetItem(xyz[2]))
        self.atm_table.setItem(row, 4, QTableWidgetItem(occ))
        self.atm_table.setItem(row, 5, QTableWidgetItem(Uiso))

    def get_scatterer(self):

        row = self.atm_table.currentRow()
        if row is not None:
            return self.get_atom_site(row)

    def get_atom_site(self, row):

        atm = self.atm_table.item(row, 0).text()
        x = self.atm_table.item(row, 1).text()
        y = self.atm_table.item(row, 2).text()
        z = self.atm_table.item(row, 3).text()
        occ = self.atm_table.item(row, 4).text()
        Uiso = self.atm_table.item(row, 5).text()
        scatterer = [atm, *[float(val) for val in [x, y, z, occ, Uiso]]]

        return scatterer

    def get_scatterers(self):

        n = self.atm_table.getRowCount()

        scatterers = []
        for row in range(n):
            scatterer = self.get_atom_site(row)
            scatterers.append(scatterer)

        return scatterers

    def set_atom(self, scatterer):

        self.atm_button.setText(scatterer[0])
        self.x_line.setText(str(scatterer[1]))
        self.y_line.setText(str(scatterer[2]))
        self.z_line.setText(str(scatterer[3]))
        self.occ_line.setText(str(scatterer[4]))
        self.Uiso_line.setText(str(scatterer[5]))

    def set_atom_table(self):

        row = self.atm_table.currentRow()

        params = self.x_line, self.y_line, self.z_line, \
                 self.occ_line, self.Uiso_line

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params and row:

            scatterer = [self.atm_button.text(), \
                         *[float(param.text()) for param in params]]

            self.set_scatterer(row, scatterer)

    def set_formula_z(self, chemical_formula, z_parameter):

        self.chem_line.setText(chemical_formula)
        self.Z_line.setText(str(z_parameter))

    def get_minimum_d_spacing(self):

        if self.dmin_line.hasAcceptableInput():

            return float(self.dmin_line.text())

    def constrain_parameters(self, const):

        params = self.a_line, self.b_line, self.c_line, \
                 self.alpha_line, self.beta_line, self.gamma_line

        for fixed, param in zip(const, params):
            param.setDisabled(fixed)

    def change_proj(self):

        if self.proj_box.isChecked():
            self.plotter.enable_parallel_projection()
        else:
            self.plotter.disable_parallel_projection()

    def reset_view(self):

        self.plotter.reset_camera()
        self.plotter.view_isometric()

    def add_atoms(self, atom_dict):

        self.plotter.clear_actors()

        T = np.eye(4)

        geoms, cmap, self.indexing = [], [], {}

        sphere = pv.Icosphere(radius=1, nsub=1)

        atm_ind = 0

        for ind, atom in enumerate(atom_dict.keys()):

            color = colors[atom]
            radius = radii[atom][0]

            coordinates, opacities, indices = atom_dict[atom]

            for i_atm, (coord, occ) in enumerate(zip(coordinates, opacities)):
                T[0,0] = T[1,1] = T[2,2] = radius
                T[:3,3] = coord
                atm = sphere.copy().transform(T)
                atm['scalars'] = np.full(sphere.n_cells, ind+1.)
                geoms.append(atm)
                self.indexing[atm_ind] = atom
                atm_ind += 1

            cmap.append(color)

        cmap = matplotlib.colors.ListedColormap(cmap)

        multiblock = pv.MultiBlock(geoms)

        _, mapper = self.plotter.add_composite(multiblock,
                                               cmap=cmap,
                                               smooth_shading=True,
                                               show_scalar_bar=False)

        self.mapper = mapper

        self.plotter.enable_block_picking(callback=self.highlight,
                                          side='left')
        self.plotter.enable_block_picking(callback=self.highlight,
                                          side='right')

        self.change_proj()

    def highlight(self, index, dataset):

        color = self.mapper.block_attr[index].color

        if color == 'pink':
            color, select = None, False
        else:
            color, select = 'pink', True

        self.mapper.block_attr[index].color = color
        
        print('atom = {}'.format(self.indexing[index]))

        return self.indexing[index], select

    def set_factors(self, hkls, ds, F2s):

        self.f2_table.setRowCount(0)
        self.f2_table.setRowCount(len(hkls))

        for row, (hkl, d, F2) in enumerate(zip(hkls, ds, F2s)):
            hkl = ['{:.0f}'.format(val) for val in hkl]
            d = '{:.4f}'.format(d)
            F2 = '{:.2f}'.format(F2)
            self.f2_table.setItem(row, 0, QTableWidgetItem(hkl[0]))
            self.f2_table.setItem(row, 1, QTableWidgetItem(hkl[1]))
            self.f2_table.setItem(row, 2, QTableWidgetItem(hkl[2]))
            self.f2_table.setItem(row, 3, QTableWidgetItem(d))
            self.f2_table.setItem(row, 4, QTableWidgetItem(F2))

    def get_hkl(self):

        params = self.h_line, self.k_line, self.l_line

        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:

            return [float(param.text()) for param in params]

    def set_equivalents(self, hkls, d, F2):

        self.f2_table.setRowCount(0)
        self.f2_table.setRowCount(len(hkls))

        d = '{:.4f}'.format(d)
        F2 = '{:.2f}'.format(F2)

        for row, hkl in enumerate(hkls):
            hkl = ['{:.0f}'.format(val) for val in hkl]
            self.f2_table.setItem(row, 0, QTableWidgetItem(hkl[0]))
            self.f2_table.setItem(row, 1, QTableWidgetItem(hkl[1]))
            self.f2_table.setItem(row, 2, QTableWidgetItem(hkl[2]))
            self.f2_table.setItem(row, 3, QTableWidgetItem(d))
            self.f2_table.setItem(row, 4, QTableWidgetItem(F2))