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

from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure

from garnet.config.atoms import colors, radii

class SampleView(QWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        sample_layout = self.__init_sample()
        viewer_layout = self.__init_viewer()
        peaks_layout = self.__init_peaks()

        vert_sep_left = QFrame()
        vert_sep_right = QFrame()

        vert_sep_left.setFrameShape(QFrame.VLine)
        vert_sep_right.setFrameShape(QFrame.VLine)

        layout = QHBoxLayout()

        layout.addLayout(sample_layout)
        layout.addWidget(vert_sep_left)
        layout.addLayout(viewer_layout)
        layout.addWidget(vert_sep_right)
        layout.addLayout(peaks_layout)

        self.setLayout(layout)

    def __init_sample(self):

        self.sample_combo = QComboBox(self)
        self.sample_combo.addItem('Sphere')
        self.sample_combo.addItem('Cylinder')
        self.sample_combo.addItem('Plate')

        self.add_sample_button = QPushButton('Add Sample', self)
        self.copy_material_button = QPushButton('Copy Material', self)

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0, 100, 5, notation=notation)

        param1_label = QLabel('Width', self)
        param2_label = QLabel('Height', self)
        param3_label = QLabel('Thickness', self)

        unit_label = QLabel('cm', self)

        self.param1_line = QLineEdit('0.05')
        self.param2_line = QLineEdit('0.05')
        self.param3_line = QLineEdit('0.05')

        self.param1_line.setValidator(validator)
        self.param2_line.setValidator(validator)
        self.param3_line.setValidator(validator)

        param_layout = QHBoxLayout()
        generate_layout = QHBoxLayout()

        generate_layout.addWidget(self.sample_combo)
        generate_layout.addWidget(self.copy_material_button)
        generate_layout.addWidget(self.add_sample_button)
    
        param_layout.addWidget(param1_label)
        param_layout.addWidget(self.param1_line)

        param_layout.addWidget(param2_label)
        param_layout.addWidget(self.param2_line)

        param_layout.addWidget(param3_label)
        param_layout.addWidget(self.param3_line)

        param_layout.addWidget(unit_label)

        material_layout = QGridLayout()

        self.chem_line = QLineEdit()
        self.Z_line = QLineEdit()
        self.V_line = QLineEdit()

        Z_label = QLabel('Z')
        V_label = QLabel('Ω')
        uc_vol_label = QLabel('Å^3')

        material_layout.addWidget(self.chem_line, 0, 0, 1, 5)
        material_layout.addWidget(Z_label, 1, 0)
        material_layout.addWidget(self.Z_line, 1, 1)
        material_layout.addWidget(V_label, 1, 2)
        material_layout.addWidget(self.V_line, 1, 3)
        material_layout.addWidget(uc_vol_label, 1, 4)

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
        V_label = QLabel('V', self)
        m_label = QLabel('m', self)

        N_unit_label = QLabel('atoms', self)
        M_unit_label = QLabel('g/mol', self)
        n_unit_label = QLabel('1/Å^3', self)
        rho_unit_label = QLabel('g/cm^3', self)
        V_unit_label = QLabel('cm^3', self)
        m_unit_label = QLabel('g', self)

        self.N_line = QLineEdit()
        self.M_line = QLineEdit()
        self.n_line = QLineEdit()
        self.rho_line = QLineEdit()
        self.V_line = QLineEdit()
        self.m_line = QLineEdit()

        self.N_line.setEnabled(False)
        self.M_line.setEnabled(False)
        self.n_line.setEnabled(False)
        self.rho_line.setEnabled(False)
        self.V_line.setEnabled(False)
        self.m_line.setEnabled(False)

        cryst_layout.addWidget(N_label, 3, 0)
        cryst_layout.addWidget(self.N_line, 3, 1, 1, 2)
        cryst_layout.addWidget(N_unit_label, 3, 3)

        cryst_layout.addWidget(M_label, 4, 0)
        cryst_layout.addWidget(self.M_line, 4, 1, 1, 2)
        cryst_layout.addWidget(M_unit_label, 4, 3)

        cryst_layout.addWidget(n_label, 5, 0)
        cryst_layout.addWidget(self.n_line, 5, 1, 1, 2)
        cryst_layout.addWidget(n_unit_label, 5, 3)

        cryst_layout.addWidget(rho_label, 6, 0)
        cryst_layout.addWidget(self.rho_line, 6, 1, 1, 2)
        cryst_layout.addWidget(rho_unit_label, 6, 3)

        cryst_layout.addWidget(V_label, 7, 0)
        cryst_layout.addWidget(self.V_line, 7, 1, 1, 2)
        cryst_layout.addWidget(V_unit_label, 7, 3)

        cryst_layout.addWidget(m_label, 8, 0)
        cryst_layout.addWidget(self.m_line, 8, 1, 1, 2)
        cryst_layout.addWidget(m_unit_label, 8, 3)

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

        sample_layout.addLayout(generate_layout)
        sample_layout.addLayout(param_layout)
        sample_layout.addStretch(1)
        sample_layout.addLayout(orient_layout)
        sample_layout.addStretch(1)
        sample_layout.addLayout(material_layout)
        sample_layout.addStretch(1)
        sample_layout.addLayout(cryst_layout)

        return sample_layout

    def __init_peaks(self):

        peaks_layout = QVBoxLayout()

        stretch = QHeaderView.Stretch

        self.peaks_table = QTableWidget()

        self.peaks_table.setRowCount(0)
        self.peaks_table.setColumnCount(6)

        self.peaks_table.horizontalHeader().setSectionResizeMode(stretch)
        self.peaks_table.setHorizontalHeaderLabels(['h','k','l','d','λ','T'])
        self.peaks_table.setEditTriggers(QTableWidget.NoEditTriggers)

        plot_layout = QVBoxLayout()

        self.canvas = FigureCanvas(Figure())

        plot_layout.addLayout(plot_layout)
        plot_layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        plot_layout.addWidget(self.canvas)

        peaks_layout.addWidget(self.peaks_table)
        peaks_layout.addLayout(plot_layout)

        return peaks_layout

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

    def change_proj(self):

        if self.proj_box.isChecked():
            self.plotter.enable_parallel_projection()
        else:
            self.plotter.disable_parallel_projection()

    def reset_view(self):

        self.plotter.reset_camera()
        self.plotter.view_isometric()

    def add_sample(self, atom_dict):

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