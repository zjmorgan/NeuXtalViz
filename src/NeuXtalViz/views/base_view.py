import sys

from qtpy.QtWidgets import (QWidget,
                            QFrame,
                            QGridLayout,
                            QHBoxLayout,
                            QVBoxLayout,
                            QPushButton,
                            QLabel,
                            QCheckBox,
                            QComboBox,
                            QLineEdit,
                            QFileDialog)

from qtpy.QtGui import QDoubleValidator
from PyQt5.QtCore import Qt

from PyQt5.QtWidgets import QApplication, QMainWindow

import numpy as np
import pyvista as pv
from pyvistaqt import QtInteractor

class NeuXtalVizWidget(QWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.proj_box = QCheckBox('Parallel Projection', self)
        self.proj_box.clicked.connect(self.change_projection)

        self.reset_button = QPushButton('Reset View', self)
        self.reset_button.clicked.connect(self.reset_view)

        self.view_combo = QComboBox(self)
        self.view_combo.addItem('[hkl]')
        self.view_combo.addItem('[uvw]')
        self.view_combo.currentIndexChanged.connect(self.update_labels)

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

        self.px_button = QPushButton('+Qx', self)
        self.py_button = QPushButton('+Qy', self)
        self.pz_button = QPushButton('+Qz', self)

        self.mx_button = QPushButton('-Qx', self)
        self.my_button = QPushButton('-Qy', self)
        self.mz_button = QPushButton('-Qz', self)

        self.px_button.clicked.connect(self.view_yz)
        self.py_button.clicked.connect(self.view_zx)
        self.pz_button.clicked.connect(self.view_xy)

        self.mx_button.clicked.connect(self.view_zy)
        self.my_button.clicked.connect(self.view_xz)
        self.mz_button.clicked.connect(self.view_yx)

        self.a_star_button = QPushButton('a*', self)
        self.b_star_button = QPushButton('b*', self)
        self.c_star_button = QPushButton('c*', self)

        self.a_button = QPushButton('a', self)
        self.b_button = QPushButton('b', self)
        self.c_button = QPushButton('c', self)

        self.recip_box = QCheckBox('Reciprocal Lattice', self)
        self.recip_box.setChecked(True)

        self.save_button = QPushButton('Save Screenshot', self)

        self.frame = QFrame()

        self.plotter = QtInteractor(self.frame)

        layout = QHBoxLayout()
        vis_layout = QVBoxLayout()
        camera_layout = QGridLayout()
        plot_layout = QHBoxLayout()

        camera_layout.addWidget(self.proj_box, 0, 0)
        camera_layout.addWidget(self.reset_button, 1, 0)

        camera_layout.addWidget(self.px_button, 0, 1)
        camera_layout.addWidget(self.py_button, 0, 2)
        camera_layout.addWidget(self.pz_button, 0, 3)
        camera_layout.addWidget(self.a_star_button, 0, 4)
        camera_layout.addWidget(self.b_star_button, 0, 5)
        camera_layout.addWidget(self.c_star_button, 0, 6)

        camera_layout.addWidget(self.mx_button, 1, 1)
        camera_layout.addWidget(self.my_button, 1, 2)
        camera_layout.addWidget(self.mz_button, 1, 3)
        camera_layout.addWidget(self.a_button, 1, 4)
        camera_layout.addWidget(self.b_button, 1, 5)
        camera_layout.addWidget(self.c_button, 1, 6)

        camera_layout.addWidget(self.axis1_label, 0, 7, Qt.AlignCenter)
        camera_layout.addWidget(self.axis2_label, 0, 8, Qt.AlignCenter)
        camera_layout.addWidget(self.axis3_label, 0, 9, Qt.AlignCenter)

        camera_layout.addWidget(self.axis1_line, 1, 7)
        camera_layout.addWidget(self.axis2_line, 1, 8)
        camera_layout.addWidget(self.axis3_line, 1, 9)

        camera_layout.addWidget(self.view_combo, 0, 10)
        camera_layout.addWidget(self.manual_button, 1, 10)

        camera_layout.addWidget(self.recip_box, 0, 11)
        camera_layout.addWidget(self.save_button, 1, 11)

        plot_layout.addWidget(self.plotter.interactor)

        vis_layout.addLayout(camera_layout)
        vis_layout.addLayout(plot_layout)

        layout.addLayout(vis_layout)

        self.setLayout(layout)

    def connect_manual_axis(self, view_manual):
        """
        Manual axis view connection.

        Parameters
        ----------
        view_manual : function
            Manual axis view handler.

        """

        self.manual_button.clicked.connect(view_manual)

    def connect_reciprocal_axes(self, view_a_star, view_b_star, view_c_star):
        """
        Reciprocal axes view connections.

        Parameters
        ----------
        view_a_star : function
            :math:`a^\ast`-axis view handler.
        view_b_star : function
            :math:`b^\ast`-axis view handler.
        view_c_star : function
            :math:`c^\ast`-axis view handler.

        """

        self.a_star_button.clicked.connect(view_a_star)
        self.b_star_button.clicked.connect(view_b_star)
        self.c_star_button.clicked.connect(view_c_star)

    def connect_real_axes(self, view_a, view_b, view_c):
        """
        Real axes view connections.

        Parameters
        ----------
        view_a : function
            :math:`a`-axis view handler.
        view_b : function
            :math:`b`-axis view handler.
        view_c : function
            :math:`c`-axis view handler.

        """

        self.a_button.clicked.connect(view_a)
        self.b_button.clicked.connect(view_b)
        self.c_button.clicked.connect(view_c)

    def connect_save_screenshot(self, save_screenshot):
        """
        Screenshot connection.

        Parameters
        ----------
        save_screenshot : function
            Screenshot handler.

        """

        self.save_button.clicked.connect(save_screenshot)

    def connect_reciprocal_real_compass(self, change_lattice):
        """
        Reciprocal/real axis compass

        Parameters
        ----------
        change_lattice : function
            Lattice handler.

        """

        self.recip_box.clicked.connect(change_lattice)

    def change_projection(self):
        """
        Enable or disable parallel projection.

        """

        if self.proj_box.isChecked():
            self.plotter.enable_parallel_projection()
        else:
            self.plotter.disable_parallel_projection()

    def reset_view(self):
        """
        Reset the camera view.

        """

        self.plotter.reset_camera()
        self.plotter.view_isometric()

    def save_screenshot(self, filename):
        """
        Save plotter screenshot.

        Parameters
        ----------
        filename : str
            Filename with .png extension.

        """

        self.plotter.screenshot(filename)

    def save_screenshot_file_dialog(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        filename, _ = QFileDialog.getSaveFileName(self,
                                                  'Save PNG file',
                                                  '',
                                                  'PNG files (*.png)',
                                                  options=options)

        return filename

    def set_transform(self, T):
        """
        Apply a transform to the axes.

        Parameters
        ----------
        T : 3x3 2d array
            Trasformation matrix.

        """

        if T is not None:

            t = pv._vtk.vtkMatrix4x4()

            for i in range(3):
                for j in range(3):
                    t.SetElement(i,j,T[i,j])

            if self.reciprocal_lattice():

                actor = self.plotter.add_axes(xlabel='a*',
                                              ylabel='b*',
                                              zlabel='c*')

            else:

                actor = self.plotter.add_axes(xlabel='a',
                                              ylabel='b',
                                              zlabel='c')

            actor.SetUserMatrix(t)

    def reciprocal_lattice(self):

        return self.recip_box.isChecked()

    def view_vector(self, vecs):
        """
        Set the camera according to given vector(s).

        Parameters
        ----------
        vecs : list of 2 or single 3 element 1d array-like
            Cameram direction and optional upward vector.

        """

        if len(vecs) == 2:
            vec = np.cross(vecs[0], vecs[1])
            self.plotter.view_vector(vecs[0], vec)
        else:
            self.plotter.view_vector(vecs)

    def update_labels(self):
        """
        Change the axes labels between Miller and fractional notation.

        """

        axes_type = self.view_combo.currentText()

        if axes_type == '[hkl]':
            self.axis1_label.setText('h')
            self.axis2_label.setText('k')
            self.axis3_label.setText('l')
        else:
            self.axis1_label.setText('u')
            self.axis2_label.setText('v')
            self.axis3_label.setText('w')

    def get_manual_axis_indices(self):
        """
        Indices of manually entered direction components.

        Returns
        -------
        axes_type : str, [hkl] or [uvw]
            Miller index or fractional coordinate.
        ind : 3-element 1d array-like
            Indices.

        """

        axes_type = self.view_combo.currentText()

        axes = [self.axis1_line, self.axis2_line, self.axis3_line]
        valid_axes = all([axis.hasAcceptableInput() for axis in axes])

        if valid_axes:

            axis1 = float(self.axis1_line.text())
            axis2 = float(self.axis2_line.text())
            axis3 = float(self.axis3_line.text())

            ind = np.array([axis1,axis2,axis3])

            return axes_type, ind

    def view_xy(self):
        """
        View :math:`xy`-plane.

        """

        self.plotter.view_xy()

    def view_yz(self):
        """
        View :math:`yz`-plane.

        """

        self.plotter.view_yz()

    def view_zx(self):
        """
        View :math:`zx`-plane.

        """

        self.plotter.view_zx()

    def view_yx(self):
        """
        View :math:`yx`-plane.

        """

        self.plotter.view_yx()

    def view_zy(self):
        """
        View :math:`zy`-plane.

        """

        self.plotter.view_zy()

    def view_xz(self):
        """
        View :math:`xz`-plane.

        """

        self.plotter.view_xz()