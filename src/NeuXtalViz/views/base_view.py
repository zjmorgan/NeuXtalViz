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
                            QLineEdit)

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
        self.proj_box.clicked.connect(self.change_proj)

        self.reset_button = QPushButton('Reset View', self)
        self.reset_button.clicked.connect(self.reset_view)

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

        self.px_button = QPushButton('+Qx', self)
        self.py_button = QPushButton('+Qy', self)
        self.pz_button = QPushButton('+Qz', self)

        self.mx_button = QPushButton('-Qx', self)
        self.my_button = QPushButton('-Qy', self)
        self.mz_button = QPushButton('-Qz', self)

        self.a_star_button = QPushButton('a*', self)
        self.b_star_button = QPushButton('b*', self)
        self.c_star_button = QPushButton('c*', self)

        self.a_button = QPushButton('a', self)
        self.b_button = QPushButton('b', self)
        self.c_button = QPushButton('c', self)

        self.frame = QFrame()

        self.plotter = QtInteractor(self.frame)

        layout = QVBoxLayout()
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

        plot_layout.addWidget(self.plotter.interactor)

        layout.addLayout(camera_layout)
        layout.addLayout(plot_layout)

        self.setLayout(layout)

    def change_proj(self):
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

    def set_transform(self, T):
        """
        Apply the crystal axis transform to the axes.

        Parameters
        ----------
        T : 3x3 2d array
            Trasformation matrix.

        """        

        if T is not None:

            b = pv._vtk.vtkMatrix4x4()
            for i in range(3):
                for j in range(3):
                    b.SetElement(i,j,T[i,j])

            actor = self.plotter.add_axes(xlabel='a*',
                                          ylabel='b*',
                                          zlabel='c*')
            actor.SetUserMatrix(b)

    def view_vector(self, vecs):
        """
        Set the camera according to given vector(s).

        Parameters
        ----------
        vecs : list of 2 or single 3 element 1d array-like 
            Cameram direction and optional upward vector.

        """
        
        if len(vecs) == 2:
            vec = np.cross(vecs[0],vecs[1])
            self.plotter.view_vector(vecs[0],vec)
        else:
            self.plotter.view_vector(vecs)

    def update_axis_labels(self):
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

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        widget = NeuXtalVizWidget()
        self.setCentralWidget(widget)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())