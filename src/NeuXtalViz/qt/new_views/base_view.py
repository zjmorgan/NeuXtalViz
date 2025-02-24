from functools import partial

from qtpy.QtWidgets import (
    QWidget,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QVBoxLayout,
    QPushButton,
    QLabel,
    QCheckBox,
    QComboBox,
    QLineEdit,
    QProgressBar,
    QStatusBar,
    QTabWidget,
    QFileDialog,
)

from qtpy.QtGui import QDoubleValidator
from PyQt5.QtCore import Qt

import numpy as np
import pyvista as pv

from pyvistaqt import QtInteractor

from NeuXtalViz.qt.views.utilities import Worker, ThreadPool

# themes = {'Default': pv.themes.Theme(),
#           'Document': pv.themes.DocumentTheme(),
#           'Dark': pv.themes.DarkTheme(),
#           'ParaView': pv.themes.ParaViewTheme()}


class NeuXtalVizWidget(QWidget):
    def __init__(self, view_model, parent=None):
        super().__init__(parent)

        self.view_model = view_model

        self.proj_box = QCheckBox("Parallel Projection", self)
        self.proj_box.setChecked(True)

        self.reset_button = QPushButton("Reset View", self)
        self.reset_button.clicked.connect(self.reset_view)

        self.camera_button = QPushButton("Reset Camera", self)
        self.camera_button.clicked.connect(self.reset_camera)

        self.recip_box = QCheckBox("Reciprocal Lattice", self)
        self.recip_box.setChecked(True)

        self.axes_box = QCheckBox("Show Axes", self)
        self.axes_box.setChecked(True)

        self.save_button = QPushButton("Save Screenshot", self)

        self.frame = QFrame()

        self.plotter = QtInteractor(self.frame)

        layout = QHBoxLayout()
        vis_layout = QVBoxLayout()

        camera_layout = QHBoxLayout()
        left_layout = QVBoxLayout()
        right_layout = QVBoxLayout()

        # left_layout.addStretch(1)
        left_layout.addWidget(self.save_button)
        left_layout.addWidget(self.reset_button)
        left_layout.addWidget(self.camera_button)

        # right_layout.addStretch(1)
        right_layout.addWidget(self.recip_box)
        right_layout.addWidget(self.axes_box)
        right_layout.addWidget(self.proj_box)

        view_tab = self.__init_view_tab()

        camera_layout.addLayout(left_layout)
        camera_layout.addWidget(view_tab)
        camera_layout.addLayout(right_layout)

        vis_layout.addLayout(camera_layout)
        vis_layout.addWidget(self.plotter.interactor)
        vis_layout.setStretch(1, 1)

        info_tab = self.__init_info_tab()

        vis_layout.addWidget(info_tab)

        self.status_bar = QStatusBar()
        self.status_bar.showMessage("Ready!")
        self.progress_bar = QProgressBar()
        self.status_bar.addPermanentWidget(self.progress_bar)

        vis_layout.addWidget(self.status_bar)

        layout.addLayout(vis_layout, stretch=1)

        self.setLayout(layout)

        self.camera_position = None
        self.T = None

        self.threadpool = ThreadPool()

        self.plotter.enable_parallel_projection()

    def __init_view_tab(self):
        view_tab = QTabWidget()

        self.view_combo = QComboBox(self)
        self.view_combo.addItem("[hkl]")
        self.view_combo.addItem("[uvw]")

        self.viewup_combo = QComboBox(self)
        self.viewup_combo.addItem("[hkl]")
        self.viewup_combo.addItem("[uvw]")

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(-100, 100, 5, notation=notation)

        self.axis1_line = QLineEdit()
        self.axis2_line = QLineEdit()
        self.axis3_line = QLineEdit()

        self.axis1_line.setValidator(validator)
        self.axis2_line.setValidator(validator)
        self.axis3_line.setValidator(validator)

        self.axis1_label = QLabel("h", self)
        self.axis2_label = QLabel("k", self)
        self.axis3_label = QLabel("l", self)

        self.axisup1_line = QLineEdit()
        self.axisup2_line = QLineEdit()
        self.axisup3_line = QLineEdit()

        self.axisup1_line.setValidator(validator)
        self.axisup2_line.setValidator(validator)
        self.axisup3_line.setValidator(validator)

        self.axisup1_label = QLabel("h", self)
        self.axisup2_label = QLabel("k", self)
        self.axisup3_label = QLabel("l", self)

        self.manual_button = QPushButton("View Axis", self)
        self.manualup_button = QPushButton("View Up Axis", self)

        self.px_button = QPushButton("+Qx", self)
        self.py_button = QPushButton("+Qy", self)
        self.pz_button = QPushButton("+Qz", self)

        self.mx_button = QPushButton("-Qx", self)
        self.my_button = QPushButton("-Qy", self)
        self.mz_button = QPushButton("-Qz", self)

        self.px_button.clicked.connect(self.view_yz)
        self.py_button.clicked.connect(self.view_zx)
        self.pz_button.clicked.connect(self.view_xy)

        self.mx_button.clicked.connect(self.view_zy)
        self.my_button.clicked.connect(self.view_xz)
        self.mz_button.clicked.connect(self.view_yx)

        self.a_star_button = QPushButton("a*", self)
        self.b_star_button = QPushButton("b*", self)
        self.c_star_button = QPushButton("c*", self)

        self.a_button = QPushButton("a", self)
        self.b_button = QPushButton("b", self)
        self.c_button = QPushButton("c", self)

        directions_layout = QGridLayout()
        manual_layout = QGridLayout()

        directions_tab = QWidget()
        manual_tab = QWidget()

        directions_layout.addWidget(self.px_button, 0, 0)
        directions_layout.addWidget(self.py_button, 0, 1)
        directions_layout.addWidget(self.pz_button, 0, 2)
        directions_layout.addWidget(self.a_star_button, 0, 3)
        directions_layout.addWidget(self.b_star_button, 0, 4)
        directions_layout.addWidget(self.c_star_button, 0, 5)

        directions_layout.addWidget(self.mx_button, 1, 0)
        directions_layout.addWidget(self.my_button, 1, 1)
        directions_layout.addWidget(self.mz_button, 1, 2)
        directions_layout.addWidget(self.a_button, 1, 3)
        directions_layout.addWidget(self.b_button, 1, 4)
        directions_layout.addWidget(self.c_button, 1, 5)

        manual_layout.addWidget(self.axis1_label, 0, 0, Qt.AlignCenter)
        manual_layout.addWidget(self.axis2_label, 0, 1, Qt.AlignCenter)
        manual_layout.addWidget(self.axis3_label, 0, 2, Qt.AlignCenter)

        manual_layout.addWidget(self.axis1_line, 1, 0)
        manual_layout.addWidget(self.axis2_line, 1, 1)
        manual_layout.addWidget(self.axis3_line, 1, 2)

        manual_layout.addWidget(self.view_combo, 0, 3)
        manual_layout.addWidget(self.manual_button, 1, 3)

        manual_layout.addWidget(self.axisup1_label, 0, 4, Qt.AlignCenter)
        manual_layout.addWidget(self.axisup2_label, 0, 5, Qt.AlignCenter)
        manual_layout.addWidget(self.axisup3_label, 0, 6, Qt.AlignCenter)

        manual_layout.addWidget(self.axisup1_line, 1, 4)
        manual_layout.addWidget(self.axisup2_line, 1, 5)
        manual_layout.addWidget(self.axisup3_line, 1, 6)

        manual_layout.addWidget(self.viewup_combo, 0, 7)
        manual_layout.addWidget(self.manualup_button, 1, 7)

        directions_tab.setLayout(directions_layout)
        manual_tab.setLayout(manual_layout)

        view_tab.addTab(directions_tab, "Direction View")
        view_tab.addTab(manual_tab, "Manual View")

        return view_tab

    def __init_info_tab(self):
        info_tab = QTabWidget()

        ub_a_label = QLabel("a:", self)
        ub_b_label = QLabel("b:", self)
        ub_c_label = QLabel("c:", self)
        ub_alpha_label = QLabel("α:", self)
        ub_beta_label = QLabel("β:", self)
        ub_gamma_label = QLabel("γ:", self)
        ub_u_label = QLabel("u:", self)
        ub_v_label = QLabel("v:", self)

        ub_angstrom_label = QLabel("Å")
        ub_degree_label = QLabel("°")

        self.ub_a_line = QLineEdit()
        self.ub_b_line = QLineEdit()
        self.ub_c_line = QLineEdit()
        self.ub_alpha_line = QLineEdit()
        self.ub_beta_line = QLineEdit()
        self.ub_gamma_line = QLineEdit()
        self.ub_u1_line = QLineEdit()
        self.ub_u2_line = QLineEdit()
        self.ub_u3_line = QLineEdit()
        self.ub_v1_line = QLineEdit()
        self.ub_v2_line = QLineEdit()
        self.ub_v3_line = QLineEdit()

        self.ub_a_line.setReadOnly(True)
        self.ub_b_line.setReadOnly(True)
        self.ub_c_line.setReadOnly(True)
        self.ub_alpha_line.setReadOnly(True)
        self.ub_beta_line.setReadOnly(True)
        self.ub_gamma_line.setReadOnly(True)
        self.ub_u1_line.setReadOnly(True)
        self.ub_u2_line.setReadOnly(True)
        self.ub_u3_line.setReadOnly(True)
        self.ub_v1_line.setReadOnly(True)
        self.ub_v2_line.setReadOnly(True)
        self.ub_v3_line.setReadOnly(True)

        lattice_layout = QGridLayout()
        orientation_layout = QGridLayout()

        lattice_tab = QWidget()
        orientation_tab = QWidget()

        lattice_layout.addWidget(ub_a_label, 0, 0)
        lattice_layout.addWidget(self.ub_a_line, 0, 1)
        lattice_layout.addWidget(ub_b_label, 0, 2)
        lattice_layout.addWidget(self.ub_b_line, 0, 3)
        lattice_layout.addWidget(ub_c_label, 0, 4)
        lattice_layout.addWidget(self.ub_c_line, 0, 5)
        lattice_layout.addWidget(ub_angstrom_label, 0, 6)

        lattice_layout.addWidget(ub_alpha_label, 1, 0)
        lattice_layout.addWidget(self.ub_alpha_line, 1, 1)
        lattice_layout.addWidget(ub_beta_label, 1, 2)
        lattice_layout.addWidget(self.ub_beta_line, 1, 3)
        lattice_layout.addWidget(ub_gamma_label, 1, 4)
        lattice_layout.addWidget(self.ub_gamma_line, 1, 5)
        lattice_layout.addWidget(ub_degree_label, 1, 6)

        orientation_layout.addWidget(ub_u_label, 0, 0)
        orientation_layout.addWidget(self.ub_u1_line, 0, 1)
        orientation_layout.addWidget(self.ub_u2_line, 0, 2)
        orientation_layout.addWidget(self.ub_u3_line, 0, 3)
        orientation_layout.addWidget(ub_v_label, 1, 0)
        orientation_layout.addWidget(self.ub_v1_line, 1, 1)
        orientation_layout.addWidget(self.ub_v2_line, 1, 2)
        orientation_layout.addWidget(self.ub_v3_line, 1, 3)

        lattice_tab.setLayout(lattice_layout)
        orientation_tab.setLayout(orientation_layout)

        info_tab.addTab(lattice_tab, "Lattice Parameters")
        info_tab.addTab(orientation_tab, "Sample Orientation")

        return info_tab

    def connect_bindings(self):
        self.view_model.show_axes_bind.connect("show_axes", self.show_axes)
        self.view_model.parallel_projection_bind.connect(
            "parallel_projection", self.change_projection
        )

        self.view_model.lattice_parameters_bind.connect(
            "oriented_lattice", self.set_oriented_lattice_parameters
        )
        self.view_model.progress_bind.connect("progress", self.set_step)
        self.view_model.status_bind.connect("status", self.set_info)
        self.view_model.up_vector_bind.connect("up_vector", self.view_up_vector)
        self.view_model.update_labels_bind.connect("update_labels", self.update_labels)
        self.view_model.vector_bind.connect("vector", self.view_vector)

    def connect_widgets(self):
        self.view_combo.currentIndexChanged.connect(self.view_model.update_axis_type)
        self.axis1_line.textChanged.connect(partial(self.view_model.update_manual_axis, 0))
        self.axis2_line.textChanged.connect(partial(self.view_model.update_manual_axis, 1))
        self.axis3_line.textChanged.connect(partial(self.view_model.update_manual_axis, 2))
        self.viewup_combo.currentIndexChanged.connect(self.view_model.update_up_axis_type)
        self.axisup1_line.textChanged.connect(
            partial(self.view_model.update_manual_up_axis, 0)
        )
        self.axisup2_line.textChanged.connect(
            partial(self.view_model.update_manual_up_axis, 1)
        )
        self.axisup3_line.textChanged.connect(
            partial(self.view_model.update_manual_up_axis, 2)
        )
        self.manual_button.clicked.connect(self.view_model.view_manual)
        self.manualup_button.clicked.connect(self.view_model.view_up_manual)
        self.a_star_button.clicked.connect(self.view_model.view_bc_star)
        self.b_star_button.clicked.connect(self.view_model.view_ca_star)
        self.c_star_button.clicked.connect(self.view_model.view_ab_star)
        self.a_button.clicked.connect(self.view_model.view_bc)
        self.b_button.clicked.connect(self.view_model.view_ca)
        self.c_button.clicked.connect(self.view_model.view_ab)
        self.save_button.clicked.connect(self.save_screenshot)
        self.recip_box.clicked.connect(self.view_model.change_lattice)
        self.axes_box.clicked.connect(self.view_model.change_axes)
        self.proj_box.clicked.connect(self.view_model.change_projection)

    def start_worker_pool(self, worker):
        """
        Create a worker pool.

        """

        self.threadpool.start_worker_pool(worker)

    def worker(self, task):
        """
        Worker task.

        """

        return Worker(task)

    def set_info(self, status):
        """
        Update status information.

        Parameters
        ----------
        status : str
            Information.

        """

        self.status_bar.showMessage(status)

    def set_step(self, progress):
        """
        Update progress step.

        Parameters
        ----------
        progress : int
            Step.

        """

        self.progress_bar.setValue(progress)

    def set_oriented_lattice_parameters(self, ol):
        """
        Update the oriented lattice paramters.

        Parameters (destructured from ol list)
        --------------------------------------
        a, b, c : float
            Lattice constants.
        alpha, beta, gamma : float
            Lattice angles.
        """
        self.ub_a_line.setText("{:.5f}".format(ol.a))
        self.ub_b_line.setText("{:.5f}".format(ol.b))
        self.ub_c_line.setText("{:.5f}".format(ol.c))
        self.ub_alpha_line.setText("{:.3f}".format(ol.alpha))
        self.ub_beta_line.setText("{:.3f}".format(ol.beta))
        self.ub_gamma_line.setText("{:.3f}".format(ol.gamma))
        self.ub_u1_line.setText("{:.4f}".format(ol.u[0]))
        self.ub_u2_line.setText("{:.4f}".format(ol.u[1]))
        self.ub_u3_line.setText("{:.4f}".format(ol.u[2]))
        self.ub_v1_line.setText("{:.4f}".format(ol.v[0]))
        self.ub_v2_line.setText("{:.4f}".format(ol.v[1]))
        self.ub_v3_line.setText("{:.4f}".format(ol.v[2]))

    def change_projection(self):
        """
        Enable or disable parallel projection.

        """

        if self.proj_box.isChecked():
            self.plotter.enable_parallel_projection()
        else:
            self.plotter.disable_parallel_projection()

    def reset_view(self, negative=False):
        """
        Reset the view.

        """

        self.plotter.reset_camera()
        self.plotter.view_isometric(negative)
        self.camera_position = self.plotter.camera_position

    def reset_camera(self):
        """
        Reset the camera.

        """

        self.plotter.reset_camera()

    def clear_scene(self):
        self.plotter.clear_plane_widgets()
        self.plotter.clear_actors()

        if self.camera_position is not None:
            self.camera_position = self.plotter.camera_position

    def reset_scene(self):
        if self.camera_position is not None:
            self.plotter.camera_position = self.camera_position
        else:
            self.reset_view()

    def save_screenshot(self, filename):
        """
        Save plotter screenshot.

        Parameters
        ----------
        filename : str
            Filename with *.png extension.

        """

        filename = self.save_screenshot_file_dialog()

        if filename:
            self.plotter.screenshot(filename)

    def save_screenshot_file_dialog(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filename, _ = file_dialog.getSaveFileName(
            self, "Save PNG file", "", "PNG files (*.png)", options=options
        )

        return filename

    def show_axes(self, data):
        T, reciprocal_lattice, show_axes = data

        if not show_axes:
            self.plotter.hide_axes()
        elif T is not None:
            t = pv._vtk.vtkMatrix4x4()
            for i in range(3):
                for j in range(3):
                    t.SetElement(i, j, T[i, j])
            if reciprocal_lattice:
                actor = self.plotter.add_axes(xlabel="a*", ylabel="b*", zlabel="c*")
            else:
                actor = self.plotter.add_axes(xlabel="a", ylabel="b", zlabel="c")
            actor.SetUserMatrix(t)

    def view_vector(self, vecs):
        """
        Set the camera according to given vector(s).

        Parameters
        ----------
        vecs : list of 2 or single 3 element 1d array-like
            Camera direction and optional upward vector.

        """

        if len(vecs) == 2:
            vec = np.cross(vecs[0], vecs[1])
            self.plotter.view_vector(vecs[0], vec)
        else:
            self.plotter.view_vector(vecs)

    def view_up_vector(self, vec):
        """
        Set the camera according to given vector(s).

        Parameters
        ----------
        vec : 3 element 1d array-like
            Camera up direction and optional upward vector.

        """

        self.plotter.set_viewup(vec)

    def update_labels(self, _):
        """
        Change the axes labels between Miller and fractional notation.

        """

        axes_type = self.view_combo.currentText()

        if axes_type == "[hkl]":
            self.axis1_label.setText("h")
            self.axis2_label.setText("k")
            self.axis3_label.setText("l")
        else:
            self.axis1_label.setText("u")
            self.axis2_label.setText("v")
            self.axis3_label.setText("w")

        axesup_type = self.viewup_combo.currentText()

        if axesup_type == "[hkl]":
            self.axisup1_label.setText("h")
            self.axisup2_label.setText("k")
            self.axisup3_label.setText("l")
        else:
            self.axisup1_label.setText("u")
            self.axisup2_label.setText("v")
            self.axisup3_label.setText("w")

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

            ind = np.array([axis1, axis2, axis3])

            return axes_type, ind

    def get_manual_axis_up_indices(self):
        """
        Indices of manually entered direction up components.

        Returns
        -------
        axes_type : str, [hkl] or [uvw]
            Miller index or fractional coordinate.
        ind : 3-element 1d array-like
            Indices.

        """

        axes_type = self.viewup_combo.currentText()

        axes = [self.axisup1_line, self.axisup2_line, self.axisup3_line]
        valid_axes = all([axis.hasAcceptableInput() for axis in axes])

        if valid_axes:
            axis1 = float(self.axisup1_line.text())
            axis2 = float(self.axisup2_line.text())
            axis3 = float(self.axisup3_line.text())

            ind = np.array([axis1, axis2, axis3])

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

    def set_position(self, pos):
        """
        Set the position.

        Parameters
        ----------
        pos : 3-element 1d array-like
            Coordinate position.

        """

        self.plotter.set_position(pos, reset=True)
