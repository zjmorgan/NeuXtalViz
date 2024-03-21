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
                            QLineEdit,
                            QFileDialog)

from qtpy.QtGui import QDoubleValidator, QIntValidator

import numpy as np
import pyvista as pv

from pyvistaqt import QtInteractor

# import matplotlib as mpl

class SatellitePeakIndexerView(QWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.proj_box = QCheckBox('Parallel Projection', self)
        self.proj_box.clicked.connect(self.change_proj)

        self.reset_button = QPushButton('Reset View', self)
        self.reset_button.clicked.connect(self.reset_view)

        self.view_combo = QComboBox(self)
        self.view_combo.addItem('[uvw]')
        self.view_combo.addItem('[hkl]')

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(-100, 100, 5, notation=notation)

        self.axis1_line = QLineEdit()
        self.axis2_line = QLineEdit()
        self.axis3_line = QLineEdit()

        self.axis1_line.setValidator(validator)
        self.axis2_line.setValidator(validator)
        self.axis3_line.setValidator(validator)

        self.manual_button = QPushButton('View Plane', self)

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

        frame = QFrame()

        self.plotter = QtInteractor(frame)
        self.plotter.add_camera_orientation_widget()

        viewer_layout = QVBoxLayout()
        camera_layout = QHBoxLayout()
        plot_layout = QHBoxLayout()
        view_layout = QGridLayout()

        camera_layout.addWidget(self.proj_box)
        camera_layout.addWidget(self.reset_button)
        camera_layout.addStretch(1)
        camera_layout.addWidget(self.axis1_line)
        camera_layout.addWidget(self.axis2_line)
        camera_layout.addWidget(self.axis3_line)
        camera_layout.addWidget(self.view_combo)
        camera_layout.addWidget(self.manual_button)

        plot_layout.addWidget(self.plotter.interactor)

        view_layout.addWidget(self.px_button, 0, 0)
        view_layout.addWidget(self.py_button, 0, 1)
        view_layout.addWidget(self.pz_button, 0, 2)
        view_layout.addWidget(self.a_star_button, 0, 3)
        view_layout.addWidget(self.b_star_button, 0, 4)
        view_layout.addWidget(self.c_star_button, 0, 5)

        view_layout.addWidget(self.mx_button, 1, 0)
        view_layout.addWidget(self.my_button, 1, 1)
        view_layout.addWidget(self.mz_button, 1, 2)
        view_layout.addWidget(self.a_button, 1, 3)
        view_layout.addWidget(self.b_button, 1, 4)
        view_layout.addWidget(self.c_button, 1, 5)

        viewer_layout.addLayout(camera_layout)
        viewer_layout.addLayout(plot_layout)
        viewer_layout.addLayout(view_layout)

        self.cluster_button = QPushButton('Cluster', self)

        dbl_validator = QDoubleValidator(0.0, 10, 5, notation=notation)
        int_validator = QIntValidator(1, 1000)

        self.param_eps_line = QLineEdit('0.025')
        self.param_min_line = QLineEdit('15')

        self.param_eps_line.setValidator(dbl_validator)
        self.param_min_line.setValidator(int_validator)

        self.table = QTableWidget()

        self.table.setRowCount(0)
        self.table.setColumnCount(3)

        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setHorizontalHeaderLabels(['h','k','l'])

        cluster_layout = QFormLayout()

        cluster_layout.addRow(self.cluster_button)
        cluster_layout.addRow('Maximum distance:', self.param_eps_line)
        cluster_layout.addRow('Minimum samples:', self.param_min_line)
        cluster_layout.addRow(self.table)

        layout = QHBoxLayout()

        layout.addLayout(viewer_layout)
        layout.addLayout(cluster_layout)

        self.setLayout(layout)

    def update_table(self, peak_info):

        centroids = peak_info['centroids'].round(3).astype(str)

        self.table.setRowCount(0)
        self.table.setRowCount(len(centroids))

        for row, centroid in enumerate(centroids):
            self.table.setItem(row, 0, QTableWidgetItem(centroid[0]))
            self.table.setItem(row, 1, QTableWidgetItem(centroid[1]))
            self.table.setItem(row, 2, QTableWidgetItem(centroid[2]))

    def get_cluster_parameters(self):

        params = [self.param_eps_line, self.param_min_line]
        valid_params = all([param.hasAcceptableInput() for param in params])

        if valid_params:

            return float(self.param_eps_line.text()), \
                   int(self.param_min_line.text())

    def add_peaks(self, peak_dict):

        self.plotter.clear()
        self.plotter.clear_actors()

        coordinates = np.array(peak_dict['coordinates'])
        clusters = np.array(peak_dict['clusters'])

        geoms, labels = [], []
        for uni in np.unique(clusters):
            coords = coordinates[clusters == uni]
            coords = np.row_stack([coords,-coords])
            points = pv.PolyData(coords)
            if uni >= 0:
                geoms.append(points)
                labels.append('C{}'.format(uni+1))
            else:
                self.plotter.add_mesh(points,
                                      color='k', 
                                      smooth_shading=True,
                                      point_size=5,
                                      render_points_as_spheres=True)

        multiblock = pv.MultiBlock(geoms)

        _, mapper = self.plotter.add_composite(multiblock,
                                               multi_colors=True,
                                               smooth_shading=True,
                                               point_size=10,
                                               render_points_as_spheres=True)

        colors = []
        for i in range(1,len(mapper.block_attr)):
            colors.append(mapper.block_attr[i].color)

        legend = [[label, color] for label, color in zip(labels, colors)]

        self.plotter.add_legend(legend,
                                loc='lower right',
                                bcolor='w',
                                face=None)

        self.plotter.enable_depth_peeling()

        self.change_proj()

    def change_proj(self):

        if self.proj_box.isChecked():
            self.plotter.enable_parallel_projection()
        else:
            self.plotter.disable_parallel_projection()

    def reset_view(self):

        self.plotter.reset_camera()
        self.plotter.view_isometric()

    def set_transform(self, T):

        if T is not None:

            b = pv._vtk.vtkMatrix4x4()
            for i in range(3):
                for j in range(3):
                    b.SetElement(i,j,T[i,j])

            actor = self.plotter.add_axes(xlabel='a*',
                                          ylabel='b*',
                                          zlabel='c*')
            actor.SetUserMatrix(b)

    def view_xy(self):

        self.plotter.view_xy()

    def view_yz(self):

        self.plotter.view_yz()

    def view_zx(self):

        self.plotter.view_zx()

    def view_yx(self):

        self.plotter.view_yx()

    def view_zy(self):

        self.plotter.view_zy()

    def view_xz(self):

        self.plotter.view_xz()

    def view_vector(self, vecs):

        if len(vecs) == 2:
            self.plotter.view_vector(*vecs)
        else:
            self.plotter.view_vector(vecs)

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