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
                            QTabWidget,
                            QFileDialog)

from qtpy.QtGui import QDoubleValidator, QIntValidator

import numpy as np
import pyvista as pv

from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure

from NeuXtalViz.views.base_view import NeuXtalVizWidget

class ModulationView(NeuXtalVizWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.tab_widget = QTabWidget(self)

        self.modulation_tab()

        self.layout().addWidget(self.tab_widget)

    def modulation_tab(self):

        mod_tab = QWidget()
        self.tab_widget.addTab(mod_tab, 'Modulation')        

        modulation_layout = QVBoxLayout()

        self.cluster_button = QPushButton('Cluster', self)
        self.load_UB_button = QPushButton('Load UB', self)
        self.load_peaks_button = QPushButton('Load Peaks', self)

        self.param_eps_line = QLineEdit('0.025')
        self.param_min_line = QLineEdit('15')

        notation = QDoubleValidator.StandardNotation

        validator = QDoubleValidator(0.0, 10, 5, notation=notation)

        self.param_eps_line.setValidator(validator)

        validator = QIntValidator(1, 1000)

        self.param_min_line.setValidator(validator)

        self.table = QTableWidget()

        self.table.setRowCount(0)
        self.table.setColumnCount(3)

        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.setHorizontalHeaderLabels(['h','k','l'])

        generate_layout = QHBoxLayout()

        generate_layout.addWidget(self.load_peaks_button)
        generate_layout.addWidget(self.load_UB_button)
        generate_layout.addWidget(self.cluster_button)

        cluster_layout = QFormLayout()

        cluster_layout.addRow('Maximum distance:', self.param_eps_line)
        cluster_layout.addRow('Minimum samples:', self.param_min_line)

        cluster_layout.addRow(self.table)

        plot_layout = QVBoxLayout()

        self.canvas = FigureCanvas(Figure(tight_layout=True))

        plot_layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        plot_layout.addWidget(self.canvas)

        fig = self.canvas.figure

        self.ax = fig.subplots(3, 1, sharex=True, sharey=True)

        for i in range(3):

            self.ax[i].set_xlim(-1,1)
            self.ax[i].set_ylim(1,100)
            self.ax[i].minorticks_on()
            self.ax[i].set_yscale('log')

        self.ax[0].set_xlabel('$[h00]$')
        self.ax[1].set_xlabel('$[0k0]$')
        self.ax[2].set_xlabel('$[00l]$')

        modulation_layout.addLayout(generate_layout)
        modulation_layout.addLayout(cluster_layout)
        modulation_layout.addLayout(plot_layout)

        mod_tab.setLayout(modulation_layout)

    def connect_cluster(self, cluster):

        self.cluster_button.clicked.connect(cluster)

    def connect_load_UB(self, load_UB):

        self.load_UB_button.clicked.connect(load_UB)

    def connect_load_peaks(self, load_peaks):

        self.load_peaks_button.clicked.connect(load_peaks)

    def load_UB_file_dialog(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        filename, _ = QFileDialog.getOpenFileName(self,
                                                  'Load UB file',
                                                  '',
                                                  'UB files (*.mat)',
                                                  options=options)

        return filename

    def load_peaks_file_dialog(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        filter_string = 'Peaks files (*.integrate *.peaks *.nxs)'

        filename, _ = QFileDialog.getOpenFileName(self,
                                                  'Load peaks file',
                                                  '',
                                                  filter_string,
                                                  options=options)

        return filename

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

        self.plotter.clear_actors()

        for i in range(3):
            self.ax[i].clear()

        bins = np.linspace(-1.025, 1.025, 42)

        coordinates = np.array(peak_dict['coordinates'])
        clusters = np.array(peak_dict['clusters'])

        vectors = peak_dict['translation']
        T = peak_dict['transform']
        T_inv = peak_dict['inverse']

        translations = np.array(np.meshgrid([-1, 0, 1],
                                            [-1, 0, 1],
                                            [-1, 0, 1])).T.reshape(-1,3)

        offsets = np.dot(translations, vectors)

        geoms, labels = [], []
        for uni in np.unique(clusters):
            coords = coordinates[clusters == uni]
            coords = np.row_stack([coords, -coords])
            coords = (coords[:,np.newaxis,:]+offsets).reshape(-1,3)
            delta = (T_inv @ coords.T).T
            mask = (np.abs(delta) < 1).all(axis=1)
            coords = coords[mask]
            delta = delta[mask]
            points = pv.PolyData(coords)
            if uni >= 0:
                geoms.append(points)
                color = 'C{}'.format(uni+1)
                labels.append(color)
                h, _ = np.histogram(delta[:,0], bins=bins)
                k, _ = np.histogram(delta[:,1], bins=bins)
                l, _ = np.histogram(delta[:,2], bins=bins)
                self.ax[0].stairs(h, bins, color=color)
                self.ax[1].stairs(k, bins, color=color)
                self.ax[2].stairs(l, bins, color=color)
            else:
                self.plotter.add_mesh(points,
                                      color='k', 
                                      smooth_shading=True,
                                      point_size=5,
                                      render_points_as_spheres=True)

        for i in range(3):
            self.ax[i].minorticks_on()
            self.ax[i].set_yscale('log')

        self.ax[0].set_xlabel('$[h00]$')
        self.ax[1].set_xlabel('$[0k0]$')
        self.ax[2].set_xlabel('$[00l]$')

        self.canvas.draw_idle()
        self.canvas.flush_events()

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

        A = np.eye(4)
        A[:3,:3] = T

        mesh = pv.Box(bounds=(-1,1,-1,1,-1,1), level=0, quads=True)
        mesh.transform(A, inplace=True)

        self.plotter.add_mesh(mesh,
                              color='k',
                              style='wireframe',
                              render_lines_as_tubes=True)

        self.plotter.add_legend(legend,
                                loc='lower right',
                                bcolor='w',
                                face=None)

        self.plotter.enable_depth_peeling()

        self.reset_view()