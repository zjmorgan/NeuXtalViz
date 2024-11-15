from qtpy.QtWidgets import (QWidget,
                            QTableWidget,
                            QTableWidgetItem,
                            QHeaderView,
                            QHBoxLayout,
                            QVBoxLayout,
                            QPushButton,
                            QCheckBox,
                            QComboBox,
                            QLineEdit,
                            QLabel,
                            QTabWidget,
                            QFileDialog)

from qtpy.QtGui import QDoubleValidator, QIntValidator

import numpy as np
import matplotlib.pyplot as plt
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

        self.layout().addWidget(self.tab_widget, stretch=1)

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

        validator = QDoubleValidator(0.0001, 10, 5, notation=notation)

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

        cluster_layout = QVBoxLayout()
        params_layout = QHBoxLayout()

        dist_label = QLabel('Maximum distance:', self)
        samp_label = QLabel('Minimum samples:', self)

        params_layout.addWidget(dist_label)
        params_layout.addWidget(self.param_eps_line)
        params_layout.addWidget(samp_label)
        params_layout.addWidget(self.param_min_line)

        cluster_layout.addLayout(params_layout)
        cluster_layout.addWidget(self.table)

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

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filename, _ = file_dialog.getOpenFileName(self,
                                                  'Load UB file',
                                                  '',
                                                  'UB files (*.mat)',
                                                  options=options)

        return filename

    def load_peaks_file_dialog(self):

        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog

        file_dialog = QFileDialog()
        file_dialog.setFileMode(QFileDialog.AnyFile)

        filter_string = 'Peaks files (*.integrate *.peaks *.nxs)'

        filename, _ = file_dialog.getOpenFileName(self,
                                                  'Load peaks file',
                                                  '',
                                                  filter_string,
                                                  options=options)

        return filename

    def update_table(self, peak_info):

        centroids = peak_info['satellites'].round(3).astype(str)

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

        multiblock = pv.MultiBlock()

        for uni in np.unique(clusters):
            coords = coordinates[clusters == uni]
            coords = (coords[:,np.newaxis,:]+offsets).reshape(-1,3)
            delta = (T_inv @ coords.T).T
            mask = (np.abs(delta) < 1).all(axis=1)
            coords = coords[mask]
            delta = delta[mask]
            points = pv.PolyData(coords)
            if uni >= 0:
                color = 'C{}'.format(uni)
                multiblock[color] = points
                if uni > 0:
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

        _, mapper = self.plotter.add_composite(multiblock,
                                               multi_colors=True,
                                               smooth_shading=True,
                                               point_size=10,
                                               render_points_as_spheres=True)

        prop_cycle = plt.rcParams['axes.prop_cycle']

        cmap = prop_cycle.by_key()['color']

        colors = []
        for i in range(1,len(mapper.block_attr)):
            colors.append(cmap[i-1])
            mapper.block_attr[i].color = cmap[i-1]

        legend = [['C{}'.format(i), color] for i, color in enumerate(colors)]

        A = np.eye(4)
        A[:3,:3] = T

        mesh = pv.Box(bounds=(-1,1,-1,1,-1,1), level=0, quads=True)
        mesh.transform(A, inplace=True)

        self.plotter.add_mesh(mesh,
                              color='k',
                              style='wireframe',
                              render_lines_as_tubes=True)

        for point in [(1,0,0), (0,1,0), (0,0,1)]:

            mesh = pv.Line(pointa=-np.array(point), pointb=point, resolution=1)
            mesh.transform(A, inplace=True)

            self.plotter.add_mesh(mesh,
                                  color='k',
                                  style='wireframe',
                                  render_lines_as_tubes=True)

        pointsa = [(-1,-1), (-1,1), (1,1), (1,-1)]
        pointsb = [(-1,1), (1,1), (1,-1), (-1,-1)]

        for i in range(4):
            
            a, b = pointsa[i], pointsb[i]

            mesh = pv.Line(pointa=(a[0],a[1],0),
                           pointb=(b[0],b[1],0),
                           resolution=1)

            mesh.transform(A, inplace=True)

            self.plotter.add_mesh(mesh,
                                  color='k',
                                  style='wireframe',
                                  render_lines_as_tubes=True)

            mesh = pv.Line(pointa=(a[0],0,a[1]),
                           pointb=(b[0],0,b[1]),
                           resolution=1)

            mesh.transform(A, inplace=True)

            self.plotter.add_mesh(mesh,
                                  color='k',
                                  style='wireframe',
                                  render_lines_as_tubes=True)

            mesh = pv.Line(pointa=(0,a[0],a[1]),
                           pointb=(0,b[0],b[1]),
                           resolution=1)

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