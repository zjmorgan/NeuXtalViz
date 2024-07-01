import sys

from qtpy.QtWidgets import (QWidget,
                            QLineEdit,
                            QLabel,
                            QPushButton,
                            QComboBox,
                            QTableWidget,
                            QTableWidgetItem,
                            QHeaderView,
                            QFrame,
                            QHBoxLayout,
                            QVBoxLayout,
                            QGridLayout)

from matplotlib.backends.backend_qtagg import FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure

from qtpy.QtGui import QDoubleValidator
from PyQt5.QtCore import Qt

from PyQt5.QtWidgets import QApplication, QMainWindow

class ReductionPlanView(QWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        experiment_layout = QHBoxLayout()
        dataset_layout = QGridLayout()

        self.instrument_combo = QComboBox(self)
        self.instrument_combo.addItem('TOPAZ')
        self.instrument_combo.addItem('MANDI')
        self.instrument_combo.addItem('CORELLI')
        self.instrument_combo.addItem('SNAP')
        self.instrument_combo.addItem('DEMAND')
        self.instrument_combo.addItem('WAND²')

        self.ipts_combo = QComboBox(self)
        self.experiment_combo = QComboBox(self)

        experiment_layout.addWidget(self.instrument_combo)
        experiment_layout.addWidget(self.ipts_combo)
        experiment_layout.addWidget(self.experiment_combo)

        self.runs_label = QLabel('Run Numbers', self)
        self.runs_line = QLineEdit(self)

        dataset_layout.addWidget(self.runs_label, 0, 0)
        dataset_layout.addWidget(self.runs_line, 0, 1)

        self.table = QTableWidget()

        header = ['Axis','Direction','Sense','Min','Max']

        self.table.setRowCount(5)
        self.table.setColumnCount(len(header))
        self.table.setHorizontalHeaderLabels(header)

        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)

        instrument_layout = QVBoxLayout()

        self.canvas = FigureCanvas(Figure())

        instrument_layout.addLayout(experiment_layout)
        instrument_layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        instrument_layout.addWidget(self.canvas)

        calibration_layout = QGridLayout()
        vanadium_layout = QGridLayout()
        grouping_layout = QGridLayout()

        self.calibration_label = QLabel('Calibration', self)
        self.detector_label = QLabel('Detector', self)
        self.tube_label = QLabel('Tube', self)
        self.vanadium_label = QLabel('Vanadium', self)
        self.counts_label = QLabel('Counts', self)
        self.spectrum_label = QLabel('Spectrum', self)
        self.grouping_label = QLabel('Grouping', self)
        self.mask_label = QLabel('Mask', self)
        self.background_label = QLabel('Background', self)

        self.calibration_combo = QComboBox(self)
        self.detector_line = QLineEdit(self)
        self.tube_line = QLineEdit(self)
        self.vanadium_combo = QComboBox(self)
        self.counts_line = QLineEdit(self)
        self.spectrum_line = QLineEdit(self)
        self.grouping_combo = QComboBox(self)
        self.mask_line = QLineEdit(self)
        self.background_line = QLineEdit(self)

        self.detector_button = QPushButton('Browse', self)
        self.tube_button = QPushButton('Browse', self)
        self.counts_button = QPushButton('Browse', self)
        self.spectrum_button = QPushButton('Browse', self)
        self.mask_button = QPushButton('Browse', self)
        self.background_button = QPushButton('Browse', self)

        calibration_layout.addWidget(self.calibration_label, 0, 0)
        calibration_layout.addWidget(self.calibration_combo, 0, 1)
        calibration_layout.addWidget(self.detector_label, 1, 0)
        calibration_layout.addWidget(self.detector_line, 1, 1)
        calibration_layout.addWidget(self.detector_button, 1, 2)
        calibration_layout.addWidget(self.tube_label, 2, 0)
        calibration_layout.addWidget(self.tube_line, 2, 1)
        calibration_layout.addWidget(self.tube_button, 2, 2)

        vanadium_layout.addWidget(self.vanadium_label, 0, 0)
        vanadium_layout.addWidget(self.vanadium_combo, 0, 1)
        vanadium_layout.addWidget(self.counts_label, 1, 0)
        vanadium_layout.addWidget(self.counts_line, 1, 1)
        vanadium_layout.addWidget(self.counts_button, 1, 2)
        vanadium_layout.addWidget(self.spectrum_label, 2, 0)
        vanadium_layout.addWidget(self.spectrum_line, 2, 1)
        vanadium_layout.addWidget(self.spectrum_button, 2, 2)

        grouping_layout.addWidget(self.grouping_label, 0, 0)
        grouping_layout.addWidget(self.grouping_combo, 0, 1)
        grouping_layout.addWidget(self.mask_label, 1, 0)
        grouping_layout.addWidget(self.mask_line, 1, 1)
        grouping_layout.addWidget(self.mask_button, 1, 2)
        grouping_layout.addWidget(self.background_label, 2, 0)
        grouping_layout.addWidget(self.background_line, 2, 1)
        grouping_layout.addWidget(self.background_button, 2, 2)

        parameters_layout = QVBoxLayout()
        parameters_layout.addLayout(dataset_layout)
        parameters_layout.addWidget(self.table)
        parameters_layout.addLayout(calibration_layout)
        parameters_layout.addLayout(vanadium_layout)
        parameters_layout.addLayout(grouping_layout)
        parameters_layout.addStretch()

        layout = QHBoxLayout()

        layout.addLayout(instrument_layout)
        layout.addLayout(parameters_layout)

        self.setLayout(layout)

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.setWindowTitle('Reduction Plan')

        widget = ReductionPlanView()
        self.setCentralWidget(widget)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())