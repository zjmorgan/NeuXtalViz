import sys
import numpy as np

from qtpy.QtWidgets import (QWidget,
                            QPushButton,
                            QLabel,
                            QHBoxLayout,
                            QGridLayout,
                            QSizePolicy)

from PyQt5.QtWidgets import QApplication, QMainWindow
from PyQt5.QtCore import Qt

from NeuXtalViz.config.atoms import indexing, groups

colors = {'Transition Metals': '#A1C9F4', # blue
          'Alkaline Earth Metals': '#FFB482', # orange
          'Nonmetals': '#8DE5A1', # green
          'Alkali Metals': '#FF9F9B', # red
          'Lanthanides': '#D0BBFF', # purple
          'Metalloids': '#DEBB9B', # brown
          'Actinides': '#FAB0E4', # pink
          'Other Metals': '#CFCFCF', # gray
          'Halogens': '#FFFEA3', # yellow
          'Noble Gases': '#B9F2F0'} # cyan

class PeriodicTable(QWidget):

    def __init__(self, parent=None):

        super().__init__(parent)

        layout = QHBoxLayout()

        table = self.__init_table()

        layout.addLayout(table)

        self.setLayout(layout)

    def __init_table(self):

        table = QGridLayout()

        for row in range(7):
            label = QLabel(str(row+1))
            table.addWidget(label, row+1, 0, Qt.AlignCenter)

        for col in range(18):
            label = QLabel(str(col+1))
            table.addWidget(label, 0, col+1, Qt.AlignCenter)

        self.atom_buttons = []

        for key in indexing.keys():
            row, col = indexing[key]
            button = QPushButton(key, self)
            button.setFixedSize(50, 50)
            group = groups.get(key)
            if group is not None:
                color = colors[group]
                button.setStyleSheet('background-color: {}'.format(color))
            self.atom_buttons.append(button)
            table.addWidget(button, row, col)

        return table

class MainWindow(QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.setWindowTitle('Periodic Table')

        widget = PeriodicTable()
        self.setCentralWidget(widget)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())