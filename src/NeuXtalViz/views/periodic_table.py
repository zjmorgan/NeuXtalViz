import sys
import numpy as np

from qtpy.QtWidgets import (QWidget,
                            QPushButton,
                            QButtonGroup,
                            QLabel,
                            QComboBox,
                            QHBoxLayout,
                            QVBoxLayout,
                            QGridLayout)

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

class PeriodicTableView(QWidget):

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

        self.atom_buttons = QButtonGroup()
        self.atom_buttons.setExclusive(True)

        for key in indexing.keys():
            row, col = indexing[key]
            button = QPushButton(key, self)
            button.setFixedSize(50, 50)
            group = groups.get(key)
            if group is not None:
                color = colors[group]
                button.setStyleSheet('background-color: {}'.format(color))
            self.atom_buttons.addButton(button)
            table.addWidget(button, row, col)

        return table

    def get_atom_view(self):

        return AtomView()

    def connect_atoms(self, atom_info):

        self.atom_info = atom_info

        self.atom_buttons.buttonClicked.connect(self.show_atom_dialog)

    def show_atom_dialog(self, button):

        return self.atom_info(button.text())

class AtomView(QWidget):

    def __init__(self):

        super().__init__()

        card = QGridLayout()

        self.z_label = QLabel('1')
        self.symbol_label = QLabel('H')
        self.name_label = QLabel('Hydrogen')

        self.isotope_combo = QComboBox(self)

        card.addWidget(self.z_label, 0, 0, Qt.AlignCenter)
        card.addWidget(self.isotope_combo, 0, 1, 1, 2)
        card.addWidget(self.symbol_label, 1, 1, 1, 2, Qt.AlignCenter)
        card.addWidget(self.name_label, 2, 1, 1, 2, Qt.AlignCenter)

        self.setLayout(card)