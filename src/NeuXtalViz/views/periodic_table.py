import re
import numpy as np

from qtpy.QtWidgets import (QWidget,
                            QPushButton,
                            QButtonGroup,
                            QLabel,
                            QComboBox,
                            QHBoxLayout,
                            QVBoxLayout,
                            QGridLayout)

from PyQt5.QtCore import Qt, pyqtSignal

from NeuXtalViz.config.atoms import indexing, groups, isotopes

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

    selection = pyqtSignal(str)

    def __init__(self, parent=None):

        super().__init__(parent)

        layout = QHBoxLayout()

        table = self.__init_table()

        layout.addLayout(table)

        self.setLayout(layout)

        self.value = 'H'

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
            if isotopes.get(key) is None:
                button.setDisabled(True)
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

    def closeEvent(self, event):

        self.selection.emit(self.value)
        event.accept()

    def connect_selected(self, value):

        self.selection.connect(value)

class AtomView(QWidget):
    
    selection = pyqtSignal(str)

    def __init__(self):

        super().__init__()

        card = QGridLayout()

        self.z_label = QLabel('1')
        self.symbol_label = QLabel('H')
        self.name_label = QLabel('Hydrogen')
        self.mass_label = QLabel('1.007825')
        self.abundance_label = QLabel('99.9885')

        self.sigma_coh_label = QLabel('')
        self.sigma_inc_label = QLabel('')
        self.sigma_tot_label = QLabel('')
        self.b_coh_label = QLabel('')
        self.b_inc_label = QLabel('')

        self.isotope_combo = QComboBox(self)

        self.select_button = QPushButton('Use Isotope', self)

        card.addWidget(self.z_label, 0, 0, Qt.AlignCenter)
        card.addWidget(self.isotope_combo, 0, 1, 1, 2)
        card.addWidget(self.symbol_label, 1, 1, 1, 2, Qt.AlignCenter)
        card.addWidget(self.name_label, 2, 1, 1, 2, Qt.AlignCenter)
        card.addWidget(self.mass_label, 3, 1, 1, 2, Qt.AlignCenter)

        card.addWidget(self.abundance_label, 0, 3)
        card.addWidget(self.sigma_tot_label, 1, 3)
        card.addWidget(self.sigma_coh_label, 2, 3)
        card.addWidget(self.sigma_inc_label, 3, 3)
        card.addWidget(self.b_coh_label, 2, 4)
        card.addWidget(self.b_inc_label, 3, 4)
        card.addWidget(self.select_button, 0, 4)

        self.setLayout(card)

    def closeEvent(self, event):

        self.selection.emit(self.get_selection())
        event.accept()

    def connect_selected(self, value):

        self.selection.connect(value)

    def connect_isotopes(self, update_info):

        self.isotope_combo.currentIndexChanged.connect(update_info)

    def connect_selection(self, use_isotope):

        self.select_button.clicked.connect(use_isotope)

    def set_symbol_name(self, symbol, name):

        self.symbol_label.setText(symbol)
        self.name_label.setText(name)

    def get_selection(self):

        isotope = self.symbol_label.text()#+self.isotope_combo.currentText()
        
        return isotope

    def set_isotope_numbers(self, numbers):

        self.isotope_combo.clear()     
        if numbers is not None:
            self.isotope_combo.addItems(np.array(numbers).astype(str).tolist())

    def get_isotope(self):

        iso = self.isotope_combo.currentText()
        if iso is not None:
            return 0 if iso == '' else int(iso)

    def set_atom_parameters(self, atom, scatt):

        self.z_label.setText(str(atom['z']))
        self.mass_label.setText(str(atom['mass']))   
        self.abundance_label.setText(str(atom['abundance']))   

        self.sigma_coh_label.setText('σ(coh) = {}'.format(scatt['sigma_coh']))
        self.sigma_inc_label.setText('σ(inc) = {}'.format(scatt['sigma_inc']))
        self.sigma_tot_label.setText('σ(tot) = {}'.format(scatt['sigma_tot']))

        self.b_coh_label.setText('b(coh) = {}+{}i'.format(scatt['b_coh_re'],
                                                          scatt['b_coh_im']))
        self.b_inc_label.setText('b(inc) = {}+{}i'.format(scatt['b_inc_re'],
                                                          scatt['b_inc_im']))