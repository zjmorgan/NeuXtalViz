from NeuXtalViz.presenters.periodic_table import PeriodicTable
from NeuXtalViz.presenters.base_presenter import NeuXtalVizPresenter

class CrystalStructure(NeuXtalVizPresenter):

    def __init__(self, view, model):

        super(CrystalStructure, self).__init__(view, model)

        self.view.connect_group_generator(self.generate_groups)
        self.view.connect_setting_generator(self.generate_settings)
        self.view.connect_F2_calculator(self.calculate_F2)
        self.view.connect_hkl_calculator(self.calculate_hkl)
        self.view.connect_row_highligter(self.highlight_row)
        self.view.connect_lattice_parameters(self.update_parameters)
        self.view.connect_atom_table(self.set_atom_table)
        self.view.connect_load_CIF(self.load_CIF)
        self.view.connect_select_isotope(self.select_isotope)

        self.generate_groups()
        self.generate_settings()

    def highlight_row(self):

        scatterer = self.view.get_scatterer()
        self.view.set_atom(scatterer)

    def set_atom_table(self):

        self.view.set_atom_table()
        self.update_atoms()

    def update_parameters(self):

        params = self.view.get_lattice_constants()
        params = self.model.update_parameters(params)
        self.model.update_lattice_parameters(*params)
        self.view.set_lattice_constants(params)
        vol = self.model.get_unit_cell_volume()
        self.view.set_unit_cell_volume(vol)

        atom_dict = self.model.generate_atom_positions()
        self.view.add_atoms(atom_dict)

        self.view.draw_cell(self.model.get_unit_cell_transform())
        self.view.set_transform(self.model.get_transform())

    def generate_groups(self):

        system = self.view.get_crystal_system()
        nos = self.model.generate_space_groups_from_crystal_system(system)
        self.view.update_space_groups(nos)

        self.generate_settings()

    def generate_settings(self):

        no = self.view.get_space_group()
        settings = self.model.generate_settings_from_space_group(no)
        self.view.update_settings(settings)

    def load_CIF(self):

        filename = self.view.load_CIF_file_dialog()

        if filename:

            self.update_processing()

            self.update_processing('Loading CIF...', 10) 

            self.model.load_CIF(filename)

            self.update_processing('Loading CIF...', 50)

            crystal_system = self.model.get_crystal_system()
            space_group = self.model.get_space_group()
            setting = self.model.get_setting()
            params = self.model.get_lattice_constants()
            scatterers = self.model.get_scatterers()

            self.view.set_crystal_system(crystal_system)
            self.generate_groups()
            self.view.set_space_group(space_group)
            self.generate_settings()
            self.view.set_setting(setting)
            self.view.set_lattice_constants(params)
            self.view.set_scatterers(scatterers)

            params = self.model.constrain_parameters()
            self.view.constrain_parameters(params)

            atom_dict = self.model.generate_atom_positions()
            self.view.add_atoms(atom_dict)

            self.update_processing('Loading CIF...', 80)

            self.view.draw_cell(self.model.get_unit_cell_transform())
            self.view.set_transform(self.model.get_transform())
            self.update_oriented_lattice()

            form, z = self.model.get_chemical_formula_z_parameter()
            self.view.set_formula_z(form, z)

            self.update_processing('Loading CIF...', 99)

            vol = self.model.get_unit_cell_volume()
            self.view.set_unit_cell_volume(vol)

            self.update_complete('CIF loaded!')

        else:

            self.update_invalid()

    def update_atoms(self):

        params = self.view.get_lattice_constants()
        setting = self.view.get_setting()
        scatterers = self.view.get_scatterers()

        self.model.set_crystal_structure(params, setting, scatterers)

        atom_dict = self.model.generate_atom_positions()
        self.view.add_atoms(atom_dict)

        form, z = self.model.get_chemical_formula_z_parameter()
        self.view.set_formula_z(form, z)

        self.view.draw_cell(self.model.get_unit_cell_transform())
        self.view.set_transform(self.model.get_transform())

    def calculate_F2(self):

        worker = self.view.worker(self.calculate_F2_process)
        worker.connect_result(self.calculate_F2_complete)
        worker.connect_finished(self.update_complete)
        worker.connect_progress(self.update_processing)

        self.view.start_worker_pool(worker)

    def calculate_F2_complete(self, result):

        if result is not None:

            self.view.set_factors(*result)

    def calculate_F2_process(self, progress):

        d_min = self.view.get_minimum_d_spacing()

        params = self.view.get_lattice_constants()

        if params is not None:

            progress('Processing...', 1)

            progress('Calculating factors...', 10)

            if d_min is None:
                d_min = min(params[0:2])*0.2

            hkls, ds, F2s = self.model.generate_F2(d_min)

            progress('Factors calculated...', 99)

            progress('Factors calculated!', 100)

            return hkls, ds, F2s

        else:

            progress('Invalid parameters.', 0)

    def calculate_hkl(self):

        worker = self.view.worker(self.calculate_hkl_process)
        worker.connect_result(self.calculate_hkl_complete)
        worker.connect_finished(self.update_complete)
        worker.connect_progress(self.update_processing)

        self.view.start_worker_pool(worker)

    def calculate_hkl_complete(self, result):

        if result is not None:

            self.view.set_equivalents(*result)

    def calculate_hkl_process(self, progress):

        hkl = self.view.get_hkl()

        if hkl is not None:

            progress('Processing...', 1)

            progress('Calculating equivalents...', 10)

            hkls, d, F2 = self.model.calculate_F2(*hkl)

            progress('Equivalents calculated...', 99)

            progress('Equivalents calculated!', 100)

            return hkls, d, F2

        else:

            progress('Invalid parameters.', 0)

    def select_isotope(self):

        atom = self.view.get_isotope()

        if atom != '':

            view = self.view.get_periodic_table()
            model = self.model.get_periodic_table(atom)

            self.periodic_table = PeriodicTable(view, model)
            self.periodic_table.view.connect_selected(self.update_selection)
            self.periodic_table.view.show()

    def update_selection(self, data):

        self.view.set_isotope(data)
        self.view.set_atom_table()
        self.update_atoms()