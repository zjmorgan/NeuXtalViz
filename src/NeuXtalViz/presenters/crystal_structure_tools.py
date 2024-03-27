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

            self.model.load_CIF(filename)

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

            self.view.draw_cell(self.model.get_unit_cell_transform())
            self.view.set_transform(self.model.get_transform())
            self.update_oriented_lattice()

            form, z = self.model.get_chemical_formula_z_parameter()
            self.view.set_formula_z(form, z)

            vol = self.model.get_unit_cell_volume()
            self.view.set_unit_cell_volume(vol)

    def update_atoms(self):

        params = self.view.get_lattice_constants()
        setting = self.view.get_setting()
        scatterers = self.view.get_scatterers()

        self.model.set_crystal_structure(params, setting, scatterers)

        atom_dict = self.model.generate_atom_positions()
        self.view.add_atoms(atom_dict)

        self.view.draw_cell(self.model.get_unit_cell_transform())
        self.view.set_transform(self.model.get_transform())

    def calculate_F2(self):

        d_min = self.view.get_minimum_d_spacing()

        params = self.view.get_lattice_constants()

        if params is not None:

            if d_min is None: 
                d_min = min(params[0:2])*0.2

            hkls, ds, F2s = self.model.generate_F2(d_min)

            self.view.set_factors(hkls, ds, F2s)

    def calculate_hkl(self):

        hkl = self.view.get_hkl()
        
        if hkl is not None:
            
            hkls, d, F2 = self.model.calculate_F2(*hkl)
            self.view.set_equivalents(hkls, d, F2)