from NeuXtalViz.presenters.base_presenter import NeuXtalVizPresenter

class Experiment(NeuXtalVizPresenter):

    def __init__(self, view, model):

        super(Experiment, self).__init__(view, model)

        self.view.connect_load_UB(self.load_UB)
        self.view.connect_switch_instrument(self.switch_instrument)
        self.view.connect_update_goniometer(self.update_goniometer)
        self.view.connect_switch_crystal(self.switch_crystal)
        self.view.connect_switch_point_group(self.switch_group)
        self.view.connect_wavelength(self.update_wavelength)
        self.view.connect_optimize(self.optimize_coverage)

        self.switch_instrument()
        self.switch_crystal()

    def load_UB(self):

        filename = self.view.load_UB_file_dialog()

        if filename:

            self.model.load_UB(filename)

            self.update_oriented_lattice()

            self.view.set_transform(self.model.get_transform())

    def switch_instrument(self):

        instrument = self.view.get_instrument()

        wavelength = self.model.get_wavelength(instrument)
        motors = self.model.get_motors(instrument)
        modes = self.model.get_modes(instrument)
        goniometers = self.model.get_goniometers(instrument, modes[0])

        self.view.set_modes(modes)
        self.view.set_wavelength(wavelength)
        self.view.update_tables(goniometers, motors)

    def switch_crystal(self):

        cs = self.view.get_crystal_system()

        point_groups = self.model.get_crystal_system_point_groups(cs)

        self.view.set_point_groups(point_groups)

        self.switch_group()

    def switch_group(self):

        pg = self.view.get_point_group()

        centerings = self.model.get_point_group_centering(pg)

        self.view.set_lattice_centerings(centerings)

    def update_goniometer(self):

        instrument = self.view.get_instrument()
        mode = self.view.get_mode()        

        goniometers = self.model.get_goniometers(instrument, mode)
        motors = self.model.get_motors(instrument)

        self.view.update_tables(goniometers, motors)

    def update_wavelength(self):

        wl_min, wl_max = self.view.get_wavelength()
        self.view.update_wavelength(wl_min)

    def optimize_coverage(self):

        instrument = self.view.get_instrument()
        wavelength = self.view.get_wavelength()
        logs = self.view.get_motors()
        instrument_name = self.model.get_instrument_name(instrument)

        self.model.generate_instrument_coverage(instrument_name,
                                                logs,
                                                wavelength)
    
        laue = self.view.get_laue_symmetry()

        self.model.apply_symmetry(laue)

        self.view.set_transform(self.model.get_transform())

        coverage_dict = self.model.get_coverage_info()
        self.view.add_coverage(coverage_dict)