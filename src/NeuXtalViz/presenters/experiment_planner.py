from NeuXtalViz.presenters.base_presenter import NeuXtalVizPresenter

class Experiment(NeuXtalVizPresenter):

    def __init__(self, view, model):

        super(Experiment, self).__init__(view, model)

        self.view.connect_load_UB(self.load_UB)
        self.view.connect_switch_instrument(self.switch_instrument)
        self.view.connect_wavelength(self.update_wavelength)
        self.view.connect_optimize(self.optimize_coverage)

        self.switch_instrument()

    def load_UB(self):

        filename = self.view.load_UB_file_dialog()

        if filename:

            self.model.load_UB(filename)

            self.update_oriented_lattice()

            self.view.set_transform(self.model.get_transform())

    def switch_instrument(self):

        instrument = self.view.get_instrument()

        wavelength = self.model.get_wavelength(instrument)
        goniometers = self.model.get_goniometers(instrument)
        motors = self.model.get_motors(instrument)

        self.view.set_wavelength(wavelength)
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

