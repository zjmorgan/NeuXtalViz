from NeuXtalViz.presenters.base_presenter import NeuXtalVizPresenter


class Experiment(NeuXtalVizPresenter):
    def __init__(self, view, model):
        super(Experiment, self).__init__(view, model)

        self.view.connect_load_UB(self.load_UB)
        self.view.connect_switch_instrument(self.switch_instrument)
        self.view.connect_update_goniometer(self.update_goniometer)
        self.view.connect_switch_crystal(self.switch_crystal)
        self.view.connect_switch_point_group(self.switch_group)
        self.view.connect_switch_lattice_centering(self.switch_centering)
        self.view.connect_wavelength(self.update_wavelength)
        self.view.connect_optimize(self.optimize_coverage)
        self.view.connect_calculate_single(self.calculate_single)
        self.view.connect_calculate_double(self.calculate_double)
        self.view.connect_add_orientation(self.add_orientation)

        self.view.connect_roi_ready(self.lookup_angle)
        self.view.connect_viz_ready(self.visualize)

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

        self.model.remove_instrument()

    def switch_crystal(self):
        cs = self.view.get_crystal_system()

        point_groups = self.model.get_crystal_system_point_groups(cs)

        self.view.set_point_groups(point_groups)

        self.switch_group()

    def switch_group(self):
        pg = self.view.get_point_group()

        centerings = self.model.get_point_group_centering(pg)

        self.view.set_lattice_centerings(centerings)

        self.visualize()

    def switch_centering(self):
        self.visualize()

    def update_goniometer(self):
        instrument = self.view.get_instrument()
        mode = self.view.get_mode()

        goniometers = self.model.get_goniometers(instrument, mode)
        motors = self.model.get_motors(instrument)

        self.view.update_tables(goniometers, motors)

    def update_wavelength(self):
        wl_min, wl_max = self.view.get_wavelength()
        self.view.update_wavelength(wl_min)

    def create_instrument(self):
        instrument = self.view.get_instrument()
        motors = self.view.get_motors()

        self.model.initialize_instrument(instrument, motors)

    def calculate_single(self):
        worker = self.view.worker(self.calculate_single_process)
        worker.connect_result(self.calculate_single_complete)
        worker.connect_finished(self.visualize)
        worker.connect_progress(self.update_processing)

        self.view.start_worker_pool(worker)

    def calculate_single_complete(self, result):
        if result is not None:
            self.view.plot_instrument(self.model.gamma, self.model.nu, *result)

    def calculate_single_process(self, progress):
        hkl_1, hkl_2 = self.view.get_input_hkls()
        wavelength = self.view.get_wavelength()

        instrument = self.view.get_instrument()
        mode = self.view.get_mode()
        axes, polarities = self.model.get_axes_polarities(instrument, mode)

        limits = self.view.get_goniometer_limits()

        if hkl_1 is not None and self.model.has_UB():
            progress("Initializing instrument", 5)

            self.create_instrument()

            progress("Instrument initialized! ", 10)

            progress("Calculating peak coverage", 15)

            gamma, nu, lamda = self.model.individual_peak(
                hkl_1, wavelength, axes, polarities, limits
            )

            progress("Peak calculated!", 0)

            return gamma, nu, lamda

        else:
            progress("Invalid parameters.", 0)

    def calculate_double(self):
        worker = self.view.worker(self.calculate_double_process)
        worker.connect_result(self.calculate_double_complete)
        worker.connect_finished(self.visualize)
        worker.connect_progress(self.update_processing)

        self.view.start_worker_pool(worker)

    def calculate_double_complete(self, result):
        if result is not None:
            self.view.plot_instrument_alternate(
                self.model.gamma, self.model.nu, *result
            )

    def calculate_double_process(self, progress):
        hkl_1, hkl_2 = self.view.get_input_hkls()
        wavelength = self.view.get_wavelength()

        instrument = self.view.get_instrument()
        mode = self.view.get_mode()
        axes, polarities = self.model.get_axes_polarities(instrument, mode)

        limits = self.view.get_goniometer_limits()

        if hkl_1 is not None and hkl_2 is not None and self.model.has_UB():
            progress("Initializing instrument", 5)

            self.create_instrument()

            progress("Instrument initialized! ", 10)

            progress("Calculating peaks coverage", 15)

            peak_1, peak_2 = self.model.simultaneous_peaks(
                hkl_1, hkl_2, wavelength, axes, polarities, limits
            )

            gamma_1, nu_1, lamda_1 = peak_1
            gamma_2, nu_2, lamda_2 = peak_2

            progress("Peaks calculated!", 0)

            return gamma_1, nu_1, lamda_1, gamma_2, nu_2, lamda_2

        else:
            progress("Invalid parameters.", 0)

    def lookup_angle(self):
        gamma = self.view.get_horizontal()
        nu = self.view.get_vertical()

        angles = self.model.get_angles(gamma, nu)
        self.view.set_angles(angles)

    def add_orientation(self):
        worker = self.view.worker(self.add_orientation_process)
        worker.connect_result(self.add_orientation_complete)
        worker.connect_finished(self.visualize)
        worker.connect_progress(self.update_processing)

        self.view.start_worker_pool(worker)

    def add_orientation_complete(self, result):
        angles, all_angles, free_angles = result

        comment = self.model.comment
        update_angles = []
        for angle, angle_name in zip(angles, all_angles):
            if angle_name in free_angles:
                update_angles.append(angle)

        self.view.add_orientation(comment, update_angles)

    def add_orientation_process(self, progress):
        angles = self.view.get_angles()
        free_angles = self.view.get_free_angles()
        all_angles = self.view.get_all_angles()

        wavelength = self.view.get_wavelength()
        d_min = self.view.get_d_min()
        # centering = self.view.get_lattice_centering()
        rows = self.view.get_number_of_orientations()

        if len(angles) > 0:
            progress("Calculation reflections", 5)

            self.model.add_orientation(angles, wavelength, d_min, rows)

            progress("Reflections calculated!", 0)

            return angles, all_angles, free_angles

        else:
            progress("Invalid parameters.", 0)

    def visualize(self):
        point_group = self.view.get_point_group()
        lattice_centering = self.view.get_lattice_centering()
        use = self.view.get_orientations_to_use()
        d_min = self.view.get_d_min()

        stats = self.model.calculate_statistics(
            point_group, lattice_centering, use, d_min
        )

        if stats is not None:
            self.view.plot_statistics(*stats)

            peak_dict = self.model.get_coverage_info(point_group)

            self.view.add_peaks(peak_dict)

    def optimize_coverage(self):
        worker = self.view.worker(self.optimize_coverage_process)
        worker.connect_result(self.optimize_coverage_complete)
        worker.connect_finished(self.visualize)
        worker.connect_progress(self.update_processing)

        self.view.start_worker_pool(worker)

    def optimize_coverage_complete(self, result):
        if result is not None:
            for angles in result:
                self.view.add_orientation("CrystalPlan", angles)

    def optimize_coverage_process(self, progress):
        point_group = self.view.get_point_group()
        lattice_centering = self.view.get_lattice_centering()
        use = self.view.get_orientations_to_use()
        opt = self.view.get_optimized_settings()
        d_min = self.view.get_d_min()
        wavelength = self.view.get_wavelength()
        n_orient = self.view.get_settings()

        n_elite = 2
        n_gener = 10
        n_indiv = 10
        mutation_rate = 0.15

        instrument = self.view.get_instrument()
        mode = self.view.get_mode()
        axes = self.model.get_goniometer_axes(instrument, mode)
        limits = self.view.get_goniometer_limits()

        if self.model.has_UB():
            progress("Initializing instrument", 5)

            self.create_instrument()

            progress("Instrument initialized! ", 10)

            cp = self.model.crystal_plan(
                use,
                opt,
                axes,
                limits,
                wavelength,
                d_min,
                point_group,
                lattice_centering,
            )

            progress("Optimizing peaks coverage", 15)

            values = cp.optimize(
                n_orient, n_indiv, n_gener, n_elite, mutation_rate
            )

            progress("Peaks coverage optimized!", 0)

            return values

        else:
            progress("Invalid parameters.", 0)
