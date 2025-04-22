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
        self.view.connect_delete_angles(self.delete_angles)
        self.view.connect_save_CSV(self.save_CSV)
        self.view.connect_save_experiment(self.save_experiment)
        self.view.connect_load_experiment(self.load_experiment)
        self.view.connect_peak_table(self.update_peaks)

        self.view.connect_roi_ready(self.lookup_angle)
        self.view.connect_viz_ready(self.visualize)

        self.view.connect_update(self.view.update_counting)
        self.view.connect_highlight_angles(self.view.highlight_angles)

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
        options = self.model.get_counting_options(instrument)

        self.view.set_modes(modes)
        self.view.set_wavelength(wavelength)
        self.view.update_tables(goniometers, motors)
        self.view.set_counting_options(options)

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

    def update_peaks(self):
        row = self.view.get_peak_list()
        if row is not None:
            peak_list = self.model.generate_table(row)
            self.view.update_peaks_table(peak_list)

    def lookup_angle(self):
        gamma = self.view.get_horizontal()
        nu = self.view.get_vertical()

        vals = self.model.get_angles(gamma, nu)
        if vals is not None:
            angles, gamma, nu = vals
            self.view.set_angles(angles)
            self.view.set_horizontal(gamma)
            self.view.set_vertical(nu)
            self.view.update_inst()

    def delete_angles(self):
        rows = self.view.delete_angles()

        if len(rows) > 0:
            self.model.delete_angles(rows)

        self.visualize()
        self.update_peaks()

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
        self.update_peaks()

    def add_orientation_process(self, progress):
        angles = self.view.get_angles()
        free_angles = self.view.get_free_angles()
        all_angles = self.view.get_all_angles()

        wavelength = self.view.get_wavelength()
        d_min = self.view.get_d_min()
        rows = self.view.get_number_of_orientations()

        if len(angles) > 0:
            progress("Calculating reflections", 5)

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

        if stats is not None and self.model.has_UB():
            self.view.plot_statistics(*stats)

            peak_dict = self.model.get_coverage_info(
                point_group, lattice_centering
            )
            if peak_dict is not None:
                peak_dict["axis_limit"] = self.view.get_d_min()

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
            self.update_peaks()

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

    def update_plan(self):
        instrument = self.view.get_instrument()
        mode = self.view.get_mode()
        settings = self.view.get_all_settings()
        comments = self.view.get_all_comments()
        counts = self.view.get_all_countings()
        values = self.view.get_all_values()
        use = self.view.get_orientations_to_use()
        names = self.view.get_free_angles()
        UB = self.model.get_UB()
        wavelength = self.view.get_wavelength()
        d_min = self.view.get_d_min()
        crysal_system = self.view.get_crystal_system()
        point_group = self.view.get_point_group()
        lattice_centering = self.view.get_lattice_centering()
        motors = self.view.get_motors()
        limits = self.view.get_goniometer_limits()
        self.model.create_plan(names, settings, comments, counts, values, use)
        self.model.create_sample(instrument, mode, UB, wavelength, d_min)
        self.model.update_sample(crysal_system, point_group, lattice_centering)
        self.model.update_goniometer_motors(limits, motors)

    def save_CSV(self):
        filename = self.view.save_CSV_file_dialog()

        if filename:
            self.update_plan()
            self.model.save_plan(filename)

    def save_experiment(self):
        filename = self.view.save_experiment_file_dialog()

        if filename:
            self.update_plan()
            self.model.save_experiment(filename)

    def load_experiment(self):
        filename = self.view.load_experiment_file_dialog()

        if filename:
            plan, config, symm = self.model.load_experiment(filename)

            settings, comments, counts, values, use = plan
            instrument, mode, wl, d_min, lims, vals = config
            cs, pg, lc = symm

            self.view.set_instrument(instrument)
            self.switch_instrument()
            self.view.set_mode(mode)
            self.update_oriented_lattice()
            self.view.set_transform(self.model.get_transform())
            self.view.set_wavelength(wl)
            self.view.set_d_min(d_min)
            self.view.set_goniometer_limits(lims)
            self.view.set_motors(vals)
            self.view.set_crystal_system(cs)
            self.switch_crystal()
            self.view.set_point_group(pg)
            self.switch_group()
            self.view.set_lattice_centering(lc)
            self.view.add_settings(settings, comments, counts, values, use)
            self.add_settings()

    def add_settings(self):
        worker = self.view.worker(self.add_settings_process)
        worker.connect_result(self.add_settings_complete)
        worker.connect_finished(self.visualize)
        worker.connect_progress(self.update_processing)

        self.view.start_worker_pool(worker)

    def add_settings_complete(self, result):
        if result is not None:
            self.update_peaks()

    def add_settings_process(self, progress):
        wavelength = self.view.get_wavelength()
        d_min = self.view.get_d_min()
        rows = self.view.get_number_of_orientations()

        instrument = self.view.get_instrument()
        mode = self.view.get_mode()
        axes, polarities = self.model.get_axes_polarities(instrument, mode)
        self.model.generate_axes(axes, polarities)
        limits = self.view.get_goniometer_limits()

        progress("Initializing instrument", 5)

        self.create_instrument()

        for row in range(rows):
            progress("Calculating settings", 90 // rows * (row + 1) + 5)

            angles = self.view.get_angle_setting(row)

            setting = self.model.get_setting(angles, limits)

            self.model.add_orientation(setting, wavelength, d_min, row)

        progress("Settings calculated!", 0)

        return rows
