from NeuXtalViz.presenters.base_presenter import NeuXtalVizPresenter


class VolumeSlicer(NeuXtalVizPresenter):
    def __init__(self, view, model):
        super(VolumeSlicer, self).__init__(view, model)

        self.view.connect_load_NXS(self.load_NXS)

        self.view.connect_slice_combo(self.redraw_data)
        self.view.connect_cut_combo(self.update_cut)

        self.view.connect_slice_thickness_line(self.update_slice)
        self.view.connect_cut_thickness_line(self.update_cut)

        self.view.connect_clim_combo(self.redraw_data)
        self.view.connect_cbar_combo(self.redraw_data)

        self.view.connect_min_slider(self.view.update_colorbar_min)
        self.view.connect_max_slider(self.view.update_colorbar_max)

        self.view.connect_vlim_combo(self.update_slice)

        self.view.connect_slice_scale_combo(self.update_slice)
        self.view.connect_cut_scale_combo(self.update_cut)

        self.view.connect_slice_line(self.redraw_data)
        self.view.connect_cut_line(self.update_cut)

        self.view.connect_slice_ready(self.update_slice)
        self.view.connect_cut_ready(self.update_cut)

        self.view.connect_vmin_line(self.update_cvals)
        self.view.connect_vmax_line(self.update_cvals)

        self.view.connect_xmin_line(self.update_lims)
        self.view.connect_xmax_line(self.update_lims)

        self.view.connect_ymin_line(self.update_lims)
        self.view.connect_ymax_line(self.update_lims)

        self.view.connect_vol_scale_combo(self.redraw_data)
        self.view.connect_opacity_combo(self.redraw_data)
        self.view.connect_range_combo(self.redraw_data)

        self.view.connect_save_slice(self.save_slice)
        self.view.connect_save_cut(self.save_cut)

    def update_lims(self):
        xmin = self.view.get_xmin_value()
        xmax = self.view.get_xmax_value()
        ymin = self.view.get_ymin_value()
        ymax = self.view.get_ymax_value()

        if (
            xmin is not None
            and xmax is not None
            and ymin is not None
            and ymax is not None
        ):
            if xmin < xmax and ymin < ymax:
                xlim = [xmin, xmax]
                ylim = [ymin, ymax]
                self.view.set_slice_lim(xlim, ylim)
                line_cut = self.view.get_cut()
                lim = xlim if line_cut == "Axis 1" else ylim
                self.view.set_cut_lim(lim)

    def update_cvals(self):
        vmin = self.view.get_vmin_value()
        vmax = self.view.get_vmax_value()
        if vmin is not None and vmax is not None:
            if vmin < vmax:
                if vmin <= 0 and self.view.get_slice_scale() == "log":
                    vmin = vmax / 10
                self.view.update_colorbar_vlims(vmin, vmax)

    def update_slice_value(self):
        self.view.update_slice_value()

        self.update_slice()

    def update_cut_value(self):
        self.view.update_cut_value()

        self.update_cut()

    def update_slice(self):
        if self.model.is_histo_loaded():
            self.slice_data()

    def update_cut(self):
        if self.model.is_histo_loaded():
            self.cut_data()

    def load_NXS(self):
        filename = self.view.load_NXS_file_dialog()

        if filename:
            self.nxs_file = filename
            worker = self.view.worker(self.load_NXS_process)
            worker.connect_result(self.load_NXS_complete)
            worker.connect_finished(self.redraw_data)
            worker.connect_progress(self.update_processing)

            self.view.start_worker_pool(worker)

    def load_NXS_complete(self, result):
        self.update_oriented_lattice()

    def load_NXS_process(self, progress):
        progress("Processing...", 1)

        progress("Loading NeXus file...", 10)

        self.model.load_md_histo_workspace(self.nxs_file)

        progress("Loading NeXus file...", 50)

        progress("Loading NeXus file...", 80)

        progress("NeXus file loaded!", 100)

    def get_normal(self):
        slice_plane = self.view.get_slice()

        if slice_plane == "Axis 1/2":
            norm = [0, 0, 1]
        elif slice_plane == "Axis 1/3":
            norm = [0, 1, 0]
        else:
            norm = [1, 0, 0]

        return norm

    def get_axis(self):
        axis = [1 if not norm else 0 for norm in self.get_normal()]
        ind = [i for i, ax in enumerate(axis) if ax == 1]

        line_cut = self.view.get_cut()

        if line_cut == "Axis 1":
            axis[ind[0]] = 0
        else:
            axis[ind[1]] = 0

        return axis

    def get_clim_method(self):
        ctype = self.view.get_clim_clip_type()

        if ctype == "μ±3×σ":
            method = "normal"
        elif ctype == "Q₃/Q₁±1.5×IQR":
            method = "boxplot"
        else:
            method = None

        return method

    def get_vlim_method(self):
        ctype = self.view.get_vlim_clip_type()

        if ctype == "μ±3×σ":
            method = "normal"
        elif ctype == "Q₃/Q₁±1.5×IQR":
            method = "boxplot"
        else:
            method = None

        return method

    def redraw_data(self):
        worker = self.view.worker(self.redraw_data_process)
        worker.connect_result(self.redraw_data_complete)
        worker.connect_finished(self.slice_data)
        worker.connect_progress(self.update_processing)

        self.view.start_worker_pool(worker)

    def redraw_data_complete(self, result):
        if result is not None:
            histo, normal, norm, value, trans = result

            self.view.add_histo(histo, normal, norm, value)

            self.view.set_transform(trans)

    def redraw_data_process(self, progress):
        if self.model.is_histo_loaded():
            progress("Processing...", 1)

            progress("Updating volume...", 20)

            norm = self.get_normal()

            histo = self.model.get_histo_info(norm)

            data = histo["signal"]

            data = self.model.calculate_clim(data, self.get_clim_method())

            progress("Updating volume...", 50)

            histo["signal"] = data

            value = self.view.get_slice_value()

            normal = -self.model.get_normal_plane(norm)

            # origin = self.model.get_normal('[hkl]', orig)

            if value is not None:
                progress("Volume drawn!", 100)

                return histo, normal, norm, value, self.model.get_transform()

            else:
                progress("Invalid parameters.", 0)

    def slice_data(self):
        worker = self.view.worker(self.slice_data_process)
        worker.connect_result(self.slice_data_complete)
        worker.connect_finished(self.cut_data)
        worker.connect_progress(self.update_processing)

        self.view.start_worker_pool(worker)

    def slice_data_complete(self, result):
        if result is not None:
            self.view.add_slice(result)

    def slice_data_process(self, progress):
        if self.model.is_histo_loaded():
            norm = self.get_normal()

            thick = self.view.get_slice_thickness()
            value = self.view.get_slice_value()

            if thick is not None:
                progress("Processing...", 1)

                progress("Updating slice...", 50)

                slice_histo = self.model.get_slice_info(norm, value, thick)

                data = slice_histo["signal"]

                data = self.model.calculate_clim(data, self.get_vlim_method())

                slice_histo["signal"] = data

                progress("Data sliced!", 100)

                return slice_histo

    def cut_data(self):
        worker = self.view.worker(self.cut_data_process)
        worker.connect_result(self.cut_data_complete)
        worker.connect_finished(self.update_complete)
        worker.connect_progress(self.update_processing)

        self.view.start_worker_pool(worker)

    def cut_data_complete(self, result):
        if result is not None:
            self.view.add_cut(result)

    def cut_data_process(self, progress):
        if self.model.is_sliced():
            value = self.view.get_cut_value()
            thick = self.view.get_cut_thickness()

            axis = self.get_axis()

            if value is not None and thick is not None:
                progress("Processing...", 1)

                progress("Updating cut...", 50)

                progress("Data cut!", 100)

                cut_histo = self.model.get_cut_info(axis, value, thick)

                return cut_histo

    def save_slice(self):
        if self.model.is_sliced():
            filename = self.view.save_file_dialog()

            if filename:
                self.model.save_slice(filename)

    def save_cut(self):
        if self.model.is_cut():
            filename = self.view.save_file_dialog()

            if filename:
                self.model.save_cut(filename)
