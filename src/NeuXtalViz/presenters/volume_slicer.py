from NeuXtalViz.presenters.base_presenter import NeuXtalVizPresenter

class VolumeSlicer(NeuXtalVizPresenter):

    def __init__(self, view, model):

        super(VolumeSlicer, self).__init__(view, model)

        self.view.connect_load_NXS(self.load_NXS)

        self.view.connect_slice_combo(self.redraw_data)
        self.view.connect_cut_combo(self.redraw_data)

        self.view.connect_slice_thickness_line(self.redraw_data)
        self.view.connect_cut_thickness_line(self.redraw_data)

        self.view.connect_clim_combo(self.redraw_data)
        self.view.connect_cbar_combo(self.redraw_data)

        self.view.connect_min_slider(self.view.update_colorbar_min)
        self.view.connect_max_slider(self.view.update_colorbar_max)

        self.view.connect_slice_slider(self.update_slice_value)
        self.view.connect_cut_slider(self.update_cut_value)

        self.view.connect_slice_scale_combo(self.update_slice)
        self.view.connect_cut_scale_combo(self.update_cut)

        self.view.connect_slice_line(self.redraw_data)
        self.view.connect_cut_line(self.redraw_data)

    def update_slice_value(self):

        self.view.update_slice_value()

        self.update_slice()

    def update_cut_value(self):

        self.view.update_cut_value()
        
        self.update_cut()

    def update_slice(self):

        if self.model.is_histo_loaded():

            self.view.update_slice_slider()

            self.slice_data()

    def update_cut(self):

        if self.model.is_histo_loaded():

            self.view.update_cut_slider()

            self.cut_data()

    def load_NXS(self):

        filename = self.view.load_NXS_file_dialog()

        if filename:

            self.update_processing()

            self.update_processing('Loading NeXus file...', 10) 

            self.model.load_md_histo_workspace(filename)

            self.update_processing('Loading NeXus file...', 50) 

            self.update_oriented_lattice()

            self.update_processing('Loading NeXus file...', 80) 

            self.redraw_data()

            self.update_complete('NeXus file loaded!')

        else:

            self.update_invalid()

    def get_normal(self):

        slice_plane = self.view.get_slice()

        if slice_plane == 'Axis 1/2':
            norm = [0,0,1]
        elif slice_plane == 'Axis 1/3':
            norm = [0,1,0]
        else:
            norm = [1,0,0]

        return norm

    def get_axis(self):

        axis = [1 if not norm else 0 for norm in self.get_normal()]
        ind = [i for i, ax in enumerate(axis) if ax == 1]

        line_cut = self.view.get_cut()

        if line_cut == 'Axis 1':
            axis[ind[0]] = 0
        else:
            axis[ind[1]] = 0

        return axis

    def get_clim_method(self):

        ctype = self.view.get_clim_clip_type()

        if ctype == 'μ±3×σ':
            method = 'normal'
        elif ctype == 'Q₃/Q₁±1.5×IQR':
            method = 'boxplot'
        else:
            method = None

        return method

    def get_slice_indices(self):

        norm = self.get_normal()

        sind = norm.index(1)

        return sind

    def get_cut_indices(self):

        axis = self.get_axis()

        cind = axis.index(1)

        return cind

    def update_slice_info(self):

        sind = self.get_slice_indices()

        self.view.update_slice_limits(self.model.shape[sind],
                                      self.model.spacing[sind],
                                      self.model.min_lim[sind])

    def update_cut_info(self):

        cind = self.get_cut_indices()

        self.view.update_cut_limits(self.model.shape[cind],
                                    self.model.spacing[cind],
                                    self.model.min_lim[cind])

    def redraw_data(self):

        if self.model.is_histo_loaded():

            self.update_processing()

            self.update_processing('Updating volume...', 20)

            histo = self.model.get_histo_info()

            data = histo['signal']

            data = self.model.calculate_clim(data, self.get_clim_method())

            self.update_processing('Updating volume...', 50)

            histo['signal'] = data

            value = self.view.get_slice_value()

            norm = self.get_normal()

            normal = -self.model.get_normal('[uvw]', norm)

            origin = norm
            origin[origin.index(1)] = value

            # origin = self.model.get_normal('[hkl]', orig)

            if value is not None:

                self.view.add_histo(histo, normal, origin)
             
                self.view.set_transform(self.model.get_transform())

                self.slice_data()

                self.update_complete('Volume drawn!')

            else:

                self.update_invalid()

    def slice_data(self):

        if self.model.is_histo_loaded():

            norm = self.get_normal()

            thick = self.view.get_slice_thickness()
            value = self.view.get_slice_value()

            if thick is not None:

                self.update_processing()

                self.update_processing('Updating slice...', 50)

                slice_histo = self.model.get_slice_info(norm, value, thick)

                data = slice_histo['signal']

                data = self.model.calculate_clim(data, self.get_clim_method())

                slice_histo['signal'] = data

                self.view.add_slice(slice_histo)

                self.update_complete('Data sliced!')

                self.update_slice_info()

                self.cut_data()

    def cut_data(self):

        if self.model.is_sliced():

            value = self.view.get_cut_value()
            thick = self.view.get_cut_thickness()

            axis = self.get_axis()

            if value is not None and thick is not None:

                self.update_processing()

                self.update_processing('Updating cut...', 50)
    
                cut_histo = self.model.get_cut_info(axis, value, thick)

                self.view.add_cut(cut_histo)

                self.update_complete('Data cut!')

                self.update_cut_info()
