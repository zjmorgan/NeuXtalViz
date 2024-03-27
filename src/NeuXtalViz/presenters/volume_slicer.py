from NeuXtalViz.presenters.base_presenter import NeuXtalVizPresenter

class VolumeSlicer(NeuXtalVizPresenter):

    def __init__(self, view, model):

        super(VolumeSlicer, self).__init__(view, model)

        self.view.connect_load_NXS(self.load_NXS)

        self.view.connect_slice(self.slice_data)
        self.view.connect_cut(self.cut_data)
        self.view.connect_redraw(self.redraw_data)

    def load_NXS(self):

        filename = self.view.load_NXS_file_dialog()

        if filename:

            self.model.load_md_histo_workspace(filename)

            self.redraw_data()

    def get_normal(self):

        slice_plane = self.view.get_slice()

        if slice_plane == 'Axis 1/2':
            norm = [0,0,1]
        elif slice_plane == 'Axis 1/3':
            norm = [0,1,0]
        else:
            norm = [1,0,0]

        return norm

    def get_clim_method(self):

        ctype = self.view.get_clim_clip_type()

        if ctype == 'μ±3×σ':
            method = 'normal'
        elif ctype == 'Q₃/Q₁±1.5×IQR':
            method = 'boxplot'
        else:
            method = None

        return method

    def redraw_data(self):

        if self.model.is_histo_loaded():

            histo = self.model.get_histo_info()

            data = histo['signal']

            data = self.model.calculate_clim(data, self.get_clim_method())

            histo['signal'] = data

            value = self.view.get_slice_value()

            norm = self.get_normal()

            normal = -self.model.get_normal('[uvw]', norm)

            if value is not None:

                self.view.add_histo(histo, normal, value)
                self.view.set_transform(self.model.get_transform())

                self.slice_data()

    def slice_data(self):

        if self.model.is_histo_loaded():

            norm = self.get_normal()

            thick = self.view.get_slice_thickness()
            value = self.view.get_slice_value()

            if thick is not None:

                slice_histo = self.model.get_slice_info(norm, value, thick)

                data = slice_histo['signal']

                data = self.model.calculate_clim(data, self.get_clim_method())

                slice_histo['signal'] = data

                self.view.add_slice(slice_histo)

            self.cut_data()

    def cut_data(self):

        if self.model.is_sliced():

            value = self.view.get_cut_value()
            line_cut = self.view.get_cut()

            axis = [1 if not norm else 0 for norm in self.model.normal]
            ind = [i for i, ax in enumerate(axis) if ax == 1]

            if line_cut == 'Axis 1':
                axis[ind[0]] = 0
            else:
                axis[ind[1]] = 0

            thick = self.view.get_cut_thickness()

            if value is not None and thick is not None:

                cut_histo = self.model.get_cut_info(axis, value, thick)
                self.view.add_cut(cut_histo)
