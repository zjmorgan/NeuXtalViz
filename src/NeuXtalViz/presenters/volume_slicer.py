from NeuXtalViz.presenters.base_presenter import NeuXtalVizPresenter

class VolumeSlicer(NeuXtalVizPresenter):

    def __init__(self, view, model):

        super(VolumeSlicer, self).__init__(view, model)

        self.view.connect_load_NXS(self.load_NXS)

        self.view.connect_slice(self.plot_data)

    def load_NXS(self):

        filename = self.view.load_NXS_file_dialog()

        if filename:

            self.model.load_md_histo_workspace(filename)

            self.plot_data()

    def plot_data(self):

        if self.model.is_histo_loaded():

            ctype = self.view.get_clim_clip_type()

            if ctype == 'μ±3×σ':
                method = 'normal'
            elif ctype == 'Q₃/Q₁±1.5×IQR':
                method = 'boxplot'
            else:
                method = None

            histo = self.model.get_histo_info()
            clim = self.model.calculate_clim(method)

            value = self.view.get_normal_value()
            slice_plane = self.view.get_slice()
    
            if slice_plane == 'Axis 1/2':
                normal = (0,0,-1)
            elif slice_plane == 'Axis 1/3':
                normal = (0,-1,0)
            else:
                normal = (-1,0,0)
    
            normal = self.model.get_normal('[uvw]', normal)
    
            if value is not None:
    
                self.view.add_histo(histo, clim, normal, value)
                self.view.set_transform(self.model.get_transform())