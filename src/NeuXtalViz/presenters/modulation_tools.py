from NeuXtalViz.presenters.base_presenter import NeuXtalVizPresenter

class Modulation(NeuXtalVizPresenter):

    def __init__(self, view, model):

        super(Modulation, self).__init__(view, model)

        self.view.connect_cluster(self.cluster)
        self.view.connect_load_UB(self.load_UB)
        self.view.connect_load_peaks(self.load_peaks)

    def cluster(self):

        params = self.view.get_cluster_parameters()

        if params is not None:

            peak_info = self.model.get_peak_info()
            if peak_info is not None:
                success = self.model.cluster_peaks(peak_info, *params)
                if success:
                    self.view.add_peaks(peak_info)
                    self.view.update_table(peak_info)

    def load_UB(self):

        filename = self.view.load_UB_file_dialog()

        if filename:

            self.model.load_UB(filename)

            self.update_oriented_lattice()

            self.view.set_transform(self.model.get_transform())

    def load_peaks(self):

        filename = self.view.load_peaks_file_dialog()

        if filename:

            self.model.load_peaks(filename)