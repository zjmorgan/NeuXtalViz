from NeuXtalViz.presenters.base_presenter import NeuXtalVizPresenter


class Modulation(NeuXtalVizPresenter):
    def __init__(self, view, model):
        super(Modulation, self).__init__(view, model)

        self.view.connect_cluster(self.cluster)
        self.view.connect_load_UB(self.load_UB)
        self.view.connect_load_peaks(self.load_peaks)

    def cluster(self):
        worker = self.view.worker(self.cluster_process)
        worker.connect_result(self.cluster_complete)
        worker.connect_progress(self.update_processing)

        self.view.start_worker_pool(worker)

    def cluster_complete(self, result):
        if result is not None:
            self.update_processing("Adding peaks.", 30)
            self.view.add_peaks(result)
            self.view.update_table(result)
            self.update_processing("Peaks added!", 0)

    def cluster_process(self, progress):
        params = self.view.get_cluster_parameters()

        if params is not None:
            progress("Invalid parameters.", 0)

            peak_info = self.model.get_peak_info()
            if peak_info is not None:
                progress("Clustering peaks.", 25)

                success = self.model.cluster_peaks(peak_info, *params)

                if success:
                    progress("Peaks clustered!", 100)

                    return peak_info

                else:
                    progress("Invalid cluster.", 0)

        else:
            progress("Invalid parameters.", 0)

    def load_UB(self):
        filename = self.view.load_UB_file_dialog()

        if filename:
            self.update_processing("Loading peaks.", 30)

            self.model.load_UB(filename)

            self.update_oriented_lattice()

            self.view.set_transform(self.model.get_transform())

            self.update_processing("UB loaded!", 0)

    def load_peaks(self):
        filename = self.view.load_peaks_file_dialog()

        if filename:
            self.update_processing("Loading peaks.", 30)

            self.model.load_peaks(filename)

            self.update_processing("Peaks loaded!", 0)
