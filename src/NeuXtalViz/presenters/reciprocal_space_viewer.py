from NeuXtalViz.presenters.base_presenter import NeuXtalVizPresenter

class ReciprocalSpaceViewer(NeuXtalVizPresenter):

    def __init__(self, view, model):

        super(ReciprocalSpaceViewer, self).__init__(view, model)

        # self.model.set_peak_workspace('peaks')
        # self.view.add_peaks(self.model.get_peak_info())
        # self.view.set_transform(self.model.get_transform())