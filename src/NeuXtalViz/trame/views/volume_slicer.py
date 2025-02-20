from nova.trame.view.layouts import GridLayout

from NeuXtalViz.trame.views.components.visualization_panel import VisualizationPanel


class VolumeSlicerView:
    def __init__(self, view_model):
        self.view_model = view_model

        self.create_ui()



    def create_ui(self):
        with GridLayout(classes="pa-2", columns=2):
            self.base_view = VisualizationPanel(self.view_model)

            # TODO: remove this temporary hack
            from trame.widgets import html

            html.Div()
