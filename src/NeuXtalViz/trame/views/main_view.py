"""Main view for NeuXtalViz."""

from nova.mvvm.trame_binding import TrameBinding
from nova.trame import ThemedApp
from trame.widgets import vuetify3 as vuetify
from trame_server.core import Server
from trame_server.state import State

from NeuXtalViz.models.volume_slicer import VolumeSlicerModel
from NeuXtalViz.trame.views.volume_slicer import VolumeSlicerView
from NeuXtalViz.view_models.volume_slicer import VolumeSlicerViewModel


class NeuXtalViz(ThemedApp):
    """Main view for NeuXtalViz."""

    def __init__(self, server: Server = None) -> None:
        self.server = server
        super().__init__(server=server)

        binding = TrameBinding(self.server.state)

        self.view_model = VolumeSlicerViewModel(VolumeSlicerModel(), binding)

        self.create_ui()

    @property
    def state(self) -> State:
        return self.server.state

    def create_ui(self) -> None:
        self.state.trame__title = "NeuXtalViz"
        self.set_theme("CompactTheme")

        with super().create_ui() as layout:
            layout.toolbar_title.set_text("NeuXtalViz")

            with layout.pre_content:
                with vuetify.VTabs(classes="pl-4", density="compact"):
                    vuetify.VTab("Volume Slicer")

            with layout.content:
                VolumeSlicerView(self.server, self.view_model)
