"""Main view for NeuXtalViz."""

from nova.mvvm.trame_binding import TrameBinding
from nova.trame import ThemedApp
from nova.trame.view.layouts import HBoxLayout
from trame.widgets import vuetify3 as vuetify
from trame_server.core import Server
from trame_server.state import State


class NeuXtalViz(ThemedApp):
    """Main view for NeuXtalViz."""

    def __init__(self, server: Server = None) -> None:
        self.server = server
        super().__init__(server=server)

        binding = TrameBinding(self.server.state)
        # self.vm = create_viewmodels(binding)
        # self.vm["volume_slicer"].config_bind.connect("config")

        # self.vm.update_view()

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
                with vuetify.VTabs(classes="pl-4"):
                    vuetify.VTab("Volume Slicer")

            with layout.content:
                with HBoxLayout(classes="pa-2"):
                    pass
