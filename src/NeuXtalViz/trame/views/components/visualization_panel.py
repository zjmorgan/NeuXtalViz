from nova.trame.view.components import InputField
from nova.trame.view.layouts import GridLayout, VBoxLayout
from pyvista import Plotter
from pyvista.trame.ui import get_viewer
from trame.widgets import html
from trame.widgets import vuetify3 as vuetify


class VisualizationPanel:
    def __init__(self, view_model):
        self.view_model = view_model
        self.view_model.controls_bind.connect("controls")

        self.plotter = self.create_plotter()

        self.create_ui()

    def create_plotter(self):
        plotter = Plotter(off_screen=True)
        plotter.background_color = "#f0f0f0"

        return plotter

    def create_ui(self):
        with vuetify.VContainer(classes="mr-2 pa-0", fluid=True):
            with GridLayout(columns=5):
                with VBoxLayout(valign="space-between"):
                    vuetify.VBtn("Save Screenshot")
                    vuetify.VBtn("Reset View")
                    vuetify.VBtn("Reset Camera")

                with VBoxLayout(column_span=3):
                    with vuetify.VTabs(v_model="controls.view_tab", classes="mb-1", update_modelValue="flushState('controls')"):
                        vuetify.VTab("Direction View", value=1)
                        vuetify.VTab("Manual View", value=2)
                    with vuetify.VWindow(v_model="controls.view_tab"):
                        with vuetify.VWindowItem(value=1):
                            with GridLayout(columns=6):
                                vuetify.VBtn("+Qx")
                                vuetify.VBtn("+Qy")
                                vuetify.VBtn("+Qz")
                                vuetify.VBtn("a*")
                                vuetify.VBtn("b*")
                                vuetify.VBtn("c*")
                                vuetify.VBtn("-Qx")
                                vuetify.VBtn("-Qy")
                                vuetify.VBtn("-Qz")
                                vuetify.VBtn("a")
                                vuetify.VBtn("b")
                                vuetify.VBtn("c")
                        with vuetify.VWindowItem(value=2):
                            with GridLayout(columns=10, halign="center"):
                                vuetify.VLabel("h")
                                vuetify.VLabel("k")
                                vuetify.VLabel("l")
                                InputField(v_model="controls.manual_axis_type", column_span=2, type="select")
                                vuetify.VLabel("h")
                                vuetify.VLabel("k")
                                vuetify.VLabel("l")
                                InputField(v_model="controls.manual_up_axis_type", column_span=2, type="select")
                                InputField(v_model="controls.manual_axes[0]")
                                InputField(v_model="controls.manual_axes[1]")
                                InputField(v_model="controls.manual_axes[2]")
                                vuetify.VBtn("View Axis", column_span=2)
                                InputField(v_model="controls.manual_up_axes[0]")
                                InputField(v_model="controls.manual_up_axes[1]")
                                InputField(v_model="controls.manual_up_axes[2]")
                                vuetify.VBtn("View Up Axis", column_span=2)

                with VBoxLayout(valign="space-between"):
                    InputField(v_model="controls.reciprocal_lattice", type="checkbox")
                    InputField(v_model="controls.show_axes", type="checkbox")
                    InputField(v_model="controls.parallel_projection", type="checkbox")

            with vuetify.VSheet(classes="mb-2", style="height: 70vh;"):
                view = get_viewer(self.plotter)
                view.ui(add_menu=False, mode="server")

            with vuetify.VProgressLinear(v_model="controls.progress", height=25, striped=True):
                html.Span("{{ controls.status }}")
