import pyvista as pv
from nova.trame.view.components import InputField
from nova.trame.view.layouts import GridLayout, VBoxLayout
from pyvista.trame.ui import get_viewer
from trame.widgets import html
from trame.widgets import vuetify3 as vuetify


class VisualizationPanel:
    def __init__(self, server, view_model):
        self.server = server
        self.view_model = view_model
        self.view_model.controls_bind.connect("controls")
        self.view_model.show_axes_bind.connect(self.show_axes)
        self.view_model.parallel_projection_bind.connect(self.change_projection)
        self.view_model.progress_bind.connect("progress")
        self.view_model.status_bind.connect("status")
        self.view_model.update_processing("Ready!", 0)

        # These are just for compatibility with the view model as it communicates with Qt slightly differently
        self.view_model.lattice_parameters_bind.connect("lattice_parameters")

        self.camera_position = None
        self.plotter = self.create_plotter()

        self.create_ui()

    @property
    def state(self):
        return self.server.state

    def create_plotter(self):
        plotter = pv.Plotter(off_screen=True)
        plotter.background_color = "#f0f0f0"

        return plotter

    def create_ui(self):
        with vuetify.VContainer(classes="mr-2 pa-0", fluid=True):
            with GridLayout(columns=5):
                with VBoxLayout(valign="space-between"):
                    vuetify.VBtn("Save Screenshot")
                    vuetify.VBtn("Reset View", click=self.reset_view)
                    vuetify.VBtn("Reset Camera", click=self.reset_camera)

                with VBoxLayout(column_span=3):
                    with vuetify.VTabs(
                        v_model="controls.camera_tab",
                        classes="mb-1",
                        update_modelValue="flushState('controls')",
                    ):
                        vuetify.VTab("Direction View", value=1)
                        vuetify.VTab("Manual View", value=2)
                    with vuetify.VWindow(v_model="controls.camera_tab"):
                        with vuetify.VWindowItem(value=1):
                            with GridLayout(columns=6):
                                vuetify.VBtn("+Qx", click=self.plotter.view_yz)
                                vuetify.VBtn("+Qy", click=self.plotter.view_zx)
                                vuetify.VBtn("+Qz", click=self.plotter.view_xy)
                                vuetify.VBtn("a*")
                                vuetify.VBtn("b*")
                                vuetify.VBtn("c*")
                                vuetify.VBtn("-Qx", click=self.plotter.view_zy)
                                vuetify.VBtn("-Qy", click=self.plotter.view_xz)
                                vuetify.VBtn("-Qz", click=self.plotter.view_yx)
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

            with vuetify.VSheet(classes="mb-2", style="height: 50vh;"):
                self.view = get_viewer(self.plotter)
                self.view.ui(add_menu=False, mode="server")

            with vuetify.VProgressLinear(
                v_model="progress",
                height=20,
                striped=("progress < 100",),
            ):
                html.Span("{{ status }}")

    def clear_scene(self):
        self.plotter.clear_plane_widgets()
        self.plotter.clear_actors()

        if self.camera_position is not None:
            self.camera_position = self.plotter.camera_position

    def reset_camera(self):
        """
        Reset the camera.

        """

        self.plotter.reset_camera()

    def reset_scene(self):
        if self.camera_position is not None:
            self.plotter.camera_position = self.camera_position
        else:
            self.reset_view()

    def reset_view(self, negative=False):
        """
        Reset the view.

        """

        self.plotter.reset_camera()
        self.plotter.view_isometric(negative)
        self.camera_position = self.plotter.camera_position

    def show_axes(self, data):
        T, reciprocal_lattice, show_axes = data

        if not show_axes:
            self.plotter.hide_axes()
        elif T is not None:
            t = pv._vtk.vtkMatrix4x4()
            for i in range(3):
                for j in range(3):
                    t.SetElement(i, j, T[i, j])
            if reciprocal_lattice:
                actor = self.plotter.add_axes(xlabel="a*", ylabel="b*", zlabel="c*")
            else:
                actor = self.plotter.add_axes(xlabel="a", ylabel="b", zlabel="c")
            actor.SetUserMatrix(t)
        self.view.update()

    def change_projection(self, parallel_projection):
        """
        Enable or disable parallel projection.

        """

        if parallel_projection:
            self.plotter.enable_parallel_projection()
        else:
            self.plotter.disable_parallel_projection()
        self.view.update()
