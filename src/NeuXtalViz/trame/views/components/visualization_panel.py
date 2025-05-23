from asyncio import ensure_future, sleep

import numpy as np
import pyvista as pv
from io import BytesIO
from nova.trame.view.components import InputField
from nova.trame.view.layouts import GridLayout, HBoxLayout, VBoxLayout
from pyvista.trame.ui import get_viewer
from trame.widgets import client, html
from trame.widgets import vuetify3 as vuetify


class VisualizationPanel:
    def __init__(self, server, view_model):
        self.server = server
        self.view_model = view_model
        self.view_model.controls_bind.connect("controls")
        self.view_model.lattice_parameters_bind.connect("lattice_parameters")
        self.view_model.show_axes_bind.connect(self.trigger_show_axes)
        self.view_model.parallel_projection_bind.connect(self.change_projection)
        self.view_model.progress_bind.connect("progress")
        self.view_model.status_bind.connect("status")
        self.view_model.up_vector_bind.connect(self.view_up_vector)
        self.view_model.update_view_bind.connect(self.update_view)
        self.view_model.vector_bind.connect(self.view_vector)

        self.view_model.update_processing("Ready!", 0)

        self.axis_data = None
        ensure_future(self.show_axes_loop())

        self.camera_position = None
        self.plotter = self.create_plotter()
        self.js_download = client.JSEval(
            exec="utils.download($event[0], $event[1], 'image/png')"
        ).exec

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
                with VBoxLayout(valign="start"):
                    vuetify.VBtn("Save Screenshot", click=self.save_screenshot)
                    vuetify.VBtn("Reset View", classes="my-1", click=self.reset_view)
                    vuetify.VBtn("Reset Camera", click=self.reset_camera)

                with VBoxLayout(column_span=3):
                    with vuetify.VTabs(
                        v_model="controls.camera_tab",
                        classes="mb-1",
                        density="compact",
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
                                vuetify.VBtn("a*", click=self.view_model.view_bc_star)
                                vuetify.VBtn("b*", click=self.view_model.view_ca_star)
                                vuetify.VBtn("c*", click=self.view_model.view_ab_star)
                                vuetify.VBtn("-Qx", click=self.plotter.view_zy)
                                vuetify.VBtn("-Qy", click=self.plotter.view_xz)
                                vuetify.VBtn("-Qz", click=self.plotter.view_yx)
                                vuetify.VBtn("a", click=self.view_model.view_bc)
                                vuetify.VBtn("b", click=self.view_model.view_ca)
                                vuetify.VBtn("c", click=self.view_model.view_ab)
                        with vuetify.VWindowItem(value=2):
                            with GridLayout(columns=12, halign="center"):
                                vuetify.VLabel("{{ controls.manual_axis_type[0] }}")
                                vuetify.VLabel("{{ controls.manual_axis_type[1] }}")
                                vuetify.VLabel("{{ controls.manual_axis_type[2] }}")
                                InputField(
                                    v_model="controls.manual_axis_type",
                                    column_span=3,
                                    type="select",
                                )
                                vuetify.VLabel("{{ controls.manual_up_axis_type[0] }}")
                                vuetify.VLabel("{{ controls.manual_up_axis_type[1] }}")
                                vuetify.VLabel("{{ controls.manual_up_axis_type[2] }}")
                                InputField(
                                    v_model="controls.manual_up_axis_type",
                                    column_span=3,
                                    type="select",
                                )
                                InputField(v_model="controls.manual_axes[0]")
                                InputField(v_model="controls.manual_axes[1]")
                                InputField(v_model="controls.manual_axes[2]")
                                vuetify.VBtn(
                                    "View Axis",
                                    column_span=3,
                                    click=self.view_model.view_manual,
                                )
                                InputField(v_model="controls.manual_up_axes[0]")
                                InputField(v_model="controls.manual_up_axes[1]")
                                InputField(v_model="controls.manual_up_axes[2]")
                                vuetify.VBtn(
                                    "View Up Axis",
                                    column_span=3,
                                    click=self.view_model.view_up_manual,
                                )

                with VBoxLayout(valign="start"):
                    InputField(
                        v_model="controls.reciprocal_lattice",
                        density="compact",
                        type="checkbox",
                    )
                    InputField(
                        v_model="controls.show_axes", density="compact", type="checkbox"
                    )
                    InputField(
                        v_model="controls.parallel_projection",
                        density="compact",
                        type="checkbox",
                    )

            with vuetify.VSheet(classes="mb-2", style="height: 40vh;"):
                self.view = get_viewer(self.plotter)
                self.view.ui(add_menu=False, mode="server")

            with vuetify.VTabs(
                v_model="controls.oriented_lattice_tab",
                classes="mb-1",
                density="compact",
                update_modelValue="flushState('controls')",
            ):
                vuetify.VTab("Lattice Parameters", value=1)
                vuetify.VTab("Sample Orientation", value=2)
            with vuetify.VWindow(v_model="controls.oriented_lattice_tab"):
                with vuetify.VWindowItem(value=1):
                    with GridLayout(columns=3):
                        InputField(v_model="lattice_parameters.a", readonly=True)
                        InputField(v_model="lattice_parameters.b", readonly=True)
                        with HBoxLayout():
                            InputField(v_model="lattice_parameters.c", readonly=True)
                            vuetify.VLabel("Å")
                        InputField(v_model="lattice_parameters.alpha", readonly=True)
                        InputField(v_model="lattice_parameters.beta", readonly=True)
                        with HBoxLayout():
                            InputField(
                                v_model="lattice_parameters.gamma", readonly=True
                            )
                            vuetify.VLabel("°")
                with vuetify.VWindowItem(value=2):
                    with HBoxLayout():
                        vuetify.VLabel("u:")
                        InputField(v_model="lattice_parameters.u[0]", readonly=True)
                        InputField(v_model="lattice_parameters.u[1]", readonly=True)
                        InputField(v_model="lattice_parameters.u[2]", readonly=True)
                    with HBoxLayout():
                        vuetify.VLabel("v:")
                        InputField(v_model="lattice_parameters.v[0]", readonly=True)
                        InputField(v_model="lattice_parameters.v[1]", readonly=True)
                        InputField(v_model="lattice_parameters.v[2]", readonly=True)

            with vuetify.VProgressLinear(
                v_model="progress",
                height=20,
                striped=("progress < 100",),
            ):
                html.Span("{{ status }}")

    def update_view(self, _):
        self.view.update()

    def save_screenshot(self):
        data = BytesIO()
        self.plotter.screenshot(data)
        data.seek(0)
        self.js_download(("neuxtalviz.png", data.read()))

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
        self.update_view(None)

    async def show_axes_loop(self):
        while True:
            if self.axis_data is not None:
                self.show_axes(self.axis_data)
                self.axis_data = None

            await sleep(0.1)

    def trigger_show_axes(self, data):
        self.axis_data = data

    def change_projection(self, parallel_projection):
        """
        Enable or disable parallel projection.

        """

        if parallel_projection:
            self.plotter.enable_parallel_projection()
        else:
            self.plotter.disable_parallel_projection()

    def view_vector(self, vecs):
        if len(vecs) == 2:
            vec = np.cross(vecs[0], vecs[1])
            self.plotter.view_vector(vecs[0], vec)
        else:
            self.plotter.view_vector(vecs)

    def view_up_vector(self, vec):
        self.plotter.set_viewup(vec)
