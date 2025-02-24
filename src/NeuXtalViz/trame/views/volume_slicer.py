from asyncio import create_task, sleep
from threading import Thread

import numpy as np
import pyvista as pv
from nova.trame.view.components import InputField
from nova.trame.view.layouts import GridLayout, HBoxLayout
from trame.widgets import html
from trame.widgets import vuetify3 as vuetify

from NeuXtalViz.trame.views.components.visualization_panel import VisualizationPanel


cmaps = {
    "Sequential": "viridis",
    "Binary": "binary",
    "Diverging": "bwr",
    "Rainbow": "turbo",
    "Modified": "modified",
}

opacities = {
    "Linear": {"Low->High": "linear", "High->Low": "linear_r"},
    "Geometric": {"Low->High": "geom", "High->Low": "geom_r"},
    "Sigmoid": {"Low->High": "sigmoid", "High->Low": "sigmoid_r"},
}


class VolumeSlicerView:
    def __init__(self, server, view_model):
        self.server = server
        self.view_model = view_model
        self.view_model.vs_controls_bind.connect("vs_controls")
        self.view_model.add_histo_bind.connect(self.add_histo)

        self.create_ui()

    @property
    def state(self):
        return self.server.state

    def create_ui(self):
        with GridLayout(classes="pa-2", columns=2):
            self.base_view = VisualizationPanel(self.server, self.view_model)

            with HBoxLayout(valign="start"):
                InputField(v_model="vs_controls.vol_scale", type="select")
                InputField(v_model="vs_controls.opacity", type="select")
                InputField(v_model="vs_controls.opacity_range", type="select")
                InputField(v_model="vs_controls.clim_clip_type", type="select")
                InputField(v_model="vs_controls.cbar", type="select")
                with html.Div():
                    # InputField(
                    #     ref="fread",
                    #     classes="d-none",
                    #     type="file",
                    #     __events=["change"],
                    #     change=(self.load_NXS, "[$event.target.files]"),
                    # )
                    # vuetify.VBtn("Load NXS", click="trame.refs.fread.click()")
                    vuetify.VBtn("Load NXS", click=self.load_NXS)

    def load_in_background(self):
        # TODO: File upload is not working well for me with the nexus file, so I've hardcoded the path for now.
        # We should switch to RemoteFileInput, anyways.
        self.view_model.load_NXS_process(
            "/home/dugganjw/mvvm_neux/tests/data/Ba3Co2O6_50K_proj_small_6_m.nxs"
        )
        self.view_model.load_NXS_complete()
        self.loading_data = False

    async def monitor_load(self):
        while self.loading_data:
            await sleep(0.1)

        self.view_model.update_processing("Loading NeXus file...", 50)
        self.view_model.update_processing("Loading NeXus file...", 80)
        self.view_model.update_processing("NeXus file loaded!", 100)
        self.redraw_data()

    def load_NXS(self):
        self.view_model.update_processing("Processing...", 1)
        self.view_model.update_processing("Loading NeXus file...", 10)

        self.loading_data = True
        thread = Thread(target=self.load_in_background, daemon=True)
        thread.start()

        create_task(self.monitor_load())

    def redraw_in_background(self):
        self.redrawing_result = self.view_model.redraw_data_process()
        self.view_model.redraw_data_complete(self.redrawing_result)
        self.redrawing = False

    async def monitor_redraw(self):
        while self.redrawing_result is None:
            await sleep(0.1)

        self.view_model.update_processing("Updating volume...", 50)

        while self.redrawing:
            await sleep(0.1)

        if self.redrawing_result is not None:
            self.view_model.update_processing("Volume drawn!", 100)
        else:
            self.view_model.update_processing("Invalid parameters.", 0)

        self.base_view.reset_view()
        # self.slice_data()

    def redraw_data(self):
        self.view_model.update_processing("Processing...", 1)
        self.view_model.update_processing("Updating volume...", 20)

        self.redrawing = True
        self.redrawing_result = None
        thread = Thread(target=self.redraw_in_background, daemon=True)
        thread.start()

        create_task(self.monitor_redraw())

    def add_histo(self, result):
        histo_dict, normal, norm, value = result

        opacity = opacities[self.view_model.get_opacity()][
            self.view_model.get_opacity_range()
        ]

        log_scale = True if self.view_model.get_vol_scale() == "Log" else False

        cmap = cmaps[self.view_model.get_colormap()]

        self.base_view.clear_scene()

        self.norm = np.array(norm).copy()
        origin = norm
        origin[origin.index(1)] = value

        signal = histo_dict["signal"]
        labels = histo_dict["labels"]

        min_lim = histo_dict["min_lim"]
        max_lim = histo_dict["max_lim"]
        spacing = histo_dict["spacing"]

        P = histo_dict["projection"]
        T = histo_dict["transform"]
        S = histo_dict["scales"]

        grid = pv.ImageData()

        grid.dimensions = np.array(signal.shape) + 1

        grid.origin = min_lim
        grid.spacing = spacing

        min_bnd = min_lim * S
        max_bnd = max_lim * S

        bounds = np.array([[min_bnd[i], max_bnd[i]] for i in [0, 1, 2]])
        limits = np.array([[min_lim[i], max_lim[i]] for i in [0, 1, 2]])

        a = pv._vtk.vtkMatrix3x3()
        b = pv._vtk.vtkMatrix4x4()
        for i in range(3):
            for j in range(3):
                a.SetElement(i, j, T[i, j])
                b.SetElement(i, j, P[i, j])

        grid.cell_data["scalars"] = signal.flatten(order="F")

        normal /= np.linalg.norm(normal)

        origin = np.dot(P, origin)

        clim = [np.nanmin(signal), np.nanmax(signal)]

        if not np.all(np.isfinite(clim)):
            clim = [0.1, 10]

        self.clip = self.base_view.plotter.add_volume_clip_plane(
            grid,
            opacity=opacity,
            log_scale=log_scale,
            clim=clim,
            normal=normal,
            origin=origin,
            origin_translation=False,
            show_scalar_bar=False,
            normal_rotation=False,
            cmap=cmap,
            user_matrix=b,
            render=False,
        )

        prop = self.clip.GetOutlineProperty()
        prop.SetOpacity(0)

        prop = self.clip.GetEdgesProperty()
        prop.SetOpacity(0)

        actor = self.base_view.plotter.show_grid(
            xtitle=labels[0],
            ytitle=labels[1],
            ztitle=labels[2],
            font_size=8,
            minor_ticks=True,
        )

        actor.SetAxisBaseForX(*T[:, 0])
        actor.SetAxisBaseForY(*T[:, 1])
        actor.SetAxisBaseForZ(*T[:, 2])

        actor.bounds = bounds.ravel()
        actor.SetXAxisRange(limits[0])
        actor.SetYAxisRange(limits[1])
        actor.SetZAxisRange(limits[2])

        axis0_args = *limits[0], actor.n_xlabels, actor.x_label_format
        axis1_args = *limits[1], actor.n_ylabels, actor.y_label_format
        axis2_args = *limits[2], actor.n_zlabels, actor.z_label_format

        axis0_label = pv.plotting.cube_axes_actor.make_axis_labels(*axis0_args)
        axis1_label = pv.plotting.cube_axes_actor.make_axis_labels(*axis1_args)
        axis2_label = pv.plotting.cube_axes_actor.make_axis_labels(*axis2_args)

        actor.SetAxisLabels(0, axis0_label)
        actor.SetAxisLabels(1, axis1_label)
        actor.SetAxisLabels(2, axis2_label)

        # self.base_view.reset_view()

        self.clip.AddObserver("InteractionEvent", self.interaction_callback)

        self.P_inv = np.linalg.inv(P)

    def interaction_callback(self, caller, event):
        orig = caller.GetOrigin()
        # norm = caller.GetNormal()

        # norm /= np.linalg.norm(norm)
        # norm = self.norm

        ind = np.array(self.norm).tolist().index(1)

        value = np.dot(self.P_inv, orig)[ind]

        self.view_model.set_number("slice_value", value)

        # self.slice_data()
