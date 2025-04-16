from asyncio import create_task, ensure_future, sleep
from functools import partial
from io import BytesIO
from tempfile import NamedTemporaryFile
from threading import Thread

import numpy as np
import pyvista as pv
from matplotlib.figure import Figure
from matplotlib.transforms import Affine2D
from nova.trame.view.components import InputField
from nova.trame.view.components.visualization import MatplotlibFigure
from nova.trame.view.layouts import GridLayout, HBoxLayout, VBoxLayout
from trame.app.file_upload import ClientFile
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
        self.view_model.slice_lim_bind.connect(self.set_slice_lim)
        self.view_model.cut_lim_bind.connect(self.set_cut_lim)
        self.view_model.colorbar_lim_bind.connect(self.update_colorbar_vlims)
        self.view_model.vs_controls_bind.connect("vs_controls")
        self.view_model.redraw_data_bind.connect(self.redraw_data)
        self.view_model.slice_data_bind.connect(self.slice_data)
        self.view_model.cut_data_bind.connect(self.cut_data)
        self.view_model.add_histo_bind.connect(self.trigger_add_histo)
        self.view_model.add_slice_bind.connect(self.add_slice)
        self.view_model.add_cut_bind.connect(self.add_cut)

        self.histo = None
        ensure_future(self.add_histo_loop())

        self.create_ui()

    @property
    def state(self):
        return self.server.state

    def create_ui(self):
        with GridLayout(classes="bg-white pa-2", columns=2):
            self.base_view = VisualizationPanel(self.server, self.view_model)

            with VBoxLayout():
                with HBoxLayout(valign="center"):
                    InputField(v_model="vs_controls.vol_scale", type="select")
                    InputField(v_model="vs_controls.opacity", type="select")
                    InputField(v_model="vs_controls.opacity_range", type="select")
                    InputField(v_model="vs_controls.clim_clip_type", type="select")
                    InputField(v_model="vs_controls.cbar", type="select")
                    with html.Div():
                        InputField(
                            ref="fread",
                            classes="d-none",
                            type="file",
                            __events=["change"],
                            change=(self.load_NXS, "[$event.target.files]"),
                        )
                        vuetify.VBtn("Load NXS", click="trame.refs.fread.click()")

                with HBoxLayout():
                    self.fig_slice = Figure(layout="constrained")
                    self.ax_slice = self.fig_slice.subplots(1, 1)
                    self.cb = None

                    self.slice_view = MatplotlibFigure(self.fig_slice, webagg=True)
                    InputField(
                        v_model="vs_controls.vmin",
                        direction="vertical",
                        max=("vs_controls.vlims[1]",),
                        min=("vs_controls.vlims[0]",),
                        type="slider",
                    )
                    InputField(
                        v_model="vs_controls.vmax",
                        direction="vertical",
                        max=("vs_controls.vlims[1]",),
                        min=("vs_controls.vlims[0]",),
                        type="slider",
                    )

                with HBoxLayout(valign="center"):
                    InputField(v_model="vs_controls.slice_plane", type="select")
                    InputField(v_model="vs_controls.slice_value")
                    InputField(v_model="vs_controls.slice_thickness")
                    vuetify.VBtn("Save Slice", click=self.save_slice)
                    InputField(v_model="vs_controls.slice_scale", type="select")

                with HBoxLayout():
                    InputField(v_model="vs_controls.xmin")
                    InputField(v_model="vs_controls.xmax")
                    InputField(v_model="vs_controls.ymin")
                    InputField(v_model="vs_controls.ymax")
                    InputField(v_model="vs_controls.vlim_clip_type", type="select")
                    InputField(v_model="vs_controls.vmin")
                    InputField(v_model="vs_controls.vmax")

                with HBoxLayout():
                    self.fig_cut = Figure(layout="constrained", figsize=(6.4, 3.2))
                    self.ax_cut = self.fig_cut.subplots(1, 1)

                    self.cut_view = MatplotlibFigure(self.fig_cut, webagg=True)

                with HBoxLayout(valign="center"):
                    InputField(v_model="vs_controls.cut_line", type="select")
                    InputField(v_model="vs_controls.cut_value")
                    InputField(v_model="vs_controls.cut_thickness")
                    vuetify.VBtn("Save Cut", click=self.save_cut)
                    InputField(v_model="vs_controls.cut_scale", type="select")

    def load_in_background(self, filename):
        self.view_model.load_NXS_process(filename)
        self.view_model.load_NXS_complete()
        self.loading_data = False

    async def monitor_load(self):
        while self.loading_data:
            await sleep(0.1)

        self.view_model.update_processing("Loading NeXus file...", 50)
        self.view_model.update_processing("Loading NeXus file...", 80)
        self.view_model.update_processing("NeXus file loaded!", 100)
        self.redraw_data()

    def load_NXS(self, files):
        if len(files) < 1:
            return

        nxs_file = ClientFile(files[0])
        with NamedTemporaryFile(suffix=".nxs", delete=False) as _file:
            _file.write(nxs_file.content)
            nxs_filename = _file.name

        self.view_model.update_processing("Processing...", 1)
        self.view_model.update_processing("Loading NeXus file...", 10)

        self.loading_data = True
        thread = Thread(
            target=partial(self.load_in_background, nxs_filename), daemon=True
        )
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

        self.slice_data()

    def redraw_data(self, _=None):
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

        self.clip.AddObserver("InteractionEvent", self.interaction_callback)

        self.P_inv = np.linalg.inv(P)

        self.base_view.reset_view()

    async def add_histo_loop(self):
        while True:
            if self.histo is not None:
                self.add_histo(self.histo)
                self.histo = None

            await sleep(0.1)

    def trigger_add_histo(self, result):
        self.histo = result

    def interaction_callback(self, caller, event):
        orig = caller.GetOrigin()

        ind = np.array(self.norm).tolist().index(1)

        value = np.dot(self.P_inv, orig)[ind]

        self.view_model.set_number("slice_value", value)

        self.slice_data()

    def slice_in_background(self):
        self.slicing_result = self.view_model.slice_data_process()
        self.view_model.slice_data_complete(self.slicing_result)
        self.slicing = False

    async def monitor_slice(self):
        while self.slicing:
            await sleep(0.1)

        self.view_model.update_processing("Data sliced!", 100)

        self.cut_data()

    def slice_data(self, _=None):
        self.view_model.update_processing("Processing...", 1)
        self.view_model.update_processing("Updating slice...", 50)

        self.slicing = True
        thread = Thread(target=self.slice_in_background, daemon=True)
        thread.start()

        create_task(self.monitor_slice())

    def __format_axis_coord(self, x, y):
        x, y, _ = np.dot(self.T_inv, [x, y, 1])
        return "x={:.3f}, y={:.3f}".format(x, y)

    def add_slice(self, slice_dict):
        cmap = cmaps[self.view_model.get_cbar()]

        x = slice_dict["x"]
        y = slice_dict["y"]

        labels = slice_dict["labels"]
        title = slice_dict["title"]
        signal = slice_dict["signal"]

        scale = self.view_model.get_slice_scale()

        vmin = np.nanmin(signal)
        vmax = np.nanmax(signal)

        self.view_model.set_vlims(vmin, vmax)

        if np.isclose(vmax, vmin) or not np.isfinite([vmin, vmax]).all():
            vmin, vmax = (0.1, 1) if scale == "log" else (0, 1)

        T = slice_dict["transform"]
        aspect = slice_dict["aspect"]

        self.T_inv = np.linalg.inv(T)

        self.ax_slice.format_coord = self.__format_axis_coord

        transform = Affine2D(T) + self.ax_slice.transData
        self.transform = transform

        xlim = np.array([x.min(), x.max()])
        ylim = np.array([y.min(), y.max()])

        if self.cb is not None:
            self.cb.remove()

        self.ax_slice.clear()

        im = self.ax_slice.pcolormesh(
            x,
            y,
            signal,
            norm=scale,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            shading="flat",
            rasterized=True,
            transform=transform,
        )

        self.im = im
        vmin, self.vmax = self.im.norm.vmin, self.im.norm.vmax

        self.view_model.set_number("vmin", self.im.norm.vmin)
        self.view_model.set_number("vmax", self.im.norm.vmax)
        self.view_model.set_number("xmin", xlim[0])
        self.view_model.set_number("xmax", xlim[1])
        self.view_model.set_number("ymin", ylim[0])
        self.view_model.set_number("ymax", ylim[1])

        self.ax_slice.set_aspect(aspect)
        self.ax_slice.set_xlabel(labels[0])
        self.ax_slice.set_ylabel(labels[1])
        self.ax_slice.set_title(title)
        self.ax_slice.minorticks_on()

        self.ax_slice.xaxis.get_major_locator().set_params(integer=True)
        self.ax_slice.yaxis.get_major_locator().set_params(integer=True)

        self.cb = self.fig_slice.colorbar(self.im, ax=self.ax_slice)
        self.cb.minorticks_on()

        self.slice_view.update(self.fig_slice)

    def cut_in_background(self):
        self.cut_result = self.view_model.cut_data_process()
        self.view_model.cut_data_complete(self.cut_result)
        self.cutting = False

    async def monitor_cut(self):
        while self.cutting:
            await sleep(0.1)

        self.view_model.update_processing("Data cut!", 100)

    def cut_data(self, _=None):
        self.view_model.update_processing("Processing...", 1)
        self.view_model.update_processing("Updating cut...", 50)

        self.cutting = True
        thread = Thread(target=self.cut_in_background, daemon=True)
        thread.start()

        create_task(self.monitor_cut())

    def add_cut(self, cut_dict):
        x = cut_dict["x"]
        y = cut_dict["y"]
        e = cut_dict["e"]

        val = cut_dict["value"]

        label = cut_dict["label"]
        title = cut_dict["title"]

        scale = self.view_model.get_cut_scale()

        line_cut = self.view_model.get_cut_value()

        lines = self.ax_slice.get_lines()
        for line in lines:
            line.remove()

        xlim = self.view_model.get_xlim()
        ylim = self.view_model.get_ylim()

        thick = self.view_model.get_cut_thickness()

        delta = 0 if thick is None else thick / 2

        if line_cut == "Axis 2":
            l0 = [val - delta, val - delta], ylim
            l1 = [val + delta, val + delta], ylim
        else:
            l0 = xlim, [val - delta, val - delta]
            l1 = xlim, [val + delta, val + delta]

        self.ax_slice.plot(*l0, "w--", linewidth=1, transform=self.transform)
        self.ax_slice.plot(*l1, "w--", linewidth=1, transform=self.transform)

        self.ax_cut.clear()

        self.ax_cut.errorbar(x, y, e)
        self.ax_cut.set_xlabel(label)
        self.ax_cut.set_yscale(scale)
        self.ax_cut.set_title(title)
        self.ax_cut.minorticks_on()

        self.ax_cut.xaxis.get_major_locator().set_params(integer=True)

        self.cut_view.update(self.fig_cut)

    def set_slice_lim(self, lims):
        xlim, ylim = lims

        if self.cb is not None:
            xmin, xmax = xlim
            ymin, ymax = ylim
            T = np.linalg.inv(self.T_inv)
            xmin, ymin, _ = np.dot(T, [xmin, ymin, 1])
            xmax, ymax, _ = np.dot(T, [xmax, ymax, 1])
            self.ax_slice.set_xlim(xmin, xmax)
            self.ax_slice.set_ylim(ymin, ymax)
            self.slice_view.update(self.fig_slice)

    def update_colorbar_vlims(self, vlims):
        vmin, vmax = vlims

        if self.cb is not None:
            self.im.set_clim(vmin=vmin, vmax=vmax)
            self.cb.update_normal(self.im)
            self.cb.minorticks_on()

            self.slice_view.update(self.fig_slice)

    def set_cut_lim(self, lim):
        if self.cb is not None:
            self.ax_cut.set_xlim(*lim)
            self.cut_view.update(self.fig_cut)

    def save_slice(self):
        data = BytesIO()
        self.fig_slice.savefig(data, format="png")
        data.seek(0)
        self.base_view.js_download(("slice.png", data.read()))

    def save_cut(self):
        data = BytesIO()
        self.fig_cut.savefig(data, format="png")
        data.seek(0)
        self.base_view.js_download(("cut.png", data.read()))
