from typing import Optional

from pydantic import BaseModel, Field

from NeuXtalViz.view_models.base_view_model import NeuXtalVizViewModel


class VolumeSlicerControls(BaseModel):
    vol_scale: str = Field(default="Linear")
    opacity: str = Field(default="Linear")
    opacity_range: str = Field("Low->High")
    clim_clip_type: str = Field(default="Q₃/Q₁±1.5×IQR")
    cbar: str = Field(default="Sequential")

    slice_plane: str = Field(default="Axis 1/2")
    slice_value: Optional[float] = Field(default=0.0, title="Slice")
    slice_thickness: Optional[float] = Field(default=0.1, title="Thickness")
    slice_scale: str = Field(default="linear")
    vlim_clip_type: str = Field(default="Q₃/Q₁±1.5×IQR")
    vmin: Optional[float] = Field(default=0.0, title="Min")
    vmax: Optional[float] = Field(default=0.0, title="Max")
    xmin: Optional[float] = Field(default=0.0, title="X Min")
    xmax: Optional[float] = Field(default=0.0, title="X Max")
    ymin: Optional[float] = Field(default=0.0, title="Y Min")
    ymax: Optional[float] = Field(default=0.0, title="Y Max")

    cut_line: str = Field(default="Axis 1")
    cut_value: Optional[float] = Field(default=0.0, title="Cut")
    cut_thickness: Optional[float] = Field(default=0.1, title="Thickness")
    cut_scale: str = Field(default="linear")


class VolumeSlicerViewModel(NeuXtalVizViewModel):
    def __init__(self, model, binding):
        super().__init__(model, binding)

        self.vs_controls = VolumeSlicerControls()
        self.draw_idle = True
        self.slice_idle = True
        self.cut_idle = True

        self.vs_controls_bind = binding.new_bind(self.vs_controls)
        self.slice_lim_bind = binding.new_bind()
        self.colorbar_lim_bind = binding.new_bind()
        self.cut_lim_bind = binding.new_bind()
        self.redraw_data_bind = binding.new_bind()
        self.slice_data_bind = binding.new_bind()
        self.cut_data_bind = binding.new_bind()
        self.add_histo_bind = binding.new_bind()
        self.add_slice_bind = binding.new_bind()
        self.add_cut_bind = binding.new_bind()

    def set_number(self, key, value):
        try:
            setattr(self.vs_controls, key, float(value))
        except Exception:
            setattr(self.vs_controls, key, None)

    def get_colormap(self):
        return self.vs_controls.cbar

    def get_vol_scale(self):
        return self.vs_controls.vol_scale.lower()

    def set_vol_scale(self, value):
        self.vs_controls.vol_scale = value
        self.redraw_data_bind.update_in_view(None)

    def get_opacity(self):
        return self.vs_controls.opacity

    def set_opacity(self, value):
        self.vs_controls.opacity = value
        self.redraw_data_bind.update_in_view(None)

    def get_opacity_range(self):
        return self.vs_controls.opacity_range

    def set_opacity_range(self, value):
        self.vs_controls.opacity_range = value
        self.redraw_data_bind.update_in_view(None)

    def set_clim_clip_type(self, value):
        self.vs_controls.clim_clip_type = value
        self.redraw_data_bind.update_in_view(None)

    def set_cbar(self, value):
        self.vs_controls.cbar = value
        self.redraw_data_bind.update_in_view(None)

    def set_slice_plane(self, value):
        self.vs_controls.slice_plane = value
        self.redraw_data_bind.update_in_view(None)

    def set_vlim_clip_type(self, value):
        self.vs_controls.vlim_clip_type = value
        self.update_slice()

    def set_slice_scale(self, value):
        self.vs_controls.slice_scale = value.lower()
        self.update_slice()

    def set_cut_line(self, value):
        self.vs_controls.cut_line = value
        self.update_cut()

    def set_cut_scale(self, value):
        self.vs_controls.cut_scale = value.lower()
        self.update_cut()

    def update_limits(self):
        xmin = self.vs_controls.xmin
        xmax = self.vs_controls.xmax
        ymin = self.vs_controls.ymin
        ymax = self.vs_controls.ymax

        if (
            xmin is not None
            and xmax is not None
            and ymin is not None
            and ymax is not None
        ):
            if xmin < xmax and ymin < ymax:
                xlim = [xmin, xmax]
                ylim = [ymin, ymax]
                self.slice_lim_bind.update_in_view((xlim, ylim))
                line_cut = self.vs_controls.cut_line
                lim = xlim if line_cut == "Axis 1" else ylim
                self.cut_lim_bind.update_in_view(lim)

    def update_cvals(self):
        vmin = self.vs_controls.vmin
        vmax = self.vs_controls.vmax
        if vmin is not None and vmax is not None:
            if vmin < vmax:
                if vmin <= 0 and self.vs_controls.slice_scale == "log":
                    vmin = vmax / 10
                self.colorbar_lim_bind.update_in_view((vmin, vmax))

    def update_slice(self):
        if self.model.is_histo_loaded():
            self.slice_data_bind.update_in_view(None)

    def update_cut(self):
        if self.model.is_histo_loaded():
            self.cut_data_bind.update_in_view(None)

    def load_NXS_complete(self, result=None):
        self.update_oriented_lattice()

    def load_NXS_process(self, filename, progress):
        progress("Processing...", 1)

        progress("Loading NeXus file...", 10)

        self.model.load_md_histo_workspace(filename)

        progress("Loading NeXus file...", 50)

        progress("Loading NeXus file...", 80)

        progress("NeXus file loaded!", 100)

    def get_normal(self):
        slice_plane = self.vs_controls.slice_plane

        if slice_plane == "Axis 1/2":
            norm = [0, 0, 1]
        elif slice_plane == "Axis 1/3":
            norm = [0, 1, 0]
        else:
            norm = [1, 0, 0]

        return norm

    def get_axis(self):
        axis = [1 if not norm else 0 for norm in self.get_normal()]
        ind = [i for i, ax in enumerate(axis) if ax == 1]

        line_cut = self.vs_controls.cut_line

        if line_cut == "Axis 1":
            axis[ind[0]] = 0
        else:
            axis[ind[1]] = 0

        return axis

    def get_clim_method(self):
        ctype = self.vs_controls.clim_clip_type

        if ctype == "μ±3×σ":
            method = "normal"
        elif ctype == "Q₃/Q₁±1.5×IQR":
            method = "boxplot"
        else:
            method = None

        return method

    def get_vlim_method(self):
        ctype = self.vs_controls.vlim_clip_type

        if ctype == "μ±3×σ":
            method = "normal"
        elif ctype == "Q₃/Q₁±1.5×IQR":
            method = "boxplot"
        else:
            method = None

        return method

    def redraw_data_complete(self, result):
        if result is not None:
            histo, normal, norm, value, trans = result

            self.add_histo_bind.update_in_view((histo, normal, norm, value))
            self.show_axes_bind.update_in_view(
                (trans, self.controls.reciprocal_lattice, self.controls.show_axes)
            )
            self.parallel_projection_bind.update_in_view(
                self.controls.parallel_projection
            )

        self.draw_idle = True

    def redraw_data_process(self, progress):
        if self.draw_idle and self.model.is_histo_loaded():
            self.draw_idle = False

            progress("Processing...", 1)

            progress("Updating volume...", 20)

            norm = self.get_normal()

            histo = self.model.get_histo_info(norm)

            data = histo["signal"]

            data = self.model.calculate_clim(data, self.get_clim_method())

            progress("Updating volume...", 50)

            histo["signal"] = data

            value = self.vs_controls.slice_value

            normal = -self.model.get_normal_plane(norm)

            # origin = self.model.get_normal('[hkl]', orig)

            if value is not None:
                progress("Volume drawn!", 100)

                return histo, normal, norm, value, self.model.get_transform()

            else:
                progress("Invalid parameters.", 0)

    def slice_data_complete(self, result):
        if result is not None:
            self.add_slice_bind.update_in_view(result)

        self.slice_idle = True

    def slice_data_process(self, progress):
        if self.slice_idle and self.model.is_histo_loaded():
            self.slice_idle = False

            norm = self.get_normal()

            thick = self.vs_controls.slice_thickness
            value = self.vs_controls.slice_value

            if thick is not None:
                progress("Processing...", 1)

                progress("Updating slice...", 50)

                slice_histo = self.model.get_slice_info(norm, value, thick)

                data = slice_histo["signal"]

                data = self.model.calculate_clim(data, self.get_vlim_method())

                slice_histo["signal"] = data

                progress("Data sliced!", 100)

                return slice_histo

    def cut_data_complete(self, result):
        if result is not None:
            self.add_cut_bind.update_in_view(result)

        self.cut_idle = True

    def cut_data_process(self, progress):
        if self.cut_idle and self.model.is_sliced():
            self.cut_idle = False

            value = self.vs_controls.cut_value
            thick = self.vs_controls.cut_thickness

            axis = self.get_axis()

            if value is not None and thick is not None:
                progress("Processing...", 1)

                progress("Updating cut...", 50)

                progress("Data cut!", 100)

                cut_histo = self.model.get_cut_info(axis, value, thick)

                return cut_histo

    def save_slice(self, filename):
        if self.model.is_sliced():
            self.model.save_slice(filename)

    def save_cut(self, filename):
        if self.model.is_cut():
            self.model.save_cut(filename)
