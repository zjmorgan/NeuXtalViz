from enum import Enum
from typing import Optional

from pydantic import BaseModel, Field

from NeuXtalViz.view_models.base_view_model import NeuXtalVizViewModel


class AxisOptions(str, Enum):
    linear = "Linear"
    log = "Log"


class OpacityOptions(str, Enum):
    linear = "Linear"
    geometric = "Geometric"
    sigmoid = "Sigmoid"


class OpacityRangeOptions(str, Enum):
    low_to_high = "Low->High"
    high_to_low = "High->Low"


class ClipTypeOptions(str, Enum):
    minmax = "Min/Max"
    normal = "μ±3×σ"
    boxplot = "Q₃/Q₁±1.5×IQR"


class ColorbarOptions(str, Enum):
    sequential = "Sequential"
    rainbow = "Rainbow"
    binary = "Binary"
    diverging = "Diverging"
    modified = "Modified"


class SlicePlaneOptions(str, Enum):
    one_half = "Axis 1/2"
    one_third = "Axis 1/3"
    two_thirds = "Axis 2/3"


class CutLineOptions(str, Enum):
    axis_one = "Axis 1"
    axis_two = "Axis 2"


class VolumeSlicerControls(BaseModel):
    vol_scale: AxisOptions = Field(default=AxisOptions.linear, title="Scale")
    opacity: OpacityOptions = Field(default=OpacityOptions.linear, title="Opacity")
    opacity_range: OpacityRangeOptions = Field(
        OpacityRangeOptions.low_to_high, title="Opacity Range"
    )
    clim_clip_type: ClipTypeOptions = Field(
        default=ClipTypeOptions.boxplot, title="Clip Type"
    )
    cbar: ColorbarOptions = Field(
        default=ColorbarOptions.sequential, title="Color Scale"
    )

    slice_plane: SlicePlaneOptions = Field(
        default=SlicePlaneOptions.one_half, title="Plane"
    )
    slice_value: Optional[float] = Field(default=0.0, title="Slice")
    slice_thickness: Optional[float] = Field(default=0.1, title="Thickness")
    slice_scale: AxisOptions = Field(default=AxisOptions.linear, tilte="Scale")
    vlim_clip_type: ClipTypeOptions = Field(
        default=ClipTypeOptions.boxplot, title="Clip Type"
    )
    vlims: list[float] = Field(default=[0.0, 0.0])
    vmin: Optional[float] = Field(default=0.0, title="Color Min")
    vmax: Optional[float] = Field(default=0.0, title="Color Max")
    xmin: Optional[float] = Field(default=0.0, title="X Min")
    xmax: Optional[float] = Field(default=0.0, title="X Max")
    ymin: Optional[float] = Field(default=0.0, title="Y Min")
    ymax: Optional[float] = Field(default=0.0, title="Y Max")

    cut_line: CutLineOptions = Field(default=CutLineOptions.axis_one, title="Line")
    cut_value: Optional[float] = Field(default=0.0, title="Cut")
    cut_thickness: Optional[float] = Field(default=0.1, title="Thickness")
    cut_scale: AxisOptions = Field(default=AxisOptions.linear, title="Scale")


class VolumeSlicerViewModel(NeuXtalVizViewModel):
    def __init__(self, model, binding):
        super().__init__(model, binding)

        self.vs_controls = VolumeSlicerControls()
        self.draw_idle = True
        self.slice_idle = True
        self.cut_idle = True

        self.vs_controls_bind = binding.new_bind(
            self.vs_controls, callback_after_update=self.process_vs_updates
        )
        self.slice_lim_bind = binding.new_bind()
        self.colorbar_lim_bind = binding.new_bind()
        self.cut_lim_bind = binding.new_bind()
        self.redraw_data_bind = binding.new_bind()
        self.slice_data_bind = binding.new_bind()
        self.cut_data_bind = binding.new_bind()
        self.add_histo_bind = binding.new_bind()
        self.add_slice_bind = binding.new_bind()
        self.add_cut_bind = binding.new_bind()

    def process_vs_updates(self, results):
        for update in results.get("updated", []):
            match update:
                case (
                    "vol_scale"
                    | "opacity"
                    | "opacity_range"
                    | "clim_clip_type"
                    | "cbar"
                    | "slice_plane"
                ):
                    pass
                    # TODO: self.redraw_data_bind.update_in_view(None)
                case (
                    "slice_value" | "slice_thickness" | "slice_scale" | "vlim_clip_type"
                ):
                    self.update_slice()
                case "xmin" | "xmax" | "ymin" | "ymax":
                    self.update_limits()
                case "vmin" | "vmax":
                    if self.vs_controls.vmin > self.vs_controls.vmax:
                        self.vs_controls.vmin = self.vs_controls.vmax
                        self.vs_controls_bind.update_in_view(self.vs_controls)
                    self.update_cvals()
                case "cut_line" | "cut_value" | "cut_thickness" | "cut_scale":
                    self.update_cut()

    def set_number(self, key, value):
        try:
            setattr(self.vs_controls, key, float(value))
        except Exception:
            setattr(self.vs_controls, key, None)
        self.vs_controls_bind.update_in_view(self.vs_controls)

    def get_colormap(self):
        return self.vs_controls.cbar.value

    def get_vol_scale(self):
        return self.vs_controls.vol_scale.value.lower()

    def set_vol_scale(self, value):
        self.vs_controls.vol_scale = AxisOptions(value)
        self.redraw_data_bind.update_in_view(None)

    def get_opacity(self):
        return self.vs_controls.opacity.value

    def set_opacity(self, value):
        self.vs_controls.opacity = OpacityOptions(value)
        self.redraw_data_bind.update_in_view(None)

    def get_opacity_range(self):
        return self.vs_controls.opacity_range.value

    def set_opacity_range(self, value):
        self.vs_controls.opacity_range = OpacityRangeOptions(value)
        self.redraw_data_bind.update_in_view(None)

    def set_clim_clip_type(self, value):
        self.vs_controls.clim_clip_type = ClipTypeOptions(value)
        self.redraw_data_bind.update_in_view(None)

    def set_cbar(self, value):
        self.vs_controls.cbar = value
        self.redraw_data_bind.update_in_view(None)

    def set_slice_plane(self, value):
        self.vs_controls.slice_plane = SlicePlaneOptions(value)
        self.redraw_data_bind.update_in_view(None)

    def set_vlims(self, vmin, vmax):
        self.vs_controls.vlims[0] = vmin
        self.vs_controls.vlims[1] = vmax

    def set_vlim_clip_type(self, value):
        self.vs_controls.vlim_clip_type = value
        self.update_slice()

    def set_slice_scale(self, value):
        self.vs_controls.slice_scale = AxisOptions(value)
        self.update_slice()

    def set_cut_line(self, value):
        self.vs_controls.cut_line = CutLineOptions(value)
        self.update_cut()

    def get_cut_scale(self):
        return self.vs_controls.cut_scale.value.lower()

    def set_cut_scale(self, value):
        self.vs_controls.cut_scale = AxisOptions(value)
        self.update_cut()

    def get_cut_value(self):
        return self.vs_controls.cut_value

    def get_cut_thickness(self):
        return self.vs_controls.cut_thickness

    def get_xlim(self):
        return [self.vs_controls.xmin, self.vs_controls.xmax]

    def get_ylim(self):
        return [self.vs_controls.ymin, self.vs_controls.ymax]

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
                lim = xlim if line_cut == CutLineOptions.axis_one else ylim
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

    def load_NXS_process(self, filename, progress=lambda status, progress: None):
        progress("Processing...", 1)

        progress("Loading NeXus file...", 10)

        self.model.load_md_histo_workspace(filename)

        progress("Loading NeXus file...", 50)

        progress("Loading NeXus file...", 80)

        progress("NeXus file loaded!", 100)

    def get_normal(self):
        slice_plane = self.vs_controls.slice_plane

        if slice_plane == SlicePlaneOptions.one_half:
            norm = [0, 0, 1]
        elif slice_plane == SlicePlaneOptions.one_third:
            norm = [0, 1, 0]
        else:
            norm = [1, 0, 0]

        return norm

    def get_axis(self):
        axis = [1 if not norm else 0 for norm in self.get_normal()]
        ind = [i for i, ax in enumerate(axis) if ax == 1]

        line_cut = self.vs_controls.cut_line

        if line_cut == CutLineOptions.axis_one:
            axis[ind[0]] = 0
        else:
            axis[ind[1]] = 0

        return axis

    def get_slice_scale(self):
        return self.vs_controls.slice_scale.value.lower()

    def get_cbar(self):
        return self.vs_controls.cbar

    def get_clim_method(self):
        ctype = self.vs_controls.clim_clip_type.value

        if ctype == ClipTypeOptions.normal:
            method = "normal"
        elif ctype == ClipTypeOptions.boxplot:
            method = "boxplot"
        else:
            method = None

        return method

    def get_vlim_method(self):
        ctype = self.vs_controls.vlim_clip_type

        if ctype == ClipTypeOptions.normal:
            method = "normal"
        elif ctype == ClipTypeOptions.boxplot:
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

    def redraw_data_process(self, progress=lambda status, progress: None):
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

    def slice_data_process(self, progress=lambda status, progress: None):
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

    def cut_data_process(self, progress=lambda status, progress: None):
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
