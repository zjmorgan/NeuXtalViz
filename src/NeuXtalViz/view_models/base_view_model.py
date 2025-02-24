from pydantic import BaseModel, Field


class Controls(BaseModel):
    camera_tab: int = Field(default=0)
    oriented_lattice_tab: int = Field(default=0)

    reciprocal_lattice: bool = Field(default=True, title="Reciprocal Lattice")
    show_axes: bool = Field(default=True, title="Show Axes")
    parallel_projection: bool = Field(default=True, title="Parallel Projection")

    manual_axis_type: str = Field(default="[hkl]")
    manual_axes: list[float] = Field(
        default=[0.0, 0.0, 0.0], title="Manual Axis Indices"
    )
    manual_up_axis_type: str = Field(default="[hkl]")
    manual_up_axes: list[float] = Field(
        default=[0.0, 0.0, 0.0], title="Manual Up Axis Indices"
    )


class LatticeParameters(BaseModel):
    a: float = Field(default=0.0)
    b: float = Field(default=0.0)
    c: float = Field(default=0.0)
    alpha: float = Field(default=0.0)
    beta: float = Field(default=0.0)
    gamma: float = Field(default=0.0)
    u: list[float] = Field(default=[0.0, 0.0, 0.0])
    v: list[float] = Field(default=[0.0, 0.0, 0.0])


class NeuXtalVizViewModel:
    def __init__(self, model, binding):
        self.model = model

        self.controls = Controls()
        self.lattice_parameters = LatticeParameters()

        self.controls_bind = binding.new_bind(
            self.controls, callback_after_update=self.process_updates
        )
        self.lattice_parameters_bind = binding.new_bind(self.lattice_parameters)
        self.show_axes_bind = binding.new_bind()
        self.parallel_projection_bind = binding.new_bind()

        self.progress_bind = binding.new_bind()
        self.status_bind = binding.new_bind()
        self.up_vector_bind = binding.new_bind()
        self.update_labels_bind = binding.new_bind()
        self.vector_bind = binding.new_bind()

    def process_updates(self, results):
        for update in results.get("updated", []):
            match update:
                case "reciprocal_lattice" | "show_axes":
                    self.update_axes()
                case "parallel_projection":
                    self.update_projection()

    def update_status(self, status):
        """
        Update status information.

        Parameters
        ----------
        status : str
            Information.

        """

        self.status_bind.update_in_view(status)

    def update_progress(self, progress):
        """
        Update progress step.

        Parameters
        ----------
        progress : int
            Step.

        """

        self.progress_bind.update_in_view(progress)

    def update_invalid(self):
        """
        Indicate invalid.

        """

        self.update_status("Invalid parameters.")
        self.update_progress(0)

    def update_complete(self, status="Complete!"):
        """
        Indicate complete.

        Parameters
        ----------
        status : str
            Information.

        """

        self.update_status(status)
        self.update_progress(0)

    def update_processing(self, status="Processing...", progress=1):
        """
        Indicate processing.

        Parameters
        ----------
        status : str
            Information.
        progress : int
            Step.

        """

        self.update_status(status)
        self.update_progress(progress)

    def update_oriented_lattice(self):
        """
        Update oriented lattice parameter display.

        """

        ol = self.model.get_oriented_lattice_parameters()
        if ol is not None:
            a, b, c, alpha, beta, gamma, u, v = ol
            self.lattice_parameters.a = a
            self.lattice_parameters.b = b
            self.lattice_parameters.c = c
            self.lattice_parameters.alpha = alpha
            self.lattice_parameters.beta = beta
            self.lattice_parameters.gamma = gamma
            self.lattice_parameters.u = u
            self.lattice_parameters.v = v

            self.lattice_parameters_bind.update_in_view(self.lattice_parameters)

    def update_axis_type(self, value):
        """
        Update the manual axis type.

        """

        self.controls.manual_axis_type = value
        self.update_labels_bind.update_in_view(None)

    def update_manual_axis(self, index, value):
        """
        Update the manual axis indices.

        """

        self.controls.manual_axes[index] = float(value)

    def update_up_axis_type(self, value):
        """
        Update the manual up axis type.

        """

        self.controls.manual_up_axis_type = value
        self.update_labels_bind.update_in_view(None)

    def update_manual_up_axis(self, index, value):
        """
        Update the manual up axis indices.

        """

        self.controls.manual_up_axes[index] = float(value)

    def update_axes(self):
        T = self.model.get_transform(self.controls.reciprocal_lattice)
        self.show_axes_bind.update_in_view(
            (T, self.controls.reciprocal_lattice, self.controls.show_axes)
        )

    def update_projection(self):
        self.parallel_projection_bind.update_in_view(self.controls.parallel_projection)
        self.update_axes()

    def change_lattice(self):
        """
        Enable or disable reciprocal lattice.

        """

        self.controls.reciprocal_lattice = not self.controls.reciprocal_lattice
        self.update_axes()

    def change_axes(self):
        self.controls.show_axes = not self.controls.show_axes
        self.update_axes()

    def change_projection(self):
        self.controls.parallel_projection = not self.controls.parallel_projection
        self.update_projection()

    def view_manual(self):
        """
        Manual axis view.

        """

        if self.controls.manual_axes:
            vec = self.model.get_vector(
                self.controls.manual_axis_type, self.controls.manual_axes
            )
            if vec is not None:
                self.vector_bind.update_in_view(vec)

    def view_up_manual(self):
        """
        Manual axis up view.

        """

        if self.controls.manual_up_axes:
            vec = self.model.get_vector(
                self.controls.manual_up_axis_type, self.controls.manual_up_axes
            )
            if vec is not None:
                self.up_vector_bind.update_in_view(vec)

    def view_ab_star(self):
        """
        :math:`c`-axis view.

        """

        vecs = self.model.ab_star_axes()
        if vecs is not None:
            self.vector_bind.update_in_view(vecs)

    def view_bc_star(self):
        """
        :math:`a`-axis view.

        """

        vecs = self.model.bc_star_axes()
        if vecs is not None:
            self.vector_bind.update_in_view(vecs)

    def view_ca_star(self):
        """
        :math:`b`-axis view.

        """

        vecs = self.model.ca_star_axes()
        if vecs is not None:
            self.vector_bind.update_in_view(vecs)

    def view_ab(self):
        """
        :math:`c^\ast`-axis view.

        """

        vecs = self.model.ab_axes()
        if vecs is not None:
            self.vector_bind.update_in_view(vecs)

    def view_bc(self):
        """
        :math:`a^\ast`-axis view.

        """

        vecs = self.model.bc_axes()
        if vecs is not None:
            self.vector_bind.update_in_view(vecs)

    def view_ca(self):
        """
        :math:`b^\ast`-axis view.

        """

        vecs = self.model.ca_axes()
        if vecs is not None:
            self.vector_bind.update_in_view(vecs)
