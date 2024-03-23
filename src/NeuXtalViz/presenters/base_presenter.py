class NeuXtalVizPresenter:

    def __init__(self, view, model):

        self.view = view
        self.model = model

        self.view.connect_manual_axis(self.view_manual)

        self.view.connect_reciprocal_axes(self.view_bc_star,
                                          self.view_ca_star,
                                          self.view_ab_star)

        self.view.connect_real_axes(self.view_bc,
                                    self.view_ca,
                                    self.view_ab)

        self.view.connect_save_screenshot(self.save_screenshot)
        self.view.connect_reciprocal_real_compass(self.change_lattice)

    def change_lattice(self):
        """
        Enable or disable reciprocal lattice.

        """

        T = self.model.get_transform(self.view.reciprocal_lattice())

        self.view.set_transform(T)

    def save_screenshot(self):
        """
        Save image.

        """        

        filename = self.view.save_screenshot_file_dialog()

        if filename:

            self.view.save_screenshot(filename)

    def view_manual(self):
        """
        Manual axis view.

        """

        indices = self.view.get_manual_axis_indices()

        if indices is not None:
            vec = self.model.get_vector(*indices)
            if vec is not None:
                self.view.view_vector(vec)

    def view_ab_star(self):
        """
        :math:`c`-axis view.

        """

        vecs = self.model.ab_star_axes()
        if vecs is not None:
            self.view.view_vector(vecs)

    def view_bc_star(self):
        """
        :math:`a`-axis view.

        """

        vecs = self.model.bc_star_axes()
        if vecs is not None:
            self.view.view_vector(vecs)

    def view_ca_star(self):
        """
        :math:`b`-axis view.

        """

        vecs = self.model.ca_star_axes()
        if vecs is not None:
            self.view.view_vector(vecs)

    def view_ab(self):
        """
        :math:`c^\ast`-axis view.

        """

        vecs = self.model.ab_axes()
        if vecs is not None:
            self.view.view_vector(vecs)

    def view_bc(self):
        """
        :math:`a^\ast`-axis view.

        """

        vecs = self.model.bc_axes()
        if vecs is not None:
            self.view.view_vector(vecs)

    def view_ca(self):
        """
        :math:`b^\ast`-axis view.

        """

        vecs = self.model.ca_axes()
        if vecs is not None:
            self.view.view_vector(vecs)