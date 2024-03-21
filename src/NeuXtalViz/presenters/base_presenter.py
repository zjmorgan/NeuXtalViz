class NeuXVizWidgetViewer:

    def __init__(self, view, model):

        self.view = view
        self.model = model

        self.view.manual_button.clicked.connect(self.view_manual)

        self.view.px_button.clicked.connect(self.view.view_yz)
        self.view.py_button.clicked.connect(self.view.view_zx)
        self.view.pz_button.clicked.connect(self.view.view_xy)

        self.view.mx_button.clicked.connect(self.view.view_zy)
        self.view.my_button.clicked.connect(self.view.view_xz)
        self.view.mz_button.clicked.connect(self.view.view_yx)

        self.view.a_star_button.clicked.connect(self.view_bc_star)
        self.view.b_star_button.clicked.connect(self.view_ca_star)
        self.view.c_star_button.clicked.connect(self.view_ab_star)

        self.view.a_button.clicked.connect(self.view_bc)
        self.view.b_button.clicked.connect(self.view_ca)
        self.view.c_button.clicked.connect(self.view_ab)

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

    def view_manual(self):
        """
        Manual axis view.

        """

        indices = self.view.get_manual_indices()

        if indices is not None:
            vec = self.model.get_vector(*indices)
            if vec is not None:
                self.view.view_vector(vec)