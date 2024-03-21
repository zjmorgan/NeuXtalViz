class Sample:

    def __init__(self, view, model):

        self.view = view
        self.model = model

        #self.view.sample_combo.currentIndexChanged.connect(self.add_shape)

        self.view.manual_button.clicked.connect(self.view_manual)

        self.view.a_star_button.clicked.connect(self.view_bc_star)
        self.view.b_star_button.clicked.connect(self.view_ca_star)
        self.view.c_star_button.clicked.connect(self.view_ab_star)

        self.view.a_button.clicked.connect(self.view_bc)
        self.view.b_button.clicked.connect(self.view_ca)
        self.view.c_button.clicked.connect(self.view_ab)

        self.view.view_combo.currentIndexChanged.connect(self.update_labels)
        self.view.proj_box.clicked.connect(self.change_proj)
        self.view.reset_button.clicked.connect(self.reset_view)

    def add_shape(self, shape_dict):

        self.plotter.clear_actors()

        self.change_proj()

    def update_labels(self):

        self.view.update_axis_labels()

    def change_proj(self):

        self.view.change_proj()

    def reset_view(self):

        self.view.reset_view()

    def view_ab_star(self):

        vecs = self.model.ab_star_axes()
        if vecs is not None:
            self.view.view_vector(vecs)

    def view_bc_star(self):

        vecs = self.model.bc_star_axes()
        if vecs is not None:
            self.view.view_vector(vecs)

    def view_ca_star(self):
        vecs = self.model.ca_star_axes()
        if vecs is not None:
            self.view.view_vector(vecs)

    def view_ab(self):

        vecs = self.model.ab_axes()
        if vecs is not None:
            self.view.view_vector(vecs)

    def view_bc(self):

        vecs = self.model.bc_axes()
        if vecs is not None:
            self.view.view_vector(vecs)

    def view_ca(self):

        vecs = self.model.ca_axes()
        if vecs is not None:
            self.view.view_vector(vecs)

    def view_manual(self):

        indices = self.view.get_manual_indices()

        if indices is not None:
            vec = self.model.get_vector(*indices)
            if vec is not None:
                self.view.view_vector(vec)