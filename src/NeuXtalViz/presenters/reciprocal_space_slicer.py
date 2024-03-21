from mantid.simpleapi import LoadMD

class ReciprocalSpaceSlicer:

    def __init__(self, view, model):

        self.view = view
        self.model = model

        self.view.manual_button.clicked.connect(self.view_manual)
        self.view.slice_button.clicked.connect(self.slice_manual)

        self.view.a_star_button.clicked.connect(self.view_bc_star)
        self.view.b_star_button.clicked.connect(self.view_ca_star)
        self.view.c_star_button.clicked.connect(self.view_ab_star)

        self.view.a_button.clicked.connect(self.view_bc)
        self.view.b_button.clicked.connect(self.view_ca)
        self.view.c_button.clicked.connect(self.view_ab)

        self.view.view_combo.currentIndexChanged.connect(self.update_labels)
        self.view.slice_combo.currentIndexChanged.connect(self.update_labels)

        self.view.replot_button.clicked.connect(self.plot_data)

        histo = LoadMD(Filename='/SNS/CORELLI/IPTS-15331/shared/normalization/Ba3Co2O6_50K_small/Ba3Co2O6_50K_small_6_m.nxs')
        #histo = LoadMD(Filename='/SNS/CORELLI/IPTS-15331/shared/normalization/Ba3Co2O6_50K_proj_small/Ba3Co2O6_50K_proj_small_6_m.nxs')
        self.model.set_md_histo_workspace('histo')
        
        self.histo = self.model.get_histo_info()
        
        self.normal = (0,0,-1)
        self.origin = 0

        self.plot_data()

    def plot_data(self):

        method = self.view.get_clim_clip_type()
        clim = self.model.calculate_clim(method)

        self.view.add_histo(self.histo, clim, self.normal, self.origin)
        self.view.set_transform(self.model.get_transform())

    def update_labels(self):

        self.view.update_axis_labels()
        self.view.update_normal_labels()

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

        indices = self.view.get_manual_axis_indices()

        if indices is not None:
            vec = self.model.get_vector(*indices)
            if vec is not None:
                self.view.view_vector(vec)

    def slice_manual(self):

        indices = self.view.get_manual_normal_indices()
        value = self.view.get_manual_normal_value()
        
        if indices is not None:
            vec = self.model.get_vector(*indices)
            if vec is not None:
                self.normal = -vec
                print(vec)
                self.origin = value if value is not None else 0
                self.plot_data()