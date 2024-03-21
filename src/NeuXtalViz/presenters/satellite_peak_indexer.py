from mantid.simpleapi import LoadIsawPeaks, LoadIsawUB

class SatellitePeakIndexer:

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

        self.view.cluster_button.clicked.connect(self.cluster)

        peaks = LoadIsawPeaks('/SNS/TOPAZ/IPTS-23996/shared/zgf/satellite_index/26079_Niggli.integrate')
        LoadIsawUB(peaks, '/SNS/TOPAZ/IPTS-23996/shared/zgf/satellite_index/26079_Niggli.mat')

        self.model.set_peak_workspace('peaks')
        self.view.set_transform(self.model.get_transform())
        
        self.cluster()

    def cluster(self):

        params = self.view.get_cluster_parameters()

        if params is not None:

            peak_info = self.model.get_peak_info()

            self.model.cluster_peaks(peak_info, *params)

            self.view.add_peaks(peak_info)
            
            self.view.update_table(peak_info)

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