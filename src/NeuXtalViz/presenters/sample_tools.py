import numpy as np

from NeuXtalViz.presenters.base_presenter import NeuXtalVizPresenter

class Sample(NeuXtalVizPresenter):

    def __init__(self, view, model):

        super(Sample, self).__init__(view, model)

        self.view.connect_row_highligter(self.highlight_row)
        self.view.connect_sample_parameters(self.update_parameters)
        self.view.connect_goniometer_table(self.set_goniometer_table)

        self.view.connect_load_UB(self.load_UB)
        self.view.connect_add_sample(self.add_sample)

    def highlight_row(self):

        goniometer = self.view.get_goniometer()
        self.view.set_angle(goniometer)

    def set_goniometer_table(self):

        self.view.set_goniometer_table()

    def load_UB(self):

        filename = self.view.load_UB_file_dialog()

        if filename:

            self.model.load_UB(filename)
            vol = self.model.get_volume()
            self.view.set_unit_cell_volume(vol)
            self.update_oriented_lattice()

    def update_parameters(self):

        params = self.view.get_sample_constants()

        shape = self.view.get_sample_shape()

        fixed = [False, False, False]

        if shape != 'Plate':
            fixed[2] = True
            params[2] = params[0]
        if shape != 'Cylinder' and shape != 'Plate':
            fixed[1] = True
            params[1] = params[0]

        self.view.set_sample_constants(params)
        self.view.constrain_size(fixed)

    def add_sample(self):

        mat_dict = shape_dict = None

        goniometers = self.view.get_goniometers()
        axes = self.model.get_goniometer_strings(goniometers)

        material_params = self.view.get_material_paremters()
        if material_params is not None:
            mat_dict = self.model.get_material_dict(*material_params)

        sample_params = self.view.get_sample_constants()

        shape = self.view.get_sample_shape()

        if sample_params is not None:
            vectors = self.view.get_face_indexing()
            angles = 0, 0, 0
            if vectors is not None:
                values = self.model.get_euler_angles(vectors[0], vectors[1])
                if values is not None:
                    angles = values
            shape_dict = self.model.get_shape_dict(shape,
                                                   sample_params,
                                                   *angles)

        if np.all([var is not None for var in (shape_dict, mat_dict, axes)]):
            self.model.set_sample(shape_dict, mat_dict, axes)

            abs_dict = self.model.get_absorption_dict()
            self.view.set_absortion_parameters(abs_dict)
            mesh = self.model.sample_mesh()
            self.view.add_sample(mesh)
            self.view.set_transform(self.model.get_transform())
