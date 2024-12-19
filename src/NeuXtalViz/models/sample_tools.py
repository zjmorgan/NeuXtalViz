from mantid.simpleapi import (
    CreateSingleValuedWorkspace,
    SetSample,
    SetGoniometer,
    LoadIsawUB,
    mtd,
)

import numpy as np
import scipy.spatial

from NeuXtalViz.models.base_model import NeuXtalVizModel


class SampleModel(NeuXtalVizModel):
    def __init__(self):
        super(SampleModel, self).__init__()

        CreateSingleValuedWorkspace(OutputWorkspace="sample")

    def load_UB(self, filename):
        LoadIsawUB(InputWorkspace="sample", Filename=filename)

        UB = mtd["sample"].sample().getOrientedLattice().getUB().copy()

        self.set_UB(UB)

    def get_volume(self):
        if self.has_UB("sample"):
            return mtd["sample"].sample().getOrientedLattice().volume()

    def get_euler_angles(self, u_vector, v_vector):
        w_vector = np.cross(u_vector, v_vector)

        if self.UB is not None and np.linalg.norm(w_vector) > 0:
            u = np.dot(self.UB, u_vector)
            v = np.dot(self.UB, v_vector)

            u /= np.linalg.norm(u)

            w = np.cross(u, v)
            w /= np.linalg.norm(w)

            v = np.cross(w, u)

            T = np.column_stack([v, w, u])
            R = scipy.spatial.transform.Rotation.from_matrix(T)

            gamma, beta, alpha = R.as_euler("ZYX", degrees=True)

            return alpha, beta, gamma

    def get_shape_dict(self, shape, params, alpha=0, beta=0, gamma=0):
        if shape == "Sphere":
            radius = params[0] / 200
            shape = ' \
            <sphere id="sphere"> \
              <radius val="{}" /> \
              <centre x="0.0" y="0.0" z="0.0" /> \
              <rotate x="{}" y="{}" z="{}" /> \
            </sphere> \
            '.format(
                radius, alpha, beta, gamma
            )
        elif shape == "Cylinder":
            radius, height = params[0] / 200, params[1] / 100
            shape = ' \
            <cylinder id="cylinder"> \
              <centre-of-bottom-base x="0.0" y="{}" z="0.0" /> \
              <axis x="0.0" y="1.0" z="0" /> \
              <radius val="{}" /> \
              <height val="{}" /> \
              <rotate x="{}" y="{}" z="{}" /> \
            </cylinder> \
            '.format(
                -height / 2, radius, height, alpha, beta, gamma
            )
        else:
            width, height, depth = (
                params[0] / 100,
                params[1] / 100,
                params[2] / 100,
            )
            shape = ' \
            <cuboid id="cuboid"> \
              <width val="{}" /> \
              <height val="{}" /> \
              <depth val="{}" /> \
              <centre x="0.0" y="0.0" z="0.0" /> \
              <rotate x="{}" y="{}" z="{}" /> \
            </cuboid> \
            '.format(
                width, height, depth, alpha, beta, gamma
            )

        return {"Shape": "CSG", "Value": shape}

    def get_material_dict(self, chemical_formula, z_parameter, volume):
        mat_dict = {
            "ChemicalFormula": chemical_formula,
            "ZParameter": z_parameter,
            "UnitCellVolume": volume,
        }

        return mat_dict

    def get_goniometer_strings(self, goniometers):
        axes = []
        for goniometer in goniometers:
            name, x, y, z, sense, angle = goniometer
            if np.linalg.norm([x, y, z]) > 0:
                axes.append("{},{},{},{},{}".format(angle, x, y, z, sense))

        if len(axes) == 3:
            return axes

    def set_sample(self, shape_dict, mat_dict, axes):
        SetGoniometer(
            Workspace="sample", Axis0=axes[0], Axis1=axes[1], Axis2=axes[2]
        )

        SetSample(
            InputWorkspace="sample", Geometry=shape_dict, Material=mat_dict
        )

    def get_absorption_dict(self):
        mat = mtd["sample"].sample().getMaterial()

        sigma_a = mat.absorbXSection()
        sigma_s = mat.totalScatterXSection()

        M = mat.relativeMolecularMass()
        n = mat.numberDensityEffective
        N = mat.totalAtoms

        V = abs(mtd["sample"].sample().getShape().volume() * 100**3)

        rho = (n / N) / 0.6022 * M
        m = rho * V

        mu_s = n * sigma_s
        mu_a = n * sigma_a

        abs_dict = {
            "sigma_a": sigma_a,  # barn
            "sigma_s": sigma_s,  # barn
            "mu_a": mu_a,  # 1/cm
            "mu_s": mu_s,  # 1/cm
            "N": N,  # atoms
            "M": M,  # g/mol
            "n": n,  # 1/A^3
            "rho": rho,  # g/cm^3
            "V": V,  # cm^3
            "m": m,
        }  # g

        return abs_dict

    def sample_mesh(self):
        shape = mtd["sample"].sample().getShape()

        return shape.getMesh() * 100
