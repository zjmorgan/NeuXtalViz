from mantid.simpleapi import HasUB

from mantid.geometry import OrientedLattice

import numpy as np


class NeuXtalVizModel:
    def __init__(self):
        self.UB = None

    def has_UB(self, ws):
        """
        Check if the oriented lattice exists on a workspace.

        Parameters
        ----------
        ws : str
            Name of workspace.

        Returns
        -------
        ol : bool
            Oriented lattice exists or not.

        """

        ol = HasUB(Workspace=ws)

        return ol

    def set_UB(self, UB):
        """
        Update the UB-matrix.

        Parameters
        ----------
        UB : 3x3 element 2d array
            UB-matrix.

        """

        self.UB = UB.copy()

    def get_oriented_lattice_parameters(self):
        """
        Obtain the oriented lattice paramters.

        Returns
        -------
        a, b, c : float
            Lattice constants.
        alpha, beta, gamma : float
            Lattice angles.

        """

        if self.UB is not None:
            ol = OrientedLattice()
            ol.setUB(self.UB)

            params = ol.a(), ol.b(), ol.c(), ol.alpha(), ol.beta(), ol.gamma()

            u, v = ol.getuVector(), ol.getvVector()

            u = list(np.array(u) / np.abs(u).max())
            v = list(np.array(v) / np.abs(v).max())

            return *params, u, v

    def orientation_matrix(self):
        return self.UB.copy()

    def get_transform(self, reciprocal=True):
        """
        Transformation matrix describing the reciprocal or crystal axes.

        Parameters
        ----------
        reciprocal : bool
            Option for the reciprocal or crystal lattice axes.
            Default is ``True``.

        Returns
        -------
        T : 3x3 element 2d array
            Normalized transformation matrix.

        """

        if self.UB is not None:
            if reciprocal:
                T = self.orientation_matrix()
            else:
                T = np.column_stack(
                    [
                        np.cross(self.UB[:, 1], self.UB[:, 2]),
                        np.cross(self.UB[:, 2], self.UB[:, 0]),
                        np.cross(self.UB[:, 0], self.UB[:, 1]),
                    ]
                )

            return T / np.linalg.norm(T, axis=0)

    def ab_star_axes(self):
        """
        :math:`c`-direction in cartesian coordinates.

        Returns
        -------
        camera : 3 element 1d array
            Cartesian camera view vector.
        updward : 3 element 1d array
            Cartesian upward view vector.

        """

        if self.UB is not None:
            return np.dot(self.orientation_matrix(), [0, 0, 1]), np.dot(
                self.orientation_matrix(), [1, 0, 0]
            )

    def bc_star_axes(self):
        """
        :math:`a`-direction in cartesian coordinates

        Returns
        -------
        camera : 3 element 1d array
            Cartesian camera view vector.
        updward : 3 element 1d array
            Cartesian upward view vector.

        """

        if self.UB is not None:
            return np.dot(self.orientation_matrix(), [1, 0, 0]), np.dot(
                self.orientation_matrix(), [0, 1, 0]
            )

    def ca_star_axes(self):
        """
        :math:`b`-direction in cartesian coordinates

        Returns
        -------
        camera : 3 element 1d array
            Cartesian camera view vector.
        updward : 3 element 1d array
            Cartesian upward view vector.

        """

        if self.UB is not None:
            return np.dot(self.orientation_matrix(), [0, 1, 0]), np.dot(
                self.orientation_matrix(), [0, 0, 1]
            )

    def ab_axes(self):
        """
        :math:`c^\ast`-direction in cartesian coordinates.

        Returns
        -------
        camera : 3 element 1d array
            Cartesian camera view vector.
        updward : 3 element 1d array
            Cartesian upward view vector.

        """

        if self.UB is not None:
            return np.cross(*self.bc_star_axes()), np.cross(
                *self.ca_star_axes()
            )

    def bc_axes(self):
        """
        :math:`a^\ast`-direction in cartesian coordinates.

        Returns
        -------
        camera : 3 element 1d array
            Cartesian camera view vector.
        updward : 3 element 1d array
            Cartesian upward view vector.

        """

        if self.UB is not None:
            return np.cross(*self.ca_star_axes()), np.cross(
                *self.ab_star_axes()
            )

    def ca_axes(self):
        """
        :math:`b^\ast`-direction in cartesian coordinates.

        Returns
        -------
        camera : 3 element 1d array
            Cartesian camera view vector.
        updward : 3 element 1d array
            Cartesian upward view vector.

        """

        if self.UB is not None:
            return np.cross(*self.ab_star_axes()), np.cross(
                *self.bc_star_axes()
            )

    def get_vector(self, axes_type, ind):
        """
        Vector corresponding to a particular crystallographic direction.

        Parameters
        ----------
        axes_type : str, [hkl] or [uvw]
            Miller index or fractional coordinate.
        ind : 3-element 1d array-like
            Indices.

        Returns
        -------
        vec : 3 element 1d array
            Cartesian vector.

        """

        if self.UB is not None:
            if axes_type == "[hkl]":
                matrix = self.orientation_matrix()
            else:
                matrix = np.cross(
                    np.dot(
                        self.orientation_matrix(), np.roll(np.eye(3), 2, 1)
                    ).T,
                    np.dot(
                        self.orientation_matrix(), np.roll(np.eye(3), 1, 1)
                    ).T,
                ).T

            vec = np.dot(matrix, ind)

            return vec
