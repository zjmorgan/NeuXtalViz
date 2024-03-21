import numpy as np

class NeuXVizWidgetModel():

    def __init__(self):

        self.UB = None

    def set_UB(self, UB):
        """
        Update the UB-matrix.

        Parameters
        ----------
        UB : 3x3 element 2d array
            UB-matrix.

        """

        self.UB = UB

    def get_transform(self):
        """
        Transformation matrix describing the crystal axes.

        Returns
        -------
        T : 3x3 element 2d array
            Normalized transformation matrix.

        """

        if self.UB is not None:

            T = self.UB/np.linalg.norm(self.UB , axis=0)

            return T

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

            return np.dot(self.UB, [0,0,1]), np.dot(self.UB, [1,0,0])

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

            return np.dot(self.UB, [1,0,0]), np.dot(self.UB, [0,1,0])

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

            return np.dot(self.UB, [0,1,0]), np.dot(self.UB, [0,0,1])

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

            return np.cross(*self.bc_star_axes()), \
                   np.cross(*self.ca_star_axes())

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

            return np.cross(*self.ca_star_axes()), \
                   np.cross(*self.ab_star_axes())

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

            return np.cross(*self.ab_star_axes()), \
                   np.cross(*self.bc_star_axes())

    def get_vector(self, axes_type, ind):
        """
        Vector corresponding to a particular crystallographic direction.

        Parameters
        ----------
        axes_type : str, [hkl] or [uvw]
            Miller index or fractional coordinate.
        ind : 3d-element 1d array-like
            Indices.

        Returns
        -------
        vec : 3 element 1d array
            Cartesian vector.

        """

        if self.UB is not None:

            if axes_type == '[hkl]':
                matrix = self.UB
            else:
                matrix = np.cross(np.dot(self.UB, np.roll(np.eye(3),2,1)).T,
                                  np.dot(self.UB, np.roll(np.eye(3),1,1)).T).T

            vec = np.dot(matrix, ind)

            return vec