from mantid.simpleapi import (CreatePeaksWorkspace,
                              HasUB,
                              mtd)

import numpy as np
import scipy.linalg

class ReciprocalSpaceViewerModel():

    def __init__(self):

        peaks_ws = '__init_rsv'

        CreatePeaksWorkspace(OutputType='LeanElasticPeak',
                             OutputWorkspace=peaks_ws)

        self.set_peak_workspace(peaks_ws)

        self.UB = None

    def set_peak_workspace(self, peaks_ws):

        self.peaks_ws = mtd[peaks_ws]

        self.set_UB()

    def set_UB(self):

        if HasUB(Workspace=self.peaks_ws):

            self.UB = self.peaks_ws.sample().getOrientedLattice().getUB().copy()

    def get_peak_info(self):

        peak_dict = {}

        Qs, Is, pk_nos, Ts = [], [], [], []

        for j, peak in enumerate(self.peaks_ws):

            T = np.zeros((4,4))

            I = peak.getIntensity()

            shape = eval(peak.getPeakShape().toJSON())

            pk_no = peak.getPeakNumber()

            Q = peak.getQSampleFrame()

            r = np.array([shape['radius0'],
                          shape['radius1'],
                          shape['radius2']])

            dir1 = np.array(shape['direction0'].split(' ')).astype(float)
            dir2 = np.array(shape['direction1'].split(' ')).astype(float)
            dir3 = np.array(shape['direction2'].split(' ')).astype(float)

            v = np.column_stack([dir1, dir2, dir3])

            P = np.dot(v, np.dot(np.diag(r), v.T))

            T[:3,:3] = P
            T[:3,-1] = Q
            T[-1,-1] = 1

            Qs.append(Q)
            Is.append(I)
            pk_nos.append(pk_no)
            Ts.append(T)

        peak_dict['coordinates'] = Qs
        peak_dict['intensities'] = Is
        peak_dict['numbers'] = pk_nos
        peak_dict['transforms'] = Ts

        return peak_dict

    def get_transform(self):

        if self.UB is not None:

            UB = self.peaks_ws.sample().getOrientedLattice().getUB()

            ub = UB/np.linalg.norm(UB, axis=0)

            return ub

    def get_peak(self, pk_no):

        cols = self.peaks_ws.getColumnNames()

        col = cols.index('PeakNumber')
        row = self.peaks_ws.column(col).index(pk_no)

        return self.peaks_ws.row(row)

    def ab_star_axes(self):

        if self.UB is not None:

            return np.dot(self.UB, [0,0,1]), np.dot(self.UB, [1,0,0])

    def bc_star_axes(self):

        if self.UB is not None:

            return np.dot(self.UB, [1,0,0]), np.dot(self.UB, [0,1,0])

    def ca_star_axes(self):

        if self.UB is not None:

            return np.dot(self.UB, [0,1,0]), np.dot(self.UB, [0,0,1])

    def ab_axes(self):

        if self.UB is not None:

            return np.cross(*self.bc_star_axes()), \
                   np.cross(*self.ca_star_axes())

    def bc_axes(self):

        if self.UB is not None:

            return np.cross(*self.ca_star_axes()), \
                   np.cross(*self.ab_star_axes())

    def ca_axes(self):

        if self.UB is not None:

            return np.cross(*self.ab_star_axes()), \
                   np.cross(*self.bc_star_axes())

    def get_vector(self, axes_type, ind):

        if self.UB is not None:

            if axes_type == '[hkl]':
                matrix = self.UB
            else:
                matrix = np.cross(np.dot(self.UB, np.roll(np.eye(3),2,1)).T,
                                  np.dot(self.UB, np.roll(np.eye(3),1,1)).T).T

            vec = np.dot(matrix, ind)

            return vec