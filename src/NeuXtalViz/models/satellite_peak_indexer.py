from mantid.simpleapi import (CreatePeaksWorkspace,
                              HasUB,
                              CalculatePeaksHKL,
                              mtd)

import numpy as np
from sklearn.cluster import DBSCAN

class SatellitePeakIndexerModel():

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

            self.UB = self.peaks_ws.sample().getOrientedLattice().getUB()

    def cluster_peaks(self, peak_info, eps=0.025, min_samples=15):

        points = np.array(peak_info['coordinates'])

        clustering = DBSCAN(eps=eps, min_samples=min_samples)

        labels = clustering.fit_predict(points)

        centroids = []
        for label in np.unique(labels):
            if label >= 0:
                centroids.append(points[labels == label].mean(axis=0))
        centroids = np.array(centroids)

        peak_info['clusters'] = labels
        peak_info['centroids'] = centroids

    def get_peak_info(self):

        peak_dict = {}

        Qs, HKLs, pk_nos = [], [], []

        UB = self.peaks_ws.sample().getOrientedLattice().getUB()

        CalculatePeaksHKL(PeaksWorkspace=self.peaks_ws, OverWrite=True)

        for j, peak in enumerate(self.peaks_ws):

            pk_no = peak.getPeakNumber()

            diff_HKL = peak.getHKL()-np.round(peak.getHKL())

            if np.sum(diff_HKL) < 0:
                diff_HKL *= -1

            Q = 2*np.pi*np.dot(UB, diff_HKL)

            Qs.append(Q)
            HKLs.append(diff_HKL)
            pk_nos.append(pk_no)

        peak_dict['coordinates'] = Qs
        peak_dict['points'] = HKLs
        peak_dict['numbers'] = pk_nos

        return peak_dict

    def get_transform(self):

        if self.UB is not None:

            UB = self.peaks_ws.sample().getOrientedLattice().getUB()

            t = UB.copy()
            t /= np.max(t, axis=1)

            return t

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