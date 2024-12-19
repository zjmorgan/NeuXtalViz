import os

from mantid.simpleapi import (
    CreatePeaksWorkspace,
    LoadIsawUB,
    LoadIsawPeaks,
    LoadNexus,
    SetUB,
    CalculatePeaksHKL,
    mtd,
)

import numpy as np
import scipy.spatial

from sklearn.cluster import DBSCAN

from NeuXtalViz.models.base_model import NeuXtalVizModel


class ModulationModel(NeuXtalVizModel):
    def __init__(self):
        super(ModulationModel, self).__init__()

        CreatePeaksWorkspace(
            OutputType="LeanElasticPeak",
            NumberOfPeaks=0,
            OutputWorkspace="peaks",
        )

    def load_UB(self, filename):
        LoadIsawUB(InputWorkspace="peaks", Filename=filename)

        CalculatePeaksHKL(PeaksWorkspace="peaks", OverWrite=True)

        self.copy_UB()

    def copy_UB(self):
        if self.has_UB("peaks"):
            UB = mtd["peaks"].sample().getOrientedLattice().getUB().copy()

            self.set_UB(UB)

    def load_peaks(self, filename):
        _, ext = os.path.splitext(filename)

        if ext != ".nxs":
            LoadIsawPeaks(Filename=filename, OutputWorkspace="peaks")
        else:
            LoadNexus(Filename=filename, OutputWorkspace="peaks")

        self.copy_UB()

        UB = self.UB

        if UB is not None:
            SetUB(Workspace="peaks", UB=UB)

        if self.has_UB("peaks"):
            CalculatePeaksHKL(PeaksWorkspace="peaks", OverWrite=True)

    def cluster_peaks(self, peak_info, eps=0.025, min_samples=15):
        T_inv = peak_info["inverse"]

        points = np.array(peak_info["coordinates"])

        clustering = DBSCAN(eps=eps, min_samples=min_samples)

        labels = clustering.fit_predict(points)

        uni_labels, inverse = np.unique(labels, return_inverse=True)

        centroids = []
        for label in uni_labels:
            if label >= 0:
                center = points[labels == label].mean(axis=0)
                centroids.append(np.dot(T_inv, center))
        centroids = np.array(centroids)

        success = False

        if centroids.shape[0] >= 0 and len(centroids.shape) == 2:
            null = np.argmin(np.linalg.norm(centroids, axis=1))

            mask = np.ones_like(centroids[:, 0], dtype=bool)
            mask[null] = False

            peaks = np.arange(mask.size)[mask]

            satellites = centroids[mask]
            nuclear = centroids[null]

            dist = scipy.spatial.distance_matrix(satellites, -satellites)

            n = dist.shape[0]

            if n > 2:
                success = True

                indices = np.column_stack(
                    [np.arange(n), np.argmin(dist, axis=0)]
                )
                indices = np.sort(indices, axis=1)
                indices = np.unique(indices, axis=0)

                clusters = labels.copy()
                clusters[labels == null] = 0

                mod = 1
                satellites = []
                for inds in indices:
                    i, j = peaks[inds[0]], peaks[inds[1]]
                    clusters[labels == i] = mod
                    clusters[labels == j] = mod
                    satellites.append(centroids[i])
                    mod += 1
                satellites = np.array(satellites)

                peak_info["clusters"] = clusters
                peak_info["nuclear"] = nuclear
                peak_info["satellites"] = satellites

        return success

    def get_peak_info(self):
        UB = self.UB

        if UB is not None:
            peak_dict = {}

            Qs, HKLs, pk_nos = [], [], []

            for j, peak in enumerate(mtd["peaks"]):
                pk_no = peak.getPeakNumber()

                diff_HKL = peak.getHKL() - np.round(peak.getHKL())

                Q = 2 * np.pi * np.dot(UB, diff_HKL)

                Qs.append(Q)
                HKLs.append(diff_HKL)
                pk_nos.append(pk_no)

                diff_HKL = peak.getHKL() - np.round(peak.getHKL())

                Q = 2 * np.pi * np.dot(UB, -diff_HKL)

                Qs.append(Q)
                HKLs.append(diff_HKL)
                pk_nos.append(-pk_no)

            peak_dict["coordinates"] = Qs
            peak_dict["points"] = HKLs
            peak_dict["numbers"] = pk_nos

            translation = (
                2 * np.pi * UB[:, 0],
                2 * np.pi * UB[:, 1],
                2 * np.pi * UB[:, 2],
            )

            peak_dict["translation"] = translation

            T = np.column_stack(translation)

            peak_dict["transform"] = T
            peak_dict["inverse"] = np.linalg.inv(T)

            return peak_dict

    def get_peak(self, pk_no):
        cols = mtd["peaks"].getColumnNames()

        col = cols.index("PeakNumber")
        row = mtd["peaks"].column(col).index(pk_no)

        return mtd["peaks"].row(row)
