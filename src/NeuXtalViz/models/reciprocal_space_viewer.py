from mantid.simpleapi import (CreatePeaksWorkspace, mtd)

import numpy as np

from NeuXtalViz.models.base_model import NeuXtalVizModel

class ReciprocalSpaceViewerModel(NeuXtalVizModel):

    def __init__(self):

        super(NeuXtalVizModel, self).__init__()

        peaks_ws = 'rsv'

        CreatePeaksWorkspace(OutputType='LeanElasticPeak',
                             OutputWorkspace=peaks_ws)

        self.set_peak_workspace(peaks_ws)

        self.UB = None

    def set_peak_workspace(self, peaks_ws):

        self.peaks_ws = mtd[peaks_ws]

        if self.has_UB(peaks_ws):

            UB = self.peaks_ws.sample().getOrientedLattice().getUB().copy()

            self.set_UB(UB)

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

    def get_peak(self, pk_no):

        cols = self.peaks_ws.getColumnNames()

        col = cols.index('PeakNumber')
        row = self.peaks_ws.column(col).index(pk_no)

        return self.peaks_ws.row(row)