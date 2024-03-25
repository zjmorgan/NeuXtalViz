from mantid.simpleapi import (LoadMD,
                              IntegrateMDHistoWorkspace,
                              mtd)

import numpy as np
import scipy.linalg

from NeuXtalViz.models.base_model import NeuXtalVizModel

class VolumeSlicerModel(NeuXtalVizModel):

    def __init__(self):

        super(VolumeSlicerModel, self).__init__()

        self.signal = None

    def load_md_histo_workspace(self, filename):

        LoadMD(Filename=filename, OutputWorkspace='histo')

        self.set_B()
        self.set_W()

    def is_histo_loaded(self):

        return mtd.doesExist('histo')

    def set_B(self):

        if self.has_UB('histo'):

            ei = mtd['histo'].getExperimentInfo(0)

            B = ei.sample().getOrientedLattice().getB().copy()

            self.set_UB(B)

    def set_W(self):

        ei = mtd['histo'].getExperimentInfo(0)

        self.W = ei.run().getLogData('W_MATRIX').value.reshape(3,3)

    def get_histo_info(self):

        histo_dict = {}

        self.signal = mtd['histo'].getSignalArray().copy()

        self.signal[self.signal <= 0] = np.nan
        self.signal[np.isinf(self.signal)] = np.nan

        histo_dict['signal'] = self.signal

        dims = [mtd['histo'].getDimension(i) for i in range(3)]

        min_lim = np.array([dim.getMinimum() for dim in dims])
        max_lim = np.array([dim.getMaximum() for dim in dims])

        spacing = np.array([dim.getX(1)-dim.getX(0) for dim in dims])

        min_lim += spacing*0.5
        max_lim -= spacing*0.5

        labels = ['{} ({})'.format(dim.name,dim.getUnits()) for dim in dims]

        histo_dict['min_lim'] = min_lim
        histo_dict['max_lim'] = max_lim
        histo_dict['spacing'] = spacing
        histo_dict['labels'] = labels

        P, T, S = self.get_transforms()

        histo_dict['transform'] = T
        histo_dict['projection'] = P
        histo_dict['scales'] = S

        return histo_dict

    def calculate_clim(self, method='normal'):

        if self.signal is not None:

            trans = np.log10(self.signal)

            vmin, vmax = np.nanmin(trans), np.nanmax(trans)

            if method == 'normal':

                mu, sigma = np.nanmean(trans), np.nanstd(trans)

                spread = 3*sigma

                cmin, cmax = mu-spread, mu+spread

            elif method == 'boxplot':

                Q1, Q3 = np.nanpercentile(trans, [25,75])

                IQR = Q3-Q1

                spread = 1.5*IQR

                cmin, cmax = Q1-spread, Q3+spread

            else:

                cmin, cmax = vmin, vmax

            clim = [cmin if cmin > vmin else vmin,
                    cmax if cmax < vmax else vmax]

            return clim

    def get_transform(self):

        if self.UB is not None:

            b = self.UB/np.linalg.norm(self.UB, axis=0)

            Bp = np.dot(self.UB, self.W)

            Q, R = scipy.linalg.qr(Bp)

            v = scipy.linalg.cholesky(np.dot(R.T, R), lower=False)

            Q = np.dot(Bp, np.linalg.inv(v))

            return np.dot(Q.T, b)

    def get_transforms(self):

        Bp = np.dot(self.UB, self.W)

        Q, R = scipy.linalg.qr(Bp)

        v = scipy.linalg.cholesky(np.dot(R.T, R), lower=False)

        s = np.linalg.norm(v, axis=0)
        t = v/s
        p = v/v[0,0]

        s = np.linalg.norm(p, axis=0) 

        return p, t, s

    def get_normal(self, axes_type, ind):

        if self.UB is not None:

            if axes_type == '[hkl]':
                matrix = self.UB
            else:
                matrix = np.cross(np.dot(self.UB, np.roll(np.eye(3),2,1)).T,
                                  np.dot(self.UB, np.roll(np.eye(3),1,1)).T).T

            vec = np.dot(matrix, ind)

            return vec
