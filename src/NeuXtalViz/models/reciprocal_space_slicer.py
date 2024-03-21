from mantid.simpleapi import (CreatePeaksWorkspace,
                              HasUB,
                              mtd)

import numpy as np
import scipy.linalg

class ReciprocalSpaceSlicerModel():

    def __init__(self):

        self.B = None
        self.signal = None

    def set_md_histo_workspace(self, md_histo_ws):

        self.md_histo_ws = mtd[md_histo_ws]

        self.set_B()
        self.set_W()

    def set_B(self):

        if HasUB(Workspace=self.md_histo_ws):

            ei = self.md_histo_ws.getExperimentInfo(0)

            self.B = ei.sample().getOrientedLattice().getB().copy()

    def set_W(self):

        ei = self.md_histo_ws.getExperimentInfo(0)

        self.W = ei.run().getLogData('W_MATRIX').value.reshape(3,3)

    def get_histo_info(self):

        histo_dict = {}

        self.signal = self.md_histo_ws.getSignalArray().copy()

        self.signal[self.signal <= 0] = np.nan
        self.signal[np.isinf(self.signal)] = np.nan

        histo_dict['signal'] = self.signal

        dims = [self.md_histo_ws.getDimension(i) for i in range(3)]

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

    def calculate_clim(self, method='3-sig'):

        if self.signal is not None:

            trans = np.log10(self.signal)

            vmin, vmax = np.nanmin(trans), np.nanmax(trans)

            if method == '3-sig':

                mu, sigma = np.nanmean(trans), np.nanstd(trans)

                spread = 3*sigma

                cmin, cmax = mu-spread, mu+spread

            elif method == '1.5-IQR':

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

        if self.B is not None:

            b = self.B/np.linalg.norm(self.B, axis=0)

            Bp = np.dot(self.B, self.W)

            Q, R = scipy.linalg.qr(Bp)

            v = scipy.linalg.cholesky(np.dot(R.T, R), lower=False)

            Q = np.dot(Bp, np.linalg.inv(v))

            return np.dot(Q.T, b)

    def get_transforms(self):

        Bp = np.dot(self.B, self.W)

        Q, R = scipy.linalg.qr(Bp)

        v = scipy.linalg.cholesky(np.dot(R.T, R), lower=False)

        s = np.linalg.norm(v, axis=0)
        t = v/s
        p = v/v[0,0]

        s = np.linalg.norm(p, axis=0) 

        return p, t, s

    def ab_star_axes(self):

        if self.B is not None:

            return np.dot(self.B, [0,0,1]), np.dot(self.B, [1,0,0])

    def bc_star_axes(self):

        if self.B is not None:

            return np.dot(self.B, [1,0,0]), np.dot(self.B, [0,1,0])

    def ca_star_axes(self):

        if self.B is not None:

            return np.dot(self.B, [0,1,0]), np.dot(self.B, [0,0,1])

    def ab_axes(self):

        if self.B is not None:

            return np.cross(*self.bc_star_axes()), \
                   np.cross(*self.ca_star_axes())

    def bc_axes(self):

        if self.B is not None:

            return np.cross(*self.ca_star_axes()), \
                   np.cross(*self.ab_star_axes())

    def ca_axes(self):

        if self.B is not None:

            return np.cross(*self.ab_star_axes()), \
                   np.cross(*self.bc_star_axes())

    def get_vector(self, axes_type, ind):

        if self.B is not None:

            if axes_type == '[hkl]':
                matrix = self.B
            else:
                matrix = np.cross(np.dot(self.B, np.roll(np.eye(3),2,1)).T,
                                  np.dot(self.B, np.roll(np.eye(3),1,1)).T).T

            vec = np.dot(matrix, ind)

            return vec