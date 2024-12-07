from mantid.simpleapi import (LoadMD,
                              IntegrateMDHistoWorkspace,
                              mtd)

import numpy as np
import scipy.linalg

import skimage.measure

from NeuXtalViz.models.base_model import NeuXtalVizModel
from NeuXtalViz.models.utilities import SaveMDToAscii

class VolumeSlicerModel(NeuXtalVizModel):

    def __init__(self):

        super(VolumeSlicerModel, self).__init__()

    def load_md_histo_workspace(self, filename):

        LoadMD(Filename=filename, OutputWorkspace='histo')

        signal = mtd['histo'].getSignalArray().copy()

        self.shape = signal.shape

        dims = [mtd['histo'].getDimension(i) for i in range(3)]

        self.min_lim = np.array([dim.getMinimum()+\
                                 dim.getBinWidth()*0.5 for dim in dims])

        self.max_lim = np.array([dim.getMaximum()-\
                                 dim.getBinWidth()*0.5 for dim in dims])

        self.labels = ['{} ({})'.format(dim.name,
                                        dim.getUnits()) for dim in dims]

        self.spacing = np.array([dim.getBinWidth() for dim in dims])

        scale = 0.25/self.spacing
        scale[scale <= 1] = 1
        scale = scale.round().astype(int)

        blocks = [(scale[0],1,1),
                  (1,scale[1],1),
                  (1,1,scale[2])]

        self.signals = []
        self.spacings = []
        for block in blocks:
            self.spacings.append(self.spacing*np.array(block))
            self.signals.append(skimage.measure.block_reduce(signal, 
                                                             block_size=block,
                                                             func=np.nanmean,
                                                             cval=np.nan))

        self.set_B()
        self.set_W()

    def save_slice(self, filename):

        SaveMDToAscii('slice', filename)

    def save_cut(self, filename):

        SaveMDToAscii('cut', filename)

    def is_histo_loaded(self):

        return mtd.doesExist('histo')

    def is_sliced(self):

        return mtd.doesExist('slice')

    def is_cut(self):

        return mtd.doesExist('cut')

    def set_B(self):

        if self.has_UB('histo'):

            ei = mtd['histo'].getExperimentInfo(0)

            B = ei.sample().getOrientedLattice().getB().copy()

            self.set_UB(B)

    def set_W(self):

        ei = mtd['histo'].getExperimentInfo(0)

        self.W = np.eye(3)

        if ei.run().hasProperty('W_MATRIX'):

            self.W = ei.run().getLogData('W_MATRIX').value.reshape(3,3)

    def get_histo_info(self, normal):

        ind = normal.index(1)

        histo_dict = {}

        histo_dict['signal'] = self.signals[ind].copy()

        histo_dict['min_lim'] = self.min_lim
        histo_dict['max_lim'] = self.max_lim
        histo_dict['spacing'] = self.spacings[ind]
        histo_dict['labels'] = self.labels

        P, T, S = self.get_transforms()

        histo_dict['transform'] = T
        histo_dict['projection'] = P
        histo_dict['scales'] = S

        return histo_dict

    def get_slice_info(self, normal, value, thickness=0.01):

        self.normal = normal

        slice_dict = {}

        integrate = [value-thickness, value+thickness]

        self.integrate = integrate

        pbin = [None if norm == 0 else integrate for norm in normal]

        IntegrateMDHistoWorkspace(InputWorkspace='histo',
                                  P1Bin=pbin[0],
                                  P2Bin=pbin[1],
                                  P3Bin=pbin[2],
                                  OutputWorkspace='slice')

        i = np.array(normal).tolist().index(1)

        form = '{} = ({:.2f},{:.2f})'

        title = form.format(mtd['slice'].getDimension(i).name, *integrate)

        dims = mtd['slice'].getNonIntegratedDimensions()

        x, y = [np.linspace(dim.getMinimum(),
                            dim.getMaximum(),
                            dim.getNBoundaries()) for dim in dims]

        labels = ['{} ({})'.format(dim.name, dim.getUnits()) for dim in dims]

        slice_dict['x'] = x
        slice_dict['y'] = y
        slice_dict['labels'] = labels

        signal = mtd['slice'].getSignalArray().T.copy().squeeze()

        signal[signal <= 0] = np.nan
        signal[np.isinf(signal)] = np.nan

        slice_dict['signal'] = signal

        Bp = np.dot(self.UB, self.W)

        Q, R = scipy.linalg.qr(Bp)

        ind = np.array(normal) != 1

        v = scipy.linalg.cholesky(np.dot(R.T, R)[ind][:,ind], lower=False)

        v /= v[0,0]

        T = np.eye(3)
        T[:2,:2] = v

        s = np.diag(T).copy()
        T[1,1] = 1

        T[0,2] = -T[0,1]*y.min()

        slice_dict['transform'] = T
        slice_dict['aspect'] = s[1]
        slice_dict['value'] = value
        slice_dict['title'] = title

        return slice_dict

    def get_cut_info(self, axis, value, thickness=0.01):

        cut_dict = {}

        integrate = [value-thickness, value+thickness]

        pbin = [None if ax == 0 else integrate for ax in axis]

        IntegrateMDHistoWorkspace(InputWorkspace='slice',
                                  P1Bin=pbin[0],
                                  P2Bin=pbin[1],
                                  P3Bin=pbin[2],
                                  OutputWorkspace='cut')

        i = np.array(self.normal).tolist().index(1)
        j = np.array(axis).tolist().index(1)

        form = '{} = ({:.2f},{:.2f})'

        title = form.format(mtd['slice'].getDimension(i).name, *self.integrate)
        title += ' / '
        title += form.format(mtd['cut'].getDimension(j).name, *integrate)

        dim = mtd['cut'].getNonIntegratedDimensions()[0]

        x = np.linspace(dim.getMinimum(), 
                        dim.getMaximum(),
                        dim.getNBoundaries())

        x = 0.5*(x[1:]+x[:-1])

        label = '{} ({})'.format(dim.name, dim.getUnits())

        cut_dict['x'] = x
        cut_dict['y'] = mtd['cut'].getSignalArray().squeeze()
        cut_dict['e'] = np.sqrt(mtd['cut'].getErrorSquaredArray().squeeze())
        cut_dict['label'] = label
        cut_dict['value'] = value
        cut_dict['title'] = title

        return cut_dict

    def calculate_clim(self, trans, method='normal'):

        trans[~np.isfinite(trans)] = np.nan

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

        trans[trans < clim[0]] = clim[0]
        trans[trans > clim[1]] = clim[1]

        return trans

    def get_transform(self, reciprocal=True):

        if self.UB is not None:

            b = self.UB/np.linalg.norm(self.UB, axis=0)

            Bp = np.dot(self.UB, self.W)

            Q, R = scipy.linalg.qr(Bp)

            v = scipy.linalg.cholesky(np.dot(R.T, R), lower=False)

            Q = np.dot(Bp, np.linalg.inv(v))

            T = np.dot(Q.T, b)

            if not reciprocal:

                T = np.column_stack([np.cross(T[:,1], T[:,2]),
                                     np.cross(T[:,2], T[:,0]),
                                     np.cross(T[:,0], T[:,1])])

            return T

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