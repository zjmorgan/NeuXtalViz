import os
import itertools

from mantid.simpleapi import (CreatePeaksWorkspace,
                              LoadIsawUB,
                              LoadEmptyInstrument,
                              LoadInstrument,
                              ExtractMonitors,
                              PreprocessDetectorsToMD,
                              AddSampleLog,
                              CloneWorkspace,
                              mtd)

from mantid.geometry import PointGroupFactory

import numpy as np
import skimage

from NeuXtalViz.models.base_model import NeuXtalVizModel
from NeuXtalViz.config.instruments import beamlines

class ExperimentModel(NeuXtalVizModel):

    def __init__(self):

        super(ExperimentModel, self).__init__()

        CreatePeaksWorkspace(OutputType='LeanElasticPeak',
                             NumberOfPeaks=0,
                             OutputWorkspace='coverage')

        self.hist = None

    def load_UB(self, filename):

        LoadIsawUB(InputWorkspace='coverage', Filename=filename)

        self.copy_UB()

    def copy_UB(self):

        if self.has_UB('coverage'):

            UB = mtd['coverage'].sample().getOrientedLattice().getUB().copy()

            self.set_UB(UB)

    def get_instrument_name(self, instrument):

        return beamlines[instrument]['Name']

    def get_goniometers(self, instrument):

        goniometers = beamlines[instrument]['Goniometer']

        return [(name, *goniometers[name][-2:]) for name in goniometers.keys()]

    def get_motors(self, instrument):

        motors = beamlines[instrument].get('Motor')

        if motors is not None:
            return [(name, motors[name]) for name in motors.keys()]
        else:
            return []

    def get_wavelength(self, instrument):

        return beamlines[instrument]['Wavelength']

    def get_symmetry_transforms(self, laue):

        if laue == 'None':
            return [np.eye(3)]
        else:
            pg = PointGroupFactory.createPointGroup(laue)

            coord = np.eye(3).astype(int)

            transforms = []
            for symop in pg.getSymmetryOperations():
                T = np.column_stack([symop.transformHKL(vec) for vec in coord])
                transforms.append(T)

            return transforms

    def generate_instrument_coverage(self, instrument, logs, wavelength):

        LoadEmptyInstrument(InstrumentName=instrument,
                            OutputWorkspace='instrument')

        for key in logs.keys():
            AddSampleLog(Workspace='instrument',
                         LogName=key,
                         LogText=str(logs[key]),
                         LogType='Number Series',
                         NumberType='Double')

        LoadInstrument(Workspace='instrument',
                       RewriteSpectraMap=False,
                       InstrumentName=instrument)

        ExtractMonitors(InputWorkspace='instrument',
                        DetectorWorkspace='instrument',
                        MonitorWorkspace='montitors')

        PreprocessDetectorsToMD(InputWorkspace='instrument',
                                OutputWorkspace='detectors',
                                GetMaskState=False)

        two_theta = np.array(mtd['detectors'].column('TwoTheta'))
        az_phi = np.array(mtd['detectors'].column('Azimuthal'))

        wl_min, wl_max = wavelength

        Q_max = 4*np.pi/wl_min*np.sin(0.5*np.max(two_theta))
        Q_min = 4*np.pi/wl_max*np.sin(0.5*np.min(two_theta))

        self.Q_min, self.Q_max = Q_min, Q_max

        n = 100
        self.Q_bins = np.linspace(-Q_max, Q_max, n+1)
        self.hist = np.zeros((n,n,n))

        kx_hat = np.sin(two_theta)*np.cos(az_phi)
        ky_hat = np.sin(two_theta)*np.sin(az_phi)
        kz_hat = np.cos(two_theta)-1

        Qx_1 = 2*np.pi/wl_max*kx_hat
        Qy_1 = 2*np.pi/wl_max*ky_hat
        Qz_1 = 2*np.pi/wl_max*kz_hat

        Qx_2 = 2*np.pi/wl_min*kx_hat
        Qy_2 = 2*np.pi/wl_min*ky_hat
        Qz_2 = 2*np.pi/wl_min*kz_hat

        Qx_1_ind = np.digitize(Qx_1, self.Q_bins, right=True)
        Qy_1_ind = np.digitize(Qy_1, self.Q_bins, right=True)
        Qz_1_ind = np.digitize(Qz_1, self.Q_bins, right=True)

        Qx_2_ind = np.digitize(Qx_2, self.Q_bins, right=True)
        Qy_2_ind = np.digitize(Qy_2, self.Q_bins, right=True)
        Qz_2_ind = np.digitize(Qz_2, self.Q_bins, right=True)

        _, indices = np.unique(np.column_stack([Qx_1_ind,Qy_1_ind,Qz_1_ind,
                                                Qx_2_ind,Qy_2_ind,Qz_2_ind]),
                               axis=0,
                               return_index=True)

        for i in indices:
            x1, y1, z1 = Qx_1_ind[i], Qy_1_ind[i], Qz_1_ind[i]
            x2, y2, z2 = Qx_2_ind[i], Qy_2_ind[i], Qz_2_ind[i]
            ix, iy, iz = skimage.draw.line_nd([x1, y1, z1],
                                              [x2, y2, z2], endpoint=False)
            self.hist[ix,iy,iz] = 1

    def apply_symmetry(self, laue):

        mask = self.hist > 0

        Q_centers = 0.5*(self.Q_bins[1:]+self.Q_bins[:-1])

        Qx_centers, Qy_centers, Qz_centers = np.meshgrid(Q_centers,
                                                         Q_centers,
                                                         Q_centers)

        x, y, z = Qx_centers[mask], Qy_centers[mask], Qz_centers[mask]

        if self.has_UB('coverage') and self.hist is not None:

            UB = mtd['coverage'].sample().getOrientedLattice().getUB().copy()
            UB_inv = np.linalg.inv(UB)

            for T in self.get_symmetry_transforms(laue):
                xp, yp, zp = np.einsum('ij,jk->ik', (UB @ T) @ UB_inv, [x,y,z])
                Qx_ind = np.digitize(xp, self.Q_bins, right=True)
                Qy_ind = np.digitize(yp, self.Q_bins, right=True)
                Qz_ind = np.digitize(zp, self.Q_bins, right=True)
                self.hist[Qx_ind,Qy_ind,Qz_ind] = 1

    def hsl_to_rgb(self, hue, saturation, lightness):
        h = np.array(hue)
        s = np.array(saturation)
        l = np.array(lightness)

        def f(h, s, l, n):
            k = (n+h/30) % 12
            a = s*np.minimum(l, 1-l)
            return l-a*np.maximum(-1, np.minimum(np.minimum(k-3, 9-k), 1))

        rgb = np.stack((f(h, s, l, 0), f(h, s, l, 8), f(h, s, l, 4)), axis=-1)
        return rgb

    def get_coverage_info(self):

        coverage_dict = {}

        coverage_dict['signal'] = self.hist
        coverage_dict['spacing'] = (self.Q_bins[1]-self.Q_bins[0],)*3

        UB = mtd['coverage'].sample().getOrientedLattice().getUB().copy()
        UB_inv = np.linalg.inv(UB)

        Q_centers = 0.5*(self.Q_bins[1:]+self.Q_bins[:-1])

        x, y, z = np.meshgrid(Q_centers, Q_centers, Q_centers)

        h, k, l = np.einsum('ij,jk->ik', UB_inv, [x.ravel(),
                                                  y.ravel(),
                                                  z.ravel()])

        r = np.sqrt(h**2+k**2+l**2)
        theta = np.arccos(l/r)
        phi = np.arctan2(k, h)

        hue = phi*180/np.pi+180
        saturation = np.ones_like(hue)
        lightness = theta/np.pi

        rgb = self.hsl_to_rgb(hue, saturation, lightness) 

        Q = np.sqrt(x**2+y**2+z**2)
        
        a = (Q-self.Q_min)/(self.Q_max-self.Q_min)*255
        a[self.hist == 0] = 0

        limits = [-self.Q_max, self.Q_max]*3

        rgb = np.array(rgb)*255
        rgba = np.column_stack([rgb, a.ravel()]).astype(np.uint8)

        coverage_dict['scalars'] = rgba

        points = list(itertools.product([-1, 0, 1], repeat=3))
        points = [point for point in points if point != (0, 0, 0)]

        labels = ['{}{}{}'.format(*point) for point in points]

        points = [np.dot(UB, point) for point in points]
        points = [self.Q_max*point/np.linalg.norm(point) for point in points]

        coverage_dict['labels'] = labels
        coverage_dict['points'] = points
        coverage_dict['limits'] = limits

        return coverage_dict