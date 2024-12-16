import os

import itertools
import csv

from mantid.simpleapi import (CreatePeaksWorkspace,
                              PredictPeaks,
                              FilterPeaks,
                              StatisticsOfPeaksWorkspace,
                              CombinePeaksWorkspaces,
                              SortPeaksWorkspace,
                              AddPeakHKL,
                              SetUB,
                              SetGoniometer,
                              LoadNexus,
                              SaveNexus,
                              LoadIsawUB,
                              LoadEmptyInstrument,
                              LoadInstrument,
                              ExtractMonitors,
                              PreprocessDetectorsToMD,
                              GroupDetectors,
                              AddSampleLog,
                              CreateEmptyTableWorkspace,
                              DeleteTableRows,
                              CloneWorkspace,
                              DeleteWorkspace,
                              RenameWorkspaces,
                              HasUB,
                              mtd)

from mantid.kernel import V3D
from mantid.geometry import PointGroupFactory

from collections import defaultdict

import numpy as np
from scipy.spatial.transform import Rotation

from NeuXtalViz.models.base_model import NeuXtalVizModel
from NeuXtalViz.models.utilities import ParallelTasks
from NeuXtalViz.config.instruments import beamlines

lattice_centering_dict = {
    'P': 'Primitive',
    'C': 'C-face centred',
    'A': 'A-face centred',
    'B': 'B-face centred',
    'I': 'Body centred',
    'F': 'All-face centred',
    'R': 'Primitive',
    'Robv': 'Rhombohedrally centred, obverse',
    'Rrev': 'Rhombohedrally centred, reverse',
}

point_group_dict = {
    '1': '1 (Triclinic)',
    '-1': '-1 (Triclinic)',
    '2': '2 (Monoclinic, unique axis b)',
    'm': 'm (Monoclinic, unique axis b)',
    '2/m': '2/m (Monoclinic, unique axis b)',
    '112': '112 (Monoclinic, unique axis c)',
    '11m': '11m (Monoclinic, unique axis c)',
    '112/m': '112/m (Monoclinic, unique axis c)',
    '222': '222 (Orthorhombic)',
    'mm2': 'mm2 (Orthorhombic)',
    'mmm': 'mmm (Orthorhombic)',
    '4': '4 (Tetragonal)',
    '-4': '-4 (Tetragonal)',
    '4/m': '4/m (Tetragonal)',
    '422': '422 (Tetragonal)',
    '4mm': '4mm (Tetragonal)',
    '-42m': '-42m (Tetragonal)',
    '-4m2': '-4m2 (Tetragonal)',
    '4/mmm': '4/mmm (Tetragonal)',
    '3 r': '3 r (Trigonal - Rhombohedral)',
    '-3 r': '-3 r (Trigonal - Rhombohedral)',
    '32 r': '32 r (Trigonal - Rhombohedral)',
    '3m r': '3m r (Trigonal - Rhombohedral)',
    '-3m r': '-3m r (Trigonal - Rhombohedral)',
    '3': '3 (Trigonal - Hexagonal)',
    '-3': '-3 (Trigonal - Hexagonal)',
    '312': '312 (Trigonal - Hexagonal)',
    '31m': '31m (Trigonal - Hexagonal)',
    '32': '32 (Trigonal - Hexagonal)',
    '321': '321 (Trigonal - Hexagonal)',
    '3m': '3m (Trigonal - Hexagonal)',
    '-31m': '-31m (Trigonal - Hexagonal)',
    '-3m': '-3m (Trigonal - Hexagonal)',
    '-3m1': '-3m1 (Trigonal - Hexagonal)',
    '6': '6 (Hexagonal)',
    '-6': '-6 (Hexagonal)',
    '6/m': '6/m (Hexagonal)',
    '622': '622 (Hexagonal)',
    '6mm': '6mm (Hexagonal)',
    '-62m': '-62m (Hexagonal)',
    '-6m2': '-6m2 (Hexagonal)',
    '6/mmm': '6/mmm (Hexagonal)',
    '23': '23 (Cubic)',
    'm-3': 'm-3 (Cubic)',
    '432': '432 (Cubic)',
    '-43m': '-43m (Cubic)',
    'm-3m': 'm-3m (Cubic)'
}

point_group_centering = {
    '1': ['P'],
    '-1': ['P'],
    '2': ['P', 'C'],
    'm': ['P', 'C'],
    '2/m': ['P', 'C'],
    '112': ['P', 'C'],
    '11m': ['P', 'C'],
    '112/m': ['P', 'C'],
    '222': ['P', 'I', 'C', 'A', 'B'],
    'mm2': ['P', 'I', 'C', 'A', 'B'],
    'mmm': ['P', 'I', 'C', 'A', 'B'],
    '4': ['P', 'I'],
    '-4': ['P', 'I'],
    '4/m': ['P', 'I'],
    '422': ['P', 'I'],
    '4mm': ['P', 'I'],
    '-42m': ['P', 'I'],
    '-4m2': ['P', 'I'],
    '4/mmm': ['P', 'I'],
    '3 r': ['R'],
    '-3 r': ['R'],
    '32 r': ['R'],
    '3m r': ['R'],
    '-3m r': ['R'],
    '3': ['Robv', 'Rrev'],
    '-3': ['Robv', 'Rrev'],
    '312': ['Robv', 'Rrev'],
    '31m': ['Robv', 'Rrev'],
    '32': ['Robv', 'Rrev'],
    '321': ['Robv', 'Rrev'],
    '3m': ['Robv', 'Rrev'],
    '-31m': ['Robv', 'Rrev'],
    '-3m': ['Robv', 'Rrev'],
    '-3m1': ['Robv', 'Rrev'],
    '6': ['P'],
    '-6': ['P'],
    '6/m': ['P'],
    '622': ['P'],
    '6mm': ['P'],
    '-62m': ['P'],
    '-6m2': ['P'],
    '6/mmm': ['P'],
    '23': ['P', 'I', 'F'],
    'm-3': ['P', 'I', 'F'],
    '432': ['P', 'I', 'F'],
    '-43m': ['P', 'I', 'F'],
    'm-3m': ['P', 'I', 'F'],
}

crystal_system_point_groups = {
    'Triclinic': ['1', '-1'],
    'Monoclinic': ['2', 'm', '2/m', '112', '11m', '112/m'],
    'Orthorhombic': ['222', 'mm2', 'mmm'],
    'Tetragonal': ['4', '-4', '4/m', '422', '4mm', '-42m', '-4m2', '4/mmm'],
    'Trigonal/Rhombohedral': ['3 r', '-3 r', '32 r', '3m r', '-3m r'],
    'Trigonal/Hexagonal': ['3', '-3', '312', '31m', '32', '321',
                           '3m', '-31m', '-3m', '-3m1'],
    'Hexagonal': ['6', '-6', '6/m', '622', '6mm', '-62m', '-6m2', '6/mmm'],
    'Cubic': ['23', 'm-3', '432', '-43m', 'm-3m']
}

class ExperimentModel(NeuXtalVizModel):

    def __init__(self):

        super(ExperimentModel, self).__init__()

        CreateEmptyTableWorkspace(OutputWorkspace='plan')

        CreatePeaksWorkspace(NumberOfPeaks=0,
                             OutputType='LeanElasticPeak',
                             OutputWorkspace='coverage')

    def initialize_instrument(self, instrument, logs):

        inst = self.get_instrument_name(instrument)

        if not mtd.doesExist('instrument'):

            LoadEmptyInstrument(InstrumentName=inst,
                                OutputWorkspace='instrument')

            for key in logs.keys():

                AddSampleLog(Workspace='instrument',
                             LogName=key,
                             LogText=str(logs[key]),
                             LogType='Number Series',
                             NumberType='Double')

            LoadInstrument(Workspace='instrument',
                           RewriteSpectraMap=False,
                           InstrumentName=inst)

            ExtractMonitors(InputWorkspace='instrument',
                            MonitorWorkspace='monitors',
                            DetectorWorkspace='instrument')

            PreprocessDetectorsToMD(InputWorkspace='instrument',
                                    OutputWorkspace='detectors')

            cols, rows = beamlines[instrument]['BankPixels']
            grouping = beamlines[instrument]['Grouping']

            c, r = [int(val) for val in grouping.split('x')]
            shape = (-1, cols, rows)

            det_map = np.array(mtd['detectors'].column(5)).reshape(*shape)

            shape = det_map.shape
            i, j, k = np.meshgrid(np.arange(shape[0]),
                                  np.arange(shape[1]),
                                  np.arange(shape[2]), indexing='ij')
            keys = np.stack((i, j // c, k // r), axis=-1)
            keys_flat = keys.reshape(-1, keys.shape[-1])
            det_map_flat = det_map.ravel().astype(str)
            grouped_ids = defaultdict(list)

            for key, detector_id in zip(map(tuple, keys_flat), det_map_flat):
                grouped_ids[key].append(detector_id)

            detector_list = ','.join('+'.join(group)\
                                     for group in grouped_ids.values())

            GroupDetectors(InputWorkspace='instrument',
                           OutputWorkspace='instrument',
                           GroupingPattern=detector_list)

            CreatePeaksWorkspace(InstrumentWorkspace='instrument',
                                 NumberOfPeaks=0,
                                 OutputType='LeanElasticPeak',
                                 OutputWorkspace='peak')

            CreatePeaksWorkspace(InstrumentWorkspace='instrument',
                                 NumberOfPeaks=0,
                                 OutputType='Peak',
                                 OutputWorkspace='peaks')

            CreatePeaksWorkspace(InstrumentWorkspace='instrument',
                                 NumberOfPeaks=0,
                                 OutputType='Peak',
                                 OutputWorkspace='combined')

            L2 = np.array(mtd['detectors'].column(1))
            tt = np.array(mtd['detectors'].column(2))
            az = np.array(mtd['detectors'].column(3))

            x = L2*np.sin(tt)*np.cos(az)
            y = L2*np.sin(tt)*np.sin(az)
            z = L2*np.cos(tt)

            self.nu = np.rad2deg(np.arcsin(y/L2))
            self.gamma = np.rad2deg(np.arctan2(x,z))

    def remove_instrument(self):

        if mtd.doesExist('instrument'):

            DeleteWorkspace(Workspace='instrument')

        if mtd.doesExist('cobmined'):

            DeleteWorkspace(Workspace='cobmined')

        if mtd.doesExist('filtered'):

            DeleteWorkspace(Workspace='filtered')

    def get_crystal_system_point_groups(self, crystal_system):

        return crystal_system_point_groups[crystal_system]

    def get_point_group_centering(self, point_group):

        return point_group_centering[point_group]

    def get_symmetry(self, point_group, centering):

        pg = point_group_dict[point_group]
        lc = lattice_centering_dict[centering]        

        return str(pg), str(lc)

    def create_plan(self, instrument, mode, angles):

        mtd['plan'].setRowCount(0)
        for col in mtd['plan'].getColumnNames():
            mtd['plan'].removeColumn(col)

        mtd['plan'].setTitle('{} {}'.format(instrument, mode))
        mtd['plan'].setComment('{} {}'.format(instrument, mode))
        mtd['plan'].addColumn('str', 'Title')

        for angle in angles:
            mtd['plan'].addColumn('float', angle)

        mtd['plan'].addColumn('str', 'comment')
        mtd['plan'].addColumn('bool', 'use')

    def load_UB(self, filename):

        LoadIsawUB(InputWorkspace='coverage', Filename=filename)

        self.copy_UB()

    def copy_UB(self):

        if self.has_UB():

            UB = mtd['coverage'].sample().getOrientedLattice().getUB().copy()

            self.set_UB(UB)

    def has_UB(self):

        if HasUB(Workspace='coverage'):
            return True
        else:
            return False

    def get_instrument_name(self, instrument):

        return beamlines[instrument]['Name']

    def get_modes(self, instrument):

        return list(beamlines[instrument]['Goniometer'].keys())

    def get_axes_polarities(self, instrument, mode):

        goniometers = beamlines[instrument]['Goniometer'][mode]

        axes = [goniometers[name][:-3] for name in goniometers.keys()]

        polarities = [goniometers[name][3] for name in goniometers.keys()]

        return axes, polarities

    def get_goniometers(self, instrument, mode):

        goniometers = beamlines[instrument]['Goniometer'][mode]

        return [(name, *goniometers[name][-2:]) for name in goniometers.keys()]

    def get_motors(self, instrument):

        motors = beamlines[instrument].get('Motor')

        if motors is not None:
            return [(name, motors[name]) for name in motors.keys()]
        else:
            return []

    def get_wavelength(self, instrument):

        return beamlines[instrument]['Wavelength']

    def save_plan(self, filename):

        plan_dict = mtd['plan'].toDict().copy()
        use_angle = plan_dict('Use')

        for key in plan_dict.keys():
            items = plan_dict[key]
            items = [item for item, use in zip(items, use_angle) if use]
            plan_dict[key] = items

        with open(filename, mode='w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=plan_dict.keys())
            writer.writeheader()
            for row in zip(*plan_dict.values()):
                writer.writerow(dict(zip(plan_dict.keys(), row)))

    def _calculate_matrices(self, axes, polarities, limits, step):

        self.axes = [None]*6

        for i, (axis, polarity) in enumerate(zip(axes, polarities)):
            self.axes[i] = '{},'
            self.axes[i] +=','.join(np.array([*axis,polarity]).astype(str))

        angular_coverage = []
        for limit in limits:
            angular_coverage.append(np.arange(limit[0], limit[1]+step, step))

        axes = np.array(axes)
        polarities = np.array(polarities)

        angle_settings = np.meshgrid(*angular_coverage, indexing='ij')
        angle_settings = np.reshape(angle_settings, (len(polarities), -1)).T

        self.angles = angle_settings.copy()

        angle_settings = angle_settings*polarities
        angle_settings = np.deg2rad(angle_settings)

        rotation_vectors = angle_settings[...,None]*axes
        rotation_vectors = rotation_vectors.reshape(-1, 3)

        all_rotations = Rotation.from_rotvec(rotation_vectors).as_matrix()
        all_rotations = all_rotations.reshape(*angle_settings.shape, 3, 3)

        Rs = []
        for i in range(all_rotations.shape[0]):
            R = np.eye(3)
            for j in range(all_rotations.shape[1]):
                R = R @ all_rotations[i,j,:,:]
            Rs.append(R)

        return Rs

    def individual_peak(self, hkl,
                              wavelength,
                              axes,
                              polarities,
                              limits,
                              step=1):

        self.comment = '('+' '.join(np.array(hkl).astype(str))+')'

        if np.isclose(wavelength[0], wavelength[1]):
            wavelength = [0.975*wavelength[0], 1.025*wavelength[1]]

        FilterPeaks(InputWorkspace='peak',
                    OutputWorkspace='peak',
                    FilterVariable='RunNumber',
                    FilterValue=-1,
                    Operator='=')

        UB = mtd['coverage'].sample().getOrientedLattice().getUB().copy()

        SetUB(Workspace='peak', UB=UB)

        AddPeakHKL(Workspace='peak', HKL=hkl)

        Q_sample = mtd['peak'].getPeak(0).getQSampleFrame()

        Q = np.sqrt(np.dot(Q_sample, Q_sample))

        Rs = self._calculate_matrices(axes, polarities, limits, step)

        mtd['peak'].run().getGoniometer().setR(np.eye(3))
        mtd['peaks'].run().getGoniometer().setR(np.eye(3))

        Q_lab = np.einsum('kij,j->ki', Rs, Q_sample)

        lamda = -4*np.pi*Q_lab[:,2]/Q**2
        mask = (lamda > wavelength[0]) & (lamda < wavelength[1])

        k = 2*np.pi/lamda

        ki = k[:,np.newaxis]*np.array([0,0,1])
        kf = Q_lab+ki

        gamma = np.rad2deg(np.arctan2(kf[:,0],kf[:,2]))[mask]
        nu = np.rad2deg(np.arcsin(kf[:,1]/k))[mask]
        lamda = lamda[mask]

        self.angles = self.angles[mask]
        self.angles_gamma = gamma.copy()
        self.angles_nu = nu.copy()

        if len(lamda) > 0:

            k = 2*np.pi/lamda
            Qx = k*np.cos(np.deg2rad(nu))*np.sin(np.deg2rad(gamma))
            Qy = k*np.sin(np.deg2rad(nu))
            Qz = k*(np.cos(np.deg2rad(nu))*np.cos(np.deg2rad(gamma))-1)

            mask = []
            for i in range(len(k)):
                peak = mtd['peaks'].createPeak(V3D(Qx[i],Qy[i],Qz[i]))
                mask.append(peak.getDetectorID() > 0)

            mask = np.array(mask)

            gamma = gamma[mask]
            nu = nu[mask]
            lamda = lamda[mask]

            self.angles = self.angles[mask]
            self.angles_gamma = gamma.copy()
            self.angles_nu = nu.copy()

        return gamma, nu, lamda

    def simultaneous_peaks(self, hkl_1,
                                 hkl_2,
                                 wavelength,
                                 axes,
                                 polarities,
                                 limits,
                                 step=1):

        self.comment = '('+' '.join(np.array(hkl_1).astype(str))+')'

        if np.isclose(wavelength[0], wavelength[1]):
            wavelength = [0.975*wavelength[0], 1.025*wavelength[1]]

        FilterPeaks(InputWorkspace='peak',
                    OutputWorkspace='peak',
                    FilterVariable='RunNumber',
                    FilterValue=-1,
                    Operator='=')

        UB = mtd['coverage'].sample().getOrientedLattice().getUB().copy()

        SetUB(Workspace='peak', UB=UB)

        AddPeakHKL(Workspace='peak', HKL=hkl_1)
        AddPeakHKL(Workspace='peak', HKL=hkl_2)

        Q0_sample = mtd['peak'].getPeak(0).getQSampleFrame()
        Q1_sample = mtd['peak'].getPeak(1).getQSampleFrame()

        Q0 = np.sqrt(np.dot(Q0_sample, Q0_sample))
        Q1 = np.sqrt(np.dot(Q1_sample, Q1_sample))

        Rs = self._calculate_matrices(axes, polarities, limits, step)

        Q0_lab = np.einsum('kij,j->ki', Rs, Q0_sample)
        Q1_lab = np.einsum('kij,j->ki', Rs, Q1_sample)

        lamda0 = -4*np.pi*Q0_lab[:,2]/Q0**2
        lamda1 = -4*np.pi*Q1_lab[:,2]/Q1**2

        mask = (lamda0 > wavelength[0]) & (lamda0 < wavelength[1])\
             & (lamda1 > wavelength[0]) & (lamda1 < wavelength[1])

        k0 = 2*np.pi/lamda0
        k1 = 2*np.pi/lamda1

        k0i = k0[:,np.newaxis]*np.array([0,0,1])
        k1i = k1[:,np.newaxis]*np.array([0,0,1])

        k0f = Q0_lab+k0i
        k1f = Q1_lab+k1i

        gamma0 = np.rad2deg(np.arctan2(k0f[:,0],k0f[:,2]))[mask]
        gamma1 = np.rad2deg(np.arctan2(k1f[:,0],k1f[:,2]))[mask]

        nu0 = np.rad2deg(np.arcsin(k0f[:,1]/k0))[mask]
        nu1 = np.rad2deg(np.arcsin(k1f[:,1]/k1))[mask]

        lamda0 = lamda0[mask]
        lamda1 = lamda1[mask]

        self.angles = self.angles[mask]
        self.angles_gamma = gamma0.copy()
        self.angles_nu = nu0.copy()

        if len(lamda0) > 0:

            k0 = 2*np.pi/lamda0
            k1 = 2*np.pi/lamda1

            Q0x = k0*np.cos(np.deg2rad(nu0))*np.sin(np.deg2rad(gamma0))
            Q1x = k1*np.cos(np.deg2rad(nu1))*np.sin(np.deg2rad(gamma1))

            Q0y = k0*np.sin(np.deg2rad(nu0))
            Q1y = k1*np.sin(np.deg2rad(nu1))

            Q0z = k0*(np.cos(np.deg2rad(nu0))*np.cos(np.deg2rad(gamma0))-1)
            Q1z = k1*(np.cos(np.deg2rad(nu1))*np.cos(np.deg2rad(gamma1))-1)

            mask = []
            for i in range(len(k1)):
                peak0 = mtd['peaks'].createPeak(V3D(Q0x[i],Q0y[i],Q0z[i]))
                peak1 = mtd['peaks'].createPeak(V3D(Q1x[i],Q1y[i],Q1z[i]))
                mask.append((peak0.getDetectorID() > 0) &\
                            (peak1.getDetectorID() > 0))

            mask = np.array(mask)

            gamma0 = gamma0[mask]
            gamma1 = gamma1[mask]

            nu0 = nu0[mask]
            nu1 = nu1[mask]

            lamda0 = lamda0[mask]
            lamda1 = lamda1[mask]

            self.angles = self.angles[mask]
            self.angles_gamma = gamma0.copy()
            self.angles_nu = nu0.copy()

        return (gamma0, nu0, lamda0), (gamma1, nu1, lamda1)

    def get_angles(self, gamma, nu):

        d2 = (self.angles_gamma-gamma)**2+(self.angles_nu-nu)**2

        i = np.argmin(d2)

        angles = self.angles[i]

        return angles

    def add_orientation(self, angles, wavelength, d_min, rows):

        if np.isclose(wavelength[0], wavelength[1]):
            wavelength = [0.975*wavelength[0], 1.025*wavelength[1]]

        axes = self.axes.copy()

        for i, angle in enumerate(angles):
            axes[i] = axes[i].format(angle)

        ol = mtd['coverage'].sample().getOrientedLattice()
        UB = ol.getUB().copy()
        
        SetUB(Workspace='instrument', UB=UB)

        SetGoniometer(Workspace='instrument',
                      Axis0=axes[0],
                      Axis1=axes[1],
                      Axis2=axes[2],
                      Axis3=axes[3],
                      Axis4=axes[4],
                      Axis5=axes[5])

        d_max = 1.1*np.max([ol.d(1,0,0),ol.d(0,1,0),ol.d(0,0,1)])

        ws = 'peaks_orientation_{}'.format(rows)

        PredictPeaks(InputWorkspace='instrument',
                     MinDSpacing=d_min,
                     MaxDSpacing=d_max,
                     WavelengthMin=wavelength[0],
                     WavelengthMax=wavelength[1],
                     ReflectionCondition='Primitive',
                     OutputWorkspace=ws)

        SortPeaksWorkspace(InputWorkspace=ws,
                           ColumnNameToSortBy='DSpacing',
                           SortAscending=False,
                           OutputWorkspace=ws)

        columns = ['l', 'k', 'h']

        for col in columns:

            SortPeaksWorkspace(InputWorkspace=ws,
                               ColumnNameToSortBy=col,
                               SortAscending=False,
                               OutputWorkspace=ws)

        for no in range(mtd[ws].getNumberPeaks()-1,0,-1):

            if (mtd[ws].getPeak(no).getHKL()-\
                mtd[ws].getPeak(no-1).getHKL()).norm2() == 0:

                DeleteTableRows(TableWorkspace=ws, Rows=no)

        SortPeaksWorkspace(InputWorkspace=ws,
                           ColumnNameToSortBy='DSpacing',
                           SortAscending=False,
                           OutputWorkspace=ws)

        for peak in mtd[ws]:
            peak.setRunNumber(rows)
            peak.setIntensity(10)
            peak.setSigmaIntensity(np.sqrt(peak.getIntensity()))

        mtd['instrument'].run().getGoniometer().setR(np.eye(3))

        SetUB(Workspace='combined', UB=UB)
        CombinePeaksWorkspaces(LHSWorkspace='combined',
                               RHSWorkspace=ws,
                               OutputWorkspace='combined')

    def calculate_statistics(self, point_group, lattice_centering, use):

        if mtd.doesExist('combined'):

            CloneWorkspace(InputWorkspace='combined',
                           OutputWorkspace='filtered')

            rows = np.arange(len(use)).tolist()

            for row in rows:
                if not use[row]:
                    FilterPeaks(InputWorkspace='filtered',
                                FilterVariable='RunNumber',
                                FilterValue=str(row),
                                Operator='!=',
                                OutputWorkspace='filtered')

            if mtd['filtered'].getNumberPeaks() > 0:

                pg, lc = self.get_symmetry(point_group, lattice_centering)

                # mantid bug
                # if lattice_centering == 'F':
                #     lc = 'F' 

                StatisticsOfPeaksWorkspace(InputWorkspace='filtered',
                                           OutputWorkspace='filtered',
                                           StatisticsTable='statistics',
                                           EquivalentsWorkspace='equivalents',
                                           PointGroup=pg,
                                           LatticeCentering=lc)

                # CloneWorkspace(InputWorkspace='tmp',
                #                OutputWorkspace='filtered')

                stats_dict = mtd['statistics'].toDict()

                # d_min = stats_dict['Resolution Min']
                # d_max = stats_dict['Resolution Max']

                shel = stats_dict['Resolution Shell']
                mult = stats_dict['Multiplicity']
                refl = stats_dict['No. of Unique Reflections']
                comp = stats_dict['Data Completeness']

                DeleteWorkspace(Workspace='statistics')
                DeleteWorkspace(Workspace='equivalents')

                return shel, comp, mult, refl

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

    def get_coverage_info(self, point_group):

        pg = PointGroupFactory.createPointGroup(point_group)

        coverage_dict = {}

        UB = mtd['coverage'].sample().getOrientedLattice().getUB().copy()
        # UB_inv = np.linalg.inv(UB)

        h = mtd['filtered'].column('h')
        k = mtd['filtered'].column('k')
        l = mtd['filtered'].column('l')

        hkls = np.array([h,k,l]).T.astype(int).tolist()

        hkl_dict = {}
        for hkl in hkls:
            equiv_hkls = pg.getEquivalents(hkl)
            for equiv_hkl in equiv_hkls:
                key = tuple(equiv_hkl)
                no = hkl_dict.get(key)
                if no is None:
                    no = 1
                else:
                    no += 1
                hkl_dict[key] = no

        nos = np.array([value for value in hkl_dict.values()])
        hkls = np.array([key for key in hkl_dict.keys()])

        r = np.sqrt(hkls[:,0]**2+hkls[:,1]**2+hkls[:,2]**2)
        theta = np.arccos(hkls[:,2]/r)
        phi = np.arctan2(hkls[:,1], hkls[:,0])

        hue = phi*180/np.pi+180
        saturation = np.ones_like(hue)
        lightness = theta/np.pi

        rgb = self.hsl_to_rgb(hue, saturation, lightness)
        coords = np.einsum('ij,nj->ni', 2*np.pi*UB, hkls)

        coverage_dict['colors'] = (rgb*255).astype(np.uint8)
        coverage_dict['sizes'] = nos/nos.max()
        coverage_dict['coords'] = coords

        return coverage_dict

class CoverageOptimizer(ParallelTasks):

    def __init__(self):

        super().__init__(None, None)

    def initialize_parameters(self, inst_name, axes, limits, UB,
                                    wl_limits, point_group, refl_cond, d_min):

        self.inst_name = inst_name
        self.axes = axes
        self.limits = limits
        self.UB = UB
        self.wl_limits = wl_limits
        self.point_group = point_group
        self.refl_cond = refl_cond
        self.d_min = d_min
        self.d_max = np.sqrt(np.diag(np.linalg.inv(UB.T @ UB))).min()*1.2

        self.best = []
        self.worst = []

        self.settings = None

    def get_coverage(self):

        return self.best, self.worst

    def get_settings(self):

        return self.settings

    def initialize_settings(self, n_orient, n_indiv, outname,
                                  outdir='/tmp/', n_proc=2):

        self.n_orient = n_orient
        self.n_indiv = n_indiv

        fname = os.path.join(outdir, outname+'_ind{}.nxs')
        fnames = [fname.format(i_indiv) for i_indiv in range(n_indiv)]

        plan = outname+'_ind{}'
        plans = [plan.format(i_indiv) for i_indiv in range(n_indiv)]

        self.function = self._generation
        self.args = (n_orient, self.inst_name, self.axes, self.limits,
                     self.UB, self.wl_limits, self.refl_cond,
                     self.d_min, self.d_max)

        self.run_tasks(fnames, n_proc)

        fitness = []
        for fname, plan in zip(fnames, plans):
            LoadNexus(Filename=fname, OutputWorkspace=plan)
            fitness.append(self._coverage(plan))

        self.fitness = np.array(fitness)
        self.fnames = np.array(fnames)
        self.plans = np.array(plans)
        self.plan = outname+'_ind{}'

        ranking = np.argsort(self.fitness)

        cov_min = self.fitness[ranking[0]]
        cov_max = self.fitness[ranking[-1]]

        self.best.append(cov_max)
        self.worst.append(cov_min)

    def _coverage(self, peaks_ws):

        StatisticsOfPeaksWorkspace(InputWorkspace=peaks_ws,
                                   OutputWorkspace=peaks_ws+'_sorted',
                                   StatisticsTable=peaks_ws+'_stats',
                                   EquivalentsWorkspace=peaks_ws+'_equiv',
                                   PointGroup=self.point_group,
                                   LatticeCentering=self.refl_cond)

        stats = mtd[peaks_ws+'_stats'].row(0)

        return stats['Data Completeness']

    def _load_instrument(self, inst_name):

        if not mtd.doesExist(inst_name):

            LoadEmptyInstrument(InstrumentName=inst_name,
                                OutputWorkspace=inst_name)

    def _set_UB(self, ws, UB):

        SetUB(Workspace=ws, UB=UB)

    def _generation(self, fnames, n_orient, inst_name, axes, limits, UB,
                          wl_limits, refl_cond, d_min, d_max, proc=1):

        if not mtd.doesExist(inst_name):

            self._load_instrument(inst_name)
            self._set_UB(inst_name, UB)

        for fname in fnames:

            CreatePeaksWorkspace(InstrumentWorkspace=inst_name,
                                 NumberOfPeaks=0,
                                 OutputWorkspace='gen_peaks')

            self._set_UB('gen_peaks', UB)

            for i_orient in range(n_orient):

                self._generate_setting(inst_name, i_orient, axes, limits, UB,
                                       wl_limits, refl_cond, d_min, d_max)

                CombinePeaksWorkspaces(LHSWorkspace='gen_peaks',
                                       RHSWorkspace='gen_peaks_ws',
                                       OutputWorkspace='gen_peaks')

                run = mtd['gen_peaks_ws'].run()

                for key in run.keys():

                    angle = run[key].value[0]

                    AddSampleLog(Workspace='gen_peaks',
                                 LogName='{}_{}'.format(key,i_orient),
                                 LogText='{}'.format(angle))

                DeleteWorkspace(Workspace='gen_peaks_ws')

            SaveNexus(InputWorkspace='gen_peaks', Filename=fname)

    def _generate_setting(self, inst_name, run, axes, limits, UB,
                                wl_limits, refl_cond, d_min, d_max):

        if not mtd.doesExist(inst_name):

            self._load_instrument(inst_name)
            self._set_UB(inst_name, UB)

        ax = []
        for axis, limit in zip(axes, limits):
            if limit is not None:
                angle = limit[0]+(limit[1]-limit[0])*np.random.random()
                ax.append(axis.format(round(angle,1)))
            else:
                ax.append(axis)
        for _ in range(6-len(axes)):
            ax.append(None)

        SetGoniometer(Workspace=inst_name,
                      Axis0=ax[0],
                      Axis1=ax[1],
                      Axis2=ax[2],
                      Axis3=ax[3],
                      Axis4=ax[4],
                      Axis5=ax[5])

        PredictPeaks(InputWorkspace=inst_name,
                     MinDSpacing=d_min,
                     MaxDSpacing=d_max,
                     WavelengthMin=wl_limits[0],
                     WavelengthMax=wl_limits[1],
                     ReflectionCondition=refl_cond,
                     OutputWorkspace='gen_peaks_ws')

        for peak in mtd['gen_peaks_ws']:
            peak.setRunNumber(run)
            peak.setIntensity(1)
            peak.setSigmaIntensity(peak.getIntensity())

    def optimize_settings(self, n_gener, elite_rate=20, mutation_rate=20):

        n_elites = int(round(elite_rate/100*self.n_indiv))

        ranking = np.argsort(self.fitness)

        for _ in range(n_gener):

            elites = self.plans[ranking[-n_elites:]]

            fraction = self.fitness/np.sum(self.fitness)

            selections = []
            while len(selections) < (self.n_indiv-n_elites+1) // 2:
                choices = np.random.choice(self.plans, size=2,
                                           p=fraction, replace=False)
                selections.append(choices)

            self._crossover(elites, selections)

            self._mutation(mutation_rate)

            self.fitness = []
            for plan in self.plans:
                self.fitness.append(self._coverage(plan))

            ranking = np.argsort(self.fitness)

            cov_min = self.fitness[ranking[0]]
            cov_max = self.fitness[ranking[-1]]

            self.best.append(cov_max)
            self.worst.append(cov_min)

    def _mutation(self, mutation_rate):

        for i_indiv, plan in enumerate(self.plans):
            for i_orient in range(self.n_orient):
                if 100*np.random.random() < mutation_rate:

                    FilterPeaks(InputWorkspace=plan,
                                FilterVariable='RunNumber',
                                FilterValue=i_orient,
                                Operator='!=',
                                OutputWorkspace=plan)

                    self._generate_setting(self.inst_name, i_orient,
                                           self.axes, self.limits, self.UB,
                                           self.wl_limits, self.refl_cond,
                                           self.d_min, self.d_max)

                    CombinePeaksWorkspaces(LHSWorkspace=plan,
                                           RHSWorkspace='gen_peaks_ws',
                                           OutputWorkspace=plan)


                    run = mtd['gen_peaks_ws'].run()

                    for key in run.keys():

                        angle = run[key].value[0]

                        AddSampleLog(Workspace=plan,
                                     LogName='{}_{}'.format(key,i_orient),
                                     LogText='{}'.format(angle))

                    DeleteWorkspace(Workspace='gen_peaks_ws')

    def _crossover(self, elites, selections):

        i_indiv = 0

        next_plan = 'next_'+self.plan
        next_plans = []

        for best in elites:

            plan = next_plan.format(i_indiv)

            CloneWorkspace(InputWorkspace=best,
                           OutputWorkspace=plan)

            next_plans.append(plan)

            i_indiv += 1

        for recombinations in selections:

            k = np.random.randint(self.n_orient)

            for parents in [recombinations, recombinations[::-1]]:

                if i_indiv < self.n_indiv:

                    plan = next_plan.format(i_indiv)

                    FilterPeaks(InputWorkspace=parents[0],
                                FilterVariable='RunNumber',
                                FilterValue=k,
                                Operator='<',
                                OutputWorkspace='cross_peaks_ws0')

                    FilterPeaks(InputWorkspace=parents[1],
                                FilterVariable='RunNumber',
                                FilterValue=k,
                                Operator='>=',
                                OutputWorkspace='cross_peaks_ws1')

                    CombinePeaksWorkspaces(LHSWorkspace='cross_peaks_ws0',
                                           RHSWorkspace='cross_peaks_ws1',
                                           OutputWorkspace=plan)

                    DeleteWorkspace(Workspace='cross_peaks_ws0')
                    DeleteWorkspace(Workspace='cross_peaks_ws1')

                    next_plans.append(plan)

                    i_indiv += 1

        RenameWorkspaces(InputWorkspaces=next_plans,
                         WorkspaceNames=self.plans)