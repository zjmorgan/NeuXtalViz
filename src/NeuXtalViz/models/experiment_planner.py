import os

import itertools
import csv

from mantid.simpleapi import (CreatePeaksWorkspace,
                              PredictPeaks,
                              FilterPeaks,
                              StatisticsOfPeaksWorkspace,
                              CombinePeaksWorkspaces,
                              SetUB,
                              SetGoniometer,
                              LoadNexus,
                              SaveNexus,
                              LoadIsawUB,
                              LoadEmptyInstrument,
                              LoadInstrument,
                              ExtractMonitors,
                              PreprocessDetectorsToMD,
                              AddSampleLog,
                              CreateEmptyTableWorkspace,
                              CloneWorkspace,
                              DeleteWorkspace,
                              RenameWorkspaces,
                              mtd)

# from mantid.geometry import PointGroupFactory

import numpy as np

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
    'R(obv)': 'Rhombohedrally centred, obverse',
    'R(rev)': 'Rhombohedrally centred, reverse',
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
    '3r': '3 r (Trigonal - Rhombohedral)',
    '-3r': '-3 r (Trigonal - Rhombohedral)',
    '32r': '32 r (Trigonal - Rhombohedral)',
    '3mr': '3m r (Trigonal - Rhombohedral)',
    '-3mr': '-3m r (Trigonal - Rhombohedral)',
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
    '3r': ['R'],
    '-3r': ['R'],
    '32r': ['R'],
    '3mr': ['R'],
    '-3mr': ['R'],
    '3': ['R(obv)', 'R(rev)'],
    '-3': ['R(obv)', 'R(rev)'],
    '312': ['R(obv)', 'R(rev)'],
    '31m': ['R(obv)', 'R(rev)'],
    '32': ['R(obv)', 'R(rev)'],
    '321': ['R(obv)', 'R(rev)'],
    '3m': ['R(obv)', 'R(rev)'],
    '-31m': ['R(obv)', 'R(rev)'],
    '-3m': ['R(obv)', 'R(rev)'],
    '-3m1': ['R(obv)', 'R(rev)'],
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
    'Trigonal/Rhombohedral': ['3r', '-3r', '32r', '3mr', '-3mr'],
    'Trigonal/Hexagonal': ['3', '-3', '312', '31m', '32', '321',
                           '3m', '-31m', '-3m', '-3m1'],
    'Hexagonal': ['6', '-6', '6/m', '622', '6mm', '-62m', '-6m2', '6/mmm'],
    'Cubic': ['23', 'm-3', '432', '-43m', 'm-3m']
}

class ExperimentModel(NeuXtalVizModel):

    def __init__(self):

        super(ExperimentModel, self).__init__()

    def get_crystal_system_point_groups(self, crystal_system):

        return crystal_system_point_groups[crystal_system]

    def get_point_group_centering(self, point_group):

        return point_group_centering[point_group]

    def get_symmetry(self, point_group, centering):

        return point_group_dict[point_group], lattice_centering_dict[centering]

    def create_plan(self, instrument, mode, angles):

        CreateEmptyTableWorkspace(OutputWorkspace='plan')

        mtd['plan'].setTitle('{} {}'.format(instrument, mode))
        mtd['plan'].setComment('{} {}'.format(instrument, mode))
        mtd['plan'].addColumn('str', 'Title')
        for angle in angles:
            mtd['plan'].addColumn('float', angle)
        mtd['plan'].addColumn('str', 'comment')
        mtd['plan'].addColumn('bool', 'use')

    def load_UB(self, filename):

        CreatePeaksWorkspace(OutputType='LeanElasticPeak',
                             NumberOfPeaks=0,
                             OutputWorkspace='coverage')        

        LoadIsawUB(InputWorkspace='coverage', Filename=filename)

        self.copy_UB()

    def copy_UB(self):

        if self.has_UB('coverage'):

            UB = mtd['coverage'].sample().getOrientedLattice().getUB().copy()

            self.set_UB(UB)

    def get_instrument_name(self, instrument):

        return beamlines[instrument]['Name']

    def get_modes(self, instrument):

        return list(beamlines[instrument]['Goniometer'].keys())

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