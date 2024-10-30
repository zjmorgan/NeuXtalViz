import os

from mantid.simpleapi import (SelectCellWithForm,
                              ShowPossibleCells,
                              TransformHKL,
                              CalculatePeaksHKL,
                              IndexPeaks,
                              FindUBUsingFFT,
                              FindUBUsingLatticeParameters,
                              FindUBUsingIndexedPeaks,
                              OptimizeLatticeForCellType,
                              CalculateUMatrix,
                              HasUB,
                              SetUB,
                              LoadIsawUB,
                              SaveIsawUB,
                              FindPeaksMD,
                              PredictPeaks,
                              PredictSatellitePeaks,
                              CentroidPeaksMD,
                              IntegratePeaksMD,
                              PeakIntensityVsRadius,
                              FilterPeaks,
                              SortPeaksWorkspace,
                              DeleteWorkspace,
                              DeleteTableRows,
                              ExtractSingleSpectrum,
                              CombinePeaksWorkspaces,
                              CreatePeaksWorkspace,
                              ConvertPeaksWorkspace,
                              ConvertQtoHKLMDHisto,
                              CompactMD,
                              CopySample,
                              CreateSampleWorkspace,
                              CloneWorkspace,
                              SaveNexus,
                              LoadNexus,
                              HB3AAdjustSampleNorm,
                              LoadWANDSCD,
                              Load,
                              Rebin,
                              SetGoniometer,
                              PreprocessDetectorsToMD,
                              ConvertToMD,
                              ConvertHFIRSCDtoMDE,
                              LoadMD,
                              SaveMD,
                              MergeMD,
                              LoadNexus,
                              LoadIsawDetCal,
                              LoadParameterFile,
                              ApplyCalibration,
                              BinMD,
                              ConvertUnits,
                              CropWorkspace,
                              mtd)

from mantid import config

config['Q.convention'] = 'Crystallography'

from mantid.geometry import PointGroupFactory, UnitCell
from mantid.kernel import V3D

import numpy as np
import scipy
import json

from NeuXtalViz.models.base_model import NeuXtalVizModel
from NeuXtalViz.config.instruments import beamlines

lattice_group = {'Triclinic': '-1',
                 'Monoclinic': '2/m',
                 'Orthorhombic': 'mmm',
                 'Tetragonal': '4/mmm',
                 'Rhombohedral': '-3m',
                 'Hexagonal': '6/mmm',
                 'Cubic': 'm-3m'}

centering_reflection = {'P': 'Primitive',
                        'I': 'Body centred',
                        'F': 'All-face centred',
                        'R': 'Primitive', # rhomb axes
                        'R(obv)': 'Rhombohderally centred, obverse', # hex axes
                        'R(rev)': 'Rhombohderally centred, reverse', # hex axes
                        'A': 'A-face centred',
                        'B': 'B-face centred',
                        'C': 'C-face centred'}

variable = {'I/σ': 'Signal/Noise',
            'I': 'Intensity',
            'd': 'DSpacing',
            'λ': 'Wavelength',
            'Q': 'QMod',
            'h^2+k^2+l^2': 'h^2+k^2+l^2',
            'm^2+n^2+p^2': 'm^2+n^2+p^2',
            'Run #': 'RunNumber'}

class UBModel(NeuXtalVizModel):

    def __init__(self):

        super(UBModel, self).__init__()

        self.Q = None
        self.table = 'ub_peaks'

        self.peak_info = None

        CreateSampleWorkspace(OutputWorkspace='ub_lattice')

    def has_Q(self):

        if self.Q is None:
            return False
        elif mtd.doesExist(self.Q):
            return True
        else:
            return False

    def has_peaks(self):

        if mtd.doesExist(self.table):
            return True
        else:
            return False

    def has_UB(self):

        if self.has_peaks():
            if HasUB(self.table):
                return True
            else:
                return False
        else:
            return False

    def get_UB(self):

        if self.has_UB():

            return mtd[self.table].sample().getOrientedLattice().getUB().copy()

    def update_UB(self):

        UB = self.get_UB()

        if UB is not None:

            self.set_UB(UB)

    def get_instrument_name(self, instrument):

        return beamlines[instrument]['Name']

    def get_goniometers(self, instrument):

        return beamlines[instrument]['Goniometers']

    def get_wavelength(self, instrument):

        return beamlines[instrument]['Wavelength']

    def get_raw_file_path(self, instrument):

        inst = beamlines[instrument]

        return os.path.join('/',
                            inst['Facility'],
                            inst['InstrumentName'],
                            'IPTS-{}',
                            inst['RawFile'])

    def get_shared_file_path(self, instrument, ipts):

        inst = beamlines[instrument]

        if ipts is not None:
            filepath = os.path.join('/',
                                    inst['Facility'],
                                    inst['InstrumentName'],
                                    'IPTS-{}'.format(ipts),
                                    'shared')
            if os.path.exists(filepath):
                return filepath

        filepath = os.path.join('/', inst['Facility'], inst['InstrumentName'])

        return filepath

    def get_calibration_file_path(self, instrument):

        inst = beamlines[instrument]

        return os.path.join('/',
                            inst['Facility'],
                            inst['InstrumentName'],
                            'shared',
                            'calibration')

    def get_vanadium_file_path(self, instrument):

        inst = beamlines[instrument]

        return os.path.join('/',
                            inst['Facility'],
                            inst['InstrumentName'],
                            'shared',
                            'Vanadium')

    def load_data(self,
                  instrument,
                  IPTS,
                  runs,
                  exp,
                  time_stop):

        filepath = self.get_raw_file_path(instrument)

        if instrument == 'DEMAND':
            filenames = [filepath.format(IPTS, exp, runs)]
            if np.all([os.path.exists(filename) for filename in filenames]):
                HB3AAdjustSampleNorm(Filename=filenames,
                                     OutputType='Detector',
                                     NormaliseBy='None',
                                     Grouping='4x4',
                                     OutputWorkspace='data')

        elif instrument == 'WAND²':
            filenames = [filepath.format(IPTS, run) for run in runs]
            if np.all([os.path.exists(filename) for filename in filenames]):
                filenames = ','.join([filename for filename in filenames])
                LoadWANDSCD(Filename=filenames,
                            Grouping='4x4',
                            OutputWorkspace='data')
        else:
            filenames = [filepath.format(IPTS, run) for run in runs]
            if np.all([os.path.exists(filename) for filename in filenames]):
                filenames = ','.join([filename for filename in filenames])
                Load(Filename=filenames,
                     FilterByTofMin=1500,
                     FilterByTofMax=16600,
                     FilterByTimeStop=time_stop,
                     OutputWorkspace='data')

    def calibrate_data(self, instrument, det_cal, tube_cal):

        filepath = self.get_raw_file_path(instrument)

        if mtd.doesExist('data'):

            goniometers = self.get_goniometers(instrument)
            while len(goniometers) < 6:
                goniometers.append(None)

            SetGoniometer(Workspace='data',
                          Axis0=goniometers[0],
                          Axis1=goniometers[1],
                          Axis2=goniometers[2],
                          Average=False if 'HFIR' in filepath else True)

            if tube_cal != '' and os.path.exists(tube_cal):
                LoadNexus(Filename=tube_cal,
                          OutputWorkspace='tube_table')
                ApplyCalibration(Workspace='data',
                                 CalibrationTable='tube_table')

            if det_cal != '' and os.path.exists(det_cal):
                if os.path.splitext(det_cal)[1] == '.xml':
                    LoadParameterFile(Workspace='data',
                                      Filename=det_cal)
                else:
                    ws = mtd['data']
                    group = ws.isGroup()
                    if group:
                        for input_ws in ws.getNames():
                            LoadIsawDetCal(InputWorkspace=input_ws,
                                           Filename=det_cal)
                    else:
                        LoadIsawDetCal(InputWorkspace='data',
                                       Filename=det_cal)

    def convert_data(self, instrument, wavelength, lorentz):

        filepath = self.get_raw_file_path(instrument)

        if mtd.doesExist('data'):

            if 'HFIR' in filepath:
                ei = mtd['data'].getExperimentInfo(0)
                two_theta = ei.run().getProperty('TwoTheta').value
                az_phi = ei.run().getProperty('Azimuthal').value
                counts = mtd['data'].getSignalArray().copy()
                counts = counts.reshape(-1, counts.shape[2])
                Q_max = 4*np.pi/wavelength[0]*np.sin(0.5*max(two_theta))
                ConvertHFIRSCDtoMDE(InputWorkspace='data',
                                    Wavelength=wavelength[0],
                                    LorentzCorrection=lorentz,
                                    MinValues=[-Q_max,-Q_max,-Q_max],
                                    MaxValues=[+Q_max,+Q_max,+Q_max],
                                    OutputWorkspace='md')

            else:
                ConvertUnits(InputWorkspace='data',
                             Target='Wavelength',
                             OutputWorkspace='data')
                CropWorkspace(InputWorkspace='data',
                              XMin=wavelength[0],
                              XMax=wavelength[1],
                              OutputWorkspace='data')
                Rebin(InputWorkspace='data',
                      OutputWorkspace='data',
                      Params=[wavelength[0], 0.01, wavelength[1]])
                group = mtd['data'].isGroup()
                input_ws = mtd['data'].getNames()[0] if group else 'data'
                PreprocessDetectorsToMD(InputWorkspace=input_ws,
                                        OutputWorkspace='detectors')
                two_theta = mtd['detectors'].column('TwoTheta')
                az_phi = mtd['detectors'].column('Azimuthal')
                counts = mtd['data'].extractY().copy()
                Q_max = 4*np.pi/min(wavelength)*np.sin(0.5*max(two_theta))
                ConvertToMD(InputWorkspace='data',
                            QDimensions='Q3D',
                            dEAnalysisMode='Elastic',
                            Q3DFrames='Q_sample',
                            LorentzCorrection=lorentz,
                            MinValues=[-Q_max,-Q_max,-Q_max],
                            MaxValues=[+Q_max,+Q_max,+Q_max],
                            PreprocDetectorsWS='detectors',
                            OutputWorkspace='md')
                if group:
                    MergeMD(InputWorkspaces='md', OutputWorkspace='md')

            self.Q = 'md'

            BinMD(InputWorkspace=self.Q,
                  AlignedDim0='Q_sample_x,{},{},768'.format(-Q_max, Q_max),
                  AlignedDim1='Q_sample_y,{},{},768'.format(-Q_max, Q_max),
                  AlignedDim2='Q_sample_z,{},{},768'.format(-Q_max, Q_max),
                  OutputWorkspace='Q3D')

            signal = mtd['Q3D'].getSignalArray().copy()

            threshold = np.nanpercentile(signal[signal > 0], 95)
            mask = signal >= threshold

            self.signal = signal[mask]

            dims = [mtd['Q3D'].getDimension(i) for i in range(3)]

            x, y, z = [np.linspace(dim.getMinimum()+dim.getBinWidth()/2,
                                   dim.getMaximum()-dim.getBinWidth()/2,
                                   dim.getNBins()) for dim in dims]

            x, y, z = np.meshgrid(x, y, z, indexing='ij')

            self.x, self.y, self.z = x[mask], y[mask], z[mask]

            self.Qx_min, self.Qx_max = self.x.min(), self.x.max()
            self.Qy_min, self.Qy_max = self.y.min(), self.y.max()
            self.Qz_min, self.Qz_max = self.z.min(), self.z.max()

            self.wavelength = wavelength
            self.counts = counts

            kf_x = np.sin(two_theta)*np.cos(az_phi)
            kf_y = np.sin(two_theta)*np.sin(az_phi)
            kf_z = np.cos(two_theta)

            self.nu = np.rad2deg(np.arcsin(kf_y))
            self.gamma = np.rad2deg(np.arctan2(kf_x, kf_z))

    def instrument_view(self):

        return self.gamma, self.nu, self.counts.sum(axis=1)

    def get_slice_info(self, U, V, W, normal, value, thickness, width):

        UB = self.get_UB()

        if self.has_UB() and self.has_Q():
    
            Wt = np.column_stack([U, V, W])        
    
            slice_dict = {}
    
            Bp = np.dot(UB, Wt)

            bp_inv = np.linalg.inv(2*np.pi*Bp)

            corners = np.array([[self.Qx_min, self.Qy_min, self.Qz_min],
                                [self.Qx_max, self.Qy_min, self.Qz_min],
                                [self.Qx_min, self.Qy_max, self.Qz_min],
                                [self.Qx_max, self.Qy_max, self.Qz_min],
                                [self.Qx_min, self.Qy_min, self.Qz_max],
                                [self.Qx_max, self.Qy_min, self.Qz_max],
                                [self.Qx_min, self.Qy_max, self.Qz_max],
                                [self.Qx_max, self.Qy_max, self.Qz_max]])

            trans_corners = np.einsum('ij,kj->ki', bp_inv, corners)

            min_values = np.ceil(np.min(trans_corners, axis=0))
            max_values = np.floor(np.max(trans_corners, axis=0))

            bin_sizes = np.round((max_values-min_values)/width).astype(int)

            min_values = min_values.tolist()
            max_values = max_values.tolist()

            bin_sizes = bin_sizes.tolist()

            extents = []
            bins = []

            integrate = [value-thickness, value+thickness]

            for ind, i in enumerate(normal):
                if i == 0:
                    extents += [min_values[ind], max_values[ind]]
                    bins += [1+bin_sizes[ind]]
                else:
                    extents += integrate
                    bins += [1]

            ConvertQtoHKLMDHisto(InputWorkspace=self.Q,
                                 PeaksWorkspace=self.table,
                                 UProj=U,
                                 VProj=V,
                                 WProj=W,
                                 Extents=extents,
                                 Bins=bins,
                                 OutputWorkspace='slice')

            CompactMD(InputWorkspace='slice', OutputWorkspace='slice')

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
    
            # signal[signal <= 0] = np.nan
            # signal[np.isinf(signal)] = np.nan
    
            slice_dict['signal'] = signal
        
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

    def get_has_Q_vol(self):

        return mtd.doesExist('Q3D')

    def get_Q_info(self):

        Q_dict = {}

        if self.get_has_Q_vol():

            Q_dict['signal'] = self.signal

            Q_dict['x'] = self.x
            Q_dict['y'] = self.y
            Q_dict['z'] = self.z

        if self.has_peaks():

            self.sort_peaks_by_hkl(self.table)

            self.sort_peaks_by_d(self.table)

            Qs, Is, inds, pk_nos, Ts = [], [], [], [], []

            for j, peak in enumerate(mtd[self.table]):

                T = np.zeros((4,4))

                I = peak.getIntensity()

                ind = (peak.getHKL().norm2() > 0)*1.0

                shape = eval(peak.getPeakShape().toJSON())

                pk_no = peak.getPeakNumber()

                Q = peak.getQSampleFrame()

                if any(['radius' in key for key in shape.keys()]):

                    if shape.get('radius0') is not None:

                        r, v = [], []
                        for i in range(3):
                            r.append(shape['radius{}'.format(i)])
                            v.append(shape['direction{}'.format(i)].split(' '))
                        r = np.array(r)
                        v = np.array(v).T.astype(float)

                    else:

                        r = np.array([shape['radius']]*3)
                        v = np.eye(3)

                    P = np.dot(v, np.dot(np.diag(r), v.T))

                else:

                    P = np.zeros((3, 3))

                T[:3,-1] = Q
                T[:3,:3] = P
                T[-1,-1] = 1

                Qs.append(Q)
                Is.append(I)
                inds.append(ind)
                pk_nos.append(pk_no)
                Ts.append(T)

            Q_dict['coordinates'] = Qs
            Q_dict['intensities'] = Is
            Q_dict['indexings'] = inds
            Q_dict['numbers'] = pk_nos
            Q_dict['transforms'] = Ts

        return Q_dict if len(Q_dict.keys()) > 0 else None

    def get_lattice_constants(self):

        if self.has_UB():

            ol = mtd[self.table].sample().getOrientedLattice()

            params = ol.a(), ol.b(), ol.c(), ol.alpha(), ol.beta(), ol.gamma()

            return params

    def simplify_vector(self, vec):

        vec = vec/np.linalg.norm(vec)

        vec *= 10

        vec = np.round(vec).astype(int)

        return vec//np.gcd.reduce(vec)

    def get_sample_directions(self):

        if self.has_UB():

            UB = mtd[self.table].sample().getOrientedLattice().getUB()

            vecs = np.linalg.inv(UB).T

            return [self.simplify_vector(vec) for vec in vecs]

    def save_UB(self, filename):
        """
        Save UB to file.

        Parameters
        ----------
        filename : str
            Name of UB file with extension .mat.

        """

        SaveIsawUB(InputWorkspace=self.table, Filename=filename)

    def load_UB(self, filename):
        """
        Load UB from file.

        Parameters
        ----------
        filename : str
            Name of UB file with extension .mat.

        """

        if self.has_peaks():

            LoadIsawUB(InputWorkspace=self.table, Filename=filename)

    def determine_UB_with_niggli_cell(self, min_d, max_d, tol=0.1):
        """
        Determine UB with primitive lattice using min/max lattice constant.

        Parameters
        ----------
        min_d : float
            Minimum lattice parameter in ansgroms.
        max_d : float
            Maximum lattice parameter in angstroms.
        tol : float, optional
            Indexing tolerance. The default is 0.1.

        """

        FindUBUsingFFT(PeaksWorkspace=self.table,
                       MinD=min_d,
                       MaxD=max_d,
                       Tolerance=tol)

    def determine_UB_with_lattice_parameters(self, a,
                                                   b,
                                                   c,
                                                   alpha,
                                                   beta,
                                                   gamma,
                                                   tol=0.1):
        """
        Determine UB with prior known lattice parameters.

        Parameters
        ----------
        a, b, c : float
            Lattice constants in angstroms.
        alpha, beta, gamma : float
            Lattice angles in degrees.
        tol : float, optional
            Indexing tolerance. The default is 0.1.

        """

        FindUBUsingLatticeParameters(PeaksWorkspace=self.table,
                                     a=a,
                                     b=b,
                                     c=c,
                                     alpha=alpha,
                                     beta=beta,
                                     gamma=gamma,
                                     Tolerance=tol)

    def refine_UB_without_constraints(self, tol=0.1, sat_tol=None):
        """
        Refine UB with unconstrained lattice parameters.

        Parameters
        ----------
        tol : float, optional
            Indexing tolerance. The default is 0.1.
        sat_tol : float, optional
            Satellite indexing tolerance. The default is None.

        """

        tol_for_sat = sat_tol if sat_tol is not None else tol

        FindUBUsingIndexedPeaks(PeaksWorkspace=self.table,
                                Tolerance=tol,
                                ToleranceForSatellite=tol_for_sat)

    def refine_UB_with_constraints(self, cell, tol=0.1):
        """
        Refine UB with constraints corresponding to lattice system.

        +----------------+---------------+----------------------+
        | Lattice system | Lengths       | Angles               |
        +================+===============+======================+
        | Cubic          | :math:`a=b=c` | :math:`α=β=γ=90`     |
        +----------------+---------------+----------------------+
        | Hexagonal      | :math:`a=b`   | :math:`α=β=90,γ=120` |
        +----------------+---------------+----------------------+
        | Rhombohedral   | :math:`a=b=c` | :math:`α=β=γ`        |
        +----------------+---------------+----------------------+
        | Tetragonal     | :math:`a=b`   | :math:`α=β=γ=90`     |
        +----------------+---------------+----------------------+
        | Orthorhombic   | None          | :math:`α=β=γ=90`     |
        +----------------+---------------+----------------------+
        | Monoclinic     | None          | :math:`α=γ=90`       |
        +----------------+---------------+----------------------+
        | Triclinic      | None          | None                 |
        +----------------+---------------+----------------------+

        Parameters
        ----------
        cell : float
            Lattice system.
        tol : float, optional
            Indexing tolerance. The default is 0.1.

        """

        OptimizeLatticeForCellType(PeaksWorkspace=self.table,
                                   CellType=cell,
                                   Apply=True,
                                   Tolerance=tol)

    def refine_U_only(self, a, b, c, alpha, beta, gamma):
        """
        Refine the U orientation only.

        Parameters
        ----------
        a, b, c : float
            Lattice constants in angstroms.
        alpha, beta, gamma : float
            Lattice angles in degrees.

        """

        CalculateUMatrix(PeaksWorkspace=self.table,
                         a=a,
                         b=b,
                         c=c,
                         alpha=alpha,
                         beta=beta,
                         gamma=gamma)

    def select_cell(self, number, tol=0.1):
        """
        Transform to conventional cell using form number.

        Parameters
        ----------
        number : int
            Form number.
        tol : float, optional
            Indexing tolerance. The default is 0.1.

        """

        SelectCellWithForm(PeaksWorkspace=self.table,
                           FormNumber=number,
                           Apply=True,
                           Tolerance=tol)

    def possible_conventional_cells(self, max_error=0.2, permutations=True):
        """
        List possible conventional cells.

        Parameters
        ----------
        max_error : float, optional
            Max scalar error to report form numbers. The default is 0.2.
        permutations : bool, optional
            Allow permutations of the lattice. The default is True.

        Returns
        -------
        vals : list
            List of form results.

        """

        result = ShowPossibleCells(PeaksWorkspace=self.table,
                                   MaxScalarError=max_error,
                                   AllowPermutations=permutations,
                                   BestOnly=False)

        vals = [json.loads(cell) for cell in result.Cells]

        cells = []
        for i, val in enumerate(vals):
            form = val['FormNumber']
            error = val['Error']
            cell = val['CellType']
            centering = val['Centering']
            bravais = cell, centering
            a = val['a']
            b = val['b']
            c = val['c']
            alpha = val['alpha']
            beta = val['beta']
            gamma = val['gamma']
            vol = val['volume']
            params = a, b, c, alpha, beta, gamma, vol
            cell = form, error, bravais, params
            cells.append(cell)

        return cells

    def transform_lattice(self, transform, tol=0.1):
        """
        Apply a cell transformation to the lattice.

        Parameters
        ----------
        transform : 3x3 array-like
            Transform to apply to hkl values.
        tol : float, optional
            Indexing tolerance. The default is 0.1.

        """

        hkl_trans = ','.join(9*['{}']).format(*transform)

        TransformHKL(PeaksWorkspace=self.table,
                     Tolerance=tol,
                     HKLTransform=hkl_trans)

    def generate_lattice_transforms(self, cell):
        """
        Obtain possible transforms compatabile with a unit cell lattice.

        Parameters
        ----------
        cell : str
            Latttice system.

        Returns
        -------
        transforms : dict
            Transform dictionary with symmetry operation as key.

        """

        symbol = lattice_group[cell]

        pg = PointGroupFactory.createPointGroup(symbol)

        coords = np.eye(3).astype(int)

        transform = {}
        for symop in pg.getSymmetryOperations():
            T = np.column_stack([symop.transformHKL(vec) for vec in coords])
            if np.linalg.det(T) > 0:
                name = '{}: '.format(symop.getOrder())+symop.getIdentifier()
                transform[name] = T

        return {key:transform[key] for key in sorted(transform.keys())}

    def index_peaks(self, tol=0.1,
                          sat_tol=None,
                          mod_vec_1=[0,0,0],
                          mod_vec_2=[0,0,0],
                          mod_vec_3=[0,0,0],
                          max_order=0,
                          cross_terms=False,
                          round_hkl=True):
        """
        Index the peaks and calculate the lattice parameter uncertainties.

        Parameters
        ----------
        tol : float, optional
            Indexing tolerance. The default is 0.1.
        sat_tol : float, optional
            Satellite indexing tolerance. The default is None.
        mod_vec_1, mod_vec_2, mod_vec_3 : list, optional
            Modulation vectors. The default is [0,0,0].
        max_order : int, optional
            Maximum order greater than zero for satellites. The default is 0.
        cross_terms : bool, optional
            Include modulation cross terms. The default is False.
        round_hkl : bool, optional
            Round integers to integer. The default is True.

        Returns
        -------
        indexing : list
            Result of indexing including number indexed and errors.

        """

        tol_for_sat = sat_tol if sat_tol is not None else tol

        save = True if max_order > 0 else False

        indexing = IndexPeaks(PeaksWorkspace=self.table,
                              Tolerance=tol,
                              ToleranceForSatellite=tol_for_sat,
                              RoundHKLs=round_hkl,
                              CommonUBForAll=True,
                              ModVector1=mod_vec_1,
                              ModVector2=mod_vec_2,
                              ModVector3=mod_vec_3,
                              MaxOrder=max_order,
                              CrossTerms=cross_terms,
                              SaveModulationInfo=save)

        return indexing

    def calculate_hkl(self):
        """
        Calculate hkl values without rounding.

        """

        CalculatePeaksHKL(PeaksWorkspace=self.table,
                          OverWrite=True)

    def find_peaks(self, min_dist,
                         density=1000,
                         max_peaks=50,
                         edge_pixels=0):
        """
        Harvest strong peak locations from Q-sample into a peaks table.

        Parameters
        ----------
        min_dist : float
            Minimum distance enforcing lower limit of peak spacing.
        density : int, optional
            Threshold density. The default is 1000.
        max_peaks : int, optional
            Maximum number of peaks to find. The default is 50.
        edge_pixels: int, optional
            Nnumber of edge pixels to exclude. The default is 0.

        """

        FindPeaksMD(InputWorkspace=self.Q,
                    PeakDistanceThreshold=min_dist,
                    MaxPeaks=max_peaks,
                    PeakFindingStrategy='VolumeNormalization',
                    DensityThresholdFactor=density,
                    EdgePixels=edge_pixels,
                    OutputWorkspace=self.table)

        self.integrate_peaks(min_dist, 1, 1, method='sphere', centroid=False)

        self.clear_intensity()

    def centroid_peaks(self, peak_radius):
        """
        Re-center peak locations using centroid within given radius

        Parameters
        ----------
        peak_radius : float
            Integration region radius.

        """

        CentroidPeaksMD(InputWorkspace=self.Q,
                        PeakRadius=peak_radius,
                        PeaksWorkspace=self.table,
                        OutputWorkspace=self.table)

    def integrate_peaks(self, peak_radius,
                              background_inner_fact=1,
                              background_outer_fact=1.5,
                              method='sphere',
                              centroid=True):
        """
        Integrate peaks using spherical or ellipsoidal regions.
        Ellipsoid integration adapts itself to the peak distribution.

        Parameters
        ----------
        peak_radius : float
            Integration region radius.
        background_inner_fact : float, optional
            Factor of peak radius for background shell. The default is 1.
        background_outer_fact : float, optional
            Factor of peak radius for background shell. The default is 1.5.
        method : str, optional
            Integration method. The default is 'sphere'.
        centroid : str, optional
            Shift peak position to centroid. The default is True.

        """

        background_inner_radius = peak_radius*background_inner_fact
        background_outer_radius = peak_radius*background_outer_fact

        if method == 'sphere' and centroid:
            self.centroid_peaks(peak_radius)

        IntegratePeaksMD(InputWorkspace=self.Q,
                         PeaksWorkspace=self.table,
                         PeakRadius=peak_radius,
                         BackgroundInnerRadius=background_inner_radius,
                         BackgroundOuterRadius=background_outer_radius,
                         Ellipsoid=True if method == 'ellipsoid' else False,
                         FixQAxis=False,
                         FixMajorAxisLength=False,
                         UseCentroid=True,
                         MaxIterations=3,
                         ReplaceIntensity=True,
                         IntegrateIfOnEdge=False,
                         AdaptiveQBackground=False,
                         MaskEdgeTubes=False,
                         OutputWorkspace=self.table)

    def clear_intensity(self):

        for peak in mtd[self.table]:
            peak.setIntensity(0)
            peak.setSigmaIntensity(0)

    def get_max_d_spacing(self, ws):
        """
        Obtain the maximum d-spacing from the oriented lattice.

        Parameters
        ----------
        ws : str
            Workspace with UB defined on oriented lattice.

        Returns
        -------
        d_max : float
            Maximum d-spacing.

        """

        if HasUB(Workspace=ws):

            if hasattr(mtd[ws], 'sample'):
                ol = mtd[ws].sample().getOrientedLattice()
            else:
                for i in range(mtd[ws].getNumExperimentInfo()):
                    sample = mtd[ws].getExperimentInfo(i).sample()
                    if sample.hasOrientedLattice():
                        ol = sample.getOrientedLattice()
                        SetUB(Workspace=ws, UB=ol.getUB())
                ol = mtd[ws].getExperimentInfo(i).sample().getOrientedLattice()

            return 1/min([ol.astar(), ol.bstar(), ol.cstar()])

    def predict_peaks(self,
                      centering,
                      d_min,
                      lamda_min,
                      lamda_max,
                      edge_pixels=0):
        """
        Predict peak Q-sample locations with UB and lattice centering.

        +--------+-----------------------+
        | Symbol | Reflection condition  |
        +========+=======================+
        | P      | None                  |
        +--------+-----------------------+
        | I      | :math:`h+k+l=2n`      |
        +--------+-----------------------+
        | F      | :math:`h,k,l` unmixed |
        +--------+-----------------------+
        | R      | None                  |
        +--------+-----------------------+
        | R(obv) | :math:`-h+k+l=3n`     |
        +--------+-----------------------+
        | R(rev) | :math:`h-k+l=3n`      |
        +--------+-----------------------+
        | A      | :math:`k+l=2n`        |
        +--------+-----------------------+
        | B      | :math:`l+h=2n`        |
        +--------+-----------------------+
        | C      | :math:`h+k=2n`        |
        +--------+-----------------------+

        Parameters
        ----------
        centering : str
            Lattice centering that provides the reflection condition.
        d_min : float
            The lower d-spacing resolution to predict peaks.
        lamda_min, lamda_max : float
            The wavelength band over which to predict peaks.
        edge_pixels: int, optional
            Nnumber of edge pixels to exclude. The default is 0.

        """

        d_max = self.get_max_d_spacing(self.table)

        PredictPeaks(InputWorkspace=self.table,
                     WavelengthMin=lamda_min,
                     WavelengthMax=lamda_max,
                     MinDSpacing=d_min,
                     MaxDSpacing=d_max*1.2,
                     ReflectionCondition=centering_reflection[centering],
                     EdgePixels=edge_pixels,
                     OutputWorkspace=self.table)

        self.integrate_peaks(0.1, 1, 1, method='sphere', centroid=False)

        self.clear_intensity()

    def predict_modulated_peaks(self, d_min,
                                      lamda_min,
                                      lamda_max,
                                      mod_vec_1=[0,0,0],
                                      mod_vec_2=[0,0,0],
                                      mod_vec_3=[0,0,0],
                                      max_order=0,
                                      cross_terms=False):
        """
        Predict the modulated peak positions based on main peaks.

        Parameters
        ----------
        centering : str
            Lattice centering that provides the reflection condition.
        d_min : float
            The lower d-spacing resolution to predict peaks.
        lamda_min, lamda_max : float
            The wavelength band over which to predict peaks.
        mod_vec_1, mod_vec_2, mod_vec_3 : list, optional
            Modulation vectors. The default is [0,0,0].
        max_order : int, optional
            Maximum order greater than zero for satellites. The default is 0.
        cross_terms : bool, optional
            Include modulation cross terms. The default is False.

        """

        d_max = self.get_max_d_spacing(self.table)

        sat_peaks = self.table+'_sat'

        PredictSatellitePeaks(Peaks=self.table,
                              SatellitePeaks=sat_peaks,
                              ModVector1=mod_vec_1,
                              ModVector2=mod_vec_2,
                              ModVector3=mod_vec_3,
                              MaxOrder=max_order,
                              CrossTerms=cross_terms,
                              IncludeIntegerHKL=False,
                              IncludeAllPeaksInRange=True,
                              WavelengthMin=lamda_min,
                              WavelengthMax=lamda_max,
                              MinDSpacing=d_min,
                              MaxDSpacing=d_max*1.2)

        CombinePeaksWorkspaces(LHSWorkspace=self.table,
                               RHSWorkspace=sat_peaks,
                               OutputWorkspace=self.table)

        DeleteWorkspace(Workspace=sat_peaks)

    def predict_satellite_peaks(self, lamda_min,
                                      lamda_max,
                                      d_min,
                                      mod_vec_1=[0,0,0],
                                      mod_vec_2=[0,0,0],
                                      mod_vec_3=[0,0,0],
                                      max_order=0,
                                      cross_terms=False):
        """
        Locate satellite peaks from goniometer angles.

        Parameters
        ----------
        d_min : float
            The lower d-spacing resolution to predict peaks.
        lamda_min : float
            Minimum wavelength.
        lamda_max : float
            Maximum wavelength.
        mod_vec_1, mod_vec_2, mod_vec_3 : list, optional
            Modulation vectors. The default is [0,0,0].
        max_order : int, optional
            Maximum order greater than zero for satellites. The default is 0.
        cross_terms : bool, optional
            Include modulation cross terms. The default is False.

        """

        Rs = self.get_all_goniometer_matrices(self.Q)

        for R in Rs:

            self.set_goniometer(self.table, R)

            self.predict_modulated_peaks(d_min,
                                         lamda_min,
                                         lamda_max,
                                         mod_vec_1,
                                         mod_vec_2,
                                         mod_vec_3,
                                         max_order,
                                         cross_terms)

            self.remove_duplicate_peaks(self.table)

    def sort_peaks_by_hkl(self, peaks):
        """
        Sort peaks table by descending hkl values.

        Parameters
        ----------
        peaks : str
            Name of peaks table.

        """

        columns = ['l', 'k', 'h']

        for col in columns:

            SortPeaksWorkspace(InputWorkspace=peaks,
                               ColumnNameToSortBy=col,
                               SortAscending=False,
                               OutputWorkspace=peaks)

    def sort_peaks_by_d(self, peaks):
        """
        Sort peaks table by descending d-spacing.

        Parameters
        ----------
        peaks : str
            Name of peaks table.

        """

        SortPeaksWorkspace(InputWorkspace=peaks,
                           ColumnNameToSortBy='DSpacing',
                           SortAscending=False,
                           OutputWorkspace=peaks)

    def remove_duplicate_peaks(self, peaks):
        """
        Omit duplicate peaks from different based on indexing.
        Table will be sorted.

        Parameters
        ----------
        peaks : str
            Name of peaks table.

        """

        self.sort_peaks_by_hkl(peaks)

        for no in range(mtd[peaks].getNumberPeaks()-1,0,-1):

            if (mtd[peaks].getPeak(no).getHKL()-\
                mtd[peaks].getPeak(no-1).getHKL()).norm2() == 0:

                DeleteTableRows(TableWorkspace=peaks, Rows=no)

    def get_all_goniometer_matrices(self, ws):
        """
        Extract all goniometer matrices.

        Parameters
        ----------
        ws : str
            Name of workspace with goniometer indexing.

        Returns
        -------
        Rs: list
            Goniometer matrices.

        """

        Rs = []

        for ei in range(mtd[ws].getNumExperimentInfo()):

            run = mtd[ws].getExperimentInfo(ei).run()

            n_gon = run.getNumGoniometers()

            Rs += [run.getGoniometer(i).getR() for i in range(n_gon)]

        return np.array(Rs)

    def renumber_runs_by_index(self, ws, peaks):
        """
        Re-label the runs by index based on goniometer setting.

        Parameters
        ----------
        ws : str
            Name of workspace with goniometer indexing.
        peaks : str
            Name of peaks table.

        """

        Rs = self.get_all_goniometer_matrices(ws)

        for no in range(mtd[peaks].getNumberPeaks()):

            peak = mtd[peaks].getPeak(no)

            R = peak.getGoniometerMatrix()

            ind = np.isclose(Rs, R).all(axis=(1,2))
            i = -1 if not np.any(ind) else ind.tolist().index(True)

            peak.setRunNumber(i+1)

    def load_Q(self, filename):
        """
        Load Q file.

        Parameters
        ----------
        filename : str
            Name of Q file with extension .nxs.

        """

        LoadMD(Filename=filename,
               OutputWorkspace=self.table)

    def save_Q(self, filename):
        """
        Save Q file.

        Parameters
        ----------
        filename : str
            Name of Q file with extension .nxs.

        """

        SaveMD(Filename=filename,
               InputWorkspace=self.table)

    def load_peaks(self, filename):
        """
        Load peaks file.

        Parameters
        ----------
        filename : str
            Name of peaks file with extension .nxs.

        """

        LoadNexus(Filename=filename,
                  OutputWorkspace=self.table)

    def save_peaks(self, filename):
        """
        Save peaks file.

        Parameters
        ----------
        filename : str
            Name of peaks file with extension .nxs.

        """

        SaveNexus(Filename=filename,
                  InputWorkspace=self.table)

    def delete_peaks(self, peaks):
        """
        Remove peaks.

        Parameters
        ----------
        peaks : str
            Name of peaks table to be added.

        """

        if mtd.doesExist(peaks):

            DeleteWorkspace(Workspace=peaks)

    def filter_peaks(self, name, operator, value):
        """
        Filter out peaks based on value and operator.

        Parameters
        ----------
        name : str
            Filter name.
        operator : float
            Filter operator.
        value : float
            The filter value.

        """

        FilterPeaks(InputWorkspace=self.table,
                    OutputWorkspace=self.table,
                    FilterVariable=variable[name],
                    FilterValue=value,
                    Operator=operator,
                    Criterion='!=',
                    BankName='None')

    def get_modulation_info(self):

        if self.has_peaks() and self.has_UB():

            ol = mtd[self.table].sample().getOrientedLattice()

            return [ol.getModVec(i) for i in range(3)]

    def get_peak_info(self):

        if self.has_peaks():

            banks = mtd[self.table].column('BankName')

            peak_info = []
            for i, peak in enumerate(mtd[self.table]):
                bank = banks[i]
                intens = peak.getIntensity()
                sigma = peak.getSigmaIntensity()
                signal_noise = peak.getIntensityOverSigma()
                hkl = list(peak.getHKL())
                int_hkl = list(peak.getIntHKL())
                int_mnp = list(peak.getIntMNP())
                lamda = peak.getWavelength()
                d = peak.getDSpacing()
                run = peak.getRunNumber()
                row = peak.getRow()
                col = peak.getCol()
                ind = peak.getHKL().norm2() > 0
                vals = hkl, d, lamda, intens, signal_noise, \
                       sigma, int_hkl, int_mnp, run, bank, row, col, ind
                peak_info.append(vals)

            self.peak_info = peak_info

            return peak_info

    def get_peak(self, i):

        if self.peak_info is not None:
            return self.peak_info[i]

    def set_peak(self, i, hkl, int_hkl, int_mnp):

        peak = self.table.getPeak(i)

        peak.setHKL(*hkl)
        peak.setIntHKL(V3D(*int_hkl))
        peak.setIntMNP(V3D(*int_mnp))

    def calculate_peaks(self, hkl_1, hkl_2, a, b, c, alpha, beta, gamma):

        uc = UnitCell(a, b, c, alpha, beta, gamma)

        d_1 = d_2 = phi_12 = None
        if hkl_1 is not None:
            d_1 = uc.d(*hkl_1)
        if hkl_2 is not None:
            d_2 = uc.d(*hkl_2)
        if hkl_1 is not None and hkl_2 is not None:
            phi_12 = uc.recAngle(*hkl_1, *hkl_2)

        return d_1, d_2, phi_12