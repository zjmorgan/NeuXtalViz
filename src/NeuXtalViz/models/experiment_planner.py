import os

import csv

from mantid.simpleapi import (
    CreatePeaksWorkspace,
    PredictPeaks,
    FilterPeaks,
    CombinePeaksWorkspaces,
    SortPeaksWorkspace,
    AddPeakHKL,
    CountReflections,
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
    MaskBTP,
    AddSampleLog,
    CreateSampleWorkspace,
    CreateEmptyTableWorkspace,
    DeleteTableRows,
    CloneWorkspace,
    DeleteWorkspace,
    RenameWorkspace,
    RenameWorkspaces,
    HasUB,
    mtd,
)

from mantid.kernel import V3D
from mantid.geometry import PointGroupFactory

from collections import defaultdict

import numpy as np
from scipy.spatial.transform import Rotation

from NeuXtalViz.models.base_model import NeuXtalVizModel
from NeuXtalViz.config.instruments import beamlines

# lattice_centering_dict = {
#     'P': 'Primitive',
#     'C': 'C-face centred',
#     'A': 'A-face centred',
#     'B': 'B-face centred',
#     'I': 'Body centred',
#     'F': 'All-face centred',
#     'R': 'Primitive',
#     'Robv': 'Rhombohedrally centred, obverse',
#     'Rrev': 'Rhombohedrally centred, reverse',
# }

# point_group_dict = {
#     '1': '1 (Triclinic)',
#     '-1': '-1 (Triclinic)',
#     '2': '2 (Monoclinic, unique axis b)',
#     'm': 'm (Monoclinic, unique axis b)',
#     '2/m': '2/m (Monoclinic, unique axis b)',
#     '112': '112 (Monoclinic, unique axis c)',
#     '11m': '11m (Monoclinic, unique axis c)',
#     '112/m': '112/m (Monoclinic, unique axis c)',
#     '222': '222 (Orthorhombic)',
#     'mm2': 'mm2 (Orthorhombic)',
#     'mmm': 'mmm (Orthorhombic)',
#     '4': '4 (Tetragonal)',
#     '-4': '-4 (Tetragonal)',
#     '4/m': '4/m (Tetragonal)',
#     '422': '422 (Tetragonal)',
#     '4mm': '4mm (Tetragonal)',
#     '-42m': '-42m (Tetragonal)',
#     '-4m2': '-4m2 (Tetragonal)',
#     '4/mmm': '4/mmm (Tetragonal)',
#     '3 r': '3 r (Trigonal - Rhombohedral)',
#     '-3 r': '-3 r (Trigonal - Rhombohedral)',
#     '32 r': '32 r (Trigonal - Rhombohedral)',
#     '3m r': '3m r (Trigonal - Rhombohedral)',
#     '-3m r': '-3m r (Trigonal - Rhombohedral)',
#     '3': '3 (Trigonal - Hexagonal)',
#     '-3': '-3 (Trigonal - Hexagonal)',
#     '312': '312 (Trigonal - Hexagonal)',
#     '31m': '31m (Trigonal - Hexagonal)',
#     '32': '32 (Trigonal - Hexagonal)',
#     '321': '321 (Trigonal - Hexagonal)',
#     '3m': '3m (Trigonal - Hexagonal)',
#     '-31m': '-31m (Trigonal - Hexagonal)',
#     '-3m': '-3m (Trigonal - Hexagonal)',
#     '-3m1': '-3m1 (Trigonal - Hexagonal)',
#     '6': '6 (Hexagonal)',
#     '-6': '-6 (Hexagonal)',
#     '6/m': '6/m (Hexagonal)',
#     '622': '622 (Hexagonal)',
#     '6mm': '6mm (Hexagonal)',
#     '-62m': '-62m (Hexagonal)',
#     '-6m2': '-6m2 (Hexagonal)',
#     '6/mmm': '6/mmm (Hexagonal)',
#     '23': '23 (Cubic)',
#     'm-3': 'm-3 (Cubic)',
#     '432': '432 (Cubic)',
#     '-43m': '-43m (Cubic)',
#     'm-3m': 'm-3m (Cubic)',
# }

point_group_centering = {
    "1": ["P"],
    "-1": ["P"],
    "2": ["P", "C"],
    "m": ["P", "C"],
    "2/m": ["P", "C"],
    "112": ["P", "C"],
    "11m": ["P", "C"],
    "112/m": ["P", "C"],
    "222": ["P", "I", "C", "A", "B"],
    "mm2": ["P", "I", "C", "A", "B"],
    "mmm": ["P", "I", "C", "A", "B"],
    "4": ["P", "I"],
    "-4": ["P", "I"],
    "4/m": ["P", "I"],
    "422": ["P", "I"],
    "4mm": ["P", "I"],
    "-42m": ["P", "I"],
    "-4m2": ["P", "I"],
    "4/mmm": ["P", "I"],
    "3 r": ["R"],
    "-3 r": ["R"],
    "32 r": ["R"],
    "3m r": ["R"],
    "-3m r": ["R"],
    "3": ["Robv", "Rrev"],
    "-3": ["Robv", "Rrev"],
    "312": ["Robv", "Rrev"],
    "31m": ["Robv", "Rrev"],
    "32": ["Robv", "Rrev"],
    "321": ["Robv", "Rrev"],
    "3m": ["Robv", "Rrev"],
    "-31m": ["Robv", "Rrev"],
    "-3m": ["Robv", "Rrev"],
    "-3m1": ["Robv", "Rrev"],
    "6": ["P"],
    "-6": ["P"],
    "6/m": ["P"],
    "622": ["P"],
    "6mm": ["P"],
    "-62m": ["P"],
    "-6m2": ["P"],
    "6/mmm": ["P"],
    "23": ["P", "I", "F"],
    "m-3": ["P", "I", "F"],
    "432": ["P", "I", "F"],
    "-43m": ["P", "I", "F"],
    "m-3m": ["P", "I", "F"],
}

crystal_system_point_groups = {
    "Triclinic": ["1", "-1"],
    "Monoclinic": ["2", "m", "2/m", "112", "11m", "112/m"],
    "Orthorhombic": ["222", "mm2", "mmm"],
    "Tetragonal": ["4", "-4", "4/m", "422", "4mm", "-42m", "-4m2", "4/mmm"],
    "Trigonal/Rhombohedral": ["3 r", "-3 r", "32 r", "3m r", "-3m r"],
    "Trigonal/Hexagonal": [
        "3",
        "-3",
        "312",
        "31m",
        "32",
        "321",
        "3m",
        "-31m",
        "-3m",
        "-3m1",
    ],
    "Hexagonal": ["6", "-6", "6/m", "622", "6mm", "-62m", "-6m2", "6/mmm"],
    "Cubic": ["23", "m-3", "432", "-43m", "m-3m"],
}

centering_conditions = {
    "P": lambda h, k, l: True,
    "I": lambda h, k, l: (h + k + l) % 2 == 0,
    "F": lambda h, k, l: (h % 2 == k % 2 == l % 2),
    "C": lambda h, k, l: k % 2 == 0 and l % 2 == 0,
    "A": lambda h, k, l: h % 2 == 0 and l % 2 == 0,
    "B": lambda h, k, l: h % 2 == 0 and k % 2 == 0,
    "R": lambda h, k, l: True,
    "Robv": lambda h, k, l: (-h + k + l) % 3 == 0,
    "Rrev": lambda h, k, l: (h - k + l) % 3 == 0,
}


class ExperimentModel(NeuXtalVizModel):
    def __init__(self):
        super(ExperimentModel, self).__init__()

        CreatePeaksWorkspace(
            NumberOfPeaks=0,
            OutputType="LeanElasticPeak",
            OutputWorkspace="coverage",
        )

    def initialize_instrument(self, instrument, logs):
        inst = self.get_instrument_name(instrument)

        if not mtd.doesExist("instrument"):
            LoadEmptyInstrument(
                InstrumentName=inst, OutputWorkspace="instrument"
            )

            for key in logs.keys():
                AddSampleLog(
                    Workspace="instrument",
                    LogName=key,
                    LogText=str(logs[key]),
                    LogType="Number Series",
                    NumberType="Double",
                )

            LoadInstrument(
                Workspace="instrument",
                RewriteSpectraMap=False,
                InstrumentName=inst,
            )

            ExtractMonitors(
                InputWorkspace="instrument",
                MonitorWorkspace="monitors",
                DetectorWorkspace="instrument",
            )

            cols, rows = beamlines[instrument]["BankPixels"]
            mask_cols, mask_rows = beamlines[instrument]["MaskEdges"]

            MaskBTP(
                Workspace="instrument",
                Instrument=inst,
                Tube="0-{},{}-{}".format(mask_cols, cols - mask_cols, cols),
            )

            MaskBTP(
                Workspace="instrument",
                Instrument=inst,
                Pixel="0-{},{}-{}".format(mask_rows, rows - mask_rows, rows),
            )

            banks = beamlines[instrument]["MaskBanks"]

            for bank in banks:
                MaskBTP(
                    Workspace="instrument",
                    Instrument=inst,
                    Bank=bank,
                )

            PreprocessDetectorsToMD(
                InputWorkspace="instrument", OutputWorkspace="detectors"
            )

            grouping = beamlines[instrument]["Grouping"]

            c, r = [int(val) for val in grouping.split("x")]
            shape = (-1, cols, rows)

            det_map = np.array(mtd["detectors"].column(5)).reshape(*shape)

            shape = det_map.shape
            i, j, k = np.meshgrid(
                np.arange(shape[0]),
                np.arange(shape[1]),
                np.arange(shape[2]),
                indexing="ij",
            )
            mask = np.array(mtd["detectors"].column(7)) == 0

            keys = np.stack((i, j // c, k // r), axis=-1)
            keys_flat = keys.reshape(-1, keys.shape[-1])
            det_map_flat = det_map.ravel()[mask].astype(str)
            grouped_ids = defaultdict(list)

            for key, detector_id in zip(map(tuple, keys_flat), det_map_flat):
                grouped_ids[key].append(detector_id)

            detector_list = ",".join(
                "+".join(group) for group in grouped_ids.values()
            )

            GroupDetectors(
                InputWorkspace="instrument",
                OutputWorkspace="instrument",
                GroupingPattern=detector_list,
            )

            CreatePeaksWorkspace(
                InstrumentWorkspace="instrument",
                NumberOfPeaks=0,
                OutputType="LeanElasticPeak",
                OutputWorkspace="peak",
            )

            CreatePeaksWorkspace(
                InstrumentWorkspace="instrument",
                NumberOfPeaks=0,
                OutputType="Peak",
                OutputWorkspace="peaks",
            )

            CreatePeaksWorkspace(
                InstrumentWorkspace="instrument",
                NumberOfPeaks=0,
                OutputType="Peak",
                OutputWorkspace="combined",
            )

            L2 = np.array(mtd["detectors"].column(1))[mask]
            tt = np.array(mtd["detectors"].column(2))[mask]
            az = np.array(mtd["detectors"].column(3))[mask]
            det_ID = np.array(mtd["detectors"].column(4))[mask]

            x = L2 * np.sin(tt) * np.cos(az)
            y = L2 * np.sin(tt) * np.sin(az)
            z = L2 * np.cos(tt)

            self.det_ID = det_ID.copy()
            self.nu = np.rad2deg(np.arcsin(y / L2))
            self.gamma = np.rad2deg(np.arctan2(x, z))

    def remove_instrument(self):
        if mtd.doesExist("instrument"):
            DeleteWorkspace(Workspace="instrument")

        if mtd.doesExist("cobmined"):
            DeleteWorkspace(Workspace="cobmined")

        if mtd.doesExist("filtered"):
            DeleteWorkspace(Workspace="filtered")

    def get_crystal_system_point_groups(self, crystal_system):
        return crystal_system_point_groups[crystal_system]

    def get_point_group_centering(self, point_group):
        return point_group_centering[point_group]

    def get_symmetry(self, point_group, centering):
        return str(point_group), str(centering)

    def create_plan(self, table):
        pv, names, titles, settings, comments, counts, values, use = table

        CreateEmptyTableWorkspace(OutputWorkspace="plan")

        mtd["plan"].addColumn("str", pv)

        for name in names:
            mtd["plan"].addColumn("float", name)

        mtd["plan"].addColumn("str", "Comment")
        mtd["plan"].addColumn("str", "Wait For")
        mtd["plan"].addColumn("float", "Value")
        mtd["plan"].addColumn("bool", "Use")

        for title, setting, comment, count, value, active in zip(
            titles, settings, comments, counts, values, use
        ):
            row = {}
            row[pv] = title
            for angle, name in zip(setting, names):
                row[name] = np.round(angle, 2)
            row["Comment"] = comment
            row["Wait For"] = count
            row["Value"] = value
            row["Use"] = active
            mtd["plan"].addRow(row)

    def create_sample(self, instrument, mode, UB, wavelength, d_min):
        CreateSampleWorkspace(OutputWorkspace="sample")

        SetUB(Workspace="sample", UB=UB)

        AddSampleLog(
            Workspace="sample",
            LogName="instrument",
            LogText=instrument,
            LogType="String",
        )

        AddSampleLog(
            Workspace="sample",
            LogName="mode",
            LogText=mode,
            LogType="String",
        )

        AddSampleLog(
            Workspace="sample",
            LogName="lamda_min",
            LogText=str(wavelength[0]),
            LogType="Number",
            NumberType="Double",
        )

        AddSampleLog(
            Workspace="sample",
            LogName="lamda_max",
            LogText=str(wavelength[1]),
            LogType="Number",
            NumberType="Double",
        )

        AddSampleLog(
            Workspace="sample",
            LogName="d_min",
            LogText=str(d_min),
            LogType="Number",
            NumberType="Double",
        )

    def update_sample(self, crytsal_system, point_group, lattice_centering):
        if mtd.doesExist("sample"):
            AddSampleLog(
                Workspace="sample",
                LogName="crystal_system",
                LogText=crytsal_system,
                LogType="String",
            )

            AddSampleLog(
                Workspace="sample",
                LogName="point_group",
                LogText=point_group,
                LogType="String",
            )

            AddSampleLog(
                Workspace="sample",
                LogName="lattice_centering",
                LogText=lattice_centering,
                LogType="String",
            )

    def update_goniometer_motors(self, limits, motors):
        if mtd.doesExist("sample"):
            mtd["sample"].run()["limits"] = np.array(limits).flatten().tolist()

            values = []
            for key in motors.keys():
                values.append(motors[key])
            if len(values) > 0:
                mtd["sample"].run()["motors"] = values

    def load_UB(self, filename):
        LoadIsawUB(InputWorkspace="coverage", Filename=filename)

        self.copy_UB()

    def get_UB(self):
        if self.has_UB():
            return mtd["coverage"].sample().getOrientedLattice().getUB().copy()

    def copy_UB(self):
        UB = self.get_UB()
        if UB is not None:
            self.set_UB(UB)

    def has_UB(self):
        if HasUB(Workspace="coverage"):
            return True
        else:
            return False

    def get_instrument_name(self, instrument):
        return beamlines[instrument]["Name"]

    def get_modes(self, instrument):
        return list(beamlines[instrument]["Goniometer"].keys())

    def get_counting_options(self, instrument):
        return beamlines[instrument]["Counting"]

    def get_scan_log(self, instrument):
        return beamlines[instrument]["Title"]

    def get_axes_polarities(self, instrument, mode):
        goniometers = beamlines[instrument]["Goniometer"][mode]

        axes = [goniometers[name][:-3] for name in goniometers.keys()]

        polarities = [goniometers[name][3] for name in goniometers.keys()]

        return axes, polarities

    def get_goniometer_axes(self, instrument, mode):
        goniometers = beamlines[instrument]["Goniometer"][mode]

        axes = [
            "{}," + ",".join(np.array(goniometers[name][:-2]).astype(str))
            for name in goniometers.keys()
        ]

        return axes

    def get_goniometers(self, instrument, mode):
        goniometers = beamlines[instrument]["Goniometer"][mode]

        return [(name, *goniometers[name][-2:]) for name in goniometers.keys()]

    def get_motors(self, instrument):
        motors = beamlines[instrument].get("Motor")

        if motors is not None:
            return [(name, motors[name]) for name in motors.keys()]
        else:
            return []

    def get_wavelength(self, instrument):
        return beamlines[instrument]["Wavelength"]

    def save_plan(self, filename):
        plan_dict = mtd["plan"].toDict().copy()
        use_angle = plan_dict["Use"]

        for key in plan_dict.keys():
            items = plan_dict[key]
            items = [item for item, use in zip(items, use_angle) if use]
            plan_dict[key] = items

        plan_dict.pop("Use")

        with open(filename, mode="w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=plan_dict.keys())
            writer.writeheader()
            for row in zip(*plan_dict.values()):
                writer.writerow(dict(zip(plan_dict.keys(), row)))

    def save_experiment(self, filename):
        if mtd.doesExist("plan"):
            SaveNexus(InputWorkspace="plan", Filename=filename)
            if mtd.doesExist("sample"):
                SaveNexus(
                    InputWorkspace="sample", Filename=filename, Append=True
                )

    def load_experiment(self, filename):
        LoadNexus(Filename=filename, OutputWorkspace="experiment")

        plan, sample = mtd["experiment"].getNames()

        UB = mtd[sample].sample().getOrientedLattice().getUB().copy()
        SetUB(Workspace="coverage", UB=UB)

        self.set_UB(UB)

        instrument = mtd[sample].run().getProperty("instrument").value
        mode = mtd[sample].run().getProperty("mode").value
        wl_min = mtd[sample].run().getProperty("lamda_min").value
        wl_max = mtd[sample].run().getProperty("lamda_max").value
        d_min = mtd[sample].run().getProperty("d_min").value
        cs = mtd[sample].run().getProperty("crystal_system").value
        pg = mtd[sample].run().getProperty("point_group").value
        lc = mtd[sample].run().getProperty("lattice_centering").value
        lims = mtd[sample].run().getProperty("limits").value
        lims = np.array(lims).reshape(-1, 2).tolist()
        vals = []
        if mtd[sample].run().hasProperty("motors"):
            vals = mtd[sample].run().getProperty("motors").value

        if np.isclose(wl_min, wl_max):
            wl = wl_min
        else:
            wl = [wl_min, wl_max]

        cols = mtd[plan].columnCount() - 5
        rows = mtd[plan].rowCount()

        titles = mtd[plan].column(0)
        comments = mtd[plan].column(cols + 1)
        counts = mtd[plan].column(cols + 2)
        values = mtd[plan].column(cols + 3)
        use = mtd[plan].column(cols + 4)

        settings = []
        for row in range(rows):
            angles = []
            for col in range(cols):
                angle = mtd[plan].cell(row, col + 1)
                angles.append(angle)
            settings.append(angles)

        plan = (titles, settings, comments, counts, values, use)
        config = (instrument, mode, wl, d_min, lims, vals)
        symm = (cs, pg, lc)

        return plan, config, symm

    def generate_axes(self, axes, polarities):
        self.axes = [None] * 6

        for i, (axis, polarity) in enumerate(zip(axes, polarities)):
            self.axes[i] = "{},"
            self.axes[i] += ",".join(np.array([*axis, polarity]).astype(str))

    def get_setting(self, free_angles, limits):
        setting = []
        col = 0
        for limit in limits:
            if np.isclose(limit[0], limit[1]):
                setting.append(limit[0])
            else:
                setting.append(free_angles[col])
                col += 1
        return setting

    def _calculate_matrices(self, axes, polarities, limits, step):
        self.generate_axes(axes, polarities)

        angular_coverage = []
        for limit in limits:
            angular_coverage.append(np.arange(limit[0], limit[1] + step, step))

        axes = np.array(axes)
        polarities = np.array(polarities)

        angle_settings = np.meshgrid(*angular_coverage, indexing="ij")
        angle_settings = np.reshape(angle_settings, (len(polarities), -1)).T

        self.angles = angle_settings.copy()

        angle_settings = angle_settings * polarities
        angle_settings = np.deg2rad(angle_settings)

        rotation_vectors = angle_settings[..., None] * axes
        rotation_vectors = rotation_vectors.reshape(-1, 3)

        all_rotations = Rotation.from_rotvec(rotation_vectors).as_matrix()
        all_rotations = all_rotations.reshape(*angle_settings.shape, 3, 3)

        Rs = []
        for i in range(all_rotations.shape[0]):
            R = np.eye(3)
            for j in range(all_rotations.shape[1]):
                R = R @ all_rotations[i, j, :, :]
            Rs.append(R)

        return Rs

    def individual_peak(
        self, hkl, wavelength, axes, polarities, limits, step=1
    ):
        self.comment = "(" + " ".join(np.array(hkl).astype(str)) + ")"

        if np.isclose(wavelength[0], wavelength[1]):
            wavelength = [0.975 * wavelength[0], 1.025 * wavelength[1]]

        FilterPeaks(
            InputWorkspace="peak",
            OutputWorkspace="peak",
            FilterVariable="RunNumber",
            FilterValue=-1,
            Operator="=",
        )

        UB = mtd["coverage"].sample().getOrientedLattice().getUB().copy()

        SetUB(Workspace="peak", UB=UB)

        AddPeakHKL(Workspace="peak", HKL=hkl)

        Q_sample = mtd["peak"].getPeak(0).getQSampleFrame()

        Q = np.sqrt(np.dot(Q_sample, Q_sample))

        Rs = self._calculate_matrices(axes, polarities, limits, step)

        mtd["peak"].run().getGoniometer().setR(np.eye(3))
        mtd["peaks"].run().getGoniometer().setR(np.eye(3))

        Q_lab = np.einsum("kij,j->ki", Rs, Q_sample)

        lamda = -4 * np.pi * Q_lab[:, 2] / Q**2
        mask = (lamda > wavelength[0]) & (lamda < wavelength[1])

        k = 2 * np.pi / lamda

        ki = k[:, np.newaxis] * np.array([0, 0, 1])
        kf = Q_lab + ki

        gamma = np.rad2deg(np.arctan2(kf[:, 0], kf[:, 2]))[mask]
        nu = np.rad2deg(np.arcsin(kf[:, 1] / k))[mask]
        lamda = lamda[mask]

        self.angles = self.angles[mask]
        self.angles_gamma = gamma.copy()
        self.angles_nu = nu.copy()

        if len(lamda) > 0:
            k = 2 * np.pi / lamda
            Qx = k * np.cos(np.deg2rad(nu)) * np.sin(np.deg2rad(gamma))
            Qy = k * np.sin(np.deg2rad(nu))
            Qz = k * (np.cos(np.deg2rad(nu)) * np.cos(np.deg2rad(gamma)) - 1)

            mask = []
            for i in range(len(k)):
                peak = mtd["peaks"].createPeak(V3D(Qx[i], Qy[i], Qz[i]))
                mask.append(peak.getDetectorID() in self.det_ID)

            mask = np.array(mask)

            gamma = gamma[mask]
            nu = nu[mask]
            lamda = lamda[mask]

            self.angles = self.angles[mask]
            self.angles_gamma = gamma.copy()
            self.angles_nu = nu.copy()

        return gamma, nu, lamda

    def simultaneous_peaks(
        self, hkl_1, hkl_2, wavelength, axes, polarities, limits, step=1
    ):
        self.comment = "(" + " ".join(np.array(hkl_1).astype(str)) + ")"

        if np.isclose(wavelength[0], wavelength[1]):
            wavelength = [0.975 * wavelength[0], 1.025 * wavelength[1]]

        FilterPeaks(
            InputWorkspace="peak",
            OutputWorkspace="peak",
            FilterVariable="RunNumber",
            FilterValue=-1,
            Operator="=",
        )

        UB = mtd["coverage"].sample().getOrientedLattice().getUB().copy()

        SetUB(Workspace="peak", UB=UB)

        AddPeakHKL(Workspace="peak", HKL=hkl_1)
        AddPeakHKL(Workspace="peak", HKL=hkl_2)

        Q0_sample = mtd["peak"].getPeak(0).getQSampleFrame()
        Q1_sample = mtd["peak"].getPeak(1).getQSampleFrame()

        Q0 = np.sqrt(np.dot(Q0_sample, Q0_sample))
        Q1 = np.sqrt(np.dot(Q1_sample, Q1_sample))

        Rs = self._calculate_matrices(axes, polarities, limits, step)

        Q0_lab = np.einsum("kij,j->ki", Rs, Q0_sample)
        Q1_lab = np.einsum("kij,j->ki", Rs, Q1_sample)

        lamda0 = -4 * np.pi * Q0_lab[:, 2] / Q0**2
        lamda1 = -4 * np.pi * Q1_lab[:, 2] / Q1**2

        mask = (
            (lamda0 > wavelength[0])
            & (lamda0 < wavelength[1])
            & (lamda1 > wavelength[0])
            & (lamda1 < wavelength[1])
        )

        k0 = 2 * np.pi / lamda0
        k1 = 2 * np.pi / lamda1

        k0i = k0[:, np.newaxis] * np.array([0, 0, 1])
        k1i = k1[:, np.newaxis] * np.array([0, 0, 1])

        k0f = Q0_lab + k0i
        k1f = Q1_lab + k1i

        gamma0 = np.rad2deg(np.arctan2(k0f[:, 0], k0f[:, 2]))[mask]
        gamma1 = np.rad2deg(np.arctan2(k1f[:, 0], k1f[:, 2]))[mask]

        nu0 = np.rad2deg(np.arcsin(k0f[:, 1] / k0))[mask]
        nu1 = np.rad2deg(np.arcsin(k1f[:, 1] / k1))[mask]

        lamda0 = lamda0[mask]
        lamda1 = lamda1[mask]

        self.angles = self.angles[mask]
        self.angles_gamma = gamma0.copy()
        self.angles_nu = nu0.copy()

        if len(lamda0) > 0:
            k0 = 2 * np.pi / lamda0
            k1 = 2 * np.pi / lamda1

            Q0x = k0 * np.cos(np.deg2rad(nu0)) * np.sin(np.deg2rad(gamma0))
            Q1x = k1 * np.cos(np.deg2rad(nu1)) * np.sin(np.deg2rad(gamma1))

            Q0y = k0 * np.sin(np.deg2rad(nu0))
            Q1y = k1 * np.sin(np.deg2rad(nu1))

            Q0z = k0 * (
                np.cos(np.deg2rad(nu0)) * np.cos(np.deg2rad(gamma0)) - 1
            )
            Q1z = k1 * (
                np.cos(np.deg2rad(nu1)) * np.cos(np.deg2rad(gamma1)) - 1
            )

            mask = []
            for i in range(len(k1)):
                peak0 = mtd["peaks"].createPeak(V3D(Q0x[i], Q0y[i], Q0z[i]))
                peak1 = mtd["peaks"].createPeak(V3D(Q1x[i], Q1y[i], Q1z[i]))
                mask.append(
                    (peak0.getDetectorID() in self.det_ID)
                    & (peak1.getDetectorID() in self.det_ID)
                )

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
        if len(self.angles_gamma) > 0:
            d2 = (self.angles_gamma - gamma) ** 2 + (self.angles_nu - nu) ** 2

            i = np.argmin(d2)

            angles = self.angles[i]

            return angles, self.angles_gamma[i], self.angles_nu[i]

    def add_orientation(self, angles, wavelength, d_min, rows):
        if np.isclose(wavelength[0], wavelength[1]):
            wavelength = [0.975 * wavelength[0], 1.025 * wavelength[1]]

        axes = self.axes.copy()

        for i, angle in enumerate(angles):
            axes[i] = axes[i].format(angle)

        ol = mtd["coverage"].sample().getOrientedLattice()
        UB = ol.getUB().copy()

        SetUB(Workspace="instrument", UB=UB)

        SetGoniometer(
            Workspace="instrument",
            Axis0=axes[0],
            Axis1=axes[1],
            Axis2=axes[2],
            Axis3=axes[3],
            Axis4=axes[4],
            Axis5=axes[5],
        )

        d_max = 1.1 * np.max([ol.d(1, 0, 0), ol.d(0, 1, 0), ol.d(0, 0, 1)])

        ws = "peaks_orientation_{}".format(rows)

        PredictPeaks(
            InputWorkspace="instrument",
            MinDSpacing=d_min,
            MaxDSpacing=d_max,
            WavelengthMin=wavelength[0],
            WavelengthMax=wavelength[1],
            ReflectionCondition="Primitive",
            OutputWorkspace=ws,
        )

        SortPeaksWorkspace(
            InputWorkspace=ws,
            ColumnNameToSortBy="DSpacing",
            SortAscending=False,
            OutputWorkspace=ws,
        )

        columns = ["l", "k", "h"]

        for col in columns:
            SortPeaksWorkspace(
                InputWorkspace=ws,
                ColumnNameToSortBy=col,
                SortAscending=False,
                OutputWorkspace=ws,
            )

        for no in range(mtd[ws].getNumberPeaks() - 1, 0, -1):
            if (
                mtd[ws].getPeak(no).getHKL() - mtd[ws].getPeak(no - 1).getHKL()
            ).norm2() == 0:
                DeleteTableRows(TableWorkspace=ws, Rows=no)

        for no in range(mtd[ws].getNumberPeaks() - 1, 0, -1):
            if mtd[ws].getPeak(no).getDetectorID() not in self.det_ID:
                DeleteTableRows(TableWorkspace=ws, Rows=no)

        SortPeaksWorkspace(
            InputWorkspace=ws,
            ColumnNameToSortBy="DSpacing",
            SortAscending=False,
            OutputWorkspace=ws,
        )

        for peak in mtd[ws]:
            peak.setRunNumber(rows)

        mtd["instrument"].run().getGoniometer().setR(np.eye(3))

        SetUB(Workspace="combined", UB=UB)
        CombinePeaksWorkspaces(
            LHSWorkspace="combined",
            RHSWorkspace=ws,
            OutputWorkspace="combined",
        )

    def generate_table(self, row):
        if row == -1 and mtd.doesExist("missing"):
            CloneWorkspace(InputWorkspace="missing", OutputWorkspace="table")
        else:
            FilterPeaks(
                InputWorkspace="combined",
                FilterVariable="RunNumber",
                FilterValue=str(row),
                Operator="=",
                OutputWorkspace="table",
            )

        SortPeaksWorkspace(
            InputWorkspace="table",
            ColumnNameToSortBy="DSpacing",
            SortAscending=False,
            OutputWorkspace="table",
        )

        columns = ["l", "k", "h"]

        for col in columns:
            SortPeaksWorkspace(
                InputWorkspace="table",
                ColumnNameToSortBy=col,
                SortAscending=False,
                OutputWorkspace="table",
            )

        SortPeaksWorkspace(
            InputWorkspace="table",
            ColumnNameToSortBy="DSpacing",
            SortAscending=False,
            OutputWorkspace="table",
        )

        h = mtd["table"].column("h")
        k = mtd["table"].column("k")
        l = mtd["table"].column("l")

        d = mtd["table"].column("DSpacing")
        lamda = mtd["table"].column("Wavelength")

        return np.array([h, k, l, d, lamda]).T.tolist()

    def calculate_statistics(self, point_group, lattice_centering, use, d_min):
        shel, comp, mult, refl = [], [], [], []
        if mtd.doesExist("combined"):
            CloneWorkspace(
                InputWorkspace="combined", OutputWorkspace="filtered"
            )

            rows = np.arange(len(use)).tolist()

            for row in rows:
                if not use[row]:
                    FilterPeaks(
                        InputWorkspace="filtered",
                        FilterVariable="RunNumber",
                        FilterValue=str(row),
                        Operator="!=",
                        OutputWorkspace="filtered",
                    )

            if mtd["filtered"].getNumberPeaks() > 0:
                ol = mtd["combined"].sample().getOrientedLattice()
                d_max = np.max([ol.d(1, 0, 0), ol.d(0, 1, 0), ol.d(0, 0, 1)])

                d = 1 / np.sqrt(np.linspace(1 / d_max**2, 1 / d_min**2, 5))

                pg, lc = self.get_symmetry(point_group, lattice_centering)

                output = CountReflections(
                    InputWorkspace="filtered",
                    PointGroup=pg,
                    LatticeCentering=lc,
                    MinDSpacing=d_min,
                    MaxDSpacing=d_max,
                    MissingReflectionsWorkspace="missing",
                )

                unique, completeness, redundancy, multiple, _ = output

                shel = ["Overall"]
                comp = [completeness * 100]
                mult = [redundancy]
                refl = [unique]

                for i in range(len(d) - 1):
                    output = CountReflections(
                        InputWorkspace="filtered",
                        PointGroup=pg,
                        LatticeCentering=lc,
                        MinDSpacing=d[i + 1],
                        MaxDSpacing=d[i],
                        MissingReflectionsWorkspace="",
                    )

                    unique, completeness, redundancy, multiple = output

                    shel.append("{:.2f}-{:.2f}".format(d[i], d[i + 1]))
                    comp.append(completeness * 100)
                    mult.append(redundancy)
                    refl.append(unique)
            else:
                return None

        return shel, comp, mult, refl

    def hsl_to_rgb(self, hue, saturation, lightness):
        h = np.array(hue)
        s = np.array(saturation)
        l = np.array(lightness)

        def f(h, s, l, n):
            k = (n + h / 30) % 12
            a = s * np.minimum(l, 1 - l)
            return l - a * np.maximum(
                -1, np.minimum(np.minimum(k - 3, 9 - k), 1)
            )

        rgb = np.stack((f(h, s, l, 0), f(h, s, l, 8), f(h, s, l, 4)), axis=-1)

        return rgb

    def delete_angles(self, rows):
        for row in rows:
            FilterPeaks(
                InputWorkspace="combined",
                FilterVariable="RunNumber",
                FilterValue=str(row),
                Operator="!=",
                OutputWorkspace="combined",
            )

        runs = mtd["combined"].column(0)

        u, new_runs = np.unique(runs, return_index=True)

        for new_run, peak in zip(new_runs.tolist(), mtd["combined"]):
            peak.setRunNumber(new_run)

    def get_coverage_info(self, point_group, lattice_centering):
        pg = PointGroupFactory.createPointGroup(point_group)

        coverage_dict = {}

        UB = mtd["coverage"].sample().getOrientedLattice().getUB().copy()
        # UB_inv = np.linalg.inv(UB)

        if mtd.doesExist("filtered"):
            h = mtd["filtered"].column("h")
            k = mtd["filtered"].column("k")
            l = mtd["filtered"].column("l")

            hkls = np.array([h, k, l]).T.astype(int).tolist()

            hkl_dict = {}
            hkl_dict[(0, 0, 0)] = 0
            for hkl in hkls:
                if centering_conditions[lattice_centering](*hkl):
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

            r = np.sqrt(hkls[:, 0] ** 2 + hkls[:, 1] ** 2 + hkls[:, 2] ** 2)
            theta = np.arccos(hkls[:, 2] / r)
            phi = np.arctan2(hkls[:, 1], hkls[:, 0])

            hue = phi * 180 / np.pi + 180
            saturation = np.ones_like(hue)
            lightness = theta / np.pi

            rgb = self.hsl_to_rgb(hue, saturation, lightness)
            coords = np.einsum("ij,nj->ni", 2 * np.pi * UB, hkls)

            coverage_dict["colors"] = (rgb * 255).astype(np.uint8)
            coverage_dict["sizes"] = nos / nos.max()
            coverage_dict["coords"] = coords

            hkls = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])

            r = np.sqrt(hkls[:, 0] ** 2 + hkls[:, 1] ** 2 + hkls[:, 2] ** 2)
            theta = np.arccos(hkls[:, 2] / r)
            phi = np.arctan2(hkls[:, 1], hkls[:, 0])

            hue = phi * 180 / np.pi + 180
            saturation = np.ones_like(hue)
            lightness = theta / np.pi

            coverage_dict["axis_coords"] = coords

            rgb = self.hsl_to_rgb(hue, saturation, lightness)
            coords = np.einsum("ij,nj->ni", 2 * np.pi * UB, hkls)

            coverage_dict["axis_colors"] = (rgb * 255).astype(np.uint8)
            coverage_dict["axis_coords"] = coords

            return coverage_dict

    def crystal_plan(self, *args):
        return CrystalPlan(*args)


class CrystalPlan:
    def __init__(
        self,
        use,
        opt,
        axes,
        limits,
        wavelength,
        d_min,
        point_group,
        lattice_centering,
    ):
        CloneWorkspace(
            InputWorkspace="combined", OutputWorkspace="crystal_plan"
        )

        self.instrument = "instrument"

        rows = np.arange(len(use)).tolist()

        for row in rows:
            if not use[row] or opt[row]:
                FilterPeaks(
                    InputWorkspace="crystal_plan",
                    FilterVariable="RunNumber",
                    FilterValue=str(row),
                    Operator="!=",
                    OutputWorkspace="crystal_plan",
                )

        if np.isclose(wavelength[0], wavelength[1]):
            wavelength = [0.975 * wavelength[0], 1.025 * wavelength[1]]

        self.wavelength = wavelength

        self.axes = axes.copy()
        self.limits = limits.copy()

        ol = mtd["coverage"].sample().getOrientedLattice()
        UB = ol.getUB().copy()

        SetUB(Workspace="instrument", UB=UB)
        SetUB(Workspace="crystal_plan", UB=UB)

        self.UB = UB.copy()
        self.d_min = d_min
        self.d_max = 1.1 * np.max(
            [ol.d(1, 0, 0), ol.d(0, 1, 0), ol.d(0, 0, 1)]
        )
        self.offset = len(use)

        self.point_group = point_group
        self.lattice_centering = lattice_centering
        self.genes = {}

        # rng seed ---------#
        np.random.seed(13)  #
        #####################

    def generation(self, i, j):
        axes = self.axes.copy()
        limits = self.limits.copy()

        ax = [None] * 6
        angles = []
        for ind, (axis, limit) in enumerate(zip(axes, limits)):
            delta = limit[1] - limit[0]
            angle = limit[0] + delta * np.random.random()
            ax[ind] = axis.format(angle)
            if not np.isclose(delta, 0):
                angles.append(angle)

        outname = "peaks_{}_{}".format(i, j)

        self.genes[outname] = angles

        SetGoniometer(
            Workspace=self.instrument,
            Axis0=ax[0],
            Axis1=ax[1],
            Axis2=ax[2],
            Axis3=ax[3],
            Axis4=ax[4],
            Axis5=ax[5],
        )

        PredictPeaks(
            InputWorkspace=self.instrument,
            CalculateStructureFactors=False,
            MinDSpacing=self.d_min,
            MaxDSpacing=self.d_max,
            WavelengthMin=self.wavelength[0],
            WavelengthMax=self.wavelength[1],
            ReflectionCondition="Primitive",
            OutputWorkspace=outname,
        )

        for pk in mtd[outname]:
            pk.setRunNumber(i + self.offset)
            pk.setIntensity(100)
            pk.setSigmaIntensity(10)

    def initialization(self, n_orient, n_indiv):
        fit = []
        for j in range(n_indiv):
            for i in range(n_orient):
                self.generation(i, j)
            self.recombination(n_orient, j)
            fit.append(self.fitness("peaks_{}".format(j)))

        return np.array(fit)

    def recombination(self, n_orient, j):
        individuals = "peaks_{}".format(j)
        for i in range(n_orient):
            genes = "peaks_{}_{}".format(i, j)
            if i == 0:
                CombinePeaksWorkspaces(
                    LHSWorkspace="crystal_plan",
                    RHSWorkspace=genes,
                    OutputWorkspace=individuals,
                )
            else:
                CombinePeaksWorkspaces(
                    LHSWorkspace=individuals,
                    RHSWorkspace=genes,
                    OutputWorkspace=individuals,
                )

    def fitness(self, peaks):
        output = CountReflections(
            InputWorkspace=peaks,
            PointGroup=self.point_group,
            LatticeCentering=self.lattice_centering,
            MinDSpacing=self.d_min,
            MaxDSpacing=self.d_max,
            MissingReflectionsWorkspace="",
        )

        unique, completeness, redundancy, multiple = output

        return completeness * 100

    def crossover(self, n_orient, n_elite, best, selection):
        j = 0

        genes = "peaks_{}_{}"
        genome = "s_{}_{}"

        workspaces = []

        for elite in best:
            for i in range(n_orient):
                CloneWorkspace(
                    InputWorkspace=genes.format(i, elite),
                    OutputWorkspace=genome.format(i, j),
                )

                workspaces.append(genome.format(i, j))

            j += 1

        for parents in selection:
            k = np.random.randint(1, n_orient)

            for i in range(k):
                CloneWorkspace(
                    InputWorkspace=genes.format(i, parents[0]),
                    OutputWorkspace=genome.format(i, j + 0),
                )

                CloneWorkspace(
                    InputWorkspace=genes.format(i, parents[1]),
                    OutputWorkspace=genome.format(i, j + 1),
                )

                workspaces.append(genome.format(i, j + 0))
                workspaces.append(genome.format(i, j + 1))

            for i in range(k, n_orient):
                CloneWorkspace(
                    InputWorkspace=genes.format(i, parents[1]),
                    OutputWorkspace=genome.format(i, j + 0),
                )

                CloneWorkspace(
                    InputWorkspace=genes.format(i, parents[0]),
                    OutputWorkspace=genome.format(i, j + 1),
                )

                workspaces.append(genome.format(i, j + 0))
                workspaces.append(genome.format(i, j + 1))

            j += 2

        RenameWorkspaces(InputWorkspaces=workspaces, Prefix="peak")

    def mutation(self, n_orient, n_indiv, mutation_rate):
        fit = []
        for j in range(n_indiv):
            for i in range(n_orient):
                if np.random.random() < mutation_rate:
                    self.generation(i, j)
            self.recombination(n_orient, j)
            fit.append(self.fitness("peaks_{}".format(j)))

        return np.array(fit)

    def optimize(self, n_orient, n_indiv, n_gener, n_elite, mutation_rate):
        fit = self.initialization(n_orient, n_indiv)

        ranking = np.argsort(fit)

        for _ in range(n_gener):
            ranking = np.argsort(fit)

            best = ranking[-n_elite:]

            fraction = fit / np.sum(fit)

            selection = []

            while len(selection) < (n_indiv - n_elite) // 2:
                selection.append(
                    np.random.choice(
                        np.arange(n_indiv), size=2, p=fraction, replace=False
                    )
                )

            self.crossover(n_orient, n_elite, best, selection)

            fit = self.mutation(n_orient, n_indiv, mutation_rate)

        ranking = np.argsort(fit)

        j = ranking[-1]

        values = []
        for i in range(n_orient):
            genes = "peaks_{}_{}".format(i, j)
            values.append(self.genes[genes])

        CloneWorkspace(
            InputWorkspace="peaks_{}".format(j), OutputWorkspace="combined"
        )

        return values
