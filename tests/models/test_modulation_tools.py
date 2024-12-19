import os

from mantid.simpleapi import mtd

from NeuXtalViz.models.modulation_tools import ModulationModel

peaks_file = os.path.join("tests/data", "26079_Niggli.integrate")
ub_file = os.path.join("tests/data", "26079_Niggli.mat")


def test_load_peaks():
    mod = ModulationModel()

    assert mtd["peaks"].getNumberPeaks() == 0

    mod.load_peaks(peaks_file)

    assert mtd["peaks"].getNumberPeaks() > 0


def test_load_UB():
    mod = ModulationModel()

    assert mod.UB == None

    mod.load_UB(ub_file)

    assert mod.UB.shape == (3, 3)


def test_cluster():
    mod = ModulationModel()

    mod.load_peaks(peaks_file)
    mod.load_UB(ub_file)

    peak_info = mod.get_peak_info()

    assert len(peak_info["coordinates"]) == len(peak_info["points"])
    assert len(peak_info["numbers"]) == len(peak_info["points"])

    mod.cluster_peaks(peak_info, eps=0.025, min_samples=15)

    assert peak_info["nuclear"].shape == (3,)
    assert peak_info["satellites"].shape == (3, 3)
