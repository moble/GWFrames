"""Generate reference results from the old SWIG build of GWFrames.

This script records results computed by the last SWIG-based build of GWFrames,
so that `tests/test_swig_reference.py` can check that the pybind11 bindings give
the same answers.  It was run with the conda-forge package `gwframes 1.0` (for
Python 2.7 on Intel macOS, under Rosetta):

    CONDA_SUBDIR=osx-64 mamba create -p ./gwf_old -c conda-forge gwframes python=2.7
    ./gwf_old/bin/python tests/reference/generate_swig_reference.py tests/data/swig_reference.npz

The input waveform is generated here and saved along with the results, so the
comparison does not depend on post-Newtonian evolution, whose tolerances and
3.5PN terms have changed since that build.

"""

from __future__ import division, print_function

import sys

import numpy as np

import GWFrames
import Quaternions


def components(Q):
    return np.array([[q[0], q[1], q[2], q[3]] for q in Q])


def main(output):
    results = {}

    # A precessing PN waveform in the inertial frame, on a uniform time grid
    W = GWFrames.PNWaveform("TaylorT4", 0.2, [0.3, 0.2, 0.4], [-0.2, 0.1, 0.3], 0.03)
    W.TransformToInertialFrame()
    W.KeepOnlyEllModes([2, 3, 4])
    T = np.linspace(W.T()[0], W.T()[0] + 0.8 * (W.T()[-1] - W.T()[0]), 400)
    W0 = W.Interpolate(T)
    results["input_T"] = W0.T()
    results["input_LM"] = np.array(W0.LM())
    results["input_Data"] = W0.Data()
    results["input_Frame"] = components(W0.Frame())
    results["input_FrameType"] = W0.FrameType()
    results["input_DataType"] = W0.DataType()

    def copy():
        return GWFrames.Waveform(W0)

    results["AngularVelocityVector"] = W0.AngularVelocityVector()
    results["LLDominantEigenvector"] = W0.LLDominantEigenvector()
    R_corot = W0.CorotatingFrame()
    results["CorotatingFrame"] = components(R_corot)
    results["Norm"] = W0.Norm()
    results["EvaluateAtPoint"] = W0.EvaluateAtPoint(0.3, 0.4)
    results["XParityViolationNormalized"] = W0.XParityViolationNormalized()
    results["ParityViolationNormalized"] = W0.ParityViolationNormalized()

    Wc = copy()
    Wc.TransformToCoprecessingFrame()
    results["Coprecessing_Data"] = Wc.Data()
    results["Coprecessing_Frame"] = components(Wc.Frame())

    Wr = copy()
    Wr.TransformToCorotatingFrame()
    results["Corotating_Data"] = Wr.Data()
    results["Corotating_Frame"] = components(Wr.Frame())

    Wa = GWFrames.Waveform(Wr)
    Wa.AlignDecompositionFrameToModes(T[len(T) // 2])
    results["Aligned_Frame"] = components(Wa.Frame())

    Wi = GWFrames.Waveform(Wr)
    Wi.TransformToInertialFrame()
    results["CorotatingToInertial_Data"] = Wi.Data()

    R = Quaternions.Quaternion(0.1, 0.2, 0.3, 0.4).normalized()
    Wb = copy()
    Wb.RotateDecompositionBasis(R)
    results["RotateDecompositionBasis_Data"] = Wb.Data()
    Wp = copy()
    Wp.RotatePhysicalSystem(R)
    results["RotatePhysicalSystem_Data"] = Wp.Data()

    T2 = np.linspace(T[10], T[-10], 300)
    results["Interpolate_T"] = T2
    results["Interpolate_Data"] = W0.Interpolate(T2).Data()

    Wd = copy()
    Wd.Differentiate()
    results["Differentiate_Data"] = Wd.Data()

    results["MinimalRotation"] = components(Quaternions.MinimalRotation(R_corot, T))
    results["FrameAngularVelocity"] = components(Quaternions.FrameAngularVelocity(R_corot, T))
    results["Squad"] = components(Quaternions.Squad(R_corot, T, T2))
    results["ApproximateMeanRotor"] = components([Quaternions.ApproximateMeanRotor(R_corot)])

    FT = GWFrames.WaveformAtAPointFT(W0, 0.5, 0.3, 0.4, 20.0)
    results["FT_F"] = np.array(FT.F())
    results["FT_Re"] = np.array(FT.Re())
    results["FT_Im"] = np.array(FT.Im())
    results["FT_SNR"] = FT.SNR()

    ellMax = 4
    rng = np.random.RandomState(1234)
    modes = rng.normal(size=(ellMax + 1) ** 2) + 1j * rng.normal(size=(ellMax + 1) ** 2)
    modes[:4] = 0.0  # Modes with ell < |s| are zero for spin weight -2
    M = GWFrames.Modes(-2, modes)
    M.SetEllMax(ellMax)
    results["Modes_input"] = modes
    results["Modes_edth"] = np.array(M.edth().Data())
    results["Modes_edthbar"] = np.array(M.edthbar().Data())
    results["Modes_bar"] = np.array(M.bar().Data())

    np.savez_compressed(output, **results)


if __name__ == "__main__":
    main(sys.argv[1])
