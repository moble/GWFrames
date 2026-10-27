"""Compare the pybind11 bindings with results from the old SWIG build.

The reference data were produced by `tests/reference/generate_swig_reference.py`;
see that file for details.  Differences come only from compilers, libraries,
and the old build's use of `-ffast-math`, so the tolerances are tight.

"""

from pathlib import Path

import numpy as np
import pytest

import GWFrames
from GWFrames import Quaternions

RTOL = 1e-8
ATOL = 1e-10


@pytest.fixture(scope="module")
def ref():
    with np.load(Path(__file__).parent / "data" / "swig_reference.npz") as data:
        return dict(data)


@pytest.fixture(scope="module")
def W0(ref):
    W = GWFrames.Waveform(ref["input_T"], ref["input_LM"], ref["input_Data"])
    W.SetFrame(ref["input_Frame"])
    W.SetFrameType(int(ref["input_FrameType"]))
    W.SetDataType(int(ref["input_DataType"]))
    return W


def components(Q):
    return np.array([[q[0], q[1], q[2], q[3]] for q in Q]).reshape(-1, 4)


def assert_close(actual, expected, rtol=RTOL, atol=ATOL):
    np.testing.assert_allclose(actual, expected, rtol=rtol, atol=atol)


def assert_rotors_close(actual, expected, atol=1e-8):
    # Rotors R and -R are equivalent, so allow either sign at each time
    actual, expected = components(actual), np.asarray(expected)
    signs = np.sign(np.sum(actual * expected, axis=1))[:, np.newaxis]
    assert_close(signs * actual, expected, atol=atol)


def test_vectors(W0, ref):
    assert_close(W0.AngularVelocityVector(), ref["AngularVelocityVector"])
    assert_close(W0.LLDominantEigenvector(), ref["LLDominantEigenvector"])
    assert_close(W0.Norm(), ref["Norm"])
    assert_close(W0.EvaluateAtPoint(0.3, 0.4), ref["EvaluateAtPoint"])
    assert_close(W0.XParityViolationNormalized(), ref["XParityViolationNormalized"])
    assert_close(W0.ParityViolationNormalized(), ref["ParityViolationNormalized"])
    assert_rotors_close(W0.CorotatingFrame(), ref["CorotatingFrame"])


def test_frame_transformations(W0, ref):
    Wc = GWFrames.Waveform(W0).TransformToCoprecessingFrame()
    assert_close(Wc.Data(), ref["Coprecessing_Data"])
    assert_rotors_close(Wc.Frame(), ref["Coprecessing_Frame"])

    Wr = GWFrames.Waveform(W0).TransformToCorotatingFrame()
    assert_close(Wr.Data(), ref["Corotating_Data"])
    assert_rotors_close(Wr.Frame(), ref["Corotating_Frame"])

    Wa = GWFrames.Waveform(Wr).AlignDecompositionFrameToModes(W0.T()[W0.NTimes() // 2])
    assert_rotors_close(Wa.Frame(), ref["Aligned_Frame"])

    Wi = GWFrames.Waveform(Wr).TransformToInertialFrame()
    assert_close(Wi.Data(), ref["CorotatingToInertial_Data"])


def test_rotations(W0, ref):
    R = Quaternions.Quaternion(0.1, 0.2, 0.3, 0.4).normalized()
    assert_close(GWFrames.Waveform(W0).RotateDecompositionBasis(R).Data(), ref["RotateDecompositionBasis_Data"])
    assert_close(GWFrames.Waveform(W0).RotatePhysicalSystem(R).Data(), ref["RotatePhysicalSystem_Data"])


def test_interpolation_and_differentiation(W0, ref):
    assert_close(W0.Interpolate(ref["Interpolate_T"]).Data(), ref["Interpolate_Data"])
    assert_close(GWFrames.Waveform(W0).Differentiate().Data(), ref["Differentiate_Data"])


def test_quaternion_arrays(W0, ref):
    T = W0.T()
    R = W0.CorotatingFrame()
    assert_rotors_close(Quaternions.MinimalRotation(R, T), ref["MinimalRotation"])
    assert_close(components(Quaternions.FrameAngularVelocity(R, T)), ref["FrameAngularVelocity"])
    assert_rotors_close(Quaternions.Squad(R, T, ref["Interpolate_T"]), ref["Squad"])
    assert_rotors_close([Quaternions.ApproximateMeanRotor(R)], ref["ApproximateMeanRotor"])


def test_fourier_transform(W0, ref):
    FT = GWFrames.WaveformAtAPointFT(W0, 0.5, 0.3, 0.4, 20.0)
    assert_close(FT.F(), ref["FT_F"])
    scale = np.max(np.abs(ref["FT_Re"]))
    assert_close(FT.Re(), ref["FT_Re"], atol=1e-10 * scale)
    assert_close(FT.Im(), ref["FT_Im"], atol=1e-10 * scale)
    assert_close(FT.SNR(), ref["FT_SNR"])


def test_scri_modes(ref):
    M = GWFrames.Modes(-2, ref["Modes_input"])
    M.SetEllMax(4)
    assert_close(M.edth().Data(), ref["Modes_edth"])
    assert_close(M.edthbar().Data(), ref["Modes_edthbar"])
    assert_close(M.bar().Data(), ref["Modes_bar"])
