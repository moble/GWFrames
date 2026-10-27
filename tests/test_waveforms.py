import copy
import pickle

import numpy as np
import pytest

import GWFrames
from GWFrames import Quaternions


@pytest.fixture(scope="module")
def pn():
    """A short precessing PN waveform, in the coorbital frame as constructed"""
    return GWFrames.PNWaveform("TaylorT4", 0.1, [0.3, 0.2, 0.4], [-0.2, 0.1, 0.3], 0.05)


@pytest.fixture
def inertial(pn):
    return GWFrames.Waveform(pn).TransformToInertialFrame()


def test_pn_waveform(pn):
    assert isinstance(pn, GWFrames.PNWaveform)
    assert isinstance(pn, GWFrames.Waveform)
    assert pn.NTimes() > 100
    assert pn.NModes() == 77
    assert pn.FrameType() == GWFrames.Coorbital
    assert pn.DataType() == GWFrames.h
    assert np.all(np.diff(pn.T()) > 0)
    assert pn.LM().shape == (77, 2)
    assert pn.Data().shape == (77, pn.NTimes())
    assert pn.chi1().shape == (pn.NTimes(), 3)
    assert pn.Omega_orbMag().shape == (pn.NTimes(),)
    assert len(pn.Frame()) == pn.NTimes()
    # The spin magnitudes are conserved by the evolution
    np.testing.assert_allclose(np.linalg.norm(pn.chi1(), axis=1), np.linalg.norm([0.3, 0.2, 0.4]), rtol=1e-6)


def test_methods_return_self(pn):
    W = GWFrames.PNWaveform(pn)
    assert W.TransformToInertialFrame() is W
    assert W.SetFrameType(GWFrames.Inertial) is W
    assert W.SetT(W.T()).SetLM(W.LM()) is W
    assert isinstance(W.SliceOfTimes(0, 100), GWFrames.Waveform)


def test_construction_from_arrays():
    T = np.linspace(0, 10, 11)
    LM = [[2, m] for m in range(-2, 3)]
    Data = np.exp(1j * np.outer(np.arange(1, 6), T))
    W = GWFrames.Waveform(T, LM, Data)
    np.testing.assert_array_equal(W.T(), T)
    np.testing.assert_array_equal(W.LM(), LM)
    np.testing.assert_array_equal(W.Data(), Data)
    np.testing.assert_array_equal(W.Data(2), Data[2])
    assert W.Data(2, 3) == Data[2, 3]
    assert W.FindModeIndex(2, 1) == 3
    with pytest.raises(IndexError):
        W.FindModeIndex(3, 0)


def test_frame_round_trip(inertial):
    W = GWFrames.Waveform(inertial).TransformToCorotatingFrame()
    assert W.FrameType() == GWFrames.Corotating
    W.TransformToInertialFrame()
    np.testing.assert_allclose(W.Data(), inertial.Data(), atol=1e-12 * np.max(np.abs(inertial.Data())))


def test_rotation_round_trip(inertial):
    R = Quaternions.Quaternion(0.1, 0.2, 0.3, 0.4).normalized()
    W = GWFrames.Waveform(inertial).RotatePhysicalSystem(R).RotatePhysicalSystem(R.inverse())
    np.testing.assert_allclose(W.Data(), inertial.Data(), atol=1e-12 * np.max(np.abs(inertial.Data())))


def test_arithmetic(inertial):
    W = inertial
    np.testing.assert_allclose((W + W).Data(), 2 * W.Data())
    np.testing.assert_allclose((2.0 * W).Data(), 2 * W.Data())
    np.testing.assert_allclose((W / 2.0).Data(), W.Data() / 2)
    np.testing.assert_allclose((W - W).Data(), 0)


def test_pickle_and_copy(inertial):
    W = inertial
    for W2 in (pickle.loads(pickle.dumps(W)), copy.copy(W), copy.deepcopy(W)):
        assert W2 is not W
        np.testing.assert_array_equal(W2.T(), W.T())
        np.testing.assert_array_equal(W2.LM(), W.LM())
        np.testing.assert_array_equal(W2.Data(), W.Data())
        assert W2.HistoryStr() == W.HistoryStr()
        assert W2.FrameType() == W.FrameType()
    assert isinstance(copy.copy(GWFrames.PNWaveform()), GWFrames.PNWaveform)


@pytest.mark.parametrize("writer, reader", [("OutputToH5", "ReadFromH5"), ("OutputToNRAR", "ReadFromNRAR")])
def test_h5_round_trip(tmp_path, pn, writer, reader):
    W = GWFrames.Waveform(pn).KeepOnlyEll2()
    getattr(W, writer)(str(tmp_path / "W.h5"))
    W2 = getattr(GWFrames, reader)(str(tmp_path / (W.GetFileNamePrefix() + "W.h5")))
    np.testing.assert_array_equal(W2.T(), W.T())
    np.testing.assert_array_equal(W2.LM(), W.LM())
    np.testing.assert_array_equal(W2.Data(), W.Data())
    assert W2.FrameType() == W.FrameType()
    assert W2.DataType() == W.DataType()
    if reader == "ReadFromH5":
        assert len(W2.Frame()) == len(W.Frame())
        assert all(a == b for a, b in zip(W2.Frame(), W.Frame()))
    assert W.HistoryStr().replace("\n", "\n#") in W2.HistoryStr()


def test_hybridize(inertial):
    W = GWFrames.Waveform(inertial).KeepOnlyEll2()
    T = W.T()
    t1, t2 = T[len(T) // 4], T[len(T) // 2]
    H = W.Hybridize(W, t1, t2)
    # Hybridizing a waveform with itself should change nothing but the time grid
    # (away from the end, where the PN data vary too rapidly for interpolation)
    T = T[(T >= H.T()[0]) & (T <= T[3 * len(T) // 4])]
    np.testing.assert_allclose(H.Interpolate(T).Data(), W.Interpolate(T).Data(), atol=1e-8)
    assert len(W.Frame()) == len(H.Frame()) == 0  # Inertial frames are empty


def test_latex_description(pn):
    assert pn.GetLaTeXDataDescription() == r"r\,h/M"
