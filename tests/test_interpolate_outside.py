"""Regression tests for `Waveform.Interpolate` with times outside the domain"""

import numpy as np
import pytest

import GWFrames
from GWFrames import Quaternions


def smooth_waveform():
    """A waveform with smooth modes and a time-dependent frame"""
    T = np.linspace(0.0, 10.0, 201)
    LM = [[2, m] for m in range(-2, 3)]
    Data = np.exp(1j * np.outer(np.arange(1, 6), 0.3 * T)) * (1.0 + 0.1 * T)
    W = GWFrames.Waveform(T, LM, Data)
    W.SetFrame(Quaternions.QuaternionArray(0.1 * T, 0.05 * T, 0.02 * T))
    W.SetFrameType(GWFrames.Corotating)
    return W


def test_new_times_before_and_after_domain():
    W = smooth_waveform()
    Inside = np.linspace(1.05, 8.95, 40)
    NewT = np.concatenate(([-3.0, -2.0, -1.0], Inside, [11.0, 12.0]))
    I = W.Interpolate(NewT, True)
    np.testing.assert_array_equal(I.T(), NewT)
    # Inside the domain, the result agrees with ordinary interpolation
    Expected = W.Interpolate(Inside)
    np.testing.assert_allclose(I.Data()[:, 3:-2], Expected.Data(), rtol=0, atol=1e-14)
    # Outside the domain, the modes are zero
    assert np.all(I.Data()[:, :3] == 0)
    assert np.all(I.Data()[:, -2:] == 0)
    # The frame matches ordinary interpolation inside the domain, and is held at
    # its value at the nearer end of the domain outside it
    F = I.Frame()
    F_inside = Expected.Frame()
    F0, F1 = W.Frame()[0], W.Frame()[-1]
    for i in range(len(Inside)):
        assert (F[3 + i] - F_inside[i]).abs() < 1e-14
    for i in range(3):
        assert (F[i] - F0).abs() == 0
    for i in range(1, 3):
        assert (F[-i] - F1).abs() == 0


@pytest.mark.parametrize("NewT", [[-3.0, -2.0, -1.0], [11.0, 12.0, 13.0], [-1.0, 11.0]])
def test_all_new_times_outside_domain(NewT):
    W = smooth_waveform()
    I = W.Interpolate(NewT, True)
    np.testing.assert_array_equal(I.T(), NewT)
    assert I.Data().shape == (W.NModes(), len(NewT))
    assert np.all(I.Data() == 0)
    F = I.Frame()
    assert len(F) == len(NewT)
    for t, f in zip(NewT, F):
        expected = W.Frame()[0] if t < W.T()[0] else W.Frame()[-1]
        assert (f - expected).abs() == 0


def test_metadata_copied():
    W = smooth_waveform()
    W.SetBoostWeight(3)
    W.SetVersionHist([("abc123", "first commit"), ("def456", "second commit")])
    for I in [W.Interpolate(np.linspace(1, 9, 17)), W.Interpolate(np.linspace(-1, 11, 17), True)]:
        assert I.BoostWeight() == 3
        assert I.VersionHist() == W.VersionHist()
        assert I.SpinWeight() == W.SpinWeight()
        assert I.FrameType() == W.FrameType()
