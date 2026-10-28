"""Regression tests: `Waveform.Compare` honors `MinTimeStep` and `MinTime`."""

import numpy as np

import GWFrames


def make_waveform(T, phase=0.0):
    LM = [[2, m] for m in range(-2, 3)]
    Data = np.exp(1j * (np.outer(np.arange(1, 6), T) + phase))
    W = GWFrames.Waveform(T, LM, Data)
    W.SetFrameType(GWFrames.Inertial)
    return W


def test_compare_default_time_step():
    T = np.linspace(0.0, 1.0, 1001)  # step 0.001, below the default 0.005
    A = make_waveform(T)
    B = make_waveform(T, phase=0.1)
    C = B.Compare(A)
    np.testing.assert_allclose(np.diff(C.T()), 0.005, rtol=1e-8)


def test_compare_min_time_step():
    T = np.linspace(0.0, 1.0, 1001)
    A = make_waveform(T)
    B = make_waveform(T, phase=0.1)
    C = B.Compare(A, MinTimeStep=0.0005)
    np.testing.assert_allclose(np.diff(C.T()), 0.001, rtol=1e-6)
    np.testing.assert_allclose(C.T(), T, atol=1e-12)


def test_compare_min_time():
    T = np.linspace(0.0, 10.0, 101)
    A = make_waveform(T)
    B = make_waveform(T, phase=0.1)
    C = B.Compare(A, MinTime=5.0)
    assert C.T()[0] == 5.0
    assert C.T()[-1] <= 10.0
    # The data are still the difference of the two waveforms
    i = np.searchsorted(T, 7.0)
    j = np.searchsorted(C.T(), 7.0)
    np.testing.assert_allclose(C.Data()[:, j], A.Data()[:, i] - B.Data()[:, i], atol=1e-6)
