"""Regression tests for Waveform.DropTimesOutside, which used to drop the first time step"""

import numpy as np
import pytest

import GWFrames
from GWFrames import Quaternions


T = np.linspace(0.0, 10.0, 11)
LM = [[2, m] for m in range(-2, 3)]
DATA = np.exp(1j * np.outer(np.arange(1, 6), T))


def waveform(frame=None):
    W = GWFrames.Waveform(T, LM, DATA)
    if frame is not None:
        W.SetFrame(frame)
    return W


def frame_rotors(n):
    return [Quaternions.Quaternion(np.cos(0.01 * i), np.sin(0.01 * i), 0.0, 0.0) for i in range(n)]


@pytest.mark.parametrize(
    "ta, tb, keep",
    [
        (-1.0, 20.0, slice(0, 11)),  # interval contains every time step
        (-1.0, 4.5, slice(0, 5)),  # ta precedes the first time step
        (0.0, 10.0, slice(0, 11)),  # endpoints equal to time steps are kept
        (2.0, 5.0, slice(2, 6)),
        (2.5, 5.5, slice(3, 6)),
        (3.0, 3.0, slice(3, 4)),
        (3.2, 3.8, slice(3, 3)),  # no time step in the interval
        (5.0, 2.0, slice(5, 5)),  # empty interval
        (-5.0, -1.0, slice(0, 0)),  # entirely before the data
        (11.0, 12.0, slice(11, 11)),  # entirely after the data
    ],
)
def test_drop_times_outside_keeps_closed_interval(ta, tb, keep):
    W = waveform()
    assert W.DropTimesOutside(ta, tb) is W
    np.testing.assert_array_equal(W.T(), T[keep])
    assert W.NTimes() == len(T[keep])
    assert W.NModes() == len(LM)
    if W.NTimes() > 0:
        np.testing.assert_array_equal(W.Data(), DATA[:, keep])


def test_drop_times_outside_trims_full_frame():
    rotors = frame_rotors(len(T))
    W = waveform(rotors).DropTimesOutside(0.0, 4.0)
    np.testing.assert_array_equal(W.T(), T[:5])
    frame = W.Frame()
    assert len(frame) == 5
    for got, expected in zip(frame, rotors[:5]):
        assert got == expected


def test_drop_times_outside_keeps_constant_frame():
    rotor = frame_rotors(2)[1]
    W = waveform([rotor]).DropTimesOutside(0.0, 4.0)
    np.testing.assert_array_equal(W.T(), T[:5])
    frame = W.Frame()
    assert len(frame) == 1
    assert frame[0] == rotor
    W = waveform().DropTimesOutside(0.0, 4.0)
    assert len(W.Frame()) == 0


def test_drop_times_outside_rejects_mismatched_frame():
    W = waveform(frame_rotors(3))
    with pytest.raises(Exception):
        W.DropTimesOutside(0.0, 4.0)
