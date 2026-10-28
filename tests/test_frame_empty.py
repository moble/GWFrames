"""Regression tests for `Waveform.Frame(TimeIndex)` when the frame is empty

An empty frame means that the frame is constant and equal to the identity, which
is the case for any Waveform constructed directly from times, modes, and data.
`Frame(TimeIndex)` used to read the first element of the empty frame, which
crashed the interpreter.
"""

import numpy as np

import GWFrames
from GWFrames import Quaternions


def components(q):
    return [q[0], q[1], q[2], q[3]]


def small_waveform():
    t = [0.0, 1.0, 2.0, 3.0]
    lm = [[2, -2], [2, -1], [2, 0], [2, 1], [2, 2]]
    data = [[complex(1 + i_m, i_t) for i_t in range(len(t))] for i_m in range(len(lm))]
    return GWFrames.Waveform(t, lm, data)


def test_empty_frame_is_identity():
    W = small_waveform()
    assert len(W.Frame()) == 0
    for i_t in range(W.NTimes()):
        assert components(W.Frame(i_t)) == [1.0, 0.0, 0.0, 0.0]


def test_constant_and_time_dependent_frames_unchanged():
    W = small_waveform()
    R = Quaternions.Quaternion(0.1, 0.2, 0.3, 0.4).normalized()
    W.SetFrame([R])
    for i_t in range(W.NTimes()):
        np.testing.assert_array_equal(components(W.Frame(i_t)), components(R))
    Rs = [Quaternions.exp(Quaternions.Quaternion(0.0, 0.0, 0.0, 0.1 * i_t)) for i_t in range(W.NTimes())]
    W.SetFrame(Rs)
    for i_t in range(W.NTimes()):
        np.testing.assert_array_equal(components(W.Frame(i_t)), components(Rs[i_t]))
