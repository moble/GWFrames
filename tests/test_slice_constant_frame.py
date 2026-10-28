"""Regression tests for the treatment of frame data by the slicing functions

Every slicing function should slice frame data with one rotor per time step,
copy a constant frame (a single rotor) unchanged, and leave an empty frame
empty.  `SliceOfTimeIndicesWithoutModes` and `SliceOfTimesWithoutModes` used
to drop a constant frame.
"""

import numpy as np
import pytest

import GWFrames
from GWFrames import Quaternions


def components(frame):
    return np.array([[R[0], R[1], R[2], R[3]] for R in frame])


def make_waveform(nframe):
    T = np.linspace(0, 10, 11)
    LM = [[2, m] for m in range(-2, 3)]
    Data = np.exp(1j * np.outer(np.arange(1, 6), T))
    W = GWFrames.Waveform(T, LM, Data)
    frame = [Quaternions.Quaternion(1.0, 0.1 * i, 0.2, 0.3).normalized() for i in range(nframe)]
    W.SetFrame(frame)
    return W, frame


index_slicers = ["SliceOfTimeIndices", "SliceOfTimeIndicesWithEll2", "SliceOfTimeIndicesWithoutModes"]
time_slicers = ["SliceOfTimes", "SliceOfTimesWithEll2", "SliceOfTimesWithoutModes"]


@pytest.mark.parametrize("slicer", index_slicers)
def test_index_slice_keeps_constant_frame(slicer):
    W, frame = make_waveform(1)
    S = getattr(W, slicer)(2, 5)
    np.testing.assert_array_equal(S.T(), W.T()[2:5])
    np.testing.assert_array_equal(components(S.Frame()), components(frame))


@pytest.mark.parametrize("slicer", time_slicers)
def test_time_slice_keeps_constant_frame(slicer):
    W, frame = make_waveform(1)
    S = getattr(W, slicer)(2.0, 5.0)
    assert S.NTimes() > 1
    np.testing.assert_array_equal(components(S.Frame()), components(frame))


@pytest.mark.parametrize("slicer", index_slicers)
def test_index_slice_of_time_dependent_frame(slicer):
    W, frame = make_waveform(11)
    S = getattr(W, slicer)(2, 5)
    np.testing.assert_array_equal(components(S.Frame()), components(frame[2:5]))


@pytest.mark.parametrize("slicer", index_slicers)
def test_index_slice_of_empty_frame(slicer):
    W, _ = make_waveform(0)
    S = getattr(W, slicer)(2, 5)
    assert len(S.Frame()) == 0


@pytest.mark.parametrize("slicer", index_slicers)
def test_index_slice_rejects_bad_frame_size(slicer):
    W, _ = make_waveform(3)
    with pytest.raises(ValueError):
        getattr(W, slicer)(2, 5)


@pytest.mark.parametrize("slicer", index_slicers)
def test_index_slice_rejects_impossible_indices(slicer):
    W, _ = make_waveform(1)
    with pytest.raises(ValueError):
        getattr(W, slicer)(5, 2)
    with pytest.raises(IndexError):
        getattr(W, slicer)(2, 12)
