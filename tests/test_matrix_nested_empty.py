"""Regression tests for converting empty nested sequences to and from C++

Before the fix, the caster for vectors of vectors located each row with
numpy's bounds-checked indexing, `array.data(i, 0)`, which raises an
IndexError whenever the rows have length zero.  So `Matrix([[]])`, and any
other argument or return value with at least one row and no columns, raised
"index 0 is out of bounds for axis 1 with size 0" instead of converting.
"""

import numpy as np
import pytest

import GWFrames


@pytest.mark.parametrize(
    "data",
    [
        [],
        [[]],
        [[], []],
        ((), (), ()),
        np.zeros((0, 0)),
        np.zeros((0, 3)),
        np.zeros((3, 0)),
        np.zeros((3, 0), dtype=int),
    ],
    ids=["[]", "[[]]", "[[],[]]", "tuples", "(0,0)", "(0,3)", "(3,0)", "int(3,0)"],
)
def test_matrix_from_empty_data(data):
    M = GWFrames.Matrix(data)
    assert (M.nrows(), M.ncols()) == (0, 0)
    assert M.numpy().shape == (0, 0)


@pytest.mark.parametrize(
    "data",
    [
        [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
        np.arange(6.0).reshape(2, 3),
        np.arange(12.0).reshape(3, 4)[:, ::2],  # not contiguous
        np.asfortranarray(np.arange(6.0).reshape(2, 3)),
    ],
)
def test_matrix_from_nonempty_data(data):
    M = GWFrames.Matrix(data)
    expected = np.asarray(data, dtype=float)
    assert (M.nrows(), M.ncols()) == expected.shape
    assert np.array_equal(M.numpy(), expected)


def test_ragged_input():
    # Ragged data are still accepted as a sequence of rows of different
    # lengths, and are returned as a list of one-dimensional arrays.
    W = GWFrames.Waveform([], [[2, 2], [2]], [])
    LM = W.LM()
    assert isinstance(LM, list)
    assert [row.tolist() for row in LM] == [[2, 2], [2]]
    # A ragged sequence whose rows are not sequences of numbers still fails to
    # convert, so that no overload matches.
    with pytest.raises(TypeError):
        GWFrames.Matrix([[1.0, 2.0], "ab"])
    with pytest.raises(TypeError):
        GWFrames.Matrix([[[1.0]], [[2.0], [3.0]]])


@pytest.mark.parametrize("data", [[[], []], np.zeros((2, 0), dtype=complex)])
def test_complex_rows_of_length_zero(data):
    # Two modes with no times: the complex caster converts the empty rows on
    # the way in, and the data come back as an array of shape (2, 0).
    W = GWFrames.Waveform([], [[2, 2], [2, -2]], data)
    assert W.NModes() == 2
    assert W.NTimes() == 0
    Data = W.Data()
    assert Data.shape == (2, 0)
    assert Data.dtype == complex
    assert np.array_equal(W.LM(), [[2, 2], [2, -2]])


@pytest.mark.parametrize("lm", [[[]], np.zeros((1, 0), dtype=int)])
def test_int_rows_of_length_zero(lm):
    # The integer caster converts rows of length zero in both directions too.
    W = GWFrames.Waveform([], lm, [])
    LM = W.LM()
    assert LM.shape == (1, 0)
    assert LM.dtype.kind == "i"


def test_nonempty_waveform_data_round_trip():
    t = np.linspace(0.0, 1.0, 5)
    lm = np.array([[2, -2], [2, 2]])
    data = np.array([np.exp(1j * t), np.exp(-1j * t)])
    W = GWFrames.Waveform(t, lm, data)
    assert np.array_equal(W.LM(), lm)
    assert np.array_equal(W.Data(), data)
