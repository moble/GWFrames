"""Regression test for Matrix construction from rows of different lengths."""

import numpy as np
import pytest

import GWFrames


def test_ragged_matrix_rows_raise():
    with pytest.raises(ValueError):
        GWFrames.Matrix([[1.0, 2.0], [3.0]])
    M = GWFrames.Matrix([[1.0, 2.0], [3.0, 4.0]])
    assert np.array_equal(M.numpy(), [[1.0, 2.0], [3.0, 4.0]])
