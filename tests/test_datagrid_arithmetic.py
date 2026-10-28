"""Regression tests for the spin weights in DataGrid arithmetic, and for scalars on the right."""

import numpy as np
import pytest

import GWFrames


def grid(spin, values):
    return GWFrames.DataGrid(spin, 3, 3, list(values))


@pytest.fixture
def values():
    rng = np.random.default_rng(1234)
    return rng.normal(size=9) + 1j * rng.normal(size=9) + 3.0


@pytest.mark.parametrize("op", ["+", "-"])
def test_sum_of_different_spins_raises(values, op):
    A, B = grid(-2, values), grid(0, values)
    with pytest.raises(ValueError):
        eval(f"A {op} B")


@pytest.mark.parametrize("op", ["+", "-"])
def test_sum_of_equal_spins(values, op):
    A, B = grid(-2, values), grid(-2, 2 * values)
    C = eval(f"A {op} B")
    assert C.Spin() == -2
    assert np.allclose(C.Data(), eval(f"values {op} 2 * values"))


@pytest.mark.parametrize("op", ["*", "/", "+", "-"])
def test_scalar_on_the_right(values, op):
    D = grid(-2, values)
    C = eval(f"D {op} 2.5")
    assert C.Spin() == -2
    assert np.allclose(C.Data(), eval(f"values {op} 2.5"))


def test_scalar_on_the_left_still_works(values):
    D = grid(-2, values)
    assert np.allclose((2.5 - D).Data(), 2.5 - values)
    assert (2.5 / D).Spin() == 2


def test_unfinished_moreschi_iteration_is_gone():
    assert not hasattr(GWFrames.SliceModes, "MoreschiIteration")
    assert hasattr(GWFrames.SuperMomenta, "MoreschiIteration")
