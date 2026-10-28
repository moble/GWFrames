"""Arithmetic between numbers and spin-weighted grids keeps track of the spin weight"""

import numpy as np
import pytest

import GWFrames


def spin_grid(s, ell_max=3, seed=4321):
    """The grid of random modes of spin weight `s`, with the modes for ell < |s| set to zero"""
    rng = np.random.default_rng(seed)
    data = rng.normal(size=(ell_max + 1) ** 2) + 1j * rng.normal(size=(ell_max + 1) ** 2)
    data[: abs(s) ** 2] = 0.0
    n = 2 * ell_max + 1
    return GWFrames.DataGrid(GWFrames.Modes(s, data), n, n)


@pytest.mark.parametrize("s", [-2, -1, 1, 2])
def test_scalar_divided_by_grid_has_negated_spin(s):
    D = spin_grid(s)
    Q = 2.5 / D
    assert Q.Spin() == -s
    np.testing.assert_allclose(Q.Data(), 2.5 / np.asarray(D.Data()))


@pytest.mark.parametrize("s", [-2, 1])
def test_scalar_division_matches_pow_and_grid_division(s):
    D = spin_grid(s)
    one = GWFrames.DataGrid(0, D.N_theta(), D.N_phi(), np.ones(len(D), dtype=complex))
    assert (1.0 / D).Spin() == D.pow(-1).Spin()
    assert (1.0 / D).Spin() == (one / D).Spin()
    assert (1.0 / D * D).Spin() == 0


@pytest.mark.parametrize("s", [-2, 0, 1])
def test_other_scalar_operators_keep_spin(s):
    D = spin_grid(s)
    assert (3.0 * D).Spin() == s
    assert (3.0 + D).Spin() == s
    assert (3.0 - D).Spin() == s
    np.testing.assert_allclose((3.0 - D).Data(), 3.0 - np.asarray(D.Data()))
