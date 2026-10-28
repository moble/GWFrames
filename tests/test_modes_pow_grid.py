"""Modes.pow evaluates negative powers on a grid as fine as for abs(p)"""

import numpy as np
import pytest

import GWFrames


def nonvanishing_modes(ell_max, seed=1234):
    """Random spin-0 modes of a function that stays well away from zero"""
    rng = np.random.default_rng(seed)
    data = 0.15 * (rng.normal(size=(ell_max + 1) ** 2) + 1j * rng.normal(size=(ell_max + 1) ** 2))
    data[0] = 2 * np.sqrt(4 * np.pi)
    return GWFrames.Modes(0, data)


def power_on_grid(M, p, n, ell_max):
    """Modes, truncated at `ell_max`, of the p-th power taken on an n-by-n grid"""
    A = GWFrames.Modes(GWFrames.DataGrid(M, n, n).pow(p))
    return np.asarray(A.Data())[: (ell_max + 1) ** 2]


@pytest.mark.parametrize("p", [-1, -2, -3])
def test_negative_power_resolution(p):
    L = 3
    M = nonvanishing_modes(L)
    Mp = M.pow(p)
    # The number of returned modes and the spin weight are as before
    assert Mp.EllMax() == L
    assert len(Mp) == (L + 1) ** 2
    assert Mp.Spin() == 0
    # A negative power is not band-limited, so compare with a very fine grid.  On the
    # old grid of 2*L+1 points the errors were about 1e-2 for each of these powers.
    reference = power_on_grid(M, p, 201, L)
    np.testing.assert_allclose(Mp.Data(), power_on_grid(M, p, 2 * L * abs(p) + 1, L), rtol=0, atol=1e-13)
    tolerance = {-1: 2e-2, -2: 1e-3, -3: 3e-5}[p]
    assert np.max(np.abs(np.asarray(Mp.Data()) - reference)) < tolerance


@pytest.mark.parametrize("p", [0, 1, 2, 3])
def test_nonnegative_power_exact(p):
    L = 3
    M = nonvanishing_modes(L)
    Mp = M.pow(p)
    assert Mp.EllMax() == max(L * p, L)
    assert Mp.Spin() == 0
    # The result is band-limited, so a finer grid changes nothing beyond roundoff
    np.testing.assert_allclose(Mp.Data(), power_on_grid(M, p, 61, Mp.EllMax()), rtol=0, atol=1e-11)
    # And the returned modes reproduce the pointwise power of the original function
    n = 2 * Mp.EllMax() + 3
    np.testing.assert_allclose(
        GWFrames.DataGrid(Mp, n, n).Data(),
        np.asarray(GWFrames.DataGrid(M, n, n).Data()) ** p,
        rtol=1e-11,
        atol=1e-11,
    )
