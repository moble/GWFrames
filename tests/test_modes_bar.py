"""Check that `Modes.bar` returns the modes of the complex-conjugate function."""

import numpy as np
import pytest

import GWFrames

ELL_MAX = 4


def random_modes(s, seed):
    rng = np.random.default_rng(seed)
    data = rng.normal(size=(ELL_MAX + 1) ** 2) + 1j * rng.normal(size=(ELL_MAX + 1) ** 2)
    data[: s * s] = 0.0  # Modes with ell < |s| vanish
    return GWFrames.Modes(s, data)


@pytest.mark.parametrize("s", [-2, -1, 0, 1, 2])
def test_bar_matches_conjugated_grid(s):
    M = random_modes(s, seed=10 + s)
    n = 2 * ELL_MAX + 1
    grid = np.array(GWFrames.DataGrid(M, n, n).Data())
    expected = GWFrames.Modes(GWFrames.DataGrid(-s, n, n, np.conj(grid)), ELL_MAX)
    Mbar = M.bar()
    assert Mbar.Spin() == -s
    np.testing.assert_allclose(Mbar.Data(), expected.Data(), rtol=0, atol=1e-12)


@pytest.mark.parametrize("s", [-2, 0, 1])
def test_bar_is_an_involution(s):
    M = random_modes(s, seed=20 + s)
    np.testing.assert_allclose(M.bar().bar().Data(), M.Data(), rtol=0, atol=1e-14)
    assert M.bar().bar().Spin() == s


def test_bar_of_real_function_is_itself():
    # A real spin-0 function has f_{l,-m} = (-1)^m conj(f_{l,m})
    M = random_modes(0, seed=30)
    data = np.array(M.Data())
    for ell in range(ELL_MAX + 1):
        for m in range(-ell, ell + 1):
            if m < 0:
                data[ell * ell + ell + m] = (-1) ** m * np.conj(data[ell * ell + ell - m])
            elif m == 0:
                data[ell * ell + ell] = data[ell * ell + ell].real
    M = GWFrames.Modes(0, data)
    np.testing.assert_allclose(M.bar().Data(), data, rtol=0, atol=1e-14)
