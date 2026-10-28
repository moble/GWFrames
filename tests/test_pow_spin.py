"""Integer powers of spin-weighted data have spin weight p times the original"""

import numpy as np
import pytest

import GWFrames


def spin_modes(s, ell_max, seed=1234):
    """Random modes of spin weight `s`, with the modes for ell < |s| set to zero"""
    rng = np.random.default_rng(seed)
    data = rng.normal(size=(ell_max + 1) ** 2) + 1j * rng.normal(size=(ell_max + 1) ** 2)
    data[: abs(s) ** 2] = 0.0
    return GWFrames.Modes(s, data)


@pytest.mark.parametrize("p", [0, 1, 2, 3, -1])
def test_datagrid_pow_spin(p):
    M = spin_modes(-2, 4)
    D = GWFrames.DataGrid(M, 9, 9)
    assert D.Spin() == -2
    Dp = D.pow(p)
    assert Dp.Spin() == p * (-2)
    np.testing.assert_allclose(Dp.Data(), np.asarray(D.Data()) ** p)


def test_datagrid_pow_matches_product():
    D = GWFrames.DataGrid(spin_modes(1, 3), 7, 7)
    assert D.pow(2).Spin() == (D * D).Spin()
    assert D.pow(3).Spin() == (D * D * D).Spin()


@pytest.mark.parametrize("s", [1, -2])
def test_modes_pow_spin_and_values(s):
    M = spin_modes(s, 3)
    M2 = M.pow(2)
    assert M2.Spin() == 2 * s
    assert M2.Spin() == (M * M).Spin()
    # The modes must represent the pointwise square of the original function
    n = 2 * M2.EllMax() + 1
    np.testing.assert_allclose(
        GWFrames.DataGrid(M2, n, n).Data(), np.asarray(GWFrames.DataGrid(M, n, n).Data()) ** 2, rtol=1e-9, atol=1e-9
    )
