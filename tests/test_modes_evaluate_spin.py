"""Evaluation of spin-weighted `Modes` at a point on the sphere

`Modes` stores every mode from ell = 0, even those with ell < |s|, which vanish
for a spin-weighted function.  `Modes.EvaluateAtPoint` once tried to evaluate
the harmonics for those modes too, and so raised an error for any nonzero spin
weight.
"""

import math

import numpy as np
import pytest

import GWFrames
from GWFrames import Quaternions

ELL_MAX = 4
SPINS = [-3, -2, -1, 0, 1, 2]
POINTS = [(theta, phi) for theta in np.linspace(0.1, 3.0, 4) for phi in np.linspace(0.2, 6.0, 4)]


def spin_modes(s, seed):
    """Random modes of spin weight `s`, with the modes for ell < |s| set to zero"""
    rng = np.random.default_rng(seed)
    data = rng.normal(size=(ELL_MAX + 1) ** 2) + 1j * rng.normal(size=(ELL_MAX + 1) ** 2)
    data[: s * s] = 0.0
    return data


def goldberg_sYlm(s, ell, m, theta, phi):
    """The spin-weighted spherical harmonic, from the explicit sum of Goldberg et al.

    This gives, for example, the familiar
    -2Y22 = sqrt(5 / (64 pi)) (1 + cos(theta))^2 exp(2 i phi).
    """
    prefactor = (-1) ** m * math.sqrt(
        math.factorial(ell + m) * math.factorial(ell - m) * (2 * ell + 1)
        / (4 * math.pi * math.factorial(ell + s) * math.factorial(ell - s))
    )
    total = 0.0
    for r in range(ell - s + 1):
        k = r + s - m
        if 0 <= k <= ell + s:
            total += (math.comb(ell - s, r) * math.comb(ell + s, k) * (-1) ** (ell - r - s)
                      / math.tan(theta / 2) ** (2 * r + s - m))
    return prefactor * math.sin(theta / 2) ** (2 * ell) * total * np.exp(1j * m * phi)


@pytest.mark.parametrize("s", SPINS)
def test_single_modes_match_explicit_harmonics(s):
    theta, phi = 0.7, 1.9
    for ell in range(abs(s), ELL_MAX + 1):
        for m in range(-ell, ell + 1):
            data = np.zeros((ELL_MAX + 1) ** 2, dtype=complex)
            data[ell * ell + ell + m] = 1.0
            M = GWFrames.Modes(s, data)
            expected = goldberg_sYlm(s, ell, m, theta, phi)
            assert M.EvaluateAtPoint(theta, phi) == pytest.approx(expected, abs=1e-13)


def test_minus_two_two_two():
    data = np.zeros((ELL_MAX + 1) ** 2, dtype=complex)
    data[2 * 2 + 2 + 2] = 1.0
    M = GWFrames.Modes(-2, data)
    for theta, phi in POINTS:
        expected = math.sqrt(5 / (64 * math.pi)) * (1 + math.cos(theta)) ** 2 * np.exp(2j * phi)
        assert M.EvaluateAtPoint(theta, phi) == pytest.approx(expected, abs=1e-14)


@pytest.mark.parametrize("s", SPINS)
def test_matches_datagrid_at_grid_points(s):
    # The spinsfast grid has theta = pi * i / (n_theta - 1) and
    # phi = 2 pi * j / n_phi, with theta varying slowest.
    M = GWFrames.Modes(s, spin_modes(s, seed=40 + s))
    n_theta, n_phi = 2 * ELL_MAX + 3, 2 * ELL_MAX + 2
    grid = np.array(GWFrames.DataGrid(M, n_theta, n_phi).Data()).reshape(n_theta, n_phi)
    values = np.array([[M.EvaluateAtPoint(np.pi * i / (n_theta - 1), 2 * np.pi * j / n_phi)
                        for j in range(n_phi)] for i in range(n_theta)])
    np.testing.assert_allclose(values, grid, rtol=0, atol=1e-12)


@pytest.mark.parametrize("s", SPINS)
def test_matches_waveform_evaluation(s):
    data = spin_modes(s, seed=50 + s)
    M = GWFrames.Modes(s, data)
    LM = [[ell, m] for ell in range(abs(s), ELL_MAX + 1) for m in range(-ell, ell + 1)]
    W = GWFrames.Waveform([0.0, 1.0], LM, np.array([[data[ell * ell + ell + m]] * 2 for ell, m in LM]))
    W.SetSpinWeight(s)
    W.SetFrameType(GWFrames.Inertial)
    for theta, phi in POINTS:
        assert M.EvaluateAtPoint(theta, phi) == pytest.approx(W.EvaluateAtPoint(theta, phi)[0], abs=1e-13)


@pytest.mark.parametrize("s", SPINS)
def test_rotor_overload_matches_angles(s):
    M = GWFrames.Modes(s, spin_modes(s, seed=60 + s))
    for theta, phi in POINTS:
        R = Quaternions.Quaternion(theta, phi)
        assert M.EvaluateAtPoint(R) == pytest.approx(M.EvaluateAtPoint(theta, phi), abs=1e-14)


@pytest.mark.parametrize("s", [-2, 1])
def test_modes_below_spin_are_skipped(s):
    data = spin_modes(s, seed=70 + s)
    junk = data.copy()
    junk[: s * s] = 1.0 + 2.0j
    M, J = GWFrames.Modes(s, data), GWFrames.Modes(s, junk)
    for theta, phi in POINTS:
        assert J.EvaluateAtPoint(theta, phi) == M.EvaluateAtPoint(theta, phi)
        R = Quaternions.Quaternion(theta, phi)
        assert J.EvaluateAtPoint(R) == M.EvaluateAtPoint(R)


@pytest.mark.parametrize("s", [-2, -1, 1])
def test_unboosted_datagrid_from_evaluation(s):
    # The boosted DataGrid constructor evaluates M with the rotor overload; with
    # zero velocity it must reproduce the ordinary transform.
    M = GWFrames.Modes(s, spin_modes(s, seed=80 + s))
    n = 2 * ELL_MAX + 1
    boosted = np.array(GWFrames.DataGrid(M, [0.0, 0.0, 0.0], n, n).Data())
    np.testing.assert_allclose(boosted, GWFrames.DataGrid(M, n, n).Data(), rtol=0, atol=1e-12)
