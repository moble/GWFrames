"""Pointwise products of Waveform objects without aliasing

The product of two waveforms is computed by transforming both to values on a
grid on the sphere, multiplying pointwise, and transforming back.  The grid
once resolved only the modes up to the smaller of the largest ell values of
the inputs, and the inputs were truncated to that ell, so the higher modes of
the product were aliased into the result.  The result must instead be the
exact product, truncated to the returned modes.
"""

import numpy as np
import pytest

import GWFrames

T = np.linspace(0.0, 1.0, 3)


def waveform(Data, s, ellMax, ellMaxPadded=None):
    """A spin-weight-s waveform with the mode data Data up to ellMax

    If ellMaxPadded is given, modes from ellMax + 1 up to ellMaxPadded are
    included with zero data, which describes the same function on the sphere.
    """
    if ellMaxPadded is None:
        ellMaxPadded = ellMax
    LM = [[ell, m] for ell in range(abs(s), ellMaxPadded + 1) for m in range(-ell, ell + 1)]
    D = np.zeros((len(LM), len(T)), dtype=complex)
    D[:len(Data)] = Data
    W = GWFrames.Waveform(T, LM, D)
    W.SetSpinWeight(s)
    W.SetFrameType(GWFrames.Inertial)
    return W


def random_data(rng, s, ellMax):
    NModes = (ellMax + 1)**2 - abs(s)**2
    return rng.normal(size=(NModes, len(T))) + 1j * rng.normal(size=(NModes, len(T)))


def values(W):
    """Values of W at a set of points that are not on the transform grid"""
    return np.array([W.EvaluateAtPoint(theta, phi)
                     for theta in np.linspace(0.1, 3.0, 7)
                     for phi in np.linspace(0.2, 6.0, 7)])


@pytest.mark.parametrize("sA, sB, ellMaxA, ellMaxB",
                         [(0, 0, 4, 4), (-2, 0, 4, 4), (-2, 2, 4, 3), (1, -1, 2, 5), (-2, 1, 6, 2)])
def test_product_is_exact(sA, sB, ellMaxA, ellMaxB):
    rng = np.random.default_rng(1000 * (sA + 3) + 100 * (sB + 3) + 10 * ellMaxA + ellMaxB)
    DA = random_data(rng, sA, ellMaxA)
    DB = random_data(rng, sB, ellMaxB)
    ellMaxC = min(ellMaxA, ellMaxB)
    C = waveform(DA, sA, ellMaxA) * waveform(DB, sB, ellMaxB)
    assert C.SpinWeight() == sA + sB
    np.testing.assert_array_equal(
        C.LM(), [[ell, m] for ell in range(abs(sA + sB), ellMaxC + 1) for m in range(-ell, ell + 1)])

    # When both factors are padded with zero modes up to the band limit of
    # the product, the returned modes hold the whole product, as its values
    # confirm.  The unpadded product must be the truncation of that.
    ellMaxFull = ellMaxA + ellMaxB
    A = waveform(DA, sA, ellMaxA, ellMaxFull)
    B = waveform(DB, sB, ellMaxB, ellMaxFull)
    Full = A * B
    np.testing.assert_allclose(values(Full), values(A) * values(B), rtol=0, atol=1e-11)
    np.testing.assert_allclose(C.Data(), Full.Data()[:C.NModes()], rtol=0, atol=1e-11)


def test_quotient_uses_every_mode_of_the_numerator():
    # P has modes up to ell = 8, but B has modes only up to ell = 4, so the
    # quotient P/B returns modes only up to ell = 4.  It must still use the
    # higher modes of P.
    rng = np.random.default_rng(31415)
    DA = random_data(rng, -2, 4)
    DB = 0.1 * random_data(rng, 0, 4)
    DB[0] += 3.0 * np.sqrt(4 * np.pi)  # Keep B far from zero
    B = waveform(DB, 0, 4)
    P = waveform(DA, -2, 4, 8) * waveform(DB, 0, 4, 8)
    assert P.EllMax() == 8
    Q = P / B
    assert Q.SpinWeight() == -2
    np.testing.assert_array_equal(Q.LM(), waveform(DA, -2, 4).LM())
    np.testing.assert_allclose(Q.Data(), DA, rtol=0, atol=1e-10)
