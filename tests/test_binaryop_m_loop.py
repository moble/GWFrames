"""Pointwise products and quotients of Waveform objects

These operators transform both waveforms to values on a grid on the sphere,
combine them pointwise, and transform back.  They once loaded only the
m = -ell coefficient of each ell into the transform, and they gave a quotient
the sum of the spin weights instead of their difference.  The modes of a
product must be exact, even when the product has higher modes than are
returned.
"""

import numpy as np
import pytest

import GWFrames

T = np.linspace(0.0, 1.0, 3)


def waveform(rng, s, ellMax, ellMaxNonzero=None):
    """A spin-weight-s waveform with every mode up to ellMax

    Modes up to ellMaxNonzero (by default ellMax) hold random data; the rest
    are zero.
    """
    if ellMaxNonzero is None:
        ellMaxNonzero = ellMax
    LM = [[ell, m] for ell in range(abs(s), ellMax + 1) for m in range(-ell, ell + 1)]
    Data = np.zeros((len(LM), len(T)), dtype=complex)
    for i, (ell, m) in enumerate(LM):
        if ell <= ellMaxNonzero:
            Data[i] = rng.normal(size=len(T)) + 1j * rng.normal(size=len(T))
    W = GWFrames.Waveform(T, LM, Data)
    W.SetSpinWeight(s)
    W.SetFrameType(GWFrames.Inertial)
    return W


def constant(c, ellMax):
    """A spin-0 waveform representing the constant function c on the sphere"""
    LM = [[ell, m] for ell in range(ellMax + 1) for m in range(-ell, ell + 1)]
    Data = np.zeros((len(LM), len(T)), dtype=complex)
    Data[0] = c * np.sqrt(4 * np.pi)
    W = GWFrames.Waveform(T, LM, Data)
    W.SetSpinWeight(0)
    W.SetFrameType(GWFrames.Inertial)
    return W


def values(W):
    """Values of W at a set of points that are not on the transform grid"""
    return np.array([W.EvaluateAtPoint(theta, phi)
                     for theta in np.linspace(0.1, 3.0, 6)
                     for phi in np.linspace(0.2, 6.0, 6)])


def expected_LM(s, ellMax):
    return [[ell, m] for ell in range(abs(s), ellMax + 1) for m in range(-ell, ell + 1)]


def test_multiply_and_divide_by_constant():
    rng = np.random.default_rng(1234)
    A = waveform(rng, -2, 4)
    for C, expected in [(A * constant(1.0, 4), A.Data()),
                        (A / constant(1.0, 4), A.Data()),
                        (A * constant(3.0, 4), 3.0 * A.Data()),
                        (A / constant(2.0, 4), A.Data() / 2.0)]:
        assert C.SpinWeight() == -2
        np.testing.assert_array_equal(C.LM(), A.LM())
        np.testing.assert_allclose(C.Data(), expected, rtol=0, atol=1e-13)


@pytest.mark.parametrize("sA, sB", [(0, 0), (-2, 0), (0, -2), (-1, 1), (1, 1), (-2, 2), (2, -1)])
def test_product_matches_pointwise_values(sA, sB):
    # Each factor has modes only up to ell = 3, so the product has modes only
    # up to ell = 6, which the result holds exactly.
    rng = np.random.default_rng(100 * (sA + 3) + sB + 3)
    A = waveform(rng, sA, 6, 3)
    B = waveform(rng, sB, 6, 3)
    C = A * B
    assert C.SpinWeight() == sA + sB
    np.testing.assert_array_equal(C.LM(), expected_LM(sA + sB, 6))
    np.testing.assert_allclose(values(C), values(A) * values(B), rtol=0, atol=1e-12)


def test_product_is_truncated_to_the_smaller_ell_max():
    rng = np.random.default_rng(5678)
    A = waveform(rng, -2, 5, 2)
    B = waveform(rng, 1, 3, 1)
    C = A * B
    assert C.SpinWeight() == -1
    np.testing.assert_array_equal(C.LM(), expected_LM(-1, 3))
    np.testing.assert_allclose(values(C), values(A) * values(B), rtol=0, atol=1e-12)


def test_truncated_product_is_exact():
    # Every mode of each factor holds data, so the product has modes up to
    # ell = 8, beyond the ell = 3 returned.  The returned modes must equal
    # those of the whole product, which is found by padding the factors with
    # zero modes up to ell = 8.
    rng = np.random.default_rng(161718)
    A = waveform(rng, -2, 5)
    B = waveform(rng, 1, 3)
    C = A * B
    np.testing.assert_array_equal(C.LM(), expected_LM(-1, 3))
    APadded = GWFrames.Waveform(T, expected_LM(-2, 8),
                                np.vstack([A.Data(), np.zeros((len(expected_LM(-2, 8)) - A.NModes(), len(T)))]))
    APadded.SetSpinWeight(-2)
    APadded.SetFrameType(GWFrames.Inertial)
    BPadded = GWFrames.Waveform(T, expected_LM(1, 8),
                                np.vstack([B.Data(), np.zeros((len(expected_LM(1, 8)) - B.NModes(), len(T)))]))
    BPadded.SetSpinWeight(1)
    BPadded.SetFrameType(GWFrames.Inertial)
    Full = APadded * BPadded
    np.testing.assert_allclose(values(Full), values(A) * values(B), rtol=0, atol=1e-11)
    np.testing.assert_allclose(C.Data(), Full.Data()[:C.NModes()], rtol=0, atol=1e-11)


def test_quotient_of_waveform_by_itself():
    rng = np.random.default_rng(91011)
    W = waveform(rng, -2, 4)
    Q = W / W
    assert Q.SpinWeight() == 0
    np.testing.assert_array_equal(Q.LM(), expected_LM(0, 4))
    np.testing.assert_allclose(Q.Data(), constant(1.0, 4).Data(), rtol=0, atol=1e-10)


def test_quotient_undoes_product():
    rng = np.random.default_rng(121314)
    A = waveform(rng, -2, 4, 2)
    B = waveform(rng, 1, 4, 1)
    Q = (A * B) / B
    assert Q.SpinWeight() == -2
    np.testing.assert_array_equal(Q.LM(), A.LM())
    np.testing.assert_allclose(Q.Data(), A.Data(), rtol=0, atol=1e-9)
