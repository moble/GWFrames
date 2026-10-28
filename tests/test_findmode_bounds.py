"""Regression tests for missing modes in FindModeIndexWithoutError and NormalizedAntisymmetry"""

import numpy as np

import GWFrames


def make_waveform(LM, spin_weight=-2):
    T = np.linspace(0, 10, 11)
    Data = np.exp(1j * np.outer(np.arange(1, len(LM) + 1), T)) * np.arange(1, len(LM) + 1)[:, np.newaxis]
    W = GWFrames.Waveform(T, LM, Data)
    W.SetSpinWeight(spin_weight)
    return W


def expected_antisymmetry(W):
    """Evaluate the normalized antisymmetry directly, treating missing modes as zero"""
    modes = {tuple(lm): W.Data(i) for i, lm in enumerate(W.LM())}
    zero = np.zeros(W.NTimes(), dtype=complex)
    diff = np.zeros(W.NTimes())
    norm = np.zeros(W.NTimes())
    for ell in range(abs(W.SpinWeight()), W.EllMax() + 1):
        for m in range(-ell, ell + 1):
            h = modes.get((ell, m), zero)
            hbar = np.conj(modes.get((ell, -m), zero))
            diff += np.abs(h - (-1) ** (ell + m) * hbar) ** 2
            norm += np.abs(h) ** 2
    return np.sqrt(diff / (4 * norm))


def test_find_mode_index_without_error_missing_mode():
    W = make_waveform([[2, m] for m in range(-2, 3)])
    assert W.FindModeIndexWithoutError(2, 1) == 3
    # (3, 0) would be expected at index 8, beyond the end of the five modes
    assert W.FindModeIndexWithoutError(3, 0) == W.NModes() + 1
    assert W.FindModeIndexWithoutError(8, 8) == W.NModes() + 1


def test_find_mode_index_without_error_out_of_order():
    LM = [[2, 2], [2, 0], [3, 1]]
    W = make_waveform(LM)
    for i, (ell, m) in enumerate(LM):
        assert W.FindModeIndexWithoutError(ell, m) == i
    assert W.FindModeIndexWithoutError(2, -2) == W.NModes() + 1


def test_normalized_antisymmetry_single_mode():
    # With only (2, 2), the missing (2, -2) mode is zero, so a = 1/sqrt(2)
    W = make_waveform([[2, 2]])
    np.testing.assert_allclose(W.NormalizedAntisymmetry(), np.full(W.NTimes(), 1 / np.sqrt(2)), rtol=1e-14)


def test_normalized_antisymmetry_missing_modes():
    LM = [[2, 2], [2, -1], [2, 0], [3, 3], [3, -2], [4, 1], [4, -1]]
    W = make_waveform(LM)
    np.testing.assert_allclose(W.NormalizedAntisymmetry(), expected_antisymmetry(W), rtol=1e-13)


def test_normalized_antisymmetry_complete_modes_unchanged():
    LM = [[ell, m] for ell in range(2, 5) for m in range(-ell, ell + 1)]
    W = make_waveform(LM)
    np.testing.assert_allclose(W.NormalizedAntisymmetry(), expected_antisymmetry(W), rtol=1e-13)
