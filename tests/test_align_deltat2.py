import numpy as np
import pytest

import GWFrames


@pytest.fixture(scope="module")
def short_inertial():
    """A short PN waveform in the inertial frame, resampled coarsely for speed"""
    pn = GWFrames.PNWaveform("TaylorT4", 0.1, [0.3, 0.2, 0.4], [-0.2, 0.1, 0.3], 0.05)
    W = GWFrames.Waveform(pn).TransformToInertialFrame()
    return W.Interpolate(np.arange(W.T()[0], W.T()[-1], 0.5))


@pytest.mark.parametrize("T0", [0.0, 1.0e4])
def test_align_waveforms_recovers_time_shift(short_inertial, T0):
    # `B` is a copy of `A` shifted later in time by 3.  The upper bound of the
    # coarse search over time shifts used to be computed from the number of
    # time steps of `B` instead of its final time, so when the times were
    # large compared to the number of steps (T0=1e4 here), the search window
    # was empty and the alignment failed.
    shift = 3.0
    A = GWFrames.Waveform(short_inertial)
    A.SetTime(A.T() + T0)
    B = GWFrames.Waveform(short_inertial)
    B.SetTime(B.T() + T0 + shift)
    GWFrames.AlignWaveforms(A, B, T0 + 40.0, T0 + 100.0)
    np.testing.assert_allclose(B.T(), A.T(), rtol=0, atol=1e-8)
