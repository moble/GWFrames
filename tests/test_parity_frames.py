import numpy as np
import pytest

import GWFrames


@pytest.fixture(scope="module")
def pn():
    """A short precessing PN waveform, in the coorbital frame as constructed"""
    return GWFrames.PNWaveform("TaylorT4", 0.1, [0.3, 0.2, 0.4], [-0.2, 0.1, 0.3], 0.05)


PARTS = [
    "XParitySymmetricPart",
    "XParityAntisymmetricPart",
    "YParitySymmetricPart",
    "YParityAntisymmetricPart",
    "ZParitySymmetricPart",
    "ZParityAntisymmetricPart",
    "ParitySymmetricPart",
    "ParityAntisymmetricPart",
    "ConjugateAntipodalSymmetricPart",
    "ConjugateAntipodalAntisymmetricPart",
]


def frame_array(W):
    return np.array([[q[0], q[1], q[2], q[3]] for q in W.Frame()])


@pytest.mark.parametrize("name", PARTS)
def test_parity_parts_keep_frame(pn, name):
    W = getattr(pn, name)()
    R = frame_array(W)
    assert R.shape == (pn.NTimes(), 4)
    np.testing.assert_array_equal(R, frame_array(pn))
    np.testing.assert_allclose(np.linalg.norm(R, axis=1), 1.0, rtol=1e-12)
    assert W.FrameType() == pn.FrameType()


@pytest.mark.parametrize("axis", ["X", "Y", "Z", ""])
def test_parity_parts_sum_to_waveform(pn, axis):
    S = getattr(pn, axis + "ParitySymmetricPart")()
    A = getattr(pn, axis + "ParityAntisymmetricPart")()
    np.testing.assert_allclose((S + A).Data(), pn.Data(), rtol=0, atol=1e-14 * np.max(np.abs(pn.Data())))
