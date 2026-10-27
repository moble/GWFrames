import copy
import pickle

import numpy as np
import pytest

from GWFrames import Quaternions
from GWFrames.Quaternions import Quaternion


def test_construction_and_access():
    q = Quaternion(1, 2, 3, 4)
    assert list(q) == [1, 2, 3, 4]
    assert len(q) == 4
    assert q[-1] == 4
    with pytest.raises(IndexError):
        q[4]
    q[0] = 5
    assert q[0] == 5
    assert Quaternion([1, 2, 3, 4]) == Quaternion(1, 2, 3, 4)
    assert Quaternion([2, 3, 4]) == Quaternion(0, 2, 3, 4)
    assert repr(Quaternion(1, 2, 3, 4)) == "Quaternion(1, 2, 3, 4)"
    assert str(Quaternion(1, 2, 3, 4)) == "[1, 2, 3, 4]"


def test_equality():
    # These were inverted in the C++ code before GWFrames 2.0
    assert Quaternion(1, 2, 3, 4) == Quaternion(1, 2, 3, 4)
    assert not (Quaternion(1, 2, 3, 4) != Quaternion(1, 2, 3, 4))
    assert Quaternion(1, 2, 3, 4) != Quaternion(1, 2, 3, 5)
    assert not (Quaternion(1, 2, 3, 4) == Quaternion(1, 2, 3, 5))


def test_algebra():
    q = Quaternion(1, 2, 3, 4)
    p = Quaternion(-0.5, 0.25, 1.5, -2.0)
    assert Quaternions.xHat * Quaternions.yHat == Quaternions.zHat
    assert 2 * q == q * 2 == Quaternion(2, 4, 6, 8)
    assert q + 1 == 1 + q == Quaternion(2, 2, 3, 4)
    assert -q == Quaternion(-1, -2, -3, -4)
    np.testing.assert_allclose(list(q * q.inverse()), [1, 0, 0, 0], atol=1e-15)
    np.testing.assert_allclose(list((q * p) / p), list(q), atol=1e-14)
    np.testing.assert_allclose(list(q.log().exp()), list(q), atol=1e-14)
    np.testing.assert_allclose(list(q**2), list(q * q), atol=1e-13)
    assert q.abs() == pytest.approx(np.sqrt(30))
    R = q.normalized()
    np.testing.assert_allclose(list(R.sqrtOfRotor() * R.sqrtOfRotor()), list(R), atol=1e-14)


def test_sequences_convert_to_quaternions():
    assert Quaternions.normalized([0, 0, 0, 2]) == Quaternions.zHat
    np.testing.assert_allclose(list(Quaternions.log([1.0, 0, 0, 0])), [0, 0, 0, 0])


def test_arrays():
    components = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0]], dtype=float)
    for Q in (components, [Quaternion(*c) for c in components], components.tolist()):
        logs = Quaternions.log(Q)
        assert isinstance(logs, np.ndarray)
        assert logs.dtype == object
        assert all(isinstance(q, Quaternion) for q in logs)
    np.testing.assert_allclose(Quaternions.abs(2 * components), [2, 2, 2])
    np.testing.assert_allclose(Quaternions.vec(components), components[:, 1:])
    assert Quaternions.log(np.empty((0, 4))).shape == (0,)


def test_pickle_and_copy():
    q = Quaternion(1, 2, 3, 4)
    assert pickle.loads(pickle.dumps(q)) == q
    assert copy.copy(q) == q
    assert copy.deepcopy(q) == q
    assert copy.copy(q) is not q


def test_errors():
    with pytest.raises(TypeError):
        Quaternion("abcd")
    with pytest.raises(TypeError):
        Quaternions.log("abcd")
