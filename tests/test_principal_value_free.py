"""Regression tests for DominantPrincipalValue on more than one matrix.

DominantPrincipalValue used to free its GSL workspace inside the loop over
matrices, so every matrix after the first was diagonalized with freed memory,
and the workspace was freed again on each pass.  The process usually aborted.
"""

import subprocess
import sys

import numpy as np

import GWFrames


def symmetric_matrices(n):
    rng = np.random.default_rng(1234)
    arrays = []
    for _ in range(n):
        a = rng.standard_normal((3, 3))
        arrays.append(a + a.T)
    return arrays


def test_dominant_principal_value_many_matrices_in_subprocess():
    # A double free aborts the interpreter, so run the call in a child process
    # to report it as a failure rather than crashing the test session.
    code = (
        "import GWFrames\n"
        "Ms = [GWFrames.Matrix([[3.0, 1.0, 0.0], [1.0, 2.0, 0.0], [0.0, 0.0, 1.0 + i]])"
        " for i in range(10)]\n"
        "print(list(GWFrames.DominantPrincipalValue(Ms)))\n"
    )
    result = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=60)
    assert result.returncode == 0, result.stderr


def test_principal_values_match_numpy():
    arrays = symmetric_matrices(20)
    expected = [np.linalg.eigvalsh(a) for a in arrays]
    dominant = GWFrames.DominantPrincipalValue([GWFrames.Matrix(a.tolist()) for a in arrays])
    subordinate = GWFrames.SubordinatePrincipalValue([GWFrames.Matrix(a.tolist()) for a in arrays])
    assert len(dominant) == len(arrays)
    assert np.allclose(dominant, [e[-1] for e in expected], rtol=1e-12, atol=1e-12)
    assert np.allclose(subordinate, [e[0] for e in expected], rtol=1e-12, atol=1e-12)
