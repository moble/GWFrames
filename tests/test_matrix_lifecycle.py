"""Regression tests for the lifecycle of an empty or cleared GWFrames.Matrix

Before the fix, `Matrix.clear` freed its storage without forgetting it, so the
destructor freed it a second time, and every query of an empty Matrix
dereferenced a null pointer.  Both crash the interpreter rather than raise, so
each check runs in a separate Python process, and a crash shows up as a failed
test instead of ending the whole test session.
"""

import os
import subprocess
import sys
import textwrap

import pytest


def run_isolated(code):
    """Run `code` in a fresh interpreter, and return its standard output"""
    preamble = "import numpy as np\nimport GWFrames\n"
    result = subprocess.run(
        [sys.executable, "-c", preamble + textwrap.dedent(code)],
        capture_output=True,
        text=True,
        timeout=60,
        env=os.environ.copy(),
    )
    assert result.returncode == 0, (
        f"subprocess exited with code {result.returncode}\n"
        f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
    )
    return result.stdout


def test_clear_twice_then_destroy():
    out = run_isolated(
        """
        M = GWFrames.Matrix(3, 3, 1.5)
        M.clear()
        M.clear()
        print(M.nrows(), M.ncols())
        del M
        import gc; gc.collect()
        print("done")
        """
    )
    assert out.split() == ["0", "0", "done"]


@pytest.mark.parametrize(
    "construct",
    [
        "GWFrames.Matrix()",
        "GWFrames.Matrix(0, 3)",
        "GWFrames.Matrix(3, 0, 2.0)",
        "GWFrames.Matrix([])",
    ],
)
def test_empty_matrix(construct):
    out = run_isolated(
        f"""
        M = {construct}
        assert (M.nrows(), M.ncols()) == (0, 0)
        C = GWFrames.Matrix(M)
        assert (C.nrows(), C.ncols()) == (0, 0)
        a = M.numpy()
        assert a.shape == (0, 0), a.shape
        assert repr(M) == "Matrix([])", repr(M)
        D = M - C
        assert (D.nrows(), D.ncols()) == (0, 0)
        try:
            M(0, 0)
        except IndexError:
            pass
        else:
            raise AssertionError("M(0, 0) did not raise")
        M.clear()
        M.resize(2, 3, 4.0)
        assert M.numpy().tolist() == [[4.0] * 3] * 2
        M.resize(0, 0)
        assert (M.nrows(), M.ncols()) == (0, 0)
        del M, C, D
        print("done")
        """
    )
    assert out.split() == ["done"]


def test_cleared_matrix_is_usable():
    out = run_isolated(
        """
        M = GWFrames.Matrix([[1.0, 2.0], [3.0, 4.0]])
        M.clear()
        assert (M.nrows(), M.ncols()) == (0, 0)
        assert M.numpy().shape == (0, 0)
        assert repr(M) == "Matrix([])"
        C = GWFrames.Matrix(M)
        assert (C.nrows(), C.ncols()) == (0, 0)
        M.resize(2, 2)
        M.set(1, 0, 7.0)
        assert M.numpy().tolist() == [[0.0, 0.0], [7.0, 0.0]]
        print("done")
        """
    )
    assert out.split() == ["done"]


def test_nonempty_matrix_unchanged():
    out = run_isolated(
        """
        M = GWFrames.Matrix([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
        C = GWFrames.Matrix(M)
        M.set(0, 0, 9.0)
        assert C.numpy().tolist() == [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        assert M.numpy().tolist() == [[9.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
        assert (M - C).numpy().tolist() == [[8.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        assert list(M * [1.0, 1.0, 1.0]) == [14.0, 15.0]
        print("done")
        """
    )
    assert out.split() == ["done"]
