"""Regression tests for spin-0 data with only the ell=0 mode.

Data with only ell=0 used to be put on a 1x1 grid, which crashed spinsfast or
gave garbage, so the crashing cases run in a subprocess.
"""

import subprocess
import sys
import textwrap

import pytest


def run(code):
    """Run code in a fresh interpreter, and return its stdout; fail if it crashes"""
    result = subprocess.run([sys.executable, "-c", textwrap.dedent(code)], capture_output=True, text=True)
    assert result.returncode == 0, f"exit code {result.returncode}\n{result.stderr[-2000:]}"
    return result.stdout


def test_product_of_ell_0_waveforms():
    # The constant function 1 has the single mode sqrt(4 pi) at (0, 0), and so
    # does its square
    out = run(
        """
        import numpy as np, GWFrames
        c = np.sqrt(4 * np.pi)
        W = GWFrames.Waveform([0.0, 1.0, 2.0], [[0, 0]], [[c, c, c]])
        W.SetSpinWeight(0)  # The default is -2, for strain
        P = W * W
        print(P.LM().tolist(), np.abs(P.Data() - c).max())
        """
    )
    lm, error = out.split("]] ")
    assert lm + "]]" == "[[0, 0]]"
    assert float(error) < 1e-12


@pytest.mark.parametrize("expression", ["M.pow(2)", "M * M"])
def test_square_of_ell_0_modes(expression):
    # The function 2 has the single mode 2 sqrt(4 pi), and its square 4 has the
    # mode 4 sqrt(4 pi)
    out = run(
        f"""
        import numpy as np, GWFrames
        c = np.sqrt(4 * np.pi)
        M = GWFrames.Modes(0, [2 * c])
        print(abs(({expression}).Data()[0] - 4 * c))
        """
    )
    assert float(out) < 1e-12

