import numpy as np

import GWFrames


def waveform():
    T = np.linspace(0, 10, 21)
    LM = [[2, m] for m in range(-2, 3)]
    Data = np.exp(1j * np.outer(np.arange(1, 6), T))
    return GWFrames.Waveform(T, LM, Data)


def test_interpolate_in_place_appends_one_line_to_history():
    W = waveform()
    W.SetHistory("# first line\n")
    before = W.HistoryStr()
    W.InterpolateInPlace(np.linspace(1, 9, 17))
    after = W.HistoryStr()
    # The old history appears exactly once, followed by a single new line
    assert after.startswith(before)
    assert after.count(before) == 1
    added = after[len(before):]
    assert added.count("\n") == 1
    assert "InterpolateInPlace" in added


def test_repeated_interpolate_in_place_grows_history_linearly():
    W = waveform()
    W.SetHistory("# first line\n")
    lengths = [len(W.HistoryStr())]
    for i in range(4):
        W.InterpolateInPlace(np.linspace(1 + 0.1 * i, 9 - 0.1 * i, 17))
        lengths.append(len(W.HistoryStr()))
    increments = np.diff(lengths)
    assert np.all(increments == increments[0])
    assert W.HistoryStr().count("# first line") == 1


def test_interpolate_in_place_data_matches_interpolate():
    # The in-place version still interpolates exactly as Interpolate does
    W = waveform()
    NewTime = np.linspace(1, 9, 17)
    W2 = W.Interpolate(NewTime)
    W.InterpolateInPlace(NewTime)
    np.testing.assert_array_equal(W.T(), W2.T())
    np.testing.assert_allclose(W.Data(), W2.Data(), rtol=0, atol=1e-14)
