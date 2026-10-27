"""Manipulate gravitational waveforms — changing frames, and so on.

This package implements the methods described in the paper "Angular
velocity of gravitational radiation from precessing binaries and the
corotating frame" <https://arxiv.org/abs/1302.2919>, along with many
others.  The main objects are `Waveform` and `PNWaveform`; quaternions
and functions of arrays of quaternions are in `GWFrames.Quaternions`.

For the methods available on `Waveform` objects, run `help(GWFrames.Waveform)`.
To construct a Waveform from a data file, use `GWFrames.ReadFromH5` or
`GWFrames.ReadFromNRAR`.

The Waveform's `frame` data records the rotors needed to rotate the
standard (x,y,z) basis into the (X,Y,Z) basis with respect to which the
Waveform's mode `data` are decomposed.

Note on Waveform types: in any system, h — being strain — should be
dimensionless.  When G=c=1, the dimensionless quantities are rMPsi4,
rhdot, and rhOverM; as are rOverM and tOverM.  When G and c are
dimensionful, the dimensionless quantities are

  - (r/c) * (M*G/c^3) * Psi4
  - (r/c) * hdot
  - (r/c) * h / (M*G/c^3)
  - (r/c) / (M*G/c^3)
  - t / (M*G/c^3)

To regain the dimensionful quantities, we simply need to remove the
relevant dimensionful elements (e.g., the r and M factors).

"""

import sys as _sys
from importlib.metadata import version as _version

from . import _GWFrames
from ._GWFrames import *  # noqa: F403
from ._GWFrames import Quaternions, Waveform, PNWaveform

# Make `import GWFrames.Quaternions` and `from GWFrames.Quaternions import ...` work
_sys.modules[__name__ + ".Quaternions"] = Quaternions
from .Quaternions import Quaternion  # noqa: E402

__version__ = _version("GWFrames")


# Pickling.  As with the old SWIG interface, a PNWaveform is pickled as a
# plain Waveform, because its post-Newtonian data cannot be restored.
def _Waveform__getstate__(self):
    return (
        self.HistoryStr(),
        self.T(),
        self.Frame(),
        self.FrameType(),
        self.DataType(),
        self.RIsScaledOut(),
        self.MIsScaledOut(),
        self.LM(),
        self.Data(),
        self.SpinWeight(),
        self.BoostWeight(),
        self.VersionHist(),
    )


def _Waveform__setstate__(self, data):
    self.SetHistory(data[0])
    self.SetTime(data[1])
    self.SetFrame(data[2])
    self.SetFrameType(data[3])
    self.SetDataType(data[4])
    self.SetRIsScaledOut(data[5])
    self.SetMIsScaledOut(data[6])
    self.SetLM(data[7])
    self.SetData(data[8])
    if len(data) > 9:
        self.SetSpinWeight(data[9])
        self.SetBoostWeight(data[10])
        self.SetVersionHist(list(data[11]))


def _Waveform__reduce__(self):
    return (Waveform, (), self.__getstate__())


Waveform.__getstate__ = _Waveform__getstate__
Waveform.__setstate__ = _Waveform__setstate__
Waveform.__reduce__ = _Waveform__reduce__


def GetFileNamePrefix(W):
    return W.DescriptorString() + "_" + W.FrameTypeString() + "_"


Waveform.GetFileNamePrefix = GetFileNamePrefix


def GetLaTeXDataDescription(W):
    LaTeXDataDescription = ""
    if W.RIsScaledOut():
        LaTeXDataDescription = r"r\,"
    if W.MIsScaledOut():
        if W.DataType() in (UnknownDataType, h, Psi2):  # noqa: F405
            LaTeXDataDescription = LaTeXDataDescription + W.DataTypeLaTeXString() + r"/M"
        elif W.DataType() in (hdot, Psi3):  # noqa: F405
            LaTeXDataDescription = LaTeXDataDescription + W.DataTypeLaTeXString()  # hdot is independent of M
        elif W.DataType() == Psi4:  # noqa: F405
            LaTeXDataDescription = LaTeXDataDescription + r"M\," + W.DataTypeLaTeXString()
        elif W.DataType() == Psi1:  # noqa: F405
            LaTeXDataDescription = LaTeXDataDescription + W.DataTypeLaTeXString() + r"/M^2"
        elif W.DataType() == Psi0:  # noqa: F405
            LaTeXDataDescription = LaTeXDataDescription + W.DataTypeLaTeXString() + r"/M^3"
    else:
        LaTeXDataDescription = LaTeXDataDescription + W.DataTypeLaTeXString()
    return LaTeXDataDescription


Waveform.GetLaTeXDataDescription = GetLaTeXDataDescription


def _AddFileNamePrefix(W, FileName):
    """Add a descriptive prefix to FileName"""
    from os.path import basename, dirname, join

    BaseName = W.GetFileNamePrefix() + basename(FileName)
    if not dirname(FileName):
        return BaseName
    return join(dirname(FileName), BaseName)


def _AsStr(value):
    """Return h5py string data as `str`; h5py 3 reads strings as `bytes`"""
    if isinstance(value, bytes):
        return value.decode("utf-8", "replace")
    return str(value)


def OutputToNRAR(W, FileName, FileWriteMode="w"):
    """
    Output the Waveform in NRAR format.

    Note that the FileName is prepended with some descriptive
    information involving the data type and the frame type, such as
    'rhOverM_Corotating_' or 'rMPsi4_Aligned_'.

    """
    import numpy as np
    from h5py import File

    Group = None
    if ".h5" in FileName and not FileName.endswith(".h5"):
        FileName, Group = FileName.split(".h5")
        FileName += ".h5"
    # Add descriptive prefix to FileName
    FileName = _AddFileNamePrefix(W, FileName)
    with File(FileName, FileWriteMode) as F:
        # If we are writing to a group within the file, create it
        G = F.create_group(Group) if Group else F
        # Now write all the data to various groups in the file
        G.attrs["OutputFormatVersion"] = "GWFrames_NRAR"
        G.create_dataset("History.txt", data=W.HistoryStr() + "OutputToNRAR(W, {0})\n".format(FileName))
        G.attrs["FrameType"] = W.FrameType()
        G.attrs["DataType"] = W.DataType()
        G.attrs["RIsScaledOut"] = int(W.RIsScaledOut())
        G.attrs["MIsScaledOut"] = int(W.MIsScaledOut())
        T = W.T()
        for i_m, (ell, m) in enumerate(W.LM()):
            d = W.Data(i_m)
            Data_m = G.create_dataset(
                "Y_l{0}_m{1}.dat".format(ell, m),
                data=np.column_stack((T, d.real, d.imag)),
                compression="gzip",
                shuffle=True,
            )
            Data_m.attrs["ell"] = ell
            Data_m.attrs["m"] = m


Waveform.OutputToNRAR = OutputToNRAR


def OutputToH5(W, FileName, FileWriteMode="w"):
    """
    Output the Waveform with all necessary information, in GWFrames format.

    This function outputs the Waveform object with *all* data, in a
    format that is not compatible with NRAR.  See also `OutputToNRAR`.

    Note that the FileName is prepended with some descriptive
    information involving the data type and the frame type, such as
    'rhOverM_Corotating_' or 'rPsi4_Aligned_'.

    """
    import numpy as np
    from h5py import File, string_dtype

    # Add descriptive prefix to FileName
    FileName = _AddFileNamePrefix(W, FileName)
    with File(FileName, FileWriteMode) as F:
        F.attrs["OutputFormatVersion"] = "GWFrames_v2"
        F.create_dataset("History", data=W.HistoryStr() + "OutputToH5(W, {0})\n".format(FileName))
        F.create_dataset("Time", data=W.T(), compression="gzip", shuffle=True)
        VersionHist = W.VersionHist()
        if len(VersionHist) > 0:
            F.create_dataset("VersionHist", data=np.array(VersionHist, dtype=object), dtype=string_dtype())
        Frame = W.Frame()
        if len(Frame) > 0:
            F.create_dataset("Frame", data=np.array([[r[0], r[1], r[2], r[3]] for r in Frame]))
        else:
            F.create_dataset("Frame", shape=(), dtype=float)
        F.attrs["FrameType"] = W.FrameType()
        F.attrs["DataType"] = W.DataType()
        F.attrs["RIsScaledOut"] = int(W.RIsScaledOut())
        F.attrs["MIsScaledOut"] = int(W.MIsScaledOut())
        Data = F.create_group("Data")
        for i_m, (ell, m) in enumerate(W.LM()):
            Data_m = Data.create_dataset(
                "l{0}_m{1:+}".format(int(ell), int(m)), data=W.Data(i_m), compression="gzip", shuffle=True
            )
            Data_m.attrs["ell"] = ell
            Data_m.attrs["m"] = m


Waveform.OutputToH5 = OutputToH5


def ReadFromH5(FileName):
    """
    Read data from an H5 file, as output by GWFrames.
    """
    import numpy as np
    from h5py import File

    with File(FileName, "r") as f:
        # Initialize the Waveform object
        W = Waveform()
        # Record the filename being read in
        W.AppendHistory("*this = GWFrames.ReadFromH5(FileName='{0}')\n".format(FileName))
        # Add the old history to the new history
        OldHistory = _AsStr(f["History"][()])
        W.AppendHistory("# <previous_history>\n#" + OldHistory.replace("\n", "\n#") + "# </previous_history>\n")
        # Get the time data
        W.SetTime(f["Time"][()])
        # Get the version history, if there is one
        if "VersionHist" in f:
            W.SetVersionHist([(_AsStr(a), _AsStr(b)) for a, b in f["VersionHist"][()]])
        # Get the frame data, if there is any
        if "Frame" in f and f["Frame"].shape != ():
            W.SetFrame(f["Frame"][()])
        # Get the descriptive items
        try:
            W.SetFrameType(int(f.attrs["FrameType"]))
            W.SetDataType(int(f.attrs["DataType"]))
            W.SetRIsScaledOut(bool(f.attrs["RIsScaledOut"]))
            W.SetMIsScaledOut(bool(f.attrs["MIsScaledOut"]))
        except KeyError:
            print(
                "\nWarning: FrameType, DataType, RIsScaledOut, and/or MIsScaledOut "
                "were not found in '{0}'.\n".format(FileName)
            )
        # Get list of data sets and the LM data, sorted by LM
        ModeData = list(f["Data"])
        LMlist = [[int(f["Data"][Data_m].attrs["ell"]), int(f["Data"][Data_m].attrs["m"])] for Data_m in ModeData]
        SortedIndices = sorted(range(len(LMlist)), key=lambda i: LMlist[i])
        Data = np.empty((len(LMlist), W.NTimes()), dtype=complex)
        for i_m, i_sorted in enumerate(SortedIndices):
            Data[i_m] = f["Data"][ModeData[i_sorted]][()]
        W.SetLM([LMlist[i] for i in SortedIndices])
        W.SetData(Data)
    return W


def MonotonicIndices(T, MinTimeStep=1.0e-5):
    """
    Given an array of times, return the indices that make the array strictly monotonic.
    """
    import numpy as np

    Ind = np.arange(len(T))
    Size = len(Ind)
    i = 1
    while i < Size:
        if T[Ind[i]] <= T[Ind[i - 1]] + MinTimeStep:
            j = 0
            while T[Ind[j]] + MinTimeStep < T[Ind[i]]:
                j += 1
            # erase data from j (inclusive) to i (exclusive)
            Ind = np.delete(Ind, range(j, i))
            Size = len(Ind)
            i = j - 1
        i += 1
    return Ind


def ReadFromNRAR(FileName):
    """
    Read data from an H5 file in NRAR format.
    """
    import re

    import h5py
    import numpy as np

    YlmRegex = re.compile(r"""Y_l(?P<L>[0-9]+)_m(?P<M>[-+0-9]+)\.dat""")
    # Initialize the Waveform object
    W = Waveform()
    # Record the filename being read in
    W.AppendHistory("*this = GWFrames.ReadFromNRAR(FileName='{0}')\n".format(FileName))
    try:
        FileName, RootGroup = FileName.rsplit(".h5", 1)
        FileName += ".h5"
    except ValueError:
        RootGroup = ""  # FileName is just a file, not a group in a file
    with h5py.File(FileName, "r") as f_h5:
        f = f_h5[RootGroup] if RootGroup else f_h5
        # Add the old history to the new history, if found
        if "History.txt" in f:
            OldHistory = _AsStr(f["History.txt"][()])
            W.AppendHistory("# <previous_history>\n#" + OldHistory.replace("\n", "\n#") + "# </previous_history>\n")
        # Get the frame data, if there is any
        if "Frame" in f:
            W.SetFrame(f["Frame"][()])
        # Get the descriptive items
        try:
            W.SetFrameType(int(f.attrs["FrameType"]))
            W.SetDataType(int(f.attrs["DataType"]))
            W.SetRIsScaledOut(bool(f.attrs["RIsScaledOut"]))
            W.SetMIsScaledOut(bool(f.attrs["MIsScaledOut"]))
        except KeyError:
            print(
                "\nWarning: FrameType, DataType, RIsScaledOut, and/or MIsScaledOut "
                "were not found in '{0}'.\n".format(FileName)
                + "Using defaults.  You may want to re-set them manually.\n\n"
            )
            W.SetFrameType(Inertial)  # noqa: F405
            W.SetRIsScaledOut(True)
            W.SetMIsScaledOut(True)
            lower = FileName.lower()
            for name, DataType in [("psi0", Psi0), ("psi1", Psi1), ("psi2", Psi2), ("psi3", Psi3)]:  # noqa: F405
                if name in lower:
                    W.SetDataType(DataType)
            if "psi4" in lower:
                W.SetDataType(Psi4)  # noqa: F405
            elif "hdot" in lower:
                W.SetDataType(hdot)  # noqa: F405
            elif "h" in lower:
                W.SetDataType(h)  # noqa: F405
        # Get the names of all the datasets in the h5 file, and check for matches
        YLMdata = [DataSet for DataSet in list(f) if YlmRegex.search(DataSet)]
        if len(YLMdata) == 0:
            raise ValueError("Couldn't understand dataset names in '{0}'.".format(FileName))
        # Sort the dataset names by increasing ell, then increasing m
        def LMof(DataSet):
            match = YlmRegex.search(DataSet)
            return [int(match.group("L")), int(match.group("M"))]

        YLMdata = sorted(YLMdata, key=LMof)
        LM = [LMof(DataSet) for DataSet in YLMdata]
        # Get the time data (assuming all are equal)
        Wdata = f[YLMdata[0]]
        NTimes = Wdata.shape[0]
        T = Wdata[:, 0]
        # Loop through, getting each mode
        Data = np.empty((len(LM), NTimes), dtype=complex)
        for i_m, DataSet in enumerate(YLMdata):
            if f[DataSet].shape[0] != NTimes:
                raise ValueError(
                    "The number of time steps in this dataset should be {0}; ".format(NTimes)
                    + "it is {0} in '{1}'.".format(f[DataSet].shape[0], DataSet)
                )
            Data[i_m] = f[DataSet][:, 1] + 1j * f[DataSet][:, 2]
    # Make sure time is monotonic and set the data
    Indices = MonotonicIndices(T)
    W.SetTime(T[Indices])
    W.SetLM(LM)
    W.SetData(Data[:, Indices])
    return W
