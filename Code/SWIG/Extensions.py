/// Add utility functions that are specific to python.  Note that
/// these are defined in the GWFrames namespace.
%insert("python") %{

# The following are necessary to allow us to attach new attributes
# (methods) to `Waveform` objects.  This is required because we are
# calling SWIG with the `-builtin` option, which makes it a little
# more complicated to add new attributes.  This defeats the whole
# point of `-builtin` for this particular object.  But it allows other
# objects (especially `Quaternion`s) to be implemented more cleanly
# and efficiently.  Note also that the `_Waveform` was given that name
# above by the line
#    %rename(_Waveform) Waveform;
# in `Waveforms.i`.

class _MetaWaveform(type(_Waveform)):
    pass

class Waveform(_Waveform):
    """Object storing data and other information for a single waveform.

    Fundamental object encapsulating waveform data, such as time,
    (l,m) information, and complex data.

    This object provides the main user interface for this collection
    of code. The various methods for this class are intended to
    provide all manipulations that might be necessary in the course of
    waveform analysis.

    For an extensive listing of all methods available for Waveform
    objects, run `help(GWFrames.Waveform)`.  Documentation is also
    available for each such method individually.

    In python, the preferred ways to construct a new Waveform from a
    data file is to use either `GWFrames.ReadFromH5` or
    `GWFrames.ReadFromNRAR`.  The copy constructor also works.

    The Waveform's `frame` data records the rotors needed to rotate
    the standard (x,y,z) basis into the (X,Y,Z) basis with respect to
    which the Waveform's mode `data` are decomposed.

    Note on Waveform Types:
    In any system, h -- being strain -- should be dimensionless.
    When G=c=1, the dimensionless quantities are rMPsi4, rhdot, and
    rhOverM; as are rOverM and tOverM.  When G and c are
    dimensionful, the dimensionless quantities are
        -  (r/c) * (M*G/c^3) * Psi4
        -  (r/c) * hdot
        -  (r/c) * h / (M*G/c^3)
        -  (r/c) / (M*G/c^3)
        -  t / (M*G/c^3)
    To regain the dimensionful quantities, we simply need to remove
    the relevant dimensionful elements (e.g., the r and M factors).

    Member data
    -----------
      int spinweight
      int boostweight
      std::stringstream history
      std::vector<double> t
      std::vector<Quaternions::Quaternion> frame
      WaveformFrameType frameType
      WaveformDataType dataType
      bool rIsScaledOut
      bool mIsScaledOut
      std::vector<std::vector<int> > lm
      MatrixC data

    """
    __metaclass__ = _MetaWaveform

_WaveformReturners = ['CopyWithoutData', 'SliceOfTimeIndices', 'SliceOfTimeIndicesWithEll2', 'SliceOfTimeIndicesWithoutModes', 'SliceOfTimes', 'SliceOfTimesWithEll2',
                      'SliceOfTimesWithoutModes', 'Interpolate', 'InterpolateInPlace', 'DropTimesOutside', 'DropEllModes', 'KeepOnlyEllModes', 'KeepOnlyEll2',
                      'SetSpinWeight', 'SetBoostWeight', 'AppendHistory', 'SetHistory', 'SetT', 'SetTime', 'SetFrame', 'SetFrameType', 'SetDataType', 'SetRIsScaledOut',
                      'SetMIsScaledOut', 'SetLM', 'SetData', 'SetData', 'ResizeData', 'Differentiate', 'RotatePhysicalSystem',
                      'RotatePhysicalSystem', 'RotateDecompositionBasis', 'RotateDecompositionBasis', 'TransformToCoprecessingFrame', 'TransformToAngularVelocityFrame',
                      'TransformToCorotatingFrame', 'TransformToInertialFrame', 'AlignDecompositionFrameToModes', 'Compare', 'Hybridize', 'Translate']

class _MetaPNWaveform(type(_PNWaveform)):
    pass

class PNWaveform(_PNWaveform):
    """Fundamental object creating precessing PN waveform.

    This object is a subclass of the GWFrames::Waveform object.  In
    addition to the data stored in a Waveform, this stores the two
    spins chi1 and chi2; the orbital angular-velocity vectors
    Omega_orb, Omega_prec; the PN angular momentum L; and the orbital
    phase Phi_orb.  All such vectors are with respect to the inertial
    frame.  We also have the usual frame data, which is initially the
    co-orbital frame.  Various methods also exist for retrieving the
    vectors, their magnitudes, and their normalized versions.

    For an extensive listing of all methods available for PNWaveform
    objects, run `help(GWFrames.PNWaveform)`.  Documentation is also
    available for each such method individually.

    Member data
    -----------
      [All of GWFrames.Waveform members and...]
      std::vector<std::vector<double> > mchi1
      std::vector<std::vector<double> > mchi2
      std::vector<std::vector<double> > mOmega_orb
      std::vector<std::vector<double> > mOmega_prec
      std::vector<std::vector<double> > mL
      std::vector<double> mPhi_orb

    Constructor
    -----------
    Parameters:
      Approximant: 'TaylorT1'|'TaylorT4'|'TaylorT5'
      delta: Normalized BH mass difference (M1-M2)/(M1+M2)
      chi1_i: Initial dimensionless spin vector of BH1
      chi2_i: Initial dimensionless spin vector of BH2
      Omega_orb_i: Initial orbital angular frequency
      Omega_orb_0: Earliest orbital angular frequency to compute (default: Omega_orb_i)
      R_frame_i: Initial rotation of the binary (default: No rotation)
      MinStepsPerOrbit: Minimum number of time steps at which to evaluate (default: 32)
      PNWaveformModeOrder: PN order at which to compute waveform modes (default: 3.5)
      PNOrbitalEvolutionOrder: PN order at which to compute orbital evolution (default: 4.0)

    [There is also a copy constructor.]

    The PN system is defined with respect to an inertial basis
    (x,y,z).  The input spin vectors must be defined with respect to
    this basis.  A new basis is (X,Y,Z) is created by rotating (x,y,z)
    by `R_frame_i`.  The black-hole positions are initialized with
    respect to this basis, with BH1 along the positive X axis and BH2
    along the negative X axis, and the orbital angular velocity along
    the positive Z axis having magnitude `Omega_orb_i`.  Note that the
    optional parameter `Omega_orb_0` may be given, in which case the
    PN system is also evolved backwards to that point.

    The TaylorTn system is first integrated to compute the dynamics of
    the binary.  The evolved spin vectors `chi1` and `chi2`, orbital
    angular-velocity vector `Omega_orb`, and co-orbital frame
    `R_frame` itself are stored, in addition to a few other niceties.
    Once this evolution is completed, the information is used to
    construct the waveform in the co-orbital frame.  (The waveform
    will be nearly co-rotating, but not quite.)

    To get the PNWaveform in an inertial frame, you must first apply
    the method TransformToInertialFrame().  To get the PNWaveform in a
    co-rotating frame, you first apply the method
    TransformToInertialFrame(), then the method
    TransformToCororatingFrame().

    """
    __metaclass__ = _MetaPNWaveform



# Now, we just make sure that any Waveform or PNWaveform member that
# returns a _Waveform or _PNWaveform gets wrapped to return the
# appropriate derived object
from functools import wraps
def _ObliterateUnderscores(func):
    @wraps(func, ('__name__', '__doc__'))
    def decorator(self, *args, **kwargs):
        def Derivedify(obj):
            typestr = type(obj).__name__
            if typestr=='_Waveform' or typestr=='_PNWaveform':
                return eval(typestr[1:])(obj)
            else:
                return obj
        ret = func(self, *args, **kwargs)
        if isinstance(ret, tuple):
            for i in range(len(ret)):
                ret[i] = Derivedify(ret[i])
            return ret
        else:
            return Derivedify(ret)
    return decorator
for _WaveformReturner in _WaveformReturners:
    setattr(Waveform, _WaveformReturner, _ObliterateUnderscores(getattr(Waveform, _WaveformReturner)))




def GetFileNamePrefix(W) :
    return W.DescriptorString() + '_' + W.FrameTypeString() + '_'
Waveform.GetFileNamePrefix = GetFileNamePrefix
PNWaveform.GetFileNamePrefix = GetFileNamePrefix

def GetLaTeXDataDescription(W) :
    from GWFrames import UnknownDataType, h, hdot, Psi4, Psi3, Psi2, Psi1, Psi0
    LaTeXDataDescription = ''
    if(W.RIsScaledOut()) :
        LaTeXDataDescription = r'r\,'
    if(W.MIsScaledOut()) :
        if(W.DataType()==UnknownDataType or W.DataType()==h or W.DataType()==Psi2) :
            LaTeXDataDescription = LaTeXDataDescription + W.DataTypeLaTeXString() + r'/M'
        elif(W.DataType()==hdot or W.DataType()==Psi3) :
            LaTeXDataDescription = LaTeXDataDescription + W.DataTypeLaTeXString() # hdot is independent of M
        elif(W.DataType()==Psi4) :
            LaTeXDataDescription = LaTeXDataDescription + r'M\,' + W.DataTypeLaTeXString()
        elif(W.DataType()==Psi1) :
            LaTeXDataDescription = LaTeXDataDescription + W.DataTypeLaTeXString() + r'/M^2'
        elif(W.DataType()==Psi0) :
            LaTeXDataDescription = LaTeXDataDescription + W.DataTypeLaTeXString() + r'/M^3'
    else :
        LaTeXDataDescription = LaTeXDataDescription + W.DataTypeLaTeXString()
    return LaTeXDataDescription
Waveform.GetLaTeXDataDescription = GetLaTeXDataDescription
PNWaveform.GetLaTeXDataDescription = GetLaTeXDataDescription

def __AddFileNamePrefix(W, FileName):
    """Add a descriptive prefix to FileName"""
    from os.path import basename, dirname
    BaseName = W.GetFileNamePrefix() + basename(FileName)
    if not dirname(FileName):
        FileName = BaseName
    else:
        FileName = dirname(FileName) + '/' + BaseName
    return FileName

def OutputToNRAR(W, FileName, FileWriteMode='w') :
    """
    Output the Waveform in NRAR format.

    Note that the FileName is prepended with some descriptive
    information involving the data type and the frame type, such as
    'rhOverM_Corotating_' or 'rMPsi4_Aligned_'.

    """
    from h5py import File
    from GWFrames import UnknownDataType, h, hdot, Psi4
    Group = None
    if('.h5' in FileName and not FileName.endswith('.h5')) :
        FileName,Group = FileName.split('.h5')
        FileName += '.h5'
    # Add descriptive prefix to FileName
    FileName = __AddFileNamePrefix(W, FileName)
    # Open the file for output
    try :
        F = File(FileName, FileWriteMode)
    except IOError : # If that did not work...
        print("OutputToNRAR was unable to open the file '{0}'.\n\n".format(FileName))
        raise # re-raise the exception after the informative message above
    try :
        # If we are writing to a group within the file, create it
        if(Group) :
            G = F.create_group(Group)
        else :
            G = F
        # Now write all the data to various groups in the file
        G.attrs['OutputFormatVersion'] = 'GWFrames_NRAR'
        G.create_dataset("History.txt", data = W.HistoryStr() + 'OutputToNRAR(W, {0})\n'.format(FileName))
        G.attrs['FrameType'] = W.FrameType()
        G.attrs['DataType'] = W.DataType()
        G.attrs['RIsScaledOut'] = int(W.RIsScaledOut())
        G.attrs['MIsScaledOut'] = int(W.MIsScaledOut())
        for i_m in range(W.NModes()) :
            ell,m = W.LM()[i_m]
            Data_m = G.create_dataset("Y_l{0}_m{1}.dat".format(ell, m), data=[[t, d.real, d.imag] for t,d in zip(W.T(),W.Data(i_m))],
                                      compression="gzip", shuffle=True)
            Data_m.attrs['ell'] = ell
            Data_m.attrs['m'] = m
    finally : # Use `finally` to make sure this happens:
        # Close the file and we are done
        F.close()
Waveform.OutputToNRAR = OutputToNRAR
PNWaveform.OutputToNRAR = OutputToNRAR

def OutputToH5(W, FileName, FileWriteMode='w') :
    """
    Output the Waveform with all necessary information, in GWFrames format.

    This function outputs the Waveform object with *all* data, in a
    format that is not compatible with NRAR.  See also `OutputToNRAR`.

    Note that the FileName is prepended with some descriptive
    information involving the data type and the frame type, such as
    'rhOverM_Corotating_' or 'rPsi4_Aligned_'.

    """
    from h5py import File, special_dtype
    from GWFrames import UnknownDataType, h, hdot, Psi4
    # Add descriptive prefix to FileName
    FileName = __AddFileNamePrefix(W, FileName)
    # Open the file for output
    try :
        F = File(FileName, FileWriteMode)
    except IOError : # If that did not work...
        print("OutputToH5 was unable to open the file '{0}'.\n\n".format(FileName))
        raise # re-raise the exception after the informative message above
    try :
        # Now write all the data to various groups in the file
        F.attrs['OutputFormatVersion'] = 'GWFrames_v2'
        F.create_dataset("History", data = W.HistoryStr() + 'OutputToH5(W, {0})\n'.format(FileName))
        F.create_dataset("Time", data=W.T().tolist(), compression="gzip", shuffle=True)
        if(len(W.VersionHist())>0) :
            F.create_dataset("VersionHist", (1,2), data = W.VersionHist(), dtype=special_dtype(vlen=bytes))
        if(len(W.Frame())>0) :
            F.create_dataset("Frame", data=[[r[0], r[1], r[2], r[3]] for r in W.Frame()])
        else :
            F.create_dataset("Frame", shape=())
        F.attrs['FrameType'] = W.FrameType()
        F.attrs['DataType'] = W.DataType()
        F.attrs['RIsScaledOut'] = int(W.RIsScaledOut())
        F.attrs['MIsScaledOut'] = int(W.MIsScaledOut())
        Data = F.create_group("Data")
        for i_m in range(W.NModes()) :
            ell,m = W.LM()[i_m]
            Data_m = Data.create_dataset("l{0}_m{1:+}".format(int(ell), int(m)), data=W.Data(i_m),
                                         compression="gzip", shuffle=True)
            Data_m.attrs['ell'] = ell
            Data_m.attrs['m'] = m
    finally : # Use `finally` to make sure this happens:
        # Close the file and we are done
        F.close()
Waveform.OutputToH5 = OutputToH5
PNWaveform.OutputToH5 = OutputToH5

def ReadFromH5(FileName) :
    """
    Read data from an H5 file, as output by GWFrames.
    """
    from h5py import File
    from GWFrames import Waveform
    from Quaternions import Quaternion
    from numpy import empty
    try :
        f = File(FileName, 'r')
    except IOError :
        print("ReadFromH5 could not open the file '{0}'\n\n".format(FileName))
        raise
    try :
        # Initialize the Waveform object
        W = Waveform()
        # Record the filename being read in
        W.AppendHistory(str("*this = GWFrames.ReadFromH5(FileName='{0}')\n".format(FileName)))
        # Add the old history to the new history
        W.AppendHistory(str("# <previous_history>\n#" + f['History'][()].replace('\n','\n#') + "# </previous_history>\n"))
        # Get the time data
        W.SetTime(f['Time'])
        # Get the frame data, converting to GWFrame.Quaternion objects
        try :
            W.SetFrame([Quaternion(r) for r in f['Frame']])
        except TypeError :
            pass # There was no frame
        # Get the descriptive items
        try :
            W.SetFrameType(int(f.attrs['FrameType']))
            W.SetDataType(int(f.attrs['DataType']))
            W.SetRIsScaledOut(bool(f.attrs['RIsScaledOut']))
            W.SetMIsScaledOut(bool(f.attrs['MIsScaledOut']))
        except KeyError :
            print("\nWarning: FrameType, DataType, RIsScaledOut, and/or MIsScaledOut were not found in '{0}'.\n".format(FileName))
        # Get list of data sets and the LM data (unsorted)
        ModeData = list(f['Data'])
        LMlist = [[f['Data'][Data_m].attrs['ell'], f['Data'][Data_m].attrs['m']] for Data_m in ModeData]
        NModes = len(LMlist)
        # Get the order of the sort by LM
        SortedIndices = sorted(range(NModes),key=lambda i : LMlist[i])
        # Initialize the data set and LM set
        Data = empty((NModes, W.NTimes()), dtype='complex128')
        LM = empty((NModes, 2), dtype='int')
        # Loop through the modes, storing them in sorted order
        for i_m in range(NModes) :
            Data[i_m] = f['Data'][ModeData[SortedIndices[i_m]]]
            LM[i_m] = LMlist[SortedIndices[i_m]]
        # Now add these data to the Waveform object
        W.SetLM(LM.tolist())
        W.SetData(Data)
    except KeyError :
        print("This H5 file appears to have not stored all the required information.\n\n")
        raise # Re-raise the exception after adding our information
    finally : # Use `finally` to make sure this happens:
        f.close()
    return W

def MonotonicIndices(T, MinTimeStep=1.e-5) :
    """
    Given an array of times, return the indices that make the array strictly monotonic.
    """
    import numpy
    Ind = range(len(T))
    Size = len(Ind)
    i=1
    while(i<Size) :
        if(T[Ind[i]]<=T[Ind[i-1]]+MinTimeStep) :
            j=0
            while(T[Ind[j]]+MinTimeStep<T[Ind[i]]) :
                j += 1
            # erase data from j (inclusive) to i (exclusive)
            Ind = numpy.delete(Ind, range(j,i))
            Size = len(Ind)
            i = j-1
        i += 1
    return Ind

def ReadFromNRAR(FileName) :
    """
    Read data from an H5 file in NRAR format.
    """
    import re
    import h5py
    from Quaternions import Quaternion
    import numpy
    YlmRegex = re.compile(r"""Y_l(?P<L>[0-9]+)_m(?P<M>[-+0-9]+)\.dat""")
    # Initialize the Waveform object
    W = Waveform()
    # Record the filename being read in
    W.AppendHistory(str("*this = GWFrames.ReadFromNRAR(FileName='{0}')\n".format(FileName)))
    try :
        FileName, RootGroup = FileName.rsplit('.h5', 1)
        FileName += '.h5'
    except ValueError :
        RootGroup = '' # FileName is just a file, not a group in a file
    try :
        f_h5 = h5py.File(FileName, 'r')
    except IOError :
        print("ReadFromNRAR could not open the file '{0}'\n\n".format(FileName))
        raise
    if(RootGroup) :
        f = f_h5[RootGroup]
    else :
        f = f_h5
    try :
        try :
            # Add the old history to the new history, if found
            OldHistory = f['History.txt'][()].encode('ascii','ignore')
            W.AppendHistory(str("# <previous_history>\n#" + OldHistory.replace('\n','\n#') + "# </previous_history>\n"))
        except KeyError :
            pass # Did not find a history
        # Get the frame data, converting to Quaternions.Quaternion objects
        try :
            W.SetFrame([Quaternion(r) for r in f['Frame']])
        except KeyError :
            pass # There was no frame data
        # Get the descriptive items
        try :
            W.SetFrameType(int(f.attrs['FrameType']))
            W.SetDataType(int(f.attrs['DataType']))
            W.SetRIsScaledOut(bool(f.attrs['RIsScaledOut']))
            W.SetMIsScaledOut(bool(f.attrs['MIsScaledOut']))
        except KeyError :
            print("\nWarning: FrameType, DataType, RIsScaledOut, and/or MIsScaledOut were not found in '{0}'.\n".format(FileName)+
                  "Using defaults.  You may want to re-set them manually.\n\n")
            W.SetFrameType(1)
            W.SetRIsScaledOut(True)
            W.SetMIsScaledOut(True)
            if('psi0' in FileName.lower()) :
                W.SetDataType(7)
            if('psi1' in FileName.lower()) :
                W.SetDataType(6)
            if('psi2' in FileName.lower()) :
                W.SetDataType(5)
            if('psi3' in FileName.lower()) :
                W.SetDataType(4)
            if('psi4' in FileName.lower()) :
                W.SetDataType(3)
            elif('hdot' in FileName.lower()) :
                W.SetDataType(2)
            elif('h' in FileName.lower()) :
                W.SetDataType(1)
        # Get the names of all the datasets in the h5 file, and check for matches
        YLMdata = [DataSet for DataSet in list(f) for m in [YlmRegex.search(DataSet)] if m]
        if(len(YLMdata)==0) :
            raise ValueError("Couldn't understand dataset names in '{0}'.".format(FileName))
        # Sort the dataset names by increasing ell, then increasing m
        YLMdata = sorted(YLMdata, key=lambda DataSet : [int(YlmRegex.search(DataSet).group('L')), int(YlmRegex.search(DataSet).group('M'))])
        # List just the ell and m numbers
        LM = sorted([[int(m.group('L')), int(m.group('M'))] for DataSet in YLMdata for m in [YlmRegex.search(DataSet)] if m])
        NModes = len(LM)
        # Get the time data (assuming all are equal)
        Wdata = f[YLMdata[0]]
        NTimes = Wdata.shape[0]
        T = Wdata[:,0]
        # Set up storage
        Re = numpy.empty((NModes, NTimes))
        Im = numpy.empty((NModes, NTimes))
        m = 0
        # Loop through, getting each mode
        for DataSet in YLMdata :
            if( not (f[DataSet].shape[0]==NTimes) ) :
                raise ValueError("The number of time steps in this dataset should be {0}; ".format(NTimes) +
                                 "it is {0} in '{1}'.".format(f[DataSet].shape[0], DataSet))
            Re[m,:] = f[DataSet][:,1]
            Im[m,:] = f[DataSet][:,2]
            m += 1
        # Make sure time is monotonic and set the data
        Indices = MonotonicIndices(T)
        BadIndices = numpy.setdiff1d(range(len(T)), Indices)
        W.SetTime(T[Indices])
        W.SetLM(LM)
        W.SetData(numpy.delete(Re, BadIndices, 1)+1j*numpy.delete(Im, BadIndices, 1))
    except KeyError :
        print("This H5 file appears to have not stored all the required information.\n\n")
        raise # Re-raise the exception after adding our information
    finally : # Use `finally` to make sure this happens:
        f_h5.close()
    return W

%}
