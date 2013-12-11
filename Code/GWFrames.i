// -*- c++ -*-

// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details


%module GWFrames

 // Quiet warnings about overloaded operators being ignored.
#pragma SWIG nowarn=362,389,401,509
%include <typemaps.i>
%include <stl.i>

// %{
//   #include "/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages/numpy/core/include/numpy/arrayobject.h"
//   #include <vector>
// %}

// %init %{
//   import_array();
// %}

// %typemap(out) std::vector<int> {
//   npy_intp result_size = $1.size();
//   npy_intp dims[1] = { result_size };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_INT);
//   int* dat = (int*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { dat[i] = $1[i]; }
//   $result = PyArray_Return(npy_arr);
// }
// %typemap(out) std::vector<int>& {
//   npy_intp result_size = $1->size();
//   npy_intp dims[1] = { result_size };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_INT);
//   int* dat = (int*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { dat[i] = (*$1)[i]; }
//   $result = PyArray_Return(npy_arr);
// }
// %typemap(out) std::vector<std::vector<int> >& {
//   npy_intp result_size = $1->size();
//   npy_intp result_size2 = (result_size>0 ? (*$1)[0].size() : 0);
//   npy_intp dims[2] = { result_size, result_size2 };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_INT);
//   int* dat = (int*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { for (size_t j = 0; j < result_size2; ++j) { dat[i*result_size2+j] = (*$1)[i][j]; } }
//   $result = PyArray_Return(npy_arr);
// }

// %typemap(out) std::vector<double> {
//   npy_intp result_size = $1.size();
//   npy_intp dims[1] = { result_size };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
//   double* dat = (double*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { dat[i] = $1[i]; }
//   $result = PyArray_Return(npy_arr);
// }
// %typemap(out) std::vector<double>& {
//   npy_intp result_size = $1->size();
//   npy_intp dims[1] = { result_size };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
//   double* dat = (double*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { dat[i] = (*$1)[i]; }
//   $result = PyArray_Return(npy_arr);
// }
// %typemap(out) std::vector<std::vector<double> >& {
//   npy_intp result_size = $1->size();
//   npy_intp result_size2 = (result_size>0 ? (*$1)[0].size() : 0);
//   npy_intp dims[2] = { result_size, result_size2 };
//   PyArrayObject* npy_arr = (PyArrayObject*)PyArray_SimpleNew(2, dims, NPY_DOUBLE);
//   double* dat = (double*) PyArray_DATA(npy_arr);
//   for (size_t i = 0; i < result_size; ++i) { for (size_t j = 0; j < result_size2; ++j) { dat[i*result_size2+j] = (*$1)[i][j]; } }
//   $result = PyArray_Return(npy_arr);
// }

%include "../Docs/GWFrames_Doc.i"

///////////////////////////////////
//// Handle exceptions cleanly ////
///////////////////////////////////
%exception {
  try {
    $action;
  } catch(int i) {
    if(i==0) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Not yet implemented.");
    } else if(i==1) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Index out of bounds.");
    } else if(i==2) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Infinitely many solutions.");
    } else if(i==3) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Input vector size mismatch.");
    } else if(i==4) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Cannot extrapolate quaternions.");
    } else if(i==5) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Matrix size mismatch.");
    } else if(i==6) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Matrix size is assumed to be 3x3 in this function.");
    } else if(i==7) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Quaternion constructor's vector size not understood; should be 3 or 4.");
    } else if(i==8) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Waveform is missing requested l,m component.");
    } else if(i==9) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Bad file name.");
    } else if(i==10) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Not enough points to take a derivative.");
    } else if(i==11) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Empty intersection requested.");
    } else if(i==12) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Failed system call.");
    } else if(i==13) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Wrong FrameType for this operation.  Maybe you forgot to `SetFrameType`?");
    } else if(i==14) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: GSL failed.");
    } else  {
      PyErr_SetString(PyExc_RuntimeError, "Unknown exception");
    }
    return NULL;
  }
}

/////////////////////////////////////////////////
//// These will be needed by the c++ wrapper ////
/////////////////////////////////////////////////
%{
  #include <iostream>
  #include <string>
  #include <sstream>
  #include <iomanip>
  #include <complex>
  #include "Quaternions.hpp"
  #include "Utilities.hpp"
  #include "Waveforms.hpp"
  #include "PNWaveforms.hpp"
%}


//////////////////////////////////////////////////////////////////////
//// The following translates between c++ and python types nicely ////
//////////////////////////////////////////////////////////////////////
//// This lets me use numpy.array in the code below
%pythoncode %{
  import numpy;
  %}
//// Make sure std::strings are dealt with appropriately
%include <std_string.i>
//// Make sure std::complex numbers are dealt with appropriately
%include <std_complex.i>
// namespace std {
//   %template(complexd) complex<double>; // Don't use this line!!!
// };
//// Make sure std::vectors are dealt with appropriately
%include <std_vector.i>
namespace GWFrames {
  class Quaternion;
 };
namespace std {
  %template(vectori) vector<int>;
  %template(vectorvectori) vector<vector<int> >;
  %template(vectord) vector<double>;
  %template(vectorvectord) vector<vector<double> >;
  %template(vectorc) vector<std::complex<double> >;
  %template(vectorvectorc) vector<vector<std::complex<double> > >;
  %template(vectorq) vector<GWFrames::Quaternion>;
  %template(vectors) vector<string>;
  %template(vectorvectors) vector<vector<std::string> >;
};


/////////////////////////////////////
//// Import the quaternion class ////
/////////////////////////////////////
%ignore GWFrames::Quaternion::operator=;
%rename(__getitem__) GWFrames::Quaternion::operator [](unsigned int const) const;
%rename(__setitem__) GWFrames::Quaternion::operator [](unsigned int const);
%include "Quaternions.hpp"
%extend GWFrames::Quaternion {
  // This function is called when printing a Quaternion object
  const char* __str__() {
    std::stringstream S;
    S << std::setprecision(14) << "["
      << $self->operator[](0) << ", "
      << $self->operator[](1) << ", "
      << $self->operator[](2) << ", " 
      << $self->operator[](3) << "]";
    const std::string& tmp = S.str();
    const char* cstr = tmp.c_str();
    return cstr;
  }
  // This prints the Quaternion nicely at the prompt and allows nucer manipulations
  %pythoncode{
    def __repr__(self):
        return 'Quaternion('+repr(self[0])+', '+repr(self[1])+', '+repr(self[2])+', '+repr(self[3])+')'
    def __pow__(self, P) :
        return self.pow(P)
    __radd__ = __add__
    def __rsub__(self, t) :
        return -self+t
    __rmul__ = __mul__
    def __rdiv__(self, t) :
        return self.inverse()*t
  };
 };



%ignore GWFrames::Matrix::operator=;
%ignore GWFrames::Matrix::operator[];
%include "Utilities.hpp"
namespace std {
  %template(vectorM) vector<GWFrames::Matrix>;
};

%extend GWFrames::Matrix {
  // Print the Matrix nicely at the prompt
  %pythoncode{
    def __repr__(self) :
        return "".join(
            [  'Matrix(['+repr([self(0,c) for c in range(self.ncols())])+',\n']
            + ['        '+repr([self(r,c) for c in range(self.ncols())])+',\n' for r in range(1,self.nrows()-1)]
            + ['        '+repr([self(r,c) for c in range(self.ncols()) for r in [self.nrows()-1]])+'])'] )
  };
 };



////////////////////////////////////////////////
//// Prepare to read in the Waveforms class ////
////////////////////////////////////////////////
//// Ignore this, as neither const nor non-const will work (copy constructor issues?)
%ignore GWFrames::Waveforms::operator[];
%rename(__getitem__) GWFrames::Waveforms::operator[] const;


////////////////////////////////////
//// Read in the Waveform class ////
////////////////////////////////////
//// Ignore things that don't translate well...
%ignore operator<<;
%ignore GWFrames::Waveform::operator=;
//// Allow us to get Quaternions returned naturally when passed by reference
%typemap(in,numinputs=0) GWFrames::Quaternion& OUTPUT (GWFrames::Quaternion temp) { $1 = &temp; }
%typemap(argout) GWFrames::Quaternion& OUTPUT {
  $result = SWIG_Python_AppendOutput($result, SWIG_NewPointerObj((new GWFrames::Quaternion(static_cast< const GWFrames::Quaternion& >(temp$argnum))), SWIGTYPE_p_GWFrames__Quaternion, SWIG_POINTER_OWN |  0 ));
}
//// These will convert the output data to numpy.ndarray for easier use
%feature("pythonappend") GWFrames::Waveform::T() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Frame() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::LM() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Re() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Im() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Abs() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Arg() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::ArgUnwrapped() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Data() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::DataDot() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::Norm() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::LdtVector() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::LLMatrix() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::OShaughnessyEtAlVector() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::SchmidtEtAlVector() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::AngularVelocityVector() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::CorotatingFrame() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::PNEquivalentOrbitalAV() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::Waveform::PNEquivalentPrecessionalAV() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
//// Allow us to extract the outputs naturally in python
%apply double& OUTPUT { double& deltat };
%apply GWFrames::Quaternion& OUTPUT { GWFrames::Quaternion& R_delta };
//// Parse the header file to generate wrappers
%include "Waveforms.hpp"
//// Make sure vectors of Waveform are understood
namespace std {
  %template(vectorW) vector<GWFrames::Waveform>;
};
//// Make any additions to the Waveform class here
%extend GWFrames::Waveform {
  //// This function is called when printing the Waveform object
  std::string __str__() {
    std::stringstream S;
    S << ($self->HistoryStr()) << "###\n"
      << "### # In python:\n"
      << "### import GWFrames\n"
      << "### print(this)" << std::endl << std::setprecision(14);
    for(unsigned int t=0; t<$self->NTimes(); ++t) {
      S << $self->T(t) << " ";
      for(unsigned int mode=0; mode<$self->NModes(); ++mode) {
  	S << $self->Re(mode, t) << " " << $self->Im(mode, t) << " ";
      }
      S << std::endl;
    }
    return S.str();
  }
  //// Allow Waveform objects to be pickled
  %insert("python") %{
    def __getstate__(self) :
      return (self.HistoryStr(),
  	      self.T(),
  	      self.Frame(),  
  	      self.FrameType(),
  	      self.DataType(),
  	      self.RIsScaledOut(),
  	      self.MIsScaledOut(),
	      self.LM(),
  	      self.Data()
  	      )
    __safe_for_unpickling__ = True
    def __reduce__(self) :
        return (Waveform, (), self.__getstate__())
    def __setstate__(self, data) :
      self.SetHistory(data[0])
      self.SetTime(data[1])
      self.SetFrame(data[2])
      self.SetFrameType(data[3])
      self.SetDataType(data[4])
      self.SetRIsScaledOut(data[5])
      self.SetMIsScaledOut(data[6])
      self.SetLM(data[7].tolist())
      self.SetData(data[8])
  %}
 };
%extend GWFrames::Waveforms {
  void __setitem__(int i, const GWFrames::Waveform& W) {
    $self->operator[](i) = W;
    return;
  }
 };


//////////////////////////////////////
//// Read in the PNWaveform class ////
//////////////////////////////////////
//// Ignore things that don't translate well...
// %ignore operator<<;
// %ignore GWFrames::Waveform::operator=;
//// These will convert the output data to numpy.ndarray for easier use
%feature("pythonappend") GWFrames::PNWaveform::chi1() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::chi2() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Omega_orb() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Omega_prec() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Omega_tot() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::L() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Omega_orbMag() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Omega_precMag() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Omega_totMag() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::LMag() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::Phi_orb() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::chiHat1() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::chiHat2() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::OmegaHat_orb() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::OmegaHat_prec() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::OmegaHat_tot() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
%feature("pythonappend") GWFrames::PNWaveform::LHat() const %{ if isinstance(val, tuple) : val = numpy.array(val) %}
//// Parse the header file to generate wrappers
%include "PNWaveforms.hpp"


/// Add utility functions that are specific to python.  Note that
/// these are defined in the GWFrames namespace.
%insert("python") %{

### This adds a new method to Waveform objects, returning unwrapped phase.
# def ArgUnwrapper(self, mode) :
#     import numpy
#     return numpy.unwrap(self.Arg(mode))
# Waveform.ArgUnwrapped = ArgUnwrapper


def GetFileNamePrefix(W) :
    return W.DescriptorString() + '_' + W.FrameTypeString() + '_'
Waveform.GetFileNamePrefix = GetFileNamePrefix


def GetLaTeXDataDescription(W) :
    from GWFrames import UnknownDataType, h, hdot, Psi4
    LaTeXDataDescription = ''
    if(W.RIsScaledOut()) :
        LaTeXDataDescription = r'r\,'
    if(W.MIsScaledOut()) :
        if(W.DataType()==UnknownDataType or W.DataType()==h) :
            LaTeXDataDescription = LaTeXDataDescription + W.DataTypeLaTeXString() + r'/M'
        elif(W.DataType()==hdot) :
            LaTeXDataDescription = LaTeXDataDescription + W.DataTypeLaTeXString() # hdot is independent of M
        elif(W.DataType()==Psi4) :
            LaTeXDataDescription = LaTeXDataDescription + r'M\,' + W.DataTypeLaTeXString()
    else :
        LaTeXDataDescription = LaTeXDataDescription + W.DataTypeLaTeXString()
    return LaTeXDataDescription
Waveform.GetLaTeXDataDescription = GetLaTeXDataDescription


def OutputToNRAR(W, FileName, FileWriteMode='w') :
    """
    Output the Waveform in NRAR format.
    
    Note that the FileName is prepended with some descriptive
    information involving the data type and the frame type, such as
    'rhOverM_' or 'rMPsi4_'.
    
    """
    from h5py import File
    from os.path import basename, dirname
    from GWFrames import UnknownDataType, h, hdot, Psi4
    Group = None
    if('.h5' in FileName and not FileName.endswith('.h5')) :
        FileName,Group = FileName.split('.h5')
        FileName += '.h5'
    # Add descriptive prefix to FileName
    if(not dirname(FileName)) :
        FileName = W.DescriptorString() + '_' + basename(FileName)
    else :
        FileName = dirname(FileName) + '/' + W.DescriptorString() + '_' + basename(FileName)
    # Open the file for output
    try :
        F = File(FileName, FileWriteMode)
    except IOError : # If that did not work...
        print("OutputToH5 was unable to open the file '{0}'.\n\n".format(FileName))
        raise # re-raise the exception after the informative message above
    try :
        # If we are writing to a group within the file, create it
        if(Group) :
            G = F.create_group(Group)
        else :
            G = F
        # Now write all the data to various groups in the file
        G.attrs['OutputFormatVersion'] = 'GWFrames_NRAR'
        G.create_dataset("History.txt", data = W.HistoryStr() + '### OutputToNRAR(W, {0})\n'.format(FileName))
	G.attrs['FrameType'] = W.FrameType()
	G.attrs['DataType'] = W.DataType()
	G.attrs['RIsScaledOut'] = int(W.RIsScaledOut())
	G.attrs['MIsScaledOut'] = int(W.MIsScaledOut())
        for i_m in range(W.NModes()) :
            ell,m = W.LM()[i_m]
            Data_m = G.create_dataset("Y_l{0}_m{1}.dat".format(ell, m), data=[[t, d.real, d.imag] for t,d in zip(W.T(),W.Data(i_m))])
            Data_m.attrs['ell'] = ell
            Data_m.attrs['m'] = m
    finally : # Use `finally` to make sure this happens:
        # Close the file and we are done
        F.close()
Waveform.OutputToNRAR = OutputToNRAR


def OutputToH5(W, FileName) :
    """
    Output the Waveform with all necessary information.
    
    Note that the FileName is prepended with some descriptive
    information involving the data type and the frame type, such as
    'rhOverM_Corotating_' or 'rPsi4_Aligned_'.
    
    """
    from h5py import File
    from os.path import basename, dirname
    from GWFrames import UnknownDataType, h, hdot, Psi4
    # Add descriptive prefix to FileName
    FileName = dirname(FileName) + '/' + W.GetFileNamePrefix() + basename(FileName)
    # Open the file for output
    try :
        F = File(FileName, 'w')
    except IOError : # If that did not work...
        print("OutputToH5 was unable to open the file '{0}'.\n\n".format(FileName))
        raise # re-raise the exception after the informative message above
    try :
        # Now write all the data to various groups in the file
        F.attrs['OutputFormatVersion'] = 'GWFrames_v2'
        F.create_dataset("History", data = W.HistoryStr() + '### OutputToH5(W, {0})\n'.format(FileName))
        F.create_dataset("Time", data=W.T().tolist())
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
            Data_m = Data.create_dataset("l{0}_m{1:+}".format(ell, m), data=W.Data(i_m))
            Data_m.attrs['ell'] = ell
            Data_m.attrs['m'] = m
    finally : # Use `finally` to make sure this happens:
        # Close the file and we are done
        F.close()
Waveform.OutputToH5 = OutputToH5


def ReadFromH5(FileName) :
    """
    Read data from an H5 file, as output by GWFrames.
    """
    from h5py import File
    from GWFrames import Waveform, Quaternion
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
        W.AppendHistory("### *this = GWFrames.ReadFromH5(FileName='{0}')\n".format(FileName))
        # Add the old history to the new history
        W.AppendHistory("##### Begin Previous History\n#" + f['History'][()].replace('\n','\n#') + "#### End Previous History\n")
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
    from GWFrames import Waveform
    from Quaternions import Quaternion
    import numpy
    YlmRegex = re.compile(r"""Y_l(?P<L>[0-9]+)_m(?P<M>[-+0-9]+)\.dat""")
    # Initialize the Waveform object
    W = Waveform()
    # Record the filename being read in
    W.AppendHistory("### *this = GWFrames.ReadFromNRAR(FileName='{0}')\n".format(FileName))
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
            OldHistory = f['History.txt'][()]
            W.AppendHistory("##### Begin Previous History\n#" + OldHistory.replace('\n','\n#') + "#### End Previous History\n")
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
            if('psi4' in FileName.lower()) :
                W.SetDataType(GWFrames.WaveformDataType[3])
            elif('hdot' in FileName.lower()) :
                W.SetDataType(GWFrames.WaveformDataType[2])
            elif('h' in FileName.lower()) :
                W.SetDataType(GWFrames.WaveformDataType[1])
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
