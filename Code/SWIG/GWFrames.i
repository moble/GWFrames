// -*- c++ -*-
// vim: filetype=cpp

// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

%module GWFrames

// Quiet warnings about overloaded operators being ignored.
//#pragma SWIG nowarn=362,389,401,509


/////////////////////////////////////////////////
//// These will be needed by the c++ wrapper ////
/////////////////////////////////////////////////
%{
  #include <iostream>
  #include <string>
  #include <sstream>
  #include <iomanip>
  #include <complex>
  #include "../Utilities.hpp"
  #include "Quaternions.hpp"
  #include "IntegrateAngularVelocity.hpp"
  #include "../Scri.hpp"
  #include "SWSHs.hpp"
  #include "../Waveforms.hpp"
  #include "../PNWaveforms.hpp"
  #include "../WaveformsAtAPointFT.hpp"
%}


////////////////////////////////////////////////////////
//// This lets me use numpy.array in the code below ////
////////////////////////////////////////////////////////
%{
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%init %{
import_array();
%}


//////////////////////////////////////////////////////////////////////
//// The following translates between c++ and python types nicely ////
//////////////////////////////////////////////////////////////////////

%import "../Quaternions/Quaternions.i"
%import "../Quaternions/Quaternions_typemaps.i"
//%import "../SphericalFunctions/SphericalFunctions.i"
%include "../Quaternions/vector_typemaps.i"

%include <typemaps.i>
%include <stl.i>
//// Make sure std::strings are dealt with appropriately
%include <std_string.i>
//// Make sure std::complex numbers are dealt with appropriately
%include <std_complex.i>


/////////////////////////////////////////////////////////////////////
//// The doxygen-generated documentation is added with this line ////
/////////////////////////////////////////////////////////////////////
%include "../Docs/GWFrames_Doc.i"


///////////////////////////////////
//// Handle error codes nicely ////
///////////////////////////////////
%include "Exceptions.i"


///////////////////////////////////////////////////////////////////////////////
//// This ensures that GWFrames has local copies of these modules imported ////
///////////////////////////////////////////////////////////////////////////////
%pythoncode %{
  ## We have to be able to import numpy
  import numpy;

  ## We have to be able to import Quaternions
  import Quaternions

  ## We do not necessarily need spinsfast
  try :
    import spinsfast
  except ImportError :
    pass
  %}


/////////////////////////////////////////////////////////////////////////////////
//// Finally, include the descriptions of the actual C++ code in this module ////
/////////////////////////////////////////////////////////////////////////////////
%include "Utilities.i"
%include "Scri.i"
%include "Waveforms.i"
%include "PNWaveforms.i"
%include "Extensions.py"
