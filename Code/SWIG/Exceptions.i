///////////////////////////////////
//// Handle exceptions cleanly ////
///////////////////////////////////

// The following will appear in the header of the `_wrap.cpp` file.
%{
  const char* const GWFramesErrors[] = {
    "This function is not yet implemented.",
    "Failed system call.",
    "Bad file name.",
    "Failed GSL call.",
    "Unknown exception",
    "Unknown exception",
    "Unknown exception",
    "Unknown exception",
    "Unknown exception",
    "Unknown exception",
    "Bad value.",
    "Bad switches; we should not have gotten here.",
    "Index out of bounds.",
    "Unknown exception",
    "Unknown exception",
    "Vector size mismatch.",
    "Matrix size mismatch.",
    "Matrix size is assumed to be 3x3 in this function.",
    "Not enough points to take a derivative.",
    "Empty intersection requested.",
    "Waveform is missing requested (ell,m) component.",
    "Wrong frame type for this operation.",
    "Bad Waveform information.",
    "Unknown GW detector name."
  };
  const int GWFramesNumberOfErrors = 24;
  PyObject* const GWFramesExceptions[] = {
    PyExc_NotImplementedError, // Not implemented
    PyExc_SystemError, // Failed system call
    PyExc_IOError, // Bad file name
    PyExc_RuntimeError, // GSL failed
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_ValueError, // Bad value
    PyExc_ValueError, // Bad switches
    PyExc_IndexError, // Index out of bounds
    PyExc_RuntimeError, // [empty]
    PyExc_RuntimeError, // [empty]
    PyExc_AssertionError, // Mismatched vector size
    PyExc_AssertionError, // Mismatched matrix size
    PyExc_AssertionError, // 3x3 matrix assumed
    PyExc_AssertionError, // Not enough points for derivative
    PyExc_AssertionError, // Empty intersection
    PyExc_IndexError, // Waveform missing ell,m
    PyExc_AssertionError, // Bad frame type
    PyExc_ValueError, // Bad Waveform information
    PyExc_ValueError, // Unknown GW detector name
  };
%}

// This will go inside every python wrapper for any function I've
// included; the code of the function itself will replace `$action`.
// It's a good idea to try to keep this part brief, just to cut down
// the size of the wrapper file.
%exception {
  try {
    $action;
  } catch(int i) {
    std::stringstream s;
    if(i>-1 && i<GWFramesNumberOfErrors) { s << "GWFrames exception: " << GWFramesErrors[i]; }
    else  { s << "GWFrames: Unknown exception number {" << i << "}"; }
    PyErr_SetString(GWFramesExceptions[i], s.str().c_str());
    return 0; // NULL;
  }
}
