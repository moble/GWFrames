// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#include "bindings.hpp"

#include <exception>
#include <sstream>

namespace {

// The C++ code reports errors by throwing the integer codes defined in the
// `Errors.hpp` files of GWFrames and its submodules (and, through custom GSL
// error handlers, GSL's own error numbers).  The codes 20, 21, and 22 mean
// different things in GWFrames and in Quaternions, so both meanings are given.
struct ErrorDescription {
  PyObject* type;
  const char* message;
};

ErrorDescription DescribeError(const int code) {
  switch (code) {
    case 0: return {PyExc_NotImplementedError, "This function is not yet implemented."};
    case 1: return {PyExc_SystemError, "Failed system call."};
    case 2: return {PyExc_IOError, "Bad file name."};
    case 3: return {PyExc_RuntimeError, "Failed GSL call."};
    case 10: return {PyExc_ValueError, "Bad value."};
    case 11: return {PyExc_ValueError, "Bad switches; we should not have gotten here."};
    case 12: return {PyExc_IndexError, "Index out of bounds."};
    case 15: return {PyExc_ValueError, "Vector size mismatch."};
    case 16: return {PyExc_ValueError, "Matrix size mismatch."};
    case 17: return {PyExc_ValueError, "Matrix size is assumed to be 3x3 in this function."};
    case 18: return {PyExc_ValueError, "Not enough points to take a derivative."};
    case 19: return {PyExc_ValueError, "Empty intersection requested."};
    case 20:
      return {PyExc_IndexError,
              "Waveform is missing requested (ell,m) component "
              "(or, from Quaternions, infinitely many solutions in quaternion algebra)."};
    case 21:
      return {PyExc_ValueError,
              "Wrong frame type for this operation "
              "(or, from Quaternions, vector size should be 3 or 4 in Quaternion constructor)."};
    case 22:
      return {PyExc_ValueError,
              "Bad Waveform information (or, from Quaternions, cannot extrapolate Quaternion arrays)."};
    case 23: return {PyExc_ValueError, "Unknown GW detector name."};
    default: return {PyExc_RuntimeError, "Unknown exception."};
  }
}

}  // namespace

PYBIND11_MODULE(_GWFrames, m) {
  namespace py = pybind11;
  m.doc() = "Compiled core of GWFrames";

  py::register_exception_translator([](std::exception_ptr p) {
    try {
      if (p) {
        std::rethrow_exception(p);
      }
    } catch (const int code) {
      const auto description = DescribeError(code);
      std::ostringstream message;
      message << "GWFrames error code " << code << ": " << description.message
              << "  Details may have been printed to stderr.";
      PyErr_SetString(description.type, message.str().c_str());
    }
  });

  auto quaternions = m.def_submodule("Quaternions", "Quaternions and functions of arrays of quaternions");
  GWFramesBindings::BindQuaternions(quaternions);
  GWFramesBindings::BindUtilities(m);
  GWFramesBindings::BindWaveforms(m);
  GWFramesBindings::BindScri(m);
}
