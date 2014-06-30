///////////////////////////////////
//// Handle exceptions cleanly ////
///////////////////////////////////

// The following will appear in the header of the `_wrap.cpp` file.
%{
  // The following allows us to elegantly fail in python from manual
  // interrupts and floating-point exceptions found in the c++ code.
  // The setjmp part of this was inspired by the post
  // <http://stackoverflow.com/a/12155582/1194883>.  The code for
  // setting the csr flags is taken from SpEC.
  #include <csetjmp>
  #include <csignal>
  #ifdef __APPLE__
    #include <xmmintrin.h>
    int fegetexcept() { return _mm_getcsr(); }
  #else
    #include <fenv.h>     // For feenableexcept. Doesn't seem to be a <cfenv>
    #ifndef __USE_GNU
      extern "C" int feenableexcept (int EXCEPTS);
    #endif
  #endif
  namespace GWFrames {
    static sigjmp_buf FloatingPointExceptionJumpBuffer;
    static sigjmp_buf InterruptExceptionJumpBuffer;
    void FloatingPointExceptionHandler(int sig) {
      siglongjmp(FloatingPointExceptionJumpBuffer, sig);
    }
    void InterruptExceptionHandler(int sig) {
      siglongjmp(InterruptExceptionJumpBuffer, sig);
    }
    class ExceptionHandlerSwitcher {
    private:
      int OriginalExceptionFlags;
      void (*OriginalFloatingPointExceptionHandler)(int);
      void (*OriginalInterruptExceptionHandler)(int);
    public:
      ExceptionHandlerSwitcher()
        : OriginalExceptionFlags(fegetexcept()),
          OriginalFloatingPointExceptionHandler(signal(SIGFPE, FloatingPointExceptionHandler)),
          OriginalInterruptExceptionHandler(signal(SIGINT, InterruptExceptionHandler))
      {
std::cout << "ExceptionHandlerSwitcher()" << std::endl;
        #ifdef __APPLE__
          _mm_setcsr( _MM_MASK_MASK &~
                     (_MM_MASK_OVERFLOW|_MM_MASK_INVALID|_MM_MASK_DIV_ZERO));
        #else
          feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW);
        #endif
      }
      ~ExceptionHandlerSwitcher() {
std::cout << "~ExceptionHandlerSwitcher()" << std::endl;
        #ifdef __APPLE__
          _mm_setcsr(OriginalExceptionFlags);
        #else
          feenableexcept(OriginalExceptionFlags);
        #endif
        signal(SIGFPE, OriginalFloatingPointExceptionHandler);
        signal(SIGINT, OriginalInterruptExceptionHandler);
      }
    };
  }

  // The following allows us to elegantly fail in python from
  // exceptions raised by the c++ code.
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

  if (!sigsetjmp(GWFrames::FloatingPointExceptionJumpBuffer, 1)) {
    if(!sigsetjmp(GWFrames::InterruptExceptionJumpBuffer, 1)) {
      try {
        const GWFrames::ExceptionHandlerSwitcher Switcher;
        $action;
      } catch(int i) {
        std::stringstream s;
        if(i>-1 && i<GWFramesNumberOfErrors) { s << "$fulldecl: " << GWFramesErrors[i]; }
        else  { s << "$fulldecl: Unknown exception number {" << i << "}"; }
        PyErr_SetString(GWFramesExceptions[i], s.str().c_str());
        return 0;
      } catch(...) {
        PyErr_SetString(PyExc_RuntimeError, "$fulldecl: Unknown exception; default handler");
        return 0;
      }
    } else {
      PyErr_SetString(PyExc_RuntimeError, "$fulldecl: Caught a manual interrupt in the c++ code.");
      return 0;
    }
  } else {
    PyErr_SetString(PyExc_RuntimeError, "$fulldecl: Caught a floating-point exception in the c++ code.");
    return 0;
  }

}
