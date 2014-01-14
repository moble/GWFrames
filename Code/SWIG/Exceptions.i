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
    } else if(i==15) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Bad Waveform information.");
    } else if(i==16) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Bad switches; we should not have gotten here.");
    } else if(i==17) {
      PyErr_SetString(PyExc_RuntimeError, "GWFrames: Bad value.");
    } else  {
      PyErr_SetString(PyExc_RuntimeError, "Unknown exception");
    }
    return NULL;
  }
}
