// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#ifndef ERRORS_HPP
#define ERRORS_HPP

// Note: These error codes are used in SWIG/Exceptions.i.
//       If you change them here, change them there.

#define GWFrames_NotYetImplemented 0
#define GWFrames_FailedSystemCall 1
#define GWFrames_BadFileName 2
#define GWFrames_FailedGSLCall 3
// #define GWFrames_ 4
// #define GWFrames_ 5
// #define GWFrames_ 6
// #define GWFrames_ 7
// #define GWFrames_ 8
// #define GWFrames_ 9
#define GWFrames_ValueError 10
#define GWFrames_BadSwitches 11
#define GWFrames_IndexOutOfBounds 12
// #define GWFrames_ 13
// #define GWFrames_ 14
#define GWFrames_VectorSizeMismatch 15
#define GWFrames_MatrixSizeMismatch 16
#define GWFrames_MatrixSizeAssumedToBeThree 17
#define GWFrames_NotEnoughPointsForDerivative 18
#define GWFrames_EmptyIntersection 19
#define GWFrames_WaveformMissingLMIndex 20
#define GWFrames_WrongFrameType 21
#define GWFrames_BadWaveformInformation 22
#define GWFrames_UnknownDetector 23

#endif // ERRORS_HPP
