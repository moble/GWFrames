// Copyright (c) 2014, Michael Boyle
// See LICENSE file for details

#ifndef ERRORS_HPP
#define ERRORS_HPP

// Note: These error codes are used in GWFrames.i

#define GWFrames_NotYetImplemented 0
#define GWFrames_IndexOutOfBounds 1
#define GWFrames_InfinitelyManySolutions 2
#define GWFrames_VectorSizeMismatch 3
#define GWFrames_CannotExtrapolateQuaternions 4
#define GWFrames_MatrixSizeMismatch 5
#define GWFrames_MatrixSizeAssumedToBeThree 6
#define GWFrames_QuaternionVectorSizeNotUnderstood 7
#define GWFrames_WaveformMissingLMIndex 8
#define GWFrames_BadFileName 9
#define GWFrames_NotEnoughPointsForDerivative 10
#define GWFrames_EmptyIntersection 11
#define GWFrames_FailedSystemCall 12
#define GWFrames_WrongFrameType 13
#define GWFrames_FailedGSLCall 14
#define GWFrames_BadWaveformInformation 15
#define GWFrames_BadSwitches 16
#define GWFrames_ValueError 17

#endif // ERRORS_HPP
