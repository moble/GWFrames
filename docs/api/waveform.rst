Waveform
========

.. currentmodule:: GWFrames

.. autoclass:: Waveform

Frame types and data types
--------------------------

The frame type and the data type of a `Waveform` are integers, given names by
these module-level constants.

.. py:data:: UnknownFrameType
             Inertial
             Coprecessing
             Coorbital
             Corotating

   The frames in which the modes of a `Waveform` may be decomposed, as
   returned by `Waveform.FrameType` and set by `Waveform.SetFrameType`.

.. py:data:: UnknownDataType
             h
             hdot
             Psi4
             Psi3
             Psi2
             Psi1
             Psi0

   The quantities that a `Waveform` may hold, as returned by
   `Waveform.DataType` and set by `Waveform.SetDataType`.

Copying, slicing, and interpolating
-----------------------------------

.. automethod:: Waveform.CopyWithoutData
.. automethod:: Waveform.SliceOfTimeIndices
.. automethod:: Waveform.SliceOfTimeIndicesWithEll2
.. automethod:: Waveform.SliceOfTimeIndicesWithoutModes
.. automethod:: Waveform.SliceOfTimes
.. automethod:: Waveform.SliceOfTimesWithEll2
.. automethod:: Waveform.SliceOfTimesWithoutModes
.. automethod:: Waveform.Interpolate
.. automethod:: Waveform.InterpolateInPlace

Accessing the data
------------------

.. automethod:: Waveform.NTimes
.. automethod:: Waveform.NModes
.. automethod:: Waveform.T
.. automethod:: Waveform.Frame
.. automethod:: Waveform.LM
.. automethod:: Waveform.Data
.. automethod:: Waveform.Re
.. automethod:: Waveform.Im
.. automethod:: Waveform.Abs
.. automethod:: Waveform.Arg
.. automethod:: Waveform.ArgUnwrapped
.. automethod:: Waveform.SpinWeight
.. automethod:: Waveform.BoostWeight
.. automethod:: Waveform.FrameType
.. automethod:: Waveform.FrameTypeString
.. automethod:: Waveform.DataType
.. automethod:: Waveform.DataTypeString
.. automethod:: Waveform.DataTypeLaTeXString
.. automethod:: Waveform.DescriptorString
.. automethod:: Waveform.GetFileNamePrefix
.. automethod:: Waveform.GetLaTeXDataDescription
.. automethod:: Waveform.RIsScaledOut
.. automethod:: Waveform.MIsScaledOut
.. automethod:: Waveform.HistoryStr
.. automethod:: Waveform.VersionHist

Changing the data
-----------------

.. automethod:: Waveform.DropTimesOutside
.. automethod:: Waveform.DropEllModes
.. automethod:: Waveform.KeepOnlyEllModes
.. automethod:: Waveform.KeepOnlyEll2
.. automethod:: Waveform.SetT
.. automethod:: Waveform.SetTime
.. automethod:: Waveform.SetFrame
.. automethod:: Waveform.SetFrameType
.. automethod:: Waveform.SetDataType
.. automethod:: Waveform.SetRIsScaledOut
.. automethod:: Waveform.SetMIsScaledOut
.. automethod:: Waveform.SetLM
.. automethod:: Waveform.SetData
.. automethod:: Waveform.ResizeData
.. automethod:: Waveform.SetSpinWeight
.. automethod:: Waveform.SetBoostWeight
.. automethod:: Waveform.AppendHistory
.. automethod:: Waveform.SetHistory
.. automethod:: Waveform.SetVersionHist
.. automethod:: Waveform.swap

Characterizing the data
-----------------------

.. automethod:: Waveform.EllMax
.. automethod:: Waveform.FindModeIndex
.. automethod:: Waveform.FindModeIndexWithoutError
.. automethod:: Waveform.DataDot
.. automethod:: Waveform.Differentiate
.. automethod:: Waveform.Norm
.. automethod:: Waveform.MaxNormIndex
.. automethod:: Waveform.MaxNormTime

Rotations and frames
--------------------

.. automethod:: Waveform.RotatePhysicalSystem
.. automethod:: Waveform.RotateDecompositionBasis
.. automethod:: Waveform.LdtVector
.. automethod:: Waveform.LLMatrix
.. automethod:: Waveform.LLDominantEigenvector
.. automethod:: Waveform.AngularVelocityVector
.. automethod:: Waveform.AngularVelocityVectorRelativeToInertial
.. automethod:: Waveform.CorotatingFrame
.. automethod:: Waveform.TransformToCoprecessingFrame
.. automethod:: Waveform.TransformToAngularVelocityFrame
.. automethod:: Waveform.TransformToCorotatingFrame
.. automethod:: Waveform.TransformToInertialFrame

Alignment, comparison, and hybridization
----------------------------------------

.. automethod:: Waveform.GetAlignmentsOfDecompositionFrameToModes
.. automethod:: Waveform.GetAlignmentOfDecompositionFrameToModes
.. automethod:: Waveform.AlignDecompositionFrameToModes
.. automethod:: Waveform.Compare
.. automethod:: Waveform.Hybridize

.. autofunction:: AlignWaveforms

Symmetries and parity
---------------------

.. automethod:: Waveform.NormalizedAntisymmetry
.. automethod:: Waveform.DipoleMoment
.. automethod:: Waveform.MinimalParityViolation
.. automethod:: Waveform.XParityInvolution
.. automethod:: Waveform.XParitySymmetricPart
.. automethod:: Waveform.XParityAntisymmetricPart
.. automethod:: Waveform.XParityViolationSquared
.. automethod:: Waveform.XParityViolationNormalized
.. automethod:: Waveform.YParityInvolution
.. automethod:: Waveform.YParitySymmetricPart
.. automethod:: Waveform.YParityAntisymmetricPart
.. automethod:: Waveform.YParityViolationSquared
.. automethod:: Waveform.YParityViolationNormalized
.. automethod:: Waveform.ZParityInvolution
.. automethod:: Waveform.ZParitySymmetricPart
.. automethod:: Waveform.ZParityAntisymmetricPart
.. automethod:: Waveform.ZParityViolationSquared
.. automethod:: Waveform.ZParityViolationNormalized
.. automethod:: Waveform.ParityInvolution
.. automethod:: Waveform.ParitySymmetricPart
.. automethod:: Waveform.ParityAntisymmetricPart
.. automethod:: Waveform.ParityViolationSquared
.. automethod:: Waveform.ParityViolationNormalized
.. automethod:: Waveform.ConjugateAntipodalEvaluation
.. automethod:: Waveform.ConjugateAntipodalSymmetricPart
.. automethod:: Waveform.ConjugateAntipodalAntisymmetricPart

Evaluation on the sphere, translations, and boosts
--------------------------------------------------

.. automethod:: Waveform.EvaluateAtPoint
.. automethod:: Waveform.InterpolateToPoint
.. automethod:: Waveform.Translate
.. automethod:: Waveform.BoostPsi4
.. automethod:: Waveform.BoostHFaked

Arithmetic
----------

Two waveforms with the same times, modes, frames, and types can be added,
subtracted, multiplied, and divided mode by mode with the usual operators, and
a waveform can be multiplied or divided by a number.

Input and output
----------------

.. automethod:: Waveform.Output
.. automethod:: Waveform.OutputToH5
.. automethod:: Waveform.OutputToNRAR

.. autofunction:: ReadFromH5
.. autofunction:: ReadFromNRAR
.. autofunction:: MonotonicIndices
