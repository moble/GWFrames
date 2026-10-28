Quaternions
===========

.. module:: GWFrames.Quaternions

.. currentmodule:: GWFrames.Quaternions

The submodule `GWFrames.Quaternions` provides the `Quaternion` class and
functions acting on single quaternions and on arrays of them.  Unit
quaternions (rotors) represent the rotations between frames.  The class is
also available as ``GWFrames.Quaternion``.

Wherever a `Quaternion` is expected, a sequence of four numbers is also
accepted, and wherever an array of quaternions is expected, a sequence of
`Quaternion` objects, a sequence of length-4 sequences, or an array of shape
``(N, 4)`` is accepted.  Arrays of quaternions are returned as
one-dimensional numpy arrays of `Quaternion` objects.

The Quaternion class
--------------------

.. autoclass:: Quaternion

.. automethod:: Quaternion.abs
.. automethod:: Quaternion.normsquared
.. automethod:: Quaternion.angle
.. automethod:: Quaternion.conjugate
.. automethod:: Quaternion.inverse
.. automethod:: Quaternion.normalized
.. automethod:: Quaternion.vec
.. automethod:: Quaternion.log
.. automethod:: Quaternion.logRotor
.. automethod:: Quaternion.exp
.. automethod:: Quaternion.sqrt
.. automethod:: Quaternion.sqrtOfRotor
.. automethod:: Quaternion.pow
.. automethod:: Quaternion.dot
.. automethod:: Quaternion.cross
.. automethod:: Quaternion.commutator
.. automethod:: Quaternion.IntrinsicDistance
.. automethod:: Quaternion.ChordalDistance
.. automethod:: Quaternion.XParityConjugateSpinor
.. automethod:: Quaternion.YParityConjugateSpinor
.. automethod:: Quaternion.ZParityConjugateSpinor
.. automethod:: Quaternion.ParityConjugateSpinor
.. automethod:: Quaternion.str

.. py:data:: Zero
             One
             xHat
             yHat
             zHat

   The quaternions 0 and 1, and the unit basis vectors :math:`\hat{x}`,
   :math:`\hat{y}`, and :math:`\hat{z}` as pure-vector quaternions.

Functions of quaternions and arrays
-----------------------------------

.. autofunction:: abs
.. autofunction:: normsquared
.. autofunction:: angle
.. autofunction:: conjugate
.. autofunction:: inverse
.. autofunction:: normalized
.. autofunction:: vec
.. autofunction:: log
.. autofunction:: logRotor
.. autofunction:: exp
.. autofunction:: sqrt
.. autofunction:: sqrtOfRotor
.. autofunction:: pow
.. autofunction:: dot
.. autofunction:: cross
.. autofunction:: commutator
.. autofunction:: IntrinsicDistance
.. autofunction:: ChordalDistance
.. autofunction:: Slerp
.. autofunction:: XParityConjugateSpinor
.. autofunction:: YParityConjugateSpinor
.. autofunction:: ZParityConjugateSpinor
.. autofunction:: ParityConjugateSpinor

Constructing and taking apart arrays
------------------------------------

.. autofunction:: QuaternionArray
.. autofunction:: Component
.. autofunction:: Component0
.. autofunction:: Component1
.. autofunction:: Component2
.. autofunction:: Component3

Calculus and frames
-------------------

.. autofunction:: FrameFromAngularVelocity
.. autofunction:: FrameFromAngularVelocity_2D
.. autofunction:: FrameFromAngularVelocity_Integrand
.. autofunction:: FrameFromAngularVelocity_2D_Integrand
.. autofunction:: FrameAngularVelocity
.. autofunction:: MinimalRotation
.. autofunction:: PrescribedRotation
.. autofunction:: FrameFromXY
.. autofunction:: FrameFromZ
.. autofunction:: FrameFromPrescribedRotation
.. autofunction:: DifferentiateRotorByLogarithm
.. autofunction:: QuaternionDerivative
.. autofunction:: IndefiniteIntegral
.. autofunction:: DefiniteIntegral
.. autofunction:: UnflipRotors
.. autofunction:: RDelta
.. autofunction:: Squad
.. autofunction:: BoostRotor

Mean rotors and alignment
-------------------------

.. autofunction:: MeanRotor
.. autofunction:: ApproximateMeanRotor
.. autofunction:: ApproximateOptimalAlignmentRotor
.. autofunction:: ApproximateOptimalAlignment
.. autofunction:: OptimalAlignment

Searching sorted arrays
-----------------------

.. autofunction:: hunt
.. autofunction:: huntRight
