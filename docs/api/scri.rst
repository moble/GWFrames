Asymptotic quantities and the BMS group
=======================================

.. currentmodule:: GWFrames

These classes represent fields on the sphere at null infinity and their
transformations under the Bondi–Metzner–Sachs (BMS) group.

Functions on the sphere
-----------------------

.. autoclass:: DataGrid

.. automethod:: DataGrid.Spin
.. automethod:: DataGrid.N_theta
.. automethod:: DataGrid.N_phi
.. automethod:: DataGrid.size
.. automethod:: DataGrid.Data
.. automethod:: DataGrid.SetSpin
.. automethod:: DataGrid.SetNTheta
.. automethod:: DataGrid.SetNPhi
.. automethod:: DataGrid.pow

.. autoclass:: Modes

.. automethod:: Modes.Spin
.. automethod:: Modes.EllMax
.. automethod:: Modes.size
.. automethod:: Modes.Data
.. automethod:: Modes.SetSpin
.. automethod:: Modes.SetEllMax
.. automethod:: Modes.EvaluateAtPoint
.. automethod:: Modes.bar
.. automethod:: Modes.edth
.. automethod:: Modes.edthbar
.. automethod:: Modes.edth2edthbar2
.. automethod:: Modes.pow

.. autoclass:: ScriFunctor
   :members:
   :special-members: __call__

Slices of null infinity
-----------------------

.. autoclass:: SliceOfScriGrid
   :members:

.. autoclass:: SliceOfScriModes
   :members:

.. autoclass:: SliceModes

.. automethod:: SliceModes.EllMax
.. automethod:: SliceModes.Mass
.. automethod:: SliceModes.FourMomentum
.. automethod:: SliceModes.SuperMomentum
.. automethod:: SliceModes.BMSTransformationOnSlice
.. automethod:: SliceModes.MoreschiIteration

.. autoclass:: Scri

.. automethod:: Scri.NTimes
.. automethod:: Scri.T
.. automethod:: Scri.BMSTransformation

.. autoclass:: SuperMomenta

.. automethod:: SuperMomenta.NTimes
.. automethod:: SuperMomenta.T
.. automethod:: SuperMomenta.BMSTransform
.. automethod:: SuperMomenta.MoreschiIteration

Boosts and conformal factors
----------------------------

.. autofunction:: Boost
.. autofunction:: ConformalFactorGrid
.. autofunction:: InverseConformalFactorGrid
.. autofunction:: InverseConformalFactorBoostedGrid
.. autofunction:: vFromOneOverK
