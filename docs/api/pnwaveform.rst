PNWaveform
==========

.. currentmodule:: GWFrames

.. note::

   For post-Newtonian waveforms, use `PostNewtonian.jl
   <https://github.com/moble/PostNewtonian.jl>`_, which is far more complete,
   accurate, and flexible than this code.

.. autoclass:: PNWaveform

`PNWaveform` is a subclass of `Waveform`, so all of the methods of `Waveform`
apply to it as well.  The methods below return the quantities of the
post-Newtonian orbital evolution.

Spins
-----

.. automethod:: PNWaveform.chi1
.. automethod:: PNWaveform.chi2
.. automethod:: PNWaveform.chi1Mag
.. automethod:: PNWaveform.chi2Mag
.. automethod:: PNWaveform.chiHat1
.. automethod:: PNWaveform.chiHat2

Orbital and precessional motion
-------------------------------

.. automethod:: PNWaveform.Omega_orb
.. automethod:: PNWaveform.Omega_orbMag
.. automethod:: PNWaveform.OmegaHat_orb
.. automethod:: PNWaveform.Omega_prec
.. automethod:: PNWaveform.Omega_precMag
.. automethod:: PNWaveform.OmegaHat_prec
.. automethod:: PNWaveform.Omega_tot
.. automethod:: PNWaveform.Omega_totMag
.. automethod:: PNWaveform.OmegaHat_tot
.. automethod:: PNWaveform.Phi_orb

Angular momentum
----------------

.. automethod:: PNWaveform.L
.. automethod:: PNWaveform.LMag
.. automethod:: PNWaveform.LHat
