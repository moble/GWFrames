#! /usr/bin/env python
# coding: utf-8

"""Hybridize a single NR waveform with PN.

This script combines the steps needed to read in an NR waveform and
its horizon data (spins, orbital frequencies, etc.), construct the
naively appropriate PN waveform, align the PN by adjusting its time
offset and attitude to match NR as well as possible, and hybridize the
two.

Note that hybridization is a very delicate matter; many things can go
wrong, and result in a bad hybrid.  At the least, you should look at
the output plot, which shows the (2,2) modes of the PN, NR, and hybrid
waveforms -- as a very basic sanity check.

"""

import sys
import os.path
import numpy as np
import h5py
import matplotlib as mpl
mpl.use('Agg')  # Must come after importing mpl, but before importing plt
import matplotlib.pyplot as plt

try : # If this matplotlib is too old, just ignore this
    from matplotlib.pyplot import tight_layout
except :
    pass

import GWFrames
import GWFrames.plot
import Quaternions

def Hybridize(PathToSystem, Waveform='rhOverM_Asymptotic_GeometricUnits.h5/Extrapolated_N2.dir',
              t1=200.0, t2=-sys.float_info.max, InitialOmega_orb=0.0, DirectAlignmentEvaluations=0, Approximant='TaylorT1',
              PlotFileName='Hybridization.pdf', HybridFileName='Hybrid.h5', MinStepsPerOrbit=32,
              PNWaveformModeOrder=3.5, PNOrbitalEvolutionOrder=4.0, Debug=False):

    print("Reading and transforming NR data"); sys.stdout.flush()

    # Read the Waveform from the file, and transform to co-rotating frame
    W_NR_corot = GWFrames.ReadFromNRAR(os.path.join(PathToSystem, Waveform))
    W_NR_corot.DropTimesOutside(0.0, W_NR_corot.T(W_NR_corot.NTimes()-1));
    if Debug:
        W_NR_orig = GWFrames.Waveform(W_NR_corot)
    W_NR_corot.TransformToCorotatingFrame();

    # Place the merger (moment of greatest waveform norm) at t=0, as
    # the only generally meaningful time.
    t0 = -W_NR_corot.MaxNormTime()
    W_NR_corot.SetT(W_NR_corot.T()+t0)
    if Debug:
        W_NR_orig.SetT(W_NR_orig.T()+t0)
    t1 = t1 + t0
    if(t2==-sys.float_info.max):
        t2 = 0.75*t1
    else:
        t2 = t2 + t0

    print("Reading and analyzing Horizons.h5 data"); sys.stdout.flush()

    # Get the BH coordinate data from the h5 file
    with h5py.File(os.path.join(PathToSystem,'Horizons.h5'), 'r') as f:
        tA = f['AhA.dir/CoordCenterInertial.dat'][:,0]+t0
        xA = f['AhA.dir/CoordCenterInertial.dat'][:,1:]
        mA = f['AhA.dir/ArealMass.dat'][:,1]
        chiA = f['AhA.dir/chiInertial.dat'][:,1:]
        tB = f['AhB.dir/CoordCenterInertial.dat'][:,0]+t0
        xB = f['AhB.dir/CoordCenterInertial.dat'][:,1:]
        mB = f['AhB.dir/ArealMass.dat'][:,1]
        chiB = f['AhB.dir/chiInertial.dat'][:,1:]

    # Some calculations with the BH coordinates
    d = xA-xB
    nHat = Quaternions.normalized(Quaternions.QuaternionArray(d))
    dnHatdt = Quaternions.QuaternionDerivative(nHat, tA)
    lambdaHat = Quaternions.normalized(dnHatdt)
    Omega = np.array([np.cross(n.vec(),ndot.vec()) for n,ndot in zip(nHat,dnHatdt)])

    # Take initial data to be simply the coordinate quantities at the relaxation time
    i_1 = abs(tA-t1).argmin()
    t_i = tA[i_1]
    ma = mA[i_1]
    mb = mB[i_1]
    delta = (ma-mb)/(ma+mb)
    chia_0 = chiA[i_1]
    chib_0 = chiB[i_1]
    Omega_orb_0 = np.sqrt(Omega[i_1,0]**2+Omega[i_1,1]**2+Omega[i_1,2]**2)
    R_frame_i = Quaternions.FrameFromXY([nHat[i_1],],[lambdaHat[i_1],])[0]
    if(InitialOmega_orb==0.0):
        InitialOmega_orb = 0.5*Omega_orb_0

    # Construct complete PN waveform
    print("Constructing PN data:");
    print("GWFrames.PNWaveform(Approximant={0}, delta={1},\n"
          "    chia_0={2}, chib_0={3},\n"
          "    Omega_orb_0={4}, InitialOmega_orb={5},\n"
          "    R_frame_i={6},\n"
          "    MinStepsPerOrbit={7}, PNWaveformModeOrder={8}, PNOrbitalEvolutionOrder={9}"
          ")".format(Approximant, delta, chia_0, chib_0, Omega_orb_0, InitialOmega_orb, R_frame_i,
                     MinStepsPerOrbit, PNWaveformModeOrder, PNOrbitalEvolutionOrder))
    sys.stdout.flush()
    W_PN_corot = GWFrames.PNWaveform(Approximant, delta, chia_0, chib_0, Omega_orb_0, InitialOmega_orb, R_frame_i,
                                     MinStepsPerOrbit, PNWaveformModeOrder, PNOrbitalEvolutionOrder)
    W_PN_corot.SetT(W_PN_corot.T()+t_i);
    if Debug:
        W_PN_orig = GWFrames.PNWaveform(W_PN_corot)
    W_PN_corot.TransformToCorotatingFrame();

    print("Aligning PN and NR waveforms"); sys.stdout.flush()

    # Minimize the distance between the PN and NR rotors in their
    # co-rotating frames.  Note that NR is kept fixed here, because
    # its merger time is set to t=0.0, which is the only physically
    # meaningful fiducial quantity present.  That quality will
    # therefore also be present in the output waveform.
    GWFrames.AlignWaveforms(W_NR_corot, W_PN_corot, t1, t2, DirectAlignmentEvaluations, [], Debug)

    print("Constructing PN-NR hybrid"); sys.stdout.flush()

    # Construct the hybrid
    W_hyb_corot = GWFrames.Waveform(W_PN_corot.Hybridize(W_NR_corot, t1, t2))

    # Optionally, make a plot showing the results.  It is strongly
    # recommended that you look at this plot to make sure things worked
    # alright.
    if(PlotFileName):
        print("Plotting to {0}".format(os.path.join(PathToSystem, PlotFileName))); sys.stdout.flush()
        plt.axvline(t1, linestyle='dotted')
        plt.axvline(t2, linestyle='dotted')
        plt.title(r'Modes of the component waveforms')
        W_hyb_corot.plot('LogAbs', [2,2], lw=3, label='Hybrid (2,2)')
        W_PN_corot.plot('LogAbs', [2,2], label='PN (2,2)')
        W_NR_corot.plot('LogAbs', [2,2], label='NR (2,2)')
        W_hyb_corot.plot('LogAbs', [2,1], lw=3, label='Hybrid (2,1)')
        W_PN_corot.plot('LogAbs', [2,1], label='PN (2,1)')
        W_NR_corot.plot('LogAbs', [2,1], label='NR (2,1)')
        if Debug:
            W_NR_orig.plot('LogAbs', [2,2], label='NR$_0$ (2,2)', ls='--')
            W_NR_orig.plot('LogAbs', [2,1], label='NR$_0$ (2,1)', ls='--')
            W_PN_orig.plot('LogAbs', [2,2], label='PN$_0$ (2,2)', ls='--')
            W_PN_orig.plot('LogAbs', [2,1], label='PN$_0$ (2,1)', ls='--')
        plt.ylim((1e-8,1.))
        plt.xlim((W_NR_corot.T(0), W_NR_corot.T()[-1]))
        plt.legend(loc='center')
        try:
            tight_layout(pad=0.1)
        except:
            pass
        plt.savefig(os.path.join(PathToSystem, PlotFileName))

    print("Transforming hybrid to inertial frame and outputting to\n\t"
          "{0}".format(os.path.join(PathToSystem, HybridFileName))); sys.stdout.flush()
    W_hyb_corot.TransformToInertialFrame()
    W_hyb_corot.OutputToNRAR(os.path.join(PathToSystem, HybridFileName))

    print("All done"); sys.stdout.flush()



if __name__ == "__main__" :
    import os
    import os.path
    import sys
    import argparse

    # Set up and run the parser
    parser = argparse.ArgumentParser(description = __doc__)
    parser.add_argument('PathToSystem', default=os.curdir,
                        help="Path (relative or absolute) to the data directory containing the system's waveform and horizon data")
    parser.add_argument('--Waveform', default='rhOverM_Asymptotic_GeometricUnits.h5/Extrapolated_N2.dir',
                        help='Path (relative to `PathToSystem`) in which to find the waveform data [default: rhOverM_Asymptotic_GeometricUnits.h5/Extrapolated_N2.dir]')
    parser.add_argument('--t1', type=float, default=200.0,
                        help='Beginning of alignment/hybridization interval [default: 200.0 from the beginning of NR]')
    parser.add_argument('--t2', type=float, default=-sys.float_info.max,
                        help='End of alignment/hybridization interval [default: 20%% closer to merger than t1]')
    parser.add_argument('--InitialOmega_orb', type=float, default=0.0,
                        help='Earliest frequency in the PN waveform [default: half the frequency at t1]')
    parser.add_argument('--DirectAlignmentEvaluations', type=int, default=0,
                        help='Number of evaluations to use in the simple first stage of alignment [default: 0]')
    parser.add_argument('--Approximant', default='TaylorT1',
                        help='PN approximant; TaylorT1|TaylorT4|TaylorT5 [default: TaylorT1]')
    parser.add_argument('--PNWaveformModeOrder', type=float, default=3.5,
                        help='PN order at which to compute waveform modes [default: 3.5]')
    parser.add_argument('--PNOrbitalEvolutionOrder', type=float, default=4.0,
                        help='PN order at which to compute orbital evolution [default: 4.0]')
    parser.add_argument('--MinStepsPerOrbit', type=int, default=32,
                        help='Minimum number of time steps at which to evaluate in each orbit [default: 32]')
    parser.add_argument('--PlotFileName', default='Hybridization.pdf',
                        help='File name to which the result is plotted for diagnostics; empty string for no output [default: Hybridization.pdf]')
    parser.add_argument('--HybridFileName', default='Hybrid.h5',
                        help='File name to which the result is plotted for diagnostics [default: Hybrid.h5]')
    parser.add_argument('--Debug', action='store_true',
                        help='Output debugging information during alignment')
    args = vars(parser.parse_args())

    if(not args):
        parser.print_help()
        sys.exit(1)

    PathToSystem = os.path.abspath(args['PathToSystem'])

    Hybridize(PathToSystem, args['Waveform'], args['t1'], args['t2'], args['InitialOmega_orb'], args['DirectAlignmentEvaluations'],
              args['Approximant'], args['PlotFileName'], args['HybridFileName'], args['MinStepsPerOrbit'],
              args['PNWaveformModeOrder'], args['PNOrbitalEvolutionOrder'], args['Debug'])


# PathToSXSCatalog = '/Users/boyle/Research/Data/SimulationAnnex/Catalog/'
# PathToSystem = PathToSXSCatalog + 'q1.0/SBBH/d19.0_q1.0_s0.5_0_0_s0_0_0/Lev5/'
