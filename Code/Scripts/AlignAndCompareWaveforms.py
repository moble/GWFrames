# List the file (and directory within each H5 file) that you want to
# compare.  The first one will be used as the fiducial waveform, and
# everything else will be shifted to best agree with it.
Files = ['/Users/boyle/Research/Data/SimulationAnnex/Catalog/SKS/d15_q1_sA_0_0_0.97_sB_0_0_0.97_ecc6e-4/Lev6/rhOverM_Asymptotic_GeometricUnits.h5/Extrapolated_N2.dir',
         '/Users/boyle/Research/Data/SimulationAnnex/Catalog/SKS/d15_q1_sA_0_0_0.97_sB_0_0_0.97_ecc6e-4/Lev6/rhOverM_Asymptotic_GeometricUnits.h5/Extrapolated_N4.dir']

# This is an automatic legend generator for each of the lines in the
# plots.  Alternatively, just write your own list of legends.  They
# should be strings, and there should be one less than the number of
# Files, because everything is being compared to the first File.
Legends = [File.split('/')[-1] for File in Files[1:]]

# These are the names of the files into which the plots will be saved
PhaseDifferencePlot = 'PhaseDifference.pdf'
AmplitudeDifferencePlot = 'AmplitudeDifference.pdf'

# The following two numbers give the range over which phase differences are minimized
t1 = 1000.
t2 = 2000.
tmid = (t1+t2)/2.

# The following two numbers just give the plot range
tmin = 0.
tmax = 6550.



## Things below should not need to change often
from numpy import pi, array, argmin, abs
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import GWFrames

# Read in the data files
Ws = [None]*len(Files)
for i,File in enumerate(Files):
    print("Reading '{0}'".format(File))
    Ws[i] = GWFrames.ReadFromNRAR(File)
print("Finished!\n")
W0 = Ws[0]

# Transform each Waveform into its corotating frame
print("Transforming to corotating frames...")
for W in Ws:
    W.TransformToCorotatingFrame()
    W.AlignDecompositionFrameToModes(tmid)
xHat = GWFrames.Quaternion(0,1,0,0)
R0_mid = GWFrames.Squad(W0.Frame(), W0.T(), [tmid])[0]
for W in Ws[1:]:
    # Rotate by pi if necessary
    Ri_mid = GWFrames.Squad(W.Frame(), W.T(), [tmid])[0]
    if((R0_mid*xHat*R0_mid.inverse()).dot(Ri_mid*xHat*Ri_mid.inverse())<0) :
        W.RotateDecompositionBasis(GWFrames.exp(GWFrames.Quaternion(0,0,0,pi/2)))
    # Align corotating frame of `W[1:]` to corotating frame of `W[0]`
    W.AlignTimeAndFrame(W0, t1, t2);
print("Finished!\n")

# First, we interpolate everything onto the common subset of time steps:
Ti = max([W.T(0) for W in Ws])
Tf = min([W.T(W.NTimes()-1) for W in Ws])
R0 = W0.Frame()[(W0.T()>=Ti) & (W0.T()<=Tf)]
T0 = W0.T()[(W0.T()>=Ti) & (W0.T()<=Tf)]
Abs0 = W0.Interpolate(T0).Abs(W0.FindModeIndex(2,2))
Imid = argmin(abs(T0-tmid))
plt.figure(0)
plt.figure(1)
for i,W in enumerate(Ws[1:]):
    print("Plotting difference {0} of {1}...".format(i, len(Ws[1:])))
    Ti = W.T()[(W.T()>=Ti) & (W.T()<=Tf)]
    Ri = GWFrames.Squad(W.Frame(), W.T(), T0)
    RDelta = GWFrames.UnflipRotors(GWFrames.RDelta(R0, Ri, Imid))
    PhiDelta = array(GWFrames.angle(RDelta))
    plt.figure(0)
    plt.semilogy(T0, 2.0*PhiDelta, label=Legends[i])
    plt.figure(1)
    Abs_i = W.Interpolate(T0).Abs(W.FindModeIndex(2,2))
    plt.semilogy(T0, abs(Abs0-Abs_i)/Abs0, label=Legends[i])
print("Finished!\n")

plt.figure(0)
plt.xlim((tmin,tmax))
plt.ylim((1e-6,1))
plt.ylabel(r'$2\Phi_\Delta$')
# ylabel(r'$\Phi_{\Delta}\ =\ 2\, \left| \log \left( R_{1}\, \bar{R}_{2} \right) \right|$')
plt.xlabel(r'$(t-r_\ast)/M$')
plt.title('Phase difference between waveforms')
plt.legend(loc='upper left')
plt.savefig(PhaseDifferencePlot)

plt.figure(1)
plt.xlim((tmin,tmax))
plt.ylim((1e-6,1))
plt.xlabel(r'$(t-r_\ast)/M$')
plt.ylabel(r'$\left(|h^{2,2}_1|-|h^{2,2}_2|\right)\,/\,|h^{2,2}_1|$')
plt.title(r'Relative amplitude differences between $(2,2)$ modes')
plt.legend(loc='upper left')
plt.savefig(AmplitudeDifferencePlot);
