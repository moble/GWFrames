## This script uses the spinsfast package <http://www.physics.miami.edu/~huffenbe/research/spinsfast/>

from numpy import loadtxt, zeros, complex, concatenate, linspace, array, exp, pi, sin, cos, sqrt
from numpy.random import normal,seed
from matplotlib.pyplot import *
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.cm import jet

import spinsfast

# This loads a dat file, though it wouldn't be hard to load an h5 file
# using h5py.
Data = zeros((154,))

# Now select the particular time that you want to be plotted.
TimeIndex = 1234
ReIm = Data[:]

# Since this is a spin-weight -2 field, SpEC doesn't bother to store
# the (trivially zero) l=0 and l=1 modes.  But spinsfast requires
# them, so we prepend those four modes to our data here, while making
# the rest of the data complex numbers.
alm = concatenate(([sqrt(4*pi),0.,0.,0.+0.*1j], [re + 1j*im for re,im in zip(ReIm[1::2], ReIm[2::2])]))

# Set up the parameters of the transformation.  These Ntheta and Nphi
# are usually overkill for plotting purposes, but don't slow things
# down noticeably.
s = 0 # Spin weight
lmax = 7
Ntheta = 2**9+1
Nphi = 2**9
Nlm = spinsfast.N_lm(lmax);

# These are the physical coordinates to be used below.
theta = linspace(0, pi, num=Ntheta)
phi = linspace(0, 2*pi, num=Nphi)

# Here, we transform the data, take its absolute value, and then
# construct a colormap from its normalized values.
f =  spinsfast.salm2map(alm,s,lmax,Ntheta,Nphi)
absf = abs(f)

# This is for 2-d plots.
pcolormesh(phi, theta, absf)
xlim((0,2*pi))
ylim((0,pi))
xlabel(r'$\phi$')
ylabel(r'$\theta$')
tight_layout()
#savefig('Waveform.png', dpi=200, transparent=True)
#show()

# This is for 3-d plots.  This is by far the slowest part of the
# process.  Increase rstride and cstride to make this go faster at the
# cost of resolution.
# x = array([[sin(t)*cos(p) for p in phi] for t in theta])
# y = array([[sin(t)*sin(p) for p in phi] for t in theta])
# z = array([[cos(t) for p in phi] for t in theta])
# Normalization = absf.max()
# Colors = jet(absf/Normalization)
# fig = figure()
# ax = fig.add_subplot(111, projection='3d')
# ax.set_xticks(())
# ax.set_yticks(())
# ax.set_zticks(())
# ax.plot_surface(x, y, z, rstride=1, cstride=1, linewidth=0, antialiased=False, facecolors=Colors)
# tight_layout()
# savefig('Waveform.png', dpi=200, transparent=True)
# show()
