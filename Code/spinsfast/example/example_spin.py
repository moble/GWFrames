#/**************************************************************************
#
#    Copyright 2010-2012  Kevin M. Huffenberger & Benjamin D. Wandelt
#
#    This file is part of spinsfast.
#
#    Spinsfast is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Spinsfast is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with spinsfast.  If not, see <http://www.gnu.org/licenses/>.
#
#***************************************************************************/
#
#/* Code revision: 104, 2012-04-13 13:00:16 -0400 (Fri, 13 Apr 2012) */

import sys
sys.path.append('lib')

from numpy import zeros, complex
from numpy.random import normal,seed
from matplotlib.pyplot import *

import spinsfast

s = 1
lmax = 127
Ntheta = 257
Nphi = 384  
# For best accuracy, have lmax < Nphi/2 and lmax < Ntheta/2
#
# For the FFT part of the code:
#      theta FFT is fastest for Ntheta = 2^n + 1
#      phi FFT is fastest for Nphi = 2^n
#
# But the FFT part is a subdominant to the overall scaling, so lmax is
# much more important to overall speed than number of pixels


Nlm = spinsfast.N_lm(lmax);
alm = zeros(Nlm,dtype=complex)


# Fill the alm with white noise
seed(3124432)
for l in range(abs(s),lmax+1) :
    for m in range(-l,l+1) :
        i = spinsfast.lm_ind(l,m,lmax)
        if (m==0) :
            alm[i] = normal()
        else :
            alm[i] = normal()/2 + 1j*normal()/2


f =  spinsfast.salm2map(alm,s,lmax,Ntheta,Nphi)
# In this pixelization, access the map with f[itheta,iphi]
# where 0 <= itheta <= Ntheta-1 and 0<= iphi <= Nphi-1
# and theta = pi*itheta/(Ntheta-1) phi = 2*pi*iphi/Nphi


alm2 = spinsfast.map2salm(f,s,lmax)

diff_max = max((alm-alm2))
print("max(alm2-alm) = "+str(diff_max))

figure()
imshow(f.real,interpolation='nearest')
colorbar()
title("Real part of f")
xlabel("iphi")
ylabel("itheta")
show()
