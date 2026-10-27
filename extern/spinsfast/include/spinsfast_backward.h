/**************************************************************************

    Copyright 2010-2012  Kevin M. Huffenberger & Benjamin D. Wandelt

    This file is part of spinsfast.

    Spinsfast is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Spinsfast is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with spinsfast.  If not, see <http://www.gnu.org/licenses/>.

***************************************************************************/

/* Code revision: 104, 2012-04-13 13:00:16 -0400 (Fri, 13 Apr 2012) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include<fftw3.h>
#include <wigner_d_halfpi.h>
#include <alm.h>

#ifndef M_PI
# define M_PI		3.14159265358979323846	/* pi */
#endif

void spinsfast_backward_transform(fftw_complex  *f, int Ntheta, int Nphi, int lmax, fftw_complex *Gmm);
void spinsfast_salm2map(fftw_complex *alm, fftw_complex *f, int s, int Ntheta, int Nphi, int lmax);

void spinsfast_backward_Gmm(const fftw_complex * restrict a, int Ntransform, const int *spins,const int lmax, fftw_complex * restrict Gmm_set, int DeltaMethod, void *Deltawork);
//void spinsfast_backward_Gmm(const fftw_complex * restrict a, const int s,const int lmax, fftw_complex * restrict Gmm, int DeltaMethod, void *Deltawork);
void spinsfast_backward_Gmm_old(const fftw_complex * restrict a, const int s, const int lmin,const int lmax, fftw_complex * restrict Gmm);
void spinsfast_backward_Gmm_preDelta(const fftw_complex * restrict a, const int s, const int lmin,const int lmax, fftw_complex * restrict Gmm,const double * restrict Delta);
void spinsfast_backward_Gmm_alm2iqu(const fftw_complex * restrict T, const fftw_complex * restrict P2,const int lmax, fftw_complex * restrict Gmm_I, fftw_complex * restrict Gmm_P, int DeltaMethod, void *Deltawork);
