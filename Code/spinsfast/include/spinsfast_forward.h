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
#include <fftw3.h>
#include <wigner_d_halfpi.h>
#include <alm.h>

#ifndef M_PI
# define M_PI		3.14159265358979323846	/* pi */
#endif

void spinsfast_map2salm(fftw_complex *f, fftw_complex *alm, int s, int Ntheta, int Nphi, int lmax);



int indx2p (int ip, int Ntheta);
int p2indx (int p, int Ntheta);
int mpm2indx( int mp, int m, int Ntheta, int Nphi );

void spinsfast_quadrature_weights(fftw_complex *W, int wsize);


void spinsfast_f_extend_MW(fftw_complex *f, fftw_complex *F, int s, int Ntheta, int Nphi);
void spinsfast_f_extend_old(fftw_complex *f, fftw_complex *F, int s, int Ntheta, int Nphi);



void spinsfast_forward_Imm (fftw_complex *f, int s, int Ntheta, int Nphi, int lmax, fftw_complex *Imm);
void spinsfast_forward_multi_Imm (fftw_complex *f_set, int *spins, int Nmap, int Ntheta, int Nphi, int lmax, fftw_complex *Imm_set);

void spinsfast_forward_Jmm (fftw_complex *f, int s, int Ntheta, int Nphi, int lmax, fftw_complex *Jmm);
void spinsfast_forward_multi_Jmm(fftw_complex *f_set, int *spins, int Nmap, int Ntheta, int Nphi, int lmax, fftw_complex *Jmm_set);

void spinsfast_forward_transform(fftw_complex * restrict a, const int Ntransform, const int *spins,const int lmax, fftw_complex * restrict Jmm_set, int DeltaMethod, void *Deltawork);
void spinsfast_forward_transform_eo(fftw_complex * restrict a, const int Ntransform, const int *spins,const int lmax, fftw_complex * restrict Jmm_set, int DeltaMethod, void *Deltawork);
//void spinsfast_forward_transform(fftw_complex * restrict a, const int s,const int lmax, fftw_complex * restrict Jmm, int pre, void *Deltawork);
void spinsfast_forward_transform_from_Imm(fftw_complex * restrict a, const int s,const int smax,const int lmax, fftw_complex * restrict Imm, int pre, void *Deltawork);
void spinsfast_forward_transform_preDelta(fftw_complex * restrict a, const int s,const int smax,const int lmax, fftw_complex * restrict Imm,const double * restrict Delta) ;
void spinsfast_forward_transform_iqu2alm(fftw_complex * restrict T,fftw_complex * P2,const int lmax, fftw_complex * restrict Jmm_I, fftw_complex * restrict Jmm_P, int DeltaMethod, void *Deltawork);
