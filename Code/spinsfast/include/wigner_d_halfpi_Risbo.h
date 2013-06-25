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
#include <string.h>

#ifndef M_SQRT1_2
# define M_SQRT1_2	0.70710678118654752440	/* 1/sqrt(2) */
#endif

typedef struct {
  double j;
  int twicej;
  int jmax, n, nmax;

  

  // double p,q,pc,qc;
  double * restrict sqt;
  
  double * restrict d;
  double * restrict dd;
} wdhp;


void wdhp_free(wdhp *wd);
wdhp *wdhp_alloc(double jmax);
void wdhp_reset(wdhp *wd);
wdhp *wdhp_init(double jmax);
void wdhp_jplushalf(wdhp *wd);
void wdhp_jplus1(wdhp *wd);
double wdhp_getj(wdhp *wd);
double wdhp_get(wdhp *wd, double m1, double m2);
double wdhp_integer_get(wdhp *wd, int m1, int m2);
double *wdhp_integer_getrow(wdhp *wd, int m1);
int wdhp_integer_idx(int l, int m1, int m2);
int wdhp_integer_N(int lmax);
double *wdhp_integer_precompute(int lmax);
