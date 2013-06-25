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
#include <chealpix.h>
#include <complex.h>

#ifndef M_PI
# define M_PI		3.14159265358979323846	/* pi */
#endif

void ecp2healpix_avg(complex *ecp, int Ntheta, int Nphi,
		     double *hp_r, double *hp_i, int nside);

void ecp2healpix_nearest(complex *ecp, int Ntheta, int Nphi,
			 double *hp_r, double *hp_i, int nside);


void healpix2ecp(double *hp_r, double *hp_i, int nside,
	    complex *ecp, int Ntheta, int Nphi);
