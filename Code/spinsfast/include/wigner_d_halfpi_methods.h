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

#define WDHP_METHOD_RISBO 0
#define WDHP_METHOD_RISBO_PRECOMPUTE 1
#define WDHP_METHOD_TN 2
#define WDHP_METHOD_TN_PLANE 3
//#define WDHP_METHOD_TN_PRECOMPUTE 4


void Delta_initialize(int DeltaMethod,void * Deltawork);
void Delta_getplane( int DeltaMethod, void * Deltawork, const double * restrict Deltal, int l);
double *Delta_getrow( int DeltaMethod, void * Deltawork, const double * restrict Deltal, int l,int twicelp1, int mp);
void Delta_increment_l( int DeltaMethod, void * Deltawork);
