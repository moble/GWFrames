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

#include<math.h>

// The following tries to deal with the `restrict` keyword used by the
// spinsfast module.
#ifndef restrict
#ifdef __restrict
#define restrict __restrict
#endif
#endif
#ifndef restrict
#if __STDC_VERSION__ < 199901L
#define restrict /* nothing */
#endif
#endif

int lm_ind(int l, int m, int lmax);
void ind_lm(int i, int *l, int *m, int lmax);
int N_lm(int lmax);
