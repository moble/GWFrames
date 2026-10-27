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

#ifndef M_SQRT1_2l
#define M_SQRT1_2l	0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
#endif

#define FPTYPE long double
//#define FPTYPE double


// For use w/ pos m1 m2


// members of this structure can be used with  wdhp_get_row_pos and wdhp_get_col_pos
// after initialization with wdhp_TN_helper_init.
typedef struct {

  int lmax;
  FPTYPE * sqt;
  FPTYPE * invsqt;
  FPTYPE * D_all_llm;
  double * Dwork;


} wdhp_TN_helper;

wdhp_TN_helper *wdhp_TN_helper_init(int lmax);
void wdhp_TN_helper_free(wdhp_TN_helper *w);

int wdhp_lmind_pos(int l,int m);


FPTYPE wdhp_get_ll0(int l);
FPTYPE wdhp_get_llm(int l,int m);
FPTYPE wdhp_get_lm1m2_pos(int l,int m1,int m2, FPTYPE *sqt, FPTYPE *invsqt); // more accurate if m2 >= m1
FPTYPE wdhp_get_lm1m2_pos2(int l,int m1,int m2, FPTYPE *sqt, FPTYPE *invsqt); // always uses shorter recurrence for accuracy

FPTYPE *wdhp_init_sqt(int lmax);
FPTYPE *wdhp_init_invsqt(int lmax, FPTYPE *sqt);


void wdhp_get_all_llm(int lmax,FPTYPE *D_all_llm);
void wdhp_get_all_llm2(int lmax,FPTYPE *D_all_llm,FPTYPE *sqt, FPTYPE *invsqt);

void wdhp_get_col_pos(int l,int m2, FPTYPE *sqt, FPTYPE *invsqt, FPTYPE *D_all_llm, double *Dcol);
void wdhp_get_row_pos(int l,int m1, FPTYPE *sqt, FPTYPE *invsqt, FPTYPE *D_all_llm, double *Drow);
void wdhp_get_quarter_plane(int l, FPTYPE *sqt, FPTYPE *invsqt, FPTYPE *D_all_llm, double *D);
void wdhp_get_quarter_plane2(int l, FPTYPE *sqt, FPTYPE *invsqt, FPTYPE *D_all_llm, double *D); // This is same as wdhp_get_row_pos for m2>m1
