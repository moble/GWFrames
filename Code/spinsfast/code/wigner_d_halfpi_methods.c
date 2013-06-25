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

#include <wigner_d_halfpi.h>

void Delta_initialize(int DeltaMethod,void * Deltawork) {
 // If Delta not precomputed, initialize it here
  if (DeltaMethod == WDHP_METHOD_RISBO) {
    wdhp_reset((wdhp *)Deltawork);
  } else if (DeltaMethod == WDHP_METHOD_TN) {
    
  } else if (DeltaMethod == WDHP_METHOD_TN_PLANE) {

  }
}



void Delta_getplane( int DeltaMethod, void * Deltawork, const double * restrict Deltal, int l) {

  // For supported methods, grab or compute the l-plane of Delta matrix
  if (DeltaMethod == WDHP_METHOD_RISBO_PRECOMPUTE) { 
    Deltal = &((double *)Deltawork)[wdhp_integer_idx(l, 0, 0)];
  } else if (DeltaMethod == WDHP_METHOD_TN_PLANE) {
    wdhp_get_quarter_plane(l, ((wdhp_TN_helper *)Deltawork)->sqt, ((wdhp_TN_helper *)Deltawork)->invsqt, ((wdhp_TN_helper *)Deltawork)->D_all_llm, ((wdhp_TN_helper *)Deltawork)->Dwork);
  }


}




double *Delta_getrow( int DeltaMethod, void * Deltawork, const double * restrict Deltal, int l,int twicelp1, int mp) {

  double *Delta_mp = NULL;

  // Grab/compute the mp row (a 1-d array) of the Wigner-d Delta matrix.
  if (DeltaMethod == WDHP_METHOD_RISBO) {
    Delta_mp = wdhp_integer_getrow((wdhp *)Deltawork,mp);     
  } else if (DeltaMethod == WDHP_METHOD_RISBO_PRECOMPUTE) {
    Delta_mp = &Deltal[mp*twicelp1];     
  } else if (DeltaMethod == WDHP_METHOD_TN) {
    wdhp_get_row_pos(l,mp, ((wdhp_TN_helper *)Deltawork)->sqt, ((wdhp_TN_helper *)Deltawork)->invsqt, ((wdhp_TN_helper *)Deltawork)->D_all_llm, ((wdhp_TN_helper *)Deltawork)->Dwork);
    Delta_mp = ((wdhp_TN_helper *)Deltawork)->Dwork;
  } else if (DeltaMethod == WDHP_METHOD_TN_PLANE) {
    Delta_mp = &(((wdhp_TN_helper *)Deltawork)->Dwork[mp*(l+1)]);
  }
  
  return(Delta_mp);
}







void Delta_increment_l( int DeltaMethod, void * Deltawork) {

  // Increment Delta to next l if Risbo not precomputed
  if ( (DeltaMethod==WDHP_METHOD_RISBO) ) {  
    wdhp_jplus1((wdhp *)Deltawork); 
  }
  
}
