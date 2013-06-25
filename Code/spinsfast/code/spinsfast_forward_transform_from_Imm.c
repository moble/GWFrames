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

#include <spinsfast_forward.h>

void spinsfast_forward_transform_from_Imm(fftw_complex * restrict a, const int s,const int smax,const int lmax, fftw_complex * restrict Imm, int pre, void *Deltawork) {
  int l,m,mp;
  
  
  int Nm = 2*lmax+1;
  double norml;
  
  fftw_complex *Itom_helper = fftw_malloc(Nm*sizeof(complex));
  fftw_complex *Itom = &Itom_helper[lmax];
  
  for (m=-lmax; m<=lmax; m++){
    Itom[m] = cpow(I,m);
  }
  
  fftw_complex negItos = cpow(-I,s);
  
  int *midx_helper = calloc(Nm,sizeof(int));
  for (m=-lmax; m<=lmax; m++){
    midx_helper[m+lmax] = (Nm + m) % Nm;
  }
  int * restrict midx = &midx_helper[lmax];
  
  int N = N_lm(lmax);
  
  for (m=0;m<N;m++) {
    a[m] = 0;
  }
  
  
  wdhp *aDelta = NULL;

  // If Delta not precomputed, initialize it here
  if (pre == 0) {
    aDelta = (wdhp *)Deltawork;
    wdhp_reset(aDelta);
    for (l=0;l<s;l++) {
      if (l<lmax) wdhp_jplus1(aDelta);
    }
  }
  
  for (l=s;l<=lmax;l++) {
    complex * restrict asl = &a[lm_ind(l,0,lmax)];

    // if precomputed, grab the correct block from the precomputed matrix
    const double * restrict Deltal = NULL;
    if (pre == 1) Deltal = &((double *)Deltawork)[wdhp_integer_idx(l, 0, 0)];


    const int twicelp1 = 2*l+1;
    norml = sqrt(twicelp1)/2./sqrt(M_PI);
    
    for (mp=-l; mp<=l; mp++){
      int mpmod = midx[mp]; 
          
      const double * restrict Delta_mp = NULL;

      // Grab the proper row of the Wigner-d Delta matrix, on-fly or precomputed
      if (pre == 0) Delta_mp = wdhp_integer_getrow(aDelta,mp);     
      else Delta_mp = &Deltal[mp*twicelp1];     


      const double Deltamps_norml = Delta_mp[-s] * norml;
      
      complex * restrict Imp = &Imm[mpmod*Nm];
      
      for (m=-l; m<=l; m++){
	const int mmod = midx[m];
	const double Delta_mpm = Delta_mp[m];
	
	const double fact = (Delta_mpm * Deltamps_norml);
	
	const complex Impm = Imp[mmod];
	const complex term = fact*Impm;
	
	asl[m] +=  term;
	
      }
      
    }
    
    // Increment Delta to next l if not precomputed
    if ( (pre==0) && (l<lmax) ) wdhp_jplus1(aDelta);
  }
  
  //wdhp_free(aDelta);
  
  // Set the phase
  for (l=s;l<=lmax;l++) {
    complex * restrict asl = &a[lm_ind(l,0,lmax)];
    for (m=-l; m<=l; m++){
      asl[m] *= Itom[m]*negItos;
    }
  }
  
  free(midx_helper);
  
  free(Itom_helper);
}
