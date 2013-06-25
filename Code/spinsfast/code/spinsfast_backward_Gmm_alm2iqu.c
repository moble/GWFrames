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

#include <spinsfast_backward.h>


//inline 
int spinsfast_backward_sign_parity(int m);


void spinsfast_backward_Gmm_alm2iqu(const fftw_complex * restrict T, const fftw_complex * restrict P2,const int lmax, fftw_complex * restrict Gmm_I, fftw_complex * restrict Gmm_P, int DeltaMethod, void *Deltawork) {
 int l,m,mp;


 int Nm = 2*lmax+1;
 double norml;
 
 //int Nlm = N_lm(lmax);
 int NGmm = Nm*Nm;
 
 for (m=0;m<NGmm;m++) {
   Gmm_I[m] = 0;
   Gmm_P[m] = 0;
 }
 
 fftw_complex *Ito_helper = fftw_malloc(Nm*sizeof(complex));
 fftw_complex *Ito = &Ito_helper[lmax];
 
 for (m=-lmax; m<=lmax; m++){
   Ito[m] = cpow(I,m);
 }
 
 int *midx_helper = calloc(Nm,sizeof(int));
 for (m=-lmax; m<=lmax; m++){
   midx_helper[m+lmax] = (Nm + m) % Nm;
 }
 int * restrict midx = &midx_helper[lmax];
 
 
  // Set up Wigner Deltas
  // If Delta not precomputed, initialize it here
  Delta_initialize(DeltaMethod,Deltawork);

  
  // Compute Gm'm for m' >= 0
  for (l=0;l<=lmax;l++) {
    
    // For supported methods, grab or compute the l-plane of Delta matrix
    const double * restrict Deltal = NULL;
    Delta_getplane(DeltaMethod, Deltawork, Deltal, l);
    
    // remember to protect spin 2 from less than l=2
    {
      // shift to correct blocks
      const complex * restrict Tl = &T[lm_ind(l,0,lmax)];
      const complex * restrict P2l = &P2[lm_ind(l,0,lmax)];
      
      
      
      
      
      const int twicelp1 = 2*l+1;
      norml = sqrt(twicelp1)/2./sqrt(M_PI);
      const int signnegm = spinsfast_backward_sign_parity(l);
      
      for (mp=0; mp<=l; mp++){
	
	const double * restrict Delta_mp = NULL;
	
	// Grab the proper row (a 1-d array) of the Wigner-d Delta matrix, 
	// however it has been computed.
	Delta_mp = Delta_getrow( DeltaMethod, Deltawork, Deltal, l,twicelp1, mp);
	
	//     if (l==17) {
	//  printf("%d %e\n",mp,Delta_mp[0]);
	// }
	
	
	
	int mpmod = midx[mp]; 
	
	complex * restrict Gmp_I = &Gmm_I[mpmod*Nm];
	complex * restrict Gmp_P = &Gmm_P[mpmod*Nm];
	//     double *Delta_minusmp = &Delta[wdhp_integer_idx(l, -mp, 0)];
	
	// Get Delta_{mp -s} from Delta_{mp |s|}
	const double Deltamp0_norml = Delta_mp[0] * norml;
	const double Deltamp2_norml = Delta_mp[2] * norml * spinsfast_backward_sign_parity(l+mp);

	Gmp_I[midx[0]] += Delta_mp[0] * Deltamp0_norml * Tl[0];
	Gmp_P[midx[0]] += Delta_mp[0] * Deltamp2_norml * P2l[0];

	for (m=1; m<=l; m++){
	  //for (m=-l; m<=l; m++){
	  const double Delta_mpm = Delta_mp[m];
	  
	  //printf("%e %e\n",Delta_minusmp[m],Delta_minusmp1[m]);
	  
	  const double fact0 = (Delta_mpm * Deltamp0_norml);
	  const double fact2 = (Delta_mpm * Deltamp2_norml);

	  Gmp_I[midx[m]] +=  fact0*Tl[m];
	  //	  Gmp_I[midx[-m]] +=  fact0*Tl[-m]*signnegm;
	  
	  Gmp_P[midx[m]] +=  fact2*P2l[m];
	  Gmp_P[midx[-m]] +=  fact2*P2l[-m]*signnegm;

	}
	
      }
      
    }
  
  
    // Increment Delta to next l if Risbo not precomputed
    if (l<lmax) {
      if ( (DeltaMethod==WDHP_METHOD_RISBO) ) {  
	Delta_increment_l(DeltaMethod, Deltawork);
      }
    }
    
  }

  // Set the phase for m' >= 0
  for (mp=0; mp<=lmax; mp++){
    int mpmod = midx[mp];
    complex *Gmp_I = &Gmm_I[mpmod*Nm];
    complex *Gmp_P = &Gmm_P[mpmod*Nm];
    
    //   for (m=0; m<=lmax; m++){
    for (m=-lmax; m<=lmax; m++){
      int mmod = midx[m];
      Gmp_I[mmod]*= Ito[m];// (-i)^m * (-i)^0
      Gmp_P[mmod]*= -Ito[m];// (-i)^m * (-i)^2
    }
    
    // really I have no clue why this is here.
    // Maybe for delta_mp-s to delta_mps repair?
    for (m=0; m<=lmax; m++){
      int mmod = midx[m];
      Gmp_I[mmod]*=spinsfast_backward_sign_parity(m);
      Gmp_P[mmod]*=spinsfast_backward_sign_parity(m);
    }
    // no clue again...
    for (m=-lmax; m<0; m++){
      int mmod = midx[m];   
      Gmp_I[mmod]*= spinsfast_backward_sign_parity(mp+m);
      Gmp_P[mmod]*= spinsfast_backward_sign_parity(mp+m);
    }
  } 
    
  /* For intensity (real field), use G_(-m')(-m) = G_m'm^*   */
  for (mp=0; mp<=lmax; mp++){
    int mpmod = midx[mp];
    int negmpmod = midx[-mp];
    
    complex *Gmp_I = &Gmm_I[mpmod*Nm];
    complex *Gnegmp_I = &Gmm_I[negmpmod*Nm];

    for (m=0; m<=lmax; m++){
      int mmod = midx[m];
      int negmmod = midx[-m];
      
      Gnegmp_I[negmmod] = creal(Gmp_I[mmod]) - I*cimag(Gmp_I[mmod]);      
    }
    
  }

    
  /*  // Use symmetry G_(-m')m = (-1)^(m+s) G_m'm */
  const int sign_helper[3] = {-1,1,-1};
  const int *sign = &sign_helper[1];
  
  for (mp=0; mp<=lmax; mp++){
    int mpmod = midx[mp];
    int negmpmod = midx[-mp];
    
    complex *Gmp_I = &Gmm_I[mpmod*Nm];
    complex *Gmp_P = &Gmm_P[mpmod*Nm];
    
    complex *Gnegmp_I = &Gmm_I[negmpmod*Nm];
    complex *Gnegmp_P = &Gmm_P[negmpmod*Nm];
    
    for (m=-lmax; m<=lmax; m++){
      int mmod = midx[m];
      
      const int signm = sign[(m)&1];
      
      if (m>=0) {
	Gnegmp_I[mmod] = signm*Gmp_I[mmod];
      } else {
	Gmp_I[mmod] = signm*Gnegmp_I[mmod];
      }



      Gnegmp_P[mmod] = signm*Gmp_P[mmod];
    }
  }
  
  
 


 
  free(midx_helper);
  
  free(Ito_helper);
}



































