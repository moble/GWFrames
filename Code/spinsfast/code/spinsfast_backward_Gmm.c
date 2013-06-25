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
int spinsfast_backward_sign_parity(int m){
  // returns  (-1)^m, 
  // = 1 if even
  // = -1 if odd.

  int eo = (m & 1);

  return( 1 - eo - eo );
}


void spinsfast_backward_Gmm(const fftw_complex * restrict a, int Ntransform, const int *spins,const int lmax, fftw_complex * restrict Gmm_set, int DeltaMethod, void *Deltawork) {
 int l,m,mp;


 int Nm = 2*lmax+1;
 double norml;
 
 int Nlm = N_lm(lmax);
 int NGmm = Nm*Nm;
 
 for (m=0;m<NGmm*Ntransform;m++) {
   Gmm_set[m] = 0;
 }
 
 int ispin, s;

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
   
   
   for (ispin = 0; ispin < Ntransform; ispin++) {
     s = spins[ispin];

     if (l>=s) {
       // shift to correct blocks
       const complex * restrict asl = &a[ispin*Nlm + lm_ind(l,0,lmax)];
       fftw_complex * restrict Gmm = &Gmm_set[ispin*NGmm];
       
       int negtol = spinsfast_backward_sign_parity(l);
       
       
       
       
       const int twicelp1 = 2*l+1;
       norml = sqrt(twicelp1)/2./sqrt(M_PI);
       
       for (mp=0; mp<=l; mp++){
	 
	 const double * restrict Delta_mp = NULL;
	 
	 // Grab the proper row (a 1-d array) of the Wigner-d Delta matrix, 
	 // however it has been computed.
	 Delta_mp = Delta_getrow( DeltaMethod, Deltawork, Deltal, l,twicelp1, mp);
	 
	 //     if (l==17) {
	 //  printf("%d %e\n",mp,Delta_mp[0]);
	 // }
	 
	 
	 
	 int mpmod = midx[mp]; 
	 
	 complex * restrict Gmp = &Gmm[mpmod*Nm];
	 //     double *Delta_minusmp = &Delta[wdhp_integer_idx(l, -mp, 0)];
	 
	 // Get Delta_{mp s} from Delta_{mp |s|}
	 const int s_sign_fudge = (s>=0) ? 1 : spinsfast_backward_sign_parity(l+mp);
	 const double Deltamps_norml = s_sign_fudge * Delta_mp[abs(s)] * norml ;
	 const double Deltamps_norml_negtol = Deltamps_norml * negtol; //  includes 1^(-1)
	 
	 // printf("l = %d -mp = %d\n",l,-mp);
	 
	 Gmp[midx[0]] += Delta_mp[0] * Deltamps_norml_negtol * asl[0];
	 //#pragma omp parallel for private(m)
	 for (m=1; m<=l; m++){
	   //      for (m=-l; m<=l; m++){
	   const int mmod = midx[m];
	   const int negmmod = midx[-m];
	   const double Delta_mpm = Delta_mp[m];
	   const complex aslm = asl[m];
	   const complex aslnegm = asl[-m];
	   
	   //printf("%e %e\n",Delta_minusmp[m],Delta_minusmp1[m]);
	   
	   
	   const double fact = (Delta_mpm * Deltamps_norml_negtol);
	   const double factneg = (Delta_mpm * Deltamps_norml);
	   
	   const complex term = fact*aslm;
	   const complex termneg = factneg*aslnegm;
	   
	   Gmp[mmod] +=  term;
	   Gmp[negmmod] +=  termneg;
	   
	   /* if (l<6) { */
	   /* 	 printf("l mp m % d % d % d | % e % e | % e % e, % e % e\n",l,mp,m, fact, factneg, creal(term), cimag(term), creal(termneg), cimag(termneg)); */
	   /* } */
	   
	 }
	 
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
 
 
 for (ispin = 0; ispin < Ntransform; ispin++) {
   s = spins[ispin];
   fftw_complex * restrict Gmm = &Gmm_set[ispin*NGmm];
   
   // Set the phase for m' >= 0
   for (mp=0; mp<=lmax; mp++){
     int mpmod = midx[mp];
     complex *Gmp = &Gmm[mpmod*Nm];
     
     //   for (m=0; m<=lmax; m++){
     for (m=-lmax; m<=lmax; m++){
       int mmod = midx[m];
       Gmp[mmod]*= Ito[s]*Ito[m];//*spinsfast_backward_sign_parity(mp+m);
     }
     
     for (m=0; m<=lmax; m++){
       int mmod = midx[m];
       Gmp[mmod]*=spinsfast_backward_sign_parity(mp+m);
     }
     for (m=-lmax; m<0; m++){
       int mmod = midx[m];   
       Gmp[mmod]*= spinsfast_backward_sign_parity(m);
     }
   }
   
   
   /*  // Use symmetry G_(-m')m = (-1)^(m+s) G_m'm */
   const int sign_helper[3] = {-1,1,-1};
   const int *sign = &sign_helper[1];
   
   for (mp=0; mp<=lmax; mp++){
     int mpmod = midx[mp];
     int negmpmod = midx[-mp];
     complex *Gmp = &Gmm[mpmod*Nm];
     complex *Gnegmp = &Gmm[negmpmod*Nm];
     
     for (m=-lmax; m<=lmax; m++){
       int mmod = midx[m];
       
       const int signm = sign[(m+s)&1];
       
       Gnegmp[mmod] = signm*Gmp[mmod];
     }
   }
   
 }
 
 free(midx_helper);
 
 free(Ito_helper);
}



































