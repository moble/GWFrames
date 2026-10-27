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

int spinsfast_forward_sign_parity(int m);


void spinsfast_forward_transform_eo(fftw_complex * restrict a, const int Ntransform, const int *spins,const int lmax, fftw_complex * restrict Jmm_set, int DeltaMethod, void *Deltawork) {
  int l,m,mp;
  
  printf("Split eo!\n");
  
  int Nm = 2*lmax+1;
  double norml;
  
 
  int s;
  int ispin;

  int *midx_helper = calloc(Nm,sizeof(int));
  for (m=-lmax; m<=lmax; m++){
    midx_helper[m+lmax] = (Nm + m) % Nm;
  }
  int * restrict midx = &midx_helper[lmax];
  
  int Nlm = N_lm(lmax);
  int NJmm = (lmax+1)*Nm;
  
  for (m=0;m<Nlm*Ntransform;m++) {
    a[m] = 0;
  }
  
  fftw_complex *Z_ = calloc(Nm,sizeof(fftw_complex));
  fftw_complex *O_ = calloc(Nm,sizeof(fftw_complex));
  fftw_complex *E_ = calloc(Nm,sizeof(fftw_complex));
  
  fftw_complex *Z = &Z_[lmax];
  fftw_complex *O = &O_[lmax];
  fftw_complex *E = &E_[lmax];


  // Set up Wigner Deltas
  // If Delta not precomputed, initialize it here
  Delta_initialize(DeltaMethod,Deltawork);


  //  Main loop over multipole l
  for (l=0;l<=lmax;l++) {
    
    // For supported methods, grab or compute the l-plane of Delta matrix
    const double * restrict Deltal = NULL;
    Delta_getplane(DeltaMethod, Deltawork, Deltal, l);
  
	for (ispin = 0; ispin < Ntransform; ispin++) {
	  s = spins[ispin];
	  
	  if (l >= abs(s)) { 
	    // shift to correct blocks
	    complex * restrict asl = &a[ispin*Nlm + lm_ind(l,0,lmax)];
	    fftw_complex * restrict Jmm = &Jmm_set[ispin*NJmm];
	    
	    for (m=-l;m<=l;m++) {
	      Z[m] = 0;
	      O[m] = 0;
	      E[m] = 0;
	    }

	    
	    const int twicelp1 = 2*l+1;
	    norml = sqrt(twicelp1)/2./sqrt(M_PI);
	    //	    const int negtol =  spinsfast_forward_sign_parity(l); // = (-1)^l
	    

	    // Zero term
	    mp = 0;
	    {
	      const int signnegm = spinsfast_forward_sign_parity(l); // = (1-)^(l+mp)
	      
	      const double * restrict Delta_mp = NULL;
	      
	      // Grab/compute the mp row (a 1-d array) of the Wigner-d Delta matrix.
	      Delta_mp = Delta_getrow( DeltaMethod, Deltawork, Deltal, l,twicelp1, mp);
	      
	      // Get Delta_{mp s} from Delta_{mp |s|}
	      const int s_sign_fudge = (s<0) ? 1 : spinsfast_forward_sign_parity(l);
	      const double Deltamps_norml = Delta_mp[abs(s)] * norml * s_sign_fudge;
	      
	      const complex * restrict Jmp = &Jmm[mp*Nm];
	      
	      if (l<5){
		for (m=-l; m<=l; m++){
		  printf("% .2e % .2e | ",creal(Jmp[midx[m]]),cimag(Jmp[midx[m]]));
		}
		printf("\n");
	      }

	      for (m=0; m<=l; m++){
	    	const double Delta_mpm = Delta_mp[m];
	    	const complex Jmpm = Jmp[midx[m]];
	    	const complex Jmpnegm = Jmp[midx[-m]];
		
       	    	const double fact = (Delta_mpm * Deltamps_norml);
				
	    	Z[m] +=  fact*Jmpm;
	    	Z[-m] += fact*Jmpnegm*signnegm;

		/* if (1) { */
		/*   		  printf("lmpm %d %d %d: Jmpm % .2e % .2e Jmp(-m) % .2e % .2e | Zm % e % e  Z(-m) % e % e\n",l,mp,m, */
		/* 	 creal(Jmpm),cimag(Jmpm), */
		/* 	 creal(Jmpnegm),cimag(Jmpnegm), */
		/* 	 creal(Z[m]),cimag(Z[m]), */
		/* 	 creal(Z[-m]),cimag(Z[-m])); */

		/* } */


	      }
	    }

	    for (mp=1; mp<=l; mp+=2){ // Odd mp > 0
	      const int signnegm = spinsfast_forward_sign_parity(l+1); // = (1-)^(l+mp)
	      
	      const double * restrict Delta_mp = NULL;
	      
	      // Grab/compute the mp row (a 1-d array) of the Wigner-d Delta matrix.
	      Delta_mp = Delta_getrow( DeltaMethod, Deltawork, Deltal, l,twicelp1, mp);
	      
	      // Get Delta_{mp s} from Delta_{mp |s|}
	      const int s_sign_fudge = (s<0) ? 1 : spinsfast_forward_sign_parity(l+1);
	      const double Deltamps_norml = Delta_mp[abs(s)] * norml * s_sign_fudge;
	      
	      const complex * restrict Jmp = &Jmm[mp*Nm];
	      
	      for (m=0; m<=l; m++){
	    	const double Delta_mpm = Delta_mp[m];
	    	const complex Jmpm = Jmp[midx[m]];
	    	const complex Jmpnegm = Jmp[midx[-m]];
		
       	    	const double fact = (Delta_mpm * Deltamps_norml);
		
	    	O[m] +=  fact*Jmpm;
	    	O[-m] += fact*Jmpnegm*signnegm;
	      }
	    }
	    
	    for (mp=2; mp<=l; mp+=2){ // Even mp > 0
	      const int signnegm = spinsfast_forward_sign_parity(l); // = (1-)^(l+mp)
	      
	      const double * restrict Delta_mp = NULL;
	      
	      // Grab/compute the mp row (a 1-d array) of the Wigner-d Delta matrix.
	      Delta_mp = Delta_getrow( DeltaMethod, Deltawork, Deltal, l,twicelp1, mp);
	      
	      // Get Delta_{mp s} from Delta_{mp |s|}
	      const int s_sign_fudge = (s<0) ? 1 : spinsfast_forward_sign_parity(l);
	      const double Deltamps_norml = Delta_mp[abs(s)] * norml * s_sign_fudge;
	      
	      const complex * restrict Jmp = &Jmm[mp*Nm];
	      
	      for (m=0; m<=l; m++){
		const double Delta_mpm = Delta_mp[m];
		const complex Jmpm = Jmp[midx[m]];
		const complex Jmpnegm = Jmp[midx[-m]];
		
       		const double fact = (Delta_mpm * Deltamps_norml);
		
		E[m] +=  fact*Jmpm;
		E[-m] += fact*Jmpnegm*signnegm;
	      }	     
	    }

	    if (1) {
	      for (m=-l;m<=l;m++) {
		printf("% d % d: ",l,m);
		printf("Z % e % e | ",creal(Z[m]),cimag(Z[m]));
		printf("O % e % e | ",creal(O[m]),cimag(O[m]));
		printf("E % e % e |\n",creal(E[m]),cimag(E[m]));
	      }
	      printf("\n");
	    }
	    

	    for (m=-l;m<=l;m++) {
	      //	      asl[m] = 0;
	      asl[m] = Z[m] + E[m] + O[m];
	    }

	  } // We are now done looping over m' and m 
	}   //
	
	// Increment Delta to next l if Risbo not precomputed
	if (l<lmax) {
	  if ( (DeltaMethod==WDHP_METHOD_RISBO) ) {  
	    Delta_increment_l(DeltaMethod, Deltawork);
	  }
	}
  }
  
  
  // Set the phases
   fftw_complex *Itom_helper = fftw_malloc(Nm*sizeof(complex));
  fftw_complex *Itom = &Itom_helper[lmax];
  
  for (m=-lmax; m<=lmax; m++){
    Itom[m] = cpow(I,m);
  }

  fftw_complex negItos[Ntransform];
  for (ispin = 0; ispin < Ntransform; ispin++) {
    s = spins[ispin];
    negItos[ispin] = cpow(-I,s);
    
    for (l=s;l<=lmax;l++) {
      complex * restrict asl = &a[ispin*Nlm + lm_ind(l,0,lmax)];
      
      asl[0] /= 2;

      for (m=-l; m<=l; m++){
	asl[m] *= Itom[m]*negItos[ispin];
      }
    }
  }
  
  free(midx_helper);
  
  free(Itom_helper);
}
