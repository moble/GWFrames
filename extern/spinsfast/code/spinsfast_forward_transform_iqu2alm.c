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

int isOdd(int l) {
  return(l&1);
}


void spinsfast_forward_transform_iqu2alm(fftw_complex * restrict T,fftw_complex * P2,const int lmax, fftw_complex * restrict Jmm_I, fftw_complex * restrict Jmm_P, int DeltaMethod, void *Deltawork) {
  int l,m,mp;
  //
  //  T, P2 are the output harmonics
  //  
  //  T is the spin=0 transform based on Jmm_I
  //  Jmm_I is Jmm computed from the intensity map
  //
  //  P2 is the spin=2 transform based on Jmm_P 
  //  Jmm_P is Jmm computed from the map of P = Q+iU
  //
  //  E/B harmonics can be computed from P2

  int Nm = 2*lmax+1;
  double norml;
  

  int *midx_helper = calloc(Nm,sizeof(int));
  for (m=-lmax; m<=lmax; m++){
    midx_helper[m+lmax] = (Nm + m) % Nm;
  }
  int * restrict midx = &midx_helper[lmax];
  
  int Nlm = N_lm(lmax);
  //  int NJmm = (lmax+1)*Nm;
  
  for (m=0;m<Nlm;m++) {
    T[m] = 0;
    P2[m] = 0;
  }
  
  fftw_complex *ZI_ = calloc(Nm,sizeof(fftw_complex));
  fftw_complex *OI_ = calloc(Nm,sizeof(fftw_complex));
  fftw_complex *EI_ = calloc(Nm,sizeof(fftw_complex));
  
  fftw_complex *ZI = &ZI_[lmax];
  fftw_complex *OI = &OI_[lmax];
  fftw_complex *EI = &EI_[lmax];

  fftw_complex *ZP_ = calloc(Nm,sizeof(fftw_complex));
  fftw_complex *OP_ = calloc(Nm,sizeof(fftw_complex));
  fftw_complex *EP_ = calloc(Nm,sizeof(fftw_complex));
  
  fftw_complex *ZP = &ZP_[lmax];
  fftw_complex *OP = &OP_[lmax];
  fftw_complex *EP = &EP_[lmax];

  // Set up Wigner Deltas
  // If Delta not precomputed, initialize it here
  Delta_initialize(DeltaMethod,Deltawork);


  //  Main loop over multipole l
  for (l=0;l<=lmax;l++) {
    
    // For supported methods, grab or compute the l-plane of Delta matrix
    const double * restrict Deltal = NULL;
    Delta_getplane(DeltaMethod, Deltawork, Deltal, l);
    
    complex * restrict Tl = &T[lm_ind(l,0,lmax)];
    complex * restrict P2l = &P2[lm_ind(l,0,lmax)];
    
    for (m=-l;m<=l;m++) {
      ZI[m] = 0;
      OI[m] = 0;
      EI[m] = 0;
      ZP[m] = 0;
      OP[m] = 0;
      EP[m] = 0;
    }

	    
    const int twicelp1 = 2*l+1;
    norml = sqrt(twicelp1)/2./sqrt(M_PI);
    //	    const int negtol =  spinsfast_forward_sign_parity(l); // = (-1)^l
	    

    // Zero term
    mp = 0;
   {
      const int signnegm = spinsfast_forward_sign_parity(l); // = (-1)^(l+mp)
      
      const double * restrict Delta_mp = NULL;
      
      // Grab/compute the mp row (a 1-d array) of the Wigner-d Delta matrix.
      Delta_mp = Delta_getrow( DeltaMethod, Deltawork, Deltal, l,twicelp1, mp);
      
      // Get Delta_{mp -s} from Delta_{mp |s|}
      const double Deltamp0_norml = Delta_mp[0] * norml;
      const double Deltamp2_norml = Delta_mp[2] * norml * spinsfast_forward_sign_parity(l);
      
      
      const complex * restrict Jmp_I = &Jmm_I[mp*Nm];
      const complex * restrict Jmp_P = &Jmm_P[mp*Nm];
      
      // Intensity, mp=0  (can skip by two)
      if (!isOdd(l)) { for (m=0; m<=l; m+=2){
	const double Delta_mpm = Delta_mp[m];
	const complex Jmpm_I = Jmp_I[midx[m]];
	//	const complex Jmpnegm_I = Jmp_I[midx[-m]];
	const double fact0 = (Delta_mpm * Deltamp0_norml);
	
	ZI[m] +=  fact0*Jmpm_I;
	//	ZI[-m] += fact0*Jmpnegm_I*signnegm;
	
      }
      }

      // Polarization, mp = 0;
      if (l>=2) 
	for (m=0; m<=l; m++){
	  const double Delta_mpm = Delta_mp[m];
	  const complex Jmpm_P = Jmp_P[midx[m]];
	  const complex Jmpnegm_P = Jmp_P[midx[-m]];
	  const double fact2 = (Delta_mpm * Deltamp2_norml);
	  
	  ZP[m] +=  fact2*Jmpm_P;
	  ZP[-m] += fact2*Jmpnegm_P*signnegm;
	}   
    }
    
   
      for (mp=1; mp<=l; mp+=2){ // Odd mp > 0
	const int signnegm = spinsfast_forward_sign_parity(l+1); // = (-1)^(l+mp)
	
      const double * restrict Delta_mp = NULL;
      
      // Grab/compute the mp row (a 1-d array) of the Wigner-d Delta matrix.
      Delta_mp = Delta_getrow( DeltaMethod, Deltawork, Deltal, l,twicelp1, mp);
      
      // Get Delta_{mp -s} from Delta_{mp |s|}
      const double Deltamp0_norml = Delta_mp[0] * norml;
      const double Deltamp2_norml = Delta_mp[2] * norml * spinsfast_forward_sign_parity(l+1);
      
      const complex * restrict Jmp_I = &Jmm_I[mp*Nm];
      const complex * restrict Jmp_P = &Jmm_P[mp*Nm];

      // I, mp odd
     if (isOdd(l)) {  for (m=0; m<=l; m++){
	const double Delta_mpm = Delta_mp[m];
	const complex Jmpm_I = Jmp_I[midx[m]];
	//	const complex Jmpnegm_I = Jmp_I[midx[-m]];

	const double fact0 = (Delta_mpm * Deltamp0_norml);
	
	OI[m] +=  fact0*Jmpm_I;
	//	OI[-m] += fact0*Jmpnegm_I*signnegm;
      }
     }
      
      // P, mp odd
      if (l>=2) 
	for (m=0; m<=l; m++){
	  const double Delta_mpm = Delta_mp[m];
	  const complex Jmpm_P = Jmp_P[midx[m]];
	  const complex Jmpnegm_P = Jmp_P[midx[-m]];
	  
	  const double fact2 = (Delta_mpm * Deltamp2_norml);
	  
	  OP[m] +=  fact2*Jmpm_P;
	  OP[-m] += fact2*Jmpnegm_P*signnegm;
	}
    }
    
	
   
      for (mp=2; mp<=l; mp+=2){ // Even mp > 0
      const int signnegm = spinsfast_forward_sign_parity(l); // = (-1)^(l+mp)
      
      const double * restrict Delta_mp = NULL;
      
      // Grab/compute the mp row (a 1-d array) of the Wigner-d Delta matrix.
      Delta_mp = Delta_getrow( DeltaMethod, Deltawork, Deltal, l,twicelp1, mp);
      
      // Get Delta_{mp -s} from Delta_{mp |s|}
      const double Deltamp0_norml = Delta_mp[0] * norml;
      const double Deltamp2_norml = Delta_mp[2] * norml * spinsfast_forward_sign_parity(l);
      
      const complex * restrict Jmp_I = &Jmm_I[mp*Nm];
      const complex * restrict Jmp_P = &Jmm_P[mp*Nm];
      
      // I, even mp > 0
    if (!isOdd(l)) {   for (m=0; m<=l; m++){
	const double Delta_mpm = Delta_mp[m];
	const complex Jmpm_I = Jmp_I[midx[m]];
	//	const complex Jmpnegm_I = Jmp_I[midx[-m]];
	const double fact0 = (Delta_mpm * Deltamp0_norml);

	
	EI[m] +=  fact0*Jmpm_I;
	//	EI[-m] += fact0*Jmpnegm_I*signnegm;
      }
    }

      // P, even mp > 0
      if (l>=2) 
	for (m=0; m<=l; m++){
	  const double Delta_mpm = Delta_mp[m];
	  const complex Jmpm_P = Jmp_P[midx[m]];
	  const complex Jmpnegm_P = Jmp_P[midx[-m]];
	  const double fact2 = (Delta_mpm * Deltamp2_norml);
	  
	  EP[m] +=  fact2*Jmpm_P;
	  EP[-m] += fact2*Jmpnegm_P*signnegm;
	}	     
    }
    

    if (0) {
      for (m=-l;m<=l;m++) {
	printf("I| % d % d: ",l,m);
	printf("Z % e % e | ",creal(ZI[m]),cimag(ZI[m]));
	printf("O % e % e | ",creal(OI[m]),cimag(OI[m]));
	printf("E % e % e |\n",creal(EI[m]),cimag(EI[m]));
      }
      printf("\n");
    }
    if (0) {
      for (m=-l;m<=l;m++) {
	printf("P| % d % d: ",l,m);
	printf("Z % e % e | ",creal(ZP[m]),cimag(ZP[m]));
	printf("O % e % e | ",creal(OP[m]),cimag(OP[m]));
	printf("E % e % e |\n",creal(EP[m]),cimag(EP[m]));
      }
      printf("\n");
    }
    
    // Collect even and odd terms
    for (m=0;m<=l;m++) {
      Tl[m] = ZI[m] + EI[m] + OI[m];
      Tl[-m] =  spinsfast_forward_sign_parity(m)*(creal(Tl[m]) - I*cimag(Tl[m]));
    }
    
    for (m=-l;m<=l;m++) {
      P2l[m] = ZP[m] + EP[m] + OP[m];
    }

    // We are now done looping over m' and m 
    //
	
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

  for (l=0;l<=lmax;l++) {
    complex * restrict Tl = &T[lm_ind(l,0,lmax)];
    //   complex * restrict El = &E[lm_ind(l,0,lmax)];
    //   complex * restrict Bl = &B[lm_ind(l,0,lmax)];
    complex * restrict P2l = &P2[lm_ind(l,0,lmax)];
    
    // Tl[0] /= 2;  // This is not a problem because we don't use [-m] terms
    P2l[0] /= 2;  // Here we add in the [-m] term, thus double counting m=0
    
    for (m=-l; m<=l; m++){
      Tl[m] *= Itom[m];  // (-i)^m * (-i)^0
      P2l[m] *= -Itom[m]; // (-i)^m * (-i)^2
    }
  }
  
  
  free(midx_helper);
  
  free(Itom_helper);
  
  free(ZI_);  free(OI_);  free(EI_); 
  free(ZP_);  free(OP_);  free(EP_); 

}
