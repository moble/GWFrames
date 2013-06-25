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
#include <complex.h>
#include <fftw3.h>

#include <alm.h>
#include <wigner_d_halfpi.h>
#include <spinsfast_forward.h>
#include <spinsfast_backward.h>

void printf_diff(complex *a, complex *a2, int lmax);



#define USAGE "%s <lmax> <Nphi> <Ntheta-1>\n",argv[0]

int main(int argc, char *argv[]) {
  
  // Read in the arguments
  if (argc!=4) {
    printf(USAGE);
    abort();
  }
  
  int iarg=1;
  int lmax = atoi(argv[iarg++]);
  int Nphi = atoi(argv[iarg++]);
  int Ntheta = 1+atoi(argv[iarg++]);


  // Sizes of various data objects
  int Npix = Nphi * Ntheta;
  int Nlm = N_lm(lmax);
  int Nm = 2*lmax+1;
  int NJmm = (lmax+1)*Nm;
  int NGmm = Nm*Nm;


  // Report the sizes
  printf("lmax = %d\n",lmax);
  printf("Nphi = %d\n",Nphi);
  printf("Ntheta = %d\n",Ntheta);
  printf("2*(Ntheta-1) = %d\n",2*(Ntheta-1));
  
  printf("Npix = %d\n", Npix);
  printf("Nm = %d\n",Nm);
  printf("Nlm = %d\n",Nlm);

  //////////////////////////////
  //
  // Allocate data objects
  //
  /////////////////////////////


  // Real space polarized function, I and P = Q+iU parts

  fftw_complex *f = calloc(2*Npix,sizeof(fftw_complex));
  fftw_complex *fI = f;
  fftw_complex *fP = &f[Npix];
 

  // Harmonic space coeffients, intensity and polarization
 
  fftw_complex *aT = calloc(Nlm,sizeof(fftw_complex));
  fftw_complex *aP = calloc(Nlm,sizeof(fftw_complex));
  
  fftw_complex *aPnew = calloc(Nlm,sizeof(fftw_complex));
  fftw_complex *aTnew = calloc(Nlm,sizeof(fftw_complex));
 
  
  // Ancillary harmonic objects, G and J

  fftw_complex *Gmm_I = fftw_malloc(NGmm*sizeof(fftw_complex));
  fftw_complex *Gmm_P = fftw_malloc(NGmm*sizeof(fftw_complex));

  fftw_complex *Jmm = fftw_malloc(2*NJmm*sizeof(fftw_complex));
  fftw_complex *Jmm_I = Jmm;
  fftw_complex *Jmm_P = &Jmm[NJmm];

  
  ////////////////////////////////
  //
  // Generate some random harmonic coefficients 
  //
  ////////////////////////////////

  int i,l,m;
  srand48(524398);

  fftw_complex tmp;

  for (l=0;l<=lmax;l++) {
    for (m=0;m<l;m++){
      tmp = drand48() + I*drand48();
      
      i = lm_ind(l, m, lmax);
      aT[i] = tmp;
      i = lm_ind(l, -m, lmax);
      aT[i] = pow(-1,m)*conj(tmp);
      
      if (m==0) 
	aT[i] = creal(tmp);
    }
  }


  for (l=2;l<=lmax;l++) {
    for (m=-l;m<l;m++){
      tmp = drand48() + I*drand48();
      
      i = lm_ind(l, m, lmax);
      aP[i] = tmp;
    }

  }
  
  printf("Sample harmonic coef. aP(l=2,m=1) = %e %e\n",creal(aP[lm_ind(2,1,lmax)]),cimag(aP[lm_ind(2,1,lmax)]));


  ////////////////////////////////////
  //
  //  Set up & execute transform:  alm -> real-space iqu  
  //
  ////////////////////////////////////


  //  Initialize the wigner Delta functions
  wdhp_TN_helper *DeltaTN = wdhp_TN_helper_init(lmax);
  
  //  Transform to Gmm (L^3 time)
  spinsfast_backward_Gmm_alm2iqu(aT,aP,lmax,Gmm_I,Gmm_P, WDHP_METHOD_TN_PLANE, (void *)DeltaTN);

  //  Transform to real space via FFTs
  spinsfast_backward_transform(fI, Ntheta, Nphi, lmax, Gmm_I);
  spinsfast_backward_transform(fP, Ntheta, Nphi, lmax, Gmm_P);
  


  ////////////////////////////////////
  //
  //  Set up & execute transform:  real-space iqu -> alm
  //
  ////////////////////////////////////

  int spins02[2] = {0,2};
  spinsfast_forward_multi_Jmm(f, spins02, 2, Ntheta, Nphi, lmax, Jmm);


 
  spinsfast_forward_transform_iqu2alm(aTnew,
				      aPnew,
				      lmax,
				      Jmm_I, 
				      Jmm_P, 
				      WDHP_METHOD_TN_PLANE,(void *)DeltaTN );

  printf("T harmonic coefficients, new - old\n");
  printf_diff(aT, aTnew, lmax);
  
  printf("P harmonic coefficients, new - old\n");
  printf_diff(aP, aPnew, lmax);




  

  fftw_free(aT);
  fftw_free(aTnew);
  
  return(0);
}

void printf_diff(complex *a, complex *a2, int lmax) {
  int i,l,m;
  complex diff = 0;
  complex maxdiff = 0;
  double fracdiff = 0;
  double fracmaxdiff = 0;
  int  maxl=0,maxm=0;
  int  fracmaxl=0,fracmaxm=0;

  double rmsdiff = 0;
  double rmsfracdiff = 0;

  for (i=0; i< N_lm(lmax); i++) {
    ind_lm(i, &l, &m, lmax);
    //    printf("% d % d % d % e % e | %e %e\n",s,l,m,creal(a2[i]),cimag(a2[i]),creal(a[i]),cimag(a[i]));
    
    diff = a2[i] - a[i];
    rmsdiff += cabs(diff)*cabs(diff);
    
    fracdiff =  cabs(diff)/cabs(a[i]);
    if (cabs(a[i])>0) rmsfracdiff += fracdiff*fracdiff;

    if ( cabs(diff) > cabs(maxdiff) ) {
      maxdiff = diff;
      maxl = l;
      maxm = m;
    }
    if (( fracdiff > fracmaxdiff )&&(cabs(a[i])>0) ){
      fracmaxdiff = fracdiff;
      fracmaxl = l;
      fracmaxm = m;
    }
  }
  
  rmsdiff /= N_lm(lmax);
  rmsdiff = sqrt(rmsdiff);

 
  rmsfracdiff /= N_lm(lmax);
  rmsfracdiff = sqrt(rmsfracdiff);
  
  printf("       max lm %6d %6d | ",maxl,maxm);
  printf("    maxdiff %+.2e %+.2e | ",creal(maxdiff),cimag(maxdiff));
  printf("    rmsdiff %+.2e\n",rmsdiff);
  printf("   fracmax lm %6d %6d | ",fracmaxl,fracmaxm);
  printf("fracmaxdiff %+.2e           | ",fracmaxdiff);
  printf("rmsfracdiff %+.2e \n",rmsfracdiff);
  
}



