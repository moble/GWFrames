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

#include <stdlib.h>
#include <stdio.h>

#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include <alm.h>
#include <wigner_d_halfpi.h>
#include <spinsfast_forward.h>
#include <spinsfast_backward.h>

void printf_diff(complex *a, complex *a2, int lmax);



#define USAGE "%s <spin s> <lmax> <Nphi> <Ntheta-1>\n",argv[0]

int main(int argc, char *argv[]) {

  if (argc!=5) {
    printf(USAGE);
    abort();
  }
  
  int iarg=1;
  int s = atoi(argv[iarg++]);
  int lmax = atoi(argv[iarg++]);
  int Nphi = atoi(argv[iarg++]);
  int Ntheta = 1+atoi(argv[iarg++]);


  // Sizes of various data objects
  int Npix = Nphi * Ntheta;
  int Nlm = N_lm(lmax);

  
  // Report the sizes
  printf("lmax = %d\n",lmax);
  printf("Nphi = %d\n",Nphi);
  printf("Ntheta = %d\n",Ntheta);
  printf("2*(Ntheta-1) = %d\n",2*(Ntheta-1));
  
  printf("Npix = %d\n", Npix);
  printf("Nm = %d\n",2*lmax+1);
  printf("Nlm = %d\n",Nlm);

  //////////////////////////////
  //
  // Allocate data objects
  //
  /////////////////////////////

  // Real space function

  fftw_complex *f = calloc(Nphi*Ntheta,sizeof(fftw_complex));
 

  // Harmonic space coeffients
  
  fftw_complex *alm = calloc(Nlm,sizeof(fftw_complex));
  fftw_complex *alm2 = calloc(Nlm,sizeof(fftw_complex));

  ////////////////////////////////
  //
  // Generate some random harmonic coefficients 
  //
  ////////////////////////////////
  
  int i,l,m;
  srand48(524398);

  for (i=0; i< Nlm; i++) {
    ind_lm(i, &l, &m, lmax);
    if ( abs(s) <= l ) {
      alm[i] = drand48() + I*drand48();
    }
  }
  
  ////////////////////////////////////
  //
  //  Set up & execute backward transform: harmonic -> real space
  //
  ////////////////////////////////////

  printf("Computing s = %d transform\n",s);
  spinsfast_salm2map(alm, f, s, Ntheta, Nphi, lmax);

  ////////////////////////////////////
  //
  //  Set up & execute foreward transform: real -> harmonic space
  //
  ////////////////////////////////////
  
  spinsfast_map2salm(f, alm2, s, Ntheta, Nphi, lmax);

  
  ////////////////////////////////////
  //
  //  Examine difference in the coefficients from backward/forward transform pair 
  //
  ////////////////////////////////////
  
  printf_diff(alm, alm2, lmax);

  // clean up

  free(alm);
  free(alm2);
  free(f);

 
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
