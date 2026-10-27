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

// closely based on algorithm of Risbo, journal of geodesy (1996) 70:383-396
#include <wigner_d_halfpi.h>


void wdhp_free(wdhp *wd) {
  free(wd->d);
  free(wd->dd);
  free(wd->sqt);
  free(wd);
}


wdhp *wdhp_alloc(double jmax) {
  int i;
  wdhp *wd = calloc(1,sizeof(wdhp));
  
  // compute size of d and aux. array dd, allocate and initialize to zero
  wd->nmax = 2*(jmax+1)+1;
  
  wd->d = calloc(wd->nmax*wd->nmax,sizeof(double));
  wd->dd = calloc(wd->nmax*wd->nmax,sizeof(double));
  
  wd->sqt = calloc(wd->nmax,sizeof(double));
  // precompute square roots
  for (i=0;i<wd->nmax;i++){
    wd->sqt[i] = sqrt(i);
  }


  return(wd);
}

void wdhp_reset(wdhp *wd) {
  int i;

  for (i=0;i<wd->nmax*wd->nmax;i++){
    wd->d[i] = 0;
  }
  for (i=0;i<wd->nmax*wd->nmax;i++){
    wd->dd[i] = 0;
  }

  // initialize j=0
  wd->d[0] = 1;
  wd->j = 0;
  wd->twicej = 0;
  wd->n = wd->twicej+1;
  
}

wdhp *wdhp_init(double jmax) {
  
  wdhp *wd = wdhp_alloc(jmax);
  wdhp_reset(wd);
  
  return(wd);
}

void wdhp_jplushalf(wdhp *wd) {
  // compute d matrix for the next half step in j
  int i,k;  
  wd->twicej += 1;  
  wd->j = wd->twicej/2.0;
  
  int twicej = wd->twicej;
  int n = wd->n;
  int n1 = twicej + 1;
  wd->n = n1;
  // use n for old d, n1 for dd and new d
  
  double * restrict dd = wd->dd;
  double * restrict d = wd->d;
  double dval;
  
  //  double p = M_SQRT1_2;
  //  double q = wd->q;
  //  double pc = wd->pc;
  //  double qc = wd->qc;

  if (twicej == 1) {
    // treat j=1/2 as special case
    d[0*n1 + 0] = M_SQRT1_2;
    d[0*n1 + 1] = M_SQRT1_2;
    d[1*n1 + 0] = -M_SQRT1_2;
    d[1*n1 + 1] = M_SQRT1_2;
  } else {
    
 
    int N = (n1)*(n1);
    // Clear out old dd array
    //memset(dd,0,N*sizeof(double));
    for (i=0;i<N;i++) {
      dd[i] = 0;
    }
    
    // Compute new dd array from old d array
    //  Because of the symmetries of the d(pi/2), we only have to loop over 1/8th of the array.
    //  Below we fill in the rest from the symmetries.
    int jp1 = twicej/2+1;
    double sqrthalf_by_twicej = M_SQRT1_2/twicej;
    double dval_sqrthalf_by_twicej;

    for (i=0;i<jp1;i++) {

      double sqrt_twicejminusi = wd->sqt[wd->twicej - i];
      double sqrt_ip1 = wd->sqt[i+1];
      int klimit = i+2;//((i+2)<jp1) ? i+2 : jp1;

      double * restrict di = &d[i*n];
      double * restrict ddi = &dd[i*n1];
      double * restrict ddip1 = &dd[(i+1)*n1];


      for (k=0;k<klimit;k++) {  // limited to below the diagonal
	double sqrt_twicejminusk = wd->sqt[twicej - k];
	int kp1 = k+1;
	double sqrt_kp1 = wd->sqt[kp1];
	dval = di[k];
	//	printf("dval = %e\n",dval);
	dval_sqrthalf_by_twicej = dval * sqrthalf_by_twicej;
	//	printf("dval sqrt(0.5) / 2j = %e\n",dval_sqrthalf_by_twicej);
	//	printf("sqrt_twicejminusk = %e\n",sqrt_twicejminusk);

	ddi[k] += ( sqrt_twicejminusi *
		    sqrt_twicejminusk *
		    dval_sqrthalf_by_twicej );
	ddip1[k] -= ( sqrt_ip1 *
		      sqrt_twicejminusk *
		      dval_sqrthalf_by_twicej );
  	ddi[kp1] += ( sqrt_twicejminusi *
		      sqrt_kp1 *
		      dval_sqrthalf_by_twicej );
    	ddip1[kp1] += ( sqrt_ip1 *
			sqrt_kp1 *
			dval_sqrthalf_by_twicej );

	//	printf("dd[0] = %e\n",dd[0]);
      }
    }
    

    //    for (i=0;i<N;i++) {
    //   d[i] = 0;
    // } 
    
    int ilimit = (twicej+2)/2;
    int twicejp1 = twicej+1;

    int twicejminusk;
    int twicejminusi;
    int signi,signimk,sign2jmk;

    // Set up sign array w/ 
    //    sign[-1] = -1
    //    sign[0] = 1
    //    sign[1] = -1
    // Then (-1)^n = sign[n%2]
    // or   (-1)^n = sign[n&1]
    const int const signhelp[3] = {-1,1,-1};
    const int const *sign = &signhelp[1];

    // Copy dd array to d, using the symmetries.
    //  This took forever to figure out.
    for (i=0;i<ilimit;i++) {
      twicejminusi = twicej-i;
      signi = sign[i & 1];

      double * restrict di = &d[i*n1];
      double * restrict ddi = &dd[i*n1];
      double * restrict d2jmi = &d[(twicejminusi)*n1];
      
      for (k=0;k<=i;k++) {
	twicejminusk = twicej - k;
	dval = ddi[k];
	signimk = sign[(i-k) & 1];
	sign2jmk = sign[twicejminusk & 1];
	int kn1 = k*n1;
	
	di[k] = dval;
	d2jmi[(twicejminusk)] = signimk * dval;
	
	d[kn1+i] = signimk * dval;
	d[(twicejminusk)*n1+(twicejminusi)] = dval;
      }
      
      if (i < twicejp1/2) {
	for (k=0;k<=i;k++) {
	  twicejminusk = twicej - k;
	  dval = ddi[k];
	  //	  signimk = sign[(i-k) & 1];
	  sign2jmk = sign[twicejminusk & 1];
	  int kn1 = k*n1;
	  
 	  // upper right quadrant
 	  d[kn1+(twicejminusi)] = signi * dval; // upper left triangle
	  di[twicejminusk] = signi * dval; // lower right triangle inc. diagonal
	  
	  // lower left quadrant
	  d2jmi[k] = sign2jmk * dval; //upper left triangle
	  d[(twicejminusk)*n1 + i] = sign2jmk * dval;// lower right triangle inc diagonal
	}
      }
    }
    
  
    
  }

}



void wdhp_jplus1(wdhp *wd){
  // compute d matrix for the next whole step in j
 
  wdhp_jplushalf(wd);
  wdhp_jplushalf(wd);
}

inline double wdhp_getj(wdhp *wd) {
  return(wd->twicej/2.0);
}

inline double wdhp_get(wdhp *wd, double m1, double m2){
  // returns d^j_{m1 m2} from precomputed wd structure
  int i,k;
  
  i = (int) (wd->j + m1);
  k = (int) (wd->j + m2);
  
  //  printf("\tik = %d %d\n",i,k);
  
  return(wd->d[i*wd->n+k]);
}


inline double wdhp_integer_get(wdhp *wd, int m1, int m2) {

  int j = wd->j;

  int i = j + m1;
  int k = j + m2;
 
  return( wd->d[i*wd->n+k] );
}

inline double *wdhp_integer_getrow(wdhp *wd, int m1) {

  int j = wd->j;

  int i = j + m1;
 
  return( &(wd->d[i*wd->n+j]) );

}



int wdhp_integer_idx(int l, int m1, int m2) {
  // returns the index for a Delta matrix, when the matrix has incremented m2s.
  
  int i = l+m1;
  int k = l+m2;
  

  // loff returns the sum of (2l+1)^2
/*   int lp1 = l+1; */
   int twicel = 2*l; 
   int twicelp1 = twicel+1; 
/*   int twicelp1_lp1 = lp1*twicelp1; */

/*   int loff = twicel* twicelp1_lp1/3 + twicelp1_lp1; */


  int lm1 = l-1;
  int twicelm1 = 2*lm1;
  int twicelm1p1 = twicelm1+1;
  int twicelm1p1_l = l*twicelm1p1;

  int loff = twicelm1* twicelm1p1_l/3 + twicelm1p1_l;
  
  // printf("l=%d loff=%d\n",l,loff);

  return(loff + i*twicelp1 + k );
}

int wdhp_integer_N(int lmax) {
  // number of entries in all matrices up to lmax

  return( wdhp_integer_idx(lmax+1, -lmax-1, -lmax-1) );

}

double *wdhp_integer_precompute(int lmax) {
  
  double *Delta = calloc(wdhp_integer_N(lmax),sizeof(double));
  wdhp *wdhp_Delta = wdhp_init(lmax);
  int l,mp,m;

  for (l=0;l<=lmax;l++) {
    
    for (mp=-l;mp<=l;mp++) {
      for (m=-l;m<=l;m++) {
	Delta[wdhp_integer_idx(l,mp,m)] = wdhp_integer_get(wdhp_Delta,mp,m);
      }
    }
    
    if (l<lmax) wdhp_jplus1(wdhp_Delta);
  }

  wdhp_free(wdhp_Delta);
  
  return(Delta);
}
