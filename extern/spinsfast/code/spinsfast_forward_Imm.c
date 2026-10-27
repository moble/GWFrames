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

int indx2p (int ip, int wsize) {
  return( (ip > wsize/2) ?  ip-wsize : ip );
}

int p2indx (int p, int wsize) {
  return ( ( p < 0 ) ? wsize + p : p);
}

void spinsfast_forward_Imm (fftw_complex *f, int s, int Ntheta, int Nphi, int lmax, fftw_complex *Imm) {
  int Nmap = 1;

  spinsfast_forward_multi_Imm (f, &s, Nmap, Ntheta, Nphi, lmax,Imm);
}


void spinsfast_quadrature_weights(fftw_complex *W, int wsize) {
  fftw_complex *w = calloc(wsize, sizeof(fftw_complex)); // fourier space weights


  // Create weights
  int ip,p;
  int eo;
  for (ip=0; ip<wsize;ip++) {
    p = indx2p(ip,wsize);
    
    eo = abs(p % 2);
    if (p == -1) {
      w[ip] = I * M_PI/2.;
    } else if (p == 1) {
      w[ip] = - I * M_PI/2.;
    } else if ( eo == 0) {
      w[ip] = 2./(1.-p*p);
    } else {
      w[ip] = 0;
    }
    
  }
  
  fftw_plan wplan = fftw_plan_dft_1d(wsize, w, W, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(wplan);
  fftw_destroy_plan(wplan);



  /*  FILE *fp = fopen("checkwW","w"); */
  /*   for (ip=0; ip<wsize;ip++) { */
  /*     p = indx2p(ip,wsize); */
  
  /*     fprintf(fp,"%d %d %d %e %e %e %e\n", ip, p, abs(p % 2), creal(w[ip]), cimag(w[ip]),creal(W[ip]), cimag(W[ip]) ); */
  /*   } */
  /*   fclose(fp); */
  
  free(w);

}


void spinsfast_f_extend_MW(fftw_complex *f, fftw_complex *F, int s, int Ntheta, int Nphi) {
  // Make a fourier transform of f extended to the whole sphere.
  //
  // Works for pixelization where first ring is north pole and the last ring is the south pole
  // ie theta = itheta*pi/(Ntheta-1)
  //    phi   = iphi*pi/Nphi
  //
  // This extends the function using the method of McEwen & Wiaux, first taking fft in phi to get f_m(theta), then extending
  // to F_m(theta) = f_m(theta) for theta <= pi and
  //    F_m(theta) = (-1)^(m+s) f_m(2pi - theta)
  //
  // This allows the function to be extended for odd values of Nphi, but I think returns numerically identical
  // results to the old method when Nphi is even.
  //

  int itheta;
  int m,im;
  int wsize = 2*(Ntheta-1);

  fftw_complex *fm = fftw_malloc(Ntheta*Nphi*sizeof(fftw_complex));
  fftw_complex *Fm = fftw_malloc(wsize*Nphi*sizeof(fftw_complex));
  fftw_complex *W = calloc(wsize, sizeof(fftw_complex)); // for real space version of w
 
  double norm =  M_PI/Nphi/(Ntheta-1);; // = 2pi/Nphi/Ntheta_extended

  spinsfast_quadrature_weights(W, wsize);
  
  // First take fft in phi for each row of theta
  int rank = 1;
  int n = Nphi;
  int howmany = Ntheta;
  int dist = Nphi;
  int stride = 1;
  
  fftw_plan fftphiplan = fftw_plan_many_dft(rank, &n, howmany,
					    f, &n,
					    stride, dist,
					    fm, &n,
					    stride, dist,
					    FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(fftphiplan);
  fftw_destroy_plan(fftphiplan);
  
  
  
  // Now loop over theta and fill in Fm to extend theta to 2pi
  int signs = pow(-1,s);
  for (itheta = 0; itheta < Ntheta; itheta++) {
    for (im=0;im<Nphi;im++) {
      
      m = (im <= Nphi/2) ? im : (im - Nphi); // convert index of m to m value
      int signm = pow(-1,m);
      
      Fm[ itheta * Nphi + im ] = creal(W[itheta]) * fm [ itheta * Nphi + im ] * norm;;
      if ( (itheta > 0) && (itheta < Ntheta) ) {
	Fm[ (wsize - itheta) * Nphi + im] = signs*signm*creal(W[wsize - itheta]) * fm[itheta * Nphi + im] * norm;
      }
      
    }
  }
  
  // Finally take fft in phi for each column of m
  rank = 1;
  n = wsize;
  howmany = Nphi;
  dist = 1;
  stride = Nphi;
  fftw_plan fftthetaplan = fftw_plan_many_dft(rank, &n, howmany,
					      Fm, &n,
					      stride, dist,
					      F, &n,
					      stride, dist,
					      FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(fftthetaplan);
  fftw_destroy_plan(fftthetaplan);
  
  free(fm);
  free(Fm);
  free(W);

}


void spinsfast_f_extend_old(fftw_complex *f, fftw_complex *F, int s, int Ntheta, int Nphi) {
  // F is complex version of f, extended to the whole sphere
  // Works for pixelization where first ring is north pole and the last ring is the south pole
  // ie theta = itheta*pi/(Ntheta-1)
  //    phi   = iphi*pi/Nphi
  //
  // The number of theta rows in F, the extended version of f, is 2*(Ntheta-1) 
  // Extend f with F(theta > pi, phi) = (-1)^s f(2pi - theta, phi + pi)

  // finally we return the fft of F

  int itheta, iphi, opp_iphi;
  int wsize = 2*(Ntheta-1);
  int signs = pow(-1,s);

  fftw_complex *W = calloc(wsize, sizeof(fftw_complex)); // for real space version of w
 
  double norm =  M_PI/Nphi/(Ntheta-1);; // = 2pi/Nphi/Ntheta_extended

  spinsfast_quadrature_weights(W, wsize);
 

  for (itheta = 0; itheta < Ntheta; itheta++) {
    for (iphi = 0; iphi < Nphi; iphi++) {
      
      F[ itheta * Nphi + iphi ] = creal(W[itheta]) * f [ itheta * Nphi + iphi ] * norm;
      //F[ itheta * Nphi + iphi ] = 0;
      //      if ( s % 2 == 0) {
      opp_iphi = (iphi + Nphi/2) % Nphi;
      //      } else {
      //opp_iphi = iphi;
      //      }
      
      
      if ( (itheta > 0) && (itheta < Ntheta) ) {
	F[ (2*(Ntheta-1) - itheta) * Nphi + opp_iphi ] = signs*creal(W[2*(Ntheta-1) - itheta]) * f[itheta * Nphi + iphi] * norm;
	//	  F[ (2*(Ntheta-1) - itheta) * Nphi + opp_iphi ] = 0;
      }
      
    }
  }

  fftw_plan fftplan = fftw_plan_dft_2d( wsize,Nphi, F, F, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(fftplan);
  fftw_destroy_plan(fftplan);

}


void spinsfast_forward_multi_Imm (fftw_complex *f_set, int *s, int Nmap, int Ntheta, int Nphi, int lmax, fftw_complex *Imm_set) {
  // This function takes the fft of the theta extended map, then zero pads and reorganizes it.
 
  int itheta, iphi, opp_iphi;
  int m,mp,im;

  int Nm = 2*lmax + 1 ;  
  int NImm = Nm*Nm;

  int wsize = 2*(Ntheta-1);

  //  fftw_complex *fm = fftw_malloc(Ntheta*Nphi*sizeof(fftw_complex));
  //  fftw_complex *Fm = fftw_malloc(wsize*Nphi*sizeof(fftw_complex));
  fftw_complex *F = fftw_malloc(wsize*Nphi*sizeof(fftw_complex));
  //  fftw_complex *F2 = fftw_malloc(wsize*Nphi*sizeof(fftw_complex));
  //  fftw_complex *W = calloc(wsize, sizeof(fftw_complex)); // for real space version of w
  


  int imap;
  double norm =  M_PI/Nphi/(Ntheta-1);; // = 2pi/Nphi/Ntheta_extended

  for (imap = 0;imap < Nmap;imap++) {
    fftw_complex *f = &f_set[imap*Ntheta*Nphi];
    
    //    spinsfast_f_extend_old(f, F, s[imap], Ntheta, Nphi);
    spinsfast_f_extend_MW(f, F, s[imap], Ntheta, Nphi);

    
    
    // copy FT to Imm
    int limit = lmax;
    fftw_complex *Imm = &Imm_set[imap*NImm];
    for (mp=0;mp<NImm;mp++) {
      Imm[mp] = 0;
    }
    
    if ( 2*limit+1 > Nphi   ) {
      printf("Imm: Nphi pixel warning\n");
      limit = (Nphi-1)/2;
    }
    if ( 2*limit+1 > 2*(Ntheta-1) ){
      printf("Imm: Ntheta pixel warning\n");
      limit = Ntheta-3;
    }
    
    for (mp=0;mp<=limit;mp++) {
      for (m=0;m<=limit;m++) {
	
	// ++
	Imm[ mp * Nm + m ] = F[ mp * Nphi + m];
	
	// +-
	if (m > 0)
	  Imm[ mp * Nm + (Nm - m)] = F[ mp * Nphi + (Nphi - m) ];
	
	// -+
	if (mp > 0) 
	  Imm[ (Nm - mp) * Nm + m] = F[ (wsize - mp) * Nphi + m ];
	
	// --
	if ( (mp > 0) && (m > 0) )
	  Imm[ (Nm - mp) * Nm + (Nm - m)] = F[ (wsize - mp) * Nphi + (Nphi - m) ];
	
      }
    }
    
  }
  
  free(F);
}


void spinsfast_forward_multi_Imm_oldextension (fftw_complex *f_set, int *s, int Nmap, int Ntheta, int Nphi, int lmax, fftw_complex *Imm_set) {
  // F is complex version of f, extended to the whole sphere
  // Works for pixelization where first ring is north pole and the last ring is the south pole
  // ie theta = itheta*pi/(Ntheta-1)
  //    phi   = iphi*pi/Nphi
  //
  // The number of theta rows in F, the extended version of f, is 2*(Ntheta-1) 
  // Extend f with F(theta > pi, phi) = (-1)^s f(2pi - theta, phi + pi)
 
  // How is Imm indexed?

  int itheta, iphi, opp_iphi;
  int m,mp;

  int Nm = 2*lmax + 1 ;  
  int NImm = Nm*Nm;
  int wsize = 2*(Ntheta-1);

  fftw_complex *F = fftw_malloc(wsize*Nphi*sizeof(fftw_complex));
  fftw_complex *W = calloc(wsize, sizeof(fftw_complex)); // for real space version of w

  spinsfast_quadrature_weights(W, wsize);

  
  
  double norm =  M_PI/Nphi/(Ntheta-1);; // = 2pi/Nphi/Ntheta_extended
  // double norm = 1.0;
  
/*   printf("==========================================\n"); */
/*   for (mp=0;mp<8;mp++) { */
/*     for (m=0;m<8;m++) { */
/*       printf("% .2e",f[ mp * Nphi + m]); */
/*     } */
/*     printf("|  %e %e \n", creal(W[mp]),cimag(W[mp]) ); */
/*   } */

  fftw_plan fftplan = fftw_plan_dft_2d( wsize,Nphi, F, F, FFTW_FORWARD, FFTW_ESTIMATE);
  
  int imap;
  for (imap = 0;imap < Nmap;imap++) {
    fftw_complex *f = &f_set[imap*Ntheta*Nphi];

    int signs = pow(-1,s[imap]);
    
    // copy f to F
    for (itheta = 0; itheta < Ntheta; itheta++) {
      for (iphi = 0; iphi < Nphi; iphi++) {
	
	F[ itheta * Nphi + iphi ] = creal(W[itheta]) * f [ itheta * Nphi + iphi ] * norm;
	
	//      if ( s % 2 == 0) {
	opp_iphi = (iphi + Nphi/2) % Nphi;
	//      } else {
    	//opp_iphi = iphi;
	//      }      
	
	
	if ( (itheta > 0) && (itheta < Ntheta) ) {
	  F[ (2*(Ntheta-1) - itheta) * Nphi + opp_iphi ] = signs*creal(W[2*(Ntheta-1) - itheta]) * f[itheta * Nphi + iphi] * norm;
	}
	
      }
    }
    
    
    /*  fp = fopen("checkF","w"); */
    /*   for (mp=0;mp<2*(Ntheta-1);mp++) { */
    /*     //  printf("%.2d ",mp); */
    /*     for (m=0;m<Nphi;m++) { */
    /*       fprintf(fp,"%d %d % .2e % .2e \n",mp,m,creal(F[ mp * Nphi + m]),cimag(F[ mp * Nphi + m])); */
    /*     } */
    /*     fprintf(fp,"\n\n"); */
    /*   } */
    /*   fclose(fp); */
    
    fftw_execute(fftplan);
    
    
    /*  printf("==========================================\n"); */
    /*   for (mp=0;mp<5;mp++) { */
    /*     for (m=0;m<5;m++) { */
    /*        printf("% e % e | ",creal(F[ mp * Nphi + m]),cimag(F[ mp * Nphi + m])); */
    /*     } */
    /*     printf("\n"); */
    /*   } */
    
    
    
    
    // copy FT to Imm
    int limit = lmax;
    fftw_complex *Imm = &Imm_set[imap*NImm];
    for (mp=0;mp<NImm;mp++) {
      Imm[mp] = 0;
    }
    
    if ( 2*limit+1 > Nphi   ) {
      printf("Imm: Nphi pixel warning\n");
      limit = (Nphi-1)/2;
    }
    if ( 2*limit+1 > 2*(Ntheta-1) ){
      printf("Imm: Ntheta pixel warning\n");
      limit = Ntheta-3;
    }
    
    for (mp=0;mp<=limit;mp++) {
      for (m=0;m<=limit;m++) {
	
	// ++
	Imm[ mp * Nm + m ] = F[ mp * Nphi + m];
	
	// +-
	if (m > 0)
	  Imm[ mp * Nm + (Nm - m)] = F[ mp * Nphi + (Nphi - m) ];
	
	// -+
	if (mp > 0) 
	  Imm[ (Nm - mp) * Nm + m] = F[ (wsize - mp) * Nphi + m ];
	
	// --
	if ( (mp > 0) && (m > 0) )
	  Imm[ (Nm - mp) * Nm + (Nm - m)] = F[ (wsize - mp) * Nphi + (Nphi - m) ];
	
      }
    }
    
  }
  
  fftw_destroy_plan(fftplan);
  
  
  free(W);
  free(F);
}
