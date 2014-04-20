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

void spinsfast_backward_transform(fftw_complex  *f, int Ntheta, int Nphi, int lmax, fftw_complex *Gmm) {
  // f is Ntheta by Nphi, access f[ itheta * Nphi + iphi ].
  // F is  version of f extended to the whole sphere
  // Works for pixelization where first ring is north pole and the last ring is the south pole
  // ie theta = itheta*pi/(Ntheta-1)
  //    phi   = iphi*pi/Nphi
  //
  // The number of theta rows in the extended version of F is 2*(Ntheta-1) 
  int itheta, iphi;
  int m,mp;
  fftw_complex *F = fftw_malloc(2*(Ntheta-1)*Nphi*sizeof(fftw_complex));
  
  int NF = 2*(Ntheta-1)*Nphi;
  for (m=0;m<NF;m++) 
    F[m] = 0;
  
  // copy Imm values into F;
  
  int Nm = 2*lmax + 1 ;  
  int limit = lmax;
  
  if ( 2*limit+1 > Nphi   ) {
    printf("backtrans Nphi warning\n");
    limit = (Nphi-1)/2;
  }
  if ( 2*limit+1 > 2*(Ntheta-1) ) {
    printf("backtrans Ntheta warning\n");
    limit = Ntheta-3;
  }  


  for (mp=0;mp<=limit;mp++) {
    for (m=0;m<=limit;m++) {
      
      // ++
      F[ mp * Nphi + m] = Gmm[ mp * Nm + m ] ;
      
      // +-
      if (m > 0)
	F[ mp * Nphi + (Nphi - m) ] =	Gmm[ mp * Nm + (Nm - m)];
      
      // -+
      if (mp > 0) 
	F[ (2*(Ntheta-1) - mp) * Nphi + m ] = Gmm[ (Nm - mp) * Nm + m];
      
      // --
      if ( (mp > 0) && (m > 0) )
	F[ (2*(Ntheta-1) - mp) * Nphi + (Nphi - m) ] = Gmm[ (Nm - mp) * Nm + (Nm - m)];

    }
  }
  
  //    printf("Ntheta Nphi = %d %d\n",Ntheta,Nphi);
  //  printf ("dumping F status = %d\n", dump_cimage(F,2*(Ntheta-1),Nphi,"!output/eff.fits"));


  fftw_plan fftplan = fftw_plan_dft_2d( 2*(Ntheta-1),Nphi, F, F, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(fftplan);
  fftw_destroy_plan(fftplan);

  
  //printf ("dumping F status = %d\n", dump_cimage(F,2*(Ntheta-1),Nphi,"!output/EFF.fits"));

  
  for (itheta = 0; itheta < Ntheta; itheta++) {
    for (iphi = 0; iphi < Nphi; iphi++) {
      
      f[ itheta * Nphi + iphi ] =  F [ itheta * Nphi + iphi ];
      
    }
  }
  
  fftw_free(F);

}



void spinsfast_salm2map(fftw_complex *alm, fftw_complex *f, int s, int Ntheta, int Nphi, int lmax){
  int Nm = 2*lmax+1;
  int NGmm = Nm*Nm;
  
  wdhp_TN_helper *DeltaTN = wdhp_TN_helper_init(lmax);
  fftw_complex *Gmm = fftw_malloc(NGmm*sizeof(fftw_complex));
  
  int Ntransform = 1;
  
  spinsfast_backward_Gmm(alm,  Ntransform, &s, lmax, Gmm, WDHP_METHOD_TN_PLANE, (void *)DeltaTN); 
  spinsfast_backward_transform(f, Ntheta, Nphi, lmax, Gmm);
  
  wdhp_TN_helper_free(DeltaTN);
  free(Gmm);
  
}
