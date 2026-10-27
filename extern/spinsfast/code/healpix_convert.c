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

// Routines to convert healpix (ring) grids to and from equiangular grids.
// Assumes the equiangular grid well-oversamples the healpix grid,
// or you can get empty pixels (w/ NaNs).

// From Kevin Huffenberger

#ifdef USE_HEALPIX

#include <healpix_convert.h>

void ecp2healpix_avg(complex *ecp, int Ntheta, int Nphi,
		     double *hp_r, double *hp_i, int nside) {
  // assign average of contributing pixels to healpix pixel

  int npix = 12*nside*nside;
  double *w = calloc(npix,sizeof(double));
  
  double theta,phi;
  double dtheta = M_PI/Ntheta;
  double dphi = 2*M_PI/Nphi;
  
  int i,j,p;
  long hp;
  
  // initialize 
  for (hp=0;hp<npix;hp++) {
    hp_r[hp] = 0.0;
    hp_i[hp] = 0.0;
    w[hp] = 0.0;
  }

  // accumulate pixel values & weight map
  for (j=0;j<Ntheta;j++) {
    theta = j*dtheta;
    
    for (i=0;i<Nphi;i++) {
      p = j * Nphi + i;
      
      phi = i * dphi;
      
      ang2pix_ring(nside, theta, phi, &hp);
      
      hp_r[hp] += creal(ecp[p]);
      hp_i[hp] += cimag(ecp[p]);
      w[hp] += 1.0;
    }
    
  }
  
  // normalize by weight map
  for (hp=0;hp<npix;hp++) {
    hp_r[hp] /=  w[hp];
    hp_i[hp] /=  w[hp];
  }



  
  free(w);
}


void ecp2healpix_nearest(complex *ecp, int Ntheta, int Nphi,
			 double *hp_r, double *hp_i, int nside) {

  int npix = 12*nside*nside;
  
  double theta,phi;
  double dtheta = M_PI/Ntheta;
  double dphi = 2*M_PI/Nphi;
  
  int i,j,p;
  long hp;
  
  for (hp=0;hp<npix;hp++) {
    pix2ang_ring(nside, hp, &theta, &phi);
    
    // shift by half pixel so that the *middle* of the pixel is at j*dtheta
    j = theta/dtheta + 0.5;  
    i = phi/dphi + 0.5;

    // correct for out-of-bounds
    if (j>Ntheta) {
      j = Ntheta-1;
      i += Nphi/2;
    }
    if (i>Nphi) {
      i = i % Nphi;
    }

    p = j * Nphi + i;

    hp_r[hp] = creal(ecp[p]);
    hp_i[hp] = cimag(ecp[p]);
 
  }


}


void healpix2ecp(double *hp_r, double *hp_i, int nside,
		 complex *ecp, int Ntheta, int Nphi) {
   // loop over pixels and grab healpix value 
  
  double theta,phi;
  double dtheta = M_PI/Ntheta;
  double dphi = 2*M_PI/Nphi;

  int i,j,p;
  long hp;


  for (j=0;j<Ntheta;j++) {
    theta = j*dtheta;
    
    for (i=0;i<Nphi;i++) {
      p = j * Nphi + i;
      
      phi = i * dphi;
      
      ang2pix_ring(nside, theta, phi, &hp);
      
      ecp[p] = hp_r[hp] + I*hp_i[hp];
    }
    
  }
  
}


#endif
