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

#ifdef USE_FITSIO


#include <complex.h>
#include <fftw3.h>
#include <cfitsio/fitsio.h>

int dump_cimage(fftw_complex *image, int Ntheta, int Nphi, const char *filename) {
  fitsfile *fptr;
  int status = 0;
  
  fits_create_file(&fptr, filename, &status);
  
  long fpixel[2] = {1,1};
  long naxes[2] = {Nphi*2, Ntheta};

  fits_create_img( fptr, DOUBLE_IMG , 2, naxes, &status);

  double *data = calloc(2*Ntheta*Nphi,sizeof(double));
  //  int i;

  int itheta, iphi;

  for (itheta=0;itheta<Ntheta;itheta++) {
    for (iphi = 0; iphi < Nphi; iphi++) {
      data[ itheta*(2*Nphi) + iphi ] = creal( image[itheta*Nphi + iphi] );
      data[ itheta*(2*Nphi) + Nphi + iphi ] = cimag( image[itheta*Nphi + iphi] );
    }
  }

  //  for (i=0; i<Ntheta*Nphi; i++) {
  //  data[i] = creal( image[i] );
  //  data[Ntheta*Nphi + i] = cimag( image[i] );
  //}

  //  printf("data[0] = %e\n",data[0]);

  fits_write_pix (fptr, TDOUBLE, fpixel, 2*Ntheta*Nphi, data, &status);
  
  printf("%s status = %d\n",filename,status);

  fits_close_file(fptr,&status);

  free(data);
  
  return(status);
  
}


#endif
