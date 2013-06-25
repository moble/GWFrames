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

void spinsfast_forward_Jmm(fftw_complex *f, int s, int Ntheta, int Nphi, int lmax, fftw_complex *Jmm) {
  int Nmap = 1;
  spinsfast_forward_multi_Jmm(f, &s, Nmap, Ntheta, Nphi, lmax, Jmm);
}


void spinsfast_forward_multi_Jmm(fftw_complex *f_set, int *spins, int Nmap, int Ntheta, int Nphi, int lmax, fftw_complex *Jmm_set) {

  int Nm = 2*lmax+1;
  int NImm = Nm*Nm;
  int NJmm = Nm*(lmax+1);
  fftw_complex *Imm_set = fftw_malloc(Nmap*Nm*Nm*sizeof(fftw_complex));

  spinsfast_forward_multi_Imm (f_set,spins, Nmap, Ntheta, Nphi, lmax, Imm_set);
 


  int imap;
  for (imap = 0; imap < Nmap; imap++) {
    int s = spins[imap];
    //   fftw_complex *f = &f_set[imap*Ntheta*Nphi];
    fftw_complex *Imm = &Imm_set[imap*NImm];
    fftw_complex *Jmm = &Jmm_set[imap*NJmm];


 
 
  
  int negtos =  ((s & 1) == 0) ? 1 : -1; // = (-1)^s
  
  int mp, m;


  //  Compute indexing helper
  int *midx_helper = calloc(Nm,sizeof(int));
  for (m=-lmax; m<=lmax; m++){
    midx_helper[m+lmax] = (Nm + m) % Nm;
  }
  int * restrict midx = &midx_helper[lmax];
  



  for (mp=0; mp<=lmax; mp++){
    const int mpmod = midx[mp]; 
    const int negmpmod = midx[-mp];
   

    for (m=-lmax; m<=lmax; m++){
      const int mmod = midx[m];
      int negtom =  ((m & 1) == 0) ? 1 : -1; // = (-1)^m


      if (mp==0) {
	//	if ((-10 < m) && (m < 10)) printf("%d %d %d\n",m, mmod,negtom);
	Jmm[mp*Nm + mmod] = Imm[mpmod*Nm + mmod];
      } else { 
	Jmm[mp*Nm + mmod] = Imm[mpmod*Nm + mmod] + negtom*negtos*Imm[negmpmod*Nm + mmod];
      }

    }
  }

  }

  free(Imm_set);

}
