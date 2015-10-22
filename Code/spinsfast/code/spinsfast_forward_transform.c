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

int spinsfast_forward_sign_parity(int m){
  // returns  (-1)^m,
  // = 1 if even
  // = -1 if odd.

  int eo = (m & 1);

  return( 1 - eo - eo );
}



void spinsfast_forward_transform(fftw_complex * restrict a, const int Ntransform, const int *spins,const int lmax, fftw_complex * restrict Jmm_set, int DeltaMethod, void *Deltawork) {
  int l,m,mp;


  int Nm = 2*lmax+1;
  double norml;

  fftw_complex *Itom_helper = fftw_malloc(Nm*sizeof(complex));
  fftw_complex *Itom = &Itom_helper[lmax];

  for (m=-lmax; m<=lmax; m++){
    Itom[m] = cpow(I,m);
  }

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


            const int twicelp1 = 2*l+1;
            norml = sqrt(twicelp1)/2./sqrt(M_PI);
            const int negtol =  spinsfast_forward_sign_parity(l); // = (-1)^l

            for (mp=0; mp<=l; mp++){
              //   int mpmod = midx[mp];
              const int negtomp =  spinsfast_forward_sign_parity(mp); // = (-1)^mp

              const int signnegm =  negtol*negtomp; // = (1-)^(l+mp)

              const double * restrict Delta_mp = NULL;

              // Grab/compute the mp row (a 1-d array) of the Wigner-d Delta matrix.
              Delta_mp = Delta_getrow( DeltaMethod, Deltawork, Deltal, l,twicelp1, mp);

              // Get Delta_{mp s} from Delta_{mp |s|}
              const int s_sign_fudge = (s>=0) ? 1 : spinsfast_forward_sign_parity(l+mp);
              const double Deltamps_norml = Delta_mp[abs(s)] * norml * s_sign_fudge;

              const complex * restrict Jmp = &Jmm[mp*Nm];

              for (m=0; m<=l; m++){
                const double Delta_mpm = Delta_mp[m];

                const int mmod = midx[m];
                const int negmmod = midx[-m];

                const double fact = (Delta_mpm * Deltamps_norml);

                const complex Jmpm = Jmp[mmod];
                const complex Jmpnegm = Jmp[negmmod];

                const complex term = fact*Jmpm*signnegm;

                const complex term_negm = fact*Jmpnegm;



                /*   if (l<5) { */
                /*     printf("% d % d % d | %.3d %.3d | % e % e | % e % e | % e % e\n",l,mp,m, mmod, negmmod, fact, fact_negm, creal(Jmpm), creal(Jmpnegm), creal(term), creal(term_negm)); */

                /*   } */

                asl[m] +=  term;
                asl[-m] += term_negm;
              }


            }

          }
        }

        // Increment Delta to next l if Risbo not precomputed
        if (l<lmax) {
          if ( (DeltaMethod==WDHP_METHOD_RISBO) ) {
            Delta_increment_l(DeltaMethod, Deltawork);
          }
        }
  }


  // Set the phases

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





void spinsfast_map2salm(fftw_complex *f, fftw_complex *alm, int s, int Ntheta, int Nphi, int lmax){
  int Ntransform = 1;
  int Nm = 2*lmax+1;
  int NJmm = (lmax+1)*Nm;

  fftw_complex *Jmm = fftw_malloc(NJmm*sizeof(fftw_complex));

  wdhp_TN_helper *DeltaTN = wdhp_TN_helper_init(lmax);

  //  Transform to Jmm via modified FFT
  spinsfast_forward_Jmm (f,s, Ntheta, Nphi, lmax, Jmm);

  // Transform Jmm to alm (L^3 time)
  spinsfast_forward_transform(alm,  Ntransform, &s, lmax, Jmm, WDHP_METHOD_TN_PLANE,(void *)DeltaTN );

  wdhp_TN_helper_free(DeltaTN);
  free(Jmm);
}
