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

#include <wigner_d_halfpi_TN.h>




//inline
int wdhp_lmind_pos(int l,int m){
  return((l*(l+1))/2 + m);
}


//inline 
int wdhp_sign_parity(int m){
  // returns  (-1)^m, 
  // = 1 if even
  // = -1 if odd.

  int eo = (m & 1);

  return( 1 - eo - eo );
}


FPTYPE wdhp_get_ll0(int l) {
  
  FPTYPE D = 1.0;
  int lp;

  for (lp = 1;lp<=l; lp++) {
    FPTYPE fact = -sqrtl( 1.0 - 1.0/(2*lp) );
    D *= fact;

    //  printf("fact = %e D = %e\n",fact,D);
  }

  return(D);

}



FPTYPE wdhp_get_llm(int l,int m) {

  FPTYPE D =  wdhp_get_ll0(l-m); // Returns Delta^(l-m)_{(l-m) 0}

  int lp, mp;

  // printf("%d %d 0 : %e\n", l-m,l-m,D);

  for (mp=1;mp<=m;mp++) {
    lp = l-m+mp;


    FPTYPE ool = 1.0/lp;
    FPTYPE mol = mp/(FPTYPE)lp;

    FPTYPE fact = sqrtl( (1-0.5*ool) / ( (1+mol)*(1+mol-ool) ));

    //  FPTYPE fact = sqrtl( (lp/2.0)*(2*lp-1)/(FPTYPE)((lp+mp)*(lp+mp-1)) );

    //     printf("%d %d | %d %d | ",l,m,lp,mp);  printf("ool = %e mol = %e | ",(double)ool,(double)mol);       printf("fact = %e\n",(double)fact);

    D *=fact;
  }

  return(D);
}


void wdhp_get_all_llm(int lmax,FPTYPE *D_all_llm) {
  int l, m, i;
  
  for (l=0;l<=lmax;l++) {
    for (m=0;m<=l;m++){
      i = wdhp_lmind_pos(l,m);
      
      D_all_llm[i] = wdhp_get_llm(l,m);
    }
  }

}


void wdhp_get_all_llm2(int lmax,FPTYPE *D_all_llm, FPTYPE *sqt, FPTYPE *invsqt) {

  int l, ll0;
  

  //////
  //
  // Fill in all Delta_^l_l0
  //
  ///////


  FPTYPE D = 1.0;
  ll0 = wdhp_lmind_pos(0,0);
  D_all_llm[ll0] = D;
  
  for (l=1;l<=lmax;l++) {
    ll0 = wdhp_lmind_pos(l,0);
    
    int twol = 2*l;
    FPTYPE fact = -sqt[twol-1]*invsqt[twol];
    D *= fact;
    
    D_all_llm[ll0] = D;
  }

  ///////
  //
  //  Trace down diagonals to fill out llm
  //
  ////////
  int lp;
  int m,llm;

  for (lp=0;lp<=lmax;lp++) {
    ll0 = wdhp_lmind_pos(lp,0);
    
    D = D_all_llm[ll0];
    
    for (m=1;lp+m<=lmax;m++) {
      l = lp + m;
      
  
      int lplusm = l+m;
      int n1 = l*(2*l-1);
      int n2 = lplusm*(lplusm-1);
      
      FPTYPE fact = sqt[n1]*invsqt[n2]*M_SQRT1_2l;
      //     printf("%d %d %d %d %e\n",l,m,n1,n2,(double)fact);

      D *= fact;
      
      llm = wdhp_lmind_pos(l,m);
      D_all_llm[llm] = D;
    }
    
  }


}


FPTYPE *wdhp_init_sqt(int lmax) {

  int N = 4*(lmax+1)*(lmax+1);

  FPTYPE *sqt = calloc(N+1,sizeof(FPTYPE));
  
  int i;

  for (i=0;i<=N;i++) {
    sqt[i] = sqrtl((FPTYPE)i);
  }

  return(sqt);
}

FPTYPE *wdhp_init_invsqt(int lmax, FPTYPE *sqt) {

  int N = 4*(lmax+1)*(lmax+1);

  FPTYPE *invsqt = calloc(N+1,sizeof(FPTYPE));
  
  int i;

  for (i=0;i<=N;i++) {
    invsqt[i] = 1.0/sqt[i];
  }

  return(invsqt);
}



//inline 
FPTYPE wdhp_rowrecurs_coef1(int l,int m1,int m2, FPTYPE *sqt, FPTYPE *invsqt){
  int n = (l-m1)*(l+m1+1);

  // FPTYPE invs = ;

  return( 2*m2*invsqt[n] );
}

//inline 
FPTYPE wdhp_rowrecurs_coef2(int l,int m1,int m2, FPTYPE *sqt, FPTYPE *invsqt){

  int lminusm1 = l-m1;
  int lplusm1 = l+m1;

  int n1 = (lminusm1 - 1)*( lplusm1 + 2);
  int n2 = lminusm1 *(lplusm1 + 1);

  return( sqt[n1]*invsqt[n2]);

}


//inline 
FPTYPE wdhp_rowrecurs(FPTYPE D1, FPTYPE D2, int l,int m1,int m2, FPTYPE *sqt, FPTYPE *invsqt){
  // D1 = D^l_{(m1+1)(m2)}
  // D2 = D^l_{(m1+2)(m2)}
  
  FPTYPE c1 = wdhp_rowrecurs_coef1(l,m1,m2,sqt,invsqt);
  FPTYPE c2 = wdhp_rowrecurs_coef2(l,m1,m2,sqt,invsqt);
  
  //  printf("m1 m2 %d %d | D1 D2 = % e % e | c1 c2 % e % e\n",m1,m2,(double)D1,(double)D2,(double)c1,(double)c2);

  return( c1*D1-c2*D2 );
}


FPTYPE wdhp_get_lm1m2_pos(int l,int m1,int m2, FPTYPE *sqt, FPTYPE *invsqt) {

  FPTYPE D, D1, D2;

  D1 =  wdhp_get_llm(l,m2);
  D2 = 0;
  
  int mp;

  D = D1;
  for (mp = l-1; mp>=m1; mp--) {

    D =  wdhp_rowrecurs(D1,D2,l,mp,m2, sqt,invsqt);

    D2 = D1;
    D1 = D;
  }

  return(D);
}


void wdhp_get_col_pos(int l,int m2, FPTYPE *sqt, FPTYPE *invsqt, FPTYPE *D_all_llm, double *Dcol) {
  // 
  //  Upon completion Dcol[m1] = D^l_{m1 m2}
  //


  FPTYPE D1, D2;

  D1 = D_all_llm[wdhp_lmind_pos(l,m2)];
  D2 = 0;

  int m1 = l;

  Dcol[l] = D1;
  for (m1 = l-1; m1>=0; m1--) {

    Dcol[m1] =  wdhp_rowrecurs(D1,D2,l,m1,m2, sqt,invsqt);

    D2 = D1;
    D1 = Dcol[m1];
  }

}





void wdhp_get_row_pos(int l,int m1, FPTYPE *sqt, FPTYPE *invsqt, FPTYPE *D_all_llm, double *Drow) {
  //
  // Upon completion Drow[m2] = D^l_{m1 m2}
  // 
  
  FPTYPE D,D1, D2;

  D1 = D_all_llm[wdhp_lmind_pos(l,m1)];
  D2 = 0;

  int m2 = l;

  // D contains the column D^l {... m1}
  // Drow contains the transposed row, D^l {m1 ...}
  D = D1;
  Drow[m2] = wdhp_sign_parity(m1+m2)*D;
  
  for (m2 = l-1; m2>=0; m2--) {
    
    D =  wdhp_rowrecurs(D1,D2,l,m2,m1, sqt,invsqt);
    Drow[m2] = wdhp_sign_parity(m1+m2)*D;
    
    D2 = D1;
    D1 = D;
  }

}

void wdhp_get_quarter_plane2(int l, FPTYPE *sqt, FPTYPE *invsqt, FPTYPE *D_all_llm, double *D) {
  // 
  //  Upon completion D[m1*(l+1)+m2] = D^l_{m1 m2}
  //  for m1,m2 >= 0.
  
  FPTYPE D0,D1, D2;

  int m1, Nm=l+1;
  for (m1 = l-1; m1>=0; m1--) {
    D1 = D_all_llm[wdhp_lmind_pos(l,m1)];
    D2 = 0;
    
    int m2 = l;

    // D contains the column D^l {... m1}
    // Drow contains the transposed row, D^l {m1 ...}
    D0 = D1;
    D[m1*Nm+m2] = wdhp_sign_parity(m1+m2)*D1;
    D[m2*Nm+m1] = D1;
    
    for (m2 = l-1; m2>=m1; m2--) {
      
      D0 =  wdhp_rowrecurs(D1,D2,l,m2,m1, sqt,invsqt);
      D[m1*Nm+m2] = wdhp_sign_parity(m1+m2)*D0;
      D[m2*Nm+m1] = D0;

      D2 = D1;
      D1 = D0;
    }
    
  }
}



void wdhp_get_quarter_plane(int l, FPTYPE *sqt, FPTYPE *invsqt, FPTYPE *D_all_llm, double *D) {
  // 
  //  Upon completion D[m1*(l+1)+m2] = D^l_{m1 m2}
  //  for m1,m2 >= 0.
    
  FPTYPE D0,D1, D2;
  int m2;
  int Nm = l + 1;
  
  for (m2 = 0; m2 <= l; m2++) {
    D1 = D_all_llm[wdhp_lmind_pos(l,m2)];
    D2 = 0;
    
    int m1 = l;
    
    D[m1*Nm+m2] = D1;
    D[m2*Nm+m1] = wdhp_sign_parity(m1+m2)*D1;

    for (m1 = l-1; m1>=m2; m1--) {
      
      D0 = wdhp_rowrecurs(D1,D2,l,m1,m2, sqt,invsqt);
      D[m1*Nm+m2] = D0;
      D[m2*Nm+m1] = wdhp_sign_parity(m1+m2)*D0;


      D2 = D1;
      D1 = D0;
    }
    
  }
}




FPTYPE wdhp_get_lm1m2_pos2(int l,int m1,int m2, FPTYPE *sqt, FPTYPE *invsqt) {

  int bigm,lilm,sign = 1;

  if (m1>m2) {
    lilm = m2;
    bigm = m1;
    sign = wdhp_sign_parity(m1+m2);
    //  printf("flippers\n");
  } else {
    lilm = m1;
    bigm = m2;
  }

  return( sign*wdhp_get_lm1m2_pos(l,lilm,bigm,sqt,invsqt) );
}



wdhp_TN_helper *wdhp_TN_helper_init(int lmax){
  wdhp_TN_helper *w = calloc(1,sizeof(wdhp_TN_helper));

  w->lmax = lmax;
  w->sqt = wdhp_init_sqt(lmax);
  w->invsqt = wdhp_init_invsqt(lmax,w->sqt);
  w->D_all_llm = calloc((lmax+1)*(lmax+2)/2,sizeof(FPTYPE));
  
  //wdhp_get_all_llm(lmax,w->D_all_llm);
  wdhp_get_all_llm2(lmax,w->D_all_llm,w->sqt,w->invsqt);

  w->Dwork = calloc((lmax+1)*(lmax+1),sizeof(FPTYPE));

  return(w);
}



void wdhp_TN_helper_free(wdhp_TN_helper *w) {
  free( w->Dwork );
  free( w->D_all_llm );
  free( w->invsqt );
  free( w->sqt );
  free( w );
}
