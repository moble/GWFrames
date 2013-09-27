// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

// Note: This header should only be included from PNWaveforms.cpp, and
// should therefore be invisible everywhere else.  You should not need
// to include this header directly.

#ifndef DOXYGEN

// Return dot product of two three-vectors
inline double dot(const double* a, const double* b) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
inline void cross(const double* a, const double* b, double* c) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
  return;
}

const GWFrames::Quaternion xHat(0,1.0,0,0);
const GWFrames::Quaternion yHat(0,0,1.0,0);
const GWFrames::Quaternion zHat(0,0,0,1.0);

// This will do the actual work of computing the right-hand side of the evolution equations.
class TaylorT1 {
private:
  const double delta;
  const double nu, pownu2, pownu3;  // various powers of nu are needed several times, and may as well be stored
  double v, powv2, powv3, powv5, powv6, powv10;  // various powers of v are needed several times, and may as well be stored
  GWFrames::Quaternion Rax, nHat, lambdaHat;
  const double *chi1, *chi2, *OmegaHat_orb; // Components in the static Cartesian basis
  double chi1n, chi2n, chi1la, chi2la, chi1l, chi2l, chi1chi1, chi1chi2, chi2chi2; // The remaining dot products
private:
  vector<double> Omega_spin(const unsigned int HoleIndex) const {
    // Eq. (4.5) of Bohé et al. (2012) <http://arxiv.org/abs/1212.5520v1>
    const double d = (HoleIndex==1 ? delta : -delta);
    const double Omega_s = powv5 * (0.75*(1-d) + 0.5*nu
				    + powv2*(0.5625*(1-d)+0.625*nu*(2+d)-0.041666666666666667*pownu2
					     + powv2*(0.84375 + d*(-0.84375 + (4.875 - 0.15625*nu)*nu) + nu*(0.1875 + (-3.28125 - 0.020833333333333332*nu)*nu))));
    vector<double> Omega_sVec(3);
    Omega_sVec[0] = Omega_s * OmegaHat_orb[0];
    Omega_sVec[1] = Omega_s * OmegaHat_orb[1];
    Omega_sVec[2] = Omega_s * OmegaHat_orb[2];
    return Omega_sVec;
  }
  double Omega_precMag() const {
    // Eqs. (4.3) and (4.4), and above Eq. (4.1), of Bohé et al.
    // (2012) <http://arxiv.org/abs/1212.5520v1>.  Note that this
    // vector is completely along nHat (by definition).
    const double gamma = 5.51146384479718e-6*pow(v, 2)*(pow(v, 2)*(-60480.0*nu + v*(15120.0*chi1l*(-6.0*delta*(delta + 1.0) + 5.0*delta*(delta + 2.0) + 5.0) + 15120.0*chi2l*(-6.0*delta*(delta - 1.0) + 5.0*delta*(delta - 2.0) + 5.0) + v*(-982800.0*nu + v*(10080.0*chi1l*(-18.0*delta*(delta + 1.0) + delta*(delta*(4.0*nu + 15.0) + 8.0*nu + 30.0) + 4.0*nu + 15.0) + 10080.0*chi2l*(-18.0*delta*(delta - 1.0) - delta*(-delta*(4.0*nu + 15.0) + 8.0*nu + 30.0) + 4.0*nu + 15.0) + v*(-nu*(-560.0*nu*(4.0*nu + 2061.0) + 541013.822520207) + 3780.0*v*(chi1l*(4.0*delta*(delta*(nu*(16.0*nu + 61.0) - 18.0) + nu*(16.0*nu + 61.0) - 18.0) + delta*(-delta*(nu*(72.0*nu + 127.0) - 60.0) - 2.0*nu*(72.0*nu + 127.0) + 120.0) - nu*(72.0*nu + 127.0) + 60.0) + chi2l*(4.0*delta*(delta*(nu*(16.0*nu + 61.0) - 18.0) - nu*(16.0*nu + 61.0) + 18.0) + delta*(-delta*(nu*(72.0*nu + 127.0) - 60.0) + 2.0*nu*(72.0*nu + 127.0) - 120.0) - nu*(72.0*nu + 127.0) + 60.0)) + 181440.0)) + 181440.0)) + 181440.0) + 181440.0);
    const double a_l_over_vcubed = 0.00694444444444444*pow(v, 4)*(36.0*chi1n*(-6.0*delta*(delta + 1.0) + 7.0*delta*(delta + 2.0) + 7.0) + 36.0*chi2n*(-6.0*delta*(delta - 1.0) + 7.0*delta*(delta - 2.0) + 7.0) + pow(v, 2)*(-12.0*chi1n*(-9.0*delta*(delta*(3.0*nu + 4.0) + 3.0*nu + 4.0) + delta*(delta*(29.0*nu + 30.0) + 58.0*nu + 60.0) + 29.0*nu + 30.0) - 12.0*chi2n*(9.0*delta*(-delta*(3.0*nu + 4.0) + 3.0*nu + 4.0) - delta*(-delta*(29.0*nu + 30.0) + 58.0*nu + 60.0) + 29.0*nu + 30.0) + pow(v, 2)*(chi1n*(-3.0*delta*(delta*(nu*(68.0*nu + 219.0) + 36.0) + nu*(68.0*nu + 219.0) + 36.0) + delta*(delta*(nu*(208.0*nu + 531.0) + 54.0) + 2.0*nu*(208.0*nu + 531.0) + 108.0) + nu*(208.0*nu + 531.0) + 54.0) + chi2n*(3.0*delta*(-delta*(nu*(68.0*nu + 219.0) + 36.0) + nu*(68.0*nu + 219.0) + 36.0) - delta*(-delta*(nu*(208.0*nu + 531.0) + 54.0) + 2.0*nu*(208.0*nu + 531.0) + 108.0) + nu*(208.0*nu + 531.0) + 54.0))));
    return gamma*a_l_over_vcubed;
  }
  double Flux() const {
    // Eqs. (C7) -- (C13) of <http://arxiv.org/abs/0810.5336v3> [NOTE VERSION NUMBER!!!]
    const double F_2 = 0.002976190476190476*(-1247. - 980.*nu);
    const double F_3 = 12.56637061435917 + chi2l*(1.375*(-1. + delta) + 1.5*nu) + chi1l*(-1.375*(1. + delta) + 1.5*nu);
    const double F_4 = 1.494791666666667*chi1l*chi1l*(1. + delta - 2.*nu) + 6.020833333333333*chi1l*chi2l*nu - 1.494791666666667*chi2l*chi2l*(-1. + delta + 2.*nu)
      + 0.005208333333333333*(-89.*chi1chi1*(1. + delta - 2.*nu) - 412.*chi1chi2*nu + 89.*chi2chi2*(-1. + delta + 2.*nu))
      + 0.0001102292768959436*(-44711. + 18.*nu*(9271. + 1820.*nu));
    const double F_5 = -0.004674989067841954*(8191. + 16324.*nu) + 0.003472222222222222*chi2l*(531.*(-1. + delta) + 4.*(908. - 701.*delta)*nu - 2512.*pownu2)
      + 0.003472222222222222*chi1l*(-531.*(1. + delta) + 4.*(908. + 701.*delta)*nu - 2512.*pownu2);
    const double F_6 = 138.3349161635852 + 0.2617993877991494*chi2l*(-65. + 65.*delta + 68.*nu) + 0.2617993877991494*chi1l*(-65.*(1. + delta) + 68.*nu)
      - 0.00001837154614932393*nu*(482923.6129757514 + 6.*nu*(283209. + 21700.*nu)) - 8.152380952380952*std::log(16.*powv2);
    const double F_7 = 0.001298608074400543*(-78168. + nu*(300643. + 154708.*nu));
    return
      6.4 * pownu2 * powv10 * (1
			       + powv2 * ( F_2
					   + v * ( F_3
						   + v * ( F_4
							   + v * ( F_5
								   + v * ( F_6
									   + v * ( F_7 )
									   )
								   )
							   )
						   )
					   )
			       );
  }
  double dMdt() const {
    // Alvi (2001) PRD 64, 104020 <http://link.aps.org/doi/10.1103/PhysRevD.64.104020>
    return
      6.4 * pownu2 * powv10 * ( powv5
				* (0.125*chi2l*(-1.*(1. + 3.*chi2chi2)*(1. - 1.*delta + (-3. + delta)*nu) + 3.*chi1chi1*(-1. + delta + 3.*nu - 1.*delta*nu)
						+ 3.*chi1chi2*(1. + delta - 1.*(3. + delta)*nu))
				   + 0.125*chi1l*(-3.*chi1chi1*(1. + delta) + 3.*chi1chi1*(3. + delta)*nu + 3.*chi1chi2*(1. - 1.*delta + (-3. + delta)*nu)
						  + (1. + 3.*chi2chi2)*(-1. - 1.*delta + (3. + delta)*nu))) );
  }
  double dEdv() const {
    // Non-spin and spin-orbit terms come from Eq. (4.6) of Bohé et
    // al. (2012) <http://arxiv.org/abs/1212.5520v1>.  Spin-spin terms
    // come from Eqs. (C1) -- (C6) of
    // <http://arxiv.org/abs/0810.5336v3> [NOTE VERSION NUMBER!!!].
    const double dEdv_2 = 0.1666666666666667*(-9. - 1.*nu);
    const double dEdv_3 = 1.666666666666667*chi1l*(2. + 2.*delta - 1.*nu) - 1.666666666666667*chi2l*(-2. + 2.*delta + nu);
    const double dEdv_4 = -2.25*chi1l*chi1l*(1. + delta - 2.*nu) - 9.*chi1l*chi2l*nu + 2.25*chi2l*chi2l*(-1. + delta + 2.*nu)
      + 0.125*(-81. + 6.*chi1chi1*(1. + delta - 2.*nu) + (57. + 24.*chi1chi2 - 1.*nu)*nu - 6.*chi2chi2*(-1. + delta + 2.*nu));
    const double dEdv_5 = 0.1944444444444444*chi1l*(72. + delta*(72. - 31.*nu) + nu*(-121. + 2.*nu))
      + 0.1944444444444444*chi2l*(72. + nu*(-121. + 2.*nu) + delta*(-72. + 31.*nu));
    const double dEdv_6 = -0.003858024691358025*(10935. + nu*(-40149.69585598816 + nu*(1674. + 7.*nu)));
    const double dEdv_7 = -0.1875*chi2l*(-324. + nu*(1119. - 2.*nu*(172. + nu)) + delta*(324. + nu*(-633. + 14.*nu)))
      + 0.1875*chi1l*(324. + nu*(-1119. + 2.*nu*(172. + nu)) + delta*(324. + nu*(-633. + 14.*nu)));
    return
      - nu * v * (1
		  + powv2 * ( dEdv_2
			      + v * ( dEdv_3
				      + v * ( dEdv_4
					      + v * ( dEdv_5
						      + v * ( dEdv_6
							      + v * ( dEdv_7 )
							      )
						      )
					      )
				      )
			      )
		  );
  }
public:
  TaylorT1(const double idelta, const vector<double>& chi1_0, const vector<double>& chi2_0, const double v_0,
	   const GWFrames::Quaternion R_0=GWFrames::Quaternion(1,0,0,0)) :
    delta(idelta),
    nu((1.0-delta*delta)/4.0), pownu2(nu*nu), pownu3(pownu2*nu),
    v(v_0), powv2(v*v), powv3(v*powv2), powv5(powv2*powv3), powv6(powv3*powv3), powv10(powv5*powv5),
    Rax(R_0), nHat(R_0 * xHat * R_0.conjugate()), lambdaHat(R_0 * yHat * R_0.conjugate()),
    // chi1, chi2, and OmegaHat_orb will point to things in the evolution system
    chi1n(Quaternion(chi1_0).dot(nHat)), chi2n(Quaternion(chi2_0).dot(nHat)), 
    chi1la(Quaternion(chi1_0).dot(lambdaHat)), chi2la(Quaternion(chi2_0).dot(lambdaHat)),
    chi1l(Quaternion(chi1_0).dot(R_0*zHat*R_0.conjugate())), chi2l(Quaternion(chi2_0).dot(R_0*zHat*R_0.conjugate())), 
    // chi1n(chi1_0[0]), chi2n(chi2_0[0]), chi1la(chi1_0[1]), chi2la(chi2_0[1]), chi1l(chi1_0[2]), chi2l(chi2_0[2]),
    chi1chi1(dot(&chi1_0[0], &chi1_0[0])), chi1chi2(dot(&chi1_0[0], &chi2_0[0])), chi2chi2(dot(&chi2_0[0], &chi2_0[0]))
  { }
  ~TaylorT1() { }
  void RecalculateValues(double t, const double* y) {
    // Compute the various stored quantities, for use in the important functions
    const double& gamma = y[11];
    const double& Phi = y[1];
    v = y[0];
    powv2 = v*v;
    powv3 = v*powv2;
    powv5 = powv2*powv3;
    powv6 = powv3*powv3;
    powv10 = powv5*powv5;
    chi1 = &y[2];
    chi2 = &y[5];
    OmegaHat_orb = &y[8];
    const GWFrames::Quaternion chi1Q(0., chi1[0], chi1[1], chi1[2]);
    const GWFrames::Quaternion chi2Q(0., chi2[0], chi2[1], chi2[2]);
    const GWFrames::Quaternion OmegaHat_orbQ
      = GWFrames::normalized(GWFrames::Quaternion(0., OmegaHat_orb[0], OmegaHat_orb[1], OmegaHat_orb[2]));
    Rax = GWFrames::sqrtOfRotor(-OmegaHat_orbQ*zHat);
    const GWFrames::Quaternion R = Rax * GWFrames::exp(((gamma+Phi)/2.)*zHat);
    nHat = R*xHat*R.conjugate();
    lambdaHat = R*yHat*R.conjugate();
    chi1n = chi1Q.dot(nHat);
    chi2n = chi2Q.dot(nHat);
    chi1la = chi1Q.dot(lambdaHat);
    chi2la = chi2Q.dot(lambdaHat);
    chi1l = chi1Q.dot(OmegaHat_orbQ);
    chi2l = chi2Q.dot(OmegaHat_orbQ);
    chi1chi1 = chi1Q.dot(chi1Q);
    chi1chi2 = chi1Q.dot(chi2Q);
    chi2chi2 = chi2Q.dot(chi2Q);
    return;
  }
  vector<double> Omega_prec() const {
    vector<double> O(3);
    const double Omega_precMagnitude = Omega_precMag();
    O[0] = Omega_precMagnitude*nHat[1];
    O[1] = Omega_precMagnitude*nHat[2];
    O[2] = Omega_precMagnitude*nHat[3];
    return O;
  }
  vector<double> L() const {
    // Eq. (4.7) of Bohé et al.(2012)
    // <http://arxiv.org/abs/1212.5520v1>.
    const GWFrames::Quaternion l =
      powv3*(chi2n*(0.5*(1. - 1.*delta) - 1.5*nu) + chi1n*(0.5*(1. + delta) - 1.5*nu)
	     + powv2*(chi1n*(1.375*(1. + delta) + 0.02083333333333333*(-227. - 29.*delta)*nu + 1.625*pownu2)
		      + chi2n*(-1.375*(-1. + delta) + 0.02083333333333333*(-227. + 29.*delta)*nu + 1.625*pownu2)
		      + powv2*(chi2n*(-3.8125*(-1.+delta) + 0.01041666666666667*(-3163.+2065.*delta)*nu + 0.02083333333333333*(2807.-8.*delta)*pownu2-0.4375*pownu3)
			       + chi1n*(3.8125*(1.+delta) + 0.01041666666666667*(-3163.-2065.*delta)*nu + 0.02083333333333333*(2807.+8.*delta)*pownu2-0.4375*pownu3)
			       )))*nHat
      + 
      powv3*(chi2la*(2.*(-1. + delta) + 5.*nu) + chi1la*(-2.*(1. + delta) + 5.*nu)
	     + powv2*(chi2la*(2.*(-1. + delta) + 0.1666666666666667*(40. - 13.*delta)*nu - 5.666666666666667*pownu2)
		      + chi1la*(-2.*(1. + delta) + 0.1666666666666667*(40. + 13.*delta)*nu - 5.666666666666667*pownu2)
		      + powv2*(chi1la*(-3.875*(1. + delta) - 0.2291666666666667*(-29. + 7.*delta)*nu + (5.5 - 1.*delta)*pownu2 + 2.666666666666667*pownu3)
			       + chi2la*(3.875*(-1. + delta) + 0.2291666666666667*(29. + 7.*delta)*nu + (5.5 + delta)*pownu2 + 2.666666666666667*pownu3)
			       )))*lambdaHat
      +
      (1 + powv2*(1.5 + 0.1666666666666667*nu
		  + v*(chi2l*(4.166666666666667*(-1. + delta) + 10.83333333333333*nu)
		       + chi1l*(-4.166666666666667*(1. + delta) + 10.83333333333333*nu)
		       + v*(3.375 - 2.375*nu + 0.04166666666666667*pownu2
			    +v*(chi2l*(6.125*(-1. + delta) - 0.04861111111111111*(-397. + 91.*delta)*nu - 11.76388888888889*pownu2)
				+ chi1l*(-6.125*(1. + delta) + 0.04861111111111111*(397. + 91.*delta)*nu - 11.76388888888889*pownu2)
				+ v*(8.4375 - 30.97970359258346*nu + 1.291666666666667*pownu2 + 0.005401234567901235*pownu3
				     + v*(chi1l*(-15.1875*(1. + delta) + 0.09375*(901. + 523.*delta)*nu + 0.0625*(-2059. - 22.*delta)*pownu2 + 3.6875*pownu3)
					  + chi2l*(15.1875*(-1. + delta) - 0.09375*(-901. + 523.*delta)*nu + 0.0625*(-2059. + 22.*delta)*pownu2 + 3.6875*pownu3)
					  )))))))*GWFrames::Quaternion(std::vector<double>(OmegaHat_orb,OmegaHat_orb+3));
    
      // (0.005208333333333333*chi2n* (-96.*(-1. + delta + 3.*nu)*powv3 -  2.*(-1. + delta)* (33. - 19.*nu + delta*(-99. + 39.*nu))*powv5 + (-1. + delta)* (-183. - 11.*(-121. + nu)*nu +  3.*delta*(183. + nu*(-933. + 7.*nu)))* powv7) + 0.005208333333333333*chi1n* (96.*(1. + delta - 3.*nu)*powv3 -  2.*(1. + delta)* (-33. + 19.*nu + delta*(-99. + 39.*nu))* powv5 + (1. + delta)* (183. + 11.*(-121. + nu)*nu +  3.*delta*(183. + nu*(-933. + 7.*nu)))* powv7))*nHat
      // + (1. + 0.1666666666666667*(9. + nu)*powv2 +  0.04166666666666667*(81. + (-57. + nu)*nu)* powv4 + 0.0007716049382716049* (10935. + nu*(-40149.69585598816 + 1674.*nu +  7.*pownu2))*powv6 +  0.003472222222222222*chi1l* (-240.*(5. + 5.*delta - 13.*nu)*powv3 +  7.*(1. + delta)* (-99. + 61.*nu + delta*(-153. + 121.*nu))* powv5 - 9.*(486. +  nu*(-2703. + 2.*(2059. - 59.*nu)*nu) +  delta*(486. + nu*(-1569. + 44.*nu)))*powv7 ) + 0.003472222222222222*chi2l* (240.*(-5. + 5.*delta + 13.*nu)*powv3 +  7.*(-1. + delta)* (99. - 61.*nu + delta*(-153. + 121.*nu))* powv5 + 9.*(-486. +  delta*(486. + nu*(-1569. + 44.*nu)) +  nu*(2703. + 2.*nu*(-2059. + 59.*nu)))* powv7))*GWFrames::Quaternion(std::vector<double>(OmegaHat_orb,OmegaHat_orb+3))
      // + (0.02083333333333333*chi1la* (-48.*(2. + 2.*delta - 5.*nu)*powv3 +  2.*(1. + delta)* (-21. - 27.*delta + 18.*nu + 34.*delta*nu)* powv5 + (-186. +  nu*(319. + 8.*nu*(33. + 16.*nu)) -  1.*delta*(186. + nu*(77. + 48.*nu)))*powv7 ) + 0.02083333333333333*chi2la* (48.*(-2. + 2.*delta + 5.*nu)*powv3 +  2.*(-1. + delta)* (21. - 18.*nu + delta*(-27. + 34.*nu))*powv5 + (-186. +  nu*(319. + 8.*nu*(33. + 16.*nu)) +  delta*(186. + nu*(77. + 48.*nu)))*powv7))* lambdaHat;
    return ((nu/v)*l).vec();
  }
  int RHS(double t, const double* y, double* dydt) {
    RecalculateValues(t, y);
    if(y[0]>=1.0) { return GSL_ETOLX; } // Stop integrating if v is greater than or equal to 1.0
    dydt[0] = - (Flux() + dMdt()) / dEdv();
    if(dydt[0]<0.0) { return GSL_ETOLF; } // Stop integrating if v is decreasing
    dydt[1] = powv3;
    const std::vector<double> Omega_spin1 = Omega_spin(1);
    const std::vector<double> Omega_spin2 = Omega_spin(2);
    const std::vector<double> Omega_precVec = Omega_prec();
    cross(&Omega_spin1[0], chi1, &dydt[2]);
    cross(&Omega_spin2[0], chi2, &dydt[5]);
    cross(&Omega_precVec[0], OmegaHat_orb, &dydt[8]);
    
    // const double k = (nu/v) * (1.0 + powv2*( 1.5+nu/6.0 + powv2*( 3.375 + nu*(-2.375+0.0416666666667*nu) ) ) );
    // dydt[8] = (-1.0/k)*(*dydt[2] + *dydt[5]);
    
    const GWFrames::Quaternion adot(0., dydt[8], dydt[9], dydt[10]);
    const GWFrames::Quaternion Raxdot = ( (-1.0/std::sqrt(2+2*y[10]))*adot*zHat - (dydt[10]/(2+2*y[10]))*Rax );
    const GWFrames::Quaternion dydt11 = 2*(Rax.conjugate() * Raxdot * zHat);
    dydt[11] = dydt11[0];
    return GSL_SUCCESS;
  }
};


#endif // DOXYGEN
