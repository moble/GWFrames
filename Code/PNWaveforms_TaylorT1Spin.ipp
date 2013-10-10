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
  double chi1_n, chi2_n, chi1_la, chi2_la, chi1_l, chi2_l, chi1chi1, chi1chi2, chi2chi2; // The remaining dot products
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
    const double S_l = 0.25*chi1_l*pow(delta + 1.0, 2) + 0.25*chi2_l*pow(-delta + 1.0, 2);
    const double S_n = 0.25*chi1_n*pow(delta + 1.0, 2) + 0.25*chi2_n*pow(-delta + 1.0, 2);
    const double Sigma_l = -0.5*chi1_l*(delta + 1.0) + 0.5*chi2_l*(-delta + 1.0);
    const double Sigma_n = -0.5*chi1_n*(delta + 1.0) + 0.5*chi2_n*(-delta + 1.0);
    const double gamma = 5.51146384479718e-6*pow(v, 2)*(pow(v, 2)*(-60480.0*nu + v*(302400.0*S_l + 181440.0*Sigma_l*delta + v*(-982800.0*nu + v*(40320.0*S_l*(4.0*nu + 15.0) + 362880.0*Sigma_l*delta + v*(-nu*(-560.0*nu*(4.0*nu + 2061.0) + 541013.822520207) - 15120.0*v*(S_l*(nu*(72.0*nu + 127.0) - 60.0) + 2.0*Sigma_l*delta*(nu*(16.0*nu + 61.0) - 18.0)) + 181440.0)) + 181440.0)) + 181440.0) + 181440.0);
    const double a_l_over_vcubed = 0.0138888888888889*pow(v, 4)*(504.0*S_n + 216.0*Sigma_n*delta + pow(v, 2)*(-24.0*S_n*(29.0*nu + 30.0) - 108.0*Sigma_n*delta*(3.0*nu + 4.0) + pow(v, 2)*(2.0*S_n*(nu*(208.0*nu + 531.0) + 54.0) + 3.0*Sigma_n*delta*(nu*(68.0*nu + 219.0) + 36.0))));
    return gamma*a_l_over_vcubed;
  }
  double Flux() const {
    // Eqs. (C7) -- (C13) of <http://arxiv.org/abs/0810.5336v3> [NOTE VERSION NUMBER!!!]
    const double F_2 = 0.002976190476190476*(-1247. - 980.*nu);
    const double F_3 = 12.56637061435917 + chi2_l*(1.375*(-1. + delta) + 1.5*nu) + chi1_l*(-1.375*(1. + delta) + 1.5*nu);
    const double F_4 = 1.494791666666667*chi1_l*chi1_l*(1. + delta - 2.*nu) + 6.020833333333333*chi1_l*chi2_l*nu - 1.494791666666667*chi2_l*chi2_l*(-1. + delta + 2.*nu)
      + 0.005208333333333333*(-89.*chi1chi1*(1. + delta - 2.*nu) - 412.*chi1chi2*nu + 89.*chi2chi2*(-1. + delta + 2.*nu))
      + 0.0001102292768959436*(-44711. + 18.*nu*(9271. + 1820.*nu));
    const double F_5 = -0.004674989067841954*(8191. + 16324.*nu) + 0.003472222222222222*chi2_l*(531.*(-1. + delta) + 4.*(908. - 701.*delta)*nu - 2512.*pownu2)
      + 0.003472222222222222*chi1_l*(-531.*(1. + delta) + 4.*(908. + 701.*delta)*nu - 2512.*pownu2);
    const double F_6 = 138.3349161635852 + 0.2617993877991494*chi2_l*(-65. + 65.*delta + 68.*nu) + 0.2617993877991494*chi1_l*(-65.*(1. + delta) + 68.*nu)
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
				* (0.125*chi2_l*(-1.*(1. + 3.*chi2chi2)*(1. - 1.*delta + (-3. + delta)*nu) + 3.*chi1chi1*(-1. + delta + 3.*nu - 1.*delta*nu)
						+ 3.*chi1chi2*(1. + delta - 1.*(3. + delta)*nu))
				   + 0.125*chi1_l*(-3.*chi1chi1*(1. + delta) + 3.*chi1chi1*(3. + delta)*nu + 3.*chi1chi2*(1. - 1.*delta + (-3. + delta)*nu)
						  + (1. + 3.*chi2chi2)*(-1. - 1.*delta + (3. + delta)*nu))) );
  }
  double dEdv() const {
    // Non-spin and spin-orbit terms come from Eq. (4.6) of Bohé et
    // al. (2012) <http://arxiv.org/abs/1212.5520v1>.  Spin-spin terms
    // come from Eqs. (C1) -- (C6) of
    // <http://arxiv.org/abs/0810.5336v3> [NOTE VERSION NUMBER!!!].
    const double dEdv_2 = 0.1666666666666667*(-9. - 1.*nu);
    const double dEdv_3 = 1.666666666666667*chi1_l*(2. + 2.*delta - 1.*nu) - 1.666666666666667*chi2_l*(-2. + 2.*delta + nu);
    const double dEdv_4 = -2.25*chi1_l*chi1_l*(1. + delta - 2.*nu) - 9.*chi1_l*chi2_l*nu + 2.25*chi2_l*chi2_l*(-1. + delta + 2.*nu)
      + 0.125*(-81. + 6.*chi1chi1*(1. + delta - 2.*nu) + (57. + 24.*chi1chi2 - 1.*nu)*nu - 6.*chi2chi2*(-1. + delta + 2.*nu));
    const double dEdv_5 = 0.1944444444444444*chi1_l*(72. + delta*(72. - 31.*nu) + nu*(-121. + 2.*nu))
      + 0.1944444444444444*chi2_l*(72. + nu*(-121. + 2.*nu) + delta*(-72. + 31.*nu));
    const double dEdv_6 = -0.003858024691358025*(10935. + nu*(-40149.69585598816 + nu*(1674. + 7.*nu)));
    const double dEdv_7 = -0.1875*chi2_l*(-324. + nu*(1119. - 2.*nu*(172. + nu)) + delta*(324. + nu*(-633. + 14.*nu)))
      + 0.1875*chi1_l*(324. + nu*(-1119. + 2.*nu*(172. + nu)) + delta*(324. + nu*(-633. + 14.*nu)));
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
    chi1_n(Quaternion(chi1_0).dot(nHat)), chi2_n(Quaternion(chi2_0).dot(nHat)),
    chi1_la(Quaternion(chi1_0).dot(lambdaHat)), chi2_la(Quaternion(chi2_0).dot(lambdaHat)),
    chi1_l(Quaternion(chi1_0).dot(R_0*zHat*R_0.conjugate())), chi2_l(Quaternion(chi2_0).dot(R_0*zHat*R_0.conjugate())),
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
    chi1_n = chi1Q.dot(nHat);
    chi2_n = chi2Q.dot(nHat);
    chi1_la = chi1Q.dot(lambdaHat);
    chi2_la = chi2Q.dot(lambdaHat);
    chi1_l = chi1Q.dot(OmegaHat_orbQ);
    chi2_l = chi2Q.dot(OmegaHat_orbQ);
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
      powv3*(chi2_n*(0.5*(1. - 1.*delta) - 1.5*nu) + chi1_n*(0.5*(1. + delta) - 1.5*nu)
	     + powv2*(chi1_n*(1.375*(1. + delta) + 0.02083333333333333*(-227. - 29.*delta)*nu + 1.625*pownu2)
		      + chi2_n*(-1.375*(-1. + delta) + 0.02083333333333333*(-227. + 29.*delta)*nu + 1.625*pownu2)
		      + powv2*(chi2_n*(-3.8125*(-1.+delta) + 0.01041666666666667*(-3163.+2065.*delta)*nu + 0.02083333333333333*(2807.-8.*delta)*pownu2-0.4375*pownu3)
			       + chi1_n*(3.8125*(1.+delta) + 0.01041666666666667*(-3163.-2065.*delta)*nu + 0.02083333333333333*(2807.+8.*delta)*pownu2-0.4375*pownu3)
			       )))*nHat
      +
      powv3*(chi2_la*(2.*(-1. + delta) + 5.*nu) + chi1_la*(-2.*(1. + delta) + 5.*nu)
	     + powv2*(chi2_la*(2.*(-1. + delta) + 0.1666666666666667*(40. - 13.*delta)*nu - 5.666666666666667*pownu2)
		      + chi1_la*(-2.*(1. + delta) + 0.1666666666666667*(40. + 13.*delta)*nu - 5.666666666666667*pownu2)
		      + powv2*(chi1_la*(-3.875*(1. + delta) - 0.2291666666666667*(-29. + 7.*delta)*nu + (5.5 - 1.*delta)*pownu2 + 2.666666666666667*pownu3)
			       + chi2_la*(3.875*(-1. + delta) + 0.2291666666666667*(29. + 7.*delta)*nu + (5.5 + delta)*pownu2 + 2.666666666666667*pownu3)
			       )))*lambdaHat
      +
      (1 + powv2*(1.5 + 0.1666666666666667*nu
		  + v*(chi2_l*(4.166666666666667*(-1. + delta) + 10.83333333333333*nu)
		       + chi1_l*(-4.166666666666667*(1. + delta) + 10.83333333333333*nu)
		       + v*(3.375 - 2.375*nu + 0.04166666666666667*pownu2
			    +v*(chi2_l*(6.125*(-1. + delta) - 0.04861111111111111*(-397. + 91.*delta)*nu - 11.76388888888889*pownu2)
				+ chi1_l*(-6.125*(1. + delta) + 0.04861111111111111*(397. + 91.*delta)*nu - 11.76388888888889*pownu2)
				+ v*(8.4375 - 30.97970359258346*nu + 1.291666666666667*pownu2 + 0.005401234567901235*pownu3
				     + v*(chi1_l*(-15.1875*(1. + delta) + 0.09375*(901. + 523.*delta)*nu + 0.0625*(-2059. - 22.*delta)*pownu2 + 3.6875*pownu3)
					  + chi2_l*(15.1875*(-1. + delta) - 0.09375*(-901. + 523.*delta)*nu + 0.0625*(-2059. + 22.*delta)*pownu2 + 3.6875*pownu3)
					  )))))))*GWFrames::Quaternion(std::vector<double>(OmegaHat_orb,OmegaHat_orb+3));
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
    const GWFrames::Quaternion adot(0., dydt[8], dydt[9], dydt[10]);
    const GWFrames::Quaternion Raxdot = ( (-1.0/std::sqrt(2+2*y[10]))*adot*zHat - (dydt[10]/(2+2*y[10]))*Rax );
    const GWFrames::Quaternion dydt11 = 2*(Rax.conjugate() * Raxdot * zHat);
    dydt[11] = dydt11[0];
    return GSL_SUCCESS;
  }

 int TaylorT1RHS(double t, const double* y, double* dydt) {
    // Stop integrating if v is greater than or equal to 1.0
    if(y[0]>=1.0) { return GSL_ETOLX; }

    // Temporary definitions
    const double m1 = (1+delta)/2.0;
    const double m2 = (1-delta)/2.0;
    const double nu__2 = pownu2;
    const double nu__3 = pownu3;

    // Fundamental atomic quantities (other than the Quaternions)
    const double v = y[0];
    const double Phi = y[1];
    const double chi1_x = y[2];
    const double chi1_y = y[3];
    const double chi1_z = y[4];
    const double chi2_x = y[5];
    const double chi2_y = y[6];
    const double chi2_z = y[7];
    const double Lhat_Nx = y[8];
    const double Lhat_Ny = y[9];
    const double Lhat_Nz = y[10];
    const double gamma = y[11];
    const GWFrames::Quaternion Rax =
	GWFrames::sqrtOfRotor(-GWFrames::normalized(GWFrames::Quaternion(0., Lhat_Nx, Lhat_Ny, Lhat_Nz))*zHat);
    const GWFrames::Quaternion R = Rax * GWFrames::exp(((gamma+Phi)/2.)*zHat);
    const GWFrames::Quaternion nhatQ = R*xHat*R.conjugate();
    const double nhat_x = nhatQ[1];
    const double nhat_y = nhatQ[2];
    const double nhat_z = nhatQ[3];

    // Non-fundamental atomic quantities
    const double chi1chi1 = pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z, 2);
    const double chi2chi2 = pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z, 2);
    const double chi1_l = Lhat_Nx*chi1_x + Lhat_Ny*chi1_y + Lhat_Nz*chi1_z;
    const double chi1_n = chi1_x*nhat_x + chi1_y*nhat_y + chi1_z*nhat_z;
    const double chi1_lambda = chi1_x*(Lhat_Ny*nhat_z - Lhat_Nz*nhat_y) + chi1_y*(-Lhat_Nx*nhat_z + Lhat_Nz*nhat_x) + chi1_z*(Lhat_Nx*nhat_y - Lhat_Ny*nhat_x);
    const double chi2_l = Lhat_Nx*chi2_x + Lhat_Ny*chi2_y + Lhat_Nz*chi2_z;
    const double chi2_n = chi2_x*nhat_x + chi2_y*nhat_y + chi2_z*nhat_z;
    const double chi2_lambda = chi2_x*(Lhat_Ny*nhat_z - Lhat_Nz*nhat_y) + chi2_y*(-Lhat_Nx*nhat_z + Lhat_Nz*nhat_x) + chi2_z*(Lhat_Nx*nhat_y - Lhat_Ny*nhat_x);
    const double sqrt1Mchi1chi1 = sqrt(-chi1chi1 + 1.0);
    const double sqrt1Mchi2chi2 = sqrt(-chi2chi2 + 1.0);
    const double S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    const double S_n = chi1_n*pow(m1, 2) + chi2_l*pow(m2, 2);
    const double Sigma_l = -chi1_l*m1 + chi2_l*m2;
    const double Sigma_n = -chi1_n*m1 + chi2_n*m2;
    const double logv = log(v);

    // Composite quantities
    const double dEnergydv = v*(-1.0*nu + pow(v, 2)*(nu*(0.166666666666667*nu + 1.5) + v*(-11.6666666666667*S_l*nu - 5.0*Sigma_l*delta*nu + v*(chi1_l*(chi1_l*(1.5*delta*nu + nu*(-3.0*nu + 1.5)) + 6.0*chi2_l*pow(nu, 2)) + chi1_lambda*(chi1_lambda*(-0.75*delta*nu + nu*(1.5*nu - 0.75)) - 3.0*chi2_lambda*pow(nu, 2)) + chi1_n*(chi1_n*(-0.75*delta*nu + nu*(1.5*nu - 0.75)) - 3.0*chi2_n*pow(nu, 2)) + pow(chi2_l, 2)*(-1.5*delta*nu + nu*(-3.0*nu + 1.5)) + pow(chi2_lambda, 2)*(0.75*delta*nu + nu*(1.5*nu - 0.75)) + pow(chi2_n, 2)*(0.75*delta*nu + nu*(1.5*nu - 0.75)) + nu*(-7.125*nu + 0.125*nu__2 + 10.125) + v*(S_l*nu*(23.7222222222222*nu - 38.5) + Sigma_l*delta*nu*(11.6666666666667*nu - 10.5) + v*(nu*(-154.898517962917*nu + 6.45833333333333*nu__2 + 0.0270061728395062*nu__3 + 42.1875) + v*(S_l*nu*(412.875*nu - 10.875*nu__2 - 151.875) + Sigma_l*delta*nu*(175.5*nu - 5.625*nu__2 - 30.375) + v*(-298.666666666667*logv*pow(nu, 2) + nu*(-769.4015*nu + nu__2*(-0.012377829218107*nu__2 + 450.663995131026) - 0.870949074074074*nu__3 + 155.0390625) + pow(v, 2)*(logv*pow(nu, 2)*(1523.2*nu + 1710.17142857143) + nu*(330.78*nu + 538.20703125) + pow(v, 2)*(16016.0*logv*pow(nu, 2) + nu*(-4116.0*nu + 1808.9736328125)))))))))));
    const double Flux = pow(v, 10)*(6.4*nu__2 + pow(v, 2)*(-18.6666666666667*nu*nu__2 - 23.752380952381*nu__2 + v*(25.6*S_l*nu__2 - 8.0*Sigma_l*delta*nu__2 + 80.4247719318987*nu__2 + v*(chi1_l*(chi1_l*(6.6*delta*nu__2 - 13.2*nu*nu__2 + 6.6*nu__2) + 24.8*chi2_l*nu*nu__2) + chi1_lambda*(chi1_lambda*(-2.96666666666667*delta*nu__2 + 5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) - 13.7333333333333*chi2_lambda*nu*nu__2) + chi1_n*(chi1_n*(-2.96666666666667*delta*nu__2 + 5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) - 13.7333333333333*chi2_n*nu*nu__2) + pow(chi2_l, 2)*(-6.6*delta*nu__2 - 13.2*nu*nu__2 + 6.6*nu__2) + pow(chi2_lambda, 2)*(2.96666666666667*delta*nu__2 + 5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) + pow(chi2_n, 2)*(2.96666666666667*delta*nu__2 + 5.93333333333333*nu*nu__2 - 2.96666666666667*nu__2) + 117.726984126984*nu*nu__2 + nu__2*(23.1111111111111*nu__2 - 31.542151675485) + v*(S_l*(193.422222222222*nu*nu__2 - 28.8*nu__2) + Sigma_l*delta*(68.8*nu*nu__2 - 5.2*nu__2) - 488.412937878093*nu*nu__2 - 245.074146910038*nu__2 + v*(-321.699087727595*S_l*nu__2 - 103.881997078702*Sigma_l*delta*nu__2 - 104.350476190476*logv*nu__2 - 56.7811420312465*nu*nu__2 + nu__2*(-199.794708994709*nu__2 - 15.3086419753086*nu__3 + 945.576192944022) - 204.89320622011*nu__2 + v*(S_l*(208.998941798942*nu*nu__2 + nu__2*(-666.074074074074*nu__2 + 448.343327454439)) + Sigma_l*delta*(93.9174603174603*nu*nu__2 + nu__2*(-266.844444444444*nu__2 + 181.619047619048)) + 2498.67153479682*nu*nu__2 + nu__2*(1285.7923710359*nu__2 - 649.661414142346) + v*(S_l*(3875.74795014869*nu*nu__2 - 729.896693184029*nu__2) + Sigma_l*delta*(1302.34474121815*nu*nu__2 - 214.316458834892*nu__2) + 337.555736961451*logv*nu__2 - 752.028100625135*nu__2 + v*(-1311.30675759439*logv*nu__2 + 4602.42139029395*nu__2 + v*(746.4952102023*logv*nu__2 - 7788.20474442907*nu__2 + v*(3031.19666031508*logv*nu__2 + 6137.18380876523*nu__2 + v*(3340.73494332217*logv*nu__2 + logv*(850.70483446712*logv*nu__2 - 15505.402624526*nu__2) + 13022.6558856344*nu__2))))))))))));
    const double Absorption = pow(v, 15)*(chi1_l*pow(m1, 3)*(-4.8*chi1chi1*nu__2 - 1.6*nu__2) + chi2_l*pow(m2, 3)*(-4.8*chi2chi2*nu__2 - 1.6*nu__2) + pow(v, 3)*(pow(m1, 4)*(chi1chi1*nu__2*(9.6*sqrt1Mchi1chi1 + 9.6) + nu__2*(3.2*sqrt1Mchi1chi1 + 3.2)) + pow(m2, 4)*(chi2chi2*nu__2*(9.6*sqrt1Mchi2chi2 + 9.6) + nu__2*(3.2*sqrt1Mchi2chi2 + 3.2))));
    const double dvdt = - (Flux + Absorption) / dEnergydv;
    const double gamma_Lhat_N = pow(v, 2)*(pow(v, 2)*(-0.333333333333333*nu + v*(1.66666666666667*S_l + Sigma_l*delta + v*(-5.41666666666667*nu + v*(S_l*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_l*delta + v*(nu*(nu*(0.0123456790123457*nu + 6.36111111111111) - 2.98177812235564) + v*(S_l*(nu*(-6.0*nu - 10.5833333333333) + 5.0) + Sigma_l*delta*(nu*(-2.66666666666667*nu - 10.1666666666667) + 3.0)) + 1.0)) + 1.0)) + 1.0) + 1.0);
    const double a_ell_overvcubed = pow(v, 4)*(7.0*S_n + 3.0*Sigma_n*delta + pow(v, 2)*(S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0) + pow(v, 2)*(S_n*(nu*(5.77777777777778*nu + 14.75) + 1.5) + Sigma_n*delta*(nu*(2.83333333333333*nu + 9.125) + 1.5))));
    const double Omega_Lhat_N = gamma_Lhat_N*a_ell_overvcubed;
    const double Omega1 = pow(v, 5)*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu - 0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75);
    const double Omega2 = pow(v, 5)*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu + 0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75);

    // Calculate derivatives of the fundamental variables
    dydt[0] = dvdt;
    if(dydt[0]<0.0) { return GSL_ETOLF; } // Stop integrating if v is decreasing
    dydt[1] = v*v*v;
    double Omega_spin1[3] = {Omega1*Lhat_Nx, Omega1*Lhat_Ny, Omega1*Lhat_Nz};
    double Omega_spin2[3] = {Omega2*Lhat_Nx, Omega2*Lhat_Ny, Omega2*Lhat_Nz};
    double Omega_prec[3]  = {Omega_Lhat_N*nhat_x, Omega_Lhat_N*nhat_y, Omega_Lhat_N*nhat_z};
    cross(&Omega_spin1[0], &y[2], &dydt[2]);
    cross(&Omega_spin2[0], &y[5], &dydt[5]);
    cross(&Omega_prec[0],  &y[8], &dydt[8]);
    const GWFrames::Quaternion adot(0., dydt[8], dydt[9], dydt[10]);
    const GWFrames::Quaternion Raxdot = ( (-1.0/std::sqrt(2+2*y[10]))*adot*zHat - (dydt[10]/(2+2*y[10]))*Rax );
    const GWFrames::Quaternion dgammadt = 2*(Rax.conjugate() * Raxdot * zHat);
    dydt[11] = dgammadt[0];

    // GSL expects this to be the returned quantity if everything went well
    return GSL_SUCCESS;
  }


  int TaylorT4RHS(double t, const double* y, double* dydt) {
    // Stop integrating if v is greater than or equal to 1.0
    if(y[0]>=1.0) { return GSL_ETOLX; }

    // Temporary definitions
    const double m1 = (1+delta)/2.0;
    const double m2 = (1-delta)/2.0;
    const double nu__2 = pownu2;
    const double nu__3 = pownu3;

    // Fundamental atomic quantities (other than the Quaternions)
    const double v = y[0];
    const double Phi = y[1];
    const double chi1_x = y[2];
    const double chi1_y = y[3];
    const double chi1_z = y[4];
    const double chi2_x = y[5];
    const double chi2_y = y[6];
    const double chi2_z = y[7];
    const double Lhat_Nx = y[8];
    const double Lhat_Ny = y[9];
    const double Lhat_Nz = y[10];
    const double gamma = y[11];
    const GWFrames::Quaternion Rax =
	GWFrames::sqrtOfRotor(-GWFrames::normalized(GWFrames::Quaternion(0., Lhat_Nx, Lhat_Ny, Lhat_Nz))*zHat);
    const GWFrames::Quaternion R = Rax * GWFrames::exp(((gamma+Phi)/2.)*zHat);
    const GWFrames::Quaternion nhatQ = R*xHat*R.conjugate();
    const double nhat_x = nhatQ[1];
    const double nhat_y = nhatQ[2];
    const double nhat_z = nhatQ[3];

    // Non-fundamental atomic quantities
    const double chi1chi1 = pow(chi1_x, 2) + pow(chi1_y, 2) + pow(chi1_z, 2);
    const double chi2chi2 = pow(chi2_x, 2) + pow(chi2_y, 2) + pow(chi2_z, 2);
    const double chi1_l = Lhat_Nx*chi1_x + Lhat_Ny*chi1_y + Lhat_Nz*chi1_z;
    const double chi1_n = chi1_x*nhat_x + chi1_y*nhat_y + chi1_z*nhat_z;
    const double chi1_lambda = chi1_x*(Lhat_Ny*nhat_z - Lhat_Nz*nhat_y) + chi1_y*(-Lhat_Nx*nhat_z + Lhat_Nz*nhat_x) + chi1_z*(Lhat_Nx*nhat_y - Lhat_Ny*nhat_x);
    const double chi2_l = Lhat_Nx*chi2_x + Lhat_Ny*chi2_y + Lhat_Nz*chi2_z;
    const double chi2_n = chi2_x*nhat_x + chi2_y*nhat_y + chi2_z*nhat_z;
    const double chi2_lambda = chi2_x*(Lhat_Ny*nhat_z - Lhat_Nz*nhat_y) + chi2_y*(-Lhat_Nx*nhat_z + Lhat_Nz*nhat_x) + chi2_z*(Lhat_Nx*nhat_y - Lhat_Ny*nhat_x);
    const double S_l = chi1_l*pow(m1, 2) + chi2_l*pow(m2, 2);
    const double S_n = chi1_n*pow(m1, 2) + chi2_l*pow(m2, 2);
    const double Sigma_l = -chi1_l*m1 + chi2_l*m2;
    const double Sigma_n = -chi1_n*m1 + chi2_n*m2;
    const double logv = log(v);

    // Composite quantities
    const double dvdt = pow(v, 9)*(pow(v, 2)*(-17.6*nu__2 + v*(v*(chi1_l*(-32.4*chi1_l*nu__2 + 63.2*chi2_l*nu__2) + chi1_lambda*(15.5333333333333*chi1_lambda*nu__2 - 32.9333333333333*chi2_lambda*nu__2) + chi1_n*(15.5333333333333*chi1_n*nu__2 - 32.9333333333333*chi2_n*nu__2) - 32.4*pow(chi2_l, 2)*nu__2 + 15.5333333333333*pow(chi2_lambda, 2)*nu__2 + 15.5333333333333*pow(chi2_n, 2)*nu__2 - 2.93333333333333*nu*nu__2 + 43.368253968254*nu__2 + v*(542.4*S_l*nu__2 + 224.8*Sigma_l*delta*nu__2 + chi1_l*pow(m1, 3)*(-4.8*chi1chi1*nu__2 - 1.6*nu__2)/nu + chi2_l*pow(m2, 3)*(-4.8*chi2chi2*nu__2 - 1.6*nu__2)/nu - 475.008809222777*nu__2 + v*(chi1_l*(chi1_l*(-23.7*delta*nu__2 + 47.4*nu*nu__2 - 29.8428571428571*nu__2) + chi2_l*(-95.0666666666667*nu*nu__2 + 9.88571428571429*nu__2)) + chi1_lambda*(chi1_lambda*(11.9055555555556*delta*nu__2 - 23.8111111111111*nu*nu__2 + 13.9769841269841*nu__2) + chi2_lambda*(47.3111111111111*nu*nu__2 - 6.94285714285714*nu__2)) + chi1_n*(chi1_n*(11.9055555555556*delta*nu__2 - 23.8111111111111*nu*nu__2 + 13.9769841269841*nu__2) + chi2_n*(47.3111111111111*nu*nu__2 - 6.94285714285714*nu__2)) + pow(chi2_l, 2)*(23.7*delta*nu__2 + 47.4*nu*nu__2 - 29.8428571428571*nu__2) + pow(chi2_lambda, 2)*(-11.9055555555556*delta*nu__2 - 23.8111111111111*nu*nu__2 + 13.9769841269841*nu__2) + pow(chi2_n, 2)*(-11.9055555555556*delta*nu__2 - 23.8111111111111*nu*nu__2 + 13.9769841269841*nu__2) + nu*(-0.488888888888889*nu*nu__2 + 128.228042328042*nu__2) + nu__2*(1.78518518518519*nu__2 - 1058.43868227317) + (S_l*(572.444444444444*S_l*nu__2 + 712.0*Sigma_l*delta*nu__2 - 1259.98809359975*nu__2) + Sigma_l*(200.0*Sigma_l*pow(delta, 2)*nu__2 - 506.005856738196*delta*nu__2) + pow(chi1_l, 2)*(3.07142857142857*delta*nu__2 + 3.07142857142857*nu__2) + pow(chi1_lambda, 2)*(-1.03571428571429*delta*nu__2 - 1.03571428571429*nu__2) + pow(chi1_n, 2)*(-1.03571428571429*delta*nu__2 - 1.03571428571429*nu__2) + pow(chi2_l, 2)*(-3.07142857142857*delta*nu__2 + 3.07142857142857*nu__2) + pow(chi2_lambda, 2)*(1.03571428571429*delta*nu__2 - 1.03571428571429*nu__2) + pow(chi2_n, 2)*(1.03571428571429*delta*nu__2 - 1.03571428571429*nu__2) - 104.350476190476*logv*nu__2 + nu__2*(-124.363756613757*nu__2 - 15.1358024691358*nu__3 + 1090.32725114508) - 204.89320622011*nu__2)/nu) + (-183.688888888889*S_l*nu__2 - 61.6380952380952*Sigma_l*delta*nu__2 - 124.43698901219*nu__2)/nu) + (pow(chi1_l, 2)*(16.2*delta*nu__2 + 16.2*nu__2) + pow(chi1_lambda, 2)*(-7.76666666666667*delta*nu__2 - 7.76666666666667*nu__2) + pow(chi1_n, 2)*(-7.76666666666667*delta*nu__2 - 7.76666666666667*nu__2) + pow(chi2_l, 2)*(-16.2*delta*nu__2 + 16.2*nu__2) + pow(chi2_lambda, 2)*(7.76666666666667*delta*nu__2 - 7.76666666666667*nu__2) + pow(chi2_n, 2)*(7.76666666666667*delta*nu__2 - 7.76666666666667*nu__2) + nu__2*(23.9111111111111*nu__2 + 12.0292768959436))/nu) + (-49.0666666666667*S_l*nu__2 - 40.0*Sigma_l*delta*nu__2 + 80.4247719318987*nu__2)/nu) - 14.152380952381*nu__2/nu) + 6.4*nu__2/nu);
    const double gamma_Lhat_N = pow(v, 2)*(pow(v, 2)*(-0.333333333333333*nu + v*(1.66666666666667*S_l + Sigma_l*delta + v*(-5.41666666666667*nu + v*(S_l*(0.888888888888889*nu + 3.33333333333333) + 2.0*Sigma_l*delta + v*(nu*(nu*(0.0123456790123457*nu + 6.36111111111111) - 2.98177812235564) + v*(S_l*(nu*(-6.0*nu - 10.5833333333333) + 5.0) + Sigma_l*delta*(nu*(-2.66666666666667*nu - 10.1666666666667) + 3.0)) + 1.0)) + 1.0)) + 1.0) + 1.0);
    const double a_ell_overvcubed = pow(v, 4)*(7.0*S_n + 3.0*Sigma_n*delta + pow(v, 2)*(S_n*(-9.66666666666667*nu - 10.0) + Sigma_n*delta*(-4.5*nu - 6.0) + pow(v, 2)*(S_n*(nu*(5.77777777777778*nu + 14.75) + 1.5) + Sigma_n*delta*(nu*(2.83333333333333*nu + 9.125) + 1.5))));
    const double Omega_Lhat_N = gamma_Lhat_N*a_ell_overvcubed;
    const double Omega1 = pow(v, 5)*(-0.75*delta + 0.5*nu + pow(v, 2)*(delta*(0.625*nu - 0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(-0.15625*nu + 4.875) - 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75);
    const double Omega2 = pow(v, 5)*(0.75*delta + 0.5*nu + pow(v, 2)*(delta*(-0.625*nu + 0.5625) + nu*(-0.0416666666666667*nu + 1.25) + pow(v, 2)*(delta*(nu*(0.15625*nu - 4.875) + 0.84375) + nu*(nu*(-0.0208333333333333*nu - 3.28125) + 0.1875) + 0.84375) + 0.5625) + 0.75);

    // Calculate derivatives of the fundamental variables
    dydt[0] = dvdt;
    if(dydt[0]<0.0) { return GSL_ETOLF; } // Stop integrating if v is decreasing
    dydt[1] = v*v*v;
    double Omega_spin1[3] = {Omega1*Lhat_Nx, Omega1*Lhat_Ny, Omega1*Lhat_Nz};
    double Omega_spin2[3] = {Omega2*Lhat_Nx, Omega2*Lhat_Ny, Omega2*Lhat_Nz};
    double Omega_prec[3]  = {Omega_Lhat_N*nhat_x, Omega_Lhat_N*nhat_y, Omega_Lhat_N*nhat_z};
    cross(&Omega_spin1[0], &y[2], &dydt[2]);
    cross(&Omega_spin2[0], &y[5], &dydt[5]);
    cross(&Omega_prec[0],  &y[8], &dydt[8]);
    const GWFrames::Quaternion adot(0., dydt[8], dydt[9], dydt[10]);
    const GWFrames::Quaternion Raxdot = ( (-1.0/std::sqrt(2+2*y[10]))*adot*zHat - (dydt[10]/(2+2*y[10]))*Rax );
    const GWFrames::Quaternion dgammadt = 2*(Rax.conjugate() * Raxdot * zHat);
    dydt[11] = dgammadt[0];

    // GSL expects this to be the returned quantity if everything went well
    return GSL_SUCCESS;
  }

};


#endif // DOXYGEN
