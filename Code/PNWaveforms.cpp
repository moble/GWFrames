// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

#include <unistd.h>
#include <sys/param.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <gsl/gsl_odeiv2.h>
#include "PNWaveforms.hpp"
#include "Quaternions.hpp"
#include "Utilities.hpp"
#include "Errors.hpp"

using std::vector;
using std::string;

// Local utility functions
std::string VectorStringForm(const std::vector<double>& V) {
  std::stringstream S;
  S << "[" << std::setprecision(16);
  for(unsigned int i=0; i<V.size()-1; ++i) {
    S << V[i] << ", ";
  }
  S << V[V.size()-1] << "]";
  return S.str();
}
inline double dotproduct(const double* a, const double* b) {
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

/// Default constructor for an empty object
GWFrames::PNWaveform::PNWaveform() :
  Waveform(), mchi1(0), mchi2(0), mOmega_orb(0), mOmega_prec(0), mL(0), mPhi_orb(0)
{
  SetFrameType(GWFrames::Coorbital);
  SetDataType(GWFrames::h);
  SetRIsScaledOut(true);
  SetMIsScaledOut(true);
  { // Overwrite the history from Waveform
    char path[MAXPATHLEN];
    getcwd(path, MAXPATHLEN);
    string pwd = path;
    char host[MAXHOSTNAMELEN];
    gethostname(host, MAXHOSTNAMELEN);
    string hostname = host;
    time_t rawtime;
    time ( &rawtime );
    string date = asctime ( localtime ( &rawtime ) );
    history.str("");
    history.clear();
    history << "### Code revision (`git rev-parse HEAD` or arXiv version) = " << CodeRevision << std::endl
	    << "### pwd = " << pwd << std::endl
	    << "### hostname = " << hostname << std::endl
	    << "### date = " << date // comes with a newline
	    << "### PNWaveform(); // empty constructor" << std::endl;
  }
}

/// Copy constructor
GWFrames::PNWaveform::PNWaveform(const PNWaveform& a) :
  Waveform(a), mchi1(a.mchi1), mchi2(a.mchi2), mOmega_orb(a.mOmega_orb), mOmega_prec(a.mOmega_prec), mL(a.mL), mPhi_orb(a.mPhi_orb)
{
  /// Simply copies all fields in the input object to the constructed
  /// object, including history
  history.seekp(0, std::ios_base::end);
}

// This file contains all the equations for evolving the dynamics of
// the PN system.
#include "PNWaveforms_TaylorT1Spin.ipp"

// This file contains a class for calculating the polarization modes
#include "PNWaveforms_PolarizationModes.ipp"

// This will be the right-hand side for the ODE integration; params
// will point to a TaylorT1 object.
int func (double t, const double y[], double dydt[], void* params) {
  TaylorT1* T1 = (TaylorT1*) params;
  return T1->RHS(t, y, dydt);
}

void WaveformModes(const double delta, const double v, const double chisl, const double chial, std::vector<std::complex<double> >& modes);

/// Constructor of PN waveform from parameters
GWFrames::PNWaveform::PNWaveform(const double delta,
				 const std::vector<double>& chi1_0, const std::vector<double>& chi2_0,
				 const double Omega_orb_0, const GWFrames::Quaternion& R_0) :
  Waveform(), mchi1(0), mchi2(0), mOmega_orb(0), mOmega_prec(0), mL(0), mPhi_orb(0)
{
  /// 
  /// \param delta Normalized BH mass difference (M1-M2)/(M1+M2)
  /// \param chi1_0 Initial dimensionless spin vector of BH1
  /// \param chi2_0 Initial dimensionless spin vector of BH2
  /// \param Omega_orb_0 Initial orbital angular frequency
  /// \param R_0 Overall rotation of the system (optional)
  /// 
  /// The PN system is initialized having the BHs along the x axis,
  /// with the orbital angular velocity along the positive z axis,
  /// having magnitude Omega_orb_0.  The input spin vectors must be
  /// defined with respect to this basis.
  /// 
  /// The TaylorT1 system is first integrated to compute the dynamics
  /// of the binary.  The evolved spin vectors chi1 and chi2, orbital
  /// angular-velocity vector Omega_orb, and orbital phase Phi_orb are
  /// stored.  Simultaneously, the minimal-rotation frame of the
  /// angular-velocity vector is computed, then rotated about the z'
  /// axis by Phi_orb, resulting in the binary's frame.  Once this
  /// step is completed, the information is used to construct the
  /// waveform in the minimal-rotation frame.  (That is, the waveform
  /// will be essentially corotating.)
  /// 
  /// Note that, to get the PNWaveform in an inertial frame, you must
  /// first apply the method TransformToCorotatingFrame().

  std::cerr << "WARNING: The precession equations for `PNWaveform` are incorrect.  Non-precessing systems"
            << "         should be correct.  Otherwise, however, take these as simply test data for other code."
            << std::endl;
  
  SetFrameType(GWFrames::Coorbital);
  SetDataType(GWFrames::h);
  SetRIsScaledOut(true);
  SetMIsScaledOut(true);
  
  { // Overwrite the history from Waveform
    char path[MAXPATHLEN];
    getcwd(path, MAXPATHLEN);
    string pwd = path;
    char host[MAXHOSTNAMELEN];
    gethostname(host, MAXHOSTNAMELEN);
    string hostname = host;
    time_t rawtime;
    time ( &rawtime );
    string date = asctime ( localtime ( &rawtime ) );
    history.str("");
    history.clear();
    history << "### Code revision (`git rev-parse HEAD` or arXiv version) = " << CodeRevision << std::endl
	    << "### pwd = " << pwd << std::endl
	    << "### hostname = " << hostname << std::endl
	    << "### date = " << date // comes with a newline
	    << "### PNWaveform(" << delta << ", " << VectorStringForm(chi1_0) << ", " << VectorStringForm(chi2_0)
	    << ", " << Omega_orb_0 << ", " << R_0 << ");" << std::endl;
  }
  
  vector<double> ystart(12);
  const GWFrames::Quaternion chi1_0Q = R_0 * GWFrames::Quaternion(chi1_0) * R_0.conjugate();
  const GWFrames::Quaternion chi2_0Q = R_0 * GWFrames::Quaternion(chi2_0) * R_0.conjugate();
  const GWFrames::Quaternion LNHatQ  = R_0 * zHat * R_0.conjugate();
  const GWFrames::Quaternion Rax = GWFrames::sqrtOfRotor(-LNHatQ*zHat);
  ystart[0] = std::pow(Omega_orb_0, 1./3.);           // v
  ystart[1] = 0.0;                                    // Phi
  ystart[2] = chi1_0Q[1];                             // chi1_x
  ystart[3] = chi1_0Q[2];                             // chi1_y
  ystart[4] = chi1_0Q[3];                             // chi1_z
  ystart[5] = chi2_0Q[1];                             // chi2_x
  ystart[6] = chi2_0Q[2];                             // chi2_y
  ystart[7] = chi2_0Q[3];                             // chi2_z
  ystart[8] = LNHatQ[1];                              // LNHat_x
  ystart[9] = LNHatQ[2];                              // LNHat_y
  ystart[10] = LNHatQ[3];                             // LNHat_z
  ystart[11] = GWFrames::angle(Rax.conjugate()*R_0);  // gamma
  
  const unsigned int MinSteps = 100000; // This is only an approximate lower limit
  const unsigned int MaxSteps = 10000000; // This is a hard upper limit
  
  const double nu = (1.0-delta*delta)/4.0;
  double time = -5.0/(256.0*nu*std::pow(ystart[0],8)); // This is the lowest-order pN time-to-merger
  double endtime = -3*time; // Give ourselves a large margin of error in case inspiral runs longer than expected
  double h = 1.0e0;
  
  TaylorT1 T1(delta, chi1_0, chi2_0, ystart[0], R_0);
  T1.RecalculateValues(time, &ystart[0]);
  const double eps_abs = 1.e-15;
  const double eps_rel = 1.e-10;
  const double hmin = 1.0e-3;
  const double hmax = (endtime-time) / (2.0*MinSteps);
  const GWFrames::Quaternion zHat(0,0,0,1);
  
  vector<double> y(ystart);
  vector<double> v;
  
  // Reserve plenty of space in the vectors
  v.reserve(MinSteps);
  mchi1.reserve(MinSteps);
  mchi2.reserve(MinSteps);
  mOmega_orb.reserve(MinSteps);
  mOmega_prec.reserve(MinSteps);
  mL.reserve(MinSteps);
  mPhi_orb.reserve(MinSteps);
  frame.reserve(MinSteps);
  t.reserve(MinSteps);
  
  // Declare and initialize the GSL ODE integrator
  const gsl_odeiv2_step_type* T = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_step* s = gsl_odeiv2_step_alloc(T, 12);
  gsl_odeiv2_control* c = gsl_odeiv2_control_y_new(eps_abs, eps_rel);
  gsl_odeiv2_evolve* e = gsl_odeiv2_evolve_alloc(12);
  gsl_odeiv2_system sys = {func, NULL, 12, (void *) &T1};
  
  // Store the data at the first step
  {
    v.push_back(y[0]);
    mchi1.push_back(std::vector<double>(&y[2],&y[5]));
    mchi2.push_back(std::vector<double>(&y[5],&y[8]));
    const double Omega_orbMag = y[0]*y[0]*y[0];
    double Omega_orb[3] = {Omega_orbMag*y[8], Omega_orbMag*y[9], Omega_orbMag*y[10]};
    mOmega_orb.push_back(std::vector<double>(Omega_orb,Omega_orb+3));
    T1.RecalculateValues(time, &y[0]);
    mOmega_prec.push_back(T1.Omega_prec());
    mL.push_back(T1.L());
    mPhi_orb.push_back(y[1]);
    frame.push_back(GWFrames::sqrtOfRotor(-GWFrames::Quaternion(mOmega_orb.back()).normalized()*zHat) * GWFrames::exp(((y[11]+y[1])/2.)*zHat));
    t.push_back(time);
  }
  
  // Run the integration
  unsigned int NSteps = 0;
  while (time < endtime) {
    // Take a step
    int status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &time, time+hmax, &h, &y[0]);
    ++NSteps;
    
    // Check if it worked and the system is still reasonable
    if(status == GSL_ETOLF) {
      std::cout << "Velocity v has become greater than 1.0" << std::endl;
      break;
    } else if(status == GSL_ETOLX) {
      std::cout << "Velocity is no longer increasing" << std::endl;
      break;
    } else if(status != GSL_SUCCESS) {
      std::cerr << "GSL odeiv2 error.  Return value=" << status << "\n" << std::endl;
      break;
    }
    
    // If so, store the data
    v.push_back(y[0]);
    mchi1.push_back(std::vector<double>(&y[2],&y[5]));
    mchi2.push_back(std::vector<double>(&y[5],&y[8]));
    const double Omega_orbMag = y[0]*y[0]*y[0];
    double Omega_orb[3] = {Omega_orbMag*y[8], Omega_orbMag*y[9], Omega_orbMag*y[10]};
    mOmega_orb.push_back(std::vector<double>(Omega_orb,Omega_orb+3));
    T1.RecalculateValues(time, &y[0]);
    mOmega_prec.push_back(T1.Omega_prec());
    mL.push_back(T1.L());
    mPhi_orb.push_back(y[1]);
    frame.push_back(GWFrames::sqrtOfRotor(-GWFrames::Quaternion(mOmega_orb.back()).normalized()*zHat) * GWFrames::exp(((y[11]+y[1])/2.)*zHat));
    t.push_back(time);
    
    // Check if we should stop because this has gone on suspiciously long
    if(time>=endtime) {
      std::cerr << "Time has gone on four times as long as expected.  This seems strange, so we'll stop."
		<< "\nNote that this is unusual.  You may have a short waveform that stops before merger." << std::endl;
      break;
    }
    
    // Check if we should stop because there have been too many steps
    if(NSteps>MaxSteps) {
      std::cerr << "\n\nThe integration has taken " << NSteps << ".  This seems excessive, so we'll stop."
		<< "\nNote that this is unusual.  You may have a short waveform that stops before merger." << std::endl;
      break;
    }
    
    // Check if we should stop because the step has gotten too small,
    // but make sure we at least take 10 steps.  [This is the
    // condition that we expect to stop us near merger.]
    if(NSteps>10 && h<hmin) { break; }
  }
  
  // Free the gsl storage
  gsl_odeiv2_evolve_free(e);
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(s);
  
  // Make the last time = 0.0
  const double tback = t.back();
  for(unsigned int i=0; i<t.size(); ++i) {
    t[i] -= tback;
  }
  
  // Set up the (ell,m) data
  // We need (2*ell+1) coefficients for each value of ell from 2 up to
  // ellMax_PNWaveforms.  That's a total of
  //   >>> from sympy import summation, symbols
  //   >>> ell, ellMax_PNWaveforms = symbols('ell ellMax_PNWaveforms', integer=True)
  //   >>> summation(2*ell+1, (ell, 2, ellMax_PNWaveforms))
  //   ellMax_PNWaveforms**2 + 2*ellMax_PNWaveforms - 3
  // The reverse process of finding the index of element (ell,m) is
  // done by taking the element of the array with index
  //   >>> summation(2*ell+1, (ell, 2, ell-1)) + ell + m
  //   ell**2 + ell + m - 4
  lm.resize(ellMax_PNWaveforms*(ellMax_PNWaveforms+2)-3, vector<int>(2,0));
  {
    unsigned int i=0;
    for(int ell=2; ell<=ellMax_PNWaveforms; ++ell) {
      for(int m=-ell; m<=ell; ++m) {
	lm[i][0] = ell;
	lm[i][1] = m;
	++i;
      }
    }
  }
  
  // Evaluate the waveform data itself, noting that we always use the
  // frame in standard position (BHs on the x axis, with angular
  // velocity along the positive z axis).  This can then be
  // transformed to a stationary frame, using the stored 'frame' data.
  data.resize(lm.size(), t.size());
  std::vector<std::complex<double> > modes(lm.size());
  PNWaveformFromDynamics rhOverM(delta);
  for(unsigned int i_t=0; i_t<t.size(); ++i_t) {
    const double v_t = v[i_t];
    const double Omega_orb = v_t*v_t*v_t;
    const double chi1l_t = dotproduct(&mchi1[i_t][0], &mOmega_orb[i_t][0]) / Omega_orb;
    const double chi2l_t = dotproduct(&mchi2[i_t][0], &mOmega_orb[i_t][0]) / Omega_orb;
    const double chisl_t = 0.5*(chi1l_t+chi2l_t);
    const double chial_t = 0.5*(chi1l_t-chi2l_t);
    rhOverM(v_t, chisl_t, chial_t, modes);
    for(unsigned int i_lm=0; i_lm<lm.size(); ++i_lm) {
      data[i_lm][i_t] = modes[i_lm];
    }
  }
  
} // end PN constructor


/// Total angular velocity of PN binary at an instant of time
std::vector<double> GWFrames::PNWaveform::Omega_tot(const unsigned int iTime) const {
  std::vector<double> Tot(3);
  Tot[0] = mOmega_orb[iTime][0]+mOmega_prec[iTime][0];
  Tot[1] = mOmega_orb[iTime][1]+mOmega_prec[iTime][1];
  Tot[2] = mOmega_orb[iTime][2]+mOmega_prec[iTime][2];
  return Tot;
}

/// Total angular velocity of PN binary at each instant of time
std::vector<std::vector<double> > GWFrames::PNWaveform::Omega_tot() const {
  const double NT = NTimes();
  std::vector<std::vector<double> > Tot(NT, std::vector<double>(3));
  for(unsigned int i_t=0; i_t<NT; ++i_t) {
    Tot[i_t] = Omega_tot(i_t);
  }
  return Tot;
}

/// Vector of magnitudes of Omega_orb at each instant of time
std::vector<double> GWFrames::PNWaveform::Omega_orbMag() const {
  const double NT = NTimes();
  std::vector<double> Mag(NT);
  for(unsigned int i=0; i<NT; ++i) {
    Mag[i] = GWFrames::abs(mOmega_orb[i]);
  }
  return Mag;
}

/// Vector of magnitudes of Omega_prec at each instant of time
std::vector<double> GWFrames::PNWaveform::Omega_precMag() const {
  const double NT = NTimes();
  std::vector<double> Mag(NT);
  for(unsigned int i=0; i<NT; ++i) {
    Mag[i] = GWFrames::abs(mOmega_prec[i]);
  }
  return Mag;
}

/// Vector of magnitudes of Omega_tot at each instant of time
std::vector<double> GWFrames::PNWaveform::Omega_totMag() const {
  const double NT = NTimes();
  std::vector<double> Mag(NT);
  for(unsigned int i=0; i<NT; ++i) {
    Mag[i] = GWFrames::abs(Omega_tot(i));
  }
  return Mag;
}

/// Vector of magnitudes of angular momentum L at each instant of time
std::vector<double> GWFrames::PNWaveform::LMag() const {
  const double NT = NTimes();
  std::vector<double> Mag(NT);
  for(unsigned int i=0; i<NT; ++i) {
    Mag[i] = GWFrames::abs(mL[i]);
  }
  return Mag;
}
