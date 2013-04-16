// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

#include <cmath>
#include <iostream>
#include <iomanip>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_odeiv2.h>

#include "Quaternions.hpp"
#include "Utilities.hpp"
#include "Errors.hpp"
using GWFrames::Quaternion;
using GWFrames::Quaternions;

// Note: Don't do 'using namespace std' because we don't want to
// confuse which log, exp, etc., is being used in any instance.
using std::cout;
using std::cerr;
using std::endl;
using std::flush;
using std::vector;

const Quaternion  one(1,0,0,0);
const Quaternion xHat(0,1,0,0);
const Quaternion yHat(0,0,1,0);
const Quaternion zHat(0,0,0,1);
inline double SQR(const double& x) { return x*x; }


#define Quaternion_Epsilon 1.0e-14


////////////////////////////////////
// Functions for the class itself //
////////////////////////////////////

/// Empty constructor -- initialized to 0s.
GWFrames::Quaternion::Quaternion()
  : w(0.0), x(0.0), y(0.0), z(0.0) { }

/// Copy constructor.
GWFrames::Quaternion::Quaternion(const Quaternion& Q)
  : w(Q.w), x(Q.x), y(Q.y), z(Q.z) { }

/// Constructor from spherical coordinates.
GWFrames::Quaternion::Quaternion(const double vartheta, const double varphi) {
  /// 
  /// \param vartheta Float representing the polar angle
  /// \param varphi Float representing the azimuthal angle
  /// 
  /// The unit Quaternion constructed in this way rotates the z axis
  /// onto the point given by the coordinates (vartheta, varphi).
  *this = GWFrames::exp((varphi/2.)*zHat) * GWFrames::exp((vartheta/2.)*yHat);
}

/// Constructor from Euler angles.
GWFrames::Quaternion::Quaternion(const double alpha, const double beta, const double gamma) {
  /// 
  /// \param alpha First Euler angle
  /// \param beta Second Euler angle
  /// \param gamma Third Euler angle
  /// 
  /// The unit Quaternion constructed in this way corresponds to a
  /// rotation by the given Euler angles.  The convention used here is
  /// the z-y-z convention.  That is, the rotations occur about the
  /// fixed axes: first a rotation by gamma about the z axis, then a
  /// rotation by beta about the y axis, and finally a rotation by
  /// alpha about the z axis.
  *this = GWFrames::exp((alpha/2.)*zHat) * GWFrames::exp((beta/2.)*yHat) * GWFrames::exp((gamma/2.)*zHat);
}

/// Constructor by components.
GWFrames::Quaternion::Quaternion(const double w0, const double x0, const double y0, const double z0)
  : w(w0), x(x0), y(y0), z(z0)
{
  ///
  /// \param w0 Scalar component of Quaternion
  /// \param x0 First vector component of Quaternion
  /// \param y0 Second vector component of Quaternion
  /// \param z0 Third vector component of Quaternion
}

/// Constructor from vector.
GWFrames::Quaternion::Quaternion(const std::vector<double>& q) {
  ///
  /// \param q Vector containing three or four components
  /// 
  /// If the input vector has three components, they are assumed to
  /// represent the vector components of the Quaternion, and the
  /// scalar component is set to zero.  If the input vector has four
  /// components, they are assumed to represent the four components of
  /// the Quaternion, with the 0 component being the scalar part.
  if(q.size()==3) {
    w = 0.0;
    x = q[0];
    y = q[1];
    z = q[2];
  } else if(q.size()==4) {
    w = q[0];
    x = q[1];
    y = q[2];
    z = q[3];
  } else {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": q.size()=" << q.size() << endl;
    throw(GWFrames_QuaternionVectorSizeNotUnderstood);
  }
}

/// Constructor from axis-angle.
GWFrames::Quaternion::Quaternion(const double angle, const std::vector<double>& axis)
  : w(std::cos(angle/2.)), x(std::sin(angle/2.)*axis[0]), y(std::sin(angle/2.)*axis[1]), z(std::sin(angle/2.)*axis[2])
{
  /// 
  /// \param angle Single number giving the rotation angle
  /// \param axis Three-component vector (assumed to be normalized) giving the axis
  /// 
  /// This constructs a rotor (assuming 'axis' is normalized)
  /// corresponding to rotation about the given axis through the given
  /// angle.
}

/// Get component of Quaternion.
double GWFrames::Quaternion::operator[](const unsigned int i) const {
  /// The 0 component is the scalar part, and the 1--3 components are
  /// the vector components.
  switch(i) {
  case 0:
    return w;
  case 1:
    return x;
  case 2:
    return y;
  case 3:
    return z;
  default:
    cerr << "\n" << __FILE__ << ":" << __LINE__ << ": i=" << i << " is not a possible quaternion index" << endl; 
    throw(GWFrames_IndexOutOfBounds);
  }
}

/// Get reference to component of Quaternion.
double& GWFrames::Quaternion::operator[](const unsigned int i) {
  /// Note: This is unavailable from python.
  switch(i) {
  case 0:
    return w;
  case 1:
    return x;
  case 2:
    return y;
  case 3:
    return z;
  default:
    cerr << "\n" << __FILE__ << ":" << __LINE__ << ": i=" << i << " is not a possible quaternion index" << endl; 
    throw(GWFrames_IndexOutOfBounds);
  }
}

/// Quaternion multiplication.
Quaternion GWFrames::Quaternion::operator*(const Quaternion& Q) const {
  return Quaternion(w*Q.w - x*Q.x - y*Q.y - z*Q.z,
		    w*Q.x + x*Q.w + y*Q.z - z*Q.y,
		    w*Q.y - x*Q.z + y*Q.w + z*Q.x,
		    w*Q.z + x*Q.y - y*Q.x + z*Q.w);
}

/// Return logarithm of Quaternion.
Quaternion GWFrames::Quaternion::log() const {
  Quaternion Result;
  const double b = std::sqrt(x*x + y*y + z*z);
  if(std::abs(b) <= Quaternion_Epsilon*std::abs(w)) {
    if(w<0.0) {
      cerr << "" << __FILE__ << ":" << __LINE__ << ": Infinitely many solutions for log of a negative scalar: w=" << w << "." << endl;
      throw(GWFrames_InfinitelyManySolutions);
    }
    Result.w = std::log(w);
  } else {
    const double v = std::atan2(b, w);
    const double f = v/b;
    Result.w = std::log(w/std::cos(v));
    Result.x = f*x;
    Result.y = f*y;
    Result.z = f*z;
  }
  return Result;
}

/// Return exponent of Quaternion.
Quaternion GWFrames::Quaternion::exp() const {
  Quaternion Result;
  const double b = std::sqrt(x*x + y*y + z*z);
  if(std::abs(b)<=Quaternion_Epsilon*std::abs(w)) {
    Result.w = std::exp(w);
  } else {
    const double e = std::exp(w);
    const double f = std::sin(b)/b; // Note: b is never 0.0 at this point
    Result.w = e*std::cos(b);
    Result.x = e*f*x;
    Result.y = e*f*y;
    Result.z = e*f*z;
  }
  return Result;
}

/// Print the quaternion nicely to stream
std::ostream& GWFrames::operator<<(std::ostream& out, const GWFrames::Quaternion& q) {
  out << "[" << q[0] << ", " << q[1] << ", " << q[2] << ", " << q[3] << "]";
  return out;
}


/////////////////////
// Array operators //
/////////////////////

/// Calculate the derivative of a rotor by the logarithm formula.
std::vector<Quaternion> GWFrames::DifferentiateRotorByLogarithm(const std::vector<Quaternion>& RIn, const std::vector<double>& tIn) {
  /// This is a much more complicated way of evaluating the derivative
  /// of a quaternion function of time, as compared to finite
  /// differencing by 'QuaternionDerivative'.  However, there may be
  /// accuracy advantages when the logarithm is smooth, and -- at the
  /// least -- this can serve as a good test of the correctness of the
  /// logarithm formula.
  const vector<Quaternion> logR = GWFrames::log(RIn);
  const vector<Quaternion> rdot = GWFrames::QuaternionDerivative(logR, tIn);
  vector<Quaternion> ROut(RIn.size());
  for(unsigned int i=0; i<logR.size(); ++i) {
    const double absquatlogR = abs(logR[i]);
    if(absquatlogR==0.0) { ROut[i] = rdot[i]; continue; }
    const double absquatlogRsquared = SQR(absquatlogR);
    const double a = SQR(std::sin(absquatlogR)/absquatlogR)/2.0;
    const double b = (absquatlogR<0.001
		      ? 0.6666666666666666 + absquatlogRsquared*(-0.13333333333333333 + absquatlogRsquared*(0.012698412698412698 + (-0.0007054673721340388 + (4*absquatlogRsquared)/155925.)*absquatlogRsquared))
		      : (absquatlogR-std::sin(absquatlogR)*std::cos(absquatlogR))/(absquatlogRsquared*absquatlogR) ) / 4.0;
    const Quaternion comm = GWFrames::commutator(logR[i],rdot[i]);
    ROut[i] = (rdot[i] + a*comm + b*GWFrames::commutator(logR[i],comm)) * RIn[i];
  }
  return ROut;
}


/// Minimal-rotation version of the input frame.
std::vector<Quaternion> GWFrames::MinimalRotation(const std::vector<Quaternion>& R, const std::vector<double>& T, const unsigned int NIterations) {
  ///
  /// \param R Vector of rotors.
  /// \param T Vector of corresponding time steps.
  /// \param NIterations Number of refinements [default: 5]
  /// 
  /// This function returns a copy of the input R, which takes the z
  /// axis to the same point as R, but adjusts the rotation about that
  /// new point by imposing the minimal-rotation condition.
  if(T.size() != R.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": T.size()=" << T.size() << " != R.size()=" << R.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  const unsigned int Size=T.size();
  const Quaternion z(0,0,0,1);
  vector<double> gammaover2dot(Size);
  vector<Quaternion> Rreturn(R);
  for(unsigned int iteration=0; iteration<NIterations; ++iteration) {
    cout << "\t\tIteration " << iteration << endl;
    const vector<Quaternion> Rdot = GWFrames::QuaternionDerivative(Rreturn, T);
    for(unsigned int i=0; i<Size; ++i) {
      gammaover2dot[i] = ( Rreturn[i].inverse() * Rdot[i] * z )[0];
    }
    const vector<double> gammaover2 = GWFrames::ScalarIntegral(gammaover2dot, T);
    for(unsigned int i=0; i<Size; ++i) {
      Rreturn[i] = Rreturn[i] * (gammaover2[i]*z).exp();
    }
  }
  cout << "\tFinished" << endl;
  return Rreturn;
}

/// Input frame with prescribed rate of rotation about Z axis.
std::vector<Quaternion> GWFrames::PrescribedRotation(const std::vector<double>& RotationRateAboutZ,
						     const std::vector<Quaternion>& R, const std::vector<double>& T, const unsigned int NIterations) {
  /// 
  /// \param RotationRateAboutZ Vector of rotation rates about the new frame's Z axis.
  /// \param R Vector of rotors.
  /// \param T Vector of corresponding time steps.
  /// \param NIterations Number of refinements [default: 5]
  /// 
  /// This function returns a copy of the input R, which takes the z
  /// axis to the same point as R, but adjusts the rotation about that
  /// new point by imposing the minimal-rotation condition, and then
  /// including an additional rotation about the new Z axis to agree
  /// with the given rotation rate.
  
  if(T.size() != R.size() || T.size() != RotationRateAboutZ.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": T.size()=" << T.size() << " != R.size()=" << R.size() << " != RotationRateAboutZ.size()=" << RotationRateAboutZ.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  const unsigned int Size=T.size();
  const Quaternion z(0,0,0,1);
  vector<double> gammaover2dot(Size);
  vector<Quaternion> Rreturn(R);
  
  // Repeat the process a few times to refine the result
  for(unsigned int iteration=0; iteration<NIterations; ++iteration) {
    cout << "\t\tIteration " << iteration << endl;
    const vector<Quaternion> Rdot = GWFrames::QuaternionDerivative(Rreturn, T);
    
    // Calculate \dot{\gamma}/2 at each instant in time
    for(unsigned int i=0; i<Size; ++i) {
      gammaover2dot[i] = 2*RotationRateAboutZ[i] + ( Rreturn[i].inverse() * Rdot[i] * z )[0];
    }
    
    // Integrate to find \gamma/2 as a function of time
    const vector<double> gammaover2 = GWFrames::ScalarIntegral(gammaover2dot, T);
    
    // Now pre-multiply the input rotor by exp(\gamma * z / 2) at each instant
    for(unsigned int i=0; i<Size; ++i) {
      Rreturn[i] = Rreturn[i] * (gammaover2[i]*z).exp();
    }
  }
  
  cout << "\tFinished" << endl;
  return Rreturn;
}

/// Construct frame given the X and Y basis vectors of that frame.
std::vector<Quaternion> GWFrames::FrameFromXY(const std::vector<Quaternion>& X, const std::vector<Quaternion>& Y) {
  ///
  /// \param X Vector of Quaternions
  /// \param Y Vector of Quaternions
  ///
  /// The input parameters are Quaternions, assumed to be pure unit
  /// vectors, representing the X and Y basis vectors of the frame at
  /// each instant of time.  The returned vector of rotors will rotate
  /// the stationary frame's (x,y,z) vectors into the new frame's
  /// (X,Y,Z) vectors.
  if(X.size() != Y.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": X.size()=" << X.size() << " != Y.size()=" << Y.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  const unsigned int Size=X.size();
  const Quaternion x(0,1,0,0);
  const Quaternion y(0,0,1,0);
  const Quaternion z(0,0,0,1);
  vector<Quaternion> R(Size);
  for(unsigned int k=0; k<Size; ++k) {
    const Quaternion Ra = sqrtOfRotor(-X[k]*x);
    const double beta = std::atan2(GWFrames::dot(Ra*z*Ra.inverse(), Y[k]),
				   GWFrames::dot(Ra*y*Ra.inverse(), Y[k]));
    R[k] = Ra * GWFrames::exp((beta/2.0)*x);
  }
  return R;
}

/// Construct minimal-rotation frame from Z basis vector of that frame.
std::vector<Quaternion> GWFrames::FrameFromZ(const std::vector<Quaternion>& Z, const std::vector<double>& T, const unsigned int NIterations) {
  ///
  /// \param Z Vector of Quaternions
  /// \param T Vector of corresponding times
  /// \param NIterations Number of refinements [default: 5]
  /// 
  /// The input vector of Quaternions, assumed to be pure unit
  /// vectors, represent the Z basis vectors of the frame at each
  /// instant of time.  The returned vector of rotors will rotate the
  /// stationary frame's (x,y,z) vectors into the new frame's (X,Y,Z)
  /// vectors.  The X and Y vectors are deduced by imposing the
  /// minimal-rotation condition.  Note that this leaves an unfixed
  /// initial rotation about z.
  if(Z.size() != T.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": Z.size()=" << Z.size() << " != T.size()=" << T.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  const unsigned int Size=Z.size();
  const Quaternion z(0,0,0,1);
  vector<Quaternion> R(Size);
  cout << "\tComputing basic frame" << endl;
  for(unsigned int k=0; k<Size; ++k) {
    R[k] = GWFrames::sqrt(-Z[k]*z);
  }
  cout << "\tRemoving rotation:" << endl;
  return GWFrames::MinimalRotation(GWFrames::UnflipRotors(R), T, NIterations);
}

/// Construct minimal-rotation frame from Z basis vector of that frame.
std::vector<Quaternion> GWFrames::FrameFromPrescribedRotation(const std::vector<Quaternion>& omega, const std::vector<double>& T, const unsigned int NIterations) {
  ///
  /// \param omega Vector of Quaternions
  /// \param T Vector of corresponding times
  /// \param NIterations Number of refinements [default: 5]
  /// 
  /// The input vector of Quaternions represent the angular-velocity
  /// vector (omega) of the frame at each instant of time.  The
  /// returned vector of rotors will rotate the stationary frame's
  /// (x,y,z) vectors into the new frame's (X,Y,Z) vectors, where Z is
  /// parallel to omega, and the X and Y vectors are deduced by
  /// enforcing the condition that the instantaneous rotation of the
  /// frame about Z is |omega|.  Note that this leaves an unfixed
  /// initial rotation in the X--Y plane.
  if(omega.size() != T.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": omega.size()=" << omega.size() << " != T.size()=" << T.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  const unsigned int Size=omega.size();
  const Quaternion z(0,0,0,1);
  vector<Quaternion> R(Size);
  cout << "\tComputing basic frame" << endl;
  for(unsigned int k=0; k<Size; ++k) {
    R[k] = GWFrames::sqrt(-omega[k].normalized()*z);
  }
  cout << "\tRemoving rotation:" << endl;
  return GWFrames::PrescribedRotation(abs(omega), GWFrames::UnflipRotors(R), T, NIterations);
}


#ifndef DOXYGEN
typedef struct {
  gsl_interp_accel* accX;
  gsl_spline* splineX;
  gsl_interp_accel* accY;
  gsl_spline* splineY;
  gsl_interp_accel* accZ;
  gsl_spline* splineZ;
} ParameterSet;
int FrameFromAngularVelocity_RHS(double t, const double ri[], double drdt[], void* mu) {
  // Interpolate and unpack some values
  ParameterSet* params = (ParameterSet*) mu;
  const double OmegaX = gsl_spline_eval(params->splineX, t, params->accX);
  const double OmegaY = gsl_spline_eval(params->splineY, t, params->accY);
  const double OmegaZ = gsl_spline_eval(params->splineZ, t, params->accZ);
  drdt[0] = OmegaX/2.0;
  drdt[1] = OmegaY/2.0;
  drdt[2] = OmegaZ/2.0;
  const double absquatlogR = std::sqrt(ri[0]*ri[0]+ri[1]*ri[1]+ri[2]*ri[2]);
  const double absOmega = std::sqrt(OmegaX*OmegaX+OmegaY*OmegaY+OmegaZ*OmegaZ);
  if(absquatlogR < Quaternion_Epsilon * absOmega) { // If the matrix is really close to the identity, return
    return GSL_SUCCESS;
  }
  if(std::abs(std::sin(absquatlogR)) < Quaternion_Epsilon) { // If the matrix is really close to singular, it's equivalent to the identity, so return
    return GSL_SUCCESS;
  }
  
  const Quaternion OmegaOver2(0., drdt[0], drdt[1], drdt[2]);
  const Quaternion rQ(0., ri[0], ri[1], ri[2]);
  const Quaternion rHat = rQ/absquatlogR;
  const Quaternion drdtQ = (OmegaOver2 - rHat*(rHat.dot(OmegaOver2)))*(absquatlogR/std::tan(absquatlogR)) + rHat*(rHat.dot(OmegaOver2)) + OmegaOver2.cross(rQ);
  drdt[0] = drdtQ[1];
  drdt[1] = drdtQ[2];
  drdt[2] = drdtQ[3];
  return GSL_SUCCESS;
}
#endif // DOXYGEN
std::vector<Quaternion> GWFrames::FrameFromAngularVelocity(const std::vector<Quaternion>& Omega, const std::vector<double>& T) {
  /// 
  /// \param Omega Vector of Quaternions.
  /// \param T Vector of corresponding times.
  /// 
  /// Note that each element of Omega should be a pure-vector
  /// Quaternion, corresponding to the angular-velocity vector at the
  /// instant of time.
  const vector<double> OmegaX = GWFrames::Component1(Omega);
  const vector<double> OmegaY = GWFrames::Component2(Omega);
  const vector<double> OmegaZ = GWFrames::Component3(Omega);
  
  // Set up the spline to interpolate Omega
  gsl_interp_accel* accX = gsl_interp_accel_alloc();
  gsl_spline* splineX = gsl_spline_alloc(gsl_interp_cspline, OmegaX.size());
  gsl_spline_init(splineX, &T[0], &OmegaX[0], OmegaX.size());
  gsl_interp_accel* accY = gsl_interp_accel_alloc();
  gsl_spline* splineY = gsl_spline_alloc(gsl_interp_cspline, OmegaY.size());
  gsl_spline_init(splineY, &T[0], &OmegaY[0], OmegaY.size());
  gsl_interp_accel* accZ = gsl_interp_accel_alloc();
  gsl_spline* splineZ = gsl_spline_alloc(gsl_interp_cspline, OmegaZ.size());
  gsl_spline_init(splineZ, &T[0], &OmegaZ[0], OmegaZ.size());
  
  // Set up the integrator
  const double hstart = (T[1]-T[0])/10.;
  const double epsabs = 1.e-9;
  const double epsrel = 1.e-9;
  ParameterSet params;
  params.accX = accX;
  params.splineX = splineX;
  params.accY = accY;
  params.splineY = splineY;
  params.accZ = accZ;
  params.splineZ = splineZ;
  gsl_odeiv2_system sys = {FrameFromAngularVelocity_RHS, NULL, 3, (void *) &params};
  // gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rkf45, hstart, epsabs, epsrel);
  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, hstart, epsabs, epsrel);
  double t = T[0];
  double r[3] = {0.0, 0.0, 0.0};
  
  // Run the integration
  vector<vector<double> > rs;
  rs.push_back(vector<double>(r, r+3));
  for(unsigned int i=1; i<T.size(); ++i) {
    int status = gsl_odeiv2_driver_apply(d, &t, T[i], r);
    if(status != GSL_SUCCESS) {
      gsl_spline_free(splineX);
      gsl_interp_accel_free(accX);
      gsl_spline_free(splineY);
      gsl_interp_accel_free(accY);
      gsl_spline_free(splineZ);
      gsl_interp_accel_free(accZ);
      gsl_odeiv2_driver_free(d);
      // cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": gsl_odeiv2_driver_apply returned an error; return value=" << status << "\n\n";
      throw(GWFrames_FailedGSLCall);
    }
    rs.push_back(vector<double>(r, r+3));
    // Reset the value of r if necessary
    const double absquatlogR = std::sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
    if(absquatlogR>M_PI/2.) {
      r[0] = (absquatlogR-M_PI)*r[0]/absquatlogR;
      r[1] = (absquatlogR-M_PI)*r[1]/absquatlogR;
      r[2] = (absquatlogR-M_PI)*r[2]/absquatlogR;
      // cerr << "Flipping r at t=" << t << endl;
    }
  }
  
  // Free the gsl storage
  gsl_spline_free(splineX);
  gsl_interp_accel_free(accX);
  gsl_spline_free(splineY);
  gsl_interp_accel_free(accY);
  gsl_spline_free(splineZ);
  gsl_interp_accel_free(accZ);
  gsl_odeiv2_driver_free(d);
  
  return UnflipRotors(exp(Quaternions(rs)));
}


/// Remove sign-ambiguity of rotors.
std::vector<Quaternion> GWFrames::UnflipRotors(const std::vector<Quaternion>& R, const double discont) {
  ///
  /// \param R Vector of rotors
  /// \param discont Acceptable discontinuity [default: sqrt(2)]
  /// 
  /// Because of the two-sided nature of quaternion rotations, the
  /// sign of a rotor may be undetermined in many cases.
  /// Discontinuous flips in that sign for rotor-valued functions of
  /// time can cause significant problems.  This function removes
  /// those flips by ensuring that the output rotors at successive
  /// instants are within 'discont' of each other.
  vector<Quaternion> Q(R.size());
  Q[0] = R[0];
  for(unsigned int i=1; i<R.size(); ++i) {
    if((R[i]-Q[i-1]).abs() > discont) {
      Q[i] = -R[i];
    } else {
      Q[i] = R[i];
    }
  }
  return Q;
}

/// Difference between frame rotors
std::vector<Quaternion> GWFrames::RDelta(const std::vector<Quaternion>& R1, const std::vector<Quaternion>& R2, const unsigned int IndexOfFiducialTime) {
  ///
  /// \param R1 Vector of rotors
  /// \param R2 Vector of rotors
  /// \param IndexOfFiducialTime Integer index of time at which
  ///        difference is set to zero [default: 0]
  if(R1.size() != R2.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": R1.size()=" << R1.size() << " != R2.size()=" << R2.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  if(R1.size()==0) {
    return vector<Quaternion>(0);
  }
  if(R1.size()<=IndexOfFiducialTime) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": R1.size()=" << R1.size() << " <= IndexOfFiducialTime=" << IndexOfFiducialTime << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  const unsigned int Size=R1.size();
  vector<Quaternion> Rd(Size);
  const Quaternion Offset = R1[IndexOfFiducialTime].inverse() * R2[IndexOfFiducialTime];
  for(unsigned int i=0; i<Size; ++i) {
    Rd[i] = R1[i] * Offset * R2[i].inverse();
  }
  return Rd;
}

/// Squad interpolation of Quaternion time series.
std::vector<Quaternion> GWFrames::Squad(const std::vector<Quaternion>& RIn, const std::vector<double>& tIn, const std::vector<double>& tOut) {
  ///
  /// \param RIn Vector of rotors
  /// \param tIn Vector of corresponding times
  /// \param tOut Vector of times to which RIn will be interpolated
  /// 
  /// This function implements a version of cubic-spline interpolation
  /// designed for unit quaternions, which delivers more accurate,
  /// smooth, and physical rotations than other forms of
  /// interpolation.
  if(RIn.size() != tIn.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": RIn.size()=" << RIn.size() << " != tIn.size()=" << tIn.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> ROut(tOut.size());
  unsigned int iIn = 0;
  unsigned int iOut = 0;
  while(iOut<tOut.size() && iIn<tIn.size() && tIn[tIn.size()-1]>=tOut[iOut]) {
    double Dtim1, Dti, Dtip1;
    Quaternion Qim1, Qi, Qip1, Qip2, Ai, Bip1;
    while(iIn+1<tIn.size() && tIn[iIn+1]<tOut[iOut]) {
      iIn += 1;
    }
    if(iIn+1==tIn.size()) {
      cerr << "\n" << __FILE__ << ":" << __LINE__ << ": Time " << tOut[iOut] << " is beyond the end of the input data (time " << tIn.back() << ")." << endl;
      throw(GWFrames_CannotExtrapolateQuaternions);
    }
    if(iIn==0) {
      Dti = tIn[iIn+1]-tIn[iIn];
      Dtim1 = Dti;
      Dtip1 = tIn[iIn+2]-tIn[iIn+1];
      Qim1 = RIn[iIn]*RIn[iIn+1].inverse()*RIn[iIn];
      Qi = RIn[iIn];
      Qip1 = RIn[iIn+1];
      Qip2 = RIn[iIn+2];
    } else if(iIn+2==tIn.size()) {
      Dtim1 = tIn[iIn]-tIn[iIn-1];
      Dti = tIn[iIn+1]-tIn[iIn];
      Dtip1 = Dti;
      Qim1 = RIn[iIn-1];
      Qi = RIn[iIn];
      Qip1 = RIn[iIn+1];
      Qip2 = RIn[iIn+1]*RIn[iIn].inverse()*RIn[iIn+1];
    } else {
      Dtim1 = tIn[iIn]-tIn[iIn-1];
      Dti = tIn[iIn+1]-tIn[iIn];
      Dtip1 = tIn[iIn+2]-tIn[iIn+1];
      Qim1 = RIn[iIn-1];
      Qi = RIn[iIn];
      Qip1 = RIn[iIn+1];
      Qip2 = RIn[iIn+2];
    }
    Ai = Qi * GWFrames::exp((
		       GWFrames::log(Qi.inverse()*Qip1)
		       +(Dti/Dtim1)*GWFrames::log(Qim1.inverse()*Qi)
		       -2*GWFrames::log(Qi.inverse()*Qip1)
		       )*0.25);
    Bip1 = Qip1 * GWFrames::exp((
			   (Dti/Dtip1)*GWFrames::log(Qip1.inverse()*Qip2)
			   +GWFrames::log(Qi.inverse()*Qip1)
			   -2*GWFrames::log(Qi.inverse()*Qip1)
			   )*-0.25);
    while(iOut<tOut.size() && tOut[iOut]<=tIn[iIn+1]) {
      const double taui = (tOut[iOut]-tIn[iIn]) / Dti;
      ROut[iOut] = GWFrames::Slerp(2*taui*(1-taui),
				   GWFrames::Slerp(taui, Qi, Qip1),
				   GWFrames::Slerp(taui, Ai, Bip1));
      iOut += 1;
    }
    iIn += 1;
  }
  return ROut;
}

#ifndef DOXYGEN
std::vector<Quaternion> GWFrames::operator+(const double a, const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]+a;
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator-(const double a, const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = a-Q[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator*(const double a, const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = a*Q[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator/(const double a, const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = a/Q[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator+(const std::vector<double>& a, const Quaternion& Q) {
  vector<Quaternion> R(a.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q+a[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator-(const std::vector<double>& a, const Quaternion& Q) {
  vector<Quaternion> R(a.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = a[i]-Q;
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator*(const std::vector<double>& a, const Quaternion& Q) {
  vector<Quaternion> R(a.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = a[i]*Q;
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator/(const std::vector<double>& a, const Quaternion& Q) {
  vector<Quaternion> R(a.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = a[i]/Q;
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator+(const std::vector<double>& a, const std::vector<Quaternion>& Q) {
  if(a.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]+a[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator-(const std::vector<double>& a, const std::vector<Quaternion>& Q) {
  if(a.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = a[i]-Q[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator*(const std::vector<double>& a, const std::vector<Quaternion>& Q) {
  if(a.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = a[i]*Q[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator/(const std::vector<double>& a, const std::vector<Quaternion>& Q) {
  if(a.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = a[i]/Q[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator+(const Quaternion& a, const std::vector<Quaternion>& Q) {
  std::vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]+a;
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator-(const Quaternion& a, const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = a-Q[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator*(const Quaternion& a, const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = a*Q[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator/(const Quaternion& a, const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = a/Q[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator+(const std::vector<Quaternion>& Q, const Quaternion& a) {
  std::vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]+a;
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator-(const std::vector<Quaternion>& Q, const Quaternion& a) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]-a;
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator*(const std::vector<Quaternion>& Q, const Quaternion& a) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]*a;
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator/(const std::vector<Quaternion>& Q, const Quaternion& a) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]/a;
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator+(const std::vector<Quaternion>& Q, const double a) {
  std::vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]+a;
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator-(const std::vector<Quaternion>& Q, const double a) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]-a;
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator*(const std::vector<Quaternion>& Q, const double a) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]*a;
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator/(const std::vector<Quaternion>& Q, const double a) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]/a;
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator+(const Quaternion& Q, const std::vector<double>& a) {
  vector<Quaternion> R(a.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q+a[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator-(const Quaternion& Q, const std::vector<double>& a) {
  vector<Quaternion> R(a.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q-a[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator*(const Quaternion& Q, const std::vector<double>& a) {
  vector<Quaternion> R(a.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q*a[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator/(const Quaternion& Q, const std::vector<double>& a) {
  vector<Quaternion> R(a.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q/a[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator+(const std::vector<Quaternion>& Q, const std::vector<double>& a) {
  if(a.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]+a[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator-(const std::vector<Quaternion>& Q, const std::vector<double>& a) {
  if(a.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]-a[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator*(const std::vector<Quaternion>& Q, const std::vector<double>& a) {
  if(a.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]*a[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator/(const std::vector<Quaternion>& Q, const std::vector<double>& a) {
  if(a.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]/a[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator+(const std::vector<Quaternion>& Q, const std::vector<Quaternion>& a) {
  if(a.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]+a[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator-(const std::vector<Quaternion>& Q, const std::vector<Quaternion>& a) {
  if(a.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]-a[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator*(const std::vector<Quaternion>& Q, const std::vector<Quaternion>& a) {
  if(a.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]*a[i];
  }
  return R;
}
std::vector<Quaternion> GWFrames::operator/(const std::vector<Quaternion>& Q, const std::vector<Quaternion>& a) {
  if(a.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": a.size()=" << a.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i]/a[i];
  }
  return R;
}


/////////////////////
// Array functions //
/////////////////////

std::vector<Quaternion> GWFrames::pow(const std::vector<Quaternion>& Q, const double t) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].pow(t);
  }
  return R;
}
std::vector<Quaternion> GWFrames::pow(const std::vector<Quaternion>& Q, const Quaternion& P) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].pow(P);
  }
  return R;
}
std::vector<Quaternion> GWFrames::pow(const Quaternion& Q, const std::vector<double>& t) {
  vector<Quaternion> R(t.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q.pow(t[i]);
  }
  return R;
}
std::vector<Quaternion> GWFrames::pow(const Quaternion& Q, const std::vector<Quaternion>& P) {
  vector<Quaternion> R(P.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q.pow(P[i]);
  }
  return R;
}
std::vector<Quaternion> GWFrames::pow(const std::vector<Quaternion>& Q, const std::vector<double>& t) {
  if(t.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": t.size()=" << t.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].pow(t[i]);
  }
  return R;
}
std::vector<Quaternion> GWFrames::pow(const std::vector<Quaternion>& Q, const std::vector<Quaternion>& P) {
  if(P.size() != Q.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": P.size()=" << P.size() << " != Q.size()=" << Q.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].pow(P[i]);
  }
  return R;
}
std::vector<double> GWFrames::abs(const std::vector<Quaternion>& Q) {
  vector<double> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].abs();
  }
  return R;
}
std::vector<Quaternion> GWFrames::log(const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].log();
  }
  return R;
}
std::vector<Quaternion> GWFrames::exp(const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].exp();
  }
  return R;
}
std::vector<Quaternion> GWFrames::sqrt(const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].sqrt();
  }
  return R;
}
std::vector<Quaternion> GWFrames::sqrtOfRotor(const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].sqrtOfRotor();
  }
  return R;
}
std::vector<double> GWFrames::angle(const std::vector<Quaternion>& Q) {
  vector<double> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].angle();
  }
  return R;
}
std::vector<Quaternion> GWFrames::inverse(const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].inverse();
  }
  return R;
}
std::vector<Quaternion> GWFrames::conjugate(const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].conjugate();
  }
  return R;
}
std::vector<Quaternion> GWFrames::normalized(const std::vector<Quaternion>& Q) {
  vector<Quaternion> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].normalized();
  }
  return R;
}
std::vector<double> GWFrames::normsquared(const std::vector<Quaternion>& Q) {
  vector<double> R(Q.size());
  for(unsigned int i=0; i<R.size(); ++i) {
    R[i] = Q[i].normsquared();
  }
  return R;
}

std::vector<Quaternion> GWFrames::Quaternions(const std::vector<double>& vartheta, const std::vector<double>& varphi) {
  vector<Quaternion> Qs(vartheta.size());
  for(unsigned int i=0; i<Qs.size(); ++i) {
    Qs[i] = Quaternion(vartheta[i], varphi[i]);
  }
  return Qs;
}

std::vector<Quaternion> GWFrames::Quaternions(const std::vector<double>& alpha, const std::vector<double>& beta, const std::vector<double>& gamma) {
  vector<Quaternion> Qs(alpha.size());
  for(unsigned int i=0; i<Qs.size(); ++i) {
    Qs[i] = Quaternion(alpha[i], beta[i], gamma[i]);
  }
  return Qs;
}

std::vector<Quaternion> GWFrames::Quaternions(const std::vector<double>& w0, const std::vector<double>& x0, const std::vector<double>& y0, const std::vector<double>& z0) {
  vector<Quaternion> Qs(w0.size());
  for(unsigned int i=0; i<Qs.size(); ++i) {
    Qs[i] = Quaternion(w0[i], x0[i], y0[i], z0[i]);
  }
  return Qs;
}

std::vector<Quaternion> GWFrames::Quaternions(const std::vector<std::vector<double> >& q) {
  vector<Quaternion> Qs(q.size());
  for(unsigned int i=0; i<Qs.size(); ++i) {
    Qs[i] = Quaternion(q[i]);
  }
  return Qs;
}

std::vector<Quaternion> GWFrames::Quaternions(const std::vector<double>& angle, const std::vector<std::vector<double> >& axis) {
  vector<Quaternion> Qs(angle.size());
  for(unsigned int i=0; i<Qs.size(); ++i) {
    Qs[i] = Quaternion(angle[i], axis[i]);
  }
  return Qs;
}

std::vector<double> GWFrames::Component(const std::vector<GWFrames::Quaternion>& Q, const unsigned int i) {
  std::vector<double> v(Q.size());
  for(unsigned int j=0; j<Q.size(); ++j) {
    v[j] = Q[j][i];
  }
  return v;
}

std::vector<double> GWFrames::Component0(const std::vector<GWFrames::Quaternion>& Q) {
  return Component(Q, 0);
}

std::vector<double> GWFrames::Component1(const std::vector<GWFrames::Quaternion>& Q) {
  return Component(Q, 1);
}

std::vector<double> GWFrames::Component2(const std::vector<GWFrames::Quaternion>& Q) {
  return Component(Q, 2);
}

std::vector<double> GWFrames::Component3(const std::vector<GWFrames::Quaternion>& Q) {
  return Component(Q, 3);
}


#endif // DOXYGEN


///////////////////////////////////////////////
// Derivative and integral utility functions //
///////////////////////////////////////////////

/// Three-point finite-differencing of vector of Quaternions.
std::vector<Quaternion> GWFrames::QuaternionDerivative(const std::vector<Quaternion>& f, const std::vector<double>& t) {
  ///
  /// \param f Vector of Quaternions.
  /// \param t Vector of corresponding time steps.
  if(f.size() != t.size()) {
    cerr << "\n\n" << __FILE__ << ":" << __LINE__ << ": f.size()=" << f.size() << " != t.size()=" << t.size() << endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  if(f.size()<3) { cerr << "\n" << __FILE__ << ":" << __LINE__ << ": size=" << f.size() << endl; throw(GWFrames_NotEnoughPointsForDerivative); }
  vector<Quaternion> D(f.size());
  const unsigned int i1 = f.size()-1;
  double hprev = t[1]-t[0];
  { // Compute first point
    const double hnext = t[2]-t[1];
    D[0] = -((2*hprev+hnext)/(hprev*(hprev+hnext)))*f[0] + ((hnext+hprev)/(hnext*hprev))*f[1] - (hprev/(hnext*(hnext+hprev)))*f[2];
  }
  for(unsigned int i=1; i<i1; ++i) { // Compute intermediate points
    const double hnext = t[i+1]-t[i];
    /// Sundquist and Veronis, Tellus XXII (1970), 1
    D[i] = (f[i+1] - f[i-1]*SQR(hnext/hprev) - f[i]*(1-SQR(hnext/hprev))) / (hnext*(1+hnext/hprev));
    hprev = hnext;
  }
  { // Compute final point
    const double hnext = t[i1]  -t[i1-1];
    const double hprev = t[i1-1]-t[i1-2];
    D[i1] = (hnext/(hprev*(hprev+hnext)))*f[i1-2] - ((hnext+hprev)/(hnext*hprev))*f[i1-1] + ((hprev+2*hnext)/(hnext*(hnext+hprev)))*f[i1];
  }
  return D;
}