// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

#include "Scri.hpp"

#include <omp.h>

#include <algorithm>


// The following are for spinsfast
#ifndef DOXYGEN
namespace GWFrames {
  #ifndef restrict
  #ifdef __restrict
  #define restrict __restrict
  #endif
  #endif
  extern "C" {
    #include <stdlib.h>
    #include <stdio.h>
    #include <math.h>
    #include <complex.h>
    #include "fftw3.h"
    #include "alm.h"
    #include "wigner_d_halfpi.h"
    #include "spinsfast_forward.h"
    #include "spinsfast_backward.h"
  }
};
#endif // DOXYGEN

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "Utilities.hpp"
#include "Quaternions.hpp"
#include "SphericalHarmonics.hpp"
#include "Waveforms.hpp"
#include "Errors.hpp"

using GWFrames::Quaternion;
// using GWFrames::StereographicCoordinate;
// using GWFrames::MobiusTransform;
using GWFrames::FourVector;
using GWFrames::DataGrid;
using GWFrames::Modes;
using GWFrames::SliceOfScri;
using GWFrames::SliceModes;
using GWFrames::Scri;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::flush;
using std::endl;
using std::complex;


int ellmax(const int maxindex) {
  int l, m;
  GWFrames::ind_lm(maxindex, &l, &m, 0);
  return l;
}

const std::complex<double> complexi(0.0,1.0);
const std::complex<double> zero(0.0,0.0);
const Quaternion zHat(0,0,0,1);


//////////////
// DataGrid //
//////////////

DataGrid::DataGrid(const int Spin, const int N_theta, const int N_phi, const std::vector<std::complex<double> >& D) 
  : s(Spin), n_theta(N_theta), n_phi(N_phi), data(D)
{
  // Check that we have the right amount of data
  if(n_theta*n_phi != int(D.size())) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: (n_theta=" << n_theta << ")*(n_phi=" << n_phi << ") != (D.size()=" << D.size() << ")"
	      << "\n       Cannot multiply data of different sizes\n"
	      << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
}

DataGrid::DataGrid(Modes M, const int N_theta, const int N_phi)
  : s(M.Spin()), n_theta(std::max(N_theta, 2*M.EllMax()+1)), n_phi(std::max(N_phi, 2*M.EllMax()+1)), data(n_phi*n_theta, zero)
{
  spinsfast_salm2map(reinterpret_cast<fftw_complex*>(&M[0]),
		     reinterpret_cast<fftw_complex*>(&data[0]),
		     M.Spin(), n_theta, n_phi, M.EllMax());
}

DataGrid::DataGrid(const Modes& M, const GWFrames::ThreeVector& v, const int N_theta, const int N_phi)
  : s(M.Spin()), n_theta(std::max(N_theta, 2*M.EllMax()+1)), n_phi(std::max(N_phi, 2*M.EllMax()+1)), data(n_phi*n_theta)
{
  const double dtheta = M_PI/double(n_theta-1); // theta should return to M_PI
  const double dphi = 2*M_PI/double(n_phi); // phi should not return to 2*M_PI
  { int i=0;
    for(int i_theta=0; i_theta<n_theta; ++i_theta) {
      for(int i_phi=0; i_phi<n_phi; ++i_phi, ++i) {
	const Quaternion Rp(dtheta*i_theta, dphi*i_phi);
	const Quaternion R_b = GWFrames::Boost(-v, (Rp*zHat*Rp.conjugate()).vec());
	data[i] = M.EvaluateAtPoint(R_b*Rp);
      }
    }
  }
}

/// Constructor on boosted grid by means of functor
DataGrid::DataGrid(const int Spin, const int N_theta, const int N_phi, const GWFrames::ThreeVector& v, const ScriFunctor& f)
  : s(Spin), n_theta(N_theta), n_phi(N_phi), data(n_phi*n_theta)
{
  /// \param Spin Integer spin weight
  /// \param N_theta Number of points in output grid in theta
  /// \param N_phi Number of points in output grid in phi
  /// \param v Three-vector velocity of boosted frame relative to current frame
  /// \param f Functor operating on a Quaternion object
  /// 
  /// The functor takes a Quaternion argument, which describes the
  /// location and orientation of the point to be evaluated.  In
  /// particular, the rotor takes the \f$\hat{z}\f$ vector into the
  /// point at which the field is to be measured, and takes \f$\hat{x}
  /// + i \hat{y}\f$ into the \f$m\f$ vector (within normalization)
  /// needed for spin-weighted fields.
  
  const double dtheta = M_PI/double(n_theta-1); // theta should return to M_PI
  const double dphi = 2*M_PI/double(n_phi); // phi should not return to 2*M_PI
  { int i=0;
    for(int i_theta=0; i_theta<n_theta; ++i_theta) {
      for(int i_phi=0; i_phi<n_phi; ++i_phi, ++i) {
	const Quaternion Rp(dtheta*i_theta, dphi*i_phi);
	const Quaternion R_b = GWFrames::Boost(-v, (Rp*zHat*Rp.conjugate()).vec());
	data[i] = f(R_b*Rp);
      }
    }
  }
}

DataGrid DataGrid::operator*(const DataGrid& A) const {
  const DataGrid& B=*this;
  
  // Check that we have the same amounts of data
  if(A.n_theta != B.n_theta) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: (A.n_theta=" << A.n_theta << ") != (B.n_theta=" << B.n_theta << ")"
	      << "\n       Cannot multiply data of different sizes\n"
	      << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  if(A.n_phi != B.n_phi) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: (A.n_phi=" << A.n_phi << ") != (B.n_phi=" << B.n_phi << ")"
	      << "\n       Cannot multiply data of different sizes\n"
	      << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  
  // Do the work
  DataGrid C(B);
  for(unsigned int i=0; i<data.size(); ++i) {
    C.data[i] *= A.data[i];
  }
  C.s = s + A.s;
  
  return C;
}

DataGrid DataGrid::operator/(const DataGrid& A) const {
  const DataGrid& B=*this;
  
  // Check that we have the same amounts of data
  if(A.n_theta != B.n_theta) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: (A.n_theta=" << A.n_theta << ") != (B.n_theta=" << B.n_theta << ")"
	      << "\n       Cannot divide data of different sizes\n"
	      << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  if(A.n_phi != B.n_phi) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: (A.n_phi=" << A.n_phi << ") != (B.n_phi=" << B.n_phi << ")"
	      << "\n       Cannot divide data of different sizes\n"
	      << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  
  // Do the work
  DataGrid C(B);
  for(unsigned int i=0; i<data.size(); ++i) {
    C.data[i] /= A.data[i];
  }
  C.s = s - A.s;
  
  return C;
}

DataGrid DataGrid::operator+(const DataGrid& A) const {
  const DataGrid& B=*this;
  
  // Check that we have the same amounts of data
  if(A.n_theta != B.n_theta) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: (A.n_theta=" << A.n_theta << ") != (B.n_theta=" << B.n_theta << ")"
	      << "\n       Cannot add data of different sizes\n"
	      << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  if(A.n_phi != B.n_phi) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: (A.n_phi=" << A.n_phi << ") != (B.n_phi=" << B.n_phi << ")"
	      << "\n       Cannot add data of different sizes\n"
	      << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  
  // Do the work
  DataGrid C(B);
  for(unsigned int i=0; i<data.size(); ++i) {
    C.data[i] += A.data[i];
  }
  
  return C;
}

DataGrid DataGrid::operator-(const DataGrid& A) const {
  const DataGrid& B=*this;
  
  // Check that we have the same amounts of data
  if(A.n_theta != B.n_theta) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: (A.n_theta=" << A.n_theta << ") != (B.n_theta=" << B.n_theta << ")"
	      << "\n       Cannot subtract data of different sizes\n"
	      << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  if(A.n_phi != B.n_phi) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: (A.n_phi=" << A.n_phi << ") != (B.n_phi=" << B.n_phi << ")"
	      << "\n       Cannot subtract data of different sizes\n"
	      << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  
  // Do the work
  DataGrid C(B);
  for(unsigned int i=0; i<data.size(); ++i) {
    C.data[i] -= A.data[i];
  }
  
  return C;
}

DataGrid DataGrid::pow(const int p) const {
  DataGrid c(*this);
  const int N = c.N_theta()*c.N_phi();
  for(int i=0; i<N; ++i) {
    c[i] = std::pow(c[i], p);
  }
  return c;
}

DataGrid GWFrames::operator*(const double& a, const DataGrid& b) {
  DataGrid c(b);
  const int N = b.N_theta()*b.N_phi();
  for(int i=0; i<N; ++i) {
    c[i] *= a;
  }
  return c;
}

DataGrid GWFrames::operator/(const double& a, const DataGrid& b) {
  DataGrid c(b);
  const int N = b.N_theta()*b.N_phi();
  for(int i=0; i<N; ++i) {
    c[i] = a/c[i];
  }
  return c;
}

DataGrid GWFrames::operator-(const double& a, const DataGrid& b) {
  DataGrid c(b);
  const int N = b.N_theta()*b.N_phi();
  for(int i=0; i<N; ++i) {
    c[i] = a-c[i];
  }
  return c;
}



// /// Construct a grid with the conformal factor at each point
// DataGrid GWFrames::ConformalFactorGrid(const GWFrames::MobiusTransform& abcd, const int n_theta, const int n_phi) {
//   std::vector<std::complex<double> > D(n_theta*n_phi);
//   const double dtheta = M_PI/double(n_theta-1); // theta should return to M_PI
//   const double dphi = 2*M_PI/double(n_phi); // phi should not return to 2*M_PI
//   int i=0;
//   for(int i_theta=0; i_theta<n_theta; ++i_theta) {
//     for(int i_phi=0; i_phi<n_phi; ++i_phi) {
//       D[i] = BoostConformalFactor(StereographicCoordinateFromAngles(dtheta*i_theta, dphi*i_phi), abcd);
//       ++i;
//     }
//   }
//   return DataGrid(0, n_theta, n_phi, D);
// }

// /// Construct a grid with the conformal factor at each point
// DataGrid GWFrames::ConformalFactorGrid(const GWFrames::ThreeVector& v, const int n_theta, const int n_phi) {
//   return GWFrames::ConformalFactorGrid(MobiusComponentsOfBoost(v), n_theta, n_phi);
// }

/// Construct a grid with the conformal factor at each point
DataGrid GWFrames::ConformalFactorGrid(const GWFrames::ThreeVector& v, const int n_theta, const int n_phi) {
  std::vector<std::complex<double> > D(n_theta*n_phi);
  const double dtheta = M_PI/double(n_theta-1); // theta should return to M_PI
  const double dphi = 2*M_PI/double(n_phi); // phi should not return to 2*M_PI
  const double gamma = 1.0/std::sqrt(1-v[0]*v[0]-v[1]*v[1]-v[2]*v[2]);
  { int i=0;
    for(int i_theta=0; i_theta<n_theta; ++i_theta) {
      for(int i_phi=0; i_phi<n_phi; ++i_phi, ++i) {
	const double theta = dtheta*i_theta;
	const double phi = dphi*i_phi;
	D[i] = 1.0/(gamma*(1
			   -v[0]*std::cos(phi)*std::sin(theta)
			   -v[1]*std::sin(phi)*std::sin(theta)
			   -v[2]*std::cos(theta)));
      }
    }
  }
  return DataGrid(0, n_theta, n_phi, D);
}

/// Construct a grid with the conformal factor at each point
DataGrid GWFrames::InverseConformalFactorGrid(const GWFrames::ThreeVector& v, const int n_theta, const int n_phi) {
  std::vector<std::complex<double> > D(n_theta*n_phi);
  const double dtheta = M_PI/double(n_theta-1); // theta should return to M_PI
  const double dphi = 2*M_PI/double(n_phi); // phi should not return to 2*M_PI
  const double gamma = 1.0/std::sqrt(1-v[0]*v[0]-v[1]*v[1]-v[2]*v[2]);
  { int i=0;
    for(int i_theta=0; i_theta<n_theta; ++i_theta) {
      for(int i_phi=0; i_phi<n_phi; ++i_phi, ++i) {
	const double theta = dtheta*i_theta;
	const double phi = dphi*i_phi;
	D[i] = gamma*(1
		      -v[0]*std::cos(phi)*std::sin(theta)
		      -v[1]*std::sin(phi)*std::sin(theta)
		      -v[2]*std::cos(theta));
      }
    }
  }
  return DataGrid(0, n_theta, n_phi, D);
}

class InverseConformalFactorFunctor : public GWFrames::ScriFunctor {
private:
  double gamma;
  Quaternion v;
public:
  InverseConformalFactorFunctor(const GWFrames::ThreeVector& vi)
    : gamma(1.0/std::sqrt(1-vi[0]*vi[0]-vi[1]*vi[1]-vi[2]*vi[2])),
      v(0., vi[0], vi[1], vi[2]) { }
  virtual double operator()(const GWFrames::Quaternion& R) const {
    return gamma*( 1 - v.dot(R*zHat*R.conjugate()) );
  }
};
/// Construct a boosted grid with the conformal factor at each point
DataGrid GWFrames::InverseConformalFactorBoostedGrid(const GWFrames::ThreeVector& v, const int n_theta, const int n_phi) {
  const InverseConformalFactorFunctor K(v);
  return DataGrid(0, n_theta, n_phi, v, K);
}




///////////
// Modes //
///////////

Modes::Modes(const int spin, const std::vector<std::complex<double> >& Data)
  : s(spin), ellMax(0), data(Data)
{
  // Find the appropriate ellMax for this data length
  for(; ellMax<=ellMax_GWFrames; ++ellMax) {
    if(N_lm(ellMax)==int(Data.size())) {
      break;
    }
  }
  
  // Make sure Data doesn't have more ell modes than ellMax_GWFrames
  if(ellMax>=ellMax_GWFrames) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: Input data has length " << Data.size() << "."
	      << "\n       This is not a recognized length for spherical-harmonic data.\n"
	      << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
}

Modes::Modes(DataGrid D)
  : s(D.Spin()), ellMax(std::min((D.N_theta()-1)/2, (D.N_phi()-1)/2)), data(N_lm(ellMax))
{
  spinsfast_map2salm(reinterpret_cast<fftw_complex*>(&D[0]),
		     reinterpret_cast<fftw_complex*>(&data[0]),
		     s, D.N_theta(), D.N_phi(), ellMax);
}

Modes Modes::bar() const {
  const Modes& A=*this;
  Modes B(A);
  for(unsigned int i=0; i<data.size(); ++i) {
    B.data[i] = std::conj(A.data[i]);
  }
  B.s = -A.s;
  return B;
}

Modes Modes::operator*(const Modes& M) const {
  const int L = std::max(EllMax(), M.EllMax());
  Modes A = Modes(DataGrid(*this,2*L+1,2*L+1) * DataGrid(M,2*L+1,2*L+1));
  A.s = s + M.s;
  return A;
}

Modes Modes::operator/(const Modes& M) const {
  const int L = std::max(EllMax(), M.EllMax());
  Modes A = Modes(DataGrid(*this,2*L+1,2*L+1) / DataGrid(M,2*L+1,2*L+1));
  A.s = s - M.s;
  return A;
}

Modes Modes::operator+(const Modes& M) const {
  // Check that both `Modes` have the same `s` values
  if(s != M.s) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: The second data set has spin M.s=" << M.s << ", which is not the same as this set's spin s=" << s << "."
	      << "\n       This addition is not defined.\n"
	      << std::endl;
    throw(GWFrames_BadWaveformInformation);
  }
  
  // Choose `A` to be the input with a larger data set
  const Modes& A = (data.size()>=M.data.size() ? *this : M);
  const Modes& B = (data.size()< M.data.size() ? *this : M);
  Modes C(A);
  
  // Just loop over the data elements in the smaller data set
  // (assuming higher modes are zero)
  for(unsigned int i=0; i<B.data.size(); ++i) {
    C.data[i] += B.data[i];
  }
  
  return C;
}

Modes Modes::edth() const {
  /// This operator is the one defined by Geroch et al. (1973).  It
  /// raises the spin weight of any field on the sphere by 1, while
  /// leaving the boost weight unchanged.
  /// 
  /// This operator is very similar to the basic Newman-Penrose edth
  /// operator, except that it preserves boost weights.  Its effect in
  /// this implementation is identical (up to a factor of
  /// \f$\sqrt{2}\f$) to the NP edth.  There is an additional term in
  /// the definition of the GHP operator, but its value is zero.  (It
  /// transforms nontrivially, though.)  In this context, we have
  /// `NPEdth() = sqrt(2)*GHPEdth()`.
  /// 
  /// The complex shear \f$\sigma\f$ has spin weight +2 and boost
  /// weight +1.  The radial coordinate \f$r\f$ has boost weight -1,
  /// and the derivative with respect to time \f$d/du\f$ has boost
  /// weight -1.  The asymptotic metric shear \f$r\, h\f$ has spin
  /// weight -2 and boost weight -1.  In particular, it seems that
  /// \f$r\, h = r^2\, \bar{\sigma}\f$.
  /// 
  /// The Newman-Penrose scalars \f$\Psi_i\f$ have spin weight and
  /// boost weight equal to \f$2-i\f$.  (E.g., \f$\Psi_4\f$ has \f$s =
  /// b = -2\f$.)  However, when these are multiplied by the
  /// appropriate factors of \f$r\f$ to find the leading-order terms,
  /// they acquire boost weights.  In particular, we need to multiply
  /// \f$\Psi_i\f$ by \f$r^{5-i}\f$ to get nonzero values at scri,
  /// which adds \f$i-5\f$ to the boost weight, so that the asymptotic
  /// NP scalars all have boost weight -3.
  
  const Modes& A=*this;
  Modes B(A);
  const int ellMin = std::abs(s+1);
  
  int i=0;
  for(int ell=0; ell<=ellMax; ++ell) {
    const double factor = (ell<ellMin ? 0.0 : std::sqrt((ell-s)*(ell+s+1.)/2.));
    for(int m=-ell; m<=ell; ++m) {
      B.data[i] *= factor;
      ++i;
    }
  }
  
  B.SetSpin(A.Spin()+1);
  return B;
}

Modes Modes::edthbar() const {
  /// This operator is the one defined by Geroch et al. (1973).  It
  /// lowers the spin weight of any field on the sphere by 1, while
  /// leaving the boost weight unchanged.
  /// 
  /// This operator is very similar to the basic Newman-Penrose edth
  /// operator, except that it preserves boost weights.  Its effect in
  /// this implementation is identical (up to a factor of
  /// \f$\sqrt{2}\f$) to the NP edth.  There is an additional term in
  /// the definition of the GHP operator, but its value is zero.  (It
  /// transforms nontrivially, though.)  In this context, we have
  /// `NPEdthBar() = sqrt(2)*GHPEdthBar()`.
  /// 
  /// The complex shear \f$\sigma\f$ has spin weight +2 and boost
  /// weight +1.  The radial coordinate \f$r\f$ has boost weight -1,
  /// and the derivative with respect to time \f$d/du\f$ has boost
  /// weight -1.  The asymptotic metric shear \f$r\, h\f$ has spin
  /// weight -2 and boost weight -1.  In particular, it seems that
  /// \f$r\, h = r^2\, \bar{\sigma}\f$.
  /// 
  /// The Newman-Penrose scalars \f$\Psi_i\f$ have spin weight and
  /// boost weight equal to \f$2-i\f$.  (E.g., \f$\Psi_4\f$ has \f$s =
  /// b = -2\f$.)  However, when these are multiplied by the
  /// appropriate factors of \f$r\f$ to find the leading-order terms,
  /// they acquire boost weights.  In particular, we need to multiply
  /// \f$\Psi_i\f$ by \f$r^{5-i}\f$ to get nonzero values at scri,
  /// which adds \f$i-5\f$ to the boost weight, so that the asymptotic
  /// NP scalars all have boost weight -3.
  
  const Modes& A=*this;
  Modes B(A);
  const int ellMin = std::abs(s-1);
  
  int i=0;
  for(int ell=0; ell<=ellMax; ++ell) {
    const double factor = (ell<ellMin ? 0.0 : -std::sqrt((ell+s)*(ell-s+1.)/2.));
    for(int m=-ell; m<=ell; ++m) {
      B.data[i] *= factor;
      ++i;
    }
  }
  
  B.SetSpin(A.Spin()-1);
  return B;
}


/// Evaluate Waveform at a particular sky location
std::complex<double> Modes::EvaluateAtPoint(const double vartheta, const double varphi) const {
  /// 
  /// \param vartheta Polar angle of detector
  /// \param varphi Azimuthal angle of detector
  /// 
  
  complex<double> d(0.0, 0.0);
  GWFrames::SWSH Y(s);
  Y.SetAngles(vartheta, varphi);
  int i=0;
  for(int ell=0; ell<=ellMax; ++ell) {
    for(int m=-ell; m<=ell; ++m) {
      d += data[i]*Y(ell,m);
      ++i;
    }
  }
  return d;
}

/// Evaluate Waveform at a particular sky location
std::complex<double> Modes::EvaluateAtPoint(const GWFrames::Quaternion& R) const {
  /// 
  /// \param R Quaternion giving point by rotation of \f$\hat{z}\f$
  /// 
  /// Note that the argument `R` might typically be thought of as the
  /// rotor taking the unit \f$z\f$ vector into a point \f$(\vartheta,
  /// \varphi)\f$.  However, more general arguments are possible; this
  /// feature is used to greatly simplify the process of constructing
  /// a `DataGrid` object from a `Modes` object with a boost.
  /// 
  
  complex<double> d(0.0, 0.0);
  GWFrames::SWSH Y(s, R);
  int i=0;
  for(int ell=0; ell<=ellMax; ++ell) {
    for(int m=-ell; m<=ell; ++m) {
      d += data[i]*Y(ell,m);
      ++i;
    }
  }
  return d;
}




/////////////////
// SliceOfScri //
/////////////////

/// Empty constructor with reserved storage
template <class D>
SliceOfScri<D>::SliceOfScri(const int size)
  : psi0(size), psi1(size), psi2(size), psi3(size), psi4(size), sigma(size), sigmadot(size)
{
  psi0.SetSpin(2);
  psi1.SetSpin(1);
  psi2.SetSpin(0);
  psi3.SetSpin(-1);
  psi4.SetSpin(-2);
  sigma.SetSpin(2);
  sigmadot.SetSpin(2);
}

// /// Constructor from data
// template <class D>
// SliceOfScri<D>::SliceOfScri(const D& Psi0, const D& Psi1, const D& Psi2,
// 			    const D& Psi3, const D& Psi4, const D& Sigma, const D& SigmaDot)
//   : psi0(Psi0), psi1(Psi1), psi2(Psi2), psi3(Psi3), psi4(Psi4), sigma(Sigma), sigmadot(SigmaDot)
// { // The following may be redundant, but it will never be wrong, and 
//   psi0.SetSpin(2);
//   psi1.SetSpin(1);
//   psi2.SetSpin(0);
//   psi3.SetSpin(-1);
//   psi4.SetSpin(-2);
//   sigma.SetSpin(2);
//   sigmadot.SetSpin(2);
// }

/// Constructor from ellMax
SliceModes::SliceModes(const int ellMax)
  : SliceOfScri<Modes>(GWFrames::lm_ind(ellMax, ellMax, ellMax)+1)
{
  psi0.SetEllMax(ellMax);
  psi1.SetEllMax(ellMax);
  psi2.SetEllMax(ellMax);
  psi3.SetEllMax(ellMax);
  psi4.SetEllMax(ellMax);
  sigma.SetEllMax(ellMax);
  sigmadot.SetEllMax(ellMax);
}

/// Find largest ell value in the data on this slice
int SliceModes::EllMax() const {
  return std::max(psi0.EllMax(),
		  std::max(psi1.EllMax(),
			   std::max(psi2.EllMax(),
				    std::max(psi3.EllMax(),
					     std::max(psi4.EllMax(),
						      std::max(sigma.EllMax(), sigmadot.EllMax())
						      )
					     )
				    )
			   )
		  );
}

/// Calculate the mass of the system from the four-momentum
double SliceModes::Mass() const {
  const FourVector p = FourMomentum();
  return std::sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]);
}

/// Calculate the four-momentum of the system from the supermomentum
GWFrames::FourVector SliceModes::FourMomentum() const {
  /// The (Bondi) four-momentum is given by the ell=0 and ell=1 modes
  /// of the supermomentum.
  const Modes Psi = SuperMomentum();
  const double sqrt4pi = std::sqrt(4*M_PI);
  const double sqrt3 = std::sqrt(3);
  const double sqrt6 = std::sqrt(6);
  FourVector p(4);
  p[0] = std::real(Psi[0])/sqrt4pi;
  p[1] = std::real((Psi[1]-Psi[3]))/(sqrt4pi*sqrt6);
  p[2] = std::real(complexi*(Psi[1]+Psi[3]))/(sqrt4pi*sqrt6);
  p[3] = std::real(Psi[2])/(sqrt4pi*sqrt3);
  return p;
}

/// Find the Moreschi supermomentum
Modes SliceModes::SuperMomentum() const {
  /// \f$\psi = \psi_2 + \sigma \dot{\bar{\sigma}} + \eth^2 \bar{\sigma}\f$
  return psi2 + sigma*sigmadot.bar() + sigma.bar().edth().edth();
}

// /// Execute a BMS transformation except for the supertranslation of points
// GWFrames::SliceGrid SliceModes::BMSTransformationOnSlice(const double u, const ThreeVector& v, const Modes& delta) const {
//   /// A full BMS transformation is only possible using information
//   /// from multiple slices due to the supertranslation moving points
//   /// "between slices".  This function simply transforms the data
//   /// within the slice by accounting for the change of grid at each
//   /// point, and the change of grid points themselves.  The returned
//   /// object is a `DataGrid` object, each point of which can then be
//   /// used to interpolate to the supertranslated time.
  
//   const int n_theta = 2*EllMax()+1;
//   const int n_phi = n_theta;
  
//   // Evaluate the functions we need on the appropriate grids
//   const DataGrid oneoverK_g = GWFrames::InverseConformalFactorGrid(v, n_theta, n_phi);
//   const DataGrid oneoverKcubed_g = oneoverK_g.pow(3);
//   const DataGrid delta_g(delta, n_theta, n_phi);
//   const DataGrid ethethdelta_g(delta.edth().edth(), n_theta, n_phi);
//   const DataGrid ethupok_g = DataGrid(Modes((u-delta_g)/oneoverK_g).edth(), n_theta, n_phi)*oneoverK_g; // (\eth u') / K
//   const DataGrid psi0_g(psi0, n_theta, n_phi);
//   const DataGrid psi1_g(psi1, n_theta, n_phi);
//   const DataGrid psi2_g(psi2, n_theta, n_phi);
//   const DataGrid psi3_g(psi3, n_theta, n_phi);
//   const DataGrid psi4_g(psi4, n_theta, n_phi);
//   const DataGrid sigma_g(sigma, n_theta, n_phi);
//   const DataGrid sigmadot_g(sigmadot, n_theta, n_phi);
  
//   // Construct new data accounting for changes of tetrad
//   const Modes psi4factor( oneoverKcubed_g*(psi4_g) );
//   const Modes psi3factor( oneoverKcubed_g*(psi3_g - ethupok_g*psi4_g) );
//   const Modes psi2factor( oneoverKcubed_g*(psi2_g - ethupok_g*(2*psi3_g - ethupok_g*psi4_g)) );
//   const Modes psi1factor( oneoverKcubed_g*(psi1_g - ethupok_g*(3*psi2_g - ethupok_g*(3*psi3_g - ethupok_g*psi4_g))) );
//   const Modes psi0factor( oneoverKcubed_g*(psi0_g - ethupok_g*(4*psi1_g - ethupok_g*(6*psi2_g - ethupok_g*(4*psi3_g - ethupok_g*psi4_g)))) );
//   const Modes sigmafactor( (sigma_g - ethethdelta_g)*oneoverK_g );
//   const Modes sigmadotfactor( sigmadot_g*(oneoverK_g.pow(2)) );
  
//   // Evaluate the primed quantities on the boosted grids
//   SliceGrid Grids;
//   Grids.psi0 = DataGrid(psi0factor, v, n_theta, n_phi);
//   Grids.psi1 = DataGrid(psi1factor, v, n_theta, n_phi);
//   Grids.psi2 = DataGrid(psi2factor, v, n_theta, n_phi);
//   Grids.psi3 = DataGrid(psi3factor, v, n_theta, n_phi);
//   Grids.psi4 = DataGrid(psi4factor, v, n_theta, n_phi);
//   Grids.sigma = DataGrid(sigmafactor, v, n_theta, n_phi);
//   Grids.sigmadot = DataGrid(sigmadotfactor, v, n_theta, n_phi);
  
//   return Grids;
// }

/// Execute a BMS transformation except for the supertranslation of points
GWFrames::SliceGrid SliceModes::BMSTransformationOnSlice(const double u, const ThreeVector& v, const Modes& delta) const {
  /// A full BMS transformation is only possible using information
  /// from multiple slices due to the supertranslation moving points
  /// "between slices".  This function simply transforms the data
  /// within the slice by accounting for the change of grid at each
  /// point, and the change of grid points themselves.  The returned
  /// object is a `DataGrid` object, each point of which can then be
  /// used to interpolate to the supertranslated time.
  
  const int n_theta = 2*EllMax()+1;
  const int n_phi = n_theta;
  
  // Evaluate the functions we need on the boosted (and appropriately spin-transformed) grid
  const DataGrid oneoverK_g = GWFrames::InverseConformalFactorBoostedGrid(v, n_theta, n_phi);
  const DataGrid oneoverKcubed_g = oneoverK_g.pow(3);
  const DataGrid ethethdelta_g(delta.edth().edth(), v, n_theta, n_phi);
  const DataGrid ethupok_g = DataGrid(Modes((u-DataGrid(delta,n_theta,n_phi))/GWFrames::InverseConformalFactorGrid(v, n_theta, n_phi)).edth(),
				      v, n_theta, n_phi)*oneoverK_g; // (\eth u') / K
  const DataGrid psi0_g(psi0, v, n_theta, n_phi);
  const DataGrid psi1_g(psi1, v, n_theta, n_phi);
  const DataGrid psi2_g(psi2, v, n_theta, n_phi);
  const DataGrid psi3_g(psi3, v, n_theta, n_phi);
  const DataGrid psi4_g(psi4, v, n_theta, n_phi);
  const DataGrid sigma_g(sigma, v, n_theta, n_phi);
  const DataGrid sigmadot_g(sigmadot, v, n_theta, n_phi);
  
  // Construct new data accounting for changes of tetrad
  SliceGrid Grids;
  Grids.psi4 = oneoverKcubed_g*(psi4_g);
  Grids.psi3 = oneoverKcubed_g*(psi3_g - ethupok_g*psi4_g);
  Grids.psi2 = oneoverKcubed_g*(psi2_g - ethupok_g*(2*psi3_g - ethupok_g*psi4_g));
  Grids.psi1 = oneoverKcubed_g*(psi1_g - ethupok_g*(3*psi2_g - ethupok_g*(3*psi3_g - ethupok_g*psi4_g)));
  Grids.psi0 = oneoverKcubed_g*(psi0_g - ethupok_g*(4*psi1_g - ethupok_g*(6*psi2_g - ethupok_g*(4*psi3_g - ethupok_g*psi4_g))));
  Grids.sigma = (sigma_g - ethethdelta_g)*oneoverK_g;
  Grids.sigmadot = sigmadot_g*(oneoverK_g.pow(2));
  
  return Grids;
}

// Explicit instantiations
template class SliceOfScri<DataGrid>;
template class SliceOfScri<Modes>;


//////////
// Scri //
//////////
Scri::Scri(const GWFrames::Waveform& psi0, const GWFrames::Waveform& psi1,
	   const GWFrames::Waveform& psi2, const GWFrames::Waveform& psi3,
	   const GWFrames::Waveform& psi4, const GWFrames::Waveform& sigma)
  : t(psi0.T()), slices(t.size(), SliceModes(psi0.EllMax()))
{
  // Check that everyone has the same NTimes().  This is a poor man's
  // way of making sure we have all the same times, and is of course
  // only necessary, not sufficient proof that the times are the same.
  if(psi0.NTimes()!=psi1.NTimes() || psi0.NTimes()!=psi2.NTimes() ||
     psi0.NTimes()!=psi3.NTimes() || psi0.NTimes()!=psi4.NTimes() || psi0.NTimes()!=sigma.NTimes()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: psi0.NTimes()=" << psi0.NTimes() << "  psi1.NTimes()=" << psi1.NTimes() << "  psi2.NTimes()=" << psi2.NTimes()
	      << "         psi3.NTimes()=" << psi3.NTimes() << "  psi4.NTimes()=" << psi4.NTimes() << "  sigma.NTimes()=" << sigma.NTimes()
	      << "\n       Cannot store data on different slices.\n"
	      << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  
  // Check that everyone has the same EllMax()
  const int ellMax = psi0.EllMax();
  if(ellMax!=psi1.EllMax() || ellMax!=psi2.EllMax() ||
     ellMax!=psi3.EllMax() || ellMax!=psi4.EllMax() || ellMax!=sigma.EllMax()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: psi0.EllMax()=" << ellMax        << "  psi1.EllMax()=" << psi1.EllMax() << "  psi2.EllMax()=" << psi2.EllMax()
	      << "         psi3.EllMax()=" << psi3.EllMax() << "  psi4.EllMax()=" << psi4.EllMax() << "  sigma.EllMax()=" << sigma.EllMax()
	      << "\n       Cannot store data with different EllMax values.\n"
	      << std::endl;
    throw(GWFrames_VectorSizeMismatch);
  }
  
  // Fill the new data
  const unsigned int ntimes = t.size();
  for(unsigned int i_t=0; i_t<ntimes; ++i_t) { // Fill everything but sigmadot
    GWFrames::SliceModes& slice = slices[i_t];
    int i_ellm=0;
    for(int ell=0; ell<=ellMax; ++ell) {
      for(int m=-ell; m<=ell; ++m, ++i_ellm) {
	slice[0][i_ellm] = psi0.Data(psi0.FindModeIndex(ell,m), i_t);
	slice[1][i_ellm] = psi1.Data(psi1.FindModeIndex(ell,m), i_t);
	slice[2][i_ellm] = psi2.Data(psi2.FindModeIndex(ell,m), i_t);
	slice[3][i_ellm] = psi3.Data(psi3.FindModeIndex(ell,m), i_t);
	slice[4][i_ellm] = psi4.Data(psi4.FindModeIndex(ell,m), i_t);
	slice[5][i_ellm] = sigma.Data(sigma.FindModeIndex(ell,m), i_t);
      }
    }
  }
  { int i_ellm=0;
    for(int ell=0; ell<=ellMax; ++ell) {
      for(int m=-ell; m<=ell; ++m, ++i_ellm) {
	const std::vector<std::complex<double> > sigmadot = sigma.DataDot(sigma.FindModeIndex(ell,m));
	for(unsigned int i_t=0; i_t<ntimes; ++i_t) { // Fill sigmadot
	  slices[i_t][6][i_ellm] = sigmadot[i_t];
	}
      }
    }
  }
}

/// Apply a (constant) BMS transformation to data on null infinity
SliceModes Scri::BMSTransformation(const double& uPrime, const ThreeVector& v, GWFrames::Modes& delta) const {
  /// \param uPrime New retarded time at which to give the data
  /// \param v Three-vector of the boost relative to the current frame
  /// \param delta Spherical-harmonic modes of the supertranslation
  /// 
  /// A general BMS transformation is expressed as a conformal
  /// transformation of the sphere (which encompasses rotations and
  /// boosts) and a supertranslation (which affects only the time
  /// coordinate on null infinity).  Here, the conformal
  /// transformation of the sphere is assumed to be a simple boost,
  /// described by the velocity vector `v`.  (The other freedom in the
  /// conformal group is a simple rotation, which we assume is zero.)
  /// The supertranslation is decomposed into (scalar) spherical
  /// harmonics.
  /// 
  /// The work to be done by this function includes (1) evaluating the
  /// data on the appropriate equi-angular grid of the final frame,
  /// while simultaneously transforming the data at each point by the
  /// appropriate spin factor, boost factor, and mixing; (2)
  /// interpolating those data to the appropriate values of retarded
  /// time of the final frame; and (3) transforming back to spectral
  /// space to store the data in their usual representation.
  /// 
  /// The relation between the new and old time coordinates is \f$u' =
  /// K(u-\delta)\f$, where \f$K\f$ is the conformal factor (which is
  /// a function of angle).  So we need to interpolate at each point
  /// to the original time coordinate \f$u = u'/K + \delta\f$, which
  /// again depends on angle.  That is, we have to interpolate to a
  /// different time for each grid point.
  
  const int n_theta = 2*slices[0].EllMax()+1;
  const int n_phi = n_theta;
  
  // (0) Find current time slices on which we need data to interpolate
  // to the new time slice
  const DataGrid u = uPrime/GWFrames::ConformalFactorGrid(v, n_theta, n_phi) + DataGrid(delta, n_theta, n_phi);
  double uMax = std::real(u[0]);
  double uMin = std::real(u[0]);
  for(int i=1; i<n_theta*n_phi; ++i) {
    const double u_i = std::real(u[i]);
    if(u_i>uMax) { uMax = u_i; }
    if(u_i<uMin) { uMin = u_i; }
  }
  if(uMin<t[0] || uMax>t.back()) {
    std::cerr << "\n\n" << __FILE__ << ":" << __LINE__
	      << "\nError: (uMin=" << uMin << ") < (t[0]=" << t[0] << ") or (uMax=" << uMax << ") > (t[-1]=" << t.back() << ")"
	      << "\n       Cannot extrapolate data.\n"
	      << std::endl;
    throw(GWFrames_ValueError);
  }
  int iMax = t.size()-1;
  while(t[iMax]>uMax && iMax>0) { --iMax; } // t[iMax] is now strictly less than uMax
  int iMin = 0;
  while(t[iMin]<uMin && iMin<int(t.size())-1) { ++iMin; } // t[iMin] is now strictly greater than uMin
  iMin = std::max(0, iMin-3);
  iMax = std::min(int(t.size())-1, std::max(iMin+7, iMax+3));
  
  // (1) Evaluate BMS-transformed data on equi-angular grids of the final frame at a series of times
  const unsigned int Nslices = iMax-iMin+1;
  vector<SliceGrid> transformedslices(Nslices);
  vector<double> u_original(Nslices);
  for(int i=iMin; i<=iMax; ++i) {
    u_original[i-iMin] = t[i];
    transformedslices[i-iMin] = slices[i].BMSTransformationOnSlice(t[i], v, delta);
  }
  
  // (2) Interpolate to new retarded time
  // Create new object to hold the data
  SliceGrid BMStransformedGrid(n_theta*n_phi);
  // Initialize the GSL interpolators for the data
  gsl_interp_accel* accRe = gsl_interp_accel_alloc();
  gsl_interp_accel* accIm = gsl_interp_accel_alloc();
  gsl_spline* splineRe = gsl_spline_alloc(gsl_interp_cspline, Nslices);
  gsl_spline* splineIm = gsl_spline_alloc(gsl_interp_cspline, Nslices);
  // Loop through, doing the work
  { int i_g=0;
    for(int i_t=0; i_t<n_theta; ++i_t) { // Loop over theta points
      for(int i_p=0; i_p<n_phi; ++i_p, ++i_g) { // Loop over phi points
	const double u_i = std::real(u[i_g]); // Interpolate the data at this point to u_i (measured in the current frame)
	for(int i_D=0; i_D<7; ++i_D) { // Loop over data types
	  // Fill the storage vectors for extrapolation
	  vector<double> re(Nslices);
	  vector<double> im(Nslices);
	  for(unsigned int i_s=0; i_s<Nslices; ++i_s) {
	    re[i_s] = std::real(transformedslices[i_s][i_D][i_g]);
	    im[i_s] = std::imag(transformedslices[i_s][i_D][i_g]);
	  }
	  // Initialize the interpolators for this data set
	  gsl_spline_init(splineRe, &(u_original)[0], &re[0], Nslices);
	  gsl_spline_init(splineIm, &(u_original)[0], &im[0], Nslices);
	  // Extrapolate real and imaginary parts and store data
	  BMStransformedGrid[i_D][i_g] = complex<double>( gsl_spline_eval(splineRe, u_i, accRe), gsl_spline_eval(splineIm, u_i, accIm) );
	}
      }
    }
  } // scope to kill i_g
  // Free the interpolators
  gsl_interp_accel_free(accRe);
  gsl_interp_accel_free(accIm);
  gsl_spline_free(splineRe);
  gsl_spline_free(splineIm);
  
  // (3) Transform back to spectral space
  SliceModes BMStransformed(slices[0].EllMax());
  for(int i_D=0; i_D<7; ++i_D) { // Loop over data types
    BMStransformed[i_D] = Modes(BMStransformedGrid[i_D]);
  }
  
  return BMStransformed;
}
