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

#include "Utilities.hpp"
#include "Quaternions.hpp"
#include "SphericalHarmonics.hpp"
#include "Waveforms.hpp"
#include "Errors.hpp"

using GWFrames::Quaternion;
using GWFrames::StereographicCoordinate;
using GWFrames::MobiusTransform;
using GWFrames::FourVector;
using GWFrames::DataGrid;
using GWFrames::Modes;
using GWFrames::SliceOfScri;
using GWFrames::Scri;

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::flush;
using std::endl;
using std::complex;



const std::complex<double> complexi(0.0,1.0);
const std::complex<double> zero(0.0,0.0);
const Quaternion zHat(0,0,0,1);


//////////////
// DataGrid //
//////////////

DataGrid::DataGrid(const DataGrid& A)
  : s(A.s), n_theta(A.n_theta), n_phi(A.n_phi), data(A.data)
{ }

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
  spinsfast_salm2map(reinterpret_cast<fftw_complex*>(&M.Data(0)),
		     reinterpret_cast<fftw_complex*>(&data[0]),
		     M.Spin(), n_theta, n_phi, M.EllMax());
}

DataGrid::DataGrid(const Modes& M, const ThreeVector& v, const int N_theta, const int N_phi)
  : s(M.Spin()), n_theta(std::max(N_theta, 2*M.EllMax()+1)), n_phi(std::max(N_phi, 2*M.EllMax()+1)), data(n_phi*n_theta, zero)
{
  const double dtheta = M_PI/double(n_theta-1); // theta should return to M_PI
  const double dphi = 2*M_PI/double(n_phi); // phi should not return to 2*M_PI
  int i=0;
  for(int i_theta=0; i_theta<n_theta; ++i_theta) {
    for(int i_phi=0; i_phi<n_phi; ++i_phi) {
      const Quaternion R(dtheta*i_theta, dphi*i_phi);
      const Quaternion Rb = GWFrames::Boost(v, (R*zHat*R.conjugate()).vec());
      data[i] = M.EvaluateAtPoint(Rb*R);
      ++i;
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
    C.data[i] = A.data[i] * B.data[i];
  }
  return C;
}



///////////
// Modes //
///////////

Modes::Modes()
  : s(0), ellMax(0), data()
{ }

Modes::Modes(const Modes& A)
  : s(A.s), ellMax(A.ellMax), data(A.data)
{ }

Modes::Modes(const int spin, const std::vector<std::complex<double> >& Data)
  : s(spin), ellMax(0), data(Data)
{
  // Find the appropriate ellMax for this data length
  for(; ellMax<=ellMax_GWFrames; ++ellMax) {
    if(N_lm(ellMax)==int(Data.size())) {
      break;
    }
  }
  
  // Make sure gamma doesn't have more ell modes than ellMax_GWFrames
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
  spinsfast_map2salm(reinterpret_cast<fftw_complex*>(&D.Data(0)),
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

/// Empty constructor
SliceOfScri::SliceOfScri()
  : u(), psi0(), psi1(), psi2(), psi3(), psi4(), sigma(), sigmadot()
{ }

/// Constructor from data
SliceOfScri::SliceOfScri(const double& U,
			 const Modes& Psi0, const Modes& Psi1, const Modes& Psi2,
			 const Modes& Psi3, const Modes& Psi4, const Modes& Sigma, const Modes& SigmaDot)
  : u(U), psi0(Psi0), psi1(Psi1), psi2(Psi2), psi3(Psi3), psi4(Psi4), sigma(Sigma), sigmadot(SigmaDot)
{ }

// /// Apply a conformal transformation to the data on the sphere
// SliceOfScri SliceOfScri::ConformalTransformation(const GWFrames::MobiusTransform& abcd) const {
//   throw(GWFrames_NotYetImplemented);
//   return SliceOfScri();
// }

/// Calculate the mass of the system from the four-momentum
double SliceOfScri::Mass() const {
  const FourVector p = FourMomentum();
  return std::sqrt(p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-p[3]*p[3]);
}

/// Calculate the four-momentum of the system from the supermomentum
GWFrames::FourVector SliceOfScri::FourMomentum() const {
  /// The (Bondi) four-momentum is given by the ell=0 and ell=1 modes
  /// of the supermomentum.
  const Modes Psi = SuperMomentum();
  const double sqrt3 = std::sqrt(3);
  const double sqrt6 = std::sqrt(6);
  FourVector p;
  p[0] = std::real(Psi.Data(0));
  p[1] = std::real((Psi.Data(1)-Psi.Data(3))/sqrt6);
  p[2] = std::real(complexi*(Psi.Data(1)+Psi.Data(3))/sqrt6);
  p[3] = std::real(Psi.Data(2)/sqrt3);
  return p;
}

Modes SliceOfScri::SuperMomentum() const {
  return psi2 + sigma*sigmadot.bar() + sigma.bar().edth().edth();
}



Scri::Scri(const GWFrames::Waveform& psi0, const GWFrames::Waveform& psi1,
	   const GWFrames::Waveform& psi2, const GWFrames::Waveform& psi3,
	   const GWFrames::Waveform& psi4, const GWFrames::Waveform& sigma)
  : t(psi0.T()), slices(t.size())
{
  throw(GWFrames_NotYetImplemented);
}

/// Apply a (constant) BMS transformation to data on null infinity
SliceOfScri Scri::BMSTransformation(const double& uPrime, const GWFrames::MobiusTransform& abcd, GWFrames::Modes& gamma) const {
  /// \param uPrime New retarded time at which to give the data
  /// \param abcd Mobius representation of the conformal transformation
  /// \param gamma Spherical-harmonic modes of the supertranslation
  /// 
  /// A general BMS transformation is expressed as a conformal
  /// transformation of the sphere (which encompasses rotations and
  /// boosts) and a supertranslation (which affects only the time
  /// coordinate on null infinity).  Here, the conformal
  /// transformation of the sphere is expressed in terms of the
  /// parameters of the Mobius transformation of the sphere's
  /// stereographic coordinates.  The supertranslation is decomposed
  /// into (scalar) spherical harmonics.
  /// 
  /// The work to be done by this function includes (1) evaluating the
  /// data on the appropriate equi-angular grid of the final frame;
  /// (2) interpolating those data to the appropriate values of
  /// retarded time of the final frame; (3) transforming the data at
  /// each point by the appropriate spin factor, boost factor, and
  /// mixing; and (4) transforming back to spectral space to store the
  /// data in their usual representation.
  /// 
  /// The relation between the new and old time coordinates is \f$u' =
  /// K(u-\gamma)\f$, where \f$K\f$ is the conformal factor (which is
  /// a function of angle).  So we need to interpolate at each point
  /// to the original time coordinate \f$u = u'/K + \gamma\f$, which
  /// again depends on angle.  That is, we have to interpolate to a
  /// different time for each grid point.
  
  throw(GWFrames_NotYetImplemented);
  
  // (1) Evaluate data on equi-angular grids of the final frame

  // (2) Interpolate to new retarded time
  
  // (3) Transform data at each point
  
  // (4) Transform back to spectral space
  
  
  return SliceOfScri();
}
