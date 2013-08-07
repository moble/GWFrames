// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

#ifndef SCRI_HPP
#define SCRI_HPP

#include <vector>
#include <complex>
#include "Utilities.hpp"
#include "Quaternions.hpp"
#include "Waveforms.hpp"

namespace GWFrames {
  
  class Modes; // Forward declaration for DataGrid constructor
  
  class DataGrid {
    /// This object holds complex spin-weighted data on the sphere in
    /// an equi-angular representation.  That is, given integers
    /// n_theta and n_phi, the data are recorded on a "rectangular"
    /// grid with n_theta*n_phi points, including n_phi points at each
    /// of the poles.  The purpose of these objects is primarily to
    /// serve as a computational tool for pointwise multiplication.
    /// Otherwise, `Modes` is expected to be the preferable
    /// representation.
  private: // Data
    int s;
    int n_theta;
    int n_phi;
    std::vector<std::complex<double> > data;
  public: // Constructors
    DataGrid(const DataGrid& A);
    DataGrid(const int Spin, const int N_theta, const int N_phi, const std::vector<std::complex<double> >& D);
    explicit DataGrid(Modes M, const int N_theta=0, const int N_phi=0); // Can't be const& because of spinsfast design
    DataGrid(const Modes& M, const GWFrames::ThreeVector& v, const int N_theta=0, const int N_phi=0);
  public: // Access and operators
    inline int Spin() const { return s; }
    inline int N_theta() const { return n_theta; }
    inline int N_phi() const { return n_phi; }
    inline std::complex<double> Data(const unsigned int i) const { return data[i]; }
    inline std::complex<double>& Data(const unsigned int i) { return data[i]; }
    inline std::vector<std::complex<double> > Data() const { return data; }
    DataGrid operator*(const DataGrid&) const;
  }; // class DataGrid
  
  
  class Modes {
    /// This object holds complex spin-weighted data on the sphere in
    /// a spin-weighted spherical-harmonic mode representation.  All
    /// modes are present, even for \f$\ell<|s|\f$.  The modes are
    /// also assumed to be stored in order, as \f$(\ell,m) = (0,0),
    /// (1,-1), (1,0), (1,1), (2,-2), \ldots\f$.
  private: // Data
    int s;
    int ellMax;
    std::vector<std::complex<double> > data;
  public: // Constructors
    Modes();
    Modes(const Modes& A);
    Modes(const int spin, const std::vector<std::complex<double> >& Data);
    explicit Modes(DataGrid D); // Can't be const& because of spinsfast design
  public: // Access
    inline int Spin() const { return s; }
    inline int EllMax() const { return ellMax; }
    inline std::complex<double> Data(const unsigned int i) const { return data[i]; }
    inline std::complex<double>& Data(const unsigned int i) { return data[i]; }
    inline std::vector<std::complex<double> > Data() const { return data; }
  public: // Operations
    Modes bar() const;
    Modes operator*(const Modes& M) const;
    Modes operator+(const Modes& M) const;
    Modes edth() const;
    Modes edthbar() const;
    std::complex<double> EvaluateAtPoint(const double vartheta, const double varphi) const;
    std::complex<double> EvaluateAtPoint(const GWFrames::Quaternion& R) const;
  }; // class Modes
  
  
  class SliceOfScri {
    /// This class holds all the necessary `Modes` objects needed to
    /// understand the geometry of a given slice of null infinity.
  private:
    double u; // retarded time of this slice
    Modes psi0, psi1, psi2, psi3, psi4, sigma, sigmadot; // complex mode data for these objects on this slice
  public:
    // Constructors
    SliceOfScri();
    SliceOfScri(const double& U,
  		const Modes& Psi0, const Modes& Psi1, const Modes& Psi2,
  		const Modes& Psi3, const Modes& Psi4, const Modes& Sigma, const Modes& SigmaDot);
    // Transformations
    // SliceOfScri ConformalTransformation(const GWFrames::MobiusTransform& abcd) const;
    // Useful quantities
    double Mass() const;
    GWFrames::FourVector FourMomentum() const;
    Modes SuperMomentum() const;
  }; // class SliceOfScri
  
  
  class Scri {
    /// A `Scri` object contains all the information needed to express
    /// the asymptotic geometry of null infinity, and symmetry
    /// transformations thereof.  This geometry is encoded in the
    /// Newman--Penrose curvature scalars \f$\psi_0, \ldots, \psi_4\f$
    /// and the complex shear \f$\sigma\f$ of outgoing null rays
    /// (strain in gravitational-wave detectors).  The general
    /// symmetry transformation is an element of the
    /// Bondi--Metzner--Sachs (BMS) group, which transforms the data
    /// contained by `Scri` among itself.
  private: // Member data
    std::vector<double> t;
    std::vector<SliceOfScri> slices;
  public: // Constructor
    Scri(const GWFrames::Waveform& psi0, const GWFrames::Waveform& psi1,
  	 const GWFrames::Waveform& psi2, const GWFrames::Waveform& psi3,
  	 const GWFrames::Waveform& psi4, const GWFrames::Waveform& sigma);
  public: // Member functions
    // Transformations
    SliceOfScri BMSTransformation(const double& uPrime, const GWFrames::MobiusTransform& abcd, GWFrames::Modes& gamma) const;
  }; // class Scri
  
  
} // namespace GWFrames

#endif // SCRI_HPP
