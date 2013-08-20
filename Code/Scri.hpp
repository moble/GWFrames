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
  
  // This header defines the basic objects `Modes` and `DataGrid`.
  // `SliceOfScri<D>` contains a set of either `Modes` or `DataGrid` objects.
  // `SliceModes` and `SliceGrid` are specific instances of this template; `SliceModes` has some extras
  // 
  
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
    DataGrid() : s(0), n_theta(0), n_phi(0), data(0) { }
    DataGrid(const DataGrid& A);
    DataGrid(const int Spin, const int N_theta, const int N_phi, const std::vector<std::complex<double> >& D);
    explicit DataGrid(Modes M, const int N_theta=0, const int N_phi=0); // Can't be const& because of spinsfast design
    DataGrid(const Modes& M, const GWFrames::ThreeVector& v, const int N_theta=0, const int N_phi=0);
  public: // Access and operators
    inline int Spin() const { return s; }
    inline int N_theta() const { return n_theta; }
    inline int N_phi() const { return n_phi; }
    inline const std::complex<double>& operator[](const unsigned int i) const { return data[i]; }
    inline std::complex<double>& operator[](const unsigned int i) { return data[i]; }
    inline std::vector<std::complex<double> > Data() const { return data; }
    DataGrid operator*(const DataGrid&) const;
    DataGrid operator/(const DataGrid&) const;
    DataGrid operator+(const DataGrid&) const;
    DataGrid operator-(const DataGrid&) const;
    DataGrid pow(const int p) const;
  }; // class DataGrid
  DataGrid operator*(const double& a, const DataGrid& b);
  DataGrid operator/(const double& a, const DataGrid& b);
  DataGrid operator-(const double& a, const DataGrid& b);
  DataGrid ConformalFactorGrid(const MobiusTransform& abcd, const int n_theta, const int n_phi);
  DataGrid ConformalFactorGrid(const ThreeVector& v, const int n_theta, const int n_phi);
  
  
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
    inline std::complex<double> operator[](const unsigned int i) const { return data[i]; }
    inline std::complex<double>& operator[](const unsigned int i) { return data[i]; }
    inline std::vector<std::complex<double> > Data() const { return data; }
  public: // Operations
    Modes bar() const;
    Modes operator*(const Modes& M) const;
    Modes operator/(const Modes& M) const;
    Modes operator+(const Modes& M) const;
    Modes edth() const;
    Modes edthbar() const;
    std::complex<double> EvaluateAtPoint(const double vartheta, const double varphi) const;
    std::complex<double> EvaluateAtPoint(const GWFrames::Quaternion& R) const;
  }; // class Modes
  
  
  template <class D>
  class SliceOfScri {
    /// This class holds all the necessary objects needed to
    /// understand the geometry of a given slice of null infinity.
  public: // Data
    double u; // retarded time of this slice
    D psi0, psi1, psi2, psi3, psi4, sigma, sigmadot; // complex mode data for these objects on this slice
  public: // Constructors
    SliceOfScri();
    SliceOfScri(const double& U,
  		const D& Psi0, const D& Psi1, const D& Psi2,
  		const D& Psi3, const D& Psi4, const D& Sigma, const D& SigmaDot);
  public: //Access
    inline const D& operator[](const unsigned int i) const {
      if(i==0) { return psi0; }
      else if(i==1) { return psi1; }
      else if(i==2) { return psi2; }
      else if(i==3) { return psi3; }
      else if(i==4) { return psi4; }
      else if(i==5) { return sigma; }
      else if(i==6) { return sigmadot; }
      else { std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << "\nError: (i=" << i << ")>6\n" << std::endl; throw(GWFrames_IndexOutOfBounds); }
    }
    inline D& operator[](const unsigned int i) {
      if(i==0) { return psi0; }
      else if(i==1) { return psi1; }
      else if(i==2) { return psi2; }
      else if(i==3) { return psi3; }
      else if(i==4) { return psi4; }
      else if(i==5) { return sigma; }
      else if(i==6) { return sigmadot; }
      else { std::cerr << "\n\n" << __FILE__ << ":" << __LINE__ << "\nError: (i=" << i << ")>6\n" << std::endl; throw(GWFrames_IndexOutOfBounds); }
    }
  }; // class SliceOfScri
  
  
  typedef SliceOfScri<DataGrid> SliceGrid;
  
  
  class SliceModes : public SliceOfScri<Modes> {
  public:
    // Useful quantities
    int EllMax() const;
    double Mass() const;
    GWFrames::FourVector FourMomentum() const;
    Modes SuperMomentum() const;
    // Transformations
    SliceGrid BMSTransformationOnSlice(const ThreeVector& v, const Modes& gamma) const;
  }; // class SliceModes
  
  
  typedef std::vector<DataGrid> SliceOfScriGrids;
  
  
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
    std::vector<SliceModes> slices;
  public: // Constructor
    Scri(const GWFrames::Waveform& psi0, const GWFrames::Waveform& psi1,
  	 const GWFrames::Waveform& psi2, const GWFrames::Waveform& psi3,
  	 const GWFrames::Waveform& psi4, const GWFrames::Waveform& sigma);
  public: // Member functions
    // Transformations
    SliceModes BMSTransformation(const double& uPrime, const ThreeVector& v, GWFrames::Modes& gamma) const;
    // SliceOfScri BMSTransformation(const DataGrid& u, const GWFrames::MobiusTransform& abcd, GWFrames::Modes& gamma, const int iMin, const iMax) const;
    // Scri BMSTransformation(const GWFrames::MobiusTransform& abcd, GWFrames::Modes& gamma) const;
  }; // class Scri
  
  
} // namespace GWFrames

#endif // SCRI_HPP
