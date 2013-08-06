// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

#ifndef SCRI_HPP
#define SCRI_HPP

#include <vector>
#include <complex>
#include "Waveforms.hpp"

namespace GWFrames {
  
  class Modes; // Forward declaration for DataGrid constructor
  
  class DataGrid {
    friend class Modes;
  protected:
    int s;
    int n_theta;
    int n_phi;
    std::vector<std::complex<double> > data;
  public:
    DataGrid(const DataGrid& A);
    DataGrid(const int Spin, const int N_theta, const int N_phi, const std::vector<std::complex<double> >& D);
    explicit DataGrid(Modes M, const int N_theta=0, const int N_phi=0); // Can't be const& because of spinsfast design
    DataGrid operator*(const DataGrid&) const;
  }; // class DataGrid
  
  
  class Modes {
    friend class DataGrid;
  protected:
    int s;
    int ellMax;
    std::vector<std::complex<double> > data;
  public:
    Modes();
    Modes(const Modes& A);
    Modes(const int spin, const std::vector<std::complex<double> >& Data);
    explicit Modes(DataGrid D); // Can't be const& because of spinsfast design
    Modes bar() const;
    Modes operator*(const Modes& M) const;
    Modes operator+(const Modes& M) const;
    Modes edth(const int s) const;
    Modes edthbar(const int s) const;
  }; // class Modes
  
  
  class SliceOfScri {
  private:
    double u; // retarded time of this slice
    GWFrames::Modes psi0, psi1, psi2, psi3, psi4, sigma, sigmadot; // complex mode data for these objects on this slice
    
  public:
    // Constructors
    SliceOfScri();
    SliceOfScri(const double& U,
		const GWFrames::Modes& Psi0, const GWFrames::Modes& Psi1, const GWFrames::Modes& Psi2,
		const GWFrames::Modes& Psi3, const GWFrames::Modes& Psi4, const GWFrames::Modes& Sigma, const GWFrames::Modes& SigmaDot);
    
    // Transformations
    SliceOfScri ConformalTransformation(const GWFrames::MobiusTransform& abcd) const;
    
    // Useful quantities
    double Mass() const;
    GWFrames::FourVector FourMomentum() const;
    GWFrames::Modes SuperMomentum() const;
    
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
    
  public: // Constructors
    Scri(const Waveform& psi0, const Waveform& psi1, const Waveform& psi2, const Waveform& psi3, const Waveform& psi4, const Waveform& sigma);
    
  public: // Member functions
    // Transformations
    SliceOfScri BMSTransformation(const double& uPrime, const GWFrames::MobiusTransform& abcd, GWFrames::Modes& gamma) const;
    
  }; // class Scri
  
  
} // namespace GWFrames

#endif // SCRI_HPP
