// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

#ifndef SPHERICALHARMONICS_HPP
#define SPHERICALHARMONICS_HPP

#include "Utilities.hpp"
#include "Quaternions.hpp"

namespace GWFrames {
  
  const int ellMax_GWFrames = 16;
  
  /// Object for pre-computing and retrieving factorials
  class FactorialFunctor {
  private:
    static const std::vector<double> FactorialTable;
  public:
    FactorialFunctor() { };
    inline double operator()(const unsigned int i) const {
      #ifdef DEBUG
      if(i>171) {
  	std::cerr << "\n\ni = " << i << "\tiMax = 171"
  		  << "\nFactorialFunctor is only implemented up to 171!; larger values overflow."
  		  << std::endl;
  	throw(GWFrames_IndexOutOfBounds);
      }
      #endif
      return FactorialTable[i];
    }
  };
  
  /// Object for pre-computing and retrieving binomials
  class BinomialCoefficientFunctor {
  private:
    static const std::vector<double> BinomialCoefficientTable;
  public:
    BinomialCoefficientFunctor() { };
    inline double operator()(const unsigned int n, const unsigned int k) const {
      #ifdef DEBUG
      if(n>2*ellMax_Utilities || k>n) {
  	std::cerr << "\n\n(n, k) = (" << n << ", " << k << ")\t2*ellMax_Utilities = " << 2*ellMax_Utilities
  		  << "\nBinomialCoefficientFunctor is only implemented up to n=2*ellMax_Utilities=" << 2*ellMax_Utilities
  		  << ".\nTo increase this bound, edit 'ellMax_Utilities' in " << __FILE__ << " and recompile." << std::endl;
  	throw(GWFrames_IndexOutOfBounds);
      }
      #endif
      return BinomialCoefficientTable[(n*(n+1))/2+k];
    }
  };
  
  /// Object for pre-computing and retrieving values of the ladder operators
  class LadderOperatorFactorFunctor {
  private:
    static const std::vector<double> FactorTable;
  public:
    LadderOperatorFactorFunctor() { };
    inline double operator()(const int ell, const int m) const {
      #ifdef DEBUG
      if(ell>ellMax_Utilities || std::abs(m)>ell) {
  	std::cerr << "\n\n(ell, m) = (" << ell << ", " << m << ")\tellMax_Utilities = " << ellMax_Utilities
  		  << "\nLadderOperatorFactorFunctor is only implemented up to ell=" << ellMax_Utilities
  		  << ".\nTo increase this bound, edit 'ellMax_Utilities' in " << __FILE__ << " and recompile." << std::endl;
  	throw(GWFrames_IndexOutOfBounds);
      }
      #endif
      return FactorTable[ell*ell+ell+m];
    }
  };
  
  /// Object for pre-computing and retrieving coefficients for the Wigner D matrices
  class WignerCoefficientFunctor {
  private:
    static const std::vector<double> CoefficientTable;
  public:
    WignerCoefficientFunctor() { };
    inline double operator()(const int ell, const int mp, const int m) const {
      #ifdef DEBUG
      if(ell>ellMax_Utilities || std::abs(mp)>ell || std::abs(m)>ell) {
  	std::cerr << "\n\n(ell, mp, m) = (" << ell << ", " << mp << ", " << m << ")\tellMax_Utilities = " << ellMax_Utilities
  		  << "\nWignerCoefficientFunctor is only implemented up to ell=" << ellMax_Utilities
  		  << ".\nTo increase this bound, edit 'ellMax_Utilities' in " << __FILE__ << " and recompile." << std::endl;
  	throw(GWFrames_IndexOutOfBounds);
      }
      #endif
      return CoefficientTable[int(ell*(ell*(1.3333333333333333*ell + 2) + 1.6666666666666667) + mp*(2*ell + 1) + m + 0.5)];
    }
  };
  
  /// Object for computing the Wigner D matrices as functions of quaternion rotors
  class WignerDMatrix {
  private:
    BinomialCoefficientFunctor BinomialCoefficient;
    WignerCoefficientFunctor WignerCoefficient;
    std::complex<double> Ra, Rb;
    double absRa, absRb, absRRatioSquared;
  public:
    WignerDMatrix(const Quaternion& iR=Quaternion(1,0,0,0));
    WignerDMatrix& SetRotation(const Quaternion& iR);
    std::complex<double> operator()(const int ell, const int mp, const int m) const;
  };
  
  /// Object for computing values of the spin-weighted spherical harmonics (SWSHs)
  class SWSH {
    /// Note that this object is a functor taking a quaternion
    /// argument.  Ordinarily, we think of this quaternion as being
    /// the rotor taking the unit $z$ vector into a point
    /// \f$(\vartheta, \varphi)\f$, which gives us the usual form of
    /// SWSHs.  However, more general arguments are possible; no
    /// checking is done to ensure that the argument has the simple
    /// form of a minimal rotation to the spherical coordinate.
  private:
    WignerDMatrix D;
    int spin;
    double sign;
  public:
    // / \@cond
    SWSH(const int s, const Quaternion& iR=Quaternion(1,0,0,0))
      : D(iR), spin(s), sign(s%2==0 ? 1.0 : -1.0)
    { }
    // / \@endcond
    inline SWSH& SetRotation(const Quaternion& iR) { D.SetRotation(iR); return *this; }
    inline SWSH& SetAngles(const double vartheta, const double varphi) { D.SetRotation(Quaternion(vartheta, varphi)); return *this; }
    inline std::complex<double> operator()(const int ell, const int m) const {
      return sign * std::sqrt((2*ell+1)/(4*M_PI)) * D(ell, m, -spin);
    }
  };
  
} // namespace GWFrames

#endif // SPHERICALHARMONICS_HPP
