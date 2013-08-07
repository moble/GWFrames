// Copyright (c) 2013, Michael Boyle
// See LICENSE file for details

#include "SphericalHarmonics.hpp"

using GWFrames::ellMax_GWFrames;
using GWFrames::FactorialFunctor;
using GWFrames::BinomialCoefficientFunctor;
using GWFrames::LadderOperatorFactorFunctor;
using GWFrames::WignerCoefficientFunctor;
using GWFrames::WignerDMatrix;
using GWFrames::StereographicCoordinate;
using GWFrames::Quaternion;
using std::vector;
using std::complex;
using std::cerr;
using std::endl;


#define SphericalHarmonics_Epsilon 1.0e-14

/// Class to create an object returning the factorial of an argument.
vector<double> FactorialTableCalculator() {
  /// Note that because a double is returned, only values up to 28!
  /// will be exact; higher values will be accurate to machine
  /// precision.  Values up to 170! only are allowed because higher
  /// values overflow.
  vector<double> FactorialTable(171);
  FactorialTable[0] = 1.0;
  for (int i=1;i<171;i++) {
    FactorialTable[i] = i*FactorialTable[i-1];
  }
  return FactorialTable;
}
const std::vector<double> FactorialFunctor::FactorialTable = FactorialTableCalculator();

vector<double> BinomialCoefficientCalculator() {
  /// We need (n+1) coefficients for each value of n from 0 (for
  /// completeness) up to 2*ellMax_GWFrames (hard coded in the header file).
  /// That's a total of
  ///   >>> from sympy import summation, symbols
  ///   >>> ellMax_GWFrames, n, k = symbols('ellMax_GWFrames n k', integer=True)
  ///   >>> summation(n+1, (n, 0, 2*ellMax_GWFrames))
  ///  2*ellMax_GWFrames**2 + 3*ellMax_GWFrames + 1
  /// With a similar calculation, we can see that the associated access
  /// operator needs element (n*(n+1)/2 + k) of the array.
  vector<double> BinomialCoefficientTable(2*ellMax_GWFrames*ellMax_GWFrames + 3*ellMax_GWFrames + 1);
  unsigned int i=0;
  FactorialFunctor Factorial;
  for(unsigned int n=0; n<=2*ellMax_GWFrames; ++n) {
    for(unsigned int k=0; k<=n; ++k) {
      BinomialCoefficientTable[i++] = std::floor(0.5+Factorial(n)/(Factorial(k)*Factorial(n-k)));
    }
  }
  return BinomialCoefficientTable;
}
const std::vector<double> BinomialCoefficientFunctor::BinomialCoefficientTable = BinomialCoefficientCalculator();

vector<double> LadderOperatorFactorCalculator() {
  /// We need (2*ell+1) coefficients for each value of ell from 0 (for
  /// completeness) up to ellMax_GWFrames (hard coded in the header file).
  /// That's a total of
  ///   >>> from sympy import summation, symbols
  ///   >>> ell, ellMax_GWFrames, m, mp = symbols('ell ellMax_GWFrames m mp', integer=True)
  ///   >>> summation(2*ell+1, (ell, 0, ellMax_GWFrames))
  ///   ellMax_GWFrames**2 + 2*ellMax_GWFrames + 1
  /// With a similar calculation, we can see that the associated access
  /// operator needs element 
  ///   >>> summation(2*ell+1, (ell, 0, ell-1)) + ell + m
  ///   ell**2 + ell + m
  std::vector<double> FactorTable(ellMax_GWFrames*ellMax_GWFrames + 2*ellMax_GWFrames + 1);
  unsigned int i=0;
  for(int ell=0; ell<=ellMax_GWFrames; ++ell) {
    for(int m=-ell; m<=ell; ++m) {
      FactorTable[i++] = std::sqrt(ell*(ell+1)-m*(m+1));
    }
  }
  return FactorTable;
}
const std::vector<double> LadderOperatorFactorFunctor::FactorTable = LadderOperatorFactorCalculator();

std::vector<double> WignerCoefficientCalculator() {
  /// We need (2*ell+1)*(2*ell+1) coefficients for each value of ell
  /// from 0 (for completenes) up to ellMax_GWFrames (hard coded in the header
  /// file).  That's a total of
  ///   >>> from sympy import summation, symbols, simplify
  ///   >>> from sympy.polys.polyfuncs import horner
  ///   >>> ell, ellMax_GWFrames, m, mp = symbols('ell ellMax_GWFrames m mp', integer=True)
  ///   >>> horner(simplify(summation((2*ell+1)**2, (ell, 0, ellMax_GWFrames))))
  ///   ellMax_GWFrames*(ellMax_GWFrames*(4*ellMax_GWFrames/3 + 4) + 11/3) + 1
  /// With a similar calculation, we can see that the associated access
  /// operator needs element
  ///   >>> horner(summation((2*ell+1)**2, (ell, 0, ell-1)) + (2*ell+1)*(ell+mp) + ell + m)
  ///   ell*(ell*(4*ell/3 + 2) + 5/3) + mp*(2*ell + 1) + m
  /// of the array.
  std::vector<double> CoefficientTable(int(ellMax_GWFrames*(ellMax_GWFrames*(1.3333333333333333*ellMax_GWFrames + 4) + 3.6666666666666667) + 1 + 0.5));
  FactorialFunctor Factorial;
  unsigned int i=0;
  for(int ell=0; ell<=ellMax_GWFrames; ++ell) {
    for(int mp=-ell; mp<=ell; ++mp) {
      for(int m=-ell; m<=ell; ++m) {
	CoefficientTable[i++] = 
	  std::sqrt( Factorial(ell+m)*Factorial(ell-m)
	  	     / double(Factorial(ell+mp)*Factorial(ell-mp)) );
      }
    }
  }
  return CoefficientTable;
}
const std::vector<double> WignerCoefficientFunctor::CoefficientTable = WignerCoefficientCalculator();


/// Construct the D matrix object given the (optional) rotor.
WignerDMatrix::WignerDMatrix(const Quaternion& R)
  : BinomialCoefficient(), WignerCoefficient(),
    Ra(R[0], R[3]), Rb(R[2], R[1]),
    absRa(abs(Ra)), absRb(abs(Rb)), absRRatioSquared(absRb*absRb/(absRa*absRa))
{ }

/// Reset the rotor for this object to the given value.
WignerDMatrix& WignerDMatrix::SetRotation(const Quaternion& R) {
  Ra = std::complex<double>(R[0], R[3]);
  Rb = std::complex<double>(R[2], R[1]);
  absRa = abs(Ra);
  absRb = abs(Rb);
  absRRatioSquared = absRb*absRb/(absRa*absRa);
  return *this;
}

/// Evaluate the D matrix element for the given (ell, mp, m) indices.
std::complex<double> WignerDMatrix::operator()(const int ell, const int mp, const int m) const {
  if(absRa < SphericalHarmonics_Epsilon) {
    return (mp!=-m ? 0.0 : ((ell+mp)%2==0 ? 1.0 : -1.0) * std::pow(Rb, 2*mp) );
  }
  if(absRb < SphericalHarmonics_Epsilon) {
    return (mp!=m ? 0.0 : std::pow(Ra, 2*mp) );
  }
  if(absRa < 1.e-3) { // Deal with NANs in certain cases
    const std::complex<double> Prefactor =
      WignerCoefficient(ell, mp, m) * std::pow(Ra, m+mp) * std::pow(Rb, m-mp);
    const int rhoMin = std::max(0,mp-m);
    const int rhoMax = std::min(ell+mp,ell-m);
    const double absRaSquared = absRa*absRa;
    const double absRbSquared = absRb*absRb;
    double Sum = 0.0;
    for(int rho=rhoMax; rho>=rhoMin; --rho) {
      const double aTerm = std::pow(absRaSquared, ell-m-rho);
      if(aTerm != aTerm || aTerm<1.e-100) { // This assumes --fast-math is off
	Sum *= absRbSquared;
	continue;
      }
      Sum = ( (rho%2==0 ? 1 : -1) * BinomialCoefficient(ell+mp,rho) * BinomialCoefficient(ell-mp, ell-rho-m) * aTerm )
	+ ( Sum * absRbSquared );
    }
    return Prefactor * Sum * std::pow(absRbSquared, rhoMin);
  }
  const std::complex<double> Prefactor =
    (WignerCoefficient(ell, mp, m) * std::pow(absRa, 2*ell-2*m))
    * std::pow(Ra, m+mp) * std::pow(Rb, m-mp);
  const int rhoMin = std::max(0,mp-m);
  const int rhoMax = std::min(ell+mp,ell-m);
  double Sum = 0.0;
  for(int rho=rhoMax; rho>=rhoMin; --rho) {
    Sum = ( (rho%2==0 ? 1 : -1) * BinomialCoefficient(ell+mp,rho) * BinomialCoefficient(ell-mp, ell-rho-m) )
      + ( Sum * absRRatioSquared );
  }
  return Prefactor * Sum * std::pow(absRRatioSquared, rhoMin);
}
